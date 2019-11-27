
#include "lusol.h"
#include "myblas.h"
#include <math.h>
#include "heap.h"
#include "matcl-blas-ext/lapack_ext/lapack_ext.h"
#include "matcl-core/utils/workspace.h"
//#include "matcl-matrep/matcl_matrep.h"

namespace lusol
{

/* ==================================================================
   lu1pq1  constructs a permutation  iperm  from the array  len.
   ------------------------------------------------------------------
   On entry:
   len(i)  holds the number of nonzeros in the i-th row (say)
           of an m by n matrix.
   num(*)  can be anything (workspace).

   On exit:
   iperm   contains a list of row numbers in the order
           rows of length 0,  rows of length 1,..., rows of length n.
   loc(nz) points to the first row containing  nz  nonzeros,
           nz = 1, n.
   inv(i)  points to the position of row i within iperm(*).
   ================================================================== */
template<class T>
void LUSOLrec<T>::LU1PQ1(INT M, INT N, INT LEN[],
            INT IPERM[], INT LOC[], INT INV[], INT NUM[])
{
    INT NZEROS, NZ, I, L;

    // Count the number of rows of each length.
    NZEROS = 0;
    for(NZ = 1; NZ <= N; NZ++) 
    {
        NUM[NZ] = 0;
        LOC[NZ] = 0;
    }
    for(I = 1; I <= M; I++) 
    {
        NZ = LEN[I];
        if(NZ==0)
            NZEROS++;
        else
            NUM[NZ]++;
    }
    // Set starting locations for each length.
    L = NZEROS+1;
    for(NZ = 1; NZ <= N; NZ++) 
    {
        LOC[NZ] = L;
        L += NUM[NZ];
        NUM[NZ] = 0;
    }
    // Form the list.
    NZEROS = 0;
    for(I = 1; I <= M; I++) 
    {
        NZ = LEN[I];
        if(NZ==0) 
        {
            NZEROS++;
            IPERM[NZEROS] = I;
        }
        else 
        {
            L = LOC[NZ]+NUM[NZ];
            IPERM[L] = I;
            NUM[NZ]++;
        }
    }
    // Define the inverse of iperm.
    for(L = 1; L <= M; L++) 
    {
        I = IPERM[L];
        INV[I] = L;
    }
};

/* ==================================================================
   lu1pq2 frees the space occupied by the pivot row,
   and updates the column permutation iq.
   Also used to free the pivot column and update the row perm ip.
   ------------------------------------------------------------------
   nzpiv   (input)    is the length of the pivot row (or column).
   nzchng  (output)   is the net change in total nonzeros.
   ------------------------------------------------------------------
   14 Apr 1989  First version.
   ================================================================== */
template<class T>
void LUSOLrec<T>::LU1PQ2(INT NZPIV, INT *NZCHNG,
            INT IND[], INT LENOLD[], INT LENNEW[], INT IXLOC[], INT IX[], INT IXINV[])
{
    INT LR, J, NZ, NZNEW, L, NEXT, LNEW, JNEW;

    *NZCHNG = 0;
    for(LR = 1; LR <= NZPIV; LR++) 
    {
        J = IND[LR];
        IND[LR] = 0;
        NZ = LENOLD[LR];
        NZNEW = LENNEW[J];
        if(NZ!=NZNEW) 
        {
            L = IXINV[J];
            *NZCHNG = (*NZCHNG+NZNEW)-NZ;
            // l above is the position of column j in iq  (so j = iq(l)).
            if(NZ<NZNEW) 
            {
                // Column  j  has to move towards the end of  iq.
              x110:
                NEXT = NZ+1;
                LNEW = IXLOC[NEXT]-1;
                if(LNEW!=L) 
                {
                    JNEW = IX[LNEW];
                    IX[L] = JNEW;
                    IXINV[JNEW] = L;
                }
                L = LNEW;
                IXLOC[NEXT] = LNEW;
                NZ = NEXT;
                if(NZ<NZNEW)
                    goto x110;
            }
            else 
            {
                // Column  j  has to move towards the front of  iq.
              x120:
                LNEW = IXLOC[NZ];
                if(LNEW!=L) 
                {
                    JNEW = IX[LNEW];
                    IX[L] = JNEW;
                    IXINV[JNEW] = L;
                }
                L = LNEW;
                IXLOC[NZ] = LNEW+1;
                NZ = NZ-1;
                if(NZ>NZNEW)
                    goto x120;
            }
            IX[LNEW] = J;
            IXINV[J] = LNEW;
        }
    }
};

/* ==================================================================
   lu1pq3  looks at the permutation  iperm(*)  and moves any entries
   to the end whose corresponding length  len(*)  is zero.
   ------------------------------------------------------------------
   09 Feb 1994: Added work array iw(*) to improve efficiency.
   ================================================================== */
template<class T>
void LUSOLrec<T>::LU1PQ3(INT MN, INT LEN[], INT IPERM[], INT IW[], INT *NRANK)
{
    INT NZEROS, K, I;

    *NRANK = 0;
    NZEROS = 0;
    for(K = 1; K <= MN; K++) 
    {
        I = IPERM[K];

        if(LEN[I] == 0) 
        {
            NZEROS++;
            IW[NZEROS] = I;
        }
        else 
        {
            (*NRANK)++;
            IPERM[*NRANK] = I;
        }
    }
    for(K = 1; K <= NZEROS; K++)
        IPERM[(*NRANK)+K] = IW[K];
}

/* ==================================================================
   lu1rec
   ------------------------------------------------------------------
   On exit:
   ltop         is the length of useful entries in ind(*), a(*).
   ind(ltop+1)  is "i" such that len(i), loc(i) belong to the last
                item in ind(*), a(*).
   ------------------------------------------------------------------
   00 Jun 1983: Original version of lu1rec followed John Reid's
                compression routine in LA05.  It recovered
                space in ind(*) and optionally a(*)
                by eliminating entries with ind(l) = 0.
                The elements of ind(*) could not be negative.
                If len(i) was positive, entry i contained
                that many elements, starting at  loc(i).
                Otherwise, entry i was eliminated.
   23 Mar 2001: Realised we could have len(i) = 0 in rare cases!
                (Mostly during TCP when the pivot row contains
                a column of length 1 that couldn't be a pivot.)
                Revised storage scheme to
                   keep        entries with       ind(l) >  0,
                   squeeze out entries with -n <= ind(l) <= 0,
                and to allow len(i) = 0.
                Empty items are moved to the end of the compressed
                ind(*) and/or a(*) arrays are given one empty space.
                Items with len(i) < 0 are still eliminated.
   27 Mar 2001: Decided to use only ind(l) > 0 and = 0 in lu1fad.
                Still have to keep entries with len(i) = 0.
   ================================================================== */
template<class T>
void LUSOLrec<T>::LU1REC(INT N, bool REALS, INT *LTOP,
                             INT IND[], INT LEN[], INT LOC[])
{
    INT  NEMPTY, I, LENI, L, K, KLAST, ILAST, LPRINT;

    NEMPTY = 0;
    for(I = 1; I <= N; I++) 
    {
        LENI = LEN[I];
        if(LENI>0) 
        {
            L = (LOC[I]+LENI)-1;
            LEN[I] = IND[L];
            IND[L] = -(N+I);
        }
        else if(LENI==0)
            NEMPTY++;
    }
    K = 0;
    // Previous k
    KLAST = 0;
    // Last entry moved.
    ILAST = 0;
    for(L = 1; L <= *LTOP; L++) 
    {
        I = IND[L];
        if(I>0) 
        {
            K++;
            IND[K] = I;
            if(REALS)
                this->a[K] = this->a[L];
        }
        else if(I<-N) 
        {
            // This is the end of entry  i.
            I = -(N+I);
            ILAST = I;
            K++;
            IND[K] = LEN[I];
            if(REALS)
                this->a[K] = this->a[L];
            LOC[I] = KLAST+1;
            LEN[I] = K-KLAST;
            KLAST = K;
        }
    }
    // Move any empty items to the end, adding 1 free entry for each.
    if(NEMPTY>0) 
    {
        for(I = 1; I <= N; I++) 
        {
            if(LEN[I]==0) 
            {
                K++;
                LOC[I] = K;
                IND[K] = 0;
                ILAST = I;
            }
        }
    }
    LPRINT = this->luparm[LUSOL_IP_PRINTLEVEL];
    if(LPRINT>=LUSOL_MSG_PIVOT)
    {
        //FIXME:
        LUSOL_report(this, 0, (char*) "lu1rec.  File compressed from %d to %d\n",
                        *LTOP,K,REALS,NEMPTY);
    };
    // ncp 
    this->luparm[LUSOL_IP_COMPRESSIONS_LU]++;
    // Return ilast in ind(ltop + 1).
    *LTOP = K;
    IND[(*LTOP)+1] = ILAST;
}

/* ==================================================================
   lu1slk  sets w(j) > 0 if column j is a unit vector.
   ------------------------------------------------------------------
   21 Nov 2000: First version.  lu1fad needs it for TCP.
                Note that w(*) is nominally an integer array,
                but the only spare space is the double array w(*).
   ================================================================== */
template<class T>
void LUSOLrec<T>::LU1SLK()
{
    INT J, LC1, LQ, LQ1, LQ2;

    for(J = 1; J <= this->n; J++) 
    {
        this->wr[J] = 0;
    }
    LQ1 = this->iqloc ? this->iqloc[1] : this->n+1;
    LQ2 = this->n;

    if(this->m>1)
        LQ2 = this->iqloc[2]-1;

    for(LQ = LQ1; LQ <= LQ2; LQ++) 
    {
        J = this->iq[LQ];
        LC1 = this->locc[J];
        if(abs(this->a[LC1])==1) 
        {
            this->wr[J] = 1;
        }
    }
};

/* ==================================================================
   lu1gau does most of the work for each step of Gaussian elimination.
   A multiple of the pivot column is added to each other column j
   in the pivot row.  The column list is fully updated.
   The row list is updated if there is room, but some fill-ins may
   remain, as indicated by ifill and jfill.
   ------------------------------------------------------------------
   Input:
      ilast    is the row    at the end of the row    list.
      jlast    is the column at the end of the column list.
      lfirst   is the first column to be processed.
      lu + 1   is the corresponding element of U in au(*).
      nfill    keeps track of pending fill-in.
      a(*)     contains the nonzeros for each column j.
      indc(*)  contains the row indices for each column j.
      al(*)    contains the new column of L.  A multiple of it is
               used to modify each column.
      mark(*)  has been set to -1, -2, -3, ... in the rows
               corresponding to nonzero 1, 2, 3, ... of the col of L.
      au(*)    contains the new row of U.  Each nonzero gives the
               required multiple of the column of L.

   Workspace:
      markl(*) marks the nonzeros of L actually used.
               (A different mark, namely j, is used for each column.)

   Output:
      ilast     New last row    in the row    list.
      jlast     New last column in the column list.
      lfirst    = 0 if all columns were completed,
                > 0 otherwise.
      lu        returns the position of the last nonzero of U
                actually used, in case we come back in again.
      nfill     keeps track of the total extra space needed in the
                row file.
      ifill(ll) counts pending fill-in for rows involved in the new
                column of L.
      jfill(lu) marks the first pending fill-in stored in columns
                involved in the new row of U.
   ------------------------------------------------------------------
   16 Apr 1989: First version of lu1gau.
   23 Apr 1989: lfirst, lu, nfill are now input and output
                to allow re-entry if elimination is interrupted.
   23 Mar 2001: Introduced ilast, jlast.
   27 Mar 2001: Allow fill-in "in situ" if there is already room
                up to but NOT INCLUDING the end of the
                row or column file.
                Seems safe way to avoid overwriting empty rows/cols
                at the end.  (May not be needed though, now that we
                have ilast and jlast.)
   ================================================================== */
template<class T>
void LUSOLrec<T>::LU1GAU(INT MELIM, INT NSPARE,
            TR SMALL, INT LPIVC1, INT LPIVC2, INT *LFIRST, INT LPIVR2,
            INT LFREE, INT MINFRE, INT ILAST, INT *JLAST, INT *LROW, INT *LCOL,
            INT *LU, INT *NFILL,
            INT MARK[],  T AL[], INT MARKL[], T AU[], INT IFILL[], INT JFILL[])
{
    bool ATEND;
    INT    LR, J, LENJ, NFREE, LC1, LC2, NDONE, NDROP, L, I, LL, K,
            LR1, LAST, LREP, L1, L2, LC, LENI;
    register T UJ;
    T       AIJ;

    for(LR = *LFIRST; LR <= LPIVR2; LR++) 
    {
        J = this->indr[LR];
        LENJ = this->lenc[J];
        NFREE = LFREE - *LCOL;
        if(NFREE<MINFRE)
            goto x900;

        // ---------------------------------------------------------------
        //   Inner loop to modify existing nonzeros in column  j.
        //   Loop 440 performs most of the arithmetic involved in the
        //   whole LU factorization.
        //   ndone counts how many multipliers were used.
        //   ndrop counts how many modified nonzeros are negligibly small.
        //   ---------------------------------------------------------------
        (*LU)++;
        UJ = AU[*LU];
        LC1 = this->locc[J];
        LC2 = (LC1+LENJ)-1;
        ATEND = (bool) (J==*JLAST);
        NDONE = 0;
        if(LENJ==0)
            goto x500;
        NDROP = 0;
        for(L = LC1; L <= LC2; L++) 
        {
            I = this->indc[L];
            LL = -MARK[I];
            if(LL>0) 
            {
                NDONE++;
                MARKL[LL] = J;
                this->a[L] += AL[LL]*UJ;
                if(abs(this->a[L])<=SMALL) 
                {
                    NDROP++;
                }
            }
        }
        //---------------------------------------------------------------
        //   Remove any negligible modified nonzeros from both
        //   the column file and the row file.
        //---------------------------------------------------------------
        if(NDROP==0)
            goto x500;
        K = LC1;
        for(L = LC1; L <= LC2; L++) 
        {
            I = this->indc[L];
            if(abs(this->a[L])<=SMALL)
                goto x460;
            this->a[K] = this->a[L];
            this->indc[K] = I;
            K++;
            continue;
            //Delete the nonzero from the row file.

          x460:
            LENJ--;
            this->lenr[I]--;
            LR1 = this->locr[I];
            LAST = LR1+this->lenr[I];
            for(LREP = LR1; LREP <= LAST; LREP++) 
            {
                if(this->indr[LREP]==J)
                    break;
            }
            this->indr[LREP] = this->indr[LAST];
            this->indr[LAST] = 0;
            if(I==ILAST)
                (*LROW)--;
        }
        // Free the deleted elements from the column file.

        memclear(this->indc+K, LC2-K+1);
        if(ATEND)
            *LCOL = K-1;
        //---------------------------------------------------------------
        //   Deal with the fill-in in column j.
        //---------------------------------------------------------------
      x500:
        if(NDONE==MELIM)
            goto x590;
        // See if column j already has room for the fill-in.
        if(ATEND)
            goto x540;
        LAST = (LC1+LENJ)-1;
        L1 = LAST+1;
        L2 = (LAST+MELIM)-NDONE;
        // 27 Mar 2001: Be sure it's not at or past end of the col file.
        if(L2>=*LCOL)
            goto x520;
        for(L = L1; L <= L2; L++) 
        {
            if(this->indc[L]!=0)
                goto x520;
        }
        goto x540;

      // We must move column j to the end of the column file.
      // First, leave some spare room at the end of the
      // current last column.
      x520:

        L1 = (*LCOL)+1;
        L2 = (*LCOL)+NSPARE;
        *LCOL = L2;
        for(L = L1; L <= L2; L++) 
        {
            // Spare space is free.
            this->indc[L] = 0;
        }
        ATEND = true;
        *JLAST = J;
        L1 = LC1;
        LC1 = (*LCOL)+1;
        this->locc[J] = LC1;
        for(L = L1; L <= LAST; L++) 
        {
            (*LCOL)++;
            this->a[*LCOL] = this->a[L];
            this->indc[*LCOL] = this->indc[L];
            // Free space.
            this->indc[L] = 0;
        }

      //---------------------------------------------------------------
      //     Inner loop for the fill-in in column j.
      //     This is usually not very expensive.
      //---------------------------------------------------------------
      x540:
        LAST = (LC1+LENJ)-1;
        LL = 0;
        for(LC = LPIVC1; LC <= LPIVC2; LC++) 
        {
            LL++;
            if(MARKL[LL]==J)
                continue;
            AIJ = AL[LL]*UJ;
            if(abs(AIJ)<=SMALL)
                continue;
            LENJ++;
            LAST++;
            this->a[LAST] = AIJ;
            I = this->indc[LC];
            this->indc[LAST] = I;
            LENI = this->lenr[I];

            //Add 1 fill-in to row i if there is already room.
            //  27 Mar 2001: Be sure it's not at or past the }
            //  of the row file.
            L = this->locr[I]+LENI;
            if(L>=*LROW)
                goto x550;
            if(this->indr[L]>0)
                goto x550;
            this->indr[L] = J;
            this->lenr[I] = LENI+1;
            continue;
          //Row i does not have room for the fill-in.
          //    Increment ifill(ll) to count how often this has
          //    happened to row i.  Also, add m to the row index
          //    indc(last) in column j to mark it as a fill-in that is
          //    still pending.
          //    If this is the first pending fill-in for row i,
          //    nfill includes the current length of row i
          //    (since the whole row has to be moved later).
          //    If this is the first pending fill-in for column j,
          //    jfill(lu) records the current length of column j
          //    (to shorten the search for pending fill-ins later). */
          x550:
            if(IFILL[LL]==0)
                (*NFILL) += LENI+NSPARE;
            if(JFILL[*LU]==0)
                JFILL[*LU] = LENJ;
            (*NFILL)++;
            IFILL[LL]++;
            this->indc[LAST] = this->m+I;
        }
        if(ATEND)
            *LCOL = LAST;
      // End loop for column  j.  Store its final length.
      x590:
        this->lenc[J] = LENJ;
    }
    // Successful completion.
    *LFIRST = 0;
    return;

  //Interruption.  We have to come back in after the
  //     column file is compressed.  Give lfirst a new value.
  //     lu and nfill will retain their current values. */
  x900:
    *LFIRST = LR;
};

/* ==================================================================
   lu1mar  uses a Markowitz criterion to select a pivot element
   for the next stage of a sparse LU factorization,
   subject to a Threshold Partial Pivoting stability criterion (TPP)
   that bounds the elements of L.
   ------------------------------------------------------------------
   gamma  is "gamma" in the tie-breaking rule TB4 in the LUSOL paper.
   ------------------------------------------------------------------
   Search cols of length nz = 1, then rows of length nz = 1,
   then   cols of length nz = 2, then rows of length nz = 2, etc.
   ------------------------------------------------------------------
   00 Jan 1986  Version documented in LUSOL paper:
                Gill, Murray, Saunders and Wright (1987),
                Maintaining LU factors of a general sparse matrix,
                Linear algebra and its applications 88/89, 239-270.
   02 Feb 1989  Following Suhl and Aittoniemi (1987), the largest
                element in each column is now kept at the start of
                the column, i.e. in position locc(j) of a and indc.
                This should speed up the Markowitz searches.
   26 Apr 1989  Both columns and rows searched during spars1 phase.
                Only columns searched during spars2 phase.
                maxtie replaced by maxcol and maxrow.
   05 Nov 1993  Initializing  "mbest = m * n"  wasn't big enough when
                m = 10, n = 3, and last column had 7 nonzeros.
   09 Feb 1994  Realised that "mbest = maxmn * maxmn" might overflow.
                Changed to    "mbest = maxmn * 1000".
   27 Apr 2000  On large example from Todd Munson,
                that allowed  "if (mbest .le. nz1**2) go to 900"
                to exit before any pivot had been found.
                Introduced kbest = mbest / nz1.
                Most pivots can be rejected with no integer multiply.
                true merit is evaluated only if it's as good as the
                best so far (or better).  There should be no danger
                of integer overflow unless A is incredibly
                large and dense.
   10 Sep 2000  TCP, aijtol added for Threshold Complete Pivoting.
   ================================================================== */
template<class T>
void LUSOLrec<T>::LU1MAR(INT MAXMN, bool TCP, TR AIJTOL, TR LTOL,
            INT MAXCOL, INT MAXROW, INT *IBEST, INT *JBEST, INT *MBEST)
{
    INT  KBEST, NCOL, NROW, NZ1, NZ, LQ1, LQ2, LQ, J, LC1, LC2, LC, I, LEN1, MERIT, LP1,
         LP2, LP, LR1, LR2, LR;
    TR ABEST, LBEST, AMAX, AIJ, CMAX;

    ABEST = 0;
    LBEST = 0;
    *IBEST = 0;
    *MBEST = -1;
    KBEST = MAXMN+1;
    NCOL = 0;
    NROW = 0;
    NZ1 = 0;
    for(NZ = 1; NZ <= MAXMN; NZ++) 
    {
        if(KBEST<=NZ1)
            goto x900;
        if(*IBEST>0) 
        {
            if(NCOL>=MAXCOL)
                goto x200;
        }
        if(NZ>this->m)
            goto x200;

        //---------------------------------------------------------------
        //   Search the set of columns of length  nz.
        //---------------------------------------------------------------
        LQ1 = this->iqloc[NZ];
        LQ2 = this->n;
        if(NZ<this->m)
            LQ2 = this->iqloc[NZ+1]-1;

        for(LQ = LQ1; LQ <= LQ2; LQ++) 
        {
            NCOL = NCOL+1;
            J = this->iq[LQ];
            LC1 = this->locc[J];
            LC2 = LC1+NZ1;
            AMAX = abs(this->a[LC1]);
            //Test all aijs in this column.
            //  amax is the largest element (the first in the column).
            //  cmax is the largest multiplier if aij becomes pivot.

            if(TCP) 
            {
                // Nothing in whole column
                if(AMAX<AIJTOL)
                    continue;
            }
            for(LC = LC1; LC <= LC2; LC++) 
            {
                I = this->indc[LC];
                LEN1 = this->lenr[I]-1;

                if(LEN1>KBEST)
                    continue;

                // aij  has a promising merit.
                // Apply the stability test.
                // We require  aij  to be sufficiently large compared to
                // all other nonzeros in column  j.  This is equivalent
                // to requiring cmax to be bounded by Ltol.

                if(LC==LC1) 
                {
                    // This is the maximum element, amax.
                    // Find the biggest element in the rest of the column
                    // and hence get cmax.  We know cmax .le. 1, but
                    // we still want it exactly in order to break ties.
                    // 27 Apr 2002: Settle for cmax = 1.
                    AIJ = AMAX;
                    CMAX = 1.;
                }
                else 
                {
                    // aij is not the biggest element, so cmax .ge. 1.
                    //Bail out if cmax will be too big.
                    AIJ = abs(this->a[LC]);
                    // Absolute test for Complete Pivoting
                    if(TCP) 
                    {
                        if(AIJ<AIJTOL)
                            continue;
                        // TPP
                    }
                    else 
                    {
                        if(AIJ*LTOL<AMAX)
                            continue;
                    }
                    CMAX = AMAX/AIJ;
                }
                // aij  is big enough.  Its maximum multiplier is cmax.
                MERIT = NZ1*LEN1;
                if(MERIT==*MBEST) 
                {
                    // Break ties.
                    // (Initializing mbest < 0 prevents getting here if
                    // nothing has been found yet.)
                    // In this version we minimize cmax
                    // but if it is already small we maximize the pivot.
                    if(LBEST<=this->parmlu[LUSOL_RP_GAMMA] && 
                            CMAX<=this->parmlu[LUSOL_RP_GAMMA]) 
                    {
                        if(ABEST>=AIJ)
                            continue;
                    }
                    else 
                    {
                        if(LBEST<=CMAX)
                            continue;
                    }
                }
                // aij  is the best pivot so far.
                *IBEST = I;
                *JBEST = J;
                KBEST = LEN1;
                *MBEST = MERIT;
                ABEST = AIJ;
                LBEST = CMAX;
                if(NZ==1)
                    goto x900;
            }
            // Finished with that column.
            if(*IBEST>0) 
            {
                if(NCOL>=MAXCOL)
                    goto x200;
            }
        }
      //---------------------------------------------------------------
      //     Search the set of rows of length  nz.
      //---------------------------------------------------------------
      x200:
        if(KBEST<=NZ)
            goto x900;
        if(*IBEST>0) 
        {
            if(NROW>=MAXROW)
                goto x290;
        }
        if(NZ>this->n)
            goto x290;
        LP1 = this->iploc[NZ];
        LP2 = this->m;
        if(NZ<this->n)
            LP2 = this->iploc[NZ+1]-1;
        for(LP = LP1; LP <= LP2; LP++) 
        {
            NROW++;
            I = this->ip[LP];
            LR1 = this->locr[I];
            LR2 = LR1+NZ1;
            for(LR = LR1; LR <= LR2; LR++) 
            {
                J = this->indr[LR];
                LEN1 = this->lenc[J]-1;
                if(LEN1>KBEST)
                    continue;
                // aij  has a promising merit.
                // Find where  aij  is in column  j.
                LC1 = this->locc[J];
                LC2 = LC1+LEN1;
                AMAX = abs(this->a[LC1]);
                for(LC = LC1; LC <= LC2; LC++) 
                {
                    if(this->indc[LC]==I)
                        break;
                }
                // Apply the same stability test as above.
                AIJ = abs(this->a[LC]);
                // Absolute test for Complete Pivoting
                if(TCP) 
                {
                    if(AIJ<AIJTOL)
                        continue;
                }
                if(LC==LC1) 
                {
                    // This is the maximum element, amax.
                    // Find the biggest element in the rest of the column
                    // and hence get cmax.  We know cmax .le. 1, but
                    // we still want it exactly in order to break ties.
                    // 27 Apr 2002: Settle for cmax = 1.
                    CMAX = 1.;
                }
                else 
                {
                    // aij is not the biggest element, so cmax .ge. 1.
                    // Bail out if cmax will be too big.
                    if(TCP) 
                    {
                        // relax
                    }
                    else 
                    {
                        if(AIJ*LTOL<AMAX)
                            continue;
                    }
                    CMAX = AMAX/AIJ;
                }
                // aij  is big enough.  Its maximum multiplier is cmax.
                MERIT = NZ1*LEN1;
                if(MERIT==*MBEST) 
                {
                    // Break ties as before.
                    // (Initializing mbest < 0 prevents getting here if
                    // nothing has been found yet.)
                    if(LBEST<=this->parmlu[LUSOL_RP_GAMMA] &&
                        CMAX<=this->parmlu[LUSOL_RP_GAMMA]) 
                    {
                        if(ABEST>=AIJ)
                            continue;
                    }
                    else 
                    {
                        if(LBEST<=CMAX)
                            continue;
                    }
                }
                // aij  is the best pivot so far.
                *IBEST = I;
                *JBEST = J;
                *MBEST = MERIT;
                KBEST = LEN1;
                ABEST = AIJ;
                LBEST = CMAX;
                if(NZ==1)
                    goto x900;
            }
            // Finished with that row.
            if(*IBEST>0) 
            {
                if(NROW>=MAXROW)
                    goto x290;
            }
        }
      // See if it's time to quit.
      x290:
        if(*IBEST>0) 
        {
            if(NROW>=MAXROW && NCOL>=MAXCOL)
                goto x900;
        }
        // Press on with next nz.
        NZ1 = NZ;
        if(*IBEST>0)
            KBEST = *MBEST/NZ1;
    }

  x900:
    ;
};


/* ==================================================================
   lu1mRP  uses a Markowitz criterion to select a pivot element
   for the next stage of a sparse LU factorization,
   subject to a Threshold Rook Pivoting stability criterion (TRP)
   that bounds the elements of L and U.
   ------------------------------------------------------------------
   11 Jun 2002: First version of lu1mRP derived from lu1mar.
   11 Jun 2002: Current version of lu1mRP.
   ================================================================== */
template<class T>
void LUSOLrec<T>::LU1MRP(INT MAXMN, TR LTOL, INT MAXCOL, INT MAXROW,
  INT *IBEST, INT *JBEST, INT *MBEST, TR AMAXR[])
{
    INT  I, J, KBEST, LC, LC1, LC2, LEN1, LP, LP1, LP2, LQ, LQ1,
         LQ2, LR, LR1, LR2, MERIT, NCOL, NROW, NZ, NZ1;
    TR ABEST, AIJ, AMAX, ATOLI, ATOLJ;

    //------------------------------------------------------------------
    //    Search cols of length nz = 1, then rows of length nz = 1,
    //    then   cols of length nz = 2, then rows of length nz = 2, etc.
    //------------------------------------------------------------------
    ABEST = 0;
    *IBEST = 0;
    KBEST = MAXMN+1;
    *MBEST = -1;
    NCOL = 0;
    NROW = 0;
    NZ1 = 0;
    for(NZ = 1; NZ <= MAXMN; NZ++) 
    {
        if(KBEST<=NZ1)
            goto x900;
        if(*IBEST>0) 
        {
            if(NCOL>=MAXCOL)
                goto x200;
        }
        if(NZ>this->m)
            goto x200;

        //---------------------------------------------------------------
        //  Search the set of columns of length  nz.
        //---------------------------------------------------------------
        LQ1 = this->iqloc[NZ];
        LQ2 = this->n;
        if(NZ<this->m)
            LQ2 = this->iqloc[NZ+1]-1;
        for(LQ = LQ1; LQ <= LQ2; LQ++) 
        {
            NCOL = NCOL+1;
            J = this->iq[LQ];
            LC1 = this->locc[J];
            LC2 = LC1+NZ1;
            AMAX = abs(this->a[LC1]);
            // Min size of pivots in col j
            ATOLJ = AMAX/LTOL;
            // Test all aijs in this column.
            for(LC = LC1; LC <= LC2; LC++) 
            {
                I = this->indc[LC];
                LEN1 = this->lenr[I]-1;
                if(LEN1>KBEST)
                    continue;

                // aij  has a promising merit.
                //  Apply the Threshold Rook Pivoting stability test.
                //  First we require aij to be sufficiently large
                //  compared to other nonzeros in column j.
                //  Then  we require aij to be sufficiently large
                //  compared to other nonzeros in row    i.
                AIJ = abs(this->a[LC]);
                if(AIJ<ATOLJ)
                    continue;
                if(AIJ*LTOL<AMAXR[I])
                    continue;
                // aij  is big enough.
                MERIT = NZ1*LEN1;
                if(MERIT==*MBEST) 
                {
                    // Break ties.
                    // (Initializing mbest < 0 prevents getting here if
                    // nothing has been found yet.)
                    if(ABEST>=AIJ)
                        continue;
                }
                // aij  is the best pivot so far.
                *IBEST = I;
                *JBEST = J;
                KBEST = LEN1;
                *MBEST = MERIT;
                ABEST = AIJ;
                if(NZ==1)
                    goto x900;
            }
            // Finished with that column.
            if(*IBEST>0) 
            {
                if(NCOL>=MAXCOL)
                    goto x200;
            }
        }
      //---------------------------------------------------------------
      //     Search the set of rows of length  nz.
      //---------------------------------------------------------------
      x200:
        if(KBEST<=NZ)
            goto x900;
        if(*IBEST>0) 
        {
            if(NROW>=MAXROW)
                goto x290;
        }
        if(NZ>this->n)
            goto x290;
        LP1 = this->iploc[NZ];
        LP2 = this->m;
        if(NZ<this->n)
            LP2 = this->iploc[NZ+1]-1;
        for(LP = LP1; LP <= LP2; LP++) 
        {
            NROW = NROW+1;
            I = this->ip[LP];
            LR1 = this->locr[I];
            LR2 = LR1+NZ1;
            // Min size of pivots in row i 
            ATOLI = AMAXR[I]/LTOL;
            for(LR = LR1; LR <= LR2; LR++) 
            {
                J = this->indr[LR];
                LEN1 = this->lenc[J]-1;
                if(LEN1>KBEST)
                    continue;

                // aij  has a promising merit.
                // Find where  aij  is in column j.
                LC1 = this->locc[J];
                LC2 = LC1+LEN1;
                AMAX = abs(this->a[LC1]);
                for(LC = LC1; LC <= LC2; LC++) 
                {
                    if(this->indc[LC]==I)
                        break;
                }
                // Apply the Threshold Rook Pivoting stability test.
                //  First we require aij to be sufficiently large
                //  compared to other nonzeros in row    i.
                //  Then  we require aij to be sufficiently large
                //  compared to other nonzeros in column j.
                AIJ = abs(this->a[LC]);
                if(AIJ<ATOLI)
                    continue;
                if(AIJ*LTOL<AMAX)
                    continue;
                // aij  is big enough.
                MERIT = NZ1*LEN1;
                if(MERIT==*MBEST) 
                {
                    // Break ties as before.
                    // (Initializing mbest < 0 prevents getting here if
                    // nothing has been found yet.)
                    if(ABEST>=AIJ)
                        continue;
                }
                // aij  is the best pivot so far.
                *IBEST = I;
                *JBEST = J;
                KBEST = LEN1;
                *MBEST = MERIT;
                ABEST = AIJ;
                if(NZ==1)
                    goto x900;
            }
            // Finished with that row.
            if(*IBEST>0) 
            {
                if(NROW>=MAXROW)
                    goto x290;
            }
        }
      // See if it's time to quit.
      x290:
        if(*IBEST>0) 
        {
            if(NROW>=MAXROW && NCOL>=MAXCOL)
                goto x900;
        }
        // Press on with next nz.
        NZ1 = NZ;
        if(*IBEST>0)
            KBEST = *MBEST/NZ1;
    }
  x900:
    ;
};

/* ==================================================================
   lu1mSP  is intended for symmetric matrices that are either
   definite or quasi-definite.
   lu1mSP  uses a Markowitz criterion to select a pivot element for
   the next stage of a sparse LU factorization of a symmetric matrix,
   subject to a Threshold Symmetric Pivoting stability criterion
   (TSP) restricted to diagonal elements to preserve symmetry.
   This bounds the elements of L and U and should have rank-revealing
   properties analogous to Threshold Rook Pivoting for unsymmetric
   matrices.
   ------------------------------------------------------------------
   14 Dec 2002: First version of lu1mSP derived from lu1mRP.
                There is no safeguard to ensure that A is symmetric.
   14 Dec 2002: Current version of lu1mSP.
   2013       : Accept only positive pivots. In case of positive 
                semidefinite matrices diagonal pivot always exits, if
                we cannot find one, then matrix is cleared.
   ================================================================== */
template<class T>
void LUSOLrec<T>::LU1MSP(INT MAXMN, TR LTOL, INT MAXCOL,
            INT *IBEST, INT *JBEST, INT *MBEST)
{
    INT  I, J, KBEST, LC, LC1, LC2, LQ, LQ1, LQ2, MERIT, NCOL, NZ, NZ1;
    TR ABEST, AIJ, AMAX, ATOLJ;

    //------------------------------------------------------------------
    //    Search cols of length nz = 1, then cols of length nz = 2, etc.
    //------------------------------------------------------------------
    ABEST = 0;
    *IBEST = 0;
    *MBEST = -1;
    KBEST = MAXMN+1;
    NCOL = 0;
    NZ1 = 0;
    for(NZ = 1; NZ <= MAXMN; NZ++) 
    {
        if(KBEST<=NZ1)
            goto x900;

        if(*IBEST>0) 
        {
            if(NCOL>=MAXCOL)
                goto x200;
        }
        if(NZ>this->m)
            goto x200;

        //---------------------------------------------------------------
        //   Search the set of columns of length  nz.
        //---------------------------------------------------------------
        LQ1 = this->iqloc[NZ];
        LQ2 = this->n;
        
        if(NZ<this->m)
            LQ2 = this->iqloc[NZ+1]-1;
        
        for(LQ = LQ1; LQ <= LQ2; LQ++) 
        {
            NCOL++;
            J = this->iq[LQ];
            LC1 = this->locc[J];
            LC2 = LC1+NZ1;
            AMAX = abs(this->a[LC1]);
            // Min size of pivots in col j
            ATOLJ = AMAX/LTOL;
            // Test all aijs in this column.
            // Ignore everything except the diagonal.
            for(LC = LC1; LC <= LC2; LC++) 
            {
                I = this->indc[LC];
                // Skip off-diagonals.
                if(I!=J)
                    continue;

                if(NZ1>KBEST)
                    continue;

                // aij  has a promising merit.
                //  Apply the Threshold Partial Pivoting stability test
                //  (which is equivalent to Threshold Rook Pivoting for
                //  symmetric matrices).
                //  We require aij to be sufficiently large
                //  compared to other nonzeros in column j.
                AIJ = matcl::lapack::real(this->a[LC]);
                if(AIJ<ATOLJ)
                    continue;
                // aij  is big enough.
                MERIT = NZ1*NZ1;
                if(MERIT==*MBEST) 
                {
                    // Break ties.
                    // (Initializing mbest < 0 prevents getting here if
                    // nothing has been found yet.)
                    if(ABEST>=AIJ)
                        continue;
                }
                // aij  is the best pivot so far.
                *IBEST = I;
                *JBEST = J;
                KBEST = NZ1;
                *MBEST = MERIT;
                ABEST = AIJ;
                if(NZ==1)
                    goto x900;
            }
            // Finished with that column.
            if(*IBEST>0) 
            {
                if(NCOL>=MAXCOL)
                    goto x200;
            }
        }
      // See if it's time to quit.
     x200:
        if(*IBEST>0) 
        {
            if(NCOL>=MAXCOL)
                goto x900;
        }
        // Press on with next nz.
        NZ1 = NZ;
        if(*IBEST>0)
            KBEST = *MBEST/NZ1;
    }
  x900:
    ;

    if (*IBEST == 0)
    {

        NZ1 = 0;
        for(NZ = 1; NZ <= MAXMN; NZ++) 
        {
            LQ1 = this->iqloc[NZ];
            LQ2 = this->n;
        
            if(NZ < this->m)
                LQ2 = this->iqloc[NZ+1]-1;
        
            for(LQ = LQ1; LQ <= LQ2; LQ++) 
            {
                J = this->iq[LQ];
                LC1 = this->locc[J];
                LC2 = LC1+NZ1;

                for(LC = LC1; LC <= LC2; LC++) 
                {
                    I = this->indc[LC];

                    this->a[LC] = 0.;
                    --this->lenr[I];
                    --this->lenc[J];
                }
            }
            NZ1 = NZ;
        };
    };
};

/* ==================================================================
   lu1mxc  moves the largest element in each of columns iq(k1:k2)
   to the top of its column.
   If k1 > k2, nothing happens.
   ------------------------------------------------------------------
   06 May 2002: (and earlier)
                All columns k1:k2 must have one or more elements.
   07 May 2002: Allow for empty columns.  The heap routines need to
                find 0.0 as the "largest element".
   29 Nov 2005: Bug fix - avoiding overwriting the next column when
                the current column is empty (i.e. LENJ==0)
                Yin Zhang <yzhang@cs.utexas.edu>
   ================================================================== */
template<class T>
void LUSOLrec<T>::LU1MXC(INT K1, INT K2, INT IX[])
{
    INT  I, J, K, L, LC, LENJ;
    T    AMAX;

    for(K = K1; K <= K2; K++) 
    {
        J = IX[K];
        LC = this->locc[J];
        LENJ = this->lenc[J];
        if(LENJ==0)  
        {}
        else 
        {
            L = amax(this->lenc[J], this->a + LC - LUSOL_ARRAYOFFSET) + LC - 1;
            if(L>LC) 
            {
                AMAX = this->a[L];
                this->a[L] = this->a[LC];
                this->a[LC] = AMAX;
                I = this->indc[L];
                this->indc[L] = this->indc[LC];
                this->indc[LC] = I;
            }
        }
    }
};

/* ==================================================================
   lu1mxr  finds the largest element in each of row ip(k1:k2)
   and stores it in Amaxr(*).  The nonzeros are stored column-wise
   in (a,indc,lenc,locc) and their structure is row-wise
   in (  indr,lenr,locr).
   ------------------------------------------------------------------
   11 Jun 2002: First version of lu1mxr.
                Allow for empty columns.
   10 Jan 2010: First f90 version.
   12 Dec 2011: Declare intent.
   03 Apr 2013: Recoded to improve efficiency.  Need new arrays
                markc(n), markr(m) and local array cols(n).
  
                First call:  mark = 0, k1 = 1, k2 = m.
                Initialize all of markc(n), markr(m), Amaxr(m).
                Columns are searched only once.
                cols(n) is not used.
  
                Later: mark := mark + 1 (greater than for previous call).
                Cols involved in rows p(k1:k2) are searched only once.
                cols(n) is local storage.
                markc(:), markr(:) are marked (= mark) in some places.
                For next call with new mark,
   06 Jun 2013: Reverted to f77 for lusol.f, trusting that the
                automatic array cols(n) is fine for f90 compilers.
   07 Jul 2013: Forgot that f2c won't like the automatic array.
                Now cols, markc, markr are passed in from
                lu1fac to lu1fad to here.

  ================================================================== */
template<class T>
void LUSOLrec<T>::LU1MXR(INT MARK, INT K1, INT K2, INT IX[], INT COLS[], INT MARKC[], 
                         INT MARKR[], TR AMAXR[])
{
    if (MARK == 0)
    {
        //First call: Find Amaxr(1:m) for original A.

        for(INT i = 1; i <= m; ++i)
        {
            MARKR[i]    = 0;
            AMAXR[i]    = 0.;
        };

        for(INT j = 1; j <= n; ++j)
        {
            MARKC[j]    = 0;
        };

        for (INT j = 1; j <= n; ++j)
        {
            INT LC1     = this->locc[j];
            INT LC2     = LC1 + this->lenc[j];

            for (INT LC = LC1; LC < LC2; ++LC)
            {
                INT i       = this->indc[LC];
                AMAXR[i]    = maximum<TR>(AMAXR[i], abs(this->a[LC]));
            };
        };
    }
    else
    {
        // Later calls: Find Amaxr(i) for rows i = p(k1:k2).

        INT ncol = 0;

        for (INT k = K1; k <= K2; ++k)
        {
            // Search rows to find which cols are involved.
            INT i           = IX[k];

            // Mark this row
            MARKR[i]        = MARK;
            AMAXR[i]        = 0;

            INT LR1         = this->locr[i];
            INT LR2         = LR1 + this->lenr[i] - 1;

            for (INT lr = LR1; lr <= LR2; ++lr)
            {
                // Mark all unmarked cols in this row.
                INT j        = this->indr[lr];
                
                // Build up a list of which ones they are.
                if (MARKC[j] != MARK)
                {
                    MARKC[j]    = MARK;
                    ncol        = ncol + 1;
                    COLS[ncol]  = j;
                };
            };
        };

        for (INT k = 1; k <= ncol; ++k)
        {         
            //Search involved columns.
            INT j               = COLS[k];
            INT lc1             = this->locc[j];
            INT lc2             = lc1 + this->lenc[j] - 1;

            for (INT lc = lc1; lc <= lc2; ++lc)
            {
                INT  i          = this->indc[lc];
                if (MARKR[i] == MARK)
                {
                    AMAXR[i]    =  maximum<TR>( AMAXR[i], abs(this->a[lc]) );
                };
            };
        };
    };
};

template<class V>
struct lusol_real_type
{
    using type = V;
};
template<class V>
struct lusol_real_type<std::complex<V>>
{
    using type = V;
};

/* ==================================================================
   lu1ful computes a dense (full) LU factorization of the
   mleft by nleft matrix that remains to be factored at the
   beginning of the nrowu-th pass through the main loop of lu1fad.
   ------------------------------------------------------------------
   02 May 1989: First version.
   05 Feb 1994: Column interchanges added to lu1DPP.
   08 Feb 1994: ipinv reconstructed, since lu1pq3 may alter ip.
   ================================================================== */
template<class T>
bool LUSOLrec<T>::LU1FUL(INT LEND, INT LU1, INT PIV,
            INT MLEFT, INT NLEFT, INT NRANK, INT NROWU,
            INT *LENL, INT *LENU, INT *NSING,
            bool KEEPLU, TR SMALL, T D[], INT IPVT[])
{
    INT  L, I, J, IPBASE, LDBASE, LQ, LC1, LC2, LC, LKK, LKN, LU, 
         IBEST, JBEST, LA, LL, NROWD, NCOLD;
    T    AI;
    T    AJ;

    using RT    = typename lusol_real_type<T>::type;

    //------------------------------------------------------------------
    //    If lu1pq3 moved any empty rows, reset ipinv = inverse of ip.
    //------------------------------------------------------------------
    if(NRANK<this->m) 
    {
        for(L = 1; L <= this->m; L++) 
        {
            I = this->ip[L];
            this->ipinv[I] = L;
        }

        if (PIV == LUSOL_PIVOT_TSP)
        {
            for(L = 1; L <= this->n; L++) 
            {
                I = this->iq[L];
                this->iqinv[I] = L;
            }
        }
    };
    //------------------------------------------------------------------
    //    Copy the remaining matrix into the dense matrix D.
    //------------------------------------------------------------------
    memclear((D+1), LEND);

    IPBASE = NROWU-1;
    LDBASE = 1-NROWU;
    INT LD;

    for(LQ = NROWU; LQ <= this->n; LQ++) 
    {
        J = this->iq[LQ];
        LC1 = this->locc[J];
        LC2 = (LC1+this->lenc[J])-1;
        for(LC = LC1; LC <= LC2; LC++) 
        {
            I = this->indc[LC];
            if (PIV == LUSOL_PIVOT_TSP)
            {
                LD = LDBASE+this->iqinv[I];
            }
            else
            {
                LD = LDBASE+this->ipinv[I];
            };
            D[LD] = this->a[LC];
        }
        LDBASE += MLEFT;
    };
    //------------------------------------------------------------------
    //    Call our favorite dense LU factorizer.
    //------------------------------------------------------------------
    {
        INT K = (MLEFT<NLEFT)?MLEFT:NLEFT;
        T* DA = D + 1;
        if(PIV == LUSOL_PIVOT_TPP)
        {
            matcl::lapack::i_type INFO = 0;
            matcl::lapack::getrf_rec(MLEFT, NLEFT, DA, MLEFT, IPVT+1, &INFO);
            *NSING = 0;
        }
        else if(PIV == LUSOL_PIVOT_TRP)
        {
            TR TOL = this->parmlu[LUSOL_RP_FACTORMAX_Lij];
            INT INFO = 0;

            T WORK_QUERY;
            matcl::lapack::getrfr(MLEFT, NLEFT, DA, MLEFT, IPVT+1, this->iq+NROWU, RT(1./TOL), RT(1./TOL), 
                                  RT(SMALL), &WORK_QUERY, -1, &INFO);

            INT LWORK       = (INT)matcl::lapack::real(WORK_QUERY);

            using VTR       = matcl::pod_type<T>;
            using workspace = matcl::pod_workspace<VTR>;
            workspace WORK  = workspace(LWORK);
            T* ptr_WORK     = reinterpret_cast<T*>(WORK.ptr());

            matcl::lapack::getrfr(MLEFT, NLEFT, DA, MLEFT, IPVT+1, this->iq+NROWU, RT(1./TOL), RT(1./TOL), 
                                  RT(SMALL), ptr_WORK, LWORK, &INFO);

            *NSING = (INT) K-INFO;
        }
        else if (PIV == LUSOL_PIVOT_TSP)
        {
            matcl::lapack::i_type rank;
            matcl::lapack::i_type INFO;
            RT* WORK = (RT*)malloc((2*MLEFT+1)*sizeof(RT));
            matcl::lapack::potfp3("upper", MLEFT, DA, MLEFT, IPVT + 1, rank, RT(SMALL), WORK, INFO);
            *NSING = K-rank;

            //permutations to row interchanges
            matcl::lapack::perm2int(MLEFT,IPVT + 1, reinterpret_cast<matcl::lapack::i_type*>(WORK), false);
            
            free(WORK);

            //change cholesky factor to LU factors
            int LD2 = MLEFT;
            for (int i = 0; i < rank; ++i)
            {
                T diag          = DA[i + i*LD2];
                RT scal         = matcl::lapack::real(diag);
                RT iscal        = RT(1.) / scal;

                DA[i + i*LD2]    = scal * DA[i + i*LD2];

                for (int j = i+1; j < MLEFT; ++j)
                {
                    DA[j + i*LD2]    = iscal * matcl::lapack::conj(DA[i + j*LD2]);
                    DA[i + j*LD2]    = scal * DA[i + j*LD2];
                };
            };

            //apply permutations to q vector
            for(int i = 0; i < MLEFT; ++i) 
            {
                int L1 = NROWU+i;
                int L2 = NROWU+IPVT[1+i] - 1;
                if(L1 != L2) 
                {
                    int I2 = this->iq[L1];
                    this->iq[L1] = this->iq[L2];
                    this->iq[L2] = I2;
                }
            };
        }
        else
        {
            matcl::lapack::i_type INFO = 0;
            matcl::lapack::getrfc(MLEFT, NLEFT, DA, MLEFT, IPVT+1, this->iq+NROWU, RT(SMALL), &INFO);            
            *NSING = K-INFO;
        }
        //L factor is stored as I - L
        for (INT I2 = 0; I2 < K; ++I2, DA += MLEFT)
        {
            for (INT J2 = I2+1; J2 < MLEFT; ++J2)
            {
                DA[J2] = -DA[J2];
            };
        };
    }

    //------------------------------------------------------------------
    //    Move D to the beginning of A,
    //    and pack L and U at the top of a, indc, indr.
    //    In the process, apply the row permutation to ip.
    //    lkk points to the diagonal of U.
    //------------------------------------------------------------------
    memcopy(this->a+1,D+1,LEND);
    #ifdef ClassicdiagU
        this->diagU = this->a + (this->lena-this->n);
    #endif
    LKK = 1;
    LKN = (LEND-MLEFT)+1;
    LU = LU1;
    
    INT* ip2 = (INT*)malloc((MLEFT+1)*sizeof(INT));
    if (ip2 == nullptr)
    {
        return true;
    };
    memcpy(ip2,this->ip+IPBASE,(MLEFT+1)*sizeof(INT));

    INT K, L1, L2;

    for(K = 1; K <= minimum(MLEFT,NLEFT); K++) 
    {
        L1 = K;
        L2 = IPVT[K];
        if(L1!=L2) 
        {
            I = ip2[L1];
            ip2[L1] = ip2[L2];
            ip2[L2] = I;
        }
    };

    for(K = 1; K <= minimum(MLEFT,NLEFT); K++) 
    {
        L1 = IPBASE+K;
        L2 = IPBASE+IPVT[K];
        if(L1!=L2) 
        {
            I = this->ip[L1];
            this->ip[L1] = this->ip[L2];
            this->ip[L2] = I;
        }
        IBEST = this->ip[L1];
        JBEST = this->iq[L1];
        if(KEEPLU) 
        {
            //===========================================================
            //  Pack the next column of L.
            //===========================================================
            LA = LKK;
            LL = LU;
            NROWD = 1;
            for(I = K+1; I <= MLEFT; I++) 
            {
                LA++;
                AI = this->a[LA];
                if(abs(AI)>SMALL) 
                {
                    NROWD = NROWD+1;
                    LL--;
                    this->a[LL] = AI;
                    this->indc[LL] = ip2[I];
                    this->indr[LL] = IBEST;
                }
            }
            //===========================================================
            //  Pack the next row of U.
            //  We go backwards through the row of D
            //  so the diagonal ends up at the front of the row of  U.
            //  Beware -- the diagonal may be zero.
            //===========================================================
            LA = LKN+MLEFT;
            LU = LL;
            NCOLD = 0;
            for(J = NLEFT; J >= K; J--) 
            {
                LA = LA-MLEFT;
                AJ = this->a[LA];
                if(abs(AJ)>SMALL || J==K) 
                {
                    NCOLD++;
                    LU--;
                    this->a[LU] = AJ;
                    this->indr[LU] = this->iq[IPBASE+J];
                }
            }
            this->lenr[IBEST] = -NCOLD;
            this->lenc[JBEST] = -NROWD;
            *LENL = ((*LENL)+NROWD)-1;
            *LENU = (*LENU)+NCOLD;
            LKN++;
        }
        else 
        {
            //===========================================================
            //  Store just the diagonal of U, in natural order.
            //===========================================================
            this->diagU[JBEST] = this->a[LKK];
        }
        LKK += MLEFT+1;
    }
    free(ip2);

    return false;
};


/* ==================================================================
   lu1or1  organizes the elements of an  m by n  matrix  A  as
   follows.  On entry, the parallel arrays   a, indc, indr,
   contain  nelem  entries of the form     aij,    i,    j,
   in any order.  nelem  must be positive.
   Entries not larger than the input parameter  small  are treated as
   zero and removed from   a, indc, indr.  The remaining entries are
   defined to be nonzero.  numnz  returns the number of such nonzeros
   and  Amax  returns the magnitude of the largest nonzero.
   The arrays  lenc, lenr  return the number of nonzeros in each
   column and row of  A.
   inform = 0  on exit, except  inform = 1  if any of the indices in
   indc, indr  imply that the element  aij  lies outside the  m by n
   dimensions of  A.
   ------------------------------------------------------------------
   xx Feb 1985: Original version.
   17 Oct 2000: a, indc, indr now have size lena to allow nelem = 0.
   ================================================================== */
template<class T>
void LUSOLrec<T>::LU1OR1(TR SMALL, TR *AMAX, INT *NUMNZ, INT *LERR, INT *INFORM)
{
    INT I, J, L, LDUMMY;

    memclear((this->lenr+1), this->m);
    memclear((this->lenc+1), this->n);

    *AMAX = 0;
    *NUMNZ = this->nelem;
    L = this->nelem+1;
    for(LDUMMY = 1; LDUMMY <= this->nelem; LDUMMY++) 
    {
        L--;
        if(abs(this->a[L])>SMALL) 
        {
            I = this->indc[L];
            J = this->indr[L];
            *AMAX = maximum<TR>(*AMAX,abs(this->a[L]));
            if(I<1 || I>this->m)
                goto x910;
            if(J<1 || J>this->n)
                goto x910;
            this->lenr[I]++;
            this->lenc[J]++;
        }
        else 
        {
            //  Replace a negligible element by last element.  Since
            //  we are going backwards, we know the last element is ok.
            this->a[L] = this->a[*NUMNZ];
            this->indc[L] = this->indc[*NUMNZ];
            this->indr[L] = this->indr[*NUMNZ];
            (*NUMNZ)--;
        }
    }
    *LERR = 0;
    *INFORM = LUSOL_INFORM_LUSUCCESS;
    return;

  x910:
    *LERR = L;
    *INFORM = LUSOL_INFORM_LUSINGULAR;
};

/* ==================================================================
   lu1or2  sorts a list of matrix elements  a(i,j)  into column
   order, given  numa  entries  a(i,j),  i,  j  in the parallel
   arrays  a, inum, jnum  respectively.  The matrix is assumed
   to have  n  columns and an arbitrary number of rows.
   On entry,  len(*)  must contain the length of each column.
   On exit,  a(*) and inum(*)  are sorted,  jnum(*) = 0,  and
   loc(j)  points to the start of column j.
   lu1or2  is derived from  mc20ad,  a routine in the Harwell
   Subroutine Library, author J. K. Reid.
   ------------------------------------------------------------------
   xx Feb 1985: Original version.
   17 Oct 2000: a, inum, jnum now have size lena to allow nelem = 0.
   ================================================================== */
template<class T>
void LUSOLrec<T>::LU1OR2()
{
    T    ACE;
    T    ACEP;
    INT  L, J, I, JCE, ICE, ICEP, JCEP, JA, JB;

    // Set  loc(j)  to point to the beginning of column  j.
    L = 1;
    for(J = 1; J <= this->n; J++) 
    {
        this->locc[J] = L;
        L += this->lenc[J];
    }
    //  Sort the elements into column order.
    //    The algorithm is an in-place sort and is of order  numa.
    for(I = 1; I <= this->nelem; I++) 
    {
        // Establish the current entry.
        JCE = this->indr[I];
        if(JCE==0)
            continue;
        ACE = this->a[I];
        ICE = this->indc[I];
        this->indr[I] = 0;
        // Chain from current entry.
        for(J = 1; J <= this->nelem; J++) 
        {
            // The current entry is not in the correct position.
            //   Determine where to store it.
            L = this->locc[JCE];
            this->locc[JCE]++;
            // Save the contents of that location.
            ACEP = this->a[L];
            ICEP = this->indc[L];
            JCEP = this->indr[L];
            // Store current entry.
            this->a[L] = ACE;
            this->indc[L] = ICE;
            this->indr[L] = 0;
            // If next current entry needs to be processed,
            //  copy it into current entry.
            if(JCEP==0)
                break;
            ACE = ACEP;
            ICE = ICEP;
            JCE = JCEP;
        }
    }
    // Reset loc(j) to point to the start of column j.
    JA = 1;
    for(J = 1; J <= this->n; J++) 
    {
        JB = this->locc[J];
        this->locc[J] = JA;
        JA = JB;
    }
};

/* ==================================================================
   lu1or3  looks for duplicate elements in an  m by n  matrix  A
   defined by the column list  indc, lenc, locc.
   iw  is used as a work vector of length  m.
   ------------------------------------------------------------------
   xx Feb 1985: Original version.
   17 Oct 2000: indc, indr now have size lena to allow nelem = 0.
   ================================================================== */
template<class T>
void LUSOLrec<T>::LU1OR3(INT *LERR, INT *INFORM)
{
    INT I, J, L1, L2, L;

    memclear((this->ip+1), this->m);

    for(J = 1; J <= this->n; J++) 
    {
        if(this->lenc[J]>0) 
        {
            L1 = this->locc[J];
            L2 = (L1+this->lenc[J])-1;
            for(L = L1; L <= L2; L++) 
            {
                I = this->indc[L];
                if(this->ip[I]==J)
                    goto x910;
                this->ip[I] = J;
            }
        }
    }
    *INFORM = LUSOL_INFORM_LUSUCCESS;
    return;
  x910:
    *LERR = L;
    *INFORM = LUSOL_INFORM_LUSINGULAR;
};

/* ==================================================================
   lu1or4 constructs a row list  indr, locr
   from a corresponding column list  indc, locc,
   given the lengths of both columns and rows in  lenc, lenr.
   ------------------------------------------------------------------
   xx Feb 1985: Original version.
   17 Oct 2000: indc, indr now have size lena to allow nelem = 0.
   ================================================================== */
template<class T>
void LUSOLrec<T>::LU1OR4()
{
    INT L, I, L2, J, JDUMMY, L1, LR;

    // Initialize  locr(i)  to point just beyond where the
    // last component of row  i  will be stored.
    L = 1;
    for(I = 1; I <= this->m; I++) 
    {
        L += this->lenr[I];
        this->locr[I] = L;
    }
    // By processing the columns backwards and decreasing  locr(i)
    //    each time it is accessed, it will end up pointing to the
    //    beginning of row  i  as required.
    L2 = this->nelem;
    J = this->n+1;
    for(JDUMMY = 1; JDUMMY <= this->n; JDUMMY++) 
    {
        J = J-1;
        if(this->lenc[J]>0) 
        {
            L1 = this->locc[J];
            for(L = L1; L <= L2; L++) 
            {
                I = this->indc[L];
                LR = this->locr[I]-1;
                this->locr[I] = LR;
                this->indr[LR] = J;
            }
            L2 = L1-1;
        }
    }
};

/* ==================================================================
   lu1pen deals with pending fill-in in the row file.
   ------------------------------------------------------------------
   ifill(ll) says if a row involved in the new column of L
             has to be updated.  If positive, it is the total
             length of the final updated row.
   jfill(lu) says if a column involved in the new row of U
             contains any pending fill-ins.  If positive, it points
             to the first fill-in in the column that has yet to be
             added to the row file.
   ------------------------------------------------------------------
   16 Apr 1989: First version of lu1pen.
   23 Mar 2001: ilast used and updated.
   ================================================================== */
template<class T>
void LUSOLrec<T>::LU1PEN(INT NSPARE, INT *ILAST,
            INT LPIVC1, INT LPIVC2, INT LPIVR1, INT LPIVR2,
            INT *LROW, INT IFILL[], INT JFILL[])
{
    INT  LL, LC, L, I, LR1, LR2, LR, LU, J, LC1, LC2, LAST;

    LL = 0;
    for(LC = LPIVC1; LC <= LPIVC2; LC++) 
    {
        LL++;
        if(IFILL[LL]==0)
            continue;
        // Another row has pending fill.
        // First, add some spare space at the }
        // of the current last row.

        LC1 = (*LROW)+1;
        LC2 = (*LROW)+NSPARE;
        *LROW = LC2;
        for(L = LC1; L <= LC2; L++) 
        {
            this->indr[L] = 0;
        }
        // Now move row i to the end of the row file.
        I = this->indc[LC];
        *ILAST = I;
        LR1 = this->locr[I];
        LR2 = (LR1+this->lenr[I])-1;
        this->locr[I] = (*LROW)+1;
        for(LR = LR1; LR <= LR2; LR++) 
        {
            (*LROW)++;
            this->indr[*LROW] = this->indr[LR];
            this->indr[LR] = 0;
        }
        (*LROW) += IFILL[LL];
    }
    // Scan all columns of  D  and insert the pending fill-in
    // into the row file.
    LU = 1;
    for(LR = LPIVR1; LR <= LPIVR2; LR++) 
    {
        LU++;
        if(JFILL[LU]==0)
            continue;
        J = this->indr[LR];
        LC1 = (this->locc[J]+JFILL[LU])-1;
        LC2 = (this->locc[J]+this->lenc[J])-1;
        for(LC = LC1; LC <= LC2; LC++) 
        {
            I = this->indc[LC]-this->m;
            if(I>0) 
            {   
                this->indc[LC] = I;
                LAST = this->locr[I]+this->lenr[I];
                this->indr[LAST] = J;
                this->lenr[I]++;
            }
        }
    }
};


/* ==================================================================
   lu1fad  is a driver for the numerical phase of lu1fac.
   At each stage it computes a column of  L  and a row of  U,
   using a Markowitz criterion to select the pivot element,
   subject to a stability criterion that bounds the elements of  L.
   ------------------------------------------------------------------
   Local variables
   ---------------
   lcol   is the length of the column file.  It points to the last
          nonzero in the column list.
   lrow   is the analogous quantity for the row file.
   lfile  is the file length (lcol or lrow) after the most recent
          compression of the column list or row list.
   nrowd  and  ncold  are the number of rows and columns in the
          matrix defined by the pivot column and row.  They are the
          dimensions of the submatrix D being altered at this stage.
   melim  and  nelim  are the number of rows and columns in the
          same matrix D, excluding the pivot column and row.
   mleft  and  nleft  are the number of rows and columns
          still left to be factored.
   nzchng is the increase in nonzeros in the matrix that remains
          to be factored after the current elimination
          (usually negative).
   nzleft is the number of nonzeros still left to be factored.
   nspare is the space we leave at the end of the last row or
          column whenever a row or column is being moved to the }
          of its file.  nspare = 1 or 2 might help reduce the
          number of file compressions when storage is tight.
   The row and column ordering permutes A into the form
                      ------------------------
                       \                     |
                        \         U1         |
                         \                   |
                          --------------------
                          |\
                          | \
                          |  \
          P A Q   =       |   \
                          |    \
                          |     --------------
                          |     |            |
                          |     |            |
                          | L1  |     A2     |
                          |     |            |
                          |     |            |
                          --------------------
   where the block A2 is factored as  A2 = L2 U2.
   The phases of the factorization are as follows.
   Utri   is true when U1 is being determined.
          Any column of length 1 is accepted immediately (if TPP).
   Ltri   is true when L1 is being determined.
          lu1mar exits as soon as an acceptable pivot is found
          in a row of length 1.
   spars1 is true while the density of the (modified) A2 is less
          than the parameter dens1 = parmlu(7) = 0.3 say.
          lu1mar searches maxcol columns and maxrow rows,
          where  maxcol = luparm(3),  maxrow = maxcol - 1.
          lu1mxc is used to keep the biggest element at the top
          of all remaining columns.
   spars2 is true while the density of the modified A2 is less
          than the parameter dens2 = parmlu(8) = 0.6 say.
          lu1mar searches maxcol columns and no rows.
          lu1mxc could fix up only the first maxcol cols (with TPP).
          22 Sep 2000:  For simplicity, lu1mxc fixes all
                        modified cols.
   dense  is true once the density of A2 reaches dens2.
          lu1mar searches only 1 column (the shortest).
          lu1mxc could fix up only the first column (with TPP).
   ------------------------------------------------------------------
   00 Jan 1986  Version documented in LUSOL paper:
                Gill, Murray, Saunders and Wright (1987),
                Maintaining LU factors of a general sparse matrix,
                Linear algebra and its applications 88/89, 239-270.
   02 Feb 1989  Following Suhl and Aittoniemi (1987), the largest
                element in each column is now kept at the start of
                the column, i.e. in position locc(j) of a and indc.
                This should speed up the Markowitz searches.
                To save time on highly triangular matrices, we wait
                until there are no further columns of length 1
                before setting and maintaining that property.
   12 Apr 1989  ipinv and iqinv added (inverses of ip and iq)
                to save searching ip and iq for rows and columns
                altered in each elimination step.  (Used in lu1pq2)
   19 Apr 1989  Code segmented to reduce its size.
                lu1gau does most of the Gaussian elimination work.
                lu1mar does just the Markowitz search.
                lu1mxc moves biggest elements to top of columns.
                lu1pen deals with pending fill-in in the row list.
                lu1pq2 updates the row and column permutations.
   26 Apr 1989  maxtie replaced by maxcol, maxrow in the Markowitz
                search.  maxcol, maxrow change as density increases.
   25 Oct 1993  keepLU implemented.
   07 Feb 1994  Exit main loop early to finish off with a dense LU.
                densLU tells lu1fad whether to do it.
   21 Dec 1994  Bug fixed.  nrank was wrong after the call to lu1ful.
   12 Nov 1999  A parallel version of dcopy gave trouble in lu1ful
                during left-shift of dense matrix D within a(*).
                Fixed this unexpected problem here in lu1fad
                by making sure the first and second D don't overlap.
   13 Sep 2000  TCP (Threshold Complete Pivoting) implemented.
                lu2max added
                (finds aijmax from biggest elems in each col).
                Utri, Ltri and Spars1 phases apply.
                No switch to Dense CP yet.  (Only TPP switches.)
   14 Sep 2000  imax needed to remember row containing aijmax.
   22 Sep 2000  For simplicity, lu1mxc always fixes all modified cols.
                (TPP spars2 used to fix just the first maxcol cols.)
   08 Nov 2000: Speed up search for aijmax.
                Don't need to search all columns if the elimination
                didn't alter the col containing the current aijmax.
   21 Nov 2000: lu1slk implemented for Utri phase with TCP
                to guard against deceptive triangular matrices.
                (Utri used to have aijtol >= 0.9999 to include
                slacks, but this allows other 1s to be accepted.)
                Utri now accepts slacks, but applies normal aijtol
                test to other pivots.
   28 Nov 2000: TCP with empty cols must call lu1mxc and lu2max
                with ( lq1, n, ... ), not just ( 1, n, ... ).
   23 Mar 2001: lu1fad bug with TCP.
                A col of length 1 might not be accepted as a pivot.
                Later it appears in a pivot row and temporarily
                has length 0 (when pivot row is removed
                but before the column is filled in).  If it is the
                last column in storage, the preceding col also thinks
                it is "last".  Trouble arises when the preceding col
                needs fill-in -- it overlaps the real "last" column.
                (Very rarely, same trouble might have happened if
                the drop tolerance caused columns to have length 0.)
                Introduced ilast to record the last row in row file,
                           jlast to record the last col in col file.
                lu1rec returns ilast = indr(lrow + 1)
                            or jlast = indc(lcol + 1).
                (Should be an output parameter, but didn't want to
                alter lu1rec's parameter list.)
                lu1rec also treats empty rows or cols safely.
                (Doesn't eliminate them!)
   26 Apr 2002: Heap routines added for TCP.
                lu2max no longer needed.
                imax, jmax used only for printing.
   01 May 2002: lu1DCP implemented (dense complete pivoting).
                Both TPP and TCP now switch to dense LU
                when density exceeds dens2.
   06 May 2002: In dense mode, store diag(U) in natural order.
   09 May 2002: lu1mCP implemented (Markowitz TCP via heap).
   11 Jun 2002: lu1mRP implemented (Markowitz TRP).
   28 Jun 2002: Fixed call to lu1mxr.
   14 Dec 2002: lu1mSP implemented (Markowitz TSP).
   15 Dec 2002: Both TPP and TSP can grab cols of length 1
                during Utri.
   19 Dec 2004: Hdelete(...) has new input argument Hlenin.
   26 Mar 2006: lu1fad returns nrank  = min( mrank, nrank )
                and ignores nsing from lu1ful
   26 Mar 2006: Allow for empty columns before calling Hbuild.
   03 Apr 2013: f90 lu1mxr recoded to improve efficiency of TRP.
   06 Jun 2013: Adapted f90 lu1mxr for use in this f77 version.
   07 Jul 2013: cols, markc, markr are new work arrays for lu1mxr.
   ================================================================== */
template<class T>
void LUSOLrec<T>::LU1FAD(INT *INFORM, INT *LENL, INT *LENU, INT *MINLEN,
            INT *MERSUM, INT *NUTRI, INT *NLTRI,
            INT *NDENS1, INT *NDENS2, INT *NRANK,
            TR *LMAX, TR *UMAX, TR *DUMAX, TR *DUMIN, TR *AKMAX)
{
    bool UTRI, LTRI, SPARS1, SPARS2, DENSE_LOC, DENSLU, KEEPLU, TCP, TPP, TRP,TSP;
    INT  HLEN, LPIV, LPRINT, MAXCOL, MAXROW, ILAST, JLAST, LFILE, LROW, LCOL,
         MINMN, MAXMN, NZLEFT, NSPARE, LU1, KK, J, LC, MLEFT, NLEFT, NROWU,
         LQ1, LQ2, JBEST, LQ, I, IBEST, MBEST, LEND, NFREE, LD, NCOLD, NROWD,
         MELIM, NELIM, JMAX, IMAX, LL1, LSAVE, LFREE, LIMIT, MINFRE, LPIVR, LPIVR1, LPIVR2,
         L, LPIVC, LPIVC1, LPIVC2, KBEST, LU, LR, LENJ, LC1, LAST, LL, LS,
         LENI, LR1, LFIRST, NFILL, NZCHNG, K, MRANK, NSING;
    TR   SMALL, USPACE, DENS1, DENS2;
    TR   LIJ, AIJMAX, AIJTOL, AMAX, DIAG, V, LTOL;
    T    ABEST;
    INT  LENA2 = this->lena;

    using RT    = typename lusol_real_type<T>::type;

    #ifdef UseTimer
        INT    eltime, mktime, ntime;
        LUSOL_timer (LUSOL, 3, "start" );
        ntime = this->n / 4;
    #endif

    AIJMAX = 0;
    AIJTOL = 0;
    HLEN   = 0;
    JBEST  = 0;
    IBEST  = 0;
    MBEST  = 0;
    LEND   = 0;
    LD     = 0;

    LPRINT = this->luparm[LUSOL_IP_PRINTLEVEL];
    MAXCOL = this->luparm[LUSOL_IP_MARKOWITZ_MAXCOL];
    LPIV   = this->luparm[LUSOL_IP_PIVOTTYPE];
    KEEPLU = (bool) (this->luparm[LUSOL_IP_KEEPLU]!=false);
    // Threshold Partial   Pivoting (normal).
    TPP = (bool) (LPIV==LUSOL_PIVOT_TPP);
    // Threshold Rook      Pivoting
    TRP = (bool) (LPIV==LUSOL_PIVOT_TRP);
    // Threshold Complete  Pivoting.
    TCP = (bool) (LPIV==LUSOL_PIVOT_TCP);
    // Threshold Symmetric Pivoting.
    TSP = (bool) (LPIV==LUSOL_PIVOT_TSP);
    DENSLU = false;
    MAXROW = MAXCOL-1;
    // Assume row m is last in the row file.
    ILAST = this->m;
    // Assume col n is last in the col file.
    JLAST = this->n;
    LFILE = this->nelem;
    LROW = this->nelem;
    LCOL = this->nelem;
    MINMN = minimum(this->m,this->n);
    MAXMN = maximum(this->m,this->n);
    NZLEFT = this->nelem;
    NSPARE = 1;

    if(KEEPLU)
        LU1 = LENA2+1;
    else 
    {
        // Store only the diagonals of U in the top of memory.
        #ifdef ClassicdiagU
            LDIAGU = LENA2-this->n;
            LU1 = LDIAGU+1;
            this->diagU = this->a+LDIAGU;
        #else
            LU1 = LENA2+1;
        #endif
    }

    LTOL    = this->parmlu[LUSOL_RP_FACTORMAX_Lij];
    SMALL   = this->parmlu[LUSOL_RP_ZEROTOLERANCE];
    USPACE  = this->parmlu[LUSOL_RP_COMPSPACE_U];
    DENS1   = this->parmlu[LUSOL_RP_MARKOWITZ_CONLY];
    DENS2   = this->parmlu[LUSOL_RP_MARKOWITZ_DENSE];
    UTRI    = true;
    LTRI    = false;
    SPARS1  = false;
    SPARS2  = false;
    DENSE_LOC   = false;
    
	// Check parameters.
    LTOL = maximum(LTOL,TR(1.0001E+0));
    DENS1 = minimum(DENS1,DENS2);
	
    // Initialize output parameters.
    //    lenL, lenU, minlen, mersum, nUtri, nLtri, ndens1, ndens2, nrank
    //    are already initialized by lu1fac.
    *LMAX  = 0;
    *UMAX  = 0;
    *DUMAX = 0;
    *DUMIN = lusol_user_params<TR>::LUSOL_BIGNUM;
    if(this->nelem==0)
        *DUMIN = 0;
		
    *AKMAX      = 0;
	INT MARK 	= 0;
	
    //  More initialization.
    //    Don't worry yet about lu1mxc.
    if(TPP || TSP) 
    {
        AIJMAX = 0;
        AIJTOL = 0;
        HLEN = 1;
        // TRP or TCP
    }
    else 
    {
        // Move biggest element to top of each column.
        // Set w(*) to mark slack columns (unit vectors).
        LU1MXC(1,this->n,this->iq);
        LU1SLK();
    }
    if(TRP)
	{
        // Find biggest element in each row.
		MARK = 0;
        LU1MXR(MARK, 1, this->m, this->ip, this->m_cols, this->m_markc, this->m_markr, this->amaxr);
	};
    if(TCP) 
    {
        //  Set Ha(1:Hlen) = biggest element in each column,
        //  Hj(1:Hlen) = corresponding column indices.
        HLEN = 0;
        for(KK = 1; KK <= this->n; KK++) 
        {
            HLEN++;
            J = this->iq[KK];
            LC = this->locc[J];
            if (this->lenc[J] == 0)
            {
                this->m_heap.set(HLEN,J,0.);
            }
            else
            {
                this->m_heap.set(HLEN,J,abs(this->a[LC]));
            };
        }
        // Build the heap, creating new Ha, Hj and setting Hk(1:Hlen).
        this->m_heap.build(HLEN);
    };

    //------------------------------------------------------------------
    //    Start of main loop.
    //------------------------------------------------------------------
    MLEFT = this->m+1;
    NLEFT = this->n+1;
    for(NROWU = 1; NROWU <= MINMN; NROWU++) 
    {
        #ifdef UseTimer
            mktime = (NROWU / ntime) + 4;
            eltime = (NROWU / ntime) + 9;
        #endif
        MLEFT--;
        NLEFT--;
        // Bail out if there are no nonzero rows left.
        if(this->iploc[1]>this->m)
            goto x900;
        // For TCP, the largest Aij is at the top of the heap.
        if(TCP) 
        {
            // Marvelously easy
            AIJMAX = this->m_heap.get_Ha()[1];
            *AKMAX = maximum(*AKMAX,AIJMAX);
            AIJTOL = AIJMAX/LTOL;
        }
        //===============================================================
        //       Find a suitable pivot element.
        //===============================================================
        if(UTRI) 
        {
            //------------------------------------------------------------
            //      So far all columns have had length 1.
            //      We are still looking for the (backward) triangular part of A
            //      that forms the first rows and columns of U.
            //------------------------------------------------------------
            LQ1 = this->iqloc[1];
            LQ2 = this->n;
            if(this->m>1)
                LQ2 = this->iqloc[2]-1;
            // There are more cols of length 1.
            if(LQ1<=LQ2) 
            {
                if(TPP) 
                {
                    // Grab the first one.
                    JBEST = this->iq[LQ1];                    
                }
                else if (TSP)
                {
                    JBEST = 0;
                    for(LQ = LQ1; LQ <= LQ2; LQ++) 
                    {
                        int tmp_JBEST   = this->iq[LQ];
                        int LC2         = this->locc[tmp_JBEST];
                        int tmp_IBEST   = this->indc[LC2];

                        if (tmp_JBEST == tmp_IBEST)
                        {
                            JBEST = tmp_JBEST;
                            break;
                        };
                    };
                }
                else 
                {
                    // Scan all columns of length 1 ... TRP or TCP
                    JBEST = 0;
                    for(LQ = LQ1; LQ <= LQ2; LQ++) 
                    {
                        J = this->iq[LQ];
                        // Accept a slack
                        if(this->wr[J]>0) 
                        {
                            JBEST = J;
                            goto x250;
                        }
                        LC = this->locc[J];
                        AMAX = abs(this->a[LC]);
                        if(TRP) 
                        {
                            I = this->indc[LC];
                            AIJTOL = this->amaxr[I]/LTOL;
                        }
                        if(AMAX>=AIJTOL) 
                        {
                            JBEST = J;
                            goto x250;
                        }
                    }
                }
              x250:
                if(JBEST>0) 
                {
                    LC = this->locc[JBEST];
                    IBEST = this->indc[LC];
                    MBEST = 0;
                    goto x300;
                }
            }
            // This is the end of the U triangle.
            // We will not return to this part of the code.
            // TPP and TSP call lu1mxc for the first time
            // (to move biggest element to top of each column).
            if(LPRINT>=LUSOL_MSG_PIVOT)
                LUSOL_report(this, 0, (char*) "Utri ended.  spars1 = true\n"); //FIXME:
            UTRI = false;
            LTRI = true;
            
            SPARS1 = true;
            *NUTRI = NROWU-1;
            
            if(TPP || TSP)
                LU1MXC(LQ1,this->n,this->iq);
        }
        if(SPARS1) 
        {
            //------------------------------------------------------------
            //      Perform a Markowitz search.
            //      Search cols of length 1, then rows of length 1,
            //      then   cols of length 2, then rows of length 2, etc.
            //------------------------------------------------------------
            #ifdef UseTimer
                LUSOL_timer (LUSOL, mktime, "start" );
            #endif
            if(TPP || TCP) 
            {
                LU1MAR(MAXMN,TCP,AIJTOL,LTOL,MAXCOL,MAXROW,&IBEST,&JBEST,&MBEST);
            }
            else if(TRP) 
            {
                LU1MRP(MAXMN,LTOL,MAXCOL,MAXROW,&IBEST,&JBEST,&MBEST,this->amaxr);
            }
            else if(TSP) 
            {
                LU1MSP(MAXMN,LTOL,MAXCOL,&IBEST,&JBEST,&MBEST);
                if(IBEST==0)
                {
                    goto x900;
                };
            }

            #ifdef UseTimer
                LUSOL_timer (LUSOL, mktime, "finish" );
            #endif
            if(LTRI) 
            {
                // So far all rows have had length 1.
                //     We are still looking for the (forward) triangle of A
                //     that forms the first rows and columns of L.
                if(MBEST>0) 
                {
                    LTRI = false;
                    *NLTRI = NROWU-1-*NUTRI;
                    if(LPRINT>=LUSOL_MSG_PIVOT)
                        LUSOL_report(this, 0, (char*) "Ltri ended.\n"); //FIXME:
                }
            }
            else 
            {
                // See if what's left is as dense as dens1.
                if(NZLEFT>=(DENS1*MLEFT)*NLEFT) 
                {
                    SPARS1 = false;
                    SPARS2 = true;
                    *NDENS1 = NLEFT;
                    MAXROW = 0;
                    if(LPRINT>=LUSOL_MSG_PIVOT)
                        LUSOL_report(this, 0, (char*) "spars1 ended.  spars2 = true\n"); //FIXME:
                }
            }
        }
        else if(SPARS2 || DENSE_LOC) 
        {
            //------------------------------------------------------------
            //     Perform a restricted Markowitz search,
            //     looking at only the first maxcol columns.  (maxrow = 0.)
            //------------------------------------------------------------
            #ifdef UseTimer
                LUSOL_timer (LUSOL, mktime, "start" );
            #endif
            if(TPP || TCP) 
            {
                LU1MAR(MAXMN,TCP,AIJTOL,LTOL,MAXCOL,MAXROW,&IBEST,&JBEST,&MBEST);
            }
            else if(TRP) 
            {
                LU1MRP(MAXMN,LTOL,MAXCOL,MAXROW,&IBEST,&JBEST,&MBEST,this->amaxr);
            }
            else if(TSP) 
            {
                LU1MSP(MAXMN,LTOL,MAXCOL,&IBEST,&JBEST,&MBEST);
                if(IBEST==0)
                {
                    goto x900;
                };
            }

            #ifdef UseTimer
                LUSOL_timer (LUSOL, mktime, "finish" );
            #endif
            // See if what's left is as dense as dens2.
            if(SPARS2) 
            {
                if(NZLEFT>=(DENS2*MLEFT)*NLEFT) 
                {
                    SPARS2      = false;
                    DENSE_LOC   = true;
                    *NDENS2     = NLEFT;
                    MAXCOL      = 1;

                    if(LPRINT>=LUSOL_MSG_PIVOT)
                        LUSOL_report(this, 0, (char*) "spars2 ended.  dense = true\n"); //FIXME:
                }
            }
        }
        //---------------------------------------------------------------
        //       See if we can finish quickly.
        //---------------------------------------------------------------
        if(DENSE_LOC) 
        {
            LEND = MLEFT*NLEFT;
            NFREE = LU1-1;
            if(NFREE>=2*LEND) 
            {
                // There is room to treat the remaining matrix as
                //     a dense matrix D.
                //     We may have to compress the column file first.
                //     12 Nov 1999: D used to be put at the
                //                  beginning of free storage (lD = lcol + 1).
                //                  Now put it at the end     (lD = lu1 - lenD)
                //                  so the left-shift in lu1ful will not
                //                  involve overlapping storage
                //                  (fatal with parallel dcopy).
                DENSLU = true;
                *NDENS2 = NLEFT;
                LD = LU1-LEND;
                if(LCOL>=LD) 
                {
                    LU1REC(this->n,true,&LCOL,this->indc,this->lenc,this->locc);
                    LFILE = LCOL;
                    JLAST = this->indc[LCOL+1];
                }
                goto x900;
            }
        }
      //===============================================================
      //       The best  aij  has been found.
      //       The pivot row  ibest  and the pivot column  jbest
      //       Define a dense matrix  D  of size  nrowd  by  ncold.
      //===============================================================
      x300:
        NCOLD = this->lenr[IBEST];
        NROWD = this->lenc[JBEST];
        MELIM = NROWD-1;
        NELIM = NCOLD-1;
        (*MERSUM) += MBEST;
        (*LENL) += MELIM;
        (*LENU) += NCOLD;
        if(LPRINT>=LUSOL_MSG_PIVOT) 
        {
            if(NROWU==1) //FIXME:
                LUSOL_report(this, 0, (char*) "lu1fad debug:\n");
            if(TPP || TRP || TSP) 
            {//FIXME:
                LUSOL_report(this, 0, (char*) "nrowu:%7d   i,jbest:%7d,%7d   nrowd,ncold:%6d,%6d\n",
                                NROWU, IBEST,JBEST, NROWD,NCOLD);
            }
            else 
            {
                JMAX = this->m_heap.get_Hj()[1];
                IMAX = this->indc[this->locc[JMAX]];//FIXME:
                LUSOL_report(this, 0, (char*) "nrowu:%7d   i,jbest:%7d,%7d   nrowd,ncold:%6d,%6d   i,jmax:%7d,%7d   aijmax:%g\n",
                                NROWU, IBEST,JBEST, NROWD,NCOLD, IMAX,JMAX, AIJMAX);
            }
        }
        //===============================================================
        //       Allocate storage for the next column of  L  and next row of  U.
        //       Initially the top of a, indc, indr are used as follows:
        //                  ncold       melim       ncold        melim
        //       a      |...........|...........|ujbest..ujn|li1......lim|
        //       indc   |...........|  lenr(i)  |  lenc(j)  |  markl(i)  |
        //       indr   |...........| iqloc(i)  |  jfill(j) |  ifill(i)  |
        //             ^           ^             ^           ^            ^
        //             lfree   lsave             lu1         ll1          oldlu1
        //       Later the correct indices are inserted:
        //       indc   |           |           |           |i1........im|
        //       indr   |           |           |jbest....jn|ibest..ibest|
        //===============================================================
        if(!KEEPLU) 
        {
            // Always point to the top spot.
            //      Only the current column of L and row of U will
            //      take up space, overwriting the previous ones.
            LU1 = LENA2+1;
        }
        // Update (left-shift) pointers to make room for the new data
        LL1 = LU1-MELIM;
        LU1 = LL1-NCOLD;
        LSAVE = LU1-NROWD;
        LFREE = LSAVE-NCOLD;

        // Check if we need to allocate more memory, and allocate if necessary
        L  = maximum(LROW, LCOL) + 2*(this->m+this->n);
        L *= LUSOL_MULT_nz_a;
        L = maximum(L, NROWD*NCOLD);

        // Do the memory expansion
        if(L > LFREE-LCOL) 
        {
            if (LUSOL_expand_a(&L, &LFREE))
            {
                LL1   += L;
                LU1   += L;
                LSAVE += L;
                #ifdef ClassicdiagU
                    this->diagU += L;
                #endif
            }
        }
        LIMIT = (INT) (USPACE*LFILE)+this->m+this->n+1000;

        // Make sure the column file has room.
        // Also force a compression if its length exceeds a certain limit.
        MINFRE = NROWD*NCOLD;
        NFREE = LFREE-LCOL;
        if(NFREE<MINFRE || LCOL>LIMIT) 
        {
            LU1REC(this->n,true,&LCOL,this->indc,this->lenc,this->locc);
            LFILE = LCOL;
            JLAST = this->indc[LCOL+1];
            NFREE = LFREE-LCOL;
            if(NFREE<MINFRE)
                goto x970;
        }
        // Make sure the row file has room.
        MINFRE = NROWD*NCOLD;
        NFREE = LFREE-LROW;
        if(NFREE<MINFRE || LROW>LIMIT) 
        {
            LU1REC(this->m,false,&LROW,this->indr,this->lenr,this->locr);
            LFILE = LROW;
            ILAST = this->indr[LROW+1];
            NFREE = LFREE-LROW;
            if(NFREE<MINFRE)
                goto x970;
        }
        //===============================================================
        //       Move the pivot element to the front of its row
        //       and to the top of its column.
        //===============================================================
        LPIVR = this->locr[IBEST];
        LPIVR1 = LPIVR+1;
        LPIVR2 = LPIVR+NELIM;
        for(L = LPIVR; L <= LPIVR2; L++) 
        {
            if(this->indr[L]==JBEST)
                break;
        }

        this->indr[L] = this->indr[LPIVR];
        this->indr[LPIVR] = JBEST;
        LPIVC = this->locc[JBEST];
        LPIVC1 = LPIVC+1;
        LPIVC2 = LPIVC+MELIM;
        for(L = LPIVC; L <= LPIVC2; L++) 
        {
            if(this->indc[L]==IBEST)
                break;
        }
        this->indc[L] = this->indc[LPIVC];
        this->indc[LPIVC] = IBEST;
        ABEST = this->a[L];
        this->a[L] = this->a[LPIVC];
        this->a[LPIVC] = ABEST;
        if(!KEEPLU)
            // Store just the diagonal of U, in natural order.
            //!!         a[ldiagU + nrowu] = abest ! This was in pivot order.
            this->diagU[JBEST] = ABEST;

        //==============================================================
        //    Delete pivot col from heap.
        //    Hk tells us where it is in the heap.
        //==============================================================
        if(TCP) 
        {
            KBEST = this->m_heap.get_Hk()[JBEST];
            this->m_heap.remove(&HLEN,KBEST);
        }
        //===============================================================
        //       Delete the pivot row from the column file
        //       and store it as the next row of  U.
        //       set  indr(lu) = 0     to initialize jfill ptrs on columns of D,
        //            indc(lu) = lenj  to save the original column lengths.
        //===============================================================
        this->a[LU1] = ABEST;
        this->indr[LU1] = JBEST;
        this->indc[LU1] = NROWD;
        LU = LU1;
        DIAG = abs(ABEST);
        *UMAX = maximum(*UMAX,DIAG);
        *DUMAX = maximum(*DUMAX,DIAG);
        *DUMIN = minimum(*DUMIN,DIAG);
        for(LR = LPIVR1; LR <= LPIVR2; LR++) 
        {
            LU++;
            J = this->indr[LR];
            LENJ = this->lenc[J];
            this->lenc[J] = LENJ-1;
            LC1 = this->locc[J];
            LAST = LC1+this->lenc[J];
            for(L = LC1; L <= LAST; L++) 
            {
                if(this->indc[L]==IBEST)
                    break;
            }
            this->a[LU] = this->a[L];
            this->indr[LU] = 0;
            this->indc[LU] = LENJ;
            *UMAX = maximum<TR>(*UMAX,abs(this->a[LU]));
            this->a[L] = this->a[LAST];
            this->indc[L] = this->indc[LAST];
            // Free entry
            this->indc[LAST] = 0;
        }
        //===============================================================
        //       Delete the pivot column from the row file
        //       and store the nonzeros of the next column of  L.
        //       Set  indc(ll) = 0     to initialize markl(*) markers,
        //            indr(ll) = 0     to initialize ifill(*) row fill-in cntrs,
        //            indc(ls) = leni  to save the original row lengths,
        //            indr(ls) = iqloc(i)    to save parts of  iqloc(*),
        //            iqloc(i) = lsave - ls  to point to the nonzeros of  L
        //                     = -1, -2, -3, ... in mark(*).
        //===============================================================
        this->indc[LSAVE] = NCOLD;
        if(MELIM==0)
        {
            goto x700;
        };
        LL = LL1-1;
        LS = LSAVE;
        ABEST = RT(1.)/ABEST;
        for(LC = LPIVC1; LC <= LPIVC2; LC++) 
        {
            LL++;
            LS++;
            I = this->indc[LC];
            LENI = this->lenr[I];
            this->lenr[I] = LENI-1;
            LR1 = this->locr[I];
            LAST = LR1+this->lenr[I];
            for(L = LR1; L <= LAST; L++) 
            {
                if(this->indr[L]==JBEST)
                    break;
            }
            this->indr[L] = this->indr[LAST];
            // Free entry
            this->indr[LAST] = 0;
            this->a[LL] = -this->a[LC]*ABEST;
            LIJ = abs(this->a[LL]);
            *LMAX = maximum(*LMAX,LIJ);
            this->indc[LL] = 0;
            this->indr[LL] = 0;
            this->indc[LS] = LENI;
            this->indr[LS] = this->iqloc[I];
            this->iqloc[I] = LSAVE-LS;
        }
        //===============================================================
        //       Do the Gaussian elimination.
        //       This involves adding a multiple of the pivot column
        //       to all other columns in the pivot row.
        //       Sometimes more than one call to lu1gau is needed to allow
        //       compression of the column file.
        //       lfirst  says which column the elimination should start with.
        //       minfre  is a bound on the storage needed for any one column.
        //       lu      points to off-diagonals of u.
        //       nfill   keeps track of pending fill-in in the row file.
        //===============================================================
        if(NELIM==0)
            goto x700;
        LFIRST = LPIVR1;
        MINFRE = MLEFT+NSPARE;
        LU = 1;
        NFILL = 0;

      x400:
        #ifdef UseTimer
            LUSOL_timer (LUSOL, eltime, "start" );
        #endif
        LU1GAU(MELIM,NSPARE,SMALL,LPIVC1,LPIVC2,&LFIRST,LPIVR2,
               LFREE,MINFRE,ILAST,&JLAST,&LROW,&LCOL,&LU,&NFILL,
               this->iqloc, this->a+LL1-LUSOL_ARRAYOFFSET,
               this->indc+LL1-LUSOL_ARRAYOFFSET, this->a+LU1-LUSOL_ARRAYOFFSET,
               this->indr+LL1-LUSOL_ARRAYOFFSET, this->indr+LU1-LUSOL_ARRAYOFFSET);
        #ifdef UseTimer
            LUSOL_timer (LUSOL, eltime, "finish" );
        #endif
        if(LFIRST>0) 
        {
            // The elimination was interrupted.
            //       Compress the column file and try again.
            //       lfirst, lu and nfill have appropriate new values.
            LU1REC(this->n,true,&LCOL, this->indc,this->lenc,this->locc);
            LFILE = LCOL;
            JLAST = this->indc[LCOL+1];
            LPIVC = this->locc[JBEST];
            LPIVC1 = LPIVC+1;
            LPIVC2 = LPIVC+MELIM;
            NFREE = LFREE-LCOL;
            if(NFREE<MINFRE) 
            {
                goto x970;
            }
            goto x400;
        }
        //===============================================================
        //       The column file has been fully updated.
        //       Deal with any pending fill-in in the row file.
        //===============================================================
        if(NFILL>0) 
        {
            //  Compress the row file if necessary.
            //      lu1gau has set nfill to be the number of pending fill-ins
            //      plus the current length of any rows that need to be moved.
            MINFRE = NFILL;
            NFREE = LFREE-LROW;
            if(NFREE<MINFRE) 
            {
                LU1REC(this->m,false,&LROW,this->indr,this->lenr,this->locr);
                LFILE = LROW;
                ILAST = this->indr[LROW+1];
                LPIVR = this->locr[IBEST];
                LPIVR1 = LPIVR+1;
                LPIVR2 = LPIVR+NELIM;
                NFREE = LFREE-LROW;
                if(NFREE<MINFRE) 
                {
                    goto x970;
                }
            }
            // Move rows that have pending fill-in to end of the row file.
            //       Then insert the fill-in.
            LU1PEN(NSPARE,&ILAST, LPIVC1,LPIVC2,LPIVR1,LPIVR2,
                 &LROW,this->indr+LL1-LUSOL_ARRAYOFFSET,this->indr+LU1-LUSOL_ARRAYOFFSET);
        }
      //===============================================================
      //       Restore the saved values of  iqloc.
      //       Insert the correct indices for the col of L and the row of U.
      //===============================================================
      x700:
        this->lenr[IBEST] = 0;
        this->lenc[JBEST] = 0;
        LL = LL1-1;
        LS = LSAVE;
        for(LC = LPIVC1; LC <= LPIVC2; LC++) 
        {
            LL++;
            LS++;
            I = this->indc[LC];
            this->iqloc[I] = this->indr[LS];
            this->indc[LL] = I;
            this->indr[LL] = IBEST;
        }
        LU = LU1-1;
        for(LR = LPIVR; LR <= LPIVR2; LR++) 
        {
            LU++;
            this->indr[LU] = this->indr[LR];
        }
        //===============================================================
        //       Free the space occupied by the pivot row
        //       and update the column permutation.
        //       Then free the space occupied by the pivot column
        //       and update the row permutation.
        //       nzchng is found in both calls to lu1pq2, but we use it only
        //       after the second.
        //===============================================================
        LU1PQ2(NCOLD, &NZCHNG, this->indr+LPIVR-LUSOL_ARRAYOFFSET,
               this->indc+LU1-LUSOL_ARRAYOFFSET, this->lenc,
               this->iqloc, this->iq, this->iqinv);
        LU1PQ2(NROWD, &NZCHNG, this->indc+LPIVC-LUSOL_ARRAYOFFSET,
               this->indc+LSAVE-LUSOL_ARRAYOFFSET, this->lenr,
               this->iploc, this->ip, this->ipinv);
        NZLEFT += NZCHNG;

        //===============================================================
        //       lu1mxr resets Amaxr(i) in each modified row i.
        //       lu1mxc moves the largest aij to the top of each modified col j.
        //       28 Jun 2002: Note that cols of L have an implicit diag of 1.0,
        //                    so lu1mxr is called with ll1, not ll1+1, whereas
        //                       lu1mxc is called with          lu1+1.
        //===============================================================
        if(UTRI && TPP) 
        {
            // Relax -- we're not keeping big elements at the top yet.
        }
        else 
        {
            if(TRP && MELIM>0)
			{
				MARK = MARK + 1;
                LU1MXR(MARK, LL1, LL, this->indc, this->m_cols, this->m_markc, 
                        this->m_markr, this->amaxr);
			}

            if(NELIM>0) 
            {
                LU1MXC(LU1+1,LU,this->indr);
                //  Update modified columns in heap
                if(TCP) 
                { 
                    for(KK = LU1+1; KK <= LU; KK++) 
                    {
                        J = this->indr[KK];
                        K = this->m_heap.get_Hk()[J];
                        // Biggest aij in column j
                        if (this->lenc[J] == 0)
                        {
                            V = 0.;
                        }
                        else
                        {
                            V = abs(this->a[this->locc[J]]);
                        };
                        this->m_heap.change(HLEN,K,V,J);
                    }
                }
            }
        }
        //==============================================================
        //       Negate lengths of pivot row and column so they will be
        //       eliminated during compressions.
        //===============================================================
        this->lenr[IBEST] = -NCOLD;
        this->lenc[JBEST] = -NROWD;

        // Test for fatal bug: row or column lists overwriting L and U. 
        if(LROW>LSAVE || LCOL>LSAVE)
            goto x980;

        // Reset the file lengths if pivot row or col was at the end.
        if(IBEST==ILAST)
            LROW = this->locr[IBEST];

        if(JBEST==JLAST)
            LCOL = this->locc[JBEST];
    }
    //------------------------------------------------------------------
    //    End of main loop.
    //------------------------------------------------------------------
    
    //------------------------------------------------------------------
    //    Normal exit.
    //    Move empty rows and cols to the end of ip, iq.
    //    Then finish with a dense LU if necessary.
    //------------------------------------------------------------------
  x900:
    *INFORM = LUSOL_INFORM_LUSUCCESS;
    LU1PQ3(this->m,this->lenr,this->ip,this->ipinv,&MRANK);
    LU1PQ3(this->n,this->lenc,this->iq,this->iqinv,NRANK);
    *NRANK = minimum(MRANK,*NRANK);
    if(DENSLU) 
    {
        #ifdef UseTimer
            LUSOL_timer (LUSOL, 17, "start" );
        #endif
        bool mem;
        mem = LU1FUL(LEND,LU1,LPIV,MLEFT,NLEFT,*NRANK,NROWU,LENL,LENU,
                &NSING,KEEPLU,SMALL,this->a+LD-LUSOL_ARRAYOFFSET,this->locr);
        if (mem)
        {
            goto x971;
        };

        *NRANK = MINMN-NSING;
        #ifdef UseTimer
            LUSOL_timer (LUSOL, 17, "finish" );
        #endif
    }
    *MINLEN = (*LENL)+(*LENU)+2*(this->m+this->n);
    goto x990;
  // Not enough space free after a compress.
  //    Set  minlen  to an estimate of the necessary value of  lena.
  x970:
    *INFORM = LUSOL_INFORM_ANEEDMEM;
    *MINLEN = LENA2+LFILE+2*(this->m+this->n);
    goto x990;
  x971:
    *INFORM = LUSOL_INFORM_ANEEDMEM;
    goto x990;
  // Fatal error.  This will never happen!
  // (Famous last words.)
  x980:
    *INFORM = LUSOL_INFORM_FATALERR;
    goto x990;

  // Fatal error with TSP.  Diagonal pivot not found.
  //x985:
  //  *INFORM = LUSOL_INFORM_NOPIVOT;

  //  Exit.
 x990:
    #ifdef UseTimer
        LUSOL_timer (LUSOL, 3 , "finish");
    #endif
    ;
};
/* ==================================================================
   lu1fac computes a factorization A = L*U, where A is a sparse
   matrix with m rows and n columns, P*L*P' is lower triangular
   and P*U*Q is upper triangular for certain permutations P, Q
   (which are returned in the arrays ip, iq).
   Stability is ensured by limiting the size of the elements of L.
   The nonzeros of A are input via the parallel arrays a, indc, indr,
   which should contain nelem entries of the form    aij,    i,    j
   in any order.  There should be no duplicate pairs         i,    j.

   ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   +        Beware !!!   The row indices i must be in indc,         +
   +              and the column indices j must be in indr.         +
   +              (Not the other way round!)                        +
   ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   It does not matter if some of the entries in a(*) are zero.
   Entries satisfying  abs( a(i) ) .le. parmlu(3)  are ignored.
   Other parameters in luparm and parmlu are described below.
   The matrix A may be singular.  On exit, nsing = luparm(11) gives
   the number of apparent singularities.  This is the number of
   "small" diagonals of the permuted factor U, as judged by
   the input tolerances Utol1 = parmlu(4) and  Utol2 = parmlu(5).
   The diagonal element diagj associated with column j of A is
   "small" if
               abs( diagj ) .le. Utol1
   or
               abs( diagj ) .le. Utol2 * max( uj ),
   where max( uj ) is the maximum element in the j-th column of U.
   The position of such elements is returned in w(*).  In general,
   w(j) = + max( uj ),  but if column j is a singularity,
   w(j) = - max( uj ).  Thus, w(j) .le. 0 if column j appears to be
   dependent on the other columns of A.
   NOTE: lu1fac (like certain other sparse LU packages) does not
   treat dense columns efficiently.  This means it will be slow
   on "arrow matrices" of the form
                A = (x       a)
                    (  x     b)
                    (    x   c)
                    (      x d)
                    (x x x x e)
   if the numerical values in the dense column allow it to be
   chosen LATE in the pivot order.
   With TPP (Threshold Partial Pivoting), the dense column is
   likely to be chosen late.
   With TCP (Threshold Complete Pivoting), if any of a,b,c,d
   is significantly larger than other elements of A, it will
   be chosen as the first pivot and the dense column will be
   eliminated, giving reasonably sparse factors.
   However, if element e is so big that TCP chooses it, the factors
   will become dense.  (It's hard to win on these examples!)
   ------------------------------------------------------------------

   Notes on the array names
   ------------------------
   During the LU factorization, the sparsity pattern of the matrix
   being factored is stored twice: in a column list and a row list.
   The column list is ( a, indc, locc, lenc )
   where
         a(*)    holds the nonzeros,
         indc(*) holds the indices for the column list,
         locc(j) points to the start of column j in a(*) and indc(*),
         lenc(j) is the number of nonzeros in column j.
   The row list is    (    indr, locr, lenr )
   where
         indr(*) holds the indices for the row list,
         locr(i) points to the start of row i in indr(*),
         lenr(i) is the number of nonzeros in row i.
   At all stages of the LU factorization, ip contains a complete
   row permutation.  At the start of stage k,  ip(1), ..., ip(k-1)
   are the first k-1 rows of the final row permutation P.
   The remaining rows are stored in an ordered list
                        ( ip, iploc, ipinv )
   where
         iploc(nz) points to the start in ip(*) of the set of rows
                   that currently contain nz nonzeros,
         ipinv(i)  points to the position of row i in ip(*).
   For example,
         iploc(1) = k   (and this is where rows of length 1 {),
         iploc(2) = k+p  if there are p rows of length 1
                        (and this is where rows of length 2 {).
   Similarly for iq, iqloc, iqinv.
   ---------------------------------------------------------------------
   INPUT PARAMETERS
   m      (not altered) is the number of rows in A.
   n      (not altered) is the number of columns in A.
   nelem  (not altered) is the number of matrix entries given in
          the arrays a, indc, indr.
   lena   (not altered) is the dimension of  a, indc, indr.
          This should be significantly larger than nelem.
          Typically one should have
             lena > max( 2*nelem, 10*m, 10*n, 10000 )
          but some applications may need more.
          On machines with virtual memory it is safe to have
          lena "far bigger than necessary", since not all of the
          arrays will be used.
   a      (overwritten) contains entries   Aij  in   a(1:nelem).
   indc   (overwritten) contains the indices i in indc(1:nelem).
   indr   (overwritten) contains the indices j in indr(1:nelem).
   luparm input parameters:                                Typical value
   luparm( 1) = nout     File number for printed messages.         6
   luparm( 2) = lprint   Print level.                              0
                    <  0 suppresses output.
                    =  0 gives error messages.
                   >= 10 gives statistics about the LU factors.
                   >= 50 gives debug output from lu1fac
                         (the pivot row and column and the
                         no. of rows and columns involved at
                         each elimination step).
   luparm( 3) = maxcol   lu1fac: maximum number of columns         5
                         searched allowed in a Markowitz-type
                         search for the next pivot element.
                         For some of the factorization, the
                         number of rows searched is
                         maxrow = maxcol - 1.
   luparm( 6) = 0    =>  TPP: Threshold Partial   Pivoting.        0
              = 1    =>  TRP: Threshold Rook      Pivoting.
              = 2    =>  TCP: Threshold Complete  Pivoting.
              = 3    =>  TSP: Threshold Symmetric Pivoting.
              = 4    =>  TDP: Threshold Diagonal  Pivoting.
                              (TDP not yet implemented).
                         TRP and TCP are more expensive than TPP but
                         more stable and better at revealing rank.
                         Take care with setting parmlu(1), especially
                         with TCP.
                         NOTE: TSP and TDP are for symmetric matrices
                         that are either definite or quasi-definite.
                         TSP is effectively TRP for symmetric matrices.
                         TDP is effectively TCP for symmetric matrices.
   luparm( 8) = keepLU   lu1fac: keepLU = 1 means the numerical    1
                         factors will be computed if possible.
                         keepLU = 0 means L and U will be discarded
                         but other information such as the row and
                         column permutations will be returned.
                         The latter option requires less storage.
   parmlu input parameters:                                Typical value
   parmlu( 1) = Ltol1    Max Lij allowed during Factor.
                                                   TPP     10.0 or 100.0
                                                   TRP      4.0 or  10.0
                                                   TCP      5.0 or  10.0
                                                   TSP      4.0 or  10.0
                         With TRP and TCP (Rook and Complete Pivoting),
                         values less than 25.0 may be expensive
                         on badly scaled data.  However,
                         values less than 10.0 may be needed
                         to obtain a reliable rank-revealing
                         factorization.
   parmlu( 2) = Ltol2    Max Lij allowed during Updates.            10.0
                         during updates.
   parmlu( 3) = small    Absolute tolerance for       eps**0.8 = 3.0d-13
                         treating reals as zero.
   parmlu( 4) = Utol1    Absolute tol for flagging    eps**0.67= 3.7d-11
                         small diagonals of U.
   parmlu( 5) = Utol2    Relative tol for flagging    eps**0.67= 3.7d-11
                         small diagonals of U.
                         (eps = machine precision)
   parmlu( 6) = Uspace   Factor limiting waste space in  U.      3.0
                         In lu1fac, the row or column lists
                         are compressed if their length
                         exceeds Uspace times the length of
                         either file after the last compression.
   parmlu( 7) = dens1    The density at which the Markowitz      0.3
                         pivot strategy should search maxcol
                         columns and no rows.
                         (Use 0.3 unless you are experimenting
                         with the pivot strategy.)
   parmlu( 8) = dens2    the density at which the Markowitz      0.5
                         strategy should search only 1 column,
                         or (if storage is available)
                         the density at which all remaining
                         rows and columns will be processed
                         by a dense LU code.
                         For example, if dens2 = 0.1 and lena is
                         large enough, a dense LU will be used
                         once more than 10 per cent of the
                         remaining matrix is nonzero.

   OUTPUT PARAMETERS
   a, indc, indr     contain the nonzero entries in the LU factors of A.
          If keepLU = 1, they are in a form suitable for use
          by other parts of the LUSOL package, such as lu6sol.
          U is stored by rows at the start of a, indr.
          L is stored by cols at the end   of a, indc.
          If keepLU = 0, only the diagonals of U are stored, at the
          end of a.
   ip, iq    are the row and column permutations defining the
          pivot order.  For example, row ip(1) and column iq(1)
          defines the first diagonal of U.
   lenc(1:numl0) contains the number of entries in nontrivial
          columns of L (in pivot order).
   lenr(1:m) contains the number of entries in each row of U
          (in original order).
   locc(1:n) = 0 (ready for the LU update routines).
   locr(1:m) points to the beginning of the rows of U in a, indr.
   iploc, iqloc, ipinv, iqinv  are undefined.
   w      indicates singularity as described above.
   inform = 0 if the LU factors were obtained successfully.
          = 1 if U appears to be singular, as judged by lu6chk.
          = 3 if some index pair indc(l), indr(l) lies outside
              the matrix dimensions 1:m , 1:n.
          = 4 if some index pair indc(l), indr(l) duplicates
              another such pair.
          = 7 if the arrays a, indc, indr were not large enough.
              Their length "lena" should be increase to at least
              the value "minlen" given in luparm(13).
          = 8 if there was some other fatal error.  (Shouldn't happen!)
          = 9 if no diagonal pivot could be found with TSP or TDP.
              The matrix must not be sufficiently definite
              or quasi-definite.
   luparm output parameters:
   luparm(10) = inform   Return code from last call to any LU routine.
   luparm(11) = nsing    No. of singularities marked in the
                         output array w(*).
   luparm(12) = jsing    Column index of last singularity.
   luparm(13) = minlen   Minimum recommended value for  lena.
   luparm(14) = maxlen   ?
   luparm(15) = nupdat   No. of updates performed by the lu8 routines.
   luparm(16) = nrank    No. of nonempty rows of U.
   luparm(17) = ndens1   No. of columns remaining when the density of
                         the matrix being factorized reached dens1.
   luparm(18) = ndens2   No. of columns remaining when the density of
                         the matrix being factorized reached dens2.
   luparm(19) = jumin    The column index associated with DUmin.
   luparm(20) = numL0    No. of columns in initial  L.
   luparm(21) = lenL0    Size of initial  L  (no. of nonzeros).
   luparm(22) = lenU0    Size of initial  U.
   luparm(23) = lenL     Size of current  L.
   luparm(24) = lenU     Size of current  U.
   luparm(25) = lrow     Length of row file.
   luparm(26) = ncp      No. of compressions of LU data structures.
   luparm(27) = mersum   lu1fac: sum of Markowitz merit counts.
   luparm(28) = nUtri    lu1fac: triangular rows in U.
   luparm(29) = nLtri    lu1fac: triangular rows in L.
   luparm(30) =
   parmlu output parameters:
   parmlu(10) = Amax     Maximum element in  A.
   parmlu(11) = Lmax     Maximum multiplier in current  L.
   parmlu(12) = Umax     Maximum element in current  U.
   parmlu(13) = DUmax    Maximum diagonal in  U.
   parmlu(14) = DUmin    Minimum diagonal in  U.
   parmlu(15) = Akmax    Maximum element generated at any stage
                         during TCP factorization.
   parmlu(16) = growth   TPP: Umax/Amax    TRP, TCP, TSP: Akmax/Amax
   parmlu(17) =
   parmlu(18) =
   parmlu(19) =
   parmlu(20) = resid    lu6sol: residual after solve with U or U'.
   ...
   parmlu(30) =
   ------------------------------------------------------------------
   00 Jun 1983  Original version.
   00 Jul 1987  nrank  saved in luparm(16).
   12 Apr 1989  ipinv, iqinv added as workspace.
   26 Apr 1989  maxtie replaced by maxcol in Markowitz search.
   16 Mar 1992  jumin  saved in luparm(19).
   10 Jun 1992  lu1fad has to move empty rows and cols to the bottom
                (via lu1pq3) before doing the dense LU.
   12 Jun 1992  Deleted dense LU (lu1ful, lu1vlu).
   25 Oct 1993  keepLU implemented.
   07 Feb 1994  Added new dense LU (lu1ful, lu1den).
   21 Dec 1994  Bugs fixed in lu1fad (nrank) and lu1ful (ipvt).
   08 Aug 1995  Use ip instead of w as parameter to lu1or3 (for F90).
   13 Sep 2000  TPP and TCP options implemented.
   17 Oct 2000  Fixed troubles due to A = empty matrix (Todd Munson).
   01 Dec 2000  Save Lmax, Umax, etc. after both lu1fad and lu6chk.
                lu1fad sets them when keepLU = false.
                lu6chk sets them otherwise, and includes items
                from the dense LU.
   11 Mar 2001  lu6chk now looks at diag(U) when keepLU = false.
   26 Apr 2002  New TCP implementation using heap routines to
                store largest element in each column.
                New workspace arrays Ha, Hj, Hk required.
                For compatibility, borrow space from a, indc, indr
                rather than adding new input parameters.
   01 May 2002  lu1den changed to lu1DPP (dense partial  pivoting).
                lu1DCP implemented       (dense complete pivoting).
                Both TPP and TCP now switch to dense mode and end.
   ================================================================== */
template<class T>
void LUSOLrec<T>::LU1FAC(INT *INFORM)
{
    bool  KEEPLU, TCP, TPP, TRP, TSP;
    INT   LPIV, NELEM0, LPRINT, MINLEN, NUML0, LENL, LENU, LROW, MERSUM,
          NUTRI, NLTRI, NDENS1, NDENS2, NRANK, NSING, JSING, JUMIN, NUMNZ, LERR,
          LU, LL, LM, LTOPL, K, I, LENUK, J, LENLK, IDUMMY, LLSAVE, NMOVE, L2, L, NCP, NBUMP;

    TR    LMAX, LTOL, SMALL, AMAX, UMAX, DUMAX, DUMIN, AKMAX, DM, DN, DELEM, DENSTY,
          AGRWTH, UGRWTH, GROWTH, CONDU, DINCR, AVGMER;

    // Grab relevant input parameters.
    NELEM0 = this->nelem;
    LPRINT = this->luparm[LUSOL_IP_PRINTLEVEL];
    LPIV   = this->luparm[LUSOL_IP_PIVOTTYPE];
    KEEPLU = (bool) (this->luparm[LUSOL_IP_KEEPLU]!=false);
    // Limit on size of Lij
    LTOL   = this->parmlu[LUSOL_RP_FACTORMAX_Lij];
    // Drop tolerance
    SMALL  = this->parmlu[LUSOL_RP_ZEROTOLERANCE];
    TPP = (bool) (LPIV==LUSOL_PIVOT_TPP);
    TRP = (bool) (LPIV==LUSOL_PIVOT_TRP);
    TCP = (bool) (LPIV==LUSOL_PIVOT_TCP);
    TSP = (bool) (LPIV==LUSOL_PIVOT_TSP);
    // Initialize output parameters.
    *INFORM = LUSOL_INFORM_LUSUCCESS;
    LERR   = 0;
    MINLEN = this->nelem + 2*(this->m+this->n);
    NUML0  = 0;
    LENL   = 0;
    LENU   = 0;
    LROW   = 0;
    MERSUM = 0;
    NUTRI  = this->m;
    NLTRI  = 0;
    NDENS1 = 0;
    NDENS2 = 0;
    NRANK  = 0;
    NSING  = 0;
    JSING  = 0;
    JUMIN  = 0;
    AMAX   = 0;
    LMAX   = 0;
    UMAX   = 0;
    DUMAX  = 0;
    DUMIN  = 0;
    AKMAX  = 0;

    // Float version of dimensions.
    DM = TR(this->m);
    DN = TR(this->n);
    DELEM = TR(this->nelem);

    // Initialize workspace parameters.
    this->luparm[LUSOL_IP_COMPRESSIONS_LU] = 0;
    if(this->lena < MINLEN) 
    {
        if(!LUSOL_realloc_a(MINLEN))
            goto x970;
    }

    //------------------------------------------------------------------
    //    Organize the  aij's  in  a, indc, indr.
    //    lu1or1  deletes small entries, tests for illegal  i,j's,
    //            and counts the nonzeros in each row and column.
    //    lu1or2  reorders the elements of  A  by columns.
    //    lu1or3  uses the column list to test for duplicate entries
    //            (same indices  i,j).
    //    lu1or4  constructs a row list from the column list.
    //------------------------------------------------------------------
    LU1OR1(SMALL,&AMAX,&NUMNZ,&LERR,INFORM);
    if(LPRINT>=LUSOL_MSG_STATISTICS) 
    {
        DENSTY = (100*DELEM)/(DM*DN);
        LUSOL_report(this, 0, (char*) "m:%6d %c n:%6d  nzcount:%9d  Amax:%g  Density:%g\n",
                           this->m, relationChar(this->m, this->n), this->n,
                           this->nelem, AMAX, DENSTY); //FIXME:
    }
    if(*INFORM!=LUSOL_INFORM_LUSUCCESS)
        goto x930;
    this->nelem = NUMNZ;
    LU1OR2();
    LU1OR3(&LERR,INFORM);
    if(*INFORM!=LUSOL_INFORM_LUSUCCESS)
        goto x940;
    LU1OR4();
    //------------------------------------------------------------------
    //    Set up lists of rows and columns with equal numbers of nonzeros,
    //    using  indc(*)  as workspace.
    //------------------------------------------------------------------
    LU1PQ1(this->m,this->n,this->lenr, this->ip,this->iploc,this->ipinv,
         this->indc+this->nelem); // LUSOL_ARRAYOFFSET implied
    LU1PQ1(this->n,this->m,this->lenc, this->iq,this->iqloc,this->iqinv,
         this->indc+this->nelem); // LUSOL_ARRAYOFFSET implied
    //------------------------------------------------------------------
    //    For TCP, Ha, Hj, Hk are allocated separately, similarly amaxr
    //    for TRP. Then compute the factorization  A = L*U.
    //------------------------------------------------------------------

    LU1FAD(INFORM,&LENL,&LENU,&MINLEN,&MERSUM,&NUTRI,&NLTRI,&NDENS1,&NDENS2,
         &NRANK,&LMAX,&UMAX,&DUMAX,&DUMIN,&AKMAX);

    this->luparm[LUSOL_IP_RANK_U]     = NRANK;
    this->luparm[LUSOL_IP_NONZEROS_L] = LENL;
    if(*INFORM==LUSOL_INFORM_ANEEDMEM)
        goto x970;
    if(*INFORM==LUSOL_INFORM_NOPIVOT)
        goto x985;
    if(*INFORM>LUSOL_INFORM_LUSUCCESS)
        goto x980;
    if(KEEPLU) 
    {
        //---------------------------------------------------------------
        //   The LU factors are at the top of  a, indc, indr,
        //   with the columns of  L  and the rows of  U  in the order
        //   ( free )   ... ( u3 ) ( l3 ) ( u2 ) ( l2 ) ( u1 ) ( l1 ).
        //   Starting with ( l1 ) and ( u1 ), move the rows of  U  to the
        //   left and the columns of  L  to the right, giving
        //   ( u1 ) ( u2 ) ( u3 ) ...   ( free )   ... ( l3 ) ( l2 ) ( l1 ).
        //   Also, set  numl0 = the number of nonempty columns of  U.
        //---------------------------------------------------------------
        LU = 0;
        LL = this->lena+1;
        LM = LL;
    
        LTOPL = LL-LENL-LENU;
        LROW = LENU;
        for(K = 1; K <= NRANK; K++) 
        {
            I = this->ip[K];
            LENUK = -this->lenr[I];
            this->lenr[I] = LENUK;
            J = this->iq[K];
            LENLK = -this->lenc[J]-1;
            if(LENLK>0) 
            {
                NUML0++;
                this->iqloc[NUML0] = LENLK;
            }
            if(LU+LENUK<LTOPL) 
            {
                //=========================================================
                // There is room to move ( uk ).  Just right-shift ( lk ).
                //=========================================================
                for(IDUMMY = 1; IDUMMY <= LENLK; IDUMMY++) 
                {
                    LL--;
                    LM--;
                    this->a[LL] = this->a[LM];
                    this->indc[LL] = this->indc[LM];
                    this->indr[LL] = this->indr[LM];
                }
            }
            else 
            {
                //=========================================================
                // There is no room for ( uk ) yet.  We have to
                // right-shift the whole of the remaining LU file.
                // Note that ( lk ) ends up in the correct place.
                //=========================================================
                LLSAVE = LL-LENLK;
                NMOVE = LM-LTOPL;
                for(IDUMMY = 1; IDUMMY <= NMOVE; IDUMMY++) 
                {
                    LL--;
                    LM--;
                    this->a[LL] = this->a[LM];
                    this->indc[LL] = this->indc[LM];
                    this->indr[LL] = this->indr[LM];
                }
                LTOPL = LL;
                LL = LLSAVE;
                LM = LL;
            }
            //======================================================
            //  Left-shift ( uk ).
            //======================================================
            this->locr[I] = LU+1;
            L2 = LM-1;
            LM = LM-LENUK;
            for(L = LM; L <= L2; L++) 
            {
                LU = LU+1;
                this->a[LU] = this->a[L];
                this->indr[LU] = this->indr[L];
            }
        }
        //---------------------------------------------------------------
        //   Save the lengths of the nonempty columns of  L,
        //   and initialize  locc(j)  for the LU update routines.
        //---------------------------------------------------------------
        for(K = 1; K <= NUML0; K++) 
        {
            this->lenc[K] = this->iqloc[K];
        }
        for(J = 1; J <= this->n; J++) 
        {
            this->locc[J] = 0;
        }
        //---------------------------------------------------------------
        //   Test for singularity.
        //   lu6chk  sets  nsing, jsing, jumin, Lmax, Umax, DUmax, DUmin
        //   (including entries from the dense LU).
        //   inform = 1  if there are singularities (nsing gt 0).
        //--------------------------------------------------------------
        LU6CHK(1,this->lena,INFORM);
        NSING = this->luparm[LUSOL_IP_SINGULARITIES];
        JSING = this->luparm[LUSOL_IP_SINGULARINDEX];
        JUMIN = this->luparm[LUSOL_IP_COLINDEX_DUMIN];
        LMAX  = this->parmlu[LUSOL_RP_MAXMULT_L];
        UMAX  = this->parmlu[LUSOL_RP_MAXELEM_U];
        DUMAX = this->parmlu[LUSOL_RP_MAXELEM_DIAGU];
        DUMIN = this->parmlu[LUSOL_RP_MINELEM_DIAGU];
    }
    else 
    {
        //---------------------------------------------------------------
        //   keepLU = 0.  L and U were not kept, just the diagonals of U.
        //   lu1fac will probably be called again soon with keepLU = .true.
        //   11 Mar 2001: lu6chk revised.  We can call it with keepLU = 0,
        //                but we want to keep Lmax, Umax from lu1fad.
        //   05 May 2002: Allow for TCP with new lu1DCP.  Diag(U) starts
        //                below lena2, not lena.  Need lena2 in next line.
        //---------------------------------------------------------------
        LU6CHK(1,this->lena,INFORM);
        NSING = this->luparm[LUSOL_IP_SINGULARITIES];
        JSING = this->luparm[LUSOL_IP_SINGULARINDEX];
        JUMIN = this->luparm[LUSOL_IP_COLINDEX_DUMIN];
        DUMAX = this->parmlu[LUSOL_RP_MAXELEM_DIAGU];
        DUMIN = this->parmlu[LUSOL_RP_MINELEM_DIAGU];
    }
    goto x990;
    //------------
    //    Error exits.
    //------------
  x930:
    *INFORM = LUSOL_INFORM_ADIMERR;
    if(LPRINT>=LUSOL_MSG_SINGULARITY)
        LUSOL_report(this, 0, (char*) "lu1fac  error...\nentry  a[%d]  has an illegal row (%d) or column (%d) index\n",
                        LERR,this->indc[LERR],this->indr[LERR]); //FIXME:
    goto x990;
  x940:
    *INFORM = LUSOL_INFORM_ADUPLICATE;
    if(LPRINT>=LUSOL_MSG_SINGULARITY)
        LUSOL_report(this, 0, (char*) "lu1fac  error...\nentry  a[%d]  is a duplicate with indeces indc=%d, indr=%d\n",
                        LERR,this->indc[LERR],this->indr[LERR]); //FIXME:
    goto x990;
  x970:
    *INFORM = LUSOL_INFORM_ANEEDMEM;
    if(LPRINT>=LUSOL_MSG_SINGULARITY)
    {
        LUSOL_report(this, 0, (char*) "lu1fac  error...\ninsufficient storage\n"); //FIXME:
        //LUSOL_report(this, 0, "lu1fac  error...\ninsufficient storage; increase  lena  from %d to at least %d\n",
        //                this->lena, MINLEN);
    };
    goto x990;
  x980:
    *INFORM = LUSOL_INFORM_FATALERR;
    if(LPRINT>=LUSOL_MSG_SINGULARITY)
        LUSOL_report(this, 0, (char*) "lu1fac  error...\nfatal bug   (sorry --- this should never happen)\n"); //FIXME:
    goto x990;
  x985:
    *INFORM = LUSOL_INFORM_NOPIVOT;
    if(LPRINT>=LUSOL_MSG_SINGULARITY)
        LUSOL_report(this, 0, (char*) "lu1fac  error...\nTSP used but diagonal pivot could not be found\n"); //FIXME:
  // Store output parameters.
  x990:
    this->nelem = NELEM0;
    this->luparm[LUSOL_IP_INFORM]          = *INFORM;
    this->luparm[LUSOL_IP_SINGULARITIES]   = NSING;
    this->luparm[LUSOL_IP_SINGULARINDEX]   = JSING;
    this->luparm[LUSOL_IP_MINIMUMLENA]     = MINLEN;
    this->luparm[LUSOL_IP_UPDATECOUNT]     = 0;
    this->luparm[LUSOL_IP_RANK_U]          = NRANK;
    this->luparm[LUSOL_IP_COLCOUNT_DENSE1] = NDENS1;
    this->luparm[LUSOL_IP_COLCOUNT_DENSE2] = NDENS2;
    this->luparm[LUSOL_IP_COLINDEX_DUMIN]  = JUMIN;
    this->luparm[LUSOL_IP_COLCOUNT_L0]     = NUML0;
    this->luparm[LUSOL_IP_NONZEROS_L0]     = LENL;
    this->luparm[LUSOL_IP_NONZEROS_U0]     = LENU;
    this->luparm[LUSOL_IP_NONZEROS_L]      = LENL;
    this->luparm[LUSOL_IP_NONZEROS_U]      = LENU;
    this->luparm[LUSOL_IP_NONZEROS_ROW]    = LROW;
    this->luparm[LUSOL_IP_MARKOWITZ_MERIT] = MERSUM;
    this->luparm[LUSOL_IP_TRIANGROWS_U]    = NUTRI;
    this->luparm[LUSOL_IP_TRIANGROWS_L]    = NLTRI;
    this->parmlu[LUSOL_RP_MAXELEM_A]       = AMAX;
    this->parmlu[LUSOL_RP_MAXMULT_L]       = LMAX;
    this->parmlu[LUSOL_RP_MAXELEM_U]       = UMAX;
    this->parmlu[LUSOL_RP_MAXELEM_DIAGU]   = DUMAX;
    this->parmlu[LUSOL_RP_MINELEM_DIAGU]   = DUMIN;
    this->parmlu[LUSOL_RP_MAXELEM_TCP]     = AKMAX;
    AGRWTH = AKMAX/(AMAX+lusol_user_params<TR>::LUSOL_SMALLNUM);
    UGRWTH = UMAX/(AMAX+lusol_user_params<TR>::LUSOL_SMALLNUM);
    if(TPP)
        GROWTH = UGRWTH;
        // TRP or TCP or TSP
    else
        GROWTH = AGRWTH;
    this->parmlu[LUSOL_RP_GROWTHRATE] = GROWTH;

    #ifdef UseRowBasedL0
        // Create a row version L0
        LU1L0T(LUSOL, &(this->L0));
        if (this->L0 == nullptr)
        {
            *INFORM = LUSOL_INFORM_ANEEDMEM;
        };
    #endif

    //------------------------------------------------------------------
    //    Print statistics for the LU factors.
    //------------------------------------------------------------------
    NCP   = this->luparm[LUSOL_IP_COMPRESSIONS_LU];
    CONDU = DUMAX/maximum(DUMIN,lusol_user_params<TR>::LUSOL_SMALLNUM);
    DINCR = TR((LENL+LENU)-this->nelem);
    DINCR = (DINCR*100)/maximum(DELEM,TR(1.));
    AVGMER = TR(MERSUM);
    AVGMER = AVGMER/DM;
    NBUMP = this->m-NUTRI-NLTRI;
    if(LPRINT>=LUSOL_MSG_STATISTICS) 
    {
        if(TPP) 
        {
            LUSOL_report(this, 0, (char*) "Merit %g %d %d %d %g %d %d %g %g %d %d %d\n",
                          AVGMER,LENL,LENL+LENU,NCP,DINCR,NUTRI,LENU,
                          LTOL,UMAX,UGRWTH,NLTRI,NDENS1,LMAX); //FIXME
        }
        else 
        {
            LUSOL_report(this, 0, (char*) "Merit %s %g %d %d %d %g %d %d %g %g %d %d %d %g %g\n",
                          LUSOL_pivotLabel(),
                          AVGMER,LENL,LENL+LENU,NCP,DINCR,NUTRI,LENU,
                          LTOL,UMAX,UGRWTH,NLTRI,NDENS1,LMAX,AKMAX,AGRWTH); //FIXME:
        }
        LUSOL_report(this, 0, (char*) "bump%9d  dense2%7d  DUmax%g DUmin%g  conDU%g\n",
                          NBUMP,NDENS2,DUMAX,DUMIN,CONDU); //FIXME:
    }
};
};
