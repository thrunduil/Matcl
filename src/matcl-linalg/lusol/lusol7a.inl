#include "lusol.h"
#include <math.h>

/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   File  lusol7a
      lu7add   lu7cyc   lu7elm   lu7for   lu7rnk   lu7zap
      Utilities for LUSOL's update routines.
      lu7for is the most important -- the forward sweep.
  01 May 2002: Derived from LUSOL's original lu7a.f file.
  01 May 2002: Current version of lusol7a.f.
   ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */

/* ==================================================================
   lu7add  inserts the first nrank elements of the vector v(*)
   as column  jadd  of  U.  We assume that  U  does not yet have any
   entries in this column.
   Elements no larger than  parmlu(3)  are treated as zero.
   klast  will be set so that the last row to be affected
   (in pivotal order) is row  ip(klast).
   ------------------------------------------------------------------
   09 May 1988: First f77 version.
   ================================================================== */
namespace lusol
{
template<class T>
void LUSOLrec<T>::LU7ADD(INT JADD, T V[], INT LENL, INT *LENU,
    INT *LROW, INT NRANK, INT *INFORM, INT *KLAST, TR *VNORM)
{
    TR SMALL;
    INT  K, I, LENI, MINFRE, NFREE, LR1, LR2, L;

    SMALL = this->parmlu[LUSOL_RP_ZEROTOLERANCE];
    *VNORM = 0;
    *KLAST = 0;
    for(K = 1; K <= NRANK; K++) 
    {
        I = this->ip[K];
        if(abs(V[I])<=SMALL)
            continue;
        *KLAST = K;
        (*VNORM) += abs(V[I]);
        LENI = this->lenr[I];
        // Compress row file if necessary.
        MINFRE = LENI+1;
        NFREE = this->lena - LENL - *LROW;
        if(NFREE<MINFRE) 
        {
            LU1REC(this->m, true,LROW,this->indr,this->lenr,this->locr);
            NFREE = this->lena - LENL - *LROW;
            if(NFREE<MINFRE) 
            {
                goto x970;
            }
        }
        //   Move row  i  to the end of the row file,
        //   unless it is already there.
        //   No need to move if there is a gap already.
        if(LENI==0)
            this->locr[I] = (*LROW) + 1;
        LR1 = this->locr[I];
        LR2 = (LR1+LENI)-1;
        if(LR2==*LROW)
            goto x150;
        if(this->indr[LR2+1]==0)
            goto x180;
        this->locr[I] = (*LROW) + 1;

        L = LR2-LR1+1;
        if(L > 0) 
        {
            LR2 = (*LROW)+1;
            mem_move(this->a+LR2,    this->a+LR1, L);
            mem_move(this->indr+LR2, this->indr+LR1, L);
            memclear(this->indr+LR1, L);
            *LROW += L;
        }

      x150:
        LR2 = *LROW;
        (*LROW)++;

      // Add the element of  v.
      x180:
        LR2++;
        this->a[LR2] = V[I];
        this->indr[LR2] = JADD;
        this->lenr[I] = LENI+1;
        (*LENU)++;
    }
    //  Normal exit.
    *INFORM = LUSOL_INFORM_LUSUCCESS;
    goto x990;

  //  Not enough storage.
  x970:
    *INFORM = LUSOL_INFORM_ANEEDMEM;
  x990:
    ;
};

/* ==================================================================
   lu7cyc performs a cyclic permutation on the row or column ordering
   stored in ip, moving entry kfirst down to klast.
   If kfirst .ge. klast, lu7cyc should not be called.
   Sometimes klast = 0 and nothing should happen.
   ------------------------------------------------------------------
   09 May 1988: First f77 version.
   ================================================================== */
template<class T>
void LUSOLrec<T>::LU7CYC(INT KFIRST, INT KLAST, INT IX[])
{
    if(KFIRST<KLAST) 
    {
        INT IFIRST, K;

        INT *IXK, *IXK1;
        IXK = IX+KFIRST;
        IFIRST = *IXK;
        for(K = KFIRST, IXK1 = IXK+1; K <= KLAST-1; K++, IXK++, IXK1++) 
        {
            *IXK = *IXK1;
        }
        *IXK = IFIRST;
    }
};

/* ==================================================================
   lu7elm  eliminates the subdiagonal elements of a vector  v(*),
   where  L*v = y  for some vector y.
   If  jelm > 0,  y  has just become column  jelm  of the matrix  A.
   lu7elm  should not be called unless  m  is greater than  nrank.
   inform = 0 if y contained no subdiagonal nonzeros to eliminate.
   inform = 1 if y contained at least one nontrivial subdiagonal.
   inform = 7 if there is insufficient storage.
   ------------------------------------------------------------------
   09 May 1988: First f77 version.
                No longer calls lu7for at end.  lu8rpc, lu8mod do so.
   ================================================================== */
template<class T>
void LUSOLrec<T>::LU7ELM(INT JELM, T V[], INT *LENL,
            INT *LROW, INT NRANK, INT *INFORM, T *DIAG)
{
    TR   VI, SMALL;
    T    VMAX;
    INT  NRANK1, MINFRE, NFREE, KMAX, L, K, I, LMAX, IMAX, L1, L2;

    LMAX = 0;

    SMALL = this->parmlu[LUSOL_RP_ZEROTOLERANCE];
    NRANK1 = NRANK+1;
    *DIAG = 0;

    TR VMAX2 = 0;

    //  Compress row file if necessary.
    MINFRE = this->m-NRANK;
    NFREE = this->lena-(*LENL)-(*LROW);
    if(NFREE>=MINFRE)
        goto x100;
    LU1REC(this->m,true,LROW,this->indr,this->lenr,this->locr);
    NFREE = this->lena-(*LENL)-(*LROW);
    if(NFREE<MINFRE) 
    {
        goto x970;
    }
    // Pack the subdiagonals of  v  into  L,  and find the largest.
  x100:
    VMAX2 = 0;
    KMAX = 0;
    L = (this->lena-(*LENL))+1;
    for(K = NRANK1; K <= this->m; K++) 
    {
        I = this->ip[K];
        VI = abs(V[I]);
        if(VI<=SMALL)
            continue;
        L--;
        this->a[L] = V[I];
        this->indc[L] = I;
        if(VMAX2>=VI)
            continue;
        VMAX2 = VI;
        KMAX = K;
        LMAX = L;
    }
    if(KMAX==0)
        goto x900;
    //------------------------------------------------------------------
    //    Remove  vmax  by overwriting it with the last packed  v(i).
    //    Then set the multipliers in  L  for the other elements.
    //------------------------------------------------------------------
    IMAX = this->ip[KMAX];
    VMAX = this->a[LMAX];
    this->a[LMAX] = this->a[L];
    this->indc[LMAX] = this->indc[L];
    L1 = L+1;
    L2 = this->lena-(*LENL);
    *LENL = ((*LENL)+L2)-L;
    for(L = L1; L <= L2; L++) 
    {
        this->a[L] /= -VMAX;
        this->indr[L] = IMAX;
    }
    // Move the row containing vmax to pivotal position nrank + 1.
    this->ip[KMAX] = this->ip[NRANK1];
    this->ip[NRANK1] = IMAX;
    *DIAG = VMAX;
    //------------------------------------------------------------------
    //    If jelm is positive, insert  vmax  into a new row of  U.
    //    This is now the only subdiagonal element.
    //------------------------------------------------------------------
    if(JELM>0) 
    {
        (*LROW)++;
        this->locr[IMAX] = *LROW;
        this->lenr[IMAX] = 1;
        this->a[*LROW] = VMAX;
        this->indr[*LROW] = JELM;
    }
    *INFORM = LUSOL_INFORM_LUSINGULAR;
    goto x990;
    //  No elements to eliminate.
  x900:
    *INFORM = LUSOL_INFORM_LUSUCCESS;
    goto x990;
  //  Not enough storage.
  x970:
    *INFORM = LUSOL_INFORM_ANEEDMEM;
  x990:
    ;
};

/* ==================================================================
   lu7for  (forward sweep) updates the LU factorization  A = L*U
   when row  iw = ip(klast)  of  U  is eliminated by a forward
   sweep of stabilized row operations, leaving  ip * U * iq  upper
   triangular.
   The row permutation  ip  is updated to preserve stability and/or
   sparsity.  The column permutation  iq  is not altered.
   kfirst  is such that row  ip(kfirst)  is the first row involved
   in eliminating row  iw.  (Hence,  kfirst  marks the first nonzero
   in row  iw  in pivotal order.)  If  kfirst  is unknown it may be
   input as  1.
   klast   is such that row  ip(klast)  is the row being eliminated.
   klast   is not altered.
   lu7for  should be called only if  kfirst .le. klast.
   If  kfirst = klast,  there are no nonzeros to eliminate, but the
   diagonal element of row  ip(klast)  may need to be moved to the
   front of the row.
   ------------------------------------------------------------------
   On entry,  locc(*)  must be zero.

   On exit:
   inform = 0  if row iw has a nonzero diagonal (could be small).
   inform = 1  if row iw has no diagonal.
   inform = 7  if there is not enough storage to finish the update.

   On a successful exit (inform le 1),  locc(*)  will again be zero.
   ------------------------------------------------------------------
      Jan 1985: Final f66 version.
   09 May 1988: First f77 version.
   ================================================================== */
template<class T>
void LUSOLrec<T>::LU7FOR(INT KFIRST, INT KLAST, INT *LENL, INT *LENU,
                     INT *LROW, INT *INFORM, T *DIAG)
{
    bool SWAPPD;
    INT    KBEGIN, IW, LENW, LW1, LW2, JFIRST, MINFRE, NFREE, L, J, KSTART, KSTOP, K,
           LFIRST, IV, LENV, LV1, JLAST, LV2, LV3, LV, JV, LW, LDIAG, LIMIT;
    TR     LTOL, USPACE, SMALL;
    T      VJ, WJ, AMULT;

    LTOL   = this->parmlu[LUSOL_RP_UPDATEMAX_Lij];
    SMALL  = this->parmlu[LUSOL_RP_ZEROTOLERANCE];
    USPACE = this->parmlu[LUSOL_RP_COMPSPACE_U];
    KBEGIN = KFIRST;
    SWAPPD = false;

    // We come back here from below if a row interchange is performed.
  x100:
    IW = this->ip[KLAST];
    LENW = this->lenr[IW];
    if(LENW==0)
        goto x910;
    LW1 = this->locr[IW];
    LW2 = (LW1+LENW)-1;
    JFIRST = this->iq[KBEGIN];
    if(KBEGIN>=KLAST)
        goto x700;
    //  Make sure there is room at the end of the row file
    //  in case row  iw  is moved there and fills in completely.
    MINFRE = this->n+1;
    NFREE = this->lena-(*LENL)-(*LROW);
    if(NFREE<MINFRE) 
    {
        LU1REC(this->m,true,LROW,this->indr,this->lenr,this->locr);
        LW1 = this->locr[IW];
        LW2 = (LW1+LENW)-1;
        NFREE = this->lena-(*LENL)-(*LROW);
        if(NFREE<MINFRE) 
        {
            goto x970;
        }
    }
    // Set markers on row  iw.
    for(L = LW1; L <= LW2; L++) 
    {
        J = this->indr[L];
        this->locc[J] = L;
    }
    //==================================================================
    //    Main elimination loop.
    //==================================================================
    KSTART = KBEGIN;
    KSTOP = minimum(KLAST,this->n);
    for(K = KSTART; K <= KSTOP; K++) 
    {
        JFIRST = this->iq[K];
        LFIRST = this->locc[JFIRST];
        if(LFIRST==0)
            goto x490;
        // Row  iw  has its first element in column  jfirst.
        WJ = this->a[LFIRST];
        if(K==KLAST)
            goto x490;
        //---------------------------------------------------------------
        //   We are about to use the first element of row  iv
        //          to eliminate the first element of row  iw.
        //   However, we may wish to interchange the rows instead,
        //   to preserve stability and/or sparsity.
        //---------------------------------------------------------------
        IV = this->ip[K];
        LENV = this->lenr[IV];
        LV1 = this->locr[IV];
        VJ = 0;
        if(LENV==0)
            goto x150;
        if(this->indr[LV1]!=JFIRST)
            goto x150;
        VJ = this->a[LV1];
        if(SWAPPD)
            goto x200;
        if(LTOL*abs(WJ)<abs(VJ))
            goto x200;
        if(LTOL*abs(VJ)<abs(WJ))
            goto x150;
        if(LENV<=LENW)
            goto x200;
        //---------------------------------------------------------------
        //   Interchange rows  iv  and  iw.
        //---------------------------------------------------------------
      x150:
        this->ip[KLAST] = IV;
        this->ip[K] = IW;
        KBEGIN = K;
        SWAPPD = true;
        goto x600;
      //---------------------------------------------------------------
      //   Delete the eliminated element from row  iw
      //   by overwriting it with the last element.
      //---------------------------------------------------------------
      x200:
        this->a[LFIRST] = this->a[LW2];
        JLAST = this->indr[LW2];
        this->indr[LFIRST] = JLAST;
        this->indr[LW2] = 0;
        this->locc[JLAST] = LFIRST;
        this->locc[JFIRST] = 0;
        LENW--;
        (*LENU)--;
        if(*LROW==LW2)
            (*LROW)--;
        LW2 = LW2-1;
        //---------------------------------------------------------------
        //   Form the multiplier and store it in the  L  file.
        //---------------------------------------------------------------
        if(abs(WJ)<=SMALL)
            goto x490;
        AMULT = -WJ/VJ;
        L = this->lena-(*LENL);
        this->a[L] = AMULT;
        this->indr[L] = IV;
        this->indc[L] = IW;
        (*LENL)++;
        //---------------------------------------------------------------
        //   Add the appropriate multiple of row  iv  to row  iw.
        //   We use two different inner loops.  The first one is for the
        //   case where row  iw  is not at the end of storage.
        //---------------------------------------------------------------
        if(LENV==1)
            goto x490;
        LV2 = LV1+1;
        LV3 = (LV1+LENV)-1;
        if(LW2==*LROW)
            goto x400;
        //...............................................................
        //   This inner loop will be interrupted only if
        //   fill-in occurs enough to bump into the next row.
        //...............................................................
        for(LV = LV2; LV <= LV3; LV++) 
        {
            JV = this->indr[LV];
            LW = this->locc[JV];
            if(LW>0) 
            {
                // No fill-in.
                this->a[LW] += AMULT*this->a[LV];
                if(abs(this->a[LW])<=SMALL) 
                {
                    // Delete small computed element.
                    this->a[LW] = this->a[LW2];
                    J = this->indr[LW2];
                    this->indr[LW] = J;
                    this->indr[LW2] = 0;
                    this->locc[J] = LW;
                    this->locc[JV] = 0;
                    (*LENU)--;
                    LENW--;
                    LW2--;
                }
            }
            else 
            {
                //   Row  iw  doesn't have an element in column  jv  yet
                //  so there is a fill-in.
                if(this->indr[LW2+1]!=0)
                    goto x360;
                (*LENU)++;
                LENW++;
                LW2++;
                this->a[LW2] = AMULT*this->a[LV];
                this->indr[LW2] = JV;
                this->locc[JV] = LW2;
            }
        }
        goto x490;
        // Fill-in interrupted the previous loop.
        //   Move row  iw  to the end of the row file.
      x360:
        LV2 = LV;
        this->locr[IW] = (*LROW)+1;

        L = LW2-LW1+1;
        if(L > 0) 
        {
            INT loci, *locp;
            for(loci = LW1, locp = this->indr+LW1; loci <= LW2; loci++, locp++) 
            {
                (*LROW)++;
                this->locc[*locp] = *LROW;
            }
            LW2 = (*LROW)-L+1;
            mem_move(this->a+LW2,    this->a+LW1, L);
            mem_move(this->indr+LW2, this->indr+LW1, L);
            memclear(this->indr+LW1, L);
        }
        LW1 = this->locr[IW];
        LW2 = *LROW;
        //...............................................................
        //   Inner loop with row  iw  at the end of storage.
        //...............................................................
      x400:
        for(LV = LV2; LV <= LV3; LV++) 
        {
            JV = this->indr[LV];
            LW = this->locc[JV];
            if(LW>0) 
            {
                //  No fill-in.
                this->a[LW] += AMULT*this->a[LV];
                if(abs(this->a[LW])<=SMALL) 
                {
                    // Delete small computed element.
                    this->a[LW] = this->a[LW2];
                    J = this->indr[LW2];
                    this->indr[LW] = J;
                    this->indr[LW2] = 0;
                    this->locc[J] = LW;
                    this->locc[JV] = 0;
                    (*LENU)--;
                    LENW--;
                    LW2--;
                }
            }
            else 
            {
                // Row  iw  doesn't have an element in column  jv  yet
                // so there is a fill-in.
                (*LENU)++;
                LENW++;
                LW2++;
                this->a[LW2] = AMULT*this->a[LV];
                this->indr[LW2] = JV;
                this->locc[JV] = LW2;
            }
        }
        *LROW = LW2;
        //  The  k-th  element of row  iw  has been processed.
        //   Reset  swappd  before looking at the next element.
      x490:
        SWAPPD = false;
    }
    //==================================================================
    //    End of main elimination loop.
    //==================================================================

  // Cancel markers on row  iw.
  x600:
    this->lenr[IW] = LENW;
    if(LENW==0)
        goto x910;
    for(L = LW1; L <= LW2; L++) 
    {
        J = this->indr[L];
        this->locc[J] = 0;
    }
    // Move the diagonal element to the front of row  iw.
    // At this stage,  lenw gt 0  and  klast le n.
  x700:
    for(L = LW1; L <= LW2; L++) 
    {
        LDIAG = L;
        if(this->indr[L]==JFIRST)
            goto x730;
    }
    goto x910;

  x730:
    *DIAG = this->a[LDIAG];
    this->a[LDIAG] = this->a[LW1];
    this->a[LW1] = *DIAG;
    this->indr[LDIAG] = this->indr[LW1];
    this->indr[LW1] = JFIRST;

    // If an interchange is needed, repeat from the beginning with the
    //  new row  iw,  knowing that the opposite interchange cannot occur.
    if(SWAPPD)
        goto x100;
    *INFORM = LUSOL_INFORM_LUSUCCESS;
    goto x950;
  // Singular.
  x910:
    *DIAG = 0;
    *INFORM = LUSOL_INFORM_LUSINGULAR;
  // Force a compression if the file for  U  is much longer than the
  //    no. of nonzeros in  U  (i.e. if  lrow  is much bigger than  lenU).
  //    This should prevent memory fragmentation when there is far more
  //    memory than necessary  (i.e. when  lena  is huge).
  x950:
    LIMIT = (INT) (USPACE*(*LENU))+this->m+this->n+1000;
    if(*LROW>LIMIT)
        LU1REC(this->m,true,LROW,this->indr,this->lenr,this->locr);
    goto x990;
  // Not enough storage.
  x970:
    *INFORM = LUSOL_INFORM_ANEEDMEM;
  //  Exit.
  x990:
    ;
};

/* ==================================================================
   lu7rnk (check rank) assumes U is currently nrank by n
   and determines if row nrank contains an acceptable pivot.
   If not, the row is deleted and nrank is decreased by 1.
   jsing is an input parameter (not altered).  If jsing is positive,
   column jsing has already been judged dependent.  A substitute
   (if any) must be some other column.
   ------------------------------------------------------------------
   -- Jul 1987: First version.
   09 May 1988: First f77 version.
   ================================================================== */
template<class T>
void LUSOLrec<T>::LU7RNK(INT JSING, INT *LENU,
            INT *LROW, INT *NRANK, INT *INFORM, T *DIAG)
{
    TR UTOL1, UMAX;
    INT  IW, LENW, L1, L2, LMAX, L, JMAX, KMAX;

    L1 = 0;
    L2 = 0;

    UTOL1 = this->parmlu[LUSOL_RP_SMALLDIAG_U];
    *DIAG = 0;
    // Find Umax, the largest element in row nrank.
    IW = this->ip[*NRANK];
    LENW = this->lenr[IW];
    if(LENW==0)
        goto x400;
    L1 = this->locr[IW];
    L2 = (L1+LENW)-1;
    UMAX = 0;
    LMAX = L1;
    for(L = L1; L <= L2; L++) 
    {
        if(UMAX<abs(this->a[L])) 
        {
            UMAX = abs(this->a[L]);
            LMAX = L;
        }
    }
    //  Find which column that guy is in (in pivotal order).
    //    Interchange him with column nrank, then move him to be
    //    the new diagonal at the front of row nrank.
    *DIAG = this->a[LMAX];
    JMAX = this->indr[LMAX];
    for(KMAX = *NRANK; KMAX <= this->n; KMAX++) 
    {
        if(this->iq[KMAX]==JMAX)
            break;
    }
    this->iq[KMAX] = this->iq[*NRANK];
    this->iq[*NRANK] = JMAX;
    this->a[LMAX] = this->a[L1];
    this->a[L1] = *DIAG;
    this->indr[LMAX] = this->indr[L1];
    this->indr[L1] = JMAX;
    // See if the new diagonal is big enough.
    if(UMAX<=UTOL1)
        goto x400;
    if(JMAX==JSING)
        goto x400;
    //------------------------------------------------------------------
    //    The rank stays the same.
    //------------------------------------------------------------------
    *INFORM = LUSOL_INFORM_LUSUCCESS;
    return;
  //------------------------------------------------------------------
  //    The rank decreases by one.
  //------------------------------------------------------------------
  x400:
    *INFORM = LUSOL_INFORM_RANKLOSS;
    (*NRANK)--;
    if(LENW>0) 
    {
        // Delete row nrank from U.
        LENU = LENU-LENW;
        this->lenr[IW] = 0;
        for(L = L1; L <= L2; L++) 
        {
            this->indr[L] = 0;
        }
        if(L2==*LROW) 
        {
            // This row was at the end of the data structure.
            //  We have to reset lrow.
            // Preceding rows might already have been deleted, so we
            //  have to be prepared to go all the way back to 1.
            for(L = 1; L <= L2; L++) 
            {
                if(this->indr[*LROW]>0)
                    goto x900;
                (*LROW)--;
            }
        }
    }
  x900:
    ;
};

/* ==================================================================
   lu7zap  eliminates all nonzeros in column  jzap  of  U.
   It also sets  kzap  to the position of  jzap  in pivotal order.
   Thus, on exit we have  iq(kzap) = jzap.
   ------------------------------------------------------------------
   -- Jul 1987: nrank added.
   10 May 1988: First f77 version.
   ================================================================== */
template<class T>
void LUSOLrec<T>::LU7ZAP(INT JZAP, INT *KZAP, INT *LENU, INT *LROW,
            INT NRANK)
{
    INT K, I, LENI, LR1, LR2, L;

    for(K = 1; K <= NRANK; K++) 
    {
        I = this->ip[K];
        LENI = this->lenr[I];
        if(LENI==0)
            goto x90;
        LR1 = this->locr[I];
        LR2 = (LR1+LENI)-1;
        for(L = LR1; L <= LR2; L++) 
        {
            if(this->indr[L]==JZAP)
                goto x60;
        }
        goto x90;
      // Delete the old element.
      x60:
        this->a[L] = this->a[LR2];
        this->indr[L] = this->indr[LR2];
        this->indr[LR2] = 0;
        this->lenr[I] = LENI-1;
        (*LENU)--;
      // Stop if we know there are no more rows containing  jzap.
      x90:
        *KZAP = K;
        if(this->iq[K]==JZAP)
            goto x800;
    }
    // nrank must be smaller than n because we haven't found kzap yet.
    for(K = NRANK+1; K <= this->n; K++) 
    {
        *KZAP = K;
        if(this->iq[K]==JZAP)
            break;
    }
  // See if we zapped the last element in the file.
  x800:
    if(*LROW>0) 
    {
        if(this->indr[*LROW]==0)
            (*LROW)--;
    }
};

};
