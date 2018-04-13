#include "lusol.h"
#include <math.h>
#include <stdlib.h>

/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   File  lusol6a
      lu6sol   lu6L     lu6Lt     lu6U     Lu6Ut   lu6LD   lu6chk
   ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   26 Apr 2002: lu6 routines put into a separate file.
   15 Dec 2002: lu6sol modularized via lu6L, lu6Lt, lu6U, lu6Ut.
                lu6LD implemented to allow solves with LDL' or L|D|L'.
   15 Dec 2002: Current version of lusol6a.f.
   ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */

/* ==================================================================
   lu6chk  looks at the LU factorization  A = L*U.
   If mode = 1, lu6chk is being called by lu1fac.
   (Other modes not yet implemented.)
   ------------------------------------------------------------------
   The important input parameters are

                  lprint = luparm(2)
                  keepLU = luparm(8)
                  Utol1  = parmlu(4)
                  Utol2  = parmlu(5)

   and the significant output parameters are

                  inform = luparm(10)
                  nsing  = luparm(11)
                  jsing  = luparm(12)
                  jumin  = luparm(19)
                  Lmax   = parmlu(11)
                  Umax   = parmlu(12)
                  DUmax  = parmlu(13)
                  DUmin  = parmlu(14)
                  and      w(*).

   Lmax  and Umax  return the largest elements in L and U.
   DUmax and DUmin return the largest and smallest diagonals of U
                   (excluding diagonals that are exactly zero).
   In general, w(j) is set to the maximum absolute element in
   the j-th column of U.  However, if the corresponding diagonal
   of U is small in absolute terms or relative to w(j)
   (as judged by the parameters Utol1, Utol2 respectively),
   then w(j) is changed to - w(j).
   Thus, if w(j) is not positive, the j-th column of A
   appears to be dependent on the other columns of A.
   The number of such columns, and the position of the last one,
   are returned as nsing and jsing.
   Note that nrank is assumed to be set already, and is not altered.
   Typically, nsing will satisfy      nrank + nsing = n,  but if
   Utol1 and Utol2 are rather large,  nsing > n - nrank   may occur.
   If keepLU = 0,
   Lmax  and Umax  are already set by lu1fac.
   The diagonals of U are in the top of A.
   Only Utol1 is used in the singularity test to set w(*).
   inform = 0  if  A  appears to have full column rank  (nsing = 0).
   inform = 1  otherwise  (nsing .gt. 0).
   ------------------------------------------------------------------
   00 Jul 1987: Early version.
   09 May 1988: f77 version.
   11 Mar 2001: Allow for keepLU = 0.
   17 Nov 2001: Briefer output for singular factors.
   05 May 2002: Comma needed in format 1100 (via Kenneth Holmstrom).
   06 May 2002: With keepLU = 0, diags of U are in natural order.
                They were not being extracted correctly.
   23 Apr 2004: TRP can judge singularity better by comparing
                all diagonals to DUmax.
   27 Jun 2004: (PEG) Allow write only if nout .gt. 0.
   ================================================================== */

namespace lusol
{
template<class T>
void LUSOLrec<T>::LU6CHK(INT /*MODE*/, INT LENA2, INT *INFORM)
{
    bool KEEPLU;
    INT    I, J, JSING, JUMIN, K, L, L1, L2, LENL, LPRINT, NDEFIC, NRANK, NSING;
    TR     AIJ, DIAG, DUMAX, DUMIN, LMAX, UMAX, UTOL1, UTOL2;

    LPRINT = this->luparm[LUSOL_IP_PRINTLEVEL];
    KEEPLU = (bool) (this->luparm[LUSOL_IP_KEEPLU]!=0);
    NRANK = this->luparm[LUSOL_IP_RANK_U];
    LENL  = this->luparm[LUSOL_IP_NONZEROS_L];
    UTOL1 = this->parmlu[LUSOL_RP_SMALLDIAG_U];
    UTOL2 = this->parmlu[LUSOL_RP_EPSDIAG_U];
    *INFORM = LUSOL_INFORM_LUSUCCESS;
    LMAX  = 0;
    UMAX  = 0;
    NSING = 0;
    JSING = 0;
    JUMIN = 0;
    DUMAX = 0;
    DUMIN = lusol_user_params<TR>::LUSOL_BIGNUM;
    
    memclear((this->wr+1), this->n);

    if(KEEPLU) 
    {
        //--------------------------------------------------------------
        //      Find  Lmax.
        //--------------------------------------------------------------
        for(L = (LENA2+1)-LENL; L <= LENA2; L++) 
        {
            LMAX = maximum<TR>(LMAX,abs(this->a[L]));
        }
        //--------------------------------------------------------------
        //  Find Umax and set w(j) = maximum element in j-th column of U.
        //--------------------------------------------------------------
        for(K = 1; K <= NRANK; K++) 
        {
            I = this->ip[K];
            L1 = this->locr[I];
            L2 = (L1+this->lenr[I])-1;
            for(L = L1; L <= L2; L++) 
            {
                J = this->indr[L];
                AIJ = abs(this->a[L]);
                this->wr[J] = maximum(this->wr[J],AIJ);
                UMAX = maximum<TR>(UMAX,AIJ);
            }
        };
        //--------------------------------------------------------------
        //Negate w(j) if the corresponding diagonal of U is
        //too small in absolute terms or relative to the other elements
        //in the same column of  U.
        //Also find DUmax and DUmin, the extreme diagonals of U.
        //--------------------------------------------------------------
        for(K = 1; K <= this->n; K++) 
        {
            J = this->iq[K];
            if(K>NRANK)
                DIAG = 0;
            else 
            {
                I = this->ip[K];
                L1 = this->locr[I];
                DIAG = abs(this->a[L1]);
                DUMAX = maximum(DUMAX,DIAG);
                if(DUMIN>DIAG) 
                {
                    DUMIN = DIAG;
                    JUMIN = J;
                }
            }
            if(DIAG<=UTOL1 || DIAG<=UTOL2*this->wr[J]) 
            {
                NSING++;
                JSING = J;
                this->wr[J] = -this->wr[J];
            }
        }
        this->parmlu[LUSOL_RP_MAXMULT_L] = LMAX;
        this->parmlu[LUSOL_RP_MAXELEM_U] = UMAX;
    }
    else 
    {
        //--------------------------------------------------------------
        //keepLU = 0.
        //Only diag(U) is stored.  Set w(*) accordingly.
        //--------------------------------------------------------------
        for(K = 1; K <= this->n; K++) 
        {
            J = this->iq[K];
            if(K>NRANK)
                DIAG = 0;
            else 
            {
                DIAG = abs(this->diagU[J]);
                this->wr[J] = DIAG;
                DUMAX = maximum(DUMAX,DIAG);
                if(DUMIN>DIAG) 
                {
                    DUMIN = DIAG;
                    JUMIN = J;
                }
            }
            if(DIAG<=UTOL1) 
            {
                NSING++;
                JSING = J;
                this->wr[J] = -this->wr[J];
            }
        }
    }
    //-----------------------------------------------------------------
    //    Set output parameters.
    //-----------------------------------------------------------------
    if(JUMIN==0)
        DUMIN = 0;
    this->luparm[LUSOL_IP_SINGULARITIES]  = NSING;
    this->luparm[LUSOL_IP_SINGULARINDEX]  = JSING;
    this->luparm[LUSOL_IP_COLINDEX_DUMIN] = JUMIN;
    this->parmlu[LUSOL_RP_MAXELEM_DIAGU]  = DUMAX;
    this->parmlu[LUSOL_RP_MINELEM_DIAGU]  = DUMIN;

    //   The matrix has been judged singular.
    if(NSING>0) 
    {
        *INFORM = LUSOL_INFORM_LUSINGULAR;
        NDEFIC = this->n-NRANK;
        if(LPRINT>=LUSOL_MSG_SINGULARITY) 
        {
            LUSOL_report(this, 0, (char*) "Singular(m%cn)  rank:%9d  n-rank:%8d  nsing:%9d\n",
                             relationChar(this->m, this->n),NRANK,NDEFIC,NSING); //FIXME:
        }
    }
    // Exit.
    this->luparm[LUSOL_IP_INFORM] = *INFORM;
};

// Create a row-based version of L0
template<class T>
void LUSOLrec<T>::LU1L0T(LUSOLmat<T> **mat)
{
    INT  K, L, L1, LEN, LENL0, NUML0, I, J, LL;
    INT  *lsumr;

    NUML0 = this->luparm[LUSOL_IP_COLCOUNT_L0];
    LENL0 = this->luparm[LUSOL_IP_NONZEROS_L0];
    if(*mat != nullptr)
    {
        LUSOLmat<T>::destroy(mat);
    }
    if(*mat == nullptr)
    {
        *mat = LUSOLmat<T>::create(this->m, LENL0);
    };
    if(*mat == nullptr)
    {
        return;
    };

    lsumr = (INT *) malloc((this->m+1)*sizeof(*lsumr));
    if (lsumr == nullptr)
    {
        LUSOLmat<T>::destroy(mat);
        return;
    };

    // Compute row non-zero counts
    memclear(((*mat)->vlen+1), this->m);

    L1 = this->lena+1;
    for(K = 1; K <= NUML0; K++) 
    {
        LEN = this->lenc[K];
        L = L1;
        L1 -= LEN;
        for(; LEN > 0; LEN--) 
        {
            L--;
            J = this->indc[L];
            (*mat)->vlen[J]++;
        }
    };
  
    // Cumulate the row counts to get vector offsets; note order!
    K = this->m;
    lsumr[K] = 1;  // Stick with Fortran array offset for consistency
    for(; K > 1; K--)
        lsumr[K-1] = lsumr[K] + (*mat)->vlen[K];
  
    // Map the matrix into row order
    L1 = this->lena+1;
    for(K = 1; K <= NUML0; K++) 
    {
        LEN = this->lenc[K];
        L = L1;
        L1 -= LEN;
        I = this->indr[L1];    // Obtain the (permuted) column index
        for(; LEN > 0; LEN--) 
        {
            L--;
            J = this->indc[L];
            LL = lsumr[J]++;
            (*mat)->a[LL] = this->a[L];
            (*mat)->indr[LL] = I;
            (*mat)->indc[LL] = J;
        }
    }
  
    // Clean up
    lusol_free(lsumr);
};

// Solve L0T v = v based on separate row-based version of L0
template<class T>
void LUSOLrec<T>::LU6L0T_v(LUSOLmat<T> *mat, T V[], INT [])
{
    INT  LEN, IPIV, K, L, L1, LENL0, NUML0;
    TR   SMALL;
    T    VPIV;

    T    *aptr;
    INT  *jptr;

    NUML0 = this->m;
    LENL0 = this->luparm[LUSOL_IP_NONZEROS_L0];
    SMALL = this->parmlu[LUSOL_RP_ZEROTOLERANCE];
    L1 = LENL0+1;
    for(K = 1; K <= NUML0; K++) 
    {
        LEN = mat->vlen[K];
        if(LEN == 0)
            continue;
        L = L1;
        L1 -= LEN;
        IPIV = mat->indc[L1];
        VPIV = V[IPIV];
        if(abs(VPIV)>SMALL) 
        {
            // ***** This loop could be coded specially.

            L--;
            for(aptr = mat->a+L, jptr = mat->indr+L; LEN > 0; LEN--, aptr--, jptr--)
                V[*jptr] += (*aptr) * VPIV;
        }
    };
};


//------------------------------------------------------------------
//   lu6L   solves   L v = v(input).
//------------------------------------------------------------------
//   15 Dec 2002: First version derived from lu6sol.
//   15 Dec 2002: Current version.
//------------------------------------------------------------------
template<class T>
void LUSOLrec<T>::LU6L(INT *INFORM, T V[], INT [])
{
    INT  JPIV, K, L, L1, LEN, LENL, LENL0, NUML, NUML0;
    TR   SMALL;
    T    VPIV;

    T    *aptr;
    INT  *iptr, *jptr;

    NUML0 = this->luparm[LUSOL_IP_COLCOUNT_L0];
    LENL0 = this->luparm[LUSOL_IP_NONZEROS_L0];
    LENL  = this->luparm[LUSOL_IP_NONZEROS_L];
    SMALL = this->parmlu[LUSOL_RP_ZEROTOLERANCE];
    *INFORM = LUSOL_INFORM_LUSUCCESS;
    L1 = this->lena+1;
    for(K = 1; K <= NUML0; K++) 
    {
        LEN = this->lenc[K];
        L = L1;
        L1 -= LEN;
        JPIV = this->indr[L1];
        VPIV = V[JPIV];
        if(abs(VPIV)>SMALL) 
        {
            //***** This loop could be coded specially.

            L--;
            for(aptr = this->a+L, iptr = this->indc+L; LEN > 0; LEN--, aptr--, iptr--)
                V[*iptr] += (*aptr) * VPIV;
        }
    }
    L = (this->lena-LENL0)+1;
    NUML = LENL-LENL0;
    // ***** This loop could be coded specially.

    L--;
    for(aptr = this->a+L, jptr = this->indr+L, iptr = this->indc+L;
        NUML > 0; NUML--, aptr--, jptr--, iptr--) 
    {
        if(abs(V[*jptr])>SMALL)
            V[*iptr] += (*aptr) * V[*jptr];
    }

    // Exit.
    this->luparm[LUSOL_IP_INFORM] = *INFORM;
};

//==================================================================
//   lu6LD  assumes lu1fac has computed factors A = LU of a
//   symmetric definite or quasi-definite matrix A,
//   using Threshold Symmetric Pivoting (TSP),   luparm(6) = 3,
//   or    Threshold Diagonal  Pivoting (TDP),   luparm(6) = 4.
//   It also assumes that no updates have been performed.
//   In such cases,  U = D L', where D = diag(U).
//   lu6LDL returns v as follows:
//
//   mode
//    1    v  solves   L D v = v(input).
//    2    v  solves   L|D|v = v(input).
//------------------------------------------------------------------
//   15 Dec 2002: First version of lu6LD.
//   15 Dec 2002: Current version.
//==================================================================
template<class T>
void LUSOLrec<T>::LU6LD(INT *INFORM, INT MODE, T V[], INT [])
{
    INT  IPIV, K, L, L1, LEN, NUML0;
    T    DIAG;
    TR   SMALL;
    T    VPIV;

    T    *aptr;
    INT  *jptr;

    //  Solve L D v(new) = v  or  L|D|v(new) = v, depending on mode.
    //  The code for L is the same as in lu6L,
    //  but when a nonzero entry of v arises, we divide by
    //  the corresponding entry of D or |D|.

    NUML0 = this->luparm[LUSOL_IP_COLCOUNT_L0];
    SMALL = this->parmlu[LUSOL_RP_ZEROTOLERANCE];
    *INFORM = LUSOL_INFORM_LUSUCCESS;
    L1 = this->lena+1;
    for(K = 1; K <= NUML0; K++) 
    {
        LEN = this->lenc[K];
        L = L1;
        L1 -= LEN;
        IPIV = this->indr[L1];
        VPIV = V[IPIV];
        if(abs(VPIV)>SMALL) 
        {
            // ***** This loop could be coded specially.

            L--;
            for(aptr = this->a+L, jptr = this->indc+L; LEN > 0; LEN--, aptr--, jptr--)
                V[*jptr] += (*aptr) * VPIV;
            // Find diag = U(ipiv,ipiv) and divide by diag or |diag|.
            L = this->locr[IPIV];
            DIAG = this->a[L];
            if(MODE==2)
                DIAG = abs(DIAG);

            V[IPIV] = VPIV/DIAG;
        };
    };
};


//==================================================================
//   lu6Lt  solves   L'v = v(input).
//------------------------------------------------------------------
//   15 Dec 2002: First version derived from lu6sol.
//   15 Dec 2002: Current version.
//==================================================================
template<class T>
void LUSOLrec<T>::LU6LT(INT *INFORM, T V[], INT [])
{
    INT  K, L, L1, L2, LEN, LENL, LENL0, NUML0;
    TR   SMALL;
    T    SUM;

    T    *aptr;
    INT  *iptr, *jptr;

    NUML0 = this->luparm[LUSOL_IP_COLCOUNT_L0];
    LENL0 = this->luparm[LUSOL_IP_NONZEROS_L0];
    LENL  = this->luparm[LUSOL_IP_NONZEROS_L];
    SMALL = this->parmlu[LUSOL_RP_ZEROTOLERANCE];
    *INFORM = LUSOL_INFORM_LUSUCCESS;
    L1 = (this->lena-LENL)+1;
    L2 = this->lena-LENL0;
  
    // ***** This loop could be coded specially.

    for(L = L1, aptr = this->a+L1, iptr = this->indr+L1, jptr = this->indc+L1;
            L <= L2; L++, aptr++, iptr++, jptr++) 
    {
        if(abs(V[*jptr])>SMALL)
            V[*iptr] += (*aptr) * V[*jptr];
    }

    #ifdef UseRowBasedL0
        LU6L0T_v(LUSOL, this->L0, V, NZidx);   
    #else

        for(K = NUML0; K >= 1; K--) 
        {
            LEN = this->lenc[K];
            SUM = 0;
            L1 = L2+1;
            L2 += LEN;
            // ***** This loop could be coded specially.

            for(L = L1, aptr = this->a+L1, jptr = this->indc+L1;
                L <= L2; L++, aptr++, jptr++)
            {
                SUM += (*aptr) * V[*jptr];
            };
            V[this->indr[L1]] += SUM;
        }  
    #endif

    // Exit.
    this->luparm[LUSOL_IP_INFORM] = *INFORM;
};

/* ==================================================================
   lu6U   solves   U w = v.          v  is not altered.
   ------------------------------------------------------------------
   15 Dec 2002: First version derived from lu6sol.
   15 Dec 2002: Current version.
   ================================================================== */
template<class T>
void LUSOLrec<T>::LU6U(INT *INFORM, T V[], T W[], INT [])
{
    INT  I, J, K, KLAST, L, L1, L2, L3, NRANK, NRANK1;
    TR   SMALL;
    T    Z;

    T    *aptr;
    INT  *jptr;

    NRANK = this->luparm[LUSOL_IP_RANK_U];
    SMALL = this->parmlu[LUSOL_RP_ZEROTOLERANCE];
    *INFORM = LUSOL_INFORM_LUSUCCESS;
    NRANK1 = NRANK+1;
    // Find the first nonzero in v(1:nrank), counting backwards.
    for(KLAST = NRANK; KLAST >= 1; KLAST--) 
    {
        I = this->ip[KLAST];
        if(abs(V[I])>SMALL)
            break;
    }
    for(K = KLAST+1; K <= this->n; K++) 
    {
        J = this->iq[K];
        W[J] = 0;
    }
    // Do the back-substitution, using rows 1:klast of U.
    for(K = KLAST; K >= 1; K--) 
    {
        I = this->ip[K];
        Z = V[I];
        L1 = this->locr[I];
        L2 = L1+1;
        L3 = (L1+this->lenr[I])-1;
        // ***** This loop could be coded specially.

        for(L = L2, aptr = this->a+L2, jptr = this->indr+L2; L <= L3; L++, aptr++, jptr++)
            Z -= (*aptr) * W[*jptr];
        J = this->iq[K];
        if(abs(Z)<=SMALL)
            Z = 0;
        else
            Z /= this->a[L1];
        W[J] = Z;
    }
    // Compute residual for overdetermined systems.
    TR Z2 = 0;
    for(K = NRANK1; K <= this->m; K++) 
    {
        I = this->ip[K];
        Z2 += abs(V[I]);
    }
    // Exit.
    if(Z2>0)
        *INFORM = LUSOL_INFORM_LUSINGULAR;
    this->luparm[LUSOL_IP_INFORM]     = *INFORM;
    this->parmlu[LUSOL_RP_RESIDUAL_U] = Z2;
};

//==================================================================
//   lu6Ut  solves   U'v = w.          w  is destroyed.
//------------------------------------------------------------------
//   15 Dec 2002: First version derived from lu6sol.
//   15 Dec 2002: Current version.
//==================================================================
template<class T>
void LUSOLrec<T>::LU6UT(INT *INFORM, T V[], T W[], INT [])
{
    INT  I, J, K, L, L1, L2, NRANK, NRANK1;
    TR   SMALL;
    T    Z;

    T    *aptr;
    INT  *jptr;

    NRANK = this->luparm[LUSOL_IP_RANK_U];
    SMALL = this->parmlu[LUSOL_RP_ZEROTOLERANCE];
    *INFORM = LUSOL_INFORM_LUSUCCESS;
    NRANK1 = NRANK+1;
    for(K = NRANK1; K <= this->m; K++) 
    {
        I = this->ip[K];
        V[I] = 0;
    }
    //  Do the forward-substitution, skipping columns of U(transpose)
    //    when the associated element of w(*) is negligible.
    for(K = 1; K <= NRANK; K++) 
    {
        I = this->ip[K];
        J = this->iq[K];
        Z = W[J];
        if(abs(Z)<=SMALL) 
        {
            V[I] = 0;
            continue;
        }
        L1 = this->locr[I];
        Z /= this->a[L1];
        V[I] = Z;
        L2 = (L1+this->lenr[I])-1;
        L1++;
        // ***** This loop could be coded specially.

        for(L = L1, aptr = this->a+L1, jptr = this->indr+L1; L <= L2; L++, aptr++, jptr++)
            W[*jptr] -= Z * (*aptr);
    }
    // Compute residual for overdetermined systems.
    TR Z2 = 0;
    for(K = NRANK1; K <= this->n; K++) 
    {
        J = this->iq[K];
        Z2 += abs(W[J]);
    }
    //  Exit.
    if(Z2>0)
        *INFORM = LUSOL_INFORM_LUSINGULAR;

    this->luparm[LUSOL_IP_INFORM]     = *INFORM;
    this->parmlu[LUSOL_RP_RESIDUAL_U] = Z2;
};

/* ==================================================================
   lu6sol  uses the factorization  A = L U  as follows:
   ------------------------------------------------------------------
   mode
    1    v  solves   L v = v(input).   w  is not touched.
    2    v  solves   L'v = v(input).   w  is not touched.
    3    w  solves   U w = v.          v  is not altered.
    4    v  solves   U'v = w.          w  is destroyed.
    5    w  solves   A w = v.          v  is altered as in 1.
    6    v  solves   A'v = w.          w  is destroyed.

   If mode = 3,4,5,6, v and w must not be the same arrays.
   If lu1fac has just been used to factorize a symmetric matrix A
   (which must be definite or quasi-definite), the factors A = L U
   may be regarded as A = LDL', where D = diag(U).  In such cases,

   mode
    7    v  solves   A v = L D L'v = v(input).   w  is not touched.
    8    v  solves       L |D| L'v = v(input).   w  is not touched.

   ip(*), iq(*)      hold row and column numbers in pivotal order.
   lenc(k)           is the length of the k-th column of initial L.
   lenr(i)           is the length of the i-th row of U.
   locc(*)           is not used.
   locr(i)           is the start  of the i-th row of U.

   U is assumed to be in upper-trapezoidal form (nrank by n).
   The first entry for each row is the diagonal element
   (according to the permutations  ip, iq).  It is stored at
   location locr(i) in a(*), indr(*).

   On exit, inform = 0 except as follows.
     if(mode = 3,4,5,6 and if U (and hence A) is singular,)
     inform = 1 if there is a nonzero residual in solving the system
     involving U.  parmlu(20) returns the norm of the residual.
   ------------------------------------------------------------------
     July 1987: Early version.
   09 May 1988: f77 version.
   27 Apr 2000: Abolished the dreaded "computed go to".
                But hard to change other "go to"s to "if then else".
   15 Dec 2002: lu6L, lu6Lt, lu6U, lu6Ut added to modularize lu6sol.
   ================================================================== */
template<class T>
void LUSOLrec<T>::LU6SOL(INT MODE, T V[], T W[], INT NZidx[], INT *INFORM)
{
    if(MODE==LUSOL_SOLVE_Lv_v) 
    {                                       //      Solve  L v(new) = v.
        LU6L(INFORM,V, NZidx);
    }
    else if(MODE==LUSOL_SOLVE_Ltv_v) 
    {                                       //      Solve  L'v(new) = v.
        LU6LT(INFORM,V, NZidx);
    }
    else if(MODE==LUSOL_SOLVE_Uw_v) 
    {                                       //      Solve  U w = v.
        LU6U(INFORM,V,W, NZidx);
    }
    else if(MODE==LUSOL_SOLVE_Utv_w) 
    {                                       //      Solve  U'v = w.
        LU6UT(INFORM,V,W, NZidx);
    }
    else if(MODE==LUSOL_SOLVE_Aw_v) 
    {                                       //      Solve  A w      = v 
        LU6L(INFORM,V, NZidx);              //      via     L v(new) = v 
        LU6U(INFORM,V,W, nullptr);             //      ... and U w = v(new). 
    }
    else if(MODE==LUSOL_SOLVE_Atv_w) 
    {                                       //      Solve  A'v = w 
        LU6UT(INFORM,V,W, NZidx);           //      via      U'v = w 
        LU6LT(INFORM,V, nullptr);              //      ... and  L'v(new) = v. 
    }
    else if(MODE==LUSOL_SOLVE_Av_v) 
    {                                       //      Solve  LDv(bar) = v 
        LU6LD(INFORM,1,V, NZidx);           //      and    L'v(new) = v(bar). 
        LU6LT(INFORM,V, nullptr);
    }
    else if(MODE==LUSOL_SOLVE_LDLtv_v) 
    {                                       //      Solve  L|D|v(bar) = v 
        LU6LD(INFORM,2,V, NZidx);           //      and    L'v(new) = v(bar). 
        LU6LT(INFORM,V, nullptr);
    }
};

};
