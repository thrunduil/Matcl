#include "lusol.h"
#include <math.h>

/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   File  lusol8a
      lu8rpc
      Sparse LU update: Replace Column
      LUSOL's sparse implementation of the Bartels-Golub update.

   01 May 2002: Derived from LUSOL's original lu8a.f file.
   01 May 2002: Current version of lusol8a.f.
   ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */

/* ==================================================================
   lu8rpc  updates the LU factorization  A = L*U  when column  jrep
   is replaced by some vector  a(new).
   lu8rpc  is an implementation of the Bartels-Golub update,
   designed for the case where A is rectangular and/or singular.
   L is a product of stabilized eliminations (m x m, nonsingular).
   P U Q is upper trapezoidal (m x n, rank nrank).

   If  mode1 = 0,  the old column is taken to be zero
                   (so it does not have to be removed from  U).
   If  mode1 = 1,  the old column need not have been zero.
   If  mode2 = 0,  the new column is taken to be zero.
                   v(*)  is not used or altered.
   If  mode2 = 1,  v(*)  must contain the new column  a(new).
                   On exit,  v(*)  will satisfy  L*v = a(new).
   If  mode2 = 2,  v(*)  must satisfy  L*v = a(new).

   The array  w(*)  is not used or altered.
   On entry, all elements of  locc  are assumed to be zero.
   On a successful exit (inform ne 7), this will again be true.
   On exit:

   inform = -1  if the rank of U decreased by 1.
   inform =  0  if the rank of U stayed the same.
   inform =  1  if the rank of U increased by 1.
   inform =  2  if the update seemed to be unstable
                (diag much bigger than vnorm).
   inform =  7  if the update was not completed (lack of storage).
   inform =  8  if jrep is not between 1 and n.
   ------------------------------------------------------------------
   -- Jan 1985: Original F66 version.
   -- Jul 1987: Modified to maintain U in trapezoidal form.
   10 May 1988: First f77 version.
   16 Oct 2000: Added test for instability (inform = 2).
   ================================================================== */
namespace lusol
{
template<class T>
void LUSOLrec<T>::LU8RPC(INT MODE1, INT MODE2,
            INT JREP, T V[], T W[],
            INT *INFORM, T *DIAG, TR *VNORM)
{
    bool SINGLR;
    INT    LPRINT, NRANK, LENL, LENU, LROW, NRANK0, KREP, KLAST, IW, L1, J1, JSING;
    TR     UTOL1, UTOL2;

    LPRINT = this->luparm[LUSOL_IP_PRINTLEVEL];
    NRANK  = this->luparm[LUSOL_IP_RANK_U];
    LENL   = this->luparm[LUSOL_IP_NONZEROS_L];
    LENU   = this->luparm[LUSOL_IP_NONZEROS_U];
    LROW   = this->luparm[LUSOL_IP_NONZEROS_ROW];
    UTOL1  = this->parmlu[LUSOL_RP_SMALLDIAG_U];
    UTOL2  = this->parmlu[LUSOL_RP_EPSDIAG_U];
    NRANK0 = NRANK;
    *DIAG  = 0;
    *VNORM = 0;
    if(JREP<1)
        goto x980;
    if(JREP>this->n)
        goto x980;
    //------------------------------------------------------------------
    //    If mode1 = 0, there are no elements to be removed from  U
    //    but we still have to set  krep  (using a backward loop).
    //    Otherwise, use lu7zap to remove column  jrep  from  U
    //    and set  krep  at the same time.
    //------------------------------------------------------------------
    if(MODE1==LUSOL_UPDATE_OLDEMPTY) 
    {
        KREP = this->n+1;
      x10:
        KREP--;
        if(this->iq[KREP]!=JREP)
            goto x10;
    }
    else
        LU7ZAP(JREP,&KREP,&LENU,&LROW,NRANK);

    //------------------------------------------------------------------
    //    Insert a new column of u and find klast.
    //------------------------------------------------------------------
    if(MODE2==LUSOL_UPDATE_NEWEMPTY) 
    {
        KLAST = 0;
    }
    else 
    {
        if(MODE2==LUSOL_UPDATE_NEWNONEMPTY) 
        {
            // Transform v = a(new) to satisfy  L*v = a(new).
            LU6SOL(LUSOL_SOLVE_Lv_v, V,W, nullptr, INFORM);
        }
        else if(V==nullptr)
        {
            // Otherwise, the V vector is taken to satisfy this already, or stored earlier.
            V=this->vLU6L;
        };
      

        //  Insert into  U  any nonzeros in the top of  v.
        //   row  ip(klast)  will contain the last nonzero in pivotal order.
        //   Note that  klast  will be in the range  ( 0, nrank ).
        LU7ADD(JREP,V,LENL,&LENU,&LROW,NRANK,INFORM,&KLAST,VNORM);
        if(*INFORM==LUSOL_INFORM_ANEEDMEM)
            goto x970;
    }
    //------------------------------------------------------------------
    //    In general, the new column causes U to look like this:
    //                krep        n                 krep  n
    //               ....a.........          ..........a...
    //                .  a        .           .        a  .
    //                 . a        .            .       a  .
    //                  .a        .             .      a  .
    //       P U Q =     a        .    or        .     a  .
    //                   b.       .               .    a  .
    //                   b .      .                .   a  .
    //                   b  .     .                 .  a  .
    //                   b   ......                  ..a...  nrank
    //                   c                             c
    //                   c                             c
    //                   c                             c     m
    //    klast points to the last nonzero "a" or "b".
    //    klast = 0 means all "a" and "b" entries are zero.
    //------------------------------------------------------------------
    if(MODE2==LUSOL_UPDATE_NEWEMPTY) 
    {
        if(KREP>NRANK)
            goto x900;
    }
    else if(NRANK<this->m) 
    {
        //  Eliminate any "c"s (in either case).
        //   Row nrank + 1 may end up containing one nonzero.
        LU7ELM(JREP,V,&LENL,&LROW,NRANK,INFORM,DIAG);
        if(*INFORM==LUSOL_INFORM_ANEEDMEM)
            goto x970;
        if(*INFORM==LUSOL_INFORM_LUSINGULAR) 
        {
            // The nonzero is apparently significant.
            //   Increase nrank by 1 and make klast point to the bottom.
            NRANK++;
            KLAST = NRANK;
        }
    }
    if(NRANK<this->n) 
    {
        //  The column rank is low.
        //   In the first case, we want the new column to end up in
        //   position nrank, so the trapezoidal columns will have a chance
        //   later on (in lu7rnk) to pivot in that position.
        //   Otherwise the new column is not part of the triangle.  We
        //   swap it into position nrank so we can judge it for singularity.
        //   lu7rnk might choose some other trapezoidal column later.
        if(KREP<NRANK)
            KLAST = NRANK;
        else 
        {
            this->iq[KREP] = this->iq[NRANK];
            this->iq[NRANK] = JREP;
            KREP = NRANK;
        }
    }
    //------------------------------------------------------------------
    //    If krep .lt. klast, there are some "b"s to eliminate:
    //                 krep
    //               ....a.........
    //                .  a        .
    //                 . a        .
    //                  .a        .
    //       P U Q =     a        .  krep
    //                   b.       .
    //                   b .      .
    //                   b  .     .
    //                   b   ......  nrank
    //    If krep .eq. klast, there are no "b"s, but the last "a" still
    //    has to be moved to the front of row krep (by lu7for).
    //------------------------------------------------------------------
    if(KREP<=KLAST) 
    {
        //  Perform a cyclic permutation on the current pivotal order,
        //   and eliminate the resulting row spike.  krep becomes klast.
        //   The final diagonal (if any) will be correctly positioned at
        //   the front of the new krep-th row.  nrank stays the same.
        LU7CYC(KREP,KLAST,this->ip);
        LU7CYC(KREP,KLAST,this->iq);
        LU7FOR(KREP,KLAST,&LENL,&LENU,&LROW,INFORM,DIAG);
        if(*INFORM==LUSOL_INFORM_ANEEDMEM)
            goto x970;
        KREP = KLAST;
        // Test for instability (diag much bigger than vnorm).
        SINGLR = (bool) ((*VNORM)<UTOL2*abs(*DIAG));
        if(SINGLR)
            goto x920;
    }
    //------------------------------------------------------------------
    //    Test for singularity in column krep (where krep .le. nrank).
    //------------------------------------------------------------------
    *DIAG = 0;
    IW = this->ip[KREP];
    SINGLR = (bool) (this->lenr[IW]==0);
    if(!SINGLR) 
    {
        L1 = this->locr[IW];
        J1 = this->indr[L1];
        SINGLR = (bool) (J1!=JREP);
        if(!SINGLR) 
        {
            *DIAG = this->a[L1];
            SINGLR = (bool) (abs(*DIAG)<=UTOL1 || abs(*DIAG)<=UTOL2*(*VNORM));
        }
    }
    if(SINGLR && KREP<NRANK) 
    {
        //  Perform cyclic permutations to move column jrep to the end.
        //   Move the corresponding row to position nrank
        //   then eliminate the resulting row spike.
        LU7CYC(KREP,NRANK,this->ip);
        LU7CYC(KREP,this->n,this->iq);
        LU7FOR(KREP,NRANK,&LENL,&LENU,&LROW,INFORM,DIAG);
        if(*INFORM==LUSOL_INFORM_ANEEDMEM)
            goto x970;
    }
    //  Find the best column to be in position nrank.
    //    If singlr, it can't be the new column, jrep.
    //    If nothing satisfactory exists, nrank will be decreased.
    if(SINGLR || NRANK<this->n) 
    {
        JSING = 0;
        if(SINGLR)
            JSING = JREP;
        LU7RNK(JSING,&LENU,&LROW,&NRANK,INFORM,DIAG);
    }
  //------------------------------------------------------------------
  //    Set inform for exit.
  //------------------------------------------------------------------
  x900:
    if(NRANK==NRANK0)
        *INFORM = LUSOL_INFORM_LUSUCCESS;
    else if(NRANK<NRANK0) 
    {
        *INFORM = LUSOL_INFORM_RANKLOSS;
        if(NRANK0==this->n) 
        {
            if(LPRINT>=LUSOL_MSG_SINGULARITY)
                LUSOL_report(this, 0, (char*) "lu8rpc  warning...\nSingularity after replacing column.    jrep=%8d    diag=%g\n",
                            JREP,DIAG); //FIXME:
        }
    }
    else
        *INFORM = LUSOL_INFORM_LUSINGULAR;
    goto x990;
  //  Instability.
  x920:
    *INFORM = LUSOL_INFORM_LUUNSTABLE;
    if(LPRINT>=LUSOL_MSG_SINGULARITY)
        LUSOL_report(this, 0, (char*) "lu8rpc  warning...\nInstability after replacing column.    jrep=%8d    diag=%g\n",
                        JREP,DIAG); //FIXME:
    goto x990;
  //  Not enough storage.
  x970:
    *INFORM = LUSOL_INFORM_ANEEDMEM;
    if(LPRINT>=LUSOL_MSG_SINGULARITY)
        LUSOL_report(this, 0, (char*) "lu8rpc  error...\nInsufficient memory.    lena=%8d\n",
                        this->lena); //FIXME:
    goto x990;
  // jrep  is out of range.
  x980:
    *INFORM = LUSOL_INFORM_FATALERR;
    if(LPRINT>=LUSOL_MSG_SINGULARITY)
        LUSOL_report(this, 0, (char*) "lu8rpc  error...\njrep  is out of range.    m=%8d    n=%8d    jrep=%8d\n",
                        this->m,this->n,JREP); //FIXME:
  // Exit.
  x990:
    this->luparm[LUSOL_IP_INFORM]       = *INFORM;
    this->luparm[LUSOL_IP_UPDATECOUNT]++;
    this->luparm[LUSOL_IP_RANK_U]       = NRANK;
    this->luparm[LUSOL_IP_NONZEROS_L]   = LENL;
    this->luparm[LUSOL_IP_NONZEROS_U]   = LENU;
    this->luparm[LUSOL_IP_NONZEROS_ROW] = LROW;
};

};
