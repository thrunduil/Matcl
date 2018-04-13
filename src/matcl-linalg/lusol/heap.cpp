#include "lusol.h"
#include "heap.h"

/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   File  lusol2 LUSOL heap management routines
   Hbuild   Hchange  Hdelete  Hdown    Hinsert  Hup
   Heap-management routines for LUSOL's lu1fac.
   May be useful for other applications.
   ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   For LUSOL, the heap structure involves three arrays of length N.
   N        is the current number of entries in the heap.
   Ha(1:N)  contains the values that the heap is partially sorting.
            For LUSOL they are double precision values -- the largest
            element in each remaining column of the updated matrix.
            The biggest entry is in Ha(1), the top of the heap.
   Hj(1:N)  contains column numbers j.
            Ha(k) is the biggest entry in column j = Hj(k).
   Hk(1:N)  contains indices within the heap.  It is the
            inverse of Hj(1:N), so  k = Hk(j)  <=>  j = Hj(k).
            Column j is entry k in the heap.
   hops     is the number of heap operations,
            i.e., the number of times an entry is moved
            (the number of "hops" up or down the heap).
   Together, Hj and Hk let us find values inside the heap
   whenever we want to change one of the values in Ha.
   For other applications, Ha may need to be some other data type,
   like the keys that sort routines operate on.
   ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
   */

namespace lusol
{

/* ==================================================================
   Hdown  updates heap by moving down tree from node k.
   ================================================================== */
template<class Val>
void heap<Val>::HDOWN(INT N, INT K)
{
    INT J, JJ, JV, N2;
    Val V;

    V   = HA[K];
    JV  = HJ[K];
    N2  = N/2;

  x100:
    
    if(K>N2)
        goto x200;
    
    J = K+K;
    
    if(J<N) 
    {
        if(HA[J]<HA[J+1])
            J++;
    }

    if(V>=HA[J])
        goto x200;

    HA[K]   = HA[J];
    JJ      = HJ[J];
    HJ[K]   = JJ;
    HK[JJ]  = K;
    K       = J;

    goto x100;

  x200:
    HA[K]   = V;
    HJ[K]   = JV;
    HK[JV]  = K;
};

/* ==================================================================
   Hup updates heap by moving up tree from node k.
   ================================================================== */
template<class Val>
void heap<Val>::HUP(INT K)
{
    INT J, JV, K2;
    Val V;

    V   = HA[K];
    JV  = HJ[K];

x100:
    if(K<2)
        goto x200;

    K2  = K/2;

    if(V<HA[K2])
        goto x200;
  
    HA[K]   = HA[K2];
    J       = HJ[K2];
    HJ[K]   = J;
    HK[J]   = K;
    K       = K2;

    goto x100;

x200:
    HA[K]   = V;
    HJ[K]   = JV;
    HK[JV]  = K;
}

/* ==================================================================
   Hinsert inserts (v,jv) into heap of length N-1
   to make heap of length N.
   ================================================================== */
template<class Val>
void heap<Val>::HINSERT(INT N, Val V, INT JV)
{
    HA[N]   = V;
    HJ[N]   = JV;
    HK[JV]  = N;
    HUP(N);
}

/* ==================================================================
   Hchange changes Ha(k) to v in heap of length N.
   ================================================================== */
template<class Val>
void heap<Val>::change(INT N, INT K, Val V, INT JV)
{
    Val V1;

    V1 = HA[K];
    HA[K] = V;
    HJ[K] = JV;
    HK[JV] = K;
    if(V1<V)
        HUP(K);
    else
        HDOWN(N,K);
}

/* ==================================================================
   Hdelete deletes Ha(k) from heap of length N.
   ================================================================== */
template<class Val>
void heap<Val>::remove(INT *N, INT K)
{
    INT JV, NX;
    Val V;

    NX  = *N;
    V   = HA[NX];
    JV  = HJ[NX];
    (*N)--;
    
    if(K<NX)
        change(NX,K,V,JV);
}

/* ==================================================================
   Hbuild initializes the heap by inserting each element of Ha.
   Input:  Ha, Hj.
   Output: Ha, Hj, Hk, hops.
   ================================================================== */
template<class Val>
void heap<Val>::build(INT N)
{
    INT JV, K, KK;
    Val V;

    for(K = 1; K <= N; K++) 
    {
        KK  = K;
        V   = HA[K];
        JV  = HJ[K];
        HINSERT(KK,V,JV);
    }
};

template<class Val>
bool heap<Val>::realloc(INT newsize,INT oldsize)
{
    HA = (Val *) clean_realloc(HA,sizeof(Val), newsize, oldsize);
    HJ = (INT *) clean_realloc(HJ,sizeof(INT), newsize, oldsize);
    HK = (INT *) clean_realloc(HK,sizeof(INT), newsize, oldsize);

    if( newsize > 0 && (HA == nullptr || HJ == nullptr || HK == nullptr ) )
    {
        return false;
    }
    else
    {
        return true;
    };
};

template<class Val>
heap<Val>::heap()
:HA(nullptr),HJ(nullptr),HK(nullptr)
{};

template<class Val>
void heap<Val>::clear(INT len)
{
    memclear(HA,  len);
    memclear(HJ,  len);
    memclear(HK,  len);
};    

template<class Val>
void heap<Val>::set(INT HLEN, INT J, Val val)
{
    HA[HLEN]    = val;
    HJ[HLEN]    = J;
    HK[J]       = HLEN;
};

template struct heap<REAL>;
template struct heap<FLOAT>;

};
