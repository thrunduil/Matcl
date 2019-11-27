/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   LUSOL routines from the Stanford Optimization Laboratory
   The parts included are:
    lusol1      Factor a given matrix A from scratch (lu1fac).
    lusol2      Heap-management routines for lu1fac.
    lusol6      Solve with the current LU factors.
    lusol7      Utilities for all update routines.
    lusol8      Replace a column (Bartels-Golub update).
   ------------------------------------------------------------------
   26 Apr 2002: TCP implemented using heap data structure.
   01 May 2002: lu1DCP implemented.
   07 May 2002: lu1mxc must put 0.0 at top of empty columns.
   09 May 2002: lu1mCP implements Markowitz with cols searched
                in heap order.
                Often faster (searching 20 or 40 cols) but more dense.
   11 Jun 2002: TRP implemented.
                lu1mRP implements Markowitz with Threshold Rook
                Pivoting.
                lu1mxc maintains max col elements  (was lu1max.)
                lu1mxr maintains max row elements.
   12 Jun 2002: lu1mCP seems too slow on big problems (e.g. memplus).
                Disabled it for the moment.  (Use lu1mar + TCP.)
   14 Dec 2002: TSP implemented.
                lu1mSP implements Markowitz with TSP.
   07 Mar 2003: character*1, character*2 changed to f90 form.
                Comments changed from * in column to ! in column 1.
                Comments kept within column 72 to avoid compiler
                warning.
   06 Mar 2004: Translation to C by Kjell Eikland with the addition
                of data wrappers, parametric constants, various 
                helper routines, and dynamic memory reallocation.
   26 May 2004: Added LUSOL_IP_UPDATELIMIT parameter and provided
                for dynamic memory expansion based on possible 
                forward requirements.
   08 Jul 2004: Revised logic in lu6chk based on new code from 
                Michael Saunders.
   01 Dec 2005: Add support for CMEX interface (disable by undef MATLAB)
                Also include various bug fixes (disable by undef YZHANG)
                Yin Zhang <yzhang@cs.utexas.edu>
   ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */

#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>
#include <string.h>
#include <float.h>
#include <math.h>
#include "lusol.h"
#include "myblas.h"

#include <assert.h>

// LUSOL Object creation and destruction

#ifdef __GNUC__
namespace {
    #define vsprintf_s vsnprintf
}
#endif

namespace lusol
{

const REAL DENSE = .5;

const REAL lusol_user_params<REAL>::LUSOL_DEFAULT_GAMMA         = 2.0;
const REAL lusol_user_params<REAL>::LUSOL_SMALLNUM              = 1.0e-20;  // IAEE doubles have precision 2.22e-16
const REAL lusol_user_params<REAL>::LUSOL_BIGNUM                = 1.0e+20;

const REAL lusol_default_params<REAL>::LUSOL_RP_ZEROTOLERANCE   = 3.0e-13;
const REAL lusol_default_params<REAL>::LUSOL_RP_SMALLDIAG_U     = 3.7e-11;
const REAL lusol_default_params<REAL>::LUSOL_RP_EPSDIAG_U       = 3.7e-11;

const FLOAT lusol_user_params<FLOAT>::LUSOL_DEFAULT_GAMMA       = 2.0f;
const FLOAT lusol_user_params<FLOAT>::LUSOL_SMALLNUM            = 1.1921e-07f * 2.0e-4f;
const FLOAT lusol_user_params<FLOAT>::LUSOL_BIGNUM              = 1.1921e+07f * 2.0e+4f;

const FLOAT lusol_default_params<FLOAT>::LUSOL_RP_ZEROTOLERANCE = 1.0e-06f;
const FLOAT lusol_default_params<FLOAT>::LUSOL_RP_SMALLDIAG_U   = 1.0e-05f;
const FLOAT lusol_default_params<FLOAT>::LUSOL_RP_EPSDIAG_U     = 1.0e-05f;

void * clean_realloc(void *oldptr, INT width, INT newsize, INT oldsize)
{
    newsize *= width;
    oldsize *= width;
    oldptr = realloc(oldptr, newsize);

    if (newsize != 0 && oldptr == nullptr)
    {
        throw std::bad_alloc();
    }

    if(newsize > oldsize)
        memset((char *)oldptr+oldsize, '\0', newsize-oldsize);
    return oldptr;
};
template<class T>
bool LUSOLrec<T>::LUSOL_realloc_a(INT newsize)
{
    int oldsize;

    if(newsize < 0)
        newsize = this->lena + maximum((INT)abs(newsize), LUSOL_MINDELTA_a);

    oldsize = this->lena;
    this->lena = newsize;
    if(newsize > 0)
        newsize++;
    if(oldsize > 0)
        oldsize++;

    this->a    = (T *)    clean_realloc(this->a,    sizeof(*(this->a)),
                                                    newsize, oldsize);
    this->indc = (INT *)  clean_realloc(this->indc, sizeof(*(this->indc)),
                                                    newsize, oldsize);
    this->indr = (INT *)  clean_realloc(this->indr, sizeof(*(this->indr)),
                                                    newsize, oldsize);
    if((newsize == 0) || (this->a != nullptr && this->indc != nullptr && this->indr != nullptr))
        return true;
    else
        return false;
};
template<class T>
bool LUSOLrec<T>::LUSOL_expand_a(INT *delta_lena, INT *right_shift)
{
    INT LENA, NFREE, LFREE;
  
    // Add expansion factor to avoid having to resize too often/too much;
    // (exponential formula suggested by Michael A. Saunders)
    LENA = this->lena;
    *delta_lena = (INT) ((*delta_lena) * pow(1.5, abs((double)(*delta_lena))/LENA));

    // XXX: the exponential formula is too aggressive for large A
    if (*delta_lena > LENA/3 && *delta_lena > 500000)
        *delta_lena = LENA/3;
  
    // Expand it!
    if(*delta_lena <= 0 || !LUSOL_realloc_a(LENA+(*delta_lena)))
        return false;

    // Make sure we return the actual memory increase of a
    *delta_lena = this->lena-LENA;

    // Shift the used memory area to the right
    LFREE = *right_shift;
    NFREE = LFREE+*delta_lena;
    LENA  -= LFREE-1;
    mem_move(this->a+NFREE,    this->a+LFREE,    LENA);
    mem_move(this->indr+NFREE, this->indr+LFREE, LENA);
    mem_move(this->indc+NFREE, this->indc+LFREE, LENA);
    
    // Also return the new starting position for the used memory area of a
    *right_shift  = NFREE;

    this->expanded_a++;
    return true;
};
template<class T>
bool LUSOLrec<T>::LUSOL_realloc_r(INT newsize)
{
    int oldsize;

    oldsize = this->maxm;
    this->maxm = newsize;
    if(newsize > 0)
        newsize++;
    if(oldsize > 0)
        oldsize++;

    this->lenr  = (INT *) clean_realloc(this->lenr,  sizeof(*(this->lenr)), 
                                                     newsize, oldsize);
    this->ip    = (INT *) clean_realloc(this->ip,    sizeof(*(this->ip)), 
                                                     newsize, oldsize);
    this->iqloc = (INT *) clean_realloc(this->iqloc, sizeof(*(this->iqloc)), 
                                                     newsize, oldsize);
    this->ipinv = (INT *) clean_realloc(this->ipinv, sizeof(*(this->ipinv)), 
                                                     newsize, oldsize);
    this->m_markr = (INT*)clean_realloc(this->m_markr, sizeof(*(this->m_markr)), 
                                                     newsize, oldsize);
    this->locr  = (INT *) clean_realloc(this->locr,  sizeof(*(this->locr)), 
                                                     newsize, oldsize);

    if(newsize == 0 || 
        (this->lenr != nullptr && this->ip != nullptr && this->iqloc != nullptr &&
        this->ipinv != nullptr && this->locr != nullptr && this->m_markr != nullptr)) 
    {

        #ifdef AlwaysSeparateHamaxR
            if(this->luparm[LUSOL_IP_PIVOTTYPE] == LUSOL_PIVOT_TRP)
        #endif
        {
            this->amaxr = (TR*) clean_realloc(this->amaxr, sizeof(*(this->amaxr)),
                                                          newsize, oldsize);
            if(newsize > 0 && this->amaxr == nullptr)
                return false;
        }
        return true;
    }
    else
        return false;
};
template<class T>
bool LUSOLrec<T>::LUSOL_realloc_c(INT newsize)
{
    int oldsize;

    oldsize = this->maxn;
    this->maxn = newsize;
    if(newsize > 0)
        newsize++;
    if(oldsize > 0)
        oldsize++;

    this->lenc  = (INT *)  clean_realloc(this->lenc,  sizeof(*(this->lenc)),
                                                      newsize, oldsize);
    this->iq    = (INT *)  clean_realloc(this->iq,    sizeof(*(this->iq)),
                                                      newsize, oldsize);
    this->iploc = (INT *)  clean_realloc(this->iploc, sizeof(*(this->iploc)),
                                                      newsize, oldsize);
    this->iqinv = (INT *)  clean_realloc(this->iqinv, sizeof(*(this->iqinv)),
                                                      newsize, oldsize);
    this->m_cols = (INT *)  clean_realloc(this->m_cols, sizeof(*(this->m_cols)),
                                                      newsize, oldsize);
    this->m_markc = (INT *)  clean_realloc(this->m_markc, sizeof(*(this->m_markc)),
                                                      newsize, oldsize);
    this->locc  = (INT *)  clean_realloc(this->locc,  sizeof(*(this->locc)),
                                                      newsize, oldsize);
    this->work  = (T *)    clean_realloc(this->work,  sizeof(*(this->work)),
                                                      newsize, oldsize);

    this->wc    = this->work;
    this->wr    = (TR*)this->work;

    #ifdef LUSOLSafeFastUpdate
        this->vLU6L = (T *) clean_realloc(this->vLU6L, sizeof(*(this->vLU6L)),
                                                      newsize, oldsize);
    #else
        this->vLU6L = this->w;
    #endif

    if(newsize == 0 ||
        (this->work != nullptr && this->lenc != nullptr &&
        this->iq != nullptr && this->iploc != nullptr &&
        this->iqinv != nullptr && this->locc != nullptr &&
        this->m_cols != nullptr && this->m_markc != nullptr)) 
    {

        if(this->luparm[LUSOL_IP_PIVOTTYPE] == LUSOL_PIVOT_TCP) 
        {
            if (!this->m_heap.realloc(newsize,oldsize))
            {
                return false;
            };
        }
        #ifndef ClassicdiagU
            if(this->luparm[LUSOL_IP_KEEPLU] == false) 
            {
                this->diagU = (T *) clean_realloc(this->diagU, sizeof(*(this->diagU)),
                                                          newsize, oldsize);
                if(newsize > 0 && this->diagU == nullptr)
                    return false;
            }
        #endif
      
        return true;
    }
    else
        return false;
};
template<class T>
lusol::LUSOLrec<T> *LUSOLrec<T>::LUSOL_create(FILE *outstream, INT msgfil, INT pivotmodel, INT updatelimit)
{
    LUSOLrec *newLU;

    newLU = (LUSOLrec *) malloc(sizeof(*newLU));    
  
    if(newLU == nullptr)
        return newLU;

    memset(newLU, 0, sizeof(*newLU));

    newLU->luparm[LUSOL_IP_SCALAR_NZA]       = LUSOL_MULT_nz_a;
    newLU->outstream                         = outstream;
    newLU->luparm[LUSOL_IP_PRINTUNIT]        = msgfil;
    newLU->luparm[LUSOL_IP_PRINTLEVEL]       = LUSOL_MSG_SINGULARITY; 

    newLU->LUSOL_setpivotmodel(pivotmodel);

    newLU->parmlu[LUSOL_RP_GAMMA]            = lusol_user_params<TR>::LUSOL_DEFAULT_GAMMA;

    newLU->parmlu[LUSOL_RP_ZEROTOLERANCE]    = lusol_default_params<TR>::LUSOL_RP_ZEROTOLERANCE;

    newLU->parmlu[LUSOL_RP_SMALLDIAG_U]      = lusol_default_params<TR>::LUSOL_RP_SMALLDIAG_U;
    newLU->parmlu[LUSOL_RP_EPSDIAG_U]        = lusol_default_params<TR>::LUSOL_RP_EPSDIAG_U;

    newLU->parmlu[LUSOL_RP_COMPSPACE_U]      = TR(3.0e+0);

    newLU->luparm[LUSOL_IP_MARKOWITZ_MAXCOL] = 5;
    newLU->parmlu[LUSOL_RP_MARKOWITZ_CONLY]  = TR(0.3e+0);
    newLU->parmlu[LUSOL_RP_MARKOWITZ_DENSE]  = DENSE;

    newLU->luparm[LUSOL_IP_KEEPLU]           = true;
    newLU->luparm[LUSOL_IP_UPDATELIMIT]      = updatelimit;
  
    return newLU;
};
template<class T>
const char *LUSOLrec<T>::LUSOL_pivotLabel()
{
    static const char *pivotText[LUSOL_PIVOT_MAX+1] = {"TPP", "TRP", "TCP", "TSP"};
    return pivotText[this->luparm[LUSOL_IP_PIVOTTYPE]];
}
template<class T>
void LUSOLrec<T>::LUSOL_setpivotmodel(INT pivotmodel)
{
    if(pivotmodel <= LUSOL_PIVOT_DEFAULT || pivotmodel > LUSOL_PIVOT_MAX)
        pivotmodel = LUSOL_PIVOT_TPP;
    this->luparm[LUSOL_IP_PIVOTTYPE] = pivotmodel;

    // Set default pivot tolerances (note that UPDATEMAX should always be <= FACTORMAX)
    if(pivotmodel == LUSOL_PIVOT_TPP) 
    {
        this->parmlu[LUSOL_RP_FACTORMAX_Lij]    = 100.0;
        this->parmlu[LUSOL_RP_UPDATEMAX_Lij]    =  10.0;
    }
    else if(pivotmodel == LUSOL_PIVOT_TRP) 
    {
        this->parmlu[LUSOL_RP_FACTORMAX_Lij]    =  5.0;
        this->parmlu[LUSOL_RP_UPDATEMAX_Lij]    =  5.0;
    }
    else 
    {
        this->parmlu[LUSOL_RP_FACTORMAX_Lij]    = 10.0;
        this->parmlu[LUSOL_RP_UPDATEMAX_Lij]    = 10.0;
    }
};
template<class T>
bool LUSOLrec<T>::LUSOL_tightenpivot()
{
    if(minimum(this->parmlu[LUSOL_RP_UPDATEMAX_Lij], this->parmlu[LUSOL_RP_UPDATEMAX_Lij]) < 1.1)
        return false;
             
    this->parmlu[LUSOL_RP_FACTORMAX_Lij] = TR(1.0) + this->parmlu[LUSOL_RP_FACTORMAX_Lij]/TR(3.0);
    this->parmlu[LUSOL_RP_UPDATEMAX_Lij] = TR(1.0) + this->parmlu[LUSOL_RP_UPDATEMAX_Lij]/TR(3.0);

    return true;
};

template<class T>
const char * LUSOLrec<T>::LUSOL_informstr(INT inform)
{
    static const char *informText[LUSOL_INFORM_MAX-LUSOL_INFORM_MIN+1] = 
    {
        "Lost rank",
        "Success",
        "Singular A",
        "Unstable factorization",
        "Row or column count exceeded",
        "Duplicate A matrix entry found",
        "",
        "",
        "Insufficient memory for factorization",
        "Fatal internal error",
        "Found no suitable pivot",
    };
    if(inform < LUSOL_INFORM_MIN || inform > LUSOL_INFORM_MAX)
        inform = this->luparm[LUSOL_IP_INFORM];

    return informText[inform-LUSOL_INFORM_MIN];
};
template<class T>
void LUSOLrec<T>::LUSOL_clear(bool nzonly)
{
    int len;

    this->nelem = 0;
    if(!nzonly) 
    {

        // lena arrays
        len = this->lena + LUSOL_ARRAYOFFSET;
        memclear(this->a,    len);
        memclear(this->indc, len);
        memclear(this->indr, len);

        // maxm arrays
        len = this->maxm + LUSOL_ARRAYOFFSET;
        memclear(this->lenr,  len);
        memclear(this->ip,    len);
        memclear(this->iqloc, len);
        memclear(this->ipinv, len);
        memclear(this->m_markc, len);
        memclear(this->locr,  len);

        if(this->amaxr != nullptr
            #ifdef AlwaysSeparateHamaxR
                && this->luparm[LUSOL_IP_PIVOTTYPE] == LUSOL_PIVOT_TRP
            #endif
            )
        {
            memclear(this->amaxr, len);
        };

        // maxn arrays
        len = this->maxn + LUSOL_ARRAYOFFSET;
        memclear(this->lenc,  len);
        memclear(this->iq,    len);
        memclear(this->iploc, len);
        memclear(this->iqinv, len);
        memclear(this->m_cols, len);
        memclear(this->m_markc, len);
        memclear(this->locc,  len);
        memclear(this->work,  len);

        if(this->luparm[LUSOL_IP_PIVOTTYPE] == LUSOL_PIVOT_TCP) 
        {
            this->m_heap.clear(len);
        }

        #ifndef ClassicdiagU
            if(this->luparm[LUSOL_IP_KEEPLU] == false) 
            {
                memclear(this->diagU, len);
            }
        #endif      
    }
};
template<class T>
bool LUSOLrec<T>::LUSOL_assign(const Integer iA[], const Integer jA[], const T Aij[], 
                               INT mmax, INT nmax, INT cols)
{
    Integer nzcount = jA[cols] - jA[0];

    // Adjust the size of the a structure
    if(!LUSOL_realloc_a(nzcount*this->luparm[LUSOL_IP_SCALAR_NZA]))
    {
        return false;
    };
    if(!LUSOL_realloc_r(mmax))
    {
        return false;
    };
    if(!LUSOL_realloc_c(nmax))
    {
        return false;
    };

    Integer pos = 1;
    for(Integer j = 0; j < cols; j++) 
    {
        for (Integer k = jA[j]; k < jA[j+1]; ++k, ++ pos)
        {
            this->indc[pos] = iA[k] + 1;
            this->indr[pos] = j + 1;
            this->a[pos]    = Aij[k];
        };
    }
    this->m = mmax;
    this->n = nmax;
    this->nelem = nzcount;
    return true;
};
template<class T>
INT LUSOLrec<T>::LUSOL_loadColumn(INT iA[], INT jA, T Aij[], INT nzcount)
{
    INT i, nz, k;

    nz = this->nelem;
    i = nz + nzcount;
    if(i > (this->lena/this->luparm[LUSOL_IP_SCALAR_NZA]) &&
        !LUSOL_realloc_a(i*this->luparm[LUSOL_IP_SCALAR_NZA]))
    {
        return -1;
    };

    k = 0;
    for(i = 1; i <= nzcount; i++) 
    {
        if(Aij[i] == T())
            continue;
        if(iA[i] <= 0 || iA[i] > this->m ||  jA <= 0 || jA > this->n) 
        {
            //FIXME:
            LUSOL_report(this,0,(char*) "Variable index outside of set bounds (r:%d/%d, c:%d/%d)\n",
                             iA[i], this->m, jA, this->n);
            continue;
        }
        k++;
        nz++;
        this->a[nz]    = Aij[i];
        this->indc[nz] = iA[i];
        this->indr[nz] = jA;
    }
    this->nelem = nz;
    return k;
}
template<class T>
void LUSOLrec<T>::LUSOL_free()
{
    LUSOL_realloc_a(0);
    LUSOL_realloc_r(0);
    LUSOL_realloc_c(0);
    #ifdef UseRowBasedL0
        if(this->L0 != nullptr)
            LUSOL_matfree(&(this->L0));
    #endif
    free(this);
};
template<class T>
void LUSOLrec<T>::LUSOL_report(LUSOLrec* LUSOL, INT msglevel, char *format, ...)
{
    va_list ap;

    va_start(ap, format);
    if(LUSOL == nullptr) 
    {
        vfprintf(stderr, format, ap);
    }
    else if(msglevel >= 0) 
    {
        if(LUSOL->writelog != nullptr) 
        {
            char buff[255];

            vsprintf_s(buff, 255, format, ap);
            LUSOL->writelog(LUSOL, LUSOL->loghandle, buff);
        }
        if(LUSOL->outstream != nullptr) 
        {
            vfprintf(LUSOL->outstream, format, ap);
            fflush(LUSOL->outstream);
        }
    }
    va_end(ap);
};
template<class T>
void LUSOLrec<T>::LUSOL_timer(INT timerid, char *text)
{
    //FIXME:
    LUSOL_report(this, -1, (char*) "TimerID %d at %s - %s\n", timerid, "", text);
};
template<class T>
INT LUSOLrec<T>::LUSOL_ftran(T b[], INT NZidx[], bool prepareupdate)
{
    INT  inform;
    T   *vector;

    if(prepareupdate)
        vector = this->vLU6L;
    else
        vector = this->wc;

    // XXX: otherwise we get seg fault when vector is empty
    memcopy((vector+1), (b+1), this->n);

    LU6SOL(LUSOL_SOLVE_Aw_v, vector, b, NZidx, &inform);
    return inform;
};
template<class T>
INT LUSOLrec<T>::LUSOL_btran(T b[], INT NZidx[])
{
    INT inform;

    // XXX: otherwise we get seg fault when this->w is empty
    memcopy((this->wc+1), (b+1), this->m);

    LU6SOL(LUSOL_SOLVE_Atv_w, b, this->wc, NZidx, &inform);
    return inform;
};

template<class T>
INT LUSOLrec<T>::LUSOL_replaceColumn(INT jcol, T v[])
{ 
    INT  inform;
    TR   VNORM;
    T    DIAG;

    LU8RPC(LUSOL_UPDATE_OLDNONEMPTY, LUSOL_UPDATE_NEWNONEMPTY,
                jcol, v, nullptr,
                &inform, &DIAG, &VNORM);
    this->replaced_c++;
    return inform;
};

// The purpose of this routine is to find the slack row/column in
//   user-index that was singular in the last unsuccessful column
//   update; zero is returned if the search was unsuccessful.
//   (Source is Michael A. Saunders; private communication to KE)
template<class T>
INT LUSOLrec<T>::LUSOL_findColumnPosition(INT jcol)
{
    INT j;

    for(j = this->m; j > 0; j--)
    {
        if(this->iq[j] == jcol)
            break;
    };
    if(j > 0)
        j = this->ip[j];
    return j;
};

template<class T>
lusol::LUSOLmat<T> *LUSOLmat<T>::create(INT dim, INT nz)
{
    LUSOLmat *newm;

    newm = (LUSOLmat *)   malloc(sizeof(*newm));
    
    if(newm != nullptr) 
    {
        memset(newm, 0, sizeof(*newm));

        newm->a    = (T *)    malloc((nz+1)*sizeof(T));
        newm->vlen = (INT *)  malloc((dim+1)*sizeof(int));
        newm->indr = (INT *)  malloc((nz+1)*sizeof(int));
        newm->indc = (INT *)  malloc((nz+1)*sizeof(int));

        if((newm->a == nullptr) || (newm->vlen == nullptr) || 
            (newm->indr == nullptr) || (newm->indc == nullptr))
        {
            destroy(&newm);
        };
    }
    return newm;
};
template<class T>
void LUSOLmat<T>::destroy(LUSOLmat **mat)
{
    lusol_free((*mat)->a);
    lusol_free((*mat)->indc);
    lusol_free((*mat)->indr);
    lusol_free((*mat)->vlen);
    lusol_free(*mat);
    *mat = nullptr;
};

};

#include "lusol1.inl"
#include "lusol6a.inl"
#include "lusol7a.inl"
#include "lusol8a.inl"

#undef I_N_T
#undef INT

template struct lusol::LUSOLmat<lusol::REAL>;
template struct lusol::LUSOLmat<lusol::FLOAT>;
template struct lusol::LUSOLmat<lusol::COMPLEX>;
template struct lusol::LUSOLmat<lusol::FCOMPLEX>;

template struct lusol::LUSOLrec<lusol::REAL>;
template struct lusol::LUSOLrec<lusol::FLOAT>;
template struct lusol::LUSOLrec<lusol::COMPLEX>;
template struct lusol::LUSOLrec<lusol::FCOMPLEX>;
