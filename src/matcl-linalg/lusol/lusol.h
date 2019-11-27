#pragma once

#include <stdio.h>
#include "commonlib.h"
#include "heap.h"

namespace lusol
{

namespace details
{
    template<class Val> struct real_type            {};
    template<>          struct real_type<REAL>      { using type = REAL; };
    template<>          struct real_type<FLOAT>     { using type = FLOAT; };
    template<>          struct real_type<COMPLEX>   { using type = REAL; };
    template<>          struct real_type<FCOMPLEX>  { using type = FLOAT; };
    
};

//template<class Ret, class In>
//Ret lusol_cast(const In&);

// Performance compiler options                                              
// -------------------------------------------------------------------------
#define LUSOLSafeFastUpdate         // Use separate array for LU6L result storage
//#define AlwaysSeparateHamaxR      // Enabled when the pivot model is fixed
//#define UseRowBasedL0             // Create a row-sorted version of L0
//#define UseTimer


// General constants and data type definitions                               
// -------------------------------------------------------------------------
static const size_t LUSOL_ARRAYOFFSET   = 1;


// User-settable default parameter values                                    
// -------------------------------------------------------------------------
template<class Val>
struct lusol_user_params {};

template<class Val>
struct lusol_default_params{};

template<>
struct lusol_user_params<REAL>
{
    static const REAL LUSOL_DEFAULT_GAMMA;
    static const REAL LUSOL_SMALLNUM;
    static const REAL LUSOL_BIGNUM;
};

template<>
struct lusol_default_params<REAL>
{
    static const REAL LUSOL_RP_ZEROTOLERANCE;
    static const REAL LUSOL_RP_SMALLDIAG_U;
    static const REAL LUSOL_RP_EPSDIAG_U;
};

template<>
struct lusol_user_params<FLOAT>
{
    static const FLOAT LUSOL_DEFAULT_GAMMA;
    static const FLOAT LUSOL_SMALLNUM;
    static const FLOAT LUSOL_BIGNUM;
};

template<>
struct lusol_default_params<FLOAT>
{
    static const FLOAT LUSOL_RP_ZEROTOLERANCE;
    static const FLOAT LUSOL_RP_SMALLDIAG_U;
    static const FLOAT LUSOL_RP_EPSDIAG_U;
};

static const INT LUSOL_MINDELTA_a       = 10000;
static const INT LUSOL_MULT_nz_a        = 2;        // Could consider 6 or 7


// Fixed system parameters (changeable only by developers)                   
// -------------------------------------------------------------------------

// parmlu INPUT parameters:
//static const INT LUSOL_RP_USERDATA_0        = 0;
static const INT LUSOL_RP_FACTORMAX_Lij     = 1;
static const INT LUSOL_RP_UPDATEMAX_Lij     = 2;
static const INT LUSOL_RP_ZEROTOLERANCE     = 3;
static const INT LUSOL_RP_SMALLDIAG_U       = 4;
static const INT LUSOL_RP_EPSDIAG_U         = 5;
static const INT LUSOL_RP_COMPSPACE_U       = 6;
static const INT LUSOL_RP_MARKOWITZ_CONLY   = 7;
static const INT LUSOL_RP_MARKOWITZ_DENSE   = 8;
static const INT LUSOL_RP_GAMMA             = 9;

// parmlu OUPUT parameters:
static const INT LUSOL_RP_MAXELEM_A         = 10;
static const INT LUSOL_RP_MAXMULT_L         = 11;
static const INT LUSOL_RP_MAXELEM_U         = 12;
static const INT LUSOL_RP_MAXELEM_DIAGU     = 13;
static const INT LUSOL_RP_MINELEM_DIAGU     = 14;
static const INT LUSOL_RP_MAXELEM_TCP       = 15;
static const INT LUSOL_RP_GROWTHRATE        = 16;
static const INT LUSOL_RP_USERDATA_1        = 17;
static const INT LUSOL_RP_USERDATA_2        = 18;
static const INT LUSOL_RP_USERDATA_3        = 19;
static const INT LUSOL_RP_RESIDUAL_U        = 20;
static const INT LUSOL_RP_LASTITEM          = LUSOL_RP_RESIDUAL_U;

// luparm INPUT parameters:
static const INT LUSOL_IP_USERDATA_0        = 0;
static const INT LUSOL_IP_PRINTUNIT         = 1;
static const INT LUSOL_IP_PRINTLEVEL        = 2;
static const INT LUSOL_IP_MARKOWITZ_MAXCOL  = 3;
static const INT LUSOL_IP_SCALAR_NZA        = 4;
static const INT LUSOL_IP_UPDATELIMIT       = 5;
static const INT LUSOL_IP_PIVOTTYPE         = 6;
static const INT LUSOL_IP_USERDATA_1        = 7;
static const INT LUSOL_IP_KEEPLU            = 8;
static const INT LUSOL_IP_USERDATA_2        = 9;

// luparm OUTPUT parameters:
static const INT LUSOL_IP_INFORM            = 10;
static const INT LUSOL_IP_SINGULARITIES     = 11;
static const INT LUSOL_IP_SINGULARINDEX     = 12;
static const INT LUSOL_IP_MINIMUMLENA       = 13;
static const INT LUSOL_IP_MAXLEN            = 14;
static const INT LUSOL_IP_UPDATECOUNT       = 15;
static const INT LUSOL_IP_RANK_U            = 16;
static const INT LUSOL_IP_COLCOUNT_DENSE1   = 17;
static const INT LUSOL_IP_COLCOUNT_DENSE2   = 18;
static const INT LUSOL_IP_COLINDEX_DUMIN    = 19;
static const INT LUSOL_IP_COLCOUNT_L0       = 20;
static const INT LUSOL_IP_NONZEROS_L0       = 21;
static const INT LUSOL_IP_NONZEROS_U0       = 22;
static const INT LUSOL_IP_NONZEROS_L        = 23;
static const INT LUSOL_IP_NONZEROS_U        = 24;
static const INT LUSOL_IP_NONZEROS_ROW      = 25;
static const INT LUSOL_IP_COMPRESSIONS_LU   = 26;
static const INT LUSOL_IP_MARKOWITZ_MERIT   = 27;
static const INT LUSOL_IP_TRIANGROWS_U      = 28;
static const INT LUSOL_IP_TRIANGROWS_L      = 29;
static const INT LUSOL_IP_USERDATA_3        = 30;
static const INT LUSOL_IP_LASTITEM          = LUSOL_IP_USERDATA_3;


// Parameter/option defines                                                  
// -------------------------------------------------------------------------
static const INT LUSOL_MSG_NONE             = -1;
static const INT LUSOL_MSG_SINGULARITY      =  0;
static const INT LUSOL_MSG_STATISTICS       = 10;
static const INT LUSOL_MSG_PIVOT            = 50;

static const INT LUSOL_PIVOT_DEFAULT        = -1;  // Set pivoting model to default
static const INT LUSOL_PIVOT_TPP            =  0;  // Threshold Partial   pivoting (normal)
static const INT LUSOL_PIVOT_TRP            =  1;  // Threshold Rook      pivoting
static const INT LUSOL_PIVOT_TCP            =  2;  // Threshold Complete  pivoting
static const INT LUSOL_PIVOT_TSP            =  3;  // Threshold Symmetric pivoting
static const INT LUSOL_PIVOT_MAX            =  LUSOL_PIVOT_TSP;

static const INT LUSOL_UPDATE_OLDEMPTY      =  0;  // No/empty current column.
static const INT LUSOL_UPDATE_OLDNONEMPTY   =  1;  // Current column need not have been empty.
static const INT LUSOL_UPDATE_NEWEMPTY      =  0;  // New column is taken to be zero.
static const INT LUSOL_UPDATE_NEWNONEMPTY   =  1;  // v(*) contains the new column;
                                                   // on exit,  v(*)  satisfies  L*v = a(new).
static const INT LUSOL_UPDATE_USEPREPARED   =  2;  // v(*)  must satisfy  L*v = a(new).

static const INT LUSOL_SOLVE_Lv_v           =  1;  // v  solves   L v = v(input). w  is not touched.
static const INT LUSOL_SOLVE_Ltv_v          =  2;  // v  solves   L'v = v(input). w  is not touched.
static const INT LUSOL_SOLVE_Uw_v           =  3;  // w  solves   U w = v.        v  is not altered.
static const INT LUSOL_SOLVE_Utv_w          =  4;  // v  solves   U'v = w.        w  is destroyed.
static const INT LUSOL_SOLVE_Aw_v           =  5;  // w  solves   A w = v.        v  is altered as in 1.
static const INT LUSOL_SOLVE_Atv_w          =  6;  // v  solves   A'v = w.        w  is destroyed.

// If mode = 3,4,5,6, v and w must not be the same arrays.
// If lu1fac has just been used to factorize a symmetric matrix A
// (which must be definite or quasi-definite), the factors A = L U
//  may be regarded as A = LDL', where D = diag(U).  In such cases,
//  the following (faster) solve codes may be used:
static const INT LUSOL_SOLVE_Av_v           =  7;  // v  solves   A v = L D L'v = v(input). w  is not touched.
static const INT LUSOL_SOLVE_LDLtv_v        =  8;  // v  solves       L |D| L'v = v(input). w  is not touched.

static const INT LUSOL_INFORM_RANKLOSS      = -1;
static const INT LUSOL_INFORM_LUSUCCESS     =  0;
static const INT LUSOL_INFORM_LUSINGULAR    =  1;
static const INT LUSOL_INFORM_LUUNSTABLE    =  2;
static const INT LUSOL_INFORM_ADIMERR       =  3;
static const INT LUSOL_INFORM_ADUPLICATE    =  4;
static const INT LUSOL_INFORM_ANEEDMEM      =  7;  // Set lena >= luparm[LUSOL_IP_MINIMUMLENA]
static const INT LUSOL_INFORM_FATALERR      =  8;
static const INT LUSOL_INFORM_NOPIVOT       =  9;  // No diagonal pivot found with TSP or TDP.

static const INT LUSOL_INFORM_MIN           =  LUSOL_INFORM_RANKLOSS;
static const INT LUSOL_INFORM_MAX           =  LUSOL_INFORM_NOPIVOT;

static const INT LUSOL_INFORM_GETLAST       = 10;  // Code for LUSOL_informstr.
static const INT LUSOL_INFORM_SERIOUS       =  LUSOL_INFORM_LUUNSTABLE;


// Prototypes for call-back functions
typedef void LUSOLlogfunc(void *lp, void *userhandle, char *buf);


// Sparse matrix data
template<class T>
struct LUSOLmat 
{
    T    *a;
    INT  *vlen, *indr, *indc;

    static LUSOLmat *  create(INT dim, INT nz);
    static void        destroy(LUSOLmat **mat);
};

// The main LUSOL data record
template<class T>
struct LUSOLrec 
{
    using TR    = typename details::real_type<T>::type;

    // General data
    FILE *          outstream;           // Output stream, initialized to STDOUT
    LUSOLlogfunc *  writelog;
    void *          loghandle;
    LUSOLlogfunc *  debuginfo;

    // Parameter storage arrays
    INT             luparm[LUSOL_IP_LASTITEM + 1];
    TR              parmlu[LUSOL_RP_LASTITEM + 1];

    // Arrays of length lena+1
    INT             lena, nelem;
    INT             *indc, *indr;
    T               *a;

    // Arrays of length maxm+1 (row storage)
    INT             maxm, m;
    INT *           ip;
    INT             *lenr, *iqloc, *ipinv, *locr;
    INT             *m_markr;

    // Arrays of length maxn+1 (column storage)
    INT             maxn, n;
    INT *           iq;
    INT             *lenc, *iploc, *iqinv, *locc;
    INT             *m_cols, *m_markc;
    T               *work;
    TR              *wr;
    T               *wc;
    T               *vLU6L;

    // Extra arrays of length n for TCP and keepLU == false
    T               *diagU;

    heap<TR>        m_heap;
  
    //              Extra arrays of length m for TRP
    TR *            amaxr;

    // Extra storage arrays for LU transposed fast solution routines
    #ifdef UseRowBasedL0
        LUSOLmat *  L0;
    #endif

    // Miscellaneous data
    INT             expanded_a;
    INT             replaced_c;
    INT             replaced_r;

    //--------------------------------------------------------------------------
    //                      lusol.cpp
    //--------------------------------------------------------------------------
    bool        LUSOL_realloc_a(INT newsize);
    bool        LUSOL_expand_a(INT *delta_lena, INT *right_shift);
    bool        LUSOL_realloc_r(INT newsize);
    bool        LUSOL_realloc_c(INT newsize);
    const char* LUSOL_pivotLabel();
    void        LUSOL_setpivotmodel(INT pivotmodel);
    bool        LUSOL_tightenpivot();
    const char* LUSOL_informstr(INT inform);
    void        LUSOL_clear(bool nzonly);
    bool        LUSOL_assign(const Integer iA[], const Integer jA[], const T Aij[], 
                            INT mmax, INT nmax, INT cols);
    INT         LUSOL_loadColumn(INT iA[], INT jA, T Aij[], INT nzcount);
    void        LUSOL_free();
    static void LUSOL_report(LUSOLrec* LUSOL, INT msglevel, char *format, ...);
    void        LUSOL_timer(INT timerid, char *text);
    INT         LUSOL_ftran(T b[], INT NZidx[], bool prepareupdate);
    INT         LUSOL_btran(T b[], INT NZidx[]);
    INT         LUSOL_replaceColumn(INT jcol, T v[]);
    INT         LUSOL_findColumnPosition(INT jcol);
    static LUSOLrec*   LUSOL_create(FILE *outstream, INT msgfil, INT pivotmodel, INT updatelimit);

    //--------------------------------------------------------------------------
    //                      lusol1.cpp
    //--------------------------------------------------------------------------
    void LU1PQ1(INT M, INT N, INT LEN[],
                INT IPERM[], INT LOC[], INT INV[], INT NUM[]);
    void LU1PQ2(INT NZPIV, INT *NZCHNG,
                INT IND[], INT LENOLD[], INT LENNEW[], INT IXLOC[], INT IX[], INT IXINV[]);
    void LU1PQ3(INT MN, INT LEN[], INT IPERM[], INT IW[], INT *NRANK);
    void LU1REC(INT N, bool REALS, INT *LTOP,
                                 INT IND[], INT LEN[], INT LOC[]);
    void LU1SLK();
    void LU1GAU(INT MELIM, INT NSPARE,
                TR SMALL, INT LPIVC1, INT LPIVC2, INT *LFIRST, INT LPIVR2,
                INT LFREE, INT MINFRE, INT ILAST, INT *JLAST, INT *LROW, INT *LCOL,
                INT *LU, INT *NFILL,
                INT MARK[],  T AL[], INT MARKL[], T AU[], INT IFILL[], INT JFILL[]);
    void LU1MAR(INT MAXMN, bool TCP, TR AIJTOL, TR LTOL,
                INT MAXCOL, INT MAXROW, INT *IBEST, INT *JBEST, INT *MBEST);
    void LU1MRP(INT MAXMN, TR LTOL, INT MAXCOL, INT MAXROW,
                INT *IBEST, INT *JBEST, INT *MBEST, TR AMAXR[]);
    void LU1MSP(INT MAXMN, TR LTOL, INT MAXCOL,
                INT *IBEST, INT *JBEST, INT *MBEST);
    void LU1MXC(INT K1, INT K2, INT IX[]);
    //void LU1MXR(INT K1, INT K2, INT IX[], TR AMAXR[]);
    void LU1MXR(INT MARK, INT K1, INT K2, INT IX[], INT COLS[], INT MARKC[], INT MARKR[], TR AMAXR[]);
    bool LU1FUL(INT LEND, INT LU1, INT PIV,
                INT MLEFT, INT NLEFT, INT NRANK, INT NROWU,
                INT *LENL, INT *LENU, INT *NSING,
                bool KEEPLU, TR SMALL, T D[], INT IPVT[]);
    void LU1OR1(TR SMALL, TR *AMAX, INT *NUMNZ, INT *LERR, INT *INFORM);
    void LU1OR2();
    void LU1OR3(INT *LERR, INT *INFORM);
    void LU1OR4();
    void LU1PEN(INT NSPARE, INT *ILAST,
                INT LPIVC1, INT LPIVC2, INT LPIVR1, INT LPIVR2,
                INT *LROW, INT IFILL[], INT JFILL[]);
    void LU1FAD(INT *INFORM, INT *LENL, INT *LENU, INT *MINLEN,
                INT *MERSUM, INT *NUTRI, INT *NLTRI,
                INT *NDENS1, INT *NDENS2, INT *NRANK,
                TR *LMAX, TR *UMAX, TR *DUMAX, TR *DUMIN, TR *AKMAX);
    void LU1FAC(INT *INFORM);

    //--------------------------------------------------------------------------
    //                      lusol6.cpp
    //--------------------------------------------------------------------------
    void LU6CHK(INT MODE, INT LENA2, INT *INFORM);
    void LU1L0T(LUSOLmat<T> **mat);
    void LU6L0T_v(LUSOLmat<T> *mat, T V[], INT NZidx[]);
    void LU6L(INT *INFORM, T V[], INT NZidx[]);
    void LU6LD(INT *INFORM, INT MODE, T V[], INT []);
    void LU6LT(INT *INFORM, T V[], INT NZidx[]);
    void LU6U(INT *INFORM, T V[], T W[], INT NZidx[]);
    void LU6UT(INT *INFORM, T V[], T W[], INT NZidx[]);
    void LU6SOL(INT MODE, T V[], T W[], INT NZidx[], INT *INFORM);

    //--------------------------------------------------------------------------
    //                      lusol7.cpp
    //--------------------------------------------------------------------------
    void LU7ADD(INT JADD, T V[], INT LENL, INT *LENU,
                INT *LROW, INT NRANK, INT *INFORM, INT *KLAST, TR *VNORM);
    void LU7CYC(INT KFIRST, INT KLAST, INT IX[]);
    void LU7ELM(INT JELM, T V[], INT *LENL,
                INT *LROW, INT NRANK, INT *INFORM, T *DIAG);
    void LU7FOR(INT KFIRST, INT KLAST, INT *LENL, INT *LENU,
                         INT *LROW, INT *INFORM, T *DIAG);
    void LU7RNK(INT JSING, INT *LENU,
                INT *LROW, INT *NRANK, INT *INFORM, T *DIAG);
    void LU7ZAP(INT JZAP, INT *KZAP, INT *LENU, INT *LROW, INT NRANK);

    //--------------------------------------------------------------------------
    //                      lusol8.cpp
    //--------------------------------------------------------------------------
    void LU8RPC(INT MODE1, INT MODE2, INT JREP, T V[], T W[],
                INT *INFORM, T *DIAG, TR *VNORM);
};

//--------------------------------------------------------------------------
//                      lusol.cpp
//--------------------------------------------------------------------------
template<class VR>
char relationChar(VR left, VR right)
{
    if(left > right)
        return '>';
    else if(left == right)
        return '=';
    else
        return '<';
};

void* clean_realloc(void *oldptr, INT width, INT newsize, INT oldsize);

};

