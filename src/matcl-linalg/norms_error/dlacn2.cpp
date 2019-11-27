#include "matcl-linalg/norms_error/norm.h"
#include "matcl-matrep/details/utils.h"
#include "matcl-matrep/lib_functions/func_unary.h"
#include "matcl-matrep/lib_functions/func_binary.h"

namespace mmlib
{

namespace md = mmlib::details;

//--------------------------------------------------------------------------------
//                      DLACN2
//--------------------------------------------------------------------------------
/*
*  Purpose
*  =======
*
*  DLACN2 estimates the 1-norm of a square, real matrix A.
*  Reverse communication is used for evaluating matrix-vector products.
*
*  Arguments
*  =========
*
*  N      (input) INTEGER
*         The order of the matrix.  N >= 1.
*
*  V      (workspace) DOUBLE PRECISION array, dimension (N)
*         On the final return, V = A*W,  where  EST = norm(V)/norm(W)
*         (W is not returned).
*
*  X      (input/output) DOUBLE PRECISION array, dimension (N)
*         On an intermediate return, X should be overwritten by
*               A * X,   if KASE=1,
*               A' * X,  if KASE=2,
*         and DLACN2 must be re-called with all the other parameters
*         unchanged.
*
*  ISGN   (workspace) INTEGER array, dimension (N)
*
*  EST    (input/output) DOUBLE PRECISION
*         On entry with KASE = 1 or 2 and ISAVE(1) = 3, EST should be
*         unchanged from the previous call to DLACN2.
*         On exit, EST is an estimate (a lower bound) for norm(A). 
*
*  KASE   (input/output) INTEGER
*         On the initial call to DLACN2, KASE should be 0.
*         On an intermediate return, KASE will be 1 or 2, indicating
*         whether X should be overwritten by A * X  or A' * X.
*         On the final return from DLACN2, KASE will again be 0.
*
*  ISAVE  (input/output) INTEGER array, dimension (3)
*         ISAVE is used to save variables between calls to DLACN2
*
*  Further Details
*  ======= =======
*
*  Contributed by Nick Higham, University of Manchester.
*  Originally named SONEST, dated March 16, 1988.
*
*  Reference: N.J. Higham, "FORTRAN codes for estimating the one-norm of
*  a real or complex matrix, with applications to condition estimation",
*  ACM Trans. Math. Soft., vol. 14, no. 4, pp. 381-396, December 1988.
*
*/

template<class Val>
typename md::real_type<Val>::type abs_sum(const Val*X, Integer N)
{
    using VR = typename md::real_type<Val>::type;
    if (N == 0)
        return VR(0);

    VR ret = abs(X[0]);
            
    for (Integer i = 1; i < N; ++i)
        ret += abs(X[i]);

    return ret;
};

template<class Val>
Integer idamax(const Val*X, Integer N)
{
    using VR = typename md::real_type<Val>::type;

    if (N == 0)
        return 0;

    Integer pos     = 1;
    VR max_abs      = abs(X[0]);

    for (Integer i = 1; i < N; ++i)
    {
        VR tmp      = abs(X[i]);
        if (tmp > max_abs)
        {
            max_abs = tmp;
            pos     = i + 1;
        };
    };

    return pos;
};

template<class Val>
void mmlib::lacn2_ex(Integer N, Val* V, Val* X, Integer* ISGN, typename md::real_type<Val>::type& EST, Integer& KASE,
           Integer* ISAVE )
{
    using VR = typename md::real_type<Val>::type;

    static const Integer ITMAX  = 15;

    if( KASE == 0 )
    {
        Val tmp         = Val(1.0) / VR( N );

        for( Integer I = 0; I < N; ++I)
            X[I]        = tmp;

        KASE            = 1;
        ISAVE[1-1]      = 1;
        return;
    };

    switch (ISAVE[1-1])
    {
        case 1:
        {
            //    FIRST ITERATION.  X HAS BEEN OVERWRITTEN BY A*X.

            if ( N == 1 )
            {
                V[1-1]      = X[1-1];
                EST         = abs( V[1-1] );
    
                // QUIT
                KASE        = 0;
                return;
            }

            EST             = abs_sum(X,N);
        
            for (Integer I = 1; I <= N; ++I)
            {
                X[I-1]      = mmlib::sign(X[I-1]);
                ISGN[I-1]   = mmlib::iround(X[I-1]);//?
            };
              
            KASE        = 2;
            ISAVE[1-1]  = 2;

            return;
        }
        case 2:
        {
            //     FIRST ITERATION.  X HAS BEEN OVERWRITTEN BY TRANSPOSE(A)*X.
            Integer pos     = idamax(X,N);

            ISAVE[2-1]      = pos;
            ISAVE[3-1]      = 2;

            //     MAIN LOOP - ITERATIONS 2,3,...,ITMAX.
          lab_50:
            for (Integer i = 0; i < N; ++i)
            {
                X[i]        = Val(0.0);
            };

            X[ISAVE[2-1]-1] = Val(1.0);
            
            KASE            = 1;
            ISAVE[1-1]      = 3;
          
            return;
        }
        case 3:
        {
            //     X HAS BEEN OVERWRITTEN BY A*X.
            for (Integer i = 0; i < N; ++i)
                V[i]        = X[i];

            VR ESTOLD       = EST;
            EST             = abs_sum(V, N);

            bool conv       = true;
            for (Integer I = 1; I <= N; ++I)
            {
                Integer sgn = isign(X[I-1]);
                
                if (sgn != ISGN[I-1])
                {
                    conv    = false;
                    break;
                };
            };

            if (conv == true)
            {
                // REPEATED SIGN VECTOR DETECTED, HENCE ALGORITHM HAS CONVERGED.
                goto lab_120;
            };

            //TEST FOR CYCLING.
            if (EST <= ESTOLD )
                goto lab_120;

            for (Integer I = 1; I <= N; ++I)
            {
                X[I-1]      = sign(X[I-1]);
                ISGN[I-1]   = isign(X[I-1]);
            };

            KASE            = 2;
            ISAVE[1-1]      = 4;
            return;
        }
        case 4:
        {
            // X HAS BEEN OVERWRITTEN BY TRANSPOSE(A)*X.
            Integer JLAST   = ISAVE[2-1];
            ISAVE[2-1]      = idamax(X, N);

            if( ( X[JLAST-1] != abs( X[ISAVE[2-1]-1] ) ) && ( ISAVE[3-1] < ITMAX ) )
            {
                ISAVE[3-1]  = ISAVE[3-1] + 1;
                goto lab_50;
            };

            goto lab_120;
        }
        case 5:
        {
            // X HAS BEEN OVERWRITTEN BY A*X.
            VR TEMP = VR(2.0) * abs_sum(X, N) / VR ( 3*N );

            if ( TEMP > EST )
            {
                for (Integer i = 0; i < N; ++i)
                    V[i] = X[i];

                EST = TEMP;
            };

            KASE    = 0;
            return;
        }
    };

    lab_120:

    // ITERATION COMPLETE.  FINAL STAGE.
    VR ALTSGN   = VR(1.0);    

    for (Integer I = 1; I <= N; ++I)
    {
        VR val  = VR(1.0) + VR(I-1)/VR(N-1);
        X[I-1]  = ALTSGN * val;
        ALTSGN  = -ALTSGN;
    };

    KASE        = 1;
    ISAVE[1-1]  = 5;
    return;
};

template void lacn2_ex<Real>(Integer N, Real* V, Real* X, Integer* ISGN, Real& EST, Integer& KASE, Integer* ISAVE );
template void lacn2_ex<Float>(Integer N, Float* V, Float* X, Integer* ISGN, Float& EST, Integer& KASE, Integer* ISAVE );
template void lacn2_ex<Complex>(Integer N, Complex* V, Complex* X, Integer* ISGN, Real& EST, Integer& KASE, 
                             Integer* ISAVE );
template void lacn2_ex<Float_complex>(Integer N, Float_complex* V, Float_complex* X, Integer* ISGN, Float& EST, 
                                      Integer& KASE, Integer* ISAVE );

};
