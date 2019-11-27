/*
 *  This file is a part of Matrix Computation Library (MATCL)
 *
 *  Copyright (c) Pawe³ Kowal 2017 - 2018
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
 */

#include "matcl-blas-lapack/blas/blas.h"
#include "matcl-blas-lapack/lapack/lapack.h"
#include "matcl-blas-ext/lapack_ext/lapack_ext.h"

namespace matcl { namespace lapack
{

template<class V> BLAS_EXT_EXPORT
typename details::enable_if_valid_real<void,V>::type
qtrtrs(const char* TRANA, i_type M, i_type N, const V* A, i_type LDA, 
            V* C, i_type LDC, i_type* INFO)
{
    using VR    = typename details::real_type<V>::type;

    // Decode and Test input parameters
    bool NOTRNA     = (TRANA[0] == 'N' || TRANA[0] == 'n');
    *INFO           = 0;

    if (NOTRNA == false && TRANA[0] != 'T' && TRANA[0] != 't' && TRANA[0] != 'C' && TRANA[0] != 'c')
        *INFO       = -1;
    else if ( M < 0 )
        *INFO       = -2;
    else if ( N < 0 )
        *INFO       = -3;
    else if (LDA < lapack::maximum(1,M) )
        *INFO       = -5;
    else if (LDC < lapack::maximum(1, M) )
        *INFO       = -7;

    if (*INFO != 0)
        return;

    //constants 
    VR SMLNUM       = lapack::lamch<VR>( "S" );
    VR EPS          = lapack::lamch<VR>( "P" );
    VR norm_A       = lapack::lange( "M", M, M, A, LDA, (typename details::real_type<V>::type*)nullptr);
    VR SMIN         = lapack::maximum( SMLNUM, EPS*norm_A);

    // Quick return if possible

    if ( M == 0 || N == 0 )
        return;

    if( NOTRNA == true)
    {
        //  Solve    A*X = C.
        //
        //  The (K,L)th block of X is determined starting from
        //  bottom-left corner column by column by
        //
        //  A(K,K)*X(K,L) = C(K,L) - R(K,L)
        //
        //  Where
        //                 M                 
        //       R(K,L) = SUM [A(K,I)*X(I,L)]
        //               I=K+1               
        //

        for(i_type L = 1; L <= N; ++L)
        {
            // Start row loop (index = K)
            // K1 (K2): row index of the first (last) row of X(K,L).

            i_type KNEXT    = M;
            for (i_type K = M; K >= 1; --K)
            {
                if ( K > KNEXT)
                    continue;

                i_type K1, K2;

                if ( K == 1 )
                {
                    K1          = K;
                    K2          = K;
                }
                else
                {
                    if( A[K-1 + (K-2)*LDA] != 0. ) 
                    {
                        K1      = K - 1;
                        K2      = K;
                        KNEXT   = K - 2;
                    }
                    else
                    {
                        K1      = K;
                        K2      = K;
                        KNEXT   = K - 1;
                    };
                };

                if ( K1 == K2 )
                {
                    V A11       = A[K1-1 + (K1-1)*LDA];

                    if( A11 == V(0))
                    {
                        *INFO   = K1;
                        return;
                    };

                    i_type pos  = lapack::minimum( K1 + 1, M);
                    V SUML      = lapack::dot( M-K1, A + K1-1 + (pos - 1) * LDA, LDA, 
                                        C + pos - 1 + (L-1) * LDC, 1 );
                    V VEC       = C[K1-1 + (L-1)*LDC] - SUML;

                    C[K1-1 + (L-1)*LDC]
                                = VEC / A11;
                }
                else if( K1 != K2 )
                {
                    i_type pos  = lapack::minimum( K2+1, M );
                    V SUML1     = lapack::dot(M-K2, A + K1-1 + (pos-1)*LDA, LDA,
                                                C + pos - 1 + (L-1)*LDC, 1 );
                    V SUML2     = lapack::dot(M-K2, A + K2-1 + (pos-1)*LDA, LDA,
                                                C + pos - 1 + (L-1)*LDC, 1 );
                    V VEC[2];
                    V X[2];

                    VEC[0]      = C[K1-1 + (L-1)*LDC] - SUML1;
                    VEC[1]      = C[K2-1 + (L-1)*LDC] - SUML2;

                    V SCALOC    = 1.;
                    V XNORM     = 0.;
                    i_type IERR = 0;

                    lapack::laln2<V>(false, 2, 1, SMIN, 1., A + K1-1 + (K1-1)*LDA, LDA, 
                                    1., 1., VEC, 2, 0.0, 0.0, X, 2, &SCALOC, &XNORM, &IERR);

                    if ( IERR != 0 || SCALOC != 1.)
                    {
                        *INFO   = K1;
                        return;
                    };

                    C[K1-1 + (L-1)*LDC]     = X[0];
                    C[K2-1 + (L-1)*LDC]     = X[1];
                };
            };
        };
    }
    else if( NOTRNA == false)
    {
        //
        //  Solve    A' *X = C.
        //
        //  The (K,L)th block of X is determined starting from
        //  upper-left corner column by column by
        //
        //  A(K,K)'*X(K,L) = C(K,L) - R(K,L)
        //
        //  Where
        //                  K-1                 
        //         R(K,L) = SUM [A(I,K)'*X(I,L)]
        //                  I=1

        for (i_type L = 1; L <= N; ++L)
        {
            // Start row loop (index = K)
            // K1 (K2): row index of the first (last) row of X(K,L)

            i_type KNEXT    = 1;
            for (i_type K = 1; K <= M; ++K)
            {
                if ( K < KNEXT )
                    continue;

                i_type K1, K2;

                if ( K == M )
                {
                    K1          = K;
                    K2          = K;
                }
                else
                {
                    if( A[K + (K-1)*LDA] != 0 )
                    {
                        K1      = K;
                        K2      = K + 1;
                        KNEXT   = K + 2;
                    }
                    else
                    {
                        K1      = K;
                        K2      = K;
                        KNEXT   = K + 1;
                    };
                };

                if ( K1 == K2 )
                {
                    V A11       = A[K1-1 + (K1-1)*LDA];
                    if( A11 == 0)
                    {
                        *INFO   = K1;
                        return;
                    };

                    V SUML      = lapack::dot( K1-1, A + (K1-1)*LDA, 1, C+ (L-1)*LDA, 1 );
                    V VEC       = C[K1-1 + (L-1)*LDC] - SUML;
                    
                    C[K1-1 + (L-1)*LDC]
                                = VEC / A11;
                }
                else
                {
                    V SUML1     = lapack::dot(K1-1, A + (K1-1)*LDA, 1, C + (L-1)*LDC, 1 );
                    V SUML2     = lapack::dot(K1-1, A + (K2-1)*LDA, 1, C + (L-1)*LDC, 1 );

                    V VEC[2];
                    V X[2];
                    V SCALOC    = 0.;
                    V XNORM     = 0.;
                    i_type IERR = 0;

                    VEC[0]      = C[K1-1 + (L-1)*LDC] - SUML1;
                    VEC[1]      = C[K2-1 + (L-1)*LDC] - SUML2;

                    lapack::laln2<V>(true, 2, 1, SMIN, 1., A + (K1-1) + (K1-1)*LDA, LDA, 1., 1., 
                                  VEC, 2, 0., 0., X, 2, &SCALOC, &XNORM, &IERR );

                    if (IERR != 0 || SCALOC != 1.)
                    {
                        *INFO   = K1;
                        return;
                    };

                    C[K1 - 1 + (L-1)*LDC]   = X[0];
                    C[K2 - 1 + (L-1)*LDC]   = X[1];
                };
            };
        };
    };

    return;
};

template BLAS_EXT_EXPORT
void qtrtrs<s_type>(const char* TRANA, i_type M, i_type N, const s_type* A, i_type LDA, 
                    s_type* C, i_type LDC, i_type* INFO);

template BLAS_EXT_EXPORT
void qtrtrs<d_type>(const char* TRANA, i_type M, i_type N, const d_type* A, i_type LDA, 
                    d_type* C, i_type LDC, i_type* INFO);

};};