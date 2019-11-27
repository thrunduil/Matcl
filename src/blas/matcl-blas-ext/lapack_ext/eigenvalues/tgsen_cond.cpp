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
#include <algorithm>

namespace matcl { namespace lapack
{

template<class V> BLAS_EXT_EXPORT
typename details::enable_if_valid<void,V>::type
lapack::tgsen_cond(i_type IJOB, i_type N, i_type M, const V *A, i_type LDA, const V *B, i_type LDB, 
        typename details::real_type<V>::type& PL, typename details::real_type<V>::type& PR,
        typename details::real_type<V>::type* DIF, V *WORK, i_type LWORK, i_type *IWORK, i_type LIWORK,
        i_type& INFO)
{
    using VR = typename details::real_type<V>::type;

    static const i_type IDIFJB = 3;

    //     Decode and test the input parameters
    INFO            = 0;
    i_type LQUERY   = ( LWORK == -1 || LIWORK == -1 );

    if ( IJOB < 1 || IJOB > 5 )
        INFO    = -1;
    else if ( N < 0 )
        INFO    = -2;
    else if ( M < 0 || M > N)
        INFO    = -3;
    else if ( LDA < std::max( 1, N ) )
        INFO    = -5;
    else if ( LDB < std::max( 1, N ) )
        INFO    = -7;

    if (INFO != 0 )
        return;

    bool WANTP  = IJOB == 1 || IJOB >= 4;
    bool WANTD1 = IJOB == 2 || IJOB == 4;
    bool WANTD2 = IJOB == 3 || IJOB == 5;
    bool WANTD  = WANTD1 || WANTD2;

    i_type LWMIN;
    i_type LIWMIN;

    if ( IJOB == 1 || IJOB == 2 || IJOB == 4 )
    {
        LWMIN   = std::max( 4*N+16, 2*M*( N-M ) );
        LIWMIN  = std::max( 1, N+6 );
    }
    else
    {
        LWMIN   = std::max( 4*N+16, 4*M*( N-M ) );
        LIWMIN  = std::max( 2*M*( N-M ), N+6 );
    };

    WORK[0]     = V(VR(LWMIN));
    IWORK[0]    = LIWMIN;

    if ( LWORK < LWMIN && !LQUERY )
        INFO    = -12;
    else if ( LIWORK < LIWMIN && !LQUERY )
        INFO    = -14;

    if (INFO != 0 )
        return;
    else if (LQUERY )
        return;

    // Quick return if possible.
    if ( M == N || M == 0 )
    {
        if( WANTP )
        {
            PL  = VR(1.0);
            PR  = VR(1.0);
        };

        if (WANTD)
        {
            VR DSCALE   = VR(0.0);
            VR DSUM     = VR(1.0);

            for (i_type I = 1; I <= N; ++I)
            {
                lapack::lassq( N, A + (I-1) * LDA, 1, DSCALE, DSUM );
                lapack::lassq( N, B + (I-1) * LDB, 1, DSCALE, DSUM );
            };

            DIF[0]  = DSCALE * std::sqrt(DSUM);
            DIF[1]  = DIF[0];
        };

        goto lab_60;
    };

    if (WANTP)
    {
        // Solve generalized Sylvester equation for R and L
        // and compute PL and PR.
        i_type N1   = M;
        i_type N2   = N - M;
        i_type I    = N1 + 1;
        i_type IJB  = 0;

        VR DSCALE;
        i_type IERR;

        lapack::lacpy("Full", N1, N2, A + (I-1)*LDA, LDA, WORK, N1 );
        lapack::lacpy("Full", N1, N2, B + (I-1)*LDB, LDB, WORK + N1*N2, N1 );
        lapack::tgsyl<V>("N", IJB, N1, N2, A, LDA, A + (I-1) + (I-1)*LDA, LDA, WORK,
                        N1, B, LDB, B + (I-1) + (I-1)*LDB, LDB, WORK + N1*N2, N1,
                        &DSCALE, &DIF[0], WORK + N1*N2*2, LWORK-2*N1*N2, IWORK, &IERR );

        // Estimate the reciprocal of norms of "projections" onto left
        // and right eigenspaces.

        VR RDSCAL   = VR(0.0);
        VR DSUM     = VR(1.0);

        lapack::lassq(N1*N2, WORK, 1, RDSCAL, DSUM );

        PL          = RDSCAL * sqrt( DSUM );

        if (PL == VR(0.0))
            PL      = VR(1.0);
        else
            PL      = DSCALE / ( sqrt( DSCALE * DSCALE / PL+PL ) * sqrt( PL ) );

        RDSCAL      = VR(0.0);
        DSUM        = VR(1.0);

        lapack::lassq(N1*N2, WORK + N1*N2, 1, RDSCAL, DSUM );
         
        PR          = RDSCAL * sqrt(DSUM);

        if ( PR == VR(0.0) )
            PR      = VR(1.0);
        else
            PR      = DSCALE / ( sqrt( DSCALE*DSCALE / PR+PR ) * sqrt( PR ) );
    };

    if (WANTD)
    {
        // Compute estimates of Difu and Difl.

        if (WANTD1)
        {
            i_type N1   = M;
            i_type N2   = N - M;
            i_type I    = N1 + 1;
            i_type IJB  = IDIFJB;
            i_type IERR;

            VR DSCALE;

            // Frobenius norm-based Difu-estimate.
            lapack::tgsyl<V>("N", IJB, N1, N2, A, LDA, A+(I-1) + (I-1)*LDA, LDA, WORK,
                    N1, B, LDB, B +(I-1) +(I-1)*LDB, LDB, WORK + N1*N2, N1, &DSCALE, &DIF[0], 
                    WORK + 2*N1*N2, LWORK-2*N1*N2, IWORK, &IERR );

            // Frobenius norm-based Difl-estimate.
            lapack::tgsyl<V>("N", IJB, N2, N1, A +(I-1) + (I-1), LDA, A, LDA, WORK,
                    N2, B + (I-1) + (I-1)*LDB, LDB, B, LDB, WORK + N1*N2, N2, &DSCALE, &DIF[1], 
                    WORK + 2*N1*N2, LWORK-2*N1*N2, IWORK, &IERR );
        }
        else
        {
            // Compute 1-norm-based estimates of Difu and Difl using
            // reversed communication with DLACN2. In each step a
            // generalized Sylvester equation or a transposed variant
            // is solved.
            
            i_type KASE     = 0;
            i_type N1       = M;
            i_type N2       = N - M;
            i_type I        = N1 + 1;
            i_type IJB      = 0;
            i_type MN2      = 2*N1*N2;

            i_type ISAVE[3];
            VR DSCALE;
            i_type IERR;

            // 1-norm-based estimate of Difu.
          lab_40:

            lapack::lacn2(MN2, WORK + MN2, WORK, IWORK, &DIF[0], &KASE, ISAVE );

            if (KASE != 0 )
            {
                if (KASE == 1 )
                {
                    // Solve generalized Sylvester equation.
                    lapack::tgsyl<V>("N", IJB, N1, N2, A, LDA, A + (I-1) + (I-1)*LDA, LDA, WORK, N1,
                            B, LDB, B+ (I-1) + (I-1)*LDB, LDB, WORK + N1*N2, N1, &DSCALE, &DIF[0],
                            WORK + 2*N1*N2, LWORK-2*N1*N2, IWORK, &IERR );
                }
                else
                {

                    // Solve the transposed variant.
                    lapack::tgsyl<V>( "T", IJB, N1, N2, A, LDA, A + (I-1) + (I-1)*LDA, LDA, WORK, N1, 
                            B, LDB, B+ (I-1) + (I-1)*LDB, LDB, WORK + N1*N2, N1, &DSCALE, &DIF[0],
                            WORK + 2*N1*N2, LWORK-2*N1*N2, IWORK, &IERR );
                };

                goto lab_40;
            };

            DIF[0]  = DSCALE / DIF[0];

            // 1-norm-based estimate of Difl.
          lab_50:

            lapack::lacn2(MN2, WORK + MN2, WORK, IWORK, &DIF[1], &KASE, ISAVE );

            if (KASE != 0 )
            {
                if (KASE == 1 )
                {
                    // Solve generalized Sylvester equation.
                    lapack::tgsyl("N", IJB, N2, N1, A + (I-1) + (I-1)*LDA, LDA, A, LDA, WORK, N2, 
                        B + (I-1) + (I-1)*LDB, LDB, B, LDB, WORK + N1*N2, N2, &DSCALE, &DIF[1],
                        WORK + 2*N1*N2, LWORK-2*N1*N2, IWORK, &IERR );
                }
                else
                {
                    // Solve the transposed variant.
                    lapack::tgsyl("T", IJB, N2, N1, A + (I-1) + (I-1)*LDA, LDA, A, LDA, WORK, N2, 
                        B + (I-1) + (I-1)*LDB, LDB, B, LDB, WORK + N1*N2, N2, &DSCALE, &DIF[1],
                        WORK + 2*N1*N2, LWORK-2*N1*N2, IWORK, &IERR );
                };

                goto lab_50;
            };

            DIF[1] = DSCALE / DIF[1];
        };        
    };

  lab_60:

    WORK[0]     = V(VR(LWMIN));
    IWORK[0]    = LIWMIN;
};

template BLAS_EXT_EXPORT
void tgsen_cond<d_type>(i_type IJOB, i_type n, i_type m, const d_type *a,i_type lda, const d_type *b, i_type ldb, 
        d_type& pl, d_type& pr, d_type* dif, d_type *work, i_type LWORK, i_type *iwork, i_type LIWORK,
        i_type& INFO);

template BLAS_EXT_EXPORT
void tgsen_cond<s_type>(i_type IJOB, i_type n, i_type m, const s_type *a,i_type lda, const s_type *b, i_type ldb, 
        s_type& pl, s_type& pr, s_type* dif, s_type *work, i_type LWORK, i_type *iwork, i_type LIWORK,
        i_type& INFO);

template BLAS_EXT_EXPORT
void tgsen_cond<c_type>(i_type IJOB, i_type n, i_type m, const c_type *a,i_type lda, const c_type *b, i_type ldb, 
        s_type& pl, s_type& pr, s_type* dif, c_type *work, i_type LWORK, i_type *iwork, i_type LIWORK,
        i_type& INFO);

template BLAS_EXT_EXPORT
void tgsen_cond<z_type>(i_type IJOB, i_type n, i_type m, const z_type *a,i_type lda, const z_type *b, i_type ldb, 
        d_type& pl, d_type& pr, d_type* dif, z_type *work, i_type LWORK, i_type *iwork, i_type LIWORK,
        i_type& INFO);

}}