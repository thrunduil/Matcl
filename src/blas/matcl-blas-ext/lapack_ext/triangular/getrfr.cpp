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

#include <cstdlib>
#include "matcl-blas-lapack/blas/blas.h"
#include "matcl-blas-lapack/lapack/lapack.h"
#include "matcl-blas-ext/lapack_ext/lapack_ext.h"
#include "blas/matcl-blas-ext/lapack_ext/utils/optim_params.h"

namespace matcl { namespace lapack
{

template<class T>
void lapack::getrfr(i_type M, i_type N, T *A, i_type LDA, i_type *IPIV, i_type *IQIV, 
            typename lapack::details::real_type<T>::type TOLC,
            typename lapack::details::real_type<T>::type TOLR, 
            typename lapack::details::real_type<T>::type TOL,
            T* WORK, i_type LWORK, i_type *INFO)
{
    --A;
    --IPIV;
    --IQIV;

    using TR = typename lapack::details::real_type<T>::type;

    i_type NB, JB, IINFO, I;
    i_type J    = 1;

    T  ONE = 1.;
    
    lapack::getf2r(M, &N, N, 1, A+1, LDA, IPIV + 1, IQIV + 1, TOLC, TOLR, TOL, WORK, -1, INFO );

    if (*INFO < 0)
        return;

    if (LWORK < 0)
        return;

    i_type LWRK = (i_type)real(WORK[0]);

    if (LWORK < LWRK)
    {
        *INFO = -11;
        return;
    };

    *INFO = 0;

    if( M == 0 || N == 0 )
        return;

    // Determine the block size for this environment.
    NB = optim_params::getrfr_block;

    i_type MN   = lapack::minimum( M, N );

    TR TOLV     = TOL;
    TR NORM_F   = TR(1.0);

    if (TOL < TR(0.0))
    {
        NORM_F  = lapack::lange<T>("F", M, N, A+1, LDA, nullptr);
        NORM_F  = NORM_F / sqrt(TR(MN));

        TR EPS  = lapack::lamch<TR>("PRECISION");
        TOL     = TR(10.0) * EPS * NORM_F;
    }
    
    if (NORM_F == TR(0.0))
    {
        N       = 0;
        goto lab_exit;
    };

    if( NB <= 1 || NB >= lapack::minimum( M, N ) )
    {
        // Use unblocked code.
         lapack::getf2r(M, &N, N, 1, A+1, LDA, IPIV + 1, IQIV + 1, TOLC, TOLR, TOLV, WORK, LWORK, INFO );
         INFO[0] = lapack::minimum( M, N );
         return;
    }    

    // Use blocked code.
    for(; J <= lapack::minimum( M, N ); J += NB)
    {
        JB = lapack::minimum( lapack::minimum( M, N )-J+1, NB );

        // Factor diagonal, subdiagonal and block row of U
        IINFO = 0;
        i_type Nptr = N - J + 1;
        lapack::getf2r(M - J + 1, &Nptr, JB, J, A+J+(J-1)*LDA, LDA, IPIV + J, IQIV + J, 
                    TOLC, TOLR, TOLV, WORK, LWORK, &IINFO );

        N = N + Nptr - (N - J + 1);

        // Adjust pivot indices.

        for(I = J; I <= lapack::minimum( M, J+JB-1 ); ++I)
            IPIV[I] = J - 1 + IPIV[I];

        // Apply interchanges to columns 1:J-1.
        lapack::laswp( J-1, A+1, LDA, J, J+JB-1, IPIV+1, 1 );

        if( J+JB <= N && J+JB <= M)
        {
            // Update trailing submatrix.
            lapack::gemm( "No transpose", "No transpose", M-J-JB+1,
                      N-J-JB+1, JB, -ONE, A+J+JB+(J-1)*LDA, LDA,
                      A+J+(J+JB-1)*LDA, LDA, ONE, A+J+JB+(J+JB-1)*LDA, LDA );
        };
    };

  lab_exit:
    MN          = lapack::minimum( M, N );

    for (I = J; I <= MN; ++I)
        IPIV[I] = I;

    INFO[0] = lapack::minimum( M, N );

    return;
};

template void BLAS_EXT_EXPORT
getrfr<s_type>(i_type M, i_type N, s_type *A, i_type LDA, i_type* IPIV, i_type* IQIV, 
                       s_type TOLC, s_type TOLR, s_type TOLV, s_type* WORK, i_type LWORK, i_type *INFO );

template void BLAS_EXT_EXPORT
getrfr<d_type>(i_type M, i_type N, d_type *A, i_type LDA, i_type* IPIV, i_type* IQIV, 
                       d_type TOLC, d_type TOLR, d_type TOLV, d_type* WORK, i_type LWORK, i_type *INFO );

template void BLAS_EXT_EXPORT
getrfr<z_type>(i_type M, i_type N, z_type *A, i_type LDA, i_type* IPIV, i_type* IQIV, 
                       d_type TOLC, d_type TOLR, d_type TOLV, z_type* WORK, i_type LWORK, i_type *INFO );

template void BLAS_EXT_EXPORT
getrfr<c_type>(i_type M, i_type N, c_type *A, i_type LDA, i_type* IPIV, i_type* IQIV, 
                       s_type TOLC, s_type TOLR, s_type TOLV, c_type* WORK, i_type LWORK, i_type *INFO );

};};