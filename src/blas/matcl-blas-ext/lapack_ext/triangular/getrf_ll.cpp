/*
 *  This file is a part of Matrix Computation Library (MATCL)
 *
 *  Copyright (c) Pawe³ Kowal 2017 - 2021
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
template<class T>
void getrf_ll(i_type M,i_type N, T* A,i_type LDA,i_type* IPIV,i_type* INFO)
{
// Purpose
// =======
//
// DGETRF computes an LU factorization of a general M-by-N matrix A
// using partial pivoting with row interchanges.
//
// The factorization has the form
//    A = P * L * U
// where P is a permutation matrix, L is lower triangular with unit
// diagonal elements (lower trapezoidal if m > n), and U is upper
// triangular (upper trapezoidal if m < n).
//
// This is the left-looking Level 3 BLAS version of the algorithm.
//
// Arguments
// =========
//
// M       (input) INTEGER
//         The number of rows of the matrix A.  M >= 0.
//
// N       (input) INTEGER
//         The number of columns of the matrix A.  N >= 0.
//
// A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
//         On entry, the M-by-N matrix to be factored.
//         On exit, the factors L and U from the factorization
//         A = P*L*U; the unit diagonal elements of L are not stored.
//
// LDA     (input) INTEGER
//         The leading dimension of the array A.  LDA >= max(1,M).
//
// IPIV    (output) INTEGER array, dimension (min(M,N))
//         The pivot indices; for 1 <= i <= min(M,N), row i of the
//         matrix was interchanged with row IPIV(i).
//
// INFO    (output) INTEGER
//         = 0:  successful exit
//         < 0:  if INFO = -i, the i-th argument had an illegal value
//         > 0:  if INFO = i, U(i,i) is exactly zero. The factorization
//               has been completed, but the factor U is exactly
//               singular, and division by zero will occur if it is used
//               to solve a system of equations.
//
// =====================================================================
    --A;
    --IPIV;

    T ONE = 1.;

    i_type I, IINFO, J, JB, K, NB;

    *INFO = 0;
    if ( M < 0 )
        *INFO = -1;
    else if ( N < 0 )
        *INFO = -2;
    else if( LDA < lapack::maximum( (i_type)1, M ) )
        *INFO = -4;

    if ( *INFO != 0 )
         return;

    // Quick return if possible
    if( M == 0 || N == 0 )
        return;

    // Determine the block size for this environment.
    NB = lapack::ilaenv( 1, "DGETRF", " ", M, N, -1, -1 );

    if( NB <= 1 || NB >= lapack::minimum( M, N ) )
    {
        // Use unblocked code.
        getf2( M, N, A+1, LDA, IPIV+1, INFO );
    }
    else
    {
        // Use blocked code.
        for(J = 1; J <= lapack::minimum( M, N ); J += NB)
        {
            JB = lapack::minimum( lapack::minimum( M, N )-J+1, NB );

            // Update before factoring the current panel
            for(K = 1; K <= J-NB; K += NB)
            {
                // Apply interchanges to rows K:K+NB-1.
                lapack::laswp( JB, A+1+(J-1)*LDA, LDA, K, K+NB-1, IPIV+1, 1 );
        
                // Compute block row of U.
                lapack::trsm( "Left", "Lower", "No transpose", "Unit",
                        NB, JB, ONE, A+K+(K-1)*LDA, LDA, A+K+(J-1)*LDA, LDA );

                // Update trailing submatrix.
                lapack::gemm( "No transpose", "No transpose", 
                      M-K-NB+1, JB, NB, -ONE, 
                      A+K+NB+(K-1)*LDA, LDA, A+K+(J-1)*LDA, LDA, ONE, A+K+NB+(J-1)*LDA, LDA );
            };

            // Factor diagonal and subdiagonal blocks and test for exact
            // singularity.

            getf2( M-J+1, JB, A+J+(J-1)*LDA, LDA, IPIV+J, &IINFO );

            // Adjust INFO and the pivot indices.
            if( *INFO == 0 && IINFO > 0 )
                *INFO = IINFO + J - 1;

            for(I = J; I <= lapack::minimum( M, J+JB-1 ); ++I)
                IPIV[I] = J - 1 + IPIV[I];
        };

        // Apply interchanges to the left-overs
        for(K = 1; K <= lapack::minimum( M, N ); K += NB)
        {
            lapack::laswp( K-1, A+1, LDA, K, 
                lapack::minimum(K+NB-1, lapack::minimum( M, N )), IPIV+1, 1 );
        };

        // Apply update to the M+1:N columns when N > M
        if ( N > M )
        {
            lapack::laswp( N-M, A+1+M*LDA, LDA, 1, M, IPIV+1, 1 );

            for(K = 1;K <= M; K += NB)
            {
                JB = lapack::minimum( M-K+1, NB );
                lapack::trsm( "Left", "Lower", "No transpose", "Unit",
                      JB, N-M, ONE, A+K+(K-1)*LDA, LDA, A+K+M*LDA, LDA );

                if ( K+NB <= M )
                {
                    lapack::gemm( "No transpose", "No transpose",
                           M-K-NB+1, N-M, NB, -ONE, 
                           A+K+NB+(K-1)*LDA, LDA, A+K+M*LDA, LDA, ONE,
                           A+K+NB+M*LDA, LDA );
                };
            };
        };
    };
    return;
};

template void BLAS_EXT_EXPORT
getrf_ll<d_type>(i_type M, i_type N, d_type *A, i_type LDA, i_type* IPIV, i_type *INFO );
template void BLAS_EXT_EXPORT
getrf_ll<z_type>(i_type M, i_type N, z_type *A, i_type LDA, i_type* IPIV, i_type *INFO );

};};