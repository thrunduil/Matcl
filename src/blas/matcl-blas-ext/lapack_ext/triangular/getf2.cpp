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
void getf2(i_type M, i_type N, T *A, i_type LDA, i_type* IPIV, i_type *INFO )
{
//  -- LAPACK routine (version 3.1) --
//     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
//     November 2006
//
//    .. Scalar Arguments ..
//    INTEGER            INFO, LDA, M, N
//    ..
//    .. Array Arguments ..
//    INTEGER            IPIV( * )
//    DOUBLE PRECISION   A( LDA, * )
//    ..
//
// Purpose
// =======
//
// DGETF2 computes an LU factorization of a general m-by-n matrix A
// using partial pivoting with row interchanges.
//
// The factorization has the form
//    A = P * L * U
// where P is a permutation matrix, L is lower triangular with unit
// diagonal elements (lower trapezoidal if m > n), and U is upper
// triangular (upper trapezoidal if m < n).
//
// This is the right-looking Level 2 BLAS version of the algorithm.
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
//         On entry, the m by n matrix to be factored.
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
//         = 0: successful exit
//         < 0: if INFO = -k, the k-th argument had an illegal value
//         > 0: if INFO = k, U(k,k) is exactly zero. The factorization
//              has been completed, but the factor U is exactly
//              singular, and division by zero will occur if it is used
//              to solve a system of equations.
//
// =====================================================================

    T  ONE     = 1.;
    T  ZERO    = 0.;
      
    i_type     J, JP;

    --A;
    --IPIV;

    //  Test the input parameters.
    *INFO = 0;
    if ( M < 0 )
        *INFO = -1;
    else if ( N < 0 )
        * INFO = -2;
    else if( LDA < lapack::maximum((i_type)1,M) )
        *INFO = -4;

    if( *INFO != 0 )
        return;

    // Quick return if possible
    if( M == 0 || N == 0 )
        return;

    for(J = 1; J <= lapack::minimum( M, N ); ++J)
    {
        // Find pivot and test for singularity.
        JP = J - 1 + lapack::amax( M-J+1, A+J+(J-1)*LDA, 1 );
        IPIV[J] = JP;

        if( A[JP+(J-1)*LDA] != ZERO )
        {
            // Apply the interchange to columns 1:N.
            if( JP != J )
                lapack::swap( N, A+J, LDA, A+JP, LDA );

            // Compute elements J+1:M of J-th column.
            if( J < M )
            {
                lapack::scal( M-J, ONE / A[J+(J-1)*LDA], A+J+1+(J-1)*LDA, 1 );
            };
        }
        else if( INFO == 0 )
        {
            *INFO = J;
        };

        if( J < lapack::minimum( M, N ) )
        {        
            // Update trailing submatrix.
            lapack::geru( M-J, N-J, -ONE, A+J+1+(J-1)*LDA, 1, A+J+J*LDA, LDA, A+J+1+J*LDA, LDA );
        };
    };
    return;
};

template void BLAS_EXT_EXPORT
getf2<d_type>(i_type M, i_type N, d_type *A, i_type LDA, i_type* IPIV, i_type *INFO );
template void BLAS_EXT_EXPORT
getf2<z_type>(i_type M, i_type N, z_type *A, i_type LDA, i_type* IPIV, i_type *INFO );

};};