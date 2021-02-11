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
void getrf_rec(i_type M, i_type N, T* A,i_type LDA,i_type* IPIV,i_type* INFO )
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
// This code implements an iterative version of Sivan Toledo's recursive
// LU algorithm[1].  For square matrices, this iterative versions should
// be within a factor of two of the optimum number of memory transfers.
//
// [1] Toledo, S. 1997. Locality of Reference in LU Decomposition with
// Partial Pivoting. SIAM J. Matrix Anal. Appl. 18, 4 (Oct. 1997),
// 1065-1081. http://dx.doi.org/10.1137/S0895479896297744
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

    T ONE = 1.;;
    T ZERO = 0.;
    T NEGONE = -1.;

    T  TMP;
    i_type  J, JP, NSTEP, NTOPIV, NPIVED, KAHEAD;
    i_type  KSTART, IPIVSTART, JPIVSTART, KCOLS = 0;

    // Test the input parameters.
    *INFO = 0;

    if( M < 0 ) 
        *INFO = -1;
    else if( N < 0 ) 
        *INFO = -2;
    else if( LDA < lapack::maximum( (i_type)1, M ) ) 
        *INFO = -4;

    if( *INFO != 0 ) 
        return;

    // Quick return if possible
    if( M == 0 || N == 0 )
        return;
        
    NSTEP = lapack::minimum( M, N );
    for(J = 1; J <= NSTEP; ++J)
    {
        KAHEAD = lapack::iand( J, -J );
        KSTART = J + 1 - KAHEAD;
        KCOLS = lapack::minimum( KAHEAD, NSTEP-J );

        // Find pivot.
        JP = J - 1 + lapack::amax( M-J+1, A+J+(J-1)*LDA, 1 );
        IPIV[J] = JP;

        // Permute just this column.
        if (JP  !=  J) 
        {
            TMP = A[J+(J-1)*LDA];
            A[J+(J-1)*LDA] = A[JP+(J-1)*LDA];
            A[JP+(J-1)*LDA] = TMP;
        };

        // Apply pending permutations to L
        NTOPIV = 1;
        IPIVSTART = J;
        JPIVSTART = J - NTOPIV;
        while ( NTOPIV  <  KAHEAD )
        {
            lapack::laswp( NTOPIV, A+1+(JPIVSTART-1)*LDA, LDA, IPIVSTART, J, IPIV+1, 1 );
            IPIVSTART = IPIVSTART - NTOPIV;
            NTOPIV = NTOPIV * 2;
            JPIVSTART = JPIVSTART - NTOPIV;
        };

        // Permute U block to match L
        lapack::laswp( KCOLS, A+1+J*LDA, LDA, KSTART, J, IPIV+1, 1 );

        // Factor the current column
        if( A[J+(J-1)*LDA] != ZERO && !lapack::isnan( A[J+(J-1)*LDA] ) ) 
        {
            lapack::scal( M-J, ONE / A[J+(J-1)*LDA], A+J+1+(J-1)*LDA, 1 );
        }
        else if( *INFO  ==  0 ) 
        {
            *INFO = J;
        };

        // Solve for U block.
        lapack::trsm( "Left", "Lower", "No transpose", "Unit", KAHEAD,
                   KCOLS, ONE, A+KSTART+(KSTART-1)*LDA, LDA, A+KSTART+J*LDA, LDA );

        // Schur complement.
        lapack::gemm( "No transpose", "No transpose", M-J,
                    KCOLS, KAHEAD, NEGONE, A+J+1+(KSTART-1)*LDA, LDA,
                    A+KSTART+J*LDA, LDA, ONE, A+J+1+J*LDA, LDA );
    };

    // Handle pivot permutations on the way out of the recursion
    NPIVED = lapack::iand( NSTEP, -NSTEP );
    J = NSTEP - NPIVED;
    while ( J > 0 )
    {
        NTOPIV = lapack::iand( J, -J );
        lapack::laswp( NTOPIV, A+1+(J-NTOPIV)*LDA, LDA, J+1, NSTEP, IPIV+1, 1 );
        J = J - NTOPIV;
    };

    // If short and wide, handle the rest of the columns.
    if ( M  <  N ) 
    {
        //KCOLS = 0;
        lapack::laswp( N-M, A+1+(M+KCOLS)*LDA, LDA, 1, M, IPIV+1, 1 );
        lapack::trsm( "Left", "Lower", "No transpose", "Unit", M,
                N-M, ONE, A+1, LDA, A+1+(M+KCOLS)*LDA, LDA );
    };

    return;
};

template void BLAS_EXT_EXPORT
getrf_rec<s_type>(i_type M, i_type N, s_type *A, i_type LDA, i_type* IPIV, i_type *INFO );

template void BLAS_EXT_EXPORT
getrf_rec<d_type>(i_type M, i_type N, d_type *A, i_type LDA, i_type* IPIV, i_type *INFO );

template void BLAS_EXT_EXPORT
getrf_rec<z_type>(i_type M, i_type N, z_type *A, i_type LDA, i_type* IPIV, i_type *INFO );

template void BLAS_EXT_EXPORT
getrf_rec<c_type>(i_type M, i_type N, c_type *A, i_type LDA, i_type* IPIV, i_type *INFO );

};};