/*
 *  This file is a part of Matrix Computation Library (MATCL)
 *
 *  Copyright (c) Pawe³ Kowal 2018 - 2021
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
#pragma once

#include "matcl-matrep/matrix/matrix.h"
#include "matcl-linalg/general/config_linalg.h"
#include "matcl-linalg/special_matrices/unitary_matrix.h"
#include "matcl-linalg/linear_eq/linsolve_object.h"

namespace matcl
{

enum class svd_algorithm
{
    qr,     // implicit zero-shift QR algorithm
    dc      // divide-conquer
};

// return type of make_bidiagonal function
using make_bidiagonal_ret = tuple<unitary_matrix,matcl::Matrix,unitary_matrix>;

// return type of gsvd function
using gsvd_ret  = tuple<Matrix,Matrix,Matrix,Matrix,Matrix,Matrix,Integer,Integer>;

// svd decomposition
//     [U,S,V] = svd(A, economy, alg)
// such that
//      U*S*V' = A, U'*U = I, V'*V = I
// A   - matrix of size M x N
// economy
//     - boolean value
// alg - algorithm used
// U   - unitary matrix of size M x M if economy = false or MxK if 
//       economy = true, where K = min(M,N)
// S   - diagonal matrix with singular values in decreasing order 
//       of size MxN if economy = false of K x K if economy = true
// V   - unitary matrix of size N x N if economy = false or NxK
//       if economy = true
//
// not available for sparse and band matrices
MATCL_LINALG_EXPORT mat_tup_3 svd(const Matrix& A, bool economy = true, svd_algorithm alg = svd_algorithm::dc);
MATCL_LINALG_EXPORT mat_tup_3 svd(Matrix&& A, bool economy = true, svd_algorithm alg = svd_algorithm::dc);

// svd decomposition computing only singular values
//     S = svd(A, alg)
// A   - matrix
// alg - algorithm used
// S   - diagonal matrix with singular values in decreasing order of size
//       K x 1, where K = min(M,N), where M, N are number of rows and columns
//       of A
//
// not available for sparse and band matrices
MATCL_LINALG_EXPORT Matrix svd1(const Matrix& A, svd_algorithm alg = svd_algorithm::dc);
MATCL_LINALG_EXPORT Matrix svd1(Matrix&& A, svd_algorithm alg = svd_algorithm::dc);

// return orthonormal basis of the null space of the matrix A; if right = true, 
// then right null space N is returned such that A * N = 0; if right = false, then
// left null space N is returned such that N * A = 0; tol is used to detect zero
// singular values, if tol < 0, then tol = max(A.rows(), A.cols()) * norm(A, 2) * eps
MATCL_LINALG_EXPORT Matrix null(const Matrix& A, bool right = true, Real tol = -1.0);
MATCL_LINALG_EXPORT Matrix null(Matrix&& A, bool right = true, Real tol = -1.0);

// return orthonormal basis of the range of A; if right = true, then the basis of 
// the right range is returned such that A * N has full rank; if right = false, then
// the basis of the left range is returned such that N * A has full rank; tol is used
// to detect zero singular values, if tol < 0, then tol = max(A.rows(), A.cols()) * 
// norm(A, 2) * eps
MATCL_LINALG_EXPORT Matrix orth(const Matrix& A, bool right = true, Real tol = -1.0);
MATCL_LINALG_EXPORT Matrix orth(Matrix&& A, bool right = true, Real tol = -1.0);

// estimate the rank of the matrix A, i.e. the number of linearly independent rows
// or columns of A; tol is used to detect zero singular values, if tol < 0, then 
// tol = max(A.rows(), A.cols()) * norm(A, 2) * eps
MATCL_LINALG_EXPORT Integer rank(const Matrix& A, Real tol = -1.0);
MATCL_LINALG_EXPORT Integer rank(Matrix&& A, Real tol = -1.0);

// compute orthonornal basis of the null space or range based on SVD decomposition
class MATCL_LINALG_EXPORT operator_spaces
{
    private:
        Matrix      m_U;
        Matrix      m_S;
        Matrix      m_V;

    public:
        // compute spaces for the scalar 1.0
        operator_spaces();

        // compute spaces for the matrix A
        operator_spaces(const Matrix& A);
        operator_spaces(Matrix&& A);

        // standard destructor
        ~operator_spaces();

        // return orthonormal basis of the right null space N of the matrix A
        // such that A * N = 0; tol is used to detect zero singular values, 
        // if tol < 0, then default_tol() is used
        Matrix      null_right(Real tol = -1.0) const;

        // return orthonormal basis of the left null space N of the matrix A
        // such that N * A = 0; tol is used to detect zero singular values, 
        // if tol < 0, then default_tol() is used
        Matrix      null_left(Real tol = -1.0) const;

        // return orthonormal basis of the right range of the matrix A such that 
        // A * N has full rank; tol is used to detect zero singular values, 
        // if tol < 0, then default_tol() is used
        Matrix      orth_right(Real tol = -1.0) const;

        // return orthonormal basis of the right range of the matrix A such that 
        // N * A has full rank; tol is used to detect zero singular values, 
        // if tol < 0, then default_tol() is used
        Matrix      orth_left(Real tol = -1.0) const;

        // estimate the rank of the matrix A, i.e. the number of linearly independent
        // rows or columns of A; tol is used to detect zero singular values, 
        // if tol < 0, then default_tol() is used
        Integer     rank(Real tol = -1.0) const;

        // return default tolerance tol = max(A.rows(), A.cols()) * norm(A, 2) * eps
        Real        default_tol() const;
};

// compute the generalized singular value decomposition (GSVD) of an M x N matrix A
// and P x N matrix B:
//         U' * A * Q	= D1 * R,       R = [Z, RR]
//         V' * B * Q	= D2 * R
//
// [U,V,Q,D1,D2,R,K,L] = gsvd(A, B)
// A:      matrix of size M x N
// B:      matrix of size P x N
// U:      unitary matrix of size  M x M
// V:      unitary matrix of size  P x P
// Q:      unitary matrix of size  N x N
// D1:     real diagonal matrix of size M * (K+L) and D1.diag(0) = [J; C] J is a vector
//         of length Kx1 with all elements 1 and C is a vector of length min(L, M-K)
// D2:     real matrix of size P * (K+L) where diagonal K is the only nonzero diagonal
//         and D2.diag(K) = [S;J] where J is a vector of length max(0, K+L-M) with
//         all elements 1 and S is a vector of length min(L, M-K);
//         additionally C.^2 + S.^2 = 1
// R:      matrix of size (K+L) * N, where RR submatrix has size (K+L) x (K+L)
//         and is nonsingular upper triangular matrix, and Z is zero matrix
// K,L:    integer scalars, K + L is the effective numerical rank of the matrix [A;B]
//
// if B is an N-by-N nonsingular matrix, then the GSVD of A and B implicitly gives the
// SVD of A*inv(B):
//                 A*inv(B) = U * (D1*inv(D2)) * V'
// if [A;B] has orthonormal columns, then the GSVD of A and B is also equal to the CS
// decomposition of A and B;
// the GSVD can be used to derive the solution of the eigenvalue problem:
//                 (A' * A) x = lambda * (B' * B) x.
// not available for sparse and band matrices
MATCL_LINALG_EXPORT gsvd_ret gsvd(const Matrix& A, const Matrix& B);
MATCL_LINALG_EXPORT gsvd_ret gsvd(Matrix&& A, const Matrix& B);
MATCL_LINALG_EXPORT gsvd_ret gsvd(const Matrix& A, Matrix&& B);
MATCL_LINALG_EXPORT gsvd_ret gsvd(Matrix&& A, Matrix&& B);

// compute generalized singular values of an M x N matrix A and P x N matrix B:
// E = gsvd1(A, B)
// A:      matrix of size M x N
// B:      matrix of size P x N
// E:      vector of size N x 1 with generalized singular values; if A and B
//         are square and B is invertible then also E = svd1(A * inv(B)) but
//         sorted differently
// not available for sparse and band matrices
MATCL_LINALG_EXPORT Matrix gsvd1(const Matrix& A, const Matrix& B);
MATCL_LINALG_EXPORT Matrix gsvd1(Matrix&& A, const Matrix& B);
MATCL_LINALG_EXPORT Matrix gsvd1(const Matrix& A, Matrix&& B);
MATCL_LINALG_EXPORT Matrix gsvd1(Matrix&& A, Matrix&& B);

// compute the singular value decomposition of a 2-by-2 matrix
//   [  A_11   A_12  ]
//   [  A_21   A_22  ]
// on return, abs(sig_max) is the larger singular value, abs(sig_min) is the
// smaller singular value, and (cos_l,sin_l) and (cos_r,sin_r) are the left and
// right singular vectors for abs(sig_max), giving the decomposition
//
//   [ cos_l  sin_l ] [ A_11  A_12  ] [ cos_r -sin_r ]  =  [ sig_max   0      ]
//   [-sin_l' cos_l ] [ A_21  A_22  ] [ sin_r' cos_r ]     [  0       sig_min ]
//
//     svd_22(A_11, A_11, A_21, A_22, cos_r, sin_r, cos_l, sin_l, 
//             sig_max, sig_min)
// A_11, A_12, A_21, A_22  - (input) elements of the 2 x 2 matrix
// cos_l, sin_l            - (output) cosine and sine of left plane rotation
// cos_r, sin_r            - (output) cosine and sine of right plane rotation
// sig_max, sig_min        - (output) defines smaller and larger singular values, 
//                           given by abs(sig_min), abs(sig_max)
//
template<class V, class Enable = typename md::enable_if_float<V,void>::type, 
        class VR = typename md::real_type<V>::type>
MATCL_LINALG_EXPORT void svd_22(const V& A_11, const V& A_12, const V& A_21, const V& A_22, 
                                VR& cos_l, V& sin_l, VR& cos_r, V& sin_r, V& sig_max, V& sig_min);

// compute the singular values of a 2-by-2 matrix
//   [  A_11   A_12  ]
//   [  A_21   A_22  ]
//
//     svd_22(A_11, A_11, A_21, A_22, sig_max, sig_min)
// A_11, A_12, A_21, A_22  - (input) elements of the 2 x 2 matrix
// sig_max, sig_min        - (output) larger and smaller singular values
template<class V, class Enable = typename md::enable_if_float<V,void>::type, 
        class VR = typename md::real_type<V>::type>
MATCL_LINALG_EXPORT void svd_22(const V& A_11, const V& A_12, const V& A_21, const V& A_22, 
                                VR& sig_max, VR& sig_min);

// reduce matrix to bidiagonal matrix by unitary thransformations
//     [U,R,V] = make_bidiagonal(A)
// such that
//     A = U * R * V'
// A   - matrix of size MxN
// U   - unitary matrix of size MxM
// R   - bidiagonal matrix, if M >= N, then R is upper bidiagonal
//       else R is lower bidiagonal of size M x N
// V   - unitary matrix of size NxN
//
// not available for sparse and band matrices
MATCL_LINALG_EXPORT make_bidiagonal_ret make_bidiagonal(const Matrix& A);
MATCL_LINALG_EXPORT make_bidiagonal_ret make_bidiagonal(Matrix&& A);

// construct linsolve_obj from SVD factors such that
//     U * S * V' = A, U' * U = I, V' * V = I
//
//     lo  = linsolve_svd(A, U, S, V, opts)
// A       - the factored matrix
// U       - unitary matrix
// S       - square diagonal matrix
// V       - unitary matrix
// opts    - options, see namespace opt::linsolve for details
// lo      - a linsolve_obj
// 
// if S is detected to be singular, then exception is thrown
MATCL_LINALG_EXPORT linsolve_obj linsolve_svd(const Matrix& A, const Matrix& U, const Matrix& S, const Matrix& V,
                                              const options& opts = options());
MATCL_LINALG_EXPORT linsolve_obj linsolve_svd(const Matrix& A, const unitary_matrix& U, const Matrix& S, 
                                              const unitary_matrix& V, const options& opts = options());

// construct linsolve_obj from bidiagonal factorization such that
//     A = U * R * V'
//
//     lo  = linsolve_bidiag(A, U, R, V, opts)
// A       - the factored matrix
// U       - unitary matrix
// R       - bidiagonal matrix, if M >= N, then R is upper bidiagonal
//           else R is lower bidiagonal of size M x N
// V       - unitary matrix
// opts    - options, see namespace opt::linsolve for details
// lo      - a linsolve_obj
// 
// if R is detected to be singular, then exception is thrown
MATCL_LINALG_EXPORT linsolve_obj linsolve_bidiagonal(const Matrix& A, const Matrix& U, 
                                        const Matrix& R, const Matrix& V, const options& opts = options());
MATCL_LINALG_EXPORT linsolve_obj linsolve_bidiagonal(const Matrix& A, const unitary_matrix& U,
                                        const Matrix& R, const unitary_matrix& V, const options& opts = options());

// construct linsolve_obj for a matrix A using SVD decomposition if required
//     lo  = linsolve_svd(A, opts)
// A       - a square matrix
// opts    - options, see namespace opt::linsolve for details
// lo      - a linsolve_obj
// 
// if A is detected to be singular, then exception is thrown
MATCL_LINALG_EXPORT linsolve_obj linsolve_svd(const Matrix& A, const options& opts = options());
MATCL_LINALG_EXPORT linsolve_obj linsolve_svd(Matrix&& A, const options& opts = options());

};