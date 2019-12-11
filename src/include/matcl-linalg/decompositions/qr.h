/*
 *  This file is a part of Matrix Computation Library (MATCL)
 *
 *  Copyright (c) Pawe³ Kowal 2018
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
#include "matcl-matrep/matrix/permvec.h"
#include "matcl-linalg/linear_eq/linsolve_object.h"

namespace matcl
{

// return type of qr_internal function
using qr_internal2_return = tuple<Matrix,Matrix>;

// return type of qr2 function
using qr2_return = tuple<unitary_matrix,Matrix>;

// return type of qr3 function
using qr3_return = tuple<unitary_matrix,Matrix,permvec>;

// internal qr decomposition
//     r = qr_internal(A);
// A   - matrix
// r   - matrix, such that R = triu(r), and QR = A, for some unitary
//       matrix Q
// 
// not available for sparse matrices
MATCL_LINALG_EXPORT Matrix qr_internal(const Matrix& A);
MATCL_LINALG_EXPORT Matrix qr_internal(Matrix&& A);

// internal qr decomposition
//     [r,tau] = qr_internal2(A);
// A   - matrix
// r   - matrix, such that R = triu(r); orthogonal matrix is stored
//       in lower triangular part of r in the form of householder
//       reflectors
// tau - the scalar factors of the Householder reflectors; one can
//       create unitary matrix Q such that Q*R = A by calling
//       Q = householder_to_unitary(r,tau)
// 
// not available for sparse matrices
MATCL_LINALG_EXPORT qr_internal2_return qr_internal2(const Matrix& A);
MATCL_LINALG_EXPORT qr_internal2_return qr_internal2(Matrix&& A);

// qr decomposition
//     R = qr(A);
// A   - matrix
// R   - upper triangular matrix, such that QR = A, for some unitary
//       matrix Q
MATCL_LINALG_EXPORT Matrix qr(const Matrix& A);
MATCL_LINALG_EXPORT Matrix qr(Matrix&& A);

// qr decomposition
//     [Q,R] = qr2(A, economy);
// such that
//     Q * R = A
// A   - matrix of size M x N
// economy
//     - boolean value
// R   - upper triangular matrix of size K x N, where K = min(M,N) if
//       economy = true and K = M otherwise
// Q   - unitary matrix of size M x K, where K is as above
MATCL_LINALG_EXPORT qr2_return qr2(const Matrix& A, bool economy = false);
MATCL_LINALG_EXPORT qr2_return qr2(Matrix&& A, bool economy = false);

// qr decomposition
//     [Q,R,p] = qr2(A, economy);
// such that
//     Q * R = A(:,p)
// A   - matrix of size M x N
// economy
//     - boolean value
// R   - upper triangular matrix of size K x N, where K = min(M,N) if
//       economy = true and K = M otherwise, elements on diagonal of R
//       are sorted decreasingly in absolute value
// Q   - unitary matrix of size M x K, where K is as above
// p   - permutation vector
// 
// not available for sparse and band matrices
MATCL_LINALG_EXPORT qr3_return qr3(const Matrix& A, bool economy = false);
MATCL_LINALG_EXPORT qr3_return qr3(Matrix&& A, bool economy = false);

// qr decomposition using givens rotarions
//     [Q,R] = qr_givens(A);
// such that
//     Q * R = A
// A   - matrix of size M x N, number of subdiagonals should be very low;
//       otherwise use qr2 function
// R   - upper triangular matrix of size M x N
// Q   - unitary matrix of size M x M
//
// not available for band and sparse matrices
MATCL_LINALG_EXPORT qr2_return qr_givens(const Matrix& A);
MATCL_LINALG_EXPORT qr2_return qr_givens(Matrix&& A);

// rq decomposition using givens rotarions
//     [Q,R] = rq_givens(A);
// such that
//     R * Q = A
// A   - matrix of size M x M, number of subdiagonals should be very low;
// R   - upper triangular matrix of size M x M
// Q   - unitary matrix of size M x M
//
// not available for band and sparse matrices
MATCL_LINALG_EXPORT qr2_return rq_givens(const Matrix& A);
MATCL_LINALG_EXPORT qr2_return rq_givens(Matrix&& A);

// lq decomposition using givens rotarions
//     [Q,L] = lq_givens(A);
// such that
//     L * Q = A
// A   - matrix of size M x M, number of superdiagonals should be very low;
// L   - lower triangular matrix of size M x M
// Q   - unitary matrix of size N x N
//
// not available for band and sparse matrices
MATCL_LINALG_EXPORT qr2_return lq_givens(const Matrix& A);
MATCL_LINALG_EXPORT qr2_return lq_givens(Matrix&& A);

// ql decomposition using givens rotarions
//     [Q,L] = ql_givens(A);
// such that
//     Q * L = A
// A   - matrix of size N x N, number of superdiagonals should be very low;
// L   - lower triangular matrix of size N x N
// Q   - unitary matrix of size N x N
//
// not available for band and sparse matrices
MATCL_LINALG_EXPORT qr2_return ql_givens(const Matrix& A);
MATCL_LINALG_EXPORT qr2_return ql_givens(Matrix&& A);

// given upper triangular matrix T, find qr factorization of a matrix 
//     B = T + U * Sigma * W', i.e. Q * R = B
// for diagonal Sigma
//     [Q, R] = qr_update(T, U, W, sigma)
// T       - upper triangular matrix of size M x N
// U       - matrix of size M x K
// W       - matrix of size N x K
// sigma   - vector of size K x 1, or scalar; Matrix Sigma is the diagonal matrix
//           constructed from the vector sigma or sigma * I if sigma is a scalar
// Q       - square unitary matrix of size M x M
// R       - upper triangular matrix of size M x N
//
// not available for band matrices, and sparse matrices
MATCL_LINALG_EXPORT qr2_return qr_update(const Matrix& R, const Matrix& U, const Matrix& W, 
                                       const Matrix& sigma);
MATCL_LINALG_EXPORT qr2_return qr_update(Matrix&& R, const Matrix& U, const Matrix& W, 
                                       const Matrix& sigma);

// given upper triangular square matrix T, find rq factorization of a matrix 
//     B = T + U * Sigma * W', i.e. R * Q = B
// for diagonal Sigma
//     [Q, R] = rq_update(T, U, W, sigma)
// T       - upper triangular matrix of size N x N
// U       - matrix of size N x K
// W       - matrix of size N x K
// sigma   - vector of size K x 1, or scalar; Matrix Sigma is the diagonal matrix
//           constructed from the vector sigma or sigma * I if sigma is a scalar
// Q       - square unitary matrix of size N x N
// R       - upper triangular matrix of size N x N
//
// not available for band matrices, and sparse matrices
MATCL_LINALG_EXPORT qr2_return rq_update(const Matrix& R, const Matrix& U, const Matrix& W, 
                                       const Matrix& sigma);
MATCL_LINALG_EXPORT qr2_return rq_update(Matrix&& R, const Matrix& U, const Matrix& W, 
                                       const Matrix& sigma);

// given QR factorization of a matrix A in the form A = Q*R, where T is upper triangular, 
// and Q is square unitary matrix find qr factorization of a matrix 
//     B = A + U * Sigma * W', i.e. Q2 * R2 = B
// for diagonal Sigma
//     [Q2, R2] = qr_update(Q, R, U, W, sigma)
// Q       - square unitary matrix of size M x M
// R       - upper triangular matrix of size M x N
// U       - matrix of size M x K
// W       - matrix of size N x K
// sigma   - vector of size K x 1, or scalar; Matrix Sigma is the diagonal matrix
//           constructed from the vector sigma or sigma * I if sigma is a scalar
// Q2      - square unitary matrix of size M x M
// R2      - upper triangular matrix of size M x N
//
// not available for band matrices, and sparse matrices
MATCL_LINALG_EXPORT qr2_return qr_update(const unitary_matrix& Q, const Matrix& R, const Matrix& U, const Matrix& W, 
                                       const Matrix& sigma);
MATCL_LINALG_EXPORT qr2_return qr_update(const unitary_matrix& Q, Matrix&& R, const Matrix& U, const Matrix& W, 
                                       const Matrix& sigma);
MATCL_LINALG_EXPORT qr2_return qr_update(const Matrix& Q, const Matrix& R, const Matrix& U, const Matrix& W, 
                                       const Matrix& sigma);
MATCL_LINALG_EXPORT qr2_return qr_update(const Matrix& Q, Matrix&& R, const Matrix& U, const Matrix& W, 
                                       const Matrix& sigma);

// given RQ factorization of a matrix A in the form A = R*Q, where T is upper triangular, 
// and Q is square unitary matrix find qr factorization of a matrix 
//     B = A + U * Sigma * W', i.e. R2 * Q2 = B
// for diagonal Sigma
//     [Q2, R2] = qr_update(Q, R, U, W, sigma)
// Q       - square unitary matrix of size M x M
// R       - upper triangular matrix of size M x N
// U       - matrix of size M x K
// W       - matrix of size N x K
// sigma   - vector of size K x 1, or scalar; Matrix Sigma is the diagonal matrix
//           constructed from the vector sigma or sigma * I if sigma is a scalar
// Q2      - square unitary matrix of size M x M
// R2      - upper triangular matrix of size M x N
//
// not available for band matrices, and sparse matrices
MATCL_LINALG_EXPORT qr2_return rq_update(const Matrix& Q, const Matrix& R, const Matrix& U, const Matrix& W, 
                                       const Matrix& sigma);
MATCL_LINALG_EXPORT qr2_return rq_update(const Matrix& Q, Matrix&& R, const Matrix& U, const Matrix& W, 
                                       const Matrix& sigma);
MATCL_LINALG_EXPORT qr2_return rq_update(const unitary_matrix& Q, const Matrix& R, const Matrix& U, const Matrix& W, 
                                       const Matrix& sigma);
MATCL_LINALG_EXPORT qr2_return rq_update(const unitary_matrix& Q, Matrix&& R, const Matrix& U, const Matrix& W, 
                                       const Matrix& sigma);

// construct linsolve_obj from qr factors such that
//     Q * R = A, or Q * R = A(:,p)
//
//     lo  = linsolve_qr(A, Q, R, opts)
//     lo  = linsolve_qr(A, Q, R, p, opts)
// A       - the factored matrix
// Q       - unitary matrix
// R       - upper triangular matrix
// p       - optional permutation vector
// opts    - options, see namespace opt::linsolve for details
// lo      - a linsolve_obj
// 
// if R is detected to be singular, then exception is thrown
MATCL_LINALG_EXPORT linsolve_obj linsolve_qr(const Matrix& A, const Matrix& Q, const Matrix& R,
                                             const options& opts = options());
MATCL_LINALG_EXPORT linsolve_obj linsolve_qr(const Matrix& A, const Matrix& Q, const Matrix& R, 
                                        const permvec& p, const options& opts = options());
MATCL_LINALG_EXPORT linsolve_obj linsolve_qr(const Matrix& A, const unitary_matrix& Q, const Matrix& R,
                                             const options& opts = options());
MATCL_LINALG_EXPORT linsolve_obj linsolve_qr(const Matrix& A, const unitary_matrix& Q, const Matrix& R, 
                                        const permvec& p, const options& opts = options());

// construct linsolve_obj from ql factors such that
//     Q * L = A
//
//     lo  = linsolve_ql(A, Q, L, opts)
// A       - the factored matrix
// Q       - unitary matrix
// L       - lower triangular matrix
// opts    - options, see namespace opt::linsolve for details
// lo      - a linsolve_obj
// 
// if R is detected to be singular, then exception is thrown
MATCL_LINALG_EXPORT linsolve_obj linsolve_ql(const Matrix& A, const Matrix& Q, const Matrix& L,
                                             const options& opts = options());
MATCL_LINALG_EXPORT linsolve_obj linsolve_ql(const Matrix& A, const unitary_matrix& Q, const Matrix& L,
                                             const options& opts = options());

// construct linsolve_obj from rq factors such that
//     R * Q = A
//
//     lo  = linsolve_rq(A, R, Q, opts)
// A       - the factored matrix
// R       - upper triangular matrix
// Q       - unitary matrix
// opts    - options, see namespace opt::linsolve for details
// lo      - a linsolve_obj
// 
// if R is detected to be singular, then exception is thrown
MATCL_LINALG_EXPORT linsolve_obj linsolve_rq(const Matrix& A, const Matrix& R, const Matrix& Q,
                                             const options& opts = options());
MATCL_LINALG_EXPORT linsolve_obj linsolve_rq(const Matrix& A, const Matrix& R, const unitary_matrix& Q,
                                             const options& opts = options());

// construct linsolve_obj from lq factors such that
//     L * Q = A
//
//     lo  = linsolve_lq(A, L, Q, opts)
// A       - the factored matrix
// L       - lower triangular matrix
// Q       - unitary matrix
// opts    - options, see namespace opt::linsolve for details
// lo      - a linsolve_obj
// 
// if R is detected to be singular, then exception is thrown
MATCL_LINALG_EXPORT linsolve_obj linsolve_lq(const Matrix& A, const Matrix& L, const Matrix& Q,
                                             const options& opts = options());
MATCL_LINALG_EXPORT linsolve_obj linsolve_lq(const Matrix& A, const Matrix& L, const unitary_matrix& Q,
                                             const options& opts = options());

// construct linsolve_obj for a matrix A using QR decomposition if required
//     lo  = linsolve_qr(A, pivot, opts)
// A       - a square matrix
// pivot   = false (default) - do not make pivoting, qr2, qr_givens, or lq_givens
//           function will be used depending on structure of the matrix A
//         = true - make pivoting; qr3 function will be used if required
// opts    - options controlling balancing; balancing is performed if option 
//           do_balancing (for pivot = false) or do_balancing_rr (for pivot = true)
//           is true using balance_gen function
// lo      - a linsolve_obj
// 
// if A is detected to be singular, then exception is thrown
MATCL_LINALG_EXPORT linsolve_obj linsolve_qr(const Matrix& A, bool pivot = false, 
                                             const options& opts = options());
MATCL_LINALG_EXPORT linsolve_obj linsolve_qr(Matrix&& A, bool pivot = false,
                                             const options& opts = options());

};