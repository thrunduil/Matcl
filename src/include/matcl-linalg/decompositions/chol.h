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
#include "matcl-matrep/matrix/permvec.h"
#include "matcl-linalg/linear_eq/linsolve_object.h"

namespace matcl
{

//---------------------------------------------------------------------------
//                      FACTORIZATION
//---------------------------------------------------------------------------

// return type of chol function
using chol_return_type      = tuple<Matrix,permvec>;

// return type of chol2 function
using chol2_return_type     = tuple<Matrix,permvec, Integer>;

// return type of chol_rr function
using chol_rr_return_type   = tuple<Matrix,permvec,Integer>;

// return type of ldl_tridiag function
using ldl_tridiag_return    = mat_tup_2;
using ldl_tridiag2_return   = tuple<Matrix,Matrix, Integer>;

// Cholesky decomposition of positive definite symmetric or hermitian matrix
// in the form
//         S'*S = A(p,p)       if upper = true,    S - upper triangular (1)
//         S * S' = A(p,p)     if upper = false,   S - lower triangular (2)
//
//     [S,p] = chol(A, upper, opt);
//
// A       - positive definite real symmetrix or complex hermitian matrix,
//           (symmetry is not checked)
// upper   - if true then factorization (1) is performed
//           otherwise factorization (2) is performed
// opt     - options controlling factorization; see namespace opt::linsolve
//           for details
// S       - upper triangular matrix such than S'*S = A(p,p) (if upper = true)
//           or lower triangular matrix with S * S' = A(p,p) (if upper = false)
// p       - permutation vector
// 
// throws error error_nonposdef if A is not positive definite;
// for sparse matrices CHOLMOD is used; CHOLMOD is not available for Float
// and Float_complex values, in this case values are converted first to Real
// or Complex
MATCL_LINALG_EXPORT chol_return_type chol(const Matrix& A, bool upper = false,
                                          const options& opts = options());
MATCL_LINALG_EXPORT chol_return_type chol(Matrix&& A, bool upper = false,
                                          const options& opts = options());

// Cholesky decomposition of positive definite symmetric or hermitian matrix
// in the form
//         S'*S = A(p,p) + E    if upper = true,    S - upper triangular (1)
//         S * S' = A(p,p) + E  if upper = false,   S - lower triangular (2)
// where E = blkdiag(Z,EE) is block diagonal matrix with the first block Z 
// being zero matrix of size NxN, and the second block EE is some matrix of 
// size (M-N)x(M-N)
//
//     [S,p,N] = chol2(A, upper, opt);
//
// A       - positive definite real symmetrix or complex hermitian matrix of 
//           size M x M, (symmetry is not checked)
// upper   - if true then upper triangular matrix is returned
//           otherwise lower triangular matrix is returned
// opt     - options controlling factorization; see namespace opt::linsolve
//           for details
// S       - upper triangular matrix of size NxM such that
//           S'*S = A(p,p) + E (if upper = true)
//           or lower triagular matrix of size MxN such that
//           S * S' = A(p,p) + E (if upper = false)
// p       - permutation vector
// N       - position at which Cholesky decomposition breaks down; if A is
//           sufficiently positive definite, then n = A.rows()
// 
// this version does not throw exception if S is not positive definite
// for sparse matrices CHOLMOD is used; CHOLMOD is not available for Float
// and Float_complex values, in this case values are converted first to Real
// or Complex
MATCL_LINALG_EXPORT chol2_return_type chol2(const Matrix& A, bool upper = false,
                                            const options& opts = options());
MATCL_LINALG_EXPORT chol2_return_type chol2(Matrix&& A, bool upper = false,
                                            const options& opts = options());

// Cholesky decomposition of positive definite symmetric or hermitian tridiagonal 
// matrix in the form
//         S' * D * S = A  if upper = true,    S - unit upper triangular
//         S * D * S' = A  if upper = false,   S - unit lower triangular
//
//     [S,D]   = ldl_tridiag(D0, D1, upper);
//     [S,D,k] = ldl_tridiag2(D0, D1, upper);
//
// ldl_tridiag throws error error_nonposdef if A is not positive definite;
// ldl_tridiag2 does not throw exception
//
// D0      - vector of size Nx1 of elements on main diagonal of A
// D1      - vector of size (N-1)x1 of elements on subdiagonal of A
// upper   - if true then factorization (1) is performed
//           otherwise factorization (2) is performed
// S       - unit upper triangular matrix if upper = true and unit lower triagular
//           otherwise
// D       - diagonal matrix
// k       - additional return of ldl_tridiag2 version; this version does not 
//           throw exception if A is not positive definite but instead return 
//           position at which Cholesky decomposition breaks down; if A is
//           sufficiently positive definite, then k = N
MATCL_LINALG_EXPORT ldl_tridiag_return  ldl_tridiag(const Matrix& D0, const Matrix& D1, 
                                                bool upper = false);
MATCL_LINALG_EXPORT ldl_tridiag2_return ldl_tridiag2(const Matrix& D0, const Matrix& D1, 
                                                bool upper = false);

// Cholesky decomposition of positive semi-definite symmetric or hermitian matrix
// with pivoting in the form
//         S' * D * S = A  if upper = true,    S - unit upper triangular
//         S * D * S' = A  if upper = false,   S - unit lower triangular
//
//     [S, p, rank] = chol_rr(A, upper, tol); 
//
// A       - positive semi-definite real symmetrix or complex hermitian matrix,
//           (symmetry is not checked)
// upper   - if true then upper triangular matrix is returned
//           otherwise lower triangular matrix is returned
// tol     - optional, tolerance used to detect singular values; if tol < 0, then 
//           default tolerance is used equal to 10.0 * eps * |A|_F / sqrt(N), where
//           |A|_F is the Frobenius norm of A; element v is considered as zero 
//           if |v| < tol
// S       - upper triangular matrix such than S'*S = A(p,p) (if upper = true)
//           or lower triangular matrix with S * S' = A(p,p) (if upper = false)
// p       - permutation vector
// rank    - estimated rank of the matrix A
// 
// chol_rr does not produce error when the matrix A appears not be positive semi-definite,
// but instead quits and returns current estimation of rank
// for sparse matrices LUSOL is used; banded matrices are converted to dense
// or sparse matrices
MATCL_LINALG_EXPORT chol_rr_return_type chol_rr(const Matrix& A, bool upper = false, Real tol = -1.);
MATCL_LINALG_EXPORT chol_rr_return_type chol_rr(Matrix&& A, bool upper = false, Real tol = -1.);

//---------------------------------------------------------------------------
//                      LINSOLVE
//---------------------------------------------------------------------------

// construct linsolve_obj from Cholesky factors
//     lo  = linsolve_chol(A, S, p, upper, opts)
// A       - the factored matrix
// S       - upper triangular matrix such than S'*S = A(p,p) (if upper = true)
//           or lower triangular matrix with S * S' = A(p,p) (if upper = false)
// p       - permutation vector
// upper   = true:     S is upper triangular
//         = false:    S is lower triangular
// opts    - options, see namespace opt::linsolve for details
// lo      - a linsolve_obj
MATCL_LINALG_EXPORT linsolve_obj linsolve_chol(const Matrix& A, const Matrix& S, 
                                    const permvec& p, bool upper = false, const options& opts = options());

// construct linsolve_obj for a matrix A using Cholesky decomposition if required
//     lo  = linsolve_chol(A, pivot, opts)
// A       - a square matrix hermitian positive definite matrix
// pivot   = false (default) - do not make pivoting, chol function will be used
//         = true - make pivoting; chol_rr function will be used if required
// opts    - options controlling balancing; balancing is performed if option 
//           do_balancing (for pivot = false) and do_balancing_rr (for pivot = true)
//           is true using balance_posdef function
// lo      - a linsolve_obj
MATCL_LINALG_EXPORT linsolve_obj linsolve_chol(const Matrix& A, bool pivot = false,
                                               const options& = options());
MATCL_LINALG_EXPORT linsolve_obj linsolve_chol(Matrix&& A, bool pivot = false,
                                               const options& = options());

// construct linsolve_obj for a tridiagonal hermitian, positive definite matrix A 
// using Cholesky decomposition
//     lo  = linsolve_chol_tridiag(D0, D1, opts)
// D0      - vector of size Nx1 of elements on main diagonal of A
// D1      - vector of size (N-1)x1 of elements on subdiagonal of A
// opts    - options, see namespace opt::linsolve for details
// lo      - a linsolve_obj
MATCL_LINALG_EXPORT linsolve_obj linsolve_chol_tridiag(const Matrix& D0, const Matrix& D1, 
                                    const options& opts = options());
MATCL_LINALG_EXPORT linsolve_obj linsolve_chol_tridiag(const Matrix& D0, Matrix&& D1, 
                                    const options& opts = options());
MATCL_LINALG_EXPORT linsolve_obj linsolve_chol_tridiag(Matrix&& D0, const Matrix& D1, 
                                    const options& opts = options());
MATCL_LINALG_EXPORT linsolve_obj linsolve_chol_tridiag(Matrix&& D0, Matrix&& D1, 
                                    const options& opts = options());

// construct linsolve_obj from factors forming the Cholesky decomposition of
// tridiagonal, hermitian, positive definite matrix
//     lo  = linsolve_chol_tridiag_fac(A, S, D, upper, opts)
// A       - the factored matrix
// S       - unit upper triangular matrix such than S'*D*S = A (if upper = true)
//           or unit lower triangular matrix with S*D*S' = A
// D       - diagonal matrix
// upper   = true:     S is upper triangular
//         = false:    S is lower triangular
// opts    - options, see namespace opt::linsolve for details
// lo      - a linsolve_obj
MATCL_LINALG_EXPORT linsolve_obj linsolve_chol_tridiag_fac(const Matrix& A, const Matrix& S, 
                                    const Matrix& D, bool upper = false, const options& opts = options());

// construct linsolve_obj from factors forming the Cholesky decomposition of
// tridiagonal, hermitian, positive definite matrix A
//     lo  = linsolve_chol_tridiag_fac(D0, D1, S, D, upper, opts)
// D0      - vector of size Nx1 of elements on main diagonal of the factored matrix A
// D1      - vector of size (N-1)x1 of elements on subdiagonal of the factored matrix A
// S       - unit upper triangular matrix such than S'*D*S = A (if upper = true)
//           or unit lower triangular matrix with S*D*S' = A
// D       - diagonal matrix
// upper   = true:     S is upper triangular
//         = false:    S is lower triangular
// lo      - a linsolve_obj
MATCL_LINALG_EXPORT linsolve_obj linsolve_chol_tridiag_fac(const Matrix& D0, const Matrix& D1,
                                    const Matrix& S, const Matrix& D, bool upper = false,
                                    const options& opts = options());

//---------------------------------------------------------------------------
//                      UPDATE
//---------------------------------------------------------------------------

// given Cholesky factorization of a matrix A, T' * T or T * T', find Cholesky
// factorization of a matrix B constructed from the matrix A by removing
// rows and columns specified by the matrix K.
//     S = chol_remove(T, upper, K)
// T       - upper or lower triangular square matrix
// upper   = true: T is upper triangular and original matrix A is equal to T' * T
//         = false: T is lower triangular and original matrix A is equal to T * T'
// K       - a matrix with indices of a row and column of A to remove; invalid 
//           positions are ignored
// S       - upper or lower triangular matrix forming the Cholesky factor of the
//           matrix B obtained by removing row and column of A
//
// not available for band matrices, and sparse matrices if upper = true
MATCL_LINALG_EXPORT Matrix chol_remove(const Matrix& T, const matcl::Matrix& K, bool upper = false);
MATCL_LINALG_EXPORT Matrix chol_remove(Matrix&& T, const matcl::Matrix& K, bool upper = false);

// given Cholesky factorization of a matrix A, T*T' or T'*T, find cholesky
// factorization of a matrix B = A + W * Sigma * W' for diagonal Sigma
//     S = chol_update(T, W, upper, sigma)
// T       - upper or lower triangular square matrix of size N x N
// W       - matrix of size N x K
// upper   = true: T is upper triangular and original matrix A is equal to T' * T
//         = false: T is lower triangular and original matrix A is equal to T * T'
// sigma   - vector of size K x 1, or scalar; Matrix Sigma is the diagonal matrix
//           constructed from the vector sigma or sigma * I if sigma is a scalar;
//           sigma must contain real values
// S       - upper or lower triangular matrix forming the Cholesky factor of the
//           matrix B
//
// not available for band matrices, and sparse matrices if upper = true
MATCL_LINALG_EXPORT Matrix chol_update(const Matrix& T, const Matrix& W, 
                                       const Matrix& sigma, bool upper = false);
MATCL_LINALG_EXPORT Matrix chol_update(Matrix&& T, const Matrix& W, 
                                       const Matrix& sigma, bool upper = false);

};
