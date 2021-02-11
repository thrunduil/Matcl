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
#include "matcl-matrep/lib_functions/matrix_gen.h"
#include "matcl-linalg/linear_eq/linsolve_object.h"

namespace matcl
{

/// return types
using hess_ret      = tuple<unitary_matrix,Matrix>;
using hess_gen_ret  = tuple<Matrix,Matrix>;
using hess_gen2_ret = tuple<Matrix,Matrix,Matrix,Matrix>;

//----------------------------------------------------------------------
//                      HESSENBERG
//----------------------------------------------------------------------

/// Hessenberg decomposition of a square matrix
///     H = hess(A);
/// A   - square matrix
/// H   - Hessenberg matrix (a matrix with one nonzero subdiagonal),
///       is A is hermitian, then H is also hermitian
/// 
/// not available for sparse and band general matrices, but available for
/// band symmetric/hermitian matrices; for sparse symmetric matrices bandwith
/// reducing permutation is applied and depending on bandwith A is converted
/// to dense or banded matrix
MATCL_LINALG_EXPORT Matrix hess(const Matrix& A);
MATCL_LINALG_EXPORT Matrix hess(Matrix&& A);

/// Hessenberg decomposition of a square matrix
///     [U,H] = hess2(A);
/// satisfying A = U*H*U', and UU' = U'U = I
/// A   - square matrix
/// U   - unitary matrix
/// H   - Hessenberg matrix (a matrix with one nonzero subdiagonal),
///       is A is hermitian, then H is also hermitian
/// 
/// not available for sparse and band general matrices, but available for
/// band symmetric/hermitian matrices; for sparse symmetric matrices bandwith
/// reducing permutation is applied and depending on bandwith A is converted
/// to dense or banded matrix
MATCL_LINALG_EXPORT hess_ret hess2(const Matrix& A);
MATCL_LINALG_EXPORT hess_ret hess2(Matrix&& A);

//----------------------------------------------------------------------
//                   GENERALIZED HESSENBERG
//----------------------------------------------------------------------

/// generalized Hessenberg decomposition of a square matrices
///     [HA,HB] = gen_hess(A,B);
/// satisfying A = Q*H_A*ctrans(Z) and B = Q*H_B*ctrans(Z) for some unitary matrices Q, Z
/// 
/// A,B - square matrices
/// H_A - upper Hessenberg matrix
/// H_B - upper triangular matrix
/// 
/// not available for sparse and band matrices
MATCL_LINALG_EXPORT hess_gen_ret hess_gen(const Matrix& A, const Matrix& B);
MATCL_LINALG_EXPORT hess_gen_ret hess_gen(const Matrix& A, Matrix&& B);
MATCL_LINALG_EXPORT hess_gen_ret hess_gen(Matrix&& A, const Matrix& B);
MATCL_LINALG_EXPORT hess_gen_ret hess_gen(Matrix&& A, Matrix&& B);

/// generalized Hessenberg decomposition of a square matrices
///     [HA,HB,Q,Z] = gen_hess(A,B);
/// satisfying A = Q*H_A*ctrans(Z) and B = Q*H_B*ctrans(Z) for some unitary matrices Q, Z
/// A,B - square matrices
/// H_A - upper Hessenberg matrix
/// H_B - upper triangular matrix
/// Q,Z - unitary matrices
///
/// not available for sparse and band matrices
MATCL_LINALG_EXPORT hess_gen2_ret hess_gen2(const Matrix& A, const Matrix& B);
MATCL_LINALG_EXPORT hess_gen2_ret hess_gen2(const Matrix& A, Matrix&& B);
MATCL_LINALG_EXPORT hess_gen2_ret hess_gen2(Matrix&& A, const Matrix& B);
MATCL_LINALG_EXPORT hess_gen2_ret hess_gen2(Matrix&& A, Matrix&& B);

//----------------------------------------------------------------------
//                   HESSENBERG EIGENVECTORS
//----------------------------------------------------------------------

/// find left eigenvectors V of Hessenberg matrix associated with specified eigenvalues
/// using inverse iteration
///     V = hess_left_eig(H, eigs, init)
/// H    - upper Hessenberg matrix
/// E    - vector of size Mx1 containing eigenvalues of H; if H is a real matrix, then
///        complex eigenvalue and its conjugation must be stored in two consecutive
///        elements of E
/// init - matrix of size N x M containing starting vector of iteration for each eigenvector
///        init matrix can also be an empty matrix, then no external guess is used;
///        not used for symmetric/hermitian problems
/// V    - matrix containing eigenvectors of size N x M, where N is size of the matrix V
///        left eigenvector y for i-th eigenvalue l_i is defined as:
///                 ctrans(y) * A = l_i * ctrans(y)
///        matrix of eigenvectors Y satisfies
///                 Y' * H  = EM * Y', where EM = diag(E)
///        if some of eigenvector fails to converge, then exception of type hess_eig_failed
///        is throws; this exception contains further details
///
/// not available for sparse and general band matrices, but available for band 
/// symmetric/hermitian tridiagonal matrices
MATCL_LINALG_EXPORT Matrix hess_left_eig(const Matrix& H, const Matrix& E, const Matrix& init = zeros(0,0));
MATCL_LINALG_EXPORT Matrix hess_left_eig(const Matrix& H, const Matrix& E, Matrix&& init);

/// find right eigenvectors V of Hessenberg matrix associated with specified eigenvalues
/// using inverse iteration
///     V = hess_right_eig(H, eigs, init)
/// H    - upper Hessenberg matrix
/// E    - vector of size Mx1 containing eigenvalues of H; if H is a real matrix, then
///        complex eigenvalue and its conjugation must be stored in two consecutive
///        elements of E
/// init - matrix of size N x M containing starting vector of iteration for each eigenvector
///        init matrix can also be an empty matrix, then no external guess is used;
///        not used for symmetric/hermitian problems
/// V    - matrix containing eigenvectors of size N x M, where N is size of the matrix V
///        right eigenvector y for i-th eigenvalue l_i is defined as:
///                 A * y = l_i * y
///        matrix of eigenvectors Y satisfies
///                 A * Y = Y * EM, where EM = diag(E)
///        if some of eigenvector fails to converge, then exception of type hess_eig_failed
///        is throws; this exception contains further details
///
/// not available for sparse and general band matrices, but available for band 
/// symmetric/hermitian tridiagonal matrices
MATCL_LINALG_EXPORT Matrix hess_right_eig(const Matrix& H, const Matrix& E, const Matrix& init = zeros(0,0));
MATCL_LINALG_EXPORT Matrix hess_right_eig(const Matrix& H, const Matrix& E, Matrix&& init);

/// find left and right eigenvectors V of Hessenberg matrix associated with specified
/// eigenvalues using inverse iteration
///     [VL, VR] = hess_eig(H, eigs, left_init, right_init)
/// this function is equivalent to separate calls
///     VL = hess_left_eig(H, eigs, left_init);
///     VR = hess_right_eig(H, eigs, right_init);
///
/// not available for sparse and general band matrices, but available for band 
/// symmetric/hermitian tridiagonal matrices
MATCL_LINALG_EXPORT mat_tup_2 hess_eig(const Matrix& H, const Matrix& E, const Matrix& left_init = zeros(0,0),
                                         const Matrix& right_init = zeros(0,0));
MATCL_LINALG_EXPORT mat_tup_2 hess_eig(const Matrix& H, const Matrix& E, const Matrix& left_init, 
                                       Matrix&& right_init);
MATCL_LINALG_EXPORT mat_tup_2 hess_eig(const Matrix& H, const Matrix& E, Matrix&& left_init,
                                       const Matrix& right_init = zeros(0,0));
MATCL_LINALG_EXPORT mat_tup_2 hess_eig(const Matrix& H, const Matrix& E, Matrix&& left_init, Matrix&& right_init);

//----------------------------------------------------------------------
//                   HESSENBERG TO TRIANGULAR
//----------------------------------------------------------------------
/// reduce upper Hessenberg matrix to upper triangular matrix using plane
/// rotations
///     [U, T] = hess2triu(H, from_left);
/// such that H = U * T     if from_left = true     (1)
/// or        H = T * U     if from_left = false    (2)
/// H           - an upper Hessenberg matrix of size M x N
/// from_left   = true: do reduction (1), this is default value
///             = false: do reduction (2), matrix H must be square
/// U           - unitary matrix
/// T           - upper triangular matrix of size M x N
/// 
/// not available for sparse and band matrices
MATCL_LINALG_EXPORT hess_ret hess2triu(const Matrix& H, bool from_left = true);
MATCL_LINALG_EXPORT hess_ret hess2triu(Matrix&& H, bool from_left = true);

/// reduce lower Hessenberg matrix to lower triangular matrix using plane
/// rotations
///     [U, T] = lhess2tril(H, from_left);
/// such that H = U * T     if from_left = true     (1)
/// or        H = T * U     if from_left = false    (2)
/// H           - a lower Hessenberg matrix of size M x N
/// from_left   = true: do reduction (1), this is default value
///             = false: do reduction (2), matrix H must be square
/// U           - unitary matrix
/// T           - lower triangular matrix of size M x N
/// 
/// not available for sparse and band matrices
MATCL_LINALG_EXPORT hess_ret lhess2tril(const Matrix& H, bool from_left = true);
MATCL_LINALG_EXPORT hess_ret lhess2tril(Matrix&& H, bool from_left = true);

/// construct linsolve_obj from Hessenberg factors such that
///     A = U * H * ctrans(U)
///
///     lo  = linsolve_hess(A, U, H, opts)
/// A       - the factored matrix
/// U       - square unitary matrix
/// H       - square upper Hessenberg matrix
/// opts    - options, see namespace opt::linsolve for details
/// lo      - a linsolve_obj
/// 
/// if H is detected to be singular, then exception is thrown
MATCL_LINALG_EXPORT linsolve_obj linsolve_hess(const Matrix& A, const Matrix& U, const Matrix& H,
                                               const options& opts = options());
MATCL_LINALG_EXPORT linsolve_obj linsolve_hess(const Matrix& A, const unitary_matrix& U, const Matrix& H,
                                               const options& opts = options());

};