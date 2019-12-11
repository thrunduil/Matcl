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
#include "matcl-linalg/general/linalg_exception.h"
#include "matcl-linalg/decompositions/schur.h"
#include "matcl-linalg/decompositions/gschur_sym.h"
#include "matcl-linalg/decompositions/speigs.h"

namespace matcl
{

// return type of pschur function
using pschur_return     = tuple<Matrix, Matrix, bool>;

// return type of eigs function
using eigs_return       = tuple<Matrix, bool>;

// return type of svds function
using svds_return       = tuple<Matrix, Matrix, Matrix, bool>;

// perform schur decomposition of the matrix A in the form
//
//          A = U * T * U'
//
//  where U is unitary matrix, T is quasi-uppertriangular
//     [U, T] = schur(A, alg);
// A   - a square matrix
// alg - algorithm used for symmetric problems (used if A has
//       symmetric/hermitian structure), see schur_sym_alg for details
// U   - unitary matrix
// T   - quasi-uppertriangular matrix
//
// this function calls methods from class schur_decomposition
MATCL_LINALG_EXPORT mat_tup_2 schur(const Matrix& A, schur_sym_alg alg = schur_sym_alg::dc);
MATCL_LINALG_EXPORT mat_tup_2 schur(Matrix&& A, schur_sym_alg alg = schur_sym_alg::dc);

// perform generalized schur decomposition of the matrix pair (A,B) 
// in the form
//
//          Q' * A * Z = TA, 
//          Q' * B * Z = TB.
//
//  where Q, Z are unitary matrices, TA is quasi-uppertriangular, 
//  and TB is uppertriangular
//     [Q, Z, TA, TB] = gschur(A, B);
// A,B - square matrices
// Q,Z - unitary matrices
// TA  - quasi-uppertriangular matrix
// TB  - uppertriangular matrix
//
// this function calls methods from class gschur_decomposition
MATCL_LINALG_EXPORT mat_tup_4 gschur(const Matrix& A, const Matrix& B);
MATCL_LINALG_EXPORT mat_tup_4 gschur(Matrix&& A, const Matrix& B);
MATCL_LINALG_EXPORT mat_tup_4 gschur(const Matrix& A, Matrix&& B);
MATCL_LINALG_EXPORT mat_tup_4 gschur(Matrix&& A, Matrix&& B);

// perform complex schur decomposition of the matrix A in the form
//
//          A = U * T * U'
//
//  where U is unitary matrix, T is uppertriangular
//     [U, T] = schur(A, alg);
// A   - a square matrix
// alg - algorithm used for symmetric problems (used if A has
//       symmetric/hermitian structure), see schur_sym_alg for details
// U   - unitary matrix
// T   - uppertriangular matrix
//
// this function calls methods from class schur_decomposition
MATCL_LINALG_EXPORT mat_tup_2 schur_compl(const Matrix& A, schur_sym_alg alg = schur_sym_alg::dc);
MATCL_LINALG_EXPORT mat_tup_2 schur_compl(Matrix&& A, schur_sym_alg alg = schur_sym_alg::dc);

// perform complex generalized schur decomposition of the matrix pair (A,B) 
// in the form
//
//          Q' * A * Z = TA, 
//          Q' * B * Z = TB.
//
//  where Q, Z are unitary matrices, TA and TB are uppertriangular
//     [Q, Z, TA, TB] = gschur(A, B);
// A,B - square matrices
// Q,Z - unitary matrices
// TA  - uppertriangular matrix
// TB  - uppertriangular matrix
//
// // this function calls methods from class gschur_decomposition
MATCL_LINALG_EXPORT mat_tup_4 gschur_compl(const Matrix& A, const Matrix& B);
MATCL_LINALG_EXPORT mat_tup_4 gschur_compl(Matrix&& A, const Matrix& B);
MATCL_LINALG_EXPORT mat_tup_4 gschur_compl(const Matrix& A, Matrix&& B);
MATCL_LINALG_EXPORT mat_tup_4 gschur_compl(Matrix&& A, Matrix&& B);

// perform generalized Schur decomposition of the matrix pair (A,B) 
// where A and B are symmetric/hermitian and B is positive definite
// in the form:
//
// type:       eigen problem:      decomposition:
// A_B:        Ax =  lambda Bx.   A * V = B * V * D    and V' * B * V = I;
// AB:         ABx = lambda x.    A * B * V = V * D    and V' * B * V = I;
// BA:         BAx = lambda x.    B * A * V = V * D    and V' * inv(B) * V = I
//
// where V is square and D is diagonal matrix
// this solver reduce this problem to the standard symmetric eigenvalue
// problem implicitly inverting the matrix B, therefore B should be 
// sufficiently positive definite,
//
//     [V, D] = gen_sym_eigen(A, B, type);
//
// A   - symmetric/hermitian matrix
// B   - positive definite symmetric/hermitian matrix
// type- type of decomposition
// V   - normalized eigenvectors
// D   - diagonal matrix with eigenvalues sorted in asceding order
//
// not available for sparse matrices; available for band matrices if decomposition
// type is A_B
//
// this function calls methods from class gschur_sym_decomposition
MATCL_LINALG_EXPORT mat_tup_2 gen_sym_eigen(const Matrix& A, const Matrix& B, 
                                            gschur_sym_type type = gschur_sym_type::A_B);
MATCL_LINALG_EXPORT mat_tup_2 gen_sym_eigen(Matrix&& A, const Matrix& B, 
                                            gschur_sym_type type = gschur_sym_type::A_B);
MATCL_LINALG_EXPORT mat_tup_2 gen_sym_eigen(const Matrix& A, Matrix&& B, 
                                            gschur_sym_type type = gschur_sym_type::A_B);
MATCL_LINALG_EXPORT mat_tup_2 gen_sym_eigen(Matrix&& A, Matrix&& B, 
                                            gschur_sym_type type = gschur_sym_type::A_B);

// reorder Schur decomposition, i.e. put selected eigen-cluster in the
// leading position.
//     [UO, TO] = ordschur(U, T, select);
// U       - an unitary matrix of size N x N
// TA      - quasi-upper triangular matrix of size N x N
// select  - selection vector of size N x 1 with 1 indicating selection
//           of an eigenvalue
// UO      - an unitary matrix
// TO      - quasi uppertriangular matrix with selected eigenvalues in
//           upper left block
//
// this function calls methods from class schur_decomposition
MATCL_LINALG_EXPORT mat_tup_2 ordschur(const Matrix& U, const Matrix& T, const Matrix& select);

// reorder generalized Schur decomposition, i.e. put selected eigen-cluster in the
// leading position.
//     [QO, Z0, TAO, TBO] = ordschur(Q, Z, TA, TB, select);
// Q,Z     - unitary matrices of size N x N
// TA      - quasi-upper triangular matrix of size N x N
// TB      - upper triangular matrix of size N x N
// select  - selection vector of size N x 1 with 1 indicating selection
//           of an eigenvalue
// QO,ZO   - nitary matrices
// TAO     - quasi-upper triangular matrix of size N x N with selected eigenvalues in
//           upper left block
// TBO     - upper triangular matrix of size N x N with selected eigenvalues in
//           upper left block
//
// // this function calls methods from class gschur_decomposition
MATCL_LINALG_EXPORT mat_tup_4 ordgschur(const Matrix& Q, const Matrix& Z,
                                        const Matrix& S, const Matrix& T, const Matrix& select);

// perform schur decomposition of the matrix A and return eigenvalues
//     E = eig(A, alg);
// A   - a square matrix
// alg - algorithm used for symmetric problems (used if A has
//       symmetric/hermitian structure), see schur_sym_alg for details
// E   - vector of eigenvalues
//
// // this function calls methods from class schur_decomposition
MATCL_LINALG_EXPORT Matrix eig(const Matrix& A, schur_sym_alg alg = schur_sym_alg::dc);
MATCL_LINALG_EXPORT Matrix eig(Matrix&& A, schur_sym_alg alg = schur_sym_alg::dc);

// perform generalized schur decomposition of the matrix pair (A,B) 
// and return generalized eigenvalues
//     E = eig(A, B);
// A,B - square matrices
// E   - vector of generalized eigenvalues
//
// this function calls methods from class gschur_decomposition
MATCL_LINALG_EXPORT Matrix eig(const Matrix& A, const Matrix& B);
MATCL_LINALG_EXPORT Matrix eig(Matrix&& A, const Matrix& B);
MATCL_LINALG_EXPORT Matrix eig(const Matrix& A, Matrix&& B);
MATCL_LINALG_EXPORT Matrix eig(Matrix&& A, Matrix&& B);

// compute all eigenvalues of a symmetric tridiagonal matrix
//     E = eig_tridiag(diag, subdiag);
// diag    - elements on diagonal
// subdiag - elements on super- and subdiagonal
// E       - eigenvalues sorted increasingly
//
// notice that tridiagonal eigensolver is also used in eig function if
// matrix has tridiagonal symmetric/hermitian structure
//
// this function calls methods from class schur_decomposition
MATCL_LINALG_EXPORT Matrix  eig_tridiag(const Matrix& diag, const Matrix& subdiag);

// compute all eigenvalues and eigenvectors of a symmetric tridiagonal matrix
//     [E, V] = eig_tridiag2(diag, subdiag);
// diag    - elements on diagonal
// subdiag - elements on subdiagonal
// E       - eigenvalues sorted increasingly
// V       - associated eigenvectors
//
// notice that tridiagonal eigensolver is also used in eig function if
// matrix has tridiagonal symmetric/hermitian structure
//
// this function calls methods from class schur_decomposition
MATCL_LINALG_EXPORT mat_tup_2 eig_tridiag2(const Matrix& diag, const Matrix& subdiag);

// compute selected eigenvalues of a symmetric tridiagonal matrix
//     E = eigsel_tridiag_range(diag, subdiag, VL, VU);
//     E = eigsel_tridiag_index(diag, subdiag, IF, IL);
// diag    - elements on diagonal
// subdiag - elements on super- and subdiagonal
// VL, VU  - find eigenvalues in half-open interval (VL,VU]
// IF, IL  - find eigenvalues with index IF through IL in set of
//           all eigenvalues
// E       - selected eigenvalues sorted increasingly
MATCL_LINALG_EXPORT Matrix  eigsel_tridiag_range(const Matrix& diag, const Matrix& subdiag, 
                                Real VL, Real VU);
MATCL_LINALG_EXPORT Matrix  eigsel_tridiag_index(const Matrix& diag, const Matrix& subdiag, 
                                Integer IF, Integer IL);

// compute selected eigenvalues and eigenvectrs of a symmetric 
// tridiagonal matrix
//     [E,V] = eigsel_tridiag_range2(diag, subdiag, VL, VU);
//     [E,V] = eigsel_tridiag_index2(diag, subdiag, IF, IL);
// diag    - elements on diagonal
// subdiag - elements on super- and subdiagonal
// VL, VU  - find eigenvalues in half-open interval (VL,VU]
// IF, IL  - find eigenvalues with index IF through IL in set of
//           all eigenvalues
// E       - selected eigenvalues sorted increasingly
// V       - eigenvectors associated with eigenvalues E; eigenvectors are
//           computed using hess_right_eig and this function may fail; see
//           hess_right_eig for details
MATCL_LINALG_EXPORT mat_tup_2  eigsel_tridiag_range2(const Matrix& diag, const Matrix& subdiag, 
                                    Real VL, Real VU);
MATCL_LINALG_EXPORT mat_tup_2  eigsel_tridiag_index2(const Matrix& diag, const Matrix& subdiag, 
                                    Integer IF, Integer IL);

// compute selected eigenvalues of a symmetric matrix
//     E = eigsel_range(A, VL, VU);
//     E = eigsel_index(A, IF, IL);
// A       - symmetric/hermitian matrix (symmetry not checked)
// VL, VU  - find eigenvalues in half-open interval (VL,VU]
// IF, IL  - find eigenvalues with index IF through IL in set of
//           all eigenvalues sorted increasingly
// E       - selected eigenvalues sorted increasingly
//
// not available for sparse matrices
MATCL_LINALG_EXPORT Matrix  eigsel_range(const Matrix& A, Real VL, Real VU);
MATCL_LINALG_EXPORT Matrix  eigsel_index(const Matrix& A, Integer IF, Integer IL);

// compute selected eigenvalues and eigenvectors of a symmetric matrix
//     [E,V] = eigsel_range2(A, VL, VU);
//     [E,V] = eigsel_index2(A, IF, IL);
// A       - symmetric/hermitian matrix (symmetry not checked)
// VL, VU  - find eigenvalues in half-open interval (VL,VU]
// IF, IL  - find eigenvalues with index IF through IL in set of
//           all eigenvalues sorted increasingly
// E       - selected eigenvalues sorted increasingly
// V       - associated eigenvectors
//
// this function reduces A to tridiagonal matrix using hess2 function,
// calculate eigenvalues and eigenvectors using eig_tridiag2 function,
// and transforms eigenvectors; eigenvectors are computed using inverse
// iteration through hess_right_eig function, which may fail; see hess_right_eig
// for details;
// not available for sparse matrices
MATCL_LINALG_EXPORT mat_tup_2  eigsel_range2(const Matrix& A, Real VL, Real VU);
MATCL_LINALG_EXPORT mat_tup_2  eigsel_index2(const Matrix& A, Integer IF, Integer IL);

// perform generalized Schur decomposition of the matrix pair (A,B) 
// where A and B are symmetric/hermitian and B is positive definite
// in the form:
//
// type:       eigen problem:      decomposition:
// A_B:        Ax =  lambda Bx.   A * V = B * V * D    and V' * B * V = I;
// AB:         ABx = lambda x.    A * B * V = V * D    and V' * B * V = I;
// BA:         BAx = lambda x.    B * A * V = V * D    and V' * inv(B) * V = I
//
// where V is square and D is diagonal matrix
// this solver reduce this problem to the standard symmetric eigenvalue
// problem implicitly inverting the matrix B, therefore B should be 
// sufficiently positive definite,
//
//     E = eig_sym(A, B, type);
//
// A   - symmetric/hermitian matrix
// B   - positive definite symmetric/hermitian matrix
// type- type of decomposition
// E   - vector of generalized eigenvalues sorted in asceding order
//
// not available for sparse matrices; available for band matrices if decomposition
// type is A_B
//
// this function calls methods from class gschur_sym_decomposition
MATCL_LINALG_EXPORT Matrix eig_sym(const Matrix& A, const Matrix& B, 
                                    gschur_sym_type type = gschur_sym_type::A_B);
MATCL_LINALG_EXPORT Matrix eig_sym(Matrix&& A, const Matrix& B, 
                                    gschur_sym_type type = gschur_sym_type::A_B);
MATCL_LINALG_EXPORT Matrix eig_sym(const Matrix& A, Matrix&& B, 
                                    gschur_sym_type type = gschur_sym_type::A_B);
MATCL_LINALG_EXPORT Matrix eig_sym(Matrix&& A, Matrix&& B, 
                                    gschur_sym_type type = gschur_sym_type::A_B);

// find a few eigenvalues and eigenvectors of an operator A using Implicitly
// Restarted Arnoldi Method; this is a wrapper around Arpack; return matrices 
// U, and T such that
//
//          A * U = U * T
//
//  where U is unitary matrix, T is quasi-upper triangular
//     [U, T, conv] = pschur(A, k, ec, opts)
// A       - the linear operator of size N x N
// k       - number of eigenvalues to find
// ec      - eigenvalue cluster to find, see cluster_type for details;
//           on default returns eigenvalues with largest absolute value
// opts    - additional options; see namespace opt::speigs for details
// U       - unitary matrix of size N x kr, kr <= r + 1
// T       - quasi-upper triangular matrix of size krxkr, T is hermitian if
//           A is hermitian and real
// conv    = true: all eigenvalues converged withing specified tolorance
//         = false: some of eigenvalue did not converged
//
// this function calls methods from class pschur_decomposition
MATCL_LINALG_EXPORT pschur_return pschur(const linear_operator& A, const Integer k, 
                        cluster_type ec = cluster_type::LM, const options& opts = options());

// perform partial Schur decomposition of a linear operator A in the form
//     A * U = U * T,  U' * B * U = I                          (1)
// where T is upper quasi triangular and B is hermitian and semi positive 
// definite. From (1) we have that U spans invariant subspace of the operator
// A and eigenvalues of T are also eigenvalues of A. The partial Schur
// decomposition can be obtained by taking qr decompostion: Q*R = U, then
//     A * Q = Q * (R * T * R^-1), Q' * Q = I
//
// if the operator A is hermitian with respect to inner product given by the
// matrix B, then T is diagonal;  A is hermitian with respect to inner product
// given by B iff
//     B * A = A' * B, or equivalently < x,Ay > = < Ax,y >     (2)
// where <z,w> = z'Bw.
//
// representation (1) may help finding eigenvalues of the matrix pair (X,Y)
// using shift and invert method by taking:
//     A = inv[X - sigma*Y] * Y, B = Y or
//     A = inv[Y] * X, B = Y
// where sigma is a shift, such that X - sigma*Y is invertible, (such shift
// always exists if the matrix pair (X, Y) is regular); however if Y can be
// factored into a Cholesky factorization Y = LL', then the pschur decomposition
// should be used for A = inv[L] * X * inv[L']
//
//     [U, T, conv] = pbschur(A, B, herm, k, ec, opts)
// A       - the linear operator of size N x N
// B       - the linear operator of size N x N, B must be hermitian and
//           positive semi definite
// herm    = true: the operator A is hermitian with respect to inner product
//                 given by the operator B
//         = false: the operator A is general
// k       - number of eigenvalues to find
// ec      - eigenvalue cluster to find, see cluster_type for details;
//           on default returns eigenvalues with largest absolute value
// opts    - additional options; see namespace opt::speigs for details
// U       - a matrix of size N x kr, kr <= r + 1
// T       - quasi-upper triangular matrix of size kr x kr, T is hermitian if
//           A is hermitian and real
// conv    = true: all eigenvalues converged withing specified tolorance
//         = false: some of eigenvalue did not converged
//
// this function calls methods from class pbschur_decomposition
MATCL_LINALG_EXPORT pschur_return pbschur(const linear_operator& A, const linear_operator& B, 
                        bool hermitian, const Integer k, cluster_type ec = cluster_type::LM, 
                        const options& opts = options());

// perform pschur decomposition of the matrix A and return eigenvalues
//     [E, conv] = eigs(A, k, ec, opts)
// A       - the linear operator of size N x N
// k       - number of eigenvalues to find
// ec      - eigenvalue cluster to find, see cluster_type for details;
//           on default returns eigenvalues with largest absolute value
// opts    - additional options; see namespace opt::speigs for details
// E       - ke x 1 vector of eigenvalues, ke <= k + 1
// conv    = true: all eigenvalues converged withing specified tolorance
//         = false: some of eigenvalues did not converged
//
// this function calls methods from class pschur_decomposition
MATCL_LINALG_EXPORT eigs_return eigs(const linear_operator& A, const Integer k, 
                        cluster_type ec = cluster_type::LM, const options& opts = options());

// perform partial Schur decomposition of a linear operator A in the form
//     A * U = U * T,  U' * B * U = I                          (1)
// where T is upper quasi triangular and B is hermitian and semi positive 
// definite; from (1) we have that U spans invariant subspace of the operator
// A and eigenvalues of T are also eigenvalues of A. The partial Schur
// decomposition can be obtained by taking qr decompostion: Q*R = U, then
//     A * Q = Q * (R * T * R^-1), Q' * Q = I
//
// if the operator A is hermitian with respect to inner product given by the
// matrix B, then T is diagonal;  A is hermitian with respect to inner product
// given by B iff
//     B * A = A' * B, or equivalently < x,Ay > = < Ax,y >     (2)
// where <z,w> = z'Bw.
//
// representation (1) may help finding eigenvalues of the matrix pair (X,Y)
// using shift and invert method by taking:
//     A = inv[X - sigma*Y] * Y, B = Y or
//     A = inv[Y] * X, B = Y
// where sigma is a shift, such that X - sigma*Y is invertible, (such shift
// always exists if the matrix pair (X, Y) is regular), especially if Y is
// nearly singular or singular; however if Y can be factored into a Cholesky
// factorization Y = LL', then the pschur decomposition should be used for
// A = inv[L] * X * inv[L']
//
//     [E, conv] = eigs(A, B, herm, k, ec, opts)
// A       - the linear operator of size N x N
// B       - the linear operator of size N x N, Bm ust be hermitian and
//           positive semi definite
// herm    = true: the operator A is hermitian with respect to inner product
//                 given by the operator B
//         = false: the operator A is general
// k       - number of eigenvalues to find
// ec      - eigenvalue cluster to find, see cluster_type for details;
//           on default returns eigenvalues with largest absolute value
// opts    - additional options; see namespace opt::speigs for details
// E       - ke x 1 vector of eigenvalues of the operator A, ke <= k + 1
// conv    = true: all eigenvalues converged withing specified tolorance
//         = false: some of eigenvalue did not converged
//
// this function calls methods from class pbschur_decomposition
MATCL_LINALG_EXPORT eigs_return eigs(const linear_operator& A, const linear_operator& B, 
                        bool hermitian, const Integer k, cluster_type ec = cluster_type::LM, 
                        const options& opts = options());

// find a few singular values of an operator A using Implicitly
// Restarted Arnoldi Method; this function calls pschur_decomposition 
// for A*A' or A'*A operator, and if A is hermitian, then compute eigenvalues
// of A (ec must be SM or LM)
//     [E, conv] = svds1(A, k, ec, opts)
// A       - the linear operator of size M x N
// k       - number of singular values to find
// ec      - eigenvalue cluster to find, see cluster_type for details;
//           on default returns singular values with largest absolute value;
//           if one seeks for smallest singular values, then convergence
//           can be very slow
// opts    - additional options; see namespace opt::speigs for details
// E       - vector of singular values of the operator A of size ke x 1,
//           with ke <= k + 1
// conv    = true: all singular values converged withing specified tolorance
//         = false: some of singular value did not converged
MATCL_LINALG_EXPORT eigs_return svds1(const linear_operator& A, const Integer k, 
                        cluster_type ec = cluster_type::LM, const options& opts = options());

// find a few singular values and singular vectors of an operator A using Implicitly
// Restarted Arnoldi Method; this function calls pschur_decomposition 
// for A*A' or A'*A operator, and if A is hermitian, then compute eigenvalues
// of A (ec must be SM or LM)
//     [U, S, V, conv] = svds(A, k, ec, opts)
// A       - the linear operator of size M x N
// k       - number of singular values to find
// ec      - eigenvalue cluster to find, see cluster_type for details;
//           on default returns singular values with largest absolute value;
//           if one seeks for smallest singular values, then convergence
//           can be very slow
// opts    - additional options; see namespace opt::speigs for details
// U       - unitary matrix of size M x ke, U'*U = I with ke <= k + 1
// S       - diagonal matrix with singular values in decreasing order 
//           of size ke x ke and S = U' * A * V
// V       - unitary matrix of size N x ke, V'*V = I
// conv    = true: all singular values converged withing specified tolorance
//         = false: some of singular value did not converged
MATCL_LINALG_EXPORT svds_return svds(const linear_operator& A, const Integer k, 
                        cluster_type ec = cluster_type::LM, const options& opts = options());

};
