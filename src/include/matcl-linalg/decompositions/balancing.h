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
#include "matcl-linalg/linear_eq/linsolve_object.h"
#include "matcl-linalg/special_matrices/unitary_matrix.h"
#include "matcl-matrep/details/enablers.h"

namespace matcl
{

namespace md = matcl::details;

//-------------------------------------------------------------------------
//                              EIGENVALUES
//-------------------------------------------------------------------------

// return diagonal scaling to improve accuracy of computation of eigenvalues of a
// matrix A
//     [AA, D] = balance_eig(A);
// A       - a square matrix
// AA      - balanced matrix
//             AA  = diag(D)^-1 * A * diag(D)
// D       - scaling vector containing values being power of 2
// 
// uses the Parlett-Reinsch Algorithm for the second norm with the correction 
// proposed by R. James, J. Langou, B.R. Lowery, "On matrix balancing and eigenvector
// computation"
MATCL_LINALG_EXPORT mat_tup_2 balance_eig(const Matrix& A);
MATCL_LINALG_EXPORT mat_tup_2 balance_eig(Matrix&& A);

// return diagonal scaling to reduce maginitude of offdiagonal elements without
// changing diagonal elements
//     [AA, D] = balance_offdiag(A);
// A       - a square matrix
// AA      - balanced matrix
//             AA  = diag(D)^-1 * A * diag(D)
// D       - scaling vector
MATCL_LINALG_EXPORT mat_tup_2 balance_offdiag(const Matrix& A);
MATCL_LINALG_EXPORT mat_tup_2 balance_offdiag(Matrix&& A);

// return diagonal scaling to improve accuracy of computation of eigenvalues of a
// matrix pair (A,B)
//     [AA, BB, Dl, Dr] = balance_eig(A, B);
// A, B    - square matrices
// AA, BB  - balanced matrices, 
//             AA = diag(Dl) * A * diag(Dr)
//             BB = diag(Dl) * B * diag(Dr)
// Dl      - row scaling vector containing values being power of 2
// Dr      - column scaling vector containing values being power of 2
// 
// use scaling technique proposed by Lemonnier and Van Dooren:
//
//   Lemonnier, Damien, and Paul Van Dooren. "Balancing regular matrix pencils." 
//   SIAM journal on matrix analysis and applications 28.1 (2006): 253-263.
MATCL_LINALG_EXPORT mat_tup_4 balance_eig(const Matrix& A, const Matrix& B);
MATCL_LINALG_EXPORT mat_tup_4 balance_eig(Matrix&& A, const Matrix& B);
MATCL_LINALG_EXPORT mat_tup_4 balance_eig(const Matrix& A, Matrix&& B);
MATCL_LINALG_EXPORT mat_tup_4 balance_eig(Matrix&& A, const Matrix&& B);

//-------------------------------------------------------------------------
//                      LINEAR EQUATIONS
//-------------------------------------------------------------------------

// compute row and column scalings intended to equilibrate a hermitian positive 
// definite matrix A and reduce its condition number (with respect to the two-norm);
// scaling factor D are given by diagonal elements of A, chosen so that the scaled
// matrix AA has ones on the diagonal; this choice of D puts the condition number 
// of A within a factor N of the smallest possible condition number over all possible
// diagonal scalings
//     [AA, D] = balance_posdef2(A, pow2, tol_sing);
//     [D]     = balance_posdef(A, pow2, tol_sing);
// A       - a square matrix, hermitian, semi positive definite matrix
// pow2    - if true, then scaling factors are restricted to powers of 2
// tol_sing- torelance used to detect singular diagonal elements, for such
//           elements scaling is set to 1
// AA      - balanced matrix if required
//             AA  = diag(D) * A * diag(D)
// D       - scaling vector
// exception is thrown if A is detected not to be semi positive definite
MATCL_LINALG_EXPORT mat_tup_2   balance_posdef2(const Matrix& A, bool pow2 = true, Real tol_sing = 0.0);
MATCL_LINALG_EXPORT mat_tup_2   balance_posdef2(Matrix&& A, bool pow2 = true, Real tol_sing = 0.0);
MATCL_LINALG_EXPORT Matrix      balance_posdef(const Matrix& A, bool pow2 = true, Real tol_sing = 0.0);

// compute diagonal scalings intended to equilibrate a symmetric/hermitian
// matrix A and reduce its condition number; scaling factor D are computed
// using Bunch's algorithm, the scaled matrix
//     B = diag(D) * A * diag(D)
// has columns and rows with infinity norm 1 possibly except columns and rows
// j with abs(A(j,j) <= min_diag, in this case infinity norm of j-th row and
// column is less or equal 1
//     [B, D]  = balance_sym2(A, pow2, tol_sing);
//     [D]     = balance_sym(A, pow2, tol_sing);
// A       - a square symmetric/hermitian matrix
// pow2    - if true, then scaling factors are restricted to powers of 2
// tol_sing- torelance used to detect singular diagonal elements, for such
//           elements scaling is set to 1
// B       - balanced matrix if required
// D       - scaling vector
MATCL_LINALG_EXPORT mat_tup_2   balance_sym2(const Matrix& A, bool pow2 = true, Real tol_sing = 0.0);
MATCL_LINALG_EXPORT mat_tup_2   balance_sym2(Matrix&& A, bool pow2 = true, Real tol_sing = 0.0);
MATCL_LINALG_EXPORT Matrix      balance_sym(const Matrix& A, bool pow2 = true, Real tol_sing = 0.0);

// compute row and column scalings intended to equilibrate am MxN matrix A and
// reduce its condition number; return row and column scale fators R, C chosen to
// try to make the largest element in each row and column of the matrix 
//         B = diag(R) * A * diag(C)
// have absolute value 1 (i.e. vector infinity norm of all columns and rows is 1); 
// if given row or column contains only zero value, then scale factor is set to 1;
// use of these scaling factors is not guaranteed to reduce the condition number 
// of A but works well in practice; this is the algorithm used in Lapack's function
// geequ
//     [B, R, C]   = balance_gen2(A, pow2, tol_sing);
//     [R, C]      = balance_gen(A, pow2, tol_sing);
// A       - a MxN matrix
// pow2    - if true, then scaling factors are restricted to powers of 2
// tol_sing- torelance used to detect singular rows and columns, for these rows and
//           columns scaling is set to 1
// B       - balanced matrix if required
//             B  = diag(R) * A * diag(C)
// R       - scaling vector of size M x 1
// C       - scaling vector of size N x 1
MATCL_LINALG_EXPORT mat_tup_3   balance_gen2(const Matrix& A, bool pow2 = true, Real tol_sing = 0.0);
MATCL_LINALG_EXPORT mat_tup_3   balance_gen2(Matrix&& A, bool pow2 = true, Real tol_sing = 0.0);
MATCL_LINALG_EXPORT mat_tup_2   balance_gen(const Matrix& A, bool pow2 = true, Real tol_sing = 0.0);

// compute row and column scalings intended to equilibrate am MxN matrix A and
// reduce its condition number; return row and column scale fators R, C chosen to
// try to make vector infinity norm of all columns and rows of the matrix
//         B = diag(R) * A * diag(C)
// equal 1; when A is symmetric/hermitian, then B is also and R = C; if given
// row or column contains only zero values, then scale factor is set to 1;
//
//     [B, R, C]   = balance_inf2(A, pow2, tol_sing, opt);
//     [R, C]      = balance_inf(A, pow2, tol_sing, opt);
// A       - a MxN matrix
// pow2    - if true, then scaling factors are restricted to powers of 2
// tol_sing- torelance used to detect singular rows and columns, for these rows
//           and columns scaling is set to 1
// opts    - options controlling convergence check; see balancing options in
//           opt::linsolve namespace; iterations are finished if infinity norm 
//           of correction to scale factors R, C is less tol or iteration limit 
//           itmax is reached; default values: tol = 1e-3, itmax = 20.
// B       - balanced matrix if required
//             B  = diag(R) * A * diag(C)
// R       - scaling vector of size M x 1
// C       - scaling vector of size N x 1
//
// this is iterative algorithm with linear convergence with an asymptotic 
// rate of 1/2
// algorithm: Ruiz, "A symmetry preserving algorthm for matrix scaling",
// Technical Report 7552, INRIA, 2001
MATCL_LINALG_EXPORT mat_tup_3   balance_inf2(const Matrix& A, bool pow2 = true, Real tol_sing = 0.0,
                                    const options& opts = options());
MATCL_LINALG_EXPORT mat_tup_3   balance_inf2(Matrix&& A, bool pow2 = true, Real tol_sing = 0.0,
                                    const options& opts = options());
MATCL_LINALG_EXPORT mat_tup_2   balance_inf(const Matrix& A, bool pow2 = true, Real tol_sing = 0.0,
                                    const options& opts = options());

};