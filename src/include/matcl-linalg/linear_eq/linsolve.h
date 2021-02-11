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
#include "matcl-matrep/matrix/permvec.h"
#include "matcl-linalg/general/config_linalg.h"
#include "matcl-linalg/options/options_linsolve.h"
#include "matcl-linalg/linear_eq/linsolve_object.h"
#include "matcl-linalg/special_matrices/unitary_matrix.h"
#include "matcl-matrep/details/enablers.h"

namespace matcl
{

namespace md = matcl::details;

//--------------------------------------------------------------------
//                  construction of linsolve object
//--------------------------------------------------------------------
/// construct linsolve_obj for a square matrix
///
///     lo  = make_linsolve_obj(A, opts)
///
/// A       - a square matrix;
///
///           factorizations will be performed only if necessary; 
///
///           Factorization will be selected according to storage and structures
///           assigned to the matrix A; considered structures:
///             all structures from predefined_struct_type,
///             unitary, posdef, semi_posdef
///           possible factorizations: lu, Cholesky, ldl
///           Cholesky factorization is always used if the matrix A has posdef 
///           or semi_posdef structure
/// opts    - options used by factorization routines, see namespace 
///           opt::linsolve for details; opts also controls factorization used:
///             use_rr      -   if true, then rank revealing factorizations will
///                             be used
///             use_ir      -   controls iterative refinement during solution phase
///             test_sol    -   controls correctness checks of computed solution
///             can_use_ldl -   if true then LDL is used for indefinite symmetric
///                             matrices
///             do_balancing-   if true and rank_revealing = false, then balancing 
///                             is performed
///             do_balancing_rr-if true and rank_revealing = true, then balancing 
///                             is performed
/// lo      - a linsolve_obj
MATCL_LINALG_EXPORT linsolve_obj make_linsolve_obj(const Matrix& A, const matcl::options& opts = matcl::options());
MATCL_LINALG_EXPORT linsolve_obj make_linsolve_obj(Matrix&& A, const matcl::options& opts = matcl::options());

/// construct linsolve_obj for a matrix A = A1 * A2 or A = A1 * A2 * A3, 
/// where A1, A2, A3 are given by linsolve_obj, when the matrix A is required
/// then A will be recreated from A1-A3; when A1-A3 are factors of some decomposition
/// then the factored matrix A should be provided, otherwise calculated backward and
/// forward errors will be underestimated
///
///     lo  = linsolve_seq(A1, A2)
///     lo  = linsolve_seq(A1, A2, A3)
///
/// A1-A3   - linsolve objects representing matrices A1, A2, A3
/// lo      - a linsolve_obj representing a matrix A1 * A2 or A1 * A2 * A3
/// 
/// if A is detected to be singular, then exception is thrown
MATCL_LINALG_EXPORT linsolve_obj linsolve_seq(const linsolve_obj& A1, const linsolve_obj& A2);
MATCL_LINALG_EXPORT linsolve_obj linsolve_seq(const linsolve_obj& A1, const linsolve_obj& A2, 
                                                 const linsolve_obj& A3);

/// construct linsolve_obj for a matrix A = A1 * A2 or A = A1 * A2 * A3, 
/// where A1, A2, A3 are given by linsolve_obj, 
///
///     lo  = linsolve_seq(A, A1, A2)
///     lo  = linsolve_seq(A, A1, A2, A3)
///
/// A       - the product of A1-A3
/// A1-A3   - linsolve objects representing matrices A1, A2, A3
/// lo      - a linsolve_obj representing a matrix A1 * A2 or A1 * A2 * A3
/// 
/// if A is detected to be singular, then exception is thrown
MATCL_LINALG_EXPORT linsolve_obj linsolve_seq(const Matrix& A, const linsolve_obj& A1, const linsolve_obj& A2);
MATCL_LINALG_EXPORT linsolve_obj linsolve_seq(const Matrix& A, const linsolve_obj& A1, const linsolve_obj& A2, 
                                                 const linsolve_obj& A3);

/// construct linsolve_obj representing a matrix with NaN values
///     lo  = linsolve_nan(N, vc)
/// N       - number of rows and columns of the matrix
/// vc      - value code of stored elements
/// lo      - a linsolve_obj
MATCL_LINALG_EXPORT linsolve_obj linsolve_nan(Integer N, value_code vc);

/// construct linsolve_obj for a unitary matrix
///     lo  = linsolve_unitary(U, opts)
/// U       - a square unitary matrix
/// lo      - a linsolve_obj
MATCL_LINALG_EXPORT linsolve_obj linsolve_unitary(const unitary_matrix& U);
MATCL_LINALG_EXPORT linsolve_obj linsolve_unitary(const Matrix& U);

/// construct linsolve_obj for a diagonal matrix
///     lo  = linsolve_diag(D, opts)
/// D       - a square diagonal matrix
/// opts    - options, see namespace opt::linsolve for details
/// lo      - a linsolve_obj
MATCL_LINALG_EXPORT linsolve_obj linsolve_diag(const Matrix& D, const options& opts = options());

/// construct linsolve_obj for a block diagonal matrix with 1x1 or 2x2 blocks
///     lo  = linsolve_diag_22(D, opts)
/// D       - a square block diagonal matrix with blocks of size 1x1 or 2x2
/// opts    - options, see namespace opt::linsolve for details
/// lo      - a linsolve_obj
MATCL_LINALG_EXPORT linsolve_obj linsolve_diag_22(const Matrix& D, const options& opts = options());

/// construct linsolve_obj for an upper triangular or lower triangular matrix
///     lo  = linsolve_triang(T, opts)
/// T       - a square lower or upper triangular matrix
/// opts    - options, see namespace opt::linsolve for details
/// lo      - a linsolve_obj
MATCL_LINALG_EXPORT linsolve_obj linsolve_triang(const Matrix& T, const options& opts = options());

/// construct linsolve_obj for a permuted upper triangular or lower 
/// triangular matrix T(p^-1, q^-1)
///     lo  = linsolve_triang(T, p, q, opts)
/// T       - a square lower or upper triangular matrix
/// p,q     - permutation vectors
/// opts    - options, see namespace opt::linsolve for details
/// lo      - a linsolve_obj
MATCL_LINALG_EXPORT linsolve_obj linsolve_triang(const Matrix& T, const permvec& p, 
                                        const permvec& q, const options& opts = options());

/// construct linsolve_obj for an upper Hessenberg matrix
///     lo  = linsolve_uhess(T, opts)
/// T       - a square upper Hessenberg matrix
/// opts    - options, see namespace opt::linsolve for details
/// lo      - a linsolve_obj
MATCL_LINALG_EXPORT linsolve_obj linsolve_uhess(const Matrix& T, const options& opts = options());

/// construct linsolve_obj for a lower Hessenberg matrix
///     lo  = linsolve_lhess(T, opts)
/// T       - a square lower Hessenberg matrix
/// opts    - options, see namespace opt::linsolve for details
/// lo      - a linsolve_obj
///
/// not available for band and sparse matrices
MATCL_LINALG_EXPORT linsolve_obj linsolve_lhess(const Matrix& T, const options& opts = options());

/// construct linsolve_obj for a matrix A balanced with diagonal
/// scaling given by vectors Dl, Dr, where balanced matrix B satisfies
///     B = diag(Dl) * A * diag(Dr)
///     lo  = linsolve_balanced(Dl, B, Dr)
///     lo  = linsolve_balanced(A, Dl, B, Dr)
/// A       - optional, original; if A is not provided and is required by 
///           any of function defined linsolve_obj, then A will be recreated
///           from Dl, B, Dr
/// Dl      - the left scale factors
/// B       - linsolve_obj for the balanced matrix B
/// Dr      - the right scale factors
/// lo      - a linsolve_obj
MATCL_LINALG_EXPORT linsolve_obj linsolve_balanced(const Matrix& Dl, const linsolve_obj& B,
                                                   const Matrix& Dr);
MATCL_LINALG_EXPORT linsolve_obj linsolve_balanced(const Matrix& A, const Matrix& Dl, const linsolve_obj& B,
                                                   const Matrix& Dr);

/// construct linsolve_obj for a matrix A balanced with symmetric diagonal
/// scaling given by vector Dlr, where balanced matrix B satisfies
///     B = diag(Dlr) * A * diag(Drr)
///     lo  = linsolve_balanced_sym(Dlr, B)
///     lo  = linsolve_balanced_sym(A, Dlr, B)
/// A       - optional, original matrix; if A is not provided and 
///           is required by any of function defined linsolve_obj, then
///           A will be recreated from Dlr, B
/// Dlr     - the left and right scale factors
/// B       - linsolve_obj for the balanced matrix B
/// lo      - a linsolve_obj
MATCL_LINALG_EXPORT linsolve_obj linsolve_balanced_sym(const Matrix& Dlr, const linsolve_obj& B);
MATCL_LINALG_EXPORT linsolve_obj linsolve_balanced_sym(const Matrix& A, const Matrix& Dlr, const linsolve_obj& B);

/// construct linsolve_obj for a matrix A from linsolve object for permuted
/// problem
///     lo  = linsolve_perm(A, Apq, p, q)
///     lo  = linsolve_perm(Apq, p, q)
/// A       - optional, original matrix; if A is not provided and 
///           is required by any of function defined linsolve_obj, then
///           A will be recreated from p, q, and Apq
/// Apq     - linsolve object for the matrix A(p,q)
/// p, q    - row and column permutation vectors
/// lo      - a linsolve_obj
MATCL_LINALG_EXPORT linsolve_obj linsolve_perm(const Matrix& A, const linsolve_obj& Apq, const permvec& p,
                                                   const permvec& q);
MATCL_LINALG_EXPORT linsolve_obj linsolve_perm(const linsolve_obj& Apq, const permvec& p, const permvec& q);

/// construct linsolve_obj for a matrix A from linsolve object for symmetrically
/// permuted problem
///     lo  = linsolve_symperm(A, App, p)
///     lo  = linsolve_symperm(App, p)
/// A       - optional, original matrix; if A is not provided and 
///           is required by any of function defined linsolve_obj, then
///           A will be recreated from p and App
/// App     - linsolve object for the matrix A(p,p)
/// p       - row and column permutation vector
/// lo      - a linsolve_obj
MATCL_LINALG_EXPORT linsolve_obj linsolve_symperm(const Matrix& A, const linsolve_obj& Apq, const permvec& p);
MATCL_LINALG_EXPORT linsolve_obj linsolve_symperm(const linsolve_obj& Apq, const permvec& p);

/// construct linsolve_obj that perform additional steps in solution phase
///     lo  = linsolve_extended_sol(lo, opts)
/// lo      - a linsolve object
/// opts    - options controlling solution phase; see namespace opt::linsolve 
///           for details, section: options controling solution phase in linsolve
///           objects
/// lo      - a linsolve_obj which performs additional steps in solution phase
MATCL_LINALG_EXPORT linsolve_obj linsolve_extended_sol(const linsolve_obj& lo, const options& opts);

//--------------------------------------------------------------------
//         solve linear equations - low level functions
//--------------------------------------------------------------------

/// solve linear equation with upper or lower triangular matrices
/// or with general matrices using LU decomposition

/// solve linear equation:
///     A X = b             if trans = trans_type::no_trans
///     trans(A) X = b      if trans = trans_type::trans
///     ctrans(A) X = b     if trans = trans_type::conj_trans
/// opts is used by LU decomposition performed internally if necessary,
/// see namespace opt::linsolve for details
MATCL_LINALG_EXPORT Matrix linsolve(const matcl::Matrix& A,const matcl::Matrix& b, 
                                    trans_type trans = trans_type::no_trans,
                                    const matcl::options& opts = matcl::options());
MATCL_LINALG_EXPORT Matrix linsolve(const matcl::Matrix& A, matcl::Matrix&& b, 
                                    trans_type trans = trans_type::no_trans,
                                    const matcl::options& opts = matcl::options());
MATCL_LINALG_EXPORT Matrix linsolve(matcl::Matrix&& A,const matcl::Matrix& b, 
                                    trans_type trans = trans_type::no_trans,
                                    const matcl::options& opts = matcl::options());
MATCL_LINALG_EXPORT Matrix linsolve(matcl::Matrix&& A, matcl::Matrix&& b, 
                                    trans_type trans = trans_type::no_trans,
                                    const matcl::options& opts = matcl::options());

/// solve linear equation:
///      A(p^-1,q^-1)          X = b   if trans = trans_type::no_trans
///      trans(A(p^-1,q^-1))   X = b   if trans = trans_type::trans
///      ctrans(A(p^-1,q^-1))  X = b   if trans = trans_type::conj_trans
///
/// where p, q are permutation matrices. This is equvaluent to
///      X = (op(A)^-1 * b(p,:))(q^-1,:);
/// where p^-1 = invperm(p); 
/// opts is used by LU decomposition performed internally if necessary,
/// see namespace opt::linsolve for details
MATCL_LINALG_EXPORT Matrix linsolve(const matcl::Matrix& A, const matcl::permvec& p, 
                                    const matcl::permvec& q, const matcl::Matrix& b, 
                                    trans_type trans = trans_type::no_trans,
                                    const matcl::options& opts = matcl::options());
MATCL_LINALG_EXPORT Matrix linsolve(const matcl::Matrix& A, const matcl::permvec& p, 
                                    const matcl::permvec& q, matcl::Matrix&& b, 
                                    trans_type trans = trans_type::no_trans,
                                    const matcl::options& opts = matcl::options());
MATCL_LINALG_EXPORT Matrix linsolve(matcl::Matrix&& A, const matcl::permvec& p, 
                                    const matcl::permvec& q, const matcl::Matrix& b, 
                                    trans_type trans = trans_type::no_trans,
                                    const matcl::options& opts = matcl::options());
MATCL_LINALG_EXPORT Matrix linsolve(matcl::Matrix&& A, const matcl::permvec& p, 
                                    const matcl::permvec& q, matcl::Matrix&& b, 
                                    trans_type trans = trans_type::no_trans,
                                    const matcl::options& opts = matcl::options());

/// solve linear equation:
///     X A = b
/// transposed versions are not available for reverted linsolve
/// opts is used by LU decomposition performed internally if necessary,
/// see namespace opt::linsolve for details
MATCL_LINALG_EXPORT Matrix linsolve_rev(const matcl::Matrix& A, const matcl::Matrix& b,
                                        const matcl::options& opts = matcl::options());
MATCL_LINALG_EXPORT Matrix linsolve_rev(const matcl::Matrix& A, matcl::Matrix&& b,
                                        const matcl::options& opts = matcl::options());
MATCL_LINALG_EXPORT Matrix linsolve_rev(matcl::Matrix&& A, const matcl::Matrix& b,
                                        const matcl::options& opts = matcl::options());
MATCL_LINALG_EXPORT Matrix linsolve_rev(matcl::Matrix&& A, matcl::Matrix&& b,
                                        const matcl::options& opts = matcl::options());

/// solve linear equation:
///         X * A(p^-1,q^-1) = b
///
/// where p, q are permutation matrices. This is equvaluent to
///     X = (b(:,q) * A^-1)(:,p^-1)
/// where p^-1 = invperm(p)
/// transposed versions are not available for reverted linsolve
/// opts is used by LU decomposition performed internally if necessary,
/// see namespace opt::linsolve for details
MATCL_LINALG_EXPORT Matrix linsolve_rev(const matcl::Matrix& A, const matcl::permvec& p, 
                                        const matcl::permvec& q, const matcl::Matrix& b,
                                        const matcl::options& opts = matcl::options());
MATCL_LINALG_EXPORT Matrix linsolve_rev(const matcl::Matrix& A, const matcl::permvec& p, 
                                        const matcl::permvec& q, matcl::Matrix&& b,
                                        const matcl::options& opts = matcl::options());
MATCL_LINALG_EXPORT Matrix linsolve_rev(matcl::Matrix&& A, const matcl::permvec& p, 
                                        const matcl::permvec& q, const matcl::Matrix& b,
                                        const matcl::options& opts = matcl::options());
MATCL_LINALG_EXPORT Matrix linsolve_rev(matcl::Matrix&& A, const matcl::permvec& p, 
                                        const matcl::permvec& q, matcl::Matrix&& b,
                                        const matcl::options& opts = matcl::options());

/// solve linear equation:
///          X * A         = b   if trans = trans_type::no_trans
///          X * trans(A)  = b   if trans = trans_type::trans
///          X * ctrans(A) = b   if trans = trans_type::conj_trans
///
/// if trans != trans_type::no_trans and A and b are sparse matrices, then A or b is transposed
/// opts is used by LU decomposition performed internally if necessary,
/// see namespace opt::linsolve for details
MATCL_LINALG_EXPORT Matrix linsolve_rev2(const matcl::Matrix& A, const matcl::Matrix& b, 
                                         trans_type trans = trans_type::no_trans,
                                         const matcl::options& opts = matcl::options());
MATCL_LINALG_EXPORT Matrix linsolve_rev2(const matcl::Matrix& A, matcl::Matrix&& b, 
                                         trans_type trans = trans_type::no_trans,
                                         const matcl::options& opts = matcl::options());
MATCL_LINALG_EXPORT Matrix linsolve_rev2(matcl::Matrix&& A, const matcl::Matrix& b, 
                                         trans_type trans = trans_type::no_trans,
                                         const matcl::options& opts = matcl::options());
MATCL_LINALG_EXPORT Matrix linsolve_rev2(matcl::Matrix&& A, matcl::Matrix&& b, 
                                         trans_type trans = trans_type::no_trans,
                                         const matcl::options& opts = matcl::options());

/// solve linear equation:
///         X * A(p^-1,q^-1)         = b   if trans = trans_type::no_trans
///         X * trans(A(p^-1,q^-1))  = b   if trans = trans_type::trans
///         X * ctrans(A(p^-1,q^-1)) = b   if trans = trans_type::conj_trans
///
/// where p, q are permutation matrices. This is equvaluent to
///     X = (b(:,q) * op(A)^-1)(:,p^-1)
/// where p^-1 = invperm(p)
///
/// if trans != trans_type::no_trans and A and b are sparse matrices, then A or b is transposed
/// opts is used by LU decomposition performed internally if necessary,
/// see namespace opt::linsolve for details
MATCL_LINALG_EXPORT Matrix linsolve_rev2(const matcl::Matrix& A, const matcl::permvec& p, 
                                         const matcl::permvec& q, const matcl::Matrix& b, 
                                         trans_type trans = trans_type::no_trans,
                                         const matcl::options& opts = matcl::options());
MATCL_LINALG_EXPORT Matrix linsolve_rev2(const matcl::Matrix& A, const matcl::permvec& p, 
                                         const matcl::permvec& q, matcl::Matrix&& b, 
                                         trans_type trans = trans_type::no_trans,
                                         const matcl::options& opts = matcl::options());
MATCL_LINALG_EXPORT Matrix linsolve_rev2(matcl::Matrix&& A, const matcl::permvec& p, 
                                         const matcl::permvec& q, const matcl::Matrix& b, 
                                         trans_type trans = trans_type::no_trans,
                                         const matcl::options& opts = matcl::options());
MATCL_LINALG_EXPORT Matrix linsolve_rev2(matcl::Matrix&& A, const matcl::permvec& p, 
                                         const matcl::permvec& q, matcl::Matrix&& b, 
                                         trans_type trans = trans_type::no_trans,
                                         const matcl::options& opts = matcl::options());

/// inverse of a square matrix A; opts is used by factorization routines 
/// performed internally if necessary; this function calls make_linsolve_object
MATCL_LINALG_EXPORT Matrix inv(const matcl::Matrix& A, const matcl::options& opts = matcl::options());
MATCL_LINALG_EXPORT Matrix inv(matcl::Matrix&& A, const matcl::options& opts = matcl::options());

/// add small perturbation to diagonal elements of a square triangular
/// matrix tol; if |d_ii| < tol_d, then i-th diagonal element is modified 
/// to sign(d_ii) * tol_d; if tol = 0 then no correction is applied; 
/// if tol < 0 then tol_d = -tol; if tol > 0, then tol_d = eps^tol * 
/// max(|T|), where eps is epsilon for given value type
///     [Td, sing, modif] = make_nonsingular(T, tol)
/// T       - square triangular matrix
/// tol     - control tolerance tol_d
/// Td      - modified matrix
/// sing    - if true then Td has exact zero on main diagonal
/// modif   - if true then Td != T
MATCL_LINALG_EXPORT tuple<Matrix,bool,bool> make_nonsingular(const Matrix& T, Real tol);
MATCL_LINALG_EXPORT tuple<Matrix,bool,bool> make_nonsingular(Matrix&& T, Real tol);

//--------------------------------------------------------------------
//                  special problems
//--------------------------------------------------------------------

/// inverse of a square matrix A of size 2 x 2:
///   [  A_11   A_12  ]
///   [  A_21   A_22  ]
/// using SVD decomposition
///
///     inv_22(A_11, A_12, A_21, A_22, sig_max, sig_min)
///
/// A_11, A_12, A_21, A_22  - on input the 2x2 matrix
///                         - on output the 2x2 inverse matrix
/// sig_max                 - maximum singular value of A
/// sig_min                 - minimum singular value of A
///
/// if the matrix A is singular, i.e. sig_min == 0, then original matrix
/// is returned and no exception is thrown
template<class V, class Enable = typename md::enable_if_float<V,void>::type,
         class VR = typename md::real_type<V>::type>
MATCL_LINALG_EXPORT void inv_22(V& A_11, V& A_12, V& A_21, V& A_22, VR& sig_max, VR& sig_min);

};