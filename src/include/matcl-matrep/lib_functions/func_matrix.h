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

#pragma once

#include "matcl-matrep/general/config.h"
#include "matcl-scalar/lib_functions/func_matrix.h"
#include "matcl-matrep/matrix/matrix.h"
#include "matcl-matrep/details/enablers.h"

namespace matcl
{

namespace md    = matcl::details;
namespace mrd   = matcl::raw::details;

// from symmetric product A * trans(A) if transpose is false
// or trans(A) * A otherwise
MATCL_MATMULT_EXPORT Matrix symprod(const Matrix& A, bool transpose = false);
MATCL_MATMULT_EXPORT Matrix symprod(Matrix&& A, bool transpose = false);

// from hermitian product A * ctrans(A) if transpose is false
// or ctrans(A) * A otherwise
MATCL_MATMULT_EXPORT Matrix herprod(const Matrix& A, bool transpose = false);
MATCL_MATMULT_EXPORT Matrix herprod(Matrix&& A, bool transpose = false);

// from symmetric sum A + trans(A)
MATCL_MATMULT_EXPORT Matrix symsum(const Matrix& A);
MATCL_MATMULT_EXPORT Matrix symsum(Matrix&& A);

// from hermitian sum A + ctrans(A)
MATCL_MATMULT_EXPORT Matrix hersum(const Matrix& A);
MATCL_MATMULT_EXPORT Matrix hersum(Matrix&& A);

// scale rows of a matrix A by scalars given in a vector D; this is
// equivalent to bdiag(D) * A; structure of sparse matrices is preserved
MATCL_MATMULT_EXPORT Matrix scale_rows(const Matrix& A, const Matrix& D);
MATCL_MATMULT_EXPORT Matrix scale_rows(Matrix&& A, const Matrix& D);

// scale columns of a matrix A by scalars given in a vector D; this is
// equivalent to A * bdiag(D); structure of sparse matrices is preserved
MATCL_MATMULT_EXPORT Matrix scale_cols(const Matrix& A, const Matrix& D);
MATCL_MATMULT_EXPORT Matrix scale_cols(Matrix&& A, const Matrix& D);

// scale rows and columns of a matrix A by scalars given in a vectors Dr
// and Dc respectively; this is equivalent to bdiag(Dr) * A * bdiag(Dc); 
// structure of sparse matrices is preserved
MATCL_MATMULT_EXPORT Matrix scale_rowscols(const Matrix& A, const Matrix& Dr, const Matrix& Dc);
MATCL_MATMULT_EXPORT Matrix scale_rowscols(Matrix&& A, const Matrix& Dr, const Matrix& Dc);

// matrix multiplication, perform op(A) x op(B), where A and B are matrices
// or scalars of conforming sizes and
// op(X) = X                               if trans_type is no_trans
// op(X) = X^T (transposition)             if trans_type is trans
// op(X) = X^C (conjuage transposition)    if trans_type is conj_trans
// transpositions are not performed explicitly, except cases A x B with
// storage: sparse x sparse if t_B != no_trans, then B is transposed
// band x sparse if t_B != no_trans and A is not diagonal, then B is 
// transposed
MATCL_MATMULT_EXPORT Matrix mmul(const Matrix& A, const Matrix& B, 
                                 trans_type t_A = trans_type::no_trans,
                                 trans_type t_B = trans_type::no_trans);
MATCL_MATMULT_EXPORT Matrix mmul(Matrix&& A, const Matrix& B,
                                 trans_type t_A = trans_type::no_trans,
                                 trans_type t_B = trans_type::no_trans);
MATCL_MATMULT_EXPORT Matrix mmul(const Matrix& A, Matrix&& B,
                                 trans_type t_A = trans_type::no_trans,
                                 trans_type t_B = trans_type::no_trans);
MATCL_MATMULT_EXPORT Matrix mmul(Matrix&& A, Matrix&& B,
                                 trans_type t_A = trans_type::no_trans,
                                 trans_type t_B = trans_type::no_trans);

// matrix multiplication perform op(A) x op(B), where A and B are matrices
// or scalars of conforming sizes and
// op(X) = X                               if trans_type is no_trans
// op(X) = X^T (transposition)             if trans_type is trans
// op(X) = X^C (conjuage transposition)    if trans_type is conj_trans
// op(X) = conj(X)                         if trans_type is conj
// if t_A or t_B is trans_type_ext::conj, then conjugate is taken on one of
// matrix argument or the returned matrix
MATCL_MATMULT_EXPORT Matrix mmul(const Matrix& A, const Matrix& B, 
                                 trans_type_ext t_A, trans_type_ext t_B);
MATCL_MATMULT_EXPORT Matrix mmul(Matrix&& A, const Matrix& B,
                                 trans_type_ext t_A, trans_type_ext t_B);
MATCL_MATMULT_EXPORT Matrix mmul(const Matrix& A, Matrix&& B,
                                 trans_type_ext t_A, trans_type_ext t_B);
MATCL_MATMULT_EXPORT Matrix mmul(Matrix&& A, Matrix&& B,
                                 trans_type_ext t_A, trans_type_ext t_B);

// matrix multiplication; equivalent to the mmul function
inline Matrix               operator*(const Matrix& A, const Matrix& B)    {return mmul(A,B); };
inline Matrix               operator*(Matrix&& A, const Matrix& B)         {return mmul(std::move(A),B); };
inline Matrix               operator*(const Matrix& A, Matrix&& B)         {return mmul(A,std::move(B)); };
inline Matrix               operator*(Matrix&& A, Matrix&& B)              {return mmul(std::move(A),std::move(B)); };

// evaluate C = alpha * trans(A, t_A) * trans(B, t_B) + beta * C
// alpha and beta must be scalars
//
// C must be a dense matrix, possibly not initialized, that stores elements
// of type v_C, where value code v_C unifies value codes of all other 
// matrix arguments (i.e. matrix C can be converted to matrix with mat_code
// of any of arguments without warning)
//
// at exist the matrix C is unique; if A or B are submatrices of C, then 
// matcl prefers to make unique A or B (since they are smaller) if this 
// helps (i.e. C, A, and B are the only references to stored data)
MATCL_MATMULT_EXPORT void   gemm(const Matrix& alpha, const Matrix& A, const Matrix& B, 
                                 trans_type t_A, trans_type t_B, const Matrix& beta, Matrix& C);

// evaluate C(sel) = alpha * trans(A, t_A) * trans(B, t_B) + beta * C(sel)
// alpha and beta must be scalars; 
//
// sel is submatrix selection that makes dense view:
//
//     C(sel)  := C(c1, c2),       c1, c2 are continuous colons
//     C(sel)  := C(c1, col),      c1 is a continuous colon and col is 
//                                 Integer
//     C(sel)  := C(row, c2),      row is Integer and c2 is a continuous
//                                 colon
//     C(sel)  := C(row, col),     row and col are Integers
//
// a continuous colon can be represented as: c = first : 1 : last for some
// first and last
//
// C must be a dense matrix, possibly not initialized, that stores elements
// of type v_C, where value code v_C unifies value codes of all other 
// matrix arguments (i.e. matrix C can be converted to matrix with mat_code
// of any of arguments without warnings)
//
// at exist the matrix C is unique; if A or B are submatrices of C, then 
// matcl prefers to make unique A or B (since they are smaller) if this 
// helps (i.e. C, A, and B are the only references to stored data)
MATCL_MATMULT_EXPORT void   gemm(const Matrix& alpha, const Matrix& A, const Matrix& B, 
                                 trans_type t_A, trans_type t_B, const Matrix& beta, sub_matrix_1&& C);
MATCL_MATMULT_EXPORT void   gemm(const Matrix& alpha, const Matrix& A, const Matrix& B, 
                                 trans_type t_A, trans_type t_B, const Matrix& beta, sub_matrix_2&& C);
MATCL_MATMULT_EXPORT void   gemm(const Matrix& alpha, const Matrix& A, const Matrix& B, 
                                 trans_type t_A, trans_type t_B, const Matrix& beta, sub_matrix&& C);

// evaluate gemm function when the matrix C is marked as effectively unique,
// call the gemm function defined for C as a matrix or as a submatrix 
// depending on argument used to create unique_matrix
MATCL_MATMULT_EXPORT void   gemm(const Matrix& alpha, const Matrix& A, const Matrix& B, 
                                 trans_type t_A, trans_type t_B, const Matrix& beta, unique_matrix&& C);

// form abs(op(A)) * abs(X) + abs(B), where op(A) = trans(A, t_A)
// X and B should be dense
MATCL_MATMULT_EXPORT Matrix mmul_abs(const Matrix& A, const Matrix& X, trans_type t_A, const Matrix& B);
MATCL_MATMULT_EXPORT Matrix mmul_abs(const Matrix& A, const Matrix& X, trans_type t_A, Matrix&& B);

// form a sequence of multiplications A1 x A2 x ... x Ak with optimal
// sequence of multiplications
Matrix                      chain_mult(const Matrix& A1, const Matrix& A2);
MATCL_MATMULT_EXPORT Matrix chain_mult(const Matrix& A1, const Matrix& A2, const Matrix& A3);
MATCL_MATMULT_EXPORT Matrix chain_mult(const Matrix& A1, const Matrix& A2, const Matrix& A3,
                                    const Matrix& A4);
MATCL_MATMULT_EXPORT Matrix chain_mult(const Matrix& A1, const Matrix& A2, const Matrix& A3,
                                    const Matrix& A4, const Matrix& A5);
MATCL_MATMULT_EXPORT Matrix chain_mult(const Matrix& A1, const Matrix& A2, const Matrix& A3,
                                       const Matrix& A4, const Matrix& A5, const Matrix& A6);

// form kronecker product of two matrices A, B, i.e a matrix
// kron(A,B) = [a_11 * B ... a_1n * B]
//             [  ...    ...   ...   ]  
//             [a_m1 * B ... a_mn * B]
MATCL_MATMULT_EXPORT Matrix kron(const Matrix& A, const Matrix& B);
MATCL_MATMULT_EXPORT Matrix kron(Matrix&& A, const Matrix& B);
MATCL_MATMULT_EXPORT Matrix kron(const Matrix& A, Matrix&& B);
MATCL_MATMULT_EXPORT Matrix kron(Matrix&& A, Matrix&& B);

};

#include "matcl-matrep/details/func_matrix.inl"