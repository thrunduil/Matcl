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
#include "matcl-linalg/special_matrices/matrix_functors.h"

namespace matcl
{

enum class basic_vector_norm
{
    norm_1,     // the first norm; sum |x_i|
    norm_2,     // the second norm; sqrt(sum (x_i)^2)
    norm_inf,   // maximum norm, max(|x|_1, ..., |x|_k)
};

// compute matrix norm:
//     E = norm(A,p)
// A - a matrix
// p - norm type; can take values 1, 2, Inf, -1 (Frobenius norm), or -2 (max(abs(A(i,j)),
//     note that this is not a norm); note that computation of the second norm of a 
//     matrix requires forming SVD decomposition and is very costly, if A is a vector, then
//     SVD is not computed
// E - norm
MATCL_LINALG_EXPORT Real norm(const Matrix& A, Real p = 1.);
MATCL_LINALG_EXPORT Real norm(const Matrix& A, basic_vector_norm p);

// compute vector norm for every row or column of the matrix A
//     E = norm_vec(A, p, d)
// A - a matrix of size M x N
// p - norm type
// d - compute norm along dimension d, if d = 1, then compute vector
//     norm for every column of A, otherwise compute vector norm for
//     every row of A
// E - computed norms; if d = 1, then E has size 1xN, otherwise E has
//     size Mx1
MATCL_LINALG_EXPORT Matrix norm_vec(const Matrix& A, basic_vector_norm p, Integer d = 1);

// compute vector norm for every elements of the matrix A; equivalent to
// norm_vec(vec(A), p))
//     E = norm_vec_all(A, p)
// A - a matrix of size M x N
// p - norm type
// E - computed norm
MATCL_LINALG_EXPORT Real norm_vec_all(const Matrix& A, basic_vector_norm p);

// estimate norm 1 of a matrix represented by the functor f; only computation of
// matrix times vector is required; uses Hager’s [1] algorithm as in Lapack;
// estimated norm g satisfies
//     g <= ||A||_1
//
// [1] William W. Hager, "Condition Estimates," SIAM J. Sci. Stat. Comput. 5, 1984,
//     311-316, 1984. 
MATCL_LINALG_EXPORT Real normest_1(const linear_operator& f);

// estimate infinity norm of a matrix represented by the functor f; equivalent to
// normest_1(ft), where ft represents transposed matrix
MATCL_LINALG_EXPORT Real normest_inf(const linear_operator& f);

// estimate norm 2 of a matrix represented by the functor f; only computation of
// matrix times vector is required; compute maximum singular value of f using
// pschur function
MATCL_LINALG_EXPORT Real normest_2(const linear_operator& f);

// return vector infinity norm | abs(A) * abs(X) |_inf, where X is a N x 1 vector and 
// A is M x N matrix given by linear operator f
MATCL_LINALG_EXPORT Real abs_normest_r_inf(const linear_operator& f, const Matrix& X);

// return vector infinity norm | abs(X) * abs(A) |_inf, where X is a 1xM vector and 
// A is M x N matrix given by linear operator f
MATCL_LINALG_EXPORT Real abs_normest_l_inf(const linear_operator& f, const Matrix& X);

// return vector infinity norm d_i = | abs(A) * abs(X(:,i)) |_inf, for each column of X,
// where X is a N x K matrix and A is M x N matrix given by linear operator f
// return 1xK matrix with elements d_i
MATCL_LINALG_EXPORT Matrix abs_normest_r_vec_inf(const linear_operator& f, const Matrix& X);

// return vector infinity norm d_i = | abs(X(:,i)) * abs(A) |_inf, for each row of X,
// where X is a K x M matrix and A is M x N matrix given by linear operator f
// return Kx1 matrix with elements d_i
MATCL_LINALG_EXPORT Matrix abs_normest_l_vec_inf(const linear_operator& f, const Matrix& X);

};