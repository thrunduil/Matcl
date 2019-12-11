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
#include "matcl-matrep/matrix/permvec.h"
#include "matcl-linalg/special_matrices/unitary_matrix.h"

namespace matcl
{

namespace md = matcl::details;

// return type of the function construct_householder
using ret_construct_householder = tuple<Matrix,Matrix,Matrix>;

// construct an elementary reflector H such that
//     H * [a] = [b], H' * H = I       (1)
//         [x]   [0]
// where a, b are scalars and x is a vector;
// H is represented in the Lapack's form:
//     H = I - tau * [1] * [1, v']
//                   [v]
// where tau is a scalar, v is a vector; if all elements of x are zero,
// then tau = 0, otherwise 1 <= tau <= 2 for real values and
// 1 <= real(tau) <= 2, abs(tau-1) <= 1 for complex values
//     [w, tau, b] = construct_householder(y)
// y   - a vector of size N x 1 (or 1 x N) representing the vector [a;x]; N >= 1
// w   - a vector of size N x 1 (or 1 x N) representing the vector [1;v]
// tau - real scalar
// b   - a scalar b as in (1)
//
// not available for sparse and band matrices
MATCL_LINALG_EXPORT ret_construct_householder construct_householder(const Matrix& y);
MATCL_LINALG_EXPORT ret_construct_householder construct_householder(Matrix&& y);

// construct unitary matrix from Householder reflectors; given K elementary
// reflectors H_1, ..., H_k, each of size N function creates an unitary
// matrix representing H_1 x ... x H_k
//     U = householder_to_unitary(T, tau)
// T   - lower triangular matrix of size N x K, where N >= 1, K <= N representing
//       elementary reflectors, i-th elementary reflector is represented by i-th
//       column of T and i-th element of tau; it is assumed that i-th column of
//       T has the form [Z; 1; v], where Z is a vector of zeroes of size (i-1)x1
//       and v is a vector representing elementary reflector; upper triangular part
//       of T is not referred and may contain any values
// tau - vector of length K; 
// cols- optional argument, number of columns of unitary matrix; 1 <= cols <= N;
//       on default cols = N
// U   - unitary matrix representing H_1 x ... x H_k of size N x N
//
// elementary reflectors must have form as returned by construct_householder
//
// not available for sparse matrices
MATCL_LINALG_EXPORT unitary_matrix householder_to_unitary(const Matrix& T, const Matrix& tau);
MATCL_LINALG_EXPORT unitary_matrix householder_to_unitary(Matrix&& T, const Matrix& tau);
MATCL_LINALG_EXPORT unitary_matrix householder_to_unitary(const Matrix& T, const Matrix& tau,
                                        Integer cols);
MATCL_LINALG_EXPORT unitary_matrix householder_to_unitary(Matrix&& T, const Matrix& tau, 
                                        Integer cols);

// create compact WY representations of a product of Householder reflectors; 
// given K elementary reflectors H_1, ..., H_k, each of size N function creates 
// matrixes V, T such that
//         H_1 x ... x H_k = I - V x Y x V'            (1)
//
//     Y = householder_WY(T, tau)
// T   - unit lower triangular matrix of size N x K, where N >= 1, K <= N representing
//       elementary reflectors, i-th elementary reflector is represented by i-th
//       column of T and i-th element of tau; it is assumed that i-th column of
//       T has the form [Z; 1; v], where Z is a vector of zeroes of size (i-1)x1
//       and v is a vector representing elementary reflector; upper triangular part
//       of T is not referred and may contain any values
// tau - vector of length K; 
// Y   - an upper triangular matrix; the matrix V in (1) is the same as the matrix T
//
// elementary reflectors must have form as returned by construct_householder
// not available for sparse and matrices
MATCL_LINALG_EXPORT Matrix householder_WY(const Matrix& T, const Matrix& tau);
MATCL_LINALG_EXPORT Matrix householder_WY(Matrix&& T, const Matrix& tau);

};
