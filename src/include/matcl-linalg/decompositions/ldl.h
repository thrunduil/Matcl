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
#include "matcl-linalg/linear_eq/linsolve_object.h"

namespace matcl
{

///return type of ldl functions
using ldl_return_type = tuple<Matrix,Matrix,permvec>;

/// Bunch-Kaufman LDL' factorization for symmetric indefinite matrices
///     [L,D,p] = ldl(A, upper);
/// such that
///     A(p,p) = L * D * trans(L)
/// A       - symmetric matrix
/// upper   - true or false
/// L       - unit lower-trianglar (if upper = false), or unit upper-triangular
///           (if upper = true) matrix
/// D       - block diagonal matrix with 1x1 or 2x2 blocks
/// p       - permutation vector
///
/// not available for sparse and band matrices
MATCL_LINALG_EXPORT ldl_return_type ldl(const Matrix& A, bool upper);
MATCL_LINALG_EXPORT ldl_return_type ldl(Matrix&& A, bool upper);

/// Bunch-Kaufman LDL' factorization for hermitian indefinite matrices
///     [L,D,p] = ldl_herm(A, upper);
/// such that 
///     A(p,p) = L * D * ctrans(L)
/// A       - hermitian matrix
/// upper   - true or false
/// L       - unit lower-trianglar (if upper = false), or upper-triangular
///           (if upper = true) matrix
/// D       - block diagonal matrix with 1x1 or 2x2 blocks
/// p       - permutation vector
///
/// not available for sparse and band matrices
MATCL_LINALG_EXPORT ldl_return_type ldl_herm(const Matrix& A, bool upper);
MATCL_LINALG_EXPORT ldl_return_type ldl_herm(Matrix&& A, bool upper);

/// construct linsolve_obj from LDL factors such that
///     A(p,p)  = L * D * trans(L)  or  (1)
///     A(p,p)  = L * D * ctrans(L)     (2)
///
///     lo  = linsolve_ldl(A, L, D, p, upper, sym, opts)
/// A       - the factored matrix
/// L       - unit lower triangular or unit upper triangular
/// D       - block diagonal matrix with 1x1 or 2x2 blocks
/// p       - permutation vector
/// sym     = true:     factorization of type (1)
///         = false:    factorization of type (2)
/// opts    - options, see namespace opt::linsolve for details
/// lo      - a linsolve_obj
/// 
/// if D is detected to be singular, then exception is thrown
MATCL_LINALG_EXPORT linsolve_obj linsolve_ldl(const Matrix& A, const Matrix& L, const Matrix& D, 
                                        const permvec& p, bool sym, const options& opts = options());

/// construct linsolve_obj for a matrix A using LDL decomposition if required
///     lo  = linsolve_ldl(A, sym, opts)
/// A       - a square matrix, complex symmetric matrix (if sym = true) or hermitian/
///           real symmetric matrix (if sym = false)
/// sym     = false:   A is real symmetric or complex hermitian
///         = true:    A is complex symmetric
/// opts    - options controlling balancing; balancing is performed if option 
///           do_balancing is true using balance_sym function
/// lo      - a linsolve_obj
/// 
/// if A is detected to be singular, then exception is thrown
/// not available for sparse and band matrices
MATCL_LINALG_EXPORT linsolve_obj linsolve_ldl(const Matrix& A, bool sym, const options& opts = options());
MATCL_LINALG_EXPORT linsolve_obj linsolve_ldl(Matrix&& A, bool sym, const options& opts = options());

};