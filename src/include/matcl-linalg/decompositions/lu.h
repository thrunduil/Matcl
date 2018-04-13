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
#include "matcl-core/options/matcl_options.h"
#include "matcl-linalg/options/options_linsolve.h"
#include "matcl-linalg/linear_eq/linsolve_object.h"

namespace matcl
{

/// return type of lu function
using lu_return_type    = tuple<Matrix,Matrix,permvec,permvec>;

/// perform LU decomposition of the matrix A in the form
///
///          A(p,q) = L * U
///
/// where L is lower triangular matrix, U is upper triangular,
/// and p, q are permutations
///     [L,U,p,q] = lu(A, opts);
/// A       - M x N matrix
/// opts    - options, see namespace opt::linsolve for details
/// L       - M x K lower triangular matrix, K <= min(M,N); using default options
///           L is not necessary unit lower triangular
/// U       - K x N upper triangular matrix, K <= min(M,N). If pivot type is
///           pivot_type::rook or pivot_type::complete, then K is estimated rank
///           of A; pivoting strategy can be set by setting options opt.
/// p,q     - permutations
MATCL_LINALG_EXPORT lu_return_type  lu(const Matrix& A, const matcl::options& opts = matcl::options());
MATCL_LINALG_EXPORT lu_return_type  lu(Matrix&& A, const matcl::options& opts = matcl::options());

/// perform incomplete LU decomposition of the matrix A in the form
///
///          A(p,q) ~ L * U
///
/// where L is lower triangular matrix, U is upper triangular,
/// and p, q are permutations
///     [L,U,p,q] = lu(A, opts);
/// A       - M x N, M >= N, if A is dense of banded then standard lu factorization
///           is performed; if A has positive_definite of semi_positive_definite flag,
///           then diagonal privoting is preferred
/// opts    - options, see namespace opt::linsolve for details
/// L       - M x N lower triangular matrix; using default options L is not 
///           necessary unit lower triangular
/// U       - N x N upper triangular matrix
/// p,q     - permutations
/// this function used ilu algorithm from superlu
MATCL_LINALG_EXPORT lu_return_type ilu(const Matrix& A, const matcl::options& opts = matcl::options());
MATCL_LINALG_EXPORT lu_return_type ilu(Matrix&& A, const matcl::options& opts = matcl::options());

/// construct linsolve_obj from LU factors of a matrix A
///     lo  = linsolve_lu(A, L, U, p, q, opts)
/// A       - the factored matrix
/// L       - lower triangular matrix
/// U       - upper triangular matrix
/// p,q     - permutation vectors
/// opts    - options, see namespace opt::linsolve for details
/// lo      - a linsolve_obj
/// 
/// if L, U are not square or L,U are detected to be singular, then exception is thrown
MATCL_LINALG_EXPORT linsolve_obj linsolve_lu(const Matrix& A, const Matrix& L, const Matrix& U, 
                                    const permvec& p, const permvec& q, const options& opts = options());

/// construct linsolve_obj for a matrix A using lu decomposition if required
///     lo  = linsolve_lu(A, opts)
/// A       - a square matrix
/// opts    - options used by LU decomposition, see namespace opt::linsolve section 
///           'options for lu decomposition' for details; additionally balancing is
///           performed if option do_balancing (for partial pivoting) or do_balancing_rr
///           (for rook and complete pivoting) is true using balance_gen function
/// lo      - a linsolve_obj
/// 
/// if A is detected to be singular, then exception is thrown
MATCL_LINALG_EXPORT linsolve_obj linsolve_lu(const Matrix& A, const matcl::options& opts = matcl::options());
MATCL_LINALG_EXPORT linsolve_obj linsolve_lu(Matrix&& A, const matcl::options& opts = matcl::options());

/// construct linsolve_obj for a tridiagonal matrix A using LU decomposition
///     lo  = linsolve_tridiag(DL, D, DU, opts)
/// DL      - vector of size (N-1)x1 of elements on subdiagonal of A
/// D       - vector of size Nx1 of elements on main diagonal of A
/// DU      - vector of size (N-1)x1 of elements on superdiagonal of A
/// opts    - options, see namespace opt::linsolve for details
/// lo      - a linsolve_obj
/// 
/// if A is detected to be singular, then exception is thrown
MATCL_LINALG_EXPORT linsolve_obj linsolve_tridiag(const Matrix& DL, const Matrix& D, const Matrix& DU,
                                                  const options& opts = options());

};