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

#include "matcl-matrep/matcl_matrep.h"
#include "matcl-linalg/general/config_linalg.h"

namespace matcl
{

//-----------------------------------------------------------------------
//                      GRAPH MANIPULATION
//-----------------------------------------------------------------------

// form symmetric graph based on A + A'
MATCL_LINALG_EXPORT Matrix symgraph_sum(const Matrix& A);

// form symmetric graph based on A * A' (trans = false) or A' * A (trans = true)
MATCL_LINALG_EXPORT Matrix symgraph_prod(const Matrix& A, bool trans = true);

// form symmetric graph based on [0 A'; A 0] (if trans = false) or [0 A; A' 0]
// (if trans = true)
MATCL_LINALG_EXPORT Matrix symgraph_bipart(const Matrix& A, bool trans = true);

// form adjacency from a symmetric matrix A as 1/2(A + A') without main diagonal; 
// only lower triangular part of the matrix A is considered
MATCL_LINALG_EXPORT Matrix make_adjancency_matrix(const Matrix& A);

// create the Laplacian matrix L for symmetric a matrix A; if use weights then 
// wighted Laplacian matrix LW with off diagonal elements LW_ij = L_ij * A_ij
// and diagonal elements LW_ii = -sum_{i != j} LW_ij; weights should be positive,
// otherwise the matrix L will not have the Laplacian matrix property, that 
// all eigenvalues of L are positive or zero
MATCL_LINALG_EXPORT Matrix symgraph_laplacian(const Matrix& A, bool use_weights = false);

// calculate first N Fieldlers vector of Laplacian matrix L, i.e. eigevectors for the
// first N smallest nonzero eigenvalue of L (including 0 eigenvalue); this function
// uses eig function and can be very costly
// [F, E] = fieldler_vector(L, N, direct)
// F       : the vector
// E       : eigenvalues of L
// direct  = true:  use eigsel_index function
//         = false: use iterative eigensolver (ARPACK) through pschur_decomposition
MATCL_LINALG_EXPORT mat_tup_2 fiedler_vectors(const Matrix& L, Integer N, bool direct);

// for given aggregation matrix aggreg (i.e. matrix of size N x L with elements
// 0 or 1, where N is number of nodes, L is number of groups, such that 
// sum(aggreg, 2) = 1) return permutation vector p of length N such that nodes
// belonging to given group k are not reordered and put after all nodes belonging 
// to groups 1:k-1
MATCL_LINALG_EXPORT permvec aggreg_to_perm(const Matrix& aggreg);

// 
MATCL_LINALG_EXPORT tuple<permvec,permvec>
biperm_to_perms(const permvec& p, Integer rows, Integer cols, bool trans = true);

};