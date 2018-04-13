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

#include "matcl-matrep/matcl_matrep.h"
#include "matcl-linalg/general/config_linalg.h"
#include "matcl-linalg/graph/graph_manip.h"
#include "matcl-linalg/graph/scotch.h"
#include "matcl-linalg/graph/metis.h"
#include "matcl-linalg/graph/dm_decomposition.h"

namespace matcl
{

/// graph algorithms consider structure of the matrix, zeros are not ignored

//-----------------------------------------------------------------------
//                          REORDERING
//-----------------------------------------------------------------------

/// nested dissection ordering of a symmetric matrix from Petsc
MATCL_LINALG_EXPORT permvec order_nested_dissection(const Matrix& A);

/// one-way dissection ordering of a symmetric matrix from Petsc
MATCL_LINALG_EXPORT permvec order_one_way_dissection(const Matrix& A);

/// reverse Cuthill-McKee ordering of a symmetric matrix from Petsc
MATCL_LINALG_EXPORT permvec order_rcm(const Matrix& A);

// the column approximate minimum degree ordering algorithm computes
// a permutation vector p such that the LU and QR factorization of A(:,p)
// tends to be sparser than that of A;  the Cholesky factorization of
// (A (:,p))'*(A (:,p)) will also tend to be sparser than that of A'*A;
MATCL_LINALG_EXPORT permvec order_colamd(const Matrix& A);

/// an approximate minimum degree ordering algorithm for Cholesky
/// factorization of symmetric matrices.
MATCL_LINALG_EXPORT permvec order_symamd(const Matrix& A);

/// multiple minimum degree ordering of a symmetric matrix
MATCL_LINALG_EXPORT permvec order_mmd(const Matrix& A);

/// find a row and column permutations p, q of a matrix A, such that 
/// A(p,q) has large entries on the diagonal and row and column
/// scaling vector such that the scaled and permuted matrix B,
/// B = scale_rowscols(A, Dr, Dc)(p,q) has maximum number of nonzeros
/// on diagonal, all nonzero diagonal elements c_ii satisfy |c_ii| = 1
/// and all off-diagonal elements c_ij satisfy |c_ij| <= 1; this function
/// uses max_match_weighted2 function
/// [p, q, Dr, Dc] = order_diag(A, bal_off)
/// A       : a matrix
/// bal_off : if true and A is square then offdiagonal elements are reduced
///           using balance_offdiag function
/// p,q     : row and column permutation vectors
/// Dr, Dc  : row and column scaling vectors
MATCL_LINALG_EXPORT tuple<permvec,permvec, Matrix, Matrix> 
order_diag(const Matrix& A, bool reduce_offdiags = true);

MATCL_LINALG_EXPORT tuple<permvec,permvec, Matrix, Matrix> 
order_diag(Matrix&& A, bool reduce_offdiags = true);

/// ordering of a symmetric matrix to reduce fill-ins for Cholesky 
/// factorization from scotch library, see scotch for details
MATCL_LINALG_EXPORT permvec order_scotch(const Matrix& A);

/// ordering of a symmetric matrix to reduce fill-ins for Cholesky 
/// factorization from metis library, see metis for details
MATCL_LINALG_EXPORT permvec order_metis(const Matrix& A);

/// spectral ordering of a symmetric matrix; this ordering tries to reorder
/// similar vertices close together using spectral information from the Laplacian
/// matrix of the graph A; when use_weigts = true, then values stored in the 
/// matrix are used to construct the Laplacian matrix, otherwise values are 
/// ignored; weigts must be positive; if direct = true then direct eigensolver
/// is used, otherwise iterative, see fiedler_vectors for details;
/// this ordering is good for graph visualization
/// but is very costly
MATCL_LINALG_EXPORT permvec order_spectral(const Matrix& A, bool use_weigts, bool direct);

/// reorder column in nondecreasing order of nonzero count
MATCL_LINALG_EXPORT permvec colperm(const Matrix& A);

/// reorder rows in nondecreasing order of nonzero count
MATCL_LINALG_EXPORT permvec rowperm(const Matrix& A);

//-----------------------------------------------------------------------
//                          PARTITIONING
//-----------------------------------------------------------------------

/// coarse a symmetric N x N matrix A using heavy-edge matching approach from Petsc;
/// this algorithm tries to finds maximal matching with the highest edge weight
/// using local search starting from the first node; initial ordering as well
/// as element values (interpreted as weights associated with edges) are important;
/// values must be positive;
/// return N x K integer matrix C with N nonzero elements equal to 1; positions of
/// nonzero elements in column i give node indices that belongs to i-th group;
/// in this way C' * A * C is coarsen matrix and C has full column rank
MATCL_LINALG_EXPORT Matrix coarser_hem(const Matrix& A);

/// coarse a symmetric N x N matrix A using maximal independent set approach from
/// Petsc; this algorithm tries to find independent sets with nodes no two of which
/// are adjacent; initial ordering is important; values of elements are not considered
/// return N x K integer matrix C with N nonzero elements equal to 1; positions of
/// nonzero elements in column i give node indices that belongs to i-th group;
/// in this way C' * A * C is coarsen matrix and C has full column rank
MATCL_LINALG_EXPORT Matrix coarser_mis(const Matrix& A);

//-----------------------------------------------------------------------
//                          MATCHING
//-----------------------------------------------------------------------

/// find maximum matching in unsymmetrix (possibly not square) matrix using
/// depth-first search to find an augmenting path from each column node to get
/// the maximum matching; 
/// algorithm: Alex Pothen and Chin-Ju Fan, Penn State University
///     [r, c] = max_match(A)
/// A:  matrix of size M x N
/// r:  vector of size M x 1 with elements 0-N, where 0 means that given row
///     is unmatched and value i > 0 means that given row is matched with
///     column i
/// c:  vector of size N x 1 with elements 0-M, where 0 means that given column
///     is unmatched and value i > 0 means that given column is matched with
///     row i
MATCL_LINALG_EXPORT mat_tup_2 max_match(const matcl::Matrix& A);

/// find maximum weighted matching in unsymmetrix (possibly not square) matrix,
/// i.e. maximum matching with highest possible sum of weights assigned to 
/// edges in the matching; this function uses the Hungarian algorithm
///     [r, c] = max_match_weighted(A)
/// A:  matrix of size M x N
/// r:  vector of size M x 1 with elements 0-N, where 0 means that given row
///     is unmatched and value i > 0 means that given row is matched with
///     column i
/// c:  vector of size N x 1 with elements 0-M, where 0 means that given column
///     is unmatched and value i > 0 means that given column is matched with
///     row i
MATCL_LINALG_EXPORT mat_tup_2 max_match_weighted(const matcl::Matrix& A);

/// find maximum weighted matching as function max_match_weighted; additionally
/// return dual row and column weights u, v such that
///     u_i + v_j =  w_ij,  for (i,j) in M
///     u_i + v_j >= w_ij,  for (i,j) not in M, i in U, j in V
///     u_i = 0, v_j = 0,   for i not in U, j not in V
/// where w_ij is the weight associated with edge (i,j) and M is the matching
/// found, U is the set of matched rows, V is the set of matched columns
///     [r, c, u, v, f, k] = max_match_weighted(A)
/// A:  matrix of size M x N
/// r:  vector of size M x 1 with elements 0-N, where 0 means that given row
///     is unmatched and value i > 0 means that given row is matched with
///     column i
/// c:  vector of size N x 1 with elements 0-M, where 0 means that given column
///     is unmatched and value i > 0 means that given column is matched with
///     row i
/// u:  dual row weights of size M x 1
/// v:  dual column weights of size N x 1
/// f:  sum of weights of edges in the maximum matching; also sum(u) + sum(v)
/// k:  cardinality of the matching; also structural rank of the matrix A
MATCL_LINALG_EXPORT tuple<Matrix,Matrix,Matrix,Matrix, Real, Integer>
max_match_weighted2(const matcl::Matrix& A);

/// structural rank of the matrix A using maximum matching algorithm; 
/// always sprank(A) >= rank(A) and there exists an assignment of values to
/// nonzero elements of A such that resulting matrix has rank sprank(A)
/// and no assignment of values can obtain higher rank
MATCL_LINALG_EXPORT Integer sprank(const Matrix& A);

/// create row and column permutations from row matching for a matrix A of
/// size M x N, puts matched rows in upper-left corner and unmatched rows in 
/// bottom-right corner;
/// [p,q] = rowsmatch_to_perms(r, N)
/// r       : row matching of size M x 1 as returned by max_match
/// N       : number of columns of A
/// p, q    : row and column permutation vectors
MATCL_LINALG_EXPORT tuple<permvec, permvec> rowsmatch_to_perms(const Matrix& r, Integer N);

/// create row and column permutations from column matching for a matrix A of
/// size M x N, puts matched columns in upper-left corner and unmatched columns
/// in bottom-right corner;
/// [p,q] = colsmatch_to_perms(c, M)
/// c       : column matching of size N x 1 as returned by max_match
/// M       : number of rows of A
/// p, q    : row and column permutation vectors
MATCL_LINALG_EXPORT tuple<permvec, permvec> colsmatch_to_perms(const Matrix& c, Integer M);

/// reorder general matrix A into block triangular form using Dulmage-Mendelsohn 
/// decomposition; see dm_decomposition for details
/// [p, q] = dmperm(A)
/// p, q:   row and column permutation vectors
MATCL_LINALG_EXPORT tuple<permvec, permvec> dmperm(const Matrix& A);

};