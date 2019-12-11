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

// construct Givens plane rotation such that
//     [c   s] [a] - [r]
//     [-s' c] [b] - [0]
//
// additionally c is a real scalar and |c|^2 + |s|^2 = 1
template<class Val, class Enable = typename md::enable_if_float<Val,void>::type>
void construct_givens(const Val& a, const Val& b, typename md::real_type<Val>::type& c, Val& s);

// construct Givens plane rotation as construct_givens; return additionally r
template<class Val, class Enable = typename md::enable_if_float<Val,void>::type>
void construct_givens(const Val& a, const Val& b, typename md::real_type<Val>::type& c, Val& s, Val& r);

// construct Givens plane rotation as construct_givens; it is more stable but slower
// version based on Lapack's routine lartg
template<class Val, class Enable = typename md::enable_if_float<Val,void>::type>
MATCL_LINALG_EXPORT
void construct_givens2(const Val& a, const Val& b, typename md::real_type<Val>::type& c, Val& s, Val& r);

// construct unitary matrix from a sequence of plane rotations;
//     U = givens_to_unitary(N, C, S, Ind)
// N       - size of the square unitary matrix U
// C       - a vector of cosines of length K containing real values
// S       - a vector of sines of length K; elements of C and S vectors should be 
//           constructed using construct_givens functions
// Ind     - integer matrix of size K x 2; each element in Ind vector should be
//           in range [1, N]
// left    = true: the sequence defined by C, S, Ind was generated from the left
//         = false: the sequence defined by C, S, Ind was generated from the right
// U       - unitary matrix representing the sequence of Givens plane rotations
//
// The i-th elements of the left (L-) or right (R-) sequence (C, S, IND1, IND2), where
// Ind = [IND1, IND2], indicates that in the i-th step of construction of the sequence
// rows or columns given by IND1 and IND2 of a matrix are multiplied by the plane rotation
// G_i constructed from C(i) and S(i), i.e. the following multiplication was performed:
//
//    [X(:,IND1) X(:,IND2)] := [X(:,IND1) X(:,IND2)] * GR_i   for R-sequences
//    [X(IND1,:);X(IND2,:)] := GL_i * [X(IND1,:);X(IND2,:)]   for L-sequences 
//
//  for some matrix X where
//       GL_i = [C(i)   S(i)],   GR_i = [C(i)  -S'(i)]
//              [-S(i)' C(i)]           [S(i)  C(i)]
//
//  in the following order:
//       GL_K * ... * GL_2 * GL_1 * X    for L-sequences
//       X  * GR_1 * GR_2 * ... * GR_K   for R-sequences
//
//  Multiplication of the unitary matrix constructed from the sequence (C, S, IND1, IND2)
//  of size N x N and a matrix Y, op(U) * Y or Y * op(U), is defined as:
//
//       op[GL_K * ... * GL_1] * Y   when L-sequence is applied from the left
//       Y * op[GR_1 * ... * GR_K]   when R-sequence is applied from the right
//       Y * op[GL_K * ... * GL_1]   when L-sequence is applied from the right
//       op[GR_1 * ... * GR_K] * X   when R-sequence is applied from the left
//
//  where op is one of the sequence operators
//
//       op[G_1 * ... * G_k] = G_1 * ... * G_k               if op = trans_type::no_trans
//       op[G_1 * ... * G_k] = trans(G_k) * ... * trans(G_1) if op = trans_type::trans
//       op[G_1 * ... * G_k] = G_k' * ... * G_1'             if op = trans_type::conj_trans
MATCL_LINALG_EXPORT unitary_matrix givens_to_unitary(Integer N, const Matrix& C, const Matrix& S, 
                                                    const Matrix& Ind, bool from_left);

};

#include "matcl-linalg/details/givens.inl"
