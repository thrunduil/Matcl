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
#include "matcl-linalg/general/config_linalg.h"
#include "matcl-matrep/matrix/permvec.h"

namespace matcl
{

// correction algorithm used by cholmod
struct correction_alg
{
    enum type   { NONE, GMW, SE };
};

// type of correction used by cholmod
struct correction_type
{
    enum type   { TYPE_I,TYPE_II };
};
        
// return type of cholmod
using cholmod_return_type = tuple<Matrix,permvec, Integer, Real>;

/**
*  Purpose
*  =======
*  compute modified Cholesky decomposition of a symmetric or hermitian matrix
*  of the form
*
*               A(P,P) + E = U' * U     if upper = true
*               A(P,P) + E = U  * U'    if upper = false
*           
*  where U is an upper triangular matrix, P is a permutation vector. E is diagonal
*  matrix. if A is sufficiently positive definite, then E = 0.  Generally 
*  ||E|| / abs (min (lambda_k)) is small and in worst case bounded by O(N), 
*  N = dim(A), where lambda_k +is k-th eigenvalue of A.
*       
*  USAGE: 
*  =======
*           [U,P,step_I,norm_E] = cholmod(A, upper, corr, type)
*  Inputs
*  =======
*  A        - symmetrix or hermitian matrix; symmetry is not checked
*  upper    = true      : return upper traingular matrix,
*           = false     : return lowe triangular matrix
*  alg      - correction algorithm
*           = NONE      : no modification, factorization stops if the matrix A 
*                         appeared to be indefinite. See notes below.
*           = GMW       : Gill, Murray, and Wright (GMW) algorithm with 2 phase
*                         strategy as in Schnabel and Eskow, proposed by Feng.
*           = SE        : Schnabel and Eskow (SE99) algorithm, default
*  corr     - correction type
*           = TYPE_I    - type I correction, see notes below.
*           = TYPE_II   - type II correction, see notes below, default.
*  tol      - small value tolerance, if tol < 0, then default value is used
*       
*  Output
*  =======
*   U       - upper triengular matrix
*   P       - permutation vector
*   step_I  - number of steps in phase I
*   norm_E  - second norm of a diagonal matrix E.
*
*  Comments
*  =======
*  This routine implements two phase modified Cholesky decomposition using 
*  BLAS 3. In the first phase the standard Cholesky decomposition is performed
*  as long as Schur complement matrix is suffiently positive definite and bounds
*  on resulting matrices U or L can be guarantied. If the matrix A is positive
*  definite, then the standard Cholesky decomposition with complete pivoting is
*  returned. If the matrix A is semi-positive definite, then step_I is estimated
*  rank of the matrix A.
*
*  If phase II diagonal elements d_k are increased to 
*
*      max (|d_k|, delta)      for type I corrections
*      max (d_k, delta)        for type II corrections
*
*  Not available for sparse and band matrices.
*
*  References
*  =======
*  [1] P. E. Gill, W. Murray, and M. H. Wright. Practical Optimization. Academic
*      Press, 1981.
*  [2] R. B. Schnabel and E. Eskow. A new modified Cholesky factorization. SIAM
*      J. Sci. Stat. Comput., 11:1136{1158, 1990.
*  [3] R. B. Schnabel and E. Eskow. A revised modified Cholesky factorization
*      algorithm. SIAM J. Optim., 9(4):1135{1148, 1999.
*  [4] H. Feng, Matrix factorizations, triadic matrices, and modified
*      Cholesky factorizations for optimization, dissertation, 2006.
*/
MATCL_LINALG_EXPORT 
cholmod_return_type  cholmod(const Matrix& A, bool upper,correction_alg::type alg = correction_alg::SE,
                         correction_type::type corr = correction_type::TYPE_II, Real tol = -1);
MATCL_LINALG_EXPORT 
cholmod_return_type  cholmod(Matrix&& A, bool upper,correction_alg::type alg = correction_alg::SE,
                         correction_type::type corr = correction_type::TYPE_II, Real tol = -1);

};