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
#include "matcl-linalg/general/linalg_exception.h"
#include "matcl-linalg/special_matrices/matrix_functors.h"

namespace matcl
{

using ret_arnoldi   = tuple<Matrix, Matrix, Matrix, Real>;

// perform Arnoldi reduction:
//     A * V = V * H + r*E_{k}^T   (1)
//     V' * V = I, V' * r = 0.
// where V is N x k matrix, H is k x k block upper Hessenberg matrix, r is a N x KB matrix
// E_k is k x KB matrix formed from last KB columns of identity matrix of size kxk; the first
// KB colums of V are spanned by the initial vector v.
//
// if the operator A is hermitian, then Lanczos reduction is performed and H is additionally
// hermitian matrix.
//
// use block householder Arnoldi method with pivoted qr to detect colinearity;
//
//     [V, H, r, nr]   = arnoldi(A, v, k)
// A   - linear operator of size N * N
// v   - initial matrix of size N * kb
// k   - number of Arnoldi vectors to generate
// V   - unitary matrix of size N * kr with kr Arnoldi vectors, kr <= k
// H   - block upper Hessenberg matrix of size kr x kr, H is hermitian if A is hermitian
// r   - residuals of size N x kbr, 0 <= kbr <= kb
// nr  - the Frobenius norm of r
MATCL_LINALG_EXPORT ret_arnoldi blk_arnoldi(const linear_operator& A, const Matrix& v, Integer k,
                                            Real tol = 0.0);

// perform Arnoldi reduction:
//     A * V = V * H + r*E_{k}^T   (1)
//     V' * V = I, V' * r = 0.
// where V is N x k matrix, H is k x k upper Hessenberg matrix, r is a N x 1 vector
// E_k is k x 1 matrix formed from last column of identity matrix of size kxk; the first
// colum of V is spanned by the initial vector v.
//
// if the operator A is hermitian, then Lanczos reduction is performed and H is additionally
// hermitian, real tridiagonal matrix.
//
// use Gram-Schmidt reorthogonalization
//     [V, H, r, nr]   = arnoldi(A, v, k)
// A   - linear operator of size N * N
// v   - initial vectpr of size N * 1
// k   - number of Arnoldi vectors to generate
// V   - unitary matrix of size N * kr with kr Arnoldi vectors, kr <= k
// H   - upper Hessenberg matrix of size kr x kr, H is hermitian and real if A is hermitian
// r   - residuals of size N x 1
// nr  - the second norm of r
MATCL_LINALG_EXPORT ret_arnoldi arnoldi(const linear_operator& A, const Matrix& v, Integer k,
                                        Real tol = 0.0);

// perform Arnoldi reduction with respect to inner product given by a matrix B:
//     A * V = V * H + r*E_{k}^T           (1)
//     V' * B * V = I, V' * B * r = 0.
// where V is N x k matrix, H is k x k upper Hessenberg matrix, r is a N x 1 vector
// E_k is k x 1 matrix formed from last column of identity matrix of size kxk; the first
// colum of V is spanned by the initial vector v. The matrix B must be hermitian and
// semi positive definite.
//
// if the operator A is hermitian with respect to inner product given by the matrix B, 
// then Lanczos reduction is performed and H is additionally hermitian, real tridiagonal
// matrix. A is hermitian with respect to inner product given by B if
//     B * A = A' * B, or equivalently < x,Ay > = < Ax,y >
// where <z,w> = z'Bw.
//
// this class implements Arnoldi methods with Gram-Schmidt reorthogonalization
//
// use Gram-Schmidt reorthogonalization
//     [V, H, r, nr]   = arnoldi(A, B, v, k)
//     [V, H, r, nr]   = lanczos(A, B, v, k)
// A   - linear operator of size N * N, if lanczos algorithm is used, then A
//       must be hermitian with respect to inner product given by B
// B   - linear operator of size N * N, B must be hermitian and semi positive 
//       definite
// v   - initial vectpr of size N * 1
// k   - number of Arnoldi vectors to generate
// V   - unitary matrix of size N * kr with kr Arnoldi vectors, kr <= k
// H   - upper Hessenberg matrix of size kr x kr, H is hermitian and real if A is hermitian
// r   - residuals of size N x 1
// nr  - the second norm of r
MATCL_LINALG_EXPORT ret_arnoldi arnoldi(const linear_operator& A, const linear_operator& B, 
                                        const Matrix& v, Integer k, Real tol = 0.0);
MATCL_LINALG_EXPORT ret_arnoldi lanczos(const linear_operator& A, const linear_operator& B, 
                                        const Matrix& v, Integer k, Real tol = 0.0);

}