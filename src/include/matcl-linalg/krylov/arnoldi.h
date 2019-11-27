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
#include "matcl-linalg/special_matrices/unitary_matrix.h"

namespace matcl { namespace details
{
    class arnoldi_impl;
    class arnoldi_b_impl;
}};

#pragma warning(push)
#pragma warning(disable:4251)  //needs to have dll-interface to be used by clients of class

namespace matcl
{

/// perform Arnoldi reduction:
///     A * V = V * H + r*E_{k}^T   (1)
///     V' * V = I, V' * r = 0.
/// where V is N x k matrix, H is k x k upper Hessenberg matrix, r is a N x 1 vector
/// E_k is k x 1 matrix formed from last column of identity matrix of size kxk; the first
/// colum of V is spanned by the initial vector v.
///
/// if the operator A is hermitian, then Lanczos reduction is performed and H is additionally
/// hermitian, real tridiagonal matrix.
///
/// this class implements Arnoldi methods with Gram-Schmidt reorthogonalization
class MATCL_LINALG_EXPORT arnoldi_iteration
{
    private:
        using impl_type = std::shared_ptr<details::arnoldi_impl>;

        friend details::arnoldi_impl;

    private:
        impl_type               m_impl;

    public:
        /// initialize Arnoldi iterations for the linear operator given by a zero scalar
        arnoldi_iteration();

        /// initialize Arnoldi iterations for a linear operator A of size NxN, preallocate
        /// space for at most init_max_k Arnoldi vectors ; one can require higher number of 
        /// Arnoldi vectors, but then additional allocations will take place.
        /// If A is hermitian, then Lanczos iteration will be used
        arnoldi_iteration(const linear_operator& A, Integer init_max_k = 12);

        /// standard destructor
        ~arnoldi_iteration();

        /// clear all results and initialize decomposition of another operator; 
        /// see constructor for details
        arnoldi_iteration&      operator()(const linear_operator &A, Integer init_max_k = 12);

        /// generate at most k Arnoldi vectors starting from the vector v of size 
        /// N x 1. Stop the iteration if norm of the residual vector r satisfies
        /// |r|_2 <= tol * |H|_2, where H is the block Hessenberg matrix associated
        /// with the Arnoldi vectors computed so far; return number of generated Arnoldi
        /// vectors
        Integer                 run(const matcl::Matrix& v, Integer k, Real tol);

        /// generate additional k Arnoldi vectors (at most); stop iteration if tolerance
        /// tol is reached; see run function for details; return number of generated 
        /// Arnoldi vectors
        Integer                 continue_run(Integer k, Real tol);

        /// generate additional k Arnoldi vectors starting from the new vector v of 
        /// size N x 1; stop iterations if tolerance tol is reached; return number of 
        /// generated vectors; this function can be used only if deflation space is already 
        /// found, otherwise (1) will not hold
        Integer                 continue_run(const matcl::Matrix& v, Integer k, Real tol);

        /// preallocate space for at most max_k Arnoldi vectors; calling this function is not
        /// required but setting proper values will avoid memory allocations
        void                    resize(Integer max_k);

        /// clear all results; momery is not released; in order to release memory one should
        /// reset linear operator
        void                    clear();

        /// maximum number of Arnoldi vectors that can be stored without additional memory
        /// allocations
        Integer                 max_k() const;

        /// number of Arnoldi vectors generated so far
        Integer                 number_vectors() const;

        /// return the operator A
        const linear_operator&  get_operator() const;

        /// get N x k matrix with orthogonal Arnoldi vectors (the matrix V), 
        /// k = number_vectors()
        matcl::Matrix           get_V() const;

        /// return the upper Hessenberg matrix H w; if the operator A is hermitian, then H 
        /// is also hermitian, real, tridiagonal.
        matcl::Matrix           get_H() const;

        /// return N x 1 vector of residuals
        matcl::Matrix           get_resid() const;

        /// return second norm of residuals
        Real                    get_norm_resid() const;

    private:
        void                    initialize(const linear_operator &A, Integer max_k);
};

/// perform Arnoldi reduction with respect to inner product given by a matrix B:
///     A * V = V * H + r*E_{k}^T           (1)
///     V' * B * V = I, V' * B * r = 0.
/// where V is N x k matrix, H is k x k upper Hessenberg matrix, r is a N x 1 vector
/// E_k is k x 1 matrix formed from last column of identity matrix of size kxk; the first
/// colum of V is spanned by the initial vector v. The matrix B must be hermitian and
/// semi positive definite.
///
/// if the operator A is hermitian with respect to inner product given by the matrix B, 
/// then Lanczos reduction is performed and H is additionally hermitian, real tridiagonal
/// matrix. A is hermitian with respect to inner product given by B if
///     B * A = A' * B, or equivalently < x,Ay > = < Ax,y >
/// where <z,w> = z'Bw.
///
/// this class implements Arnoldi methods with Gram-Schmidt reorthogonalization
class MATCL_LINALG_EXPORT arnoldi_b_iteration
{
    private:
        using impl_type = std::shared_ptr<details::arnoldi_b_impl>;

        friend details::arnoldi_b_impl;

    private:
        impl_type               m_impl;

    public:
        /// initialize Arnoldi iterations for the linear operator given by a zero scalar
        /// and B = 1.0
        arnoldi_b_iteration();

        /// initialize Arnoldi iterations for a linear operators A, B of size NxN, preallocate
        /// space for at most init_max_k Arnoldi vectors ; one can require higher number of 
        /// Arnoldi vectors, but then additional allocations will take place.
        /// If lanczos = true, then assume that A is hermitian with respect to inner product
        /// given by B and use then Lanczos iteration
        arnoldi_b_iteration(const linear_operator& A, const linear_operator& B, 
                            bool lanczos = false, Integer init_max_k = 12);

        /// standard destructor
        ~arnoldi_b_iteration();

        /// clear all results and initialize decomposition of another operators; 
        /// see constructor for details
        arnoldi_b_iteration&    operator()(const linear_operator &A, const linear_operator &B, 
                                           bool lanczos = false, Integer init_max_k = 12);

        /// generate at most k Arnoldi vectors starting from the vector v of size 
        /// N x 1. Stop the iteration if norm of the residual vector r satisfies
        /// |r|_2 <= tol * |H|_2, where H is the block Hessenberg matrix associated
        /// with the Arnoldi vectors computed so far; return number of generated Arnoldi
        /// vectors
        Integer                 run(const matcl::Matrix& v, Integer k, Real tol);

        /// generate additional k Arnoldi vectors (at most); stop iteration if tolerance
        /// tol is reached; see run function for details; return number of generated 
        /// Arnoldi vectors
        Integer                 continue_run(Integer k, Real tol);

        /// generate additional k Arnoldi vectors starting from the new vector v of 
        /// size N x 1; stop iterations if tolerance tol is reached; return number of 
        /// generated vectors; this function can be used only if deflation space is already 
        /// found, otherwise (1) will not hold
        Integer                 continue_run(const matcl::Matrix& v, Integer k, Real tol);

        /// preallocate space for at most max_k Arnoldi vectors; calling this function is not
        /// required but setting proper values will avoid memory allocations
        void                    resize(Integer max_k);

        /// clear all results; momery is not released; in order to release memory one should
        /// reset linear operator
        void                    clear();

        /// maximum number of Arnoldi vectors that can be stored without additional memory
        /// allocations
        Integer                 max_k() const;

        /// number of Arnoldi vectors generated so far
        Integer                 number_vectors() const;

        /// return the operator A
        const linear_operator&  get_operator() const;

        /// return the operator B
        const linear_operator&  get_operator_B() const;

        /// get N x k matrix with orthogonal Arnoldi vectors (the matrix V), 
        /// k = number_vectors()
        matcl::Matrix           get_V() const;

        /// return the upper Hessenberg matrix H w; if the operator A is hermitian, then H 
        /// is also hermitian, real, tridiagonal.
        matcl::Matrix           get_H() const;

        /// return N x 1 vector of residuals
        matcl::Matrix           get_resid() const;

        /// return second norm of residuals
        Real                    get_norm_resid() const;

    private:
        void                    initialize(const linear_operator &A, const linear_operator &B, 
                                           bool lanczos, Integer max_k);
};

};

#pragma warning(pop)