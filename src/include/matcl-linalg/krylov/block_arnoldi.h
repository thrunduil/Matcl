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
    class block_arnoldi_impl;
}};

#pragma warning(push)
#pragma warning(disable:4251)  //needs to have dll-interface to be used by clients of class

namespace matcl
{

/// perform Arnoldi reduction:
///     A * V = V * H + r*E_{k}^T   (1)
///     V' * V = I, V' * r = 0.
/// where V is N x k matrix, H is k x k block upper Hessenberg matrix, r is a N x KB matrix
/// E_k is k x KB matrix formed from last KB columns of identity matrix of size kxk; the first
/// KB colums of V are spanned by the initial vector v.
///
/// if the operator A is hermitian, then Lanczos reduction is performed and H is additionally
/// hermitian matrix.
///
/// this class implements block householder Arnoldi method with pivoted qr to detect colinearity;
/// see for example
///     Baglama, James. "Augmented block householder Arnoldi method." 
///     Linear algebra and its applications 429.10 (2008): 2315-2334.
class MATCL_LINALG_EXPORT block_arnoldi_iteration
{
    private:
        using impl_type = std::shared_ptr<details::block_arnoldi_impl>;

        friend details::block_arnoldi_impl;

    private:
        impl_type               m_impl;

    public:
        /// initialize Arnoldi iterations for the linear operator given by a zero scalar
        block_arnoldi_iteration();

        /// initialize Arnoldi iterations for a linear operator A of size NxN, preallocate
        /// space for at most init_max_k Arnoldi vectors and iterations with block size at
        /// most init_max_kb; one can require higher number of Arnoldi vectors or use
        /// block size higher than init_max_kb, but then additional allocations will
        /// take place.
        /// If A is hermitian, then Lanczos iteration will be used
        block_arnoldi_iteration(const linear_operator& A, Integer init_max_k = 12, 
                                Integer init_max_kb = 4);

        /// standard destructor
        ~block_arnoldi_iteration();

        /// clear all results and initialize decomposition of another operator; 
        /// see constructor for details
        block_arnoldi_iteration&
                                operator()(const linear_operator &A, Integer init_max_k = 12, 
                                        Integer init_max_kb = 4);

        /// generate at most k Arnoldi vectors starting from the matrix v of size 
        /// N x KB. Stop the iteration if norm of the residual vector r satisfies
        /// |r|_2 <= tol * |H|_2, where H is the block Hessenberg matrix associated
        /// with the Arnoldi vectors computed so far; 
        /// return number of generated Arnoldi vectors
        Integer                 run(const matcl::Matrix& v, Integer k, Real tol);

        /// generate additional k Arnoldi vectors (at most); stop iteration if tolerance
        /// tol is reached; see run function for details; return number of generated 
        /// Arnoldi vectors
        Integer                 continue_run(Integer k, Real tol);

        /// generate additional k Arnoldi vectors starting from the new Matrix v of 
        /// size N x KB, KB can be different than previously used; stop iterations if
        /// tolerance tol is reached; return number of generated vectors; this function
        /// can be used only if deflation space is already found, otherwise (1) will not
        /// hold
        Integer                 continue_run(const matcl::Matrix& v, Integer k, Real tol);

        /// preallocate space for at most max_k Arnoldi vectors and block size at most
        /// max_kb; calling this function is not required but setting proper values will
        /// avoid memory allocations
        void                    resize(Integer max_k, Integer max_kb);

        /// clear all results; momery is not released; in order to release memory one should
        /// reset linear operator
        void                    clear();

        /// maximum number of Arnoldi vectors that can be stored without additional memory
        /// allocations
        Integer                 max_k() const;

        /// maximum block size that can be used withot additional memory allocations
        Integer                 max_kb() const;

        /// number of Arnoldi vectors generated so far
        Integer                 number_vectors() const;

        /// current block size of the residual vectors, i.e. r is a matrix of size 
        /// N x current_block_size(), block size can change due to call to continue_run
        /// starting from another matrix or if colinearity of residual vectors is detected
        Integer                 current_block_size() const;

        /// return the operator A
        const linear_operator&  get_operator() const;

        /// get N x k matrix with orthogonal Arnoldi vectors (the matrix V), 
        /// k = number_vectors()
        matcl::Matrix           get_V() const;

        /// the matrix V is stored as a product of elementary reflectors,
        /// [U,tau] = get_V_as_reflectors()
        /// U - unit lower triangular matrix of size N x k with elementary reflectors 
        ///     (k = number_vectors() )
        /// tau - vector of scalar factors of elementary reflectors of size k
        /// the matrix V can be created by calling householder_to_unitary(U,tau,k)
        mat_tup_2               get_V_as_reflectors() const;

        /// return the block upper Hessenberg matrix H with at most HLD subdiagonals, 
        /// where HLD = get_H_ldiags(); if the operator A is hermitian, then H is also
        /// hermitian.
        matcl::Matrix           get_H() const;

        /// number of subdiagonals of the block Hessenberg matrix H, this is at most
        /// 2*KB - 1, where KB is the block size of the initial vector v; 
        Integer                 get_H_ldiags() const;

        /// return N x kb vector of residuals, kb = current_block_size()
        matcl::Matrix           get_resid() const;

        /// return second norm of residuals
        Real                    get_norm_resid() const;

    private:
        void                    initialize(const linear_operator &A, Integer max_k, Integer max_kb);
};

};

#pragma warning(pop)