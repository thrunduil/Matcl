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

#include "matcl-linalg/decompositions/eig_functions.h"
#include "arnoldi_process.h"
#include "matcl-linalg/decompositions/eig/arpack_wrapper.h"
#include "matcl-linalg/special_matrices/matrix_functors.h"

namespace matcl { namespace raw
{

namespace md = matcl::details;

// perform Arnoldi reduction:
//     A * V = V * H + r*e_{k}^T
//     V' * B * V = I, V' * B* r = 0.
// where V is N x k matrix, H is k x k upper Hessenberg matrix, r is a vector
// e_k is k x 1 vector with zeroes on all positions except last and 1 on the last
// position; B is semi positive definite matrix
//
// if lanczos = true, then Lanczos reduction is performed and H is hermitian tridiagonal
// real matrix. The matrix A must be symmetric/hermitian with respect to inner product
// defined by the matrix B, i.e.
//     B * A = A' * B, or equivalently < x,A y > = < A x,y >  
// where < z,w > = z'Bw (for example A = S * B, where S is symmetric/hermitian)
template<class T>
class arnoldi_impl
{
    private:        
        using TR            = typename md::real_type<T>::type;
        using Mat_D         = raw::Matrix<T, struct_dense>;
        using Mat_DR        = raw::Matrix<TR, struct_dense>;

    private:
        arnoldi_impl(const arnoldi_impl&) = delete;
        arnoldi_impl& operator=(const arnoldi_impl&) = delete;

    public:
        // prepare Arnoldi process for matrix A; allow for storing at most max_k 
        // Arnoldi vectors; if lanczos = true, then call Lanczos version; in this
        // case A must be hermitian 
        arnoldi_impl(const linear_operator& A, bool lanczos, Integer max_k);
        // version with B != I
        arnoldi_impl(const linear_operator& A, const linear_operator& B, 
                     bool lanczos, Integer max_k);

        ~arnoldi_impl();

        // change maximum number of allowed Arnoldi vectors
        void                resize(Integer max_k);

        // discard all previously computer Arnoldi vectors; capacity is not changed
        void                clear();

        // perform Arnoldi process for the vector v; perform at most min(k, max_k)
        // steps; finish iterations if norm of residual vector is less than tol;
        // return number of steps performed; 
        // previous computations are discarded
        Integer             run(const matcl::Matrix& v, Integer k, Real tol);

        // perform additional k steps of the Arnoldi process
        // return number of steps performed
        Integer             continue_run(Integer k, Real tol);

        // perform additional k steps of the Arnoldi process with new initial vector v
        // at the beginning the vector v is orthogonalized againts get_V()
        // return number of steps performed
        Integer             continue_run(const matcl::Matrix& v, Integer k, Real tol);

        // return number of rows of the operator A;
        Integer             size() const;

        // return maximum number of Arnoldi vectors that can be computed
        Integer             max_k() const;

        // return number of computed Arnoldi vectors
        Integer             number_vectors() const;

        // get N x k matrix with Arnoldi vectors
        matcl::Matrix       get_V() const;

        // get k x k upper hessenberg matrix, this matrix is symmetric if Lanczos
        // argorithm is used
        matcl::Matrix       get_H() const;

        // return N x 1 vector of residuals
        matcl::Matrix       get_resid() const;

        // return second norm of residuals
        TR                  get_norm_resid() const;

        const linear_operator& 
                            get_operator() const;
        const linear_operator& 
                            get_operator_B() const;

    private:
        void                init_arrays();
        void                init(const linear_operator& A, bool lanczos, Integer max_k);        
        Integer             run_impl(Integer k, Integer K, Real tol);
        void                init_work_v(const matcl::Matrix& v);
        value_code          get_value_code() const;
        bool                orthogonalize_resid();
        void                assign_to_work(Integer start, Integer N, const matcl::Matrix& val);

    private:
        Integer             m_N;
        Integer             m_k;
        Integer             m_max_k;

        linear_operator     m_A;
        linear_operator     m_B;

        bool                m_has_B;
        bool                m_lanczos;

        matcl::Matrix       m_arnoldi_vec;
        T*                  m_arnoldi_vec_ptr;        
        Integer             m_arnoldi_vec_ld;

        matcl::Matrix       m_hess;
        Integer             m_hess_ld;

        TR                  m_resid_norm;
        TR                  m_scale;
        TR                  m_sumsq;

        matcl::Matrix       m_work;
        T*                  m_work_ptr;
};

}};