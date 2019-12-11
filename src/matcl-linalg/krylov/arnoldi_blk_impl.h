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

#include "matcl-linalg/decompositions/eig_functions.h"
#include "arnoldi_process.h"
#include "matcl-linalg/decompositions/eig/arpack_wrapper.h"
#include "matcl-linalg/special_matrices/matrix_functors.h"

namespace matcl { namespace raw
{

namespace md = matcl::details;

// perform Arnoldi reduction:
//     A * V = V * H + r*E_{k}^T
//     V' * V = I, V' * r = 0.
// where V is N x k matrix, H is k x k block upper Hessenberg matrix, r is a N x KB matrix
// E_k is k x KB matrix formed from last KB columns of identity matrix of size kxk
//
// if lanczos = true, then Lanczos reduction is performed and H is additionally hermitian matrix. 
// The matrix A must be symmetric/hermitian.
template<class T>
class arnoldi_blk_impl
{
    private:        
        using TR            = typename md::real_type<T>::type;
        using Mat_D         = raw::Matrix<T, struct_dense>;
        using Mat_DR        = raw::Matrix<TR, struct_dense>;

    private:
        arnoldi_blk_impl(const arnoldi_blk_impl&) = delete;
        arnoldi_blk_impl& operator=(const arnoldi_blk_impl&) = delete;

    public:
        arnoldi_blk_impl(const linear_operator& A, bool lanczos, Integer max_k, Integer max_kb);

        ~arnoldi_blk_impl();

        void                resize(Integer max_k, Integer max_kb);
        void                clear();
        Integer             run(const matcl::Matrix& v, Integer k, Real tol);

        Integer             continue_run(Integer k, Real tol);
        Integer             continue_run(const matcl::Matrix& v, Integer k, Real tol);

        // return number of rows of the operator A;
        Integer             size() const;

        // return maximum number of Arnoldi vectors that can be computed
        Integer             max_k() const;

        // return maximum block size allowed
        Integer             max_kb() const;

        // return number of computed Arnoldi vectors
        Integer             number_vectors() const;

        // return current block size
        Integer             current_block_size() const;

        // return block size of the initial vecto
        Integer             get_H_ldiags() const;

        // get N x k matrix with Arnoldi vectors
        matcl::Matrix       get_V() const;

        mat_tup_2           get_V_as_reflectors() const;

        // get k x k upper hessenberg matrix, this matrix is symmetric if Lanczos
        // argorithm is used
        matcl::Matrix       get_H() const;

        // return N x 1 vector of residuals
        matcl::Matrix       get_resid() const;

        // return second norm of residuals
        TR                  get_norm_resid() const;

        const linear_operator& 
                            get_operator() const;

    private:
        void                init(const linear_operator& A, bool lanczos, Integer max_k, Integer max_kb);
        value_code          get_value_code() const;
        void                init_v(const matcl::Matrix& v, bool cont);
        Integer             run_impl(bool need_init, Integer k, Integer K, Real tol);
        void                init_arrays();
        void                init_workspace(Integer N, Integer max_k, Integer KB);

    private:
        Integer             m_N;
        Integer             m_k;
        Integer             m_KB;
        Integer             m_max_k;
        Integer             m_max_kb;
        linear_operator     m_A;
        bool                m_lanczos;

        matcl::Matrix       m_arnoldi_vec;
        T*                  m_arnoldi_vec_ptr;        
        Integer             m_arnoldi_vec_ld;
        matcl::Matrix       m_tau;
        T*                  m_tau_ptr;

        matcl::Matrix       m_hess;
        Integer             m_hess_ld;
        bool                m_hess_struct;
        Integer             m_hess_ldiags;
        TR                  m_scale;
        TR                  m_sumsq;

        matcl::Matrix       m_x;
        matcl::Matrix       m_res;
        TR                  m_norm_res;

        matcl::Matrix       m_work;
        T*                  m_work_ptr;
        Integer             m_work_size;

        matcl::Matrix       m_T;
        T*                  m_T_ptr;
        Integer             m_T_ld;
        Integer             m_NT;
        Integer             m_NB;
};

}};