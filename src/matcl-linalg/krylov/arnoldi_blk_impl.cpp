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

#include "arnoldi_blk_impl.h"
#include "matcl-matrep/matcl_matrep.h"
#include "matcl-linalg/norms_error/norm.h"
#include "arnoldi_process.h"
#include "matcl-blas-ext/lapack_ext/lapack_ext.h"
#include "matcl-internals/func/lapack_utils.h"
#include "matcl-blas-lapack/level1/level1.h"

namespace matcl { namespace raw
{

template<class T>
arnoldi_blk_impl<T>::arnoldi_blk_impl(const linear_operator& A, bool lanczos, Integer max_k, Integer max_kb)    
{
    init(A, lanczos, max_k, max_kb);
}

template<class T>
arnoldi_blk_impl<T>::~arnoldi_blk_impl()
{};

template<class T>
void arnoldi_blk_impl<T>::init(const linear_operator& A, bool lanczos, Integer max_k, Integer max_kb)
{
    m_lanczos           = lanczos;
    m_A                 = A;
    m_N                 = A.rows();
    m_k                 = 0;
    m_KB                = 0;
    m_max_k             = std::max(0,max_k);
    m_max_kb            = std::max(1,max_kb);

    m_arnoldi_vec_ptr   = nullptr;
    m_arnoldi_vec_ld    = 1;
    m_tau_ptr           = nullptr;

    m_hess_ld           = 1;
    m_hess_struct       = false;
    m_hess_ldiags       = 1;

    m_work_ptr          = nullptr;
    m_work_size         = 0;

    m_T_ptr             = nullptr;
    m_T_ld              = 0;
    m_NT                = 0;

    m_scale             = TR(0.0);
    m_sumsq             = TR(1.0);
    m_norm_res          = TR(0.0);

    init_arrays();
};

template<class T>
void arnoldi_blk_impl<T>::init_arrays()
{
    if (m_lanczos == false)
    {
        m_hess          = zeros(m_max_k, m_max_k, get_value_code());
        m_hess_ld       = std::max(1,m_max_k);
    }
    else
    {
        Integer max_kb  = std::min(2*m_max_kb, m_max_k);
        m_hess          = zeros(max_kb, m_max_k, get_value_code());
        m_hess_ld       = std::max(1,max_kb);
    };

    m_arnoldi_vec       = make_dense_noinit(m_N, m_max_k, get_value_code());
    m_arnoldi_vec_ptr   = m_arnoldi_vec.get_array_unique<T>();
    m_arnoldi_vec_ld    = std::max(m_N,1);
    
    m_tau               = make_dense_noinit(m_max_k, 1, get_value_code());
    m_tau_ptr           = m_tau.get_array_unique<T>();

    m_NB                = baitr_block_size<T>();
    Integer NTB         = m_max_k/m_NB;
    m_T                 = make_dense_noinit(m_NB*m_NB, NTB, get_value_code());
    m_T_ld              = std::max(1,m_NB*m_NB);
    m_T_ptr             = m_T.get_array_unique<T>();
};

template<class T>
void arnoldi_blk_impl<T>::init_workspace(Integer N, Integer max_k, Integer KB)
{
    T work_query;
    Integer info;
    Integer NT = 0;
    Integer N1  = std::max(N,1);
    baitr<T>(m_lanczos, N, 0, max_k, KB, TR(0.0), nullptr, N1, true, nullptr, N1, m_norm_res, nullptr, N1, 
             nullptr, nullptr, nullptr, N1, m_scale, m_sumsq, nullptr, m_NB*m_NB, NT, &work_query, 
             -1, nullptr, info);

    Integer work_size   = (Integer)real(work_query);

    if (m_work_size == 0)
    {
        m_work          = make_dense_noinit(work_size, 1, get_value_code());
        m_work_ptr      = m_work.get_array_unique<T>(); 
        m_work_size     = work_size;
    }
    else
    {
        m_work_size     = std::max(m_work_size, work_size);
        m_work.resize(m_work_size, 1);
        m_work_ptr      = m_work.get_array_unique<T>(); 
    };
};

template<class T>
void arnoldi_blk_impl<T>::clear()
{
    m_k                 = 0;
    m_KB                = 0;
    m_hess_ldiags       = 1;
    m_hess_struct       = false;
    m_NT                = 0;
};

template<class T>
void arnoldi_blk_impl<T>::resize(Integer max_k, Integer max_kb)
{
    if (max_k <= m_max_k && max_kb <= m_max_kb)
        return;

    max_k               = std::max(max_k, 0);
    max_kb              = std::max(max_kb, m_max_kb);

    if (m_lanczos == false)
    {
        m_hess.resize(max_k, max_k);
        m_hess_ld       = m_hess.impl<Mat_D>().ld();
    }
    else
    {
        Integer kb      = std::min(2*max_kb, max_k);
        m_hess.resize(kb, max_k);
        m_hess_ld       = m_hess.impl<Mat_D>().ld();
    }

    m_arnoldi_vec.resize(m_N, max_k);
    m_arnoldi_vec_ptr   = m_arnoldi_vec.get_array_unique<T>();
    m_arnoldi_vec_ld    = m_arnoldi_vec.impl<Mat_D>().ld();

    m_tau.resize(max_k,1);
    m_tau_ptr           = m_tau.get_array_unique<T>();

    m_k                 = std::min(m_k, max_k);
    m_max_k             = max_k;
    m_max_kb            = max_kb;

    Integer NTB         = max_k/m_NB;
    m_T.resize(m_NB*m_NB, NTB);
    m_T_ptr             = m_T.get_array_unique<T>();
};

template<class T>
void arnoldi_blk_impl<T>::init_v(const matcl::Matrix& v, bool cont)
{
    m_x                 = v;
    m_res               = v;
    m_KB                = v.cols();
    m_hess_struct       = (cont == false)? (m_KB == 1) : (m_hess_struct && m_KB == 1);
    init_workspace(m_N, m_max_k, m_max_kb);
};

template<class T>
Integer arnoldi_blk_impl<T>::run(const matcl::Matrix& v, Integer k, Real tol)
{
    init_v(v, false);
    return run_impl(true, 0, std::min(k, m_max_k), tol);
};

template<class T>
Integer arnoldi_blk_impl<T>::continue_run(Integer k, Real tol)
{
    if (m_KB == 0)
        return 0;

    Integer K       = std::min(m_k + k, m_max_k);
    return run_impl(false, m_k, K, tol);
};

template<class T>
Integer arnoldi_blk_impl<T>::continue_run(const matcl::Matrix& v, Integer k, Real tol)
{
    init_v(v, true);    
    Integer K   = std::min(m_k + k, m_max_k);    

    return run_impl(true, m_k, K, tol);
};

template<class T>
Integer arnoldi_blk_impl<T>::run_impl(bool need_init, Integer k, Integer K, Real tol)
{ 
    Integer info    = 0;

    struct callback_impl : callback_baitr<T>
    {    
        arnoldi_blk_impl* owner;

        callback_impl(arnoldi_blk_impl* o)
            :owner(o)
        {};

        virtual void eval(Integer M, Integer NB, T* X, Integer LDX, T* Y, Integer LDY) override
        {
            //eval Y = OP * X 
            matcl::Matrix x     = make_dense_foreign(M, NB, X, LDX);
            matcl::Matrix y     = make_dense_foreign(M, NB, Y, LDY);
            owner->m_A.mmul_right(x, trans_type::no_trans, y);
        }
    };

    callback_impl callback(this);

    Integer np          = std::max(K - k, 0);

    T* H_S_ptr          = (m_lanczos == true) ? m_hess.get_array_unique<T>() : nullptr;
    T* H_NS_ptr         = (m_lanczos == false) ? m_hess.get_array_unique<T>() : nullptr;
    Mat_D mat_v         = m_x.impl_unique<Mat_D>();
    Mat_D mat_res       = m_res.impl_unique<Mat_D>();
    Integer KB_old      = m_KB;

    baitr<T>(m_lanczos, m_N, k, np, m_KB, (TR)tol, mat_v.ptr(), mat_v.ld(), need_init, mat_res.ptr(), mat_res.ld(),
             m_norm_res, m_arnoldi_vec_ptr, m_arnoldi_vec_ld, m_tau_ptr, H_NS_ptr, H_S_ptr, m_hess_ld, m_scale, m_sumsq,
             m_T_ptr, m_T_ld, m_NT, m_work_ptr, m_work_size, &callback, info);

    if (info < 0)
    {
        std::ostringstream msg;
        msg << "baitr returned code: " << info;
        throw error::error_general(msg.str());
    }        

    Integer k_new       = info > 0 ? info : K;
    Integer n_steps     = k_new - m_k;

    m_k                 = k_new;

    if (n_steps > 0 && m_lanczos == true)
    {
        m_hess_ldiags       = std::max(2*KB_old, m_hess_ldiags);
        m_hess_ldiags       = std::min(m_hess_ldiags, m_max_k);
    };

    return n_steps;
};

template<class T>
Integer arnoldi_blk_impl<T>::size() const
{
    return m_N;
};

template<class T>
Integer arnoldi_blk_impl<T>::max_k() const
{
    return m_max_k;
};

template<class T>
Integer arnoldi_blk_impl<T>::max_kb() const
{
    return m_max_kb;
};

template<class T>
Integer arnoldi_blk_impl<T>::number_vectors() const
{
    return m_k;
};

template<class T>
matcl::Matrix arnoldi_blk_impl<T>::get_V() const
{
    using md::lap;

    Mat_D V(ti::ti_empty(), T(0.0), m_N, m_k);
    T* ptr_V        = V.ptr();
    Integer V_ld    = V.ld();

    for (Integer i = 0; i < m_k; ++i)
        ptr_V[i + i * V_ld] = T(1.0);

    Integer NT  = m_NT;

    apply_U<T>("N", m_N, m_k, m_k, m_arnoldi_vec_ptr, m_arnoldi_vec_ld, m_tau_ptr,
            ptr_V, V_ld, m_T_ptr, m_T_ld, NT, m_NB, m_work_ptr, m_work_size);

    return matcl::Matrix(V,false);
}

template<class T>
mat_tup_2 arnoldi_blk_impl<T>::get_V_as_reflectors() const
{
    matcl::Matrix V    = m_arnoldi_vec(colon(), colon(1,m_k));
    matcl::Matrix tau  = m_tau(colon(1,m_k));
    return mat_tup_2(V,tau);
};

template<class T>
matcl::Matrix arnoldi_blk_impl<T>::get_H() const
{
    if (m_lanczos == true)
    {        
        value_code vc       = get_value_code();
        matcl::Matrix ret   = make_band_noinit(m_k, m_k, -m_hess_ldiags+1, m_hess_ldiags-1, vc);

        using Mat_B         = raw::Matrix<T,struct_banded>;
        using Mat_D         = raw::Matrix<T,struct_dense>;
        Mat_B mat_ret       = ret.get_impl_unique<Mat_B>();
        const Mat_D& hess   = m_hess.get_impl<Mat_D>();

        T* ptr_ret          = mat_ret.rep_ptr();
        Integer ret_ld      = mat_ret.ld();
        const T* ptr_hess   = hess.ptr();
        Integer hess_ld     = hess.ld();

        matcl_assert(mat_ret.has_diag(0) == true, "matrix must have main diagonal");

        Integer ret_diags   = mat_ret.number_subdiagonals();

        Integer N           = m_k;

        //main diagonal
        T* ptr_ret2         = ptr_ret + mat_ret.first_elem_diag(0);
        const T* ptr_hess2  = ptr_hess + 0;
        for (Integer j = 0; j < N; ++j)
        {
            ptr_ret2[0]     = ptr_hess2[0];
            ptr_ret2        += ret_ld;
            ptr_hess2       += hess_ld;
        }

        //subdiagonals        
        for (Integer i = 1; i <= ret_diags; ++i)
        {
            ptr_ret2        = ptr_ret + mat_ret.first_elem_diag(-i);
            T* ptr_ret3     = ptr_ret + mat_ret.first_elem_diag(i);
            ptr_hess2       = ptr_hess + i;

            for (Integer j = 0; j < N-i; ++j)
            {
                ptr_ret2[0]     = ptr_hess2[0];
                ptr_ret3[0]     = conj(ptr_hess2[0]);
                ptr_ret2        += ret_ld;
                ptr_ret3        += ret_ld;
                ptr_hess2       += hess_ld;
            }
        };

        ret.add_struct(predefined_struct_type::her);
        return ret;
    }
    else
    {
        matcl::Matrix ret   = m_hess(colon(1,m_k), colon(1,m_k));

        if (m_hess_struct)
            ret.set_struct(predefined_struct_type::hessu);
        else
            ret.set_struct(struct_flag());

        return ret;
    };
};

template<class T>
Integer arnoldi_blk_impl<T>::get_H_ldiags() const
{
    return m_hess_ldiags - 1;
};

template<class T>
matcl::Matrix arnoldi_blk_impl<T>::get_resid() const
{
    return m_res(colon(), colon(1,m_KB));
};

template<class T>
Integer arnoldi_blk_impl<T>::current_block_size() const
{
    return m_KB;
};

template<class T>
typename arnoldi_blk_impl<T>::TR 
arnoldi_blk_impl<T>::get_norm_resid() const
{
    return m_norm_res;
};

template<class T>
const linear_operator& arnoldi_blk_impl<T>::get_operator() const
{
    return m_A;
};

template<class T>
value_code arnoldi_blk_impl<T>::get_value_code() const
{
    return matrix_traits::value_code<T>::value;
}

template class arnoldi_blk_impl<Real>;
template class arnoldi_blk_impl<Float>;
template class arnoldi_blk_impl<Complex>;
template class arnoldi_blk_impl<Float_complex>;

}};