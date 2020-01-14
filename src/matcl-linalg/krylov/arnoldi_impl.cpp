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

#include "arnoldi_impl.h"
#include "matcl-matrep/matcl_matrep.h"
#include "matcl-linalg/norms_error/norm.h"
#include "arnoldi_process.h"
#include "matcl-blas-ext/lapack_ext/lapack_ext.h"
#include "matcl-internals/func/lapack_utils.h"

namespace matcl { namespace raw
{

template<class T>
arnoldi_impl<T>::arnoldi_impl(const linear_operator& A, bool lanczos, Integer max_k)    
{
    init(A, lanczos, max_k);
}
template<class T>
arnoldi_impl<T>::arnoldi_impl(const linear_operator& A, const linear_operator& B,
                              bool lanczos, Integer max_k)
{
    init(A, lanczos, max_k);
    m_B                 = B;
    m_has_B             = true;
};

template<class T>
void arnoldi_impl<T>::init(const linear_operator& A, bool lanczos, Integer max_k)
{
    m_has_B             = false;
    m_lanczos           = lanczos;
    m_A                 = A;
    m_N                 = A.rows();
    m_k                 = 0;
    m_max_k             = std::max(0,max_k);

    m_resid_norm        = TR(0.0);
    m_scale             = TR(0.0);
    m_sumsq             = TR(1.0);
    m_arnoldi_vec_ptr   = nullptr;
    m_arnoldi_vec_ld    = 1;

    m_hess_ld           = 1;
    m_work_ptr          = nullptr;

    init_arrays();
};

template<class T>
void arnoldi_impl<T>::init_arrays()
{
    if (m_lanczos)
    {
        value_code vcr  = matrix_traits::value_code<TR>::value;
        m_hess          = zeros(m_max_k, 2, vcr);        
        m_hess_ld       = m_max_k;
    }
    else
    {
        m_hess          = zeros(m_max_k, m_max_k, get_value_code());
        m_hess_ld       = m_max_k;
    };

    m_work              = make_dense_noinit(m_N, 4, get_value_code());
    m_work_ptr          = m_work.get_array_unique<T>();

    m_arnoldi_vec       = make_dense_noinit(m_N, m_max_k, get_value_code());
    m_arnoldi_vec_ptr   = m_arnoldi_vec.get_array_unique<T>();
    m_arnoldi_vec_ld    = m_N;
};

template<class T>
void arnoldi_impl<T>::clear()
{
    m_k                 = 0;
};

template<class T>
void arnoldi_impl<T>::resize(Integer max_k)
{
    max_k               = std::max(max_k, 0);

    if (m_lanczos)
    {
        m_hess.resize(max_k, 2);        
        m_hess_ld       = convert(m_hess, Mat_DR::matrix_code).get_impl<Mat_DR>().ld();
    }
    else
    {
        m_hess.resize(max_k, max_k);
        m_hess_ld       = convert(m_hess, Mat_D::matrix_code).get_impl<Mat_D>().ld();
    };

    m_arnoldi_vec.resize(m_N, max_k);
    m_arnoldi_vec_ptr   = m_arnoldi_vec.get_array_unique<T>();
    m_arnoldi_vec_ld    = convert(m_arnoldi_vec, Mat_D::matrix_code).get_impl<Mat_D>().ld();

    m_k                 = std::min(m_k, max_k);
    m_max_k             = max_k;
};

template<class T>
void arnoldi_impl<T>::init_work_v(const matcl::Matrix& v)
{
    using md::lap;

    const T* v_ptr      = v.get_array<T>();
    lapack::copy(m_N, lap(v_ptr), 1, lap(m_work_ptr), 1);

    m_resid_norm        = TR(0.0);
};

template<class T>
arnoldi_impl<T>::~arnoldi_impl()
{};

template<class T>
Integer arnoldi_impl<T>::run(const matcl::Matrix& v, Integer k, Real tol)
{
    init_work_v(v);
    return run_impl(0, std::min(k, m_max_k), tol);
};

template<class T>
Integer arnoldi_impl<T>::continue_run(Integer k, Real tol)
{
    Integer K       = std::min(m_k + k, m_max_k);
    m_resid_norm    = TR(0.0);
    return run_impl(m_k, K, tol);
};

template<class T>
Integer arnoldi_impl<T>::continue_run(const matcl::Matrix& v, Integer k, Real tol)
{
    init_work_v(v);
    bool suc = orthogonalize_resid();
    
    if (suc == false)
        return 0;

    Integer K       = std::min(m_k + k, m_max_k);    

    return run_impl(m_k, K, tol);
};

template<class T>
Integer arnoldi_impl<T>::run_impl(Integer k, Integer K, Real tol)
{    
    Integer info    = 0;

    struct callback_impl : callback_aitr<T>
    {    
        arnoldi_impl* owner;

        callback_impl(arnoldi_impl* o)
            :owner(o)
        {};

        virtual void eval(Integer type, Integer N, T* X, T* Y, T* BX) override
        {
            (void)BX;

            if (type == 1)
            {
                //eval Y = OP * X 
                matcl::Matrix x     = make_dense_foreign(N, 1, X, N);
                matcl::Matrix y     = make_dense_foreign(N, 1, Y, N);
                owner->m_A.mmul_right(x, trans_type::no_trans, y);
            }
            else if (type == 2)
            {
                //eval Y = B * X
                matcl::Matrix x     = make_dense_foreign(N, 1, X, N);                
                matcl::Matrix y     = make_dense_foreign(N, 1, Y, N);
                owner->m_B.mmul_right(x, trans_type::no_trans, y);
            };
        }
    };

    callback_impl callback(this);

    Integer np          = std::max(K - k, 0);

    TR* H_S_ptr         = (m_lanczos == true) ? m_hess.get_array_unique<TR>() : nullptr;
    T* H_NS_ptr         = (m_lanczos == false) ? m_hess.get_array_unique<T>() : nullptr;

    aitr<T>(m_lanczos, m_has_B, m_N, k, np, (TR)tol, m_resid_norm, m_arnoldi_vec_ptr,
                    m_arnoldi_vec_ld, H_NS_ptr, H_S_ptr, m_hess_ld, m_scale, m_sumsq, 
                    m_work_ptr, &callback, info);

    if (info < 0)
    {
        std::ostringstream msg;
        msg << "aitr returned code: " << info;
        throw error::error_general(msg.str());
    }        

    Integer k_new       = info > 0 ? info : K;
    Integer n_steps     = k_new - m_k;

    m_k                 = k_new;

    return n_steps;
};

template<class T>
void arnoldi_impl<T>::assign_to_work(Integer start, Integer N, const matcl::Matrix& val)
{
    const T* ptr    = val.get_array<T>();
    T* dest         = m_work_ptr + start - 1;

    for (Integer i = 0; i < N; ++i)
        dest[i]     = ptr[i];
};

template<class T>
bool arnoldi_impl<T>::orthogonalize_resid()
{
    using md::lap;

    Integer POS_RES = 0;
    Integer POS_W1  = m_N;
    Integer POS_W2  = m_N*2;
    Integer POS_B   = m_N*3;

    TR RNORM;
    Integer INFO2;
    
    if (m_has_B)
    {
        struct callback
        {
            using TL            = typename md::lapack_value_type<T>::type;

            arnoldi_impl*   m_arnoldi;

            static void eval_gsg(void* CTX, Integer N, TL* X, TL* Y)
            { 
                callback* owner     = reinterpret_cast<callback*>(CTX);
                owner->eval_impl(N, (T*)X, (T*)Y);
            };

            void eval_impl(Integer N, T* X, T* Y)
            {
                //eval Y = B * X
                matcl::Matrix x     = make_dense_foreign(N, 1, X, N);
                matcl::Matrix y     = make_dense_foreign(N, 1, Y, N);
                m_arnoldi->m_B.mmul_right(x, trans_type::no_trans, y);
            };
        };

        callback ctx{this};

        lapack::gsg_ort(&callback::eval_gsg, &ctx, "N", m_N, m_k, 
            lap(m_arnoldi_vec_ptr), m_arnoldi_vec_ld, lap(m_work_ptr + POS_RES), 1,
            lap(m_work_ptr + POS_W1), 1, RNORM, lap(m_work_ptr + POS_B), 
            lap(m_work_ptr + POS_W2), INFO2);
    }
    else
    {
        lapack::gs_ort("N", m_N, m_k, lap(m_arnoldi_vec_ptr), m_arnoldi_vec_ld, 
                lap(m_work_ptr + POS_RES), 1, lap(m_work_ptr + POS_W1), 1, RNORM, 
                lap(m_work_ptr + POS_W2), INFO2);
    };

    if (INFO2 == 0)
        return true;
    else
        return false;
};

template<class T>
Integer arnoldi_impl<T>::size() const
{
    return m_N;
};

template<class T>
const linear_operator& arnoldi_impl<T>::get_operator() const
{
    return m_A;
};
template<class T>
const linear_operator& arnoldi_impl<T>::get_operator_B() const
{
    return m_B;
};

template<class T>
Integer arnoldi_impl<T>::max_k() const
{
    return m_max_k;
};

template<class T>
Integer arnoldi_impl<T>::number_vectors() const
{
    return m_k;
};

template<class T>
matcl::Matrix arnoldi_impl<T>::get_V() const
{
    matcl::Matrix ret = m_arnoldi_vec(colon(), colon(1,m_k));
    return ret.clone();
}

template<class T>
matcl::Matrix arnoldi_impl<T>::get_H() const
{
    if (m_lanczos == true)
    {
        matcl::Matrix d     = (mat_row(), 0, 1, -1);

        value_code vcr      = matrix_traits::value_code<TR>::value;
        matcl::Matrix ret   = make_band_noinit(m_k, m_k, -1, 1, vcr);

        if (m_k == 0)
            return ret;

        ret.diag(0)         = m_hess(colon(1,m_k), 2);

        if (m_k > 1)
        {
            ret.diag(1)         = m_hess(colon(2,m_k), 1);
            ret.diag(-1)        = m_hess(colon(2,m_k), 1);
        };

        ret.add_struct(predefined_struct_type::her);
        return ret;
    }
    else
    {
        matcl::Matrix ret   = m_hess(colon(1,m_k), colon(1,m_k));

        ret.set_struct(predefined_struct_type::hessu);
        return ret.clone();
    };
};

template<class T>
matcl::Matrix arnoldi_impl<T>::get_resid() const
{
    return m_work(colon(), 1).clone();
};

template<class T>
typename arnoldi_impl<T>::TR 
arnoldi_impl<T>::get_norm_resid() const
{
    return m_resid_norm;
};

template<class T>
value_code arnoldi_impl<T>::get_value_code() const
{
    return matrix_traits::value_code<T>::value;
}

template class arnoldi_impl<Real>;
template class arnoldi_impl<Float>;
template class arnoldi_impl<Complex>;
template class arnoldi_impl<Float_complex>;

}};