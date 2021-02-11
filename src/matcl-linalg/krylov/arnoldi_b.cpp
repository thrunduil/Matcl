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

#include "matcl-linalg/krylov/arnoldi.h"
#include "matcl-core/matrix/matrix_traits.h"
#include "matcl-matrep/lib_functions/matrix_gen.h"
#include "matcl-matrep/details/extract_type_switch.h" 
#include "arnoldi_impl.h"

#include "matcl-internals/container/mat_s.h"
#include "matcl-internals/container/mat_d.h"
#include "matcl-internals/container/mat_b.h"
#include "matcl-internals/func/converter.h"

namespace matcl { namespace details
{

class arnoldi_b_impl
{
    public:
        using impl_type = std::shared_ptr<details::arnoldi_b_impl>;

    public:
        static void set_impl(arnoldi_b_iteration& ret, const impl_type& impl)
        {
            ret.m_impl = impl;
        };

        virtual ~arnoldi_b_impl(){};

        virtual Integer         run(const matcl::Matrix& v, Integer k, Real tol) = 0;

        virtual Integer         continue_run(Integer k, Real tol) = 0;
        virtual Integer         continue_run(const matcl::Matrix& v, Integer k, Real tol) = 0;
        virtual void            resize(Integer max_k) = 0;
        virtual void            clear() = 0;

        virtual Integer         max_k() const = 0;
        virtual Integer         number_vectors() const = 0;
        virtual const linear_operator&
                                get_operator() const = 0;
        virtual const linear_operator&
                                get_operator_B() const = 0;
        virtual matcl::Matrix   get_V() const = 0;
        virtual matcl::Matrix   get_H() const = 0;
        virtual matcl::Matrix   get_resid() const = 0;
        virtual Real            get_norm_resid() const = 0;
};

template<class Val>
class arnoldi_b_impl_val : public arnoldi_b_impl
{
    private:
        using raw_impl_type = matcl::raw::arnoldi_impl<Val>;

    private:
        raw_impl_type   m_impl;
        bool            m_has_run;

    public:
        arnoldi_b_impl_val(const linear_operator& A, const linear_operator& B, bool lanczos, Integer max_k)
            :m_impl(A, B, lanczos, max_k), m_has_run(false)
        {};

        virtual ~arnoldi_b_impl_val() override{};

        static void construct(const linear_operator& A, const linear_operator& B, arnoldi_b_iteration& ret,
                              bool lanczos, Integer max_k)
        {            
            impl_type impl(new arnoldi_b_impl_val(A, B, lanczos, max_k));

            set_impl(ret, impl);
        };

        virtual Integer run(const matcl::Matrix& v, Integer np, Real tol) override
        {
            Integer M   = v.rows();
            Integer N   = v.cols();

            if (m_impl.size() != M || N != 1)
                throw error::invalid_size2(M, N, m_impl.size(), 1);

            Integer max_k   = std::max(np, m_impl.max_k());            
            max_k           = std::min(max_k, M);
            Integer k       = std::min(np, M);
            k               = std::max(k,0);

            m_impl.resize(max_k);

            Integer ret_k   = m_impl.run(v, k, tol);
            m_has_run       = true;

            return ret_k;
        };

        virtual Integer continue_run(Integer np, Real tol) override
        {
            if (m_has_run == false)
                throw error::unable_to_continue_arnoldi();

            Integer M       = m_impl.size();
            Integer k       = m_impl.number_vectors();
            Integer max_k   = std::max(k+np, m_impl.max_k());

            max_k           = std::min(max_k, M);

            m_impl.resize(max_k);

            Integer ret_k   = m_impl.continue_run(std::max(0,np), tol);
            return ret_k;
        };

        virtual Integer continue_run(const matcl::Matrix& v, Integer np, Real tol) override
        {
            if (m_has_run == false)
                return run(v,np,tol);

            Integer M   = v.rows();
            Integer N   = v.cols();

            if (m_impl.size() != M || N != 1)
                throw error::invalid_size2(M, N, m_impl.size(), 1);

            Integer k       = m_impl.number_vectors();
            Integer max_k   = std::max(k+np, m_impl.max_k());

            max_k           = std::min(max_k, M);

            m_impl.resize(max_k);

            Integer ret_k   = m_impl.continue_run(v, std::max(0,np), tol);
            return ret_k;
        };

        virtual void resize(Integer max_k) override
        {
            m_impl.resize(max_k);
        };
        virtual void clear() override
        {
            m_impl.clear();
            m_has_run = false;
        };

        virtual Integer max_k() const override
        {
            return m_impl.max_k();
        };
        virtual Integer number_vectors() const override
        {
            return m_impl.number_vectors();
        };
        virtual const linear_operator& get_operator() const override
        {
            return m_impl.get_operator();
        };
        virtual const linear_operator& get_operator_B() const override
        {
            return m_impl.get_operator_B();
        };

        virtual matcl::Matrix get_V() const override
        {
            return m_impl.get_V();
        };
        virtual matcl::Matrix get_H() const override
        {
            return m_impl.get_H();
        };
        virtual matcl::Matrix get_resid() const override
        {
            return m_impl.get_resid();
        };
        virtual Real get_norm_resid() const override
        {
            return m_impl.get_norm_resid();
        };

    private:
        bool    get_lanczos(const linear_operator& A) const { return A.is_hermitian(); };
};

template<>
struct arnoldi_b_impl_val<Integer>
{
    static void construct(const linear_operator& A, const linear_operator& B, arnoldi_b_iteration& ret,
                          bool lanczos, Integer max_k)
    {
        linear_operator Ac  = A.convert(value_code::v_real);
        linear_operator Bc  = B.convert(value_code::v_real);
        return arnoldi_b_impl_val<Real>::construct(Ac,Bc,ret,lanczos,max_k);
    };
};

template<>
struct arnoldi_b_impl_val<Object>
{
    static void construct(const linear_operator&, const linear_operator&, arnoldi_b_iteration&, bool, Integer)
    {
        throw error::object_value_type_not_allowed("arnoldi_b_iteration");
    };
};

struct vis_arnoldi_cons : public extract_type_switch<void,vis_arnoldi_cons,true>
{
    template<class T>
    static void eval(const Matrix&, const T&, const linear_operator& A, arnoldi_b_iteration& ret,
                     Integer max_k)
    {
        using V = typename T::value_type;
        return arnoldi_b_impl_val<V>::construct(A, ret, max_k);
    };

    template<class T>
    static void eval_scalar(const Matrix&, const T&, const linear_operator& A, arnoldi_b_iteration& ret,
                            Integer max_k)
    {
        return arnoldi_b_impl_val<T>::construct(A, ret, max_k);
    };
};

struct vis_arnoldi_b_cons : public extract_type_switch<void,vis_arnoldi_b_cons,true>
{
    template<class T>
    static void eval(const Matrix&, const T&, const linear_operator& A, const linear_operator& B, 
                     arnoldi_b_iteration& ret, bool lanczos, Integer max_k)
    {
        using V = typename T::value_type;
        return arnoldi_b_impl_val<V>::construct(A, B, ret, lanczos, max_k);
    };

    template<class T>
    static void eval_scalar(const Matrix&, const T&, const linear_operator& A, const linear_operator& B,
                            arnoldi_b_iteration& ret, bool lanczos, Integer max_k)
    {
        return arnoldi_b_impl_val<T>::construct(A, B, ret, lanczos, max_k);
    };
};

}};

namespace matcl
{

arnoldi_b_iteration::arnoldi_b_iteration()    
{
    initialize(Matrix(),Matrix(1.0), false, 1);
};

arnoldi_b_iteration::arnoldi_b_iteration(const linear_operator& A, const linear_operator& B, 
                                         bool lanczos, Integer init_max_k)
{
    initialize(A,B,lanczos, init_max_k);
};

void arnoldi_b_iteration::initialize(const linear_operator &A, const linear_operator& B, bool lanczos, 
                                     Integer max_k)
{
    if (A.rows() != A.cols())
        throw error::square_matrix_required(A.rows(), A.cols());

    if (B.rows() != B.cols())
        throw error::square_matrix_required(B.rows(), B.cols());

    if (A.rows() != B.rows())
        throw error::invalid_size2(B.rows(), B.cols(), A.rows(), A.rows());

    value_code vc_A = A.get_value_code();
    value_code vc_B = B.get_value_code();
    value_code vc0  = matrix_traits::unify_value_types(vc_A, vc_B);
    value_code vc   = matrix_traits::unify_value_types(vc0, value_code::v_float);

    linear_operator Ac  = (vc == vc_A)? A : A.convert(vc);
    linear_operator Bc  = (vc == vc_B)? B : B.convert(vc);

    Matrix dum      = zeros(0,0,vc);
    max_k           = std::min(max_k, A.rows());
    return details::vis_arnoldi_b_cons::make<const Matrix&>(dum, Ac, Bc, *this, lanczos, max_k);
};

arnoldi_b_iteration::~arnoldi_b_iteration()
{};

arnoldi_b_iteration& arnoldi_b_iteration::operator()(const linear_operator &A, const linear_operator &B,
                                                     bool lanczos, Integer init_max_k)
{
    initialize(A,B,lanczos,init_max_k);
    return *this;
};

Integer arnoldi_b_iteration::run(const matcl::Matrix& v, Integer k, Real tol)
{
    return m_impl->run(v,k,tol);
};

Integer arnoldi_b_iteration::continue_run(Integer k, Real tol)
{
    return m_impl->continue_run(k,tol);
};
Integer arnoldi_b_iteration::continue_run(const matcl::Matrix& v, Integer k, Real tol)
{
    return m_impl->continue_run(v,k,tol);
};

void arnoldi_b_iteration::resize(Integer max_k)
{
    m_impl->resize(max_k);
};
void arnoldi_b_iteration::clear()
{
    m_impl->clear();
};

Integer arnoldi_b_iteration::max_k() const
{
    return m_impl->max_k();
};

Integer arnoldi_b_iteration::number_vectors() const
{
    return m_impl->number_vectors();
};

const linear_operator& arnoldi_b_iteration::get_operator() const
{
    return m_impl->get_operator();
};
const linear_operator& arnoldi_b_iteration::get_operator_B() const
{
    return m_impl->get_operator_B();
};

matcl::Matrix arnoldi_b_iteration::get_V() const
{
    return m_impl->get_V();
};

matcl::Matrix arnoldi_b_iteration::get_H() const
{
    return m_impl->get_H();
};

matcl::Matrix arnoldi_b_iteration::get_resid() const
{
    return m_impl->get_resid();
};

Real arnoldi_b_iteration::get_norm_resid() const
{
    return m_impl->get_norm_resid();
};

};