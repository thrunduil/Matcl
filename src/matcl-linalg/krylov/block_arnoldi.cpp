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

#include "matcl-linalg/krylov/block_arnoldi.h"
#include "matcl-core/matrix/matrix_traits.h"
#include "matcl-matrep/lib_functions/matrix_gen.h"
#include "matcl-matrep/details/extract_type_switch.h" 
#include "arnoldi_blk_impl.h"
#include "matcl-linalg/decompositions/householder.h"

#include "matcl-internals/container/mat_s.h"
#include "matcl-internals/container/mat_d.h"
#include "matcl-internals/container/mat_b.h"
#include "matcl-internals/func/converter.h"

namespace matcl { namespace details
{

class block_arnoldi_impl
{
    public:
        using impl_type = std::shared_ptr<details::block_arnoldi_impl>;

    public:
        static void set_impl(block_arnoldi_iteration& ret, const impl_type& impl)
        {
            ret.m_impl = impl;
        };

        virtual ~block_arnoldi_impl(){};

        virtual Integer         run(const matcl::Matrix& v, Integer k, Real tol) = 0;

        virtual Integer         continue_run(Integer k, Real tol) = 0;
        virtual Integer         continue_run(const matcl::Matrix& v, Integer k, Real tol) = 0;
        virtual void            resize(Integer max_k, Integer max_kb) = 0;
        virtual void            clear() = 0;

        virtual Integer         max_k() const = 0;
        virtual Integer         max_kb() const = 0;
        virtual Integer         number_vectors() const = 0;
        virtual Integer         current_block_size() const = 0;
        virtual Integer         get_H_ldiags() const = 0;
        virtual const linear_operator&
                                get_operator() const = 0;
        virtual matcl::Matrix   get_V() const = 0;
        virtual mat_tup_2       get_V_as_reflectors() const = 0;
        virtual unitary_matrix  get_V_as_unitary_matrix() const = 0;
        virtual matcl::Matrix   get_H() const = 0;
        virtual matcl::Matrix   get_resid() const = 0;
        virtual Real            get_norm_resid() const = 0;
};

template<class Val>
class block_arnoldi_impl_val : public block_arnoldi_impl
{
    private:
        using raw_impl_type = matcl::raw::arnoldi_blk_impl<Val>;

    private:
        raw_impl_type   m_impl;
        bool            m_has_run;

    public:
        block_arnoldi_impl_val(const linear_operator& A, Integer max_k, Integer max_kb)
            :m_impl(A, get_lanczos(A), max_k, max_kb), m_has_run(false)
        {};

        virtual ~block_arnoldi_impl_val() override{};

        static void construct(const linear_operator& A, block_arnoldi_iteration& ret,
                              Integer max_k, Integer max_kb)
        {            
            impl_type impl(new block_arnoldi_impl_val(A, max_k, max_kb));

            set_impl(ret, impl);
        };

        virtual Integer run(const matcl::Matrix& v, Integer np, Real tol) override
        {
            Integer M   = v.rows();
            Integer N   = v.cols();

            if (m_impl.size() != M || N < 1)
                throw error::invalid_size2(M, N, m_impl.size(), std::max(N,1));

            Integer max_k   = std::max(np, m_impl.max_k());            
            Integer max_kb  = std::max(N, m_impl.max_kb());
            max_k           = std::min(max_k, M);
            max_kb          = std::min(max_kb, M);
            Integer k       = std::min(np, M);
            k               = std::max(k,0);

            m_impl.resize(max_k, max_kb);

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
            Integer max_kb  = m_impl.max_kb();

            max_k           = std::min(max_k, M);
            max_kb          = std::min(max_kb, M);

            m_impl.resize(max_k, max_kb);

            Integer ret_k   = m_impl.continue_run(std::max(0,np), tol);
            return ret_k;
        };

        virtual Integer continue_run(const matcl::Matrix& v, Integer np, Real tol) override
        {
            if (m_has_run == false)
                return run(v,np,tol);

            Integer M   = v.rows();
            Integer N   = v.cols();

            if (m_impl.size() != M || N < 1)
                throw error::invalid_size2(M, N, m_impl.size(), std::max(N,1));

            Integer k       = m_impl.number_vectors();
            Integer max_k   = std::max(k+np, m_impl.max_k());
            Integer max_kb  = std::max(N, m_impl.max_kb());

            max_k           = std::min(max_k, M);
            max_kb          = std::min(max_kb, M);

            m_impl.resize(max_k, max_kb);

            Integer ret_k   = m_impl.continue_run(v, std::max(0,np), tol);
            return ret_k;
        };

        virtual void resize(Integer max_k, Integer max_kb) override
        {
            m_impl.resize(max_k, max_kb);
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
        virtual Integer max_kb() const override
        {
            return m_impl.max_kb();
        };
        virtual Integer number_vectors() const override
        {
            return m_impl.number_vectors();
        };
        virtual Integer current_block_size() const override
        {
            return m_impl.current_block_size();
        };
        virtual Integer get_H_ldiags() const override
        {
            return m_impl.get_H_ldiags();
        };
        virtual const linear_operator& get_operator() const override
        {
            return m_impl.get_operator();
        };
        virtual matcl::Matrix get_V() const override
        {
            return m_impl.get_V();
        };
        virtual mat_tup_2 get_V_as_reflectors() const override
        {
            return m_impl.get_V_as_reflectors();
        };
        virtual unitary_matrix  get_V_as_unitary_matrix() const override
        {
            Matrix V, tau;
            tie(V,tau) = get_V_as_reflectors();

            return householder_to_unitary(V, tau, V.cols());
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
struct block_arnoldi_impl_val<Integer>
{
    static void construct(const linear_operator& A, block_arnoldi_iteration& ret,
                          Integer max_k, Integer max_kb)
    {
        linear_operator Ac  = A.convert(value_code::v_real);
        return block_arnoldi_impl_val<Real>::construct(Ac,ret,max_k,max_kb);
    };
};

template<>
struct block_arnoldi_impl_val<Object>
{
    static void construct(const linear_operator&, block_arnoldi_iteration&,
                          Integer, Integer)
    {
        throw error::object_value_type_not_allowed("block_arnoldi_iteration");
    };
};

struct vis_block_arnoldi_cons : public extract_type_switch<void,vis_block_arnoldi_cons,true>
{
    template<class T>
    static void eval(const Matrix&, const T&, const linear_operator& A, block_arnoldi_iteration& ret,
                     Integer max_k, Integer max_kb)
    {
        using V = typename T::value_type;
        return block_arnoldi_impl_val<V>::construct(A, ret, max_k, max_kb);
    };

    template<class T>
    static void eval_scalar(const Matrix&, const T&, const linear_operator& A, block_arnoldi_iteration& ret,
                            Integer max_k, Integer max_kb)
    {
        return block_arnoldi_impl_val<T>::construct(A, ret, max_k, max_kb);
    };
};

}};

namespace matcl
{

block_arnoldi_iteration::block_arnoldi_iteration()    
{
    initialize(Matrix(),1,1);
};

block_arnoldi_iteration::block_arnoldi_iteration(const linear_operator& A, Integer init_max_k, 
                                                 Integer init_max_kb)
{
    initialize(A,init_max_k,init_max_kb);
};

void block_arnoldi_iteration::initialize(const linear_operator &A,Integer max_k, Integer max_kb)
{
    if (A.rows() != A.cols())
        throw error::square_matrix_required(A.rows(), A.cols());

    value_code vc_A = A.get_value_code();
    value_code vc   = matrix_traits::unify_value_types(vc_A, value_code::v_float);

    linear_operator Ac  = (vc == vc_A)? A : A.convert(vc);

    Matrix dum      = zeros(0,0,vc);
    max_k           = std::min(max_k, A.rows());
    max_kb          = std::min(max_kb, A.rows());
    return details::vis_block_arnoldi_cons::make<const Matrix&>(dum, Ac, *this, max_k, max_kb);
};

block_arnoldi_iteration::~block_arnoldi_iteration()
{};

block_arnoldi_iteration& block_arnoldi_iteration::operator()(const linear_operator &A, Integer init_max_k, 
                                Integer init_max_kb)
{
    initialize(A,init_max_k,init_max_kb);
    return *this;
};

Integer block_arnoldi_iteration::run(const matcl::Matrix& v, Integer k, Real tol)
{
    return m_impl->run(v,k,tol);
};

Integer block_arnoldi_iteration::continue_run(Integer k, Real tol)
{
    return m_impl->continue_run(k,tol);
};
Integer block_arnoldi_iteration::continue_run(const matcl::Matrix& v, Integer k, Real tol)
{
    return m_impl->continue_run(v,k,tol);
};

void block_arnoldi_iteration::resize(Integer max_k, Integer max_kb)
{
    m_impl->resize(max_k, max_kb);
};
void block_arnoldi_iteration::clear()
{
    m_impl->clear();
};

Integer block_arnoldi_iteration::max_k() const
{
    return m_impl->max_k();
};
Integer block_arnoldi_iteration::max_kb() const
{
    return m_impl->max_kb();
};

Integer block_arnoldi_iteration::number_vectors() const
{
    return m_impl->number_vectors();
};
Integer block_arnoldi_iteration::current_block_size() const
{
    return m_impl->current_block_size();
};

Integer block_arnoldi_iteration::get_H_ldiags() const
{
    return m_impl->get_H_ldiags();
};

const linear_operator& block_arnoldi_iteration::get_operator() const
{
    return m_impl->get_operator();
};

matcl::Matrix block_arnoldi_iteration::get_V() const
{
    return m_impl->get_V();
};

mat_tup_2 block_arnoldi_iteration::get_V_as_reflectors() const
{
    return m_impl->get_V_as_reflectors();
};

matcl::Matrix block_arnoldi_iteration::get_H() const
{
    return m_impl->get_H();
};

matcl::Matrix block_arnoldi_iteration::get_resid() const
{
    return m_impl->get_resid();
};

Real block_arnoldi_iteration::get_norm_resid() const
{
    return m_impl->get_norm_resid();
};

};