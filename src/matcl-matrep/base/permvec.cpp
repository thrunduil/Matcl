/*
 *  This file is a part of Matrix Computation Library (MATCL)
 *
 *  Copyright (c) Pawe³ Kowal 2017 - 2018
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

#include "matcl-matrep/matrix/permvec.h"
#include "matcl-matrep/lib_functions/matrix_gen.h"
#include "matcl-matrep/lib_functions/manip.h"
#include "matcl-matrep/matrix/colon.h"
#include "matcl-matrep/details/matrix_details_subs.h"
#include "matcl-internals/error/error_check_basic.h"
#include "matcl-core/utils/workspace.h"

namespace matcl
{

//---------------------------------------------------------------------------
//                          pv_impl
//---------------------------------------------------------------------------
class pv_impl
{
    private:
        using abool         = matcl::atomic<bool>;

    public:        
        ~pv_impl();

        Integer             length() const;
        bool                is_id() const;       
		
        Matrix              to_matrix() const;
        Matrix              to_interchanges_matrix() const;
        const Integer*      to_array() const;

		permvec             invperm(const permvec& owner) const;
        permvec             apply(const pv_impl& p) const;

        static permvec      identity(Integer length);
        static permvec      from_matrix(const Matrix& p, bool check);

    private:
        pv_impl(const Matrix& p, const Matrix& pinv);
        pv_impl(const Matrix& p);
        pv_impl(Integer length);        //make identity

        void                set_int_matrix(const Matrix& p_int) const;
        void                set_inv_matrix(const Matrix& p_inv) const;
        void                set_id_matrix(const Matrix& p) const;

    private:
        Integer             m_length;
        bool                m_is_identity;

        mutable Matrix      m_p;
        mutable Matrix      m_pinv;
        mutable Matrix      m_pint;

        mutable abool       m_has_p;
        mutable abool       m_has_inv;
        mutable abool       m_has_int;        
};

pv_impl::pv_impl(const Matrix& p, const Matrix& pinv)
    :m_has_inv(1), m_has_int(0), m_has_p(1), m_is_identity(false), m_length(p.length())
    , m_p(p), m_pinv(pinv)
{};

pv_impl::pv_impl(const Matrix& p)
    :m_has_inv(0), m_has_int(0), m_has_p(1), m_is_identity(false), m_length(p.length())
    , m_p(p)
{};

pv_impl::pv_impl(Integer length)
    :m_has_inv(0), m_has_int(0), m_has_p(0), m_is_identity(true), m_length(length)
{};

pv_impl::~pv_impl()
{};

Integer pv_impl::length() const
{
    return m_length;
};

bool pv_impl::is_id() const
{
    return m_is_identity;
};

Matrix pv_impl::to_matrix() const
{
    if (this->is_id() == false)
        return vec(m_p);

    if (this->m_has_p != 0)
    {
        std::atomic_thread_fence(std::memory_order_acquire);

        return vec(m_p);
    };

    Matrix p = matcl::irange(1,this->m_length);
    set_id_matrix(p);
    
    return vec(p);
};

//taken from blas
static void perm2int(Integer M, Integer *IPIV, Integer* WORK, bool ZB)
{
    if (M<=0)
        return;
    
    if (ZB == true)
    {
        for(Integer i=0; i<M;i++)
            WORK[IPIV[i]]	= i;

        for (Integer i=0;i<M;i++)
        {
            IPIV[WORK[i]]	= IPIV[i];
            WORK[IPIV[i]]	= WORK[i];
        };
    }
    else
    {
        for(Integer i = 0; i < M; i++)
            WORK[IPIV[i] - 1]	= i + 1;

        for (Integer i = 0; i < M; i++)
        {
            IPIV[WORK[i] - 1]	= IPIV[i];
            WORK[IPIV[i] - 1]	= WORK[i];
        };
    }
};

Matrix pv_impl::to_interchanges_matrix() const
{
    if (m_has_int != 0)
    {
        std::atomic_thread_fence(std::memory_order_acquire);
        return m_pint;
    }

    if (this->is_id() == true)
    {
        Matrix p = matcl::irange(1,this->m_length);
        set_id_matrix(p);
    
        return p;
    };

    Matrix p_int                        = m_p;
    matcl::pod_workspace<Integer> work(this->m_length,1);

    Integer* p_ptr  = p_int.get_array_unique<Integer>();
    Integer* w_ptr  = work.ptr();

    perm2int(this->m_length, p_ptr, w_ptr, false);

    set_int_matrix(p_int);
    return m_pint;
};

const Integer* pv_impl::to_array() const
{
    if (this->is_id() == false)
        return m_p.get_array<Integer>();

    if (this->m_has_p != 0)
    {
        std::atomic_thread_fence(std::memory_order_acquire);

        return m_p.get_array<Integer>();
    };

    Matrix p = matcl::irange(1,this->m_length);
    set_id_matrix(p);
    
    return m_p.get_array<Integer>();
};

permvec pv_impl::invperm(const permvec& owner) const
{
    if (this->m_is_identity == true)
        return owner;

    if (this->m_has_inv == 0)
    {
        Matrix pinv     = izeros(1,this->m_length);
        pinv(this->m_p) = irange(1,this->m_length);

        set_inv_matrix(pinv);
    };

    std::atomic_thread_fence(std::memory_order_acquire);

    std::shared_ptr<pv_impl> impl(new pv_impl(this->m_pinv, this->m_p));
    return permvec(impl);
};

permvec pv_impl::apply(const pv_impl& p) const
{
    //this and p is not identity

    Matrix out = this->m_p(matcl::colon(p.m_p));
    return pv_impl::from_matrix(out,false);
}

permvec pv_impl::identity(Integer length)
{
    std::shared_ptr<pv_impl> impl(new pv_impl(length));
    return permvec(impl);
};

permvec pv_impl::from_matrix(const Matrix& p, bool check)
{
    if (check == true)
    {
        bool trivial = error::check_permutation_vector(p, p.length());
        if (trivial == true)
        {
            std::shared_ptr<pv_impl> impl(new pv_impl(p.length()));
            return permvec(impl);
        };
    };

    std::shared_ptr<pv_impl> impl(new pv_impl(p));
    return permvec(impl);
};

void pv_impl::set_int_matrix(const Matrix& p_int) const
{
    this->m_pint = p_int;

    std::atomic_thread_fence(std::memory_order_release);
    this->m_has_int = true;
};

void pv_impl::set_inv_matrix(const Matrix& p_inv) const
{
    this->m_pinv = p_inv;

    std::atomic_thread_fence(std::memory_order_release);
    this->m_has_inv = true;
};

void pv_impl::set_id_matrix(const Matrix& p) const
{
    this->m_p       = p;
    this->m_pinv    = p;
    this->m_pint    = p;

    std::atomic_thread_fence(std::memory_order_release);

    this->m_has_inv = true;
    this->m_has_int = true;
    this->m_has_p = true;
};

//---------------------------------------------------------------------------
//                          permvec
//---------------------------------------------------------------------------
permvec::permvec()
{};

permvec::permvec(const impl_type& impl)
    :m_data(impl)
{};

permvec::~permvec()
{};

Integer permvec::length() const
{
    if (!m_data)
        return 0;
    else
        return m_data->length();
};

bool permvec::is_id() const
{
    if (!m_data)
        return true;
    else
        return m_data->is_id();
};

permvec permvec::invperm() const
{
    if (!m_data)
        return *this;
    else
        return m_data->invperm(*this);
};

permvec permvec::operator()(const permvec& p) const
{
    if (!m_data)
    {
        if (p.length() != 0)
            throw error::invalid_perm_vectors_composition(0,p.length());
        else
            return *this;
    };

    if (p.length() == 0)
        throw error::invalid_perm_vectors_composition(this->length(),p.length());

    if (p.is_id() == true)
        return *this;

    if (this->is_id() == true)
        return p;

    return this->m_data->apply(*p.m_data);
};

Matrix permvec::to_matrix() const
{
    if (!m_data)
        return matcl::izeros(0,1);
    else
        return m_data->to_matrix();
};

Matrix permvec::to_interchanges_matrix() const
{
    if (!m_data)
        return matcl::izeros(0,1);
    else
        return m_data->to_interchanges_matrix();
};

const Integer* permvec::to_array() const
{
    if (!m_data)
        return nullptr;
    else
        return m_data->to_array();
};

permvec permvec::identity(Integer length)
{
    if (length == 0)
        return permvec();

    return pv_impl::identity(length);
};

permvec permvec::from_matrix(const Matrix& p)
{
    if (p.length() == 0)
        return permvec();

    return pv_impl::from_matrix(p, true);
};

permvec permvec::from_matrix_nocheck(const Matrix& p)
{
    if (p.length() == 0)
        return permvec();

    return pv_impl::from_matrix(p, false);
};

};
