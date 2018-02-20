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

#include "matcl-internals/container/mat_base.h"
#include "matcl-internals/error/error_check_basic.h"
#include "matcl-matrep/base/alloc.h"
#include "matcl-matrep/base/serialize.h"
#include "matcl-matrep/details/debug.h"
#include "matcl-core/details/integer.h"
#include "matcl-scalar/details/scalfunc_helpers.h"

namespace matcl { namespace raw 
{

namespace mrd = matcl::raw::details;

namespace details
{

template<class Val>
root_ptr<Val>::root_ptr()
    :m_ptr(nullptr), m_foreign(false)
{};

template<class Val>
void root_ptr<Val>::destroy(Integer max_rows, Integer max_cols)
{
    if (m_ptr && m_foreign == false)
    {
        matcl::details::allocator<Val>::free(m_ptr, max_rows*max_cols); 
        m_ptr  = nullptr;
    };
};

template<class Val>
void root_ptr<Val>::set(Val* arr)
{
    m_ptr       = arr;
    m_foreign   = false;
};

template<class Val>
void root_ptr<Val>::set_foreign(Val* arr)
{
    m_ptr       = arr;
    m_foreign   = true;
};

template<class Val>
ptrdiff_t root_ptr<Val>::offset(const Val* ptr) const
{
    return ptr - m_ptr;
};

template<class Val>
root_ptr<Val>::root_ptr(const root_ptr& other)
    :m_ptr(other.m_ptr), m_foreign(other.m_foreign)
{};

template<class Val>
root_ptr<Val>& root_ptr<Val>::operator=(const root_ptr& other)
{
    m_ptr       = other.m_ptr;
    m_foreign   = other.m_foreign;
    return *this;
};

template<class val_type>
mat_data<val_type>::mat_data(type_info ti)
    :m_rows(0),m_cols(0),m_max_rows(0),m_max_cols(0),m_ld(1),m_size(0)
    ,m_root_ptr(root_ptr<val_type>()), m_ptr(nullptr)
    ,m_ti(ti),m_refcount(refcount_str::create(1)),m_effective_unique(false)
{};

template<class val_type>
void mat_data<val_type>::mark_unique(bool unique)
{
    m_effective_unique = unique;
};

template<class val_type>
bool mat_data<val_type>::is_unique() const
{
    return m_refcount->get_count() == 1
            || m_effective_unique == true;
};

template<class val_type>
mat_data<val_type>::mat_data(const mat_data& other)
    : m_rows(other.m_rows)
    , m_cols(other.m_cols)
    , m_max_rows(other.m_max_rows)
    , m_max_cols(other.m_max_cols)
    , m_ld(other.m_ld)
    , m_size(other.m_size)
    , m_root_ptr(other.m_root_ptr)
    , m_ptr(other.m_ptr)
    , m_ti(other.m_ti)
    , m_flag(other.m_flag)
    , m_refcount(other.m_refcount)
    , m_effective_unique(false)
{
    m_refcount  = other.m_refcount->increase();
};

template<class val_type>
mat_data<val_type>::mat_data(mat_data&& other)
    : m_rows(other.m_rows)
    , m_cols(other.m_cols)
    , m_max_rows(other.m_max_rows)
    , m_max_cols(other.m_max_cols)
    , m_ld(other.m_ld)
    , m_size(other.m_size)
    , m_root_ptr(other.m_root_ptr)
    , m_ptr(other.m_ptr)
    , m_ti(other.m_ti)
    , m_flag(other.m_flag)
    , m_refcount(other.m_refcount)
    , m_effective_unique(false)
{
    other.m_refcount = nullptr;
};

template<class val_type>
void mat_data<val_type>::reset_unique()
{
    destroy_data();
    m_rows      = m_cols = m_size = m_max_rows = m_max_cols = 0;
    m_ld        = 1;
    m_ptr       = nullptr;
    m_root_ptr.set(nullptr);
    m_flag      = struct_flag();
};

template<class val_type>
void mat_data<val_type>::alloc(Integer n)
{
    if (n > 0)
    {
        val_type* arr   = matcl::details::allocator<val_type>::malloc(m_ti,n);
        error::check_alloc(arr, n, sizeof(val_type));

        m_root_ptr.set(arr);
        m_ptr       = arr;        
    }
    else
    {
        if (n < 0)
            error::check_alloc(nullptr, n, sizeof(val_type));

        m_root_ptr.set(nullptr);
        m_ptr       = nullptr;
    };
};

template<class val_type>
void mat_data<val_type>::take_foreign(val_type* arr)
{
    m_root_ptr.set_foreign(arr);
    m_ptr       = arr;
};

template<class val_type>
mat_data<val_type>::mat_data(type_info ti, Integer r, Integer c)
    :m_rows(r), m_cols(c), m_max_rows(r), m_max_cols(c), m_ti(ti)
{
    error::check_size(r, c);
    m_ld   = std::max(r,1);
    m_size = imult(m_rows,m_cols);

    alloc(m_size);

    m_refcount = refcount_str::create(1);
    m_effective_unique = false;
};

//does not take ownership of arr
template<class val_type>
mat_data<val_type>::mat_data(type_info ti, Integer r, Integer c, val_type* arr, Integer ld)
    :m_rows(r), m_cols(c), m_max_rows(r), m_max_cols(c), m_ti(ti)
{
    error::check_size(r, c);
    m_ld   = std::max(r,ld);
    m_size = imult(m_rows,m_cols);

    take_foreign(arr);

    m_refcount = refcount_str::create(1);
    m_effective_unique = false;
};

template<class val_type>
void mat_data<val_type>::destroy_data()
{
    m_root_ptr.destroy(m_max_rows, m_max_cols);
    m_ptr       = nullptr;
};

template<class val_type>
void mat_data<val_type>::reset_unique(Integer r, Integer c)
{
    destroy_data();
    error::check_size(r, c);

    m_rows      = r;
    m_cols      = c;
    m_max_rows  = r;
    m_max_cols  = c;
    m_ld        = std::max(r,1);
    m_size      = imult(m_rows,m_cols);

    m_flag.reset();
    alloc(m_size);
};

template<class val_type>
mat_data<val_type>& mat_data<val_type>::assign_to_fresh(const mat_data& other)
{
    other.m_refcount->increase();
 
    if (m_refcount->decrease())
    {
        destroy_data();
        m_refcount->destroy();
    };

    m_refcount  = other.m_refcount;
    m_rows      = other.m_rows;
    m_cols      = other.m_cols;
    m_max_rows  = other.m_max_rows;
    m_max_cols  = other.m_max_cols;
    m_ld        = other.m_ld;
    m_size      = other.m_size;
    m_root_ptr  = other.m_root_ptr;
    m_ptr       = other.m_ptr;
    m_ti        = other.m_ti;
    m_flag      = other.m_flag;
    m_effective_unique = false;

    return *this;
};

template<class val_type>
mat_data<val_type>& mat_data<val_type>::assign_to_fresh(mat_data&& other)
{
    std::swap(m_refcount, other.m_refcount);
    std::swap(m_max_cols, other.m_max_cols);
    std::swap(m_max_rows, other.m_max_rows);
    std::swap(m_root_ptr, other.m_root_ptr);

    m_rows      = other.m_rows;
    m_cols      = other.m_cols;
    m_ld        = other.m_ld;
    m_size      = other.m_size;
    m_ptr       = other.m_ptr;
    m_ti        = other.m_ti;
    m_flag      = other.m_flag;
    m_effective_unique = false;

    return *this;
};

template<class val_type>
mat_data<val_type>::~mat_data()
{
    if (m_refcount && m_refcount->decrease())
    {       
        destroy_data();
        m_refcount->destroy();
    }
};

template<class value_type> 
inline void mat_data<value_type>::serialize(oarchive_impl & ar, const unsigned int v) const
{
    Integer m_ld_s = std::max(m_rows,1);
    ar << m_size;
    ar << m_rows;
    ar << m_cols;
    ar << m_ld_s;
    matcl::details::serialize_save(ar,m_ti,v);
    matcl::details::serialize_save(ar,m_flag,v);

    const value_type* loc_ptr = ptr();

    for (Integer j = 0; j < m_cols; ++j)
    {
        for (Integer i = 0; i < m_rows; ++i)
            matcl::details::serialize_save(ar,loc_ptr[i],v);

        loc_ptr += m_ld;
    };
};

template<class value_type> 
void mat_data<value_type>::serialize(iarchive_impl & ar, const unsigned int v)
{  
    ar >> m_size;
    ar >> m_rows;    
    ar >> m_cols;
    ar >> m_ld;
    matcl::details::serialize_load(ar,m_ti,v);
    matcl::details::serialize_load(ar,m_flag,v);

    m_max_rows = m_rows;
    m_max_cols = m_cols;

    error::check_size(m_rows, m_cols);
    alloc(m_size);

    try
    {
        for (Integer i = 0; i < m_size; ++i)
            matcl::details::serialize_load(ar,m_ptr[i],v);
    }
    catch(...)
    {
        destroy_data();
        throw;
    };
};

};

template <class value_type>
inline dense_matrix_base<value_type>::dense_matrix_base(tinfo ti)
:m_data(ti)
{};

template <class value_type>
Integer dense_matrix_base<value_type>::nnz() const
{
    struct_flag::diag_type ld = m_data.get_struct().get_ldiags();
    struct_flag::diag_type ud = m_data.get_struct().get_udiags();

    if (ld == struct_flag::general && ud == struct_flag::general)
        return m_data.m_size;

    Integer nz_d    = std::min(m_data.m_rows,m_data.m_cols);
    Integer nz_l    = 0;
    Integer nz_u    = 0;

    if (nz_d == 0)
        return 0;

    switch (ld)
    {
        case struct_flag::zero:
            break;
        case struct_flag::one:
        case struct_flag::qtriang:
            nz_l        = std::max(std::min(m_data.m_rows-1,m_data.m_cols),0);
            break;
        default:
            if (m_data.m_rows-1 <= m_data.m_cols)
            {
                nz_l    = ((m_data.m_rows-1) * m_data.m_rows) / 2;
            }
            else
            {
                nz_l    = (m_data.m_cols * (m_data.m_cols + 1)) / 2
                        + m_data.m_cols * (m_data.m_rows-1 - m_data.m_cols);
            };
    };

    switch (ud)
    {
        case struct_flag::zero:
            break;
        case struct_flag::one:
        case struct_flag::qtriang:
            nz_u        = std::max(std::min(m_data.m_rows,m_data.m_cols-1),0);
            break;
        default:
            if (m_data.m_cols-1 <= m_data.m_rows)
            {
                nz_u    = ((m_data.m_cols-1) * m_data.m_cols) / 2;
            }
            else
            {
                nz_u    = (m_data.m_rows * (m_data.m_rows + 1)) / 2
                        + m_data.m_rows * (m_data.m_cols-1 - m_data.m_rows);
            };
    };

    return nz_d + nz_l + nz_u;
}

template <class value_type>
dense_matrix_base<value_type>::dense_matrix_base(tinfo ti, Integer r, Integer c)
:m_data(ti,r,c)
{
};

template <class value_type>
dense_matrix_base<value_type>::dense_matrix_base(tinfo ti, Integer r, Integer c, value_type* arr, Integer ld)
:m_data(ti,r,c, arr, ld)
{
};

template <class value_type>
dense_matrix_base<value_type>::dense_matrix_base(const dense_matrix_base &mat) 
: m_data(mat.m_data)
{}

template <class value_type>
dense_matrix_base<value_type>::dense_matrix_base(dense_matrix_base &&mat) 
: m_data(std::move(mat.m_data))
{}

template <class value_type>
dense_matrix_base<value_type>&
dense_matrix_base<value_type>::reset_unique()
{
    m_data.reset_unique();
    return *this;
}

template <class value_type>
dense_matrix_base<value_type>&
dense_matrix_base<value_type>::reset_unique(Integer r, Integer c)
{
    m_data.reset_unique(r,c);
    return *this;
}

template <class value_type>
ti::ti_type<value_type> dense_matrix_base<value_type>::get_type() const
{
    return m_data.get_type();
};

template <class value_type>
dense_matrix_base<value_type>&
dense_matrix_base<value_type>::assign_to_fresh(const dense_matrix_base &mat)
{
    m_data.assign_to_fresh(mat.m_data);
    return *this;
}

template <class value_type>
dense_matrix_base<value_type>&
dense_matrix_base<value_type>::assign_to_fresh(dense_matrix_base &&mat)
{
    m_data.assign_to_fresh(std::move(mat.m_data));
    return *this;
}

template <class value_type>
void dense_matrix_base<value_type>::read_from(const value_type* ptr, Integer ld)
{
    if (m_data.m_size == 0)
        return;

    value_type* this_ptr = m_data.ptr();
    Integer this_ld     = this->ld();

    for(Integer i = 0; i < m_data.m_cols; ++i)
    {
        for (Integer j = 0; j < m_data.m_rows; ++j)
            mrd::reset_helper(this_ptr[j], ptr[j]);

        this_ptr    += this_ld;
        ptr         += ld;
    };
};

template<class value_type> 
inline void dense_matrix_base<value_type>::serialize(oarchive_impl & ar, const unsigned int ) const
{
    ar << m_data;
};

template<class value_type> 
void dense_matrix_base<value_type>::serialize(iarchive_impl & ar, const unsigned int )
{  
    ar >> m_data;
};

};};

template class matcl::raw::details::mat_data<matcl::Integer>;
template class matcl::raw::details::mat_data<matcl::Real>;
template class matcl::raw::details::mat_data<matcl::Complex>;
template class matcl::raw::details::mat_data<matcl::Object>;
template class matcl::raw::details::mat_data<matcl::Float>;
template class matcl::raw::details::mat_data<matcl::Float_complex>;

template class matcl::raw::details::root_ptr<matcl::Integer>;
template class matcl::raw::details::root_ptr<matcl::Real>;
template class matcl::raw::details::root_ptr<matcl::Complex>;
template class matcl::raw::details::root_ptr<matcl::Object>;
template class matcl::raw::details::root_ptr<matcl::Float>;
template class matcl::raw::details::root_ptr<matcl::Float_complex>;

template class matcl::raw::dense_matrix_base<matcl::Integer>;
template class matcl::raw::dense_matrix_base<matcl::Real>;
template class matcl::raw::dense_matrix_base<matcl::Complex>;
template class matcl::raw::dense_matrix_base<matcl::Object>;
template class matcl::raw::dense_matrix_base<matcl::Float>;
template class matcl::raw::dense_matrix_base<matcl::Float_complex>;
