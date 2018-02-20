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

#include "matcl-internals/container/sparse_ccs.h"
#include "matcl-internals/error/error_check_basic.h"
#include "matcl-matrep/base/alloc.h"
#include "matcl-internals/base/sort.h"
#include "matcl-internals/base/utils.h"
#include "matcl-matrep/base/serialize.h"
#include "matcl-scalar/details/matfunc_helpers.h"
#include "matcl-matrep/utils/workspace.h"
#include "matcl-scalar/details/scalfunc_helpers.h"

namespace matcl { namespace raw { namespace details
{

namespace md = matcl::details;
namespace mrd = matcl::raw::details;

template <class value_type>
sparse_ccs<value_type>::~sparse_ccs()
{
    if (m_refcount && m_refcount->decrease())
    {        
        destroy_data();
        m_refcount->destroy();
    };
}

template <class value_type>
void sparse_ccs<value_type>::mark_unique(bool is_unique)
{
    m_effective_unique = is_unique;
};

template <class value_type>
bool sparse_ccs<value_type>::is_unique() const
{
    return m_refcount->get_count() == 1 || m_effective_unique == true;
};

template <class value_type>
sparse_ccs<value_type>::sparse_ccs(tinfo ti)
:m_ti(ti),m_refcount(refcount_str::create(1))
{
    alloc(0, 0, 0);
}

template <class value_type>
sparse_ccs<value_type>::sparse_ccs(tinfo ti,Integer r, Integer c, Integer nzmax)
:m_ti(ti),m_refcount(refcount_str::create(1))
{
    nzmax = std::max(nzmax,0);
    error::check_size_sp(r, c);
    alloc(r, c, nzmax);
}

template <class value_type>
inline sparse_ccs<value_type>::sparse_ccs(const sparse_ccs &other)
:m_ti(other.m_ti)
{	
    m_refcount  = other.m_refcount->increase();
    m_nzmax	    = other.m_nzmax;
    m_rows      = other.m_rows;
    m_cols      = other.m_cols;
    m_max_cols  = other.m_max_cols;
    m_offset    = other.m_offset;
    m_c 	    = other.m_c;
    m_c_root    = other.m_c_root;
    m_r 	    = other.m_r;
    m_x 	    = other.m_x;
    m_flag      = other.m_flag;

    m_effective_unique  = false;
}

template <class value_type>
inline sparse_ccs<value_type>::sparse_ccs(sparse_ccs &&other)
:m_ti(other.m_ti)
{	
    m_refcount  = other.m_refcount;
    m_nzmax	    = other.m_nzmax;
    m_rows      = other.m_rows;
    m_cols      = other.m_cols;
    m_max_cols  = other.m_max_cols;
    m_offset    = other.m_offset;
    m_c 	    = other.m_c;
    m_c_root    = other.m_c_root;
    m_r 	    = other.m_r;
    m_x 	    = other.m_x;
    m_flag      = other.m_flag;

    m_effective_unique  = false;
    other.m_refcount    = nullptr;
}

template <class value_type>
void sparse_ccs<value_type>::assign_to_fresh(const sparse_ccs &other)
{	
    other.m_refcount->increase();

    if (m_refcount->decrease())
    {
        destroy_data();
        m_refcount->destroy();
    };

    m_refcount  = other.m_refcount;
    m_nzmax	    = other.m_nzmax;
    m_rows      = other.m_rows;
    m_cols      = other.m_cols;
    m_max_cols  = other.m_max_cols;
    m_offset    = other.m_offset;
    m_c 	    = other.m_c;
    m_c_root    = other.m_c_root;
    m_r 	    = other.m_r;
    m_x 	    = other.m_x;
    m_flag      = other.m_flag;

    m_effective_unique  = false;
}

template <class value_type>
void sparse_ccs<value_type>::assign_to_fresh(sparse_ccs &&other)
{	
    std::swap(m_refcount, other.m_refcount);
    std::swap(m_r, other.m_r);
    std::swap(m_x, other.m_x);
    std::swap(m_nzmax, other.m_nzmax);
    std::swap(m_c_root, other.m_c_root);
    std::swap(m_max_cols, other.m_max_cols);

    m_rows      = other.m_rows;
    m_cols      = other.m_cols;
    m_offset    = other.m_offset;
    m_c 	    = other.m_c;
    m_flag      = other.m_flag;

    m_effective_unique  = false;
}

template <class value_type>
bool sparse_ccs<value_type>::has_element(Integer ri, Integer ci, Integer &k) const
{    
    if (nnz() == 0)
    {
        k = offset();
        return false;
    };

    const Integer* d_c = ptr_c();
    const Integer* d_r = ptr_r();

    Integer k_median;
    Integer k_high;

    k       = d_c[ci];
    k_high  = d_c[ci+1] - 1;

    if (k == k_high+1)
        return false;
    if (d_r[k] == ri)
        return true;
    if (d_r[k] > ri)
        return false;
    
    if (d_r[k_high] < ri) 
    { 
        k = k_high + 1; 
        return false; 
    } 

    if (d_r[k_high] == ri) 
    {
        k = k_high;
        return true; 
    }

    for(;;) 
    {
        k_median = k + (k_high - k) / 2;

        if (k == k_median) 
        { 
            k = k_high; 
            return false; 
        }
        if (d_r[k_median] == ri) 
        { 
            k = k_median; 
            return true; 
        }
        if (d_r[k_median] > ri)
        {
            k_high = k_median; 
        }
        else
        {
            k = k_median;
        };
    };
}

template <class value_type>
void sparse_ccs<value_type>::add_element(Integer ri, Integer ci, Integer k)
{
    matcl_assert(offset() == 0,"");

    Integer i;
    const Integer nz = nnz();

    if (nz == nzmax())
        add_memory();

    Integer * ptr_c     = this->ptr_c() + ci + 1;
    Integer * ptr_r     = this->ptr_r() + nz;
    value_type * ptr_x  = this->ptr_x() + nz;
    Integer c           = cols();

    for (i = ci + 1; i <= c; ++i)
        *(ptr_c++)	+= 1;

    for (i = nz; i > k; --i)
    {
        ptr_r[0] = ptr_r[-1];
        mrd::reset_helper(ptr_x[0],ptr_x[-1]);

        --ptr_r;
        --ptr_x;
    }

    *ptr_r = ri;
    mrd::reset_helper(*ptr_x,matcl::details::default_value<value_type>(get_type())); 
}

template <class value_type>
value_type sparse_ccs<value_type>::get_elem(Integer i, Integer j) const
{
    Integer k;

    if (!has_element(i, j, k))
        return matcl::details::default_value<value_type>(get_type());

    return value_type(ptr_x()[k]);
};

template <class value_type>
void sparse_ccs<value_type>::sort()
{
    Integer N = m_cols;

    if (nnz() == 0)
        return;

    for (Integer i = 0; i < N; ++i)
    {
        Integer k = ptr_c()[i];
        utils::sort_q(ptr_r()+k, ptr_x()+k, ptr_c()[i + 1] - k);
    }
}

template <class value_type>
bool sparse_ccs<value_type>::is_sorted() const
{
    Integer N = m_cols;
    if (nnz() == 0)
        return true;

    const Integer* d_c  = ptr_c();
    const Integer* d_r  = ptr_r();

    for (Integer i = 0; i < N; ++i)
    {
        Integer start = d_r[d_c[i]];
        for(Integer j = d_c[i]+1; j < d_c[i+1]; ++j)
        {
            Integer pos = d_r[j];
            if (pos <= start)
                return false;

            start = pos;
        };
    }
    return true;
};

template <class value_type>
void sparse_ccs<value_type>::add_duplications()
{    
    matcl_assert(offset() == 0,"");

    Integer N = cols();

    if (nnz() == 0)
        return;

    matcl::pod_workspace<Integer> ci(N + 1, 0);

    Integer* d_c    = ptr_c();
    Integer* d_r    = ptr_r();
    value_type* d_x = ptr_x();

    Integer sav_k   = -1;

    for (Integer i = 0; i < N; ++i)
    {
        Integer sav_row = -1;
        ci[i + 1]       = ci[i];

        for (Integer k = d_c[i]; k < d_c[i + 1]; ++k)
        {
            Integer cur_row = d_r[k];

            if (cur_row == sav_row)
            {
                mrd::reset_helper(d_x[sav_k],plus_helper<value_type,value_type>
                                    ::eval(d_x[sav_k],d_x[k]));
            }
            else
            {
                ++ci[i + 1];
                d_r[++sav_k]    = cur_row;
                sav_row         = cur_row;
                mrd::reset_helper(d_x[sav_k],d_x[k]);                
            }
        }
    }

    for (Integer i = 0; i <= N; ++i)
        d_c[i] = ci[i];
}

template <class value_type>
void sparse_ccs<value_type>::remove_element(Integer , Integer ci, Integer k)
{
    Integer nz = nnz();

    if (nz == 0)
        return;

    Integer * ptr_c     = this->ptr_c() + ci + 1;
    Integer * ptr_r     = this->ptr_r() + k;
    value_type * ptr_x  = this->ptr_x() + k;
    Integer c           = cols();

    for (Integer i = ci + 1; i <= c; ++i)
        *(ptr_c++)	-= 1;

    for (Integer i = k+1; i < nz; ++i)
    {
        ptr_r[0] = ptr_r[1];
        mrd::reset_helper(ptr_x[0],ptr_x[1]);

        ++ptr_r;
        ++ptr_x;
    }

    get_struct().reset();
}

template <class value_type>
void sparse_ccs<value_type>::alloc(Integer r, Integer c, Integer nzmax)
{
    nzmax       = std::max(nzmax,Integer());

    m_nzmax     = 0;
    m_c         = 0;
    m_r         = 0;
    m_x         = 0;
    m_c_root    = 0;
    m_offset    = 0;

    m_effective_unique  = false;

    try
    {
        m_nzmax = matcl::details::allocator<Integer>::malloc(ti::ti_empty(), 1);
        error::check_alloc(m_nzmax, 1, sizeof(Integer));	

        m_rows      = r; 
        m_cols      = c; 
        m_max_cols  = c;
        *m_nzmax    = nzmax;

        m_r     = matcl::details::allocator<Integer*>::malloc(1);
        m_x     = matcl::details::allocator<value_type*>::malloc(1);
        m_c_root= matcl::details::allocator<Integer*>::malloc(1);

        if (m_r)
            m_r[0] = 0;			

        if (m_x)
            m_x[0] = 0;

        if (m_c_root)
            m_c_root[0] = 0;

        error::check_alloc((void*) m_r, 1, sizeof(Integer*));
        error::check_alloc((void*) m_x, 1, sizeof(value_type*));
        error::check_alloc((void*) m_c_root, 1, sizeof(Integer*));

        if (c && r)
        {
            Integer* c_ptr  = matcl::details::allocator<Integer>::malloc(ti::ti_int(),c + 1);
            error::check_alloc(c_ptr, c + 1,sizeof(Integer));

            m_c         = c_ptr;
            m_c_root[0] = m_c;            

            ::memset(c_ptr, 0, (c+1) * sizeof(Integer));
        };

        if (nzmax && c && r)
        {
            Integer* r_ptr      = matcl::details::allocator<Integer>::malloc(ti::ti_int(),nzmax);
            error::check_alloc(r_ptr, nzmax,sizeof(Integer));

            value_type * x_ptr  = matcl::details::allocator<value_type>::malloc(m_ti,nzmax);
            error::check_alloc(x_ptr, nzmax,sizeof(value_type));

            *m_r    = r_ptr; 									
            *m_x    = x_ptr;
        };
    }
    catch(...)
    {
        destroy_data();
        throw;
    };
};

template <class value_type>
void sparse_ccs<value_type>::resize(Integer r, Integer c, Integer nzmax)
{
    matcl_assert(offset() == 0,"");

    nzmax       = std::max(nzmax,Integer());

    try
    {
        m_rows      = r; 
        m_cols      = c; 
        m_max_cols  = c;
        m_nzmax[0]  = nzmax;

        if (c && r)
        {
            m_c         = matcl::details::allocator<Integer>::malloc(ti::ti_empty(), c + 1);            
            error::check_alloc(m_c, c + 1,sizeof(Integer));

            m_c_root[0] = m_c;

            ::memset(m_c, 0, (c+1) * sizeof(Integer));
        };

        if (m_r)
            m_r[0] = 0;			

        if (m_x)
            m_x[0] = 0;

        if (nzmax && c && r)
        {
            Integer* r_ptr      = matcl::details::allocator<Integer>::malloc(ti::ti_int(),nzmax);
            error::check_alloc(r_ptr, nzmax,sizeof(Integer));

            value_type * x_ptr  = matcl::details::allocator<value_type>::malloc(m_ti,nzmax);
            error::check_alloc(x_ptr, nzmax,sizeof(value_type));

            *m_r = r_ptr; 									
            *m_x = x_ptr;
            m_flag.reset();
        }
        else
        {
            m_flag = struct_flag(predefined_struct_type::diag);
        };
    }
    catch(...)
    {
        destroy_data();
        throw;
    };
};

template <class value_type>
void sparse_ccs<value_type>::add_memory(Integer s)
{
    Integer off     = this->offset();    

    if (m_rows == 0 || m_cols == 0)
        return;

    Integer nzmax_old = m_nzmax[0];

    if (s < 0)
    {
        if (m_nzmax[0] == off) 
            return;
        else 
            m_nzmax[0] = off + nnz();
    }    

    if (s == 0)
    {
        if (m_nzmax[0] == off)
            m_nzmax[0] = 1 + off;
        else 
            m_nzmax[0] = imult_c(m_nzmax[0],2);
    }

    if (s > 0) 
        m_nzmax[0] += s;

    if (m_nzmax[0])
    {
        *m_r = matcl::details::allocator<Integer>::realloc(ti::ti_int(),*m_r, nzmax_old, m_nzmax[0]);
        error::check_alloc((void*) *m_r, m_nzmax[0],sizeof(Integer));

        value_type* ptr = m_x[0];
        ptr = matcl::details::allocator<value_type>::realloc(m_ti,ptr, nzmax_old, m_nzmax[0]);
        error::check_alloc((void*) ptr, m_nzmax[0],sizeof(value_type));
        *m_x = ptr;
    }
    else
    {
        matcl::details::allocator<Integer>::free(*m_r,nzmax_old);
        matcl::details::allocator<value_type>::free(*m_x,nzmax_old);
        *m_r = nullptr;
        *m_x = nullptr;
    }
}

template <class value_type>
void sparse_ccs<value_type>::add_columns(Integer nc)
{
    if (!m_c)
    {
        Integer new_c   = m_cols + nc;

        m_c_root= matcl::details::allocator<Integer*>::malloc(1);

        if (m_c_root)
            m_c_root[0] = 0;			

        error::check_alloc((void*) m_c_root, 1, sizeof(Integer*));

        Integer* c_ptr  = matcl::details::allocator<Integer>::malloc(ti::ti_empty(), new_c + 1);
        error::check_alloc(c_ptr, new_c + 1,sizeof(Integer));

        m_c_root[0]     = c_ptr;
        
        this->m_c       = c_ptr;
        this->m_max_cols= new_c;
    }
    else
    {
        Integer dc      = Integer(m_c - *m_c_root);
        Integer new_c   = m_max_cols + nc;

        *m_c_root       = matcl::details::allocator<Integer>
                            ::realloc(ti::ti_int(),*m_c_root, m_max_cols + 1, new_c + 1);

        error::check_alloc((void*) *m_c_root, new_c + 1, sizeof(Integer));

        m_c             = m_c_root[0] + dc;
        this->m_max_cols= new_c;
    }
}

template<class value_type> 
void sparse_ccs<value_type>::destroy_data()
{
    if (m_r)		matcl::details::allocator<Integer>::free(*m_r, *m_nzmax);
    if (m_x)		matcl::details::allocator<value_type>::free(*m_x, *m_nzmax);
    if (m_c_root)   matcl::details::allocator<Integer>::free(*m_c_root, m_max_cols + 1);

    matcl::details::allocator<Integer*>::free(m_c_root, 1);
    matcl::details::allocator<Integer*>::free(m_r, 1);
    matcl::details::allocator<value_type*>::free(m_x, 1);
    matcl::details::allocator<Integer>::free(m_nzmax, 1);
};

template <class value_type>
inline sparse_ccs<value_type> sparse_ccs<value_type>::make_unique(bool keep_maxcol) const
{
    if (this->is_unique())
        return *this;
    else
        return this->copy(keep_maxcol);
};

template <class value_type>
inline sparse_ccs<value_type> sparse_ccs<value_type>::copy(bool keep_maxcol) const
{
    Integer c = (keep_maxcol == true)? max_cols() : cols();
    Integer nz = nnz();

    sparse_ccs<value_type> tmp(get_type(),rows(), c, nz);
    tmp.m_flag = m_flag;
    tmp.m_cols = cols();

    if (nnz() == 0)
        return tmp;

    Integer* ptr_c      = tmp.ptr_c();
    Integer* ptr_r      = tmp.ptr_r();
    value_type* ptr_x   = tmp.ptr_x();

    const Integer* ptr_this_c   = sparse_ccs::ptr_c();
    const Integer* ptr_this_r   = sparse_ccs::ptr_r()+offset();
    const value_type* ptr_this_x = sparse_ccs::ptr_x()+offset();

    for(Integer i = 0; i <= c; ++i)
        ptr_c[i] = ptr_this_c[i] - offset();

    memcpy(ptr_r,ptr_this_r,imult_c(nz,sizeof(Integer)));

    //for Object assignment must be called, ptr_x for objects is 
    //initialized with default values.
    for(Integer i = 0; i < nz; ++i)
        mrd::reset_helper(ptr_x[i],ptr_this_x[i]);

    return tmp;
}

template <class value_type>
inline sparse_ccs<value_type> sparse_ccs<value_type>::clone(bool keep_maxcol) const
{
    Integer c = (keep_maxcol == true)? max_cols() : cols();
    Integer nz = nnz();

    sparse_ccs<value_type> tmp(get_type(),rows(), c, nz);
    tmp.m_flag = m_flag;
    tmp.m_cols = cols();

    if (nnz() == 0)
        return tmp;

    Integer* ptr_c      = tmp.ptr_c();
    Integer* ptr_r      = tmp.ptr_r();
    value_type* ptr_x   = tmp.ptr_x();

    const Integer* ptr_this_c   = sparse_ccs::ptr_c();
    const Integer* ptr_this_r   = sparse_ccs::ptr_r()+offset();
    const value_type* ptr_this_x = sparse_ccs::ptr_x()+offset();

    for(Integer i = 0; i <= c; ++i)
        ptr_c[i] = ptr_this_c[i] - offset();

    memcpy(ptr_r,ptr_this_r,imult_c(nz,sizeof(Integer)));

    //for Object copy constructor must be called
    for(Integer i = 0; i < nz; ++i)
        mrd::reset_helper(ptr_x[i],matcl::raw::details::clone_helper<value_type>::eval(ptr_this_x[i]));

    return tmp;
}

template <class value_type>
void sparse_ccs<value_type>::serialize(oarchive_impl & ar, const unsigned int v) const
{
    tinfo ti = get_type();
    Integer nz = nnz();
    Integer m_offset_s = 0;

    ar << m_rows;
    ar << m_cols;
    ar << m_offset_s;
    ar << m_nzmax[0];
    ar << nz;
    matcl::details::serialize_save(ar,ti,v);
    matcl::details::serialize_save(ar,m_flag,v);

    if (m_c)
    {
        if (nz == 0)
            return;
        
        for (Integer i = 0; i <= m_cols; ++i)
        {
            Integer c = m_c[i] - m_offset;
            matcl::details::serialize_save(ar,c,v);
        };
        
        const Integer* ptr_r = *m_r + m_offset;
        const value_type* ptr_x = *m_x + m_offset;

        for (Integer i = 0; i < nz; ++i)
            matcl::details::serialize_save(ar,*(ptr_r++),v);

        for (Integer i = 0; i < nz; ++i)
            matcl::details::serialize_save(ar,*(ptr_x++),v);
    };
};

template <class value_type>
void sparse_ccs<value_type>::serialize(iarchive_impl & ar, const unsigned int v)
{
    Integer r, c, nzmax, nz;

    ar >> r;
    ar >> c;
    ar >> m_offset;
    ar >> nzmax;
    ar >> nz;
    matcl::details::serialize_load(ar,m_ti,v);
    matcl::details::serialize_load(ar,m_flag,v);

    m_max_cols = c;

    resize(r,c,nzmax);

    if (nz == 0)
        return;

    try
    {
        Integer* ptr_c = m_c;
        for (Integer i = 0; i <= c; ++i)
            matcl::details::serialize_load(ar,*(ptr_c++),v);

        nz = m_c[c];

        Integer* ptr_r      = *m_r;
        value_type* ptr_x   = *m_x;

        for (Integer i = 0; i < nz; ++i)
            matcl::details::serialize_load(ar,*(ptr_r++),v);

        for (Integer i = 0; i < nz; ++i)
            matcl::details::serialize_load(ar,*(ptr_x++),v);
    }
    catch(...)
    {
        destroy_data();
        throw;
    };
};

};};};

template class matcl::raw::details::sparse_ccs<matcl::Integer>;
template class matcl::raw::details::sparse_ccs<matcl::Real>;
template class matcl::raw::details::sparse_ccs<matcl::Complex>;
template class matcl::raw::details::sparse_ccs<matcl::Object>;
template class matcl::raw::details::sparse_ccs<matcl::Float>;
template class matcl::raw::details::sparse_ccs<matcl::Float_complex>;
