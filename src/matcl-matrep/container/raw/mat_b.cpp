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

#include "matcl-internals/container/mat_b.h"
#include "matcl-internals/container/mat_d.h"
#include "matcl-internals/container/mat_s.h"
#include "matcl-matrep/container/matrix2.inl"
#include "matcl-matrep/func/raw/raw_manip.h"
#include "matcl-scalar/details/scalfunc_helpers.h"
#include "matcl-matrep/details/struct_flag_predefined.h"

namespace matcl { namespace raw 
{

namespace gr = matcl::raw;
namespace mrd = matcl::raw::details;
namespace md = matcl::details;

template <class value_type>
Matrix<value_type,struct_banded>::Matrix(tinfo ti) 
: base_type(ti)
{
    m_rows = m_cols = 0;
}

template<class value_type>
void Matrix<value_type,struct_banded>::destroy_data()
{
    Matrix<value_type,struct_banded>::base_type::m_data.destroy_data();
};

template <class value_type>
Matrix<value_type,struct_banded>::Matrix(tinfo ti, Integer r, Integer c, Integer fd, Integer ld)
: base_type(ti)
{
    Integer l   = std::max(-fd, 0);
    Integer u   = std::max(ld, 0);
    l           = std::min(l,r-1);
    u           = std::min(u,c-1);
    l           = std::max(l,0);
    u           = std::max(u,0);

    error::check_size_band(r, c, l, u);

    m_rows          = r;
    m_cols          = c;
    m_ldiags        = l;
    m_udiags        = u;

    if (r > 0)
        base_type::reset_unique(l + u + 1, c);
    else
        base_type::reset_unique(0, c);

    this->update_struct();
}

template <class value_type>
Matrix<value_type,struct_banded>::Matrix(tinfo ti,const value_type &val, 
                                         Integer r, Integer c, Integer fd, Integer ld)
    :Matrix(ti,r,c,fd,ld)
{
    value_type* ptr = this->rep_ptr();
    Integer s       = this->impl_size();

    for (Integer i = 0; i < s; ++i)
        ptr[i]      = val;
}

template <class value_type>
Matrix<value_type,struct_banded>::Matrix(const Matrix &m) 
    : base_type(m)
{
    m_rows          = m.m_rows;
    m_cols          = m.m_cols;
    m_ldiags        = m.m_ldiags;
    m_udiags        = m.m_udiags;
}

template <class value_type>
Matrix<value_type,struct_banded>::Matrix(Matrix &&m) 
    : base_type(std::move(m))
{
    m_rows          = m.m_rows;
    m_cols          = m.m_cols;
    m_ldiags        = m.m_ldiags;
    m_udiags        = m.m_udiags;
}

template <class value_type>
bool Matrix<value_type,struct_banded>::is_same_matrix(const Matrix &m) const
{
    if (this->rep_ptr() == m.rep_ptr() && this->rows() == m.rows() && this->cols() == m.cols())
        return true;
    else
        return false;
}

template <class value_type>
Matrix<value_type,struct_banded>&
Matrix<value_type,struct_banded>::reset_unique()
{
    m_rows = m_cols = 0;
    base_type::reset_unique();

    return *this;
}

template <class value_type>
Matrix<value_type,struct_banded>&
Matrix<value_type,struct_banded>::reset_unique(Integer r, Integer c, Integer fd, Integer ld)
{
    Integer l           = -fd;
    Integer u           = ld;

    l           = std::min(l,r-1);
    u           = std::min(u,c-1);
    l           = std::max(l,0);
    u           = std::max(u,0);
    error::check_size_band(r, c, l, u);

    m_rows      = r;
    m_cols      = c;
    m_ldiags    = l;
    m_udiags    = u;

    if (r > 0)
        base_type::reset_unique(l + u + 1, c);
    else
        base_type::reset_unique(0, c);

    return *this;
}

template <class value_type>
void Matrix<value_type,struct_banded>::assign_to_fresh(const Matrix &m)
{
    base_type::assign_to_fresh(m);

    m_rows      = m.m_rows;
    m_cols      = m.m_cols;
    m_ldiags    = m.m_ldiags;
    m_udiags    = m.m_udiags;
}

template <class value_type>
void Matrix<value_type,struct_banded>::assign_to_fresh(Matrix &&m)
{
    base_type::assign_to_fresh(std::move(m));

    m_rows      = m.m_rows;
    m_cols      = m.m_cols;
    m_ldiags    = m.m_ldiags;
    m_udiags    = m.m_udiags;
}

template <class value_type>
typename Matrix<value_type,struct_banded>::DenseMatrix
Matrix<value_type,struct_banded>::get_diag(Integer d) const
{
    using cur_matr_type = Matrix<value_type,struct_banded>;
    error::check_diag(d, m_rows, m_cols);

    Integer s;
    if (d > 0)
        s = std::min(m_cols - d, m_rows);
    else
        s = std::min(m_rows + d, m_cols);

    if (d > m_udiags || -d > m_ldiags)
    {
        value_type Z = matcl::details::default_value<value_type>(cur_matr_type::get_type());
        return DenseMatrix(cur_matr_type::get_type(),Z,s,1);
    };

    Integer st;
    if (d > 0)
        st = imult(d,ld()) + m_udiags - d;
    else
        st = m_udiags - d;

    DenseMatrix res(cur_matr_type::get_type(), s, 1);

    if (s == 0)
        return res;

    value_type * ptr           	= res.ptr();
    const value_type * ptr_this = Matrix::ptr() + st;
    Integer ld                  = this->ld();

    for (Integer i = 0; i < s; ++i)
    {
        mrd::reset_helper(*(ptr++),*ptr_this);
        ptr_this += ld;
    };

    return res;
};


template <class value_type>
Integer Matrix<value_type,struct_banded>::nnz() const
{
    Integer out = 0;
    for (Integer i = 1; i <= m_ldiags; ++i)
        out += std::min(m_rows - i, m_cols);

    for (Integer i = 1; i <= m_udiags; ++i)
        out += std::min(m_rows, m_cols - i);

    out += std::min(m_rows, m_cols);
    return out;
};

template <class value_type>
Matrix<value_type,struct_banded> Matrix<value_type,struct_banded>::make_unique(bool keep_bufor) const
{
    if (this->is_unique())
        return *this;
    else
        return this->copy(keep_bufor);
};

template <class value_type>
Matrix<value_type,struct_banded> Matrix<value_type,struct_banded>::copy(bool keep_bufor) const
{
    using cur_matr_type = Matrix<value_type,struct_banded>;

    Integer c   = (keep_bufor == true)? cur_matr_type::max_cols() : cols();
    Integer r   = (keep_bufor == true)? cur_matr_type::max_rows() : cur_matr_type::m_data.m_rows;

    Matrix out(cur_matr_type::get_type(), r, c, -std::max(r-1,0), 0);

    out.m_data.m_flag   = cur_matr_type::m_data.m_flag;
    out.m_rows          = rows();
    out.m_cols          = cols();
    out.m_ldiags        = m_ldiags;
    out.m_udiags        = m_udiags;

    if (keep_bufor == true)
    {
        out.m_data.m_rows   = cur_matr_type::m_data.m_rows;
        out.m_data.m_cols   = cur_matr_type::m_data.m_cols;
        out.m_data.m_size   = cur_matr_type::m_data.m_size;
    };

    Integer ld              = this->ld();

    if (keep_bufor == true)
    {
        Integer off         = cur_matr_type::m_data.m_root_ptr.offset(cur_matr_type::m_data.m_ptr) % ld;
        out.m_data.m_ptr    += off;
    };

    const value_type* ptr_this = ptr();
    value_type* ptr_out     = out.ptr();
    Integer out_ld          = out.ld();
    Integer this_ld         = this->ld();

    for (Integer j = 0; j < m_cols; ++j)
    {
        Integer fr      = first_row(j);
        Integer lr      = last_row(j);
        Integer pos     = first_elem_pos(j);

        for(Integer i = fr; i <= lr; ++i, ++pos)
            mrd::reset_helper(ptr_out[pos],ptr_this[pos]);

        ptr_out         += out_ld;
        ptr_this        += this_ld;
    };

    return out;
};

template <class value_type>
Matrix<value_type,struct_banded> Matrix<value_type,struct_banded>::clone(bool keep_bufor) const
{
    using cur_matr_type = Matrix<value_type,struct_banded>;

    Integer c   = (keep_bufor == true)? cur_matr_type::max_cols() : cols();
    Integer r   = (keep_bufor == true)? cur_matr_type::max_rows() : cur_matr_type::m_data.m_rows;

    Matrix out(cur_matr_type::get_type(), r, c, -std::max(r-1,0), 0);

    out.m_data.m_flag   = cur_matr_type::m_data.m_flag;
    out.m_rows          = rows();
    out.m_cols          = cols();
    out.m_ldiags        = m_ldiags;
    out.m_udiags        = m_udiags;
    Integer this_ld     = this->ld();

    if (keep_bufor == true)
    {
        out.m_data.m_rows   = cur_matr_type::m_data.m_rows;
        out.m_data.m_cols   = cur_matr_type::m_data.m_cols;
        out.m_data.m_size   = cur_matr_type::m_data.m_size;
    };

    if (keep_bufor == true)
    {
        Integer off         = cur_matr_type::m_data.m_root_ptr.offset(cur_matr_type::m_data.m_ptr)%this_ld;
        out.m_data.m_ptr    += off;
    };

    const value_type* ptr_this = ptr();
    value_type* ptr_out = out.ptr();
    Integer out_ld      = out.ld();

    for (Integer j = 0; j < m_cols; ++j)
    {
        Integer fr      = first_row(j);
        Integer lr      = last_row(j);
        Integer pos     = first_elem_pos(j);

        for(Integer i = fr; i <= lr; ++i, ++pos)
            mrd::reset_helper(ptr_out[pos],mrd::clone_helper<value_type>::eval(ptr_this[pos]));

        ptr_out         += out_ld;
        ptr_this        += this_ld;
    };

    return out;
};

template <class val_type>
matcl::Matrix Matrix<val_type,struct_banded>::fast_optim() const
{
    using cur_matr_type = Matrix<val_type,struct_banded>;
    if (rows() == 1 && cols() == 1)
        return operator()(1,1);

    if (cur_matr_type::get_struct().is_diag() == true)
        return matcl::Matrix(get_diag_band(),false);

    if (Real(nnz())/(rows()+1.0)/(cols()+1.0) > optim_params::max_sparse_density_max)
    {
        using full_matrix   = Matrix<val_type,struct_dense>;
        full_matrix tmp     = raw::converter<full_matrix,Matrix>::eval(*this);

        return matcl::Matrix(tmp,false);
    }

    return matcl::Matrix(*this,false);
};

template <class val_type>
Matrix<val_type,struct_banded>
Matrix<val_type,struct_banded>::get_diag_band() const
{
    using cur_matr_type = Matrix<val_type,struct_banded>;

    if (m_ldiags == 0 && m_udiags == 0)
        return *this;

    Matrix out(cur_matr_type::get_type(), m_rows, m_cols, 0, 0);
    Integer s = std::min(m_cols, m_rows);

    val_type * ptr              = out.ptr();
    const val_type * ptr_this	= Matrix::ptr() + m_udiags;
    Integer this_ld             = this->ld();

    for (Integer i = 0; i < s; ++i)
    {
        mrd::reset_helper(*(ptr++),*ptr_this);
        ptr_this += this_ld;
    };

    return out;
};

template <class val_type>
Matrix<val_type,struct_banded> 
Matrix<val_type,struct_banded>::reserve(Integer r, Integer c) const
{
    return reserve(r,c,0,0);
};

template <class val_type>
const Matrix<val_type,struct_banded> 
Matrix<val_type,struct_banded>::reserve(Integer r, Integer c, Integer fd, Integer ld) const
{
    using cur_matr_type = Matrix<val_type,struct_banded>;

    Integer l           = -fd;
    Integer u           = ld;

    if (c <= cur_matr_type::m_data.m_max_cols && l <= max_ldiags() && u <= max_udiags())
        return *this;

    Integer m_c = cols();
    Integer m_r = rows();
    Integer rm  = std::max(r,m_r);    
    Integer cm  = std::max(c,cur_matr_type::max_cols());
    Integer um  = std::max(m_udiags,u);
    Integer lm  = std::max(m_ldiags, l);    
    um          = std::max(std::min(cm-1,um), 0);
    lm          = std::max(std::min(rm-1,lm), 0);

    Matrix out(cur_matr_type::get_type(), rm, cm, -lm, um);

    out.m_rows          = m_r;
    out.m_cols          = m_c;
    out.m_ldiags        = m_ldiags;
    out.m_udiags        = m_udiags;
    out.m_data.m_rows   = m_ldiags + m_udiags + 1;
    out.m_data.m_size   = imult(out.m_data.m_rows,out.m_data.m_cols);
    out.m_data.m_ptr    = out.m_data.m_ptr + um - m_udiags;

    out.set_struct(cur_matr_type::get_struct());
    out.update_struct();

    const value_type* this_ptr = ptr();
    value_type* ptr     = out.ptr();
    Integer out_ld      = out.ld();
    Integer this_ld     = this->ld();

    Integer m_r_rep     = m_ldiags + m_udiags + 1;

    if (m_r > 0)
    {
        for (Integer j = 0; j < m_c; ++j)
        {
            for (Integer i = 0; i < m_r_rep; ++i)
                mrd::reset_helper(ptr[i],this_ptr[i]);

            ptr         += out_ld;
            this_ptr    += this_ld;
        }; 
    };

    return out;
};

template <class val_type>
const Matrix<val_type,struct_banded>
Matrix<val_type,struct_banded>::resize(Integer r, Integer c) const
{
    return resize(r, c, -m_ldiags, m_udiags);
};

template <class val_type>
Matrix<val_type,struct_banded> Matrix<val_type,struct_banded>::resize(Integer r, Integer c)
{
    return resize(r, c, -m_ldiags, m_udiags);
};

template <class val_type>
const Matrix<val_type,struct_banded> 
Matrix<val_type,struct_banded>::resize(Integer r, Integer c, Integer fd, Integer ld) const
{
    using cur_matr_type = Matrix<val_type,struct_banded>;

    error::check_resize(r,c);

    Integer l           = -fd;
    Integer u           = ld;

    if (c <= cols() && l <= m_ldiags && u <= m_udiags)
        return make_view(1, r, c, -l, u);

    if (c <= cur_matr_type::m_data.m_max_cols && l <= max_ldiags() && u <= max_udiags())
    {
        Matrix out = copy(true);
        return out.resize(r,c,-l,u);
    };

    Matrix out = reserve(r,c,-l,u);
    return out.resize(r,c,-l,u);
};

template <class val_type>
Matrix<val_type,struct_banded> 
Matrix<val_type,struct_banded>::resize(Integer r, Integer c, Integer fd, Integer ld)
{
    error::check_resize(r,c);

    Integer l           = -fd;
    Integer u           = ld;

    l   = std::max(std::min(r-1,l),0);
    u   = std::max(std::min(c-1,u),0);

    if (c <= cols() && l <= m_ldiags && u <= m_udiags)
        return make_view(1,r,c,-l,u);

    Matrix out      = reserve(r,c,-l,u);
    Integer m_c     = cols();
    value_type Z    = matcl::details::default_value<value_type>(Matrix<val_type,struct_banded>::get_type());
    Integer out_ld  = out.ld();

    if (u > m_udiags)
    {
        Integer du          = u - m_udiags;
        out.m_data.m_ptr    = out.m_data.m_ptr - du;
        value_type* ptr     = out.ptr();

        for (Integer j = 0; j < c; ++j)
        {
            for (Integer i = 0; i < du; ++i)
                mrd::reset_helper(ptr[i],Z);

            ptr += out_ld;
        };
    }
    else if (u < m_udiags)
    {
        Integer du          = m_udiags - u;
        out.m_data.m_ptr    = out.m_data.m_ptr + du;
    };

    out.m_udiags            = u;
    out.m_ldiags            = l;
    out.m_rows              = r;
    out.m_cols              = c;
    out.m_data.m_rows       = l + u + 1;
    out.m_data.m_cols       = c;
    out.m_data.m_size       = imult(out.m_data.m_rows,out.m_data.m_cols);

    value_type* ptr         = out.ptr(); 

    for (Integer j = 0; j < m_c; ++j)
    {
        Integer fr          = out.first_row(j);
        Integer pos         = out.first_elem_pos(j) - fr;
        Integer lr1         = pos + this->last_row(j);
        Integer lr2         = pos + out.last_row(j);

        for (Integer i = lr1+1; i <= lr2; ++i)
            mrd::reset_helper(ptr[i],Z);

        ptr += out_ld;
    };
    
    for (Integer j = m_c; j < c; ++j)
    {
        Integer fr          = out.first_row(j);
        Integer lr          = out.last_row(j);
        Integer pos         = out.first_elem_pos(j);

        for (Integer i = fr; i <= lr; ++i, ++pos)
            mrd::reset_helper(ptr[pos],Z);

        ptr += out_ld;
    }; 

    bool is_sym = (r == c) && (out.m_ldiags == out.m_udiags);
    out.m_data.get_struct() = md::predefined_struct::get_resize(out.m_data.get_struct(), is_sym,
                                                                is_real_matrix(out));
    out.update_struct();

    return out;
};

template <class val_type>
const Matrix<val_type,struct_banded>
Matrix<val_type,struct_banded>::make_view(Integer rc0, Integer r, Integer c) const
{
    return make_view(rc0, r,c, -m_ldiags, m_udiags);
};

template <class val_type>
const Matrix<val_type,struct_banded>
Matrix<val_type,struct_banded>::make_view(Integer rc0, Integer re, Integer ce, Integer fd, Integer ld) const
{
    Matrix out(*this);    

    Integer  r          = re - rc0 + 1;
    Integer  c          = ce - rc0 + 1;
    Integer out_ld      = out.ld();

    Integer l           = -fd;
    Integer u           = ld;
    l                   = std::max(std::min(r-1,l),0);
    u                   = std::max(std::min(c-1,u),0);

    out.m_rows          = r;
    out.m_cols          = c;
    out.m_ldiags        = l;
    out.m_udiags        = u;

    Integer du          = m_udiags - u;

    out.m_data.m_rows   = l+u+1;
    out.m_data.m_cols   = c;
    out.m_data.m_size   = imult(out.m_data.m_rows,out.m_data.m_cols);
    out.m_data.m_ptr    = const_cast<value_type*>(ptr()) + du + (rc0-1)*out_ld;

    bool is_sym             = (r == c) && (out.m_ldiags == out.m_udiags);
    out.m_data.get_struct() = md::predefined_struct::get_rectangle_view(out.m_data.get_struct(), is_sym);
    out.update_struct();

    return out;
};

template<class value_type>
bool Matrix<value_type,struct_banded>::all_finite() const
{
    return mr::all_finite_helper<Matrix>::eval(*this);
};

template<class value_type>
void Matrix<value_type,struct_banded>::serialize(oarchive_impl & ar, const unsigned int ) const
{
    ar << boost::serialization::base_object<base_type>(const_cast<Matrix&>(*this));

    ar << m_ldiags;
    ar << m_udiags;
    ar << m_rows;
    ar << m_cols;
};

template<class value_type>
void Matrix<value_type,struct_banded>::serialize(iarchive_impl & ar, const unsigned int )
{
    ar >> boost::serialization::base_object<base_type>(*this);

    ar >> m_ldiags;
    ar >> m_udiags;
    ar >> m_rows;
    ar >> m_cols;
};

template<class value_type>
void Matrix<value_type,struct_banded>::update_struct() const
{
    if (m_ldiags <= 1)
    {
        if (m_ldiags == 0)
            const_cast<struct_flag&>(get_struct()).add_ldiags(struct_flag::zero);
        else
            const_cast<struct_flag&>(get_struct()).add_ldiags(struct_flag::one);
    }
    if (m_udiags <= 1)
    {
        if (m_udiags == 0)
            const_cast<struct_flag&>(get_struct()).add_udiags(struct_flag::zero);
        else
            const_cast<struct_flag&>(get_struct()).add_udiags(struct_flag::one);
    }
};

};};

template class matcl::raw::Matrix<matcl::Integer,matcl::struct_banded>;
template class matcl::raw::Matrix<matcl::Real,matcl::struct_banded>;
template class matcl::raw::Matrix<matcl::Complex,matcl::struct_banded>;
template class matcl::raw::Matrix<matcl::Object,matcl::struct_banded>;
template class matcl::raw::Matrix<matcl::Float,matcl::struct_banded>;
template class matcl::raw::Matrix<matcl::Float_complex,matcl::struct_banded>;
