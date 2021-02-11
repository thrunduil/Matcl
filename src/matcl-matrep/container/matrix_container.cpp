/*
 *  This file is a part of Matrix Computation Library (MATCL)
 *
 *  Copyright (c) Pawe³ Kowal 2017 - 2021
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

#include "matcl-matrep/container/matrix_container.inl"
#include "matcl-core/details/integer.h"

namespace matcl { namespace details
{

struct eval_is_same_matrix : public eval_type<bool,eval_is_same_matrix>
{
    template<class T>
    static bool eval(T* mat, const matrix_container_base* other)
    { 
        return mat->is_same_matrix_impl(other); 
    };
};

bool matrix_container_base::is_same_matrix(const matrix_container_base& other) const
{
    if (this->m_type != other.m_type)
        return false;

    return eval_is_same_matrix::make<const matrix_container_base*>(this, &other);
};

Integer matrix_container_base::rows() const
{
    return eval_type_rows::make(this);
}; 

Integer matrix_container_base::cols() const
{
    return eval_type_cols::make(this);
}; 

Integer matrix_container_base::structural_nnz() const
{
    return eval_type_struct_nnz::make(this);
}; 

Integer matrix_container_base::structural_ldiags(bool use_flags) const
{
    return eval_type_struct_ldiags::make<bool>(this, use_flags);
}; 

Integer matrix_container_base::structural_udiags(bool use_flags) const
{
    return eval_type_struct_udiags::make<bool>(this, use_flags);
}; 

int_tup_2 matrix_container_base::size() const
{
    return eval_type_size::make(this);
}; 

bool matrix_container_base::is_scalar_true() const
{
    return eval_type_is_scalar_true::make(this);
}; 

ti::ti_object matrix_container_base::get_type_info() const
{
    return eval_type_get_ti::make(this);
}; 

void matrix_container_base::disp(const disp_stream_ptr& os) const
{
    return eval_type_disp::make<const disp_stream_ptr&>(this,os);
}; 

struct_flag& matrix_container_base::get_struct()
{
    return eval_type_get_struct::make(this);
};

const struct_flag& matrix_container_base::get_struct() const
{
    return eval_type_get_struct_const::make(this);
};

void matrix_container_base::set_struct(const struct_flag& fl) const
{
    return eval_type_set_struct::make<const struct_flag&>(this,fl);
};

void matrix_container_base::add_struct(const struct_flag& fl) const
{
    return eval_type_add_struct::make<const struct_flag&>(this,fl);
};

matcl::Matrix matrix_container_base::get_scalar() const
{
    return eval_type_get_scalar::make(this);
};

matcl::Matrix matrix_container_base::reserve(Integer r, Integer c) const
{ 
    return eval_type_reserve::make<Integer,Integer>(this,r,c);
};

matcl::Matrix matrix_container_base::resize(Integer r, Integer c) const
{
    return eval_type_resize::make<Integer,Integer>(this,r,c);
};

struct int_4
{
    Integer r,c,fd,ld;
    int_4(Integer r_, Integer c_, Integer fd_, Integer ld_)
        :r(r_),c(c_),fd(fd_),ld(ld_)
    {};
};

matcl::Matrix matrix_container_base::reserve_band(Integer r, Integer c, Integer fd, Integer ld) const
{
    return eval_type_reserve_band::make<const int_4&>(this,int_4(r,c,fd,ld));
};

matcl::Matrix matrix_container_base::resize_band(Integer r, Integer c, Integer fd, Integer ld) const
{
    return eval_type_resize_band::make<const int_4&>(this,int_4(r,c,fd,ld));
};

};};