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

#include "matcl-matrep/container/matrix2.inl"
#include "matcl-matrep/container/matrix_container.inl"
#include "matcl-matrep/matrix/matrix_concat.h"
#include "matcl-core/error/exception_classes.h"
#include "matcl-core/matrix/matrix_traits.h"
#include "matcl-matrep/lib_functions/matrix_gen.h"
#include "matcl-matrep/matrix/colon.h"
#include "matcl-core/lib_functions/constants.h"
#include "matcl-matrep/lib_functions/manip.h"
#include "matcl-matrep/details/extract_type_switch.h" 
#include "matcl-matrep/IO/matrix_io.h"
#include "matcl-matrep/matrix/unique_matrix.h"

namespace matcl { namespace details
{

namespace mr = matcl :: raw;

struct all_finite_vis : public extract_type_switch<bool,all_finite_vis,true>
{
    template<class T>
    static bool eval(const Matrix&, const T& mat)
    {
        return mr::all_finite_helper<T>::eval(mat);
    };

    template<class T>
    static bool eval_scalar(const Matrix&, const T& v)
    {
        return (bool)mrd::isfinite_helper<T>::eval(v);
    };
};

}};

namespace matcl
{

Matrix::Matrix()
{
    details::constructor_helper<double>::eval(0.,true,*this);
};

Matrix::Matrix(bool val)
{
    details::constructor_helper<Integer>::eval(Integer(val),false,*this);
};

Matrix::Matrix(const Matrix& mat)
    :matrix_base(mat)
{};

Matrix::Matrix(Matrix&& mat)
    :matrix_base(std::move(mat))
{};

Matrix::Matrix(const mat_row& mr)
    :matrix_base(mr.to_matrix())
{};

Matrix::Matrix(const mat_col& mc)
    :matrix_base(mc.to_matrix())
{};

Matrix::~Matrix()
{};

Integer Matrix::rows() const
{
    if (m_type < mat_code::integer_dense)
        return 1;

    return m_value.m_mat.m_mat_ptr->rows();
}; 

Integer	Matrix::cols() const
{
    if (m_type < mat_code::integer_dense)
        return 1;

    return m_value.m_mat.m_mat_ptr->cols();
}; 

Integer	Matrix::structural_nnz() const
{
    switch(m_type)
    {
        case mat_code::integer_scalar:
        case mat_code::real_scalar:
        case mat_code::float_scalar:
        case mat_code::complex_scalar:
        case mat_code::float_complex_scalar:
        case mat_code::object_scalar:
        {
            return 1;
        }
        default:
        {
            return m_value.m_mat.m_mat_ptr->structural_nnz();
        }
    };
}; 

Integer	Matrix::structural_ldiags(bool use_flags) const
{
    if (m_type < mat_code::integer_dense)
        return 0;

    return m_value.m_mat.m_mat_ptr->structural_ldiags(use_flags);
};

Integer Matrix::structural_udiags(bool use_flags) const
{
    if (m_type < mat_code::integer_dense)
        return 0;

    return m_value.m_mat.m_mat_ptr->structural_udiags(use_flags);
};

Integer	Matrix::length() const
{
    Integer m   = rows();
    Integer n   = cols();
    
    if (m == 0 || n == 0)
        return 0;

    return (m > n ? m : n);
};

Real Matrix::numel() const
{
    Integer m   = rows();
    Integer n   = cols();

    return Real(m)*Real(n);
};

bool Matrix::is_empty() const
{
    Integer m   = rows();
    Integer n   = cols();
    return (m == 0 || n == 0);
};

bool Matrix::is_scalar() const
{
    Integer m   = rows();
    Integer n   = cols();
    return (m == 1 && n == 1);
};

bool Matrix::is_matrix_type() const
{
    if (m_type < mat_code::integer_dense)
        return false;

    return true;
};

bool Matrix::is_scalar_type() const
{
    if (m_type < mat_code::integer_dense)
        return true;

    return false;
};

bool Matrix::is_square() const
{
    Integer m   = rows();
    Integer n   = cols();
    return (m == n);
};

bool Matrix::is_vector() const
{
    Integer m   = rows();
    Integer n   = cols();
    return (m == 1 || n == 1 || m == 0 || n == 0);
};

bool Matrix::is_unique() const
{
    if (m_type < mat_code::integer_dense)
        return true;

    return this->is_effective_unique_mat();
};

void Matrix::mark_unique(bool unique)
{
    if (m_type < mat_code::integer_dense)
        return;

    m_value.m_mat.m_mat_ptr->mark_unique(unique);
};

matcl::value_code Matrix::get_value_code() const
{
    return matrix_traits::get_value_type(m_type);
};

matcl::struct_code Matrix::get_struct_code() const
{
    return matrix_traits::get_struct_type(m_type);
};

Matrix& Matrix::operator=(const matcl::sub_matrix& mat) &
{
    matrix_base::operator=(mat.to_matrix());
    return *this;
}

Matrix& Matrix::operator=(const matcl::sub_matrix_1& mat) &
{
    matrix_base::operator=(mat.to_matrix());
    return *this;
}

Matrix& Matrix::operator=(const matcl::sub_matrix_2& mat) &
{
    matrix_base::operator=(mat.to_matrix());
    return *this;
}

Matrix& Matrix::operator=(const Matrix& mat) &
{
    matrix_base::operator =(mat);
    return *this;
};

Matrix& Matrix::operator=(Matrix&& mat) &
{
    matrix_base::operator =(std::move(mat));
    return *this;
};

struct_flag& Matrix::get_struct() const
{
    static struct_flag struct_scalar;

    if (m_type < mat_code::integer_dense)
        return struct_scalar;

    struct_flag& ret = m_value.m_mat.m_mat_ptr->get_struct();
    return ret;
};

void Matrix::set_struct(const struct_flag& fl) const
{
    if (m_type < mat_code::integer_dense)
        return;

    //this is probably thread safe
    //struct is always proper
    //if different thread does not see the change, then all operations are still valid
    //unless struct is checked many times during one operation
    //struct change must be atomic
    m_value.m_mat.m_mat_ptr->set_struct(fl);
};

void Matrix::add_struct(const struct_flag& fl) const
{
    if (m_type < mat_code::integer_dense)
        return;

    //this is probably thread safe
    //struct is always proper
    //if different thread does not see the change, then all operations are still valid
    //unless struct is checked many times during one operation
    //struct change must be atomic
    m_value.m_mat.m_mat_ptr->add_struct(fl);
};

const details::matrix_container_base* Matrix::get_mat_ptr_prv() const
{
    switch(m_type)
    {
        case mat_code::integer_scalar:
        case mat_code::real_scalar:
        case mat_code::float_scalar:
        case mat_code::complex_scalar:
        case mat_code::float_complex_scalar:
        case mat_code::object_scalar:
        {
            throw error::invalid_return_type();
        }
        default:
        {
            return m_value.m_mat.m_mat_ptr;
        }
    };
};

Matrix::operator bool() const
{
    switch(m_type)
    {
        case mat_code::integer_scalar:
        {
            return (m_value.val_int != 0);
        }
        case mat_code::real_scalar:
        {
            return (m_value.val_real != 0.);
        }
        case mat_code::float_scalar:
        {
            return (m_value.val_float != 0.f);
        }
        case mat_code::complex_scalar:
        {
            return (m_value.val_complex[0] != 0. || m_value.val_complex[1] != 0.);
        }
        case mat_code::float_complex_scalar:
        {
            return (m_value.val_fcomplex[0] != 0.f || m_value.val_fcomplex[1] != 0.f);
        }
        case mat_code::object_scalar:
        {
            return get_object().is_zero() == false;
        }
        default:
        {
            return m_value.m_mat.m_mat_ptr->is_scalar_true();
        }
    };
};

Matrix::Matrix(const matcl::sub_matrix& sm)
    :matrix_base(sm.to_matrix())
{};

Matrix::Matrix(const matcl::sub_matrix_1& sm)
    :matrix_base(sm.to_matrix())
{};

Matrix::Matrix(const matcl::sub_matrix_2& sm)
    :matrix_base(sm.to_matrix())
{};

Matrix::Matrix(unique_matrix&& umat)
    :matrix_base(std::move(umat).to_matrix())
{
    mark_unique(true);
};

ti::ti_object Matrix::get_type() const
{
    switch (m_type)
    {
        case mat_code::integer_scalar:      return ti::predefined::get_ti_int();
        case mat_code::real_scalar:         return ti::predefined::get_ti_real();
        case mat_code::float_scalar:        return ti::predefined::get_ti_float();
        case mat_code::complex_scalar:      return ti::predefined::get_ti_complex();
        case mat_code::float_complex_scalar:return ti::predefined::get_ti_float_complex();
        case mat_code::object_scalar:       return base_type::get_object().get_type();
    };

    return m_value.m_mat.m_mat_ptr->get_type_info();
};

matcl::sub_matrix Matrix::diag(Integer d)
{
    matcl::sub_matrix ret = matcl::sub_matrix(d,this);
    return ret;
};

const Matrix Matrix::diag(Integer d) const
{
    return get_diag(*this,d);
};

void Matrix::resize(Integer r, Integer c)
{
    if (m_type < mat_code::integer_dense)
    {
        *this = full(*this);
        *this = m_value.m_mat.m_mat_ptr->resize(r,c);
        return;
    };
    *this = m_value.m_mat.m_mat_ptr->resize(r,c);
};

void Matrix::reserve(Integer r, Integer c)
{
    if (m_type < mat_code::integer_dense)
    {
        *this = full(*this);
        *this = m_value.m_mat.m_mat_ptr->reserve(r,c);
        return;
    };
    *this = m_value.m_mat.m_mat_ptr->reserve(r,c);
};

void Matrix::resize_band(Integer r, Integer c, Integer fd, Integer ld)
{
    if (m_type < mat_code::integer_dense)
    {
        *this = full(*this);
        *this = m_value.m_mat.m_mat_ptr->resize_band(r,c,fd,ld);
        return;
    };
    *this = m_value.m_mat.m_mat_ptr->resize_band(r,c,fd,ld);
};

void Matrix::reserve_band(Integer r, Integer c, Integer fd, Integer ld)
{
    if (m_type < mat_code::integer_dense)
    {
        *this = full(*this);
        *this = m_value.m_mat.m_mat_ptr->reserve_band(r,c,fd,ld);
        return;
    };
    *this = m_value.m_mat.m_mat_ptr->reserve_band(r,c,fd,ld);
};

const Matrix& Matrix::make_unique() const
{
    switch(m_type)
    {
        case mat_code::integer_scalar:
        case mat_code::real_scalar:
        case mat_code::float_scalar:
        case mat_code::complex_scalar:
        case mat_code::float_complex_scalar:
            return *this;
        case mat_code::object_scalar:
        {
            return *this;
        }
        default:
        {
            if (m_value.m_mat.m_refcount->is_unique())
            {
                return *this;
            };
        }
    };

    if (this->is_effective_unique_mat() == true)
        return *this;

    Matrix m2                   = this->clone();

    //this is thread safe since instances cannot be shared between threads
    const_cast<Matrix&>(*this)  = m2;

    return *this;
};

void details::matrix_data_accesser::assign_to_const_mat(const Matrix& old_m,const Matrix& new_m)
{
    old_m.matrix_base::assign(new_m);
};

bool Matrix::all_finite() const
{
    return details::all_finite_vis::make<const Matrix&>(*this);
};

};

namespace matcl { namespace details 
{
    void matrix_base::decrease_refcount() const
    {
        if (m_type < mat_code::integer_dense)
            return;

        if(m_value.m_mat.m_refcount->decrease())
        {                
            destroy_container(m_value.m_mat.m_mat_ptr);
            m_value.m_mat.m_refcount->destroy();
        }
        else if (m_value.m_mat.m_mat_ptr->decrease())
        {   
            free_container(m_value.m_mat.m_mat_ptr);
        };
    };

    void matrix_base::increase_refcount() const
    {
        if (m_type < mat_code::integer_dense)
            return;

        m_value.m_mat.m_refcount->increase();
        m_value.m_mat.m_mat_ptr->increase();
    };
}}
