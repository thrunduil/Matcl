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

#pragma once

#include "matcl-matrep/matrix/matrix.h"
#include "matcl-core/error/exception_classes.h"
#include "matcl-matrep/details/matrix_details_subs.h"
#include "matcl-matrep/matrix/matrix_concat.h"
#include "matcl-core/matrix/matrix_traits.h"
#include "matcl-matrep/lib_functions/manip.h"

namespace matcl
{

namespace details
{
    template<class T>
    struct get_rep_functor
    {
        static MATCL_MATREP_EXPORT const T* eval_array_c(const Matrix& m);
        static MATCL_MATREP_EXPORT T*       eval_array_nc(Matrix& m);
    };
};

template<class T> 
inline Matrix::Matrix(T&& val, typename details::enable_if_is_conv_to_mat<T, true>::type)
{
    using T0 = typename std::decay<T>::type;
    using T2 = typename matcl::details::hide_object_type<T0>::type;

    details::constructor_helper<T2>::eval(std::forward<T>(val),true,*this);
};

template<class T> 
inline Matrix::Matrix(const object_type<T>& val)
{
    using T2 = Object;
    details::constructor_helper<T2>::eval(Object(val),true,*this);
};

template<class T> 
inline Matrix::Matrix(object_type<T>&& val)
{
    using T2 = Object;
    details::constructor_helper<T2>::eval(Object(std::move(val)),true,*this);
};

template<class T> 
inline Matrix::Matrix(T&& val,bool allow_conversions,
                typename details::enable_if_is_conv_to_mat<T, false>::type)
{
    using T0 = typename std::decay<T>::type;
    using T2 = typename matcl::details::hide_object_type<T0>::type;
    details::constructor_helper<T2>::eval(std::forward<T>(val),allow_conversions,*this);
};

template<class T, bool Is_safe>
inline Matrix::Matrix(const dense_matrix<T,Is_safe>& mat)
    :matrix_base(mat.to_matrix())
{};

template<class T, bool Is_safe>
inline Matrix::Matrix(dense_matrix<T,Is_safe>&& mat)
    :matrix_base(mat.to_matrix())
{};

template<class T, bool Is_const>
inline Matrix::Matrix(const sparse_matrix<T,Is_const>& mat)
    :matrix_base(mat.to_matrix())
{};

template<class T, bool Is_const>
inline Matrix::Matrix(sparse_matrix<T,Is_const>&& mat)
    :matrix_base(mat.to_matrix())
{};

template<class T, bool Is_const>
inline Matrix::Matrix(const band_matrix<T,Is_const>& mat)
    :matrix_base(mat.to_matrix())
{};

template<class T, bool Is_const>
inline Matrix::Matrix(band_matrix<T,Is_const>&& mat)
    :matrix_base(mat.to_matrix())
{};

template<class T>
inline Matrix::Matrix(const object_matrix<T>& mat)
    :matrix_base(mat.to_matrix())
{};

template<class T>
inline Matrix::Matrix(object_matrix<T>&& mat)
    :matrix_base(mat.to_matrix())
{};

template<class T>
inline Matrix::Matrix(const sub_dense_matrix<T>& mat)
    :Matrix(mat.to_matrix())
{};

template<class T>
inline Matrix::Matrix(const dense_row<T>& mat)
    :Matrix(mat.to_matrix())
{};

template<class T>
inline Matrix::Matrix(const dense_col<T>& mat)
    :Matrix(mat.to_matrix())
{};

template<class T>
inline Matrix::Matrix(const sub_sparse_matrix<T>& mat)
    :Matrix(mat.to_matrix())
{};

template<class T>
inline Matrix::Matrix(const sub_sparse_matrix_1<T>& mat)
    :Matrix(mat.to_matrix())
{};

template<class T>
inline Matrix::Matrix(const sub_sparse_matrix_2<T>& mat)
    :Matrix(mat.to_matrix())
{};

template<class T>
inline Matrix::Matrix(const sparse_row<T>& mat)
    :Matrix(mat.to_matrix())
{};

template<class T>
inline Matrix::Matrix(const sparse_col<T>& mat)
    :Matrix(mat.to_matrix())
{};

template<class T>
inline Matrix::Matrix(const sub_band_matrix<T>& mat)
    :Matrix(mat.to_matrix())
{};

template<class T>
inline Matrix::Matrix(const sub_band_matrix_1<T>& mat)
    :Matrix(mat.to_matrix())
{};

template<class T>
inline Matrix::Matrix(const sub_band_matrix_2<T>& mat)
    :Matrix(mat.to_matrix())
{};

template<class T>
inline Matrix::Matrix(const object_row<T>& mat)
    :Matrix(mat.to_matrix())
{};

template<class T>
inline Matrix::Matrix(const object_col<T>& mat)
    :Matrix(mat.to_matrix())
{};

template<class T>
inline Matrix::Matrix(const sub_object_matrix<T>& mat)
    :Matrix(mat.to_matrix())
{};

template<class T>
inline Matrix::Matrix(const sub_object_matrix_1<T>& mat)
    :Matrix(mat.to_matrix())
{};

template<class T>
inline Matrix::Matrix(const sub_object_matrix_2<T>& mat)
    :Matrix(mat.to_matrix())
{};

template<class M>  
inline const M&	Matrix::get_impl() const
{ 
    return details::get_functor<M>::eval(*this);
};

template<class M>  
inline details::rvalue_holder<M> Matrix::move_impl()
{ 
    return details::get_functor<M>::eval_move(*this);
};

template<class M>		 
inline M& Matrix::get_impl_unique()			
{ 
    make_unique();
    return details::get_functor<M>::eval(*this);
};

template<class M> 
inline const M Matrix::impl() const
{
    matcl::mat_code mat_type = matrix_traits::mat_type_info_type<M>::matrix_code;
    Matrix tmp = convert(*this,mat_type);
    return tmp.get_impl<M>();
};

template<class M> 
inline M& Matrix::impl_unique()
{
    matcl::mat_code mat_type = matrix_traits::mat_type_info_type<M>::matrix_code;
    if (mat_type == this->get_matrix_code())
    {
        return this->get_impl_unique<M>();
    }
    else
    {
        *this = convert(*this,mat_type);
        return this->get_impl_unique<M>();
    };
};

force_inline matcl::sub_matrix_1 Matrix::operator()(Integer i)
{
    matcl::sub_matrix_1 ret = matcl::sub_matrix_1(this,i);
    return ret;
};

force_inline matcl::sub_matrix_2 Matrix::operator()(Integer i,Integer j)
{
    matcl::sub_matrix_2 ret = matcl::sub_matrix_2(this,i,j);
    return ret;
};

force_inline matcl::sub_matrix Matrix::operator()(const colon& c1)
{
    matcl::sub_matrix ret = matcl::sub_matrix(this,c1);
    return ret;
};

force_inline matcl::sub_matrix Matrix::operator()(const colon& c1,const colon& c2)
{
    matcl::sub_matrix ret = matcl::sub_matrix(this,c1,c2);
    return ret;
};

template<class S>
inline typename details::enable_if<details::is_scalar<S>::value,mat_row&>::type
mat_row::add(S&& val)
{
    make_unique();
    using val_type = typename details::promote_scalar<S>::type;
    return add_scalar<val_type>(std::forward<S>(val));
};

template<class S>
inline typename details::enable_if<details::is_scalar<S>::value,mat_col&>::type
mat_col::add(S&& val)
{
    make_unique();
    using val_type = typename details::promote_scalar<S>::type;
    return add_scalar<val_type>(std::forward<S>(val));
};

template<class S>
inline mat_row& mat_row::add_scalar(const S& val)
{
    add_scalar_impl(val);
    return *this;
};

template<class S>
inline mat_row& mat_row::add_scalar(typename std::decay<S>::type&& val)
{
    S tmp(std::move(val));
    add_scalar_impl(tmp);
    return *this;
};

template<class S>
inline mat_col& mat_col::add_scalar(const S& val)
{
    add_scalar_impl(val);
    return *this;
};

template<class S>
inline mat_col& mat_col::add_scalar(typename std::decay<S>::type&& val)
{
    S tmp(std::move(val));
    add_scalar_impl(tmp);
    return *this;
};

template<class V, class Enable> 
inline V Matrix::get_scalar() const
{
    return details::get_scalar_functor<V>::eval_get(*this);
};

template<class V, class Enable> 
inline V& Matrix::get_scalar_unique()
{
    return details::get_scalar_functor<V>::eval_get_unique(*this);
};

template<class T, class Enable>
const T* Matrix::get_array() const
{
    return details::get_rep_functor<T>::eval_array_c(*this);
};

template<class T, class Enable>
T* Matrix::get_array_unique()
{
    return details::get_rep_functor<T>::eval_array_nc(*this);
};

namespace details
{
    struct matrix_data_accesser
    {
        using mat_ptr       = matrix_container_base;
        using refcount_str  = refcount_str_default;

        static mat_ptr*	get_ptr(Matrix& A)
        {
            return A.m_value.m_mat.m_mat_ptr;
        };

        static const mat_ptr* get_ptr(const Matrix& A)
        {
            return A.m_value.m_mat.m_mat_ptr;
        };

        static const refcount_str* get_refcount(const Matrix& A)
        {
            return A.m_value.m_mat.m_refcount;
        }

        static const Integer& get_val_int(const Matrix& A)
        {
            return A.m_value.val_int;
        };

        static Integer& get_val_int(Matrix& A)
        {
            return A.m_value.val_int;
        };

        static const Real& get_val_real(const Matrix& A)
        {
            return A.m_value.val_real;
        };

        static Real& get_val_real(Matrix& A)
        {
            return A.m_value.val_real;
        };

        static const Float& get_val_float(const Matrix& A)
        {
            return A.m_value.val_float;
        };

        static Float& get_val_float(Matrix& A)
        {
            return A.m_value.val_float;
        };

        static const Real& get_val_complex(const Matrix& A, size_t pos)
        {
            return A.m_value.val_complex[pos];
        };

        static Real& get_val_complex(Matrix& A, size_t pos)
        {
            return A.m_value.val_complex[pos];
        };

        static const Float& get_val_fcomplex(const Matrix& A, size_t pos)
        {
            return A.m_value.val_fcomplex[pos];
        };

        static Float& get_val_fcomplex(Matrix& A, size_t pos)
        {
            return A.m_value.val_fcomplex[pos];
        };

        static const Complex& get_val_complex(const Matrix& A)
        {
            return *reinterpret_cast<const Complex*>(A.m_value.val_complex);
        };

        static Complex&	get_val_complex(Matrix& A)
        {
            return *reinterpret_cast<Complex*>(A.m_value.val_complex);
        };

        static const Float_complex& get_val_fcomplex(const Matrix& A)
        {
            return *reinterpret_cast<const Float_complex*>(A.m_value.val_fcomplex);
        };

        static Float_complex& get_val_fcomplex(Matrix& A)
        {
            return *reinterpret_cast<Float_complex*>(A.m_value.val_fcomplex);
        };

        static const Object& get_val_object(const Matrix& A)
        {
            return A.get_object();
        };

        static Object& get_val_object(Matrix& A)
        {
            return A.get_object();
        };
        
        static matrix_base& get_base(Matrix& A)
        {
            return A;
        };
        
        static const matrix_base& get_base(const Matrix& A)
        {
            return A;
        };
        
        static void assign_to_const_mat(const Matrix& old_m,const Matrix& new_m);
    };
};

};
