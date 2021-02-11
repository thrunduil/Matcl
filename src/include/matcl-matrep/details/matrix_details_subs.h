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

#pragma once

#include "matcl-matrep/matrix/matrix_sub.h"

namespace matcl
{

inline matcl::sub_matrix::sub_matrix(const sub_matrix& m)
:m_matrix(m.m_matrix),m_colon_1(m.m_colon_1),m_colon_2(m.m_colon_2),m_d(m.m_d)
{};

inline matcl::sub_matrix_1::sub_matrix_1(const sub_matrix_1& m)
:m_matrix(m.m_matrix),m_ind_1(m.m_ind_1)
{};

inline matcl::sub_matrix_2::sub_matrix_2(const sub_matrix_2& m)
:m_matrix(m.m_matrix),m_ind_1(m.m_ind_1),m_ind_2(m.m_ind_2)
{};

inline matcl::sub_matrix::sub_matrix(sub_matrix&& m)
:m_matrix(m.m_matrix),m_colon_1(m.m_colon_1),m_colon_2(m.m_colon_2),m_d(m.m_d)
{};

inline matcl::sub_matrix_1::sub_matrix_1(sub_matrix_1&& m)
:m_matrix(m.m_matrix),m_ind_1(m.m_ind_1)
{};

inline matcl::sub_matrix_2::sub_matrix_2(sub_matrix_2&& m)
:m_matrix(m.m_matrix),m_ind_1(m.m_ind_1),m_ind_2(m.m_ind_2)
{};

inline matcl::sub_matrix_1::sub_matrix_1(Matrix* m, Integer i)                        
    :m_matrix(m),m_ind_1(i)
{};

inline matcl::sub_matrix_2::sub_matrix_2(Matrix* m, Integer i, Integer j)            
    :m_matrix(m),m_ind_1(i),m_ind_2(j)
{};

inline matcl::sub_matrix::sub_matrix(Matrix* m, const colon& c1)                
    :m_matrix(m), m_colon_1(&c1),m_colon_2(nullptr),m_d(0)
{};

inline matcl::sub_matrix::sub_matrix(Matrix* m, const colon& c1, const colon& c2)
    :m_matrix(m), m_colon_1(&c1),m_colon_2(&c2),m_d(0)
{};

inline matcl::sub_matrix::sub_matrix(Integer d, Matrix* m)
    :m_matrix(m), m_colon_1(nullptr),m_colon_2(nullptr),m_d(d)
{};

namespace details
{

struct assign_scal_1_functor : public details::extract_type_switch<Matrix&,assign_scal_1_functor,false>
{
    template<class T, class S>
    static Matrix& eval(Matrix& handle, T& mat,Integer ind, const S& val)
    {
        return details::assign_functor<T,S>::eval(handle,mat,val,ind);
    };

    template<class T, class S>
    static Matrix& eval_scalar(Matrix& handle, T&,Integer ind, const S& val)
    {
        return details::assign_functor_scal<T,S>::eval(handle,val,ind);
    };
};

struct assign_scal_2_functor : public details::extract_type_switch<Matrix&,assign_scal_2_functor,false>
{
    template<class T, class S>
    static Matrix& eval(Matrix& handle, T& mat,Integer ind1, Integer ind2, const S& val)
    {
        return details::assign_functor<T,S>::eval(handle,mat,val,ind1,ind2);
    };

    template<class T, class S>
    static Matrix& eval_scalar(Matrix& handle, T&, Integer ind1,Integer ind2, const S& val)
    {
        return details::assign_functor_scal<T,S>::eval(handle,val,ind1,ind2);
    };
};

struct assign_scal_3_functor : public details::extract_type_switch<Matrix&,assign_scal_3_functor,false>
{
    template<class T, class S>
    static Matrix& eval(Matrix& handle, T& mat,const colon& c1, const S& val)
    {
        return details::assign_functor<T,S>::eval(handle,mat,val,c1);
    };

    template<class T, class S>
    static Matrix& eval_scalar(Matrix& handle, T&, const colon& c1, const S& val)
    {
        return details::assign_functor_scal<T,S>::eval(handle,val,c1);
    };
};

struct assign_scal_4_functor : public details::extract_type_switch<Matrix&,assign_scal_4_functor,false>
{
    template<class T, class S>
    static Matrix& eval(Matrix& handle, T& mat,const colon& c1,const colon& c2, const S& val)
    {
        return details::assign_functor<T,S>::eval(handle,mat,val,c1,c2);
    };

    template<class T, class S>
    static Matrix& eval_scalar(Matrix& handle, T&, const colon& c1, const colon& c2, const S& val)
    {
        return details::assign_functor_scal<T,S>::eval(handle,val,c1,c2);
    };
};

struct assign_diag_functor : public details::extract_type_switch<Matrix&,assign_diag_functor,false>
{
    template<class T, class S>
    static Matrix& eval(Matrix& handle, T& mat,Integer d, const S& val)
    {
        return details::assign_functor<T,S>::eval_diag(handle,mat,val,d);
    };

    template<class T, class S>
    static Matrix& eval_scalar(Matrix& handle, T&, Integer d, const S& val)
    {
        if (d != 0)
            throw error::invalid_diag(d, 1, 1);

        handle = val;
        return handle;
    };
};

};

template<class S>
inline Matrix& sub_matrix_1::assign_scalar(const S& val) const
{
    Matrix& tmp = *m_matrix;
    return details::assign_scal_1_functor::make<Matrix&>(tmp,m_ind_1,val);
};

template<class S>
inline Matrix& sub_matrix_2::assign_scalar(const S& val) const
{
    return details::assign_scal_2_functor::make<Matrix&>(*m_matrix,m_ind_1,m_ind_2,val);
};

template<class S>
inline Matrix& sub_matrix::assign_scalar(const S& val) const
{
    if (m_colon_2)
        return details::assign_scal_4_functor::make<Matrix&>(*m_matrix,*m_colon_1,*m_colon_2,val);
    else if(m_colon_1)
        return details::assign_scal_3_functor::make<Matrix&,const colon&,const S&>(*m_matrix,*m_colon_1,val);
    else
        return details::assign_diag_functor::make<Matrix&,const Integer&,const S&>(*m_matrix,m_d,val);
};

template<class S>
typename details::enable_if<details::is_scalar<S>::value,Matrix&>::type
inline sub_matrix::operator=(const S& val) const &&
{
    using val_type = typename matcl::details::promote_scalar<S>::type;
    return assign_scalar<val_type>(val);
};

template<class S>
typename details::enable_if<details::is_scalar<S>::value,Matrix&>::type
inline sub_matrix_1::operator=(const S& val) const &&
{
    using val_type = typename matcl::details::promote_scalar<S>::type;
    return assign_scalar<val_type>(val);
};

template<class S>
typename details::enable_if<details::is_scalar<S>::value,Matrix&>::type
inline sub_matrix_2::operator=(const S& val) const &&
{
    using val_type = typename matcl::details::promote_scalar<S>::type;
    return assign_scalar<val_type>(val);
};

template<class V> 
inline V sub_matrix::get_scalar() const
{
    return Matrix(*this).get_scalar<V>();
};
 
inline Integer sub_matrix::rows() const
{
    return Matrix(*this).rows();
};

inline Integer sub_matrix::cols() const
{
    return Matrix(*this).cols();
};

inline Integer sub_matrix::length() const
{
    return Matrix(*this).length();
};

inline Real sub_matrix::numel() const
{
    return Matrix(*this).numel();
};

inline bool sub_matrix::all_finite() const
{
    return Matrix(*this).all_finite();
};

inline sub_matrix::operator bool() const
{
    return Matrix(*this).operator bool();
};

template<class V> 
inline V sub_matrix_1::get_scalar() const
{
    return Matrix(*this).get_scalar<V>();
};

inline Integer sub_matrix_1::rows() const
{
    return Matrix(*this).rows();
};

inline Integer sub_matrix_1::cols() const
{
    return Matrix(*this).cols();
};

inline Integer sub_matrix_1::length() const
{
    return Matrix(*this).length();
};

inline Real sub_matrix_1::numel() const
{
    return Matrix(*this).numel();
};

inline bool sub_matrix_1::all_finite() const
{
    return Matrix(*this).all_finite();
};

inline sub_matrix_1::operator bool() const
{
    return Matrix(*this).operator bool();
};

template<class V> 
inline V sub_matrix_2::get_scalar() const
{
    return Matrix(*this).get_scalar<V>();
};

inline Integer sub_matrix_2::rows() const
{
    return Matrix(*this).rows();
};

inline Integer sub_matrix_2::cols() const
{
    return Matrix(*this).cols();
};

inline Integer sub_matrix_2::length() const
{
    return Matrix(*this).length();
};

inline Real sub_matrix_2::numel() const
{
    return Matrix(*this).numel();
};

inline bool sub_matrix_2::all_finite() const
{
    return Matrix(*this).all_finite();
};

inline sub_matrix_2::operator bool() const
{
    return Matrix(*this).operator bool();
};

};
