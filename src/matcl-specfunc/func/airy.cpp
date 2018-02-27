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

#include "matcl-specfunc/lib_functions/airy.h"
#include "matcl-matrep/lib_functions/func_unary.h"
#include "matcl-specfunc/objects/object_func.h"
#include "matcl-internals/func/converter.h"
#include "matcl-matrep/lib_functions/manip.h"

#pragma warning(push)
#pragma warning(disable:4127)   //conditional expression is constant
#include <boost/math/special_functions.hpp>
#pragma warning(pop)

namespace matcl { namespace details
{

namespace md = matcl::details;

namespace ignore_errors 
{
    using boost__ignore_policy
        = ::boost::math::policies::policy 
          <
            ::boost::math::policies::domain_error<::boost::math::policies::ignore_error>,
            ::boost::math::policies::pole_error<::boost::math::policies::ignore_error>,
            ::boost::math::policies::overflow_error<::boost::math::policies::ignore_error>,
            ::boost::math::policies::underflow_error<::boost::math::policies::ignore_error>,
            ::boost::math::policies::rounding_error<::boost::math::policies::ignore_error>,
            ::boost::math::policies::denorm_error<::boost::math::policies::ignore_error>,
            ::boost::math::policies::indeterminate_result_error<::boost::math::policies::ignore_error>,
            ::boost::math::policies::evaluation_error<::boost::math::policies::ignore_error>
         >;
    BOOST_MATH_DECLARE_SPECIAL_FUNCTIONS(boost__ignore_policy)
}

//-------------------------------------------------------------------------------
//                                  REAL
//-------------------------------------------------------------------------------
template<>
typename airy_helper<Real>::return_type
airy_helper<Real>::eval_airy_ai(const Real& x)
{
    return ignore_errors::airy_ai(x);
};

template<>
typename airy_helper<Real>::return_type
airy_helper<Real>::eval_airy_bi(const Real& x)
{
    return ignore_errors::airy_bi(x);
};

template<>
typename airy_helper<Real>::return_type
airy_helper<Real>::eval_airy_ai_dif(const Real& x)
{
    return ignore_errors::airy_ai_prime(x);
};

template<>
typename airy_helper<Real>::return_type
airy_helper<Real>::eval_airy_bi_dif(const Real& x)
{
    return ignore_errors::airy_bi_prime(x);
};

template<>
typename airy_helper<Real>::return_type
airy_helper<Real>::eval_airy_ai_zero(Integer x)
{
    return ignore_errors::airy_ai_zero<Real>(x);
};

template<>
typename airy_helper<Real>::return_type
airy_helper<Real>::eval_airy_bi_zero(Integer x)
{
    return ignore_errors::airy_bi_zero<Real>(x);
};

//-------------------------------------------------------------------------------
//                                  FLOAT
//-------------------------------------------------------------------------------
template<>
typename airy_helper<Float>::return_type
airy_helper<Float>::eval_airy_ai(const Float& x)
{
    return ignore_errors::airy_ai(x);
};

template<>
typename airy_helper<Float>::return_type
airy_helper<Float>::eval_airy_bi(const Float& x)
{
    return ignore_errors::airy_bi(x);
};

template<>
typename airy_helper<Float>::return_type
airy_helper<Float>::eval_airy_ai_dif(const Float& x)
{
    return ignore_errors::airy_ai_prime(x);
};

template<>
typename airy_helper<Float>::return_type
airy_helper<Float>::eval_airy_bi_dif(const Float& x)
{
    return ignore_errors::airy_bi_prime(x);
};

template<>
typename airy_helper<Float>::return_type
airy_helper<Float>::eval_airy_ai_zero(Integer x)
{
    return ignore_errors::airy_ai_zero<Float>(x);
};

template<>
typename airy_helper<Float>::return_type
airy_helper<Float>::eval_airy_bi_zero(Integer x)
{
    return ignore_errors::airy_bi_zero<Float>(x);
};

//-------------------------------------------------------------------------------
//                                  COMPLEX
//-------------------------------------------------------------------------------
template<>
typename airy_helper<Complex>::return_type
airy_helper<Complex>::eval_airy_ai(const Complex& x)
{
    if (imag(x) == 0.0)
        return airy_helper<Real>::eval_airy_ai(real(x));
    else
        throw matcl::error::function_not_defined_for_complex("airy_ai"); 
};

template<>
typename airy_helper<Complex>::return_type
airy_helper<Complex>::eval_airy_bi(const Complex& x)
{
    if (imag(x) == 0.0)
        return airy_helper<Real>::eval_airy_bi(real(x));
    else
        throw matcl::error::function_not_defined_for_complex("airy_bi"); 
};

template<>
typename airy_helper<Complex>::return_type
airy_helper<Complex>::eval_airy_ai_dif(const Complex& x)
{
    if (imag(x) == 0.0)
        return airy_helper<Real>::eval_airy_ai_dif(real(x));
    else
        throw matcl::error::function_not_defined_for_complex("airy_ai_dif"); 
};

template<>
typename airy_helper<Complex>::return_type
airy_helper<Complex>::eval_airy_bi_dif(const Complex& x)
{
    if (imag(x) == 0.0)
        return airy_helper<Real>::eval_airy_bi_dif(real(x));
    else
        throw matcl::error::function_not_defined_for_complex("airy_bi_dif"); 
};

template<>
typename airy_helper<Complex>::return_type
airy_helper<Complex>::eval_airy_ai_zero(Integer x)
{
    return ignore_errors::airy_ai_zero<Real>(x);
};

template<>
typename airy_helper<Complex>::return_type
airy_helper<Complex>::eval_airy_bi_zero(Integer x)
{
    return ignore_errors::airy_bi_zero<Real>(x);
};

//-------------------------------------------------------------------------------
//                              FLOAT_COMPLEX
//-------------------------------------------------------------------------------
template<>
typename airy_helper<Float_complex>::return_type
airy_helper<Float_complex>::eval_airy_ai(const Float_complex& x)
{
    if (imag(x) == 0.0)
        return airy_helper<Float>::eval_airy_ai(real(x));
    else
        throw matcl::error::function_not_defined_for_complex("airy_ai"); 
};

template<>
typename airy_helper<Float_complex>::return_type
airy_helper<Float_complex>::eval_airy_bi(const Float_complex& x)
{
    if (imag(x) == 0.0)
        return airy_helper<Float>::eval_airy_bi(real(x));
    else
        throw matcl::error::function_not_defined_for_complex("airy_bi"); 
};

template<>
typename airy_helper<Float_complex>::return_type
airy_helper<Float_complex>::eval_airy_ai_dif(const Float_complex& x)
{
    if (imag(x) == 0.0)
        return airy_helper<Float>::eval_airy_ai_dif(real(x));
    else
        throw matcl::error::function_not_defined_for_complex("airy_ai_dif"); 
};
template<>
typename airy_helper<Float_complex>::return_type
airy_helper<Float_complex>::eval_airy_bi_dif(const Float_complex& x)
{
    if (imag(x) == 0.0)
        return airy_helper<Float>::eval_airy_bi_dif(real(x));
    else
        throw matcl::error::function_not_defined_for_complex("airy_bi_dif"); 
};
template<>
typename airy_helper<Float_complex>::return_type
airy_helper<Float_complex>::eval_airy_ai_zero(Integer x)
{
    return ignore_errors::airy_ai_zero<Float>(x);
};
template<>
typename airy_helper<Float_complex>::return_type
airy_helper<Float_complex>::eval_airy_bi_zero(Integer x)
{
    return ignore_errors::airy_bi_zero<Float>(x);
};

//-------------------------------------------------------------------------------
//                              INTEGER
//-------------------------------------------------------------------------------
template<>
typename airy_helper<Integer>::return_type
airy_helper<Integer>::eval_airy_ai(const Integer& x)
{
    return airy_helper<Real>::eval_airy_ai((Real)x);
};
template<>
typename airy_helper<Integer>::return_type
airy_helper<Integer>::eval_airy_bi(const Integer& x)
{
    return airy_helper<Real>::eval_airy_bi((Real)x);
};
template<>
typename airy_helper<Integer>::return_type
airy_helper<Integer>::eval_airy_ai_dif(const Integer& x)
{
    return airy_helper<Real>::eval_airy_ai_dif((Real)x);
};
template<>
typename airy_helper<Integer>::return_type
airy_helper<Integer>::eval_airy_bi_dif(const Integer& x)
{
    return airy_helper<Real>::eval_airy_bi_dif((Real)x);
};
template<>
typename airy_helper<Integer>::return_type
airy_helper<Integer>::eval_airy_ai_zero(Integer x)
{
    return airy_helper<Real>::eval_airy_ai_zero(x);
};
template<>
typename airy_helper<Integer>::return_type
airy_helper<Integer>::eval_airy_bi_zero(Integer x)
{
    return airy_helper<Real>::eval_airy_bi_zero(x);
};

//-------------------------------------------------------------------------------
//                              OBJECT
//-------------------------------------------------------------------------------
template<>
typename airy_helper<Object>::return_type
airy_helper<Object>::eval_airy_ai(const Object& x)
{
    return object_func::airy_ai(x);
};
template<>
typename airy_helper<Object>::return_type
airy_helper<Object>::eval_airy_bi(const Object& x)
{
    return object_func::airy_bi(x);
};
template<>
typename airy_helper<Object>::return_type
airy_helper<Object>::eval_airy_ai_dif(const Object& x)
{
    return object_func::airy_ai_dif(x);
};
template<>
typename airy_helper<Object>::return_type
airy_helper<Object>::eval_airy_bi_dif(const Object& x)
{
    return object_func::airy_bi_dif(x);
};
template<>
typename airy_helper<Object>::return_type
airy_helper<Object>::eval_airy_ai_zero(Integer x)
{
    (void)x;
    //TODO: impl for object
    throw matcl::error::object_value_type_not_allowed("airy_ai_zero");
};
template<>
typename airy_helper<Object>::return_type
airy_helper<Object>::eval_airy_bi_zero(Integer x)
{
    (void)x;
    //TODO: impl for object
    throw matcl::error::object_value_type_not_allowed("airy_bi_zero");
};


template struct airy_helper<Integer>;
template struct airy_helper<Real>;
template struct airy_helper<Float>;
template struct airy_helper<Complex>;
template struct airy_helper<Float_complex>;
template struct airy_helper<Object>;

struct func_airy_ai
{
    template<class T>
    auto eval(const T& v) const -> decltype(airy_helper<T>::eval_airy_ai(v))
    {
        return airy_helper<T>::eval_airy_ai(v);
    }
};
struct func_airy_bi
{
    template<class T>
    auto eval(const T& v) const -> decltype(airy_helper<T>::eval_airy_bi(v))
    {
        return airy_helper<T>::eval_airy_bi(v);
    }
};
struct func_airy_ai_dif
{
    template<class T>
    auto eval(const T& v) const -> decltype(airy_helper<T>::eval_airy_ai_dif(v))
    {
        return airy_helper<T>::eval_airy_ai_dif(v);
    }
};
struct func_airy_bi_dif
{
    template<class T>
    auto eval(const T& v) const -> decltype(airy_helper<T>::eval_airy_bi_dif(v))
    {
        return airy_helper<T>::eval_airy_bi_dif(v);
    }
};

struct func_airy_ai_zero
{
    template<class T>
    auto eval(const T& v) const -> Real
    {
        Integer m = mr::converter<Integer,T>::eval(v);
        return airy_helper<Real>::eval_airy_ai_zero(m);
    }
};

struct func_airy_bi_zero
{
    template<class T>
    auto eval(const T& v) const -> Real
    {
        Integer m = mr::converter<Integer,T>::eval(v);
        return airy_helper<Real>::eval_airy_bi_zero(m);
    }
};

};};

namespace matcl
{

Matrix matcl::airy_ai(const Matrix &A)
{
    return eval_scalar_func(A, details::func_airy_ai());
};

Matrix matcl::airy_bi(const Matrix &A)
{
    return eval_scalar_func(A, details::func_airy_bi());
};

Matrix matcl::airy_ai_dif(const Matrix &A)
{
    return eval_scalar_func(A, details::func_airy_ai_dif());
};

Matrix matcl::airy_bi_dif(const Matrix &A)
{
    Matrix m = matcl::convert(A,mat_code::integer_dense);
    return eval_scalar_func(m, details::func_airy_ai_zero());
};
Matrix matcl::airy_ai_zero(const Matrix &A)
{
    Matrix m = matcl::convert(A,mat_code::integer_dense);
    return eval_scalar_func(m, details::func_airy_ai_zero());
};
Matrix matcl::airy_bi_zero(const Matrix &A)
{
    Matrix m = matcl::convert(A,mat_code::integer_dense);
    return eval_scalar_func(m, details::func_airy_bi_zero());
};

};
