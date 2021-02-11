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

#include "matcl-specfunc/lib_functions/error.h"
#include "matcl-matrep/lib_functions/func_unary.h"
#include "matcl-specfunc/func/raw/complex_func.h"
#include "matcl-specfunc/objects/object_func.h"

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

template<>
typename error_helper<Integer>::return_type error_helper<Integer>::eval_erf(const Integer& x)
{
    return ignore_errors::erf((Real)x);
};
template<>
typename error_helper<Integer>::return_type error_helper<Integer>::eval_erfc(const Integer& x)
{
    return ignore_errors::erfc((Real)x);
};
template<>
typename error_helper<Integer>::return_type error_helper<Integer>::eval_erf_inv(const Integer& x)
{
    return ignore_errors::erf_inv((Real)x);
};
template<>
typename error_helper<Integer>::return_type error_helper<Integer>::eval_erfc_inv(const Integer& x)
{
    return ignore_errors::erfc_inv((Real)x);
};

template<>
typename error_helper<Real>::return_type error_helper<Real>::eval_erf(const Real& x)
{
    return ignore_errors::erf(x);
};
template<>
typename error_helper<Real>::return_type error_helper<Real>::eval_erfc(const Real& x)
{
    return ignore_errors::erfc(x);
};
template<>
typename error_helper<Real>::return_type error_helper<Real>::eval_erf_inv(const Real& x)
{
    return ignore_errors::erf_inv(x);
};
template<>
typename error_helper<Real>::return_type error_helper<Real>::eval_erfc_inv(const Real& x)
{
    return ignore_errors::erfc_inv(x);
};

template<>
typename error_helper<Float>::return_type error_helper<Float>::eval_erf(const Float& x)
{
    return ignore_errors::erf(x);
};
template<>
typename error_helper<Float>::return_type error_helper<Float>::eval_erfc(const Float& x)
{
    return ignore_errors::erfc(x);
};
template<>
typename error_helper<Float>::return_type error_helper<Float>::eval_erf_inv(const Float& x)
{
    return ignore_errors::erf_inv(x);
};
template<>
typename error_helper<Float>::return_type error_helper<Float>::eval_erfc_inv(const Float& x)
{
    return ignore_errors::erfc_inv(x);
};

template<>
typename error_helper<Complex>::return_type 
error_helper<Complex>::eval_erf(const Complex& x)
{
    return scal_func::impl_erf(x);
};
template<>
typename error_helper<Complex>::return_type 
error_helper<Complex>::eval_erfc(const Complex& x)
{
    return scal_func::impl_erfc(x);
};
template<>
typename error_helper<Complex>::return_type 
error_helper<Complex>::eval_erf_inv(const Complex& x)
{
    if (imag(x) == 0.0)
        return error_helper<Real>::eval_erf_inv(real(x));
    else
        throw matcl::error::function_not_defined_for_complex("erf_inv"); 
};
template<>
typename error_helper<Complex>::return_type 
error_helper<Complex>::eval_erfc_inv(const Complex& x)
{
    if (imag(x) == 0.0)
        return error_helper<Real>::eval_erfc_inv(real(x));
    else
        throw matcl::error::function_not_defined_for_complex("erfc_inv"); 
};

template<>
typename error_helper<Float_complex>::return_type 
error_helper<Float_complex>::eval_erf(const Float_complex& x)
{
    return scal_func::impl_erf(x);
};
template<>
typename error_helper<Float_complex>::return_type 
error_helper<Float_complex>::eval_erfc(const Float_complex& x)
{
    return scal_func::impl_erfc(x);
};
template<>
typename error_helper<Float_complex>::return_type 
error_helper<Float_complex>::eval_erf_inv(const Float_complex& x)
{
    if (imag(x) == 0.0f)
        return error_helper<Float>::eval_erf_inv(real(x));
    else
        throw matcl::error::function_not_defined_for_complex("erf_inv"); 
};
template<>
typename error_helper<Float_complex>::return_type 
error_helper<Float_complex>::eval_erfc_inv(const Float_complex& x)
{
    if (imag(x) == 0.0)
        return error_helper<Float>::eval_erfc_inv(real(x));
    else
        throw matcl::error::function_not_defined_for_complex("erfc_inv"); 
};

template<>
typename error_helper<Object>::return_type error_helper<Object>::eval_erf(const Object& x)
{
    return object_func::erf(x);
};
template<>
typename error_helper<Object>::return_type error_helper<Object>::eval_erfc(const Object& x)
{
    return object_func::erfc(x);
};
template<>
typename error_helper<Object>::return_type error_helper<Object>::eval_erf_inv(const Object& x)
{
    return object_func::erf_inv(x);
};
template<>
typename error_helper<Object>::return_type error_helper<Object>::eval_erfc_inv(const Object& x)
{
    return object_func::erfc_inv(x);
};

template struct error_helper<Integer>;
template struct error_helper<Real>;
template struct error_helper<Float>;
template struct error_helper<Complex>;
template struct error_helper<Float_complex>;
template struct error_helper<Object>;

struct func_erf
{
    template<class T>
    auto eval(const T& v) const -> decltype(error_helper<T>::eval_erf(v))
    {
        return error_helper<T>::eval_erf(v);
    }
};

struct func_erfc
{
    template<class T>
    auto eval(const T& v) const -> decltype(error_helper<T>::eval_erfc(v))
    {
        return error_helper<T>::eval_erfc(v);
    }
};

struct func_erf_inv
{
    template<class T>
    auto eval(const T& v) const -> decltype(error_helper<T>::eval_erf_inv(v))
    {
        return error_helper<T>::eval_erf_inv(v);
    }
};

struct func_erfc_inv
{
    template<class T>
    auto eval(const T& v) const -> decltype(error_helper<T>::eval_erfc_inv(v))
    {
        return error_helper<T>::eval_erfc_inv(v);
    }
};

};};

namespace matcl
{

Matrix matcl::erf(const Matrix &A)
{
    return eval_scalar_func(A, details::func_erf());
};
Matrix matcl::erfc(const Matrix &A)
{
    return eval_scalar_func(A, details::func_erfc());
};
Matrix matcl::erf_inv(const Matrix &A)
{
    return eval_scalar_func(A, details::func_erf_inv());
};
Matrix matcl::erfc_inv(const Matrix &A)
{
    return eval_scalar_func(A, details::func_erfc_inv());
};

};
