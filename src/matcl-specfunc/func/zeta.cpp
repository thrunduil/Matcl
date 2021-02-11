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

#include "matcl-specfunc/lib_functions/zeta.h"
#include "matcl-matrep/lib_functions/func_unary.h"
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
typename zeta_helper<Integer>::return_type zeta_helper<Integer>::eval_zeta(const Integer& x)
{
    return ignore_errors::zeta((Real)x);
};

template<>
typename zeta_helper<Real>::return_type zeta_helper<Real>::eval_zeta(const Real& x)
{
    return ignore_errors::zeta(x);
};

template<>
typename zeta_helper<Float>::return_type zeta_helper<Float>::eval_zeta(const Float& x)
{
    return ignore_errors::zeta(x);
};

template<>
typename zeta_helper<Complex>::return_type 
zeta_helper<Complex>::eval_zeta(const Complex& x)
{
    if (imag(x) == 0.0)
        return zeta_helper<Real>::eval_zeta(real(x));
    else
        throw matcl::error::function_not_defined_for_complex("zeta"); 
};

template<>
typename zeta_helper<Float_complex>::return_type 
zeta_helper<Float_complex>::eval_zeta(const Float_complex& x)
{
    if (imag(x) == 0.0f)
        return zeta_helper<Float>::eval_zeta(real(x));
    else
        throw matcl::error::function_not_defined_for_complex("zeta"); 
};

template<>
typename zeta_helper<Object>::return_type zeta_helper<Object>::eval_zeta(const Object& x)
{
    return object_func::zeta(x);
};

template struct zeta_helper<Integer>;
template struct zeta_helper<Real>;
template struct zeta_helper<Float>;
template struct zeta_helper<Complex>;
template struct zeta_helper<Float_complex>;
template struct zeta_helper<Object>;

struct func_zeta
{
    template<class T>
    auto eval(const T& v) const -> decltype(zeta_helper<T>::eval_zeta(v))
    {
        return zeta_helper<T>::eval_zeta(v);
    }
};

};};

namespace matcl
{

Matrix matcl::zeta(const Matrix &A)
{
    return eval_scalar_func(A, details::func_zeta());
};

};
