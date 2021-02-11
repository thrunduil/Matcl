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

#include "matcl-scalar/lib_functions/func_binary.h"
#include "matcl-core/details/scalfunc_complex.h"

#pragma warning(push)
#pragma warning(disable:4127)   //conditional expression is constant
#include <boost/math/special_functions.hpp>
#include <boost/math/special_functions/next.hpp>
#pragma warning(pop)

namespace matcl { namespace raw { namespace details { namespace scal_func
{
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
};};}}

namespace matcl { namespace details
{

namespace md = matcl::details;

template<>
Real basic_binfunc_helper<Real>::eval_powm1(const Real& x, const Real& y)
{
    //powm1 is broken for x == 0
    if (x == 0.0)
    {
        if (is_nan(y))
            return constants::nan();
        else if (y == 0.0)
            return 0.0;
        else if (y > 0.0)
            return -1.0;
        else if (signbit(x))
            return -constants::inf();
        else
            return constants::inf();
    };

    // powm1 is broken for x == nan, inf
    if (is_finite(x) == false)
        return matcl::pow(x, y) - 1.0;

    if (is_finite(y) == false)
        return matcl::pow(x, y) - 1.0;

    return mrd::scal_func::ignore_errors::powm1(x,y);
}

template<>
Float basic_binfunc_helper<Float>::eval_powm1(const Float& x, const Float& y)
{
    //powm1 is broken for x == 0
    if (x == 0.0f)
    {
        if (is_nan(y))
            return constants::f_nan();
        else if (y == 0.0f)
            return 0.0;
        else if (y > 0.0f)
            return -1.0f;
        else if (signbit(x))
            return -constants::f_inf();
        else
            return constants::f_inf();
    };

    // powm1 is broken for x == nan, inf
    if (is_finite(x) == false)
        return matcl::pow(x, y) - 1.0f;

    if (is_finite(y) == false)
        return matcl::pow(x, y) - 1.0f;

    return mrd::scal_func::ignore_errors::powm1(x,y);
}

template<>
Complex basic_binfunc_helper<Complex>::eval_powm1(const Complex& x, const Complex& y)
{
    if (imag(x) == 0 && imag(y) == 0)
        return basic_binfunc_helper<Real>::eval_powm1(real(x), real(y));

    throw error::function_not_defined_for_complex("powm1");
}

template<>
Float_complex basic_binfunc_helper<Float_complex>::eval_powm1(const Float_complex& x, const Float_complex& y)
{
    if (imag(x) == 0 && imag(y) == 0)
        return basic_binfunc_helper<Float>::eval_powm1(real(x), real(y));

    throw error::function_not_defined_for_complex("powm1");
}

template<>
Real basic_binfunc_helper<Integer>::eval_powm1(const Integer& x, const Integer& y)
{
    return basic_binfunc_helper<Real>::eval_powm1(Real(x),Real(y));
}

template<>
Object basic_binfunc_helper<Object>::eval_powm1(const Object& x, const Object& y)
{
    return dynamic::powm1(x,y);
}

template struct basic_binfunc_helper<Real>;
template struct basic_binfunc_helper<Float>;
template struct basic_binfunc_helper<Complex>;
template struct basic_binfunc_helper<Float_complex>;
template struct basic_binfunc_helper<Integer>;
template struct basic_binfunc_helper<Object>;

};};
