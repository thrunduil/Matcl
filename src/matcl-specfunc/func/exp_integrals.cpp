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

#include "matcl-specfunc/lib_functions/exp_integrals.h"
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

//-------------------------------------------------------------------------------
//                                  REAL
//-------------------------------------------------------------------------------
template<>
typename expint_helper<Real>::return_type 
expint_helper<Real>::eval_expint(const Real& x, Integer n0)
{
    if (n0 < 0)
        return constants::nan();

    unsigned n = unsigned(n0);
    return return_type(ignore_errors::expint(n,x));
}
template<>
typename expint_helper<Real>::return_type 
expint_helper<Real>::eval_expint(const Real& x)
{
    return return_type(ignore_errors::expint(x));
}

//-------------------------------------------------------------------------------
//                                  FLOAT
//-------------------------------------------------------------------------------
template<>
typename expint_helper<Float>::return_type 
expint_helper<Float>::eval_expint(const Float& x, Integer n0)
{
    if (n0 < 0)
        return constants::f_nan();

    unsigned n = unsigned(n0);
    return return_type(ignore_errors::expint(n,x));
}
template<>
typename expint_helper<Float>::return_type 
expint_helper<Float>::eval_expint(const Float& x)
{
    return return_type(ignore_errors::expint(x));
}

//-------------------------------------------------------------------------------
//                                  COMPLEX
//-------------------------------------------------------------------------------
template<>
typename expint_helper<Complex>::return_type 
expint_helper<Complex>::eval_expint(const Complex& x, Integer n0)
{
    if (imag(x) == 0.0)
        return expint_helper<Real>::eval_expint(real(x),n0);
    else
        throw matcl::error::function_not_defined_for_complex("expint"); 
}
template<>
typename expint_helper<Complex>::return_type 
expint_helper<Complex>::eval_expint(const Complex& x)
{
    if (imag(x) == 0.0)
        return expint_helper<Real>::eval_expint(real(x));
    else
        throw matcl::error::function_not_defined_for_complex("expint"); 
}

//-------------------------------------------------------------------------------
//                                  FLOAT_COMPLEX
//-------------------------------------------------------------------------------
template<>
typename expint_helper<Float_complex>::return_type 
expint_helper<Float_complex>::eval_expint(const Float_complex& x, Integer n0)
{
    if (imag(x) == 0.0)
        return expint_helper<Float>::eval_expint(real(x),n0);
    else
        throw matcl::error::function_not_defined_for_complex("expint"); 
}
template<>
typename expint_helper<Float_complex>::return_type 
expint_helper<Float_complex>::eval_expint(const Float_complex& x)
{
    if (imag(x) == 0.0)
        return expint_helper<Float>::eval_expint(real(x));
    else
        throw matcl::error::function_not_defined_for_complex("expint"); 
}

//-------------------------------------------------------------------------------
//                                  INTEGER
//-------------------------------------------------------------------------------
template<>
typename expint_helper<Integer>::return_type 
expint_helper<Integer>::eval_expint(const Integer& x, Integer n0)
{
    return expint_helper<Real>::eval_expint(Real(x),n0);
}
template<>
typename expint_helper<Integer>::return_type 
expint_helper<Integer>::eval_expint(const Integer& x)
{
    return expint_helper<Real>::eval_expint(Real(x));
}

//-------------------------------------------------------------------------------
//                                  OBJECT
//-------------------------------------------------------------------------------
template<>
typename expint_helper<Object>::return_type 
expint_helper<Object>::eval_expint(const Object& x, Integer n0)
{
    return object_func::expint(x,n0);
}
template<>
typename expint_helper<Object>::return_type 
expint_helper<Object>::eval_expint(const Object& x)
{
    return object_func::expint(x);
}

template struct expint_helper<Real>;
template struct expint_helper<Float>;
template struct expint_helper<Complex>;
template struct expint_helper<Float_complex>;
template struct expint_helper<Object>;
template struct expint_helper<Integer>;

struct func_expint
{
    template<class T>
    auto eval(const T& v) const -> decltype(expint_helper<T>::eval_expint(v))
    {
        return expint_helper<T>::eval_expint(v);
    }
};

};};

namespace matcl
{

Matrix matcl::expint(const Matrix &A)
{
    return eval_scalar_func(A, details::func_expint());
};

void test_expint_inst()
{
    Real x = 0.0;
    Integer n = 0;

    expint(x,n);
    expint(x);
};

};
