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

#include "matcl-specfunc/lib_functions/owens_t.h"
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
typename owens_t_helper<Real>::return_type 
owens_t_helper<Real>::eval_owens_t(const Real& h, const Real& a)
{
    return return_type(ignore_errors::owens_t(h,a));
}

//-------------------------------------------------------------------------------
//                                  FLOAT
//-------------------------------------------------------------------------------
template<>
typename owens_t_helper<Float>::return_type 
owens_t_helper<Float>::eval_owens_t(const Float& h, const Float& a)
{
    return return_type(ignore_errors::owens_t(h,a));
}

//-------------------------------------------------------------------------------
//                                  COMPLEX
//-------------------------------------------------------------------------------
template<>
typename owens_t_helper<Complex>::return_type 
owens_t_helper<Complex>::eval_owens_t(const Complex& h, const Complex& a)
{
    if (imag(h) == 0.0 && imag(a) == 0.0)
        return owens_t_helper<Real>::eval_owens_t(real(h),real(a));
    else
        throw matcl::error::function_not_defined_for_complex("owens_t"); 
}

//-------------------------------------------------------------------------------
//                                  FLOAT_COMPLEX
//-------------------------------------------------------------------------------
template<>
typename owens_t_helper<Float_complex>::return_type 
owens_t_helper<Float_complex>::eval_owens_t(const Float_complex& h, const Float_complex& a)
{
    if (imag(h) == 0.0 && imag(a) == 0.0)
        return owens_t_helper<Float>::eval_owens_t(real(h),real(a));
    else
        throw matcl::error::function_not_defined_for_complex("owens_t"); 
}

//-------------------------------------------------------------------------------
//                                  OBJECT
//-------------------------------------------------------------------------------
template<>
typename owens_t_helper<Object>::return_type 
owens_t_helper<Object>::eval_owens_t(const Object& h, const Object& a)
{
    return object_func::owens_t(h,a);
}

template struct owens_t_helper<Real>;
template struct owens_t_helper<Float>;
template struct owens_t_helper<Complex>;
template struct owens_t_helper<Float_complex>;
template struct owens_t_helper<Object>;

};};

namespace matcl
{

void test_owens_t_inst()
{
    Real h = 0.0;
    Real a = 0.0;

    owens_t(h,a);
};

};
