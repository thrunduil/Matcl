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

#include "matcl-specfunc/lib_functions/hankel.h"
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
typename hankel_helper<Real>::return_type 
hankel_helper<Real>::eval_cyl_hankel_1(const Real& v, const Real& x)
{
    return return_type(ignore_errors::cyl_hankel_1(v,x));
}

template<>
typename hankel_helper<Real>::return_type 
hankel_helper<Real>::eval_cyl_hankel_2(const Real& v, const Real& x)
{
    return return_type(ignore_errors::cyl_hankel_2(v,x));
}

template<>
typename hankel_helper<Real>::return_type 
hankel_helper<Real>::eval_sph_hankel_1(const Real& v, const Real& x)
{
    return return_type(ignore_errors::sph_hankel_1(v,x));
}

template<>
typename hankel_helper<Real>::return_type 
hankel_helper<Real>::eval_sph_hankel_2(const Real& v, const Real& x)
{
    return return_type(ignore_errors::sph_hankel_2(v,x));
}

//-------------------------------------------------------------------------------
//                                  FLOAT
//-------------------------------------------------------------------------------
template<>
typename hankel_helper<Float>::return_type 
hankel_helper<Float>::eval_cyl_hankel_1(const Float& v, const Float& x)
{
    return return_type(ignore_errors::cyl_hankel_1(v,x));
}

template<>
typename hankel_helper<Float>::return_type 
hankel_helper<Float>::eval_cyl_hankel_2(const Float& v, const Float& x)
{
    return return_type(ignore_errors::cyl_hankel_2(v,x));
}

template<>
typename hankel_helper<Float>::return_type 
hankel_helper<Float>::eval_sph_hankel_1(const Float& v, const Float& x)
{
    return return_type(ignore_errors::sph_hankel_1(v,x));
}

template<>
typename hankel_helper<Float>::return_type 
hankel_helper<Float>::eval_sph_hankel_2(const Float& v, const Float& x)
{
    return return_type(ignore_errors::sph_hankel_2(v,x));
}

//-------------------------------------------------------------------------------
//                                  COMPLEX
//-------------------------------------------------------------------------------
template<>
typename hankel_helper<Complex>::return_type 
hankel_helper<Complex>::eval_cyl_hankel_1(const Complex& v, const Complex& x)
{
    if (imag(v) == 0.0 && imag(x) == 0.0)
        return hankel_helper<Real>::eval_cyl_hankel_1(real(v),real(x));
    else
        throw matcl::error::function_not_defined_for_complex("cyl_hankel_1"); 
}

template<>
typename hankel_helper<Complex>::return_type 
hankel_helper<Complex>::eval_cyl_hankel_2(const Complex& v, const Complex& x)
{
    if (imag(v) == 0.0 && imag(x) == 0.0)
        return hankel_helper<Real>::eval_cyl_hankel_2(real(v),real(x));
    else
        throw matcl::error::function_not_defined_for_complex("cyl_hankel_2"); 
}

template<>
typename hankel_helper<Complex>::return_type 
hankel_helper<Complex>::eval_sph_hankel_1(const Complex& v, const Complex& x)
{
    if (imag(v) == 0.0 && imag(x) == 0.0)
        return hankel_helper<Real>::eval_sph_hankel_1(real(v),real(x));
    else
        throw matcl::error::function_not_defined_for_complex("sph_hankel_1"); 
}

template<>
typename hankel_helper<Complex>::return_type 
hankel_helper<Complex>::eval_sph_hankel_2(const Complex& v, const Complex& x)
{
    if (imag(v) == 0.0 && imag(x) == 0.0)
        return hankel_helper<Real>::eval_sph_hankel_2(real(v),real(x));
    else
        throw matcl::error::function_not_defined_for_complex("sph_hankel_2"); 
}

//-------------------------------------------------------------------------------
//                          FLOAT_COMPLEX
//-------------------------------------------------------------------------------
template<>
typename hankel_helper<Float_complex>::return_type 
hankel_helper<Float_complex>::eval_cyl_hankel_1(const Float_complex& v, const Float_complex& x)
{
    if (imag(v) == 0.0 && imag(x) == 0.0)
        return hankel_helper<Float>::eval_cyl_hankel_1(real(v),real(x));
    else
        throw matcl::error::function_not_defined_for_complex("cyl_hankel_1"); 
}

template<>
typename hankel_helper<Float_complex>::return_type 
hankel_helper<Float_complex>::eval_cyl_hankel_2(const Float_complex& v, const Float_complex& x)
{
    if (imag(v) == 0.0 && imag(x) == 0.0)
        return hankel_helper<Float>::eval_cyl_hankel_2(real(v),real(x));
    else
        throw matcl::error::function_not_defined_for_complex("cyl_hankel_2"); 
}

template<>
typename hankel_helper<Float_complex>::return_type 
hankel_helper<Float_complex>::eval_sph_hankel_1(const Float_complex& v, const Float_complex& x)
{
    if (imag(v) == 0.0 && imag(x) == 0.0)
        return hankel_helper<Float>::eval_sph_hankel_1(real(v),real(x));
    else
        throw matcl::error::function_not_defined_for_complex("sph_hankel_1"); 
}

template<>
typename hankel_helper<Float_complex>::return_type 
hankel_helper<Float_complex>::eval_sph_hankel_2(const Float_complex& v, const Float_complex& x)
{
    if (imag(v) == 0.0 && imag(x) == 0.0)
        return hankel_helper<Float>::eval_sph_hankel_2(real(v),real(x));
    else
        throw matcl::error::function_not_defined_for_complex("sph_hankel_2"); 
}

//-------------------------------------------------------------------------------
//                          INTEGER
//-------------------------------------------------------------------------------
template<>
typename hankel_helper<Integer>::return_type 
hankel_helper<Integer>::eval_cyl_hankel_1(const Integer& v, const Integer& x)
{
    return hankel_helper<Real>::eval_cyl_hankel_1((Real)v,(Real)x);
}

template<>
typename hankel_helper<Integer>::return_type 
hankel_helper<Integer>::eval_cyl_hankel_2(const Integer& v, const Integer& x)
{
    return hankel_helper<Real>::eval_cyl_hankel_2((Real)v,(Real)x);
}

template<>
typename hankel_helper<Integer>::return_type 
hankel_helper<Integer>::eval_sph_hankel_1(const Integer& v, const Integer& x)
{
    return hankel_helper<Real>::eval_sph_hankel_1((Real)v,(Real)x);
}

template<>
typename hankel_helper<Integer>::return_type 
hankel_helper<Integer>::eval_sph_hankel_2(const Integer& v, const Integer& x)
{
    return hankel_helper<Real>::eval_sph_hankel_2((Real)v,(Real)x);
}

//-------------------------------------------------------------------------------
//                          OBJECT
//-------------------------------------------------------------------------------
template<>
typename hankel_helper<Object>::return_type 
hankel_helper<Object>::eval_cyl_hankel_1(const Object& v, const Object& x)
{
    return object_func::cyl_hankel_1(v,x);
}

template<>
typename hankel_helper<Object>::return_type 
hankel_helper<Object>::eval_cyl_hankel_2(const Object& v, const Object& x)
{
    return object_func::cyl_hankel_2(v,x);
}

template<>
typename hankel_helper<Object>::return_type 
hankel_helper<Object>::eval_sph_hankel_1(const Object& v, const Object& x)
{
    return object_func::sph_hankel_1(v,x);
}

template<>
typename hankel_helper<Object>::return_type 
hankel_helper<Object>::eval_sph_hankel_2(const Object& v, const Object& x)
{
    return object_func::sph_hankel_2(v,x);
}

template struct hankel_helper<Integer>;
template struct hankel_helper<Real>;
template struct hankel_helper<Float>;
template struct hankel_helper<Complex>;
template struct hankel_helper<Float_complex>;
template struct hankel_helper<Object>;

};};

namespace matcl
{

void test_hankel_inst()
{
    Real v = 0.0;
    Real x = 0.0;

    cyl_hankel_1(v,x);
    cyl_hankel_2(v,x);
    sph_hankel_1(v,x);
    sph_hankel_2(v,x);
};

};
