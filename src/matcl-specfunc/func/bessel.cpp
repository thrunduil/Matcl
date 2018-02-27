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

#include "matcl-specfunc/lib_functions/bessel.h"
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
typename bessel_helper<Real>::return_type
bessel_helper<Real>::eval_cyl_bessel_j(const Real& v, const Real& x)
{
    return ignore_errors::cyl_bessel_j(v,x);
}

template<>
typename bessel_helper<Real>::return_type
bessel_helper<Real>::eval_cyl_neumann(const Real& v, const Real& x)
{
    return ignore_errors::cyl_neumann(v,x);
}

template<>
typename bessel_helper<Real>::return_type 
bessel_helper<Real>::eval_cyl_bessel_j_zero(const Real& v, Integer m)
{
    return ignore_errors::cyl_bessel_j_zero(v,m);
}

template<>
typename bessel_helper<Real>::return_type 
bessel_helper<Real>::eval_cyl_neumann_zero(const Real& v, Integer m)
{
    return ignore_errors::cyl_neumann_zero(v,m);
}

template<>
typename bessel_helper<Real>::return_type
bessel_helper<Real>::eval_cyl_bessel_j_dif(const Real& v, const Real& x)
{
    return ignore_errors::cyl_bessel_j_prime(v,x);
}

template<>
typename bessel_helper<Real>::return_type
bessel_helper<Real>::eval_cyl_neumann_dif(const Real& v, const Real& x)
{
    return ignore_errors::cyl_neumann_prime(v,x);
}

template<>
typename bessel_helper<Real>::return_type 
bessel_helper<Real>::eval_cyl_bessel_i(const Real& v, const Real& x)
{
    return ignore_errors::cyl_bessel_i(v,x);
}

template<>
typename bessel_helper<Real>::return_type
bessel_helper<Real>::eval_cyl_bessel_k(const Real& v, const Real& x)
{
    return ignore_errors::cyl_bessel_k(v,x);
}

template<>
typename bessel_helper<Real>::return_type
bessel_helper<Real>::eval_cyl_bessel_i_dif(const Real& v, const Real& x)
{
    return ignore_errors::cyl_bessel_i_prime(v,x);
}

template<>
typename bessel_helper<Real>::return_type
bessel_helper<Real>::eval_cyl_bessel_k_dif(const Real& v, const Real& x)
{
    return ignore_errors::cyl_bessel_k_prime(v,x);
}

template<>
typename bessel_helper<Real>::return_type
bessel_helper<Real>::eval_sph_bessel(Integer n, const Real& x)
{
    //convert negative index to equivalent positive index
    unsigned nu = (n >= 0) ? (unsigned)n : unsigned(-n - 1);
    return ignore_errors::sph_bessel(nu,x);
}

template<>
typename bessel_helper<Real>::return_type
bessel_helper<Real>::eval_sph_neumann(Integer n, const Real& x)
{
    //convert negative index to equivalent positive index
    unsigned nu = (n >= 0) ? (unsigned)n : unsigned(-n - 1);
    return ignore_errors::sph_neumann(nu,x);
}

template<>
typename bessel_helper<Real>::return_type
bessel_helper<Real>::eval_sph_bessel_dif(Integer n, const Real& x)
{
    //convert negative index to equivalent positive index
    unsigned nu = (n >= 0) ? (unsigned)n : unsigned(-n - 1);
    return ignore_errors::sph_bessel_prime(nu,x);
}

template<>
typename bessel_helper<Real>::return_type
bessel_helper<Real>::eval_sph_neumann_dif(Integer n, const Real& x)
{
    //convert negative index to equivalent positive index
    unsigned nu = (n >= 0) ? (unsigned)n : unsigned(-n - 1);
    return ignore_errors::sph_neumann_prime(nu,x);
}

//-------------------------------------------------------------------------------
//                                  FLOAT
//-------------------------------------------------------------------------------
template<>
typename bessel_helper<Float>::return_type
bessel_helper<Float>::eval_cyl_bessel_j(const Float& v, const Float& x)
{
    return ignore_errors::cyl_bessel_j(v,x);
}

template<>
typename bessel_helper<Float>::return_type
bessel_helper<Float>::eval_cyl_neumann(const Float& v, const Float& x)
{
    return ignore_errors::cyl_neumann(v,x);
}

template<>
typename bessel_helper<Float>::return_type 
bessel_helper<Float>::eval_cyl_bessel_j_zero(const Float& v, Integer m)
{
    return ignore_errors::cyl_bessel_j_zero(v,m);
}

template<>
typename bessel_helper<Float>::return_type 
bessel_helper<Float>::eval_cyl_neumann_zero(const Float& v, Integer m)
{
    return ignore_errors::cyl_neumann_zero(v,m);
}

template<>
typename bessel_helper<Float>::return_type
bessel_helper<Float>::eval_cyl_bessel_j_dif(const Float& v, const Float& x)
{
    return ignore_errors::cyl_bessel_j_prime(v,x);
}

template<>
typename bessel_helper<Float>::return_type
bessel_helper<Float>::eval_cyl_neumann_dif(const Float& v, const Float& x)
{
    return ignore_errors::cyl_neumann_prime(v,x);
}

template<>
typename bessel_helper<Float>::return_type 
bessel_helper<Float>::eval_cyl_bessel_i(const Float& v, const Float& x)
{
    return ignore_errors::cyl_bessel_i(v,x);
}

template<>
typename bessel_helper<Float>::return_type
bessel_helper<Float>::eval_cyl_bessel_k(const Float& v, const Float& x)
{
    return ignore_errors::cyl_bessel_k(v,x);
}

template<>
typename bessel_helper<Float>::return_type
bessel_helper<Float>::eval_cyl_bessel_i_dif(const Float& v, const Float& x)
{
    return ignore_errors::cyl_bessel_i_prime(v,x);
}

template<>
typename bessel_helper<Float>::return_type
bessel_helper<Float>::eval_cyl_bessel_k_dif(const Float& v, const Float& x)
{
    return ignore_errors::cyl_bessel_k_prime(v,x);
}

template<>
typename bessel_helper<Float>::return_type
bessel_helper<Float>::eval_sph_bessel(Integer n, const Float& x)
{
    //convert negative index to equivalent positive index
    unsigned nu = (n >= 0) ? (unsigned)n : unsigned(-n - 1);
    return ignore_errors::sph_bessel(nu,x);
}

template<>
typename bessel_helper<Float>::return_type
bessel_helper<Float>::eval_sph_neumann(Integer n, const Float& x)
{
    //convert negative index to equivalent positive index
    unsigned nu = (n >= 0) ? (unsigned)n : unsigned(-n - 1);
    return ignore_errors::sph_neumann(nu,x);
}

template<>
typename bessel_helper<Float>::return_type
bessel_helper<Float>::eval_sph_bessel_dif(Integer n, const Float& x)
{
    //convert negative index to equivalent positive index
    unsigned nu = (n >= 0) ? (unsigned)n : unsigned(-n - 1);
    return ignore_errors::sph_bessel_prime(nu,x);
}

template<>
typename bessel_helper<Float>::return_type
bessel_helper<Float>::eval_sph_neumann_dif(Integer n, const Float& x)
{
    //convert negative index to equivalent positive index
    unsigned nu = (n >= 0) ? (unsigned)n : unsigned(-n - 1);
    return ignore_errors::sph_neumann_prime(nu,x);
}

//-------------------------------------------------------------------------------
//                                  COMPLEX
//-------------------------------------------------------------------------------
template<>
typename bessel_helper<Complex>::return_type
bessel_helper<Complex>::eval_cyl_bessel_j(const Complex& v, const Complex& x)
{
    if (imag(v) == 0.0 && imag(x) == 0.0)
        return bessel_helper<Real>::eval_cyl_bessel_j(real(v),real(x));
    else
        throw matcl::error::function_not_defined_for_complex("cyl_bessel_j"); 
}

template<>
typename bessel_helper<Complex>::return_type
bessel_helper<Complex>::eval_cyl_neumann(const Complex& v, const Complex& x)
{
    if (imag(v) == 0.0 && imag(x) == 0.0)
        return bessel_helper<Real>::eval_cyl_neumann(real(v),real(x));
    else
        throw matcl::error::function_not_defined_for_complex("cyl_neumann"); 
}

template<>
typename bessel_helper<Complex>::return_type 
bessel_helper<Complex>::eval_cyl_bessel_j_zero(const Complex& v, Integer m)
{
    if (imag(v) == 0.0)
        return bessel_helper<Real>::eval_cyl_bessel_j_zero(real(v),m);
    else
        throw matcl::error::function_not_defined_for_complex("cyl_bessel_j_zero"); 
}

template<>
typename bessel_helper<Complex>::return_type 
bessel_helper<Complex>::eval_cyl_neumann_zero(const Complex& v, Integer m)
{
    if (imag(v) == 0.0)
        return bessel_helper<Real>::eval_cyl_neumann_zero(real(v),m);
    else
        throw matcl::error::function_not_defined_for_complex("cyl_neumann_zero"); 
}

template<>
typename bessel_helper<Complex>::return_type
bessel_helper<Complex>::eval_cyl_bessel_j_dif(const Complex& v, const Complex& x)
{
    if (imag(v) == 0.0 && imag(x) == 0.0)
        return bessel_helper<Real>::eval_cyl_bessel_j_dif(real(v),real(x));
    else
        throw matcl::error::function_not_defined_for_complex("cyl_bessel_j_dif"); 
}

template<>
typename bessel_helper<Complex>::return_type
bessel_helper<Complex>::eval_cyl_neumann_dif(const Complex& v, const Complex& x)
{
    if (imag(v) == 0.0 && imag(x) == 0.0)
        return bessel_helper<Real>::eval_cyl_neumann_dif(real(v),real(x));
    else
        throw matcl::error::function_not_defined_for_complex("cyl_neumann_dif"); 
}

template<>
typename bessel_helper<Complex>::return_type 
bessel_helper<Complex>::eval_cyl_bessel_i(const Complex& v, const Complex& x)
{
    if (imag(v) == 0.0 && imag(x) == 0.0)
        return bessel_helper<Real>::eval_cyl_bessel_i(real(v),real(x));
    else
        throw matcl::error::function_not_defined_for_complex("cyl_bessel_i"); 
}

template<>
typename bessel_helper<Complex>::return_type
bessel_helper<Complex>::eval_cyl_bessel_k(const Complex& v, const Complex& x)
{
    if (imag(v) == 0.0 && imag(x) == 0.0)
        return bessel_helper<Real>::eval_cyl_bessel_k(real(v),real(x));
    else
        throw matcl::error::function_not_defined_for_complex("cyl_bessel_k"); 
}

template<>
typename bessel_helper<Complex>::return_type
bessel_helper<Complex>::eval_cyl_bessel_i_dif(const Complex& v, const Complex& x)
{
    if (imag(v) == 0.0 && imag(x) == 0.0)
        return bessel_helper<Real>::eval_cyl_bessel_i_dif(real(v),real(x));
    else
        throw matcl::error::function_not_defined_for_complex("cyl_bessel_i_dif"); 
}

template<>
typename bessel_helper<Complex>::return_type
bessel_helper<Complex>::eval_cyl_bessel_k_dif(const Complex& v, const Complex& x)
{
    if (imag(v) == 0.0 && imag(x) == 0.0)
        return bessel_helper<Real>::eval_cyl_bessel_k_dif(real(v),real(x));
    else
        throw matcl::error::function_not_defined_for_complex("cyl_bessel_k_dif"); 
}

template<>
typename bessel_helper<Complex>::return_type
bessel_helper<Complex>::eval_sph_bessel(Integer n, const Complex& x)
{
    if (imag(x) == 0.0)
        return bessel_helper<Real>::eval_sph_bessel(n,real(x));
    else
        throw matcl::error::function_not_defined_for_complex("sph_bessel"); 
}

template<>
typename bessel_helper<Complex>::return_type
bessel_helper<Complex>::eval_sph_neumann(Integer n, const Complex& x)
{
    if (imag(x) == 0.0)
        return bessel_helper<Real>::eval_sph_neumann(n,real(x));
    else
        throw matcl::error::function_not_defined_for_complex("sph_neumann"); 
}

template<>
typename bessel_helper<Complex>::return_type
bessel_helper<Complex>::eval_sph_bessel_dif(Integer n, const Complex& x)
{
    if (imag(x) == 0.0)
        return bessel_helper<Real>::eval_sph_bessel_dif(n,real(x));
    else
        throw matcl::error::function_not_defined_for_complex("sph_bessel_dif"); 
}

template<>
typename bessel_helper<Complex>::return_type
bessel_helper<Complex>::eval_sph_neumann_dif(Integer n, const Complex& x)
{
    if (imag(x) == 0.0)
        return bessel_helper<Real>::eval_sph_neumann_dif(n,real(x));
    else
        throw matcl::error::function_not_defined_for_complex("sph_neumann_dif"); 
}

//-------------------------------------------------------------------------------
//                          FLOAT_COMPLEX
//-------------------------------------------------------------------------------
template<>
typename bessel_helper<Float_complex>::return_type
bessel_helper<Float_complex>::eval_cyl_bessel_j(const Float_complex& v, const Float_complex& x)
{
    if (imag(v) == 0.0 && imag(x) == 0.0)
        return bessel_helper<Float>::eval_cyl_bessel_j(real(v),real(x));
    else
        throw matcl::error::function_not_defined_for_complex("cyl_bessel_j"); 
}

template<>
typename bessel_helper<Float_complex>::return_type
bessel_helper<Float_complex>::eval_cyl_neumann(const Float_complex& v, const Float_complex& x)
{
    if (imag(v) == 0.0 && imag(x) == 0.0)
        return bessel_helper<Float>::eval_cyl_neumann(real(v),real(x));
    else
        throw matcl::error::function_not_defined_for_complex("cyl_neumann"); 
}

template<>
typename bessel_helper<Float_complex>::return_type 
bessel_helper<Float_complex>::eval_cyl_bessel_j_zero(const Float_complex& v, Integer m)
{
    if (imag(v) == 0.0)
        return bessel_helper<Float>::eval_cyl_bessel_j_zero(real(v),m);
    else
        throw matcl::error::function_not_defined_for_complex("cyl_bessel_j_zero"); 
}

template<>
typename bessel_helper<Float_complex>::return_type 
bessel_helper<Float_complex>::eval_cyl_neumann_zero(const Float_complex& v, Integer m)
{
    if (imag(v) == 0.0)
        return bessel_helper<Float>::eval_cyl_neumann_zero(real(v),m);
    else
        throw matcl::error::function_not_defined_for_complex("cyl_neumann_zero"); 
}

template<>
typename bessel_helper<Float_complex>::return_type
bessel_helper<Float_complex>::eval_cyl_bessel_j_dif(const Float_complex& v, const Float_complex& x)
{
    if (imag(v) == 0.0 && imag(x) == 0.0)
        return bessel_helper<Float>::eval_cyl_bessel_j_dif(real(v),real(x));
    else
        throw matcl::error::function_not_defined_for_complex("cyl_bessel_j_dif"); 
}

template<>
typename bessel_helper<Float_complex>::return_type
bessel_helper<Float_complex>::eval_cyl_neumann_dif(const Float_complex& v, const Float_complex& x)
{
    if (imag(v) == 0.0 && imag(x) == 0.0)
        return bessel_helper<Float>::eval_cyl_neumann_dif(real(v),real(x));
    else
        throw matcl::error::function_not_defined_for_complex("cyl_neumann_dif"); 
}

template<>
typename bessel_helper<Float_complex>::return_type 
bessel_helper<Float_complex>::eval_cyl_bessel_i(const Float_complex& v, const Float_complex& x)
{
    if (imag(v) == 0.0 && imag(x) == 0.0)
        return bessel_helper<Float>::eval_cyl_bessel_i(real(v),real(x));
    else
        throw matcl::error::function_not_defined_for_complex("cyl_bessel_i"); 
}

template<>
typename bessel_helper<Float_complex>::return_type
bessel_helper<Float_complex>::eval_cyl_bessel_k(const Float_complex& v, const Float_complex& x)
{
    if (imag(v) == 0.0 && imag(x) == 0.0)
        return bessel_helper<Float>::eval_cyl_bessel_k(real(v),real(x));
    else
        throw matcl::error::function_not_defined_for_complex("cyl_bessel_k"); 
}

template<>
typename bessel_helper<Float_complex>::return_type
bessel_helper<Float_complex>::eval_cyl_bessel_i_dif(const Float_complex& v, const Float_complex& x)
{
    if (imag(v) == 0.0 && imag(x) == 0.0)
        return bessel_helper<Float>::eval_cyl_bessel_i_dif(real(v),real(x));
    else
        throw matcl::error::function_not_defined_for_complex("cyl_bessel_i_dif"); 
}

template<>
typename bessel_helper<Float_complex>::return_type
bessel_helper<Float_complex>::eval_cyl_bessel_k_dif(const Float_complex& v, const Float_complex& x)
{
    if (imag(v) == 0.0 && imag(x) == 0.0)
        return bessel_helper<Float>::eval_cyl_bessel_k_dif(real(v),real(x));
    else
        throw matcl::error::function_not_defined_for_complex("cyl_bessel_k_dif"); 
}

template<>
typename bessel_helper<Float_complex>::return_type
bessel_helper<Float_complex>::eval_sph_bessel(Integer n, const Float_complex& x)
{
    if (imag(x) == 0.0)
        return bessel_helper<Float>::eval_sph_bessel(n,real(x));
    else
        throw matcl::error::function_not_defined_for_complex("sph_bessel"); 
}

template<>
typename bessel_helper<Float_complex>::return_type
bessel_helper<Float_complex>::eval_sph_neumann(Integer n, const Float_complex& x)
{
    if (imag(x) == 0.0)
        return bessel_helper<Float>::eval_sph_neumann(n,real(x));
    else
        throw matcl::error::function_not_defined_for_complex("sph_neumann"); 
}

template<>
typename bessel_helper<Float_complex>::return_type
bessel_helper<Float_complex>::eval_sph_bessel_dif(Integer n, const Float_complex& x)
{
    if (imag(x) == 0.0)
        return bessel_helper<Float>::eval_sph_bessel_dif(n,real(x));
    else
        throw matcl::error::function_not_defined_for_complex("sph_bessel_dif"); 
}

template<>
typename bessel_helper<Float_complex>::return_type
bessel_helper<Float_complex>::eval_sph_neumann_dif(Integer n, const Float_complex& x)
{
    if (imag(x) == 0.0)
        return bessel_helper<Float>::eval_sph_neumann_dif(n,real(x));
    else
        throw matcl::error::function_not_defined_for_complex("sph_neumann_dif"); 
}

//-------------------------------------------------------------------------------
//                                  INTEGER
//-------------------------------------------------------------------------------
template<>
typename bessel_helper<Integer>::return_type
bessel_helper<Integer>::eval_cyl_bessel_j(const Integer& v, const Integer& x)
{
    return bessel_helper<Real>::eval_cyl_bessel_j((Real)v,(Real)x);
}

template<>
typename bessel_helper<Integer>::return_type
bessel_helper<Integer>::eval_cyl_neumann(const Integer& v, const Integer& x)
{
    return bessel_helper<Real>::eval_cyl_neumann((Real)v,(Real)x);
}

template<>
typename bessel_helper<Integer>::return_type 
bessel_helper<Integer>::eval_cyl_bessel_j_zero(const Integer& v, Integer m)
{
    return bessel_helper<Real>::eval_cyl_bessel_j_zero((Real)v,m);
}

template<>
typename bessel_helper<Integer>::return_type 
bessel_helper<Integer>::eval_cyl_neumann_zero(const Integer& v, Integer m)
{
    return bessel_helper<Real>::eval_cyl_neumann_zero((Real)v,m);
}

template<>
typename bessel_helper<Integer>::return_type
bessel_helper<Integer>::eval_cyl_bessel_j_dif(const Integer& v, const Integer& x)
{
    return bessel_helper<Real>::eval_cyl_bessel_j_dif((Real)v,(Real)x);
}

template<>
typename bessel_helper<Integer>::return_type
bessel_helper<Integer>::eval_cyl_neumann_dif(const Integer& v, const Integer& x)
{
    return bessel_helper<Real>::eval_cyl_neumann_dif((Real)v,(Real)x);
}

template<>
typename bessel_helper<Integer>::return_type 
bessel_helper<Integer>::eval_cyl_bessel_i(const Integer& v, const Integer& x)
{
    return bessel_helper<Real>::eval_cyl_bessel_i((Real)v,(Real)x);
}

template<>
typename bessel_helper<Integer>::return_type
bessel_helper<Integer>::eval_cyl_bessel_k(const Integer& v, const Integer& x)
{
    return bessel_helper<Real>::eval_cyl_bessel_k((Real)v,(Real)x);
}

template<>
typename bessel_helper<Integer>::return_type
bessel_helper<Integer>::eval_cyl_bessel_i_dif(const Integer& v, const Integer& x)
{
    return bessel_helper<Real>::eval_cyl_bessel_i_dif((Real)v,(Real)x);
}

template<>
typename bessel_helper<Integer>::return_type
bessel_helper<Integer>::eval_cyl_bessel_k_dif(const Integer& v, const Integer& x)
{
    return bessel_helper<Real>::eval_cyl_bessel_k_dif((Real)v,(Real)x);
}

template<>
typename bessel_helper<Integer>::return_type
bessel_helper<Integer>::eval_sph_bessel(Integer n, const Integer& x)
{
    return bessel_helper<Real>::eval_sph_bessel(n,(Real)x);
}

template<>
typename bessel_helper<Integer>::return_type
bessel_helper<Integer>::eval_sph_neumann(Integer n, const Integer& x)
{
    return bessel_helper<Real>::eval_sph_neumann(n,(Real)x);
}

template<>
typename bessel_helper<Integer>::return_type
bessel_helper<Integer>::eval_sph_bessel_dif(Integer n, const Integer& x)
{
    return bessel_helper<Real>::eval_sph_bessel_dif(n,(Real)x);
}

template<>
typename bessel_helper<Integer>::return_type
bessel_helper<Integer>::eval_sph_neumann_dif(Integer n, const Integer& x)
{
    return bessel_helper<Real>::eval_sph_neumann_dif(n,(Real)x);
}

//-------------------------------------------------------------------------------
//                                  OBJECT
//-------------------------------------------------------------------------------
template<>
typename bessel_helper<Object>::return_type
bessel_helper<Object>::eval_cyl_bessel_j(const Object& v, const Object& x)
{
    return object_func::cyl_bessel_j(v,x);
}

template<>
typename bessel_helper<Object>::return_type
bessel_helper<Object>::eval_cyl_neumann(const Object& v, const Object& x)
{
    return object_func::cyl_neumann(v,x);
}

template<>
typename bessel_helper<Object>::return_type 
bessel_helper<Object>::eval_cyl_bessel_j_zero(const Object& v, Integer m)
{
    return object_func::cyl_bessel_j_zero(v,m);
}

template<>
typename bessel_helper<Object>::return_type 
bessel_helper<Object>::eval_cyl_neumann_zero(const Object& v, Integer m)
{
    return object_func::cyl_neumann_zero(v,m);
}

template<>
typename bessel_helper<Object>::return_type
bessel_helper<Object>::eval_cyl_bessel_j_dif(const Object& v, const Object& x)
{
    return object_func::cyl_bessel_j_dif(v,x);
}

template<>
typename bessel_helper<Object>::return_type
bessel_helper<Object>::eval_cyl_neumann_dif(const Object& v, const Object& x)
{
    return object_func::cyl_neumann_dif(v,x);
}

template<>
typename bessel_helper<Object>::return_type 
bessel_helper<Object>::eval_cyl_bessel_i(const Object& v, const Object& x)
{
    return object_func::cyl_bessel_i(v,x);
}

template<>
typename bessel_helper<Object>::return_type
bessel_helper<Object>::eval_cyl_bessel_k(const Object& v, const Object& x)
{
    return object_func::cyl_bessel_k(v,x);
}

template<>
typename bessel_helper<Object>::return_type
bessel_helper<Object>::eval_cyl_bessel_i_dif(const Object& v, const Object& x)
{
    return object_func::cyl_bessel_i_dif(v,x);
}

template<>
typename bessel_helper<Object>::return_type
bessel_helper<Object>::eval_cyl_bessel_k_dif(const Object& v, const Object& x)
{
    return object_func::cyl_bessel_k_dif(v,x);
}

template<>
typename bessel_helper<Object>::return_type
bessel_helper<Object>::eval_sph_bessel(Integer n, const Object& x)
{
    //convert negative index to equivalent positive index
    unsigned nu = (n >= 0) ? (unsigned)n : unsigned(-n - 1);
    return object_func::sph_bessel(nu,x);
}

template<>
typename bessel_helper<Object>::return_type
bessel_helper<Object>::eval_sph_neumann(Integer n, const Object& x)
{
    //convert negative index to equivalent positive index
    unsigned nu = (n >= 0) ? (unsigned)n : unsigned(-n - 1);
    return object_func::sph_neumann(nu,x);
}

template<>
typename bessel_helper<Object>::return_type
bessel_helper<Object>::eval_sph_bessel_dif(Integer n, const Object& x)
{
    //convert negative index to equivalent positive index
    unsigned nu = (n >= 0) ? (unsigned)n : unsigned(-n - 1);
    return object_func::sph_bessel_dif(nu,x);
}

template<>
typename bessel_helper<Object>::return_type
bessel_helper<Object>::eval_sph_neumann_dif(Integer n, const Object& x)
{
    //convert negative index to equivalent positive index
    unsigned nu = (n >= 0) ? (unsigned)n : unsigned(-n - 1);
    return object_func::sph_neumann_dif(nu,x);
}

template struct bessel_helper<Integer>;
template struct bessel_helper<Real>;
template struct bessel_helper<Float>;
template struct bessel_helper<Complex>;
template struct bessel_helper<Float_complex>;
template struct bessel_helper<Object>;

};};

namespace matcl
{

void test_bessel_inst()
{
    Real v = 0.0;
    Real x = 0.0;
    Integer n = 1;

    cyl_bessel_j(v,x);
    cyl_neumann(v,x);
    cyl_bessel_j_zero(v,1);
    cyl_neumann_zero(v,1);
    cyl_bessel_j_dif(v,x);
    cyl_neumann_dif(v,x);

    cyl_bessel_i(v,x);
    cyl_bessel_k(v,x);
    cyl_bessel_i_dif(v,x);
    cyl_bessel_k_dif(v,x);

    sph_bessel(n,x);
    sph_neumann(n,x);
    sph_bessel_dif(n,x);
    sph_neumann_dif(n,x);
};

};
