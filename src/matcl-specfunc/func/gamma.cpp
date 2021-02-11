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

#include "matcl-specfunc/lib_functions/gamma.h"
#include "matcl-matrep/lib_functions/func_unary.h"
#include "matcl-specfunc/func/raw/complex_func.h"
#include "matcl-specfunc/objects/object_func.h"

#pragma warning(push)
#pragma warning(disable:4127)   //conditional expression is constant
#pragma warning(disable:4244)   // :  '=': conversion from 'double' to 'T', possible loss of data
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
typename gamma_helper<Real>::return_type gamma_helper<Real>::eval_gamma(const Real& x)
{
    if(mrd::scal_func::isnan(x)) 
        return x;
    else if(mrd::scal_func::isinf(x))
        return x > 0. ? constants::inf() : constants::nan();

    return ignore_errors::tgamma(x);
};
template<>
typename gamma_helper<Real>::return_type gamma_helper<Real>::eval_gammaln(const Real& x)
{
    if(mrd::scal_func::isnan(x))
        return x;
    else if(mrd::scal_func::isinf(x))
        return x > 0. ? constants::inf() : constants::nan();

    Integer sign;
    Real val = ignore_errors::lgamma(x,&sign);

    if (sign < 0)
        return constants::nan();
    else
        return val;
};

template<>
typename gamma_helper<Real>::return_type 
gamma_helper<Real>::eval_gamma1pm1(const Real& x)
{
    return ignore_errors::tgamma1pm1(x);
};
template<>
typename gamma_helper<Real>::return_type 
gamma_helper<Real>::eval_digamma(const Real& x)
{
    return ignore_errors::digamma(x);
};
template<>
typename gamma_helper<Real>::return_type 
gamma_helper<Real>::eval_trigamma(const Real& x)
{
    return ignore_errors::trigamma(x);
};
template<>
typename gamma_helper<Real>::return_type 
gamma_helper<Real>::eval_gammaln(const Real& x, int& sign)
{
    return ignore_errors::lgamma(x,&sign);
};
template<>
typename gamma_helper<Real>::return_type 
gamma_helper<Real>::eval_polygamma(const Real& x, Integer n)
{
    return ignore_errors::polygamma(n,x);
};
template<>
typename gamma_helper<Real>::return_type 
gamma_helper<Real>::eval_gamma_ratio(const Real& x, const Real& y)
{
    return ignore_errors::tgamma_ratio(x,y);
};
template<>
typename gamma_helper<Real>::return_type 
gamma_helper<Real>::eval_gamma_delta_ratio(const Real& x, const Real& y)
{
    return ignore_errors::tgamma_delta_ratio(x,y);
};
template<>
typename gamma_helper<Real>::return_type 
gamma_helper<Real>::eval_igamma_lower(const Real& x, const Real& y)
{
    return ignore_errors::tgamma_lower(x,y);
};
template<>
typename gamma_helper<Real>::return_type 
gamma_helper<Real>::eval_igamma_upper(const Real& x, const Real& y)
{
    return ignore_errors::tgamma(x,y);
};
template<>
typename gamma_helper<Real>::return_type 
gamma_helper<Real>::eval_igamma_lower_norm(const Real& x, const Real& y)
{
    return ignore_errors::gamma_p(x,y);
};
template<>
typename gamma_helper<Real>::return_type 
gamma_helper<Real>::eval_igamma_upper_norm(const Real& x, const Real& y)
{
    return ignore_errors::gamma_q(x,y);
};
template<>
typename gamma_helper<Real>::return_type 
gamma_helper<Real>::eval_igamma_lower_inv(const Real& x, const Real& y)
{
    return ignore_errors::gamma_p_inv(x,y);
};
template<>
typename gamma_helper<Real>::return_type 
gamma_helper<Real>::eval_igamma_upper_inv(const Real& x, const Real& y)
{
    return ignore_errors::gamma_q_inv(x,y);
};
template<>
typename gamma_helper<Real>::return_type 
gamma_helper<Real>::eval_igamma_lower_inva(const Real& x, const Real& y)
{
    return ignore_errors::gamma_p_inva(x,y);
};
template<>
typename gamma_helper<Real>::return_type 
gamma_helper<Real>::eval_igamma_upper_inva(const Real& x, const Real& y)
{
    return ignore_errors::gamma_q_inva(x,y);
};
template<>
typename gamma_helper<Real>::return_type 
gamma_helper<Real>::eval_igamma_lower_dif(const Real& x, const Real& y)
{
    return ignore_errors::gamma_p_derivative(x,y);
};

//-------------------------------------------------------------------------------
//                                  FLOAT
//-------------------------------------------------------------------------------
template<>
typename gamma_helper<Float>::return_type gamma_helper<Float>::eval_gamma(const Float& x)
{
    if(mrd::scal_func::isnan(x)) 
        return x;
    else if(mrd::scal_func::isinf(x))
        return x > 0.f ? constants::f_inf() : constants::f_nan();

    return ignore_errors::tgamma(x);
};
template<>
typename gamma_helper<Float>::return_type gamma_helper<Float>::eval_gammaln(const Float& x)
{
    if(mrd::scal_func::isnan(x)) 
        return x;
    else if(mrd::scal_func::isinf(x))
        return x > 0.f ? constants::f_inf() : constants::f_nan();

    Integer sign;
    Float val = ignore_errors::lgamma(x,&sign);

    if (sign < 0)
        return constants::f_nan();
    else
        return val;
};

template<>
typename gamma_helper<Float>::return_type 
gamma_helper<Float>::eval_gamma1pm1(const Float& x)
{
    return ignore_errors::tgamma1pm1(x);
};
template<>
typename gamma_helper<Float>::return_type 
gamma_helper<Float>::eval_digamma(const Float& x)
{
    return ignore_errors::digamma(x);
};
template<>
typename gamma_helper<Float>::return_type 
gamma_helper<Float>::eval_trigamma(const Float& x)
{
    return ignore_errors::trigamma(x);
};
template<>
typename gamma_helper<Float>::return_type 
gamma_helper<Float>::eval_gammaln(const Float& x, int& sign)
{
    return ignore_errors::lgamma(x,&sign);
};
template<>
typename gamma_helper<Float>::return_type 
gamma_helper<Float>::eval_polygamma(const Float& x, Integer n)
{
    return ignore_errors::polygamma(n,x);
};
template<>
typename gamma_helper<Float>::return_type 
gamma_helper<Float>::eval_gamma_ratio(const Float& x, const Float& y)
{
    return ignore_errors::tgamma_ratio(x,y);
};
template<>
typename gamma_helper<Float>::return_type 
gamma_helper<Float>::eval_gamma_delta_ratio(const Float& x, const Float& y)
{
    return ignore_errors::tgamma_delta_ratio(x,y);
};
template<>
typename gamma_helper<Float>::return_type 
gamma_helper<Float>::eval_igamma_lower(const Float& x, const Float& y)
{
    return ignore_errors::tgamma_lower(x,y);
};
template<>
typename gamma_helper<Float>::return_type 
gamma_helper<Float>::eval_igamma_upper(const Float& x, const Float& y)
{
    return ignore_errors::tgamma(x,y);
};
template<>
typename gamma_helper<Float>::return_type 
gamma_helper<Float>::eval_igamma_lower_norm(const Float& x, const Float& y)
{
    return ignore_errors::gamma_p(x,y);
};
template<>
typename gamma_helper<Float>::return_type 
gamma_helper<Float>::eval_igamma_upper_norm(const Float& x, const Float& y)
{
    return ignore_errors::gamma_q(x,y);
};
template<>
typename gamma_helper<Float>::return_type 
gamma_helper<Float>::eval_igamma_lower_inv(const Float& x, const Float& y)
{
    return ignore_errors::gamma_p_inv(x,y);
};
template<>
typename gamma_helper<Float>::return_type 
gamma_helper<Float>::eval_igamma_upper_inv(const Float& x, const Float& y)
{
    return ignore_errors::gamma_q_inv(x,y);
};
template<>
typename gamma_helper<Float>::return_type 
gamma_helper<Float>::eval_igamma_lower_inva(const Float& x, const Float& y)
{
    return ignore_errors::gamma_p_inva(x,y);
};
template<>
typename gamma_helper<Float>::return_type 
gamma_helper<Float>::eval_igamma_upper_inva(const Float& x, const Float& y)
{
    return ignore_errors::gamma_q_inva(x,y);
};
template<>
typename gamma_helper<Float>::return_type 
gamma_helper<Float>::eval_igamma_lower_dif(const Float& x, const Float& y)
{
    return ignore_errors::gamma_p_derivative(x,y);
};

//-------------------------------------------------------------------------------
//                                  Complex
//-------------------------------------------------------------------------------
template<>
typename gamma_helper<Complex>::return_type 
gamma_helper<Complex>::eval_gamma(const Complex& x)
{
    return scal_func::impl_gamma(x);
};
template<>
typename gamma_helper<Complex>::return_type 
gamma_helper<Complex>::eval_gammaln(const Complex& x)
{
    return scal_func::impl_gammaln(x);
};
template<>
typename gamma_helper<Complex>::return_type 
gamma_helper<Complex>::eval_gamma1pm1(const Complex& x)
{
    if (imag(x) == 0.0)
        return gamma_helper<Real>::eval_gamma1pm1(real(x));
    else
        throw matcl::error::function_not_defined_for_complex("gamma1pm1"); 
};
template<>
typename gamma_helper<Complex>::return_type 
gamma_helper<Complex>::eval_digamma(const Complex& x)
{
    if (imag(x) == 0.0)
        return gamma_helper<Real>::eval_digamma(real(x));
    else
        throw matcl::error::function_not_defined_for_complex("digamma"); 
};
template<>
typename gamma_helper<Complex>::return_type 
gamma_helper<Complex>::eval_trigamma(const Complex& x)
{
    if (imag(x) == 0.0)
        return gamma_helper<Real>::eval_trigamma(real(x));
    else
        throw matcl::error::function_not_defined_for_complex("trigamma"); 
};
template<>
typename gamma_helper<Complex>::return_type 
gamma_helper<Complex>::eval_gammaln(const Complex& x, int& sign)
{
    if (imag(x) == 0.0)
        return gamma_helper<Real>::eval_gammaln(real(x),sign);
    else
        throw matcl::error::function_not_defined_for_complex("gammaln"); 
};
template<>
typename gamma_helper<Complex>::return_type 
gamma_helper<Complex>::eval_polygamma(const Complex& x, Integer n)
{
    if (imag(x) == 0.0)
        return gamma_helper<Real>::eval_polygamma(real(x),n);
    else
        throw matcl::error::function_not_defined_for_complex("polygamma"); 
};
template<>
typename gamma_helper<Complex>::return_type 
gamma_helper<Complex>::eval_gamma_ratio(const Complex& x, const Complex& y)
{
    if (imag(x) == 0.0 && imag(y) == 0.0)
        return gamma_helper<Real>::eval_gamma_ratio(real(x),real(y));
    else
        throw matcl::error::function_not_defined_for_complex("gamma_ratio"); 
};
template<>
typename gamma_helper<Complex>::return_type 
gamma_helper<Complex>::eval_gamma_delta_ratio(const Complex& x, const Complex& y)
{
    if (imag(x) == 0.0 && imag(y) == 0.0)
        return gamma_helper<Real>::eval_gamma_delta_ratio(real(x),real(y));
    else
        throw matcl::error::function_not_defined_for_complex("gamma_delta_ratio"); 
};
template<>
typename gamma_helper<Complex>::return_type 
gamma_helper<Complex>::eval_igamma_lower(const Complex& x, const Complex& y)
{
    if (imag(x) == 0.0 && imag(y) == 0.0)
        return gamma_helper<Real>::eval_igamma_lower(real(x),real(y));
    else
        throw matcl::error::function_not_defined_for_complex("igamma_lower"); 
};
template<>
typename gamma_helper<Complex>::return_type 
gamma_helper<Complex>::eval_igamma_upper(const Complex& x, const Complex& y)
{
    if (imag(x) == 0.0 && imag(y) == 0.0)
        return gamma_helper<Real>::eval_igamma_upper(real(x),real(y));
    else
        throw matcl::error::function_not_defined_for_complex("igamma_upper"); 
};
template<>
typename gamma_helper<Complex>::return_type 
gamma_helper<Complex>::eval_igamma_lower_norm(const Complex& x, const Complex& y)
{
    if (imag(x) == 0.0 && imag(y) == 0.0)
        return gamma_helper<Real>::eval_igamma_lower_norm(real(x),real(y));
    else
        throw matcl::error::function_not_defined_for_complex("igamma_lower_norm"); 
};
template<>
typename gamma_helper<Complex>::return_type 
gamma_helper<Complex>::eval_igamma_upper_norm(const Complex& x, const Complex& y)
{
    if (imag(x) == 0.0 && imag(y) == 0.0)
        return gamma_helper<Real>::eval_igamma_upper_norm(real(x),real(y));
    else
        throw matcl::error::function_not_defined_for_complex("igamma_upper_norm"); 
};
template<>
typename gamma_helper<Complex>::return_type 
gamma_helper<Complex>::eval_igamma_lower_inv(const Complex& x, const Complex& y)
{
    if (imag(x) == 0.0 && imag(y) == 0.0)
        return gamma_helper<Real>::eval_igamma_lower_inv(real(x),real(y));
    else
        throw matcl::error::function_not_defined_for_complex("igamma_lower_inv"); 
};
template<>
typename gamma_helper<Complex>::return_type 
gamma_helper<Complex>::eval_igamma_upper_inv(const Complex& x, const Complex& y)
{
    if (imag(x) == 0.0 && imag(y) == 0.0)
        return gamma_helper<Real>::eval_igamma_upper_inv(real(x),real(y));
    else
        throw matcl::error::function_not_defined_for_complex("igamma_upper_inv"); 
};
template<>
typename gamma_helper<Complex>::return_type 
gamma_helper<Complex>::eval_igamma_lower_inva(const Complex& x, const Complex& y)
{
    if (imag(x) == 0.0 && imag(y) == 0.0)
        return gamma_helper<Real>::eval_igamma_lower_inva(real(x),real(y));
    else
        throw matcl::error::function_not_defined_for_complex("igamma_lower_inva"); 
};
template<>
typename gamma_helper<Complex>::return_type 
gamma_helper<Complex>::eval_igamma_upper_inva(const Complex& x, const Complex& y)
{
    if (imag(x) == 0.0 && imag(y) == 0.0)
        return gamma_helper<Real>::eval_igamma_upper_inva(real(x),real(y));
    else
        throw matcl::error::function_not_defined_for_complex("igamma_upper_inva"); 
};
template<>
typename gamma_helper<Complex>::return_type 
gamma_helper<Complex>::eval_igamma_lower_dif(const Complex& x, const Complex& y)
{
    if (imag(x) == 0.0 && imag(y) == 0.0)
        return gamma_helper<Real>::eval_igamma_lower_dif(real(x),real(y));
    else
        throw matcl::error::function_not_defined_for_complex("igamma_lower_dif"); 
};

//-------------------------------------------------------------------------------
//                                  Float_complex
//-------------------------------------------------------------------------------
template<>
typename gamma_helper<Float_complex>::return_type 
gamma_helper<Float_complex>::eval_gamma(const Float_complex& x)
{
    return scal_func::impl_gamma(x);
};
template<>
typename gamma_helper<Float_complex>::return_type 
gamma_helper<Float_complex>::eval_gammaln(const Float_complex& x)
{
    return scal_func::impl_gammaln(x);
};
template<>
typename gamma_helper<Float_complex>::return_type 
gamma_helper<Float_complex>::eval_gamma1pm1(const Float_complex& x)
{
    if (imag(x) == 0.0)
        return gamma_helper<Float>::eval_gamma1pm1(real(x));
    else
        throw matcl::error::function_not_defined_for_complex("gamma1pm1"); 
};
template<>
typename gamma_helper<Float_complex>::return_type 
gamma_helper<Float_complex>::eval_digamma(const Float_complex& x)
{
    if (imag(x) == 0.0)
        return gamma_helper<Float>::eval_digamma(real(x));
    else
        throw matcl::error::function_not_defined_for_complex("digamma"); 
};
template<>
typename gamma_helper<Float_complex>::return_type 
gamma_helper<Float_complex>::eval_trigamma(const Float_complex& x)
{
    if (imag(x) == 0.0)
        return gamma_helper<Float>::eval_trigamma(real(x));
    else
        throw matcl::error::function_not_defined_for_complex("trigamma"); 
};
template<>
typename gamma_helper<Float_complex>::return_type 
gamma_helper<Float_complex>::eval_gammaln(const Float_complex& x, int& sign)
{
    if (imag(x) == 0.0)
        return gamma_helper<Float>::eval_gammaln(real(x),sign);
    else
        throw matcl::error::function_not_defined_for_complex("gammaln"); 
};
template<>
typename gamma_helper<Float_complex>::return_type 
gamma_helper<Float_complex>::eval_polygamma(const Float_complex& x, Integer n)
{
    if (imag(x) == 0.0)
        return gamma_helper<Float>::eval_polygamma(real(x),n);
    else
        throw matcl::error::function_not_defined_for_complex("polygamma"); 
};
template<>
typename gamma_helper<Float_complex>::return_type 
gamma_helper<Float_complex>::eval_gamma_ratio(const Float_complex& x, const Float_complex& y)
{
    if (imag(x) == 0.0 && imag(y) == 0.0)
        return gamma_helper<Float>::eval_gamma_ratio(real(x),real(y));
    else
        throw matcl::error::function_not_defined_for_complex("gamma_ratio"); 
};
template<>
typename gamma_helper<Float_complex>::return_type 
gamma_helper<Float_complex>::eval_gamma_delta_ratio(const Float_complex& x, const Float_complex& y)
{
    if (imag(x) == 0.0 && imag(y) == 0.0)
        return gamma_helper<Float>::eval_gamma_delta_ratio(real(x),real(y));
    else
        throw matcl::error::function_not_defined_for_complex("gamma_delta_ratio"); 
};
template<>
typename gamma_helper<Float_complex>::return_type 
gamma_helper<Float_complex>::eval_igamma_lower(const Float_complex& x, const Float_complex& y)
{
    if (imag(x) == 0.0 && imag(y) == 0.0)
        return gamma_helper<Float>::eval_igamma_lower(real(x),real(y));
    else
        throw matcl::error::function_not_defined_for_complex("igamma_lower"); 
};
template<>
typename gamma_helper<Float_complex>::return_type 
gamma_helper<Float_complex>::eval_igamma_upper(const Float_complex& x, const Float_complex& y)
{
    if (imag(x) == 0.0 && imag(y) == 0.0)
        return gamma_helper<Float>::eval_igamma_upper(real(x),real(y));
    else
        throw matcl::error::function_not_defined_for_complex("igamma_upper"); 
};
template<>
typename gamma_helper<Float_complex>::return_type 
gamma_helper<Float_complex>::eval_igamma_lower_norm(const Float_complex& x, const Float_complex& y)
{
    if (imag(x) == 0.0 && imag(y) == 0.0)
        return gamma_helper<Float>::eval_igamma_lower_norm(real(x),real(y));
    else
        throw matcl::error::function_not_defined_for_complex("igamma_lower_norm"); 
};
template<>
typename gamma_helper<Float_complex>::return_type 
gamma_helper<Float_complex>::eval_igamma_upper_norm(const Float_complex& x, const Float_complex& y)
{
    if (imag(x) == 0.0 && imag(y) == 0.0)
        return gamma_helper<Float>::eval_igamma_upper_norm(real(x),real(y));
    else
        throw matcl::error::function_not_defined_for_complex("igamma_upper_norm"); 
};
template<>
typename gamma_helper<Float_complex>::return_type 
gamma_helper<Float_complex>::eval_igamma_lower_inv(const Float_complex& x, const Float_complex& y)
{
    if (imag(x) == 0.0 && imag(y) == 0.0)
        return gamma_helper<Float>::eval_igamma_lower_inv(real(x),real(y));
    else
        throw matcl::error::function_not_defined_for_complex("igamma_lower_inv"); 
};
template<>
typename gamma_helper<Float_complex>::return_type 
gamma_helper<Float_complex>::eval_igamma_upper_inv(const Float_complex& x, const Float_complex& y)
{
    if (imag(x) == 0.0 && imag(y) == 0.0)
        return gamma_helper<Float>::eval_igamma_upper_inv(real(x),real(y));
    else
        throw matcl::error::function_not_defined_for_complex("igamma_upper_inv"); 
};
template<>
typename gamma_helper<Float_complex>::return_type 
gamma_helper<Float_complex>::eval_igamma_lower_inva(const Float_complex& x, const Float_complex& y)
{
    if (imag(x) == 0.0 && imag(y) == 0.0)
        return gamma_helper<Float>::eval_igamma_lower_inva(real(x),real(y));
    else
        throw matcl::error::function_not_defined_for_complex("igamma_lower_inva"); 
};
template<>
typename gamma_helper<Float_complex>::return_type 
gamma_helper<Float_complex>::eval_igamma_upper_inva(const Float_complex& x, const Float_complex& y)
{
    if (imag(x) == 0.0 && imag(y) == 0.0)
        return gamma_helper<Float>::eval_igamma_upper_inva(real(x),real(y));
    else
        throw matcl::error::function_not_defined_for_complex("igamma_upper_inva"); 
};
template<>
typename gamma_helper<Float_complex>::return_type 
gamma_helper<Float_complex>::eval_igamma_lower_dif(const Float_complex& x, const Float_complex& y)
{
    if (imag(x) == 0.0 && imag(y) == 0.0)
        return gamma_helper<Float>::eval_igamma_lower_dif(real(x),real(y));
    else
        throw matcl::error::function_not_defined_for_complex("igamma_lower_dif"); 
};

//-------------------------------------------------------------------------------
//                                  Integer
//-------------------------------------------------------------------------------
template<>
typename gamma_helper<Integer>::return_type gamma_helper<Integer>::eval_gamma(const Integer& x)
{
    return gamma_helper<Real>::eval_gamma((Real)x);
};
template<>
typename gamma_helper<Integer>::return_type gamma_helper<Integer>::eval_gammaln(const Integer& x)
{
    return gamma_helper<Real>::eval_gammaln((Real)x);
};
template<>
typename gamma_helper<Integer>::return_type 
gamma_helper<Integer>::eval_gamma1pm1(const Integer& x)
{
    return gamma_helper<Real>::eval_gamma1pm1((Real)x);
};
template<>
typename gamma_helper<Integer>::return_type 
gamma_helper<Integer>::eval_digamma(const Integer& x)
{
    return gamma_helper<Real>::eval_digamma((Real)x);
};
template<>
typename gamma_helper<Integer>::return_type 
gamma_helper<Integer>::eval_trigamma(const Integer& x)
{
    return gamma_helper<Real>::eval_trigamma((Real)x);
};
template<>
typename gamma_helper<Integer>::return_type 
gamma_helper<Integer>::eval_gammaln(const Integer& x, int& sign)
{
    return gamma_helper<Real>::eval_gammaln((Real)x,sign);
};
template<>
typename gamma_helper<Integer>::return_type 
gamma_helper<Integer>::eval_polygamma(const Integer& x, Integer n)
{
    return gamma_helper<Real>::eval_polygamma((Real)x,n);
};
template<>
typename gamma_helper<Integer>::return_type 
gamma_helper<Integer>::eval_gamma_ratio(const Integer& x, const Integer& y)
{
    return gamma_helper<Real>::eval_gamma_ratio((Real)x,(Real)y);
};
template<>
typename gamma_helper<Integer>::return_type 
gamma_helper<Integer>::eval_gamma_delta_ratio(const Integer& x, const Integer& y)
{
    return gamma_helper<Real>::eval_gamma_delta_ratio((Real)x,(Real)y);
};
template<>
typename gamma_helper<Integer>::return_type 
gamma_helper<Integer>::eval_igamma_lower(const Integer& x, const Integer& y)
{
    return gamma_helper<Real>::eval_igamma_lower((Real)x,(Real)y);
};
template<>
typename gamma_helper<Integer>::return_type 
gamma_helper<Integer>::eval_igamma_upper(const Integer& x, const Integer& y)
{
    return gamma_helper<Real>::eval_igamma_upper((Real)x,(Real)y);
};
template<>
typename gamma_helper<Integer>::return_type 
gamma_helper<Integer>::eval_igamma_lower_norm(const Integer& x, const Integer& y)
{
    return gamma_helper<Real>::eval_igamma_lower_norm((Real)x,(Real)y);
};
template<>
typename gamma_helper<Integer>::return_type 
gamma_helper<Integer>::eval_igamma_upper_norm(const Integer& x, const Integer& y)
{
    return gamma_helper<Real>::eval_igamma_upper_norm((Real)x,(Real)y);
};
template<>
typename gamma_helper<Integer>::return_type 
gamma_helper<Integer>::eval_igamma_lower_inv(const Integer& x, const Integer& y)
{
    return gamma_helper<Real>::eval_igamma_lower_inv((Real)x,(Real)y);
};
template<>
typename gamma_helper<Integer>::return_type 
gamma_helper<Integer>::eval_igamma_upper_inv(const Integer& x, const Integer& y)
{
    return gamma_helper<Real>::eval_igamma_upper_inv((Real)x,(Real)y);
};
template<>
typename gamma_helper<Integer>::return_type 
gamma_helper<Integer>::eval_igamma_lower_inva(const Integer& x, const Integer& y)
{
    return gamma_helper<Real>::eval_igamma_lower_inva((Real)x,(Real)y);
};
template<>
typename gamma_helper<Integer>::return_type 
gamma_helper<Integer>::eval_igamma_upper_inva(const Integer& x, const Integer& y)
{
    return gamma_helper<Real>::eval_igamma_upper_inva((Real)x,(Real)y);
};
template<>
typename gamma_helper<Integer>::return_type 
gamma_helper<Integer>::eval_igamma_lower_dif(const Integer& x, const Integer& y)
{
    return gamma_helper<Real>::eval_igamma_lower_dif((Real)x,(Real)y);
};

//-------------------------------------------------------------------------------
//                                  Object
//-------------------------------------------------------------------------------
template<>
typename gamma_helper<Object>::return_type gamma_helper<Object>::eval_gamma(const Object& x)
{
    return object_func::gamma(x);
};
template<>
typename gamma_helper<Object>::return_type gamma_helper<Object>::eval_gammaln(const Object& x)
{
    return object_func::gammaln(x);
};
template<>
typename gamma_helper<Object>::return_type 
gamma_helper<Object>::eval_gamma1pm1(const Object& x)
{
    return object_func::gamma1pm1(x);
};
template<>
typename gamma_helper<Object>::return_type 
gamma_helper<Object>::eval_digamma(const Object& x)
{
    return object_func::digamma(x);
};
template<>
typename gamma_helper<Object>::return_type 
gamma_helper<Object>::eval_trigamma(const Object& x)
{
    return object_func::trigamma(x);
};

template<>
typename gamma_helper<Object>::return_type 
gamma_helper<Object>::eval_gammaln(const Object& x, int& sign)
{
    (void)x;
    (void)sign;
    //TODO: impl for object
    throw matcl::error::object_value_type_not_allowed("gammaln");
};

template<>
typename gamma_helper<Object>::return_type 
gamma_helper<Object>::eval_polygamma(const Object& x, Integer n)
{
    return object_func::polygamma(x,n);
};
template<>
typename gamma_helper<Object>::return_type 
gamma_helper<Object>::eval_gamma_ratio(const Object& x, const Object& y)
{
    return object_func::gamma_ratio(x,y);
};
template<>
typename gamma_helper<Object>::return_type 
gamma_helper<Object>::eval_gamma_delta_ratio(const Object& x, const Object& y)
{
    return object_func::gamma_delta_ratio(x,y);
};
template<>
typename gamma_helper<Object>::return_type 
gamma_helper<Object>::eval_igamma_lower(const Object& x, const Object& y)
{
    return object_func::igamma_lower(x,y);
};
template<>
typename gamma_helper<Object>::return_type 
gamma_helper<Object>::eval_igamma_upper(const Object& x, const Object& y)
{
    return object_func::igamma_upper(x,y);
};
template<>
typename gamma_helper<Object>::return_type 
gamma_helper<Object>::eval_igamma_lower_norm(const Object& x, const Object& y)
{
    return object_func::igamma_lower_norm(x,y);
};
template<>
typename gamma_helper<Object>::return_type 
gamma_helper<Object>::eval_igamma_upper_norm(const Object& x, const Object& y)
{
    return object_func::igamma_upper_norm(x,y);
};
template<>
typename gamma_helper<Object>::return_type 
gamma_helper<Object>::eval_igamma_lower_inv(const Object& x, const Object& y)
{
    return object_func::igamma_lower_inv(x,y);
};
template<>
typename gamma_helper<Object>::return_type 
gamma_helper<Object>::eval_igamma_upper_inv(const Object& x, const Object& y)
{
    return object_func::igamma_upper_inv(x,y);
};
template<>
typename gamma_helper<Object>::return_type 
gamma_helper<Object>::eval_igamma_lower_inva(const Object& x, const Object& y)
{
    return object_func::igamma_lower_inva(x,y);
};
template<>
typename gamma_helper<Object>::return_type 
gamma_helper<Object>::eval_igamma_upper_inva(const Object& x, const Object& y)
{
    return object_func::igamma_upper_inva(x,y);
};
template<>
typename gamma_helper<Object>::return_type 
gamma_helper<Object>::eval_igamma_lower_dif(const Object& x, const Object& y)
{
    return object_func::igamma_lower_dif(x,y);
};

template struct gamma_helper<Integer>;
template struct gamma_helper<Real>;
template struct gamma_helper<Float>;
template struct gamma_helper<Complex>;
template struct gamma_helper<Float_complex>;
template struct gamma_helper<Object>;

//-------------------------------------------------------------------------------
//                                  scalar functors
//-------------------------------------------------------------------------------
struct func_gamma
{
    template<class T>
    auto eval(const T& v) const -> decltype(md::gamma_helper<T>::eval_gamma(v))
    {
        return md::gamma_helper<T>::eval_gamma(v);
    }
};

struct func_gammaln
{
    template<class T>
    auto eval(const T& v) const -> decltype(md::gamma_helper<T>::eval_gammaln(v))
    {
        return md::gamma_helper<T>::eval_gammaln(v);
    }
};

struct func_digamma
{
    template<class T>
    auto eval(const T& v) const -> decltype(md::gamma_helper<T>::eval_digamma(v))
    {
        return md::gamma_helper<T>::eval_digamma(v);
    }
};

struct func_trigamma
{
    template<class T>
    auto eval(const T& v) const -> decltype(md::gamma_helper<T>::eval_trigamma(v))
    {
        return md::gamma_helper<T>::eval_trigamma(v);
    }
};

struct func_polygamma
{
    Integer m_order;

    func_polygamma(Integer n)   : m_order(n){};

    template<class T>
    auto eval(const T& v) const -> decltype(md::gamma_helper<T>::eval_polygamma(v,m_order))
    {
        return md::gamma_helper<T>::eval_polygamma(v,m_order);
    }
};

};};

namespace matcl
{

Matrix matcl::gamma(const Matrix &A)
{
    return eval_scalar_func(A, details::func_gamma());
};
Matrix matcl::gammaln(const Matrix &A)
{
    return eval_scalar_func(A, details::func_gammaln());
};
Matrix matcl::digamma(const Matrix &A)
{
    return eval_scalar_func(A, details::func_digamma());
};
Matrix matcl::trigamma(const Matrix &A)
{
    return eval_scalar_func(A, details::func_trigamma());
};

Matrix matcl::polygamma(const Matrix &A, Integer n)
{
    return eval_scalar_func(A, details::func_polygamma(n));
};

};
