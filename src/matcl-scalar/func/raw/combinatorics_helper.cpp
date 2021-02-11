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

#include "matcl-scalar/details/combinatorics.inl"

#pragma warning(push)
#pragma warning(disable:4127)   //conditional expression is constant
#include <boost/math/special_functions.hpp>
#include <boost/math/special_functions/prime.hpp>
#include <boost/math/special_functions/bernoulli.hpp>
#include <boost/math/special_functions/factorials.hpp>
#include <boost/math/special_functions/binomial.hpp>
#pragma warning(pop)

namespace matcl { namespace details
{

namespace md = matcl::details;

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

size_t primes_helper::eval_prime(Integer n)
{
    return boost::math::prime(n,boost__ignore_policy());
};
Integer primes_helper::eval_prime_max_count()
{
    return boost::math::max_prime;
};

unsigned error_unsigned(Integer, const std::string& func_name)
{
    throw error::integer_value_type_not_allowed(func_name);
};

//-------------------------------------------------------------------------------
//                                  INTEGER
//-------------------------------------------------------------------------------
template<>
Real combinatorics_helper<Integer>::eval_bernoulli_b2n(Integer n)
{
    return boost::math::bernoulli_b2n<Real>(n,boost__ignore_policy());
}

template<>
Integer combinatorics_helper<Integer>::max_bernoulli_b2n()
{
    return boost::math::max_bernoulli_b2n<Real>::value;
}

template<>
Real combinatorics_helper<Integer>::eval_factorial(Integer i)
{
    if (i < 0)
        return constants::nan();

    return boost::math::factorial<Real>(i,boost__ignore_policy());
}

template<>
Real combinatorics_helper<Integer>::eval_double_factorial(Integer i)
{
    if (i < 0)
        return constants::nan();

    return boost::math::double_factorial<Real>(i,boost__ignore_policy());
}

template<>
Real combinatorics_helper<Integer>::eval_binomial_coefficient(Integer n, Integer k)
{
    if (n < 0 || k < 0 || k > n)
        return constants::nan();

    return boost::math::binomial_coefficient<Real>(n,k,boost__ignore_policy());
}

template<>
Real combinatorics_helper<Integer>::eval_rising_factorial(const Integer& x, Integer i)
{
    return boost::math::rising_factorial<Real>(x,i,boost__ignore_policy());
}

template<>
Real combinatorics_helper<Integer>::eval_falling_factorial(const Integer& x, Integer i)
{
    if (i < 0)
        return constants::nan();

    return boost::math::falling_factorial<Real>(x,i,boost__ignore_policy());
}

//-------------------------------------------------------------------------------
//                                  REAL
//-------------------------------------------------------------------------------
template<>
Real combinatorics_helper<Real>::eval_bernoulli_b2n(Integer n)
{
    return boost::math::bernoulli_b2n<Real>(n,boost__ignore_policy());
}

template<>
Integer combinatorics_helper<Real>::max_bernoulli_b2n()
{
    return boost::math::max_bernoulli_b2n<Real>::value;
}

template<>
Real combinatorics_helper<Real>::eval_factorial(Integer i)
{
    if (i < 0)
        return constants::nan();

    return boost::math::factorial<Real>(i,boost__ignore_policy());
}

template<>
Real combinatorics_helper<Real>::eval_double_factorial(Integer i)
{
    if (i < 0)
        return constants::nan();

    return boost::math::double_factorial<Real>(i,boost__ignore_policy());
}

template<>
Real combinatorics_helper<Real>::eval_binomial_coefficient(Integer n, Integer k)
{
    if (n < 0 || k < 0 || k > n)
        return constants::nan();

    return boost::math::binomial_coefficient<Real>(n,k,boost__ignore_policy());
}

template<>
Real combinatorics_helper<Real>::eval_rising_factorial(const Real& x, Integer i)
{
    return boost::math::rising_factorial(x,i,boost__ignore_policy());
}

template<>
Real combinatorics_helper<Real>::eval_falling_factorial(const Real& x, Integer i)
{
    if (i < 0)
        return constants::nan();

    return boost::math::falling_factorial(x,i,boost__ignore_policy());
}

//-------------------------------------------------------------------------------
//                                  FLOAT
//-------------------------------------------------------------------------------
template<>
Float combinatorics_helper<Float>::eval_bernoulli_b2n(Integer n)
{
    return boost::math::bernoulli_b2n<Float>(n,boost__ignore_policy());
}

template<>
Integer combinatorics_helper<Float>::max_bernoulli_b2n()
{
    return boost::math::max_bernoulli_b2n<Float>::value;
}

template<>
Float combinatorics_helper<Float>::eval_factorial(Integer i)
{
    if (i < 0)
        return constants::f_nan();

    return boost::math::factorial<Float>(i,boost__ignore_policy());
}

template<>
Float combinatorics_helper<Float>::eval_double_factorial(Integer i)
{
    if (i < 0)
        return constants::f_nan();

    return boost::math::double_factorial<Float>(i,boost__ignore_policy());
}

template<>
Float combinatorics_helper<Float>::eval_binomial_coefficient(Integer n, Integer k)
{
    if (n < 0 || k < 0 || k > n)
        return constants::f_nan();

    return boost::math::binomial_coefficient<Float>(n,k,boost__ignore_policy());
}

template<>
Float combinatorics_helper<Float>::eval_rising_factorial(const Float& x, Integer i)
{
    return boost::math::rising_factorial(x,i,boost__ignore_policy());
}

template<>
Float combinatorics_helper<Float>::eval_falling_factorial(const Float& x, Integer i)
{
    if (i < 0)
        return constants::f_nan();

    return boost::math::falling_factorial(x,i,boost__ignore_policy());
}

//-------------------------------------------------------------------------------
//                                  COMPLEX
//-------------------------------------------------------------------------------
template<>
Complex combinatorics_helper<Complex>::eval_bernoulli_b2n(Integer n)
{
    return boost::math::bernoulli_b2n<Real>(n,boost__ignore_policy());
}

template<>
Integer combinatorics_helper<Complex>::max_bernoulli_b2n()
{
    return boost::math::max_bernoulli_b2n<Real>::value;
}

template<>
Complex combinatorics_helper<Complex>::eval_factorial(Integer i)
{
    return boost::math::factorial<Real>(i,boost__ignore_policy());
}

template<>
Complex combinatorics_helper<Complex>::eval_double_factorial(Integer i)
{
    return boost::math::double_factorial<Real>(i,boost__ignore_policy());
}

template<>
Complex combinatorics_helper<Complex>::eval_binomial_coefficient(Integer n, Integer k)
{
    return boost::math::binomial_coefficient<Real>(n,k,boost__ignore_policy());
}

template<>
Complex combinatorics_helper<Complex>::eval_rising_factorial(const Complex& x, Integer i)
{
    if (imag(x) == 0)
        return combinatorics_helper<Real>::eval_rising_factorial(real(x),i);
    else
        throw matcl::error::function_not_defined_for_complex("rising_factorial"); 
}

template<>
Complex combinatorics_helper<Complex>::eval_falling_factorial(const Complex& x, Integer i)
{
    if (imag(x) == 0)
        return combinatorics_helper<Real>::eval_falling_factorial(real(x),i);
    else
        throw matcl::error::function_not_defined_for_complex("falling_factorial"); 
}

//-------------------------------------------------------------------------------
//                                  FLOAT_COMPLEX
//-------------------------------------------------------------------------------
template<>
Float_complex combinatorics_helper<Float_complex>::eval_bernoulli_b2n(Integer n)
{
    return boost::math::bernoulli_b2n<Float>(n,boost__ignore_policy());
}

template<>
Integer combinatorics_helper<Float_complex>::max_bernoulli_b2n()
{
    return boost::math::max_bernoulli_b2n<Float>::value;
}

template<>
Float_complex combinatorics_helper<Float_complex>::eval_factorial(Integer i)
{
    return boost::math::factorial<Float>(i,boost__ignore_policy());
}

template<>
Float_complex combinatorics_helper<Float_complex>::eval_double_factorial(Integer i)
{
    return boost::math::double_factorial<Float>(i,boost__ignore_policy());
}

template<>
Float_complex combinatorics_helper<Float_complex>::eval_binomial_coefficient(Integer n, Integer k)
{
    return boost::math::binomial_coefficient<Float>(n,k,boost__ignore_policy());
}

template<>
Float_complex combinatorics_helper<Float_complex>::eval_rising_factorial(const Float_complex& x, Integer i)
{
    if (imag(x) == 0)
        return combinatorics_helper<Float>::eval_rising_factorial(real(x),i);
    else
        throw matcl::error::function_not_defined_for_complex("rising_factorial"); 
}

template<>
Float_complex combinatorics_helper<Float_complex>::eval_falling_factorial(const Float_complex& x, Integer i)
{
    if (imag(x) == 0)
        return combinatorics_helper<Float>::eval_falling_factorial(real(x),i);
    else
        throw matcl::error::function_not_defined_for_complex("falling_factorial"); 
}

//-------------------------------------------------------------------------------
//                                  OBJECT
//-------------------------------------------------------------------------------
template<>
dynamic::object combinatorics_helper<dynamic::object>::eval_bernoulli_b2n(Integer)
{
    throw error::object_value_type_not_allowed("bernoulli_b2n");
}

template<>
Integer combinatorics_helper<dynamic::object>::max_bernoulli_b2n()
{
    throw error::object_value_type_not_allowed("max_bernoulli_b2n");
}

template<>
dynamic::object combinatorics_helper<dynamic::object>::eval_factorial(Integer)
{
    throw error::object_value_type_not_allowed("factorial");
}

template<>
dynamic::object combinatorics_helper<dynamic::object>::eval_double_factorial(Integer)
{
    throw error::object_value_type_not_allowed("double_factorial");
}

template<>
dynamic::object combinatorics_helper<dynamic::object>::eval_binomial_coefficient(Integer, Integer)
{
    throw error::object_value_type_not_allowed("binomial_coefficient");
}

template<>
dynamic::object 
combinatorics_helper<dynamic::object>::eval_rising_factorial(const dynamic::object& x, Integer i)
{
    return dynamic::rising_factorial(x,i);
}

template<>
dynamic::object 
combinatorics_helper<dynamic::object>::eval_falling_factorial(const dynamic::object& x, Integer i)
{
    return dynamic::falling_factorial(x,i);
}

template struct combinatorics_helper<Integer>;
template struct combinatorics_helper<Real>;
template struct combinatorics_helper<Float>;
template struct combinatorics_helper<Complex>;
template struct combinatorics_helper<Float_complex>;
template struct combinatorics_helper<dynamic::object>;

void test_combinatorics_inst()
{
    Integer n   = 0;
    Integer i   = 0;
    Integer k   = 0;
    Real x      = 0;
    Integer xi  = 0;

    bernoulli_b2n<Integer>(n);
    bernoulli_b2n<Real>(n);
    bernoulli_b2n(n);
    fbernoulli_b2n(n);

    max_bernoulli_b2n<Real>();
    max_bernoulli_b2n<Integer>();
    max_bernoulli_b2n();
    fmax_bernoulli_b2n();

    prime(n);
    prime_max_count();

    factorial<Real>(i);
    factorial<Integer>(i);
    factorial(i);
    ffactorial(i);

    double_factorial<Real>(i);
    double_factorial<Integer>(i);
    double_factorial(i);
    fdouble_factorial(i);

    rising_factorial(x, i);
    rising_factorial(xi, i);

    falling_factorial(x, i);
    falling_factorial(xi, i);
                        
    binomial_coefficient<Real>(n, k);
    binomial_coefficient<Integer>(n, k);
    binomial_coefficient(n, k);
    fbinomial_coefficient(n, k);
};

};};
