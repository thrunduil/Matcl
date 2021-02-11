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

#pragma once

#include "matcl-scalar/lib_functions/func_unary.h"

namespace matcl { namespace details
{

template<class S>
struct MATCL_SCALAR_EXPORT combinatorics_helper
{
    using SF = typename unify_types<S,Float>::type;

    static SF       eval_bernoulli_b2n(Integer n);
    static Integer  max_bernoulli_b2n();
    static SF       eval_factorial(Integer i);
    static SF       eval_double_factorial(Integer i);
    static SF       eval_binomial_coefficient(Integer n, Integer k);

    static SF       eval_rising_factorial(const S& x, Integer i);
    static SF       eval_falling_factorial(const S& x, Integer i);
};

struct MATCL_SCALAR_EXPORT primes_helper
{
    static size_t   eval_prime(Integer n);
    static Integer  eval_prime_max_count();
};

};};

namespace matcl
{

template<class S1, class Enable>
typename md::unify_types_promote<S1,Float>::type
matcl::bernoulli_b2n(Integer n)
{
    using S = typename md::unify_types_promote<S1,Float>::type;
    return details::combinatorics_helper<S>::eval_bernoulli_b2n(n);
};

inline Real matcl::bernoulli_b2n(Integer n)
{
    return details::combinatorics_helper<Real>::eval_bernoulli_b2n(n);
};
inline Float matcl::fbernoulli_b2n(Integer n)
{
    return details::combinatorics_helper<Float>::eval_bernoulli_b2n(n);
};

template<class S1, class Enable>
Integer matcl::max_bernoulli_b2n()
{
    using S = typename md::promote_scalar<S1>::type;
    return details::combinatorics_helper<S>::max_bernoulli_b2n();
};

inline Integer matcl::max_bernoulli_b2n()
{
    return details::combinatorics_helper<Real>::max_bernoulli_b2n();
};

inline Integer matcl::fmax_bernoulli_b2n()
{
    return details::combinatorics_helper<Float>::max_bernoulli_b2n();
};

inline size_t matcl::prime(Integer n)
{
    return details::primes_helper::eval_prime(n);
};

inline Integer matcl::prime_max_count()
{
    return details::primes_helper::eval_prime_max_count();
};

template<class S1, class Enable>
typename md::unify_types_promote<S1,Float>::type
matcl::factorial(Integer i)
{
    using S = typename md::unify_types_promote<S1,Float>::type;
    return details::combinatorics_helper<S>::eval_factorial(i);
};

inline Real matcl::factorial(Integer i)
{
    return details::combinatorics_helper<Real>::eval_factorial(i);
};
inline Float matcl::ffactorial(Integer i)
{
    return details::combinatorics_helper<Float>::eval_factorial(i);
};

template<class S1, class Enable>
typename md::unify_types_promote<S1,Float>::type
matcl::double_factorial(Integer i)
{
    using S = typename md::unify_types_promote<S1,Float>::type;
    return details::combinatorics_helper<S>::eval_double_factorial(i);
};

inline Real matcl::double_factorial(Integer i)
{
    return details::combinatorics_helper<Real>::eval_double_factorial(i);
};
inline Float matcl::fdouble_factorial(Integer i)
{
    return details::combinatorics_helper<Float>::eval_double_factorial(i);
};

template<class S1, class Enable>
typename md::unify_types_promote<S1,Float>::type
matcl::rising_factorial(const S1& x, Integer i)
{
    static_assert(md::is_complex<S1>::value == false, "complex scalars not allowed");

    using S = typename md::unify_types_promote<S1,Float>::type;
    return details::combinatorics_helper<S>::eval_rising_factorial(S(x),i);
};

template<class S1, class Enable>
typename md::unify_types_promote<S1,Float>::type
matcl::falling_factorial(const S1& x, Integer i)
{
    static_assert(md::is_complex<S1>::value == false, "complex scalars not allowed");

    using S = typename md::unify_types_promote<S1,Float>::type;
    return details::combinatorics_helper<S>::eval_falling_factorial(S(x),i);
};

template<class S1, class Enable>
typename md::unify_types_promote<S1,Float>::type
matcl::binomial_coefficient(Integer n, Integer k)
{
    using S = typename md::unify_types_promote<S1,Float>::type;
    return details::combinatorics_helper<S>::eval_binomial_coefficient(n,k);
};

inline Real matcl::binomial_coefficient(Integer n, Integer k)
{
    return details::combinatorics_helper<Real>::eval_binomial_coefficient(n,k);
};
inline Float matcl::fbinomial_coefficient(Integer n, Integer k)
{
    return details::combinatorics_helper<Float>::eval_binomial_coefficient(n,k);
};

};