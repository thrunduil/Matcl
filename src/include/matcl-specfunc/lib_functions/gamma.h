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

#pragma once

#include "matcl-specfunc/general/config.h"
#include "matcl-matrep/matrix/matrix.h"

namespace matcl
{

namespace md = matcl::details;

// gamma function of the elements of given matrix;
// gamma function for scalar x is defined as
//     gamma(x) = int_0^inf t^(x-1) * exp(-t) dt
MATCL_SF_EXPORT Matrix	    gamma(const Matrix &m);

template<class S1, class Enable = typename md::enable_if_scalar<S1,void>::type>
typename md::unify_types_promote<S1,Float>::type
                            gamma(const S1& x);

// natural logarithm of the gamma function at the elements of given matrix m.
//     gammaln(x) = log gamma(x)
MATCL_SF_EXPORT Matrix	    gammaln(const Matrix &m);

template<class S1, class Enable = typename md::enable_if_scalar<S1,void>::type>
typename md::unify_types_promote<S1,Float>::type
                            gammaln(const S1& x);

// version of the function gammaln, which returns additionally sign of gamma
// function evaluated at given point:
//     sign = isign(gamma(x))
// not available for complex values
template<class S1, class Enable = typename md::enable_if_scalar<S1,void>::type>
typename md::unify_types_promote<S1,Float>::type
                            gammaln(const S1& x, int& sign);

// evaluate function:
//     gamma1pm1(x) = gamma(x + 1) - 1
template<class S1, class Enable = typename md::enable_if_scalar<S1,void>::type>
typename md::unify_types_promote<S1,Float>::type
                            gamma1pm1(const S1& x);

// digamma (or phi) function of the elements of given matrix;
// digamma is defined as the logarithmic derivative of the gamma function:
//     digamma(x) = d/dx[ln(gamma(x))] = d/dx[gamma(x)] / gamma(x)
// not available for complex values
MATCL_SF_EXPORT Matrix	    digamma(const Matrix &m);

template<class S1, class Enable = typename md::enable_if_scalar<S1,void>::type>
typename md::unify_types_promote<S1,Float>::type
                            digamma(const S1& x);

// trigamma function of the elements of given matrix;
// trigamma is defined as the derivative of the digamma function:
//     trigamma(x) = d/dx[digamma(x)]
// not available for complex values
MATCL_SF_EXPORT Matrix	    trigamma(const Matrix &m);

template<class S1, class Enable = typename md::enable_if_scalar<S1,void>::type>
typename md::unify_types_promote<S1,Float>::type
                            trigamma(const S1& x);

// polygamma function;
// polygamma is defined as the n-th derivative of the digamma function:
//     polygamma(x,n) = d^n/dx_n[digamma(x)]
// not available for complex values
MATCL_SF_EXPORT Matrix	    polygamma(const Matrix &m, Integer n);

template<class S1, class Enable = typename md::enable_if_scalar<S1,void>::type>
typename md::unify_types_promote<S1,Float>::type
                            polygamma(const S1& x, Integer n);

// return the ratio of gamma functions:
//     gamma_ratio(x,y) = gamma(x) / gamma(y)
// not available for complex values
template<class S1, class S2, class Enable = typename md::enable_if_scalar2<S1,S2,void>::type>
typename md::unify_types2_promote<S1,S2,Float>::type
                            gamma_ratio(const S1& x, const S2& y);

// return the ratio of gamma functions:
//     gamma_ratio(x,delta) = gamma(x) / gamma(x + delta)
// not available for complex values
template<class S1, class S2, class Enable = typename md::enable_if_scalar2<S1,S2,void>::type>
typename md::unify_types2_promote<S1,S2,Float>::type
                            gamma_delta_ratio(const S1& x, const S2& delta);

// return lower incomplete gamma function defined as:
//     igamma_lower(a,z) = int_0^z t^(a-1) exp(-t) dt
// not available for complex values
template<class S1, class S2, class Enable = typename md::enable_if_scalar2<S1,S2,void>::type>
typename md::unify_types2_promote<S1,S2,Float>::type
                            igamma_lower(const S1& a, const S2& z);

// return upper incomplete gamma function defined as:
//     igamma_upper(a,z) = int_z^infty t^(a-1) exp(-t) dt
// not available for complex values
template<class S1, class S2, class Enable = typename md::enable_if_scalar2<S1,S2,void>::type>
typename md::unify_types2_promote<S1,S2,Float>::type
                            igamma_upper(const S1& a, const S2& z);

// return normalized lower incomplete gamma function defined as:
//     igamma_lower_norm(a,z) = igamma_lower(a,z)/gamma(a)
// not available for complex values
template<class S1, class S2, class Enable = typename md::enable_if_scalar2<S1,S2,void>::type>
typename md::unify_types2_promote<S1,S2,Float>::type
                            igamma_lower_norm(const S1& a, const S2& z);

// return normalized upper incomplete gamma function defined as:
//     igamma_upper_norm(a,z) = igamma_upper(a,z)/gamma(a)
// not available for complex values
template<class S1, class S2, class Enable = typename md::enable_if_scalar2<S1,S2,void>::type>
typename md::unify_types2_promote<S1,S2,Float>::type
                            igamma_upper_norm(const S1& a, const S2& z);

// return inverse of normalized lower incomplete gamma function defined as:
//     x = igamma_lower_inv(a,q) <=> q = igamma_lower_norm(a,x)
// not available for complex values
template<class S1, class S2, class Enable = typename md::enable_if_scalar2<S1,S2,void>::type>
typename md::unify_types2_promote<S1,S2,Float>::type
                            igamma_lower_inv(const S1& a, const S2& q);

// return inverse of normalized upper incomplete gamma function defined as:
//     x = igamma_upper_inv(a,q) <=> q = igamma_upper_norm(a,x)
// not available for complex values
template<class S1, class S2, class Enable = typename md::enable_if_scalar2<S1,S2,void>::type>
typename md::unify_types2_promote<S1,S2,Float>::type
                            igamma_upper_inv(const S1& a, const S2& q);

// return inverse of normalized lower incomplete gamma function with respect
// to parameter 'a' defined as:
//     a = igamma_lower_inva(x,q) <=> q = igamma_lower_norm(a,x)
// not available for complex values
template<class S1, class S2, class Enable = typename md::enable_if_scalar2<S1,S2,void>::type>
typename md::unify_types2_promote<S1,S2,Float>::type
                            igamma_lower_inva(const S1& x, const S2& q);

// return inverse of normalized upper incomplete gamma function with respect
// to parameter 'a' defined as:
//     a = igamma_upper_inva(x,q) <=> q = igamma_upper_norm(a,x)
// not available for complex values
template<class S1, class S2, class Enable = typename md::enable_if_scalar2<S1,S2,void>::type>
typename md::unify_types2_promote<S1,S2,Float>::type
                            igamma_upper_inva(const S1& x, const S2& q);

// return partial derivative of normalized lower incomplete gamma function defined as:
//     igamma_lower_dif(a,x) = d/dx[igamma_lower_norm(a,x)] = x^(a-1)exp(-x)/gamma(a)
// not available for complex values
template<class S1, class S2, class Enable = typename md::enable_if_scalar2<S1,S2,void>::type>
typename md::unify_types2_promote<S1,S2,Float>::type
                            igamma_lower_dif(const S1& a, const S2& x);

};

#include "matcl-specfunc/details/gamma.inl"