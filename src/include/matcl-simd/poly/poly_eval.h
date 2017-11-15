/*
 *  This file is a part of Matrix Computation Library (MATCL)
 *
 *  Copyright (c) Pawe³ Kowal 2017
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

#include "matcl-simd/simd.h"
#include "matcl-core/float/twofold.h"

namespace matcl { namespace simd
{

//-----------------------------------------------------------------------
//                      EVALUATION OF POLYNOMIALS
//-----------------------------------------------------------------------
// evaluation of a polynomial P(x) = sum_{i = 0}^{N - 1} a_i * x^i

// Template arguments:
//  N           - number of coefficients of the polynomial
//  Arg_type    - floating point type (float, double, or a simd type)
//  Coef_type   - floating point type of polynomial coefficients,
//                  this type must be convertible to Arg_type

// evaluate a polynomial at point x using the Horner's scheme
//      res = small_horner(x, a0, a1, ..., a_N-1)
// where ai is the i-th coefficientof the polynomial
//
// Template arguments:
//  Arg_type    - floating point type (float, double, or a simd type)
//  Args...     - floating point type of a polynomial coefficient,
//                  this type must be convertible to Arg_type
template<class Arg_type, class ... Args>
force_inline
Arg_type    small_horner(const Arg_type& x, const Args& ... coef);

// evaluate a polynomial at point x using the Horner's scheme
// polynomial is represented as an array of size N:
//      poly = {a_0, a_1, ..., a_{N-1}}
template<int N, class Arg_type, class Coef_type>
force_inline
Arg_type    horner(const Arg_type& x, const Coef_type* poly);

// evaluate a polynomial at point x using the Horner's scheme
// polynomial is represented as an array of size N:
//      poly = {a_0, a_1, ..., a_{N-1}}
template<class Arg_type, class Coef_type>
Arg_type    horner(const Arg_type& x, int N, const Coef_type* poly);

// evaluate a polynomial at point x using the Estrin's scheme
// polynomial is represented as an array of size N:
//      poly = {a_0, a_1, ..., a_{N-1}}
//
// Note:
//  this function is generally faster than horner function for 
//  sufficiently large polynomials (N > 10), but is less accurate
template<int N, class Arg_type, class Coef_type>
force_inline
Arg_type    estrin(const Arg_type& x, const Coef_type* poly);

// evaluate a polynomial at point x using the Horner's scheme
// polynomial is represented as an array of size N:
//      poly = {a_0, a_1, ..., a_{N-1}}
//
// Note:
//  this function can be must faster than horner function for 
//  sufficiently large polynomials (N > 100), but is less accurate.
template<class Arg_type, class Coef_type>
Arg_type    estrin(const Arg_type& x, int N, const Coef_type* poly);

// evaluate a polynomial at point x using the compensated Horner's scheme
// polynomial is represented as an array of size N:
//      poly = {a_0, a_1, ..., a_{N-1}}
//
// Note:
//  this function is evaluates a polynomial with high accuracy:
//      |p(x) - p_ap(x)|/|p(x)| <= u + phi * u^2
// where p(x) is the true value, p_ap(x) is computed, u is the unit
// roundoff, phi = b(n)*cond(p,x) * u^2, b(n) ~ 4n^2, and cond(p,x) is the
// condition number of p as returned by horner_apriori_cond function.
// Thus, if p is not too ill-conditioned, then p_app(x) is calculated with
// full full accuracy (0.5 ulp error). This function is however more costly,
// than the horner function (3-8 times slower)
//
// References:
//  [1]. Faithful Polynomial Evaluation with Compensated Horner Algorithm,
//      P. Langlois, N. Louvet, 2006
template<class Arg_type, class Coef_type>
Arg_type    compensated_horner(const Arg_type& x, int N, const Coef_type* poly);

// return the condition number of a polynomial at point x evaluated using
// the Horner's method. 
// The relative accuracy of computed polynomial p_ap(x) is given by
//      |p(x) - p_ap(x)| / |p(x)| <= a(2N) * cond(p, x) * u + O(u^2)
// where p(x) is the true value, cond(p, x) is returned condition number
// u is the unit roundoff (half of the machine epsilon), and a(2N) = 2N/(1-2Nu) ~ 2N.
//
// The condition number cond(p, x) is given by:
//      cond(p, x) = sum_{i = 0}^{N - 1} |a_i| * |x|^i / |p(x)| 
//                := |p|(|x|) / p(x)
//
// Note:
//  this is an a priori bound and so takes no account of the actual rounding
//  errors that occur, therefore can be too pessimistic.
template<class Arg_type, class Coef_type>
Arg_type    horner_apriori_cond(const Arg_type& x, int N, const Coef_type* poly);

// return the condition number of a polynomial at point x evaluated using
// the Horner's method as in function horner_apriori_cond; this function
// additionally returns |p|(|x|) in 'val_abs' argument and p(x) in 'val'
// argument
template<class Arg_type, class Coef_type>
Arg_type    horner_apriori_cond(const Arg_type& x, int N, const Coef_type* poly,
                                Arg_type& val, Arg_type& val_abs);

// evaluate a polynomial at point x using the Horner's scheme
// polynomial is represented as an array of size N:
//      poly = {a_0, a_1, ..., a_{N-1}}
// this function also returns aposteriori absolute forward error
// estimator 'error', such that:
//      |res - p(x) | < 'error'
// where res is the computed value, p(x) is true value
template<class Arg_type, class Coef_type>
Arg_type    horner_and_error(const Arg_type& x, int N, const Coef_type* poly,
                Arg_type& error);

// evaluate a polynomial at point x using the compensated Horner's scheme
// as in case of the function compensated_horner function, but slightly
// more accurately; polynomial is represented as an array of size N:
//      poly = {a_0, a_1, ..., a_{N-1}}
// this function also returns aposteriori absolute forward error
// estimator 'error', such that:
//      |res - p(x) | < 'error'
// where res is the computed value, p(x) is true value
// if returned value of 'is_exactly_rounded' argument is true, then the
// result is proven to be correct up to 1/2 ulp.
//
// References:
//  [1]. Faithful Polynomial Evaluation with Compensated Horner Algorithm,
//      P. Langlois, N. Louvet, 2006
template<class Arg_type, class Coef_type>
Arg_type    compensated_horner_and_error(const Arg_type& x, int N, const Coef_type* poly,
                Arg_type& error, bool& is_exactly_rounded);

}};

#include "matcl-simd/details/poly/poly_eval.inl"
