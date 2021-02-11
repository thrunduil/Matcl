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

#include "matcl-scalar/config.h"
#include "matcl-scalar/object.h"
#include "matcl-core/matrix/scalar_types.h"
#include "matcl-scalar/details/utils.h"
#include "matcl-scalar/details/enablers.h"

namespace matcl
{

namespace md = matcl::details;

// complex conversion function:
// for a function fun defined for real numbers in a set D and defined
// for complex numbers, the complex conversion function (denoted by func_c)
// operating on a matrix A is defined as follows:
//     if each element of A is in the domain of D, then function func
//         is called on every element of A
//     otherwise A is converted to a complex matrix and the function func_c
//         is caled on every element of the converted matrix

//---------------------------------------------------------------------------
//                      VALUE CLASSIFICATION
//---------------------------------------------------------------------------

// return type of floating point number; not defined for complex values
template<class S, class Enable = typename md::enable_if_scalar_ntobj<S,void>::type>
fp_type                 fpclassify(const S& x);

// check if value is finite
template<class S, class Enable = typename md::enable_if_scalar_ntobj<S,void>::type>
typename md::bool_or_object<S>::type
                        is_finite(const S& x);

// check if value is infinite
template<class S, class Enable = typename md::enable_if_scalar_ntobj<S,void>::type>
typename md::bool_or_object<S>::type
                        is_nan(const S& x);

// check if value is infinite
template<class S, class Enable = typename md::enable_if_scalar_ntobj<S,void>::type>
typename md::bool_or_object<S>::type
                        is_inf(const S& x);

// check if value is regular (i.e., neither NaN, nor an infinity nor zero)
// for a complex value false is retured if real or imaginary part is
// NaN or infinity or if both real and imaginary part is zero
template<class S, class Enable = typename md::enable_if_scalar_ntobj<S,void>::type>
typename md::bool_or_object<S>::type
                        is_regular(const S& x);

// check if value is regular (i.e., neither NaN, nor an infinity nor zero
// nor subnormal)
template<class S, class Enable = typename md::enable_if_scalar_ntobj<S,void>::type>
typename md::bool_or_object<S>::type
                        is_normal(const S& x);

// check if value is an integer
template<class S, class Enable = typename md::enable_if_scalar_ntobj<S,void>::type>
typename md::bool_or_object<S>::type
                        is_int(const S& x);

// check if value is real; for a complex value return true if the imaginary
// part is zero
template<class S, class Enable = typename md::enable_if_scalar_ntobj<S,void>::type>
typename md::bool_or_object<S>::type
                        is_real(const S& x);

// check if value is zero; not defined for a matrix
template<class S, class Enable = typename md::enable_if_scalar_ntobj<S,void>::type>
bool                    is_zero(const S& x);

// check if value is one; not defined for a matrix
template<class S, class Enable = typename md::enable_if_scalar_ntobj<S,void>::type>
bool                    is_one(const S& x);

//---------------------------------------------------------------------------
//                      LOGICAL OPERATIONS
//---------------------------------------------------------------------------

template<class S, class Enable = typename md::enable_if_scalar_ntobj<S,void>::type>
bool                    cast_bool(const S& x);

template<class S, class Enable = typename md::enable_if_scalar_ntobj<S,void>::type>
bool                    op_not(const S& x);

// return logical negation of x; this is a different name for the op_not
// function; 
template<class S, class Enable = typename md::enable_if_scalar_ntobj<S,void>::type>
bool                    operator!(const S& x)           { return op_not(x); };

// call the element-wise negation operator (i.e. ~x); for scalar types this
// function is equivalent to !x (even for Integer types)
template<class S, class Enable = typename md::enable_if_scalar_ntobj<S,void>::type>
typename md::bool_or_object<S>::type
                        neg(const S& x);

// call the element-wise negation operator (i.e. ~x); this is a different 
// name for the neg function
template<class S, class Enable = typename md::enable_if_scalar_ntobj<S,void>::type>
typename md::bool_or_object<S>::type
                        is_false(const S& x)            { return neg(x); };

// call the element-wise cast bool operator; for scalar types this function 
// is equivalent to (bool)x
template<class S, class Enable = typename md::enable_if_scalar_ntobj<S,void>::type>
typename md::bool_or_object<S>::type
                        is_true(const S& x);

//---------------------------------------------------------------------------
//                      ARITHMETIC FUNCTIONS
//---------------------------------------------------------------------------

// unary minus, i.e. -A; different name of the operator- function
template<class S, class Enable = typename md::enable_if_scalar_ntobj<S,void>::type>
typename md::promote_scalar<S>::type
                        uminus(const S& x);

// return unary minus, i.e -x;
template<class S, class Enable = typename md::enable_if_scalar_ntobj<S,void>::type> inline
typename md::promote_scalar<S>::type
                        operator-(const S& x)           { return uminus(x); };

// unary plus operator; always return the input value x
template<class S, class Enable = typename md::enable_if_scalar_ntobj<S,void>::type> inline
typename md::promote_scalar<S>::type
                        operator+(const S& x)           { return x; };

// inverse of x, i.e. 1/x; if no errors occure, then mul(x, invs(x)) = 1
template<class S, class Enable = typename md::enable_if_scalar_ntobj<S,void>::type>
typename md::unify_types_promote<S, Float>::type
                        invs(const S& x);

// inverse of x, i.e. 1/x; if no errors occur, then mmul(x, inv(x)) = 1;
// for all scalar types and most of objects this function is equivalent to
// invs (and mul is equivalent to mmul); however for OMatrix this is a 
// matrix inversion
template<class S, class Enable = typename md::enable_if_scalar_ntobj<S,void>::type>
typename md::unify_types_promote<S, Float>::type
                        inv(const S& x);

// return imaginary part of a value x
template<class S, class Enable = typename md::enable_if_scalar_ntobj<S,void>::type>
typename md::real_type_promote<S>::type
                        imag(const S& x);

// return real part of a value x
template<class S, class Enable = typename md::enable_if_scalar_ntobj<S,void>::type>
typename md::real_type_promote<S>::type
                        real(const S& x);

// return the absolute value of x
template<class S, class Enable = typename md::enable_if_scalar_ntobj<S,void>::type>
typename md::real_type_promote<S>::type
                        abs(const S& x);

// return the absolute value squared of x
template<class S, class Enable = typename md::enable_if_scalar_ntobj<S,void>::type>
typename md::real_type_promote<S>::type
                        abs2(const S& x);

// return argument (polar angle) of a complex number
template<class S, class Enable = typename md::enable_if_scalar_ntobj<S,void>::type>
typename md::real_type_promote_unify<S,Float>::type
                        arg(const S& x);

// different name of the arg function
template<class S, class Enable = typename md::enable_if_scalar_ntobj<S,void>::type>
typename md::real_type_promote_unify<S,Float>::type
                        angle(const S& x);

// return the conjugate transposition of x
template<class S, class Enable = typename md::enable_if_scalar_ntobj<S,void>::type>
typename md::promote_scalar<S>::type
                        conj(const S& x);

//---------------------------------------------------------------------------
//               EXPONENTIAL, LOGARITHMIC AND POWER FUNCTIONS
//---------------------------------------------------------------------------

// return the square root of x
template<class S, class Enable = typename md::enable_if_scalar_ntobj<S,void>::type>
typename md::unify_types_promote<S,Float>::type
                         sqrt(const S& x);

// return the square root of x converted to a complex value
template<class S, class Enable = typename md::enable_if_scalar_ntobj<S,void>::type>
typename md::unify_types_promote<S,Float_complex>::type
                        sqrt_c(const S& x);

// return sqrt(1+x) - 1; more accurate than sqrt for x close to 0
template<class S, class Enable = typename md::enable_if_scalar_ntobj<S,void>::type>
typename md::unify_types_promote<S,Float>::type
                        sqrt1pm1(const S& x);

// return the value of sqrt1pm1 function for x converted to a complex value
template<class S, class Enable = typename md::enable_if_scalar_ntobj<S,void>::type>
typename md::unify_types_promote<S,Float_complex>::type
                        sqrt1pm1_c(const S& x);

// return the cubic root of x, for negative values result
// is also negative; not available for complex values
template<class S, class Enable = typename md::enable_if_scalar_ntobj<S,void>::type>
typename md::unify_types_promote<S,Float>::type
                        cbrt(const S& x);

// return the exponential function, i.e. the e (Euler's number, 2.7182818)
// raised to the given power x; 
// for complex z = x + i*y, exp(z) = exp(x) * (cos(y) + i*sin(Y))
template<class S, class Enable = typename md::enable_if_scalar_ntobj<S,void>::type>
typename md::unify_types_promote<S,Float>::type
                        exp(const S& x);

// return 2^x
template<class S, class Enable = typename md::enable_if_scalar_ntobj<S,void>::type>
typename md::unify_types_promote<S,Float>::type
                        exp2(const S& x);

// different name of the exp2 function
template<class S, class Enable = typename md::enable_if_scalar_ntobj<S,void>::type>
typename md::unify_types_promote<S,Float>::type
                        pow2(const S& x)        { return exp2(x); }

// return 10^x
template<class S, class Enable = typename md::enable_if_scalar_ntobj<S,void>::type>
typename md::unify_types_promote<S,Float>::type
                        exp10(const S& x);

// different name of the exp10 function
template<class S, class Enable = typename md::enable_if_scalar_ntobj<S,void>::type>
typename md::unify_types_promote<S,Float>::type
                        pow10(const S& x)       { return exp10(x); };

// return exp(x) - 1; this function is more accurate for x near zero
template<class S, class Enable = typename md::enable_if_scalar_ntobj<S,void>::type>
typename md::unify_types_promote<S,Float>::type
                        expm1(const S& x);

// return exp(x * 1i), where 1i is an imaginary unit
template<class S, class Enable = typename md::enable_if_scalar_ntobj<S,void>::type>
typename md::unify_types_promote<S,Float_complex>::type
                        expi(const S& x);

// scale x by 2 raised to the power of n, returning ldexp(x,n) = x * 2^n
template<class S1, class Enable = typename md::enable_if_scalar_ntobj<S1,void>::type>
typename md::unify_types_promote<S1,Float>::type
                        ldexp(const S1& x, Integer exp);

// function equivalent to scalbn
template<class S1, class Enable = typename md::enable_if_scalar_ntobj<S1,void>::type>
typename md::unify_types_promote<S1,Float>::type
                        scalbn(const S1& x, Integer n);

// computes the natural (base e, Euler's number, 2.7182818) logarithm 
// of x; for complex argument z = x + i*y compute logarithm of a complex 
// value z with a branch cut along the negative real axis; if no errors
// occur, the complex natural logarithm of z is returned, in the range of
// a strip in the interval [-i pi, +i pi] along the imaginary axis and 
// mathematically unbounded along the real axis;
// Log(z) = log(|z|) + i Arg(z) = log( hypot(x,y) ) + i atan2(y,x)
template<class S, class Enable = typename md::enable_if_scalar_ntobj<S,void>::type>
typename md::unify_types_promote<S,Float>::type
                        log(const S& x);

// return the value of log function for x converted to a complex value
template<class S, class Enable = typename md::enable_if_scalar_ntobj<S,void>::type>
typename md::unify_types_promote<S,Float_complex>::type
                        log_c(const S& x);

// compute log(1+x); more accurate for m close to 0 
template<class S, class Enable = typename md::enable_if_scalar_ntobj<S,void>::type>
typename md::unify_types_promote<S,Float>::type
                        log1p(const S& x);

// return the value of log1p function for x converted to a complex value
template<class S, class Enable = typename md::enable_if_scalar_ntobj<S,void>::type>
typename md::unify_types_promote<S,Float_complex>::type
                        log1p_c(const S& x);

// computes the base 2 logarithm of x
template<class S, class Enable = typename md::enable_if_scalar_ntobj<S,void>::type>
typename md::unify_types_promote<S,Float>::type
                        log2(const S& x);

// return the value of log2 function for x converted to a complex value
template<class S, class Enable = typename md::enable_if_scalar_ntobj<S,void>::type>
typename md::unify_types_promote<S,Float_complex>::type
                        log2_c(const S& x);

// computes the base 10 logarithm of x
template<class S, class Enable = typename md::enable_if_scalar_ntobj<S,void>::type>
typename md::unify_types_promote<S,Float>::type
                        log10(const S& x);

// return the value of log10 function for x converted to a complex value
template<class S, class Enable = typename md::enable_if_scalar_ntobj<S,void>::type>
typename md::unify_types_promote<S,Float_complex>::type
                        log10_c(const S& x);

// extract the value of the unbiased exponent from the floating-point argument
// x, and returns it as a floating-point value; equivalent to log2(|x|) rounded
// to integer number toward -INF; additionally for regular values 
// logb(x) = exp - 1, where exp is the exponent returned by frexp, 
// logb(0) = -inf, logb(+- inf) = inf, logb(nan) = nan; for complex argument 
// logb(|x|) is returned
template<class S, class Enable = typename md::enable_if_scalar_ntobj<S,void>::type>
typename md::real_unify_types_promote<S,Float>::type
                        logb(const S& x);

// cast result of logb to Integer
template<class S, class Enable = typename md::enable_if_scalar_ntobj<S,void>::type>
typename md::integer_or_object<S>::type
                        ilogb(const S& x);

//---------------------------------------------------------------------------
//                      TRIGONOMETRIC FUNCTIONS
//---------------------------------------------------------------------------

// return the sine function of x in radians;
// for complex z = x + i*y, sin(z) = sin(x) * cosh(y) + i cos(x) * sinh(y)
template<class S, class Enable = typename md::enable_if_scalar_ntobj<S,void>::type>
typename md::unify_types_promote<S,Float>::type
                        sin(const S& x);

// return the cosine function of x in radians;
// for complex z = x + i*y, cos(z) = cos(x) * cosh(y) - i sin(x) * sinh(y)
template<class S, class Enable = typename md::enable_if_scalar_ntobj<S,void>::type>
typename md::unify_types_promote<S,Float>::type
                        cos(const S& x);

// return the tangent function of x in radians;
// for complex z = x + i*y, tan(z) = -i * tanh(i*z)
template<class S, class Enable = typename md::enable_if_scalar_ntobj<S,void>::type>
typename md::unify_types_promote<S,Float>::type
                        tan(const S& x);

// return the cotangent function of x in radians;
// for complex z = x + i*y, cot(z) = i * coth(i*z)
template<class S, class Enable = typename md::enable_if_scalar_ntobj<S,void>::type>
typename md::unify_types_promote<S,Float>::type
                        cot(const S& x);

// return the secant function of x in radians;
// for complex z = x + i*y, sec(z) = 1 / cos(z)
template<class S, class Enable = typename md::enable_if_scalar_ntobj<S,void>::type>
typename md::unify_types_promote<S,Float>::type
                        sec(const S& x);

// return the cosecant function of x in radians;
// for complex z = x + i*y, csc(z) = 1 / sin(z)
template<class S, class Enable = typename md::enable_if_scalar_ntobj<S,void>::type>
typename md::unify_types_promote<S,Float>::type
                        csc(const S& x);

// return the hyperbolic sine function of x in radians;
// for complex z = x + i*y, sinh(z) = sinh(x) * cos(y) + i cosh(x) * sin(y)
template<class S, class Enable = typename md::enable_if_scalar_ntobj<S,void>::type>
typename md::unify_types_promote<S,Float>::type
                        sinh(const S& x);

// return the hyperbolic cosine function of x in radians;
// for complex z = x + i*y, cosh(z) = cosh(x) * cos(y) + i sinh(x) * sin(y)
template<class S, class Enable = typename md::enable_if_scalar_ntobj<S,void>::type>
typename md::unify_types_promote<S,Float>::type
                        cosh(const S& x);

// return the hyperbolic tangent function of x in radians;
// for complex z = x + i*y, tanh(z) = sinh(z) / cosh(z)
template<class S, class Enable = typename md::enable_if_scalar_ntobj<S,void>::type>
typename md::unify_types_promote<S,Float>::type
                        tanh(const S& x);

// return the hyperbolic cotangent function of x in radians;
// for complex z = x + i*y, coth(z) = cosh(z) / sinh(z)
template<class S, class Enable = typename md::enable_if_scalar_ntobj<S,void>::type>
typename md::unify_types_promote<S,Float>::type
                        coth(const S& x);

// return the hyperbolic secant function of x in radians;
// for complex z = x + i*y, sech(z) = 1 / cosh(z)
template<class S, class Enable = typename md::enable_if_scalar_ntobj<S,void>::type>
typename md::unify_types_promote<S,Float>::type
                        sech(const S& x);

// return the hyperbolic cosecant function of x in radians;
// for complex z = x + i*y, csch(z) = 1 / sinh(z)
template<class S, class Enable = typename md::enable_if_scalar_ntobj<S,void>::type>
typename md::unify_types_promote<S,Float>::type
                        csch(const S& x);

// return the inverse sine function of x, result in radians;
// for complex z:
// asin z = 1/i log(iz + sqrt(|1 - z^2|) * exp(i/2 arg(1-z^2)))
template<class S, class Enable = typename md::enable_if_scalar_ntobj<S,void>::type>
typename md::unify_types_promote<S,Float>::type
                        asin(const S& x);

// return the value of asin function for x converted to a complex value
template<class S, class Enable = typename md::enable_if_scalar_ntobj<S,void>::type>
typename md::unify_types_promote<S,Float_complex>::type
                        asin_c(const S& x);

// return the inverse cosine function of x, result in radians
// for complex z:
// acos z = 1/i log(z + i * sqrt(|1 - z^2|) * exp(i/2 arg(1-z^2)))
template<class S, class Enable = typename md::enable_if_scalar_ntobj<S,void>::type>
typename md::unify_types_promote<S,Float>::type
                        acos(const S& x);

// return the value of acos function for x converted to a complex value
template<class S, class Enable = typename md::enable_if_scalar_ntobj<S,void>::type>
typename md::unify_types_promote<S,Float_complex>::type
                        acos_c(const S& x);

// return the inverse tangent function of x, result in radians;
// for a complex z, atan(z) = -i atanh(iz)
template<class S, class Enable = typename md::enable_if_scalar_ntobj<S,void>::type>
typename md::unify_types_promote<S,Float>::type
                        atan(const S& x);

// return the inverse cotangent function of x, result in radians;
// for a complex z, acot(z) = atan(1/z)
template<class S, class Enable = typename md::enable_if_scalar_ntobj<S,void>::type>
typename md::unify_types_promote<S,Float>::type
                        acot(const S& x);

// return the inverse secant function of x, result in radians;
// for a complex z, asec(z) = acos(1/z)
template<class S, class Enable = typename md::enable_if_scalar_ntobj<S,void>::type>
typename md::unify_types_promote<S,Float>::type
                        asec(const S& x);

// return the value of asec function for x converted to a complex value
template<class S, class Enable = typename md::enable_if_scalar_ntobj<S,void>::type>
typename md::unify_types_promote<S,Float_complex>::type
                        asec_c(const S& x);

// return the inverse cosecant function of x, result in radians;
// for a complex z, acsc(z) = asin(1/z)
template<class S, class Enable = typename md::enable_if_scalar_ntobj<S,void>::type>
typename md::unify_types_promote<S,Float>::type
                        acsc(const S& x);

// return the value of acsc function for x converted to a complex value
template<class S, class Enable = typename md::enable_if_scalar_ntobj<S,void>::type>
typename md::unify_types_promote<S,Float_complex>::type
                        acsc_c(const S& x);

// return the inverse hyperbolic sine function of x, result in radians
// for a complex z, asinh(z) = i asin(-i z)
template<class S, class Enable = typename md::enable_if_scalar_ntobj<S,void>::type>
typename md::unify_types_promote<S,Float>::type
                        asinh(const S& x);

// return the inverse hyperbolic cosine function of x, result in radians
// for a complex z, acosh(z) = +-i acos(z) with posive sign if imag(acos(z)) < 0,
// thus real( acosh(z) ) >= 0
template<class S, class Enable = typename md::enable_if_scalar_ntobj<S,void>::type>
typename md::unify_types_promote<S,Float>::type
                        acosh(const S& x);

// return the value of acosh function for x converted to a complex value
template<class S, class Enable = typename md::enable_if_scalar_ntobj<S,void>::type>
typename md::unify_types_promote<S,Float_complex>::type
                        acosh_c(const S& x);

// return the inverse hyperbolic tangent function of x, result in radians
// for complex z = x + i*y, atanh(z) = re + im *i, where
// re = log1p(4 * x / ( (x - 1) * (x - 1) + y^2))
// im = atan(2y, (1 - x) * (1 + x) - y^2 ) / 2
template<class S, class Enable = typename md::enable_if_scalar_ntobj<S,void>::type>
typename md::unify_types_promote<S,Float>::type
                        atanh(const S& x);

// return the value of atanh function for x converted to a complex value
template<class S, class Enable = typename md::enable_if_scalar_ntobj<S,void>::type>
typename md::unify_types_promote<S,Float_complex>::type
                        atanh_c(const S& x);

// return the inverse hyperbolic cotangent function of x, result in radians;
// for a complex z, acoth(z) = atanh(1/z)
template<class S, class Enable = typename md::enable_if_scalar_ntobj<S,void>::type>
typename md::unify_types_promote<S,Float>::type
                        acoth(const S& x);

// return the value of acoth function for x converted to a complex value
template<class S, class Enable = typename md::enable_if_scalar_ntobj<S,void>::type>
typename md::unify_types_promote<S,Float_complex>::type
                        acoth_c(const S& x);

// return the inverse hyperbolic secant function of x, result in radians;
// for a complex z, asech(z) = acosh(1/z)
template<class S, class Enable = typename md::enable_if_scalar_ntobj<S,void>::type>
typename md::unify_types_promote<S,Float>::type
                        asech(const S& x);

// return the value of asech function for x converted to a complex value
template<class S, class Enable = typename md::enable_if_scalar_ntobj<S,void>::type>
typename md::unify_types_promote<S,Float_complex>::type
                        asech_c(const S& x);

// return the inverse hyperbolic cosecant function of x, result in radians;
// for a complex z, acsch(z) = asinh(1/z)
template<class S, class Enable = typename md::enable_if_scalar_ntobj<S,void>::type>
typename md::unify_types_promote<S,Float>::type
                        acsch(const S& x);

//---------------------------------------------------------------------------
//                      ROUNDING FUNCTIONS
//---------------------------------------------------------------------------

// return epsilon value eps, i.e positive distance from abs(x) to the next
// larger in magnitude floating point number of the same precision as x
template<class S, class Enable = typename md::enable_if_scalar_ntobj<S,void>::type>
typename md::real_unify_types_promote<S,Float>::type
                        eps(const S& x);

// round to nearest integer towards -INF; for complex numbers
// real and imaginary part are rounded separately
template<class S, class Enable = typename md::enable_if_scalar_ntobj<S,void>::type>
typename md::promote_scalar<S>::type
                        floor(const S& x);

// round to nearest integer towards +INF; for complex numbers
// real and imaginary part are rounded separately
template<class S, class Enable = typename md::enable_if_scalar_ntobj<S,void>::type>
typename md::promote_scalar<S>::type
                        ceil(const S& x);

// round to nearest integer, rounding halfway cases to nearest even integer;
// for complex numbers real and imaginary part are rounded separately
template<class S, class Enable = typename md::enable_if_scalar_ntobj<S,void>::type>
typename md::promote_scalar<S>::type
                        round(const S& x);

// for each element in the matrix m of size M x N call the trunc function; 
// return an M x N matrix 
template<class S, class Enable = typename md::enable_if_scalar_ntobj<S,void>::type>
typename md::promote_scalar<S>::type
                        trunc(const S& x);

// different name of the trunc function
template<class S, class Enable = typename md::enable_if_scalar_ntobj<S,void>::type>
typename md::promote_scalar<S>::type
                        fix(const S& x)         { return trunc(x); };

// round to nearest integer towards -INF and cast to integer;
// not available for complex numbers
template<class S, class Enable = typename md::enable_if_scalar_ntobj<S,void>::type>
typename md::integer_or_object<S>::type
                        ifloor(const S& x);

// round to nearest integer towards +INF and cast to integer;
// not available for complex numbers
template<class S, class Enable = typename md::enable_if_scalar_ntobj<S,void>::type>
typename md::integer_or_object<S>::type
                        iceil(const S& x);

// round to nearest integer, rounding halfway cases away from zero
// not available for complex numbers
template<class S, class Enable = typename md::enable_if_scalar_ntobj<S,void>::type>
typename md::integer_or_object<S>::type
                        iround(const S& x);

// round to nearest integer towards zero and cast to integer; 
// not available for complex numbers
template<class S, class Enable = typename md::enable_if_scalar_ntobj<S,void>::type>
typename md::integer_or_object<S>::type
                        itrunc(const S& x);

// different name of the itrunc function
template<class S, class Enable = typename md::enable_if_scalar_ntobj<S,void>::type>
typename md::integer_or_object<S>::type
                        ifix(const S& x)        { return itrunc(x); };

//---------------------------------------------------------------------------
//                 VALUE DECOMPOSITION FUNCTIONS
//---------------------------------------------------------------------------

// the sign of x (-1 for negative, 1 for positive, 0 for zero, NaN if x is
// NaN); for complex arguments return x / abs(x), i.e. x normalized to the
// unit circle
template<class S, class Enable = typename md::enable_if_scalar_ntobj<S,void>::type>
typename md::promote_scalar<S>::type
                        sign(const S& x);

// the sign of x returned as integer (-1 for negative, 1 for positive, 
// 0 for plus or minus zero); not available for complex numbers
template<class S, class Enable = typename md::enable_if_scalar_ntobj<S,void>::type>
typename md::integer_or_object<S>::type
                        isign(const S& x);

// check whether the sign of x is negative; not available for complex numbers
template<class S, class Enable = typename md::enable_if_scalar_ntobj<S,void>::type>
typename md::bool_or_object<S>::type
                        signbit(const S& x);

// decompose given floating point value arg into a normalized fraction and an 
// integral power of two, i.e. x is represented as x = d x 2 ^ exp, where 
// 0.5 <= |d| < 1; if x is zero, then d = 0 and exp = 0; if is Inf or NaN, then
// d = x; not available for complex numbers
template<class S1, class Enable = typename md::enable_if_scalar_ntobj<S1,void>::type>
typename md::real_unify_types_promote<S1,Float>::type
                        frexp(const S1& x, Integer& exp);

// decompose given floating point value x into integral and fractional parts, 
// each having the same type and sign as x, the integral part is stored in 
// int_part argument and fractional part is returned; not available for complex
// numbers
template<class S1, class Enable = typename md::enable_if_scalar_ntobj<S1,void>::type,
        class Return = typename md::real_unify_types_promote<S1,Float>::type>
Return                  modf(const S1& x, Return& int_part);

//---------------------------------------------------------------------------
//                      COMBINATORICS
//---------------------------------------------------------------------------

// return the (2 * n)-th Bernoulli number B_{2n}, n >= 0
Real                    bernoulli_b2n(Integer n);
Float                   fbernoulli_b2n(Integer n);

// return the (2 * n)-th Bernoulli number B_{2n}, n >= 0; select single or double
// precision using the template argument
template<class S1, class Enable = typename md::enable_if_scalar_ntobj<S1,void>::type>
typename md::unify_types_promote<S1,Float>::type
                        bernoulli_b2n(Integer n);

// return largest value, such that bernoulli_b2n does not overflow
Integer                 max_bernoulli_b2n();
Integer                 fmax_bernoulli_b2n();

// return largest value, such that bernoulli_b2n<T> does not overflow; select
// single or double precision using the template argument
template<class S1, class Enable = typename md::enable_if_scalar_ntobj<S1,void>::type>
Integer                 max_bernoulli_b2n();

// return n-th prime numbers based on table lookup, where 0 <= n < N, and 
// N is given by the function prime_max_index()
size_t                  prime(Integer n);

// return the number of precomputed prime numbers
Integer                 prime_max_count();

// return the factorial i!, i.e. prod_{k=1}^i k
Real                    factorial(Integer i);
Float                   ffactorial(Integer i);

// return the factorial i!, i.e. prod_{k=1}^i k; select single or double precision
// using the template argument
template<class S1, class Enable = typename md::enable_if_scalar_ntobj<S1,void>::type>
typename md::unify_types_promote<S1,Float>::type
                        factorial(Integer i);

// return double factorial i!!:
//     i!! = i * (i-2) * ... * 1   for i odd
//     i!! = i * (i-2) * ... * 2   for i even
//     i!! = 1                     for i = 0, -1
Real                    double_factorial(Integer i);
Float                   fdouble_factorial(Integer i);

// return double factorial i!!; select single or double precision using the 
// template argument
template<class S1, class Enable = typename md::enable_if_scalar_ntobj<S1,void>::type>
typename md::unify_types_promote<S1,Float>::type
                        double_factorial(Integer i);

// return the rising factorial of x and i:
//     rising_factorial(x, i)  = x*(x+1)*...*(x+i-1)
// both x and i can be positive or negative; not defined for complex values
template<class S1, class Enable = typename md::enable_if_scalar_ntobj<S1,void>::type>
typename md::unify_types_promote<S1,Float>::type
                        rising_factorial(const S1& x, Integer i);

// return the falling factorial of x and i:
//     falling_factorial(x, i) = x*(x-1)*...*(x-i+1)
// x can be positive or negative, i must be positive
// not defined for complex values
template<class S1, class Enable = typename md::enable_if_scalar_ntobj<S1,void>::type>
typename md::unify_types_promote<S1,Float>::type
                        falling_factorial(const S1& x, Integer i);

// returns the binomial coefficient:
//     binomial_coefficient(n,k)   = n! / k! / (n-k)!
// requires k <= n, and k,n >= 0
Real                    binomial_coefficient(Integer n, Integer k);
Float                   fbinomial_coefficient(Integer n, Integer k);

// returns the binomial coefficient; select single or double precision using the 
// template argument
template<class S1, class Enable = typename md::enable_if_scalar_ntobj<S1,void>::type>
typename md::unify_types_promote<S1,Float>::type
                        binomial_coefficient(Integer n, Integer k);

};

#include "matcl-scalar/details/func_unary.inl"
#include "matcl-scalar/details/combinatorics.inl"