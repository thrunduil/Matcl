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

#include "matcl-scalar/config.h"
#include "matcl-dynamic/function_name.h"
#include "matcl-dynamic/predefined_functions_names.h"
#include "matcl-core/matrix/enums.h"
#include "matcl-core/matrix/scalar_types.h"

namespace matcl { namespace dynamic { namespace functions
{

//----------------------------------------------------------
//          already defined functions for any object
//----------------------------------------------------------

//these functions are defined for objects in different part
//of matcl; here we define functions names for dynamic calls

// disp(object); defined for any object type; function with
// this name must return Unit type
struct MATCL_SCALAR_EXPORT disp        { static function_name eval(); };

// to_string(object); defined for any object type; function with
// this name must return String
struct MATCL_SCALAR_EXPORT to_string   { static function_name eval(); };

}}};

namespace matcl { namespace dynamic
{

//--------------------------------------------------------------------
//               functions for object
//--------------------------------------------------------------------

// call the element-wise negation operator (i.e. ~x); for scalar types
// this function should be equivalent to op_not; there is a registered
// function bool op_neg(Any), that calls the op_not function
MATCL_SCALAR_EXPORT object     op_neg(const object& x);

// call the element-wise cast bool operator; for scalar types this 
// function should be equivalent to cast_bool; there is a registered
// function bool op_true(Any), that calls the cast_bool function
MATCL_SCALAR_EXPORT object     op_true(const object& x);

// argument (polar angle) of a complex number
MATCL_SCALAR_EXPORT object     arg(const object& x);

// conjugate transposition
MATCL_SCALAR_EXPORT object     conj(const object& x);

// absolute value
MATCL_SCALAR_EXPORT object     abs(const object& x);

// absolute value squared
MATCL_SCALAR_EXPORT object     abs2(const object& x);

// check if value is nan
MATCL_SCALAR_EXPORT object     is_nan(const object& x);

// check if value is finite
MATCL_SCALAR_EXPORT object     is_finite(const object& x);

// check if value is infinite
MATCL_SCALAR_EXPORT object     is_inf(const object& x);

// check if value is regular (i.e., neither NaN, nor an infinity nor zero); 
MATCL_SCALAR_EXPORT object     is_regular(const object& x);

// check if value is normal (i.e., neither NaN, nor an infinity nor 
// zero nor subnormal); 
MATCL_SCALAR_EXPORT object     is_normal(const object& x);

// check if value is an integer
MATCL_SCALAR_EXPORT object     is_int(const object& x);

// check if value is real (i.e. integer scalar, real floating point
// scalar or complex scalar with zero imaginary part)
MATCL_SCALAR_EXPORT object     is_real(const object& x);

// return value class of floating point number
MATCL_SCALAR_EXPORT fp_type    fpclassify(const object& x);

// scalar inverse invs(object); if no errors occure, then mul(x, invs(x)) = 1;
MATCL_SCALAR_EXPORT object     invs(const object& x);

// matrix inverse inv(object); if no errors occure, then mmul(x, inv(x)) = 1
MATCL_SCALAR_EXPORT object     inv(const object& x);

// element by element multiplication
MATCL_SCALAR_EXPORT object     elem_mul(const object& x, const object& y);

// four quadrant arctangent of x and y; -pi <= atan2(x, y) <= pi
MATCL_SCALAR_EXPORT object     atan2(const object& x, const object& y);

// calculate sqrt(|x|^2 + |y|^2) avoiding owerflow and underflow at
// intermediate stages of computation
MATCL_SCALAR_EXPORT object     hypot(const object& x, const object& y);

// return is x - n * y where n = floor(x / y) if y != 0 and 0 otherwise;
MATCL_SCALAR_EXPORT object     mod(const object& x, const object& y);

// return is x - n * y where n = trunc(x / y) if y != 0 and 0 otherwise
MATCL_SCALAR_EXPORT object     rem(const object& x, const object& y);

// copy sign from the second operand to the first operand
MATCL_SCALAR_EXPORT object     copysign(const object& x, const object& y);

// return the next representable value after x in the direction of y
MATCL_SCALAR_EXPORT object     nextafter(const object& x, const object& y);

// return the next representable value after x in the direction +INF
MATCL_SCALAR_EXPORT object     nextabove(const object& x);

// return the next representable value after x in the direction -INF
MATCL_SCALAR_EXPORT object     nextbelow(const object& x);

// check whether the sign of x is negative
MATCL_SCALAR_EXPORT object     signbit(const object& x);

// the sign of x (-1 for negative, 1 for positive, 0 for zero, 
// NaN if x is NaN); for complex arguments return x / abs(x)
MATCL_SCALAR_EXPORT object     sign(const object& x);

// the sign of x returned as integer (-1 for negative, 1 for positive, 
// 0 for plus or minus zero); for complex elements return sign of the
// real part
MATCL_SCALAR_EXPORT object     isign(const object& x);

// scale x by 2 raised to the power of n, returning ldexp(x,n) = x * 2^n
MATCL_SCALAR_EXPORT object     ldexp(const object& x, const object& n);

// function equivalent to ldexp
MATCL_SCALAR_EXPORT object     scalbn(const object& x, const object& n);

// decompose given floating point value arg into a normalized fraction
// and an integral power of two, i.e. x is represented as x = d x 2 ^ exp,
// where 0.5 <= |d| < 1; if x is zero, then d = 0 and exp = 0;
MATCL_SCALAR_EXPORT object     frexp(const object& x, object& exp);

// decompose given floating point value x into integral and fractional
// parts, each having the same type and sign as x, the integral part is 
// stored in the second argument and fractional part is returned;
MATCL_SCALAR_EXPORT object     modf(const object& x, object& n);

// return max(x - y, 0) or NaN if x or y is NaN
MATCL_SCALAR_EXPORT object     fdim(const object& x, const object& n);

// extract the value of the unbiased exponent from the floating-point
// argument x, and returns it as a floating-point value; equivalent to
// log2(|x|) rounded to integer number toward -INF
MATCL_SCALAR_EXPORT object     logb(const object& x);

// cast result of logb to integer
MATCL_SCALAR_EXPORT object     ilogb(const object& x);

// return fused multiply add (a*b + c)
MATCL_SCALAR_EXPORT object     fma(const object& x, const object& y, const object& c);

// return fused multiply subtract (a*b - c)
MATCL_SCALAR_EXPORT object     fms(const object& x, const object& y, const object& c);

// calculate a * b + c * d accurately
MATCL_SCALAR_EXPORT object     dot2_ac(const object& x, const object& y, const object& c, 
                                const object& d);

// return epsilon value eps, i.e positive distance from abs(x) to the next larger in
// magnitude floating point number of the same precision as x
MATCL_SCALAR_EXPORT object     eps(const object& x);

// return the square root of x
MATCL_SCALAR_EXPORT object     sqrt(const object& x);

// compute the complex conversion version of the sqrt function
MATCL_SCALAR_EXPORT object     sqrt_c(const object& x);

// return sqrt(1 + x) - 1
MATCL_SCALAR_EXPORT object     sqrt1pm1(const object& x);

// compute complex conversion of the sqrt function
MATCL_SCALAR_EXPORT object     sqrt1pm1_c(const object& x);

// return the cubic root of x, for negative values result is also negative
MATCL_SCALAR_EXPORT object     cbrt(const object& x);

// return the exponential function, i.e. the e (Euler's number, 2.7182818) 
// raised to the given power x; 
MATCL_SCALAR_EXPORT object     exp(const object& x);

// return exp(x) - 1
MATCL_SCALAR_EXPORT object     expm1(const object& x);

// return exp(x * 1i), where 1i is an imaginary unit
MATCL_SCALAR_EXPORT object     expi(const object& x);

// return 2^x, i.e two raised to the given power x
MATCL_SCALAR_EXPORT object     exp2(const object& x);

// return 10^x, i.e ten raised to the given power x
MATCL_SCALAR_EXPORT object     exp10(const object& x);

// compute the natural logarithm of x;
MATCL_SCALAR_EXPORT object     log(const object& x);

// compute the complex conversion version of the natural logarithm of x
MATCL_SCALAR_EXPORT object     log_c(const object& x);

// compute log(1 + x)
MATCL_SCALAR_EXPORT object     log1p(const object& x);

// compute the complex conversion version of the log1p function
MATCL_SCALAR_EXPORT object     log1p_c(const object& x);

// compute the base 2 logarithm of x;
MATCL_SCALAR_EXPORT object     log2(const object& x);

// compute the complex conversion version of the base 2 logarithm of x
MATCL_SCALAR_EXPORT object     log2_c(const object& x);

// compute the base 10 logarithm of x;
MATCL_SCALAR_EXPORT object     log10(const object& x);

// compute the complex conversion version of the base 10 logarithm of x
MATCL_SCALAR_EXPORT object     log10_c(const object& x);

// return the sine function of x in radians.
MATCL_SCALAR_EXPORT object     sin(const object& x);

// return the cosine function of x in radians.
MATCL_SCALAR_EXPORT object     cos(const object& x);

// return the tangent function of x in radians.
MATCL_SCALAR_EXPORT object     tan(const object& x);

// return the cotangent function of x in radians.
MATCL_SCALAR_EXPORT object     cot(const object& x);

// return the secant function of x in radians.
MATCL_SCALAR_EXPORT object     sec(const object& x);

// return the cosecant function of x in radians.
MATCL_SCALAR_EXPORT object     csc(const object& x);

// return the hyperbolic sine function of x in radians.
MATCL_SCALAR_EXPORT object     sinh(const object& x);

// return the hyperbolic cosine function of x in radians.
MATCL_SCALAR_EXPORT object     cosh(const object& x);

// return the  hyperbolic tangent function of x in radians.
MATCL_SCALAR_EXPORT object     tanh(const object& x);

// return the  hyperbolic cotangent function of x in radians.
MATCL_SCALAR_EXPORT object     coth(const object& x);

// return the hyperbolic secant function of x in radians.
MATCL_SCALAR_EXPORT object     sech(const object& x);

// return the hyperbolic cosecant function of x in radians.
MATCL_SCALAR_EXPORT object     csch(const object& x);

// return the inverse sine function
MATCL_SCALAR_EXPORT object     asin(const object& x);

// compute complex conversion of the asin function
MATCL_SCALAR_EXPORT object     asin_c(const object& x);

// return the inverse cosine function
MATCL_SCALAR_EXPORT object     acos(const object& x);

// compute complex conversion of the acos function
MATCL_SCALAR_EXPORT object     acos_c(const object& x);

// return the inverse tangent function
MATCL_SCALAR_EXPORT object     atan(const object& x);

// return the inverse cotangent function
MATCL_SCALAR_EXPORT object     acot(const object& x);

// return the inverse secant function
MATCL_SCALAR_EXPORT object     asec(const object& x);

// compute complex conversion of the asec function
MATCL_SCALAR_EXPORT object     asec_c(const object& x);

// return the inverse cosecant function
MATCL_SCALAR_EXPORT object     acsc(const object& x);

// compute complex conversion of the acsc function
MATCL_SCALAR_EXPORT object     acsc_c(const object& x);

// return the inverse hyperbolic sine function
MATCL_SCALAR_EXPORT object     asinh(const object& x);

// return the inverse hyperbolic cosine function
MATCL_SCALAR_EXPORT object     acosh(const object& x);

// compute complex conversion of the acosh function
MATCL_SCALAR_EXPORT object     acosh_c(const object& x);

// return the inverse hyperbolic tangent function
MATCL_SCALAR_EXPORT object     atanh(const object& x);

// compute complex conversion of the atanh function
MATCL_SCALAR_EXPORT object     atanh_c(const object& x);

// return the inverse hyperbolic cotangent function
MATCL_SCALAR_EXPORT object     acoth(const object& x);

// compute complex conversion of the acoth function
MATCL_SCALAR_EXPORT object     acoth_c(const object& x);

// return the inverse hyperbolic secant function
MATCL_SCALAR_EXPORT object     asech(const object& x);

// compute complex conversion of the asech function
MATCL_SCALAR_EXPORT object     asech_c(const object& x);

// return the inverse hyperbolic cosecant function
MATCL_SCALAR_EXPORT object     acsch(const object& x);

// round to nearest integer towards -INF
MATCL_SCALAR_EXPORT object     floor(const object& x);

// round to nearest integer towards +INF
MATCL_SCALAR_EXPORT object     ceil(const object& x);

// round to nearest integer, rounding halfway cases away from zero; 
MATCL_SCALAR_EXPORT object     round(const object& x);

// round to nearest integer towards zero
MATCL_SCALAR_EXPORT object     trunc(const object& x);

// round to nearest integer towards -INF and cast to integer
MATCL_SCALAR_EXPORT object     ifloor(const object& x);

// round to nearest integer towards +INF and cast to integer
MATCL_SCALAR_EXPORT object     iceil(const object& x);

// round to nearest integer, rounding halfway cases away from zero
// and cast to integer
MATCL_SCALAR_EXPORT object     iround(const object& x);

// round to nearest integer towards zero and cast to integer
MATCL_SCALAR_EXPORT object     itrunc(const object& x);

// equivalent to a/b
MATCL_SCALAR_EXPORT object     div(const object& x, const object& y);

// equivalent to a/b but 0/0 gives 0
MATCL_SCALAR_EXPORT object     div_0(const object& x, const object& y);

// equivalent to a/b but 0/0 gives 1
MATCL_SCALAR_EXPORT object     div_1(const object& x, const object& y);

// compute power, i.e. x^y
MATCL_SCALAR_EXPORT object     pow(const object& x, const object& y);

// compute complex conversion of the pow function
MATCL_SCALAR_EXPORT object     pow_c(const object& x, const object& y);

// return x^y - 1
MATCL_SCALAR_EXPORT object     powm1(const object& x, const object& y);

// return x && y
MATCL_SCALAR_EXPORT bool       op_and(const object& x, const object& y);

// return x || y
MATCL_SCALAR_EXPORT bool       op_or(const object& x, const object& y);

// return x ^ y
MATCL_SCALAR_EXPORT bool       op_xor(const object& x, const object& y);

// return x & y
MATCL_SCALAR_EXPORT object     elem_and(const object& x, const object& y);

// return x | y
MATCL_SCALAR_EXPORT object     elem_or(const object& x, const object& y);

// return x ^ y
MATCL_SCALAR_EXPORT object     elem_xor(const object& x, const object& y);

// return max(x, y)
MATCL_SCALAR_EXPORT object     max(const object& x, const object& y);

// return min(x, y)
MATCL_SCALAR_EXPORT object     min(const object& x, const object& y);

// return x == y assuming that nan values are equal
MATCL_SCALAR_EXPORT object     eeq_nan(const object& x, const object& y);

// return x != y assuming that nan values are equal
MATCL_SCALAR_EXPORT object     neq_nan(const object& x, const object& y);

// return the rising factorial of x and integer i:
// rising_factorial(x, i)  = x*(x+1)*...*(x+i-1)
MATCL_SCALAR_EXPORT object     rising_factorial(const object& x, Integer i);

// return the falling factorial of x and integer i:
// falling_factorial(x, i) = x*(x-1)*...*(x-i+1)
MATCL_SCALAR_EXPORT object     falling_factorial(const object& x, Integer i);

//--------------------------------------------------------------------
//               template functions for object
//--------------------------------------------------------------------
// return the (2 * n)-th Bernoulli number B_{2n}, n >= 0
MATCL_SCALAR_EXPORT object     bernoulli_b2n(Type t, Integer n);

// return largest value, such that bernoulli_b2n does not overflow
MATCL_SCALAR_EXPORT Integer    max_bernoulli_b2n(Type t);

// return the factorial i!, i.e. prod_{k=1}^i k
MATCL_SCALAR_EXPORT object     factorial(Type t, Integer i);

// return double factorial i!!:
//     i!! = i * (i-2) * ... * 1   for i odd
//     i!! = i * (i-2) * ... * 2   for i even
//     i!! = 1                     for i = 0, -1
MATCL_SCALAR_EXPORT object     double_factorial(Type t, Integer i);

// returns the binomial coefficient:
//     binomial_coefficient(n,k)   = n! / k! / (n-k)!
// requires k <= n, and k,n >= 0
MATCL_SCALAR_EXPORT object     binomial_coefficient(Type t, Integer n, Integer k);

};};
    