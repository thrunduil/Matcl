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

#include "matcl-dynamic/function_name.h"

namespace matcl { namespace dynamic { namespace functions
{

// function names used in matcl; these functions are not
// registered in matcl-dynamic

// call the element-wise negation operator (i.e. ~x); for scalar types
// this function should be equivalent to op_not
struct MATCL_DYN_EXPORT op_neg          { static const function_name& eval(); };

// call the element-wise cast bool operator; for scalar types this 
// function should be equivalent to cast_bool
struct MATCL_DYN_EXPORT op_true         { static const function_name& eval(); };

// scalar inverse invs(object); if no errors occure, then mul(x, invs(x)) = 1
struct MATCL_DYN_EXPORT invs            { static const function_name& eval(); };

// matrix inverse inv(object); if no errors occure, then mmul(x, inv(x)) = 1
struct MATCL_DYN_EXPORT inv             { static const function_name& eval(); };

// the power function x^y
struct MATCL_DYN_EXPORT pow             { static const function_name& eval(); };

// compute complex conversion of the pow function
struct MATCL_DYN_EXPORT pow_c           { static const function_name& eval(); };

// return x^y - 1
struct MATCL_DYN_EXPORT powm1           { static const function_name& eval(); };

// argument (polar angle) of a complex number
struct MATCL_DYN_EXPORT arg             { static const function_name& eval(); };    

// conjugate transposition
struct MATCL_DYN_EXPORT conj            { static const function_name& eval(); };

// absolute value
struct MATCL_DYN_EXPORT abs             { static const function_name& eval(); };

// absolute value squared
struct MATCL_DYN_EXPORT abs2            { static const function_name& eval(); };

// check if value is nan
struct MATCL_DYN_EXPORT is_nan          { static const function_name& eval(); };

// check if value is finite
struct MATCL_DYN_EXPORT is_finite       { static const function_name& eval(); };

// check if value is infinite
struct MATCL_DYN_EXPORT is_inf          { static const function_name& eval(); };

// check if value is regular (i.e., neither NaN, nor an infinity nor zero); 
struct MATCL_DYN_EXPORT is_regular      { static const function_name& eval(); };

// check if value is regular (i.e., neither NaN, nor an infinity nor zero
struct MATCL_DYN_EXPORT is_normal       { static const function_name& eval(); };

// check if value is integer
struct MATCL_DYN_EXPORT is_int          { static const function_name& eval(); };

// check if value is real (i.e. integer scalar, real floating point
// scalar or complex scalar with zero imaginary part)
struct MATCL_DYN_EXPORT is_real         { static const function_name& eval(); };

// return value class of floating point number, must return OInteger
struct MATCL_DYN_EXPORT fpclassify      { static const function_name& eval(); };    

// return fused multiply add (a*b + c)
struct MATCL_DYN_EXPORT fma             { static const function_name& eval(); };

// return fused multiply subtract (a*b - c)
struct MATCL_DYN_EXPORT fms             { static const function_name& eval(); };

// form a*b + c*d accurately
struct MATCL_DYN_EXPORT dot2_ac         { static const function_name& eval(); };

// element by element multiplication
struct MATCL_DYN_EXPORT elem_mul        { static const function_name& eval(); }; 

// division x / y; integer division gives Real result; 0/0 division gives 0
struct MATCL_DYN_EXPORT div_0           { static const function_name& eval(); };

// division x / y; integer division gives Real result; 0/0 division gives 1
struct MATCL_DYN_EXPORT div_1           { static const function_name& eval(); };

// name of atan2 binary function
struct MATCL_DYN_EXPORT atan2           { static const function_name& eval(); }; 

// name of hypot binary function
struct MATCL_DYN_EXPORT hypot           { static const function_name& eval(); }; 

// return is x - n * y where n = floor(x / y) if y != 0 and 0 otherwise;
struct MATCL_DYN_EXPORT mod             { static const function_name& eval(); }; 

// return is x - n * y where n = trunc(x / y) if y != 0 and 0 otherwise
struct MATCL_DYN_EXPORT rem             { static const function_name& eval(); }; 

// copy sign from the second operand to the first operand
struct MATCL_DYN_EXPORT copysign        { static const function_name& eval(); };

// return the next representable value after x in the direction of y
struct MATCL_DYN_EXPORT nextafter       { static const function_name& eval(); };

// return the next representable value after x in the direction of -INF
struct MATCL_DYN_EXPORT nextbelow       { static const function_name& eval(); };

// return the next representable value after x in the direction of +INF
struct MATCL_DYN_EXPORT nextabove       { static const function_name& eval(); };

// return number of distinct representations between x and y;
// if one of arguments is NaN, then NaN is returned;
struct MATCL_DYN_EXPORT float_distance  { static const function_name& eval(); };

// determines if the given value is negative
struct MATCL_DYN_EXPORT signbit         { static const function_name& eval(); };    

// the sign of x (-1 for negative, 1 for positive, 0 for zero, 
// NaN if x is NaN); for complex arguments return x / abs(x)
struct MATCL_DYN_EXPORT sign            { static const function_name& eval(); };

// the sign of x returned as integer (-1 for negative, 1 for positive, 
// 0 for plus or minus zero); for complex elements return sign of the
// real part
struct MATCL_DYN_EXPORT isign           { static const function_name& eval(); };

// scale x by 2 raised to the power of n, returning ldexp(x,n) = x * 2^n
struct MATCL_DYN_EXPORT ldexp           { static const function_name& eval(); };    

// decompose given floating point value arg into a normalized fraction
// and an integral power of two, i.e. x is represented as x = d x 2 ^ exp;
// for functions with this name the second argument must have a reference type
struct MATCL_DYN_EXPORT frexp           { static const function_name& eval(); };    

// decompose given floating point value x into integral and fractional
// parts, each having the same type and sign as x, the integral part is 
// stored in the second argument and fractional part is returned;
// for function with this name the second argument must have a reference type
struct MATCL_DYN_EXPORT modf            { static const function_name& eval(); };    

// return max(x - y, 0) or NaN if x or y is NaN
struct MATCL_DYN_EXPORT fdim            { static const function_name& eval(); };    

// extract the value of the unbiased exponent from the floating-point
// argument x, and returns it as a floating-point value; equivalent to
// log2(|x|) rounded to integer number toward -INF
struct MATCL_DYN_EXPORT logb            { static const function_name& eval(); };    

// cast result of logb to integer
struct MATCL_DYN_EXPORT ilogb           { static const function_name& eval(); };

// return epsilon value eps, i.e positive distance from abs(x) to the
// next larger in magnitude floating point number of the same precision as x
struct MATCL_DYN_EXPORT eps             { static const function_name& eval(); };    

// return the square root of x
struct MATCL_DYN_EXPORT sqrt            { static const function_name& eval(); };    

// compute complex conversion of the sqrt function
struct MATCL_DYN_EXPORT sqrt_c          { static const function_name& eval(); };    

// return sqrt(1 + x) - 1
struct MATCL_DYN_EXPORT sqrt1pm1        { static const function_name& eval(); };

// compute complex conversion of the sqrt function
struct MATCL_DYN_EXPORT sqrt1pm1_c      { static const function_name& eval(); };

// return the cubic root of x, for negative values result is also negative
struct MATCL_DYN_EXPORT cbrt            { static const function_name& eval(); };    

// return the exponential function, i.e. the e (Euler's number, 2.7182818) 
// raised to the given power x; 
struct MATCL_DYN_EXPORT exp             { static const function_name& eval(); };    

// return exp(x) - 1
struct MATCL_DYN_EXPORT expm1           { static const function_name& eval(); };    

// return exp(x * 1i), where 1i is an imaginary unit
struct MATCL_DYN_EXPORT expi            { static const function_name& eval(); };    

// return 2^x, i.e two raised to the given power x
struct MATCL_DYN_EXPORT exp2            { static const function_name& eval(); };    

// return 10^x, i.e ten raised to the given power x
struct MATCL_DYN_EXPORT exp10           { static const function_name& eval(); };    

// return the natural logarithm fuction
struct MATCL_DYN_EXPORT log             { static const function_name& eval(); };    

// compute complex conversion of the natural logarithm fuction
struct MATCL_DYN_EXPORT log_c           { static const function_name& eval(); };    

// return log(1+x)
struct MATCL_DYN_EXPORT log1p           { static const function_name& eval(); };    

// compute complex conversion of the log1p fuction
struct MATCL_DYN_EXPORT log1p_c         { static const function_name& eval(); };    

// return the base 2 logarithm fuction
struct MATCL_DYN_EXPORT log2            { static const function_name& eval(); };    

// compute complex conversion of the base 2 logarithm fuction
struct MATCL_DYN_EXPORT log2_c          { static const function_name& eval(); };    

// return the base 10 logarithm fuction
struct MATCL_DYN_EXPORT log10           { static const function_name& eval(); };    

// compute complex conversion of the base 10 logarithm fuction
struct MATCL_DYN_EXPORT log10_c         { static const function_name& eval(); };    

// return the sine function of x in radians.
struct MATCL_DYN_EXPORT sin             { static const function_name& eval(); };

// return the cosine function of x in radians.
struct MATCL_DYN_EXPORT cos             { static const function_name& eval(); };

// return the tangent function of x in radians.
struct MATCL_DYN_EXPORT tan             { static const function_name& eval(); };

// return the cotangent function of x in radians.
struct MATCL_DYN_EXPORT cot             { static const function_name& eval(); };

// return the secant function of x in radians.
struct MATCL_DYN_EXPORT sec             { static const function_name& eval(); };

// return the cosecant function of x in radians.
struct MATCL_DYN_EXPORT csc             { static const function_name& eval(); };

// return the hyperbolic sine function of x in radians.
struct MATCL_DYN_EXPORT sinh            { static const function_name& eval(); };

// return the hyperbolic cosine function of x in radians.
struct MATCL_DYN_EXPORT cosh            { static const function_name& eval(); };

// return the hyperbolic tangent function of x in radians.
struct MATCL_DYN_EXPORT tanh            { static const function_name& eval(); };

// return the hyperbolic cotangent function of x in radians.
struct MATCL_DYN_EXPORT coth            { static const function_name& eval(); };

// return the hyperbolic secant function of x in radians.
struct MATCL_DYN_EXPORT sech            { static const function_name& eval(); };

// return the hyperbolic cosecant function of x in radians.
struct MATCL_DYN_EXPORT csch            { static const function_name& eval(); };

// return the inverse sine function
struct MATCL_DYN_EXPORT asin            { static const function_name& eval(); };

// compute complex conversion of the asin function
struct MATCL_DYN_EXPORT asin_c          { static const function_name& eval(); };

// return the inverse cosine function
struct MATCL_DYN_EXPORT acos            { static const function_name& eval(); };

// compute complex conversion of the acos function
struct MATCL_DYN_EXPORT acos_c          { static const function_name& eval(); };

// return the inverse tangent function
struct MATCL_DYN_EXPORT atan            { static const function_name& eval(); };

// return the inverse cotangent function
struct MATCL_DYN_EXPORT acot            { static const function_name& eval(); };

// return the inverse secant function
struct MATCL_DYN_EXPORT asec            { static const function_name& eval(); };

// compute complex conversion of the asec function
struct MATCL_DYN_EXPORT asec_c          { static const function_name& eval(); };

// return the inverse cosecant function
struct MATCL_DYN_EXPORT acsc            { static const function_name& eval(); };

// compute complex conversion of the acsc function
struct MATCL_DYN_EXPORT acsc_c          { static const function_name& eval(); };

// return the inverse hyperbolic sine function
struct MATCL_DYN_EXPORT asinh           { static const function_name& eval(); };

// return the inverse hyperbolic cosine function
struct MATCL_DYN_EXPORT acosh           { static const function_name& eval(); };

// compute complex conversion of the acosh function
struct MATCL_DYN_EXPORT acosh_c         { static const function_name& eval(); };

// return the inverse hyperbolic tangent function
struct MATCL_DYN_EXPORT atanh           { static const function_name& eval(); };

// return the inverse hyperbolic cotangent function
struct MATCL_DYN_EXPORT acoth           { static const function_name& eval(); };

// compute complex conversion of the acoth function
struct MATCL_DYN_EXPORT acoth_c         { static const function_name& eval(); };

// return the inverse hyperbolic secant function
struct MATCL_DYN_EXPORT asech           { static const function_name& eval(); };

// compute complex conversion of the asech function
struct MATCL_DYN_EXPORT asech_c         { static const function_name& eval(); };

// return the inverse hyperbolic cosecant function
struct MATCL_DYN_EXPORT acsch           { static const function_name& eval(); };

// compute complex conversion of the atanh function
struct MATCL_DYN_EXPORT atanh_c         { static const function_name& eval(); };

// round to nearest integer towards -INF
struct MATCL_DYN_EXPORT floor           { static const function_name& eval(); };

// round to nearest integer towards +INF
struct MATCL_DYN_EXPORT ceil            { static const function_name& eval(); };

// round to nearest integer towards zero
struct MATCL_DYN_EXPORT trunc           { static const function_name& eval(); };

// round to nearest integer, rounding halfway cases away from zero; 
struct MATCL_DYN_EXPORT round           { static const function_name& eval(); };

// round to nearest integer towards -INF and cast to integer
struct MATCL_DYN_EXPORT ifloor          { static const function_name& eval(); };

// round to nearest integer towards +INF and cast to integer
struct MATCL_DYN_EXPORT iceil           { static const function_name& eval(); };

// round to nearest integer towards zero and cast to integer
struct MATCL_DYN_EXPORT itrunc          { static const function_name& eval(); };

// round to nearest integer, rounding halfway cases away from zero
// and cast to integer
struct MATCL_DYN_EXPORT iround          { static const function_name& eval(); };

// return x && y; functions with this name must return OBool
struct MATCL_DYN_EXPORT op_and          { static const function_name& eval(); };

// return x || y; functions with this name must return OBool
struct MATCL_DYN_EXPORT op_or           { static const function_name& eval(); };

// return x ^ y; functions with this name must return OBool
struct MATCL_DYN_EXPORT op_xor          { static const function_name& eval(); };

// return x & y
struct MATCL_DYN_EXPORT elem_and        { static const function_name& eval(); };

// return x | y
struct MATCL_DYN_EXPORT elem_or         { static const function_name& eval(); };

// return x ^ y
struct MATCL_DYN_EXPORT elem_xor        { static const function_name& eval(); };

// return x == y assuming that nan values are equal
struct MATCL_DYN_EXPORT eeq_nan         { static const function_name& eval(); };

// return x != y assuming that nan values are equal
struct MATCL_DYN_EXPORT neq_nan         { static const function_name& eval(); };

// return max(x, y)
struct MATCL_DYN_EXPORT max             { static const function_name& eval(); };

// return min(x, y)
struct MATCL_DYN_EXPORT min             { static const function_name& eval(); };

// return the rising factorial of x and integer i:
// rising_factorial(x, i)  = x*(x+1)*...*(x+i-1)
struct MATCL_DYN_EXPORT rising_factorial    { static const function_name& eval(); };

// return the falling factorial of x and integer i:
// falling_factorial(x, i) = x*(x-1)*...*(x-i+1)
struct MATCL_DYN_EXPORT falling_factorial   { static const function_name& eval(); };

// return the (2 * n)-th Bernoulli number B_{2n}, n >= 0
struct MATCL_DYN_EXPORT bernoulli_b2n       { static const function_name& eval(); };

// return largest value, such that bernoulli_b2n does not overflow;
// function with this name must return OInteger
struct MATCL_DYN_EXPORT max_bernoulli_b2n   { static const function_name& eval(); };

// return the factorial i!, i.e. prod_{k=1}^i k
struct MATCL_DYN_EXPORT factorial           { static const function_name& eval(); };

// return double factorial i!!
struct MATCL_DYN_EXPORT double_factorial    { static const function_name& eval(); };

// returns the binomial coefficient:
//     binomial_coefficient(n,k)   = n! / k! / (n-k)!
// requires k <= n, and k,n >= 0
struct MATCL_DYN_EXPORT binomial_coefficient{ static const function_name& eval(); };

};};};
