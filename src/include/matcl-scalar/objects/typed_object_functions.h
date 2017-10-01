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

#include "matcl-core/matrix/enums.h"
#include "matcl-dynamic/result_of.h"

#include "matcl-dynamic/function_name.h"
#include "matcl-dynamic/predefined_functions_names.h"

#include "matcl-scalar/object.h"
#include "matcl-dynamic/typed_object_functions.h"

namespace matcl { namespace dynamic
{

//--------------------------------------------------------------------
//               predefined functions for object_type
//--------------------------------------------------------------------
// a unary function func is enabled for object_type<T> if a function func(a)
// can be found by unqualified call in the namespace matcl for arguments of 
// type T unless otherwise stated
//
// similarly a binary function func is enabled for arguments of type 
// func(object_type<T>, object_type<S>), func(object_type<T>, S), func(T, object_type<S>)
// if a function func(a,b) can be found by unqualified call in the namespace 
// matcl for arguments of type T and S unless otherwise stated
//
// notice, that some functions can be automatically defined in the namespace
// matcl, see func_forwarding.h

//---------------------------------------------------------------------------
//                      VALUE CLASSIFICATION
//---------------------------------------------------------------------------
// return value class of floating point number
template<class T, class Enable = typename result_of::result_of_fpclassify<T>::type_object>
fp_type                 fpclassify(const object_type<T>& x);

// check if value is nan
template<class T, class Ret = typename result_of::result_of_is_nan<T>::type_object>
Ret                     is_nan(const object_type<T>& x);

// check if value is finite
template<class T, class Ret = typename result_of::result_of_is_finite<T>::type_object>
Ret                     is_finite(const object_type<T>& x);

// check if value is infinite
template<class T, class Ret = typename result_of::result_of_is_inf<T>::type_object>
Ret                     is_inf(const object_type<T>& x);

// check if value is regular (i.e., neither NaN, nor an infinity nor zero); 
template<class T, class Ret = typename result_of::result_of_is_regular<T>::type_object>
Ret                     is_regular(const object_type<T>& x);

// check if value is regular (i.e., neither NaN, nor an infinity nor zero
// nor subnormal)
template<class T, class Ret = typename result_of::result_of_is_normal<T>::type_object>
Ret                     is_normal(const object_type<T>& x);

// check if value is integer
template<class T, class Ret = typename result_of::result_of_is_int<T>::type_object>
Ret                     is_int(const object_type<T>& x);

// check if value is real (i.e. integer scalar, real floating point
template<class T, class Ret = typename result_of::result_of_is_real<T>::type_object>
Ret                     is_real(const object_type<T>& x);

//is_zero and is_one are defined in matcl-dynamic

//---------------------------------------------------------------------------
//                      LOGICAL OPERATIONS
//---------------------------------------------------------------------------

// logical negation; equivalent to operator!
template<class T, class Res = typename result_of::result_of_neg<T>::type_object>
Res                     neg(const object_type<T>& x);

// logical negation; equivalent to operator!
template<class T, class Res = typename result_of::result_of_neg<T>::type_object>
Res                     is_false(const object_type<T>& x)   { return neg(x); };

// cast to boolean value; equivalent to cast_bool
template<class T, class Res = typename result_of::result_of_is_true<T>::type_object>
Res                     is_true(const object_type<T>& x);

//operator! and cast_bool is defined in matcl-dynamic

//---------------------------------------------------------------------------
//                      ARITHMETIC FUNCTIONS
//---------------------------------------------------------------------------

// unary minus; equivalent to operator-
template<class T>
typename result_of::result_of_uminus<T>::type_object
                        uminus(const object_type<T>& x)     { return -x; };

// unary plus
// this function is always enabled; return value is always equal to
// the argument
template<class T>
object_type<T>          operator+(const object_type<T>& x)  { return x; };

// matrix inverse inv(object); if no errors occure, then mmul(x, inv(x)) = 1
template<class T>
typename result_of::result_of_inv<T>::type_object
                        inv(const object_type<T>& x);

// scalar inverse invs(object); if no errors occure, then mul(x, invs(x)) = 1
template<class T>
typename result_of::result_of_invs<T>::type_object
                        invs(const object_type<T>& x);

// absolute value
template<class T>
typename result_of::result_of_abs<T>::type_object
                        abs(const object_type<T>& x);

// absolute value squared
template<class T>
typename result_of::result_of_abs2<T>::type_object
                        abs2(const object_type<T>& x);

// argument (polar angle) of a complex number
template<class T>
typename result_of::result_of_arg<T>::type_object
                        arg(const object_type<T>& x);

// argument (polar angle) of a complex number; equivalent to arg
template<class T>
typename result_of::result_of_arg<T>::type_object
                        angle(const object_type<T>& x);

// conjugate transposition
template<class T>
typename result_of::result_of_conj<T>::type_object
                        conj(const object_type<T>& x);

// operator-, real, imag are defined in matcl-dynamic

//---------------------------------------------------------------------------
//               EXPONENTIAL, LOGARITHMIC AND POWER FUNCTIONS
//---------------------------------------------------------------------------

// return the square root function of x
template<class T>
typename result_of::result_of_sqrt<T>::type_object
                        sqrt(const object_type<T>& x);

// compute the complex conversion version of the square root fuction
template<class T>
typename result_of::result_of_sqrt_c<T>::type_object
                        sqrt_c(const object_type<T>& x);

// return sqrt(1+x) - 1
template<class T>
typename result_of::result_of_sqrt1pm1<T>::type_object
                        sqrt1pm1(const object_type<T>& x);

// return a complex conversion of the sqrt1pm1 function
template<class T>
typename result_of::result_of_sqrt1pm1_c<T>::type_object
                        sqrt1pm1_c(const object_type<T>& x);

// return the cubic root of x, for negative values result is also negative
template<class T>
typename result_of::result_of_cbrt<T>::type_object
                        cbrt(const object_type<T>& x);

// return the exponential function, i.e. the e (Euler's number, 2.7182818) 
// raised to the given power x; 
template<class T>
typename result_of::result_of_exp<T>::type_object
                        exp(const object_type<T>& x);

// return 2^x, i.e two raised to the given power x
template<class T>
typename result_of::result_of_exp2<T>::type_object
                        exp2(const object_type<T>& x);

// return 2^x, equivalent to exp2; 
template<class T>
typename result_of::result_of_exp2<T>::type_object
                        pow2(const object_type<T>& x)   { return exp2(x); };

// return 10^x, i.e ten raised to the given power x
template<class T>
typename result_of::result_of_exp10<T>::type_object
                        exp10(const object_type<T>& x);

// return 10^x; equivalent to exp10
template<class T>
typename result_of::result_of_exp10<T>::type_object
                        pow10(const object_type<T>& x)  { return exp10(x); };

// return exp(x) - 1
template<class T>
typename result_of::result_of_expm1<T>::type_object
                        expm1(const object_type<T>& x);

// return exp(x * 1i), where 1i is an imaginary unit
template<class T>
typename result_of::result_of_expi<T>::type_object
                        expi(const object_type<T>& x);

// scale x by 2 raised to the power of n, returning ldexp(x,n) = x * 2^n
template<class T, class Ret = typename result_of::result_of_ldexp<T,Integer>::type_object>
Ret                     ldexp(const object_type<T>& x, Integer exp);

template<class T, class Ret = typename result_of::result_of_ldexp<T,Integer>::type_object>
Ret                     ldexp(const object_type<T>& x, const object_type<Integer>& exp);

// scale x by 2 raised to the power of n, returning ldexp(x,n) = x * 2^n
template<class T, class Ret = typename result_of::result_of_scalbn<T,Integer>::type_object>
Ret                     scalbn(const object_type<T>& x, Integer exp);

template<class T, class Ret = typename result_of::result_of_scalbn<T,Integer>::type_object>
Ret                     scalbn(const object_type<T>& x, const object_type<Integer>& exp);

// compute the natural logarithm of x
template<class T>
typename result_of::result_of_log<T>::type_object
                        log(const object_type<T>& x);

// compute the complex conversion version of the natural logarithm fuction
template<class T>
typename result_of::result_of_log_c<T>::type_object
                        log_c(const object_type<T>& x);

// compute the natural logarithm of x + 1
template<class T>
typename result_of::result_of_log1p<T>::type_object
                        log1p(const object_type<T>& x);

// compute the complex conversion version of the log1p fuction
template<class T>
typename result_of::result_of_log1p_c<T>::type_object
                        log1p_c(const object_type<T>& x);

// compute the base 2 logarithm of x
template<class T>
typename result_of::result_of_log2<T>::type_object
                        log2(const object_type<T>& x);

// compute the complex conversion version of the base 2 logarithm fuction
template<class T>
typename result_of::result_of_log2_c<T>::type_object
                        log2_c(const object_type<T>& x);

// compute the base 10 logarithm of x
template<class T>
typename result_of::result_of_log10<T>::type_object
                        log10(const object_type<T>& x);

// compute the complex conversion version of the base 10 logarithm fuction
template<class T>
typename result_of::result_of_log10_c<T>::type_object
                        log10_c(const object_type<T>& x);

// extract the value of the unbiased exponent from the floating-point
// argument x, and returns it as a floating-point value; equivalent to
// log2(|x|) rounded to integer number toward -INF
template<class T, class Ret = typename result_of::result_of_logb<T>::type_object>
Ret                     logb(const object_type<T>& x);

// cast result of logb to integer
template<class T, class Ret = typename result_of::result_of_ilogb<T>::type_object>
Ret                     ilogb(const object_type<T>& x);

//---------------------------------------------------------------------------
//                      TRIGONOMETRIC FUNCTIONS
//---------------------------------------------------------------------------
// return the sine function of x in radians
template<class T>
typename result_of::result_of_sin<T>::type_object
                        sin(const object_type<T>& x);

// return the cosine function of x in radians
template<class T>
typename result_of::result_of_cos<T>::type_object
                        cos(const object_type<T>& x);

// return the tangent function of x in radians
template<class T>
typename result_of::result_of_tan<T>::type_object
                        tan(const object_type<T>& x);

// return the cotangent function of x in radians
template<class T>
typename result_of::result_of_cot<T>::type_object
                        cot(const object_type<T>& x);

// return the secant function of x in radians
template<class T>
typename result_of::result_of_sec<T>::type_object
                        sec(const object_type<T>& x);

// return the cosecant function of x in radians
template<class T>
typename result_of::result_of_csc<T>::type_object
                        csc(const object_type<T>& x);

// return the hyperbolic sine function of x in radians
template<class T>
typename result_of::result_of_sinh<T>::type_object
                        sinh(const object_type<T>& x);

// return the hyperbolic cosine function of x in radians
template<class T>
typename result_of::result_of_cosh<T>::type_object
                        cosh(const object_type<T>& x);

// return the hyperbolic tangent function of x in radians
template<class T>
typename result_of::result_of_tanh<T>::type_object
                        tanh(const object_type<T>& x);

// return the hyperbolic cotangent function of x in radians
template<class T>
typename result_of::result_of_coth<T>::type_object
                        coth(const object_type<T>& x);

// return the hyperbolic secant function of x in radians
template<class T>
typename result_of::result_of_sech<T>::type_object
                        sech(const object_type<T>& x);

// return the hyperbolic cosecant function of x in radians
template<class T>
typename result_of::result_of_csch<T>::type_object
                        csch(const object_type<T>& x);

// return the inverse sine function
template<class T>
typename result_of::result_of_asin<T>::type_object
                        asin(const object_type<T>& x);

// compute complex conversion of the asin function
template<class T>
typename result_of::result_of_asin_c<T>::type_object
                        asin_c(const object_type<T>& x);

// return the inverse cosine function
template<class T>
typename result_of::result_of_acos<T>::type_object
                        acos(const object_type<T>& x);

// compute complex conversion of the acos function
template<class T>
typename result_of::result_of_acos_c<T>::type_object
                        acos_c(const object_type<T>& x);

// return the inverse tangent function of x in radians
template<class T>
typename result_of::result_of_atan<T>::type_object
                        atan(const object_type<T>& x);

// return the inverse cotangent function of x in radians
template<class T>
typename result_of::result_of_acot<T>::type_object
                        acot(const object_type<T>& x);

// return the inverse secant function of x in radians
template<class T>
typename result_of::result_of_asec<T>::type_object
                        asec(const object_type<T>& x);

// compute complex conversion of the asec function
template<class T>
typename result_of::result_of_asec_c<T>::type_object
                        asec_c(const object_type<T>& x);

// return the inverse cosecant function of x in radians
template<class T>
typename result_of::result_of_acsc<T>::type_object
                        acsc(const object_type<T>& x);

// compute complex conversion of the acsc function
template<class T>
typename result_of::result_of_acsc_c<T>::type_object
                        acsc_c(const object_type<T>& x);

// return the inverse hyperbolic sine function
template<class T>
typename result_of::result_of_asinh<T>::type_object
                        asinh(const object_type<T>& x);

// return the inverse hyperbolic cosine function
template<class T>
typename result_of::result_of_acosh<T>::type_object
                        acosh(const object_type<T>& x);

// compute complex conversion of the acosh function
template<class T>
typename result_of::result_of_acosh_c<T>::type_object
                        acosh_c(const object_type<T>& x);

// return the inverse hyperbolic tangent function of x in radians
template<class T>
typename result_of::result_of_atanh<T>::type_object
                        atanh(const object_type<T>& x);

// compute complex conversion of the atanh function
template<class T>
typename result_of::result_of_atanh_c<T>::type_object
                        atanh_c(const object_type<T>& x);

// return the inverse hyperbolic cotangent function of x in radians
template<class T>
typename result_of::result_of_acoth<T>::type_object
                        acoth(const object_type<T>& x);

// compute complex conversion of the acoth function
template<class T>
typename result_of::result_of_acoth_c<T>::type_object
                        acoth_c(const object_type<T>& x);

// return the inverse hyperbolic secant function of x in radians
template<class T>
typename result_of::result_of_asech<T>::type_object
                        asech(const object_type<T>& x);

// compute complex conversion of the asech function
template<class T>
typename result_of::result_of_asech_c<T>::type_object
                        asech_c(const object_type<T>& x);

// return the inverse hyperbolic cosecant function of x in radians
template<class T>
typename result_of::result_of_acsch<T>::type_object
                        acsch(const object_type<T>& x);

//---------------------------------------------------------------------------
//                      ROUNDING FUNCTIONS
//---------------------------------------------------------------------------
// return epsilon value eps, i.e positive distance from abs(x) to the next larger in
// magnitude floating point number of the same precision as x
template<class T>
typename result_of::result_of_eps<T>::type_object
                        eps(const object_type<T>& x);

// round to nearest integer towards -INF
template<class T>
typename result_of::result_of_floor<T>::type_object
                        floor(const object_type<T>& x);

// round to nearest integer towards +INF
template<class T>
typename result_of::result_of_ceil<T>::type_object
                        ceil(const object_type<T>& x);

// round to nearest integer, rounding halfway cases away from zero
template<class T>
typename result_of::result_of_round<T>::type_object
                        round(const object_type<T>& x);

// round to nearest integer towards zero
template<class T>
typename result_of::result_of_trunc<T>::type_object
                        trunc(const object_type<T>& x);

// round to nearest integer towards zero; equivalent to trunc
template<class T>
typename result_of::result_of_trunc<T>::type_object
                        fix(const object_type<T>& x)    { return trunc(x); };

// round to nearest integer towards -INF and cast to integer
template<class T>
typename result_of::result_of_ifloor<T>::type_object
                        ifloor(const object_type<T>& x);

// round to nearest integer towards +INF and cast to integer
template<class T>
typename result_of::result_of_iceil<T>::type_object
                        iceil(const object_type<T>& x);

// round to nearest integer, rounding halfway cases away from zero
// and cast to integer
template<class T>
typename result_of::result_of_iround<T>::type_object
                        iround(const object_type<T>& x);

// round to nearest integer towards zero and cast to integer
template<class T>
typename result_of::result_of_itrunc<T>::type_object
                        itrunc(const object_type<T>& x);

// round to nearest integer towards zero and cast to integer
// equivalent to itrunc
template<class T>
typename result_of::result_of_itrunc<T>::type_object
                        ifix(const object_type<T>& x)   { return itrunc(x); };

//---------------------------------------------------------------------------
//                 VALUE DECOMPOSITION FUNCTIONS
//---------------------------------------------------------------------------
// the sign of x (-1 for negative, 1 for positive, 0 for zero, 
// NaN if x is NaN); for complex arguments return x / abs(x)
template<class T, class Ret = typename result_of::result_of_sign<T>::type_object>
Ret                     sign(const object_type<T>& x);

// the sign of x returned as integer (-1 for negative, 1 for positive, 
// 0 for plus or minus zero); for complex elements return sign of the
// real part
template<class T, class Ret = typename result_of::result_of_isign<T>::type_object>
Ret                     isign(const object_type<T>& x);

// check whether the sign of x is negative
template<class T, class Ret = typename result_of::result_of_signbit<T>::type_object>
Ret                     signbit(const object_type<T>& x);

// decompose given floating point value arg into a normalized fraction
// and an integral power of two, i.e. x is represented as x = d x 2 ^ exp,
// where 0.5 <= |d| < 1; if x is zero, then d = 0 and exp = 0
template<class T, class S, class Ret = typename result_of::result_of_frexp<T,S&>::type_object>
Ret                     frexp(const object_type<T>& x, object_type<S>& exp);
template<class T, class S, class Ret = typename result_of::result_of_frexp<T,S&>::type_object>
Ret                     frexp(const object_type<T>& x, S& exp);
template<class T, class S, class Ret = typename result_of::result_of_frexp<T,S&>::type_object>
Ret                     frexp(const T& x, object_type<S>& exp);

// decompose given floating point value x into integral and fractional
// parts, each having the same type and sign as x, the integral part is 
// stored in the second argument and fractional part is returned
template<class T, class S, class Ret = typename result_of::result_of_modf<T,S&>::type_object>
Ret                     modf(const object_type<T>& x, object_type<S>& int_part);
template<class T, class S, class Ret = typename result_of::result_of_modf<T,S&>::type_object>
Ret                     modf(const object_type<T>& x, S& int_part);
template<class T, class S, class Ret = typename result_of::result_of_modf<T,S&>::type_object>
Ret                     modf(const T& x, object_type<S>& int_part);

//---------------------------------------------------------------------------
//                     BINARY ARITHMETIC FUNCTIONS
//---------------------------------------------------------------------------

// binary plus
// this function is enabled if a function T + S can be found
// by unqualified call in the namespace matcl
template<class T, class S, class Ret = typename result_of::result_of_plus<T,S>::type_object>
Ret                     plus(const object_type<T>& x, const object_type<S>& y);
template<class T, class S, class Ret = typename result_of::result_of_plus<T,S>::type_object>
Ret                     plus(const object_type<T>& x, const S& y);
template<class T, class S, class Ret = typename result_of::result_of_plus<T,S>::type_object>
Ret                     plus(const T& x, const object_type<S>& y);

// binary minus
// this function is enabled if a function T - S can be found
// by unqualified call in the namespace matcl
template<class T, class S, class Ret = typename result_of::result_of_minus<T,S>::type_object>
Ret                     minus(const object_type<T>& x, const object_type<S>& y);
template<class T, class S, class Ret = typename result_of::result_of_minus<T,S>::type_object>
Ret                     minus(const object_type<T>& x, const S& y);
template<class T, class S, class Ret = typename result_of::result_of_minus<T,S>::type_object>
Ret                     minus(const T& x, const object_type<S>& y);

// multiplication
// this function is enabled if a function mul(T,S) can be found
// by unqualified call in the namespace matcl
template<class T, class S, class Ret = typename result_of::result_of_mul<T,S>::type_object>
Ret                     mul(const object_type<T>& x, const object_type<S>& y);
template<class T, class S, class Ret = typename result_of::result_of_mul<T,S>::type_object>
Ret                     mul(const object_type<T>& x, const S& y);
template<class T, class S, class Ret = typename result_of::result_of_mul<T,S>::type_object>
Ret                     mul(const T& x, const object_type<S>& y);

// division; equivalent to a/b
// this function is enabled if a function T / S can be found
// by unqualified call in the namespace matcl
template<class T, class S, class Ret = typename result_of::result_of_div<T,S>::type_object>
Ret                     div(const object_type<T>& x, const object_type<S>& y);
template<class T, class S, class Ret = typename result_of::result_of_div<T,S>::type_object>
Ret                     div(const object_type<T>& x, const S& y);
template<class T, class S, class Ret = typename result_of::result_of_div<T,S>::type_object>
Ret                     div(const T& x, const object_type<S>& y);

// equivalent to a/b but 0/0 gives 0
// this function is enabled if a function div_0(T, S) can be found
// by unqualified call in the namespace matcl
template<class T, class S, class Ret = typename result_of::result_of_div_0<T,S>::type_object>
Ret                     div_0(const object_type<T>& x, const object_type<S>& y);
template<class T, class S, class Ret = typename result_of::result_of_div_0<T,S>::type_object>
Ret                     div_0(const object_type<T>& x, const S& y);
template<class T, class S, class Ret = typename result_of::result_of_div_0<T,S>::type_object>
Ret                     div_0(const T& x, const object_type<S>& y);

// equivalent to a/b but 0/0 gives 1
// this function is enabled if a function div_1(T, S) can be found
// by unqualified call in the namespace matcl and the return type
// has defined one value
template<class T, class S, class Ret = typename result_of::result_of_div_1<T,S>::type_object>
Ret                     div_1(const object_type<T>& x, const object_type<S>& y);
template<class T, class S, class Ret = typename result_of::result_of_div_1<T,S>::type_object>
Ret                     div_1(const object_type<T>& x, const S& y);
template<class T, class S, class Ret = typename result_of::result_of_div_1<T,S>::type_object>
Ret                     div_1(const T& x, const object_type<S>& y);

// return is x - n * y where n = floor(x / y) if y != 0 and 0 otherwise;
// is mod(x,y) is nonzero, then has the same sign as y; if signs of x and y
// are the same, then mod(x,y) = rem(x,y)
template<class T, class S, class Ret = result_of::result_of_mod<T,S>::type_object>
Ret                     mod(const object_type<T>& x, const object_type<S>& y);
template<class T, class S, class Ret = result_of::result_of_mod<T,S>::type_object>
Ret                     mod(const object_type<T>& x, const S& y);
template<class T, class S, class Ret = result_of::result_of_mod<T,S>::type_object>
Ret                     mod(const T& x, const object_type<S>& y);

// return is x - n * y where n = floor(x / y) if y != 0 and 0 otherwise;
// is mod(x,y) is nonzero, then has the same sign as y; if signs of x and y
// are the same, then mod(x,y) = rem(x,y)
template<class T, class S, class Ret = result_of::result_of_rem<T,S>::type_object>
Ret                     rem(const object_type<T>& x, const object_type<S>& y);
template<class T, class S, class Ret = result_of::result_of_rem<T,S>::type_object>
Ret                     rem(const object_type<T>& x, const S& y);
template<class T, class S, class Ret = result_of::result_of_rem<T,S>::type_object>
Ret                     rem(const T& x, const object_type<S>& y);

// four quadrant arctangent of x and y; -pi <= atan2(x, y) <= pi;
template<class T, class S, class Ret = typename result_of::result_of_atan2<T,S>::type_object>
Ret                     atan2(const object_type<T>& x, const object_type<S>& y);
template<class T, class S, class Ret = typename result_of::result_of_atan2<T,S>::type_object>
Ret                     atan2(const object_type<T>& x, const S& y);
template<class T, class S, class Ret = typename result_of::result_of_atan2<T,S>::type_object>
Ret                     atan2(const T& x, const object_type<S>& y);

// calculate sqrt(|x|^2 + |y|^2) avoiding owerflow and underflow at
// intermediate stages of computation
template<class T, class S, class Ret = result_of::result_of_hypot<T,S>::type_object>
Ret                     hypot(const object_type<T>& x, const object_type<S>& y);
template<class T, class S, class Ret = result_of::result_of_hypot<T,S>::type_object>
Ret                     hypot(const object_type<T>& x, const S& y);
template<class T, class S, class Ret = result_of::result_of_hypot<T,S>::type_object>
Ret                     hypot(const T& x, const object_type<S>& y);

// the power function x^y
template<class T, class S, class Ret = typename result_of::result_of_pow<T,S>::type_object>
Ret                     pow(const object_type<T>& x, const object_type<S>& y);
template<class T, class S, class Ret = typename result_of::result_of_pow<T,S>::type_object>
Ret                     pow(const object_type<T>& x, const S& y);
template<class T, class S, class Ret = typename result_of::result_of_pow<T,S>::type_object>
Ret                     pow(const T& x, const object_type<S>& y);

// compute complex conversion of the pow function
template<class T, class S, class Ret = typename result_of::result_of_pow_c<T,S>::type_object>
Ret                     pow_c(const object_type<T>& x, const object_type<S>& y);
template<class T, class S, class Ret = typename result_of::result_of_pow_c<T,S>::type_object>
Ret                     pow_c(const object_type<T>& x, const S& y);
template<class T, class S, class Ret = typename result_of::result_of_pow_c<T,S>::type_object>
Ret                     pow_c(const T& x, const object_type<S>& y);

// return x^y - 1
template<class T, class S, class Ret = typename result_of::result_of_powm1<T,S>::type_object>
Ret                     powm1(const object_type<T>& x, const object_type<S>& y);
template<class T, class S, class Ret = typename result_of::result_of_powm1<T,S>::type_object>
Ret                     powm1(const object_type<T>& x, const S& y);
template<class T, class S, class Ret = typename result_of::result_of_powm1<T,S>::type_object>
Ret                     powm1(const T& x, const object_type<S>& y);

// operators +, -, *, / and idiv are defined in matcl-dynamic
//---------------------------------------------------------------------------
//                      MATRIX FUNCTIONS
//---------------------------------------------------------------------------

// equivalent to x * y; this function is enabled if a function
// T * S can be found by unqualified call in the namespace matcl
template<class T, class S, class Ret = typename result_of::result_of_op_mul<T,S>::type_object>
Ret                     mmul(const object_type<T>& x, const object_type<S>& y);
template<class T, class S, class Ret = typename result_of::result_of_op_mul<T,S>::type_object>
Ret                     mmul(const object_type<T>& x, const S& y);
template<class T, class S, class Ret = typename result_of::result_of_op_mul<T,S>::type_object>
Ret                     mmul(const T& x, const object_type<S>& y);

// equivalent to x * y; this function is enabled if a function
// T * S can be found by unqualified call in the namespace matcl
template<class T, class S, class Ret = typename result_of::result_of_op_mul<T,S>::type_object>
Ret                     kron(const object_type<T>& x, const object_type<S>& y);
template<class T, class S, class Ret = typename result_of::result_of_op_mul<T,S>::type_object>
Ret                     kron(const object_type<T>& x, const S& y);
template<class T, class S, class Ret = typename result_of::result_of_op_mul<T,S>::type_object>
Ret                     kron(const T& x, const object_type<S>& y);

//---------------------------------------------------------------------------
//                      LOGICAL FUNCTIONS
//---------------------------------------------------------------------------
// logical and x && y; this function is enabled if a function
// T && S can be found by unqualified call in the namespace matcl
template<class T, class S, class Ret = typename result_of::result_of_op_and<T,S>::type_object>
bool                    operator&&(const object_type<T>& x, const object_type<S>& y);
template<class T, class S, class Ret = typename result_of::result_of_op_and<T,S>::type_object>
bool                    operator&&(const object_type<T>& x, const S& y);
template<class T, class S, class Ret = typename result_of::result_of_op_and<T,S>::type_object>
bool                    operator&&(const T& x, const object_type<S>& y);

// logical or x || y; this function is enabled if a function
// T || S can be found by unqualified call in the namespace matcl
template<class T, class S, class Ret = typename result_of::result_of_op_or<T,S>::type_object>
bool                    operator||(const object_type<T>& x, const object_type<S>& y);
template<class T, class S, class Ret = typename result_of::result_of_op_or<T,S>::type_object>
bool                    operator||(const object_type<T>& x, const S& y);
template<class T, class S, class Ret = typename result_of::result_of_op_or<T,S>::type_object>
bool                    operator||(const T& x, const object_type<S>& y);

// equivalent to x && y; this function is enabled if a function
// T && S can be found by unqualified call in the namespace matcl
template<class T, class S, class Ret = typename result_of::result_of_op_and<T,S>::type_object>
bool                    op_and(const object_type<T>& x, const object_type<S>& y);
template<class T, class S, class Ret = typename result_of::result_of_op_and<T,S>::type_object>
bool                    op_and(const object_type<T>& x, const S& y);
template<class T, class S, class Ret = typename result_of::result_of_op_and<T,S>::type_object>
bool                    op_and(const T& x, const object_type<S>& y);

// equivalent to x || y; this function is enabled if a function
// T || S can be found by unqualified call in the namespace matcl
template<class T, class S, class Ret = typename result_of::result_of_op_or<T,S>::type_object>
bool                    op_or(const object_type<T>& x, const object_type<S>& y);
template<class T, class S, class Ret = typename result_of::result_of_op_or<T,S>::type_object>
bool                    op_or(const object_type<T>& x, const S& y);
template<class T, class S, class Ret = typename result_of::result_of_op_or<T,S>::type_object>
bool                    op_or(const T& x, const object_type<S>& y);

// logical symmetric difference; both x and y are converted to boolean
// value first; this function is enabled if a function op_xor(T,S)
// can be found by unqualified call in the namespace matcl
template<class T, class S, class Ret = typename result_of::result_of_op_xor<T,S>::type_object>
bool                    op_xor(const object_type<T>& x, const object_type<S>& y);
template<class T, class S, class Ret = typename result_of::result_of_op_xor<T,S>::type_object>
bool                    op_xor(const object_type<T>& x, const S& y);
template<class T, class S, class Ret = typename result_of::result_of_op_xor<T,S>::type_object>
bool                    op_xor(const T& x, const object_type<S>& y);

// element-wise and x & y; this function is enabled if a function
// elem_and(T, S) can be found by unqualified call in the namespace matcl
template<class T, class S, class Ret = typename result_of::result_of_elem_and<T,S>::type_object>
Ret                     operator&(const object_type<T>& x, const object_type<S>& y);
template<class T, class S, class Ret = typename result_of::result_of_elem_and<T,S>::type_object>
Ret                     operator&(const object_type<T>& x, const S& y);
template<class T, class S, class Ret = typename result_of::result_of_elem_and<T,S>::type_object>
Ret                     operator&(const T& x, const object_type<S>& y);

// element-wise or x | y; this function is enabled if a function
// elem_or(T, S) can be found by unqualified call in the namespace matcl
template<class T, class S, class Ret = typename result_of::result_of_elem_or<T,S>::type_object>
Ret                     operator|(const object_type<T>& x, const object_type<S>& y);
template<class T, class S, class Ret = typename result_of::result_of_elem_or<T,S>::type_object>
Ret                     operator|(const object_type<T>& x, const S& y);
template<class T, class S, class Ret = typename result_of::result_of_elem_or<T,S>::type_object>
Ret                     operator|(const T& x, const object_type<S>& y);

// equivalent to x & y; this function is enabled if a function
// elem_and(T, S) can be found by unqualified call in the namespace matcl
template<class T, class S, class Ret = typename result_of::result_of_elem_and<T,S>::type_object>
Ret                     elem_and(const object_type<T>& x, const object_type<S>& y);
template<class T, class S, class Ret = typename result_of::result_of_elem_and<T,S>::type_object>
Ret                     elem_and(const object_type<T>& x, const S& y);
template<class T, class S, class Ret = typename result_of::result_of_elem_and<T,S>::type_object>
Ret                     elem_and(const T& x, const object_type<S>& y);

// equivalent to x | y; this function is enabled if a function
// elem_or(T, S) can be found by unqualified call in the namespace matcl
template<class T, class S, class Ret = typename result_of::result_of_elem_or<T,S>::type_object>
Ret                     elem_or(const object_type<T>& x, const object_type<S>& y);
template<class T, class S, class Ret = typename result_of::result_of_elem_or<T,S>::type_object>
Ret                     elem_or(const object_type<T>& x, const S& y);
template<class T, class S, class Ret = typename result_of::result_of_elem_or<T,S>::type_object>
Ret                     elem_or(const T& x, const object_type<S>& y);

// element-wise symmetric difference; this function is enabled if 
// a function elem_xor(T,S) can be found by unqualified call in the
// namespace matcl
template<class T, class S, class Ret = typename result_of::result_of_elem_xor<T,S>::type_object>
Ret                     elem_xor(const object_type<T>& x, const object_type<S>& y);
template<class T, class S, class Ret = typename result_of::result_of_elem_xor<T,S>::type_object>
Ret                     elem_xor(const object_type<T>& x, const S& y);
template<class T, class S, class Ret = typename result_of::result_of_elem_xor<T,S>::type_object>
Ret                     elem_xor(const T& x, const object_type<S>& y);

//---------------------------------------------------------------------------
//                      COMPARISON FUNCTIONS
//---------------------------------------------------------------------------

// equivalent to x == y; this function is enabled if a function
// T == S can be found by unqualified call in the namespace matcl
template<class T, class S, class Ret = typename result_of::result_of_eeq<T,S>::type_object>
Ret                     eeq(const object_type<T>& x, const object_type<S>& y);
template<class T, class S, class Ret = typename result_of::result_of_eeq<T,S>::type_object>
Ret                     eeq(const object_type<T>& x, const S& y);
template<class T, class S, class Ret = typename result_of::result_of_eeq<T,S>::type_object>
Ret                     eeq(const T& x, const object_type<S>& y);

// return x == y; nan values are considered as equal; this function
// is enabled if a function eeq_nan(T,S) can be found by unqualified
// call in the namespace matcl
template<class T, class S, class Ret = typename result_of::result_of_eeq_nan<T,S>::type_object>
Ret                     eeq_nan(const object_type<T>& x, const object_type<S>& y);
template<class T, class S, class Ret = typename result_of::result_of_eeq_nan<T,S>::type_object>
Ret                     eeq_nan(const object_type<T>& x, const S& y);
template<class T, class S, class Ret = typename result_of::result_of_eeq_nan<T,S>::type_object>
Ret                     eeq_nan(const T& x, const object_type<S>& y);

// equivalent to x != y; this function is enabled if a function
// T != S can be found by unqualified call in the namespace matcl
template<class T, class S, class Ret = typename result_of::result_of_neq<T,S>::type_object>
Ret                     neq(const object_type<T>& x, const object_type<S>& y);
template<class T, class S, class Ret = typename result_of::result_of_neq<T,S>::type_object>
Ret                     neq(const object_type<T>& x, const S& y);
template<class T, class S, class Ret = typename result_of::result_of_neq<T,S>::type_object>
Ret                     neq(const T& x, const object_type<S>& y);

// return x != y; nan values are considered as equal; this function
// is enabled if a function neq_nan(T,S) can be found by unqualified
// call in the namespace matcl
template<class T, class S, class Ret = typename result_of::result_of_neq_nan<T,S>::type_object>
Ret                     neq_nan(const object_type<T>& x, const object_type<S>& y);
template<class T, class S, class Ret = typename result_of::result_of_neq_nan<T,S>::type_object>
Ret                     neq_nan(const object_type<T>& x, const S& y);
template<class T, class S, class Ret = typename result_of::result_of_neq_nan<T,S>::type_object>
Ret                     neq_nan(const T& x, const object_type<S>& y);

// equivalent to x >= y; this function is enabled if a function
// T >= S can be found by unqualified call in the namespace matcl
template<class T, class S, class Ret = typename result_of::result_of_geq<T,S>::type_object>
Ret                     geq(const object_type<T>& x, const object_type<S>& y);
template<class T, class S, class Ret = typename result_of::result_of_geq<T,S>::type_object>
Ret                     geq(const object_type<T>& x, const S& y);
template<class T, class S, class Ret = typename result_of::result_of_geq<T,S>::type_object>
Ret                     geq(const T& x, const object_type<S>& y);

// equivalent to x <= y; this function is enabled if a function
// T <= S can be found by unqualified call in the namespace matcl
template<class T, class S, class Ret = typename result_of::result_of_leq<T,S>::type_object>
Ret                     leq(const object_type<T>& x, const object_type<S>& y);
template<class T, class S, class Ret = typename result_of::result_of_leq<T,S>::type_object>
Ret                     leq(const object_type<T>& x, const S& y);
template<class T, class S, class Ret = typename result_of::result_of_leq<T,S>::type_object>
Ret                     leq(const T& x, const object_type<S>& y);

// equivalent to x > y; this function is enabled if a function
// T > S can be found by unqualified call in the namespace matcl
template<class T, class S, class Ret = typename result_of::result_of_gt<T,S>::type_object>
Ret                     gt(const object_type<T>& x, const object_type<S>& y);
template<class T, class S, class Ret = typename result_of::result_of_gt<T,S>::type_object>
Ret                     gt(const object_type<T>& x, const S& y);
template<class T, class S, class Ret = typename result_of::result_of_gt<T,S>::type_object>
Ret                     gt(const T& x, const object_type<S>& y);

// equivalent to x < y; this function is enabled if a function
// T < S can be found by unqualified call in the namespace matcl
template<class T, class S, class Ret = typename result_of::result_of_lt<T,S>::type_object>
Ret                     lt(const object_type<T>& x, const object_type<S>& y);
template<class T, class S, class Ret = typename result_of::result_of_lt<T,S>::type_object>
Ret                     lt(const object_type<T>& x, const S& y);
template<class T, class S, class Ret = typename result_of::result_of_lt<T,S>::type_object>
Ret                     lt(const T& x, const object_type<S>& y);

// maximum of x and y; this function is enabled if a function
// max(T, S) can be found by unqualified call in the namespace matcl
template<class T, class S, class Ret = typename result_of::result_of_max<T,S>::type_object>
Ret                     max(const object_type<T>& x, const object_type<S>& y);
template<class T, class Ret = typename result_of::result_of_max<T,T>::type_object>
Ret                     max(const object_type<T>& x, const object_type<T>& y);
template<class T, class S, class Ret = typename result_of::result_of_max<T,S>::type_object>
Ret                     max(const object_type<T>& x, const S& y);
template<class T, class S, class Ret = typename result_of::result_of_max<T,S>::type_object>
Ret                     max(const T& x, const object_type<S>& y);

// minimum of x and y; this function is enabled if a function
// min(T, S) can be found by unqualified call in the namespace matcl
template<class T, class S, class Ret = typename result_of::result_of_min<T,S>::type_object>
Ret                     min(const object_type<T>& x, const object_type<S>& y);
template<class T, class Ret = typename result_of::result_of_min<T,T>::type_object>
Ret                     min(const object_type<T>& x, const object_type<T>& y);
template<class T, class S, class Ret = typename result_of::result_of_min<T,S>::type_object>
Ret                     min(const object_type<T>& x, const S& y);
template<class T, class S, class Ret = typename result_of::result_of_min<T,S>::type_object>
Ret                     min(const T& x, const object_type<S>& y);

//---------------------------------------------------------------------------
//                      MISCELLANEOUS BINARY FUNCTIONS
//---------------------------------------------------------------------------

// copy sign from the second operand to the first operand
template<class T, class S, class Ret = result_of::result_of_copysign<T,S>::type_object>
Ret                     copysign(const object_type<T>& x, const object_type<S>& y);
template<class T, class S, class Ret = result_of::result_of_copysign<T,S>::type_object>
Ret                     copysign(const object_type<T>& x, const S& y);
template<class T, class S, class Ret = result_of::result_of_copysign<T,S>::type_object>
Ret                     copysign(const T& x, const object_type<S>& y);

// return the next representable value after x in the direction of y
template<class T, class S, class Ret = result_of::result_of_nextafter<T,S>::type_object>
Ret                     nextafter(const object_type<T>& x, const object_type<S>& y);
template<class T, class S, class Ret = result_of::result_of_nextafter<T,S>::type_object>
Ret                     nextafter(const object_type<T>& x, const S& y);
template<class T, class S, class Ret = result_of::result_of_nextafter<T,S>::type_object>
Ret                     nextafter(const T& x, const object_type<S>& y);

// return the next representable value after x in the direction +INF
template<class T>
typename result_of::result_of_nextabove<T>::type_object
                        nextabove(const object_type<T>& x);

// return the next representable value after x in the direction -INF
template<class T>
typename result_of::result_of_nextbelow<T>::type_object
                        nextbelow(const object_type<T>& x);

// return max(x - y, 0) or NaN if x or y is NaN
template<class T, class S, class Ret = result_of::result_of_fdim<T,S>::type_object>
Ret                     fdim(const object_type<T>& x, const object_type<S>& y);
template<class T, class S, class Ret = result_of::result_of_fdim<T,S>::type_object>
Ret                     fdim(const object_type<T>& x, const S& y);
template<class T, class S, class Ret = result_of::result_of_fdim<T,S>::type_object>
Ret                     fdim(const T& x, const object_type<S>& y);

// return the rising factorial of x and integer i:
// rising_factorial(x, i)  = x*(x+1)*...*(x+i-1)
template<class T, class Ret = result_of::result_of_rising_factorial<T,Integer>::type_object>
Ret                     rising_factorial(const object_type<T>& x, Integer i);

// return the falling factorial of x and integer i:
// falling_factorial(x, i) = x*(x-1)*...*(x-i+1)
template<class T, class Ret = result_of::result_of_falling_factorial<T,Integer>::type_object>
Ret                     falling_factorial(const object_type<T>& x, Integer i);

// return the (2 * n)-th Bernoulli number B_{2n}, n >= 0
template<class T, class TB = typename details::base_object_type<T>::type,
        class Ret = result_of::result_of_bernoulli_b2n<TB,Integer>::type_object>
Ret                     bernoulli_b2n(Integer n);

// return largest value, such that bernoulli_b2n does not overflow
template<class T, class TB = typename T::value_type, 
        class Ret = result_of::result_of_max_bernoulli_b2n<TB>::type_object>
Integer                 max_bernoulli_b2n();

// return the factorial i!, i.e. prod_{k=1}^i k
template<class T, class TB = typename T::value_type, 
        class Ret = result_of::result_of_factorial<TB,Integer>::type_object>
Ret                     factorial(Integer i);

// return double factorial i!!:
//     i!! = i * (i-2) * ... * 1   for i odd
//     i!! = i * (i-2) * ... * 2   for i even
//     i!! = 1                     for i = 0, -1
template<class T, class TB = typename T::value_type, 
        class Ret = result_of::result_of_double_factorial<TB,Integer>::type_object>
Ret                     double_factorial(Integer i);

// returns the binomial coefficient:
//     binomial_coefficient(n,k)   = n! / k! / (n-k)!
// requires k <= n, and k,n >= 0
template<class T, class TB = typename T::value_type, 
        class Ret = result_of::result_of_binomial_coefficient<TB,Integer,Integer>::type_object>
Ret                     binomial_coefficient(Integer n, Integer k);

};};

#include "matcl-scalar/details/typed_object_functions.inl"
