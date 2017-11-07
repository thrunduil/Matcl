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

#include "matcl-mp/details/initializer.h"
#include "matcl-mp/mp_int.h"
#include "matcl-mp/mp_float.h"
#include "matcl-mp/mp_rational.h"
#include "matcl-mp/mp_complex.h"
#include "matcl-mp/details/enablers.h"

namespace matcl
{

// precision argument gives required precision of the function result
// if this argument is zero, then function result will have the same
// precision as maximum precision of input arguments

// functions for real arguments are accurate up to 0.5 ulp (units in the
// last place) unless documented otherwise; functions for complex arguments
// are accurate up to 1 ulp, but can be more accurate (this fact is 
// documented); additionally:
//     - unary minus (-x) is accurate, but uminus(x, p) is accurate up to
//       0.5 ulp if required precision p is less than precision of x (i.e.
//       rounding occures)
//     - complex and real rounding is accurate up to 0.5 ulp

// convert a scalar of type From to a scalar of type To;
// warnings are not printed if this conversion leads to precision lost
template<class To, class From>
inline To                       convert_scalar(const From& s, 
                                    typename mp::details::enable_mp_bin<To,From,void*>::type = 0);

// unary minus
MATCL_MP_EXPORT mp_int          operator-(const mp_int& a);
MATCL_MP_EXPORT mp_float        operator-(const mp_float& a);
MATCL_MP_EXPORT mp_rational     operator-(const mp_rational& a);
MATCL_MP_EXPORT mp_complex      operator-(const mp_complex& a);

// equivaluent to operator-
MATCL_MP_EXPORT mp_int          uminus(const mp_int& a, precision p = precision());
MATCL_MP_EXPORT mp_float        uminus(const mp_float& a, precision p = precision());
MATCL_MP_EXPORT mp_rational     uminus(const mp_rational& a, precision p = precision());
MATCL_MP_EXPORT mp_complex      uminus(const mp_complex& a, precision p = precision());

// return the inverse, i.e. 1/x; returned value has the same precision as x
// if x is integer or rational zero, then zero is returned and integer
// overflow flag is set
MATCL_MP_EXPORT mp_rational     inv(const mp_int& x, precision p = precision());
MATCL_MP_EXPORT mp_float        inv(const mp_float& x, precision p = precision());
MATCL_MP_EXPORT mp_rational     inv(const mp_rational& x, precision p = precision());
MATCL_MP_EXPORT mp_complex      inv(const mp_complex& x, precision p = precision());

// return real part; this function is accurate
inline const mp_int&            real(const mp_int& a)       { return a; }; 
inline const mp_float&          real(const mp_float& a)     { return a; }; 
inline const mp_rational&       real(const mp_rational& a)  { return a; }; 
inline const mp_float&          real(const mp_complex& a)   { return a.real(); }; 

// return imaginary part; this function is accurate
inline mp_int                   imag(const mp_int&)         { return mp_int(); };
inline mp_float                 imag(const mp_float& a)     { return mp_float(a.get_precision()); };
inline mp_rational              imag(const mp_rational&)    { return mp_rational(); };
inline const mp_float&          imag(const mp_complex& a)   { return a.imag(); };

// conjugate transposition; this function is accurate
inline const mp_int&            conj(const mp_int& a)       { return a; };
inline const mp_float&          conj(const mp_float& a)     { return a; };
inline const mp_rational&       conj(const mp_rational& a)  { return a; };
MATCL_MP_EXPORT mp_complex      conj(const mp_complex& a);

// check if value is nan
inline bool                     is_nan(const mp_int&)           { return false; };
MATCL_MP_EXPORT bool            is_nan(const mp_float& val);
inline bool                     is_nan(const mp_rational&)      { return false; };
MATCL_MP_EXPORT bool            is_nan(const mp_complex& val);

// check if value is infinite
inline bool                     is_inf(const mp_int&)           { return false; };
MATCL_MP_EXPORT bool            is_inf(const mp_float& val);
inline bool                     is_inf(const mp_rational&)      { return false; };
MATCL_MP_EXPORT bool            is_inf(const mp_complex& val);

// check if value is finite
inline bool                     is_finite(const mp_int&)        { return true; };
MATCL_MP_EXPORT bool            is_finite(const mp_float& val);
inline bool                     is_finite(const mp_rational&)   { return true; };
MATCL_MP_EXPORT bool            is_finite(const mp_complex& val);

// check if value is regular (i.e., neither NaN, nor an infinity nor zero)
// for a complex value false is retured if real or imaginary part is
// NaN or infinity or if both real and imaginary part is zero
MATCL_MP_EXPORT bool            is_regular(const mp_int& val);
MATCL_MP_EXPORT bool            is_regular(const mp_float& val);
MATCL_MP_EXPORT bool            is_regular(const mp_rational& val);
MATCL_MP_EXPORT bool            is_regular(const mp_complex& val);

// check if value is regular (i.e., neither NaN, nor an infinity 
// nor zero nor subnormal); equivalent to is_regular since multiprecision
// floating point cannot be subnormal
MATCL_MP_EXPORT bool            is_normal(const mp_int& val);
MATCL_MP_EXPORT bool            is_normal(const mp_float& val);
MATCL_MP_EXPORT bool            is_normal(const mp_rational& val);
MATCL_MP_EXPORT bool            is_normal(const mp_complex& val);

// check if value is an integer
inline bool                     is_int(const mp_int&)           { return true; }
MATCL_MP_EXPORT bool            is_int(const mp_float& val);
MATCL_MP_EXPORT bool            is_int(const mp_rational&);
MATCL_MP_EXPORT bool            is_int(const mp_complex& val);

// check if value is real
inline bool                     is_real(const mp_int&)          { return true; };
inline bool                     is_real(const mp_float&)        { return true; };
inline bool                     is_real(const mp_rational&)     { return true; };
MATCL_MP_EXPORT bool            is_real(const mp_complex&);

// check if value is zero
MATCL_MP_EXPORT bool            is_zero(const mp_int& val);
MATCL_MP_EXPORT bool            is_zero(const mp_float& val);
MATCL_MP_EXPORT bool            is_zero(const mp_rational& val);
MATCL_MP_EXPORT bool            is_zero(const mp_complex& val);

// check if value is one
MATCL_MP_EXPORT bool            is_one(const mp_int& val);
MATCL_MP_EXPORT bool            is_one(const mp_float& val);
MATCL_MP_EXPORT bool            is_one(const mp_rational& val);
MATCL_MP_EXPORT bool            is_one(const mp_complex& val);

// return value class of floating point number; see also fp_type
// enumeration; not defined for complex values
MATCL_MP_EXPORT fp_type         fpclassify(const mp_int& val);
MATCL_MP_EXPORT fp_type         fpclassify(const mp_float& val);
MATCL_MP_EXPORT fp_type         fpclassify(const mp_rational& val);

// absolute value; real and complex abs is accurate up to
// 0.5 ulp; abs for real arguments is exact if no rounding is required
MATCL_MP_EXPORT mp_int          abs(const mp_int& val, precision p = precision());
MATCL_MP_EXPORT mp_float        abs(const mp_float& val, precision p = precision());
MATCL_MP_EXPORT mp_rational     abs(const mp_rational& val, precision p = precision());
MATCL_MP_EXPORT mp_float        abs(const mp_complex& val, precision p = precision());

// absolute value squared
MATCL_MP_EXPORT mp_int          abs2(const mp_int& val, precision p = precision());
MATCL_MP_EXPORT mp_float        abs2(const mp_float& val, precision p = precision());
MATCL_MP_EXPORT mp_rational     abs2(const mp_rational& val, precision p = precision());
MATCL_MP_EXPORT mp_float        abs2(const mp_complex& val, precision p = precision());

// argument (polar angle) of a complex number
MATCL_MP_EXPORT mp_float        arg(const mp_int& val, precision p = precision());
MATCL_MP_EXPORT mp_float        arg(const mp_float& val, precision p = precision());
MATCL_MP_EXPORT mp_float        arg(const mp_rational& val, precision p = precision());
MATCL_MP_EXPORT mp_float        arg(const mp_complex& val, precision p = precision());

// return the next representable value after x in the direction +INF
// if x is a complex value, then imaginary part is ignored
MATCL_MP_EXPORT mp_float        nextabove(const mp_int& x);
MATCL_MP_EXPORT mp_float        nextabove(const mp_float& x);
MATCL_MP_EXPORT mp_float        nextabove(const mp_rational& x);

// return the next representable value after x in the direction -INF
// if x is a complex value, then imaginary part is ignored
MATCL_MP_EXPORT mp_float        nextbelow(const mp_int& x);
MATCL_MP_EXPORT mp_float        nextbelow(const mp_float& x);
MATCL_MP_EXPORT mp_float        nextbelow(const mp_rational& x);

// the sign of x (-1 for negative, 1 for positive, 0 for zero, 
// NaN if x is NaN); for complex arguments return x / abs(x),
// i.e. x normalized to the unit circle
MATCL_MP_EXPORT mp_int          sign(const mp_int& x, precision p = precision());
MATCL_MP_EXPORT mp_float        sign(const mp_float& x, precision p = precision());
MATCL_MP_EXPORT mp_rational     sign(const mp_rational& x, precision p = precision());
MATCL_MP_EXPORT mp_complex      sign(const mp_complex& x, precision p = precision());

// the sign of x returned as integer (-1 for negative, 1 for positive, 
// 0 for plus or minus zero); not available for complex numbers
MATCL_MP_EXPORT Integer         isign(const mp_int& x);
MATCL_MP_EXPORT Integer         isign(const mp_float& x);
MATCL_MP_EXPORT Integer         isign(const mp_rational& x);

// check whether the sign of x is negative; not available for complex numbers
MATCL_MP_EXPORT bool            signbit(const mp_int& x);
MATCL_MP_EXPORT bool            signbit(const mp_float& x);
MATCL_MP_EXPORT bool            signbit(const mp_rational& x);

// scale x by 2 raised to the power of n, returning ldexp(x,n) = x * 2^n
// returned result is exact
MATCL_MP_EXPORT mp_rational     ldexp(const mp_int& x, Integer n);
MATCL_MP_EXPORT mp_float        ldexp(const mp_float& x, Integer n);
MATCL_MP_EXPORT mp_rational     ldexp(const mp_rational& x, Integer n);
MATCL_MP_EXPORT mp_complex      ldexp(const mp_complex& x, Integer n);

// function equivalent to ldexp
MATCL_MP_EXPORT mp_rational     scalbn(const mp_int& x, Integer n);
MATCL_MP_EXPORT mp_float        scalbn(const mp_float& x, Integer n);
MATCL_MP_EXPORT mp_rational     scalbn(const mp_rational& x, Integer n);
MATCL_MP_EXPORT mp_complex      scalbn(const mp_complex& x, Integer n);

// decompose given floating point value arg into a normalized fraction
// and an integral power of two, i.e. x is represented as x = d x 2 ^ exp,
// where 0.5 <= |d| < 1; if x is zero, then d = 0 and exp = 0; if is Inf
// or NaN, then d = x; returned result is exact
// not available for complex numbers
MATCL_MP_EXPORT mp_float        frexp(const mp_int& x, Integer& exp);
MATCL_MP_EXPORT mp_float        frexp(const mp_float& x, Integer& exp);
MATCL_MP_EXPORT mp_float        frexp(const mp_rational& x, Integer& exp);

// decompose given floating point value x into integral and fractional
// parts, each having the same type and sign as x, the integral part is 
// stored in int_part argument and fractional part is returned;
// not available for complex numbers
MATCL_MP_EXPORT mp_float        modf(const mp_int& x, mp_float& int_part);
MATCL_MP_EXPORT mp_float        modf(const mp_float& x, mp_float& int_part);
MATCL_MP_EXPORT mp_float        modf(const mp_rational& x, mp_float& int_part);

// extract the value of the unbiased exponent from the floating-point
// argument x, and returns it as a floating-point value; equivalent to
// log2(|x|) rounded to integer number toward -INF; additionally for 
// regular values logb(x) = exp - 1, where exp is the exponent returned 
// by frexp, logb(0) = -inf, logb(+- inf) = inf, logb(nan) = nan;
// for complex argument logb(|x|) is returned
MATCL_MP_EXPORT mp_float        logb(const mp_int& x);
MATCL_MP_EXPORT mp_float        logb(const mp_float& x);
MATCL_MP_EXPORT mp_float        logb(const mp_rational& x);
MATCL_MP_EXPORT mp_float        logb(const mp_complex& x);

// cast result of logb(x) to integer
MATCL_MP_EXPORT Integer         ilogb(const mp_int& x);
MATCL_MP_EXPORT Integer         ilogb(const mp_float& x);
MATCL_MP_EXPORT Integer         ilogb(const mp_rational& x);
MATCL_MP_EXPORT Integer         ilogb(const mp_complex& x);

// return epsilon value eps, i.e positive distance from abs(x)
// to the next larger in magnitude floating point number of the
// same precision as x
MATCL_MP_EXPORT mp_float        eps(const mp_int& x);
MATCL_MP_EXPORT mp_float        eps(const mp_float& x);
MATCL_MP_EXPORT mp_float        eps(const mp_rational& x);
MATCL_MP_EXPORT mp_float        eps(const mp_complex& x);

// return the square root of x
MATCL_MP_EXPORT mp_float        sqrt(const mp_int& x, precision p = precision());
MATCL_MP_EXPORT mp_float        sqrt(const mp_float& x, precision p = precision());
MATCL_MP_EXPORT mp_float        sqrt(const mp_rational& x, precision p = precision());
MATCL_MP_EXPORT mp_complex      sqrt(const mp_complex& x, precision p = precision());

// return the cubic root of x, for negative values result
// is also negative; not available for complex values
MATCL_MP_EXPORT mp_float        cbrt(const mp_int& x, precision p = precision());
MATCL_MP_EXPORT mp_float        cbrt(const mp_float& x, precision p = precision());
MATCL_MP_EXPORT mp_float        cbrt(const mp_rational& x, precision p = precision());

// return the exponential function, i.e. the e (Euler's number,
// 2.7182818) raised to the given power x; for complex z = x + i*y, 
// exp(z) = exp(x) * (cos(y) + i*sin(Y))
MATCL_MP_EXPORT mp_float        exp(const mp_int& x, precision p = precision());
MATCL_MP_EXPORT mp_float        exp(const mp_float& x, precision p = precision());
MATCL_MP_EXPORT mp_float        exp(const mp_rational& x, precision p = precision());
MATCL_MP_EXPORT mp_complex      exp(const mp_complex& x, precision p = precision());

// return 2^x
MATCL_MP_EXPORT mp_float        exp2(const mp_int& x, precision p = precision());
MATCL_MP_EXPORT mp_float        exp2(const mp_float& x, precision p = precision());
MATCL_MP_EXPORT mp_float        exp2(const mp_rational& x, precision p = precision());
MATCL_MP_EXPORT mp_complex      exp2(const mp_complex& x, precision p = precision());

// return 10^x
MATCL_MP_EXPORT mp_float        exp10(const mp_int& x, precision p = precision());
MATCL_MP_EXPORT mp_float        exp10(const mp_float& x, precision p = precision());
MATCL_MP_EXPORT mp_float        exp10(const mp_rational& x, precision p = precision());
MATCL_MP_EXPORT mp_complex      exp10(const mp_complex& x, precision p = precision());

// return exp(x) - 1
MATCL_MP_EXPORT mp_float        expm1(const mp_int& x, precision p = precision());
MATCL_MP_EXPORT mp_float        expm1(const mp_float& x, precision p = precision());
MATCL_MP_EXPORT mp_float        expm1(const mp_rational& x, precision p = precision());
MATCL_MP_EXPORT mp_complex      expm1(const mp_complex& x, precision p = precision());

// return exp(x * 1i), where 1i is an imaginary unit
MATCL_MP_EXPORT mp_complex      expi(const mp_int& x, precision p = precision());
MATCL_MP_EXPORT mp_complex      expi(const mp_float& x, precision p = precision());
MATCL_MP_EXPORT mp_complex      expi(const mp_rational& x, precision p = precision());
MATCL_MP_EXPORT mp_complex      expi(const mp_complex& x, precision p = precision());

// computes the natural (base e, Euler's number, 2.7182818) logarithm 
// of x; for complex argument z = x + i*y compute logarithm of a complex 
// value z with a branch cut along the negative real axis; if no errors occur,
// the complex natural logarithm of z is returned, in the range of a strip in
// the interval [-i pi, +i pi] along the imaginary axis and mathematically 
// unbounded along the real axis;
// Log(z) = log(|z|) + i Arg(z) = log( hypot(x,y) ) + i atan2(y,x)
MATCL_MP_EXPORT mp_float        log(const mp_int& x, precision p = precision());
MATCL_MP_EXPORT mp_float        log(const mp_float& x, precision p = precision());
MATCL_MP_EXPORT mp_float        log(const mp_rational& x, precision p = precision());
MATCL_MP_EXPORT mp_complex      log(const mp_complex& x, precision p = precision());

// compute log(1+x); 
MATCL_MP_EXPORT mp_float        log1p(const mp_int& x, precision p = precision());
MATCL_MP_EXPORT mp_float        log1p(const mp_float& x, precision p = precision());
MATCL_MP_EXPORT mp_float        log1p(const mp_rational& x, precision p = precision());
MATCL_MP_EXPORT mp_complex      log1p(const mp_complex& x, precision p = precision());

// computes the base 2 logarithm of x; see also log
MATCL_MP_EXPORT mp_float        log2(const mp_int& x, precision p = precision());
MATCL_MP_EXPORT mp_float        log2(const mp_float& x, precision p = precision());
MATCL_MP_EXPORT mp_float        log2(const mp_rational& x, precision p = precision());
MATCL_MP_EXPORT mp_complex      log2(const mp_complex& x, precision p = precision());

// computes the base 10 logarithm of x; see also log
MATCL_MP_EXPORT mp_float        log10(const mp_int& x, precision p = precision());
MATCL_MP_EXPORT mp_float        log10(const mp_float& x, precision p = precision());
MATCL_MP_EXPORT mp_float        log10(const mp_rational& x, precision p = precision());
MATCL_MP_EXPORT mp_complex      log10(const mp_complex& x, precision p = precision());

// compute simultaneously the sine of x and the cosine of x; not defined for
// a complex argument
MATCL_MP_EXPORT void            sin_cos(const mp_int& x, mp_float& sin, mp_float& cos, 
                                        precision p = precision());
MATCL_MP_EXPORT void            sin_cos(const mp_float& x, mp_float& sin, mp_float& cos, 
                                        precision p = precision());
MATCL_MP_EXPORT void            sin_cos(const mp_rational& x, mp_float& sin, mp_float& cos, 
                                        precision p = precision());

// return the sine function of x in radians.
// for complex z = x + i*y, sin(z) = sin(x) * cosh(y) + i cos(x) * sinh(y)
MATCL_MP_EXPORT mp_float        sin(const mp_int& x, precision p = precision());
MATCL_MP_EXPORT mp_float        sin(const mp_float& x, precision p = precision());
MATCL_MP_EXPORT mp_float        sin(const mp_rational& x, precision p = precision());
MATCL_MP_EXPORT mp_complex      sin(const mp_complex& x, precision p = precision());

// return the cosine function of x in radians.
// for complex z = x + i*y, cos(z) = cos(x) * cosh(y) - i sin(x) * sinh(y)
MATCL_MP_EXPORT mp_float        cos(const mp_int& x, precision p = precision());
MATCL_MP_EXPORT mp_float        cos(const mp_float& x, precision p = precision());
MATCL_MP_EXPORT mp_float        cos(const mp_rational& x, precision p = precision());
MATCL_MP_EXPORT mp_complex      cos(const mp_complex& x, precision p = precision());

// return the tangent function of x in radians.
// for complex z = x + i*y, tan(z) = -i * tanh(i*z)
MATCL_MP_EXPORT mp_float        tan(const mp_int& x, precision p = precision());
MATCL_MP_EXPORT mp_float        tan(const mp_float& x, precision p = precision());
MATCL_MP_EXPORT mp_float        tan(const mp_rational& x, precision p = precision());
MATCL_MP_EXPORT mp_complex      tan(const mp_complex& x, precision p = precision());

// return the cotangent function of x in radians.
// for complex z = x + i*y, cot(z) = i * coth(i*z)
MATCL_MP_EXPORT mp_float        cot(const mp_int& x, precision p = precision());
MATCL_MP_EXPORT mp_float        cot(const mp_float& x, precision p = precision());
MATCL_MP_EXPORT mp_float        cot(const mp_rational& x, precision p = precision());
MATCL_MP_EXPORT mp_complex      cot(const mp_complex& x, precision p = precision());

// return the secant function of x in radians.
// for complex z = x + i*y, sec(z) = 1 / cos(z)
MATCL_MP_EXPORT mp_float        sec(const mp_int& x, precision p = precision());
MATCL_MP_EXPORT mp_float        sec(const mp_float& x, precision p = precision());
MATCL_MP_EXPORT mp_float        sec(const mp_rational& x, precision p = precision());
MATCL_MP_EXPORT mp_complex      sec(const mp_complex& x, precision p = precision());

// return the cosecant function of x in radians.
// for complex z = x + i*y, csc(z) = 1 / sin(z)
MATCL_MP_EXPORT mp_float        csc(const mp_int& x, precision p = precision());
MATCL_MP_EXPORT mp_float        csc(const mp_float& x, precision p = precision());
MATCL_MP_EXPORT mp_float        csc(const mp_rational& x, precision p = precision());
MATCL_MP_EXPORT mp_complex      csc(const mp_complex& x, precision p = precision());

// compute simultaneously the hyperbolic sine of x and the hyperbolic cosine of x;
// not defined for a complex argument
MATCL_MP_EXPORT void            sinh_cosh(const mp_int& x, mp_float& sinh, mp_float& cosh, 
                                        precision p = precision());
MATCL_MP_EXPORT void            sinh_cosh(const mp_float& x, mp_float& sinh, mp_float& cosh, 
                                        precision p = precision());
MATCL_MP_EXPORT void            sinh_cosh(const mp_rational& x, mp_float& sinh, mp_float& cosh, 
                                        precision p = precision());

// return the hyperbolic sine function of x in radians.
// for complex z = x + i*y, sinh(z) = sinh(x) * cos(y) + i cosh(x) * sin(y)
MATCL_MP_EXPORT mp_float        sinh(const mp_int& x, precision p = precision());
MATCL_MP_EXPORT mp_float        sinh(const mp_float& x, precision p = precision());
MATCL_MP_EXPORT mp_float        sinh(const mp_rational& x, precision p = precision());
MATCL_MP_EXPORT mp_complex      sinh(const mp_complex& x, precision p = precision());

// return the hyperbolic cosine function of x in radians.
// for complex z = x + i*y, cosh(z) = cosh(x) * cos(y) + i sinh(x) * sin(y)
MATCL_MP_EXPORT mp_float        cosh(const mp_int& x, precision p = precision());
MATCL_MP_EXPORT mp_float        cosh(const mp_float& x, precision p = precision());
MATCL_MP_EXPORT mp_float        cosh(const mp_rational& x, precision p = precision());
MATCL_MP_EXPORT mp_complex      cosh(const mp_complex& x, precision p = precision());

// return the hyperbolic tangent function of x in radians.
// for complex z = x + i*y, tanh(z) = sinh(z) / cosh(z)
MATCL_MP_EXPORT mp_float        tanh(const mp_int& x, precision p = precision());
MATCL_MP_EXPORT mp_float        tanh(const mp_float& x, precision p = precision());
MATCL_MP_EXPORT mp_float        tanh(const mp_rational& x, precision p = precision());
MATCL_MP_EXPORT mp_complex      tanh(const mp_complex& x, precision p = precision());

// return the hyperbolic cotangent function of x in radians.
// for complex z = x + i*y, coth(z) = cosh(z) / sinh(z)
MATCL_MP_EXPORT mp_float        coth(const mp_int& x, precision p = precision());
MATCL_MP_EXPORT mp_float        coth(const mp_float& x, precision p = precision());
MATCL_MP_EXPORT mp_float        coth(const mp_rational& x, precision p = precision());
MATCL_MP_EXPORT mp_complex      coth(const mp_complex& x, precision p = precision());

// return the hyperbolic secant function of x in radians.
// for complex z = x + i*y, sech(z) = 1 / cosh(z)
MATCL_MP_EXPORT mp_float        sech(const mp_int& x, precision p = precision());
MATCL_MP_EXPORT mp_float        sech(const mp_float& x, precision p = precision());
MATCL_MP_EXPORT mp_float        sech(const mp_rational& x, precision p = precision());
MATCL_MP_EXPORT mp_complex      sech(const mp_complex& x, precision p = precision());

// return the hyperbolic cosecant function of x in radians.
// for complex z = x + i*y, csch(z) = 1 / sinh(z)
MATCL_MP_EXPORT mp_float        csch(const mp_int& x, precision p = precision());
MATCL_MP_EXPORT mp_float        csch(const mp_float& x, precision p = precision());
MATCL_MP_EXPORT mp_float        csch(const mp_rational& x, precision p = precision());
MATCL_MP_EXPORT mp_complex      csch(const mp_complex& x, precision p = precision());

// return the inverse sine function of x, result in radians;
// for complex z:
// asin z = 1/i log(iz + sqrt(|1 - z^2|) * exp(i/2 arg(1-z^2)))
MATCL_MP_EXPORT mp_float        asin(const mp_int& x, precision p = precision());
MATCL_MP_EXPORT mp_float        asin(const mp_float& x, precision p = precision());
MATCL_MP_EXPORT mp_float        asin(const mp_rational& x, precision p = precision());
MATCL_MP_EXPORT mp_complex      asin(const mp_complex& x, precision p = precision());

// return the inverse cosine function of x, result in radians
// for complex z:
// acos z = 1/i log(z + i * sqrt(|1 - z^2|) * exp(i/2 arg(1-z^2)))
MATCL_MP_EXPORT mp_float        acos(const mp_int& x, precision p = precision());
MATCL_MP_EXPORT mp_float        acos(const mp_float& x, precision p = precision());
MATCL_MP_EXPORT mp_float        acos(const mp_rational& x, precision p = precision());
MATCL_MP_EXPORT mp_complex      acos(const mp_complex& x, precision p = precision());

// return the inverse tangent function of x, result in radians;
// for a complex z, atan(z) = -i atanh(iz)
MATCL_MP_EXPORT mp_float        atan(const mp_int& x, precision p = precision());
MATCL_MP_EXPORT mp_float        atan(const mp_float& x, precision p = precision());
MATCL_MP_EXPORT mp_float        atan(const mp_rational& x, precision p = precision());
MATCL_MP_EXPORT mp_complex      atan(const mp_complex& x, precision p = precision());

// return the inverse hyperbolic sine function of x, result in radians
// for a complex z, asinh(z) = i asin(-i z)
MATCL_MP_EXPORT mp_float        asinh(const mp_int& x, precision p = precision());
MATCL_MP_EXPORT mp_float        asinh(const mp_float& x, precision p = precision());
MATCL_MP_EXPORT mp_float        asinh(const mp_rational& x, precision p = precision());
MATCL_MP_EXPORT mp_complex      asinh(const mp_complex& x, precision p = precision());

// return the inverse hyperbolic cosine function of x, result in radians
// for a complex z, acosh(z) = +-i acos(z) with posive sign if imag(acos(z)) < 0,
// thus real( acosh(z) ) >= 0
MATCL_MP_EXPORT mp_float        acosh(const mp_int& x, precision p = precision());
MATCL_MP_EXPORT mp_float        acosh(const mp_float& x, precision p = precision());
MATCL_MP_EXPORT mp_float        acosh(const mp_rational& x, precision p = precision());
MATCL_MP_EXPORT mp_complex      acosh(const mp_complex& x, precision p = precision());

// return the inverse hyperbolic tangent function of x, result in radians
// for complex z = x + i*y, atanh(z) = re + im *i, where
// re = log1p(4 * x / ( (x - 1) * (x - 1) + y^2))
// im = atan(2y, (1 - x) * (1 + x) - y^2 ) / 2
MATCL_MP_EXPORT mp_float        atanh(const mp_int& x, precision p = precision());
MATCL_MP_EXPORT mp_float        atanh(const mp_float& x, precision p = precision());
MATCL_MP_EXPORT mp_float        atanh(const mp_rational& x, precision p = precision());
MATCL_MP_EXPORT mp_complex      atanh(const mp_complex& x, precision p = precision());

// round to nearest integer towards -INF; for complex numbers
// real and imaginary part are rounded separately
MATCL_MP_EXPORT mp_int          floor(const mp_int& x);
MATCL_MP_EXPORT mp_float        floor(const mp_float& x);
MATCL_MP_EXPORT mp_float        floor(const mp_rational& x);
MATCL_MP_EXPORT mp_complex      floor(const mp_complex& x);

// round to nearest integer towards +INF; for complex numbers
// real and imaginary part are rounded separately
MATCL_MP_EXPORT mp_int          ceil(const mp_int& x);
MATCL_MP_EXPORT mp_float        ceil(const mp_float& x);
MATCL_MP_EXPORT mp_float        ceil(const mp_rational& x);
MATCL_MP_EXPORT mp_complex      ceil(const mp_complex& x);

// round to nearest integer, rounding halfway cases to nearest even integer;
// for complex numbers real and imaginary part are rounded separately
MATCL_MP_EXPORT mp_int          round(const mp_int& x);
MATCL_MP_EXPORT mp_float        round(const mp_float& x);
MATCL_MP_EXPORT mp_float        round(const mp_rational& x);
MATCL_MP_EXPORT mp_complex      round(const mp_complex& x);

// round to nearest integer towards zero; for complex numbers
// real and imaginary part are rounded separately
MATCL_MP_EXPORT mp_int          trunc(const mp_int& x);
MATCL_MP_EXPORT mp_float        trunc(const mp_float& x);
MATCL_MP_EXPORT mp_float        trunc(const mp_rational& x);
MATCL_MP_EXPORT mp_complex      trunc(const mp_complex& x);

// round to nearest integer towards -INF and cast to integer;
// not available for complex numbers
MATCL_MP_EXPORT mp_int          ifloor(const mp_int& x);
MATCL_MP_EXPORT mp_int          ifloor(const mp_float& x);
MATCL_MP_EXPORT mp_int          ifloor(const mp_rational& x);

// round to nearest integer towards +INF and cast to integer;
// not available for complex numbers
MATCL_MP_EXPORT mp_int          iceil(const mp_int& x);
MATCL_MP_EXPORT mp_int          iceil(const mp_float& x);
MATCL_MP_EXPORT mp_int          iceil(const mp_rational& x);

// round to nearest integer, rounding halfway cases away from zero
// not available for complex numbers
MATCL_MP_EXPORT mp_int          iround(const mp_int& x);
MATCL_MP_EXPORT mp_int          iround(const mp_float& x);
MATCL_MP_EXPORT mp_int          iround(const mp_rational& x);

// round to nearest integer towards zero and cast to integer; 
// not available for complex numbers
MATCL_MP_EXPORT mp_int          itrunc(const mp_int& x);
MATCL_MP_EXPORT mp_int          itrunc(const mp_float& x);
MATCL_MP_EXPORT mp_int          itrunc(const mp_rational& x);

// return the factorial i!, i.e. prod_{k=1}^i k
MATCL_MP_EXPORT mp_int          mpi_factorial(Integer i, precision p = precision());
MATCL_MP_EXPORT mp_float        mpf_factorial(Integer i, precision p = precision());

template<class S1, class Enable = typename mp::details::enable_mp<S1>::type>
MATCL_MP_EXPORT S1              factorial(Integer i, precision p = precision());

// return double factorial i!!:
//     i!! = i * (i-2) * ... * 1   for i odd
//     i!! = i * (i-2) * ... * 2   for i even
//     i!! = 1                     for i = 0, -1
MATCL_MP_EXPORT mp_int          mpi_double_factorial(Integer i, precision p = precision());
MATCL_MP_EXPORT mp_float        mpf_double_factorial(Integer i, precision p = precision());

template<class S1, class Enable = typename mp::details::enable_mp<S1>::type>
MATCL_MP_EXPORT S1              double_factorial(Integer i, precision p = precision());

// returns the binomial coefficient:
//     binomial_coefficient(n,k)   = n! / k! / (n-k)!
// requires k <= n, and k,n >= 0
MATCL_MP_EXPORT mp_int          mpi_binomial_coefficient(Integer n, Integer k, precision p = precision());
MATCL_MP_EXPORT mp_float        mpf_binomial_coefficient(Integer n, Integer k, precision p = precision());

template<class S1, class Enable = typename mp::details::enable_mp<S1>::type>
MATCL_MP_EXPORT S1              binomial_coefficient(Integer n, Integer k, precision p = precision());

};

#include "matcl-mp/details/func_unary.inl"