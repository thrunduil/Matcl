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
#include "matcl-mp/details/arithmetic.h"
#include "matcl-mp/details/comparison.h"

namespace matcl
{

// generally functions for real arguments are accurate up to 0.5 ulp 
// (units in the last place); functions for complex arguments are accurate
// up to 1 ulp; additionally functions plus, minus (and associaced operators)
// for complex arguments are accurate up to 0.5 ulp; accuracy can be different
// if this is stated in description of given function

// equality comparison, x == y
template<class T1, class T2, 
        class Enable = typename mp::details::enable_mp_bin<T1,T2,void>::type>
bool            operator==(const T1& x, const T2& y);

// equality comparison, x == y; nan values are equal
template<class T1, class T2, 
        class Enable = typename mp::details::enable_mp_bin<T1,T2,void>::type>
bool            eeq_nan(const T1& x, const T2& y);

// inequality comparison, x != y
template<class T1, class T2, 
        class Enable = typename mp::details::enable_mp_bin<T1,T2,void>::type>
bool            operator!=(const T1& x, const T2& y);

// inequality comparison, x != y; nan values are equal
template<class T1, class T2, 
        class Enable = typename mp::details::enable_mp_bin<T1,T2,void>::type>
bool            neq_nan(const T1& x, const T2& y);

// greater or equal comparison, x >= y
template<class T1, class T2, 
        class Enable = typename mp::details::enable_mp_bin<T1,T2,void>::type>
bool            operator>=(const T1& x, const T2& y);

// less or equal comparison, x <= y
template<class T1, class T2, 
        class Enable = typename mp::details::enable_mp_bin<T1,T2,void>::type>
bool            operator<=(const T1& x, const T2& y);

// greater than comparison, x > y
template<class T1, class T2, 
        class Enable = typename mp::details::enable_mp_bin<T1,T2,void>::type>
bool            operator>(const T1& x, const T2& y);

// less than comparison, x < y
template<class T1, class T2, 
        class Enable = typename mp::details::enable_mp_bin<T1,T2,void>::type>
bool            operator<(const T1& x, const T2& y);

// addition x + y
template<class T1, class T2, 
        class Enable = typename mp::details::enable_mp_bin<T1,T2,void>::type,
        class Result = decltype(mp::details::plus_impl(T1(),T2(),precision()))>
Result          operator+(const T1& x, const T2& y);

// subtraction x - y
template<class T1, class T2, 
        class Enable = typename mp::details::enable_mp_bin<T1,T2,void>::type,
        class Result = decltype(mp::details::minus_impl(T1(),T2(),precision()))>
Result          operator-(const T1& x, const T2& y);

// multiplication x * y
template<class T1, class T2, 
        class Enable = typename mp::details::enable_mp_bin<T1,T2,void>::type,
        class Result = decltype(mp::details::mul_impl(T1(),T2(),precision()))>
Result          operator*(const T1& x, const T2& y);

// division x / y; integer division gives rational result
// if x and y are integer or rational values and y is zero, then
// integer overflow flag is set and result is zero
template<class T1, class T2, 
        class Enable = typename mp::details::enable_mp_bin<T1,T2,void>::type,
        class Result = decltype(mp::details::div_impl(T1(),T2(),precision()))>
Result          operator/(const T1& x, const T2& y);

// plus function is equivalent to x+y
template<class T1, class T2, 
        class Enable = typename mp::details::enable_mp_bin<T1,T2,void>::type,
        class Result = decltype(mp::details::plus_impl(T1(),T2(),precision()))>
Result          plus(const T1& x, const T2& y, precision p = precision());

// minus function is equivalent to x-y
template<class T1, class T2, 
        class Enable = typename mp::details::enable_mp_bin<T1,T2,void>::type,
        class Result = decltype(mp::details::minus_impl(T1(),T2(),precision()))>
Result          minus(const T1& x, const T2& y, precision p = precision());

// mul function is equivalent to x*y
template<class T1, class T2, 
        class Enable = typename mp::details::enable_mp_bin<T1,T2,void>::type,
        class Result = decltype(mp::details::mul_impl(T1(),T2(),precision()))>
Result          mul(const T1& x, const T2& y, precision p = precision());

// div function is equivalent to x/y
template<class T1, class T2, 
        class Enable = typename mp::details::enable_mp_bin<T1,T2,void>::type,
        class Result = decltype(mp::details::div_impl(T1(),T2(),precision()))>
Result          div(const T1& x, const T2& y, precision p = precision());

// idiv function; if one of argument has floating point type or rational type, then
// result is equivalent to div function; however if both arguments are integers, then
// integer division is performed; if x and y are integer or rational values and y is
// zero, then integer overflow flag is set and result is zero
template<class T1, class T2, 
        class Enable = typename mp::details::enable_mp_bin<T1,T2,void>::type,
        class Result = decltype(mp::details::idiv_impl(T1(),T2(),precision()))>
Result          idiv(const T1& x, const T2& y, precision p = precision());

// division function is equivalent to x/y but 0/0 gives 0
template<class T1, class T2, 
        class Enable = typename mp::details::enable_mp_bin<T1,T2,void>::type,
        class Result = decltype(mp::details::div_impl(T1(),T2(),precision()))>
Result          div_0(const T1& x, const T2& y, precision p = precision());

// division function is equivalent to x/y but 0/0 gives 1
template<class T1, class T2, 
        class Enable = typename mp::details::enable_mp_bin<T1,T2,void>::type,
        class Result = decltype(mp::details::div_impl(T1(),T2(),precision()))>
Result          div_1(const T1& x, const T2& y, precision p = precision());

// calculate power x^y
template<class T1, class T2, 
        class Enable = typename mp::details::enable_mp_bin<T1,T2,void>::type,
        class Result = typename mp::details::result_of_pow<T1,T2>::type>
Result          pow(const T1& x, const T2& y, precision p = precision());

// calculate sqrt(|x|^2 + |y|^2) avoiding owerflow and underflow at
// intermediate stages of computation
template<class T1, class T2, 
        class Enable = typename mp::details::enable_mp_bin<T1,T2,void>::type>
mp_float        hypot(const T1& x, const T2& y, precision p = precision());

// four quadrant arctangent of x and y; -pi <= atan2(x, y) <= pi;
// not available for complex arguments
template<class T1, class T2, 
        class Enable = typename mp::details::enable_mp_bin<T1,T2,void>::type>
mp_float        atan2(const T1& x, const T2& y, precision p = precision());

// return the next representable value after x in the direction of y;
// return mp_float value with the same precision as the the first argument
// (possible after conversion to mp_float)
// not available for complex arguments
template<class T1, class T2, 
        class Enable = typename mp::details::enable_mp_bin<T1,T2,void>::type>
mp_float        nextafter(const T1& x, const T2& y);

// returns a value with the magnitude of x and the sign of y
// not available for complex arguments
template<class T1, class T2, 
        class Enable = typename mp::details::enable_mp_bin<T1,T2,void>::type>
mp_float        copysign(const T1& x, const T2& y);

// return is x - n * y where n = floor(x / y) if y != 0 and 0 otherwise;
// is mod(x,y) is nonzero, then has the same sign as y; if signs of x and y
// are the same, then mod(x,y) = rem(x,y); not available for complex arguments; 
// this function is accurate up to 1 ulp even for real arguments
template<class T1, class T2, 
        class Enable = typename mp::details::enable_mp_bin<T1,T2,void>::type,
        class Ret    = typename md::unify_types<T1,T2>::type>
Ret             mod(const T1& x, const T2& y, precision p = precision());

// return is x - n * y where n = trunc(x / y) if y != 0 and 0 otherwise
// is rem(x,y) is nonzero, then has the same sign as x; if signs of x and y
// are the same, then mod(x,y) = rem(x,y); not available for complex arguments
template<class T1, class T2, 
        class Enable = typename mp::details::enable_mp_bin<T1,T2,void>::type,
        class Ret    = typename md::unify_types<T1,T2>::type>
Ret             rem(const T1& x, const T2& y, precision p = precision());

// return max(x - y, 0) or NaN if x or y is NaN; not available for complex
// numbers
template<class T1, class T2, 
        class Enable = typename mp::details::enable_mp_bin<T1,T2,void>::type,
        class Ret    = typename matcl::details::unify_types2<T1,T2, Real>::type>
Ret             fdim(const T1& x, const T2& y, precision p = precision());

// compute a x b + c with one rounding at the end to precision p
MATCL_MP_EXPORT 
mp_float        fma(const mp_float& a, const mp_float& b, const mp_float& c,
                        precision p = precision());

// compute a x b - c with one rounding at the end to precision p
MATCL_MP_EXPORT 
mp_float        fms(const mp_float& a, const mp_float& b, const mp_float& c,
                        precision p = precision());

// compute a * b + c * d with high accuracy using Kahan's algorithm;
// this function is accurate up to 1 ulp
MATCL_MP_EXPORT
mp_float        dot2_ac(const mp_float& a, const mp_float& b, const mp_float& c,
                        const mp_float& d, precision p = precision());

};

#include "matcl-mp/details/func_binary.inl"