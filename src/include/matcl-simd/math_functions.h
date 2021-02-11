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

#include "matcl-simd/simd.h"
#include "matcl-simd/details/utils.h"

namespace matcl { namespace simd
{

//-----------------------------------------------------------------------
//                   MATHEMATICAL FUNCTIONS
//-----------------------------------------------------------------------

// return a vector of floating point values containing 2^k; k must represent an
// integer and min_k < k < max_k, where
//      for single precision: min_k = -127,  max_k = 128
//      for double precision: min_k = -1023, max_k = 1024
// additionally pow2k(min_k) = 0, pow2k(max_k) = inf, and for min_k < l < max_k
// pow2k(l) is a regular number; no checks are performed
template<class Val, int Bits, class Simd_tag>
simd<Val, Bits, Simd_tag> 
pow2k(const simd<Val, Bits, Simd_tag>& k);

// return a vector of floats containing 2^k; k is a vector of int32_t values; 
// this function requires, that -127 < k < 128, additionally pow2k(-127) = 0, 
// pow2k(128) = inf, and for -127 < l < 128 pow2k(l) is a regular number; 
// no checks are performed
template<int Bits, class Simd_tag>
simd<float, Bits, Simd_tag> 
pow2ki(const simd<int32_t, Bits, Simd_tag>& k);

// return a vector of floats containing 2^k; k is a vector of int64_t values; 
// this function requires, that -1023 < k < 1024, where
// additionally pow2k(-1023) = 0, pow2k(1024) = inf, and for -1023 < l < 1024
// pow2k(l) is a regular number; no checks are performed
template<int Bits, class Simd_tag>
simd<double, Bits, Simd_tag> 
pow2ki(const simd<int64_t, Bits, Simd_tag>& k);

// return the fraction part f of a floating point number x: 
//      x = f * 2^exp, 0.5 <= |f| < 1 
// this functions returns correct valuest only for regular numbers (i.e. not 0, 
// +- inf, NaN and denormal numbers) and for irregular numbers:
//      +- inf  -> +- 0.5
//      +- 0    -> +- 0.5
//      other   -> v,       where 0.5 <= |v| < 1 
//  
// no checks are performed
template<class Val, int Bits, class Simd_tag>
simd<Val, Bits, Simd_tag> 
fraction(const simd<Val, Bits, Simd_tag>& x);

// return the exponent part exp of a floating point number x represented as integer: 
//      x = f * 2^exp, 0.5 <= |f| < 1 
// this functions returns correct valuest only for regular numbers (i.e. not 0, 
// +- inf, NaN and denormal numbers)
// no checks are performed
template<class Val, int Bits, class Simd_tag>
simd<typename details::integer_type<Val>::type, Bits, Simd_tag> 
iexponent(const simd<Val, Bits, Simd_tag>& x);

// return the exponent part exp of a floating point number x represented as floating
// point number: 
//      x = f * 2^exp, 0.5 <= |f| < 1 
// this functions returns correct valuest only for regular numbers (i.e. not 0, 
// +- inf, NaN and denormal numbers) and for irregular numbers:
//      inf, nan     -> (1025,   129)    (double exponent, float exponent)
//      0, denormals -> (-1022, -126)
// for regular values: 
//      double precision: exp in [-1021, 1024]
//      single precision: exp in [-125, 128]
// no checks are performed
template<class Val, int Bits, class Simd_tag>
simd<Val, Bits, Simd_tag> 
exponent(const simd<Val, Bits, Simd_tag>& x);

// composes a floating point value with the magnitude of x and the sign of y
template<class Val, int Bits, class Simd_tag>
simd<Val, Bits, Simd_tag> 
copysign(const simd<Val, Bits, Simd_tag>& x, const simd<Val, Bits, Simd_tag>& y);

// return the exponential function, i.e. the e (Euler's number, 2.7182818)
// raised to the given power x; 
template<class Val, int Bits, class Simd_tag>
simd<Val, Bits, Simd_tag> 
exp(const simd<Val, Bits, Simd_tag>& x);

// computes the natural (base e, Euler's number, 2.7182818) logarithm of x
template<class Val, int Bits, class Simd_tag>
simd<Val, Bits, Simd_tag> 
log(const simd<Val, Bits, Simd_tag>& x);

// return the sine function of x in radians
template<class Val, int Bits, class Simd_tag>
simd<Val, Bits, Simd_tag> 
sin(const simd<Val, Bits, Simd_tag>& x);

// return the cosine function of x in radians
template<class Val, int Bits, class Simd_tag>
simd<Val, Bits, Simd_tag> 
cos(const simd<Val, Bits, Simd_tag>& x);

// return the tangent function of x in radians
template<class Val, int Bits, class Simd_tag>
simd<Val, Bits, Simd_tag> 
tan(const simd<Val, Bits, Simd_tag>& x);

// return the cotangent function of x in radians
template<class Val, int Bits, class Simd_tag>
simd<Val, Bits, Simd_tag> 
cot(const simd<Val, Bits, Simd_tag>& x);

}}
