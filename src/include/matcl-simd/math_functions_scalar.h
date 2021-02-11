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
//              MATHEMATICAL FUNCTIONS FOR SCALAR ARGUMENTS
//-----------------------------------------------------------------------

// return 2^k; k must represent an integer and min_k < k < max_k, where
//      for single precision: min_k = -127,  max_k = 128
//      for double precision: min_k = -1023, max_k = 1024
// additionally pow2k(min_k) = 0, pow2k(max_k) = inf, and for min_k < l < max_k
// pow2k(l) is a regular number; no checks are performed
double  pow2k(double x);
float   pow2k(float x);

// return 2^k; k is a int32_t value;  this function requires, that -127 < k < 128, 
// additionally pow2k(-127) = 0, pow2k(128) = inf, and for -127 < l < 128
// pow2k(l) is a regular number; no checks are performed
float   pow2ki(int32_t k);

// return 2^k; k is a int64_t value;  this function requires, that -1023 < k < 1024,
// additionally pow2k(-1023) = 0, pow2k(1024) = inf, and for -1023 < l < 1024
// pow2k(l) is a regular number; no checks are performed
double  pow2ki(int64_t k);

// return the fraction part f of a floating point number x: 
//      x = f * 2^exp, 0.5 <= |f| < 1 
// this functions returns correct valuest only for regular numbers (i.e. not 0, 
// +- inf, NaN and denormal numbers) and for irregular numbers:
//      +- inf  -> +- 0.5
//      +- 0    -> +- 0.5
//      other   -> v,       where 0.5 <= |v| < 1 
//  
// no checks are performed
double  fraction(double x);
float   fraction(float x);

// return the exponent part exp of a floating point number x represented as integer: 
//      x = f * 2^exp, 0.5 <= |f| < 1 
// this functions returns correct valuest only for regular numbers (i.e. not 0, 
// +- inf, NaN and denormal numbers)
// no checks are performed
int32_t iexponent(float x);
int64_t iexponent(double x);

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
float   exponent(float x);
double  exponent(double x);

// composes a floating point value with the magnitude of x and the sign of y
float   copysign(float x, float y);
double  copysign(double x, double y);

// return the exponential function, i.e. the e (Euler's number, 2.7182818)
// raised to the given power x; 
double  exp(double x);
float   exp(float x);

// computes the natural (base e, Euler's number, 2.7182818) logarithm of x
double  log(double x);
float   log(float x);

// return the sine function of x in radians
double  sin(double x);
float   sin(float x);

// return the cosine function of x in radians
double  cos(double x);
float   cos(float x);

// return the tangent function of x in radians
double  tan(double x);
float   tan(float x);

// return the cotangent function of x in radians
double  cot(double x);
float   cot(float x);

}}
