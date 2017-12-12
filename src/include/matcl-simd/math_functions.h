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
// this function requires, that -127 < k < 128, where
// additionally pow2k(-127) = 0, pow2k(128) = inf, and for -127 < l < 128
// pow2k(l) is a regular number; no checks are performed
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

// return the exponential function, i.e. the e (Euler's number, 2.7182818)
// raised to the given power x; 
template<class Val, int Bits, class Simd_tag>
simd<Val, Bits, Simd_tag> 
exp(const simd<Val, Bits, Simd_tag>& x);

// return the exponential function, i.e. the e (Euler's number, 2.7182818)
// raised to the given power x, where x is a floating point scalar; 
double  exp(double x);
float   exp(float x);

}}
