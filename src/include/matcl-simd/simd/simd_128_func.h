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

#include "matcl-simd/arch/reg_128/simd_128.h"

namespace matcl { namespace simd
{

// vector of elements in reverse order
template<class Val>
simd<Val,reg_128> reverse(const simd<Val,reg_128>& x);

// vector multiply x * y
template<class Val>
simd<Val,reg_128> operator*(const simd<Val,reg_128>& x, const simd<Val,reg_128>& y);

// vector division x / y
template<class Val>
simd<Val,reg_128> operator/(const simd<Val,reg_128>& x, const simd<Val,reg_128>& y);

// vector add x + y
template<class Val>
simd<Val,reg_128> operator+(const simd<Val,reg_128>& x, const simd<Val,reg_128>& y);

// vector subtract x - y
template<class Val>
simd<Val,reg_128> operator-(const simd<Val,reg_128>& x, const simd<Val,reg_128>& y);

// vector unary minus x
template<class Val>
simd<Val,reg_128> operator-(const simd<Val,reg_128>& x);

// vector unary minus x
template<class Val>
simd<Val,reg_128> operator-(const simd<Val,reg_128>& x);

// absolute value
template<class Val>
simd<Val,reg_128> abs(const simd<Val,reg_128>& x);

//res = x*y + z; only one rounding at the end
template<class Val>
simd<Val,reg_128> fma(const simd<Val,reg_128>& x, const simd<Val,reg_128>& y, 
                       const simd<Val,reg_128>& z);

//res = x*y - z; only one rounding at the end
template<class Val>
simd<Val,reg_128> fms(const simd<Val,reg_128>& x, const simd<Val,reg_128>& y, 
                       const simd<Val,reg_128>& z);

// sum of all elements stored in the vector x
template<class Val>
Val sum_all(const simd<Val,reg_128>& x);

// element by element maximum
template<class Val>
simd<Val,reg_128>  max(const simd<Val,reg_128>& x, const simd<Val,reg_128>& y);

// element by element minumum
template<class Val>
simd<Val,reg_128>  min(const simd<Val,reg_128>& x, const simd<Val,reg_128>& y);

// test x == y 
template<class Val>
simd<Val,reg_128>  eeq(const simd<Val,reg_128>& x, const simd<Val,reg_128>& y);

// test x != y 
template<class Val>
simd<Val,reg_128>  neq(const simd<Val,reg_128>& x, const simd<Val,reg_128>& y);

// test x < y 
template<class Val>
simd<Val,reg_128>  lt(const simd<Val,reg_128>& x, const simd<Val,reg_128>& y);

// test x == y 
template<class Val>
simd<Val,reg_128>  gt(const simd<Val,reg_128>& x, const simd<Val,reg_128>& y);

// test x <= y 
template<class Val>
simd<Val,reg_128>  leq(const simd<Val,reg_128>& x, const simd<Val,reg_128>& y);

// test x >= y 
template<class Val>
simd<Val,reg_128>  geq(const simd<Val,reg_128>& x, const simd<Val,reg_128>& y);

// square root
template<class Val>
simd<Val,reg_128>  sqrt(const simd<Val,reg_128>& x);

// round toward nearest integer
template<class Val>
simd<Val,reg_128>  round(const simd<Val,reg_128>& x);

// round toward -INF
template<class Val>
simd<Val,reg_128>  floor(const simd<Val,reg_128>& x);

// round toward +INF
template<class Val>
simd<Val,reg_128>  ceil(const simd<Val,reg_128>& x);

// round toward zero
template<class Val>
simd<Val,reg_128>  trunc(const simd<Val,reg_128>& x);

}}
