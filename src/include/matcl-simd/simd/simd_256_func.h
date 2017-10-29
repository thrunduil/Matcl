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

#include "matcl-simd/arch/reg_256/simd_256.h"

namespace matcl { namespace simd
{

// vector multiply x * y
template<class Val>
simd<Val,reg_256>  operator*(const simd<Val,reg_256>& x, const simd<Val,reg_256>& y);

// vector division x / y
template<class Val>
simd<Val,reg_256>  operator/(const simd<Val,reg_256>& x, const simd<Val,reg_256>& y);

// vector add x + y
template<class Val>
simd<Val,reg_256>  operator+(const simd<Val,reg_256>& x, const simd<Val,reg_256>& y);

// vector subtract x - y
template<class Val>
simd<Val,reg_256>  operator-(const simd<Val,reg_256>& x, const simd<Val,reg_256>& y);

// vector unary minus x
template<class Val>
simd<Val,reg_256>  operator-(const simd<Val,reg_256>& x);

// absolute value
template<class Val>
simd<Val,reg_256> abs(const simd<Val,reg_256>& x);

// vector of elements in reverse order
template<class Val>
simd<Val,reg_256>  reverse(const simd<Val,reg_256>& x);

// compute x*y + z; only one rounding at the end
template<class Val>
simd<Val,reg_256>  fma(const simd<Val,reg_256>& x, const simd<Val,reg_256>& y, 
                        const simd<Val,reg_256>& z);
// compute x*y - z; only one rounding at the end
template<class Val>
simd<Val,reg_256>  fms(const simd<Val,reg_256>& x, const simd<Val,reg_256>& y, 
                        const simd<Val,reg_256>& z);

// sum of all elements stored in the vector x
template<class Val>
Val sum_all(const simd<Val,reg_256>& x);

// element by element maximum
template<class Val>
simd<Val,reg_256>  max(const simd<Val,reg_256>& x, const simd<Val,reg_256>& y);

// element by element minumum
template<class Val>
simd<Val,reg_256>  min(const simd<Val,reg_256>& x, const simd<Val,reg_256>& y);

// test x == y 
template<class Val>
simd<Val,reg_256>  eeq(const simd<Val,reg_256>& x, const simd<Val,reg_256>& y);

// test x != y 
template<class Val>
simd<Val,reg_256>  neq(const simd<Val,reg_256>& x, const simd<Val,reg_256>& y);

// test x < y 
template<class Val>
simd<Val,reg_256>  lt(const simd<Val,reg_256>& x, const simd<Val,reg_256>& y);

// test x == y 
template<class Val>
simd<Val,reg_256>  gt(const simd<Val,reg_256>& x, const simd<Val,reg_256>& y);

// test x <= y 
template<class Val>
simd<Val,reg_256>  leq(const simd<Val,reg_256>& x, const simd<Val,reg_256>& y);

// test x >= y 
template<class Val>
simd<Val,reg_256>  geq(const simd<Val,reg_256>& x, const simd<Val,reg_256>& y);

// square root
template<class Val>
simd<Val,reg_256>  sqrt(const simd<Val,reg_256>& x);

// round toward nearest integer
template<class Val>
simd<Val,reg_256>  round(const simd<Val,reg_256>& x);

// round toward -INF
template<class Val>
simd<Val,reg_256>  floor(const simd<Val,reg_256>& x);

// round toward +INF
template<class Val>
simd<Val,reg_256>  ceil(const simd<Val,reg_256>& x);

// round toward zero
template<class Val>
simd<Val,reg_256>  trunc(const simd<Val,reg_256>& x);

}}
