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

#include "matcl-simd/simd_general.h"

namespace matcl { namespace simd
{

// vector of elements in reverse order
template<class Val, int Bits, class Simd_tag>
simd<Val, Bits, Simd_tag> 
reverse(const simd<Val, Bits, Simd_tag>& x);

// vector multiply x * y
template<class Val, int Bits, class Simd_tag>
simd<Val, Bits, Simd_tag> 
operator*(const simd<Val, Bits, Simd_tag>& x, const simd<Val, Bits, Simd_tag>& y);

// vector division x / y
template<class Val, int Bits, class Simd_tag>
simd<Val, Bits, Simd_tag> 
operator/(const simd<Val, Bits, Simd_tag>& x, const simd<Val, Bits, Simd_tag>& y);

// vector add x + y
template<class Val, int Bits, class Simd_tag>
simd<Val, Bits, Simd_tag> 
operator+(const simd<Val, Bits, Simd_tag>& x, const simd<Val, Bits, Simd_tag>& y);

// vector subtract x - y
template<class Val, int Bits, class Simd_tag>
simd<Val, Bits, Simd_tag> 
operator-(const simd<Val, Bits, Simd_tag>& x, const simd<Val, Bits, Simd_tag>& y);

// vector unary minus x
template<class Val, int Bits, class Simd_tag>
simd<Val, Bits, Simd_tag> 
operator-(const simd<Val, Bits, Simd_tag>& x);

// absolute value
template<class Val, int Bits, class Simd_tag>
simd<Val, Bits, Simd_tag> 
abs(const simd<Val, Bits, Simd_tag>& x);

// alternatively subtract and add elements in x and y
// i.e. form [x[0] - y[0], x[1] + y[1], x[2] - y[2], x[3] + y[3], ...]
template<class Val, int Bits, class Simd_tag>
simd<Val, Bits, Simd_tag> 
sub_add(const simd<Val, Bits, Simd_tag>& x, const simd<Val, Bits, Simd_tag>& y);

//res = x*y + z; only one rounding at the end
template<class Val, int Bits, class Simd_tag>
simd<Val, Bits, Simd_tag> 
fma(const simd<Val, Bits, Simd_tag>& x, const simd<Val, Bits, Simd_tag>& y, 
                       const simd<Val, Bits, Simd_tag>& z);

//res = x*y - z; only one rounding at the end
template<class Val, int Bits, class Simd_tag>
simd<Val, Bits, Simd_tag> 
fms(const simd<Val, Bits, Simd_tag>& x, const simd<Val, Bits, Simd_tag>& y, 
                       const simd<Val, Bits, Simd_tag>& z);

// sum of all elements stored in the vector x
template<class Val, int Bits, class Simd_tag>
Val sum_all(const simd<Val, Bits, Simd_tag>& x);

// element by element maximum
template<class Val, int Bits, class Simd_tag>
simd<Val, Bits, Simd_tag>  
max(const simd<Val, Bits, Simd_tag>& x, const simd<Val, Bits, Simd_tag>& y);

// element by element minumum
template<class Val, int Bits, class Simd_tag>
simd<Val, Bits, Simd_tag>  
min(const simd<Val, Bits, Simd_tag>& x, const simd<Val, Bits, Simd_tag>& y);

// test x == y 
template<class Val, int Bits, class Simd_tag>
simd<Val, Bits, Simd_tag>  
eeq(const simd<Val, Bits, Simd_tag>& x, const simd<Val, Bits, Simd_tag>& y);

// test x != y 
template<class Val, int Bits, class Simd_tag>
simd<Val, Bits, Simd_tag>  
neq(const simd<Val, Bits, Simd_tag>& x, const simd<Val, Bits, Simd_tag>& y);

// test x < y 
template<class Val, int Bits, class Simd_tag>
simd<Val, Bits, Simd_tag>  
lt(const simd<Val, Bits, Simd_tag>& x, const simd<Val, Bits, Simd_tag>& y);

// test x == y 
template<class Val, int Bits, class Simd_tag>
simd<Val, Bits, Simd_tag>  
gt(const simd<Val, Bits, Simd_tag>& x, const simd<Val, Bits, Simd_tag>& y);

// test x <= y 
template<class Val, int Bits, class Simd_tag>
simd<Val, Bits, Simd_tag> 
leq(const simd<Val, Bits, Simd_tag>& x, const simd<Val, Bits, Simd_tag>& y);

// test x >= y 
template<class Val, int Bits, class Simd_tag>
simd<Val, Bits, Simd_tag>  
geq(const simd<Val, Bits, Simd_tag>& x, const simd<Val, Bits, Simd_tag>& y);

// square root
template<class Val, int Bits, class Simd_tag>
simd<Val, Bits, Simd_tag>  
sqrt(const simd<Val, Bits, Simd_tag>& x);

// round toward nearest integer
template<class Val, int Bits, class Simd_tag>
simd<Val, Bits, Simd_tag>  
round(const simd<Val, Bits, Simd_tag>& x);

// round toward -INF
template<class Val, int Bits, class Simd_tag>
simd<Val, Bits, Simd_tag>  
floor(const simd<Val, Bits, Simd_tag>& x);

// round toward +INF
template<class Val, int Bits, class Simd_tag>
simd<Val, Bits, Simd_tag>  
ceil(const simd<Val, Bits, Simd_tag>& x);

// round toward zero
template<class Val, int Bits, class Simd_tag>
simd<Val, Bits, Simd_tag>  
trunc(const simd<Val, Bits, Simd_tag>& x);

}}
