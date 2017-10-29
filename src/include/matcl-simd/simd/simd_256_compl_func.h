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

#include "matcl-simd/arch/reg_256/simd_256_compl.h"

namespace matcl { namespace simd
{

// vector multiply x * y
template<class Val> simd_compl<Val,reg_256>
operator*(const simd_compl<Val,reg_256>& x, const simd_compl<Val,reg_256>& y);

// vector multiply x * y, where x is a real vector
template<class Val> simd_compl<Val,reg_256>
operator*(const simd<Val,reg_256>& x, const simd_compl<Val,reg_256>& y);

// vector multiply x * y, where y is a real vector
template<class Val> simd_compl<Val,reg_256>
operator*(const simd_compl<Val,reg_256>& x, const simd<Val,reg_256>& y);

// vector division x / y
template<class Val>
simd_compl<Val,reg_256>
operator/(const simd_compl<Val,reg_256>& x, const simd_compl<Val,reg_256>& y);

// vector add x + y
template<class Val>
simd_compl<Val,reg_256>
operator+(const simd_compl<Val,reg_256>& x, const simd_compl<Val,reg_256>& y);

// vector subtract x - y
template<class Val>
simd_compl<Val,reg_256>
operator-(const simd_compl<Val,reg_256>& x, const simd_compl<Val,reg_256>& y);

// vector unary minus x
template<class Val>
simd_compl<Val,reg_256>
operator-(const simd_compl<Val,reg_256>& x);

// vector of elements in reverse order
template<class Val>
simd_compl<Val,reg_256> reverse(const simd_compl<Val,reg_256>& x);

//res = x*y + z; rounding different, than required by the fma function
template<class Val>
simd_compl<Val,reg_256>
fma(const simd_compl<Val,reg_256>& x, const simd_compl<Val,reg_256>& y, 
                         const simd_compl<Val,reg_256>& z);
template<class Val>
simd_compl<Val,reg_256>
fma(const simd<Val,reg_256>& x, const simd_compl<Val,reg_256>& y, 
                         const simd_compl<Val,reg_256>& z);

//res = x*y - z; rounding different, than required by the fms function
template<class Val>
simd_compl<Val,reg_256>
fms(const simd_compl<Val,reg_256>& x, const simd_compl<Val,reg_256>& y, 
                         const simd_compl<Val,reg_256>& z);
template<class Val>
simd_compl<Val,reg_256>
fms(const simd<Val,reg_256>& x, const simd_compl<Val,reg_256>& y, 
                         const simd_compl<Val,reg_256>& z);

// sum of all elements stored in the vector x
template<class Val>
typename simd_compl<Val,reg_256>::value_type
sum_all(const simd_compl<Val,reg_256>& x);

}}

