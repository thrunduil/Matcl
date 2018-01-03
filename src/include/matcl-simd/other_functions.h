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

//-----------------------------------------------------------------------
//                   BIT MANIPULATION
//-----------------------------------------------------------------------

// count the number of leading zero bits in unsigned integer x
uint32_t number_leading_zeros(uint32_t x);
uint64_t number_leading_zeros(uint64_t x);

// count the number of trailing zero bits in unsigned integer x
uint32_t number_trailing_zeros(uint32_t x);
uint64_t number_trailing_zeros(uint64_t x);

// count the number of bits set to 1 in unsigned integer x
uint32_t number_bits_set(uint32_t x);
uint64_t number_bits_set(uint64_t x);

// extract the lowest set bit from unsigned integer x and set the corresponding bit 
// in the result res; all other bits in res are zeroed, and all bits are zeroed if no
// bits are set in x.
uint32_t least_significant_bit(uint32_t x);
uint64_t least_significant_bit(uint64_t x);

//-----------------------------------------------------------------------
//                   ARITHMETIC
//-----------------------------------------------------------------------

// return the high 64 bits of the 128-bit result of the multiplication
uint64_t mulh(uint64_t x, uint64_t y);

}}
