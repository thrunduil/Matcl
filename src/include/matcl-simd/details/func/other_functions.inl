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

#include "matcl-simd/other_functions.h"

#if MATCL_ARCHITECTURE_HAS_LZCNT || MATCL_ARCHITECTURE_HAS_BMI1
    #include "immintrin.h"
#endif

#if MATCL_ARCHITECTURE_HAS_POPCNT
    #include "nmmintrin.h"
#endif

#if MATCL_ARCHITECTURE_64 && defined (_MSC_VER)
    #include <intrin.h>  
#endif

namespace matcl
{

//-----------------------------------------------------------------------
//                   number_leading_zeros
//-----------------------------------------------------------------------

#if MATCL_ARCHITECTURE_HAS_LZCNT

    force_inline
    uint32_t simd::number_leading_zeros(uint32_t x)
    {
        return _lzcnt_u32(x);
    };

#else
    
    force_inline
    uint32_t simd::number_leading_zeros(uint32_t x)
    {
        //taken from The Aggregate Magic Algorithms
        x       |= (x >> 1);
        x       |= (x >> 2);
        x       |= (x >> 4);
        x       |= (x >> 8);
        x       |= (x >> 16);

        return 32 - number_bits_set(x);
    };

#endif

#if MATCL_ARCHITECTURE_HAS_LZCNT && MATCL_ARCHITECTURE_64

    force_inline
    uint64_t simd::number_leading_zeros(uint64_t x)
    {
        return _lzcnt_u64(x);
    };

#else
    
    force_inline
    uint64_t simd::number_leading_zeros(uint64_t x)
    {
	    //taken from The Aggregate Magic Algorithms
	    x |= (x >> 1);
	    x |= (x >> 2);
	    x |= (x >> 4);
	    x |= (x >> 8);
	    x |= (x >> 16);
	    x |= (x >> 32);

	    return 64 - number_bits_set(x);
    };

#endif

//-----------------------------------------------------------------------
//                   number_trailing_zeros
//-----------------------------------------------------------------------

#if MATCL_ARCHITECTURE_HAS_BMI1
    
    force_inline
    uint32_t simd::number_trailing_zeros(uint32_t x)
    {
        return _tzcnt_u32(x);
    };

#else

    force_inline
    uint32_t simd::number_trailing_zeros(uint32_t x)
    {
        //taken from The Aggregate Magic Algorithms

        int32_t xi = int32_t(x);
        return number_bits_set(uint32_t((xi & -xi) - 1));
    };

#endif

#if MATCL_ARCHITECTURE_HAS_BMI1 && MATCL_ARCHITECTURE_64
    
    force_inline
    uint64_t simd::number_trailing_zeros(uint64_t x)
    {
        return _tzcnt_u64(x);
    };

#else

    force_inline
    uint64_t simd::number_trailing_zeros(uint64_t x)
    {
        //taken from The Aggregate Magic Algorithms

        int64_t xi = int64_t(x);
        return number_bits_set(uint64_t((xi & -xi) - 1));
    };

#endif

//-----------------------------------------------------------------------
//                   least_significant_bit
//-----------------------------------------------------------------------

#if MATCL_ARCHITECTURE_HAS_BMI1
    
    force_inline
    uint32_t simd::least_significant_bit(uint32_t x)
    {
        return _blsi_u32(x);
    };

#else

    #pragma warning(push)
    #pragma warning(disable : 4146) //unary minus operator applied to unsigned type

    force_inline
    uint32_t simd::least_significant_bit(uint32_t x)
    {
        //taken from The Aggregate Magic Algorithms
        return x & -x;
    };

    #pragma warning(pop)
#endif

#if MATCL_ARCHITECTURE_HAS_BMI1 && MATCL_ARCHITECTURE_64
    
    force_inline
    uint64_t simd::least_significant_bit(uint64_t x)
    {
        return _blsi_u64(x);
    }

#else

    #pragma warning(push)
    #pragma warning(disable : 4146) //unary minus operator applied to unsigned type

    force_inline
    uint64_t simd::least_significant_bit(uint64_t x)
    {
        //taken from The Aggregate Magic Algorithms
        return x & -x;
    }

    #pragma warning(pop)
#endif

//-----------------------------------------------------------------------
//                   number_bits_set
//-----------------------------------------------------------------------

#if MATCL_ARCHITECTURE_HAS_POPCNT

    force_inline
    uint32_t simd::number_bits_set(uint32_t x)
    {
        return _mm_popcnt_u32(x);
    };

#else

    force_inline
    uint32_t simd::number_bits_set(uint32_t bits)
    {
        //taken from The Aggregate Magic Algorithms
        bits = bits - ((bits >> size_t(1)) & 0x55555555);
        bits = (bits & 0x33333333) + ((bits >> 2) & 0x33333333);
        bits = (((bits + (bits >> 4)) & 0x0F0F0F0F) * 0x01010101) >> 24;

        return bits;
    };

#endif

#if MATCL_ARCHITECTURE_HAS_POPCNT && MATCL_ARCHITECTURE_64

    force_inline
    uint64_t simd::number_bits_set(uint64_t x)
    {
        return _mm_popcnt_u64(x);
    };

#else

    force_inline
    uint64_t simd::number_bits_set(uint64_t x)
    {
        // Hamming weight algorithm with fast multiplication
        // http://en.wikipedia.org/wiki/Hamming_weight

	    static const uint64_t m1  = 0x5555555555555555; //binary: 0101...
	    static const uint64_t m2  = 0x3333333333333333; //binary: 00110011..
	    static const uint64_t m4  = 0x0f0f0f0f0f0f0f0f; //binary:  4 zeros,  4 ones ...
	    static const uint64_t h01 = 0x0101010101010101; //the sum of 256 to the power of 0,1,2,3...

	    x -= (x >> 1) & m1;             //put count of each 2 bits into those 2 bits
	    x = (x & m2) + ((x >> 2) & m2); //put count of each 4 bits into those 4 bits 
	    x = (x + (x >> 4)) & m4;        //put count of each 8 bits into those 8 bits 
	    return (x * h01) >> 56;         //returns left 8 bits of x + (x<<8) + (x<<16) + (x<<24) + ... 
    };

#endif

//-----------------------------------------------------------------------
//                   mulh
//-----------------------------------------------------------------------

force_inline
uint64_t simd::mulh(uint64_t x, uint64_t y)
{
    #if MATCL_ARCHITECTURE_64 && defined (_MSC_VER)

        return __umulh(x, y);

    #elif MATCL_ARCHITECTURE_HAS_BMI2 && MATCL_ARCHITECTURE_64

        uint64_t lo, hi;
        lo = _mulx_u64 (x, y, &hi);
        return hi;

    #else
    
        static const
        uint64_t mask_lo    = uint32_t(-1);

        uint64_t L1         = (uint32_t)x;
        uint64_t H1         = x >> 32;
        uint64_t L2         = (uint32_t)y;
        uint64_t H2         = y >> 32;

        uint64_t HH         = H1 * H2;
        uint64_t HL         = H1 * L2;
        uint64_t LH         = L1 * H2;
        uint64_t LL         = L1 * L2;

        uint64_t HL_L       = HL & mask_lo;
        uint64_t LH_L       = LH & mask_lo;
        uint64_t LL_H       = LL >> 32;

        uint64_t carry      = (HL_L + LH_L + LL_H) >> 32;
        uint64_t res        = HH + (HL >> 32) + (LH >> 32) + carry;

        return res;
    #endif
};

}
