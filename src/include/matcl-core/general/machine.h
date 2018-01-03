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

//--------------------------------------------------------------------------
//                      Machine dependent parameters
//--------------------------------------------------------------------------
// cache line size in bytes; this is optimization parameters, invalid value
// may have negative impact on performance
#ifndef MATCL_CACHE_LINE_SIZE
    #define MATCL_CACHE_LINE_SIZE   64
#endif

// set value of this macro to 1 if SSE2 instruction set is available and 0
// otherwise
#ifndef MATCL_ARCHITECTURE_HAS_SSE2
    #define MATCL_ARCHITECTURE_HAS_SSE2 1
#endif

// set value of this macro to 1 if SSE3 instruction set is available and 0
// otherwise
#ifndef MATCL_ARCHITECTURE_HAS_SSE3
    #define MATCL_ARCHITECTURE_HAS_SSE3 1
#endif

// set value of this macro to 1 if SSE4.1 instruction set is available and 0
// otherwise
#ifndef MATCL_ARCHITECTURE_HAS_SSE41
    #define MATCL_ARCHITECTURE_HAS_SSE41 1
#endif

// set value of this macro to 1 if SSE4.2 instruction set is available and 0
// otherwise
#ifndef MATCL_ARCHITECTURE_HAS_SSE42
    #define MATCL_ARCHITECTURE_HAS_SSE42 1
#endif

// set value of this macro to 1 if AVX instruction set is available and 0
// otherwise
#ifndef MATCL_ARCHITECTURE_HAS_AVX
    #define MATCL_ARCHITECTURE_HAS_AVX 1
#endif

// set value of this macro to 1 if AVX2 instruction set is available and 0
// otherwise
#ifndef MATCL_ARCHITECTURE_HAS_AVX2
    #define MATCL_ARCHITECTURE_HAS_AVX2 1
#endif

// set value of this macro to 1 if FMA instruction set is available and 0
// otherwise
#ifndef MATCL_ARCHITECTURE_HAS_FMA
    #define MATCL_ARCHITECTURE_HAS_FMA 1
#endif

// set value of this macro to 1 if POPCNT instruction is available and 0
// otherwise
#ifndef MATCL_ARCHITECTURE_HAS_POPCNT
    #define MATCL_ARCHITECTURE_HAS_POPCNT 1
#endif

// set value of this macro to 1 if LZCNT instruction is available and 0
// otherwise
#ifndef MATCL_ARCHITECTURE_HAS_LZCNT
    #define MATCL_ARCHITECTURE_HAS_LZCNT 1
#endif

// set value of this macro to 1 if BMI1 instruction set is available and 0
// otherwise
#ifndef MATCL_ARCHITECTURE_HAS_BMI1
    #define MATCL_ARCHITECTURE_HAS_BMI1 1
#endif

// set value of this macro to 1 if BMI2 instruction set is available and 0
// otherwise
#ifndef MATCL_ARCHITECTURE_HAS_BMI2
    #define MATCL_ARCHITECTURE_HAS_BMI2 1
#endif

// alignment for SIMD related types; this is optimization parameters, invalid
// value may have negative impact on performance
#ifndef MATCL_SIMD_ALIGNMENT

    #if MATCL_ARCHITECTURE_HAS_AVX
        #define MATCL_SIMD_ALIGNMENT    32
    #elif MATCL_ARCHITECTURE_HAS_SSE2
        #define MATCL_SIMD_ALIGNMENT    16
    #else
        #define MATCL_SIMD_ALIGNMENT    8
    #endif

#endif

// set value of this macro to 1 if SIMD registers are available
#ifndef MATCL_ARCHITECTURE_HAS_SIMD

    #if MATCL_ARCHITECTURE_HAS_AVX
        #define MATCL_ARCHITECTURE_HAS_SIMD 1
    #elif MATCL_ARCHITECTURE_HAS_SSE2
        #define MATCL_ARCHITECTURE_HAS_SIMD 1
    #else
        #define MATCL_ARCHITECTURE_HAS_SIMD 0
    #endif

#endif

// set MATCL_ARCHITECTURE_32 to 1 if compiling on 32-bit architecture
// set MATCL_ARCHITECTURE_64 to 1 if compiling on 64-bit architecture
#if defined (_M_X64)
    #define MATCL_ARCHITECTURE_32 0
    #define MATCL_ARCHITECTURE_64 1
#else
    #define MATCL_ARCHITECTURE_32 1
    #define MATCL_ARCHITECTURE_64 0
#endif
