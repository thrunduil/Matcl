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

#ifdef MATCL_TEST_SIMD_FMA
    #define MATCL_ARCHITECTURE_HAS_SSE2 1
    #define MATCL_ARCHITECTURE_HAS_SSE3 1
    #define MATCL_ARCHITECTURE_HAS_SSE41 1
    #define MATCL_ARCHITECTURE_HAS_SSE42 1
    #define MATCL_ARCHITECTURE_HAS_AVX 1
    #define MATCL_ARCHITECTURE_HAS_AVX2 1
    #define MATCL_ARCHITECTURE_HAS_FMA 1

    #define MATCL_TEST_SIMD_TAG "fma"

#elif defined MATCL_TEST_SIMD_AVX2
    #define MATCL_ARCHITECTURE_HAS_SSE2 1
    #define MATCL_ARCHITECTURE_HAS_SSE3 1
    #define MATCL_ARCHITECTURE_HAS_SSE41 1
    #define MATCL_ARCHITECTURE_HAS_SSE42 1
    #define MATCL_ARCHITECTURE_HAS_AVX 1
    #define MATCL_ARCHITECTURE_HAS_AVX2 1
    #define MATCL_ARCHITECTURE_HAS_FMA 0

    #define MATCL_TEST_SIMD_TAG "avx2"

#elif defined MATCL_TEST_SIMD_AVX
    #define MATCL_ARCHITECTURE_HAS_SSE2 1
    #define MATCL_ARCHITECTURE_HAS_SSE3 1
    #define MATCL_ARCHITECTURE_HAS_SSE41 1
    #define MATCL_ARCHITECTURE_HAS_SSE42 1
    #define MATCL_ARCHITECTURE_HAS_AVX 1
    #define MATCL_ARCHITECTURE_HAS_AVX2 0
    #define MATCL_ARCHITECTURE_HAS_FMA 0

    #define MATCL_TEST_SIMD_TAG "avx"

#elif defined MATCL_TEST_SIMD_NO_SSE42
    #define MATCL_ARCHITECTURE_HAS_SSE2 1
    #define MATCL_ARCHITECTURE_HAS_SSE3 1
    #define MATCL_ARCHITECTURE_HAS_SSE41 1
    #define MATCL_ARCHITECTURE_HAS_SSE42 0
    #define MATCL_ARCHITECTURE_HAS_AVX 0
    #define MATCL_ARCHITECTURE_HAS_AVX2 0
    #define MATCL_ARCHITECTURE_HAS_FMA 0

    #define MATCL_TEST_SIMD_TAG "no_sse42"

#elif defined MATCL_TEST_SIMD_NO_SSE41
    #define MATCL_ARCHITECTURE_HAS_SSE2 1
    #define MATCL_ARCHITECTURE_HAS_SSE3 1
    #define MATCL_ARCHITECTURE_HAS_SSE41 0
    #define MATCL_ARCHITECTURE_HAS_SSE42 0
    #define MATCL_ARCHITECTURE_HAS_AVX 0
    #define MATCL_ARCHITECTURE_HAS_AVX2 0
    #define MATCL_ARCHITECTURE_HAS_FMA 0

    #define MATCL_TEST_SIMD_TAG "no_sse41"

#elif defined MATCL_TEST_SIMD_NO_SSE3
    #define MATCL_ARCHITECTURE_HAS_SSE2 1
    #define MATCL_ARCHITECTURE_HAS_SSE3 0
    #define MATCL_ARCHITECTURE_HAS_SSE41 0
    #define MATCL_ARCHITECTURE_HAS_SSE42 0
    #define MATCL_ARCHITECTURE_HAS_AVX 0
    #define MATCL_ARCHITECTURE_HAS_AVX2 0
    #define MATCL_ARCHITECTURE_HAS_FMA 0

    #define MATCL_TEST_SIMD_TAG "no_sse3"

#elif defined MATCL_TEST_SIMD_DEFAULT

    #define MATCL_TEST_SIMD_TAG "default"

#else
    #error unknown simd type
#endif