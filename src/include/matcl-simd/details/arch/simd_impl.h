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

#include "matcl-simd/config.h"

#include "matcl-simd/arch/nosimd/simd_float_128.h"
#include "matcl-simd/arch/nosimd/simd_double_128.h"
#include "matcl-simd/arch/nosimd/simd_int32_128.h"
#include "matcl-simd/arch/nosimd/simd_int64_128.h"

#include "matcl-simd/arch/nosimd/simd_scalar_float_128.h"
#include "matcl-simd/arch/nosimd/simd_scalar_double_128.h"
#include "matcl-simd/arch/nosimd/simd_scalar_int32_128.h"
#include "matcl-simd/arch/nosimd/simd_scalar_int64_128.h"

#include "matcl-simd/arch/nosimd/simd_float_256.h"
#include "matcl-simd/arch/nosimd/simd_double_256.h"
#include "matcl-simd/arch/nosimd/simd_int32_256.h"
#include "matcl-simd/arch/nosimd/simd_int64_256.h"

#if MATCL_ARCHITECTURE_HAS_SSE2
    #include "matcl-simd/arch/sse/simd_double_128.h"
    #include "matcl-simd/arch/sse/simd_float_128.h"
    #include "matcl-simd/arch/sse/simd_int32_128.h"
    #include "matcl-simd/arch/sse/simd_int64_128.h"
    
    #include "matcl-simd/arch/sse/simd_double_256.h"
    #include "matcl-simd/arch/sse/simd_float_256.h"
    #include "matcl-simd/arch/sse/simd_int32_256.h"
    #include "matcl-simd/arch/sse/simd_int64_256.h"

    #include "matcl-simd/arch/sse/simd_scalar_double_128.h"
    #include "matcl-simd/arch/sse/simd_scalar_float_128.h"
    #include "matcl-simd/arch/sse/simd_scalar_int32_128.h"
    #include "matcl-simd/arch/sse/simd_scalar_int64_128.h"
#endif

#if MATCL_ARCHITECTURE_HAS_AVX
    #include "matcl-simd/arch/avx/simd_double_256.h"
    #include "matcl-simd/arch/avx/simd_float_256.h"
    #include "matcl-simd/arch/avx/simd_int32_256.h"
    #include "matcl-simd/arch/avx/simd_int64_256.h"
#endif