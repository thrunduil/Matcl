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

#include "matcl-core/config.h"
#include "matcl-simd/simd_fwd.h"

namespace matcl
{

// implementation of fma function using dekker algorithm

// evaluate a * b using the Veltkamp/Dekker algorithm; result is exact; 
// this function is equivalent to twofold_mult and is called by twofold_mult
// and other functions, when FMA instruction is not available
void twofold_mult_dekker(double a, double b, double& val, double& err);

// evaluate x * y + z with one rounding using the Veltkamp/Dekker algorithm
// for double precision arguments and double arithmetics for single precision
// arguments; 
// this function is equivalent to FMA instruction and is called by twofold's
// function, when FMA instruction is not available
float   fma_dekker(float x, float y, float z);

MATCL_CORE_EXPORT
double  fma_dekker(double x, double y, double z);

// evaluate x * y + z with one rounding using the Veltkamp/Dekker algorithm
// for simd types
#if MATCL_ARCHITECTURE_HAS_SSE2

    simd::simd<float, 128, simd::sse_tag>
    fma_dekker_simd(const simd::simd<float, 128, simd::sse_tag>& x, 
                    const simd::simd<float, 128, simd::sse_tag>& y, 
                    const simd::simd<float, 128, simd::sse_tag>& z);

    MATCL_CORE_EXPORT
    simd::simd<double, 128, simd::sse_tag>
    fma_dekker_simd(const simd::simd<double, 128, simd::sse_tag>& x, 
                    const simd::simd<double, 128, simd::sse_tag>& y, 
                    const simd::simd<double, 128, simd::sse_tag>& z);
#endif

#if MATCL_ARCHITECTURE_HAS_AVX

    simd::simd<float, 256, simd::avx_tag>
    fma_dekker_simd(const simd::simd<float, 256, simd::avx_tag>& x, 
                    const simd::simd<float, 256, simd::avx_tag>& y, 
                    const simd::simd<float, 256, simd::avx_tag>& z);

    MATCL_CORE_EXPORT
    simd::simd<double, 256, simd::avx_tag>
    fma_dekker_simd(const simd::simd<double, 256, simd::avx_tag>& x, 
                    const simd::simd<double, 256, simd::avx_tag>& y, 
                    const simd::simd<double, 256, simd::avx_tag>& z);
#endif

}

#include "matcl-core/details/float/fma_dekker.inl"

