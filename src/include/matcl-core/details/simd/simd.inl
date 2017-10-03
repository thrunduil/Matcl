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

#include "matcl-core/details/simd/simd.h"
#include <immintrin.h>

namespace matcl
{

inline double simd::fma(const double& x, const double& y, const double& z)
{
    __m128d xs  = _mm_set1_pd(x);
    __m128d ys  = _mm_set1_pd(y);
    __m128d zs  = _mm_set1_pd(z);

    __m128d ret = _mm_fmadd_pd(xs, ys, zs);

    return ret.m128d_f64[0];
};

inline float simd::fma(const float& x, const float& y, const float& z)
{
    __m128 xs   = _mm_set1_ps(x);
    __m128 ys   = _mm_set1_ps(y);
    __m128 zs   = _mm_set1_ps(z);

    __m128 ret  = _mm_fmadd_ps(xs, ys, zs);

    return ret.m128_f32[0];
}

inline double simd::fms(const double& x, const double& y, const double& z)
{
    __m128d xs  = _mm_set1_pd(x);
    __m128d ys  = _mm_set1_pd(y);
    __m128d zs  = _mm_set1_pd(z);

    __m128d ret = _mm_fmsub_pd(xs, ys, zs);

    return ret.m128d_f64[0];
}

inline float  simd::fms(const float& x, const float& y, const float& z)
{
    __m128 xs   = _mm_set1_ps(x);
    __m128 ys   = _mm_set1_ps(y);
    __m128 zs   = _mm_set1_ps(z);

    __m128 ret  = _mm_fmsub_ps(xs, ys, zs);

    return ret.m128_f32[0];
}

}