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

#include "test_simd_config.h"
#include "matcl-simd/simd.h"
#include "matcl-simd/simd_complex.h"

namespace matcl { namespace simd
{

template simd<float, 128, nosimd_tag>;
template simd<double, 128, nosimd_tag>;
template simd<int32_t, 128, nosimd_tag>;
template simd<int64_t, 128, nosimd_tag>;

template simd<float, 128, scalar_nosimd_tag>;
template simd<double, 128, scalar_nosimd_tag>;
template simd<int32_t, 128, scalar_nosimd_tag>;
template simd<int64_t, 128, scalar_nosimd_tag>;

template simd_compl<float, 128, nosimd_tag>;
template simd_compl<double, 128, nosimd_tag>;

template simd<float, 256, nosimd_tag>;
template simd<double, 256, nosimd_tag>;
template simd<int32_t, 256, nosimd_tag>;
template simd<int64_t, 256, nosimd_tag>;

template simd_compl<float, 256, nosimd_tag>;
template simd_compl<double, 256, nosimd_tag>;

template void simd<float, 128, nosimd_tag>::scatter<2>(float*) const;
template void simd<double, 128, nosimd_tag>::scatter<2>(double*) const;
template void simd<int32_t, 128, nosimd_tag>::scatter<2>(int32_t*) const;
template void simd<int64_t, 128, nosimd_tag>::scatter<2>(int64_t*) const;

template simd<float, 128, nosimd_tag> simd<float, 128, nosimd_tag>::select<3,2,1,0>() const;
template simd<double, 128, nosimd_tag> simd<double, 128, nosimd_tag>::select<1,1>() const;
template simd<int32_t, 128, nosimd_tag> simd<int32_t, 128, nosimd_tag>::select<3,2,1,0>() const;
template simd<int64_t, 128, nosimd_tag> simd<int64_t, 128, nosimd_tag>::select<1,1>() const;

template void simd<float, 128, scalar_nosimd_tag>::scatter<2>(float*) const;
template void simd<double, 128, scalar_nosimd_tag>::scatter<2>(double*) const;
template void simd<int32_t, 128, scalar_nosimd_tag>::scatter<2>(int32_t*) const;
template void simd<int64_t, 128, scalar_nosimd_tag>::scatter<2>(int64_t*) const;

template void simd<float, 256, nosimd_tag>::scatter<2>(float*) const;
template void simd<double, 256, nosimd_tag>::scatter<2>(double*) const;
template void simd<int32_t, 256, nosimd_tag>::scatter<2>(int32_t*) const;
template void simd<int64_t, 256, nosimd_tag>::scatter<2>(int64_t*) const;

template simd<float, 256, nosimd_tag> simd<float, 256, nosimd_tag>::select<3,2,1,0,0,1,2,3>() const;
template simd<double, 256, nosimd_tag> simd<double, 256, nosimd_tag>::select<3,2,1,0>() const;
template simd<int32_t, 256, nosimd_tag> simd<int32_t, 256, nosimd_tag>::select<3,2,1,0,0,1,2,3>() const;
template simd<int64_t, 256, nosimd_tag> simd<int64_t, 256, nosimd_tag>::select<3,2,1,0>() const;

template void simd_compl<float, 128, nosimd_tag>::scatter<2>(simd_single_complex*) const;
template void simd_compl<double, 128, nosimd_tag>::scatter<2>(simd_double_complex*) const;

template void simd_compl<float, 256, nosimd_tag>::scatter<2>(simd_single_complex*) const;
template void simd_compl<double, 256, nosimd_tag>::scatter<2>(simd_double_complex*) const;

#if MATCL_ARCHITECTURE_HAS_SSE2
    template simd<float, 128, sse_tag>;
    template simd<double, 128, sse_tag>;
    template simd<int32_t, 128, sse_tag>;
    template simd<int64_t, 128, sse_tag>;

    template simd<float, 128, scalar_sse_tag>;
    template simd<double, 128, scalar_sse_tag>;
    template simd<int32_t, 128, scalar_sse_tag>;
    template simd<int64_t, 128, scalar_sse_tag>;

    template simd_compl<float, 128, sse_tag>;
    template simd_compl<double, 128, sse_tag>;

    template simd<float, 256, sse_tag>;
    template simd<double, 256, sse_tag>;
    template simd<int32_t, 256, sse_tag>;
    template simd<int64_t, 256, sse_tag>;

    template simd_compl<float, 256, sse_tag>;
    template simd_compl<double, 256, sse_tag>;

    template void simd<float, 128, sse_tag>::scatter<2>(float*) const;
    template void simd<double, 128, sse_tag>::scatter<2>(double*) const;
    template void simd<int32_t, 128, sse_tag>::scatter<2>(int32_t*) const;
    template void simd<int64_t, 128, sse_tag>::scatter<2>(int64_t*) const;

    template simd<float, 128, sse_tag> simd<float, 128, sse_tag>::select<3,2,1,0>() const;
    template simd<double, 128, sse_tag> simd<double, 128, sse_tag>::select<1,1>() const;
    template simd<int32_t, 128, sse_tag> simd<int32_t, 128, sse_tag>::select<3,2,1,0>() const;
    template simd<int64_t, 128, sse_tag> simd<int64_t, 128, sse_tag>::select<1,1>() const;

    template void simd<float, 128, scalar_sse_tag>::scatter<2>(float*) const;
    template void simd<double, 128, scalar_sse_tag>::scatter<2>(double*) const;
    template void simd<int32_t, 128, scalar_sse_tag>::scatter<2>(int32_t*) const;
    template void simd<int64_t, 128, scalar_sse_tag>::scatter<2>(int64_t*) const;

    template void simd<float, 256, sse_tag>::scatter<2>(float*) const;
    template void simd<double, 256, sse_tag>::scatter<2>(double*) const;
    template void simd<int32_t, 256, sse_tag>::scatter<2>(int32_t*) const;
    template void simd<int64_t, 256, sse_tag>::scatter<2>(int64_t*) const;

    template simd<float, 256, sse_tag> simd<float, 256, sse_tag>::select<3,2,1,0,0,1,2,3>() const;
    template simd<double, 256, sse_tag> simd<double, 256, sse_tag>::select<3,2,1,0>() const;
    template simd<int32_t, 256, sse_tag> simd<int32_t, 256, sse_tag>::select<3,2,1,0,0,1,2,3>() const;
    template simd<int64_t, 256, sse_tag> simd<int64_t, 256, sse_tag>::select<3,2,1,0>() const;

    template void simd_compl<float, 128, sse_tag>::scatter<2>(simd_single_complex*) const;
    template void simd_compl<double, 128, sse_tag>::scatter<2>(simd_double_complex*) const;
    template void simd_compl<float, 256, sse_tag>::scatter<2>(simd_single_complex*) const;
    template void simd_compl<double, 256, sse_tag>::scatter<2>(simd_double_complex*) const;
#endif

#if MATCL_ARCHITECTURE_HAS_AVX
    template simd<float, 256, avx_tag>;
    template simd<double, 256, avx_tag>;
    template simd<int32_t, 256, avx_tag>;
    template simd<int64_t, 256, avx_tag>;

    template simd_compl<float, 256, avx_tag>;
    template simd_compl<double, 256, avx_tag>;

    template void simd<float, 256, avx_tag>::scatter<2>(float*) const;
    template void simd<double, 256, avx_tag>::scatter<2>(double*) const;
    template void simd<int32_t, 256, avx_tag>::scatter<2>(int32_t*) const;
    template void simd<int64_t, 256, avx_tag>::scatter<2>(int64_t*) const;

    template void simd_compl<float, 256, avx_tag>::scatter<2>(simd_single_complex*) const;
    template void simd_compl<double, 256, avx_tag>::scatter<2>(simd_double_complex*) const;

    template simd<float, 256, avx_tag> simd<float, 256, avx_tag>::select<3,2,1,0,0,1,2,3>() const;
    template simd<double, 256, avx_tag> simd<double, 256, avx_tag>::select<3,2,1,0>() const;
    template simd<int32_t, 256, avx_tag> simd<int32_t, 256, avx_tag>::select<3,2,1,0,0,1,2,3>() const;
    template simd<int64_t, 256, avx_tag> simd<int64_t, 256, avx_tag>::select<3,2,1,0>() const;
#endif

}}
