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

#include "matcl-core/float/fma_dekker.h"
#include "matcl-core/float/twofold.h"
#include "matcl-simd/simd.h"

namespace matcl { namespace details
{

template<class T>
struct nan_recovery
{};

template<>
struct nan_recovery<float>
{
    static bool test(float x)
    {
        namespace mrds  = matcl::raw::details::scal_func;
        return mrds::isnan(x);
    }

    static float eval(float res, float x, float y, float z)
    {
        (void)res;
        return x * y + z;
    }
};

template<>
struct nan_recovery<double>
{
    static bool test(double x)
    {
        namespace mrds  = matcl::raw::details::scal_func;
        return mrds::isnan(x);
    }

    static double eval(double res, double x, double y, double z)
    {
        (void)res;
        return x * y + z;
    }
};

template<class T, int Bits, class Tag>
struct nan_recovery<simd::simd<T, Bits, Tag>>
{
    using simd_type = simd::simd<T, Bits, Tag>;

    static bool test(const simd_type& x)
    {
        return simd::any_nan(x);
    }

    static simd_type eval(const simd_type& res, const simd_type& x, 
                          const simd_type& y, const simd_type& z)
    {
        using nosimd_type = simd::simd<T, Bits, simd::nosimd_tag>;

        nosimd_type x0(x);
        nosimd_type y0(y);
        nosimd_type z0(z);

        (void)res;
        return simd_type(simd::fma_f(x0, y0, z0));
    }
};

template<class T>
force_inline
T fma_dekker_double_impl(const T& x, const T& y, const T& z)
{
    T xy_val, xy_err;

    // form x * y = xy_val + xy_err exactly
    twofold_mult_dekker_double(x, y, xy_val, xy_err);

    // form xy_val + z = s_val + s_err
    twofold<T> s   = twofold_plus(xy_val, z);

    T res       = xy_err + s.error;
    res         = res + s.value;

    if (nan_recovery<T>::test(res) == false)
        return res;
    else
        return nan_recovery<T>::eval(res, x, y, z);
};

}}

namespace matcl
{

double matcl::fma_dekker(double x, double y, double z)
{
    return details::fma_dekker_double_impl(x, y, z);
};

#if MATCL_ARCHITECTURE_HAS_SSE2

    MATCL_CORE_EXPORT
    simd::simd<double, 128, simd::sse_tag>
    matcl::fma_dekker_simd(const simd::simd<double, 128, simd::sse_tag>& x, 
                    const simd::simd<double, 128, simd::sse_tag>& y, 
                    const simd::simd<double, 128, simd::sse_tag>& z)
    {
        return details::fma_dekker_double_impl(x, y, z);
    };

#endif

#if MATCL_ARCHITECTURE_HAS_AVX

    MATCL_CORE_EXPORT
    simd::simd<double, 256, simd::avx_tag>
    matcl::fma_dekker_simd(const simd::simd<double, 256, simd::avx_tag>& x, 
                    const simd::simd<double, 256, simd::avx_tag>& y, 
                    const simd::simd<double, 256, simd::avx_tag>& z)
    {
        return details::fma_dekker_double_impl(x, y, z);
    };

#endif

}