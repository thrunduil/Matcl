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

#include "matcl-core/details/simd/simd.h"
#include "matcl-simd/simd.h"

namespace matcl { namespace simd
{

#if MATCL_ARCHITECTURE_HAS_FMA

force_inline
double details::fma_f(double x, double y, double z)
{
    using simd_type = matcl::simd::simd<double,128,sse_tag>;

    simd_type xs(x);
    simd_type ys(y);
    simd_type zs(z);

    simd_type ret   = simd_type(_mm_fmadd_sd( xs.data, ys.data, zs.data));
    return ret.first();
};

force_inline
float details::fma_f(float x, float y, float z)
{
    using simd_type = matcl::simd::simd<float,128,sse_tag>;

    simd_type xs(x);
    simd_type ys(y);
    simd_type zs(z);

    simd_type ret   = simd_type(_mm_fmadd_ss( xs.data, ys.data, zs.data));
    return ret.first();
}

force_inline
double details::fms_f(double x, double y, double z)
{
    using simd_type = matcl::simd::simd<double,128,sse_tag>;

    simd_type xs(x);
    simd_type ys(y);
    simd_type zs(z);

    simd_type ret   = simd_type(_mm_fmsub_sd( xs.data, ys.data, zs.data));
    return ret.first();
}

force_inline
float details::fms_f(float x, float y, float z)
{
    using simd_type = matcl::simd::simd<float,128,sse_tag>;

    simd_type xs(x);
    simd_type ys(y);
    simd_type zs(z);

    simd_type ret   = simd_type(_mm_fmsub_ss( xs.data, ys.data, zs.data));
    return ret.first();
}

force_inline
double details::fnma_f(double x, double y, double z)
{
    using simd_type = matcl::simd::simd<double,128,sse_tag>;

    simd_type xs(x);
    simd_type ys(y);
    simd_type zs(z);

    simd_type ret   = simd_type(_mm_fnmadd_pd( xs.data, ys.data, zs.data));
    return ret.first();
};

force_inline
float details::fnma_f(float x, float y, float z)
{
    using simd_type = matcl::simd::simd<float,128,sse_tag>;

    simd_type xs(x);
    simd_type ys(y);
    simd_type zs(z);

    simd_type ret   = simd_type(_mm_fnmadd_ss( xs.data, ys.data, zs.data));
    return ret.first();
}

force_inline
double details::fnms_f(double x, double y, double z)
{
    using simd_type = matcl::simd::simd<double,128,sse_tag>;

    simd_type xs(x);
    simd_type ys(y);
    simd_type zs(z);

    simd_type ret   = simd_type(_mm_fnmsub_sd( xs.data, ys.data, zs.data));
    return ret.first();
}

force_inline
float details::fnms_f(float x, float y, float z)
{
    using simd_type = matcl::simd::simd<float,128,sse_tag>;

    simd_type xs(x);
    simd_type ys(y);
    simd_type zs(z);

    simd_type ret   = simd_type(_mm_fnmsub_ss( xs.data, ys.data, zs.data));
    return ret.first();
}

#endif

}}