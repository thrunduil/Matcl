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
#include "matcl-simd/simd.h"

namespace matcl
{

force_inline
double simd::fma_f(const double& x, const double& y, const double& z)
{
    using simd_type = matcl::simd::default_simd_type<double>::type;

    simd_type xs(x);
    simd_type ys(y);
    simd_type zs(z);

    simd_type ret = fma_f(xs, ys, zs);
    return ret.first();
};

force_inline
float simd::fma_f(const float& x, const float& y, const float& z)
{
    using simd_type = matcl::simd::default_simd_type<float>::type;

    simd_type xs(x);
    simd_type ys(y);
    simd_type zs(z);

    simd_type ret = fma_f(xs, ys, zs);
    return ret.first();
}

force_inline
double simd::fms_f(const double& x, const double& y, const double& z)
{
    using simd_type = matcl::simd::default_simd_type<double>::type;

    simd_type xs(x);
    simd_type ys(y);
    simd_type zs(z);

    simd_type ret = fms_f(xs, ys, zs);
    return ret.first();
}

force_inline
float simd::fms_f(const float& x, const float& y, const float& z)
{
    using simd_type = matcl::simd::default_simd_type<float>::type;

    simd_type xs(x);
    simd_type ys(y);
    simd_type zs(z);

    simd_type ret = fms_f(xs, ys, zs);
    return ret.first();
}

}