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

#include "matcl-simd/default_simd.h"
#include "matcl-simd/math_functions.h"

namespace matcl { namespace simd
{

//-----------------------------------------------------------------------
//                   MATHEMATICAL FUNCTIONS
//-----------------------------------------------------------------------
force_inline double
ms::exp(double x)
{
    using simd_type     = default_scalar_simd_type<double>::type;

    simd_type xs        = simd_type(x);
    simd_type res       = ms::exp(xs);

    return res.first();
};

force_inline float
ms::exp(float x)
{
    using simd_type     = default_scalar_simd_type<float>::type;

    simd_type xs        = simd_type(x);
    simd_type res       = ms::exp(xs);

    return res.first();
};

}}
