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

force_inline double
ms::log(double x)
{
    using simd_type     = default_scalar_simd_type<double>::type;

    simd_type xs        = simd_type(x);
    simd_type res       = ms::log(xs);

    return res.first();
};

force_inline float
ms::log(float x)
{
    using simd_type     = default_scalar_simd_type<float>::type;

    simd_type xs        = simd_type(x);
    simd_type res       = ms::log(xs);

    return res.first();
};

force_inline double ms::pow2k(double x)
{
    // do not use SIMD-based scalars
    using simd_type = simd<double, 128, scalar_nosimd_tag>;
    simd_type res   = pow2k(simd_type(x));
    return res.first();
}

force_inline float ms::pow2k(float x)
{
    // do not use SIMD-based scalars
    using simd_type = simd<float, 128, scalar_nosimd_tag>;
    simd_type res   = pow2k(simd_type(x));
    return res.first();
}

force_inline double ms::fraction(double x)
{
    // do not use SIMD-based scalars
    using simd_type = simd<double, 128, scalar_nosimd_tag>;
    simd_type res   = fraction(simd_type(x));
    return res.first();
}

force_inline float ms::fraction(float x)
{
    // do not use SIMD-based scalars
    using simd_type = simd<float, 128, scalar_nosimd_tag>;
    simd_type res   = fraction(simd_type(x));
    return res.first();
}

force_inline double ms::exponent(double x)
{
    // do not use SIMD-based scalars
    using simd_type = simd<double, 128, scalar_nosimd_tag>;
    simd_type res   = exponent(simd_type(x));
    return res.first();
}

force_inline float ms::exponent(float x)
{
    // do not use SIMD-based scalars
    using simd_type = simd<float, 128, scalar_nosimd_tag>;
    simd_type res   = exponent(simd_type(x));
    return res.first();
}

force_inline int32_t ms::iexponent(float x)
{
    // do not use SIMD-based scalars
    using simd_type = simd<float, 128, scalar_nosimd_tag>;
    auto res        = iexponent(simd_type(x));
    return res.first();
}

force_inline int64_t ms::iexponent(double x)
{
    // do not use SIMD-based scalars
    using simd_type = simd<double, 128, scalar_nosimd_tag>;
    auto res        = iexponent(simd_type(x));
    return res.first();
}

force_inline float ms::copysign(float x, float y)
{
    using simd_type     = default_scalar_simd_type<float>::type;
    simd_type res       = copysign(simd_type(x), simd_type(y));
    return res.first();
}

force_inline double ms::copysign(double x, double y)
{
    using simd_type     = default_scalar_simd_type<double>::type;
    simd_type res       = copysign(simd_type(x), simd_type(y));
    return res.first();
}

force_inline float ms::pow2ki(int32_t k)
{
    // do not use SIMD-based scalars
    using simd_type = simd<int32_t, 128, scalar_nosimd_tag>;
    auto res        = pow2ki(simd_type(k));
    return res.first();
}

force_inline double ms::pow2ki(int64_t k)
{
    // do not use SIMD-based scalars
    using simd_type = simd<int64_t, 128, scalar_nosimd_tag>;
    auto res        = pow2ki(simd_type(k));
    return res.first();
};

force_inline double
ms::sin(double x)
{
    using simd_type     = default_scalar_simd_type<double>::type;

    simd_type xs        = simd_type(x);
    simd_type res       = ms::sin(xs);

    return res.first();
};

force_inline float
ms::sin(float x)
{
    using simd_type     = default_scalar_simd_type<float>::type;

    simd_type xs        = simd_type(x);
    simd_type res       = ms::sin(xs);

    return res.first();
};

force_inline double
ms::cos(double x)
{
    using simd_type     = default_scalar_simd_type<double>::type;

    simd_type xs        = simd_type(x);
    simd_type res       = ms::cos(xs);

    return res.first();
};

force_inline float
ms::cos(float x)
{
    using simd_type     = default_scalar_simd_type<float>::type;

    simd_type xs        = simd_type(x);
    simd_type res       = ms::cos(xs);

    return res.first();
};

force_inline double
ms::tan(double x)
{
    using simd_type     = default_scalar_simd_type<double>::type;

    simd_type xs        = simd_type(x);
    simd_type res       = ms::tan(xs);

    return res.first();
};

force_inline float
ms::tan(float x)
{
    using simd_type     = default_scalar_simd_type<float>::type;

    simd_type xs        = simd_type(x);
    simd_type res       = ms::tan(xs);

    return res.first();
};

force_inline double
ms::cot(double x)
{
    using simd_type     = default_scalar_simd_type<double>::type;

    simd_type xs        = simd_type(x);
    simd_type res       = ms::cot(xs);

    return res.first();
};

force_inline float
ms::cot(float x)
{
    using simd_type     = default_scalar_simd_type<float>::type;

    simd_type xs        = simd_type(x);
    simd_type res       = ms::cot(xs);

    return res.first();
};

}}
