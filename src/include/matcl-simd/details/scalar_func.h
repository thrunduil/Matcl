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
#include <cmath>

namespace matcl { namespace simd { namespace scalar_func
{

template<class T> struct integer_type{};
template<>        struct integer_type<float>    { using type = int32_t;};
template<>        struct integer_type<double>   { using type = int64_t;};

template<class T>
struct is_signed_impl
{
    using int_type  = typename integer_type<T>::type;

    static bool eval(const T& x)
    {
        int_type xi     = *reinterpret_cast<const int_type*>(&x);
        return xi < 0;
    }
};

template<class T>
struct round_impl
{
    force_inline
    static T eval(T x)
    {
        // std::round has wrong tie rule (we require round to even in case
        // of ties); std::nearbyint is very slow;

        const T einv = details::eps_inv<T>();
        
        T xa        = std::abs(x);
        T tmp       = xa + einv;
        tmp         = tmp - einv;
        T resa      = (xa < einv) ? tmp : xa;
        bool sign   = is_signed_impl<T>::eval(x);
        return sign ? -resa : resa;
    };
};

// implementation of missing simd functions

force_inline float floor(float f)     { return std::floor(f); };
force_inline float ceil(float f)      { return std::ceil(f); };
force_inline float trunc(float f)     { return std::trunc(f); };
force_inline float round(float f)     { return round_impl<float>::eval(f); };

force_inline double floor(double f)   { return std::floor(f); };
force_inline double ceil(double f)    { return std::ceil(f); };
force_inline double trunc(double f)   { return std::trunc(f); };
force_inline double round(double f)   { return round_impl<double>::eval(f); };

force_inline bool is_nan(double f)   { return f != f; }
force_inline bool is_nan(float f)    { return f != f; }

force_inline float fma_f(float x, float y, float z)
{
    return x * y + z;
};

force_inline float fms_f(float x, float y, float z)
{
    return x * y - z;
};

force_inline float fnma_f(float x, float y, float z)
{
    return z - x * y;
};

force_inline float fnms_f(float x, float y, float z)
{
    return -(x * y + z);
};

force_inline double fma_f(double x, double y, double z)
{
    return x * y + z;
};

force_inline double fms_f(double x, double y, double z)
{
    return x * y - z;
};

force_inline double fnma_f(double x, double y, double z)
{
    return z - x * y;
};

force_inline double fnms_f(double x, double y, double z)
{
    return -(x * y + z);
};

force_inline float fma_a(float x, float y, float z)
{
    return fma_dekker(x, y, z);
}

force_inline float fms_a(float x, float y, float z)
{
    return fma_dekker(x, y, -z);
}

force_inline float fnma_a(float x, float y, float z)
{
    return fma_dekker(-x, y, z);
}

force_inline float fnms_a(float x, float y, float z)
{
    return -fma_dekker(x, y, z);
}

force_inline double fma_a(double x, double y, double z)
{
    return fma_dekker(x, y, z);
}

force_inline double fms_a(double x, double y, double z)
{
    return fma_dekker(x, y, -z);
}

force_inline double fnma_a(double x, double y, double z)
{
    return fma_dekker(-x, y, z);
}

force_inline double fnms_a(double x, double y, double z)
{
    return -fma_dekker(x, y, z);
}

}}}
