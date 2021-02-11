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

#include "matcl-simd/arch/nosimd/simd_double_256.h"

#pragma warning(push)
#pragma warning(disable : 4127) // conditional expression is constant

namespace matcl { namespace simd
{

//-------------------------------------------------------------------
//                          GENERAL DOUBLE
//-------------------------------------------------------------------

force_inline
simd<double, 256, nosimd_tag>::simd(int32_t val)
    : simd(double(val)) 
{}

force_inline
simd<double, 256, nosimd_tag>::simd(int64_t val)
    : simd(double(val)) 
{}

force_inline
simd<double, 256, nosimd_tag>::simd(float val)
    : simd(double(val))
{}

force_inline
simd<double, 256, nosimd_tag>::simd(double val)
{
    data[0] = val;
    data[1] = val;
    data[2] = val;
    data[3] = val;
}

force_inline
simd<double, 256, nosimd_tag>::simd(double v1, double v2, double v3, double v4) 
{
    data[0] = v1;
    data[1] = v2;
    data[2] = v3;
    data[3] = v4;
}

force_inline
simd<double, 256, nosimd_tag>::simd(const simd_half& lo, const simd_half& hi)
{
    data[0] = lo.data[0];
    data[1] = lo.data[1];
    data[2] = hi.data[0];
    data[3] = hi.data[1];
}

force_inline
simd<double, 256, nosimd_tag>::simd(const simd_half& lo_hi)
{
    data[0] = lo_hi.data[0];
    data[1] = lo_hi.data[1];
    data[2] = lo_hi.data[0];
    data[3] = lo_hi.data[1];
}

#if MATCL_ARCHITECTURE_HAS_SSE2
    force_inline
    simd<double, 256, nosimd_tag>::simd(const simd<double, 256, sse_tag>& s)
    {
        s.store(data, std::true_type());
    }

    force_inline
    simd<double, 256, nosimd_tag>::simd(const simd<double, 128, scalar_sse_tag>& s)
        :simd(s.first())
    {}
#endif

#if MATCL_ARCHITECTURE_HAS_AVX
    force_inline
    simd<double, 256, nosimd_tag>::simd(const simd<double, 256, avx_tag>& s)
    {
        s.store(data, std::true_type());
    }
#endif

force_inline
simd<double, 256, nosimd_tag>::simd(const simd<double, 128, scalar_nosimd_tag>& s)
    :simd(s.first())
{}

force_inline
double simd<double, 256, nosimd_tag>::get(int pos) const
{ 
    return data[pos] ; 
};

force_inline
double simd<double, 256, nosimd_tag>::first() const
{ 
    return data[0]; 
};

force_inline
void simd<double, 256, nosimd_tag>::set(int pos, double val)
{ 
    data[pos] = val; 
};

force_inline
const double* simd<double, 256, nosimd_tag>::get_raw_ptr() const
{ 
    return data; 
};

force_inline
double* simd<double, 256, nosimd_tag>::get_raw_ptr()
{ 
    return data; 
};

force_inline simd<double, 256, nosimd_tag>::simd_half
simd<double, 256, nosimd_tag>::extract_low() const
{
    return simd_half(data[0], data[1]);
}

force_inline simd<double, 256, nosimd_tag>::simd_half
simd<double, 256, nosimd_tag>::extract_high() const
{
    return simd_half(data[2], data[3]);
}

template<int I1, int I2, int I3, int I4>
force_inline simd<double, 256, nosimd_tag>
simd<double, 256, nosimd_tag>::select() const
{
    static_assert(I1 >= 0 && I1 <= 3, "invalid index in select function");
    static_assert(I2 >= 0 && I2 <= 3, "invalid index in select function");
    static_assert(I3 >= 0 && I3 <= 3, "invalid index in select function");
    static_assert(I4 >= 0 && I4 <= 3, "invalid index in select function");

    if (I1 == 0 && I2 == 1 && I3 == 2 && I4 == 3) 
        return *this;

    simd<double, 256, nosimd_tag> ret;
    ret.data[0] = data[I1];
    ret.data[1] = data[I2];
    ret.data[2] = data[I3];
    ret.data[3] = data[I4];

    return ret;
};

force_inline
simd<double, 256, nosimd_tag> simd<double, 256, nosimd_tag>::zero()
{
    simd<double, 256, nosimd_tag> ret;
    ret.data[0] = 0.0;
    ret.data[1] = 0.0;
    ret.data[2] = 0.0;
    ret.data[3] = 0.0;

    return ret;
}

force_inline
simd<double, 256, nosimd_tag> simd<double, 256, nosimd_tag>::minus_zero()
{
    simd<double, 256, nosimd_tag> ret;
    ret.data[0] = -0.0;
    ret.data[1] = -0.0;
    ret.data[2] = -0.0;
    ret.data[3] = -0.0;

    return ret;
}

force_inline
simd<double, 256, nosimd_tag> simd<double, 256, nosimd_tag>::one()
{
    simd<double, 256, nosimd_tag> ret;
    ret.data[0] = 1.0;
    ret.data[1] = 1.0;
    ret.data[2] = 1.0;
    ret.data[3] = 1.0;

    return ret;
}

force_inline
simd<double, 256, nosimd_tag> simd<double, 256, nosimd_tag>::minus_one()
{
    simd<double, 256, nosimd_tag> ret;
    ret.data[0] = -1.0;
    ret.data[1] = -1.0;
    ret.data[2] = -1.0;
    ret.data[3] = -1.0;

    return ret;
}

force_inline simd<double, 256, nosimd_tag> 
simd<double, 256, nosimd_tag>::load(const double* arr, std::true_type aligned)
{
    (void)aligned;
    simd<double, 256, nosimd_tag> ret;
    ret.data[0] = arr[0];
    ret.data[1] = arr[1];
    ret.data[2] = arr[2];
    ret.data[3] = arr[3];

    return ret;
};

force_inline simd<double, 256, nosimd_tag> 
simd<double, 256, nosimd_tag>::load(const double* arr, std::false_type not_aligned)
{
    (void)not_aligned;
    simd<double, 256, nosimd_tag> ret;
    ret.data[0] = arr[0];
    ret.data[1] = arr[1];
    ret.data[2] = arr[2];
    ret.data[3] = arr[3];

    return ret;
};

force_inline simd<double, 256, nosimd_tag> 
simd<double, 256, nosimd_tag>::gather(const double* arr, const simd_int32_half& ind)
{
    simd<double, 256, nosimd_tag> ret;
    ret.data[0] = arr[ind.data[0]];
    ret.data[1] = arr[ind.data[1]];
    ret.data[2] = arr[ind.data[2]];
    ret.data[3] = arr[ind.data[3]];

    return ret;
}

force_inline simd<double, 256, nosimd_tag> 
simd<double, 256, nosimd_tag>::gather(const double* arr, const simd_int32& ind)
{
    simd<double, 256, nosimd_tag> ret;
    ret.data[0] = arr[ind.data[0]];
    ret.data[1] = arr[ind.data[1]];
    ret.data[2] = arr[ind.data[2]];
    ret.data[3] = arr[ind.data[3]];

    return ret;
}

force_inline simd<double, 256, nosimd_tag> 
simd<double, 256, nosimd_tag>::gather(const double* arr, const simd_int64& ind)
{
    simd<double, 256, nosimd_tag> ret;
    ret.data[0] = arr[ind.data[0]];
    ret.data[1] = arr[ind.data[1]];
    ret.data[2] = arr[ind.data[2]];
    ret.data[3] = arr[ind.data[3]];

    return ret;
}

force_inline simd<double, 256, nosimd_tag> 
simd<double, 256, nosimd_tag>::broadcast(const double* arr)
{
    simd<double, 256, nosimd_tag> ret;
    ret.data[0] = arr[0];
    ret.data[1] = arr[0];
    ret.data[2] = arr[0];
    ret.data[3] = arr[0];

    return ret;
};

force_inline simd<double, 256, nosimd_tag> 
simd<double, 256, nosimd_tag>::broadcast(const double& arr)
{
    simd<double, 256, nosimd_tag> ret;
    ret.data[0] = arr;
    ret.data[1] = arr;
    ret.data[2] = arr;
    ret.data[3] = arr;

    return ret;
};

force_inline void 
simd<double, 256, nosimd_tag>::store(double* arr, std::true_type aligned) const
{
    (void)aligned;
    arr[0]  = data[0];
    arr[1]  = data[1];
    arr[2]  = data[2];
    arr[3]  = data[3];
};

force_inline void 
simd<double, 256, nosimd_tag>::store(double* arr, std::false_type not_aligned) const
{
    (void)not_aligned;
    arr[0]  = data[0];
    arr[1]  = data[1];
    arr[2]  = data[2];
    arr[3]  = data[3];
};

force_inline
simd<float, 128, nosimd_tag>
simd<double, 256, nosimd_tag>::convert_to_float() const
{
    simd<float, 128, nosimd_tag> res;

    res.data[0] = scalar_func::convert_double_float(data[0]);
    res.data[1] = scalar_func::convert_double_float(data[1]);
    res.data[2] = scalar_func::convert_double_float(data[2]);
    res.data[3] = scalar_func::convert_double_float(data[3]);

    return res;
};

force_inline
simd<int32_t, 128, nosimd_tag>
simd<double, 256, nosimd_tag>::convert_to_int32() const
{
    simd<int32_t, 128, nosimd_tag> res;

    res.data[0] = scalar_func::convert_double_int32(data[0]);
    res.data[1] = scalar_func::convert_double_int32(data[1]);
    res.data[2] = scalar_func::convert_double_int32(data[2]);
    res.data[3] = scalar_func::convert_double_int32(data[3]);

    return res;
};

force_inline simd<int64_t, 256, nosimd_tag> 
simd<double, 256, nosimd_tag>::convert_to_int64() const
{
    simd<int64_t, 256, nosimd_tag> res;

    res.data[0] = scalar_func::convert_double_int64(data[0]);
    res.data[1] = scalar_func::convert_double_int64(data[1]);
    res.data[2] = scalar_func::convert_double_int64(data[2]);
    res.data[3] = scalar_func::convert_double_int64(data[3]);

    return res;
};

force_inline simd<int64_t, 256, nosimd_tag>
simd<double, 256, nosimd_tag>::reinterpret_as_int64() const
{
    return *reinterpret_cast<const simd_int64*>(this);
}

force_inline simd<float, 256, nosimd_tag>
simd<double, 256, nosimd_tag>::reinterpret_as_float() const
{
    return *reinterpret_cast<const simd_float*>(this);
}

force_inline simd<int32_t, 256, nosimd_tag>
simd<double, 256, nosimd_tag>::reinterpret_as_int32() const
{
    return *reinterpret_cast<const simd_int32*>(this);
}

template<int Step>
force_inline
void simd<double, 256, nosimd_tag>::scatter(double* arr) const
{
    arr[0*Step] = data[0];
    arr[1*Step] = data[1];
    arr[2*Step] = data[2];
    arr[3*Step] = data[3];
};

force_inline simd<double, 256, nosimd_tag>& 
simd<double, 256, nosimd_tag>::operator+=(const simd& x)
{
    *this = *this + x;
    return *this;
}

force_inline simd<double, 256, nosimd_tag>& 
simd<double, 256, nosimd_tag>::operator-=(const simd& x)
{
    *this = *this - x;
    return *this;
}

force_inline simd<double, 256, nosimd_tag>& 
simd<double, 256, nosimd_tag>::operator*=(const simd& x)
{
    *this = *this * x;
    return *this;
}

force_inline simd<double, 256, nosimd_tag>& 
simd<double, 256, nosimd_tag>::operator/=(const simd& x)
{
    *this = *this / x;
    return *this;
}

}}

#pragma warning(pop)