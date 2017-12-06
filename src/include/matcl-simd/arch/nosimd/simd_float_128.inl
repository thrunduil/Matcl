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

#include "matcl-simd/arch/nosimd/simd_float_128.h"

namespace matcl { namespace simd
{

//-------------------------------------------------------------------
//                          GENERIC SINGLE
//-------------------------------------------------------------------

force_inline
simd<float, 128, nosimd_tag>::simd(float val)
{
    data[0] = val;
    data[1] = val;
    data[2] = val;
    data[3] = val;
}

force_inline
simd<float, 128, nosimd_tag>::simd(float v1, float v2, float v3, float v4) 
{
    data[0] = v1;
    data[1] = v2;
    data[2] = v3;
    data[3] = v4;
}

force_inline
simd<float, 128, nosimd_tag>::simd(const simd& lo, const simd& hi)
{
    data[0] = lo.data[0];
    data[1] = lo.data[1];
    data[2] = hi.data[0];
    data[3] = hi.data[1];
}

force_inline
simd<float, 128, nosimd_tag>::simd(const impl_type& v)
{
    data[0] = v[0];
    data[1] = v[1];
    data[2] = v[2];
    data[3] = v[3];
};

force_inline
simd<float, 128, nosimd_tag>::simd(const simd<float, 128, sse_tag>& s)
{
    s.store(data, std::true_type());
}

force_inline
simd<float, 128, nosimd_tag> simd<float, 128, nosimd_tag>::broadcast(const float* arr)
{
    simd<float, 128, nosimd_tag> ret;

    ret.data[0] = arr[0];
    ret.data[1] = arr[0];
    ret.data[2] = arr[0];
    ret.data[3] = arr[0];

    return ret;
};

force_inline
simd<float, 128, nosimd_tag> simd<float, 128, nosimd_tag>::broadcast(const float& arr)
{ 
    simd<float, 128, nosimd_tag> ret;

    ret.data[0] = arr;
    ret.data[1] = arr;
    ret.data[2] = arr;
    ret.data[3] = arr;

    return ret;
};

force_inline
float simd<float, 128, nosimd_tag>::get(int pos) const
{ 
    return data[pos] ; 
};

force_inline
float simd<float, 128, nosimd_tag>::first() const
{ 
    return data[0]; 
};

template<int Pos>
force_inline
float simd<float, 128, nosimd_tag>::get() const
{ 
    return data[Pos]; 
};

force_inline
void simd<float, 128, nosimd_tag>::set(int pos, float val)
{ 
    data[pos] = val; 
};

template<int Pos>
force_inline
void simd<float, 128, nosimd_tag>::set(float val)
{ 
    data[Pos] = val; 
};

force_inline
const float* simd<float, 128, nosimd_tag>::get_raw_ptr() const
{ 
    return data; 
};

force_inline
float* simd<float, 128, nosimd_tag>::get_raw_ptr()
{ 
    return data; 
};

force_inline
simd<float, 128, nosimd_tag> simd<float, 128, nosimd_tag>::zero()
{
    simd<float, 128, nosimd_tag> ret;

    ret.data[0] = 0.0f;
    ret.data[1] = 0.0f;
    ret.data[2] = 0.0f;
    ret.data[3] = 0.0f;

    return ret;
}

force_inline
simd<float, 128, nosimd_tag> simd<float, 128, nosimd_tag>::minus_zero()
{
    simd<float, 128, nosimd_tag> ret;

    ret.data[0] = -0.0f;
    ret.data[1] = -0.0f;
    ret.data[2] = -0.0f;
    ret.data[3] = -0.0f;

    return ret;
}

force_inline
simd<float, 128, nosimd_tag> simd<float, 128, nosimd_tag>::one()
{
    simd<float, 128, nosimd_tag> ret;

    ret.data[0] = 1.0f;
    ret.data[1] = 1.0f;
    ret.data[2] = 1.0f;
    ret.data[3] = 1.0f;

    return ret;
}

force_inline
simd<float, 128, nosimd_tag> simd<float, 128, nosimd_tag>::minus_one()
{
    simd<float, 128, nosimd_tag> ret;

    ret.data[0] = -1.0f;
    ret.data[1] = -1.0f;
    ret.data[2] = -1.0f;
    ret.data[3] = -1.0f;

    return ret;
}

force_inline simd<float, 128, nosimd_tag> 
simd<float, 128, nosimd_tag>::load(const float* arr, std::true_type aligned)
{
    (void)aligned;
    simd<float, 128, nosimd_tag> ret;

    ret.data[0] = arr[0];
    ret.data[1] = arr[1];
    ret.data[2] = arr[2];
    ret.data[3] = arr[3];

    return ret;
};

force_inline simd<float, 128, nosimd_tag> 
simd<float, 128, nosimd_tag>::load(const float* arr, std::false_type not_aligned)
{
    (void)not_aligned;
    simd<float, 128, nosimd_tag> ret;

    ret.data[0] = arr[0];
    ret.data[1] = arr[1];
    ret.data[2] = arr[2];
    ret.data[3] = arr[3];

    return ret;
};

force_inline simd<float, 128, nosimd_tag>
simd<float, 128, nosimd_tag>::gather(const float* arr, const simd_128_int32& ind)
{
    simd<float, 128, nosimd_tag> ret;

    ret.data[0] = arr[ind.data[0]];
    ret.data[1] = arr[ind.data[1]];
    ret.data[2] = arr[ind.data[2]];
    ret.data[3] = arr[ind.data[3]];

    return ret;
}

force_inline simd<float, 128, nosimd_tag> 
simd<float, 128, nosimd_tag>::gather(const float* arr, const simd_128_int64& ind)
{
    simd<float, 128, nosimd_tag> ret;

    ret.data[0] = arr[ind.data[0]];
    ret.data[1] = arr[ind.data[1]];
    ret.data[2] = ret.data[0];
    ret.data[3] = ret.data[1];

    return ret;
}

force_inline simd<float, 128, nosimd_tag> 
simd<float, 128, nosimd_tag>::set_lower(float v)
{
    simd<float, 128, nosimd_tag> ret;

    ret.data[0] = v;
    ret.data[1] = 0.0f;
    ret.data[2] = 0.0f;
    ret.data[3] = 0.0f;

    return ret;
};

force_inline simd<double, 128, nosimd_tag>
simd<float, 128, nosimd_tag>::convert_low_to_double() const
{
    simd<double, 128, nosimd_tag> ret;

    ret.data[0] = data[0];
    ret.data[1] = data[1];

    return ret;
};

force_inline simd<double, 128, nosimd_tag>
simd<float, 128, nosimd_tag>::convert_high_to_double() const
{
    simd<double, 128, nosimd_tag> ret;

    ret.data[0] = data[2];
    ret.data[1] = data[3];

    return ret;
};

force_inline
simd<double, 256, nosimd_tag>
simd<float, 128, nosimd_tag>::convert_to_double() const
{
    simd<double, 256, nosimd_tag> ret;

    ret.data[0] = data[0];
    ret.data[1] = data[1];
    ret.data[2] = data[2];
    ret.data[3] = data[3];

    return ret;
};

force_inline
simd<int32_t, 128, nosimd_tag>
simd<float, 128, nosimd_tag>::convert_to_int32() const
{
    simd<int32_t, 128, nosimd_tag> ret;

    ret.data[0] = int32_t(data[0]);
    ret.data[1] = int32_t(data[1]);
    ret.data[2] = int32_t(data[2]);
    ret.data[3] = int32_t(data[3]);

    return ret;
};

force_inline simd<int32_t, 128, nosimd_tag>
simd<float, 128, nosimd_tag>::reinterpret_as_int32() const
{
    return *reinterpret_cast<const simd_128_int32*>(this);
}

force_inline simd<int64_t, 128, nosimd_tag>
simd<float, 128, nosimd_tag>::reinterpret_as_int64() const
{
    return *reinterpret_cast<const simd_128_int64*>(this);
}

force_inline simd<double, 128, nosimd_tag>
simd<float, 128, nosimd_tag>::reinterpret_as_double() const
{
    return *reinterpret_cast<const simd_128_double*>(this);
}

force_inline simd<float, 128, nosimd_tag> 
simd<float, 128, nosimd_tag>::extract_low() const
{
    return *this;
}

force_inline simd<float, 128, nosimd_tag> 
simd<float, 128, nosimd_tag>::extract_high() const
{
    simd<float, 128, nosimd_tag> ret;
    ret.data[0] = data[2];
    ret.data[1] = data[3];
    ret.data[2] = data[2];
    ret.data[3] = data[3];

    return ret;
}

force_inline void 
simd<float, 128, nosimd_tag>::store(float* arr, std::true_type aligned) const
{
    (void)aligned;
    arr[0] = data[0];
    arr[1] = data[1];
    arr[2] = data[2];
    arr[3] = data[3];
};


force_inline void 
simd<float, 128, nosimd_tag>::store(float* arr, std::false_type not_aligned) const
{
    (void)not_aligned;
    arr[0] = data[0];
    arr[1] = data[1];
    arr[2] = data[2];
    arr[3] = data[3];
};

template<int Step>
force_inline
void simd<float, 128, nosimd_tag>::scatter(float* arr) const
{
    //no scatter intrinsic
    arr[0*Step] = data[0];
    arr[1*Step] = data[1];
    arr[2*Step] = data[2];
    arr[3*Step] = data[3];
};

force_inline simd<float, 128, nosimd_tag>& 
simd<float, 128, nosimd_tag>::operator+=(const simd& x)
{
    *this = *this + x;
    return *this;
}

force_inline simd<float, 128, nosimd_tag>& 
simd<float, 128, nosimd_tag>::operator-=(const simd& x)
{
    *this = *this - x;
    return *this;
}

force_inline simd<float, 128, nosimd_tag>& 
simd<float, 128, nosimd_tag>::operator*=(const simd& x)
{
    *this = *this * x;
    return *this;
}

force_inline simd<float, 128, nosimd_tag>& 
simd<float, 128, nosimd_tag>::operator/=(const simd& x)
{
    *this = *this / x;
    return *this;
}

}}
