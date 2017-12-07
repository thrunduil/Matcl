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

#include "matcl-simd/arch/nosimd/simd_float_256.h"

namespace matcl { namespace simd
{

//-------------------------------------------------------------------
//                          GENERAL SINGLE
//-------------------------------------------------------------------

force_inline
simd<float, 256, nosimd_tag>::simd(float val) 
{
    data[0] = val;
    data[1] = val;
    data[2] = val;
    data[3] = val;
    data[4] = val;
    data[5] = val;
    data[6] = val;
    data[7] = val;
}

force_inline
simd<float, 256, nosimd_tag>::simd(const simd_half& lo_hi)
{
    data[0] = lo_hi.data[0];
    data[1] = lo_hi.data[1];
    data[2] = lo_hi.data[2];
    data[3] = lo_hi.data[3];
    data[4] = lo_hi.data[0];
    data[5] = lo_hi.data[1];
    data[6] = lo_hi.data[2];
    data[7] = lo_hi.data[3];
}

force_inline
simd<float, 256, nosimd_tag>::simd(const simd_half& lo, const simd_half& hi)
{
    data[0] = lo.data[0];
    data[1] = lo.data[1];
    data[2] = lo.data[2];
    data[3] = lo.data[3];
    data[4] = hi.data[0];
    data[5] = hi.data[1];
    data[6] = hi.data[2];
    data[7] = hi.data[3];
}

force_inline
simd<float, 256, nosimd_tag>::simd(float v1, float v2, float v3, float v4, float v5, float v6, float v7, float v8) 
{
    data[0] = v1;
    data[1] = v2;
    data[2] = v3;
    data[3] = v4;
    data[4] = v5;
    data[5] = v6;
    data[6] = v7;
    data[7] = v8;
}

force_inline
simd<float, 256, nosimd_tag>::simd(const impl_type& v)
{
    data[0] = v[0];
    data[1] = v[1];
    data[2] = v[2];
    data[3] = v[3];
    data[4] = v[4];
    data[5] = v[5];
    data[6] = v[6];
    data[7] = v[7];
};

force_inline
simd<float, 256, nosimd_tag>::simd(const simd<float, 256, sse_tag>& s)
{
    s.store(data, std::true_type());
}

force_inline
simd<float, 256, nosimd_tag>::simd(const simd<float, 256, avx_tag>& s)
{
    s.store(data, std::true_type());
}

force_inline
float simd<float, 256, nosimd_tag>::get(int pos) const  
{ 
    return data[pos] ; 
};

force_inline
float simd<float, 256, nosimd_tag>::first() const  
{ 
    return data[0]; 
};

template<int Pos>
force_inline
float simd<float, 256, nosimd_tag>::get() const  
{ 
    return data[Pos]; 
};

force_inline
void simd<float, 256, nosimd_tag>::set(int pos, float val)
{ 
    data[pos] = val; 
};

template<int Pos>
force_inline
void simd<float, 256, nosimd_tag>::set(float val)
{ 
    data[Pos] = val; 
};

force_inline
const float* simd<float, 256, nosimd_tag>::get_raw_ptr() const
{ 
    return data; 
};

force_inline
float* simd<float, 256, nosimd_tag>::get_raw_ptr()
{ 
    return data; 
};

force_inline simd<float, 256, nosimd_tag>::simd_half
simd<float, 256, nosimd_tag>::extract_low() const
{
    return simd_half(data[0], data[1], data[2], data[3]);
}

force_inline simd<float, 256, nosimd_tag>::simd_half
simd<float, 256, nosimd_tag>::extract_high() const
{
    return simd_half(data[4], data[5], data[6], data[7]);
}

force_inline
simd<float, 256, nosimd_tag> simd<float, 256, nosimd_tag>::zero()
{
    simd<float, 256, nosimd_tag> ret;

    ret.data[0] = 0.0f;
    ret.data[1] = 0.0f;
    ret.data[2] = 0.0f;
    ret.data[3] = 0.0f;
    ret.data[4] = 0.0f;
    ret.data[5] = 0.0f;
    ret.data[6] = 0.0f;
    ret.data[7] = 0.0f;

    return ret;
}

force_inline
simd<float, 256, nosimd_tag> simd<float, 256, nosimd_tag>::minus_zero()
{
    simd<float, 256, nosimd_tag> ret;

    ret.data[0] = -0.0f;
    ret.data[1] = -0.0f;
    ret.data[2] = -0.0f;
    ret.data[3] = -0.0f;
    ret.data[4] = -0.0f;
    ret.data[5] = -0.0f;
    ret.data[6] = -0.0f;
    ret.data[7] = -0.0f;

    return ret;
}

force_inline
simd<float, 256, nosimd_tag> simd<float, 256, nosimd_tag>::one()
{
    simd<float, 256, nosimd_tag> ret;

    ret.data[0] = 1.0f;
    ret.data[1] = 1.0f;
    ret.data[2] = 1.0f;
    ret.data[3] = 1.0f;
    ret.data[4] = 1.0f;
    ret.data[5] = 1.0f;
    ret.data[6] = 1.0f;
    ret.data[7] = 1.0f;

    return ret;
}

force_inline
simd<float, 256, nosimd_tag> simd<float, 256, nosimd_tag>::minus_one()
{
    simd<float, 256, nosimd_tag> ret;

    ret.data[0] = -1.0f;
    ret.data[1] = -1.0f;
    ret.data[2] = -1.0f;
    ret.data[3] = -1.0f;
    ret.data[4] = -1.0f;
    ret.data[5] = -1.0f;
    ret.data[6] = -1.0f;
    ret.data[7] = -1.0f;

    return ret;
}

force_inline simd<float, 256, nosimd_tag> 
simd<float, 256, nosimd_tag>::load(const float* arr, std::true_type aligned)
{
    (void)aligned;
    simd<float, 256, nosimd_tag> ret;

    ret.data[0] = arr[0];
    ret.data[1] = arr[1];
    ret.data[2] = arr[2];
    ret.data[3] = arr[3];
    ret.data[4] = arr[4];
    ret.data[5] = arr[5];
    ret.data[6] = arr[6];
    ret.data[7] = arr[7];

    return ret;
};

force_inline simd<float, 256, nosimd_tag> 
simd<float, 256, nosimd_tag>::load(const float* arr, std::false_type not_aligned)
{
    (void)not_aligned;
    simd<float, 256, nosimd_tag> ret;

    ret.data[0] = arr[0];
    ret.data[1] = arr[1];
    ret.data[2] = arr[2];
    ret.data[3] = arr[3];
    ret.data[4] = arr[4];
    ret.data[5] = arr[5];
    ret.data[6] = arr[6];
    ret.data[7] = arr[7];

    return ret;
};

force_inline simd<float, 256, nosimd_tag> 
simd<float, 256, nosimd_tag>::gather(const float* arr, const simd_256_int32& ind)
{
    simd<float, 256, nosimd_tag> ret;

    ret.data[0] = arr[ind.data[0]];
    ret.data[1] = arr[ind.data[1]];
    ret.data[2] = arr[ind.data[2]];
    ret.data[3] = arr[ind.data[3]];
    ret.data[4] = arr[ind.data[4]];
    ret.data[5] = arr[ind.data[5]];
    ret.data[6] = arr[ind.data[6]];
    ret.data[7] = arr[ind.data[7]];

    return ret;
}

force_inline simd<float, 256, nosimd_tag> 
simd<float, 256, nosimd_tag>::gather(const float* arr, const simd_256_int64& ind)
{
    simd<float, 256, nosimd_tag> ret;

    ret.data[0] = arr[ind.data[0]];
    ret.data[1] = arr[ind.data[1]];
    ret.data[2] = arr[ind.data[2]];
    ret.data[3] = arr[ind.data[3]];

    ret.data[4] = ret.data[0];
    ret.data[5] = ret.data[1];
    ret.data[6] = ret.data[2];
    ret.data[7] = ret.data[3];

    return ret;
}

force_inline simd<float, 256, nosimd_tag> 
simd<float, 256, nosimd_tag>::broadcast(const float* arr)
{
    simd<float, 256, nosimd_tag> ret;

    ret.data[0] = arr[0];
    ret.data[1] = arr[0];
    ret.data[2] = arr[0];
    ret.data[3] = arr[0];
    ret.data[4] = arr[0];
    ret.data[5] = arr[0];
    ret.data[6] = arr[0];
    ret.data[7] = arr[0];

    return ret;
};

force_inline simd<float, 256, nosimd_tag> 
simd<float, 256, nosimd_tag>::broadcast(const float& arr)
{
    simd<float, 256, nosimd_tag> ret;

    ret.data[0] = arr;
    ret.data[1] = arr;
    ret.data[2] = arr;
    ret.data[3] = arr;
    ret.data[4] = arr;
    ret.data[5] = arr;
    ret.data[6] = arr;
    ret.data[7] = arr;

    return ret;
};

force_inline void 
simd<float, 256, nosimd_tag>::store(float* arr, std::true_type aligned) const
{
    (void)aligned;
    arr[0] = data[0];
    arr[1] = data[1];
    arr[2] = data[2];
    arr[3] = data[3];
    arr[4] = data[4];
    arr[5] = data[5];
    arr[6] = data[6];
    arr[7] = data[7];
};

force_inline void 
simd<float, 256, nosimd_tag>::store(float* arr, std::false_type not_aligned) const
{
    (void)not_aligned;
    arr[0] = data[0];
    arr[1] = data[1];
    arr[2] = data[2];
    arr[3] = data[3];
    arr[4] = data[4];
    arr[5] = data[5];
    arr[6] = data[6];
    arr[7] = data[7];
};

force_inline simd<double, 256, nosimd_tag>
simd<float, 256, nosimd_tag>::convert_low_to_double() const
{
    simd<double, 256, nosimd_tag> res;
    res.data[0] = data[0];
    res.data[1] = data[1];
    res.data[2] = data[2];
    res.data[3] = data[3];

    return res;
};

force_inline simd<double, 256, nosimd_tag>
simd<float, 256, nosimd_tag>::convert_high_to_double() const
{
    simd<double, 256, nosimd_tag> res;
    res.data[0] = data[4];
    res.data[1] = data[5];
    res.data[2] = data[6];
    res.data[3] = data[7];

    return res;
}

force_inline simd<int32_t, 256, nosimd_tag>
simd<float, 256, nosimd_tag>::convert_to_int32() const
{
    simd<int32_t, 256, nosimd_tag> res;
    res.data[0] = int32_t(data[0]);
    res.data[1] = int32_t(data[1]);
    res.data[2] = int32_t(data[2]);
    res.data[3] = int32_t(data[3]);
    res.data[4] = int32_t(data[4]);
    res.data[5] = int32_t(data[5]);
    res.data[6] = int32_t(data[6]);
    res.data[7] = int32_t(data[7]);

    return res;
}

force_inline simd<int64_t, 256, nosimd_tag>
simd<float, 256, nosimd_tag>::reinterpret_as_int64() const
{
    return *reinterpret_cast<const simd_256_int64*>(this);
}

force_inline simd<double, 256, nosimd_tag>
simd<float, 256, nosimd_tag>::reinterpret_as_double() const
{
    return *reinterpret_cast<const simd_256_double*>(this);
}

force_inline simd<int32_t, 256, nosimd_tag>
simd<float, 256, nosimd_tag>::reinterpret_as_int32() const
{
    return *reinterpret_cast<const simd_256_int32*>(this);
}

template<int Step>
force_inline
void simd<float, 256, nosimd_tag>::scatter(float* arr) const
{
    arr[0*Step] = data[0];
    arr[1*Step] = data[1];
    arr[2*Step] = data[2];
    arr[3*Step] = data[3];
    arr[4*Step] = data[4];
    arr[5*Step] = data[5];
    arr[6*Step] = data[6];
    arr[7*Step] = data[7];
};

force_inline simd<float, 256, nosimd_tag>& 
simd<float, 256, nosimd_tag>::operator+=(const simd& x)
{
    *this = *this + x;
    return *this;
}

force_inline simd<float, 256, nosimd_tag>& 
simd<float, 256, nosimd_tag>::operator-=(const simd& x)
{
    *this = *this - x;
    return *this;
}

force_inline simd<float, 256, nosimd_tag>& 
simd<float, 256, nosimd_tag>::operator*=(const simd& x)
{
    *this = *this * x;
    return *this;
}

force_inline simd<float, 256, nosimd_tag>& 
simd<float, 256, nosimd_tag>::operator/=(const simd& x)
{
    *this = *this / x;
    return *this;
}

}}
