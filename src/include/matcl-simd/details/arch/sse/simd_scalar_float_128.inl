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

#include "matcl-simd/arch/sse/simd_scalar_float_128.h"

namespace matcl { namespace simd
{

//-------------------------------------------------------------------
//                          SSE SCALAR SINGLE
//-------------------------------------------------------------------

force_inline
simd<float, 128, scalar_sse_tag>::simd(float val)
    : data(_mm_set1_ps(val)) 
{}

force_inline
simd<float, 128, scalar_sse_tag>::simd(const impl_type& v)
    : data(v) 
{};

force_inline
simd<float, 128, scalar_sse_tag>::simd(const simd<float, 128, nosimd_tag>& s)
    : data(_mm_load_ps(s.data))
{};

force_inline
simd<float, 128, scalar_sse_tag>::simd(const simd<float, 128, sse_tag>& s)
    : data(s.data)
{};

force_inline
simd<float, 128, scalar_sse_tag>::simd(const simd<float, 128, scalar_nosimd_tag>& s)
    : simd(broadcast(s.data))
{};

force_inline simd<float, 128, scalar_sse_tag> 
simd<float, 128, scalar_sse_tag>::broadcast(const float* arr)
{ 
    return _mm_load1_ps(arr);
};

force_inline simd<float, 128, scalar_sse_tag> 
simd<float, 128, scalar_sse_tag>::broadcast(const float& arr)
{ 
    return _mm_load1_ps(&arr);
};

force_inline
float simd<float, 128, scalar_sse_tag>::get(int pos) const
{ 
    (void)pos;
    return first();
};

force_inline
float simd<float, 128, scalar_sse_tag>::first() const
{ 
    return _mm_cvtss_f32(data); 
};

force_inline
void simd<float, 128, scalar_sse_tag>::set(int pos, float val)
{ 
    (void)pos;
    data = _mm_set1_ps(val);
};

force_inline
const float* simd<float, 128, scalar_sse_tag>::get_raw_ptr() const
{ 
    return reinterpret_cast<const float*>(&data); 
};

force_inline
float* simd<float, 128, scalar_sse_tag>::get_raw_ptr()
{ 
    return reinterpret_cast<float*>(&data); 
};

force_inline
simd<float, 128, scalar_sse_tag> simd<float, 128, scalar_sse_tag>::zero()
{
    impl_type data   = _mm_setzero_ps();
    return data;
}

force_inline
simd<float, 128, scalar_sse_tag> simd<float, 128, scalar_sse_tag>::minus_zero()
{
    return simd(-0.0f);
}

force_inline
simd<float, 128, scalar_sse_tag> simd<float, 128, scalar_sse_tag>::one()
{
    return simd(1.0f);
}

force_inline
simd<float, 128, scalar_sse_tag> simd<float, 128, scalar_sse_tag>::minus_one()
{
    return simd(-1.0f);
}

force_inline simd<float, 128, scalar_sse_tag> 
simd<float, 128, scalar_sse_tag>::load(const float* arr, std::true_type aligned)
{
    (void)aligned;
    return _mm_load1_ps(arr);
};

force_inline simd<float, 128, scalar_sse_tag> 
simd<float, 128, scalar_sse_tag>::load(const float* arr, std::false_type not_aligned)
{
    (void)not_aligned;
    return _mm_load1_ps(arr);
};

force_inline simd<float, 128, scalar_sse_tag>
simd<float, 128, scalar_sse_tag>::gather(const float* arr, const simd_int32& ind)
{
    return simd(arr[ind.first()]);
}

force_inline simd<float, 128, scalar_sse_tag> 
simd<float, 128, scalar_sse_tag>::gather(const float* arr, const simd_int64& ind)
{
    return simd(arr[ind.first()]);
}

force_inline simd<float, 128, scalar_sse_tag> 
simd<float, 128, scalar_sse_tag>::set_lower(float v)
{
    return _mm_set_ss(v);
};

force_inline void simd<float, 128, scalar_sse_tag>::store(float* arr, std::true_type aligned) const
{
    (void)aligned;
    arr[0] = first();
};


force_inline void 
simd<float, 128, scalar_sse_tag>::store(float* arr, std::false_type not_aligned) const
{
    (void)not_aligned;
    arr[0] = first();
};

force_inline
typename simd<double, 128, scalar_sse_tag>
simd<float, 128, scalar_sse_tag>::convert_to_double() const
{
    return _mm_cvtps_pd(data);  
};

force_inline simd<double, 128, scalar_sse_tag>
simd<float, 128, scalar_sse_tag>::reinterpret_as_double() const
{
    return _mm_castps_pd(data);
}

force_inline simd<int32_t, 128, scalar_sse_tag>
simd<float, 128, scalar_sse_tag>::reinterpret_as_int32() const
{
    return _mm_castps_si128(data);
}

force_inline simd<int64_t, 128, scalar_sse_tag>
simd<float, 128, scalar_sse_tag>::reinterpret_as_int64() const
{
    return _mm_castps_si128(data);
}

force_inline simd<int32_t, 128, scalar_sse_tag>
simd<float, 128, scalar_sse_tag>::convert_to_int32() const
{
    return _mm_cvtps_epi32(data);
};

force_inline
simd<int64_t, 128, scalar_sse_tag>
simd<float, 128, scalar_sse_tag>::convert_to_int64() const
{
    return simd<int64_t, 128, scalar_sse_tag>(scalar_func::convert_float_int64(first()));
};

force_inline simd<float, 128, sse_tag>
simd<float, 128, scalar_sse_tag>::as_vector() const
{
    return data;
}

template<int Step>
force_inline
void simd<float, 128, scalar_sse_tag>::scatter(float* arr) const
{
    arr[0] = first();
};

force_inline simd<float, 128, scalar_sse_tag>& 
simd<float, 128, scalar_sse_tag>::operator+=(const simd& x)
{
    *this = *this + x;
    return *this;
}

force_inline simd<float, 128, scalar_sse_tag>& 
simd<float, 128, scalar_sse_tag>::operator-=(const simd& x)
{
    *this = *this - x;
    return *this;
}

force_inline simd<float, 128, scalar_sse_tag>& 
simd<float, 128, scalar_sse_tag>::operator*=(const simd& x)
{
    *this = *this * x;
    return *this;
}

force_inline simd<float, 128, scalar_sse_tag>& 
simd<float, 128, scalar_sse_tag>::operator/=(const simd& x)
{
    *this = *this / x;
    return *this;
}

}}
