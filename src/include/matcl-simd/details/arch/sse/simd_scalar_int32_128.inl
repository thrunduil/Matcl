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

#include "matcl-simd/arch/sse/simd_scalar_int32_128.h"
#include "matcl-simd/details/arch/sse/func/missing_intrinsics.h"

namespace matcl { namespace simd
{

//-------------------------------------------------------------------
//                          SSE SCALAR int32_t
//-------------------------------------------------------------------

force_inline
simd<int32_t, 128, scalar_sse_tag>::simd(int32_t val)
    : data(_mm_set1_epi32(val)) 
{}

force_inline
simd<int32_t, 128, scalar_sse_tag>::simd(const impl_type& v)
    : data(v) 
{};

force_inline
simd<int32_t, 128, scalar_sse_tag>::simd(const simd<int32_t, 128, nosimd_tag>& s)
    : data(_mm_load_si128((const __m128i*)s.data))
{};

force_inline
simd<int32_t, 128, scalar_sse_tag>::simd(const simd<int32_t, 128, sse_tag>& s)
    : data(s.data)
{};

force_inline
simd<int32_t, 128, scalar_sse_tag>::simd(const simd<int32_t, 128, scalar_nosimd_tag>& s)
    : simd(broadcast(s.data))
{};

force_inline simd<int32_t, 128, scalar_sse_tag> 
simd<int32_t, 128, scalar_sse_tag>::broadcast(const int32_t* arr)
{ 
    return _mm_set1_epi32(arr[0]);
};

force_inline simd<int32_t, 128, scalar_sse_tag> 
simd<int32_t, 128, scalar_sse_tag>::broadcast(const int32_t& arr)
{ 
    return _mm_set1_epi32(arr);
};

force_inline
int32_t simd<int32_t, 128, scalar_sse_tag>::get(int pos) const
{ 
    (void)pos;
    return first();
};

force_inline
int32_t simd<int32_t, 128, scalar_sse_tag>::first() const
{ 
    return _mm_cvtsi128_si32(data); 
};

force_inline
void simd<int32_t, 128, scalar_sse_tag>::set(int pos, int32_t val)
{ 
    (void)pos;
    data = _mm_set1_epi32(val); 
};

force_inline
const int32_t* simd<int32_t, 128, scalar_sse_tag>::get_raw_ptr() const
{ 
    return reinterpret_cast<const int32_t*>(&data); 
};

force_inline
int32_t* simd<int32_t, 128, scalar_sse_tag>::get_raw_ptr()
{ 
    return reinterpret_cast<int32_t*>(&data); 
};

force_inline
simd<int32_t, 128, scalar_sse_tag> simd<int32_t, 128, scalar_sse_tag>::zero()
{
    impl_type data   = _mm_setzero_si128();
    return data;
}

force_inline
simd<int32_t, 128, scalar_sse_tag> simd<int32_t, 128, scalar_sse_tag>::one()
{
    return simd(1);
}

force_inline
simd<int32_t, 128, scalar_sse_tag> simd<int32_t, 128, scalar_sse_tag>::minus_one()
{
    return simd(-1);
}

force_inline simd<int32_t, 128, scalar_sse_tag> 
simd<int32_t, 128, scalar_sse_tag>::load(const int32_t* arr, std::true_type aligned)
{
    (void)aligned;
    return _mm_set1_epi32(arr[0]);
};

force_inline simd<int32_t, 128, scalar_sse_tag> 
simd<int32_t, 128, scalar_sse_tag>::load(const int32_t* arr, std::false_type not_aligned)
{
    (void)not_aligned;
    return _mm_set1_epi32(arr[0]);
};

force_inline simd<int32_t, 128, scalar_sse_tag>
simd<int32_t, 128, scalar_sse_tag>::gather(const int32_t* arr, const simd_int32& ind)
{
    return simd(arr[ind.first()]);
}

force_inline simd<int32_t, 128, scalar_sse_tag>
simd<int32_t, 128, scalar_sse_tag>::gather(const int32_t* arr, const simd_int64& ind)
{
    return simd(arr[ind.first()]);
}

force_inline simd<int32_t, 128, scalar_sse_tag> 
simd<int32_t, 128, scalar_sse_tag>::set_lower(int32_t v)
{
    return _mm_cvtsi32_si128(v);
};

force_inline void 
simd<int32_t, 128, scalar_sse_tag>::store(int32_t* arr, std::true_type aligned) const
{
    (void)aligned;
    arr[0] = first();
};


force_inline void 
simd<int32_t, 128, scalar_sse_tag>::store(int32_t* arr, std::false_type not_aligned) const
{
    (void)not_aligned;
    arr[0] = first();
};

force_inline simd<int64_t, 128, scalar_sse_tag>
simd<int32_t, 128, scalar_sse_tag>::convert_to_int64() const
{
    return simd<int64_t, 128, scalar_sse_tag>(int64_t(first()));
};

force_inline simd<float, 128, scalar_sse_tag>
simd<int32_t, 128, scalar_sse_tag>::convert_to_float() const
{
    return _mm_cvtepi32_ps(data);
};

force_inline simd<double, 128, scalar_sse_tag>
simd<int32_t, 128, scalar_sse_tag>::convert_to_double() const
{
    return _mm_cvtepi32_pd(data);
};

force_inline simd<double, 128, scalar_sse_tag>
simd<int32_t, 128, scalar_sse_tag>::reinterpret_as_double() const
{
    return _mm_castsi128_pd(data);
}

force_inline simd<float, 128, scalar_sse_tag>
simd<int32_t, 128, scalar_sse_tag>::reinterpret_as_float() const
{
    return _mm_castsi128_ps(data);
}

force_inline simd<int64_t, 128, scalar_sse_tag>
simd<int32_t, 128, scalar_sse_tag>::reinterpret_as_int64() const
{
    return data;
}

force_inline simd<int64_t, 128, sse_tag>
simd<int64_t, 128, scalar_sse_tag>::as_vector() const
{
    return data;
}

template<int Step>
force_inline
void simd<int32_t, 128, scalar_sse_tag>::scatter(int32_t* arr) const
{
    arr[0] = first();
};

force_inline simd<int32_t, 128, scalar_sse_tag>& 
simd<int32_t, 128, scalar_sse_tag>::operator+=(const simd& x)
{
    *this = *this + x;
    return *this;
}

force_inline simd<int32_t, 128, scalar_sse_tag>& 
simd<int32_t, 128, scalar_sse_tag>::operator-=(const simd& x)
{
    *this = *this - x;
    return *this;
}

force_inline simd<int32_t, 128, scalar_sse_tag>& 
simd<int32_t, 128, scalar_sse_tag>::operator*=(const simd& x)
{
    *this = *this * x;
    return *this;
}

}}
