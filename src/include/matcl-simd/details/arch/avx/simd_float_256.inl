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

#include "matcl-simd/arch/avx/simd_float_256.h"

namespace matcl { namespace simd
{

//-------------------------------------------------------------------
//                          AVX SINGLE
//-------------------------------------------------------------------

force_inline
simd<float, 256, avx_tag>::simd(float val) 
    : data(_mm256_set1_ps(val)) 
{}

force_inline
simd<float, 256, avx_tag>::simd(const simd_half& lo_hi)
    : data(_mm256_setr_m128(lo_hi.data, lo_hi.data))
{}

force_inline
simd<float, 256, avx_tag>::simd(float v1, float v2, float v3, float v4, float v5, float v6, float v7, float v8) 
    : data(_mm256_setr_ps(v1, v2, v3, v4, v5, v6, v7, v8)) 
{}

force_inline
simd<float, 256, avx_tag>::simd(const impl_type& v)
    : data(v) 
{};

force_inline
simd<float, 256, avx_tag>::simd(const simd_half& lo, const simd_half& hi)
    : data(_mm256_setr_m128(lo.data, hi.data))
{}

force_inline
simd<float, 256, avx_tag>::simd(const simd<float, 256, nosimd_tag>& s)
    : data(_mm256_load_ps(s.data))
{}

force_inline
simd<float, 256, avx_tag>::simd(const simd<float, 256, sse_tag>& s)
    : simd(s.data[0], s.data[1])
{}

force_inline
float simd<float, 256, avx_tag>::get(int pos) const  
{ 
    return get_raw_ptr()[pos] ; 
};

force_inline
float simd<float, 256, avx_tag>::first() const
{ 
    __m128 ds  = _mm256_castps256_ps128(data);
    return _mm_cvtss_f32(ds);
};

force_inline
void simd<float, 256, avx_tag>::set(int pos, float val)
{ 
    get_raw_ptr()[pos] = val; 
};

force_inline
const float* simd<float, 256, avx_tag>::get_raw_ptr() const
{ 
    return reinterpret_cast<const float*>(&data); 
};

force_inline
float* simd<float, 256, avx_tag>::get_raw_ptr()
{ 
    return reinterpret_cast<float*>(&data); 
};

force_inline simd<float, 256, avx_tag>::simd_half
simd<float, 256, avx_tag>::extract_low() const
{
    return _mm256_castps256_ps128(data);
}

force_inline simd<float, 256, avx_tag>::simd_half
simd<float, 256, avx_tag>::extract_high() const
{
    return _mm256_extractf128_ps(data, 1);
}

force_inline
simd<float, 256, avx_tag> simd<float, 256, avx_tag>::zero()
{
    impl_type data  = _mm256_setzero_ps();
    return data;
}

force_inline
simd<float, 256, avx_tag> simd<float, 256, avx_tag>::minus_zero()
{
    return simd(-0.0f);
}

force_inline
simd<float, 256, avx_tag> simd<float, 256, avx_tag>::one()
{
    return simd(1.0f);
}

force_inline
simd<float, 256, avx_tag> simd<float, 256, avx_tag>::minus_one()
{
    return simd(-1.0f);
}

force_inline simd<float, 256, avx_tag> 
simd<float, 256, avx_tag>::load(const float* arr, std::true_type aligned)
{
    (void)aligned;
    return _mm256_load_ps(arr);
};

force_inline simd<float, 256, avx_tag> 
simd<float, 256, avx_tag>::load(const float* arr, std::false_type not_aligned)
{
    (void)not_aligned;
    return _mm256_loadu_ps(arr);
};

force_inline simd<float, 256, avx_tag> 
simd<float, 256, avx_tag>::gather(const float* arr, const simd_256_int32& ind)
{
    #if MATCL_ARCHITECTURE_HAS_AVX2
        return _mm256_i32gather_ps(arr, ind.data, 4);
    #else
        simd res;
        float* res_ptr          = res.get_raw_ptr();
        const int32_t* ind_ptr  = ind.get_raw_ptr();

        for (int i = 0; i < vector_size; ++i)
            res_ptr[i]  = arr[ind_ptr[i]];

        return res;
    #endif    
}

force_inline simd<float, 256, avx_tag> 
simd<float, 256, avx_tag>::gather(const float* arr, const simd_256_int64& ind)
{
    #if MATCL_ARCHITECTURE_HAS_AVX2
        return simd(_mm256_i64gather_ps(arr, ind.data, 4));
    #else
        simd_half res;
        float* res_ptr          = res.get_raw_ptr();
        const int64_t* ind_ptr  = ind.get_raw_ptr();

        for (int i = 0; i < vector_size/2; ++i)
            res_ptr[i]  = arr[ind_ptr[i]];

        return simd(res, res);
    #endif    
}

force_inline simd<float, 256, avx_tag> 
simd<float, 256, avx_tag>::broadcast(const float* arr)
{
    return _mm256_broadcast_ss(arr);
};

force_inline simd<float, 256, avx_tag> 
simd<float, 256, avx_tag>::broadcast(const float& arr)
{
    return _mm256_broadcast_ss(&arr);
};

force_inline void simd<float, 256, avx_tag>::store(float* arr, std::true_type aligned) const
{
    (void)aligned;
    _mm256_store_ps(arr, data);
};

force_inline void simd<float, 256, avx_tag>::store(float* arr, std::false_type not_aligned) const
{
    (void)not_aligned;
    _mm256_storeu_ps(arr, data);
};

force_inline simd<double, 256, avx_tag>
simd<float, 256, avx_tag>::convert_low_to_double() const
{
    simd_half lo    = this->extract_low();
    return _mm256_cvtps_pd(lo.data);
};

force_inline
simd<double, 256, avx_tag>
simd<float, 256, avx_tag>::convert_high_to_double() const
{
    simd_half hi    = this->extract_high();
    return _mm256_cvtps_pd(hi.data);
};

force_inline simd<int32_t, 256, avx_tag>
simd<float, 256, avx_tag>::convert_to_int32() const
{
    return _mm256_cvtps_epi32(data);
};

force_inline simd<double, 256, avx_tag>
simd<float, 256, avx_tag>::reinterpret_as_double() const
{
    return _mm256_castps_pd(data);
};

force_inline simd<int32_t, 256, avx_tag>
simd<float, 256, avx_tag>::reinterpret_as_int32() const
{
    return _mm256_castps_si256(data);
};

force_inline simd<int64_t, 256, avx_tag>
simd<float, 256, avx_tag>::reinterpret_as_int64() const
{
    return _mm256_castps_si256(data);
};

template<int Step>
force_inline
void simd<float, 256, avx_tag>::scatter(float* arr) const
{
    const float* ptr = get_raw_ptr();

    //no scatter intrinsic
    arr[0*Step] = ptr[0];
    arr[1*Step] = ptr[1];
    arr[2*Step] = ptr[2];
    arr[3*Step] = ptr[3];
    arr[4*Step] = ptr[4];
    arr[5*Step] = ptr[5];
    arr[6*Step] = ptr[6];
    arr[7*Step] = ptr[7];
};

force_inline simd<float, 256, avx_tag>& 
simd<float, 256, avx_tag>::operator+=(const simd& x)
{
    *this = *this + x;
    return *this;
}

force_inline simd<float, 256, avx_tag>& 
simd<float, 256, avx_tag>::operator-=(const simd& x)
{
    *this = *this - x;
    return *this;
}

force_inline simd<float, 256, avx_tag>& 
simd<float, 256, avx_tag>::operator*=(const simd& x)
{
    *this = *this * x;
    return *this;
}

force_inline simd<float, 256, avx_tag>& 
simd<float, 256, avx_tag>::operator/=(const simd& x)
{
    *this = *this / x;
    return *this;
}

}}
