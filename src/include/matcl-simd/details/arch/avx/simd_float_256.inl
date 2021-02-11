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

#include "matcl-simd/arch/avx/simd_float_256.h"

#pragma warning(push)
#pragma warning(disable : 4127) // conditional expression is constant

namespace matcl { namespace simd { namespace details
{

template<int I1, int I2>
struct is_continuous
{
    static const bool value = (I1 % 2 == 0) && I2 == I1 + 1;
};

}}}

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
simd<float, 256, avx_tag>::simd(const simd<float, 128, scalar_sse_tag>& s)
{
    #if MATCL_ARCHITECTURE_HAS_AVX2
        data = _mm256_broadcastss_ps(s.data);
    #else
        data = _mm256_set1_ps(s.first());
    #endif
}

force_inline
simd<float, 256, avx_tag>::simd(const simd<float, 128, scalar_nosimd_tag>& s)
    :simd(s.first())
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

template<int I1, int I2, int I3, int I4, int I5, int I6, int I7, int I8>
force_inline simd<float, 256, avx_tag>
simd<float, 256, avx_tag>::select() const
{
    static_assert(I1 >= 0 && I1 <= 7, "invalid index in select function");
    static_assert(I2 >= 0 && I2 <= 7, "invalid index in select function");
    static_assert(I3 >= 0 && I3 <= 7, "invalid index in select function");
    static_assert(I4 >= 0 && I4 <= 7, "invalid index in select function");
    static_assert(I5 >= 0 && I5 <= 7, "invalid index in select function");
    static_assert(I6 >= 0 && I6 <= 7, "invalid index in select function");
    static_assert(I7 >= 0 && I7 <= 7, "invalid index in select function");
    static_assert(I8 >= 0 && I8 <= 7, "invalid index in select function");

    if (I1 == 0 && I2 == 1 && I3 == 2 && I4 == 3 && I5 == 4 && I6 == 5 && I7 == 6 && I8 == 7) 
        return data;

    if (details::is_continuous<I1, I2>::value && details::is_continuous<I3, I4>::value &&
        details::is_continuous<I5, I6>::value && details::is_continuous<I7, I8>::value)
    {
        //64-bit select can be used
        return this->reinterpret_as_double().select<I1/2, I3/2, I5/2, I7/2>().reinterpret_as_float();
    }

    if (I1 <= 3 && I2 <= 3 && I3 <= 3 && I4 <= 3 &&
        I5 > 3  && I6 > 3  && I7 > 3  && I8 > 3) 
    {
        // no exchange of data between low and high half

        if (I5 == I1 + 4 && I6 == I2 + 4 && I7 == I3 + 4 && I8 == I4 + 4)
        {
            // the same pattern in low and high half

            static const int ind0 = I1 + (I2 << 2) + (I3 << 4) + (I4 << 6);
            static const int ind  = std::min(std::max(ind0, 0), 255);
            return _mm256_shuffle_ps(data, data, ind);
        }

        __m256i ind = details::vector_8_int<I1, I2, I3, I4, I5, I6, I7, I8>();
        return _mm256_permutevar_ps(data, ind);
    }

    // general case

    #if MATCL_ARCHITECTURE_HAS_AVX2
        static const __m256i ind = details::vector_8_int<I1, I2, I3, I4, I5, I6, I7, I8>();
        return _mm256_permutevar8x32_ps(data, ind);
    #else
        simd_half lo        = this->extract_low();
        simd_half hi        = this->extract_high();

        simd_half res_lo    = simd_half::combine<I1, I2, I3, I4>(lo, hi);
        simd_half res_hi    = simd_half::combine<I5, I6, I7, I8>(lo, hi);

        return simd(res_lo, res_hi);
    #endif
};

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
simd<float, 256, avx_tag>::gather(const float* arr, const simd_int32& ind)
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
simd<float, 256, avx_tag>::gather(const float* arr, const simd_int64& ind)
{
    #if MATCL_ARCHITECTURE_HAS_AVX2
        return simd(simd_half(_mm256_i64gather_ps(arr, ind.data, 4)), simd_half::zero());
    #else
        simd_half res;
        float* res_ptr          = res.get_raw_ptr();
        const int64_t* ind_ptr  = ind.get_raw_ptr();

        for (int i = 0; i < vector_size/2; ++i)
            res_ptr[i]  = arr[ind_ptr[i]];

        return simd(res, simd_half::zero());
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

force_inline simd<int64_t, 256, avx_tag>
simd<float, 256, avx_tag>::convert_low_to_int64() const
{
    // intrinsic is not available
    
    simd_int64 ret;

    const float* ptr    = this->get_raw_ptr();
    int64_t* ptr_ret    = ret.get_raw_ptr();

    ptr_ret[0]          = scalar_func::convert_float_int64(ptr[0]);
    ptr_ret[1]          = scalar_func::convert_float_int64(ptr[1]);
    ptr_ret[2]          = scalar_func::convert_float_int64(ptr[2]);
    ptr_ret[3]          = scalar_func::convert_float_int64(ptr[3]);

    return ret;
};

force_inline
simd<int64_t, 256, avx_tag>
simd<float, 256, avx_tag>::convert_high_to_int64() const
{
    // intrinsic is not available
    
    simd_int64 ret;

    const float* ptr    = this->get_raw_ptr();
    int64_t* ptr_ret    = ret.get_raw_ptr();

    ptr_ret[0]          = scalar_func::convert_float_int64(ptr[4]);
    ptr_ret[1]          = scalar_func::convert_float_int64(ptr[5]);
    ptr_ret[2]          = scalar_func::convert_float_int64(ptr[6]);
    ptr_ret[3]          = scalar_func::convert_float_int64(ptr[7]);

    return ret;
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

#pragma warning(pop)