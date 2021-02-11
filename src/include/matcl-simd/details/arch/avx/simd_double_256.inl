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

#include "matcl-simd/arch/avx/simd_double_256.h"

#pragma warning(push)
#pragma warning(disable : 4127) // conditional expression is constant

namespace matcl { namespace simd
{

//-------------------------------------------------------------------
//                          AVX DOUBLE
//-------------------------------------------------------------------

force_inline
simd<double, 256, avx_tag>::simd(int32_t val)
    : data(_mm256_set1_pd(val)) 
{}

force_inline
simd<double, 256, avx_tag>::simd(int64_t val)
    : data(_mm256_set1_pd(double(val))) 
{}

force_inline
simd<double, 256, avx_tag>::simd(float val)
    : data(_mm256_set1_pd(val)) 
{}

force_inline
simd<double, 256, avx_tag>::simd(double val)
    : data(_mm256_set1_pd(val)) 
{}

force_inline
simd<double, 256, avx_tag>::simd(double v1, double v2, double v3, double v4) 
    : data(_mm256_setr_pd(v1, v2, v3, v4)) 
{}

force_inline
simd<double, 256, avx_tag>::simd(const simd_half& lo, const simd_half& hi)
    : data(_mm256_setr_m128d(lo.data, hi.data))
{}

force_inline
simd<double, 256, avx_tag>::simd(const simd_half& lo_hi)
    : data(_mm256_setr_m128d(lo_hi.data, lo_hi.data))
{}

force_inline
simd<double, 256, avx_tag>::simd(const impl_type& v)
    : data(v)
{};

force_inline
simd<double, 256, avx_tag>::simd(const simd<double, 256, nosimd_tag>& s)
    : data(_mm256_load_pd(s.data))
{}

force_inline
simd<double, 256, avx_tag>::simd(const simd<double, 256, sse_tag>& s)
    : simd(s.data[0], s.data[1])
{}

force_inline
simd<double, 256, avx_tag>::simd(const simd<double, 128, scalar_sse_tag>& s)
{
    #if MATCL_ARCHITECTURE_HAS_AVX2
        data = _mm256_broadcastsd_pd(s.data);
    #else
        data = _mm256_set1_pd(s.first());
    #endif
}

force_inline
simd<double, 256, avx_tag>::simd(const simd<double, 128, scalar_nosimd_tag>& s)
    :simd(s.first())
{}

force_inline
double simd<double, 256, avx_tag>::get(int pos) const
{ 
    return get_raw_ptr()[pos] ; 
};

force_inline
double simd<double, 256, avx_tag>::first() const
{ 
    __m128d ds  = _mm256_castpd256_pd128(data);
    return _mm_cvtsd_f64(ds);
};

force_inline
void simd<double, 256, avx_tag>::set(int pos, double val)
{ 
    get_raw_ptr()[pos] = val; 
};

force_inline
const double* simd<double, 256, avx_tag>::get_raw_ptr() const
{ 
    return reinterpret_cast<const double*>(&data); 
};

force_inline
double* simd<double, 256, avx_tag>::get_raw_ptr()
{ 
    return reinterpret_cast<double*>(&data); 
};

force_inline simd<double, 256, avx_tag>::simd_half
simd<double, 256, avx_tag>::extract_low() const
{
    return _mm256_castpd256_pd128(data);
}

force_inline simd<double, 256, avx_tag>::simd_half
simd<double, 256, avx_tag>::extract_high() const
{
    return _mm256_extractf128_pd(data, 1);
}

template<int I1, int I2, int I3, int I4>
force_inline simd<double, 256, avx_tag>
simd<double, 256, avx_tag>::select() const
{
    static_assert(I1 >= 0 && I1 <= 3, "invalid index in select function");
    static_assert(I2 >= 0 && I2 <= 3, "invalid index in select function");
    static_assert(I3 >= 0 && I3 <= 3, "invalid index in select function");
    static_assert(I4 >= 0 && I4 <= 3, "invalid index in select function");

    if (I1 == 0 && I2 == 1 && I3 == 2 && I4 == 3) 
        return data;

    if (I1 <= 1 && I2 <= 1 && I3 > 1 && I4 > 1) 
    {
        // no exchange of data between low and high half

        static const int ind0   = I1 + (I2 << 1) + ((I3 - 2) << 2) + ((I4 - 2) << 3);
        static const int ind    = std::min(std::max(ind0, 0), 15);
        return _mm256_shuffle_pd (data, data, ind);
    }

    // general case

    #if MATCL_ARCHITECTURE_HAS_AVX2
        static const int ind   = I1 + (I2 << 2) + (I3 << 4) + (I4 << 6);

        __m256i datai   = _mm256_castpd_si256(data);
        __m256i reti    = _mm256_permute4x64_epi64(datai, ind);
        __m256d ret     = _mm256_castsi256_pd(reti);
        return ret;
    #else
        simd_half lo        = this->extract_low();
        simd_half hi        = this->extract_high();

        simd_half res_lo    = simd_half::combine<I1, I2>(lo, hi);
        simd_half res_hi    = simd_half::combine<I3, I4>(lo, hi);

        return simd(res_lo, res_hi);
    #endif
};

force_inline
simd<double, 256, avx_tag> simd<double, 256, avx_tag>::zero()
{
    impl_type data  = _mm256_setzero_pd();
    return data;
}

force_inline
simd<double, 256, avx_tag> simd<double, 256, avx_tag>::minus_zero()
{
    return simd(-0.0);
}

force_inline
simd<double, 256, avx_tag> simd<double, 256, avx_tag>::one()
{
    return simd(1.0);
}

force_inline
simd<double, 256, avx_tag> simd<double, 256, avx_tag>::minus_one()
{
    return simd(-1.0);
}

force_inline simd<double, 256, avx_tag> 
simd<double, 256, avx_tag>::load(const double* arr, std::true_type aligned)
{
    (void)aligned;
    return _mm256_load_pd(arr);
};

force_inline simd<double, 256, avx_tag> 
simd<double, 256, avx_tag>::load(const double* arr, std::false_type not_aligned)
{
    (void)not_aligned;
    return _mm256_loadu_pd(arr);
};

force_inline simd<double, 256, avx_tag> 
simd<double, 256, avx_tag>::gather(const double* arr, const simd_int32_half& ind)
{
    #if MATCL_ARCHITECTURE_HAS_AVX2
        return _mm256_i32gather_pd(arr, ind.data, 8);
    #else
        simd<double, 256, avx_tag> res;
        double* res_ptr         = res.get_raw_ptr();
        const int32_t* ind_ptr  = ind.get_raw_ptr();

        for (int i = 0; i < vector_size; ++i)
            res_ptr[i]  = arr[ind_ptr[i]];

        return res;
    #endif
}

force_inline simd<double, 256, avx_tag> 
simd<double, 256, avx_tag>::gather(const double* arr, const simd_int32& ind)
{
    return gather(arr, ind.extract_low());
}

force_inline simd<double, 256, avx_tag> 
simd<double, 256, avx_tag>::gather(const double* arr, const simd_int64& ind)
{
    #if MATCL_ARCHITECTURE_HAS_AVX2
        return _mm256_i64gather_pd(arr, ind.data, 8);
    #else
        simd<double, 256, avx_tag> res;
        double* res_ptr         = res.get_raw_ptr();
        const int64_t* ind_ptr  = ind.get_raw_ptr();

        for (int i = 0; i < vector_size; ++i)
            res_ptr[i]  = arr[ind_ptr[i]];

        return res;
    #endif
}

force_inline simd<double, 256, avx_tag> 
simd<double, 256, avx_tag>::broadcast(const double* arr)
{
    return _mm256_broadcast_sd(arr);
};

force_inline simd<double, 256, avx_tag> 
simd<double, 256, avx_tag>::broadcast(const double& arr)
{
    return _mm256_broadcast_sd(&arr);
};

force_inline void 
simd<double, 256, avx_tag>::store(double* arr, std::true_type aligned) const
{
    (void)aligned;
    _mm256_store_pd(arr, data);
};

force_inline void 
simd<double, 256, avx_tag>::store(double* arr, std::false_type not_aligned) const
{
    (void)not_aligned;
    _mm256_storeu_pd(arr, data);
};

force_inline simd<float, 128, sse_tag> 
simd<double, 256, avx_tag>::convert_to_float() const
{
    return _mm256_cvtpd_ps(data);
};

force_inline simd<int32_t, 128, sse_tag> 
simd<double, 256, avx_tag>::convert_to_int32() const
{
    return _mm256_cvtpd_epi32(data);
};

force_inline simd<int64_t, 256, avx_tag> 
simd<double, 256, avx_tag>::convert_to_int64() const
{
    simd<int64_t, 256, avx_tag> ret;

    int64_t* ret_ptr    = ret.get_raw_ptr();
    const double* ptr   = this->get_raw_ptr();

    for (int i = 0; i < vector_size; ++i)
        ret_ptr[i]      = scalar_func::convert_double_int64(ptr[i]);
    
    return ret;
};

force_inline simd<float, 256, avx_tag> 
simd<double, 256, avx_tag>::reinterpret_as_float() const
{
    return _mm256_castpd_ps(data);
}

force_inline simd<int32_t, 256, avx_tag> 
simd<double, 256, avx_tag>::reinterpret_as_int32() const
{
    return _mm256_castpd_si256(data);
};

force_inline simd<int64_t, 256, avx_tag> 
simd<double, 256, avx_tag>::reinterpret_as_int64() const
{
    return _mm256_castpd_si256(data);
}

template<int Step>
force_inline
void simd<double, 256, avx_tag>::scatter(double* arr) const
{
    const double* ptr = get_raw_ptr();

    //no scatter intrinsic
    arr[0*Step] = ptr[0];
    arr[1*Step] = ptr[1];
    arr[2*Step] = ptr[2];
    arr[3*Step] = ptr[3];
};

force_inline simd<double, 256, avx_tag>& 
simd<double, 256, avx_tag>::operator+=(const simd& x)
{
    *this = *this + x;
    return *this;
}

force_inline simd<double, 256, avx_tag>& 
simd<double, 256, avx_tag>::operator-=(const simd& x)
{
    *this = *this - x;
    return *this;
}

force_inline simd<double, 256, avx_tag>& 
simd<double, 256, avx_tag>::operator*=(const simd& x)
{
    *this = *this * x;
    return *this;
}

force_inline simd<double, 256, avx_tag>& 
simd<double, 256, avx_tag>::operator/=(const simd& x)
{
    *this = *this / x;
    return *this;
}

}}

#pragma warning(pop)