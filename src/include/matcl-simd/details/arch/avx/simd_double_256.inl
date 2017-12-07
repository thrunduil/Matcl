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

#include "matcl-simd/arch/avx/simd_double_256.h"

namespace matcl { namespace simd
{

//-------------------------------------------------------------------
//                          AVX DOUBLE
//-------------------------------------------------------------------

force_inline
simd<double, 256, avx_tag>::simd(Integer val)
    : data(_mm256_set1_pd(val)) 
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
simd<double, 256, avx_tag>::gather(const double* arr, const simd_128_int32& ind)
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
simd<double, 256, avx_tag>::gather(const double* arr, const simd_256_int64& ind)
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
