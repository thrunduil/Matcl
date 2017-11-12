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

#include "matcl-simd/arch/avx/simd_256.h"

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
simd<double, 256, avx_tag>::broadcast(const double* arr)
{
    return _mm256_broadcast_sd(arr);
};

force_inline simd<double, 256, avx_tag> 
simd<double, 256, avx_tag>::broadcast(const double& arr)
{
    return _mm256_broadcast_sd(&arr);
};

force_inline void simd<double, 256, avx_tag>::store(double* arr, std::true_type aligned) const
{
    (void)aligned;
    _mm256_store_pd(arr, data);
};

force_inline void simd<double, 256, avx_tag>::store(double* arr, std::false_type not_aligned) const
{
    (void)not_aligned;
    _mm256_storeu_pd(arr, data);
};

force_inline simd<float, 128, sse_tag> 
simd<double, 256, avx_tag>::cast_to_float() const
{
    return _mm256_cvtpd_ps(data);
};

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
simd<float, 256, avx_tag>::cast_low_to_double() const
{
    simd_half lo    = this->extract_low();
    return _mm256_cvtps_pd(lo.data);
};

force_inline
simd<double, 256, avx_tag>
simd<float, 256, avx_tag>::cast_high_to_double() const
{
    simd_half hi    = this->extract_high();
    return _mm256_cvtps_pd(hi.data);
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
