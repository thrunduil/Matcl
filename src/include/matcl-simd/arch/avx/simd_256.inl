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
simd<double, 256, avx_tag>::simd(const double& val)
    : data(_mm256_broadcast_sd(&val)) 
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
    : data(_mm256_broadcast_pd(&lo_hi.data))
{}

force_inline
simd<double, 256, avx_tag>::simd(const impl_type& v)
    : data(v)
{};

force_inline
double simd<double, 256, avx_tag>::get(int pos) const
{ 
    return data.m256d_f64[pos] ; 
};

template<int Pos>
force_inline
double simd<double, 256, avx_tag>::get() const
{ 
    return data.m256d_f64[Pos]; 
};

force_inline
void simd<double, 256, avx_tag>::set(int pos, double val)
{ 
    data.m256d_f64[pos] = val; 
};

template<int Pos>
force_inline
void simd<double, 256, avx_tag>::set(double val)
{ 
    data.m256d_f64[Pos] = val; 
};

force_inline simd<double, 256, avx_tag>::simd_half
simd<double, 256, avx_tag>::extract_low() const
{
    return _mm256_extractf128_pd(data, 0);
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

template<int Step>
force_inline
void simd<double, 256, avx_tag>::scatter(double* arr) const
{
    //no scatter intrinsic
    arr[0*Step] = data.m256d_f64[0];
    arr[1*Step] = data.m256d_f64[1];
    arr[2*Step] = data.m256d_f64[2];
    arr[3*Step] = data.m256d_f64[3];
};

//-------------------------------------------------------------------
//                          AVX SINGLE
//-------------------------------------------------------------------
force_inline
simd<float, 256, avx_tag>::simd(const float& val) 
    : data(_mm256_broadcast_ss(&val)) 
{}

force_inline
simd<float, 256, avx_tag>::simd(const simd_half& lo_hi)
    : data(_mm256_broadcast_ps(&lo_hi.data))
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
float simd<float, 256, avx_tag>::get(int pos) const  
{ 
    return data.m256_f32[pos] ; 
};

template<int Pos>
force_inline
float simd<float, 256, avx_tag>::get() const
{ 
    return data.m256_f32[Pos]; 
};

force_inline
void simd<float, 256, avx_tag>::set(int pos, float val)
{ 
    data.m256_f32[pos] = val; 
};

template<int Pos>
force_inline
void simd<float, 256, avx_tag>::set(float val)
{ 
    data.m256_f32[Pos] = val; 
};

force_inline simd<float, 256, avx_tag>::simd_half
simd<float, 256, avx_tag>::extract_low() const
{
    return _mm256_extractf128_ps(data, 0);
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

template<int Step>
force_inline
void simd<float, 256, avx_tag>::scatter(float* arr) const
{
    //no scatter intrinsic
    arr[0*Step] = data.m256_f32[0];
    arr[1*Step] = data.m256_f32[1];
    arr[2*Step] = data.m256_f32[2];
    arr[3*Step] = data.m256_f32[3];
    arr[4*Step] = data.m256_f32[4];
    arr[5*Step] = data.m256_f32[5];
    arr[6*Step] = data.m256_f32[6];
    arr[7*Step] = data.m256_f32[7];
};

}}
