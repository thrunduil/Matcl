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

#include "matcl-simd/arch/reg_256/simd_256.h"
#include "matcl-simd/simd/simd_256_func.h"

namespace matcl { namespace simd
{

//-------------------------------------------------------------------
//                          AVX DOUBLE
//-------------------------------------------------------------------

force_inline
simd<double, reg_256>::simd(Integer val)
    : data(_mm256_set1_pd(val)) 
{}

force_inline
simd<double, reg_256>::simd(float val)
    : data(_mm256_set1_pd(val)) 
{}

force_inline
simd<double, reg_256>::simd(const double& val)
    : data(_mm256_broadcast_sd(&val)) 
{}

force_inline
simd<double, reg_256>::simd(double v1, double v2, double v3, double v4) 
    : data(_mm256_setr_pd(v1, v2, v3, v4)) 
{}

force_inline
simd<double, reg_256>::simd(simd_half lo, simd_half hi)
    : data(_mm256_setr_m128d(lo.data, hi.data))
{}

force_inline
simd<double, reg_256>::simd(const simd_half& lo_hi)
    : data(_mm256_broadcast_pd(&lo_hi.data))
{}

force_inline
simd<double, reg_256>::simd(const impl_type& v)
    : data(v)
{};

force_inline
double simd<double, reg_256>::get(int pos) const
{ 
    return data.m256d_f64[pos] ; 
};

force_inline simd<double, reg_256>::simd_half
simd<double, reg_256>::extract_low() const
{
    return _mm256_extractf128_pd(data, 0);
}

force_inline simd<double, reg_256>::simd_half
simd<double, reg_256>::extract_high() const
{
    return _mm256_extractf128_pd(data, 1);
}

force_inline
simd<double, reg_256> simd<double, reg_256>::zero()
{
    impl_type data  = _mm256_setzero_pd();
    return data;
}

force_inline simd<double,reg_256> 
simd<double,reg_256>::load(const double* arr, std::true_type aligned)
{
    (void)aligned;
    return _mm256_load_pd(arr);
};

force_inline simd<double,reg_256> 
simd<double,reg_256>::load(const double* arr, std::false_type not_aligned)
{
    (void)not_aligned;
    return _mm256_loadu_pd(arr);
};

force_inline simd<double,reg_256> 
simd<double,reg_256>::broadcast(const double* arr)
{
    return _mm256_broadcast_sd(arr);
};

force_inline simd<double,reg_256> 
simd<double,reg_256>::broadcast(const double& arr)
{
    return _mm256_broadcast_sd(&arr);
};

force_inline simd<double,reg_256> 
simd<double,reg_256>::load_reverse(const double* arr, std::true_type aligned)
{
    static const int off        = vector_size - 1;
    simd<double,reg_256> ret    = load(arr - off, aligned);
    return reverse(ret);
}

force_inline simd<double,reg_256> 
simd<double,reg_256>::load_reverse(const double* arr, std::false_type not_aligned)
{
    static const int off        = vector_size - 1;
    simd<double,reg_256> ret    = load(arr - off, not_aligned);
    return reverse(ret);
};

force_inline void simd<double,reg_256>::store(double* arr, std::true_type aligned) const
{
    (void)aligned;
    _mm256_store_pd(arr, data);
};

force_inline void simd<double,reg_256>::store(double* arr, std::false_type not_aligned) const
{
    (void)not_aligned;
    _mm256_storeu_pd(arr, data);
};

force_inline void simd<double,reg_256>::store_reverse(double* arr, std::true_type aligned) const
{
    static const int off        = vector_size - 1;
    simd<double,reg_256> rev    = reverse(*this);
    rev.store(arr - off, aligned);
}

force_inline void simd<double,reg_256>::store_reverse(double* arr, std::false_type not_aligned) const
{
    static const int off        = vector_size - 1;
    simd<double,reg_256> rev    = reverse(*this);
    rev.store(arr - off, not_aligned);
};

template<int Step>
force_inline
void simd<double,reg_256>::scatter(double* arr) const
{
    //no scatter intrinsic

    __m128d a = extract_low();
    __m128d b = extract_high();

    _mm_storel_pd(arr + 0*Step, a);
    _mm_storeh_pd(arr + 1*Step, a);
    _mm_storel_pd(arr + 2*Step, b);
    _mm_storeh_pd(arr + 3*Step, b);
};

//-------------------------------------------------------------------
//                          AVX SINGLE
//-------------------------------------------------------------------
force_inline
simd<float, reg_256>::simd(const float& val) 
    : data(_mm256_broadcast_ss(&val)) 
{}

force_inline
simd<float, reg_256>::simd(const simd_half& lo_hi)
    : data(_mm256_broadcast_ps(&lo_hi.data))
{}

force_inline
simd<float, reg_256>::simd(float v1, float v2, float v3, float v4, float v5, float v6, float v7, float v8) 
    : data(_mm256_setr_ps(v1, v2, v3, v4, v5, v6, v7, v8)) 
{}

force_inline
simd<float, reg_256>::simd(const impl_type& v)
    : data(v) 
{};

force_inline
simd<float, reg_256>::simd(simd_half lo, simd_half hi)
    : data(_mm256_setr_m128(lo.data, hi.data))
{}

force_inline
float simd<float, reg_256>::get(int pos) const  
{ 
    return data.m256_f32[pos] ; 
};

force_inline simd<float, reg_256>::simd_half
simd<float, reg_256>::extract_low() const
{
    return _mm256_extractf128_ps(data, 0);
}

force_inline simd<float, reg_256>::simd_half
simd<float, reg_256>::extract_high() const
{
    return _mm256_extractf128_ps(data, 1);
}

force_inline
simd<float, reg_256> simd<float, reg_256>::zero()
{
    impl_type data  = _mm256_setzero_ps();
    return data;
}

force_inline simd<float,reg_256> 
simd<float,reg_256>::load(const float* arr, std::true_type aligned)
{
    (void)aligned;
    return _mm256_load_ps(arr);
};

force_inline simd<float,reg_256> 
simd<float,reg_256>::load(const float* arr, std::false_type not_aligned)
{
    (void)not_aligned;
    return _mm256_loadu_ps(arr);
};

force_inline simd<float,reg_256> 
simd<float,reg_256>::broadcast(const float* arr)
{
    return _mm256_broadcast_ss(arr);
};

force_inline simd<float,reg_256> 
simd<float,reg_256>::broadcast(const float& arr)
{
    return _mm256_broadcast_ss(&arr);
};

force_inline simd<float,reg_256> 
simd<float,reg_256>::load_reverse(const float* arr, std::true_type aligned)
{
    static const int off        = vector_size - 1;
    simd<float,reg_256> ret     = load(arr - off, aligned);
    return reverse(ret);
}

force_inline simd<float,reg_256> 
simd<float,reg_256>::load_reverse(const float* arr, std::false_type not_aligned)
{
    static const int off        = vector_size - 1;
    simd<float,reg_256> ret     = load(arr - off, not_aligned);
    return reverse(ret);
};

force_inline void simd<float,reg_256>::store(float* arr, std::true_type aligned) const
{
    (void)aligned;
    _mm256_store_ps(arr, data);
};

force_inline void simd<float,reg_256>::store(float* arr, std::false_type not_aligned) const
{
    (void)not_aligned;
    _mm256_storeu_ps(arr, data);
};

force_inline void simd<float,reg_256>::store_reverse(float* arr, std::true_type aligned) const
{
    static const int off        = vector_size - 1;
    simd<float,reg_256> rev     = reverse(*this);
    rev.store(arr - off, aligned);
}

force_inline void simd<float,reg_256>::store_reverse(float* arr, std::false_type not_aligned) const
{
    static const int off        = vector_size - 1;
    simd<float,reg_256> rev     = reverse(*this);
    rev.store(arr - off, not_aligned);
};

template<int Step>
force_inline
void simd<float,reg_256>::scatter(float* arr) const
{
    //no scatter intrinsic
    arr + 0*Step    = data.m256_f32[0];
    arr + 1*Step    = data.m256_f32[1];
    arr + 2*Step    = data.m256_f32[2];
    arr + 3*Step    = data.m256_f32[3];
    arr + 4*Step    = data.m256_f32[4];
    arr + 5*Step    = data.m256_f32[5];
    arr + 6*Step    = data.m256_f32[6];
    arr + 7*Step    = data.m256_f32[7];
};

}}
