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

#include "matcl-simd/arch/reg_128/simd_128.h"
#include "matcl-simd/simd/simd_128_func.h"

namespace matcl { namespace simd
{

//-------------------------------------------------------------------
//                          SSE2 DOUBLE
//-------------------------------------------------------------------

force_inline
simd<double, reg_128>::simd(Integer val)
    : data(_mm_set1_pd(val)) 
{}

force_inline
simd<double, reg_128>::simd(float val)
    : data(_mm_set1_pd(val)) 
{}

force_inline
simd<double, reg_128>::simd(double val)
    : data(_mm_set1_pd(val)) 
{}

force_inline
simd<double, reg_128>::simd(double v1, double v2)
    : data(_mm_setr_pd(v1, v2)) 
{}

force_inline
simd<double, reg_128>::simd(const impl_type& v)
    : data(v) 
{};

force_inline
simd<double, reg_128> simd<double, reg_128>::broadcast(const double* arr)
{ 
    return _mm_load1_pd(arr);
};

force_inline
simd<double, reg_128> simd<double, reg_128>::broadcast(const double& arr)
{ 
    return _mm_load1_pd(&arr);
};

force_inline
simd<double, reg_128> simd<double, reg_128>::set_lower(double v)
{
    return _mm_set_sd(v);
};

force_inline
double simd<double, reg_128>::get(int pos) const
{ 
    return data.m128d_f64[pos] ; 
};

force_inline
simd<double, reg_128> simd<double, reg_128>::zero()
{
    impl_type data  = _mm_setzero_pd();
    return data;
}

force_inline simd<double,reg_128> 
simd<double,reg_128>::load(const double* arr, std::true_type aligned)
{
    (void)aligned;
    return _mm_load_pd(arr);
};

force_inline simd<double,reg_128> 
simd<double,reg_128>::load(const double* arr, std::false_type not_aligned)
{
    (void)not_aligned;
    return _mm_loadu_pd(arr);
};

force_inline simd<double,reg_128> 
simd<double,reg_128>::load_reverse(const double* arr, std::true_type aligned)
{
    (void)aligned;
    static const int off        = vector_size - 1;
    simd<double,reg_128> ret    = _mm_loadr_pd(arr-off);
    return ret;
}

force_inline simd<double,reg_128> 
simd<double,reg_128>::load_reverse(const double* arr, std::false_type not_aligned)
{
    static const int off        = vector_size - 1;
    simd<double,reg_128> ret    = load(arr - off, not_aligned);
    return reverse(ret);
};

force_inline void simd<double,reg_128>::store(double* arr, std::true_type aligned) const
{
    (void)aligned;
    _mm_store_pd(arr, data);
};

force_inline void simd<double,reg_128>::store(double* arr, std::false_type not_aligned) const
{
    (void)not_aligned;
    _mm_storeu_pd(arr, data);
};

force_inline void simd<double,reg_128>::store_reverse(double* arr, std::true_type aligned) const
{
    (void)aligned;
    static const int off = vector_size - 1;

    _mm_storer_pd(arr - off, data);
}

force_inline void simd<double,reg_128>::store_reverse(double* arr, std::false_type not_aligned) const
{
    (void)not_aligned;
    arr[-1]     = data.m128d_f64[1];
    arr[0]      = data.m128d_f64[0];
};

template<int Step>
force_inline
void simd<double,reg_128>::scatter(double* arr) const
{
    //no scatter intrinsic
    _mm_storel_pd(arr + 0*Step, data);
    _mm_storeh_pd(arr + 1*Step, data);
};

//-------------------------------------------------------------------
//                          SSE2 SINGLE
//-------------------------------------------------------------------
force_inline
simd<float, reg_128>::simd(float val)
    : data(_mm_set1_ps(val)) 
{}

force_inline
simd<float, reg_128>::simd(float v1, float v2, float v3, float v4) 
    : data(_mm_setr_ps(v1, v2, v3, v4)) 
{}

force_inline
simd<float, reg_128>::simd(const impl_type& v)
    : data(v) 
{};

force_inline
simd<float, reg_128> simd<float, reg_128>::broadcast(const float* arr)
{ 
    return _mm_load1_ps(arr);
};

force_inline
simd<float, reg_128> simd<float, reg_128>::broadcast(const float& arr)
{ 
    return _mm_load1_ps(&arr);
};

force_inline
float simd<float, reg_128>::get(int pos) const
{ 
    return data.m128_f32[pos] ; 
};

force_inline
simd<float, reg_128> simd<float, reg_128>::zero()
{
    impl_type data   = _mm_setzero_ps();
    return data;
}

force_inline simd<float,reg_128> 
simd<float,reg_128>::load(const float* arr, std::true_type aligned)
{
    (void)aligned;
    return _mm_load_ps(arr);
};

force_inline simd<float,reg_128> 
simd<float,reg_128>::load(const float* arr, std::false_type not_aligned)
{
    (void)not_aligned;
    return _mm_loadu_ps(arr);
};

force_inline simd<float,reg_128> 
simd<float,reg_128>::load_reverse(const float* arr, std::true_type aligned)
{
    (void)aligned;
    static const int off        = vector_size - 1;
    simd<float,reg_128> ret     = _mm_loadr_ps(arr-off);
    return ret;
}

force_inline simd<float,reg_128> 
simd<float,reg_128>::load_reverse(const float* arr, std::false_type not_aligned)
{
    static const int off    = vector_size - 1;
    simd<float,reg_128> ret = load(arr - off, not_aligned);
    return reverse(ret);
};

force_inline simd<float,reg_128> 
simd<float,reg_128>::set_lower(float v)
{
    return _mm_set_ss(v);
};

force_inline void simd<float,reg_128>::store(float* arr, std::true_type aligned) const
{
    (void)aligned;
    _mm_store_ps(arr, data);
};


force_inline void simd<float,reg_128>::store(float* arr, std::false_type not_aligned) const
{
    (void)not_aligned;
    _mm_storeu_ps(arr, data);
};

force_inline void simd<float,reg_128>::store_reverse(float* arr, std::true_type aligned) const
{
    (void)aligned;
    static const int off = vector_size - 1;

    _mm_storer_ps(arr - off, data);
}

force_inline void simd<float,reg_128>::store_reverse(float* arr, std::false_type not_aligned) const
{
    static const int off        = vector_size - 1;
    simd<float,reg_128> rev     = reverse(*this);
    rev.store(arr - off, not_aligned);
};

template<int Step>
force_inline
void simd<float,reg_128>::scatter(float* arr) const
{
    //no scatter intrinsic
    arr + 0*Step    = data.m128_f32[0];
    arr + 1*Step    = data.m128_f32[1];
    arr + 2*Step    = data.m128_f32[2];
    arr + 3*Step    = data.m128_f32[3];
};

}}
