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

#include "matcl-simd/arch/sse/simd_128.h"

namespace matcl { namespace simd { namespace details
{

template<class Ret>
struct cast_to_double_impl
{
    using simd_type = simd<float,128, matcl::simd::sse_tag>;

    force_inline
    static Ret eval(const simd_type& s)
    {
        return Ret(s.cast_low_to_double(), s.cast_high_to_double());
    };
};

template<>
struct cast_to_double_impl<simd<double,256, matcl::simd::avx_tag>>
{
    using Ret       = simd<double,256, matcl::simd::avx_tag>;
    using simd_type = simd<float,128, matcl::simd::sse_tag>;

    force_inline
    static Ret eval(const simd_type& s)
    {
        return _mm256_cvtps_pd(s.data);
    };
};

}}};

namespace matcl { namespace simd
{

//-------------------------------------------------------------------
//                          SSE2 DOUBLE
//-------------------------------------------------------------------

force_inline
simd<double, 128, sse_tag>::simd(Integer val)
    : data(_mm_set1_pd(val)) 
{}

force_inline
simd<double, 128, sse_tag>::simd(float val)
    : data(_mm_set1_pd(val)) 
{}

force_inline
simd<double, 128, sse_tag>::simd(double val)
    : data(_mm_set1_pd(val)) 
{}

force_inline
simd<double, 128, sse_tag>::simd(double v1, double v2)
    : data(_mm_setr_pd(v1, v2)) 
{}

force_inline
simd<double, 128, sse_tag>::simd(const simd& lo, const simd& hi)
    : data(_mm_shuffle_pd(lo.data, hi.data, 0)) 
{};

force_inline
simd<double, 128, sse_tag>::simd(const impl_type& v)
    : data(v) 
{};

force_inline
simd<double, 128, sse_tag>::simd(const simd<double, 128, nosimd_tag>& s)
    : data(_mm_load_pd(s.data))
{};

force_inline
simd<double, 128, sse_tag> simd<double, 128, sse_tag>::broadcast(const double* arr)
{ 
    return _mm_load1_pd(arr);
};

force_inline
simd<double, 128, sse_tag> simd<double, 128, sse_tag>::broadcast(const double& arr)
{ 
    return _mm_load1_pd(&arr);
};

force_inline
simd<double, 128, sse_tag> simd<double, 128, sse_tag>::set_lower(double v)
{
    return _mm_set_sd(v);
};

force_inline
simd<float, 128, sse_tag>
simd<double, 128, sse_tag>::cast_to_float() const
{
    return _mm_cvtpd_ps(data);
};

force_inline
double simd<double, 128, sse_tag>::get(int pos) const
{ 
    return get_raw_ptr()[pos] ; 
};

force_inline
double simd<double, 128, sse_tag>::first() const
{ 
    return _mm_cvtsd_f64(data); 
};

force_inline
void simd<double, 128, sse_tag>::set(int pos, double val)
{ 
    get_raw_ptr()[pos] = val; 
};

force_inline
const double* simd<double, 128, sse_tag>::get_raw_ptr() const
{ 
    return reinterpret_cast<const double*>(&data); 
};

force_inline
double* simd<double, 128, sse_tag>::get_raw_ptr()
{ 
    return reinterpret_cast<double*>(&data); 
};

force_inline
simd<double, 128, sse_tag> simd<double, 128, sse_tag>::zero()
{
    impl_type data  = _mm_setzero_pd();
    return data;
}

force_inline
simd<double, 128, sse_tag> simd<double, 128, sse_tag>::minus_zero()
{
    return simd(-0.0);
}

force_inline
simd<double, 128, sse_tag> simd<double, 128, sse_tag>::one()
{
    return simd(1.0);
}

force_inline
simd<double, 128, sse_tag> simd<double, 128, sse_tag>::minus_one()
{
    return simd(-1.0);
}

force_inline simd<double, 128, sse_tag> 
simd<double, 128, sse_tag>::load(const double* arr, std::true_type aligned)
{
    (void)aligned;
    return _mm_load_pd(arr);
};

force_inline simd<double, 128, sse_tag> 
simd<double, 128, sse_tag>::load(const double* arr, std::false_type not_aligned)
{
    (void)not_aligned;
    return _mm_loadu_pd(arr);
};

force_inline void simd<double, 128, sse_tag>::store(double* arr, std::true_type aligned) const
{
    (void)aligned;
    _mm_store_pd(arr, data);
};

force_inline void simd<double, 128, sse_tag>::store(double* arr, std::false_type not_aligned) const
{
    (void)not_aligned;
    _mm_storeu_pd(arr, data);
};

template<int Step>
force_inline
void simd<double, 128, sse_tag>::scatter(double* arr) const
{
    //no scatter intrinsic
    _mm_storel_pd(arr + 0*Step, data);
    _mm_storeh_pd(arr + 1*Step, data);
};

force_inline simd<double, 128, sse_tag>& 
simd<double, 128, sse_tag>::operator+=(const simd& x)
{
    *this = *this + x;
    return *this;
}

force_inline simd<double, 128, sse_tag>& 
simd<double, 128, sse_tag>::operator-=(const simd& x)
{
    *this = *this - x;
    return *this;
}

force_inline simd<double, 128, sse_tag>& 
simd<double, 128, sse_tag>::operator*=(const simd& x)
{
    *this = *this * x;
    return *this;
}

force_inline simd<double, 128, sse_tag>& 
simd<double, 128, sse_tag>::operator/=(const simd& x)
{
    *this = *this / x;
    return *this;
}

//-------------------------------------------------------------------
//                          SSE2 SINGLE
//-------------------------------------------------------------------

force_inline
simd<float, 128, sse_tag>::simd(float val)
    : data(_mm_set1_ps(val)) 
{}

force_inline
simd<float, 128, sse_tag>::simd(float v1, float v2, float v3, float v4) 
    : data(_mm_setr_ps(v1, v2, v3, v4)) 
{}

force_inline
simd<float, 128, sse_tag>::simd(const impl_type& v)
    : data(v) 
{};

force_inline
simd<float, 128, sse_tag>::simd(const simd& lo, const simd& hi)
    : data(_mm_movelh_ps(lo.data, hi.data))
{};

force_inline
simd<float, 128, sse_tag>::simd(const simd<float, 128, nosimd_tag>& s)
    : data(_mm_load_ps(s.data))
{};

force_inline
simd<float, 128, sse_tag> simd<float, 128, sse_tag>::broadcast(const float* arr)
{ 
    return _mm_load1_ps(arr);
};

force_inline
simd<float, 128, sse_tag> simd<float, 128, sse_tag>::broadcast(const float& arr)
{ 
    return _mm_load1_ps(&arr);
};

force_inline
float simd<float, 128, sse_tag>::get(int pos) const
{ 
    return get_raw_ptr()[pos] ; 
};

force_inline
float simd<float, 128, sse_tag>::first() const
{ 
    return _mm_cvtss_f32(data); 
};

force_inline
void simd<float, 128, sse_tag>::set(int pos, float val)
{ 
    get_raw_ptr()[pos] = val; 
};

force_inline
const float* simd<float, 128, sse_tag>::get_raw_ptr() const
{ 
    return reinterpret_cast<const float*>(&data); 
};

force_inline
float* simd<float, 128, sse_tag>::get_raw_ptr()
{ 
    return reinterpret_cast<float*>(&data); 
};

force_inline
simd<float, 128, sse_tag> simd<float, 128, sse_tag>::zero()
{
    impl_type data   = _mm_setzero_ps();
    return data;
}

force_inline
simd<float, 128, sse_tag> simd<float, 128, sse_tag>::minus_zero()
{
    return simd(-0.0f);
}

force_inline
simd<float, 128, sse_tag> simd<float, 128, sse_tag>::one()
{
    return simd(1.0f);
}

force_inline
simd<float, 128, sse_tag> simd<float, 128, sse_tag>::minus_one()
{
    return simd(-1.0f);
}

force_inline simd<float, 128, sse_tag> 
simd<float, 128, sse_tag>::load(const float* arr, std::true_type aligned)
{
    (void)aligned;
    return _mm_load_ps(arr);
};

force_inline simd<float, 128, sse_tag> 
simd<float, 128, sse_tag>::load(const float* arr, std::false_type not_aligned)
{
    (void)not_aligned;
    return _mm_loadu_ps(arr);
};

force_inline simd<float, 128, sse_tag> 
simd<float, 128, sse_tag>::set_lower(float v)
{
    return _mm_set_ss(v);
};

force_inline void simd<float, 128, sse_tag>::store(float* arr, std::true_type aligned) const
{
    (void)aligned;
    _mm_store_ps(arr, data);
};


force_inline void 
simd<float, 128, sse_tag>::store(float* arr, std::false_type not_aligned) const
{
    (void)not_aligned;
    _mm_storeu_ps(arr, data);
};

force_inline simd<double, 128, sse_tag>
simd<float, 128, sse_tag>::cast_low_to_double() const
{
    return _mm_cvtps_pd(data);
};

force_inline
simd<double, 128, sse_tag>
simd<float, 128, sse_tag>::cast_high_to_double() const
{
    __m128 hi = _mm_movehl_ps(data, data);
    return _mm_cvtps_pd(hi);
};

force_inline
typename simd<float, 128, sse_tag>::simd_double_2
simd<float, 128, sse_tag>::cast_to_double() const
{
    return details::cast_to_double_impl<simd_double_2>::eval(data);    
};

template<int Step>
force_inline
void simd<float, 128, sse_tag>::scatter(float* arr) const
{
    const float* this_data  = get_raw_ptr();

    //no scatter intrinsic
    arr[0*Step]    = this_data[0];
    arr[1*Step]    = this_data[1];
    arr[2*Step]    = this_data[2];
    arr[3*Step]    = this_data[3];
};

force_inline simd<float, 128, sse_tag>& 
simd<float, 128, sse_tag>::operator+=(const simd& x)
{
    *this = *this + x;
    return *this;
}

force_inline simd<float, 128, sse_tag>& 
simd<float, 128, sse_tag>::operator-=(const simd& x)
{
    *this = *this - x;
    return *this;
}

force_inline simd<float, 128, sse_tag>& 
simd<float, 128, sse_tag>::operator*=(const simd& x)
{
    *this = *this * x;
    return *this;
}

force_inline simd<float, 128, sse_tag>& 
simd<float, 128, sse_tag>::operator/=(const simd& x)
{
    *this = *this / x;
    return *this;
}

}}
