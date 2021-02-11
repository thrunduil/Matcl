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

#include "matcl-simd/arch/sse/simd_double_128.h"
#include "matcl-simd/details/arch/sse/func/missing_intrinsics.h"

#include <algorithm>

#pragma warning(push)
#pragma warning(disable : 4127) // conditional expression is constant

namespace matcl { namespace simd
{

//-------------------------------------------------------------------
//                          SSE2 DOUBLE
//-------------------------------------------------------------------

force_inline
simd<double, 128, sse_tag>::simd(int32_t val)
    : data(_mm_set1_pd(val)) 
{}

force_inline
simd<double, 128, sse_tag>::simd(int64_t val)
    : data(_mm_set1_pd(double(val))) 
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
simd<double, 128, sse_tag>::simd(const simd<double, 128, scalar_sse_tag>& s)
    :data(_mm_shuffle_pd(s.data, s.data, 0))
{};

force_inline
simd<double, 128, sse_tag>::simd(const simd<double, 128, scalar_nosimd_tag>& s)
    : simd(s.first())
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
simd<double, 128, sse_tag>::convert_to_float() const
{
    return _mm_cvtpd_ps(data);
};

force_inline simd<int32_t, 128, sse_tag>
simd<double, 128, sse_tag>::convert_to_int32() const
{
    return _mm_cvtpd_epi32(data);
};

force_inline simd<int64_t, 128, sse_tag> 
simd<double, 128, sse_tag>::convert_to_int64() const
{
    simd<int64_t, 128, sse_tag> ret;

    int64_t* ret_ptr    = ret.get_raw_ptr();
    const double* ptr   = this->get_raw_ptr();

    for (int i = 0; i < vector_size; ++i)
        ret_ptr[i]      = scalar_func::convert_double_int64(ptr[i]);
    
    return ret;
};

force_inline simd<float, 128, sse_tag>
simd<double, 128, sse_tag>::reinterpret_as_float() const
{
    return _mm_castpd_ps(data);
}

force_inline simd<int32_t, 128, sse_tag>
simd<double, 128, sse_tag>::reinterpret_as_int32() const
{
    return _mm_castpd_si128(data);
}

force_inline simd<int64_t, 128, sse_tag>
simd<double, 128, sse_tag>::reinterpret_as_int64() const
{
    return _mm_castpd_si128(data);
}

force_inline simd<double, 128, sse_tag> 
simd<double, 128, sse_tag>::extract_low() const
{
    return *this;
}

force_inline simd<double, 128, sse_tag> 
simd<double, 128, sse_tag>::extract_high() const
{
    return missing::mm_movehl_pd(data, data);
}

template<int I1, int I2>
force_inline simd<double, 128, sse_tag>  
simd<double, 128, sse_tag>::select() const
{
    static_assert(I1 >= 0 && I1 <= 1, "invalid index in select function");
    static_assert(I2 >= 0 && I2 <= 1, "invalid index in select function");

    if (I1 == 0 && I2 == 1)
        return data;

    static const int ind0   = I1 + (I2 << 1);
    static const int ind    = (ind0 < 0) ? 0 : ind0;

    return _mm_shuffle_pd(data, data, ind);
}

template<int I1, int I2>
force_inline simd<double, 128, sse_tag>  
simd<double, 128, sse_tag>::combine(const simd& x, const simd& y)
{
    static_assert(I1 >= 0 && I1 <= 3, "invalid index in select function");
    static_assert(I2 >= 0 && I2 <= 3, "invalid index in select function");

    if (I1 <= 1 && I2 <= 1)
    {
        //only elements from x vector

        static const int I1_s   = std::min(I1, 1);
        static const int I2_s   = std::min(I2, 1);

        return x.select<I1_s, I2_s>();
    }

    if (I1 >= 2 && I2 >= 2)
    {
        //only elements from y vector

        static const int I1_s   = std::max(I1 - 2, 0);
        static const int I2_s   = std::max(I2 - 2, 0);

        return y.select<I1_s, I2_s>();
    }

    if (I1 <= 1)
    {
        // first element from x, second from y

        static const int I1_s   = std::min(I1, 1);
        static const int I2_s   = std::max(I2, 2);

        static const int ind   = I1_s + ((I2_s - 2) << 1);
        return _mm_shuffle_pd(x.data, y.data, ind);
    }
    else
    {
        // first element from y, second from x

        static const int I1_s   = std::max(I1, 2);
        static const int I2_s   = std::min(I2, 1);

        static const int ind    = (I1_s - 2) + (I2_s << 1);
        return _mm_shuffle_pd(y.data, x.data, ind);
    };
}

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

force_inline simd<double, 128, sse_tag> 
simd<double, 128, sse_tag>::gather(const double* arr, const simd_int32& ind)
{
    #if MATCL_ARCHITECTURE_HAS_AVX2
        return _mm_i32gather_pd(arr, ind.data, 8);
    #else
        simd res;
        double* res_ptr         = res.get_raw_ptr();
        const int32_t* ind_ptr  = ind.get_raw_ptr();

        for (int i = 0; i < vector_size; ++i)
            res_ptr[i]  = arr[ind_ptr[i]];

        return res;
    #endif    
}

force_inline simd<double, 128, sse_tag> 
simd<double, 128, sse_tag>::gather(const double* arr, const simd_int64& ind)
{
    #if MATCL_ARCHITECTURE_HAS_AVX2
        return _mm_i64gather_pd(arr, ind.data, 8);
    #else
        simd res;
        double* res_ptr         = res.get_raw_ptr();
        const int64_t* ind_ptr  = ind.get_raw_ptr();

        for (int i = 0; i < vector_size; ++i)
            res_ptr[i]  = arr[ind_ptr[i]];

        return res;
    #endif
}

force_inline void 
simd<double, 128, sse_tag>::store(double* arr, std::true_type aligned) const
{
    (void)aligned;
    _mm_store_pd(arr, data);
};

force_inline void 
simd<double, 128, sse_tag>::store(double* arr, std::false_type not_aligned) const
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

}}

#pragma warning(pop)