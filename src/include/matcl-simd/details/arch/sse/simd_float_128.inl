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

#include "matcl-simd/arch/sse/simd_float_128.h"
#include "matcl-simd/details/arch/sse/helpers.h"

#pragma warning(push)
#pragma warning(disable : 4127) // conditional expression is constant

namespace matcl { namespace simd { namespace details
{

template<class Ret>
struct cast_to_double_impl
{
    using simd_type = simd<float,128, matcl::simd::sse_tag>;

    force_inline
    static Ret eval(const simd_type& s)
    {
        return Ret(s.convert_low_to_double(), s.convert_high_to_double());
    };
};

#if MATCL_ARCHITECTURE_HAS_AVX
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
#endif

}}};

namespace matcl { namespace simd
{

//-------------------------------------------------------------------
//                          SSE SINGLE
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
simd<float, 128, sse_tag>::simd(const simd<float, 128, scalar_sse_tag>& s)
    : data(_mm_shuffle_ps(s.data, s.data, 0))
{};

force_inline
simd<float, 128, sse_tag>::simd(const simd<float, 128, scalar_nosimd_tag>& s)
    : simd(s.first())
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
simd<float, 128, sse_tag>::gather(const float* arr, const simd_int32& ind)
{
    #if MATCL_ARCHITECTURE_HAS_AVX2
        return _mm_i32gather_ps(arr, ind.data, 4);
    #else
        simd res;
        float* res_ptr          = res.get_raw_ptr();
        const int32_t* ind_ptr  = ind.get_raw_ptr();

        for (int i = 0; i < vector_size; ++i)
            res_ptr[i]  = arr[ind_ptr[i]];

        return res;
    #endif    
}

force_inline simd<float, 128, sse_tag> 
simd<float, 128, sse_tag>::gather(const float* arr, const simd_int64& ind)
{
    #if MATCL_ARCHITECTURE_HAS_AVX2
        return simd(_mm_i64gather_ps(arr, ind.data, 4));
    #else
        simd res;
        float* res_ptr          = res.get_raw_ptr();
        const int64_t* ind_ptr  = ind.get_raw_ptr();

        for (int i = 0; i < vector_size/2; ++i)
            res_ptr[i]  = arr[ind_ptr[i]];

        return simd(res, simd::zero());
    #endif        
}

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
simd<float, 128, sse_tag>::convert_low_to_double() const
{
    return _mm_cvtps_pd(data);
};

force_inline simd<double, 128, sse_tag>
simd<float, 128, sse_tag>::convert_high_to_double() const
{
    __m128 hi = _mm_movehl_ps(data, data);
    return _mm_cvtps_pd(hi);
};

force_inline
typename simd<float, 128, sse_tag>::simd_double_2
simd<float, 128, sse_tag>::convert_to_double() const
{
    return details::cast_to_double_impl<simd_double_2>::eval(data);    
};

force_inline
simd<int64_t, 128, sse_tag>
simd<float, 128, sse_tag>::convert_low_to_int64() const
{
    //SIMD intrinsic not available
    simd<int64_t, 128, sse_tag> ret;

    const float* ptr    = this->get_raw_ptr();
    int64_t* ret_ptr    = ret.get_raw_ptr();

    ret_ptr[0]  = scalar_func::convert_float_int64(ptr[0]);
    ret_ptr[1]  = scalar_func::convert_float_int64(ptr[1]);

    return ret;
};

force_inline
simd<int64_t, 128, sse_tag>
simd<float, 128, sse_tag>::convert_high_to_int64() const
{
    //SIMD intrinsic not available
    simd<int64_t, 128, sse_tag> ret;

    const float* ptr    = this->get_raw_ptr();
    int64_t* ret_ptr    = ret.get_raw_ptr();

    ret_ptr[0]  = scalar_func::convert_float_int64(ptr[2]);
    ret_ptr[1]  = scalar_func::convert_float_int64(ptr[3]);

    return ret;

};

force_inline
typename simd<float, 128, sse_tag>::simd_int64_2
simd<float, 128, sse_tag>::convert_to_int64() const
{
    //SIMD intrinsic not available
    simd_int64_2 ret;

    const float* ptr    = this->get_raw_ptr();
    int64_t* ret_ptr    = ret.get_raw_ptr();

    ret_ptr[0]  = scalar_func::convert_float_int64(ptr[0]);
    ret_ptr[1]  = scalar_func::convert_float_int64(ptr[1]);
    ret_ptr[2]  = scalar_func::convert_float_int64(ptr[2]);
    ret_ptr[3]  = scalar_func::convert_float_int64(ptr[3]);

    return ret;
};

force_inline simd<double, 128, sse_tag>
simd<float, 128, sse_tag>::reinterpret_as_double() const
{
    return _mm_castps_pd(data);
}

force_inline simd<int32_t, 128, sse_tag>
simd<float, 128, sse_tag>::reinterpret_as_int32() const
{
    return _mm_castps_si128(data);
}

force_inline simd<int64_t, 128, sse_tag>
simd<float, 128, sse_tag>::reinterpret_as_int64() const
{
    return _mm_castps_si128(data);
}

force_inline simd<int32_t, 128, sse_tag>
simd<float, 128, sse_tag>::convert_to_int32() const
{
    return _mm_cvtps_epi32(data);
};

force_inline simd<float, 128, sse_tag> 
simd<float, 128, sse_tag>::extract_low() const
{
    return *this;
}

force_inline simd<float, 128, sse_tag> 
simd<float, 128, sse_tag>::extract_high() const
{
    return _mm_movehl_ps(data, data);
}

template<int I1, int I2, int I3, int I4>
force_inline simd<float, 128, sse_tag>  
simd<float, 128, sse_tag>::select() const
{
    static_assert(I1 >= 0 && I1 <= 3, "invalid index in select function");
    static_assert(I2 >= 0 && I2 <= 3, "invalid index in select function");
    static_assert(I3 >= 0 && I3 <= 3, "invalid index in select function");
    static_assert(I4 >= 0 && I4 <= 3, "invalid index in select function");

    if (I1 == 0 && I2 == 1 && I3 == 2 && I4 == 3)
        return data;

    static const int ind    = I1 + (I2 << 2) + (I3 << 4) + (I4 << 6);
    return _mm_shuffle_ps(data, data, ind);
}

template<int I1, int I2, int I3, int I4>
force_inline simd<float, 128, sse_tag>  
simd<float, 128, sse_tag>::combine(const simd& x, const simd& y)
{
    static_assert(I1 >= 0 && I1 <= 7, "invalid index in select function");
    static_assert(I2 >= 0 && I2 <= 7, "invalid index in select function");
    static_assert(I3 >= 0 && I3 <= 7, "invalid index in select function");
    static_assert(I4 >= 0 && I4 <= 7, "invalid index in select function");

    if (I1 <= 3 && I2 <= 3 && I3 <= 3 && I4 <= 3)
    {
        //only elements from x vector

        static const int I1_s   = std::min(I1, 3);
        static const int I2_s   = std::min(I2, 3);
        static const int I3_s   = std::min(I3, 3);
        static const int I4_s   = std::min(I4, 3);

        return x.select<I1_s, I2_s, I3_s, I4_s>();
    }

    if (I1 >= 4 && I2 >= 4 && I3 >= 4 && I4 >= 4)
    {
        //only elements from y vector

        static const int I1_s   = std::max(I1 - 4, 0);
        static const int I2_s   = std::max(I2 - 4, 0);
        static const int I3_s   = std::max(I3 - 4, 0);
        static const int I4_s   = std::max(I4 - 4, 0);

        return y.select<I1_s, I2_s, I3_s, I4_s>();
    }

    if (I1 == 0 && I2 == 4 && I3 == 1 && I4 == 5)
        return _mm_unpacklo_ps (x.data, y.data);

    if (I1 == 4 && I2 == 0 && I3 == 5 && I4 == 1)
        return _mm_unpacklo_ps (y.data, x.data);

    if (I1 == 2 && I2 == 6 && I3 == 3 && I4 == 7)
        return _mm_unpackhi_ps (x.data, y.data);

    if (I1 == 6 && I2 == 2 && I3 == 7 && I4 == 3)
        return _mm_unpackhi_ps (y.data, x.data);

    if (I1 <= 3 && I2 <= 3 && I3 >= 4 && I4 >= 4)
    {
        // first two elements from x, last two from y

        static const int I3_s   = std::max(I3 - 4, 0);
        static const int I4_s   = std::max(I4 - 4, 0);

        static const int ind0   = I1 + (I2 << 2) + (I3_s << 4) + (I4_s << 6);
        static const int ind    = std::min(std::max(ind0, 0), 255);
        return _mm_shuffle_ps(x.data, y.data, ind);
    }

    if (I1 >= 4 && I2 >= 4 && I3 <= 3 && I4 <= 3)
    {
        // first two elements from y, last two from x

        static const int I1_s   = std::max(I1 - 4, 0);
        static const int I2_s   = std::max(I2 - 4, 0);

        static const int ind0   = I1_s + (I2_s << 2) + (I3 << 4) + (I4 << 6);
        static const int ind    = std::min(std::max(ind0, 0), 255);

        return _mm_shuffle_ps(y.data, x.data, ind);
    }

    static const bool is_cont_x = ((I1 == 0 || I1 >= 4) && (I2 == 1 || I2 >= 4) && (I3 == 2 || I3 >= 4) &&
                                    (I4 == 3 || I4 >= 4));
    static const bool is_cont_y = ((I1 <= 3 || I1 == 4) && (I2 <= 3 || I2 == 5) && (I3 <= 3 || I3 == 6) &&
                                    (I4 <= 3 || I4 == 7));

    simd x_perm = is_cont_x ? x : x.select<std::min(I1, 3), std::min(I2, 3), std::min(I3, 3), std::min(I4,3)>();
    simd y_perm = is_cont_y ? y : y.select<std::max(I1-4, 0), std::max(I2-4, 0), std::max(I3-4, 0), 
                                           std::max(I4-4,0)>();


    // use if_then_else

    using simd_int  = simd<int32_t, 128, sse_tag>;

    simd_int cond   = details::vector_4_int<(I1 <= 3) ? -1 : 0, (I2 <= 3) ? -1 : 0,
                                            (I3 <= 3) ? -1 : 0, (I4 <= 3) ? -1 : 0>();        
         
    return if_then_else(cond.reinterpret_as_float(), x_perm, y_perm);
}

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

#pragma warning(pop)