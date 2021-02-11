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

#include "matcl-simd/arch/sse/simd_int32_128.h"
#include "matcl-simd/details/arch/sse/func/missing_intrinsics.h"

#pragma warning(push)
#pragma warning(disable : 4127) // conditional expression is constant

namespace matcl { namespace simd { namespace details
{

template<class Ret>
struct cast_int32_int64_impl
{
    using simd_type = simd<int32_t, 128, matcl::simd::sse_tag>;

    force_inline
    static Ret eval(const simd_type& s)
    {
        return Ret(s.convert_low_to_int64(), s.convert_high_to_int64());
    };
};

#if MATCL_ARCHITECTURE_HAS_AVX
    template<>
    struct cast_int32_int64_impl<simd<int64_t, 256, matcl::simd::avx_tag>>
    {
        using Ret       = simd<int64_t,256, matcl::simd::avx_tag>;
        using simd_type = simd<int32_t, 128, matcl::simd::sse_tag>;

        force_inline
        static Ret eval(const simd_type& s)
        {
            #if MATCL_ARCHITECTURE_HAS_AVX2
                return _mm256_cvtepi32_epi64(s.data);
            #else      
                return Ret(s.convert_low_to_int64(), s.convert_high_to_int64());
            #endif
        };
    };
#endif

template<class Ret>
struct cast_int32_double_impl
{
    using simd_type = simd<int32_t, 128, matcl::simd::sse_tag>;

    force_inline
    static Ret eval(const simd_type& s)
    {
        return Ret(s.convert_low_to_double(), s.convert_high_to_double());
    };
};

#if MATCL_ARCHITECTURE_HAS_AVX
    template<>
    struct cast_int32_double_impl<simd<double, 256, matcl::simd::avx_tag>>
    {
        using Ret       = simd<double,256, matcl::simd::avx_tag>;
        using simd_type = simd<int32_t, 128, matcl::simd::sse_tag>;

        force_inline
        static Ret eval(const simd_type& s)
        {
            #if MATCL_ARCHITECTURE_HAS_AVX2
                return _mm256_cvtepi32_pd(s.data);
            #else
                return Ret(s.convert_low_to_double(), s.convert_high_to_double());
            #endif
        };
    };
#endif

}}};

namespace matcl { namespace simd
{

//-------------------------------------------------------------------
//                          SSE int32_t
//-------------------------------------------------------------------

force_inline
simd<int32_t, 128, sse_tag>::simd(int32_t val)
    : data(_mm_set1_epi32(val)) 
{}

force_inline
simd<int32_t, 128, sse_tag>::simd(int32_t v1, int32_t v2, int32_t v3, int32_t v4) 
    : data(_mm_setr_epi32(v1, v2, v3, v4)) 
{}

force_inline
simd<int32_t, 128, sse_tag>::simd(const impl_type& v)
    : data(v) 
{};

force_inline
simd<int32_t, 128, sse_tag>::simd(const simd& lo, const simd& hi)
    : data(_mm_unpacklo_epi64(lo.data, hi.data))
{};

force_inline
simd<int32_t, 128, sse_tag>::simd(const simd<int32_t, 128, nosimd_tag>& s)
    : data(_mm_load_si128((const __m128i*)s.data))
{};

force_inline
simd<int32_t, 128, sse_tag>::simd(const simd<int32_t, 128, scalar_sse_tag>& s)
    : data(_mm_shuffle_epi32(s.data, 0))
{};

force_inline
simd<int32_t, 128, sse_tag>::simd(const simd<int32_t, 128, scalar_nosimd_tag>& s)
    : simd(s.first())
{};

force_inline simd<int32_t, 128, sse_tag> 
simd<int32_t, 128, sse_tag>::broadcast(const int32_t* arr)
{ 
    return _mm_set1_epi32(arr[0]);
};

force_inline simd<int32_t, 128, sse_tag> 
simd<int32_t, 128, sse_tag>::broadcast(const int32_t& arr)
{ 
    return _mm_set1_epi32(arr);
};

force_inline
int32_t simd<int32_t, 128, sse_tag>::get(int pos) const
{ 
    return get_raw_ptr()[pos] ; 
};

force_inline
int32_t simd<int32_t, 128, sse_tag>::first() const
{ 
    return _mm_cvtsi128_si32(data); 
};

force_inline
void simd<int32_t, 128, sse_tag>::set(int pos, int32_t val)
{ 
    get_raw_ptr()[pos] = val; 
};

force_inline
const int32_t* simd<int32_t, 128, sse_tag>::get_raw_ptr() const
{ 
    return reinterpret_cast<const int32_t*>(&data); 
};

force_inline
int32_t* simd<int32_t, 128, sse_tag>::get_raw_ptr()
{ 
    return reinterpret_cast<int32_t*>(&data); 
};

force_inline
simd<int32_t, 128, sse_tag> simd<int32_t, 128, sse_tag>::zero()
{
    impl_type data   = _mm_setzero_si128();
    return data;
}

force_inline
simd<int32_t, 128, sse_tag> simd<int32_t, 128, sse_tag>::one()
{
    return simd(1);
}

force_inline
simd<int32_t, 128, sse_tag> simd<int32_t, 128, sse_tag>::minus_one()
{
    return simd(-1);
}

force_inline simd<int32_t, 128, sse_tag> 
simd<int32_t, 128, sse_tag>::load(const int32_t* arr, std::true_type aligned)
{
    (void)aligned;
    return _mm_load_si128 ((const __m128i*)arr);
};

force_inline simd<int32_t, 128, sse_tag> 
simd<int32_t, 128, sse_tag>::load(const int32_t* arr, std::false_type not_aligned)
{
    (void)not_aligned;
    return _mm_loadu_si128 ((const __m128i*)arr);
};

force_inline simd<int32_t, 128, sse_tag>
simd<int32_t, 128, sse_tag>::gather(const int32_t* arr, const simd_int32& ind)
{
    #if MATCL_ARCHITECTURE_HAS_AVX2
        return _mm_i32gather_epi32(arr, ind.data, 4);
    #else
        simd res;
        int32_t* res_ptr        = res.get_raw_ptr();
        const int32_t* ind_ptr  = ind.get_raw_ptr();

        for (int i = 0; i < vector_size; ++i)
            res_ptr[i]  = arr[ind_ptr[i]];

        return res;
    #endif
}

force_inline simd<int32_t, 128, sse_tag>
simd<int32_t, 128, sse_tag>::gather(const int32_t* arr, const simd_int64& ind)
{
    #if MATCL_ARCHITECTURE_HAS_AVX2
        return simd(_mm_i64gather_epi32(arr, ind.data, 4));
    #else
        simd res;
        int32_t* res_ptr        = res.get_raw_ptr();
        const int64_t* ind_ptr  = ind.get_raw_ptr();

        for (int i = 0; i < vector_size/2; ++i)
            res_ptr[i]  = arr[ind_ptr[i]];

        return simd(res, simd::zero());
    #endif
}

force_inline simd<int32_t, 128, sse_tag> 
simd<int32_t, 128, sse_tag>::set_lower(int32_t v)
{
    return _mm_cvtsi32_si128(v);
};

force_inline void 
simd<int32_t, 128, sse_tag>::store(int32_t* arr, std::true_type aligned) const
{
    (void)aligned;
    _mm_store_si128((__m128i*)arr, data);
};


force_inline void 
simd<int32_t, 128, sse_tag>::store(int32_t* arr, std::false_type not_aligned) const
{
    (void)not_aligned;
    _mm_storeu_si128((__m128i*)arr, data);
};

force_inline simd<int64_t, 128, sse_tag>
simd<int32_t, 128, sse_tag>::convert_low_to_int64() const
{
    #if MATCL_ARCHITECTURE_HAS_SSE41
        return _mm_cvtepi32_epi64(data);
    #else
        // sign bit
        __m128i sign = _mm_srai_epi32(data, 31);

        // interleave with sign extensions
        return _mm_unpacklo_epi32(data, sign);
    #endif
};

force_inline
simd<int64_t, 128, sse_tag>
simd<int32_t, 128, sse_tag>::convert_high_to_int64() const
{    
    #if MATCL_ARCHITECTURE_HAS_SSE41
        __m128i hi = missing::mm_movehl_epi32(data, data);
        return _mm_cvtepi32_epi64(hi);
    #else
        // sign bit
        __m128i sign = _mm_srai_epi32(data, 31);

        // interleave with sign extensions
        return _mm_unpackhi_epi32(data, sign);
    #endif
};

force_inline
typename simd<int32_t, 128, sse_tag>::simd_int64_2
simd<int32_t, 128, sse_tag>::convert_to_int64() const
{
    return details::cast_int32_int64_impl<simd_int64_2>::eval(data);    
};

force_inline simd<float, 128, sse_tag>
simd<int32_t, 128, sse_tag>::convert_to_float() const
{
    return _mm_cvtepi32_ps(data);
};

force_inline simd<double, 128, sse_tag>
simd<int32_t, 128, sse_tag>::convert_low_to_double() const
{
    return _mm_cvtepi32_pd(data);
};

force_inline simd<double, 128, sse_tag>
simd<int32_t, 128, sse_tag>::convert_high_to_double() const
{
    simd hi = this->extract_high();
    return _mm_cvtepi32_pd(hi.data);
};

force_inline typename simd<int32_t, 128, sse_tag>::simd_double_2
simd<int32_t, 128, sse_tag>::convert_to_double() const
{
    return details::cast_int32_double_impl<simd_double_2>::eval(data);
};

force_inline simd<double, 128, sse_tag>
simd<int32_t, 128, sse_tag>::reinterpret_as_double() const
{
    return _mm_castsi128_pd(data);
}

force_inline simd<float, 128, sse_tag>
simd<int32_t, 128, sse_tag>::reinterpret_as_float() const
{
    return _mm_castsi128_ps(data);
}

force_inline simd<int64_t, 128, sse_tag>
simd<int32_t, 128, sse_tag>::reinterpret_as_int64() const
{
    return data;
}

force_inline simd<int32_t, 128, sse_tag> 
simd<int32_t, 128, sse_tag>::extract_low() const
{
    return *this;
}

force_inline simd<int32_t, 128, sse_tag> 
simd<int32_t, 128, sse_tag>::extract_high() const
{
    return missing::mm_movehl_epi32(data, data);
}

template<int I1, int I2, int I3, int I4>
force_inline simd<int32_t, 128, sse_tag>  
simd<int32_t, 128, sse_tag>::select() const
{
    static_assert(I1 >= 0 && I1 <= 3, "invalid index in select function");
    static_assert(I2 >= 0 && I2 <= 3, "invalid index in select function");
    static_assert(I3 >= 0 && I3 <= 3, "invalid index in select function");
    static_assert(I4 >= 0 && I4 <= 3, "invalid index in select function");

    if (I1 == 0 && I2 == 1 && I3 == 2 && I4 == 3)
        return data;

    static const int ind    = I1 + (I2 << 2) + (I3 << 4) + (I4 << 6);
    return _mm_shuffle_epi32(data, ind);
}

template<int I1, int I2, int I3, int I4>
force_inline simd<int32_t, 128, sse_tag>  
simd<int32_t, 128, sse_tag>::combine(const simd& x, const simd& y)
{
    using simd_float = simd<float, 128, sse_tag>;
    return simd_float::combine<I1,I2,I3,I4>(x.reinterpret_as_float(),
                    y.reinterpret_as_float()).reinterpret_as_int32();
}

template<int Step>
force_inline
void simd<int32_t, 128, sse_tag>::scatter(int32_t* arr) const
{
    const int32_t* this_data  = get_raw_ptr();

    //no scatter intrinsic
    arr[0*Step]    = this_data[0];
    arr[1*Step]    = this_data[1];
    arr[2*Step]    = this_data[2];
    arr[3*Step]    = this_data[3];
};

force_inline simd<int32_t, 128, sse_tag>& 
simd<int32_t, 128, sse_tag>::operator+=(const simd& x)
{
    *this = *this + x;
    return *this;
}

force_inline simd<int32_t, 128, sse_tag>& 
simd<int32_t, 128, sse_tag>::operator-=(const simd& x)
{
    *this = *this - x;
    return *this;
}

force_inline simd<int32_t, 128, sse_tag>& 
simd<int32_t, 128, sse_tag>::operator*=(const simd& x)
{
    *this = *this * x;
    return *this;
}

}}

#pragma warning(pop)