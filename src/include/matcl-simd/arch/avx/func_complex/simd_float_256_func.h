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

#include "matcl-simd/arch/simd_impl.h"
#include "matcl-simd/complex/simd_complex_impl.h"
#include "matcl-simd/func/simd_func_complex_def.h"
#include "matcl-simd/arch/avx/helpers.h"

namespace matcl { namespace simd
{

template<>
struct simd_compl_reverse<float, 256, avx_tag>
{
    using simd_type         = simd_compl<float, 256, avx_tag>;

    force_inline
    static simd_type eval(const simd_type& x)
    {
        #if MATCL_ARCHITECTURE_HAS_AVX2 && !MATCL_TEST_MISSING
            static const __m256i control = details::vector_8_int<1,0,3,2,5,4,7,6>();
            return _mm256_permutevar8x32_ps(x.data.data, control);
        #else
            return simd_type(reverse(x.extract_high()), reverse(x.extract_low()));
        #endif
    };
};

template<>
struct simd_compl_mult<float, 256, avx_tag>
{
    using simd_type         = simd_compl<float, 256, avx_tag>;
    using simd_real_type    = simd<float, 256, avx_tag>;

    // (x.re * y.re - x.im * y.im,  x.re * y.im + x.re * y.im)
    force_inline
    static simd_type eval(const simd_type& x, const simd_type& y)
    {
        __m256 y_flip = _mm256_shuffle_ps(y.data.data, y.data.data, 0xB1);  // swap y.re and y.im
        __m256 x_im   = _mm256_shuffle_ps(x.data.data, x.data.data, 0xF5);  // imag of x in both
        __m256 x_re   = _mm256_shuffle_ps(x.data.data, x.data.data, 0xA0);  // real of x in both
        __m256 x_imy  = _mm256_mul_ps(x_im, y_flip);                        // (x.im*y.im, x.im*y.re)

        #if MATCL_ARCHITECTURE_HAS_FMA && !MATCL_TEST_MISSING
            return  _mm256_fmaddsub_ps(x_re, y.data.data, x_imy);           // a_re * y -/+ x_imy
        #else
            __m256 x_rey = _mm256_mul_ps(x_re, y.data.data);                // a_re * y
            simd_real_type xv_rey(x_rey);
            simd_real_type xv_imy(x_imy);
            return sub_add(xv_rey, xv_imy);                                 // a_re * y +/- x_imy
        #endif
    };

    force_inline
    static simd_type eval(const simd_real_type& x, const simd_type& y)
    {
        __m256 res  = _mm256_mul_ps(x.data, y.data.data);
        return res;
    }

    force_inline
    static simd_type eval(const simd_type& x, const simd_real_type& y)
    {
        __m256 res  = _mm256_mul_ps(x.data.data, y.data);
        return res;
    }
};

template<>
struct simd_compl_div<float, 256, avx_tag>
{
    using simd_type = simd_compl<float, 256, avx_tag>;
    using simd_real = simd<float, 256, avx_tag>;

    // (a.re * b.re + a.im * b.im, b.re * a.im - a.re * b.im) / (b.re * b.re + b.im * b.im)
    force_inline
    static simd_type eval(const simd_type& x, const simd_type& y)
    {
        __m256 y_flip = _mm256_shuffle_ps(y.data.data, y.data.data, 0xB1);  // swap y.re and y.im
        __m256 x_im   = _mm256_shuffle_ps(x.data.data, x.data.data, 0xF5);  // imag of x in both
        __m256 x_re   = _mm256_shuffle_ps(x.data.data, x.data.data, 0xA0);  // real of x in both
        __m256 x_rey  = _mm256_mul_ps(x_re, y.data.data);                   // (x.re*b.re, x.re*b.im)  
        __m256 yy     = _mm256_mul_ps(y.data.data, y.data.data);            // (y.re*y.re, y.im*y.im)

        #if MATCL_ARCHITECTURE_HAS_FMA && !MATCL_TEST_MISSING
            __m256 n      = _mm256_fmsubadd_ps(x_im, y_flip, x_rey);        // x_re * y +/- x_imy
        #else
            __m256 x_imy    = _mm256_mul_ps(x_im, y_flip);                  // (x_im * y_im, x_im * y_re)
            simd_real xv_imy(x_imy);
            simd_real xv_rey(x_rey);
            __m256 n        = sub_add(xv_imy, -xv_rey).data;                // x_re * y +/- x_imy
        #endif
        
        __m256 yy2    = _mm256_shuffle_ps(yy,yy,0xB1);                      // Swap yy.re and yy.im
        __m256 yy3    = _mm256_add_ps(yy,yy2);                              // (y.re*y.re + y.im*y.im) dublicated

        return _mm256_div_ps(n, yy3);
    };
};

template<>
struct simd_compl_plus<float, 256, avx_tag>
{
    using simd_type         = simd_compl<float, 256, avx_tag>;

    force_inline
    static simd_type eval(const simd_type& x, const simd_type& y)
    {
        return _mm256_add_ps( x.data.data, y.data.data );
    };
};

template<>
struct simd_compl_minus<float, 256, avx_tag>
{
    using simd_type         = simd_compl<float, 256, avx_tag>;

    force_inline
    static simd_type eval(const simd_type& x, const simd_type& y)
    {
        return _mm256_sub_ps( x.data.data, y.data.data );
    };
};

template<>
struct simd_compl_uminus<float, 256, avx_tag>
{
    using simd_type         = simd_compl<float, 256, avx_tag>;

    force_inline
    static simd_type eval(const simd_type& x)
    {
        float Z = -0.0f;
        return _mm256_xor_ps(x.data.data, _mm256_broadcast_ss(&Z));
    };
};

template<>
struct simd_compl_fma<float, 256, avx_tag>
{
    using simd_type         = simd_compl<float, 256, avx_tag>;
    using simd_real_type    = simd<float, 256, avx_tag>;

    force_inline
    static simd_type eval(const simd_type& x, const simd_type& y, const simd_type& z)
    {
        return x * y + z;
    };

    force_inline
    static simd_type eval(const simd_real_type& x, const simd_type& y, const simd_type& z)
    {
        return x * y + z;
    };
};

template<>
struct simd_compl_fms<float, 256, avx_tag>
{
    using simd_type         = simd_compl<float, 256, avx_tag>;
    using simd_real_type    = simd<float, 256, avx_tag>;

    force_inline
    static simd_type eval(const simd_type& x, const simd_type& y, const simd_type& z)
    {
        return x * y - z;
    };

    force_inline
    static simd_type eval(const simd_real_type& x, const simd_type& y, const simd_type& z)
    {
        return x * y - z;
    };
};

template<>
struct simd_compl_sum_all<float, 256, avx_tag>
{
    using simd_type         = simd_compl<float, 256, avx_tag>;

    force_inline
    static simd_single_complex eval(const simd_type& x)
    {
        using simd_compl = simd_compl<float, 256, avx_tag>;

        __m128 lo   = _mm256_extractf128_ps(x.data.data, 0);
        __m128 hi   = _mm256_extractf128_ps(x.data.data, 1);

        __m128 sum  = _mm_add_ps(lo, hi);
        __m128 rev  = _mm_shuffle_ps(sum, sum, _MM_SHUFFLE(1,0,3,2));
        __m128 sum2 = _mm_add_ps(sum, rev);

        return simd_single_complex(sum2.m128_f32[0], sum2.m128_f32[1]);
    };
};

}}
