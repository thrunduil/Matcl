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

#include "matcl-simd/details/arch/simd_impl.h"
#include "matcl-simd/details/complex/simd_complex_impl.h"
#include "matcl-simd/details/func/simd_func_complex_def.h"
#include "matcl-simd/details/arch/avx/helpers.h"
#include "matcl-simd/details/complex/recover_nan.h"

namespace matcl { namespace simd { namespace details
{

template<>
struct simd_compl_reverse<float, 256, avx_tag>
{
    using simd_type         = simd_compl<float, 256, avx_tag>;

    force_inline
    static simd_type eval(const simd_type& x)
    {
        #if MATCL_ARCHITECTURE_HAS_AVX2
            static const __m256i control = details::vector_8_int<6,7,4,5,2,3,0,1>();
            return _mm256_permutevar8x32_ps(x.data.data, control);
        #else
            __m256 xp   = _mm256_permute2f128_ps(x.data.data, x.data.data, 1);
            xp          = _mm256_permute_ps(xp, _MM_SHUFFLE(1,0,3,2));
            return xp;
        #endif
    };
};

template<>
struct simd_compl_conj<float, 256, avx_tag>
{
    using simd_type         = simd_compl<float, 256, avx_tag>;

    force_inline
    static simd_type eval(const simd_type& x)
    {
        const __m256 mask = _mm256_setr_ps(0.0f, -0.0f, 0.0f, -0.0f, 0.0f, -0.0f, 0.0f, -0.0f);
        __m256 res        = _mm256_xor_ps(x.data.data, mask);
        return res;
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

        #if MATCL_ARCHITECTURE_HAS_FMA
            __m256 res  = _mm256_fmaddsub_ps(x_re, y.data.data, x_imy);     // a_re * y -/+ x_imy
        #else
            __m256 x_rey = _mm256_mul_ps(x_re, y.data.data);                // a_re * y
            simd_real_type xv_rey(x_rey);
            simd_real_type xv_imy(x_imy);
            __m256 res = sub_add(xv_rey, xv_imy).data;                      // a_re * y +/- x_imy
        #endif

        // check for NaN
        __m256 nt       = _mm256_cmp_ps(res, res, _CMP_NEQ_UQ);
        int have_nan    = _mm256_movemask_ps(nt);

        if (have_nan != 0)
            return recover_nan(x, y, simd_type(res));
        else
            return res;
    };

    static simd_type recover_nan(const simd_type& x, const simd_type& y, const simd_type& xy)
    {
        using value_type            = typename simd_type::value_type;
        using mult_impl             = details::recover_nan_mul<float>;
        static const int vec_size   = simd_type::vector_size;

        simd_type res;

        value_type* res_ptr         = res.get_raw_ptr();
        const value_type* xy_ptr    = xy.get_raw_ptr();
        const value_type* x_ptr     = x.get_raw_ptr();
        const value_type* y_ptr     = y.get_raw_ptr();

        for (int i = 0; i < vec_size; ++i)
        {
            float r_re      = real(xy_ptr[i]);
            float r_im      = imag(xy_ptr[i]);

            value_type res2 = mult_impl::eval(x_ptr[i], y_ptr[i], r_re, r_im);
            res_ptr[i]      = res2;
        };

        return res;
    };

    force_inline
    static simd_type eval_rc(const simd_real_type& x, const simd_type& y)
    {
        __m256 res  = _mm256_mul_ps(x.data, y.data.data);
        return res;
    }

    force_inline
    static simd_type eval_cr(const simd_type& x, const simd_real_type& y)
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
        bool o1         = check_overflow(x);
        bool o2         = check_overflow(y);

        __m256 y_flip = _mm256_shuffle_ps(y.data.data, y.data.data, 0xB1);  // swap y.re and y.im
        __m256 x_im   = _mm256_shuffle_ps(x.data.data, x.data.data, 0xF5);  // imag of x in both
        __m256 x_re   = _mm256_shuffle_ps(x.data.data, x.data.data, 0xA0);  // real of x in both
        __m256 x_rey  = _mm256_mul_ps(x_re, y.data.data);                   // (x.re*b.re, x.re*b.im)  
        __m256 yy     = _mm256_mul_ps(y.data.data, y.data.data);            // (y.re*y.re, y.im*y.im)

        __m256 yy2    = _mm256_shuffle_ps(yy,yy,0xB1);                      // Swap yy.re and yy.im
        __m256 yy3    = _mm256_add_ps(yy,yy2);                              // (y.re*y.re + y.im*y.im) dublicated

        #if MATCL_ARCHITECTURE_HAS_FMA
            __m256 n      = _mm256_fmsubadd_ps(x_im, y_flip, x_rey);        // x_re * y +/- x_imy
        #else
            __m256 x_imy    = _mm256_mul_ps(x_im, y_flip);                  // (x_im * y_im, x_im * y_re)
            simd_real xv_imy(x_imy);
            simd_real xv_rey(x_rey);
            __m256 n        = sub_add(xv_imy, -xv_rey).data;                // x_re * y +/- x_imy
        #endif       

        __m256 res      = _mm256_div_ps(n, yy3);

        if (o1 || o2)
            return recover_nan(x, y);
        else
            return res;
    };

    // (a * b.re, - a * b.im) / (b.re * b.re + b.im * b.im)
    force_inline
    static simd_type eval_rc(const simd_real& x, const simd_type& y)
    {
        __m256 x_re     = _mm256_shuffle_ps(x.data, x.data, 0xA0);            // real of x in both
        bool o1         = check_overflow(simd_real(x_re));
        bool o2         = check_overflow(y);

        const __m256 mask = _mm256_setr_ps(0.0f, -0.0f, 0.0f, -0.0f, 0.0f, -0.0f, 0.0f, -0.0f);

        __m256 y_flip   = _mm256_shuffle_ps(y.data.data, y.data.data, 0xB1);  // swap y.re and y.im        
        __m256 x_rey    = _mm256_mul_ps(x_re, y.data.data);                   // (x.re*b.re, x.re*b.im)  

        __m256 yy       = _mm256_mul_ps(y.data.data, y.data.data);            // (y.re*y.re, y.im*y.im)
        __m256 yy2      = _mm256_shuffle_ps(yy,yy,0xB1);                      // Swap yy.re and yy.im
        __m256 yy3      = _mm256_add_ps(yy,yy2);                              // (y.re*y.re + y.im*y.im) dublicated

        __m256 n        = _mm256_xor_ps(x_rey, mask);                       // +/- x_rey       
        __m256 res      = _mm256_div_ps(n, yy3);

        if (o1 || o2)
            return recover_nan_rc(x, y);
        else
            return res;
    };

    force_inline
    static simd_type eval_cr(const simd_type& x, const simd_real& y)
    {
        __m256 res  = _mm256_div_ps(x.data.data, y.data);
        return res;
    }

    force_inline
    static bool check_overflow(const simd_type& x)
    {
        //MIN = 1.175494351e-38F
        //MAX = 3.402823466e+38F

        simd_real xa        = abs(x.data);
        const __m256 max    = _mm256_set1_ps(1.30e+19f); // max < sqrt(MAX/2)
        const __m256 min    = _mm256_set1_ps(4.45e-16f); // min > sqrt(MIN/eps*2)

        bool res1           = !all(gt(xa, simd_real(min)));
        bool res2           = !all(lt(xa, simd_real(max)));

        return res1 || res2;
    };

    static simd_type recover_nan(const simd_type& x, const simd_type& y)
    {
        using value_type            = typename simd_type::value_type;
        using div_impl              = details::recover_nan_div<float>;
        static const int vec_size   = simd_type::vector_size;

        simd_type res;

        value_type* res_ptr         = res.get_raw_ptr();
        const value_type* x_ptr     = x.get_raw_ptr();
        const value_type* y_ptr     = y.get_raw_ptr();

        for (int i = 0; i < vec_size; ++i)
        {
            value_type res2 = div_impl::eval(x_ptr[i], y_ptr[i]);
            res_ptr[i]      = res2;
        };

        return res;
    };

    static simd_type recover_nan_rc(const simd_real& x, const simd_type& y)
    {
        using value_type            = typename simd_type::value_type;
        using div_impl              = details::recover_nan_div_rc<float>;
        static const int vec_size   = simd_type::vector_size;

        simd_type res;

        value_type* res_ptr         = res.get_raw_ptr();
        const float* x_ptr          = x.get_raw_ptr();
        const value_type* y_ptr     = y.get_raw_ptr();

        for (int i = 0; i < vec_size; ++i)
        {
            value_type res2 = div_impl::eval(x_ptr[2*i], y_ptr[i]);
            res_ptr[i]      = res2;
        };

        return res;
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
    using simd_real         = simd<float, 256, avx_tag>;

    force_inline
    static simd_type eval(const simd_type& x)
    {
        const simd_real mzero = simd_real::minus_zero();
        return _mm256_xor_ps(x.data.data, mzero.data);
    };
};

template<>
struct simd_compl_sum_all<float, 256, avx_tag>
{
    using simd_type         = simd_compl<float, 256, avx_tag>;
    using complex_type      = simd_single_complex;

    force_inline
    static complex_type eval(const simd_type& x)
    {
        using simd_compl = simd_compl<float, 256, avx_tag>;

        __m128 lo   = _mm256_extractf128_ps(x.data.data, 0);
        __m128 hi   = _mm256_extractf128_ps(x.data.data, 1);

        __m128 sum  = _mm_add_ps(lo, hi);
        __m128 rev  = _mm_shuffle_ps(sum, sum, _MM_SHUFFLE(1,0,3,2));
        __m128 sum2 = _mm_add_ps(sum, rev);

        const complex_type* tmp = reinterpret_cast<const complex_type*>(&sum2);

        return complex_type(tmp[0]);
    };
};

}}}
