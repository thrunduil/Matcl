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

namespace matcl { namespace simd
{

template<>
struct simd_compl_reverse<double, 256, avx_tag>
{
    using simd_type         = simd_compl<double, 256, avx_tag>;

    force_inline
    static simd_type eval(const simd_type& x)
    {
        #if MATCL_ARCHITECTURE_HAS_AVX2 && !MATCL_TEST_MISSING
            static const int control = 2 + (3 << 2) + (0 << 4) + (1 << 6);
            return _mm256_permute4x64_pd(x.data.data, control);
        #else
            return simd_type(x.extract_high(), x.extract_low());
        #endif
    };
};

template<>
struct simd_compl_mult<double, 256, avx_tag>
{
    using simd_type         = simd_compl<double, 256, avx_tag>;
    using simd_real_type    = simd<double, 256, avx_tag>;

    // (x.re * y.re - x.im * y.im,  x.re * y.im + x.re * y.im)
    force_inline
    static simd_type eval(const simd_type& x, const simd_type& y)
    {
        __m256d y_flip = _mm256_shuffle_pd(y.data.data, y.data.data, 5);    // swap y.re and y.im
        __m256d x_im   = _mm256_shuffle_pd(x.data.data, x.data.data, 0xF);  // imag of x in both
        __m256d x_re   = _mm256_shuffle_pd(x.data.data, x.data.data, 0);    // real of x in both
        __m256d x_imy  = _mm256_mul_pd(x_im, y_flip);                       // (x.im*y.im, x.im*y.re)

        #if MATCL_ARCHITECTURE_HAS_FMA && !MATCL_TEST_MISSING
            return  _mm256_fmaddsub_pd(x_re, y.data.data, x_imy);           // a_re * y -/+ x_imy
        #else
            __m256d x_rey = _mm256_mul_pd(x_re, y.data.data);               // a_re * y
            simd_real_type xv_rey(x_rey);
            simd_real_type xv_imy(x_imy);
            return sub_add(xv_rey, xv_imy);                                 // a_re * y -/+ x_imy
        #endif
    };

    force_inline
    static simd_type eval(const simd_real_type& x, const simd_type& y)
    {
        __m256d res  = _mm256_mul_pd(x.data, y.data.data);
        return res;
    }

    force_inline
    static simd_type eval(const simd_type& x, const simd_real_type& y)
    {
        __m256d res  = _mm256_mul_pd(x.data.data, y.data);
        return res;
    }
};

template<>
struct simd_compl_div<double, 256, avx_tag>
{
    using simd_type = simd_compl<double, 256, avx_tag>;
    using simd_real = simd<double, 256, avx_tag>;

    // (a.re * b.re + a.im * b.im, b.re * a.im - a.re * b.im) / (b.re * b.re + b.im * b.im)
    force_inline
    static simd_type eval(const simd_type& x, const simd_type& y)
    {
        __m256d y_flip = _mm256_shuffle_pd(y.data.data, y.data.data, 5);    // swap y.re and y.im
        __m256d x_im   = _mm256_shuffle_pd(x.data.data, x.data.data, 0xF);  // imag of x in both
        __m256d x_re   = _mm256_shuffle_pd(x.data.data, x.data.data, 0);    // real of x in both
        __m256d x_rey  = _mm256_mul_pd(x_re, y.data.data);                  // (x.re*b.re, x.re*b.im)  
        simd_real yy   = _mm256_mul_pd(y.data.data, y.data.data);           // (y.re*y.re, y.im*y.im)

        #if MATCL_ARCHITECTURE_HAS_FMA && !MATCL_TEST_MISSING
            __m256d n      = _mm256_fmsubadd_pd(x_im, y_flip, x_rey);       // x_re * y +/- x_imy
        #else
            __m256d x_imy   = _mm256_mul_pd(x_im, y_flip);                  // (x_im * y_im, x_im * y_re)
            simd_real xv_imy(x_imy);
            simd_real xv_rey(x_rey);
            __m256d n       = sub_add(xv_imy, -xv_rey).data;                // x_re * y +/- x_imy
        #endif
        
        __m256d yy2    = horizontal_add(yy,yy).data;                        // (y.re*y.re + y.im*y.im) 

        return _mm256_div_pd(n, yy2);
    };
};

template<>
struct simd_compl_plus<double, 256, avx_tag>
{
    using simd_type         = simd_compl<double, 256, avx_tag>;

    force_inline
    static simd_type eval(const simd_type& x, const simd_type& y)
    {
        return _mm256_add_pd( x.data.data, y.data.data );
    };
};

template<>
struct simd_compl_minus<double, 256, avx_tag>
{
    using simd_type         = simd_compl<double, 256, avx_tag>;

    force_inline
    static simd_type eval(const simd_type& x, const simd_type& y)
    {
        return _mm256_sub_pd( x.data.data, y.data.data );
    };
};

template<>
struct simd_compl_uminus<double, 256, avx_tag>
{
    using simd_type         = simd_compl<double, 256, avx_tag>;

    force_inline
    static simd_type eval(const simd_type& x)
    {
        double Z = -0.0;
        return _mm256_xor_pd(x.data.data, _mm256_broadcast_sd(&Z));
    };
};

template<>
struct simd_compl_fma<double, 256, avx_tag>
{
    using simd_type         = simd_compl<double, 256, avx_tag>;
    using simd_real_type    = simd<double, 256, avx_tag>;

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
struct simd_compl_fms<double, 256, avx_tag>
{
    using simd_type         = simd_compl<double, 256, avx_tag>;
    using simd_real_type    = simd<double, 256, avx_tag>;

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
struct simd_compl_sum_all<double, 256, avx_tag>
{
    using simd_type = simd_compl<double, 256, avx_tag>;

    force_inline
    static simd_double_complex eval(const simd_type& x)
    {
        using simd_half     = simd_type::simd_half;

        simd_half x1    = x.extract_low();
        simd_half x2    = x.extract_high();
        simd_half res   = x1 + x2;

        return res.get<0>();
    };
};

}}
