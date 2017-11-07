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
#include "matcl-simd/complex/recover_nan.h"

namespace matcl { namespace simd
{

template<>
struct simd_compl_reverse<double, 256, avx_tag>
{
    using simd_type         = simd_compl<double, 256, avx_tag>;

    force_inline
    static simd_type eval(const simd_type& x)
    {
        #if MATCL_ARCHITECTURE_HAS_AVX2
            static const int control = 2 + (3 << 2) + (0 << 4) + (1 << 6);
            return _mm256_permute4x64_pd(x.data.data, control);
        #else
            return simd_type(x.extract_high(), x.extract_low());
        #endif
    };
};

template<>
struct simd_compl_conj<double, 256, avx_tag>
{
    using simd_type         = simd_compl<double, 256, avx_tag>;

    force_inline
    static simd_type eval(const simd_type& x)
    {
        const __m256d mask = _mm256_setr_pd(0.0, -0.0, 0.0, -0.0);
        __m256d res        = _mm256_xor_pd(x.data.data, mask);
        return res;
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

        #if MATCL_ARCHITECTURE_HAS_FMA
            __m256d res = _mm256_fmaddsub_pd(x_re, y.data.data, x_imy);     // a_re * y -/+ x_imy
        #else
            __m256d x_rey = _mm256_mul_pd(x_re, y.data.data);               // a_re * y
            simd_real_type xv_rey(x_rey);
            simd_real_type xv_imy(x_imy);
            __m256d res = sub_add(xv_rey, xv_imy).data;                     // a_re * y -/+ x_imy
        #endif

        // check for NaN
        __m256d nt      = _mm256_cmp_pd(res, res, _CMP_NEQ_UQ);
        int have_nan    = _mm256_movemask_pd(nt);

        if (have_nan != 0)
            return recover_nan(x, y, simd_type(res));
        else
            return res;
    };

    static simd_type recover_nan(const simd_type& x, const simd_type& y, const simd_type& xy)
    {
        using value_type            = typename simd_type::value_type;
        using mult_impl             = details::recover_nan_mul<double>;
        static const int vec_size   = simd_type::vector_size;

        simd_type res;

        for (int i = 0; i < vec_size; ++i)
        {
            double r_re     = real(xy.get(i));
            double r_im     = imag(xy.get(i));
            value_type res2 = mult_impl::eval(x.get(i), y.get(i), r_re, r_im);
            res.set(i, res2);
        };

        return res;
    };

    force_inline
    static simd_type eval_rc(const simd_real_type& x, const simd_type& y)
    {
        __m256d res  = _mm256_mul_pd(x.data, y.data.data);
        return res;
    }

    force_inline
    static simd_type eval_cr(const simd_type& x, const simd_real_type& y)
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
        bool o1         = check_overflow(x);
        bool o2         = check_overflow(y);

        __m256d y_flip = _mm256_shuffle_pd(y.data.data, y.data.data, 5);    // swap y.re and y.im
        __m256d x_im   = _mm256_shuffle_pd(x.data.data, x.data.data, 0xF);  // imag of x in both
        __m256d x_re   = _mm256_shuffle_pd(x.data.data, x.data.data, 0);    // real of x in both

        __m256d x_rey  = _mm256_mul_pd(x_re, y.data.data);                  // (x.re*b.re, x.re*b.im)  
        __m256d yy     = _mm256_mul_pd(y.data.data, y.data.data);           // (y.re*y.re, y.im*y.im)
        __m256d yy2    = _mm256_hadd_pd(yy,yy);                             // (y.re*y.re + y.im*y.im)

        #if MATCL_ARCHITECTURE_HAS_FMA
            __m256d n      = _mm256_fmsubadd_pd(x_im, y_flip, x_rey);       // x_re * y +/- x_imy
        #else
            __m256d x_imy   = _mm256_mul_pd(x_im, y_flip);                  // (x_im * y_im, x_im * y_re)
            simd_real xv_imy(x_imy);
            simd_real xv_rey(x_rey);
            __m256d n       = sub_add(xv_imy, -xv_rey).data;                // x_re * y +/- x_imy
        #endif
        
        __m256d res     = _mm256_div_pd(n, yy2);

        if (o1 || o2)
            return recover_nan(x, y);
        else
            return res;
    };

    // (a * b.re, - a * b.im) / (b.re * b.re + b.im * b.im)
    force_inline
    static simd_type eval_rc(const simd_real& x, const simd_type& y)
    {
        __m256d x_re   = _mm256_shuffle_pd(x.data, x.data, 0);              // real of x in both

        bool o1         = check_overflow(simd_real(x_re));
        bool o2         = check_overflow(y);

        const __m256d mask = _mm256_setr_pd(0.0, -0.0, 0.0, -0.0);

        __m256d y_flip = _mm256_shuffle_pd(y.data.data, y.data.data, 5);    // swap y.re and y.im        
        __m256d x_rey  = _mm256_mul_pd(x_re, y.data.data);                  // (x.re*b.re, x.re*b.im)  
        __m256d yy     = _mm256_mul_pd(y.data.data, y.data.data);           // (y.re*y.re, y.im*y.im)
        __m256d yy2    = _mm256_hadd_pd(yy,yy);                             // (y.re*y.re + y.im*y.im) 

        __m256d n      = _mm256_xor_pd(x_rey, mask);                        // +/- x_rey
        __m256d res    = _mm256_div_pd(n, yy2);

        if (o1 || o2)
            return recover_nan_rc(x, y);
        else
            return res;
    };

    force_inline
    static simd_type eval_cr(const simd_type& x, const simd_real& y)
    {
        __m256d res  = _mm256_div_pd(x.data.data, y.data);
        return res;
    }

    force_inline
    static bool check_overflow(const simd_type& x)
    {
        //MIN = 2.2250738585072014e-308
        //MAX = 1.7976931348623158e+308

        simd_real xa        = abs(x.data);
        const __m256d max   = _mm256_set1_pd(9.48e+153);    // max < sqrt(MAX/2)
        const __m256d min   = _mm256_set1_pd(1.42e-146);    // min > sqrt(MIN/eps*2)

        bool res1           = !all(gt(xa, simd_real(min)));
        bool res2           = !all(lt(xa, simd_real(max)));

        return res1 || res2;
    };

    static simd_type recover_nan(const simd_type& x, const simd_type& y)
    {
        using value_type            = typename simd_type::value_type;
        using div_impl              = details::recover_nan_div<double>;
        static const int vec_size   = simd_type::vector_size;

        simd_type res;

        for (int i = 0; i < vec_size; ++i)
        {
            value_type res2 = div_impl::eval(x.get(i), y.get(i));
            res.set(i, res2);
        };

        return res;
    };

    static simd_type recover_nan_rc(const simd_real& x, const simd_type& y)
    {
        using value_type            = typename simd_type::value_type;
        using div_impl              = details::recover_nan_div_rc<double>;
        static const int vec_size   = simd_type::vector_size;

        simd_type res;

        for (int i = 0; i < vec_size; ++i)
        {
            value_type res2 = div_impl::eval(x.get(2*i), y.get(i));
            res.set(i, res2);
        };

        return res;
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
