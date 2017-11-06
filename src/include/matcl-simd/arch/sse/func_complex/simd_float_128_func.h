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
struct simd_compl_reverse<float, 128, sse_tag>
{
    using simd_type         = simd_compl<float, 128, sse_tag>;

    force_inline
    static simd_type eval(const simd_type& x)
    {
        return _mm_shuffle_ps(x.data.data, x.data.data, _MM_SHUFFLE(1,0,3,2));
    };
};

template<>
struct simd_compl_conj<float, 128, sse_tag>
{
    using simd_type         = simd_compl<float, 128, sse_tag>;

    force_inline
    static simd_type eval(const simd_type& x)
    {
        const __m128 mask = _mm_setr_ps(0.0f, -0.0f, 0.0f, -0.0f);
        __m128 res        = _mm_xor_ps(x.data.data, mask);
        return res;
    };
};

template<>
struct simd_compl_mult<float, 128, sse_tag>
{
    using simd_type         = simd_compl<float, 128, sse_tag>;
    using simd_real_type    = simd<float, 128, sse_tag>;

    // (x.re * y.re - x.im * y.im,  x.re * y.im + x.re * y.im)
    force_inline
    static simd_type eval(const simd_type& x, const simd_type& y)
    {        
        __m128 y_flip = _mm_shuffle_ps(y.data.data, y.data.data, 0xB1); // swap y.re and y.im
        __m128 x_im   = _mm_shuffle_ps(x.data.data, x.data.data, 0xF5); // imag of x in both
        __m128 x_re   = _mm_shuffle_ps(x.data.data, x.data.data, 0xA0); // real of x in both
        __m128 x_imy  = _mm_mul_ps(x_im, y_flip);                       // (x.im*y.im, x.im*y.re)

        #if MATCL_ARCHITECTURE_HAS_FMA
            __m128 res  = _mm_fmaddsub_ps(x_re, y.data.data, x_imy);    // a_re * y -/+ x_imy
        #else
            __m128 x_rey = _mm_mul_ps(x_re, y.data.data);               // a_re * y
            simd_real_type xv_rey(x_rey);
            simd_real_type xv_imy(x_imy);
            __m128 res  = sub_add(xv_rey, xv_imy).data;                 // a_re * y -/+ x_imy
        #endif

        // check for NaN
        __m128 nt       = _mm_cmp_ps(res, res, _CMP_NEQ_UQ);
        int have_nan    = _mm_movemask_ps(nt);

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

        for (int i = 0; i < vec_size; ++i)
        {
            float r_re      = real(xy.get(i));
            float r_im      = imag(xy.get(i));
            value_type res2 = mult_impl::eval(x.get(i), y.get(i), r_re, r_im);
            res.set(i, res2);
        };

        return res;
    };

    force_inline
    static simd_type eval_rc(const simd_real_type& x, const simd_type& y)
    {
        __m128 res  = _mm_mul_ps(x.data, y.data.data);
        return res;
    }

    force_inline
    static simd_type eval_cr(const simd_type& x, const simd_real_type& y)
    {
        __m128 res  = _mm_mul_ps(x.data.data, y.data);
        return res;
    }
};

template<>
struct simd_compl_div<float, 128, sse_tag>
{
    using simd_type = simd_compl<float, 128, sse_tag>;
    using simd_real = simd<float, 128, sse_tag>;

    // (a.re * b.re + a.im * b.im, b.re * a.im - a.re * b.im) / (b.re * b.re + b.im * b.im)
    force_inline
    static simd_type eval(const simd_type& x, const simd_type& y)
    {
        __m128 y_flip   = _mm_shuffle_ps(y.data.data, y.data.data, 0xB1);   // swap y.re and y.im
        __m128 x_im     = _mm_shuffle_ps(x.data.data, x.data.data, 0xF5);   // imag of x in both
        __m128 x_re     = _mm_shuffle_ps(x.data.data, x.data.data, 0xA0);   // real of x in both
        simd_real x_rey = _mm_mul_ps(x_re, y.data.data);                    // (x.re*b.re, x.re*b.im)  
        __m128 yy       = _mm_mul_ps(y.data.data, y.data.data);             // (y.re*y.re, y.im*y.im)

        __m128 yy2    = _mm_shuffle_ps(yy, yy, 0xB1);                       // swap yy.re and yy.im
        __m128 yy3    = _mm_add_ps(yy,yy2);                                 // add pairwise

        bool overflow1 = any_inf(simd_real(yy3));

        #if MATCL_ARCHITECTURE_HAS_FMA
            __m128 n    = _mm_fmsubadd_ps(x_im, y_flip, x_rey.data);        // (x_im * y_im, x_im * y_re) +/- x_rey
        #else
            __m128 x_imy    = _mm_mul_ps(x_im, y_flip);                     // (x_im * y_im, x_im * y_re)
            simd_real xv_imy(x_imy);
            simd_real xv_rey(x_rey);
            __m128 n        = sub_add(xv_imy, -xv_rey).data;                // x_re * y +/- x_imy
        #endif
        
        bool overflow2  = any_inf(simd_real(n));
        __m128 res      = _mm_div_ps(n, yy3);
        bool nan        = any_nan(simd_real(res));

        if (nan || overflow1 || overflow2)
            return recover_nan(x, y, simd_type(res));
        else
            return res;
    };

    // (a * b.re, - a * b.im) / (b.re * b.re + b.im * b.im)
    force_inline
    static simd_type eval_rc(const simd_real& x, const simd_type& y)
    {
        const __m128 mask = _mm_setr_ps(0.0f, -0.0f, 0.0f, -0.0f);

        __m128 y_flip   = _mm_shuffle_ps(y.data.data, y.data.data, 0xB1);   // swap y.re and y.im
        __m128 x_re     = _mm_shuffle_ps(x.data, x.data, 0xA0);             // real of x in both
        __m128 x_rey    = _mm_mul_ps(x_re, y.data.data);                    // (x.re*b.re, x.re*b.im)  
        __m128 yy       = _mm_mul_ps(y.data.data, y.data.data);             // (y.re*y.re, y.im*y.im)        
        
        __m128 yy2      = _mm_shuffle_ps(yy, yy, 0xB1);       // swap yy.re and yy.im
        __m128 yy3      = _mm_add_ps(yy,yy2);                 // add pairwise
        bool overflow1  = any_inf(simd_real(yy3));

        __m128 n        = _mm_xor_ps(x_rey, mask);                          // +/- x_rey
        bool overflow2  = any_inf(simd_real(n));

        __m128 res      = _mm_div_ps(n, yy3);
        bool nan        = any_nan(simd_real(res));

        if (nan || overflow1 || overflow2)
            return recover_nan_rc(x, y, simd_type(res));
        else
            return res;
    };

    force_inline
    static simd_type eval_cr(const simd_type& x, const simd_real& y)
    {
        __m128 res  = _mm_div_ps(x.data.data, y.data);
        return res;
    }

    static simd_type recover_nan(const simd_type& x, const simd_type& y, const simd_type& xy)
    {
        using value_type            = typename simd_type::value_type;
        using div_impl              = details::recover_nan_div<float>;
        static const int vec_size   = simd_type::vector_size;

        simd_type res;

        for (int i = 0; i < vec_size; ++i)
        {
            value_type res2 = div_impl::eval(x.get(i), y.get(i));
            res.set(i, res2);
        };

        return res;
    };

    static simd_type recover_nan_rc(const simd_real& x, const simd_type& y, const simd_type& xy)
    {
        using value_type            = typename simd_type::value_type;
        using div_impl              = details::recover_nan_div_rc<float>;;
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
struct simd_compl_plus<float, 128, sse_tag>
{
    using simd_type         = simd_compl<float, 128, sse_tag>;

    force_inline
    static simd_type eval(const simd_type& x, const simd_type& y)
    {
        return _mm_add_ps( x.data.data, y.data.data );
    };
};

template<>
struct simd_compl_minus<float, 128, sse_tag>
{
    using simd_type         = simd_compl<float, 128, sse_tag>;

    force_inline
    static simd_type eval(const simd_type& x, const simd_type& y)
    {
        return _mm_sub_ps( x.data.data, y.data.data );
    };
};

template<>
struct simd_compl_uminus<float, 128, sse_tag>
{
    using simd_type         = simd_compl<float, 128, sse_tag>;

    force_inline
    static simd_type eval(const simd_type& x)
    {
        return _mm_xor_ps( x.data.data, _mm_set1_ps(-0.0f) );
    };
};

template<>
struct simd_compl_sum_all<float, 128, sse_tag>
{
    using simd_type         = simd_compl<float, 128, sse_tag>;

    force_inline
    static simd_single_complex eval(const simd_type& x)
    {
        simd_type res = x + reverse(x);
        return res.get<0>();
    };
};

}}
