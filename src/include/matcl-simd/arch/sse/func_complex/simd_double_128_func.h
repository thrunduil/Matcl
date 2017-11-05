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
struct simd_compl_reverse<double, 128, sse_tag>
{
    using simd_type         = simd_compl<double, 128, sse_tag>;

    force_inline
    static simd_type eval(const simd_type& x)
    {
        return _mm_shuffle_pd(x.data.data, x.data.data, _MM_SHUFFLE2(0,1));
    };
};

template<>
struct simd_compl_mult<double, 128, sse_tag>
{
    using simd_type         = simd_compl<double, 128, sse_tag>;
    using simd_real_type    = simd<double, 128, sse_tag>;

    // (x.re * y.re - x.im * y.im,  x.re * y.im + x.re * y.im)
    force_inline
    static simd_type eval(const simd_type& x, const simd_type& y)
    {        
        __m128d y_flip = _mm_shuffle_pd(y.data.data, y.data.data, 1);   // swap y.re and y.im
        __m128d x_im   = _mm_shuffle_pd(x.data.data, x.data.data, 3);   // imag of x in both
        __m128d x_re   = _mm_shuffle_pd(x.data.data, x.data.data, 0);   // real of x in both
        __m128d x_imy  = _mm_mul_pd(x_im, y_flip);                      // (x.im*y.im, x.im*y.re)

        #if MATCL_ARCHITECTURE_HAS_FMA
            return  _mm_fmaddsub_pd(x_re, y.data.data, x_imy);          // a_re * y -/+ x_imy
        #else
            __m128d x_rey = _mm_mul_pd(x_re, y.data.data);              // a_re * y
            simd_real_type xv_rey(x_rey);
            simd_real_type xv_imy(x_imy);
            return sub_add(xv_rey, xv_imy);                             // a_re * y -/+ x_imy
        #endif
    };

    force_inline
    static simd_type eval(const simd_real_type& x, const simd_type& y)
    {
        __m128d res  = _mm_mul_pd(x.data, y.data.data);
        return res;
    }

    force_inline
    static simd_type eval(const simd_type& x, const simd_real_type& y)
    {
        __m128d res  = _mm_mul_pd(x.data.data, y.data);
        return res;
    }
};

template<>
struct simd_compl_div<double, 128, sse_tag>
{
    using simd_type         = simd_compl<double, 128, sse_tag>;
    using simd_real         = simd<double, 128, sse_tag>;

    // (a.re * b.re + a.im * b.im, b.re * a.im - a.re * b.im) / (b.re * b.re + b.im * b.im)
    force_inline
    static simd_type eval(const simd_type& x, const simd_type& y)
    {        
        __m128d y_flip = _mm_shuffle_pd(y.data.data, y.data.data, 1);   // swap y.re and y.im
        __m128d x_im   = _mm_shuffle_pd(x.data.data, x.data.data, 3);   // imag of x in both
        __m128d x_re   = _mm_shuffle_pd(x.data.data, x.data.data, 0);   // real of x in both
        __m128d x_rey  = _mm_mul_pd(x_re, y.data.data);                 // (x.re*b.re, x.re*b.im)  
        __m128d yy     = _mm_mul_pd(y.data.data, y.data.data);          // (y.re*y.re, y.im*y.im)

        #if MATCL_ARCHITECTURE_HAS_FMA
            __m128d n      = _mm_fmsubadd_pd(x_im, y_flip, x_rey);      // (x_im * y_im, x_im * y_re) +/- x_rey
        #else
            __m128d x_imy  = _mm_mul_pd(x_im, y_flip);                  // (x_im * y_im, x_im * y_re)
            simd_real xv_imy(x_imy);
            simd_real xv_rey(x_rey);
            __m128d n      = sub_add(xv_imy, -xv_rey).data;             // x_re * y +/- x_imy
        #endif

        #if MATCL_ARCHITECTURE_HAS_SSE3
            __m128d yy2 = _mm_hadd_pd(yy,yy);                           // (y.re*y.re + y.im*y.im) 
        #else
            double s    = yy.m128d_f64[0] + yy.m128d_f64[1];
            __m128d yy2 = _mm_set1_pd(s);
        #endif

        return _mm_div_pd(n, yy2);
    };
};

template<>
struct simd_compl_plus<double, 128, sse_tag>
{
    using simd_type         = simd_compl<double, 128, sse_tag>;

    force_inline
    static simd_type eval(const simd_type& x, const simd_type& y)
    {
        return _mm_add_pd( x.data.data, y.data.data );
    };
};

template<>
struct simd_compl_minus<double, 128, sse_tag>
{
    using simd_type         = simd_compl<double, 128, sse_tag>;

    force_inline
    static simd_type eval(const simd_type& x, const simd_type& y)
    {
        return _mm_sub_pd( x.data.data, y.data.data );
    };
};

template<>
struct simd_compl_uminus<double, 128, sse_tag>
{
    using simd_type         = simd_compl<double, 128, sse_tag>;

    force_inline
    static simd_type eval(const simd_type& x)
    {
        return _mm_xor_pd( x.data.data, _mm_set1_pd(-0.0) );
    };
};

template<>
struct simd_compl_fma<double, 128, sse_tag>
{
    using simd_type         = simd_compl<double, 128, sse_tag>;
    using simd_real_type    = simd<double, 128, sse_tag>;

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
struct simd_compl_fms<double, 128, sse_tag>
{
    using simd_type         = simd_compl<double, 128, sse_tag>;
    using simd_real_type    = simd<double, 128, sse_tag>;

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
struct simd_compl_sum_all<double, 128, sse_tag>
{
    using simd_type         = simd_compl<double, 128, sse_tag>;

    force_inline
    static simd_double_complex eval(const simd_type& x)
    {
        return x.get<0>();
    };
};

}}
