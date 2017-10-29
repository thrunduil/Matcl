
/*
 *  This file is a part of Morfa Matrix Lib.
 *
 *  Copyright (c) Pawe³ Kowal 2015-2016
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

#include "matcl-simd/arch/reg_128/simd_128_compl.h"
#include "matcl-simd/simd/simd_128_compl_func.h"
#include "matcl-simd/impl/simd_helpers.inl"

namespace matcl { namespace simd
{

template<>
force_inline
simd_compl<double,reg_128> reverse(const simd_compl<double,reg_128>& x)
{
    return x;
};
template<>
force_inline
simd_compl<float,reg_128> reverse(const simd_compl<float,reg_128>& x)
{
    return _mm_shuffle_ps(x.data.data, x.data.data, _MM_SHUFFLE(1,0,3,2));
};

template<>
force_inline
simd_compl<float,reg_128> operator+(const simd_compl<float,reg_128>& x, 
                                       const simd_compl<float,reg_128>& y)
{
    return _mm_add_ps( x.data.data, y.data.data );
}

template<>
force_inline
simd_compl<double,reg_128> operator+(const simd_compl<double,reg_128>& x, 
                                        const simd_compl<double,reg_128>& y)
{
    //TODO
    return _mm_add_pd( x.data.data, y.data.data );
}

template<>
force_inline
simd_compl<double,reg_128> operator-(const simd_compl<double,reg_128>& x, 
                                        const simd_compl<double,reg_128>& y)
{
    //TODO
    return _mm_sub_pd( x.data.data, y.data.data );
}

template<>
force_inline
simd_compl<double,reg_128> operator-(const simd_compl<double,reg_128>& x)
{
    //TODO
    return _mm_xor_pd( x.data.data, _mm_set1_pd(-0.0) );
}

// (x.re * y.re - x.im * y.im,  x.re * y.im + x.re * y.im)
template<>
force_inline
simd_compl<double,reg_128> operator*(const simd_compl<double,reg_128>& x, 
                                        const simd_compl<double,reg_128>& y)
{
    //TODO
    __m128d y_flip = _mm_shuffle_pd(y.data.data, y.data.data, 1);   // swap y.re and y.im
    __m128d x_im   = _mm_shuffle_pd(x.data.data, x.data.data, 3);   // imag of x in both
    __m128d x_re   = _mm_shuffle_pd(x.data.data, x.data.data, 0);   // real of x in both
    __m128d x_imy  = _mm_mul_pd(x_im, y_flip);            // (x.im*y.im, x.im*y.re)

    return  _mm_fmaddsub_pd(x_re, y.data.data, x_imy);         // a_re * y +/- x_imy
}

template<>
force_inline
simd_compl<float,reg_128> operator*(const simd<float,reg_128>& x, 
                                       const simd_compl<float,reg_128>& y)
{
    __m128 res  = _mm_mul_ps(x.data, y.data.data);
    return res;
}

template<>
force_inline
simd_compl<float,reg_128> operator*(const simd_compl<float,reg_128>& x,
                                    const simd<float,reg_128>& y)
{
    __m128 res  = _mm_mul_ps(x.data.data, y.data);
    return res;
}

template<>
force_inline
simd_compl<double,reg_128> operator*(const simd<double,reg_128>& x, 
                                       const simd_compl<double,reg_128>& y)
{
    //TODO
    __m128d res  = _mm_mul_pd(x.data, y.data.data);
    return res;
}

template<>
force_inline
simd_compl<double,reg_128> operator*(const simd_compl<double,reg_128>& x,
                                    const simd<double,reg_128>& y)
{
    //TODO
    __m128d res  = _mm_mul_pd(x.data.data, y.data);
    return res;
}

// (a.re * b.re + a.im * b.im, b.re * a.im - a.re * b.im) / (b.re * b.re + b.im * b.im)
template<>
force_inline
simd_compl<double,reg_128> operator/(const simd_compl<double,reg_128>& x, 
                                        const simd_compl<double,reg_128>& y)
{
    //TODO
    __m128d y_flip = _mm_shuffle_pd(y.data.data, y.data.data, 1);      // swap y.re and y.im
    __m128d x_im   = _mm_shuffle_pd(x.data.data, x.data.data, 3);      // imag of x in both
    __m128d x_re   = _mm_shuffle_pd(x.data.data, x.data.data, 0);      // real of x in both
    __m128d x_rey  = _mm_mul_pd(x_re, y.data.data);               // (x.re*b.re, x.re*b.im)  

    __m128d n      = _mm_fmsubadd_pd(x_im, y_flip, x_rey);   // x_re * y +/- x_imy
    __m128d yy     = _mm_mul_pd(y.data.data, y.data.data);             // (y.re*y.re, y.im*y.im)
    __m128d yy2    = _mm_hadd_pd(yy,yy);                     // (y.re*y.re + y.im*y.im) 

    return _mm_div_pd(n, yy2);
}

template<>
force_inline
simd_compl<float,reg_128> operator-(const simd_compl<float,reg_128>& x,
                                       const simd_compl<float,reg_128>& y)
{
    return _mm_sub_ps( x.data.data, y.data.data );
}

template<>
force_inline
simd_compl<float,reg_128> operator-(const simd_compl<float,reg_128>& x)
{
    return _mm_xor_ps( x.data.data, _mm_set1_ps(-0.0f) );
}

// (x.re * y.re - x.im * y.im,  x.re * y.im + x.re * y.im)
template<>
force_inline
simd_compl<float,reg_128> operator*(const simd_compl<float,reg_128>& x, 
                                       const simd_compl<float,reg_128>& y)
{
    __m128 y_flip = _mm_shuffle_ps(y.data.data, y.data.data, 0xB1);      // swap y.re and y.im
    __m128 x_im   = _mm_shuffle_ps(x.data.data, x.data.data, 0xF5);      // imag of x in both
    __m128 x_re   = _mm_shuffle_ps(x.data.data, x.data.data, 0xA0);      // real of x in both
    __m128 x_imy  = _mm_mul_ps(x_im, y_flip);                  // (x.im*y.im, x.im*y.re)

    return  _mm_fmaddsub_ps(x_re, y.data.data, x_imy);              // a_re * y +/- x_imy
}

// (a.re * b.re + a.im * b.im, b.re * a.im - a.re * b.im) / (b.re * b.re + b.im * b.im)
template<>
force_inline
simd_compl<float,reg_128> operator/(const simd_compl<float,reg_128>& x, 
                                       const simd_compl<float,reg_128>& y)
{
    __m128 y_flip = _mm_shuffle_ps(y.data.data, y.data.data, 0xB1);      // swap y.re and y.im
    __m128 x_im   = _mm_shuffle_ps(x.data.data, x.data.data, 0xF5);      // imag of x in both
    __m128 x_re   = _mm_shuffle_ps(x.data.data, x.data.data, 0xA0);      // real of x in both
    __m128 x_rey  = _mm_mul_ps(x_re, y.data.data);                  // (x.re*b.re, x.re*b.im)  

    __m128 n      = _mm_fmsubadd_ps(x_im, y_flip, x_rey);   // x_re * y +/- x_imy
    __m128 yy     = _mm_mul_ps(y.data.data, y.data.data);             // (y.re*y.re, y.im*y.im)
    __m128 yy2    = _mm_shuffle_ps(yy,yy, 0xB1);            // swap yy.re and yy.im
    __m128 yy3    = _mm_add_ps(yy,yy2);                     // add pairwise

    return _mm_div_ps(n, yy3);
}

//res = x*y + z
template<class Val>
force_inline
simd_compl<Val,reg_128>
fma(const simd_compl<Val,reg_128>& x, const simd_compl<Val,reg_128>& y, 
                       const simd_compl<Val,reg_128>& z)
{
    return x * y + z;
}

template<class Val>
force_inline
simd_compl<Val,reg_128>
fma(const simd<Val,reg_128>& x, const simd_compl<Val,reg_128>& y, 
                       const simd_compl<Val,reg_128>& z)
{
    return x * y + z;
}

//res = x*y - z
template<class Val>
force_inline
simd_compl<Val,reg_128>
fms(const simd_compl<Val,reg_128>& x, const simd_compl<Val,reg_128>& y, 
                       const simd_compl<Val,reg_128>& z)
{
    return x * y - z;
}

template<class Val>
force_inline
simd_compl<Val,reg_128>
fms(const simd<Val,reg_128>& x, const simd_compl<Val,reg_128>& y, 
                       const simd_compl<Val,reg_128>& z)
{
    return x * y - z;
}

template<>
force_inline
simd_double_complex sum_all(const simd_compl<double, reg_128>& x)
{
    return simd_double_complex(x.data.data.m128d_f64[0], x.data.data.m128d_f64[1]);
};

template<>
force_inline
simd_single_complex sum_all(const simd_compl<float, reg_128>& x)
{
    using simd_compl = simd_compl<float, reg_128>;

    simd_compl res = x + reverse(x);
    return simd_single_complex(res.data.data.m128_f32[0], res.data.data.m128_f32[1]);
};

}}
