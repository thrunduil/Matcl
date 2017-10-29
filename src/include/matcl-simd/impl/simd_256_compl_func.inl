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

#include "matcl-simd/arch/reg_256/simd_256_compl.h"
#include "matcl-simd/simd/simd_256_compl_func.h"
#include "matcl-simd/impl/simd_helpers.inl"

namespace matcl { namespace simd
{

template<>
force_inline
simd_compl<double,reg_256> reverse(const simd_compl<double,reg_256>& x)
{
    static const int control = 2 + (3 << 2) + (0 << 4) + (1 << 6);
    return _mm256_permute4x64_pd(x.data.data, control);
};

template<>
force_inline
simd_compl<float,reg_256> reverse(const simd_compl<float,reg_256>& x)
{
    static const __m256i control = vector_8_int<1,0,3,2,5,4,7,6>();
    return _mm256_permutevar8x32_ps(x.data.data, control);
};

template<>
force_inline
simd_compl<double,reg_256> operator+(const simd_compl<double,reg_256>& x, 
                                        const simd_compl<double,reg_256>& y)
{
    return _mm256_add_pd( x.data.data, y.data.data );
}

template<>
force_inline
simd_compl<double,reg_256> operator-(const simd_compl<double,reg_256>& x,
                                        const simd_compl<double,reg_256>& y)
{
    return _mm256_sub_pd( x.data.data, y.data.data );
}
template<>
force_inline
simd_compl<double,reg_256> operator-(const simd_compl<double,reg_256>& x)
{
    double Z = -0.0;
    return _mm256_xor_pd(x.data.data, _mm256_broadcast_sd(&Z));
}

// (x.re * y.re - x.im * y.im,  x.re * y.im + x.re * y.im)
template<>
force_inline
simd_compl<double,reg_256> operator*(const simd_compl<double,reg_256>& x, 
                                        const simd_compl<double,reg_256>& y)
{
    __m256d y_flip = _mm256_shuffle_pd(y.data.data, y.data.data, 5);   // swap y.re and y.im
    __m256d x_im   = _mm256_shuffle_pd(x.data.data, x.data.data, 0xF); // imag of x in both
    __m256d x_re   = _mm256_shuffle_pd(x.data.data, x.data.data, 0);   // real of x in both
    __m256d x_imy  = _mm256_mul_pd(x_im, y_flip);            // (x.im*y.im, x.im*y.re)

    return  _mm256_fmaddsub_pd(x_re, y.data.data, x_imy);         // a_re * y +/- x_imy
}

// (a.re * b.re + a.im * b.im, b.re * a.im - a.re * b.im) / (b.re * b.re + b.im * b.im)
template<>
force_inline
simd_compl<double,reg_256> operator/(const simd_compl<double,reg_256>& x,
                                        const simd_compl<double,reg_256>& y)
{
    __m256d y_flip = _mm256_shuffle_pd(y.data.data, y.data.data, 5);      // swap y.re and y.im
    __m256d x_im   = _mm256_shuffle_pd(x.data.data, x.data.data, 0xF);    // imag of x in both
    __m256d x_re   = _mm256_shuffle_pd(x.data.data, x.data.data, 0);      // real of x in both
    __m256d x_rey  = _mm256_mul_pd(x_re, y.data.data);               // (x.re*b.re, x.re*b.im)  

    __m256d n      = _mm256_fmsubadd_pd(x_im, y_flip, x_rey);   // x_re * y +/- x_imy
    __m256d yy     = _mm256_mul_pd(y.data.data, y.data.data);             // (y.re*y.re, y.im*y.im)
    __m256d yy2    = _mm256_hadd_pd(yy,yy);                     // (y.re*y.re + y.im*y.im) 

    return _mm256_div_pd(n, yy2);
}

//-------------------------------------------------------------------
//                          AVX FLOAT COMPLEX
//-------------------------------------------------------------------
template<>
force_inline
simd_compl<float,reg_256> operator+(const simd_compl<float,reg_256>& x,
                                       const simd_compl<float,reg_256>& y)
{
    return _mm256_add_ps( x.data.data, y.data.data );
}

template<>
force_inline
simd_compl<float,reg_256> operator-(const simd_compl<float,reg_256>& x,
                                       const simd_compl<float,reg_256>& y)
{
    return _mm256_sub_ps( x.data.data, y.data.data );
}

template<>
force_inline
simd_compl<float,reg_256> operator-(const simd_compl<float,reg_256>& x)
{
    float Z = -0.0f;
    return _mm256_xor_ps(x.data.data, _mm256_broadcast_ss(&Z));
}

// (x.re * y.re - x.im * y.im,  x.re * y.im + x.re * y.im)
template<>
force_inline
simd_compl<float,reg_256> operator*(const simd_compl<float,reg_256>& x, 
                                       const simd_compl<float,reg_256>& y)
{
    __m256 y_flip = _mm256_shuffle_ps(y.data.data, y.data.data, 0xB1);   // swap y.re and y.im
    __m256 x_im   = _mm256_shuffle_ps(x.data.data, x.data.data, 0xF5);   // imag of x in both
    __m256 x_re   = _mm256_shuffle_ps(x.data.data, x.data.data, 0xA0);   // real of x in both
    __m256 x_imy  = _mm256_mul_ps(x_im, y_flip);               // (x.im*y.im, x.im*y.re)

    return  _mm256_fmaddsub_ps(x_re, y.data.data, x_imy);           // a_re * y +/- x_imy
}

template<>
force_inline
simd_compl<float,reg_256> operator*(const simd<float,reg_256>& x, 
                                       const simd_compl<float,reg_256>& y)
{
    __m256 res  = _mm256_mul_ps(x.data, y.data.data);
    return res;
}

template<>
force_inline
simd_compl<float,reg_256> operator*(const simd_compl<float,reg_256>& x,
                                    const simd<float,reg_256>& y)
{
    __m256 res  = _mm256_mul_ps(x.data.data, y.data);
    return res;
}

template<>
force_inline
simd_compl<double,reg_256> operator*(const simd<double,reg_256>& x, 
                                       const simd_compl<double,reg_256>& y)
{
    __m256d res  = _mm256_mul_pd(x.data, y.data.data);
    return res;
}

template<>
force_inline
simd_compl<double,reg_256> operator*(const simd_compl<double,reg_256>& x,
                                    const simd<double,reg_256>& y)
{
    __m256d res  = _mm256_mul_pd(x.data.data, y.data);
    return res;
}

// (a.re * b.re + a.im * b.im, b.re * a.im - a.re * b.im) / (b.re * b.re + b.im * b.im)
template<>
force_inline
simd_compl<float,reg_256> operator/(const simd_compl<float,reg_256>& x,
                                       const simd_compl<float,reg_256>& y)
{
    __m256 y_flip = _mm256_shuffle_ps(y.data.data, y.data.data, 0xB1);   // swap y.re and y.im
    __m256 x_im   = _mm256_shuffle_ps(x.data.data, x.data.data, 0xF5);   // imag of x in both
    __m256 x_re   = _mm256_shuffle_ps(x.data.data, x.data.data, 0xA0);   // real of x in both
    __m256 x_rey  = _mm256_mul_ps(x_re, y.data.data);               // (x.re*b.re, x.re*b.im)  

    __m256 n      = _mm256_fmsubadd_ps(x_im, y_flip, x_rey);   // x_re * y +/- x_imy
    __m256 yy     = _mm256_mul_ps(y.data.data, y.data.data);             // (y.re*y.re, y.im*y.im)
    __m256 yy2    = _mm256_shuffle_ps(yy,yy,0xB1);             // Swap yy.re and yy.im
    __m256 yy3    = _mm256_add_ps(yy,yy2);                     // (y.re*y.re + y.im*y.im) dublicated

    return _mm256_div_ps(n, yy3);
}

//res = x*y + z
template<class Val>
force_inline
simd_compl<Val,reg_256>
fma(const simd_compl<Val,reg_256>& x, const simd_compl<Val,reg_256>& y, 
                         const simd_compl<Val,reg_256>& z)
{
    return x*y + z;
}
template<class Val>
force_inline
simd_compl<Val,reg_256>
fma(const simd<Val,reg_256>& x, const simd_compl<Val,reg_256>& y, 
                         const simd_compl<Val,reg_256>& z)
{
    return x*y + z;
}

//res = x*y - z
template<class Val>
force_inline
simd_compl<Val,reg_256>
fms(const simd_compl<Val,reg_256>& x, const simd_compl<Val,reg_256>& y, 
                         const simd_compl<Val,reg_256>& z)
{
    return x*y - z;
}

template<class Val>
force_inline
simd_compl<Val,reg_256>
fms(const simd<Val,reg_256>& x, const simd_compl<Val,reg_256>& y, 
                         const simd_compl<Val,reg_256>& z)
{
    return x*y - z;
}

template<>
force_inline
simd_double_complex sum_all(const simd_compl<double, reg_256>& x)
{
    using simd_compl = simd_compl<double, reg_256>;
    simd_compl res = x + reverse(x);
    return simd_double_complex(res.data.data.m256d_f64[0], res.data.data.m256d_f64[1]);
};

template<>
force_inline
simd_single_complex sum_all(const simd_compl<float, reg_256>& x)
{
    using simd_compl = simd_compl<float, reg_256>;

    __m128 lo   = _mm256_extractf128_ps(x.data.data, 0);
    __m128 hi   = _mm256_extractf128_ps(x.data.data, 1);

    __m128 sum  = _mm_add_ps(lo, hi);
    __m128 rev  = _mm_shuffle_ps(sum, sum, _MM_SHUFFLE(1,0,3,2));
    __m128 sum2 = _mm_add_ps(sum, rev);

    return simd_single_complex(sum2.m128_f32[0], sum2.m128_f32[1]);
};

}}
