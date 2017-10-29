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

#pragma once

#include "matcl-simd/arch/reg_128/simd_128.h"
#include "matcl-simd/simd/simd_128_func.h"
#include "matcl-simd/impl/simd_helpers.inl"

namespace matcl { namespace simd
{

//-------------------------------------------------------------------
//                          AVX 1 DOUBLE
//-------------------------------------------------------------------

template<>
force_inline
simd<double,reg_128> reverse(const simd<double,reg_128>& x)
{
    return _mm_shuffle_pd(x.data, x.data, _MM_SHUFFLE2(0,1));
};

template<>
force_inline
simd<double,reg_128> fma(const simd<double,reg_128>& x, const simd<double,reg_128>& y, 
                          const simd<double,reg_128>& z)
{
    return _mm_fmadd_pd( x.data, y.data, z.data);
};

template<>
force_inline
simd<double,reg_128> fms(const simd<double,reg_128>& x, const simd<double,reg_128>& y, 
                          const simd<double,reg_128>& z)
{
    return _mm_fmsub_pd( x.data, y.data, z.data);
};

template<>
force_inline
double sum_all(const simd<double,reg_128>& x)
{
    __m128d s   = _mm_hadd_pd(x.data, x.data);
    double* arr = (double*)&s;
    return arr[0];
};

template<>
force_inline
simd<double,reg_128> operator*(const simd<double,reg_128>& x, const simd<double,reg_128>& y)
{
    return _mm_mul_pd( x.data, y.data );
}

template<>
force_inline
simd<double,reg_128> operator/(const simd<double,reg_128>& x, const simd<double,reg_128>& y)
{
    return _mm_div_pd( x.data, y.data );
};

template<>
force_inline
simd<double,reg_128> operator+(const simd<double,reg_128>& x, const simd<double,reg_128>& y)
{
    return _mm_add_pd( x.data, y.data );
}

template<>
force_inline
simd<double,reg_128> operator-(const simd<double,reg_128>& x, const simd<double,reg_128>& y)
{
    return _mm_sub_pd( x.data, y.data );
}

template<>
force_inline
simd<double,reg_128> operator-(const simd<double,reg_128>& x)
{
    return _mm_xor_pd( x.data, _mm_set1_pd(-0.0) );
}

//-------------------------------------------------------------------
//                          AVX 1 SINGLE
//-------------------------------------------------------------------

template<>
force_inline
simd<float,reg_128> reverse(const simd<float,reg_128>& x)
{
    return _mm_shuffle_ps(x.data, x.data, _MM_SHUFFLE(0,1,2,3));
};

template<>
force_inline
simd<float,reg_128> fma(const simd<float,reg_128>& x, const simd<float,reg_128>& y, 
                          const simd<float,reg_128>& z)
{
    return _mm_fmadd_ps( x.data, y.data, z.data);
};

template<>
force_inline
simd<float,reg_128> fms(const simd<float,reg_128>& x, const simd<float,reg_128>& y, 
                          const simd<float,reg_128>& z)
{
    return _mm_fmsub_ps( x.data, y.data, z.data);
};

template<>
force_inline
float sum_all(const simd<float,reg_128>& x)
{
    __m128 s1   = _mm_hadd_ps(x.data, x.data);
    __m128 s2   = _mm_hadd_ps(s1,s1);

    float* arr  = (float*)&s2;
    return arr[0];
};

template<>
force_inline
simd<float,reg_128> operator*(const simd<float,reg_128>& x, const simd<float,reg_128>& y)
{
    return _mm_mul_ps( x.data, y.data );
}

template<>
force_inline
simd<float,reg_128> operator/(const simd<float,reg_128>& x, const simd<float,reg_128>& y)
{
    return _mm_div_ps( x.data, y.data );
};

template<>
force_inline
simd<float,reg_128> operator+(const simd<float,reg_128>& x, const simd<float,reg_128>& y)
{
    return _mm_add_ps( x.data, y.data );
}

template<>
force_inline
simd<float,reg_128> operator-(const simd<float,reg_128>& x, const simd<float,reg_128>& y)
{
    return _mm_sub_ps( x.data, y.data );
}

template<>
force_inline
simd<float,reg_128> operator-(const simd<float,reg_128>& x)
{
    return _mm_xor_ps( x.data, _mm_set1_ps(-0.0f) );
}

template<>
force_inline
simd<float,reg_128> abs(const simd<float,reg_128>& x)
{
    const __m128 sign_mask = _mm_set1_ps(-0.f); // -0.f = 1 << 31
    return _mm_andnot_ps(sign_mask, x.data);
}

template<>
force_inline
simd<double,reg_128> abs(const simd<double,reg_128>& x)
{
    const __m128d sign_mask = _mm_set1_pd(-0.); // -0. = 1 << 63
    return _mm_andnot_pd(sign_mask, x.data);
}

template<>
force_inline
simd<float,reg_128> max(const simd<float,reg_128>& x, const simd<float,reg_128>& y)
{
    return _mm_max_ps(x.data, y.data);
}

template<>
force_inline
simd<double,reg_128> max(const simd<double,reg_128>& x, const simd<double,reg_128>& y)
{
    return _mm_max_pd(x.data, y.data);
}

template<>
force_inline
simd<float,reg_128> min(const simd<float,reg_128>& x, const simd<float,reg_128>& y)
{
    return _mm_min_ps(x.data, y.data);
}

template<>
force_inline
simd<double,reg_128> min(const simd<double,reg_128>& x, const simd<double,reg_128>& y)
{
    return _mm_min_pd(x.data, y.data);
}

template<>
force_inline
simd<float,reg_128> sqrt(const simd<float,reg_128>& x)
{
    return _mm_sqrt_ps(x.data);
}

template<>
force_inline
simd<double,reg_128> sqrt(const simd<double,reg_128>& x)
{
    return _mm_sqrt_pd(x.data);
}

template<>
force_inline
simd<float,reg_128> round(const simd<float,reg_128>& x)
{
    return _mm_round_ps(x.data, _MM_FROUND_TO_NEAREST_INT);
}

template<>
force_inline
simd<double,reg_128> round(const simd<double,reg_128>& x)
{
    return _mm_round_pd(x.data, _MM_FROUND_TO_NEAREST_INT);
}

template<>
force_inline
simd<float,reg_128> floor(const simd<float,reg_128>& x)
{
    return _mm_round_ps(x.data, _MM_FROUND_TO_NEG_INF);
}

template<>
force_inline
simd<double,reg_128> floor(const simd<double,reg_128>& x)
{
    return _mm_round_pd(x.data, _MM_FROUND_TO_NEG_INF);
}

template<>
force_inline
simd<float,reg_128> ceil(const simd<float,reg_128>& x)
{
    return _mm_round_ps(x.data, _MM_FROUND_TO_POS_INF);
}

template<>
force_inline
simd<double,reg_128> ceil(const simd<double,reg_128>& x)
{
    return _mm_round_pd(x.data, _MM_FROUND_TO_POS_INF);
}

template<>
force_inline
simd<float,reg_128> trunc(const simd<float,reg_128>& x)
{
    return _mm_round_ps(x.data, _MM_FROUND_TO_ZERO);
}

template<>
force_inline
simd<double,reg_128> trunc(const simd<double,reg_128>& x)
{
    return _mm_round_pd(x.data, _MM_FROUND_TO_ZERO);
}


template<>
force_inline
simd<double,reg_128> eeq(const simd<double,reg_128>& x, const simd<double,reg_128>& y)
{
    return _mm_cmp_pd(x.data, y.data, _CMP_EQ_OQ);
}

template<>
force_inline
simd<double,reg_128> neq(const simd<double,reg_128>& x, const simd<double,reg_128>& y)
{
    return _mm_cmp_pd(x.data, y.data, _CMP_NEQ_OQ);
}

template<>
force_inline
simd<double,reg_128> lt(const simd<double,reg_128>& x, const simd<double,reg_128>& y)
{
    return _mm_cmp_pd(x.data, y.data, _CMP_LT_OQ);
}

template<>
force_inline
simd<double,reg_128> gt(const simd<double,reg_128>& x, const simd<double,reg_128>& y)
{
    return _mm_cmp_pd(x.data, y.data, _CMP_GT_OQ);
}

template<>
force_inline
simd<double,reg_128> leq(const simd<double,reg_128>& x, const simd<double,reg_128>& y)
{
    return _mm_cmp_pd(x.data, y.data, _CMP_LE_OQ);
}

template<>
force_inline
simd<double,reg_128> geq(const simd<double,reg_128>& x, const simd<double,reg_128>& y)
{
    return _mm_cmp_pd(x.data, y.data, _CMP_GE_OQ);
}

template<>
force_inline
simd<float,reg_128> eeq(const simd<float,reg_128>& x, const simd<float,reg_128>& y)
{
    return _mm_cmp_ps(x.data, y.data, _CMP_EQ_OQ);
}

template<>
force_inline
simd<float,reg_128> neq(const simd<float,reg_128>& x, const simd<float,reg_128>& y)
{
    return _mm_cmp_ps(x.data, y.data, _CMP_NEQ_OQ);
}

template<>
force_inline
simd<float,reg_128> lt(const simd<float,reg_128>& x, const simd<float,reg_128>& y)
{
    return _mm_cmp_ps(x.data, y.data, _CMP_LT_OQ);
}

template<>
force_inline
simd<float,reg_128> gt(const simd<float,reg_128>& x, const simd<float,reg_128>& y)
{
    return _mm_cmp_ps(x.data, y.data, _CMP_GT_OQ);
}

template<>
force_inline
simd<float,reg_128> leq(const simd<float,reg_128>& x, const simd<float,reg_128>& y)
{
    return _mm_cmp_ps(x.data, y.data, _CMP_LE_OQ);
}

template<>
force_inline
simd<float,reg_128> geq(const simd<float,reg_128>& x, const simd<float,reg_128>& y)
{
    return _mm_cmp_ps(x.data, y.data, _CMP_GE_OQ);
}

}}
