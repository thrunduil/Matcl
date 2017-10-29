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

#include "matcl-simd/arch/reg_256/simd_256.h"
#include "matcl-simd/simd/simd_256_func.h"
#include "matcl-simd/impl/simd_helpers.inl"

namespace matcl { namespace simd
{

//-------------------------------------------------------------------
//                          AVX 2 DOUBLE
//-------------------------------------------------------------------

template<>
force_inline
simd<double,reg_256> reverse(const simd<double,reg_256>& x)
{
    static const int control = 3 + (2 << 2) + (1 << 4);
    return _mm256_permute4x64_pd(x.data, control);
};

template<>
force_inline
simd<double,reg_256> operator*(const simd<double,reg_256>& x, const simd<double,reg_256>& y)
{
    return _mm256_mul_pd( x.data, y.data );
}

template<>
force_inline
simd<double,reg_256> operator/(const simd<double,reg_256>& x, const simd<double,reg_256>& y)
{
    return _mm256_div_pd( x.data, y.data );
};

template<>
force_inline
simd<double,reg_256> operator+(const simd<double,reg_256>& x, const simd<double,reg_256>& y)
{
    return _mm256_add_pd( x.data, y.data );
}

template<>
force_inline
simd<double,reg_256> operator-(const simd<double,reg_256>& x, const simd<double,reg_256>& y)
{
    return _mm256_sub_pd( x.data, y.data );
}

template<>
force_inline
simd<double,reg_256> operator-(const simd<double,reg_256>& x)
{
    double Z = -0.0;
    return _mm256_xor_pd(x.data, _mm256_broadcast_sd(&Z));
}

template<>
force_inline
simd<double,reg_256> fma(const simd<double,reg_256>& x, const simd<double,reg_256>& y, 
                          const simd<double,reg_256>& z)
{
    return _mm256_fmadd_pd( x.data, y.data, z.data);
};

template<>
force_inline
simd<double,reg_256> fms(const simd<double,reg_256>& x, const simd<double,reg_256>& y, 
                          const simd<double,reg_256>& z)
{
    return _mm256_fmsub_pd( x.data, y.data, z.data);
};

template<>
force_inline
double sum_all(const simd<double,reg_256>& x)
{
    __m256d sum     = _mm256_hadd_pd(x.data, x.data);
    double* arr     = (double*)&sum;
    return arr[0] + arr[2];
};

//-------------------------------------------------------------------
//                          AVX 2 SINGLE
//-------------------------------------------------------------------

template<>
force_inline
simd<float,reg_256> reverse(const simd<float,reg_256>& x)
{
    static const __m256i control = vector_8_int<0,1,2,3,4,5,6,7>();
    return _mm256_permutevar8x32_ps(x.data, control);
};

template<>
force_inline
simd<float,reg_256> operator*(const simd<float,reg_256>& x, const simd<float,reg_256>& y)
{
    return _mm256_mul_ps( x.data, y.data );
}

template<>
force_inline
simd<float,reg_256> operator/(const simd<float,reg_256>& x, const simd<float,reg_256>& y)
{
    return _mm256_div_ps( x.data, y.data );
};

template<>
force_inline
simd<float,reg_256> operator+(const simd<float,reg_256>& x, const simd<float,reg_256>& y)
{
    return _mm256_add_ps( x.data, y.data );
}

template<>
force_inline
simd<float,reg_256> operator-(const simd<float,reg_256>& x, const simd<float,reg_256>& y)
{
    return _mm256_sub_ps( x.data, y.data );
}

template<>
force_inline
simd<float,reg_256> operator-(const simd<float,reg_256>& x)
{
    float Z = -0.0f;
    return _mm256_xor_ps(x.data, _mm256_broadcast_ss(&Z));
}

template<>
force_inline
simd<float,reg_256> fma(const simd<float,reg_256>& x, const simd<float,reg_256>& y, 
                          const simd<float,reg_256>& z)
{
    return _mm256_fmadd_ps( x.data, y.data, z.data);
};

template<>
force_inline
simd<float,reg_256> fms(const simd<float,reg_256>& x, const simd<float,reg_256>& y, 
                          const simd<float,reg_256>& z)
{
    return _mm256_fmsub_ps( x.data, y.data, z.data);
};

template<>
force_inline
float sum_all(const simd<float,reg_256>& x)
{
    __m256 s    = _mm256_hadd_ps(x.data, x.data);
    __m256 s2   = _mm256_hadd_ps(s,s);

    float* arr  = (float*)&s2;
    return arr[0] + arr[4];
};

template<>
force_inline
simd<float,reg_256> abs(const simd<float,reg_256>& x)
{
    const __m256 sign_mask = _mm256_set1_ps(-0.f); // -0.f = 1 << 31
    return _mm256_andnot_ps(sign_mask, x.data);
}

template<>
force_inline
simd<double,reg_256> abs(const simd<double,reg_256>& x)
{
    const __m256d sign_mask = _mm256_set1_pd(-0.); // -0. = 1 << 63
    return _mm256_andnot_pd(sign_mask, x.data);
}

template<>
force_inline
simd<float,reg_256> max(const simd<float,reg_256>& x, const simd<float,reg_256>& y)
{
    return _mm256_max_ps(x.data, y.data);
}

template<>
force_inline
simd<double,reg_256> max(const simd<double,reg_256>& x, const simd<double,reg_256>& y)
{
    return _mm256_max_pd(x.data, y.data);
}

template<>
force_inline
simd<float,reg_256> min(const simd<float,reg_256>& x, const simd<float,reg_256>& y)
{
    return _mm256_min_ps(x.data, y.data);
}

template<>
force_inline
simd<double,reg_256> min(const simd<double,reg_256>& x, const simd<double,reg_256>& y)
{
    return _mm256_min_pd(x.data, y.data);
}

template<>
force_inline
simd<float,reg_256> sqrt(const simd<float,reg_256>& x)
{
    return _mm256_sqrt_ps(x.data);
}

template<>
force_inline
simd<double,reg_256> sqrt(const simd<double,reg_256>& x)
{
    return _mm256_sqrt_pd(x.data);
}

template<>
force_inline
simd<float,reg_256> round(const simd<float,reg_256>& x)
{
    return _mm256_round_ps(x.data, _MM_FROUND_TO_NEAREST_INT);
}

template<>
force_inline
simd<double,reg_256> round(const simd<double,reg_256>& x)
{
    return _mm256_round_pd(x.data, _MM_FROUND_TO_NEAREST_INT);
}

template<>
force_inline
simd<float,reg_256> floor(const simd<float,reg_256>& x)
{
    return _mm256_round_ps(x.data, _MM_FROUND_TO_NEG_INF);
}

template<>
force_inline
simd<double,reg_256> floor(const simd<double,reg_256>& x)
{
    return _mm256_round_pd(x.data, _MM_FROUND_TO_NEG_INF);
}

template<>
force_inline
simd<float,reg_256> ceil(const simd<float,reg_256>& x)
{
    return _mm256_round_ps(x.data, _MM_FROUND_TO_POS_INF);
}

template<>
force_inline
simd<double,reg_256> ceil(const simd<double,reg_256>& x)
{
    return _mm256_round_pd(x.data, _MM_FROUND_TO_POS_INF);
}

template<>
force_inline
simd<float,reg_256> trunc(const simd<float,reg_256>& x)
{
    return _mm256_round_ps(x.data, _MM_FROUND_TO_ZERO);
}

template<>
force_inline
simd<double,reg_256> trunc(const simd<double,reg_256>& x)
{
    return _mm256_round_pd(x.data, _MM_FROUND_TO_ZERO);
}

template<>
force_inline
simd<double,reg_256> eeq(const simd<double,reg_256>& x, const simd<double,reg_256>& y)
{
    return _mm256_cmp_pd(x.data, y.data, _CMP_EQ_OQ);
}

template<>
force_inline
simd<double,reg_256> neq(const simd<double,reg_256>& x, const simd<double,reg_256>& y)
{
    return _mm256_cmp_pd(x.data, y.data, _CMP_NEQ_OQ);
}

template<>
force_inline
simd<double,reg_256> lt(const simd<double,reg_256>& x, const simd<double,reg_256>& y)
{
    return _mm256_cmp_pd(x.data, y.data, _CMP_LT_OQ);
}

template<>
force_inline
simd<double,reg_256> gt(const simd<double,reg_256>& x, const simd<double,reg_256>& y)
{
    return _mm256_cmp_pd(x.data, y.data, _CMP_GT_OQ);
}

template<>
force_inline
simd<double,reg_256> leq(const simd<double,reg_256>& x, const simd<double,reg_256>& y)
{
    return _mm256_cmp_pd(x.data, y.data, _CMP_LE_OQ);
}

template<>
force_inline
simd<double,reg_256> geq(const simd<double,reg_256>& x, const simd<double,reg_256>& y)
{
    return _mm256_cmp_pd(x.data, y.data, _CMP_GE_OQ);
}

template<>
force_inline
simd<float,reg_256> eeq(const simd<float,reg_256>& x, const simd<float,reg_256>& y)
{
    return _mm256_cmp_ps(x.data, y.data, _CMP_EQ_OQ);
}

template<>
force_inline
simd<float,reg_256> neq(const simd<float,reg_256>& x, const simd<float,reg_256>& y)
{
    return _mm256_cmp_ps(x.data, y.data, _CMP_NEQ_OQ);
}

template<>
force_inline
simd<float,reg_256> lt(const simd<float,reg_256>& x, const simd<float,reg_256>& y)
{
    return _mm256_cmp_ps(x.data, y.data, _CMP_LT_OQ);
}

template<>
force_inline
simd<float,reg_256> gt(const simd<float,reg_256>& x, const simd<float,reg_256>& y)
{
    return _mm256_cmp_ps(x.data, y.data, _CMP_GT_OQ);
}

template<>
force_inline
simd<float,reg_256> leq(const simd<float,reg_256>& x, const simd<float,reg_256>& y)
{
    return _mm256_cmp_ps(x.data, y.data, _CMP_LE_OQ);
}

template<>
force_inline
simd<float,reg_256> geq(const simd<float,reg_256>& x, const simd<float,reg_256>& y)
{
    return _mm256_cmp_ps(x.data, y.data, _CMP_GE_OQ);
}

}}
