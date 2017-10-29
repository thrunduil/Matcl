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

#pragma once

#include "matcl-simd/arch/reg_256/simd_256_compl.h"
#include "matcl-simd/simd/simd_256_compl_func.h"

namespace matcl { namespace simd
{

//-------------------------------------------------------------------
//                          AVX DOUBLE COMPLEX
//-------------------------------------------------------------------
force_inline
simd_compl<double, reg_256>::simd_compl(const impl_type& v)
    : data(v) 
{};

force_inline
simd_compl<double, reg_256>::simd_compl(Integer re)
    :simd_compl(double(re))
{};

force_inline
simd_compl<double, reg_256>::simd_compl(float re)
    :simd_compl(double(re))
{};

force_inline
simd_compl<double, reg_256>::simd_compl(double re)
    :simd_compl(simd_half(re))
{};

force_inline
simd_compl<double, reg_256>::simd_compl(const simd_half& lo_hi)
    : data(lo_hi.data)
{}

force_inline
simd_compl<double, reg_256>::simd_compl(double re, double im)
    :simd_compl(simd_half(re, im))
{};

force_inline
simd_compl<double, reg_256>::simd_compl(const simd_double_complex& val)
    :simd_compl(simd_half(val))
{};

force_inline
simd_compl<double, reg_256>::simd_compl(const simd_single_complex& val)
    :simd_compl(simd_half(val))
{};

force_inline
simd_compl<double, reg_256>::simd_compl(const simd_double_complex& v0, const simd_double_complex& v1)
    : simd_compl(simd_half(v0), simd_half(v1))
{};

force_inline
simd_compl<double, reg_256>::simd_compl(double re_0, double im_0, double re_1, double im_1)
    : data(re_0, im_0, re_1, im_1)
{};

force_inline
simd_compl<double, reg_256> simd_compl<double, reg_256>::zero()
{
    return impl_type::zero();
}

force_inline
simd_compl<double, reg_256> simd_compl<double, reg_256>::broadcast(const simd_double_complex* arr)
{
    return simd_compl(arr[0]);
}

force_inline
simd_compl<double, reg_256> simd_compl<double, reg_256>::broadcast(const simd_double_complex& arr)
{
    return simd_compl(arr);
}

force_inline
simd_compl<double, reg_256> simd_compl<double, reg_256>::broadcast(const double* arr)
{
    return simd_compl(arr[0]);
}

force_inline
simd_compl<double, reg_256> simd_compl<double, reg_256>::broadcast(const double& arr)
{
    return simd_compl(arr);
}

force_inline simd_compl<double,reg_256> 
simd_compl<double,reg_256>::load(const simd_double_complex* arr, std::true_type aligned)
{
    return impl_type::load((const double*)arr, aligned);
};

force_inline simd_compl<double,reg_256> 
simd_compl<double,reg_256>::load(const simd_double_complex* arr, std::false_type not_aligned)
{
    return impl_type::load((const double*)arr, not_aligned);
};

force_inline void 
simd_compl<double,reg_256>::store(simd_double_complex* arr, std::true_type aligned) const
{
    data.store((double*)arr, aligned);
};

force_inline void 
simd_compl<double,reg_256>::store(simd_double_complex* arr, std::false_type not_aligned) const
{
    data.store((double*)arr, not_aligned);
};

force_inline simd_compl<double,reg_256> 
simd_compl<double,reg_256>::load_reverse(const simd_double_complex* arr, std::true_type aligned)
{
    static const int off            = vector_size - 1;
    simd_compl<double,reg_256> ret  = load(arr - off, aligned);
    return reverse(ret);
}

force_inline simd_compl<double,reg_256> 
simd_compl<double,reg_256>::load_reverse(const simd_double_complex* arr, std::false_type not_aligned)
{
    static const int off            = vector_size - 1;
    simd_compl<double,reg_256> ret  = load(arr - off, not_aligned);
    return reverse(ret);
};

force_inline void 
simd_compl<double,reg_256>::store_reverse(simd_double_complex* arr, std::true_type aligned) const
{
    static const int off            = vector_size - 1;
    simd_compl<double,reg_256> rev  = reverse(*this);
    rev.store(arr - off, aligned);
}
force_inline void 
simd_compl<double,reg_256>::store_reverse(simd_double_complex* arr, std::false_type not_aligned) const
{
    static const int off            = vector_size - 1;
    simd_compl<double,reg_256> rev  = reverse(*this);
    rev.store(arr - off, not_aligned);
};

force_inline
simd_double_complex simd_compl<double, reg_256>::get(int pos) const
{
    return simd_double_complex(data.data.m256d_f64[pos*2], data.data.m256d_f64[pos*2 + 1]);
};

template<int Step>
force_inline
void simd_compl<double,reg_256>::scatter(simd_double_complex* arr0) const
{
    //no scatter intrinsic
    double* arr = (double*) arr0;

    arr + 0*Step*2      = data.data.m256d_f64[0];
    arr + 0*Step*2+1    = data.data.m256d_f64[1];
    arr + 1*Step*2      = data.data.m256d_f64[2];
    arr + 1*Step*2+1    = data.data.m256d_f64[3];
};

force_inline simd_compl<double,reg_256>::simd_half
simd_compl<double,reg_256>::extract_low() const
{
    return data.extract_low();
}

force_inline simd_compl<double,reg_256>::simd_half 
simd_compl<double,reg_256>::extract_high() const
{
    return data.extract_high();
}

//-------------------------------------------------------------------
//                          AVX FLOAT COMPLEX
//-------------------------------------------------------------------

force_inline
simd_compl<float, reg_256>::simd_compl(const impl_type& v)  
    : data(v) 
{};

force_inline
simd_compl<float, reg_256>::simd_compl(float re)
    :simd_compl(simd_half(re))
{}

force_inline
simd_compl<float, reg_256>::simd_compl(const simd_half& lo_hi)
    : data(lo_hi.data)
{}

force_inline
simd_compl<float, reg_256>::simd_compl(float re, float im)
    : simd_compl(simd_half(re, im))
{}

force_inline
simd_compl<float, reg_256>::simd_compl(const simd_single_complex& val)
    : simd_compl(simd_half(val))
{}

force_inline
simd_compl<float, reg_256>::simd_compl(const simd_single_complex& v0, const simd_single_complex& v1,
                   const simd_single_complex& v2, const simd_single_complex& v3)
    : simd_compl(simd_half(v0, v1), simd_half(v2, v3))
{}

force_inline
simd_compl<float, reg_256>::simd_compl(float re_0, float im_0, float re_1, float im_1,
                   float re_2, float im_2, float re_3, float im_3)
    : data(re_0, im_0, re_1, im_1, re_2, im_2, re_3, im_3)
{};

force_inline
simd_compl<float, reg_256> simd_compl<float, reg_256>::zero()
{
    return impl_type::zero();
}

force_inline
simd_compl<float, reg_256> simd_compl<float, reg_256>::broadcast(const simd_single_complex* arr)
{
    return simd_compl(arr[0]);
}

force_inline
simd_compl<float, reg_256> simd_compl<float, reg_256>::broadcast(const simd_single_complex& arr)
{
    return simd_compl(arr);
};

force_inline
simd_compl<float, reg_256> simd_compl<float, reg_256>::broadcast(const float* arr)
{
    return simd_compl(arr[0]);
}

force_inline
simd_compl<float, reg_256> simd_compl<float, reg_256>::broadcast(const float& arr)
{
    return simd_compl(arr);
};

force_inline simd_compl<float,reg_256> 
simd_compl<float,reg_256>::load(const simd_single_complex* arr, std::true_type aligned)
{
    return impl_type::load((const float*)arr, aligned);
};

force_inline simd_compl<float,reg_256> 
simd_compl<float,reg_256>::load(const simd_single_complex* arr, std::false_type not_aligned)
{
    return impl_type::load((const float*)arr, not_aligned);
};

force_inline void 
simd_compl<float,reg_256>::store(simd_single_complex* arr, std::true_type aligned) const
{
    data.store((float*)arr, aligned);
};

force_inline void 
simd_compl<float,reg_256>::store(simd_single_complex* arr, std::false_type not_aligned) const
{
    data.store((float*)arr, not_aligned);
};

force_inline simd_compl<float,reg_256> 
simd_compl<float,reg_256>::load_reverse(const simd_single_complex* arr, std::true_type aligned)
{
    static const int off            = vector_size - 1;
    simd_compl<float,reg_256> ret   = load(arr - off, aligned);
    return reverse(ret);
}

force_inline simd_compl<float,reg_256> 
simd_compl<float,reg_256>::load_reverse(const simd_single_complex* arr, std::false_type not_aligned)
{
    static const int off            = (vector_size - 1)*2;
    simd_compl<float,reg_256> ret   = load(arr - off, not_aligned);
    return reverse(ret);
};

force_inline void 
simd_compl<float,reg_256>::store_reverse(simd_single_complex* arr, std::true_type aligned) const
{
    static const int off            = vector_size - 1;
    simd_compl<float,reg_256> rev   = reverse(*this);
    rev.store(arr - off, aligned);
}

force_inline void 
simd_compl<float,reg_256>::store_reverse(simd_single_complex* arr, std::false_type not_aligned) const
{
    static const int off            = vector_size - 1;
    simd_compl<float,reg_256> rev   = reverse(*this);
    rev.store(arr - off, not_aligned);
};

force_inline
simd_single_complex simd_compl<float, reg_256>::get(int pos) const
{
    return simd_single_complex(data.data.m256_f32[pos*2], data.data.m256_f32[pos*2 + 1]);
};

template<int Step>
force_inline
void simd_compl<float,reg_256>::scatter(simd_single_complex* arr0) const
{
    float* arr      = (float*)arr0;

    //no scatter intrinsic
    arr + 0*Step*2      = data.m256_f32[0];
    arr + 0*Step*2+1    = data.m256_f32[1];
    arr + 1*Step*2      = data.m256_f32[2];
    arr + 1*Step*2+1    = data.m256_f32[3];
    arr + 2*Step*2      = data.m256_f32[4];
    arr + 2*Step*2+1    = data.m256_f32[5];
    arr + 3*Step*2      = data.m256_f32[6];
    arr + 3*Step*2+1    = data.m256_f32[7];
};

force_inline simd_compl<float,reg_256>::simd_half
simd_compl<float,reg_256>::extract_low() const
{
    return data.extract_low();
}

force_inline simd_compl<float,reg_256>::simd_half 
simd_compl<float,reg_256>::extract_high() const
{
    return data.extract_high();
}

}}
