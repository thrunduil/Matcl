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

#include "matcl-simd/arch/reg_128/simd_128_compl.h"
#include "matcl-simd/simd/simd_128_compl_func.h"

namespace matcl { namespace simd
{

//-------------------------------------------------------------------
//                          SSE2 DOUBLE COMPLEX
//-------------------------------------------------------------------

force_inline
simd_compl<double, reg_128>::simd_compl(const impl_type& v) 
    : data(v) 
{};

force_inline
simd_compl<double, reg_128>::simd_compl(Integer re)
    : data(impl_type::set_lower(double(re)))
{};

force_inline
simd_compl<double, reg_128>::simd_compl(float re)
    : data(impl_type::set_lower(double(re)))
{};

force_inline
simd_compl<double, reg_128>::simd_compl(double re)
    : data(impl_type::set_lower(re))
{}

force_inline
simd_compl<double, reg_128>::simd_compl(double re, double im)
    : data(re, im)
{}

force_inline
simd_compl<double, reg_128>::simd_compl(const simd_double_complex& val)
{
    data = impl_type::load((const double*)&val, std::false_type());
}

force_inline
simd_compl<double, reg_128>::simd_compl(const simd_single_complex& val)
    : data(real(val), imag(val))
{}

force_inline
simd_compl<double, reg_128> simd_compl<double, reg_128>::zero()
{
    return impl_type::zero();
}

force_inline
simd_compl<double, reg_128> simd_compl<double, reg_128>::broadcast(const simd_double_complex* arr)
{
    return impl_type::load((const double*)arr, std::false_type());
}

force_inline
simd_compl<double, reg_128> simd_compl<double, reg_128>::broadcast(const simd_double_complex& arr)
{
    return impl_type::load((const double*)&arr, std::false_type());
}

force_inline
simd_compl<double, reg_128> simd_compl<double, reg_128>::broadcast(const double* arr)
{
    return simd_compl(arr[0]);
}

force_inline
simd_compl<double, reg_128> simd_compl<double, reg_128>::broadcast(const double& arr)
{
    return simd_compl(arr);
}

force_inline simd_compl<double,reg_128> 
simd_compl<double,reg_128>::load(const simd_double_complex* arr, std::true_type aligned)
{
    return impl_type::load((const double*)arr, aligned);
};

force_inline simd_compl<double,reg_128> 
simd_compl<double,reg_128>::load(const simd_double_complex* arr, std::false_type not_aligned)
{
    return impl_type::load((const double*)arr, not_aligned);
};

force_inline void 
simd_compl<double,reg_128>::store(simd_double_complex* arr, std::true_type aligned) const
{
    data.store((double*)arr, aligned);
};

force_inline void 
simd_compl<double,reg_128>::store(simd_double_complex* arr, std::false_type not_aligned) const
{
    data.store((double*)arr, not_aligned);
};

force_inline simd_compl<double,reg_128> 
simd_compl<double,reg_128>::load_reverse(const simd_double_complex* arr, std::true_type aligned)
{
    return load(arr, aligned);
}

force_inline simd_compl<double,reg_128> 
simd_compl<double,reg_128>::load_reverse(const simd_double_complex* arr, std::false_type not_aligned)
{
    return load(arr, not_aligned);
};

force_inline void 
simd_compl<double,reg_128>::store_reverse(simd_double_complex* arr, std::true_type aligned) const
{
    store(arr, aligned);
}

force_inline void 
simd_compl<double,reg_128>::store_reverse(simd_double_complex* arr, std::false_type not_aligned) const
{
    store(arr, not_aligned);
};

force_inline
simd_double_complex simd_compl<double, reg_128>::get(Integer pos) const
{
    (void)pos;
    return simd_double_complex(data.data.m128d_f64[0], data.data.m128d_f64[1]);
};

template<int Step>
force_inline
void simd_compl<double,reg_128>::scatter(simd_double_complex* arr) const
{
    return store(arr, std::false_type());
};

//-------------------------------------------------------------------
//                          SSE2 FLOAT COMPLEX
//-------------------------------------------------------------------
force_inline
simd_compl<float, reg_128>::simd_compl(const impl_type& v) 
    : data(v) 
{};

force_inline
simd_compl<float, reg_128>::simd_compl(float re)
    : data( re, 0.0f, re, 0.0f)
{};

force_inline
simd_compl<float, reg_128>::simd_compl(float re, float im)
    : data( re, im, re, im)
{};

force_inline
simd_compl<float, reg_128>::simd_compl(float re_0, float im_0, float re_1, float im_1)
    : data( re_0, im_0, re_1, im_1)
{};

force_inline
simd_compl<float, reg_128>::simd_compl(const simd_single_complex& v0, const simd_single_complex& v1)
    : data( real(v0), imag(v0), real(v1), imag(v1))
{};

force_inline
simd_compl<float, reg_128>::simd_compl(const simd_single_complex& val)
    : data( real(val), imag(val), real(val), imag(val))
{};

force_inline
simd_compl<float, reg_128> simd_compl<float, reg_128>::zero()
{
    return impl_type::zero();
}

force_inline
simd_compl<float, reg_128> simd_compl<float, reg_128>::broadcast(const simd_single_complex* arr)
{
    return simd_compl(arr[0]);
}

force_inline
simd_compl<float, reg_128> simd_compl<float, reg_128>::broadcast(const simd_single_complex& arr)
{
    return simd_compl(arr);
}

force_inline
simd_compl<float, reg_128> simd_compl<float, reg_128>::broadcast(const float* arr)
{
    return simd_compl(arr[0]);
}

force_inline
simd_compl<float, reg_128> simd_compl<float, reg_128>::broadcast(const float& arr)
{
    return simd_compl(arr);
}

force_inline simd_compl<float,reg_128> 
simd_compl<float,reg_128>::load(const simd_single_complex* arr, std::true_type aligned)
{
    return impl_type::load((const float*)arr, aligned);
};

force_inline simd_compl<float,reg_128> 
simd_compl<float,reg_128>::load(const simd_single_complex* arr, std::false_type not_aligned)
{
    return impl_type::load((const float*)arr, not_aligned);
};

force_inline void 
simd_compl<float,reg_128>::store(simd_single_complex* arr, std::false_type not_aligned) const
{
    data.store((float*)arr, not_aligned);
};

force_inline void 
simd_compl<float,reg_128>::store(simd_single_complex* arr, std::true_type aligned) const
{
    data.store((float*)arr, aligned);
};

force_inline simd_compl<float,reg_128> 
simd_compl<float,reg_128>::load_reverse(const simd_single_complex* arr, std::true_type aligned)
{
    static const int off            = vector_size - 1;
    simd_compl<float,reg_128> ret   = load(arr - off, aligned);
    return reverse(ret);
}

force_inline simd_compl<float,reg_128> 
simd_compl<float,reg_128>::load_reverse(const simd_single_complex* arr, std::false_type not_aligned)
{
    static const int off            = vector_size - 1;
    simd_compl<float,reg_128> ret   = load(arr - off, not_aligned);
    return reverse(ret);
};

force_inline void 
simd_compl<float,reg_128>::store_reverse(simd_single_complex* arr, std::true_type aligned) const
{
    static const int off            = vector_size - 1;
    simd_compl<float,reg_128> rev   = reverse(*this);
    rev.store(arr - off, aligned);
}

force_inline void 
simd_compl<float,reg_128>::store_reverse(simd_single_complex* arr, std::false_type not_aligned) const
{
    static const int off            = vector_size - 1;
    simd_compl<float,reg_128> rev   = reverse(*this);
    rev.store(arr - off, not_aligned);
};

force_inline
simd_single_complex simd_compl<float, reg_128>::get(int pos) const
{
    return simd_single_complex(data.get(pos*2), data.get(pos*2 + 1));
};

template<int Step>
force_inline
void simd_compl<float,reg_128>::scatter(simd_single_complex* arr0) const
{
    float* arr  = (float*)arr0;

    //no scatter intrinsic
    arr + 0*Step*2      = data.data.m128_f32[0];
    arr + 0*Step*2+1    = data.data.m128_f32[1];
    arr + 1*Step*2      = data.data.m128_f32[2];
    arr + 1*Step*2+1    = data.data.m128_f32[3];
};

}}
