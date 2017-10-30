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

#include "matcl-simd/complex/simd_128_compl.h"
#include "matcl-simd/arch/simd_impl.inl"

namespace matcl { namespace simd
{

//-------------------------------------------------------------------
//                          DOUBLE COMPLEX
//-------------------------------------------------------------------

template<class Simd_tag>
force_inline
simd_compl<double, 128, Simd_tag>::simd_compl(const impl_type& v) 
    : data(v) 
{};

template<class Simd_tag>
force_inline
simd_compl<double, 128, Simd_tag>::simd_compl(Integer re)
    : data(impl_type::set_lower(double(re)))
{};

template<class Simd_tag>
force_inline
simd_compl<double, 128, Simd_tag>::simd_compl(float re)
    : data(impl_type::set_lower(double(re)))
{};

template<class Simd_tag>
force_inline
simd_compl<double, 128, Simd_tag>::simd_compl(double re)
    : data(impl_type::set_lower(re))
{}

template<class Simd_tag>
force_inline
simd_compl<double, 128, Simd_tag>::simd_compl(double re, double im)
    : data(re, im)
{}

template<class Simd_tag>
force_inline
simd_compl<double, 128, Simd_tag>::simd_compl(const simd_double_complex& val)
{
    data = impl_type::load((const double*)&val, std::false_type());
}

template<class Simd_tag>
force_inline
simd_compl<double, 128, Simd_tag>::simd_compl(const simd_single_complex& val)
    : data(real(val), imag(val))
{}

template<class Simd_tag>
force_inline
simd_compl<double, 128, Simd_tag> simd_compl<double, 128, Simd_tag>::zero()
{
    return impl_type::zero();
}

template<class Simd_tag>
force_inline
simd_compl<double, 128, Simd_tag> simd_compl<double, 128, Simd_tag>::broadcast(const simd_double_complex* arr)
{
    return impl_type::load((const double*)arr, std::false_type());
}

template<class Simd_tag>
force_inline
simd_compl<double, 128, Simd_tag> simd_compl<double, 128, Simd_tag>::broadcast(const simd_double_complex& arr)
{
    return impl_type::load((const double*)&arr, std::false_type());
}

template<class Simd_tag>
force_inline
simd_compl<double, 128, Simd_tag> simd_compl<double, 128, Simd_tag>::broadcast(const double* arr)
{
    return simd_compl(arr[0]);
}

template<class Simd_tag>
force_inline
simd_compl<double, 128, Simd_tag> simd_compl<double, 128, Simd_tag>::broadcast(const double& arr)
{
    return simd_compl(arr);
}

template<class Simd_tag>
force_inline simd_compl<double, 128, Simd_tag> 
simd_compl<double, 128, Simd_tag>::load(const simd_double_complex* arr, std::true_type aligned)
{
    return impl_type::load((const double*)arr, aligned);
};

template<class Simd_tag>
force_inline simd_compl<double, 128, Simd_tag> 
simd_compl<double, 128, Simd_tag>::load(const simd_double_complex* arr, std::false_type not_aligned)
{
    return impl_type::load((const double*)arr, not_aligned);
};

template<class Simd_tag>
force_inline void 
simd_compl<double, 128, Simd_tag>::store(simd_double_complex* arr, std::true_type aligned) const
{
    data.store((double*)arr, aligned);
};

template<class Simd_tag>
force_inline void 
simd_compl<double, 128, Simd_tag>::store(simd_double_complex* arr, std::false_type not_aligned) const
{
    data.store((double*)arr, not_aligned);
};

template<class Simd_tag>
force_inline
simd_double_complex simd_compl<double, 128, Simd_tag>::get(Integer pos) const
{
    (void)pos;
    return simd_double_complex(data.get<0>(), data.get<1>());
};

template<class Simd_tag>
template<int Pos>
force_inline
simd_double_complex simd_compl<double, 128, Simd_tag>::get() const
{
    return simd_double_complex(data.get<0>(), data.get<1>());
};

template<class Simd_tag>
template<int Step>
force_inline
void simd_compl<double, 128, Simd_tag>::scatter(simd_double_complex* arr) const
{
    return store(arr, std::false_type());
};

//-------------------------------------------------------------------
//                          FLOAT COMPLEX
//-------------------------------------------------------------------
template<class Simd_tag>
force_inline
simd_compl<float, 128, Simd_tag>::simd_compl(const impl_type& v) 
    : data(v) 
{};

template<class Simd_tag>
force_inline
simd_compl<float, 128, Simd_tag>::simd_compl(float re)
    : data( re, 0.0f, re, 0.0f)
{};

template<class Simd_tag>
force_inline
simd_compl<float, 128, Simd_tag>::simd_compl(float re, float im)
    : data( re, im, re, im)
{};

template<class Simd_tag>
force_inline
simd_compl<float, 128, Simd_tag>::simd_compl(float re_0, float im_0, float re_1, float im_1)
    : data( re_0, im_0, re_1, im_1)
{};

template<class Simd_tag>
force_inline
simd_compl<float, 128, Simd_tag>::simd_compl(const simd_single_complex& v0, const simd_single_complex& v1)
    : data( real(v0), imag(v0), real(v1), imag(v1))
{};

template<class Simd_tag>
force_inline
simd_compl<float, 128, Simd_tag>::simd_compl(const simd_single_complex& val)
    : data( real(val), imag(val), real(val), imag(val))
{};

template<class Simd_tag>
force_inline
simd_compl<float, 128, Simd_tag> simd_compl<float, 128, Simd_tag>::zero()
{
    return impl_type::zero();
}

template<class Simd_tag>
force_inline
simd_compl<float, 128, Simd_tag> simd_compl<float, 128, Simd_tag>::broadcast(const simd_single_complex* arr)
{
    return simd_compl(arr[0]);
}

template<class Simd_tag>
force_inline
simd_compl<float, 128, Simd_tag> simd_compl<float, 128, Simd_tag>::broadcast(const simd_single_complex& arr)
{
    return simd_compl(arr);
}

template<class Simd_tag>
force_inline
simd_compl<float, 128, Simd_tag> simd_compl<float, 128, Simd_tag>::broadcast(const float* arr)
{
    return simd_compl(arr[0]);
}

template<class Simd_tag>
force_inline
simd_compl<float, 128, Simd_tag> simd_compl<float, 128, Simd_tag>::broadcast(const float& arr)
{
    return simd_compl(arr);
}

template<class Simd_tag>
force_inline simd_compl<float, 128, Simd_tag> 
simd_compl<float, 128, Simd_tag>::load(const simd_single_complex* arr, std::true_type aligned)
{
    return impl_type::load((const float*)arr, aligned);
};

template<class Simd_tag>
force_inline simd_compl<float, 128, Simd_tag> 
simd_compl<float, 128, Simd_tag>::load(const simd_single_complex* arr, std::false_type not_aligned)
{
    return impl_type::load((const float*)arr, not_aligned);
};

template<class Simd_tag>
force_inline void 
simd_compl<float, 128, Simd_tag>::store(simd_single_complex* arr, std::false_type not_aligned) const
{
    data.store((float*)arr, not_aligned);
};

template<class Simd_tag>
force_inline void 
simd_compl<float, 128, Simd_tag>::store(simd_single_complex* arr, std::true_type aligned) const
{
    data.store((float*)arr, aligned);
};

template<class Simd_tag>
force_inline
simd_single_complex simd_compl<float, 128, Simd_tag>::get(int pos) const
{
    return simd_single_complex(data.get(pos*2), data.get(pos*2 + 1));
};

template<class Simd_tag>
template<int Pos>
force_inline
simd_single_complex simd_compl<float, 128, Simd_tag>::get() const
{
    return simd_single_complex(data.get<Pos*2>(), data.get<Pos*2 + 1>());
};

template<class Simd_tag>
template<int Step>
force_inline
void simd_compl<float, 128, Simd_tag>::scatter(simd_single_complex* arr0) const
{
    float* arr  = (float*)arr0;

    //no scatter intrinsic
    arr[0*Step*2]      = data.get<0>();
    arr[0*Step*2+1]    = data.get<1>();
    arr[1*Step*2]      = data.get<2>();
    arr[1*Step*2+1]    = data.get<3>();
};

}}
