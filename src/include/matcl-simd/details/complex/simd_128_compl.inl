/*
 *  This file is a part of Matrix Computation Library (MATCL)
 *
 *  Copyright (c) Pawe³ Kowal 2017 - 2021
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
#include "matcl-simd/details/arch/simd_impl.inl"

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
simd_compl<double, 128, Simd_tag>::simd_compl(int32_t re)
    : data(impl_type::set_lower(double(re)))
{};

template<class Simd_tag>
force_inline
simd_compl<double, 128, Simd_tag>::simd_compl(int64_t re)
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
template<class Tag>
force_inline
simd_compl<double, 128, Simd_tag>::simd_compl(const simd_compl<double, 128, Tag>& s)
    : data(s.data)
{}

template<class Simd_tag>
force_inline
simd_compl<double, 128, Simd_tag> simd_compl<double, 128, Simd_tag>::zero()
{
    return impl_type::zero();
}

template<class Simd_tag>
force_inline simd_compl<double, 128, Simd_tag> 
simd_compl<double, 128, Simd_tag>::broadcast(const simd_double_complex* arr)
{
    return impl_type::load((const double*)arr, std::false_type());
}

template<class Simd_tag>
force_inline simd_compl<double, 128, Simd_tag> 
simd_compl<double, 128, Simd_tag>::broadcast(const simd_double_complex& arr)
{
    return impl_type::load((const double*)&arr, std::false_type());
}

template<class Simd_tag>
force_inline simd_compl<double, 128, Simd_tag> 
simd_compl<double, 128, Simd_tag>::broadcast(const double* arr)
{
    return simd_compl(arr[0]);
}

template<class Simd_tag>
force_inline simd_compl<double, 128, Simd_tag> 
simd_compl<double, 128, Simd_tag>::broadcast(const double& arr)
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
simd_double_complex simd_compl<double, 128, Simd_tag>::get(int pos) const
{
    (void)pos;
    const simd_double_complex* ptr = get_raw_ptr();
    return simd_double_complex(*ptr);
};

template<class Simd_tag>
force_inline
simd_double_complex simd_compl<double, 128, Simd_tag>::first() const
{
    return simd_double_complex(get_raw_ptr()[0]);
};

template<class Simd_tag>
force_inline
void simd_compl<double, 128, Simd_tag>::set(int pos, const simd_double_complex& val)
{
    simd_double_complex* ptr    = this->get_raw_ptr();
    ptr[pos]                    = val;
};

template<class Simd_tag>
force_inline
const simd_double_complex*
simd_compl<double, 128, Simd_tag>::get_raw_ptr() const
{
    return reinterpret_cast<const simd_double_complex*>(data.get_raw_ptr());
}

template<class Simd_tag>
force_inline
simd_double_complex*
simd_compl<double, 128, Simd_tag>::get_raw_ptr()
{
    return reinterpret_cast<simd_double_complex*>(data.get_raw_ptr());
}

template<class Simd_tag>
force_inline
simd_compl<float, 128, Simd_tag>
simd_compl<double, 128, Simd_tag>::convert_to_float() const
{
    return data.convert_to_float();
};

template<class Simd_tag>
template<int Step>
force_inline
void simd_compl<double, 128, Simd_tag>::scatter(simd_double_complex* arr) const
{
    return store(arr, std::false_type());
};

template<class Simd_tag>
force_inline simd_compl<double, 128, Simd_tag>& 
simd_compl<double, 128, Simd_tag>::operator+=(const simd_compl& x)
{
    *this = *this + x;
    return *this;
}

template<class Simd_tag>
force_inline simd_compl<double, 128, Simd_tag>& 
simd_compl<double, 128, Simd_tag>::operator-=(const simd_compl& x)
{
    *this = *this - x;
    return *this;
}

template<class Simd_tag>
force_inline simd_compl<double, 128, Simd_tag>& 
simd_compl<double, 128, Simd_tag>::operator*=(const simd_compl& x)
{
    *this = *this * x;
    return *this;
}

template<class Simd_tag>
force_inline simd_compl<double, 128, Simd_tag>& 
simd_compl<double, 128, Simd_tag>::operator/=(const simd_compl& x)
{
    *this = *this / x;
    return *this;
}

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
simd_compl<float, 128, Simd_tag>::simd_compl(const simd_compl& lo, const simd_compl& hi)
    :data(lo.data, hi.data)
{};

template<class Simd_tag>
force_inline
simd_compl<float, 128, Simd_tag>::simd_compl(const simd_single_complex& val)
    : data( real(val), imag(val), real(val), imag(val))
{};

template<class Simd_tag>
template<class Tag>
force_inline
simd_compl<float, 128, Simd_tag>::simd_compl(const simd_compl<float, 128, Tag>& s)
    : data(s.data)
{}

template<class Simd_tag>
force_inline
simd_compl<float, 128, Simd_tag> simd_compl<float, 128, Simd_tag>::zero()
{
    return impl_type::zero();
}

template<class Simd_tag>
force_inline simd_compl<float, 128, Simd_tag> 
simd_compl<float, 128, Simd_tag>::broadcast(const simd_single_complex* arr)
{
    return simd_compl(arr[0]);
}

template<class Simd_tag>
force_inline simd_compl<float, 128, Simd_tag>
simd_compl<float, 128, Simd_tag>::broadcast(const simd_single_complex& arr)
{
    return simd_compl(arr);
}

template<class Simd_tag>
force_inline simd_compl<float, 128, Simd_tag> 
simd_compl<float, 128, Simd_tag>::broadcast(const float* arr)
{
    return simd_compl(arr[0]);
}

template<class Simd_tag>
force_inline simd_compl<float, 128, Simd_tag>
simd_compl<float, 128, Simd_tag>::broadcast(const float& arr)
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
    const simd_single_complex* ptr  = this->get_raw_ptr();
    return simd_single_complex(ptr[pos]);
};

template<class Simd_tag>
force_inline
simd_single_complex simd_compl<float, 128, Simd_tag>::first() const
{
    return simd_single_complex(get_raw_ptr()[0]);
};

template<class Simd_tag>
force_inline
void simd_compl<float, 128, Simd_tag>::set(int pos, const simd_single_complex& val)
{
    simd_single_complex* ptr    = this->get_raw_ptr();
    ptr[pos]                    = val;
};

template<class Simd_tag>
force_inline
const simd_single_complex*
simd_compl<float, 128, Simd_tag>::get_raw_ptr() const
{
    return reinterpret_cast<const simd_single_complex*>(data.get_raw_ptr());
}

template<class Simd_tag>
force_inline
simd_single_complex*
simd_compl<float, 128, Simd_tag>::get_raw_ptr()
{
    return reinterpret_cast<simd_single_complex*>(data.get_raw_ptr());
}

template<class Simd_tag>
force_inline simd_compl<double, 128, Simd_tag>
simd_compl<float, 128, Simd_tag>::convert_low_to_double() const
{
    return data.convert_low_to_double();
}

template<class Simd_tag>
force_inline simd_compl<double, 128, Simd_tag>
simd_compl<float, 128, Simd_tag>::convert_high_to_double() const
{
    return data.convert_high_to_double();
}

template<class Simd_tag>
force_inline typename simd_compl<float, 128, Simd_tag>::simd_double_2
simd_compl<float, 128, Simd_tag>::convert_to_double() const
{
    return data.convert_to_double();
}

template<class Simd_tag>
template<int Step>
force_inline
void simd_compl<float, 128, Simd_tag>::scatter(simd_single_complex* arr0) const
{
    float* arr          = (float*)arr0;
    const float* ptr    = data.get_raw_ptr();

    //no scatter intrinsic
    arr[0*Step*2]      =ptr[0];
    arr[0*Step*2+1]    =ptr[1];
    arr[1*Step*2]      =ptr[2];
    arr[1*Step*2+1]    =ptr[3];
};

template<class Simd_tag>
force_inline simd_compl<float, 128, Simd_tag>& 
simd_compl<float, 128, Simd_tag>::operator+=(const simd_compl& x)
{
    *this = *this + x;
    return *this;
}

template<class Simd_tag>
force_inline simd_compl<float, 128, Simd_tag>& 
simd_compl<float, 128, Simd_tag>::operator-=(const simd_compl& x)
{
    *this = *this - x;
    return *this;
}

template<class Simd_tag>
force_inline simd_compl<float, 128, Simd_tag>& 
simd_compl<float, 128, Simd_tag>::operator*=(const simd_compl& x)
{
    *this = *this * x;
    return *this;
}

template<class Simd_tag>
force_inline simd_compl<float, 128, Simd_tag>& 
simd_compl<float, 128, Simd_tag>::operator/=(const simd_compl& x)
{
    *this = *this / x;
    return *this;
}

}}
