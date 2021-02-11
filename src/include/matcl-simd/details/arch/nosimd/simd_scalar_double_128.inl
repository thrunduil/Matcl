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

#include "matcl-simd/arch/nosimd/simd_double_128.h"

namespace matcl { namespace simd
{

//-------------------------------------------------------------------
//                          GENERIC DOUBLE
//-------------------------------------------------------------------

force_inline
simd<double, 128, scalar_nosimd_tag>::simd(int32_t val)
    : data(double(val)) 
{}

force_inline
simd<double, 128, scalar_nosimd_tag>::simd(int64_t val)
    : data(double(val)) 
{}

force_inline
simd<double, 128, scalar_nosimd_tag>::simd(float val)
    : data(double(val)) 
{}

force_inline
simd<double, 128, scalar_nosimd_tag>::simd(double val)
{
    data = val;
}

#if MATCL_ARCHITECTURE_HAS_SSE2
    force_inline
    simd<double, 128, scalar_nosimd_tag>::simd(const simd<double, 128, sse_tag>& s)
        :data(s.first())
    {}

    force_inline
    simd<double, 128, scalar_nosimd_tag>::simd(const simd<double, 128, scalar_sse_tag>& s)
        :data(s.first())
    {}
#endif

force_inline
simd<double, 128, scalar_nosimd_tag>::simd(const simd<double, 128, nosimd_tag>& s)
    :data(s.first())
{}

force_inline
simd<double, 128, scalar_nosimd_tag> simd<double, 128, scalar_nosimd_tag>::broadcast(const double* arr)
{ 
    return simd(arr[0]);
};

force_inline
simd<double, 128, scalar_nosimd_tag> simd<double, 128, scalar_nosimd_tag>::broadcast(const double& arr)
{ 
    return simd(arr);
};

force_inline
simd<double, 128, scalar_nosimd_tag> simd<double, 128, scalar_nosimd_tag>::set_lower(double v)
{
    return simd(v);
};

force_inline
simd<float, 128, scalar_nosimd_tag>
simd<double, 128, scalar_nosimd_tag>::convert_to_float() const
{
    return simd<float, 128, scalar_nosimd_tag>(scalar_func::convert_double_float(data));
};

force_inline
simd<int32_t, 128, scalar_nosimd_tag>
simd<double, 128, scalar_nosimd_tag>::convert_to_int32() const
{
    return simd<int32_t, 128, scalar_nosimd_tag>(scalar_func::convert_double_int32(data));
};

force_inline
simd<int64_t, 128, scalar_nosimd_tag>
simd<double, 128, scalar_nosimd_tag>::convert_to_int64() const
{
    return simd<int64_t, 128, scalar_nosimd_tag>(scalar_func::convert_double_int64(data));
};

force_inline simd<int32_t, 128, scalar_nosimd_tag>
simd<double, 128, scalar_nosimd_tag>::reinterpret_as_int32() const
{
    return *reinterpret_cast<const simd_int32*>(this);
}

force_inline simd<int64_t, 128, scalar_nosimd_tag>
simd<double, 128, scalar_nosimd_tag>::reinterpret_as_int64() const
{
    return *reinterpret_cast<const simd_int64*>(this);
}

force_inline simd<float, 128, scalar_nosimd_tag>
simd<double, 128, scalar_nosimd_tag>::reinterpret_as_float() const
{
    return *reinterpret_cast<const simd_float*>(this);
}

force_inline simd<double, 128, nosimd_tag>
simd<double, 128, scalar_nosimd_tag>::as_vector() const
{
    return simd<double, 128, nosimd_tag>(data);
}

force_inline
double simd<double, 128, scalar_nosimd_tag>::get(int pos) const
{ 
    (void)pos;
    return data; 
};

force_inline
double simd<double, 128, scalar_nosimd_tag>::first() const
{ 
    return data; 
};

force_inline
void simd<double, 128, scalar_nosimd_tag>::set(int pos, double val)
{ 
    (void)pos;
    data = val; 
};

force_inline
const double* simd<double, 128, scalar_nosimd_tag>::get_raw_ptr() const
{ 
    return &data; 
};

force_inline
double* simd<double, 128, scalar_nosimd_tag>::get_raw_ptr()
{ 
    return &data; 
};

force_inline
simd<double, 128, scalar_nosimd_tag> simd<double, 128, scalar_nosimd_tag>::zero()
{
    return simd(0.0);
}

force_inline
simd<double, 128, scalar_nosimd_tag> simd<double, 128, scalar_nosimd_tag>::minus_zero()
{
    return simd(-0.0);
}

force_inline
simd<double, 128, scalar_nosimd_tag> simd<double, 128, scalar_nosimd_tag>::one()
{
    return simd(1.0);
}

force_inline
simd<double, 128, scalar_nosimd_tag> simd<double, 128, scalar_nosimd_tag>::minus_one()
{
    return simd(-1.0);
}

force_inline simd<double, 128, scalar_nosimd_tag> 
simd<double, 128, scalar_nosimd_tag>::load(const double* arr, std::true_type aligned)
{
    (void)aligned;
    return simd(arr[0]);
};

force_inline simd<double, 128, scalar_nosimd_tag> 
simd<double, 128, scalar_nosimd_tag>::load(const double* arr, std::false_type not_aligned)
{
    (void)not_aligned;
    return simd(arr[0]);
};

force_inline simd<double, 128, scalar_nosimd_tag> 
simd<double, 128, scalar_nosimd_tag>::gather(const double* arr, const simd_int32& ind)
{
    return simd(arr[ind.first()]);
}

force_inline simd<double, 128, scalar_nosimd_tag> 
simd<double, 128, scalar_nosimd_tag>::gather(const double* arr, const simd_int64& ind)
{
    return simd(arr[ind.first()]);
}

force_inline void 
simd<double, 128, scalar_nosimd_tag>::store(double* arr, std::true_type aligned) const
{
    (void)aligned;
    arr[0] = data;
};

force_inline void 
simd<double, 128, scalar_nosimd_tag>::store(double* arr, std::false_type not_aligned) const
{
    (void)not_aligned;
    arr[0] = data;
};

template<int Step>
force_inline
void simd<double, 128, scalar_nosimd_tag>::scatter(double* arr) const
{
    arr[0*Step] = data;
};

force_inline simd<double, 128, scalar_nosimd_tag>& 
simd<double, 128, scalar_nosimd_tag>::operator+=(const simd& x)
{
    *this = *this + x;
    return *this;
}

force_inline simd<double, 128, scalar_nosimd_tag>& 
simd<double, 128, scalar_nosimd_tag>::operator-=(const simd& x)
{
    *this = *this - x;
    return *this;
}

force_inline simd<double, 128, scalar_nosimd_tag>& 
simd<double, 128, scalar_nosimd_tag>::operator*=(const simd& x)
{
    *this = *this * x;
    return *this;
}

force_inline simd<double, 128, scalar_nosimd_tag>& 
simd<double, 128, scalar_nosimd_tag>::operator/=(const simd& x)
{
    *this = *this / x;
    return *this;
}

}}
