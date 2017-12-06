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

#include "matcl-simd/arch/nosimd/simd_int64_128.h"

namespace matcl { namespace simd
{

//-------------------------------------------------------------------
//                          GENERIC INT64_T
//-------------------------------------------------------------------

force_inline
simd<int64_t, 128, nosimd_tag>::simd(int32_t val)
    : simd(int64_t(val)) 
{}

force_inline
simd<int64_t, 128, nosimd_tag>::simd(int64_t val)
{
    data[0] = val;
    data[1] = val;
}

force_inline
simd<int64_t, 128, nosimd_tag>::simd(int64_t v1, int64_t v2)
{
    data[0] = v1;
    data[1] = v2;
}

force_inline
simd<int64_t, 128, nosimd_tag>::simd(const simd& lo, const simd& hi)
{
    data[0] = lo.data[0];
    data[1] = hi.data[0];
}

force_inline
simd<int64_t, 128, nosimd_tag>::simd(const impl_type& v) 
{
    data[0] = v[0];
    data[1] = v[1];
};

force_inline
simd<int64_t, 128, nosimd_tag>::simd(const simd<int64_t, 128, sse_tag>& s)
{
    s.store(data, std::true_type());
}

force_inline simd<int64_t, 128, nosimd_tag> 
simd<int64_t, 128, nosimd_tag>::broadcast(const int64_t* arr)
{ 
    simd<int64_t, 128, nosimd_tag> ret;
    ret.data[0] = arr[0];
    ret.data[1] = arr[0];

    return ret;
};

force_inline simd<int64_t, 128, nosimd_tag>
simd<int64_t, 128, nosimd_tag>::broadcast(const int64_t& arr)
{ 
    simd<int64_t, 128, nosimd_tag> ret;
    ret.data[0] = arr;
    ret.data[1] = arr;

    return ret;
};

force_inline
simd<int64_t, 128, nosimd_tag> simd<int64_t, 128, nosimd_tag>::set_lower(int64_t v)
{
    simd<int64_t, 128, nosimd_tag> ret;
    ret.data[0] = v;
    ret.data[1] = 0;

    return ret;
};

force_inline
simd<int32_t, 128, nosimd_tag>
simd<int64_t, 128, nosimd_tag>::convert_to_int32() const
{
    simd<int32_t, 128, nosimd_tag> ret;
    ret.data[0] = (int32_t)data[0];
    ret.data[1] = (int32_t)data[1];
    ret.data[2] = 0;
    ret.data[3] = 0;

    return ret;
};

force_inline simd<float, 128, nosimd_tag>
simd<int64_t, 128, nosimd_tag>::reinterpret_as_float() const
{
    return *reinterpret_cast<const simd_128_float*>(this);
}

force_inline simd<int32_t, 128, nosimd_tag>
simd<int64_t, 128, nosimd_tag>::reinterpret_as_int32() const
{
    return *reinterpret_cast<const simd_128_int32*>(this);
}

force_inline simd<double, 128, nosimd_tag>
simd<int64_t, 128, nosimd_tag>::reinterpret_as_double() const
{
    return *reinterpret_cast<const simd_128_double*>(this);
}

force_inline simd<int64_t, 128, nosimd_tag> 
simd<int64_t, 128, nosimd_tag>::extract_low() const
{
    return *this;
}

force_inline simd<int64_t, 128, nosimd_tag> 
simd<int64_t, 128, nosimd_tag>::extract_high() const
{
    simd<int64_t, 128, nosimd_tag> ret;
    ret.data[0] = data[1];
    ret.data[1] = data[1];

    return ret;
}

force_inline
int64_t simd<int64_t, 128, nosimd_tag>::get(int pos) const
{ 
    return data[pos] ; 
};

force_inline
int64_t simd<int64_t, 128, nosimd_tag>::first() const
{ 
    return data[0]; 
};

template<int Pos>
force_inline
int64_t simd<int64_t, 128, nosimd_tag>::get() const
{ 
    return data[Pos] ; 
};

force_inline
void simd<int64_t, 128, nosimd_tag>::set(int pos, int64_t val)
{ 
    data[pos] = val; 
};

template<int Pos>
force_inline
void simd<int64_t, 128, nosimd_tag>::set(int64_t val)
{ 
    data[Pos] = val; 
};

force_inline
const int64_t* simd<int64_t, 128, nosimd_tag>::get_raw_ptr() const
{ 
    return data; 
};

force_inline
int64_t* simd<int64_t, 128, nosimd_tag>::get_raw_ptr()
{ 
    return data; 
};

force_inline
simd<int64_t, 128, nosimd_tag> simd<int64_t, 128, nosimd_tag>::zero()
{
    simd<int64_t, 128, nosimd_tag> ret;
    ret.data[0] = 0;
    ret.data[1] = 0;

    return ret;
}

force_inline
simd<int64_t, 128, nosimd_tag> simd<int64_t, 128, nosimd_tag>::one()
{
    simd<int64_t, 128, nosimd_tag> ret;
    ret.data[0] = 1;
    ret.data[1] = 1;

    return ret;
}

force_inline
simd<int64_t, 128, nosimd_tag> simd<int64_t, 128, nosimd_tag>::minus_one()
{
    simd<int64_t, 128, nosimd_tag> ret;
    ret.data[0] = -1;
    ret.data[1] = -1;

    return ret;
}

force_inline simd<int64_t, 128, nosimd_tag> 
simd<int64_t, 128, nosimd_tag>::load(const int64_t* arr, std::true_type aligned)
{
    (void)aligned;
    simd<int64_t, 128, nosimd_tag> ret;
    ret.data[0] = arr[0];
    ret.data[1] = arr[1];

    return ret;
};

force_inline simd<int64_t, 128, nosimd_tag> 
simd<int64_t, 128, nosimd_tag>::load(const int64_t* arr, std::false_type not_aligned)
{
    (void)not_aligned;
    simd<int64_t, 128, nosimd_tag> ret;
    ret.data[0] = arr[0];
    ret.data[1] = arr[1];

    return ret;
};

force_inline simd<int64_t, 128, nosimd_tag> 
simd<int64_t, 128, nosimd_tag>::gather(const int64_t* arr, const simd_128_int32& ind)
{
    simd<int64_t, 128, nosimd_tag> ret;
    ret.data[0] = arr[ind.data[0]];
    ret.data[1] = arr[ind.data[1]];

    return ret;
}

force_inline simd<int64_t, 128, nosimd_tag> 
simd<int64_t, 128, nosimd_tag>::gather(const int64_t* arr, const simd_128_int64& ind)
{
    simd<int64_t, 128, nosimd_tag> ret;
    ret.data[0] = arr[ind.data[0]];
    ret.data[1] = arr[ind.data[1]];

    return ret;
}

force_inline void 
simd<int64_t, 128, nosimd_tag>::store(int64_t* arr, std::true_type aligned) const
{
    (void)aligned;
    arr[0] = data[0];
    arr[1] = data[1];
};

force_inline void 
simd<int64_t, 128, nosimd_tag>::store(int64_t* arr, std::false_type not_aligned) const
{
    (void)not_aligned;
    arr[0] = data[0];
    arr[1] = data[1];
};

template<int Step>
force_inline
void simd<int64_t, 128, nosimd_tag>::scatter(int64_t* arr) const
{
    arr[0*Step] = data[0];
    arr[1*Step] = data[1];
};

}}
