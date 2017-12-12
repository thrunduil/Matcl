/*
 *  This file is a part of Matrix Computation Library (MATCL)
 *
 *  Copyright (c) Pawe� Kowal 2017
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

#include "matcl-simd/arch/nosimd/simd_int32_128.h"

namespace matcl { namespace simd
{

//-------------------------------------------------------------------
//                          GENERIC int32_t
//-------------------------------------------------------------------

force_inline
simd<int32_t, 128, nosimd_tag>::simd(int32_t val)
{
    data[0] = val;
    data[1] = val;
    data[2] = val;
    data[3] = val;
}

force_inline
simd<int32_t, 128, nosimd_tag>::simd(int32_t v1, int32_t v2, int32_t v3, int32_t v4) 
{
    data[0] = v1;
    data[1] = v2;
    data[2] = v3;
    data[3] = v4;
}

force_inline
simd<int32_t, 128, nosimd_tag>::simd(const simd& lo, const simd& hi)
{
    data[0] = lo.data[0];
    data[1] = lo.data[1];
    data[2] = hi.data[0];
    data[3] = hi.data[1];
}

force_inline
simd<int32_t, 128, nosimd_tag>::simd(const simd<int32_t, 128, sse_tag>& s)
{
    s.store(data, std::true_type());
}

force_inline
simd<int32_t, 128, nosimd_tag>::simd(const simd<int32_t, 128, scalar_sse_tag>& s)
    :simd(s.first())
{};

force_inline
simd<int32_t, 128, nosimd_tag>::simd(const simd<int32_t, 128, scalar_nosimd_tag>& s)
    :simd(s.first())
{}

force_inline simd<int32_t, 128, nosimd_tag> 
simd<int32_t, 128, nosimd_tag>::broadcast(const int32_t* arr)
{
    simd<int32_t, 128, nosimd_tag> ret;

    ret.data[0] = arr[0];
    ret.data[1] = arr[0];
    ret.data[2] = arr[0];
    ret.data[3] = arr[0];

    return ret;
};

force_inline simd<int32_t, 128, nosimd_tag> 
simd<int32_t, 128, nosimd_tag>::broadcast(const int32_t& arr)
{ 
    simd<int32_t, 128, nosimd_tag> ret;

    ret.data[0] = arr;
    ret.data[1] = arr;
    ret.data[2] = arr;
    ret.data[3] = arr;

    return ret;
};

force_inline
int32_t simd<int32_t, 128, nosimd_tag>::get(int pos) const
{ 
    return data[pos] ; 
};

force_inline
int32_t simd<int32_t, 128, nosimd_tag>::first() const
{ 
    return data[0]; 
};

force_inline
void simd<int32_t, 128, nosimd_tag>::set(int pos, int32_t val)
{ 
    data[pos] = val; 
};

force_inline
const int32_t* simd<int32_t, 128, nosimd_tag>::get_raw_ptr() const
{ 
    return data; 
};

force_inline
int32_t* simd<int32_t, 128, nosimd_tag>::get_raw_ptr()
{ 
    return data; 
};

force_inline
simd<int32_t, 128, nosimd_tag> simd<int32_t, 128, nosimd_tag>::zero()
{
    simd<int32_t, 128, nosimd_tag> ret;

    ret.data[0] = 0;
    ret.data[1] = 0;
    ret.data[2] = 0;
    ret.data[3] = 0;

    return ret;
}

force_inline
simd<int32_t, 128, nosimd_tag> simd<int32_t, 128, nosimd_tag>::one()
{
    simd<int32_t, 128, nosimd_tag> ret;

    ret.data[0] = 1;
    ret.data[1] = 1;
    ret.data[2] = 1;
    ret.data[3] = 1;

    return ret;
}

force_inline
simd<int32_t, 128, nosimd_tag> simd<int32_t, 128, nosimd_tag>::minus_one()
{
    simd<int32_t, 128, nosimd_tag> ret;

    ret.data[0] = -1;
    ret.data[1] = -1;
    ret.data[2] = -1;
    ret.data[3] = -1;

    return ret;
}

force_inline simd<int32_t, 128, nosimd_tag> 
simd<int32_t, 128, nosimd_tag>::load(const int32_t* arr, std::true_type aligned)
{
    (void)aligned;
    simd<int32_t, 128, nosimd_tag> ret;

    ret.data[0] = arr[0];
    ret.data[1] = arr[1];
    ret.data[2] = arr[2];
    ret.data[3] = arr[3];

    return ret;
};

force_inline simd<int32_t, 128, nosimd_tag> 
simd<int32_t, 128, nosimd_tag>::load(const int32_t* arr, std::false_type not_aligned)
{
    (void)not_aligned;
    simd<int32_t, 128, nosimd_tag> ret;

    ret.data[0] = arr[0];
    ret.data[1] = arr[1];
    ret.data[2] = arr[2];
    ret.data[3] = arr[3];

    return ret;
};

force_inline simd<int32_t, 128, nosimd_tag>
simd<int32_t, 128, nosimd_tag>::gather(const int32_t* arr, const simd_int32& ind)
{
    simd<int32_t, 128, nosimd_tag> ret;

    ret.data[0] = arr[ind.data[0]];
    ret.data[1] = arr[ind.data[1]];
    ret.data[2] = arr[ind.data[2]];
    ret.data[3] = arr[ind.data[3]];

    return ret;
}

force_inline simd<int32_t, 128, nosimd_tag> 
simd<int32_t, 128, nosimd_tag>::gather(const int32_t* arr, const simd_int64& ind)
{
    simd<int32_t, 128, nosimd_tag> ret;

    ret.data[0] = arr[ind.data[0]];
    ret.data[1] = arr[ind.data[1]];

    ret.data[2] = 0;
    ret.data[3] = 0;

    return ret;
}

force_inline simd<int32_t, 128, nosimd_tag> 
simd<int32_t, 128, nosimd_tag>::set_lower(int32_t v)
{
    simd<int32_t, 128, nosimd_tag> ret;

    ret.data[0] = v;
    ret.data[1] = 0;
    ret.data[2] = 0;
    ret.data[3] = 0;

    return ret;
};

force_inline simd<int64_t, 128, nosimd_tag>
simd<int32_t, 128, nosimd_tag>::convert_low_to_int64() const
{
    simd<int64_t, 128, nosimd_tag> ret;

    ret.data[0] = scalar_func::convert_int32_int64(data[0]);
    ret.data[1] = scalar_func::convert_int32_int64(data[1]);

    return ret;
};

force_inline simd<int64_t, 128, nosimd_tag>
simd<int32_t, 128, nosimd_tag>::convert_high_to_int64() const
{
    simd<int64_t, 128, nosimd_tag> ret;

    ret.data[0] = scalar_func::convert_int32_int64(data[2]);
    ret.data[1] = scalar_func::convert_int32_int64(data[3]);

    return ret;
};

force_inline
simd<int64_t, 256, nosimd_tag>
simd<int32_t, 128, nosimd_tag>::convert_to_int64() const
{
    simd<int64_t, 256, nosimd_tag> ret;

    ret.data[0] = scalar_func::convert_int32_int64(data[0]);
    ret.data[1] = scalar_func::convert_int32_int64(data[1]);
    ret.data[2] = scalar_func::convert_int32_int64(data[2]);
    ret.data[3] = scalar_func::convert_int32_int64(data[3]);

    return ret;
};

force_inline simd<double, 128, nosimd_tag>
simd<int32_t, 128, nosimd_tag>::convert_low_to_double() const
{
    simd<double, 128, nosimd_tag> ret;

    ret.data[0] = scalar_func::convert_int32_double(data[0]);
    ret.data[1] = scalar_func::convert_int32_double(data[1]);

    return ret;
};

force_inline simd<double, 128, nosimd_tag>
simd<int32_t, 128, nosimd_tag>::convert_high_to_double() const
{
    simd<double, 128, nosimd_tag> ret;

    ret.data[0] = scalar_func::convert_int32_double(data[2]);
    ret.data[1] = scalar_func::convert_int32_double(data[3]);

    return ret;
};

force_inline
simd<double, 256, nosimd_tag>
simd<int32_t, 128, nosimd_tag>::convert_to_double() const
{
    simd<double, 256, nosimd_tag> ret;

    ret.data[0] = scalar_func::convert_int32_double(data[0]);
    ret.data[1] = scalar_func::convert_int32_double(data[1]);
    ret.data[2] = scalar_func::convert_int32_double(data[2]);
    ret.data[3] = scalar_func::convert_int32_double(data[3]);

    return ret;
};

force_inline
simd<float, 128, nosimd_tag>
simd<int32_t, 128, nosimd_tag>::convert_to_float() const
{
    simd<float, 128, nosimd_tag> ret;

    ret.data[0] = scalar_func::convert_int32_float(data[0]);
    ret.data[1] = scalar_func::convert_int32_float(data[1]);
    ret.data[2] = scalar_func::convert_int32_float(data[2]);
    ret.data[3] = scalar_func::convert_int32_float(data[3]);

    return ret;
};

force_inline simd<float, 128, nosimd_tag>
simd<int32_t, 128, nosimd_tag>::reinterpret_as_float() const
{
    return *reinterpret_cast<const simd_float*>(this);
}

force_inline simd<int64_t, 128, nosimd_tag>
simd<int32_t, 128, nosimd_tag>::reinterpret_as_int64() const
{
    return *reinterpret_cast<const simd_int64*>(this);
}

force_inline simd<double, 128, nosimd_tag>
simd<int32_t, 128, nosimd_tag>::reinterpret_as_double() const
{
    return *reinterpret_cast<const simd_double*>(this);
}

force_inline simd<int32_t, 128, nosimd_tag> 
simd<int32_t, 128, nosimd_tag>::extract_low() const
{
    return *this;
}

force_inline simd<int32_t, 128, nosimd_tag> 
simd<int32_t, 128, nosimd_tag>::extract_high() const
{
    simd<int32_t, 128, nosimd_tag> ret;
    ret.data[0] = data[2];
    ret.data[1] = data[3];
    ret.data[2] = data[2];
    ret.data[3] = data[3];

    return ret;
}

force_inline void 
simd<int32_t, 128, nosimd_tag>::store(int32_t* arr, std::true_type aligned) const
{
    (void)aligned;
    arr[0] = data[0];
    arr[1] = data[1];
    arr[2] = data[2];
    arr[3] = data[3];
};

force_inline void 
simd<int32_t, 128, nosimd_tag>::store(int32_t* arr, std::false_type not_aligned) const
{
    (void)not_aligned;
    arr[0] = data[0];
    arr[1] = data[1];
    arr[2] = data[2];
    arr[3] = data[3];
};

template<int Step>
force_inline
void simd<int32_t, 128, nosimd_tag>::scatter(int32_t* arr) const
{
    //no scatter intrinsic
    arr[0*Step] = data[0];
    arr[1*Step] = data[1];
    arr[2*Step] = data[2];
    arr[3*Step] = data[3];
};

force_inline simd<int32_t, 128, nosimd_tag>& 
simd<int32_t, 128, nosimd_tag>::operator+=(const simd& x)
{
    *this = *this + x;
    return *this;
}

force_inline simd<int32_t, 128, nosimd_tag>& 
simd<int32_t, 128, nosimd_tag>::operator-=(const simd& x)
{
    *this = *this - x;
    return *this;
}

force_inline simd<int32_t, 128, nosimd_tag>& 
simd<int32_t, 128, nosimd_tag>::operator*=(const simd& x)
{
    *this = *this * x;
    return *this;
}

}}
