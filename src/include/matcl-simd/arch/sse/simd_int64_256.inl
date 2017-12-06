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

#include "matcl-simd/arch/sse/simd_int64_256.h"

namespace matcl { namespace simd
{

//-------------------------------------------------------------------
//                          SSE2 INT64_T
//-------------------------------------------------------------------

force_inline
simd<int64_t, 256, sse_tag>::simd(int32_t val)
    : simd(int64_t(val))
{}

force_inline
simd<int64_t, 256, sse_tag>::simd(int64_t val)    
{
    data[0] = simd_half(val); 
    data[1] = simd_half(val);
}

force_inline
simd<int64_t, 256, sse_tag>::simd(int64_t v1, int64_t v2, int64_t v3, int64_t v4) 
{
    data[0] = simd_half(v1, v2); 
    data[1] = simd_half(v3, v4);
}

force_inline
simd<int64_t, 256, sse_tag>::simd(const simd_half& lo, const simd_half& hi)
{
    data[0] = lo; 
    data[1] = hi;
}

force_inline
simd<int64_t, 256, sse_tag>::simd(const simd_half& lo_hi)
{
    data[0] = lo_hi; 
    data[1] = lo_hi;
}

force_inline
simd<int64_t, 256, sse_tag>::simd(const impl_type& v)
{
    data[0] = v[0];
    data[1] = v[1];
};

force_inline
simd<int64_t, 256, sse_tag>::simd(const simd<int64_t, 256, nosimd_tag>& s)
{
    std::true_type aligned;

    data[0] = simd_half::load(s.data + 0, aligned);
    data[1] = simd_half::load(s.data + 1, aligned);
}

force_inline
simd<int64_t, 256, sse_tag>::simd(const simd<int64_t, 256, avx_tag>& s)
{
    data[0] = s.extract_low();
    data[1] = s.extract_high();
}

force_inline
int64_t simd<int64_t, 256, sse_tag>::get(int pos) const
{ 
    return get_raw_ptr()[pos] ; 
};

force_inline
int64_t simd<int64_t, 256, sse_tag>::first() const
{ 
    return data[0].first();
};

force_inline
void simd<int64_t, 256, sse_tag>::set(int pos, int64_t val)
{ 
    get_raw_ptr()[pos] = val; 
};

force_inline
const int64_t* simd<int64_t, 256, sse_tag>::get_raw_ptr() const
{ 
    return reinterpret_cast<const int64_t*>(&data); 
};

force_inline
int64_t* simd<int64_t, 256, sse_tag>::get_raw_ptr()
{ 
    return reinterpret_cast<int64_t*>(&data); 
};

force_inline simd<int64_t, 256, sse_tag>::simd_half
simd<int64_t, 256, sse_tag>::extract_low() const
{
    return data[0];
}

force_inline simd<int64_t, 256, sse_tag>::simd_half
simd<int64_t, 256, sse_tag>::extract_high() const
{
    return data[1];
}

force_inline
simd<int64_t, 256, sse_tag> simd<int64_t, 256, sse_tag>::zero()
{
    simd<int64_t, 256, sse_tag> ret;
    ret.data[0] = simd_half::zero();
    ret.data[1] = simd_half::zero();
    return ret;
}

force_inline
simd<int64_t, 256, sse_tag> simd<int64_t, 256, sse_tag>::one()
{
    simd<int64_t, 256, sse_tag> ret;
    ret.data[0] = simd_half::one();
    ret.data[1] = simd_half::one();
    return ret;
}

force_inline
simd<int64_t, 256, sse_tag> simd<int64_t, 256, sse_tag>::minus_one()
{
    simd<int64_t, 256, sse_tag> ret;
    ret.data[0] = simd_half::minus_one();
    ret.data[1] = simd_half::minus_one();
    return ret;
}

force_inline simd<int64_t, 256, sse_tag> 
simd<int64_t, 256, sse_tag>::load(const int64_t* arr, std::true_type aligned)
{
    simd<int64_t, 256, sse_tag> ret;
    ret.data[0] = simd_half::load(arr, aligned);
    ret.data[1] = simd_half::load(arr + 2, aligned);
    return ret;
};

force_inline simd<int64_t, 256, sse_tag> 
simd<int64_t, 256, sse_tag>::load(const int64_t* arr, std::false_type not_aligned)
{
    simd<int64_t, 256, sse_tag> ret;
    ret.data[0] = simd_half::load(arr, not_aligned);
    ret.data[1] = simd_half::load(arr + 2, not_aligned);
    return ret;
};

force_inline simd<int64_t, 256, sse_tag> 
simd<int64_t, 256, sse_tag>::gather(const int64_t* arr, const simd_128_int32& ind)
{
    simd<int64_t, 256, sse_tag> ret;
    ret.data[0] = simd_half::gather(arr, ind);
    ret.data[1] = simd_half::gather(arr, ind.extract_high());
    return ret;
}

force_inline simd<int64_t, 256, sse_tag> 
simd<int64_t, 256, sse_tag>::gather(const int64_t* arr, const simd_256_int64& ind)
{
    simd<int64_t, 128, sse_tag> ret_lo, ret_hi;
    ret_lo  = simd_half::gather(arr, ind.extract_low());
    ret_hi  = simd_half::gather(arr, ind.extract_high());

    simd<int64_t, 128, sse_tag> ret_1(ret_lo, ret_hi);
    return simd(ret_1, ret_1);
}

force_inline simd<int64_t, 256, sse_tag> 
simd<int64_t, 256, sse_tag>::broadcast(const int64_t* arr)
{
    simd<int64_t, 256, sse_tag> ret;
    ret.data[0] = simd_half::broadcast(arr);
    ret.data[1] = simd_half::broadcast(arr);
    return ret;
};

force_inline simd<int64_t, 256, sse_tag> 
simd<int64_t, 256, sse_tag>::broadcast(const int64_t& arr)
{
    simd<int64_t, 256, sse_tag> ret;
    ret.data[0] = simd_half::broadcast(arr);
    ret.data[1] = simd_half::broadcast(arr);
    return ret;
};

force_inline void 
simd<int64_t, 256, sse_tag>::store(int64_t* arr, std::true_type aligned) const
{
    data[0].store(arr, aligned);
    data[1].store(arr + 2, aligned);
};

force_inline void 
simd<int64_t, 256, sse_tag>::store(int64_t* arr, std::false_type not_aligned) const
{
    data[0].store(arr, not_aligned);
    data[1].store(arr + 2, not_aligned);
};

force_inline simd<int32_t, 128, sse_tag>
simd<int64_t, 256, sse_tag>::convert_to_int32() const
{
    return simd_128_int32(data[0].convert_to_int32(), data[1].convert_to_int32());
}

force_inline simd<double, 256, sse_tag>
simd<int64_t, 256, sse_tag>::reinterpret_as_double() const
{
    return *reinterpret_cast<const simd_256_double*>(this);
}

force_inline simd<float, 256, sse_tag>
simd<int64_t, 256, sse_tag>::reinterpret_as_float() const
{
    return *reinterpret_cast<const simd_256_float*>(this);
}

force_inline simd<int32_t, 256, sse_tag>
simd<int64_t, 256, sse_tag>::reinterpret_as_int32() const
{
    return *reinterpret_cast<const simd_256_int32*>(this);
}

template<int Step>
force_inline
void simd<int64_t, 256, sse_tag>::scatter(int64_t* arr) const
{
    const int64_t* this_data  = get_raw_ptr();

    //no scatter intrinsic
    arr[0*Step] = this_data[0];
    arr[1*Step] = this_data[1];
    arr[2*Step] = this_data[2];
    arr[3*Step] = this_data[3];
};

}}
