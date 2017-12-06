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

#include "matcl-simd/arch/sse/simd_double_256.h"

namespace matcl { namespace simd
{

//-------------------------------------------------------------------
//                          SSE2 DOUBLE
//-------------------------------------------------------------------

force_inline
simd<double, 256, sse_tag>::simd(Integer val)
    : simd(double(val))
{}

force_inline
simd<double, 256, sse_tag>::simd(float val)
    : simd(double(val))
{}

force_inline
simd<double, 256, sse_tag>::simd(double val)    
{
    data[0] = simd_half(val); 
    data[1] = simd_half(val);
}

force_inline
simd<double, 256, sse_tag>::simd(double v1, double v2, double v3, double v4) 
{
    data[0] = simd_half(v1, v2); 
    data[1] = simd_half(v3, v4);
}

force_inline
simd<double, 256, sse_tag>::simd(const simd_half& lo, const simd_half& hi)
{
    data[0] = lo; 
    data[1] = hi;
}

force_inline
simd<double, 256, sse_tag>::simd(const simd_half& lo_hi)
{
    data[0] = lo_hi; 
    data[1] = lo_hi;
}

force_inline
simd<double, 256, sse_tag>::simd(const impl_type& v)
{
    data[0] = v[0];
    data[1] = v[1];
};

force_inline
simd<double, 256, sse_tag>::simd(const simd<double, 256, nosimd_tag>& s)
{
    std::true_type aligned;

    data[0] = simd_half::load(s.data + 0, aligned);
    data[1] = simd_half::load(s.data + 1, aligned);
}

force_inline
simd<double, 256, sse_tag>::simd(const simd<double, 256, avx_tag>& s)
{
    data[0] = s.extract_low();
    data[1] = s.extract_high();
}

force_inline
double simd<double, 256, sse_tag>::get(int pos) const
{ 
    return get_raw_ptr()[pos] ; 
};

force_inline
double simd<double, 256, sse_tag>::first() const
{ 
    return data[0].first();
};

force_inline
void simd<double, 256, sse_tag>::set(int pos, double val)
{ 
    get_raw_ptr()[pos] = val; 
};

force_inline
const double* simd<double, 256, sse_tag>::get_raw_ptr() const
{ 
    return reinterpret_cast<const double*>(&data); 
};

force_inline
double* simd<double, 256, sse_tag>::get_raw_ptr()
{ 
    return reinterpret_cast<double*>(&data); 
};

force_inline simd<double, 256, sse_tag>::simd_half
simd<double, 256, sse_tag>::extract_low() const
{
    return data[0];
}

force_inline simd<double, 256, sse_tag>::simd_half
simd<double, 256, sse_tag>::extract_high() const
{
    return data[1];
}

force_inline
simd<double, 256, sse_tag> simd<double, 256, sse_tag>::zero()
{
    simd<double, 256, sse_tag> ret;
    ret.data[0] = simd_half::zero();
    ret.data[1] = simd_half::zero();
    return ret;
}

force_inline
simd<double, 256, sse_tag> simd<double, 256, sse_tag>::minus_zero()
{
    simd<double, 256, sse_tag> ret;
    ret.data[0] = simd_half::minus_zero();
    ret.data[1] = simd_half::minus_zero();
    return ret;
}

force_inline
simd<double, 256, sse_tag> simd<double, 256, sse_tag>::one()
{
    simd<double, 256, sse_tag> ret;
    ret.data[0] = simd_half::one();
    ret.data[1] = simd_half::one();
    return ret;
}

force_inline
simd<double, 256, sse_tag> simd<double, 256, sse_tag>::minus_one()
{
    simd<double, 256, sse_tag> ret;
    ret.data[0] = simd_half::minus_one();
    ret.data[1] = simd_half::minus_one();
    return ret;
}

force_inline simd<double, 256, sse_tag> 
simd<double, 256, sse_tag>::load(const double* arr, std::true_type aligned)
{
    simd<double, 256, sse_tag> ret;
    ret.data[0] = simd_half::load(arr, aligned);
    ret.data[1] = simd_half::load(arr + 2, aligned);
    return ret;
};

force_inline simd<double, 256, sse_tag> 
simd<double, 256, sse_tag>::load(const double* arr, std::false_type not_aligned)
{
    simd<double, 256, sse_tag> ret;
    ret.data[0] = simd_half::load(arr, not_aligned);
    ret.data[1] = simd_half::load(arr + 2, not_aligned);
    return ret;
};

force_inline simd<double, 256, sse_tag> 
simd<double, 256, sse_tag>::gather(const double* arr, const simd_128_int32& ind)
{
    simd<double, 256, sse_tag> ret;
    ret.data[0] = simd_half::gather(arr, ind);
    ret.data[1] = simd_half::gather(arr, ind.extract_high());
    return ret;
}

force_inline simd<double, 256, sse_tag> 
simd<double, 256, sse_tag>::gather(const double* arr, const simd_256_int64& ind)
{
    simd<double, 128, sse_tag> ret_lo, ret_hi;
    ret_lo  = simd_half::gather(arr, ind.extract_low());
    ret_hi  = simd_half::gather(arr, ind.extract_high());

    simd<double, 128, sse_tag> ret_1(ret_lo, ret_hi);
    return simd(ret_1, ret_1);
}

force_inline simd<double, 256, sse_tag> 
simd<double, 256, sse_tag>::broadcast(const double* arr)
{
    simd<double, 256, sse_tag> ret;
    ret.data[0] = simd_half::broadcast(arr);
    ret.data[1] = simd_half::broadcast(arr);
    return ret;
};

force_inline simd<double, 256, sse_tag> 
simd<double, 256, sse_tag>::broadcast(const double& arr)
{
    simd<double, 256, sse_tag> ret;
    ret.data[0] = simd_half::broadcast(arr);
    ret.data[1] = simd_half::broadcast(arr);
    return ret;
};

force_inline void 
simd<double, 256, sse_tag>::store(double* arr, std::true_type aligned) const
{
    data[0].store(arr, aligned);
    data[1].store(arr + 2, aligned);
};

force_inline void 
simd<double, 256, sse_tag>::store(double* arr, std::false_type not_aligned) const
{
    data[0].store(arr, not_aligned);
    data[1].store(arr + 2, not_aligned);
};

force_inline simd<float, 128, sse_tag>
simd<double, 256, sse_tag>::convert_to_float() const
{
    return simd_128_float(data[0].convert_to_float(), data[1].convert_to_float());
}

force_inline simd<int32_t, 128, sse_tag>
simd<double, 256, sse_tag>::convert_to_int32() const
{
    return simd_128_int32(data[0].convert_to_int32(), data[1].convert_to_int32());
}

force_inline simd<float, 256, sse_tag>
simd<double, 256, sse_tag>::reinterpret_as_float() const
{
    return *reinterpret_cast<const simd_256_float*>(this);
}

force_inline simd<int32_t, 256, sse_tag>
simd<double, 256, sse_tag>::reinterpret_as_int32() const
{
    return *reinterpret_cast<const simd_256_int32*>(this);
}

force_inline simd<int64_t, 256, sse_tag>
simd<double, 256, sse_tag>::reinterpret_as_int64() const
{
    return *reinterpret_cast<const simd_256_int64*>(this);
}

template<int Step>
force_inline
void simd<double, 256, sse_tag>::scatter(double* arr) const
{
    const double* this_data  = get_raw_ptr();

    //no scatter intrinsic
    arr[0*Step] = this_data[0];
    arr[1*Step] = this_data[1];
    arr[2*Step] = this_data[2];
    arr[3*Step] = this_data[3];
};

force_inline simd<double, 256, sse_tag>& 
simd<double, 256, sse_tag>::operator+=(const simd& x)
{
    *this = *this + x;
    return *this;
}

force_inline simd<double, 256, sse_tag>& 
simd<double, 256, sse_tag>::operator-=(const simd& x)
{
    *this = *this - x;
    return *this;
}

force_inline simd<double, 256, sse_tag>& 
simd<double, 256, sse_tag>::operator*=(const simd& x)
{
    *this = *this * x;
    return *this;
}

force_inline simd<double, 256, sse_tag>& 
simd<double, 256, sse_tag>::operator/=(const simd& x)
{
    *this = *this / x;
    return *this;
}

}}
