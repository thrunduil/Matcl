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

#include "matcl-simd/arch/sse/simd_256.h"

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
simd<double, 256, sse_tag>::simd(const double& val)    
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
double simd<double, 256, sse_tag>::get(int pos) const
{ 
    return data[0].data.m128d_f64[pos] ; 
};

template<int Pos>
force_inline
double simd<double, 256, sse_tag>::get() const
{ 
    return data[0].data.m128d_f64[Pos]; 
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

force_inline void simd<double, 256, sse_tag>::store(double* arr, std::true_type aligned) const
{
    data[0].store(arr, aligned);
    data[1].store(arr + 2, aligned);
};

force_inline void simd<double, 256, sse_tag>::store(double* arr, std::false_type not_aligned) const
{
    data[0].store(arr, not_aligned);
    data[1].store(arr + 2, not_aligned);
};

template<int Step>
force_inline
void simd<double, 256, sse_tag>::scatter(double* arr) const
{
    //no scatter intrinsic
    arr[0*Step] = get<0>();
    arr[1*Step] = get<1>();
    arr[2*Step] = get<2>();
    arr[3*Step] = get<3>();
};

//-------------------------------------------------------------------
//                          SSE2 SINGLE
//-------------------------------------------------------------------
force_inline
simd<float, 256, sse_tag>::simd(const float& val) 
{
    data[0] = simd_half(val); 
    data[1] = simd_half(val);
}

force_inline
simd<float, 256, sse_tag>::simd(const simd_half& lo_hi)
{
    data[0] = lo_hi; 
    data[1] = lo_hi;
}

force_inline
simd<float, 256, sse_tag>::simd(float v1, float v2, float v3, float v4, float v5, float v6, float v7, float v8) 
{
    data[0] = simd_half(v1, v2, v3, v4);
    data[1] = simd_half(v5, v6, v7, v8);
}

force_inline
simd<float, 256, sse_tag>::simd(const impl_type& v)
{
    data[0] = v[0];
    data[1] = v[1];
};

force_inline
simd<float, 256, sse_tag>::simd(const simd_half& lo, const simd_half& hi)
{
    data[0] = lo;
    data[1] = hi;
}

force_inline
float simd<float, 256, sse_tag>::get(int pos) const  
{ 
    return data[0].data.m128_f32[pos] ; 
};

template<int Pos>
force_inline
float simd<float, 256, sse_tag>::get() const
{ 
    return data[0].data.m128_f32[Pos]; 
};

force_inline simd<float, 256, sse_tag>::simd_half
simd<float, 256, sse_tag>::extract_low() const
{
    return data[0];
}

force_inline simd<float, 256, sse_tag>::simd_half
simd<float, 256, sse_tag>::extract_high() const
{
    return data[1];
}

force_inline
simd<float, 256, sse_tag> simd<float, 256, sse_tag>::zero()
{
    simd<float, 256, sse_tag> ret;
    ret.data[0] = simd_half::zero();
    ret.data[1] = simd_half::zero();
    return ret;
}

force_inline simd<float, 256, sse_tag> 
simd<float, 256, sse_tag>::load(const float* arr, std::true_type aligned)
{
    simd<float, 256, sse_tag> ret;
    ret.data[0] = simd_half::load(arr, aligned);
    ret.data[1] = simd_half::load(arr + 4, aligned);
    return ret;

};

force_inline simd<float, 256, sse_tag> 
simd<float, 256, sse_tag>::load(const float* arr, std::false_type not_aligned)
{
    simd<float, 256, sse_tag> ret;
    ret.data[0] = simd_half::load(arr, not_aligned);
    ret.data[1] = simd_half::load(arr + 4, not_aligned);
    return ret;
};

force_inline simd<float, 256, sse_tag> 
simd<float, 256, sse_tag>::broadcast(const float* arr)
{
    simd<float, 256, sse_tag> ret;
    ret.data[0] = simd_half::broadcast(arr);
    ret.data[1] = simd_half::broadcast(arr);
    return ret;
};

force_inline simd<float, 256, sse_tag> 
simd<float, 256, sse_tag>::broadcast(const float& arr)
{
    simd<float, 256, sse_tag> ret;
    ret.data[0] = simd_half::broadcast(arr);
    ret.data[1] = simd_half::broadcast(arr);
    return ret;
};

force_inline void simd<float, 256, sse_tag>::store(float* arr, std::true_type aligned) const
{
    data[0].store(arr, aligned);
    data[1].store(arr + 4, aligned);
};

force_inline void simd<float, 256, sse_tag>::store(float* arr, std::false_type not_aligned) const
{
    data[0].store(arr, not_aligned);
    data[1].store(arr + 4, not_aligned);
};

template<int Step>
force_inline
void simd<float, 256, sse_tag>::scatter(float* arr) const
{
    //no scatter intrinsic
    arr[0*Step] = get<0>();
    arr[1*Step] = get<1>();
    arr[2*Step] = get<2>();
    arr[3*Step] = get<3>();
    arr[4*Step] = get<4>();
    arr[5*Step] = get<5>();
    arr[6*Step] = get<6>();
    arr[7*Step] = get<7>();
};

}}
