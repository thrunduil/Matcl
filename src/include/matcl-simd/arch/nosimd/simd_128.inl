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

#include "matcl-simd/arch/nosimd/simd_128.h"

namespace matcl { namespace simd
{

//-------------------------------------------------------------------
//                          GENERIC DOUBLE
//-------------------------------------------------------------------

force_inline
simd<double, 128, nosimd_tag>::simd(Integer val)
    : simd(double(val)) 
{}

force_inline
simd<double, 128, nosimd_tag>::simd(float val)
    : simd(double(val)) 
{}

force_inline
simd<double, 128, nosimd_tag>::simd(double val)
{
    data[0] = val;
    data[1] = val;
}

force_inline
simd<double, 128, nosimd_tag>::simd(double v1, double v2)
{
    data[0] = v1;
    data[1] = v2;
}

force_inline
simd<double, 128, nosimd_tag>::simd(const impl_type& v) 
{
    data[0] = v[0];
    data[1] = v[1];
};

force_inline
simd<double, 128, nosimd_tag> simd<double, 128, nosimd_tag>::broadcast(const double* arr)
{ 
    simd<double, 128, nosimd_tag> ret;
    ret.data[0] = arr[0];
    ret.data[1] = arr[0];

    return ret;
};

force_inline
simd<double, 128, nosimd_tag> simd<double, 128, nosimd_tag>::broadcast(const double& arr)
{ 
    simd<double, 128, nosimd_tag> ret;
    ret.data[0] = arr;
    ret.data[1] = arr;

    return ret;
};

force_inline
simd<double, 128, nosimd_tag> simd<double, 128, nosimd_tag>::set_lower(double v)
{
    simd<double, 128, nosimd_tag> ret;
    ret.data[0] = v;
    ret.data[1] = 0.0;

    return ret;
};

force_inline
double simd<double, 128, nosimd_tag>::get(int pos) const
{ 
    return data[pos] ; 
};

template<int Pos>
force_inline
double simd<double, 128, nosimd_tag>::get() const
{ 
    return data[Pos] ; 
};

force_inline
simd<double, 128, nosimd_tag> simd<double, 128, nosimd_tag>::zero()
{
    simd<double, 128, nosimd_tag> ret;
    ret.data[0] = 0.0;
    ret.data[1] = 0.0;

    return ret;
}

force_inline simd<double, 128, nosimd_tag> 
simd<double, 128, nosimd_tag>::load(const double* arr, std::true_type aligned)
{
    (void)aligned;
    simd<double, 128, nosimd_tag> ret;
    ret.data[0] = arr[0];
    ret.data[1] = arr[1];

    return ret;
};

force_inline simd<double, 128, nosimd_tag> 
simd<double, 128, nosimd_tag>::load(const double* arr, std::false_type not_aligned)
{
    (void)not_aligned;
    simd<double, 128, nosimd_tag> ret;
    ret.data[0] = arr[0];
    ret.data[1] = arr[1];

    return ret;
};

force_inline void simd<double, 128, nosimd_tag>::store(double* arr, std::true_type aligned) const
{
    (void)aligned;
    arr[0] = data[0];
    arr[1] = data[1];
};

force_inline void simd<double, 128, nosimd_tag>::store(double* arr, std::false_type not_aligned) const
{
    (void)not_aligned;
    arr[0] = data[0];
    arr[1] = data[1];
};

template<int Step>
force_inline
void simd<double, 128, nosimd_tag>::scatter(double* arr) const
{
    arr[0*Step] = data[0];
    arr[1*Step] = data[1];
};

//-------------------------------------------------------------------
//                          GENERIC SINGLE
//-------------------------------------------------------------------
force_inline
simd<float, 128, nosimd_tag>::simd(float val)
{
    data[0] = val;
    data[1] = val;
    data[2] = val;
    data[3] = val;
}

force_inline
simd<float, 128, nosimd_tag>::simd(float v1, float v2, float v3, float v4) 
{
    data[0] = v1;
    data[1] = v2;
    data[2] = v3;
    data[3] = v4;
}

force_inline
simd<float, 128, nosimd_tag>::simd(const impl_type& v)
{
    data[0] = v[0];
    data[1] = v[1];
    data[2] = v[2];
    data[3] = v[3];
};

force_inline
simd<float, 128, nosimd_tag> simd<float, 128, nosimd_tag>::broadcast(const float* arr)
{
    simd<float, 128, nosimd_tag> ret;

    ret.data[0] = arr[0];
    ret.data[1] = arr[0];
    ret.data[2] = arr[0];
    ret.data[3] = arr[0];

    return ret;
};

force_inline
simd<float, 128, nosimd_tag> simd<float, 128, nosimd_tag>::broadcast(const float& arr)
{ 
    simd<float, 128, nosimd_tag> ret;

    ret.data[0] = arr;
    ret.data[1] = arr;
    ret.data[2] = arr;
    ret.data[3] = arr;

    return ret;
};

force_inline
float simd<float, 128, nosimd_tag>::get(int pos) const
{ 
    return data[pos] ; 
};

template<int Pos>
force_inline
float simd<float, 128, nosimd_tag>::get() const
{ 
    return data[Pos]; 
};

force_inline
simd<float, 128, nosimd_tag> simd<float, 128, nosimd_tag>::zero()
{
    simd<float, 128, nosimd_tag> ret;

    ret.data[0] = 0.0f;
    ret.data[1] = 0.0f;
    ret.data[2] = 0.0f;
    ret.data[3] = 0.0f;

    return ret;
}

force_inline simd<float, 128, nosimd_tag> 
simd<float, 128, nosimd_tag>::load(const float* arr, std::true_type aligned)
{
    (void)aligned;
    simd<float, 128, nosimd_tag> ret;

    ret.data[0] = arr[0];
    ret.data[1] = arr[1];
    ret.data[2] = arr[2];
    ret.data[3] = arr[3];

    return ret;
};

force_inline simd<float, 128, nosimd_tag> 
simd<float, 128, nosimd_tag>::load(const float* arr, std::false_type not_aligned)
{
    (void)not_aligned;
    simd<float, 128, nosimd_tag> ret;

    ret.data[0] = arr[0];
    ret.data[1] = arr[1];
    ret.data[2] = arr[2];
    ret.data[3] = arr[3];

    return ret;
};

force_inline simd<float, 128, nosimd_tag> 
simd<float, 128, nosimd_tag>::set_lower(float v)
{
    simd<float, 128, nosimd_tag> ret;

    ret.data[0] = v;
    ret.data[1] = 0.0f;
    ret.data[2] = 0.0f;
    ret.data[3] = 0.0f;

    return ret;
};

force_inline void simd<float, 128, nosimd_tag>::store(float* arr, std::true_type aligned) const
{
    (void)aligned;
    arr[0] = data[0];
    arr[1] = data[1];
    arr[2] = data[2];
    arr[3] = data[3];
};


force_inline void simd<float, 128, nosimd_tag>::store(float* arr, std::false_type not_aligned) const
{
    (void)not_aligned;
    arr[0] = data[0];
    arr[1] = data[1];
    arr[2] = data[2];
    arr[3] = data[3];
};

template<int Step>
force_inline
void simd<float, 128, nosimd_tag>::scatter(float* arr) const
{
    //no scatter intrinsic
    arr[0*Step] = data[0];
    arr[1*Step] = data[1];
    arr[2*Step] = data[2];
    arr[3*Step] = data[3];
};

}}
