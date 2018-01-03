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

#include "matcl-simd/simd.h"
#include "matcl-core/details/scalfunc_real.h"

namespace matcl { namespace simd { namespace details
{

constexpr int eval_log2(int Val)
{
    return (Val < 2) ? 0 : 1 + eval_log2(Val/2);
};

constexpr int eval_pow2(int Val)
{
    return (Val <= 0)? 1 : 2 * eval_pow2(Val - 1);
};

template<class Arg_type, class Coef_type>
struct eval_fma
{
    force_inline
    static Arg_type eval(const Arg_type& x, const Arg_type& y, const Coef_type& z)
    {
        return fma_f(x, y, Arg_type(z));
        //return x * y + Arg_type(z);
    }
};

template<>
struct eval_fma<float, float>
{
    force_inline
    static float eval(float x, float y, float z)
    {
        namespace mrds = matcl::raw::details::scal_func;
        return mrds::fma_f(x, y, z);
    }
};

template<>
struct eval_fma<double, double>
{
    force_inline
    static double eval(double x, double y, double z)
    {
        namespace mrds = matcl::raw::details::scal_func;
        return mrds::fma_f(x, y, z);
    }
};

template<class T, class TZ>
force_inline
T fma(const T& x, const T& y, const TZ& z)
{
    return eval_fma<T, TZ>::eval(x, y, z);
};

template<class T>
struct eval_abs
{
    force_inline
    static T eval(const T& x)
    {
        return abs(x);
    }
};

template<>
struct eval_abs<float>
{
    force_inline
    static float eval(float x)
    {
        namespace mrds = matcl::raw::details::scal_func;
        return mrds::abs(x);
    }
};

template<>
struct eval_abs<double>
{
    force_inline
    static double eval(double x)
    {
        namespace mrds = matcl::raw::details::scal_func;
        return mrds::abs(x);
    }
};

struct trans_id
{
    template<class T>
    force_inline
    static T eval(const T& arg)
    {
        return arg;
    }
};

struct trans_abs
{
    template<class T>
    force_inline
    static T eval(const T& arg)
    {
        return eval_abs<T>::eval(arg);
    }
};

template<class T>
struct eval_eps
{
    force_inline
    static T eval()
    {
        return std::numeric_limits<T>::epsilon();
    }
};

template<class T, int Bits, class Tag>
struct eval_eps<matcl::simd::simd<T, Bits, Tag>>
{
    using simd_type = matcl::simd::simd<T, Bits, Tag>;

    force_inline
    static simd_type eval()
    {
        T e = eval_eps<T>::eval();
        return simd_type(e);
    }
};

template<class T>
struct eval_lt
{
    force_inline
    static bool eval(T x, T y)
    {
        return x < y;
    }
};

template<class T, int Bits, class Tag>
struct eval_lt<matcl::simd::simd<T, Bits, Tag>>
{
    using simd_type = matcl::simd::simd<T, Bits, Tag>;

    static bool eval(const simd_type& x, const simd_type& y)
    {
        simd_type res = matcl::simd::lt(x, y);
        return all(res);
    };
};

template<class Arg>
struct broadcast
{
    template<class Coef>
    force_inline
    static Arg eval(const Coef* arr)
    {
        return Arg::broadcast(arr);
    }

    template<class Coef>
    force_inline
    static Arg eval(const Coef& arr)
    {
        return Arg(arr);
    }

    force_inline
    static Arg eval(const Arg* arr)
    {
        return arr[0];
    }

    force_inline
    static Arg eval(const Arg& arr)
    {
        return arr;
    }
};

template<>
struct broadcast<double>
{
    template<class Coef>
    force_inline
    static double eval(const Coef* arr)
    {
        return double(arr[0]);
    }

    template<class Coef>
    force_inline
    static double eval(const Coef& arr)
    {
        return double(arr);
    }

    force_inline
    static double eval(const double& arr)
    {
        return arr;
    }

    force_inline
    static double eval(const double* arr)
    {
        return arr[0];
    }
};

template<>
struct broadcast<float>
{
    template<class Coef>
    force_inline
    static float eval(const Coef* arr)
    {
        return float(arr[0]);
    }

    template<class Coef>
    force_inline
    static float eval(const Coef& arr)
    {
        return float(arr);
    }

    force_inline
    static float eval(const float& arr)
    {
        return arr;
    }

    force_inline
    static float eval(const float* arr)
    {
        return arr[0];
    }
};

}}}

