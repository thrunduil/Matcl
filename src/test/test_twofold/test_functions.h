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

#include "matcl-core/config.h"
#include "matcl-mp/matcl_mp.h"
#include "matcl-simd/simd.h"

#include <iomanip>

namespace test_functions
{

namespace ms = matcl::simd;

using matcl::mp_float;
using matcl::twofold;

inline 
double get_value(double x)
{
    return x;
};

inline 
float get_value(float x)
{
    return x;
};

inline 
mp_float get_value(const mp_float& x)
{
    return x;
};

template<class Float>
Float get_value(const twofold<Float>& x)
{
    return x.value;
};

template<class TB, int Bits, class Tag>
matcl::simd::simd<TB, Bits, Tag> 
get_value(const matcl::simd::simd<TB, Bits, Tag>& x)
{
    return x;
};


// missing functions
inline mp_float twofold_sum(const mp_float& a, const mp_float& b)          { return a + b; };
inline mp_float twofold_sum_sorted(const mp_float& a, const mp_float& b)   { return a + b; };
inline mp_float twofold_minus(const mp_float& a, const mp_float& b)        { return a - b; };
inline mp_float twofold_mult(const mp_float& a, const mp_float& b)         { return a * b; };
inline mp_float twofold_mult_dekker(const mp_float& a, const mp_float& b)  { return a * b; };
inline mp_float twofold_div(const mp_float& a, const mp_float& b)          { return a / b; };
inline mp_float twofold_sqrt(const mp_float& a)                            { return sqrt(a); };

template<class T> struct get_prec_sum{};
template<> struct get_prec_sum<double>
{
    static size_t eval()
    {
        return matcl::precision::precision_double();
    }
};

template<> struct get_prec_sum<float>
{
    static size_t eval()
    {
        return matcl::precision::precision_float();
    }
};

template<class T>
mp_float twofold_sum(const mp_float& a)
{
    matcl::precision prec = matcl::precision(get_prec_sum<T>::eval());
    return mp_float(a, prec);
};

inline
mp_float twofold_minus_sorted(const mp_float& a, const mp_float& b)
{
    return a - b; 
};

// double
inline double twofold_sum(double a, double b)          { return a + b; };
inline double twofold_sum_sorted(double a, double b)   { return a + b; };
inline double twofold_minus(double a, double b)        { return a - b; };
inline double twofold_mult(double a, double b)         { return a * b; };
inline double twofold_mult_dekker(double a, double b)  { return a * b; };
inline double twofold_div(double a, double b)          { return a / b; };
inline double twofold_sqrt(double a)                   { return matcl::sqrt(a); };

template<class T>
double twofold_sum(double a)                    { return a; };

//simd
template<class TB, int Bits, class Tag>
matcl::simd::simd<TB, Bits, Tag>
twofold_sum(const matcl::simd::simd<TB, Bits, Tag>& a)
{ 
    return a; 
};

template<class T, int Bits, class Tag>
matcl::simd::simd<T, Bits, Tag>
twofold_sqrt(const matcl::simd::simd<T, Bits, Tag>& a)
{ 
    return matcl::simd::sqrt(a); 
};

template<class T, int Bits, class Tag>
matcl::simd::simd<T, Bits, Tag>
twofold_mult(const matcl::simd::simd<T, Bits, Tag>& a, const matcl::simd::simd<T, Bits, Tag>& b)
{ 
    return a * b; 
};

template<class T, int Bits, class Tag>
matcl::simd::simd<T, Bits, Tag>
twofold_mult_dekker(const matcl::simd::simd<T, Bits, Tag>& a, const matcl::simd::simd<T, Bits, Tag>& b)
{ 
    return a * b; 
};

template<class T, int Bits, class Tag>
matcl::simd::simd<T, Bits, Tag>
twofold_sum(const matcl::simd::simd<T, Bits, Tag>& a, const matcl::simd::simd<T, Bits, Tag>& b)
{ 
    return a + b; 
};

template<class T, int Bits, class Tag>
matcl::simd::simd<T, Bits, Tag>
twofold_sum_sorted(const matcl::simd::simd<T, Bits, Tag>& a, const matcl::simd::simd<T, Bits, Tag>& b)
{ 
    return a + b; 
};

template<class T, int Bits, class Tag>
matcl::simd::simd<T, Bits, Tag>
twofold_minus(const matcl::simd::simd<T, Bits, Tag>& a, const matcl::simd::simd<T, Bits, Tag>& b)
{ 
    return a - b; 
};

template<class T, int Bits, class Tag>
matcl::simd::simd<T, Bits, Tag>
twofold_minus_sorted(const matcl::simd::simd<T, Bits, Tag>& a, const matcl::simd::simd<T, Bits, Tag>& b)
{ 
    return a - b; 
};

template<class T, int Bits, class Tag>
matcl::simd::simd<T, Bits, Tag>
twofold_div(const matcl::simd::simd<T, Bits, Tag>& a, const matcl::simd::simd<T, Bits, Tag>& b)
{ 
    return a / b; 
};

inline
double twofold_minus_sorted(double a, double b)
{
    return a - b; 
};

//float
inline float twofold_sum(float a, float b)         { return a + b; };
inline float twofold_sum_sorted(float a, float b)  { return a + b; };
inline float twofold_minus(float a, float b)        { return a - b; };
inline float twofold_mult(float a, float b)         { return a * b; };
inline float twofold_mult_dekker(float a, float b)  { return a * b; };
inline float twofold_div(float a, float b)          { return a / b; };
inline float twofold_sqrt(float a)                  { return matcl::sqrt(a); };

template<class T>
float twofold_sum(float a)                   { return a; };

inline float twofold_minus_sorted(float a, float b)
{
    return a - b; 
};

template<class Float>
force_inline
twofold<Float> twofold_sum(const twofold<Float>& a, const twofold<Float>& b)
{ 
    return matcl::twofold_sum(a.value, b.value); 
};

template<class Float>
force_inline
twofold<Float> twofold_sum_sorted(const twofold<Float>& a, const twofold<Float>& b)
{ 
    return matcl::twofold_sum_sorted(a.value, b.value); 
};

template<class T, int Bits, class Tag>
force_inline
twofold<matcl::simd::simd<T, Bits, Tag>> 
twofold_sum_sorted(const twofold<matcl::simd::simd<T, Bits, Tag>>& a, 
                     const twofold<matcl::simd::simd<T, Bits, Tag>>& b)
{ 
    return matcl::twofold_sum_sorted(a.value, b.value); 
};

template<class Float>
force_inline
twofold<Float> twofold_minus(const twofold<Float>& a, const twofold<Float>& b)
{ 
    return matcl::twofold_minus(a.value, b.value); 
};

template<class Float>
force_inline
twofold<Float> twofold_minus_sorted(const twofold<Float>& a, const twofold<Float>& b)
{ 
    return matcl::twofold_minus_sorted(a.value, b.value); 
};

template<class T, int Bits, class Tag>
force_inline
twofold<matcl::simd::simd<T, Bits, Tag>> 
twofold_minus_sorted(const twofold<matcl::simd::simd<T, Bits, Tag>>& a, 
                     const twofold<matcl::simd::simd<T, Bits, Tag>>& b)
{ 
    return matcl::twofold_minus_sorted(a.value, b.value); 
};

template<class Float>
force_inline
twofold<Float> twofold_mult(const twofold<Float>& a, const twofold<Float>& b)
{ 
    return matcl::twofold_mult(a.value, b.value); 
};

template<class Float>
force_inline
twofold<Float> twofold_mult_dekker(const twofold<Float>& a, const twofold<Float>& b)
{ 
    twofold<Float> ret = twofold<Float>(twofold<Float>::uninitialized());

    matcl::twofold_mult_dekker(a.value, b.value, ret.value, ret.error); 
    return ret;
};

template<class Float>
force_inline
twofold<Float> twofold_div(const twofold<Float>& a, const twofold<Float>& b)
{ 
    return matcl::twofold_div(a.value, b.value); 
};

template<class Float>
force_inline
twofold<Float> twofold_sqrt(const twofold<Float>& a)
{ 
    return matcl::twofold_sqrt(a.value); 
};

template<class T, class Float>
force_inline
twofold<Float> twofold_sum(const twofold<Float>& a)
{ 
    return twofold<Float>(a.sum());
};

struct Func_sqrt_1
{
    template<class T>    
    static T eval(const T& x1)
    {
        return twofold_sqrt(x1);
    }

    template<class T2>
    static T2 convert_arg_1(const T2& x)
    {
        return T2(x.value);
    };

    static std::string name()
    { 
        return "sqrt_1"; 
    };
};

struct Func_abs
{
    template<class T>    
    static T eval(const T& x1)
    {
        return abs(x1);
    }

    template<class T2>
    static T2 convert_arg_1(const T2& x)
    {
        return x;
    };

    static std::string name()
    { 
        return "abs"; 
    };
};

template<class T_base>
struct Func_sum
{
    template<class T>    
    static T eval(const T& x1)
    {
        return twofold_sum<T_base>(x1);
    }

    template<class T2>
    static T2 convert_arg_1(const T2& x)
    {
        return x;
    };

    static std::string name()
    { 
        return "sum"; 
    };
};

struct Func_sqrt_2
{
    template<class T>    
    static T eval(const T& x1)
    {
        return sqrt(x1);
    }

    template<class T2>
    static T2 convert_arg_1(const T2& x)
    {
        return x;
    };

    static std::string name()
    { 
        return "sqrt_2"; 
    };
};

struct Func_uminus_2
{
    template<class T>    
    static T eval(const T& x1)
    {
        return -x1;
    }

    template<class T2>
    static T2 convert_arg_1(const T2& x)
    {
        return x;
    };

    static std::string name()
    { 
        return "uminus_2"; 
    };
};

struct Func_plus_11
{
    template<class T>    
    static T eval(const T& x1, const T& x2)
    {
        return twofold_sum(x1, x2);
    }

    template<class T2>
    static T2 convert_arg_1(const T2& x)
    {
        return T2(x.value);
    };

    template<class T2>
    static T2 convert_arg_2(const T2& x)
    {
        return T2(x.value);
    };

    static std::string name()
    { 
        return "plus_11"; 
    };
};

struct Func_plus_11_sort
{
    template<class T>    
    static T eval(const T& x1, const T& x2)
    {
        return twofold_sum_sorted(x1, x2);
    }

    template<class T2>
    static T2 convert_arg_1(const T2& x)
    {
        return T2(x.value);
    };

    template<class T2>
    static T2 convert_arg_2(const T2& x)
    {
        return T2(x.value);
    };

    static std::string name()
    { 
        return "plus_11_sort"; 
    };
};

struct Func_minus_11
{
    template<class T>    
    static T eval(const T& x1, const T& x2)
    {
        return twofold_minus(x1, x2);
    }

    template<class T2>
    static T2 convert_arg_1(const T2& x)
    {
        return T2(x.value);
    };

    template<class T2>
    static T2 convert_arg_2(const T2& x)
    {
        return T2(x.value);
    };

    static std::string name()
    { 
        return "minus_11"; 
    };
};

struct Func_minus_11_sort
{
    template<class T>    
    static T eval(const T& x1, const T& x2)
    {
        return twofold_minus_sorted(x1, x2);
    }

    template<class T2>
    static T2 convert_arg_1(const T2& x)
    {
        return T2(x.value);
    };

    template<class T2>
    static T2 convert_arg_2(const T2& x)
    {
        return T2(x.value);
    };

    static std::string name()
    { 
        return "minus_11_sort"; 
    };
};

struct Func_mult_11
{
    template<class T>    
    static T eval(const T& x1, const T& x2)
    {
        return twofold_mult(x1, x2);
    }

    template<class T2>
    static T2 convert_arg_1(const T2& x)
    {
        return T2(x.value);
    };

    template<class T2>
    static T2 convert_arg_2(const T2& x)
    {
        return T2(x.value);
    };

    static std::string name()
    { 
        return "mult_11"; 
    };
};

template<class T_base>
struct Func_mult_dekker
{
    template<class T>    
    static T eval(const T& x1, const T& x2)
    {
        (void)x1;
        (void)x2;
        return T(0.0);
    }

    template<class T2>
    static T2 convert_arg_1(const T2& x)
    {
        return T2(x.value);
    };

    template<class T2>
    static T2 convert_arg_2(const T2& x)
    {
        return T2(x.value);
    };

    static std::string name()
    { 
        return "mult_dekker"; 
    };
};

template<>
struct Func_mult_dekker<double>
{
    template<class T>    
    static T eval(const T& x1, const T& x2)
    {
        return twofold_mult_dekker(x1, x2);
    }

    template<class T2>
    static T2 convert_arg_1(const T2& x)
    {
        return T2(x.value);
    };

    template<class T2>
    static T2 convert_arg_2(const T2& x)
    {
        return T2(x.value);
    };

    static std::string name()
    { 
        return "mult_dekker"; 
    };
};

struct Func_div_11
{
    template<class T>    
    static T eval(const T& x1, const T& x2)
    {
        return twofold_div(x1, x2);
    }

    template<class T2>
    static T2 convert_arg_1(const T2& x)
    {
        return T2(x.value);
    };

    template<class T2>
    static T2 convert_arg_2(const T2& x)
    {
        return T2(x.value);
    };

    static std::string name()
    { 
        return "div_11"; 
    };
};

struct Func_plus_12
{
    template<class T>    
    static T eval(const T& x1, const T& x2)
    {
        return get_value(x1) + x2;
    }

    template<class T2>
    static T2 convert_arg_1(const T2& x)
    {
        return T2(x.value);
    };

    template<class T2>
    static T2 convert_arg_2(const T2& x)
    {
        return x;
    };

    static std::string name()
    { 
        return "plus_12"; 
    };
};

struct Func_plus_21
{
    template<class T>    
    static T eval(const T& x1, const T& x2)
    {
        return x1 + get_value(x2);
    }

    template<class T2>
    static T2 convert_arg_1(const T2& x)
    {
        return x;
    };

    template<class T2>
    static T2 convert_arg_2(const T2& x)
    {
        return T2(x.value);
    };

    static std::string name()
    { 
        return "plus_21"; 
    };
};

struct Func_plus_22
{
    template<class T>    
    static T eval(const T& x1, const T& x2)
    {
        return x1 + x2;
    }

    template<class T2>
    static T2 convert_arg_1(const T2& x)
    {
        return x;
    };

    template<class T2>
    static T2 convert_arg_2(const T2& x)
    {
        return x;
    };

    static std::string name()
    { 
        return "plus_22"; 
    };
};

//
struct Func_minus_12
{
    template<class T>    
    static T eval(const T& x1, const T& x2)
    {
        return get_value(x1) - x2;
    }

    template<class T2>
    static T2 convert_arg_1(const T2& x)
    {
        return T2(x.value);
    };

    template<class T2>
    static T2 convert_arg_2(const T2& x)
    {
        return x;
    };

    static std::string name()
    { 
        return "minus_12"; 
    };
};

struct Func_minus_21
{
    template<class T>    
    static T eval(const T& x1, const T& x2)
    {
        return x1 - get_value(x2);
    }

    template<class T2>
    static T2 convert_arg_1(const T2& x)
    {
        return x;
    };

    template<class T2>
    static T2 convert_arg_2(const T2& x)
    {
        return T2(x.value);
    };

    static std::string name()
    { 
        return "minus_21"; 
    };
};

struct Func_minus_22
{
    template<class T>    
    static T eval(const T& x1, const T& x2)
    {
        return x1 - x2;
    }

    template<class T2>
    static T2 convert_arg_1(const T2& x)
    {
        return x;
    };

    template<class T2>
    static T2 convert_arg_2(const T2& x)
    {
        return x;
    };

    static std::string name()
    { 
        return "minus_22"; 
    };
};

//
struct Func_mult_12
{
    template<class T>    
    static T eval(const T& x1, const T& x2)
    {
        return get_value(x1) * x2;
    }

    template<class T2>
    static T2 convert_arg_1(const T2& x)
    {
        return T2(x.value);
    };

    template<class T2>
    static T2 convert_arg_2(const T2& x)
    {
        return x;
    };

    static std::string name()
    { 
        return "mult_12"; 
    };
};

struct Func_mult_21
{
    template<class T>    
    static T eval(const T& x1, const T& x2)
    {
        return x1 * get_value(x2);
    }

    template<class T2>
    static T2 convert_arg_1(const T2& x)
    {
        return x;
    };

    template<class T2>
    static T2 convert_arg_2(const T2& x)
    {
        return T2(x.value);
    };

    static std::string name()
    { 
        return "mult_21"; 
    };
};

struct Func_mult_22
{
    template<class T>    
    static T eval(const T& x1, const T& x2)
    {
        return x1 * x2;
    }

    template<class T2>
    static T2 convert_arg_1(const T2& x)
    {
        return x;
    };

    template<class T2>
    static T2 convert_arg_2(const T2& x)
    {
        return x;
    };

    static std::string name()
    { 
        return "mult_22"; 
    };
};

//
struct Func_div_12
{
    template<class T>    
    static T eval(const T& x1, const T& x2)
    {
        return get_value(x1) / x2;
    }

    template<class T2>
    static T2 convert_arg_1(const T2& x)
    {
        return T2(x.value);
    };

    template<class T2>
    static T2 convert_arg_2(const T2& x)
    {
        return x;
    };

    static std::string name()
    { 
        return "div_12"; 
    };
};

struct Func_div_21
{
    template<class T>    
    static T eval(const T& x1, const T& x2)
    {
        return x1 / get_value(x2);
    }

    template<class T2>
    static T2 convert_arg_1(const T2& x)
    {
        return x;
    };

    template<class T2>
    static T2 convert_arg_2(const T2& x)
    {
        return T2(x.value);
    };

    static std::string name()
    { 
        return "div_21"; 
    };
};

struct Func_div_22
{
    template<class T>    
    static T eval(const T& x1, const T& x2)
    {
        return x1 / x2;
    }

    template<class T2>
    static T2 convert_arg_1(const T2& x)
    {
        return x;
    };

    template<class T2>
    static T2 convert_arg_2(const T2& x)
    {
        return x;
    };

    static std::string name()
    { 
        return "div_22"; 
    };
};

struct Func_fma
{
    template<class T>    
    static T eval(const T& x1, const T& x2, const T& x3)
    {
        return matcl::fma_f(x1, x2, x3);
    }

    template<class T>    
    static T eval_twofold(const T& x1, const T& x2, const T& x3)
    {
        return matcl::fma_dekker(x1, x2, x3);
    }

    static std::string name()
    { 
        return "fma"; 
    };
};

struct Func_save_load
{
    template<class T>    
    static T eval(const T& x1)
    {
        std::ostringstream os;
        os << std::setprecision(16) << x1;

        std::istringstream is(os.str());
        T x2;
        is >> x2;

        return x2;
    }

    static std::string name()
    { 
        return "save_load"; 
    };
};

}
