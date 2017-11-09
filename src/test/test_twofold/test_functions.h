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

#include <iomanip>

namespace test_functions
{

namespace ms = matcl::simd;

using matcl::mp_float;
using matcl::twofold;

double get_value(double x)
{
    return x;
};

mp_float get_value(const mp_float& x)
{
    return x;
};

double get_value(const twofold& x)
{
    return x.value;
};


// missing functions
mp_float twofold_plus(const mp_float& a, const mp_float& b)         { return a + b; };
mp_float twofold_plus_sorted(const mp_float& a, const mp_float& b)  { return a + b; };
mp_float twofold_minus(const mp_float& a, const mp_float& b)        { return a - b; };
mp_float twofold_mult(const mp_float& a, const mp_float& b)         { return a * b; };
mp_float twofold_mult_dekker(const mp_float& a, const mp_float& b)  { return a * b; };
mp_float twofold_div(const mp_float& a, const mp_float& b)          { return a / b; };
mp_float twofold_sqrt(const mp_float& a)                            { return sqrt(a); };
mp_float twofold_sum(const mp_float& a)                             { return a.cast_float(); };

mp_float twofold_minus_sorted(const mp_float& a, const mp_float& b)
{
    if (abs(a) >= abs(b))
        return a - b; 
    else
        return b - a; 
};

double twofold_plus(double a, double b)         { return a + b; };
double twofold_plus_sorted(double a, double b)  { return a + b; };
double twofold_minus(double a, double b)        { return a - b; };
double twofold_mult(double a, double b)         { return a * b; };
double twofold_mult_dekker(double a, double b)  { return a * b; };
double twofold_div(double a, double b)          { return a / b; };
double twofold_sqrt(double a)                   { return matcl::sqrt(a); };
double twofold_sum(double a)                    { return a; };

double twofold_minus_sorted(double a, double b)
{
    if (abs(a) >= abs(b))
        return a - b; 
    else
        return b - a; 
};

force_inline
twofold twofold_plus(const twofold& a, const twofold& b)
{ 
    return matcl::twofold_plus(a.value, b.value); 
};

force_inline
twofold twofold_plus_sorted(const twofold& a, const twofold& b)
{ 
    if (abs(a.value) >= abs(b.value))
        return matcl::twofold_plus_sorted(a.value, b.value); 
    else
        return matcl::twofold_plus_sorted(b.value, a.value); 
};

force_inline
twofold twofold_minus(const twofold& a, const twofold& b)
{ 
    return matcl::twofold_minus(a.value, b.value); 
};

force_inline
twofold twofold_minus_sorted(const twofold& a, const twofold& b)
{ 
    if (abs(a.value) >= abs(b.value))
        return matcl::twofold_minus_sorted(a.value, b.value); 
    else
        return matcl::twofold_minus_sorted(b.value, a.value); 
};

force_inline
twofold twofold_mult(const twofold& a, const twofold& b)
{ 
    return matcl::twofold_mult(a.value, b.value); 
};

force_inline
twofold twofold_mult_dekker(const twofold& a, const twofold& b)
{ 
    return matcl::twofold_mult_dekker(a.value, b.value); 
};

force_inline
twofold twofold_div(const twofold& a, const twofold& b)
{ 
    return matcl::twofold_div(a.value, b.value); 
};

force_inline
twofold twofold_sqrt(const twofold& a)
{ 
    return matcl::twofold_sqrt(a.value); 
};

force_inline
twofold twofold_sum(const twofold& a)
{ 
    return twofold(a.sum());
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

struct Func_sum
{
    template<class T>    
    static T eval(const T& x1)
    {
        return twofold_sum(x1);
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
        return twofold_plus(x1, x2);
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
        return twofold_plus_sorted(x1, x2);
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

struct Func_mult_dekker
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
