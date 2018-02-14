/*
 *  This file is a part of Matrix Computation Library (MATCL)
 *
 *  Copyright (c) Pawe³ Kowal 2018
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
#include "matcl-scalar/IO/formatted_disp.h"
#include "matcl-mp/matcl_mp.h"
#include "matcl-scalar/lib_functions/sequence/limest.h"
#include "matcl-scalar/lib_functions/func_binary.h"

#include <string>

namespace matcl { namespace test
{

template<class T>
using lim_func  = std::function<T (const T&)>;

struct lim_func_1
{
    template<class T>
    static lim_func<T> get()
    {
        return [](const T& x) -> T { return sqrt(x) * log(x); };
    }
};

struct lim_func_2
{
    template<class T>
    static lim_func<T> get()
    {
        return [](const T& x) -> T { return (sqrt(x) + 1.0)/x; };
    }
};

struct lim_func_3
{
    template<class T>
    static lim_func<T> get()
    {
        return [](const T& x) -> T { return sin(x)/x; };
    }
};

struct lim_func_4
{
    template<class T>
    static lim_func<T> get()
    {
        return [](const T& x) -> T { return x / (1-exp(2*x)); }; 
    }
};

struct lim_func_5
{
    template<class T>
    static lim_func<T> get()
    {
        return [](const T& x) -> T { return (exp(x)-1) / x; }; 
    }
};

struct lim_func_6
{
    template<class T>
    static lim_func<T> get()
    {
        return [](const T& x) -> T { return (x * exp(x) - exp(x) + 1) / (x*x); }; 
    }
};

struct lim_func_7
{
    template<class T>
    static lim_func<T> get()
    {
        return [](const T& x) -> T { return (exp(x) - 1 - x - x*x/2) / (x*x*x); };
    }
};

struct lim_func_8
{
    template<class T>
    static lim_func<T> get(const T& s)
    {
        return [s](const T& x) -> T { return (matcl::pow(1.0 + x, s) - 1)/ x; };
    }
};

struct lim_func_9
{
    template<class T>
    static lim_func<T> get(const T& m, const T& n)
    {
        return [m, n](const T& x) -> T 
                            { return (matcl::pow(x, 1/n) - 1)/(matcl::pow(x, 1/m) - 1); };
    }
};

struct lim_func_10
{
    template<class T>
    static lim_func<T> get()
    {
        return [](const T& x) -> T { return x/(x-1.0) - 1.0/log(x); };
    }
};

struct lim_func_11
{
    template<class T>
    static lim_func<T> get()
    {
        return [](const T& x) -> T { return matcl::pow(abs(x),x);};
    }
};

struct lim_func_12
{
    template<class T>
    static lim_func<T> get()
    {
        return [](const T& x) -> T { return sqrt(log(x+ 1)) - sqrt(log(x));};
    }
};

struct lim_func_13
{
    template<class T>
    static lim_func<T> get()
    {
        return [](const T& x) -> T { return (x - sin(x)) /(x+sin(x)); };
    }
};

struct lim_func_14
{
    template<class T>
    static lim_func<T> get()
    {
        return [](const T& x) -> T { return matcl::pow(x,5)/exp(x);};
    }
};

struct lim_func_15
{
    template<class T>
    static lim_func<T> get()
    {
        return [](const T& x) -> T { return x*log(x); };
    }
};

struct lim_func_16
{
    template<class T>
    static lim_func<T> get()
    {
        return [](const T& x) -> T { return ((log(log(x)+log(log(x)))-log(log(x))) 
                                                            / log(log(x) +log(log(log(x)))))*log(x); };
    }
};

struct lim_func_17
{
    template<class T>
    static lim_func<T> get()
    {
        return [](const T& x) -> T { return 1/cos(x) * exp(-x);};
    }
};

struct lim_func_18
{
    template<class T>
    static lim_func<T> get()
    {
        return [](const T& x) -> T { return sin(1.0/x - exp(-x)) - sin(1.0/x);};
    }
};

struct lim_func_19
{
    template<class T>
    static lim_func<T> get()
    {
        return [](const T& x) -> T { return sqrt(abs(x)); };
    }
};

struct lim_func_20
{
    template<class T>
    static lim_func<T> get()
    {
        return [](const T& x) -> T { return x; };
    }
};

struct lim_func_21
{
    template<class T>
    static lim_func<T> get()
    {
        return [](const T& x) -> T { return 1.0/matcl::sqrt(matcl::abs(x)); };
    }
};

struct lim_func_22
{
    template<class T>
    static lim_func<T> get()
    {
        return [](const T& x) -> T 
        { 
            return (x < 0 ? -1.0 : 1.0) 
                    * matcl::cos(1.0/seq_helpers::cast_double<T>::eval(x)); 
        };
    }
};

struct lim_func_23
{
    template<class T>
    static lim_func<T> get()
    {
        return [](const T& x) -> T { return (x < 0 ? -1.0 : 1.0) 
                                        * matcl::pow(matcl::abs(x), 1.0/32.0); };
    }
};

struct lim_func_24
{
    template<class T>
    static lim_func<T> get()
    {
        return [](const T& x) -> T { return x + ((matcl::abs(x) < 1.0e-4) 
                                                ? matcl::pow(matcl::abs(x), -1.0/32.0) : T(0.0)); };
    }
};

}}
