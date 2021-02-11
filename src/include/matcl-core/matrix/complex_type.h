/*
 *  This file is a part of Matrix Computation Library (MATCL)
 *
 *  Copyright (c) Pawe³ Kowal 2017 - 2021
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
#include "matcl-core/general/fwd_decls.h"

#include <complex>

namespace matcl { namespace details
{

template<class From, class To, bool no_loss>
struct enable_complex_conv{};

template<>
struct enable_complex_conv<float,double,true>
{ 
    using type = void*;
};

template<>
struct enable_complex_conv<double,float,false>
{ 
    using type = void*;
};

}};

namespace matcl
{

// Type for complex numbers used in matcl. 
//Some functions with the same name as in std namespace must be changed
//in order to satisfy the embedness property: f_R(x) == f_C(Complex(x)) for x in R
// and for any function f defined for real numbers (f_R) and complex numbers (f_C).
template<class T>
struct complex
{
    using real_type = T;

    std::complex<T>             value;

    // standard constructors
    complex()                   : value(T(),T()){};
    complex(T re)               : value(re,T()) {};
    complex(T re, T im)         : value(re,im)  {};

    complex(const complex<T>& other)
        : value(other.value)
    {};

    // conversion from std::complex must be explicit
    template<class S>
    explicit complex(const std::complex<S>& other)
        : value(other)
    {};

    // implicit conversion from complex<S>
    template<class S>
    complex(const complex<S>& other)
        : value(other.value)
    {};

    // explicit cast to boolean value
    explicit operator bool() const
    { 
        return std::real(value) != T() || std::imag(value) != T(); 
    };

    // save and load using boost::serialization library.
    void        serialize(oarchive_impl & ar, unsigned int version);
    void        serialize(iarchive_impl & ar, unsigned int version);
};

// get real and imaginary part of complex number
template<class T> 
T   imag(const complex<T>& val)     { return std::imag(val.value); };

template<class T>
T   real(const complex<T>& val)     { return std::real(val.value); };

//-----------------------------------------------------------------
//          operators +=, -=, *=, /= 
//          defined only for complex numbers
//-----------------------------------------------------------------
template<class T> 
complex<T>&         operator+=(complex<T>& arg, const complex<T>& val);

template<class T>
complex<T>&         operator-=(complex<T>& arg, const complex<T>& val);

template<class T> 
complex<T>&         operator*=(complex<T>& arg, const complex<T>& val);

template<class T> 
complex<T>&         operator/=(complex<T>& arg, const complex<T>& val);

//-----------------------------------------------------------------
//                  comparison operators
//-----------------------------------------------------------------
template<class T> 
bool                operator==(const complex<T>& arg, const complex<T>& val);

template<class T> 
bool                operator!=(const complex<T>& arg, const complex<T>& val);

template<class T> 
bool                operator>=(const complex<T>& arg, const complex<T>& val);

template<class T>
bool                operator<=(const complex<T>& arg, const complex<T>& val);

template<class T>
bool                operator<(const complex<T>& arg, const complex<T>& val);

template<class T>
bool                operator>(const complex<T>& arg, const complex<T>& val);

//-----------------------------------------------------------------
//                  arithmetic operators
//-----------------------------------------------------------------
template<class T> 
complex<T>          operator+(const complex<T>& arg1, const complex<T>& arg2);

template<class T> 
complex<T>          operator-(const complex<T>& arg1, const complex<T>& arg2);

template<class T> 
complex<T>          operator/(const complex<T>& arg1, const complex<T>& arg2);

template<class T> 
complex<T>          operator*(const complex<T>& arg1, const complex<T>& arg2);

template<class T> 
complex<T>          operator-(const complex<T>& arg1);

template<class T> 
complex<T>          operator+(const complex<T>& arg1, const T& arg2);

template<class T> 
complex<T>          operator-(const complex<T>& arg1, const T& arg2);

template<class T> 
complex<T>          operator/(const complex<T>& arg1, const T& arg2);

template<class T> 
complex<T>          operator*(const complex<T>& arg1, const T& arg2);

template<class T> 
complex<T>          operator+(const T& arg1, const complex<T>& arg2);

template<class T> 
complex<T>          operator-(const T& arg1, const complex<T>& arg2);

template<class T>
complex<T>          operator/(const T& arg1, const complex<T>& arg2);

template<class T> 
complex<T>          operator*(const T& arg1, const complex<T>& arg2);

// other operators and functions, which are defined also for matrices, 
// are specializations of generic functions and can be found elsewhere.

};
