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

#include "matcl-simd/config.h"
#include "matcl-simd/complex/simd_complex.h"

namespace matcl { namespace level1
{

using simd_single_complex = simd::simd_single_complex;
using simd_double_complex = simd::simd_double_complex;

}};

namespace matcl { namespace level1 { namespace details
{

//helper classes
struct true_t{};
struct false_t{};

//check if T is a complex type
template<class T>
struct is_complex   { static const bool value = false; };

template<> struct is_complex<simd_double_complex>   { static const bool value = true; };
template<> struct is_complex<simd_single_complex>   { static const bool value = true; };

//check if T is an arithmetic type (integral types, floating point types and scalar types)
template<class T>
struct is_arithmetic_type
{
    static const bool value = std::is_arithmetic<T>::value || is_complex<T>::value;
};

//check if T is a single precision floating point type
template<class T>
struct is_single_precision  { static const bool value = false; };

template<> struct is_single_precision<float>                { static const bool value = true; };
template<> struct is_single_precision<simd_single_complex>  { static const bool value = true; };

//check if T is a double precision floating point type
template<class T>
struct is_double_precision  { static const bool value = false; };

template<> struct is_double_precision<double>               { static const bool value = true; };
template<> struct is_double_precision<simd_double_complex>  { static const bool value = true; };

//promote a scalar T to floating point type of the same 
//precision as Ret; if Ret is are not a floating point type
//or T is not an arithmetic type, then do nothing
template<class T, class Ret, 
        bool Is_scal = is_arithmetic_type<T>::value>
struct promote_scalar;

template<class T, class Ret>
struct promote_scalar<T, Ret, false>
{
    using type  = T;
};

template<class T, bool single_prec, bool double_prec,
        bool is_compl = is_complex<T>::value>
struct promote_scalar_impl;

template<class T, bool is_compl>
struct promote_scalar_impl<T, false, false, is_compl>
{
    using type = T;
};

template<class T>
struct promote_scalar_impl<T, true, false, false>
{
    using type = float;
};
template<class T>
struct promote_scalar_impl<T, false, true, false>
{
    using type = double;
};

template<class T>
struct promote_scalar_impl<T, true, false, true>
{
    using type = simd_single_complex;
};
template<class T>
struct promote_scalar_impl<T, false, true, true>
{
    using type = simd_double_complex;
};

template<class T, class Ret>
struct promote_scalar<T, Ret, true>
{
    static const bool is_single = is_single_precision<Ret>::value;
    static const bool is_double = is_double_precision<Ret>::value;

    using type  = typename promote_scalar_impl<T, is_single, is_double>::type;
};

//check if simd can be enabled for a given type
template<class T>
struct allow_simd                       { static const bool value = false; };

template<>
struct allow_simd<double>               { static const bool value = true; };

template<>
struct allow_simd<float>                { static const bool value = true; };

template<>
struct allow_simd<simd_double_complex>  { static const bool value = true; };

template<>
struct allow_simd<simd_single_complex>  { static const bool value = true; };

//check if all argumets Args are equal to T
template<class T, class ... Args>
struct equal_args;

template<class T>
struct equal_args<T>
{
    static const bool value = true;
};

template<class T, class S, class ... Args>
struct equal_args<T,S,Args...>
{
    static const bool value = std::is_same<T,S>::value
                                && equal_args<T, Args...>::value;
};

//if check_simd<Args...>::value is true, then vectorized version can be used
//Args are types of arrays
template<class Ret, class ... Args>
struct check_simd
{ 
    static const bool value = allow_simd<Ret>::value
                            && equal_args<Ret, Args...>::value; 
};

//if check_simd_scal<A, Args...>::value is true, then vectorized version can be used
//Args are types of arrays, A is a type of a scalar
template<class A, class Ret, class ... Args>
struct check_simd_scal
{ 
    static const bool value = allow_simd<Ret>::value
                            && equal_args<Ret, Args...>::value
                            && is_arithmetic_type<A>::value;
};

//if check_simd_scal2<A, B, Args...>::value is true, then vectorized version can be used
//Args are types of arrays, A, B are types of scalars
template<class A, class B, class Ret, class ... Args>
struct check_simd_scal2
{ 
    static const bool value = allow_simd<Ret>::value
                            && equal_args<Ret, Args...>::value
                            && is_arithmetic_type<A>::value
                            && is_arithmetic_type<B>::value;
};

//return real type
template<class T> struct real_type              { using type = T; };
template<> struct real_type<simd_double_complex>{ using type = double; };
template<> struct real_type<simd_single_complex>{ using type = float; };

//check if given scalar has a complex type but is real
template<class T, bool is_compl = is_complex<T>::value>
struct is_real;

template<class T>
struct is_real<T, false>
{
    static bool eval(const T&)  { return false; };
};

template<class T>
struct is_real<T, true>
{
    using TR    = typename real_type<T>::type;

    static bool eval(const T& v) 
    { 
        return imag(v) == TR(0.0); 
    };
};

//eval real 
template<class T> struct eval_real
{ 
    static T eval(const T& val) { return val; };
};
template<> struct eval_real<simd_double_complex>
{ 
    using T     = simd_double_complex;
    using TR    = double;

    static TR eval(const T& val) { return real(val); };
};
template<> struct eval_real<simd_single_complex>
{ 
    using T     = simd_single_complex;
    using TR    = float;

    static TR eval(const T& val) { return real(val); };
};

}}}