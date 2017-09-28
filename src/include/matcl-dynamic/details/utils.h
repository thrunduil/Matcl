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

#include "matcl-dynamic/config.h"
#include "matcl-core/config.h"
#include "matcl-dynamic/details/fwd_decls.h"

namespace matcl { namespace dynamic { namespace details
{

template<class T> struct make_type          { using type = const typename std::decay<T>::type&; };
template<class T> struct make_type<T&>      { using type = T&; };
template<class T> struct make_type<const T&>{ using type = const typename std::decay<T>::type&; };

template<class T> struct base_object_type   {};

template<class T> struct base_object_type<object_type<T>>
{
    using type = T;
};

template<class T>
struct is_object                        { static const bool value = std::is_base_of<object,T>::value; };

template<>
struct is_object<dynamic::object>       { static const bool value = true; };

template<class T>
struct is_object<object_type<T>>        { static const bool value = true; };

template<class T>
struct is_typed_object                  { static const bool value = false; };

template<class T>
struct is_typed_object<object_type<T>>  { static const bool value = true; };

template<class T>
struct is_any { static const bool value = std::is_same<T,any_type>::value; };

//is_assignable scalar; remove assignments that are treated as explicit
//see predefined_conversions
template<class To, class From>
struct is_assignable_scalar_impl
{
    static const bool value = true;
};
template<> struct is_assignable_scalar_impl<Integer, Float>
{
    static const bool value = false;
};
template<> struct is_assignable_scalar_impl<Integer, Real>
{
    static const bool value = false;
};
template<> struct is_assignable_scalar_impl<Integer, Complex>
{
    static const bool value = false;
};
template<> struct is_assignable_scalar_impl<Integer, Float_complex>
{
    static const bool value = false;
};
template<> struct is_assignable_scalar_impl<Float, Complex>
{
    static const bool value = false;
};
template<> struct is_assignable_scalar_impl<Real, Complex>
{
    static const bool value = false;
};
template<> struct is_assignable_scalar_impl<Float, Float_complex>
{
    static const bool value = false;
};
template<> struct is_assignable_scalar_impl<Real, Float_complex>
{
    static const bool value = false;
};

//is_constructible_scalar_impl scalar; remove conversions that are treated as casts
//see predefined_conversions
template<class To, class From>
struct is_constructible_scalar_impl
    :is_assignable_scalar_impl<To,From>
{};

template<class To, class From>
struct is_assignable_scalar
{
    using To_base   = typename std::decay<To>::type;
    using From_base = typename std::decay<From>::type;

    static const bool value = is_assignable_scalar_impl<To_base, From_base>::value;
};

//correct default implicit conversions for basic scalars
template<class To, class From>
struct is_convertible_scalar
{
    using To_base   = typename std::decay<To>::type;
    using From_base = typename std::decay<From>::type;

    static const bool value = is_assignable_scalar_impl<To_base, From_base>::value;
};

//correct default explicit conversions for basic scalars
template<class To, class From>
struct is_constructible_scalar
{
    using To_base   = typename std::decay<To>::type;
    using From_base = typename std::decay<From>::type;

    static const bool value = is_constructible_scalar_impl<To_base, From_base>::value;
};

template<class T, class Ret>
struct enable_if_nonobject 
    : std::enable_if<is_object<typename std::decay<T>::type>::value == false, Ret>
{};

template<class T, class Ret>
struct enable_if_nonobject_any 
    : std::enable_if<is_object<typename std::decay<T>::type>::value == false 
            && is_any<typename std::decay<T>::type>::value == false, Ret>
{};

template<class From, class To, bool Implicit>
struct is_convertible;

template<class From, class To>
struct is_convertible<From,To,true>
{
    static const bool value = std::is_convertible<From,To>::value
                            && is_convertible_scalar<To,From>::value;
};

template<class From, class To>
struct is_convertible<From,To,false>
{
    static const bool value = is_convertible<From,To,true>::value == false
                            && std::is_constructible<To,From>::value
                            && is_constructible_scalar<To,From>::value;
};

template<class From, class To>
struct is_assignable
{
    using To_ref    = typename std::add_lvalue_reference<To>::type;
    static const bool value = std::is_assignable<To_ref,From>::value
                            && is_assignable_scalar<To, From>::value;
};

template<class T, class S, class Ret>
struct enable_if_different : std::enable_if<std::is_same<T,S>::value == false, Ret>
{};

template<class From, class To, bool Implicit, class Type>
struct enable_if_different_conv
    : std::enable_if<std::is_same<From,To>::value == false 
            && is_convertible<From,To,Implicit>::value, Type>
{};

template<class From, class To, bool Implicit, class Type>
struct enable_if_nonobject_conv
    : std::enable_if<is_object<typename std::decay<From>::type>::value == false
            && is_convertible<From,To,Implicit>::value, Type>
{};

template<class From, class To, bool Implicit, class Type>
struct enable_if_conv
    : std::enable_if<is_convertible<From,To,Implicit>::value, Type>
{};

template<class Type, class Test, class Result>
struct enable_if_not
    : std::enable_if<std::is_same<typename std::decay<Type>::type, Test>::value == false, Type>
{};

template<class From, class To>
struct is_convertible_any
{
    static const bool value = is_convertible<From,To,true>::value 
                            || is_convertible<From,To,false>::value;
};

template<class From, class To, class Type>
struct enable_convert
    : std::enable_if<is_convertible_any<From, To>::value, Type>
{};

template<class From, class To, class Type>
struct enable_if_nonobject_assign
    : std::enable_if<is_object<typename std::decay<From>::type>::value == false
            && is_assignable<From,To>::value, Type>
{};

template<class From, class To, class Type>
struct enable_if_assign
    : std::enable_if<is_assignable<From,To>::value, Type>
{};

};};};
