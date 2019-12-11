/*
 *  This file is a part of Matrix Computation Library (MATCL)
 *
 *  Copyright (c) Pawe� Kowal 2019
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

#include "mkgen/details/mkgen_fwd.h"
#include <type_traits>

// call has_member_function_x, where x is some identifier generates function
// that check, whether class C has member function
//          Ret x(Arg1 arg1, ..., Argk argk)

#define has_member_function_x(name)                                         \
template<class C, class T>                                                  \
struct has_member_function_##name                                           \
{                                                                           \
    static_assert(details::dependent_false<C>::value,                       \
            "second template parameter needs to be of function type.");     \
};                                                                          \
template<class C, class Ret, class... Args>                                 \
struct has_member_function_##name<C, Ret (Args...)>                         \
{                                                                           \
    /* attempt to call it and see if the return type is correct */          \
    template<typename T>                                                    \
    static constexpr auto check(T*)                                         \
        -> typename std::is_same<Ret,                                       \
            decltype( std::declval<T>().name( std::declval<Args>()... ) )   \
                        >::type;                                            \
                                                                            \
    template<typename>                                                      \
    static constexpr std::false_type check(...);                            \
                                                                            \
    using type  = decltype(check<C>(0));                                    \
                                                                            \
    static const bool value     = type::value;                              \
};

// call has_static_member_function_x, where x is some identifier generates
// function that check, whether class C has member function
//          static Ret x(Arg1 arg1, ..., Argk argk)

#define has_static_member_function_x(name)                                  \
template<class C, class T>                                                  \
    struct has_static_member_function_##name                                \
{                                                                           \
    static_assert(details::dependent_false<C>::value,                       \
            "second template parameter needs to be of function type.");     \
};                                                                          \
template<class C, class Ret, class... Args>                                 \
struct has_static_member_function_##name<C, Ret (Args...)>                  \
{                                                                           \
    /* attempt to call it and see if the return type is correct */          \
    template<typename T>                                                    \
    static constexpr auto check(T*)                                         \
        -> typename std::is_same<Ret,                                       \
                            decltype( T::name( std::declval<Args>()... ) )  \
                        >::type;                                            \
                                                                            \
    template<typename>                                                      \
    static constexpr std::false_type check(...);                            \
                                                                            \
    using type  = decltype(check<C>(0));                                    \
                                                                            \
    static const bool value     = type::value;                              \
};

// call has_static_template_function_x, where x is some identifier generates
// function that check, whether class C has member function
//          static Ret x<TArgs>(Arg1 arg1, ..., Argk argk)

#define has_static_member_template_function_x(name)                         \
template<class C, class TArgs, class T>                                     \
struct has_static_member_template_function_##name                           \
{                                                                           \
    static_assert(details::dependent_false<C>::value,                       \
            "third template parameter needs to be of function type.");      \
};                                                                          \
template<class C, class TArgs, class Ret, class... Args>                    \
struct has_static_member_template_function_##name<C, TArgs, Ret (Args...)>  \
{                                                                           \
    /* attempt to call it and see if the return type is correct */          \
    template<typename T>                                                    \
    static constexpr auto check(T*)                                         \
        -> typename std::is_same<Ret,                                       \
             decltype( T::template name<TArgs>( std::declval<Args>()... ))  \
                        >::type;                                            \
                                                                            \
    template<typename>                                                      \
    static constexpr std::false_type check(...);                            \
                                                                            \
    using type  = decltype(check<C>(0));                                    \
                                                                            \
    static const bool value     = type::value;                              \
};