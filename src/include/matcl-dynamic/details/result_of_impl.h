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

#include "matcl-dynamic/object_type.h"
#include "matcl-dynamic/details/utils.h"

#define matcl_FUN_FUNCTOR(fun)                          \
namespace matcl { namespace impl                        \
{                                                       \
struct fun##_functor                                    \
{                                                       \
    template<class T>                                   \
    auto operator()(T arg_) -> decltype(fun(arg_));     \
                                                        \
    auto operator()(...) -> false_type;                 \
};                                                      \
}}                                                      \

#define matcl_OP_FUN_FUNCTOR(fun,op)                    \
namespace matcl { namespace impl                        \
{                                                       \
struct fun##_functor                                    \
{                                                       \
    template<class T>                                   \
    auto operator()(T arg1) -> decltype(op(arg1));      \
                                                        \
    auto operator()(...) -> false_type;                 \
};                                                      \
}}                                                      \

#define matcl_BIN_FUN_FUNCTOR(fun)                      \
namespace matcl { namespace impl                        \
{                                                       \
struct fun##_functor                                    \
{                                                       \
    template<class T, class S>                          \
    auto operator()(T a, S b) -> decltype(fun(a,b));    \
                                                        \
    auto operator()(...) -> false_type;                 \
};                                                      \
}}                                                      \

#define matcl_BINOP_FUN_FUNCTOR(fun,op)                 \
namespace matcl { namespace impl                        \
{                                                       \
struct fun##_functor                                    \
{                                                       \
    template<class T, class S>                          \
    auto operator()(T a, S b) -> decltype(a op b);      \
                                                        \
    auto operator()(...) -> false_type;                 \
};                                                      \
}}                                                      \

#define matcl_TEMPL_1_1_FUN_FUNCTOR(fun)                \
namespace matcl { namespace impl                        \
{                                                       \
struct fun##_functor                                    \
{                                                       \
    template<class T, class S>                          \
    auto operator()(T, S b) -> decltype(fun<T>(b));     \
                                                        \
    auto operator()(...) -> false_type;                 \
};                                                      \
}}                                                      \

#define matcl_TEMPL_1_0_FUN_FUNCTOR(fun)                \
namespace matcl { namespace impl                        \
{                                                       \
struct fun##_functor                                    \
{                                                       \
    template<class T>                                   \
    auto operator()(T) -> decltype(matcl::template fun<T>());\
                                                        \
    auto operator()(...) -> false_type;                 \
};                                                      \
}}                                                      \

#define matcl_TEMPL_1_2_FUN_FUNCTOR(fun)                \
namespace matcl { namespace impl                        \
{                                                       \
struct fun##_functor                                    \
{                                                       \
    template<class T, class S1, class S2>               \
    auto operator()(T, S1 b, S2 c)                      \
        -> decltype(fun<T>(b,c));                       \
                                                        \
    auto operator()(...) -> false_type;                 \
};                                                      \
}}                                                      \

#define MATCL_RESULT_OF_UNARY_FUN_IMPL(fun)                                                     \
namespace matcl { namespace result_of {                                                         \
template<class T>                                                                               \
struct has_##fun                                                                                \
{                                                                                               \
    using TR    = typename matcl::dynamic::details::make_type<T>::type;                         \
    using type  = typename std::result_of<matcl::impl::fun##_functor (TR)>::type;               \
    static const bool value = (std::is_same<type, matcl::impl::false_type>::value == false);    \
};                                                                                              \
template<class T, bool Has = has_##fun<T>::value>                                               \
struct result_of_##fun                                                                          \
{                                                                                               \
    using type          = typename has_##fun<T>::type;                                          \
    using type_object   = dynamic::object_type<typename std::decay<type>::type>;                \
};                                                                                              \
template<class T>                                                                               \
struct result_of_##fun<T,false>                                                                 \
{};                                                                                             \
}}                                                                                              \

#define MATCL_RESULT_OF_BINARY_FUN_IMPL(fun)                                                    \
namespace matcl { namespace result_of {                                                         \
template<class T, class S>                                                                      \
struct has_##fun                                                                                \
{                                                                                               \
    using TR = typename matcl::dynamic::details::make_type<T>::type;                            \
    using SR = typename matcl::dynamic::details::make_type<S>::type;                            \
    using type  = typename std::result_of<matcl::impl::fun##_functor (TR,SR)>::type;            \
    static const bool value = (std::is_same<type, matcl::impl::false_type>::value == false);    \
};                                                                                              \
template<class T, class S, bool Has = has_##fun<T,S>::value>                                    \
struct result_of_##fun                                                                          \
{                                                                                               \
    using type          = typename has_##fun<T,S>::type;                                        \
    using type_object   = dynamic::object_type<typename std::decay<type>::type>;                \
};                                                                                              \
template<class T, class S>                                                                      \
struct result_of_##fun<T,S,false>                                                               \
{};                                                                                             \
}}                                                                                              \

#define MATCL_RESULT_OF_BINARY_FUN_ADD_INT(fun)                                                 \
namespace matcl { namespace result_of {                                                         \
template<>                                                                                      \
struct result_of_##fun<Integer,Float, true>                                                     \
{                                                                                               \
    using type          = Real;                                                                 \
    using type_object   = dynamic::object_type<typename std::decay<type>::type>;                \
};                                                                                              \
template<>                                                                                      \
struct result_of_##fun<Float, Integer,true>                                                     \
{                                                                                               \
    using type          = Real;                                                                 \
    using type_object   = dynamic::object_type<typename std::decay<type>::type>;                \
};                                                                                              \
}}                                                                                              \

#define MATCL_RESULT_OF_BINARY_FUN_ADD_INT2(fun)                                                \
MATCL_RESULT_OF_BINARY_FUN_ADD_INT(fun)                                                         \
namespace matcl { namespace result_of {                                                         \
template<>                                                                                      \
struct result_of_##fun<Integer, Integer,true>                                                   \
{                                                                                               \
    using type          = Real;                                                                 \
    using type_object   = dynamic::object_type<typename std::decay<type>::type>;                \
};                                                                                              \
}}                                                                                              \

#define MATCL_RESULT_OF_TEMPL_1_1_FUN_IMPL(fun)                                                 \
namespace matcl { namespace result_of {                                                         \
template<class T, class S>                                                                      \
struct has_##fun                                                                                \
{                                                                                               \
    using TR = typename matcl::dynamic::details::make_type<T>::type;                            \
    using SR = typename matcl::dynamic::details::make_type<S>::type;                            \
    using type  = typename std::result_of<matcl::impl::fun##_functor (TR,SR)>::type;            \
    static const bool value = (std::is_same<type, matcl::impl::false_type>::value == false);    \
};                                                                                              \
template<class T, class S, bool Has = has_##fun<T,S>::value>                                    \
struct result_of_##fun                                                                          \
{                                                                                               \
    using type          = typename has_##fun<T,S>::type;                                        \
    using type_object   = dynamic::object_type<typename std::decay<type>::type>;                \
};                                                                                              \
template<class T, class S>                                                                      \
struct result_of_##fun<T,S,false>                                                               \
{};                                                                                             \
}}                                                                                              \

#define MATCL_RESULT_OF_TEMPL_1_2_FUN_IMPL(fun)                                                 \
namespace matcl { namespace result_of {                                                         \
template<class T, class S1, class S2>                                                           \
struct has_##fun                                                                                \
{                                                                                               \
    using TR = typename matcl::dynamic::details::make_type<T>::type;                            \
    using S1R = typename matcl::dynamic::details::make_type<S1>::type;                          \
    using S2R = typename matcl::dynamic::details::make_type<S2>::type;                          \
    using type  = typename std::result_of<matcl::impl::fun##_functor (TR,S1R,S2R)>::type;       \
    static const bool value = (std::is_same<type, matcl::impl::false_type>::value == false);    \
};                                                                                              \
template<class T, class S1, class S2, bool Has = has_##fun<T,S1,S2>::value>                     \
struct result_of_##fun                                                                          \
{                                                                                               \
    using type          = typename has_##fun<T,S1,S2>::type;                                    \
    using type_object   = dynamic::object_type<typename std::decay<type>::type>;                \
};                                                                                              \
template<class T, class S1, class S2>                                                           \
struct result_of_##fun<T,S1,S2,false>                                                           \
{};                                                                                             \
}}                                                                                              \

#define MATCL_RESULT_OF_TEMPL_1_0_FUN_IMPL(fun)                                                 \
namespace matcl { namespace result_of {                                                         \
template<class T>                                                                               \
struct has_##fun                                                                                \
{                                                                                               \
    using TR = typename matcl::dynamic::details::make_type<T>::type;                            \
    using type  = typename std::result_of<matcl::impl::fun##_functor (TR)>::type;               \
    static const bool value = (std::is_same<type, matcl::impl::false_type>::value == false);    \
};                                                                                              \
template<class T, bool Has = has_##fun<T>::value>                                               \
struct result_of_##fun                                                                          \
{                                                                                               \
    using type          = typename has_##fun<T>::type;                                          \
    using type_object   = dynamic::object_type<typename std::decay<type>::type>;                \
};                                                                                              \
template<class T>                                                                               \
struct result_of_##fun<T,false>                                                                 \
{};                                                                                             \
}}                                                                                              \

#define MATCL_RESULT_OF_UNARY(fun)          \
matcl_FUN_FUNCTOR(fun)                      \
MATCL_RESULT_OF_UNARY_FUN_IMPL(fun)         \

#define MATCL_RESULT_OF_BINARY(fun)         \
matcl_BIN_FUN_FUNCTOR(fun)                  \
MATCL_RESULT_OF_BINARY_FUN_IMPL(fun)        \

#define MATCL_RESULT_OF_UNARY_OP(fun,op)    \
matcl_OP_FUN_FUNCTOR(fun, op)               \
MATCL_RESULT_OF_UNARY_FUN_IMPL(fun)         \

#define MATCL_RESULT_OF_BINARY_OP(fun,op)   \
matcl_BINOP_FUN_FUNCTOR(fun,op)             \
MATCL_RESULT_OF_BINARY_FUN_IMPL(fun)        \

#define MATCL_RESULT_OF_TEMPL_1_0(fun)      \
matcl_TEMPL_1_0_FUN_FUNCTOR(fun)            \
MATCL_RESULT_OF_TEMPL_1_0_FUN_IMPL(fun)     \

#define MATCL_RESULT_OF_TEMPL_1_1(fun)      \
matcl_TEMPL_1_1_FUN_FUNCTOR(fun)            \
MATCL_RESULT_OF_TEMPL_1_1_FUN_IMPL(fun)     \

#define MATCL_RESULT_OF_TEMPL_1_2(fun)      \
matcl_TEMPL_1_2_FUN_FUNCTOR(fun)            \
MATCL_RESULT_OF_TEMPL_1_2_FUN_IMPL(fun)     \

