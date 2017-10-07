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
#include "matcl-dynamic/register_function.h"
#include "matcl-dynamic/details/utils.h"

namespace matcl
{

#define MATCL_SELECT_OVERLOAD_1(unique_name,func,T1)                \
struct unique_name                                                  \
{                                                                   \
    using T1R = matcl::dynamic::details::make_type<T1>::type;       \
    static auto eval(T1R arg_) -> decltype(func(arg_))              \
    {                                                               \
        return func(arg_);                                          \
    };                                                              \
};                                                                  \

#define MATCL_SELECT_OVERLOAD_2(unique_name,func,T1,T2)             \
struct unique_name                                                  \
{                                                                   \
    using T1R = matcl::dynamic::details::make_type<T1>::type;       \
    using T2R = matcl::dynamic::details::make_type<T2>::type;       \
    static auto eval(T1R a1, T2R a2)                                \
                -> decltype(func(a1, a2))                           \
    {                                                               \
        return func(a1,a2);                                         \
    };                                                              \
};                                                                  \

#define MATCL_SELECT_OVERLOAD_OPERATOR(unique_name,op,T1,T2)        \
struct unique_name                                                  \
{                                                                   \
    using T1R = matcl::dynamic::details::make_type<T1>::type;       \
    using T2R = matcl::dynamic::details::make_type<T2>::type;       \
    static auto eval(T1R a1, T2R a2) -> decltype(a1 op a2)          \
    {                                                               \
        return a1 op a2;                                            \
    };                                                              \
};                                                                  \

#define MATCL_SELECT_OVERLOAD_3(unique_name,func,T1,T2,T3)          \
struct unique_name                                                  \
{                                                                   \
    using T1R = matcl::dynamic::details::make_type<T1>::type;       \
    using T2R = matcl::dynamic::details::make_type<T2>::type;       \
    using T3R = matcl::dynamic::details::make_type<T3>::type;       \
    static auto eval(T1R a1, T2R a2, T3R a3)                        \
                -> decltype(func(a1, a2, a3))                       \
    {                                                               \
        return func(a1,a2,a3);                                      \
    };                                                              \
};                                                                  \

#define MATCL_SELECT_OVERLOAD_4(unique_name,func,T1,T2,T3,T4)       \
struct unique_name                                                  \
{                                                                   \
    using T1R = matcl::dynamic::details::make_type<T1>::type;       \
    using T2R = matcl::dynamic::details::make_type<T2>::type;       \
    using T3R = matcl::dynamic::details::make_type<T3>::type;       \
    using T4R = matcl::dynamic::details::make_type<T4>::type;       \
    static auto eval(T1R a1, T2R a2, T3R a3, T4R a4)                \
                -> decltype(func(a1,a2,a3,a4))                      \
    {                                                               \
        return func(a1,a2,a3,a4);                                   \
    };                                                              \
};                                                                  \

#define MATCL_SELECT_OVERLOAD_TEMPL_0(unique_name,func,templ,of)                        \
struct unique_name                                                                      \
    : dynamic::register_function_template<unique_name,of,dynamic::object_type<templ>>   \
{                                                                                       \
    static auto eval()                                                                  \
        -> dynamic::object_type<decltype(func<templ>())>                                \
    {                                                                                   \
        using Ret = dynamic::object_type<decltype(func<templ>())>;                      \
        return Ret(func<templ>());                                                      \
    };                                                                                  \
};                                                                                      \


#define MATCL_SELECT_OVERLOAD_TEMPL_1(unique_name,func,templ,T1,of)                     \
struct unique_name                                                                      \
    : dynamic::register_function_template<unique_name,of,dynamic::object_type<templ>>   \
{                                                                                       \
    using T1_O  = matcl::dynamic::object_type<T1>;                                      \
    using T1_R  = matcl::dynamic::details::make_type<T1_O>::type;                       \
    static auto eval(T1_R arg_)                                                         \
        -> dynamic::object_type<decltype(func<templ>(arg_.get()))>                      \
    {                                                                                   \
        using Ret = dynamic::object_type<decltype(func<templ>(arg_.get()))>;            \
        return Ret(func<templ>(arg_.get()));                                            \
    };                                                                                  \
};                                                                                      \


#define MATCL_SELECT_OVERLOAD_TEMPL_2(unique_name,func,templ,T1,T2,of)                  \
struct unique_name                                                                      \
    : dynamic::register_function_template<unique_name,of,dynamic::object_type<templ>>   \
{                                                                                       \
    using T1_O  = matcl::dynamic::object_type<T1>;                                      \
    using T2_O  = matcl::dynamic::object_type<T2>;                                      \
    using T1_R  = matcl::dynamic::details::make_type<T1_O>::type;                       \
    using T2_R  = matcl::dynamic::details::make_type<T2_O>::type;                       \
    static auto eval(T1_R a, T2_R b)                                                    \
        -> dynamic::object_type<decltype(func<templ>(a.get(), b.get()))>                \
    {                                                                                   \
        using Ret = dynamic::object_type<decltype(func<templ>(a.get(),b.get()))>;       \
        return Ret(func<templ>(a.get(),b.get()));                                       \
    };                                                                                  \
};                                                                                      \


#define MATCL_REGISTER_SCALAR_FUNC(unique_name,func,type,object_func_name)  \
MATCL_SELECT_OVERLOAD_1(unique_name,func,type)                              \
template struct dynamic::register_function_ptr<decltype(&unique_name::eval),\
                    &unique_name::eval,object_func_name>;                   \

#define MATCL_REGISTER_BIN_FUNC(unique_name,func,T1,T2,object_func_name)    \
MATCL_SELECT_OVERLOAD_2(unique_name,func,T1,T2)                             \
template struct dynamic::register_function_ptr<decltype(&unique_name::eval),\
                    &unique_name::eval,object_func_name>;                   \

#define MATCL_REGISTER_OPERATOR(unique_name,op,T1,T2,object_func_name)      \
MATCL_SELECT_OVERLOAD_OPERATOR(unique_name,op,T1,T2)                        \
template struct dynamic::register_function_ptr<decltype(&unique_name::eval),\
                    &unique_name::eval,object_func_name>;                   \


#define MATCL_REGISTER_3_FUNC(unique_name,func,T1,T2,T3,object_func_name)   \
MATCL_SELECT_OVERLOAD_3(unique_name,func,T1,T2,T3)                          \
template struct dynamic::register_function_ptr<decltype(&unique_name::eval),\
                    &unique_name::eval,object_func_name>;                   \

#define MATCL_REGISTER_4_FUNC(unique_name,func,T1,T2,T3,T4,object_func_name)\
MATCL_SELECT_OVERLOAD_4(unique_name,func,T1,T2,T3,T4)                       \
template struct dynamic::register_function_ptr<decltype(&unique_name::eval),\
                    &unique_name::eval,object_func_name>;                   \

#define MATCL_REGISTER_TEMPL_0_FUNC(unique_name,func,templ,object_func_name)\
MATCL_SELECT_OVERLOAD_TEMPL_0(unique_name,func,templ,object_func_name)      \

#define MATCL_REGISTER_TEMPL_1_FUNC(unique_name,func,templ,type,object_func_name) \
MATCL_SELECT_OVERLOAD_TEMPL_1(unique_name,func,templ,type,object_func_name) \

#define MATCL_REGISTER_TEMPL_2_FUNC(unique_name,func,templ,type1,type2,object_func_name)\
MATCL_SELECT_OVERLOAD_TEMPL_2(unique_name,func,templ,type1,type2,object_func_name)      \

// utility macro to define a structure, that returns function name;
// such template is required by register_function class and associated
// macroes
#define MATCL_DEFINE_FUNCTION_NAME(name)        \
struct matcl_dynfunc_name_##name                \
{                                               \
    static function_name eval()                 \
    {                                           \
        static function_name func(#name);       \
        return func;                            \
    }                                           \
};

// return structure name generated by MATCL_DEFINE_FUNCTION_NAME
#define MATCL_FUNCTION_NAME(name)               \
matcl_dynfunc_name_##name

}