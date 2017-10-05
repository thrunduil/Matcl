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

#include <utility>

namespace matcl { namespace dynamic { namespace details
{

template<int n_arg, bool void_return, class data_constructor,
            class function_traits>
struct function_evaler
{};

template<int n_arg, bool void_return, class data_constructor,
            class function_traits, class cl_type>
struct member_evaler
{};

//-------------------------------------------------------------------
//                  helpers
//-------------------------------------------------------------------
template<int n_args, class List, class... Build>
struct make_arg_pos_list
{
    using type  =  typename make_arg_pos_list
                <   n_args- 1, 
                    List, 
                    std::pair
                            <
                                std::integral_constant< int, n_args - 1>,
                                typename get_elem<List, n_args - 1>::type
                            >,
                    Build...
                > :: type;
};

template<class List, class... Build>
struct make_arg_pos_list<0, List, Build...>
{
    using type = make_list<Build...>;
};

template<class data_constructor, class Type>
struct get_arg
{
    using Pos = typename Type::first_type;
    using Arg = typename Type::second_type;

    static Arg eval(data_constructor& dh)
    {        
        return dh.template get<Pos::value, Arg>();
    };
};

template<class return_type, class data_constructor, class Function, 
            class Arguments>
struct make_call
{};

template<class cl_type, class return_type, class data_constructor,
            class Function, class Arguments>
struct make_call_member
{};

template<class return_type, class data_constructor, 
         class Function, class ... A>
struct make_call<return_type, data_constructor, Function, make_list<A...>>
{
    force_inline
    static void eval(data_constructor& dh,Function fun)
    {
        (void)dh;
        return_type ret = fun(get_arg<data_constructor, A>::eval(dh) ...);
        dh.set_return(std::move(ret));
    };
};

template<class data_constructor, 
         class Function, class ... A>
struct make_call<void, data_constructor, Function, make_list<A...>>
{
    force_inline
    static void eval(data_constructor& dh,Function fun)
    {
        (void)dh;
        fun(get_arg<data_constructor, A>::eval(dh) ...);
    };
};

template<class cl_type, class return_type, class data_constructor, 
         class Function, class ... A>
struct make_call_member<cl_type, return_type, data_constructor, 
                        Function, make_list<A...>>
{
    force_inline
    static void eval(data_constructor& dh,Function fun, cl_type obj)
    {
        return_type ret = (obj.*fun)(get_arg<data_constructor, A>::eval(dh) ...);
        dh.set_return(std::move(val));
    };
};

template<class cl_type, class data_constructor, 
         class Function, class ... A>
struct make_call_member<cl_type, void, data_constructor, 
                        Function, make_list<A...>>
{
    force_inline
    static void eval(data_constructor& dh,Function fun, cl_type obj)
    {
        (obj.*fun)(get_arg<data_constructor, A>::eval(dh) ...);
    };
};

//===================================================================
//                  function, void return
//===================================================================
template<int n_args, class data_constructor,class function_traits>
struct function_evaler<n_args,true,data_constructor,function_traits>
{
    using Function = typename function_traits::function;

    force_inline
    static void eval(data_constructor& dh,Function fun)
    {
        using in_type           = typename function_traits::input_type;
        using return_type       = typename function_traits::return_type;
        using arguments_type    = typename make_arg_pos_list<n_args, in_type>::type;
        
        make_call<return_type, data_constructor, Function, arguments_type>
                            ::eval(dh,fun);
    };
};

//===================================================================
//                  function, value return
//===================================================================
template<int n_args, class data_constructor,class function_traits>
struct function_evaler<n_args,false,data_constructor,function_traits>
{
    using Function = typename function_traits::function;

    force_inline
    static void eval(data_constructor& dh,Function fun)
    {
        using in_type           = typename function_traits::input_type;        
        using arguments_type    = typename make_arg_pos_list<n_args, in_type>::type;
        using return_type       = typename function_traits::return_type;

        make_call<return_type, data_constructor, Function, arguments_type>
                            ::eval(dh, fun);
    };
};

//===================================================================
//                  member, void return
//===================================================================
template<int n_args, class data_constructor,class function_traits, class cl_type>
struct member_evaler<n_args,true,data_constructor,function_traits,cl_type>
{
    using Function = typename function_traits::function;

    force_inline
    static void eval(data_constructor ds,Function fun, cl_type obj)
    {
        using in_type           = typename function_traits::input_type;
        using return_type       = typename function_traits::return_type;
        using arguments_type    = typename make_arg_pos_list<n_args, in_type>::type;

        make_call_member<cl_type, return_type, data_constructor, Function, arguments_type>
                            ::eval(ds,fun, obj);
    };
};

//===================================================================
//                  member, value return
//===================================================================
template<int n_args, class data_constructor,class function_traits, class cl_type>
struct member_evaler<n_args,false,data_constructor,function_traits,cl_type>
{
    using Function = typename function_traits::function;

    force_inline
    static void eval(data_constructor ds,Function fun, cl_type obj)
    {
        using in_type           = typename function_traits::input_type;
        using return_type       = typename function_traits::return_type;
        using arguments_type    = typename make_arg_pos_list<n_args, in_type>::type;

        make_call_member<cl_type, return_type, data_constructor, Function, 
                                arguments_type>::eval(ds,fun, obj);
    };
};

};};};