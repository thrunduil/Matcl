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

#include "type_list.h"

namespace matcl { namespace dynamic { namespace details
{

//=============================================================================
//                  function traits
//=============================================================================
struct function_t{};
struct function_ptr_t{};
struct member_function_ptr_t{};
struct const_member_function_ptr_t{};
struct void_class_type{};

template<class Fun> struct funtion_traits_impl {};

template<class Ret, class ... A> 
struct funtion_traits_impl<Ret (A...)> 
{ 
    using function_type = function_t;
    using return_type   = Ret;
    using class_type    = void_class_type;
    using input_type    = make_list<A...>;
};

//-----------------------------------------------------------------------------
//                  function pointer
//-----------------------------------------------------------------------------
template<class Ret, class ... A> 
struct funtion_traits_impl<Ret (*)(A...)> 
{ 
    using function_type = function_ptr_t;
    using return_type   = Ret;
    using class_type    = void_class_type;
    using input_type    = make_list<A...>;
};

//-----------------------------------------------------------------------------
//                  member function pointer
//-----------------------------------------------------------------------------
template<class Ret, class X, class ... A> 
struct funtion_traits_impl<Ret (X::*)(A...)> 
{ 
    using function_type = member_function_ptr_t;
    using return_type   = Ret;
    using class_type    = X;
    using input_type    = make_list<A...>;
};

//-----------------------------------------------------------------------------
//                  constant member function pointer
//-----------------------------------------------------------------------------
template<class Ret, class X, class ... A> 
struct funtion_traits_impl<Ret (X::*)(A...) const> 
{ 
    using function_type = const_member_function_ptr_t;
    using return_type   = Ret;
    using class_type    = X;
    using input_type    = make_list<A...>;
};

template<class Fun, class fun_t>
struct proper_func_type                     { using type = Fun; };

template<class Fun>
struct proper_func_type<Fun,function_t>     { using type = Fun*; };

template<class Function> 
struct function_traits : public details::funtion_traits_impl<Function>
{
    private:
        using base_type = funtion_traits_impl<Function>;

    public:
        using function          = Function;
        using function_type     = typename base_type::function_type;
        using return_type       = typename base_type::return_type;
        using class_type        = typename base_type::class_type;
        using input_type        = typename base_type::input_type;

        static const int n_inputs = details::size<typename base_type::input_type>::value;
};

};};};