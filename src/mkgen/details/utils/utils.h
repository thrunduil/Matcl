/*
 *  This file is a part of Matrix Computation Library (MATCL)
 *
 *  Copyright (c) Pawe³ Kowal 2019 - 2021
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

//TODO
namespace matcl { namespace mkgen { namespace details
{

//-----------------------------------------------------------------------
//                         array_item
//-----------------------------------------------------------------------
template<class Elem, Integer Step, class Type>
struct array_item
{
    using elem                  = Elem;
    using elem_type             = Type;
    static const Integer step   = Step;
};

// array item types
struct array_item_scalar{};
struct array_item_temp{};
struct array_item_extern{};

//-----------------------------------------------------------------------
//                         isa functions
//-----------------------------------------------------------------------
// return true if T is expr_plus_sd<...>
template<class T> 
struct is_plus_expr                                 {static const bool value = false; };

template<bool Flag, class ... T> 
struct is_plus_expr<expr_plus_sd<Flag, T...>>       {static const bool value = true; };

// return true if T is expr_mult_sd<...>
template<class T> 
struct is_mult_expr                                 {static const bool value = false; };

template<bool Flag, class... T> 
struct is_mult_expr<expr_mult_sd<Flag, T...>>       {static const bool value = true; };

//-----------------------------------------------------------------------
//                         static if
//-----------------------------------------------------------------------
template<bool Cond, class If_Expr_Type, class Else_Type>
struct static_if
{
    using type = If_Expr_Type;
};

template<class If_Expr_Type, class Else_Type>
struct static_if<false, If_Expr_Type, Else_Type>
{
    using type = Else_Type;
};

}}}