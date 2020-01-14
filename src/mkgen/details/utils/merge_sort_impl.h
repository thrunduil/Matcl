/*
 *  This file is a part of Matrix Computation Library (MATCL)
 *
 *  Copyright (c) Pawe³ Kowal 2019
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

#include "mkgen/mkgen_fwd.h"
#include "mkgen/utils/list.h"

#include <type_traits>

namespace merge_sort_impl
{

namespace mk        = matcl::mkgen;
namespace mkl       = mk::list;

//----------------------------------------------------------------------------------
//                              split
//----------------------------------------------------------------------------------
template <int K, class L>
struct split;

template <int K, class ... Elems>
struct split<K, mkl::list<Elems ...>> 
{
    private:
        using first         = split<K / 2, mkl::list<Elems ...>>;
        using second        = split<K - K / 2, typename first::right>;

    public:
        using left          = typename mkl::concat<typename first :: left, 
                                        typename second :: left> :: type;
        using right         = typename second::right;
};

template<class Elem, class ... Elems>
struct split<0, mkl::list<Elem, Elems...>> 
{
    using left  = mkl::list<>;
    using right = mkl::list<Elem, Elems ...>;
};

template<class Elem, class ... Elems>
struct split<1, mkl::list<Elem, Elems...>> 
{
    using left  = mkl::list<Elem>;
    using right = mkl::list<Elems...>;
};

template<int K>
struct split<K, mkl::list<>> 
{
    using left  = mkl::list<>;
    using right = mkl::list<>;
};

//----------------------------------------------------------------------------------
//                              subdivide
//----------------------------------------------------------------------------------
// split a list into two roughly equal lists
template <class TL>
struct subdivide : split<mkl::size<TL>::value / 2, TL> 
{};

//----------------------------------------------------------------------------------
//                              merge
//----------------------------------------------------------------------------------
template<class Compare, class A, class B>
struct merge;

template<bool Cond, class Compare, class L1, class L2>
struct merge_impl;

struct merge_nil{};

template<class Elem>
struct merge_elem
{
    using head  = Elem;
    using tail  = merge_nil;
};

template<class Compare, class L11, class L12, class L21, class L22>
struct merge<Compare, merge<Compare, L11, L12>, merge<Compare, L21, L22>> 
{
    using A         = merge<Compare, L11, L12>;
    using B         = merge<Compare, L21, L22>;

    static const bool cond  = Compare::template value<typename A :: head, typename B :: head>;

    using merger    = merge_impl<cond, Compare, A, B>;

    using head      = typename merger::head;
    using tail      = typename merger::tail;
};

template<class Compare, class L11, class L12, class E2>
struct merge<Compare, merge<Compare, L11, L12>, merge_elem<E2>> 
{
    using A         = merge<Compare, L11, L12>;
    using B         = merge_elem<E2>;

    static const bool cond  = Compare::template value<typename A :: head, E2>;
    using merger    = merge_impl<cond, Compare, A, B>;

    using head      = typename merger::head;
    using tail      = typename merger::tail;
};

template<class Compare, class E1, class L21, class L22>
struct merge<Compare, merge_elem<E1>, merge<Compare, L21, L22>> 
{
    using A         = merge_elem<E1>;
    using B         = merge<Compare, L21, L22>;

    static const bool cond  = Compare::template value<E1, typename B :: head>;
    using merger    = merge_impl<cond, Compare, A, B>;

    using head      = typename merger::head;
    using tail      = typename merger::tail;
};

template<class Compare, class E1, class E2>
struct merge<Compare, merge_elem<E1>, merge_elem<E2>> 
{
    using A         = merge_elem<E1>;
    using B         = merge_elem<E2>;

    static const bool cond  = Compare::template value<E1, E2>;
    using merger    = merge_impl<cond, Compare, A, B>;

    using head      = typename merger::head;
    using tail      = typename merger::tail;
};

template<class Compare, class L2>
struct merge<Compare, merge_nil, L2> 
{
    using head      = typename L2::head;
    using tail      = typename L2::tail;
};

template<class Compare, class L1>
struct merge<Compare, L1, merge_nil> 
{
    using head      = typename L1::head;
    using tail      = typename L1::tail;
};

template<class Compare, class L1, class L2>
struct merge_impl<true, Compare, L1, L2> 
{
    // a < b

    using head  = typename L1 :: head;
    using tail  = merge<Compare, typename L1 :: tail, L2>;
};

template<class Compare, class L1, class L2>
struct merge_impl<false, Compare, L1, L2> 
{
    // !(a < b)

    using head  = typename L2 :: head;
    using tail  = merge<Compare, L1, typename L2 :: tail>;
};

//----------------------------------------------------------------------------------
//                              merge_sort_impl
//----------------------------------------------------------------------------------
template<class Compare, class TL>
struct merge_sort_impl;

template <class Compare>
struct merge_sort_impl<Compare, mkl::list<>> 
{
    using type = merge_nil;
};

template<class Compare, class Elem>
struct merge_sort_impl<Compare, mkl::list<Elem>> 
{
    using type  = merge_elem<Elem>;    
};

template <class Compare, class Elem, class... Elems>
struct merge_sort_impl<Compare, mkl::list<Elem, Elems...>>
{
    using input_t       = mkl::list<Elem, Elems ... >;
    using subdivide_t   = subdivide<input_t>;

    using left          = typename subdivide_t::left;
    using right         = typename subdivide_t::right;

    using left_sort_t   = typename merge_sort_impl<Compare, left>::type;
    using right_sort_t  = typename merge_sort_impl<Compare, right>::type;

    using type          = merge<Compare, left_sort_t, right_sort_t>;
};

//----------------------------------------------------------------------------------
//                              expand
//----------------------------------------------------------------------------------
// convert the first Depth elements of a linked list List to list
// and return remaining elements unmodified
template<int Depth, class List>
struct expand
{
    static_assert(dependent_false<List>::value);
};

template<int Depth, class Compare, class L1, class L2>
struct expand<Depth, merge<Compare, L1, L2>>
{
    using merge_t   = merge<Compare, L1, L2>;

    using tail_m    = typename merge_t :: tail;
    using expander  = expand<Depth - 1, tail_m>;

    using head_e    = typename expander::head;
    using head_m    = typename merge_t :: head;

    using head      = typename mkl::push_front<head_e, head_m> :: type;
    using tail      = typename expander::tail;
};

template<class Compare, class L1, class L2>
struct expand<0, merge<Compare, L1, L2>>
{
    using merge_t   = merge<Compare, L1, L2>;
    using head      = mkl::list<typename merge_t::head>;
    using tail      = typename merge_t::tail;
};

template<int Depth>
struct expand<Depth, merge_nil>
{
    using head  = mkl::list<>;
    using tail  = merge_nil;
};

template<>
struct expand<0, merge_nil>
{
    using head  = mkl::list<>;
    using tail  = merge_nil;
};

//----------------------------------------------------------------------------------
//                              full_expand
//----------------------------------------------------------------------------------
// convert all elements in a linked_list L to list<...>
template<class L>
struct full_expand
{
    // number of elements converted in one step; cannot be too high due to recursion
    // depth limit
    static const int expand_depth   = 100;

    using expander1 = expand<expand_depth, L>;
    using tail      = typename expander1::tail;

    using tail_e    = typename full_expand<tail>::type;
    using head      = typename expander1::head;   

    using type      = typename mkl::concat<head, tail_e> :: type;
};

template<>
struct full_expand<merge_nil>
{
    using type      = mkl::list<>;
};

//----------------------------------------------------------------------------------
//                              merge_sort
//----------------------------------------------------------------------------------
template<class Compare, class TL>
struct merge_sort
{
    static_assert(mkl::is_list<TL>::value, "list<...> type required");

    using type_l    = typename merge_sort_impl<Compare, TL>::type;
    using type      = typename full_expand<type_l>::type;
};

}
