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

#include "matcl-core/details/mpl.h"

namespace matcl { namespace mkgen { namespace list { namespace details
{

template<class Elem, class List, bool Is_Member>
struct insert_new_elem;

}}}}

namespace matcl { namespace mkgen { namespace list
{

//-------------------------------------------------------------------------
//                      list
//-------------------------------------------------------------------------

// compile time list with elements Elem
template<class ... Elems>
struct list{};

// empty list
using nil = list<>;

//-------------------------------------------------------------------------
//                      is_list
//-------------------------------------------------------------------------
// is_list<List>::value == true if List is a list type
template<class List>
struct is_list
{
    static const bool value = false;
};

template<class ... Elem>
struct is_list<list<Elem...>>
{
    static const bool value = true;
};

//-------------------------------------------------------------------------
//                      push_back
//-------------------------------------------------------------------------

// add element Elem at the end of the list List
template<class List, class Elem>
struct push_back
{
    static_assert(details::dependent_false<List>::value,
                  "List must be list<...> type");
};

template<class... Elems, class Elem>
struct push_back<list<Elems...>, Elem>
{
    using type = list<Elems...,Elem>;
};

//-------------------------------------------------------------------------
//                      push_front
//-------------------------------------------------------------------------

// add element Elem at the beginning of the list List
template<class List, class Elem>
struct push_front
{
    static_assert(details::dependent_false<List>::value,
                  "List must be list<...> type");
};

template<class... Elems, class Elem>
struct push_front<list<Elems...>,Elem>
{
    using type = list<Elem, Elems...>;
};

//-------------------------------------------------------------------------
//                      size
//-------------------------------------------------------------------------
//size of the list
template<class List>
struct size
{
    static_assert(details::dependent_false<List>::value,
                  "List must be list<...> type");
};

template<class ... Args>
struct size<list<Args...>>
{
    static const Integer value = sizeof... (Args);
};

//-------------------------------------------------------------------------
//                      is_member_tuple
//-------------------------------------------------------------------------
// check if type Elem is equal to one of types Elems
template<class Elem, class ... Elems>
struct is_member_tuple
{
    static_assert(details::dependent_false<Elem>::value, 
                "this type should not be instantiated");
};

template<class Elem> 
struct is_member_tuple<Elem>
{ 
    static const bool value = false;
};

template<class Elem, class ... E> 
struct is_member_tuple<Elem,Elem, E...> 
{ 
    static const bool value = true;
};

template<class Elem, class E1, class ... E>
struct is_member_tuple<Elem, E1, E...>
{
    static const bool value = is_member_tuple<Elem, E...>::value;
};

//-------------------------------------------------------------------------
//                      is_member
//-------------------------------------------------------------------------
// check if Elem is in the container Container
template<class Elem, class Container>
struct is_member
{
    static_assert(details::dependent_false<Elem>::value, 
                "this type should not be instantiated");
};

template<class Elem, template<class ...Args> class Container, class ... Elems>
struct is_member<Elem, Container<Elems...>>
    : is_member_tuple<Elem, Elems...>
{};

//-------------------------------------------------------------------------
//                      unique_list
//-------------------------------------------------------------------------
// remove duplicates from the list List
template<class List>
struct unique_list
{
    static_assert(details::dependent_false<List>::value,
                  "List must be list<...> type");
};

template<class Elem, class ...Elems>
struct unique_list<list<Elem, Elems...> >
{    
    using type_1                = typename unique_list<list::list<Elems...>>::type;
    static const bool found     = is_member<Elem, type_1>::value;
    using type                  = typename details::insert_new_elem<Elem, type_1, found>::type;    
};

template<>
struct unique_list<list<>>
{
    using type                  = list<>;
};

template<class Elem>
struct unique_list<list<Elem>>
{
    using type                  = list<Elem>;
};

//-------------------------------------------------------------------------
//                      elem_at_pos
//-------------------------------------------------------------------------
// get element from the list List at position Pos (0-based)
template<class List, Integer Pos>
struct elem_at_pos
{
    static_assert(details::dependent_false<List>::value,
                  "List must be list<...> type");
};

template<class Elem, class ... Elems, Integer Pos>
struct elem_at_pos<list<Elem, Elems...>, Pos>
{
    using type = typename elem_at_pos<list<Elems...>, Pos-1>::type;
};

template<class Elem, class ... Elems>
struct elem_at_pos<list<Elem, Elems...>, 0>
{
    using type = Elem;
};

template<Integer Pos>
struct elem_at_pos<list<>, Pos>
{
    static_assert(details::dependent_value_false<Integer, Pos>::value,
                  "invalid Pos");
};

//-------------------------------------------------------------------------
//                      elem_pos
//-------------------------------------------------------------------------
// return position of the element Elem in the list List (position is 0-based)
template<class List, class Elem, Integer Pos = 0>
struct elem_pos
{
    static_assert(details::dependent_false<List>::value,
                  "this type should not be instantiated");
};

template<class Elem, class... T, Integer Pos>
struct elem_pos<list<Elem, T...>, Elem, Pos> 
{
    static const Integer value = Pos;
};

template<class Elem1, class... T, class Elem2, Integer Pos>
struct elem_pos<list<Elem1, T...>, Elem2, Pos> 
{
    static const Integer value = elem_pos<list<T...>, Elem2, Pos+1>::value;
};

template<class Elem, Integer Pos>
struct elem_pos<list<>, Elem, Pos> 
{
    static_assert(details::dependent_false<Elem>::value,
                  "element not found in list");
};

}}}

namespace matcl { namespace mkgen { namespace list { namespace details
{

template<class Elem, class List, bool Is_Member>
struct insert_new_elem
{
    using type = List;
};

template<class Elem, class ... Elems>
struct insert_new_elem<Elem, list<Elems...>, false>
{
    using type = list::list<Elem,Elems...>;
};

}}}}