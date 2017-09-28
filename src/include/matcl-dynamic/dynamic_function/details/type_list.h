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

namespace matcl { namespace dynamic { namespace details
{

//-----------------------------------------------------------------------------
//                  type list
//-----------------------------------------------------------------------------

template<class ... A1>
struct make_list
{};

using nil = make_list<>;

//-----------------------------------------------------------------------------
//                  size
//-----------------------------------------------------------------------------
template<class list>
struct size
{
    //static_assert(false, "make_list type is required");
};
template<class ... A>
struct size<make_list<A...>>
{
    static const int value = sizeof ...(A);
};

//-----------------------------------------------------------------------------
//                  get_elem
//-----------------------------------------------------------------------------
template<class list,int elem>
struct get_elem
{
    //static_assert(false, "make_list type is required");
};
template<int elem>
struct get_elem<nil, elem>
{
    //static_assert(false, "list is empty");
};
template<class A, class ... B>
struct get_elem<make_list<A, B...>, 0>
{
    using type = A;
};
template<int elem, class A, class ... B>
struct get_elem<make_list<A, B...>, elem>
{
    using type = typename get_elem<make_list<B...>, elem-1>::type;
};

//-----------------------------------------------------------------------------
//                  tail
//-----------------------------------------------------------------------------
template<class list>
struct tail
{
    //static_assert(false, "make_list type is required");
};
template<>
struct tail<nil>
{
    //static_assert(false, "list is empty");
};

template<class A, class ... B>
struct tail<make_list<A,B...>>
{
    using type = make_list<B...>;
};

};};};