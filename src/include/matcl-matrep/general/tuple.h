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

#include <cassert>
#include <tuple>

namespace matcl
{

// tuple class extending std::tuple
template<class ... Args>
struct tuple : public std::tuple<Args...>
{
    using standard_tuple_type   = std::tuple<Args...>;

    // make tuple or copy constructor
    template<class ... Args2>
    tuple(Args2&& ... args) 
        : standard_tuple_type(std::forward<Args2>(args) ...) 
    {};

    // assignment
    template<class Arg>
    tuple& operator=(Arg&& arg)
    {
        standard_tuple_type::operator=(std::forward<Arg>(arg));
        return *this;
    };

    // get N-th element of tuple (1-based indexing)
    template<int N>
    auto get() const -> typename std::tuple_element<N-1, standard_tuple_type>::type
    {
        return std::get<N-1>(*this);
    };

    // convert matcl::tuple to std::tuple
    const standard_tuple_type&  to_standard_tuple() const & 
    { 
        return *this; 
    };
    
    // convert matcl::tuple to std::tuple
    standard_tuple_type&& to_standard_tuple() &&
    { 
        return std::move(*this); 
    };
};

// assign to many objects from tuple;
// for example: Matrix a,b; tie(a,b) = tuple<Matrix,Matrix>()
template<class ... Args>
tuple<Args& ...> tie(Args& ... args)
{ 
    return tuple<Args& ...>(args...); 
};

// get type of N-th element in the tuple (1-based indexing)
template<int N, class Tuple>
struct tuple_element : std::tuple_element<N-1, typename Tuple::standard_tuple_type>
{};

// construct matcl::tuple from std::tuple
template<class ...Args>
auto from_standard_tuple(const std::tuple<Args...>& t) -> tuple<Args...>
{ 
    return tuple<Args...>(t); 
};

// construct matcl::tuple from std::tuple
template<class ...Args>
auto from_standard_tuple(std::tuple<Args...>&& t) -> tuple<Args...>
{ 
    return tuple<Args...>(std::move(t)); 
};

};