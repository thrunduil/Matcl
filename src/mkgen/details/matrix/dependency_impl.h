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

#include "mkgen/matrix/dependency.h"
#include "mkgen/utils/list.h"

namespace matcl { namespace mkgen { namespace details
{

//-----------------------------------------------------------------------
//                      checks dependency
//-----------------------------------------------------------------------
template<class Dep>
struct check_valid_dep
{
    static_assert(md::dependent_false<Deps>::value, "type is not dep<>");
};

template<class Temp_Tag, Integer Size, class Type>
struct check_valid_dep<dep<Temp_Tag, Size, Type>>
{
    // arguments are already checked
    using type = void;
};

template<class... T>
struct check_dps_argument
{
    using type = dps<typename check_valid_dep<T>::type ...>;
};

template<class Temp_Tag, Integer Size, class Type>
struct check_dep_argument
{
    using type1 = typename check_valid_dep_tag<Temp_Tag>::type;
    using type2 = typename check_valid_dep_size<Size>::type;
    using type3 = typename check_valid_dep_type<Type>::type;

    using type  = void;
};

template<class Tag>
struct check_valid_dep_tag
{
    //TODO
    using type = void;
};

template<Integer Size>
struct check_valid_dep_size
{
    //TODO
    using type = void;
};

template<class Type>
struct check_valid_dep_type
{
    static const bool value =   std::is_same<Type, dep_scalar>::value
                            ||  std::is_same<Type, dep_extern>::value
                            ||  std::is_same<Type, dep_temp>::value
                            ||  std::is_same<Type, dep_return>::value
                            ||  std::is_same<Type, dep_computation>::value
                            ;

    static_assert(value == true, "invalid type argument passed to dep<>");

    using type  = void;
};

//-----------------------------------------------------------------------
//                      insert_new_dep
//-----------------------------------------------------------------------
// if Is_Member == false, then add New_Elem to the list of dependencies List
template<class New_Elem, class List, bool Is_Member>
struct insert_new_dep
{
    using type = List;
};

template<class New_Elem, class ... Deps>
struct insert_new_dep<New_Elem, dps<Deps...>, false>
{
    using type = dps<New_Elem, Deps...>;
};

//-----------------------------------------------------------------------
//                      get_new_dps
//-----------------------------------------------------------------------
// get sublist of dependencies from the set Deps1, which does not belong to the set Deps2
template<class Deps1, class Deps2>
struct get_new_dps
{};

template<class Deps2>
struct get_new_dps<dps<>, Deps2>
{
    using type  = dps<>;
};

template<class Dep1, class Deps2>
struct get_new_dps<dps<Dep1>, Deps2>
{
    static const bool is_member = list::is_member<Dep1, Deps2>::value;
    using type                  = typename std::conditional<is_member, dps<>, dps<Dep1>> :: type;
};

template<class Dep1, class ...Deps1, class Deps2>
struct get_new_dps<dps<Dep1, Deps1...>, Deps2>
{    
    using type_1                = typename get_new_dps<dps<Deps1...>, Deps2>::type;
    static const bool is_member = list::is_member<Dep1, Deps2>::value;
    using type                  = typename insert_new_dep<Dep1, type_1, is_member>::type;    
};

//-----------------------------------------------------------------------
//                      merge_deps_set
//-----------------------------------------------------------------------
// merge two sets of dependencies
template<class Dep1, class Dep2>
struct merge_deps_set
{};

template<class ... Dep1, class ... Dep2>
struct merge_deps_set<dps<Dep1...>, dps<Dep2...>>
{
    using type = dps<Dep1..., Dep2...>;
};

//-----------------------------------------------------------------------
//                      link_deps_impl
//-----------------------------------------------------------------------
template<class Deps_1, class Deps_2>
struct link_deps_impl
{
    static_assert(md::dependent_false<v>::value, 
                "this type should not be instantiated");
};

template<>
struct link_deps_impl<dps<>, dps<>>
{
    using type  = dps<>;
};

template<class Deps1>
struct link_deps_impl<Deps1, dps<>>
{
    using type  = Deps1;
};

template<class Deps2>
struct link_deps<dps<>, Deps2>
{
    using type  = Deps2;
};

template<class Dep1>
struct link_deps_impl<Dep1, Dep1>
{
    using type = Dep1;
};

// this specialization is required in order to avoid ambiguities between
// link_deps_impl<Dep1,Dep1>
// link_deps_impl<dps<Dep1, Deps1 ... >, Dep2> 
template<class Dep1, class... Deps1>
struct link_deps_impl<dps<Dep1, Deps1 ... >, dps<Dep1, Deps1 ... >>
{
    using type      = dps<Dep1, Deps1 ... >;
};

template<class Dep1, class... Deps1, class Dep2>
struct link_deps_impl<dps<Dep1, Deps1 ... >, Dep2>
{
    using dps_1     = dps<Dep1, Deps1 ... >;

    // get sublist of dependencies from the set dps_1, which does not belong to the set Dep2
    using unique    = typename get_new_dps<dps_1, Dep2>::type;

    // merge sets, these sets are disjoint
    using type      = typename merge_deps_set<unique, Dep2>::type;
};

template<class Dep1, class... Deps1>
struct link_deps_impl<dps<Dep1, Deps1 ... >, dps<>>
{
    using type = dps<Dep1, Deps1 ... >;
};

}}};

namespace matcl { namespace mkgen 
{

namespace mkd = matcl::mkgen::details;

//-----------------------------------------------------------------------
//                      link_deps
//-----------------------------------------------------------------------
template<class Deps_1, class Deps_2>
struct link_deps
{
    static const bool is_deps_1 = is_dps<Deps_1>::value;
    static const bool is_deps_2 = is_dps<Deps_2>::value;

    static_assert(is_deps_1 == true && is_deps_2 == true, "arguments must have type dps<...>");

    using type = typename mkd::link_deps_impl<Deps_1, Deps_2>::type;
};

//------------------------------------------------------------------------
//                     isa functions
//------------------------------------------------------------------------
template<class T>
struct is_dps
{
    static const bool value = false;
};

template<class ... Args>
struct is_dps<dps<Args...>>
{
    static const bool value = true;
};

template<class T>
struct is_dep
{
    static const bool value = false;
};

template<class Temp_Tag, Integer Size, class Type>
struct is_dep<dep<Temp_Tag, Size, Type>>
{
    static const bool value = true;
};


}};