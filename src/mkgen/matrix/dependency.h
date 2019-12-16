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

#include "mkgen/details/mkgen_fwd.h"

namespace matcl { namespace mkgen 
{

//TODO
//------------------------------------------------------------------------
//                     single dependency
//------------------------------------------------------------------------

// type of dependency
struct dep_scalar;      // dependency on computed scalar
struct dep_extern;
struct dep_temp;
struct dep_return;
struct dep_computation;

// single dependency to data represented by tag Temp_Tag
template<class Temp_Tag, Integer Size, class Type>
struct dep
{
    // check arguments
    using check1    = typename details::check_dep_argument<Temp_Tag, Size, Type>::type;

    using tag                   = Temp_Tag;
    static const Integer size   = Size;
    using dep_type              = Type;
};

template<class Tag>
using extern_dep = dep<Tag, 0, dep_extern>;

template<class Tag>
using scalar_dep = dep<Tag, 0, dep_scalar>;

template<class Tag, Integer R, Integer C>
using temp_dep = dep<Tag, R * C, dep_temp>;

template<class Tag, Integer R, Integer C>
using return_dep = dep<Tag, R * C, dep_return>;

//------------------------------------------------------------------------
//                     list of dependencies
//------------------------------------------------------------------------

// list of dependencies; all items should be unique
template<class... T>
struct dps
{
    // check arguments
    using check1    = typename details::check_dps_argument<T...>::type;
};

// empty list of dependencies
using empty_deps = dps<>;

template<class Tag>
using extern_deps = dps<extern_dep<Tag>>;

template<class Tag, Integer R, Integer C>
using temp_deps = dps<temp_dep<Tag, R, C>>;

template<class Tag, Integer R, Integer C>
using return_deps = dps<return_dep<Tag, R, C>>;

//------------------------------------------------------------------------
//                     isa functions
//------------------------------------------------------------------------
// return true if T is dps<...> type
template<class T>
struct is_dps;

// return true if T is dep<...> type
template<class T>
struct is_dep;

//------------------------------------------------------------------------
//                     process dependencies
//------------------------------------------------------------------------

// construct list of unique dependencies dps<...> from two sets
// of dependencies dps<...>
template<class Deps_1, class Deps_2>
struct link_deps;

}}

#include "mkgen/details/matrix/dependency_impl.h"