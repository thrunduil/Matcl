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

namespace matcl { namespace mkgen { namespace details
{

//-----------------------------------------------------------------------
//                      checks dependency
//-----------------------------------------------------------------------
template<class Deps>
struct check_valid_dps
{
    static_assert(details::dependent_false<Deps>::value, "type is not dps<>");
};

template<class ...T>
struct check_valid_dps<dps<T...>>
{
    // arguments T are already checked
    using type = void;
};

template<class Dep>
struct check_valid_dep
{
    static_assert(details::dependent_false<Deps>::value, "type is not dep<>");
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

}}};