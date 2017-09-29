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

//#include "matcl-scalar/details/mpl.h"
#include "matcl-core/details/type_codes.h"
#include "matcl-core/details/val_struct_codes.h"
//#include "matcl-scalar/objects/object_type.h"
#include "matcl-core/details/utils.h"

namespace matcl { namespace details
{

template<class S, bool is_obj = md::is_object<S>::value>
struct integer_or_object
{
    using type = Integer;
};
template<class S>
struct integer_or_object<S,true>
{
    using type = Object;
};

template<class S, bool is_obj = md::is_object<S>::value>
struct bool_or_object
{
    using type = bool;
};
template<class S>
struct bool_or_object<S,true>
{
    using type = Object;
};

template<class S1, class S2, bool is_obj = is_object<S1>::value || is_object<S2>::value>
struct bool_or_object2
{
    using type = bool;
};
template<class S1, class S2>
struct bool_or_object2<S1,S2,true>
{
    using type = Object;
};

template<class Value_1, class Value_2, class Value_3>
struct unify_types2
{
    using type_1            = typename unify_types<Value_1, Value_2>::type;
    using type              = typename unify_types<type_1, Value_3>::type;
};

template<class Value_1, class Value_2, class Value_3, class Value_4>
struct unify_types3
{
    using type_1            = typename unify_types<Value_1, Value_2>::type;
    using type_2            = typename unify_types<Value_3, Value_4>::type;
    using type              = typename unify_types<type_1, type_2>::type;
};

template<class Value_1, class Value_2, class Value_3, class Value_4, class Value_5>
struct unify_types4
{
    using type_1            = typename unify_types<Value_1, Value_2>::type;
    using type_2            = typename unify_types<Value_3, Value_4>::type;
    using type_3            = typename unify_types<type_1, type_2>::type;
    using type              = typename unify_types<type_3, Value_5>::type;
};

template <class T> 
struct real_type_promote : real_type<typename promote_scalar<T>::type>
{};

template <class T, class S> 
struct real_type_promote_unify : unify_types<S, typename real_type<typename promote_scalar<T>::type>::type>
{};

template<class T1, class T2>
struct unify_types_promote : unify_types<typename promote_scalar<T1>::type, 
                                        typename promote_scalar<T2>::type>
{};

template<class T1, class T2>
struct real_unify_types_promote 
    : real_type<typename unify_types_promote <T1, T2>::type>
{};

template<class T1, class T2, class T3>
struct unify_types2_promote : unify_types2<typename promote_scalar<T1>::type, 
                                           typename promote_scalar<T2>::type,
                                           typename promote_scalar<T3>::type>
{};

template<class T1, class T2, class T3>
struct real_unify_types2_promote 
    : real_type<typename unify_types2_promote <T1, T2, T3>::type>
{};

template<class T1, class T2, class T3, class T4>
struct unify_types3_promote : unify_types3<typename promote_scalar<T1>::type, 
                                           typename promote_scalar<T2>::type,
                                           typename promote_scalar<T3>::type,
                                           typename promote_scalar<T4>::type>
{};

template<class T1, class T2, class T3, class T4, class T5>
struct unify_types4_promote : unify_types4<typename promote_scalar<T1>::type, 
                                           typename promote_scalar<T2>::type,
                                           typename promote_scalar<T3>::type,
                                           typename promote_scalar<T4>::type,
                                           typename promote_scalar<T5>::type>
{};

};};