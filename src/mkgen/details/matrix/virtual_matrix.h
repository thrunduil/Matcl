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

#include "mkgen/matrix/matrix.h"
#include "mkgen/TODO/expression/mat_assign.h"

namespace matcl { namespace mkgen { namespace details
{

//----------------------------------------------------------------------------------
//                              forward declarations
//----------------------------------------------------------------------------------

template<Integer M, Integer N, class Array1, class Deps, class ... Assign_List>
struct expand_virtual_impl;

template<Integer M, Integer N, class Deps, class ... Assign_List>
struct expand_virtual_impl2;

//----------------------------------------------------------------------------------
//                              make_colon_assignment
//----------------------------------------------------------------------------------

template<Integer M, Integer N, class Array1, class Deps, class Assign_Info>
struct make_colon_assignment
{
    static_assert(details::dependent_false<Array1>::value, 
                "this type should not be instantiated");
};

template<Integer M, Integer N, class Deps, class Assign_Info>
struct make_colon_assignment2
{
    static_assert(details::dependent_false<Deps>::value, 
                "this type should not be instantiated");
};

template<Integer M, Integer N, class Array1, class Deps, Integer M2, Integer N2, class Array2, class Colon_1>
struct make_colon_assignment<M, N, Array1, Deps, assign_item<M, N, M2, N2, Array2, Colon_1>>
{
    using array_type    = mkd::mat_assign_array_colon<M, N, Array1, Colon_1, M2, N2, Array2>;
    using type          = ct_matrix<M2, N2, array_type, Deps>;
};

template<Integer M, Integer N, class Deps, Integer M2, Integer N2, class Array2, class Colon_1>
struct make_colon_assignment2<M, N, Deps, assign_item<M, N, M2, N2, Array2, Colon_1>>
{
    using type          = list::list<Colon_1, ct_matrix<M2,N2,Array2, Deps>>;
};

//TODO
//----------------------------------------------------------------------------------
//                              expand_virtual_matrix
//----------------------------------------------------------------------------------
template<class Matrix>
struct expand_virtual_matrix
{
    static_assert(details::dependent_false<Matrix>::value, 
                "this type should not be instantiated");
};

template<Integer M, Integer N, class Array, class Deps>
struct expand_virtual_matrix<ct_matrix<M, N, Array, Deps>>
{
    using type = list::list<ct_matrix<M, N, Array, Deps>>;
};

template<Integer M, Integer N, class Array1, class Virt_Tag, class ... Assign_List, class Deps>
struct expand_virtual_matrix<ct_matrix<M, N, 
        mkd::mat_assign_array<M, N, Array1, mkd::virtual_array<Virt_Tag, Assign_List...>>, Deps>>
{
    using type = typename expand_virtual_impl<M, N, Array1, Deps, Assign_List...>::type;
};

template<Integer M, Integer N, class Array1, class Deps, class ... Assign_List>
struct expand_virtual_impl
{
    using type = list::list<typename make_colon_assignment<M, N, Array1, Deps, Assign_List>::type ... >;
};

//----------------------------------------------------------------------------------
//                              expand_virtual_matrix2
//----------------------------------------------------------------------------------
template<class Matrix>
struct expand_virtual_matrix2
{
    static_assert(details::dependent_false<Matrix>::value, 
                "this type should not be instantiated");
};

template<Integer M, Integer N, class Array, class Deps>
struct expand_virtual_matrix2<ct_matrix<M, N, Array, Deps>>
{
    using type = list::list<list::list<colon_all, ct_matrix<M, N, Array, Deps>>>;
};

template<Integer M, Integer N, class Virt_Tag, class ... Assign_List, class Deps>
struct expand_virtual_matrix2<ct_matrix<M, N, mkd::virtual_array<Virt_Tag, Assign_List...>, Deps>>
{
    using type = typename expand_virtual_impl2<M,N,Deps,Assign_List...>::type;
};

template<Integer M, Integer N, class Deps, class ... Assign_List>
struct expand_virtual_impl2
{
    using type = list::list<typename make_colon_assignment2<M, N, Deps, Assign_List>::type ... >;
};

}}}
