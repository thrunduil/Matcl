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

#include "mkgen/matrix/matrix.h"

namespace matcl { namespace mkgen { namespace details
{

//----------------------------------------------------------------------------------
//                              forward declarations
//----------------------------------------------------------------------------------

template<Integer M, Integer N, class Array1, class Deps, class ... Assign_List>
struct expand_virtual_impl;

template<Integer M, Integer N, class Deps, class ... Assign_List>
struct expand_virtual_impl2;

template<Integer M1, Integer N1, class Colon_1, Integer M2, Integer N2>
struct is_assignment_valid;

template<Integer M1, Integer N1, class Colon_1>
struct is_assignment_valid_scalar;

template<class Array, class Item>
struct add_to_virtual_array;

template<bool Found, Integer Row, Integer Col, class Item, class ... Items>
struct get_assignment_item;

template<Integer Row, Integer Col, class Item>
struct get_assignment_item_impl;

template<Integer Row, Integer Col, class Item>
struct is_member_assign;

//----------------------------------------------------------------------------------
//                              virtual_assign_item
//----------------------------------------------------------------------------------

// make assignment A(Colon_1) = B; where A is a matrix of size M1 x N1, B is a matrix
// of size M2 x N2 with array Array2
template<Integer M1, Integer N1, Integer M2, Integer N2, class Array2, class Colon_1>
struct virtual_assign_item
{
    static_assert(is_assignment_valid<M1, N1, Colon_1, M2, N2>::value, "invalid assignment");
};

// make assignment A(Colon_1) = Scalar; where A is a matrix of size M1 x N1
template<Integer M1, Integer N1, class Scalar, class Colon_1>
struct virtual_assign_item_scalar
{
    static_assert(is_assignment_valid_scalar<M1, N1, Colon_1>::value, "invalid assignment");
    static_assert(is_scalar<Scalar>::value, "scalar required");
};

//----------------------------------------------------------------------------------
//                              is_assignment_valid
//----------------------------------------------------------------------------------
// check if assignment A(Colon_1) = B is valid; where A is M1 x N1 matrix and B is
// M2 x N2 matrix
template<Integer M1, Integer N1, class Colon_1, Integer M2, Integer N2>
struct is_assignment_valid
{
    using colon                     = Colon_1;
    static const Integer A_size     = M1 * N1;
    static const Integer colon_size = colon_func::size<colon, A_size>::value;

    static const bool value         = (M2 == 1 || N2 == 1) 
                                    && (M2 * N2 == colon_size);
};

// check if assignment A(Colon_1) = B is valid; where A is M1 x N1 matrix and B is
// scalar
template<Integer M1, Integer N1, class Colon_1>
struct is_assignment_valid_scalar
{
    using colon                     = Colon_1;
    static const Integer A_size     = M1 * N1;
    static const Integer colon_size = colon_func::size<colon, A_size>::value;

    static const bool value         = true;
};

//----------------------------------------------------------------------------------
//                              mat_virtual_assign_1
//----------------------------------------------------------------------------------
// make virtual matrix assignment A(Colon_1) = B
template<Matrix A, Mat_or_scalar B, Colon Colon_1>
struct mat_virtual_assign_1
{
    static_assert(md::dependent_false<A>::value, "invalid arguments");
};

template<Integer M1, Integer N1, Mat_array Array_1, DPS Deps_1, 
         Integer M2, Integer N2, Mat_array Array_2, DPS Deps_2, 
         Colon Colon_1>
struct mat_virtual_assign_1<ct_matrix<M1, N1, Array_1, Deps_1>,
                            ct_matrix<M2, N2, Array_2, Deps_2>, Colon_1>
{
    using A             = ct_matrix<M1, N1, Array_1, Deps_1>;
    using B             = ct_matrix<M2, N2, Array_2, Deps_2>;

    static const bool is_vA = is_virtual_matrix<A>::value;
    static_assert(is_vA == true, "virtual matrix is required");
    
    //array_1 has type mkd::virtual_array

    using item          = virtual_assign_item<M1, N1, M2, N2, Array_2, Colon_1>;
    using array_type    = typename add_to_virtual_array<Array_1, item>::type;
    using deps          = typename link_deps<Deps_1, Deps_2>::type;

    // output
    using type          = ct_matrix<M1, N1, array_type, deps>;
};

template<Integer M1, Integer N1, Mat_array Array_1, DPS Deps_1, 
         Scal_data Array_2, DPS Deps_2, Colon Colon_1>
struct mat_virtual_assign_1<ct_matrix<M1, N1, Array_1, Deps_1>,
                            ct_scalar<Array_2,Deps_2>, Colon_1>
{
    using A                 = ct_matrix<M1, N1, Array_1, Deps_1>;
    using B                 = ct_scalar<Array_2,Deps_2>;

    static const bool is_vA = is_virtual_matrix<A>::value;
    static_assert(is_vA == true, "virtual matrix is required");
    
    //array_1 has type mkd::virtual_array

    using item          = virtual_assign_item_scalar<M1, N1, B, Colon_1>;
    using array_type    = typename add_to_virtual_array<Array_1, item>::type;
    using deps          = typename link_deps<Deps_1, Deps_2>::type;

    // output
    using type          = ct_matrix<M1, N1, array_type, deps>;
};

//----------------------------------------------------------------------------------
//                              add_to_virtual_array
//----------------------------------------------------------------------------------
// insert Item to virtual_array assignments
template<class Array, class Item>
struct add_to_virtual_array
{
    static_assert(md::dependent_false<Array>::value, "virtual_array required");
};

template<class Tag, class ... Args, class Item>
struct add_to_virtual_array<mkd::virtual_array<Tag, Args...>, Item>
{
    using type = mkd::virtual_array<Tag, Item, Args...>;
};

//----------------------------------------------------------------------------------
//                              get_virtual_array_assignment
//----------------------------------------------------------------------------------
// implements virtual_array::get_element
template<Integer Row, Integer Col, class... Items>
struct get_virtual_array_assignment;

template<Integer Row, Integer Col, class Item, class ... Items>
struct get_virtual_array_assignment<Row, Col, Item, Items...>
{
    static const bool found = is_member_assign<Row,Col, Item>::value;
    using type = typename get_assignment_item<found, Row, Col, Item, Items...>::type;
};

template<Integer Row, Integer Col>
struct get_virtual_array_assignment<Row, Col>
{
    static_assert(md::dependent_value_false<Integer,Row>::value,
                  "element is not assigned to given virtual matrix");
};

//----------------------------------------------------------------------------------
//                              get_assignment_item
//----------------------------------------------------------------------------------
// if Found == true, then get element from given assignment Item, 
// otherwise search in remaining assignments
template<bool Found, Integer Row, Integer Col, class Item, class ... Items>
struct get_assignment_item
{
    using type = typename get_virtual_array_assignment<Row, Col, Items...>::type;
};

template<Integer Row, Integer Col, class Item, class ... Items>
struct get_assignment_item<true, Row, Col, Item, Items...>
{
    using type = typename get_assignment_item_impl<Row,Col, Item>::type;
};

//----------------------------------------------------------------------------------
//                              get_assignment_item_impl
//----------------------------------------------------------------------------------
template<Integer Row, Integer Col, class Item>
struct get_assignment_item_impl
{
    static_assert(md::dependent_false<Item>::value, "invalid arguments");
};

template<Integer Row, Integer Col, 
        Integer M1, Integer N1, Integer M2, Integer N2, class Array2, class Colon_1>
struct get_assignment_item_impl<Row, Col, 
                virtual_assign_item<M1, N1, M2, N2, Array2, Colon_1>>
{
    // we assume, that RHS is assignment is always a vector
    static_assert(M2 == 1 || N2 == 1, "rhs of assignment should be a vector");

    // matrix has size M1 x N1, this (Row, Col) is at position pos:
    static const Integer pos    = (Col - 1) * M1 + Row;

    // which colon element points to position pos?
    static const Integer rel_pos= get_relative_pos<pos, Colon_1>::value;
    
    // size of RHS
    static const Integer size_B = M2 * N2;
    
    // rel_pos must be a valid position in RHS matrix
    static_assert(rel_pos >= 1 && rel_pos <= size_B, "invalid element");

    // convert position to (row, col) pair
    static const Integer col2   = (rel_pos - 1) / M2 + 1;
    static const Integer row2   = rel_pos - (col2 - 1) * M2;

    using type = typename Array2::template get_element<row2, col2>::type;
};

template<Integer Row, Integer Col, 
        Integer M1, Integer N1, class Scalar, class Colon_1>
struct get_assignment_item_impl<Row, Col, 
                    virtual_assign_item_scalar<M1, N1, Scalar, Colon_1>>
{
    // TODO: check requirements for type
    using type = Scalar;
};

//----------------------------------------------------------------------------------
//                              is_member_assign
//----------------------------------------------------------------------------------
// return true if assignment given by Item defines element at (Row, Col)
template<Integer Row, Integer Col, class Item>
struct is_member_assign
{
    static_assert(md::dependent_false<Item>::value, "invalid item type, this should not happened");
};

template<Integer Row, Integer Col, 
        Integer M1, Integer N1, Integer M2, Integer N2, class Array2, class Colon_1>
struct is_member_assign<Row, Col, virtual_assign_item<M1, N1, M2, N2, Array2, Colon_1>>
{
    static const bool value = is_in_colon<Row, Col, Colon_1, M1, N1>::value;
};

template<Integer Row, Integer Col, Integer M1, Integer N1, class Scalar, class Colon_1>
struct is_member_assign<Row,Col,virtual_assign_item_scalar<M1, N1, Scalar, Colon_1>>
{
    static const bool value = is_in_colon<Row, Col, Colon_1, M1, N1>::value;
};

//----------------------------------------------------------------------------------
//                              is_in_colon
//----------------------------------------------------------------------------------
// return true if element (Row, Col) is one of elements selected by M(Colon_1),
// where M is Mat_Row x Mat_Col matrix
template<Integer Row, Integer Col, class Colon_1, Integer Mat_Row, Integer Mat_Col>
struct is_in_colon
{
    static_assert(md::dependent_false<Colon_1>::value, 
                  "invalid item type, this should not happened");
};

template<Integer Row, Integer Col, Integer Mat_Row, Integer Mat_Col>
struct is_in_colon<Row, Col, colon_all, Mat_Row, Mat_Col>
{
    static_assert(Col <= Mat_Col && Col >= 1 && Row <= Mat_Row && Row >= 1, 
                    "position out of matrix range");

    static const bool value = true;
};

template<Integer Row, Integer Col, Integer Pos, Integer Mat_Row, Integer Mat_Col>
struct is_in_colon<Row, Col, colon<Pos>, Mat_Row, Mat_Col>
{
    static_assert(Col <= Mat_Col && Col >= 1 && Row <= Mat_Row && Row >= 1, 
                    "position out of matrix range");

    // convert Pos to (row, col) pair
    static const Integer pos_col = (Pos - 1) / Mat_Row + 1;
    static const Integer pos_row = Pos - (pos_col - 1) * Mat_Row;

    static const bool value = (pos_col == Col && pos_row == Row);
};

template<Integer Row, Integer Col, Integer Start, Integer End, 
        Integer Mat_Row, Integer Mat_Col>
struct is_in_colon<Row, Col, colon2<Start, End>, Mat_Row, Mat_Col>
{
    static_assert(Col <= Mat_Col && Col >= 1 && Row <= Mat_Row && Row >= 1, 
                  "position out of matrix range");

    // convert (Row, Col) to position
    static const Integer pos = (Col - 1) * Mat_Row + Row;

    static const bool value = (pos >= Start && pos <= End);
};

template<Integer Row, Integer Col, Integer Start, Integer Step, Integer End, 
                    Integer Mat_Row, Integer Mat_Col>
struct is_in_colon<Row, Col, colon3<Start, Step, End>, Mat_Row, Mat_Col>
{
    static_assert(Col <= Mat_Col && Col >= 1 && Row <= Mat_Row && Row >= 1, 
                    "position out of matrix range");

    // convert (Row, Col) to position
    static const Integer pos    = (Col - 1) * Mat_Row + Row;
    static const Integer dif    = pos - Start;


    static const bool value     
        = dif == 0 
        || (dif > 0 && Step > 0 && dif % Step == 0 && pos >= Start && pos <= End)
        || (dif < 0 && Step < 0 && ((-dif) % (-Step) == 0) && pos >= End && pos <= Start);
};

//----------------------------------------------------------------------------------
//                              get_relative_pos
//----------------------------------------------------------------------------------
// return which element selected by Colon points to element at position Pos
template<Integer Pos, class Colon>
struct get_relative_pos
{
    static_assert(md::dependent_false<Colon>::value, 
                  "invalid item type, this should not happened");
};

template<Integer Pos>
struct get_relative_pos<Pos, colon_all>
{
    static const Integer value = Pos;
};

template<Integer Pos, Integer Sel>
struct get_relative_pos<Pos, colon<Sel>>
{
    static_assert(Pos == Sel, "invalid element");
    static const Integer value = 1;
};

template<Integer Pos, Integer Start, Integer End>
struct get_relative_pos<Pos, colon2<Start, End>>
{
    static_assert(Pos >= Start && Pos <= End, "invalid element");
    static const Integer value = Pos - Start + 1;
};

template<Integer Pos, Integer Start, Integer Step, Integer End>
struct get_relative_pos<Pos, colon3<Start, Step, End>>
{
    static_assert(Step > 0 && (Pos >= Start && Pos <= End)
                  || Step < 0 && (Pos >= End && Pos <= Start) , "invalid element");

    static const Integer dif    = Pos - Start;

    static_assert(dif == 0 
                  || dif > 0 && Step > 0 && dif % Step == 0 
                  || dif < 0 && Step < 0 && (-dif) % (-Step) == 0, "invalid element");
   
    static const Integer steps  = dif / Step;
    static const Integer value  = steps + 1;
};

//----------------------------------------------------------------------------------
//                              make_colon_assignment
//----------------------------------------------------------------------------------
//TODO
template<Integer M, Integer N, class Array1, class Deps, class Assign_Info>
struct make_colon_assignment
{
    static_assert(md::dependent_false<Array1>::value, 
                "this type should not be instantiated");
};

template<Integer M, Integer N, class Deps, class Assign_Info>
struct make_colon_assignment2
{
    static_assert(md::dependent_false<Deps>::value, 
                "this type should not be instantiated");
};

template<Integer M, Integer N, class Array1, DPS Deps, Integer M2, Integer N2, class Array2, class Colon_1>
struct make_colon_assignment<M, N, Array1, Deps, virtual_assign_item<M, N, M2, N2, Array2, Colon_1>>
{
    using array_type    = mkd::mat_assign_array_colon<M, N, Array1, Colon_1, M2, N2, Array2>;
    using type          = ct_matrix<M2, N2, array_type, Deps>;
};

template<Integer M, Integer N, DPS Deps, Integer M2, Integer N2, class Array2, class Colon_1>
struct make_colon_assignment2<M, N, Deps, virtual_assign_item<M, N, M2, N2, Array2, Colon_1>>
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
    static_assert(md::dependent_false<Matrix>::value, 
                "this type should not be instantiated");
};

template<Integer M, Integer N, Mat_array Array, DPS Deps>
struct expand_virtual_matrix<ct_matrix<M, N, Array, Deps>>
{
    using type = list::list<ct_matrix<M, N, Array, Deps>>;
};

template<Integer M, Integer N, class Array1, class Virt_Tag, class ... Assign_List, DPS Deps>
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
    static_assert(md::dependent_false<Matrix>::value, 
                "this type should not be instantiated");
};

template<Integer M, Integer N, Mat_array Array, DPS Deps>
struct expand_virtual_matrix2<ct_matrix<M, N, Array, Deps>>
{
    using type = list::list<list::list<colon_all, ct_matrix<M, N, Array, Deps>>>;
};

template<Integer M, Integer N, class Virt_Tag, class ... Assign_List, DPS Deps>
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
