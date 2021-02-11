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

#include "mkgen/expression/expressions.h"
#include "mkgen/details/expressions/expr_mult.h"
#include "mkgen/details/expressions/expr_dot_scalar_data.h"

namespace matcl { namespace mkgen { namespace details
{

//----------------------------------------------------------------------------------
//                              forward declarations
//----------------------------------------------------------------------------------

// get element (Row, Col) from mat_mult_array
template<class Array, Integer Row, Integer Col>
struct mat_mult_array_get_elem;

// get element (Row, Col) from mat_scal_mult_array
template<class Array, Integer Row, Integer Col>
struct mat_scal_mult_array_get_elem;

template<class Array, Integer Row, Integer Col>
struct mult_rows_array_get_elem;

template<class Array, Integer Row, Integer Col>
struct mult_cols_array_get_elem;

template<class Array, Integer Row, Integer Col>
struct mult_array_get_elem;

template<class Array, Integer Row, Integer Col>
struct div_array_get_elem;

template<class Array, Integer Row, Integer Col>
struct div_array_mat_scal_get_elem;

template<class Array, Integer Row, Integer Col>
struct div_array_scal_mat_get_elem;

// make element (Row, Col) from multiplication of matrices Array_1 x Array_2
// i.e. K-element dot product
template<Integer Row, Integer Col, Integer K_start, Integer K_end, Integer Length, 
    class Array_1, class Array_2>
struct make_array_mat_mult;

template<class List_1, class List_2>
struct merge_dots;

//----------------------------------------------------------------------------------
//                              matrix arrays
//----------------------------------------------------------------------------------
template<Integer K, class Array1, class Array2>
struct mat_mult_array : public matrix_array<mat_mult_array<K, Array1, Array2>>
{
    using this_type = mat_mult_array<K, Array1, Array2>;

    template<Integer Row, Integer Col>
    using get_element_impl  = mat_mult_array_get_elem<this_type, Row, Col>;
};

template<class Array1, class Scalar2>
struct mat_scal_mult_array : public matrix_array<mat_scal_mult_array<Array1, Scalar2>>
{
    using this_type = mat_scal_mult_array<Array1, Scalar2>;

    template<Integer Row, Integer Col>
    using get_element_impl  = mat_scal_mult_array_get_elem<this_type, Row, Col>;
};

template<class Array1, class Array2>
struct mult_rows_array : public matrix_array<mult_rows_array<Array1, Array2>>
{
    using this_type = mult_rows_array<Array1, Array2>;

    template<Integer Row, Integer Col>
    using get_element_impl  = mult_rows_array_get_elem<this_type, Row, Col>;
};

template<class Array1, class Array2>
struct mult_cols_array : public matrix_array<mult_cols_array<Array1, Array2>>
{
    using this_type = mult_cols_array<Array1, Array2>;

    template<Integer Row, Integer Col>
    using get_element_impl  = mult_cols_array_get_elem<this_type, Row, Col>;
};

template<class Array1, class Array2>
struct mult_array : public matrix_array<mult_array<Array1, Array2>>
{
    using this_type = mult_array<Array1, Array2>;

    template<Integer Row, Integer Col>
    using get_element_impl  = mult_array_get_elem<this_type, Row, Col>;
};

template<class Array1, class Array2>
struct div_array : public matrix_array<div_array<Array1, Array2>>
{
    using this_type = div_array<Array1, Array2>;

    template<Integer Row, Integer Col>
    using get_element_impl  = div_array_get_elem<this_type, Row, Col>;
};

template<class Array1, class Scal2>
struct div_array_mat_scal : public matrix_array<div_array_mat_scal<Array1, Scal2>>
{
    using this_type = div_array_mat_scal<Array1, Scal2>;

    template<Integer Row, Integer Col>
    using get_element_impl  = div_array_mat_scal_get_elem<this_type, Row, Col>;
};

template<class Array2, class Scal1>
struct div_array_scal_mat : public matrix_array<div_array_scal_mat<Array2, Scal1>>
{
    using this_type = div_array_scal_mat<Array2, Scal1>;

    template<Integer Row, Integer Col>
    using get_element_impl  = div_array_scal_mat_get_elem<this_type, Row, Col>;
};

//----------------------------------------------------------------------------------
//                              get_elem impl
//----------------------------------------------------------------------------------

template<class Array, Integer Row, Integer Col>
struct mat_mult_array_get_elem
{};

template<Integer K, class Array1, class Array2, Integer Row, Integer Col>
struct mat_mult_array_get_elem<mkd::mat_mult_array<K, Array1, Array2>, Row, Col>
{
    using type = typename make_array_mat_mult<Row,Col,1,K,K,Array1,Array2>::type;
};

template<class Array1, class Array2, Integer Row, Integer Col>
struct mat_mult_array_get_elem<mkd::mat_mult_array<0, Array1, Array2>, Row, Col>
{
    using type  = zero_sd;
};

template<class Array1, class Array2, Integer Row, Integer Col>
struct mat_mult_array_get_elem<mkd::mat_mult_array<1, Array1, Array2>, Row, Col>
{
    using elem_1 = typename Array1 :: template get_element<Row, 1>::type;
    using elem_2 = typename Array2 :: template get_element<1, Col>::type;

    //TODO: why element_step?
    using elem_s = mkd::element_step<elem_2,0>;

    using type  = typename mkd::make_mult_root<elem_1,elem_s>::type;
};

template<class Array, Integer Row, Integer Col>
struct mat_scal_mult_array_get_elem
{};

template<class Array1, class Array2, Integer Row, Integer Col>
struct mat_scal_mult_array_get_elem<mkd::mat_scal_mult_array<Array1, Array2>, Row, Col>
{
    using elem_1    = typename Array1 :: template get_element<Row,Col>::type;
    using elem_2    = Array2;

    using type      = typename mkd::make_mult_root<elem_1,elem_2>::type;
};

template<class Array, Integer Row, Integer Col>
struct mult_rows_array_get_elem
{};

template<class Array1, class Array2, Integer Row, Integer Col>
struct mult_rows_array_get_elem<mkd::mult_rows_array<Array1,Array2>, Row, Col>
{
    using elem_1    = typename Array1 :: template get_element<Row, Col>::type;
    using elem_2    = typename Array2 :: template get_element<Row, 1>::type;

    using type      = typename mkd::make_mult_root<elem_1,elem_2>::type;
};

template<class Array, Integer Row, Integer Col>
struct mult_cols_array_get_elem
{};

template<class Array1, class Array2, Integer Row, Integer Col>
struct mult_cols_array_get_elem<mkd::mult_cols_array<Array1,Array2>, Row, Col>
{
    using elem_1    = typename Array1 :: template get_element<Row, Col>::type;
    using elem_2    = typename Array2 :: template get_element<1, Col>::type;

    using type      = typename mkd::make_mult_root<elem_1,elem_2>::type;
};

template<class Array, Integer Row, Integer Col>
struct mult_array_get_elem
{};

template<class Array1, class Array2, Integer Row, Integer Col>
struct mult_array_get_elem<mkd::mult_array<Array1,Array2>, Row, Col>
{
    using elem_1    = typename Array1 :: template get_element<Row, Col>::type;
    using elem_2    = typename Array2 :: template get_elemend<Row, Col>::type;

    using type      = typename mkd::make_mult_root<elem_1,elem_2>::type;
};

template<class Array, Integer Row, Integer Col>
struct div_array_get_elem
{};

template<class Array1, class Array2, Integer Row, Integer Col>
struct div_array_get_elem<mkd::div_array<Array1, Array2>, Row, Col>
{
    using elem_1    = typename Array1 :: template get_element<Row, Col>::type;
    using elem_2    = typename Array2 :: template get_element<Row, Col>::type;

    using type      = typename mkd::make_div_root<elem_1,elem_2>::type;
};

template<class Array, Integer Row, Integer Col>
struct div_array_mat_scal_get_elem
{};

template<class Array1, class Scal2, Integer Row, Integer Col>
struct div_array_mat_scal_get_elem<mkd::div_array_mat_scal<Array1, Scal2>, Row, Col>
{
    using elem_1    = typename Array1 :: template get_element<Row, Col>::type;
    using elem_2    = Scal2;

    using type      = typename mkd::make_div_root<elem_1,elem_2>::type;
};

template<class Array, Integer Row, Integer Col>
struct div_array_scal_mat_get_elem
{};

template<class Array2, class Scal1, Integer Row, Integer Col>
struct div_array_scal_mat_get_elem<mkd::div_array_scal_mat<Array2, Scal1>, Row, Col>
{
    using elem_1    = Scal1;
    using elem_2    = typename Array2 :: template get_element<Row, Col>::type;    
    
    using type      = typename mkd::make_div_root<elem_1,elem_2>::type;
};

//----------------------------------------------------------------------------------
//                              mat_mult_impl
//----------------------------------------------------------------------------------
template<class M1, class M2>
struct mat_mult_impl
{
    static_assert(md::dependent_false<M1>::value, "M1, M2 must be ct_matrix or ct_scalar");
};

template<Integer M1, Integer N1, Mat_array Array1, DPS Deps1,
        Integer M2, Integer N2, Mat_array Array2, DPS Deps2>
struct mat_mult_impl<ct_matrix<M1, N1, Array1, Deps1>, ct_matrix<M2, N2, Array2, Deps2>>
{
    static_assert(N1 == M2, "conconformant matrix sizes");

    static const Integer N1_M2  = N1;

    using array_type    = mkd::mat_mult_array<N1_M2, Array1, Array2>;
    using deps          = typename link_deps<Deps1, Deps2>::type;
    using type          = ct_matrix<M1, N2, array_type, deps>;
};

template<Integer M1, Integer N1, Mat_array Array1, DPS Deps1, 
         Scal_data Array2, DPS Deps2>
struct mat_mult_impl<ct_matrix<M1, N1, Array1, Deps1>, ct_scalar<Array2,Deps2> >
{
    using array_type    = mkd::mat_scal_mult_array<Array1, Array2>;
    using deps          = typename link_deps<Deps1, Deps2>::type;
    using type          = ct_matrix<M1, N1, array_type, deps>;
};

template<Scal_data Array1, DPS Deps1, 
         Integer M2, Integer N2, Mat_array Array2, DPS Deps2>
struct mat_mult_impl<ct_scalar<Array1,Deps1>, ct_matrix<M2, N2, Array2, Deps2>>
{
    using array_type    = mkd::mat_scal_mult_array<Array2, Array1>;
    using deps          = typename link_deps<Deps1, Deps2>::type;
    using type          = ct_matrix<M2, N2, array_type,deps>;
};

template<Scal_data Array1, DPS Deps1, Scal_data Array2, DPS Deps2>
struct mat_mult_impl<ct_scalar<Array1, Deps1>, ct_scalar<Array2, Deps2>>
{
    using scal_1    = Array1;
    using scal_2    = Array2;

    using mult_type = typename mkd::make_mult_root<scal_1, scal_2>::type;
    using deps      = typename link_deps<Deps1, Deps2>::type;
    using type      = ct_scalar<mult_type, deps>;
};

//----------------------------------------------------------------------------------
//                              mult_rows_impl
//----------------------------------------------------------------------------------
template<class M1, class M2>
struct mult_rows_impl
{
    static_assert(md::dependent_false<M1>::value, "M1, M2 must be ct_matrix or ct_scalar");
};

template<Integer M1, Integer N1, Mat_array Array1, DPS Deps1, 
         Integer M2, Integer N2, Mat_array Array2, DPS Deps2>
struct mult_rows_impl<ct_matrix<M1, N1, Array1, Deps1>, ct_matrix<M2,N2,Array2,Deps2>>
{
    static_assert(M1 == M2, "conconformant matrix sizes");
    static_assert(N2 == 1, "vector expected");

    static const Integer M1_M2  = M1;

    using array_type    = mkd::mult_rows_array<Array1, Array2>;
    using deps          = typename link_deps<Deps1, Deps2>::type;
    using type          = ct_matrix<M1_M2, N1, array_type,deps>;
};

template<Integer M1, Integer N1, Mat_array Array1, DPS Deps1, 
         Scal_data Array2, DPS Deps2>
struct mult_rows_impl<ct_matrix<M1, N1, Array1, Deps1>, ct_scalar<Array2,Deps2> >
{
    static_assert(M1 == 1, "vector expected");

    using A             = ct_matrix<M1, N1, Array1, Deps1>;
    using B             = ct_scalar<Array2,Deps2>;

    // equivalent to A * B
    using type          = typename mat_mult_impl<A, B>::type;
};

template<Scal_data Array1, DPS Deps1, 
         Integer M2, Integer N2, Mat_array Array2, DPS Deps2>
struct mult_rows_impl<ct_scalar<Array1,Deps1>, ct_matrix<M2, N2, Array2, Deps2>>
{
    static_assert(M2 == 1 && N2 == 1, "1x1 matrix expected");

    using A             = ct_scalar<Array1,Deps1>;
    using B             = ct_matrix<M2, N2, Array2, Deps2>;

    // equivalent to A * B
    using type          = typename mat_mult_impl<A, B>::type;
};

template<Scal_data Array1, DPS Deps1, Scal_data Array2, DPS Deps2>
struct mult_rows_impl<ct_scalar<Array1, Deps1>, ct_scalar<Array2, Deps2>>
{
    using scal_1    = Array1;
    using scal_2    = Array2;

    using mult_type = typename mkd::make_mult_root<scal_1, scal_2>::type;
    using deps      = typename link_deps<Deps1, Deps2>::type;
    using type      = ct_scalar<mult_type, deps>;
};

//----------------------------------------------------------------------------------
//                              mult_cols_impl
//----------------------------------------------------------------------------------
template<class M1, class M2>
struct mult_cols_impl
{
    static_assert(md::dependent_false<M1>::value, "M1, M2 must be ct_matrix or ct_scalar");
};

template<Integer M1, Integer N1, Mat_array Array1, DPS Deps1, 
         Integer M2, Integer N2, Mat_array Array2, DPS Deps2>
struct mult_cols_impl<ct_matrix<M1, N1, Array1, Deps1>, ct_matrix<M2, N2, Array2, Deps2>>
{
    static_assert(N1 == M2, "conconformant matrix sizes");
    static_assert(N2 == 1, "vector expected");

    static const Integer N1_M2  = N1;

    using array_type    = mkd::mult_cols_array<Array1, Array2>;
    using deps          = typename link_deps<Deps1, Deps2>::type;
    using type          = ct_matrix<M1, N1_M2, array_type,deps>;
};

template<Integer M1, Integer N1, Mat_array Array1, DPS Deps1, 
         Scal_data Array2, DPS Deps2>
struct mult_cols_impl<ct_matrix<M1, N1, Array1, Deps1>, ct_scalar<Array2,Deps2> >
{
    static_assert(N1 == 1, "vector expected");

    using A             = ct_matrix<M1, N1, Array1, Deps1>;
    using B             = ct_scalar<Array2,Deps2>;

    // equivalent to A * B
    using type          = typename mat_mult_impl<A, B>::type;
};

template<Scal_data Array1, DPS Deps1, 
         Integer M2, Integer N2, Mat_array Array2, DPS Deps2>
struct mult_cols_impl<ct_scalar<Array1, Deps1>, ct_matrix<M2, N2, Array2, Deps2>>
{
    static_assert(M2 == 1 && N2 == 1, "1x1 matrix expected");

    using A             = ct_scalar<Array1,Deps1>;
    using B             = ct_matrix<M2, N2, Array2, Deps2>;

    // equivalent to A * B
    using type          = typename mat_mult_impl<A, B>::type;
};

template<Scal_data Array1, DPS Deps1, Scal_data Array2, DPS Deps2>
struct mult_cols_impl<ct_scalar<Array1, Deps1>, ct_scalar<Array2, Deps2>>
{
    using scal_1    = Array1;
    using scal_2    = Array2;

    using mult_type = typename mkd::make_mult_root<scal_1, scal_2>::type;
    using deps      = typename link_deps<Deps1, Deps2>::type;
    using type      = ct_scalar<mult_type, deps>;
};

//----------------------------------------------------------------------------------
//                              mul_impl
//----------------------------------------------------------------------------------
template<class M1, class M2>
struct mul_impl
{
    static_assert(md::dependent_false<M1>::value, "M1, M2 must be ct_matrix or ct_scalar");
};

template<Integer M1, Integer N1, Mat_array Array1, DPS Deps1,
         Integer M2, Integer N2, Mat_array Array2, DPS Deps2>
struct mul_impl<ct_matrix<M1, N1, Array1, Deps1>, ct_matrix<M2, N2, Array2, Deps2>>
{
    static_assert(M1 == M2, "conconformant matrix sizes");
    static_assert(N1 == N2, "conconformant matrix sizes");

    static const Integer M1_M2  = M1;
    static const Integer N1_N2  = N1;

    using array_type    = mkd::mult_array<Array1, Array2>;
    using deps          = typename link_deps<Deps1, Deps2>::type;
    using type          = ct_matrix<M1_M2, N1_N2, array_type,deps>;
};

template<Integer M1, Integer N1, Mat_array Array1, DPS Deps1, 
         Scal_data Array2, DPS Deps2>
struct mul_impl<ct_matrix<M1, N1, Array1, Deps1>, ct_scalar<Array2,Deps2> >
{
    using A             = ct_matrix<M1, N1, Array1, Deps1>;
    using B             = ct_scalar<Array2,Deps2>;

    // equivalent to A * B
    using type          = typename mat_mult_impl<A, B>::type;
};

template<Scal_data Array1, DPS Deps1, 
         Integer M2, Integer N2, Mat_array Array2, DPS Deps2>
struct mul_impl<ct_scalar<Array1,Deps1>, ct_matrix<M2, N2, Array2, Deps2>>
{
    using A             = ct_scalar<Array1,Deps1>;
    using B             = ct_matrix<M2, N2, Array2, Deps2>;

    // equivalent to A * B
    using type          = typename mat_mult_impl<A, B>::type;
};

template<Scal_data Array1, DPS Deps1, Scal_data Array2, DPS Deps2>
struct mul_impl<ct_scalar<Array1, Deps1>, ct_scalar<Array2, Deps2>>
{
    using scal_1    = Array1;
    using scal_2    = Array2;

    using mult_type = typename mkd::make_mult_root<scal_1, scal_2>::type;
    using deps      = typename link_deps<Deps1, Deps2>::type;
    using type      = ct_scalar<mult_type, deps>;
};

//----------------------------------------------------------------------------------
//                              mat_div_impl
//----------------------------------------------------------------------------------
template<class M1, class M2>
struct mat_div_impl
{
    static_assert(md::dependent_false<M1>::value, "M1, M2 must be ct_matrix or ct_scalar");
};

template<Integer M1, Integer N1, Mat_array Array1, DPS Deps1, 
         Integer M2, Integer N2, Mat_array Array2, DPS Deps2>
struct mat_div_impl<ct_matrix<M1, N1, Array1, Deps1>, ct_matrix<M2, N2, Array2, Deps2>>
{
    static_assert(M1 == M2, "conconformant matrix sizes");
    static_assert(N1 == N2, "conconformant matrix sizes");

    static const Integer M1_M2  = M1;
    static const Integer N1_N2  = N1;

    using array_type    = mkd::div_array<Array1, Array2>;
    using deps          = typename link_deps<Deps1, Deps2>::type;
    using type          = ct_matrix<M1_M2, N1_N2, array_type,deps>;
};

template<Integer M1, Integer N1, Mat_array Array1, DPS Deps1, 
        Scal_data Array2, DPS Deps2>
struct mat_div_impl<ct_matrix<M1, N1, Array1, Deps1>, ct_scalar<Array2,Deps2>>
{
    using array_type    = mkd::div_array_mat_scal<Array1, Array2>;
    using deps          = typename link_deps<Deps1, Deps2>::type;
    using type          = ct_matrix<M1, N1, array_type,deps>;
};

template<Scal_data Array1, DPS Deps1, 
        Integer M2, Integer N2, Mat_array Array2, DPS Deps2>
struct mat_div_impl<ct_scalar<Array1,Deps1>, ct_matrix<M2, N2, Array2, Deps2>>
{
    using array_type    = mkd::div_array_scal_mat<Array2, Array1>;
    using deps          = typename link_deps<Deps1, Deps2>::type;
    using type          = ct_matrix<M2, N2, array_type,deps>;
};

template<Scal_data Array1, DPS Deps1, Scal_data Array2, DPS Deps2>
struct mat_div_impl<ct_scalar<Array1, Deps1>, ct_scalar<Array2, Deps2>>
{
    using div_type      = typename mkd::make_div_root<Array1, Array2>::type;
    using deps          = typename link_deps<Deps1, Deps2>::type;
    using type          = ct_scalar<div_type,deps>;
};

//----------------------------------------------------------------------------------
//                              make_array_mat_mult
//----------------------------------------------------------------------------------
template<Integer Row, Integer Col, Integer K_start, Integer K_end, Integer Length, 
    class Array_1, class Array_2>
struct make_array_mat_mult 
{
    static const Integer K_end_1    = K_start + Length/2 - 1;

    using list_1    = typename make_array_mat_mult<Row,Col,K_start, K_end_1, Length/2, 
                                Array_1, Array_2>::type;
    using list_2    = typename make_array_mat_mult<Row,Col,K_end_1+1, K_end, Length - Length/2, 
                                Array_1, Array_2>::type;

    using type      = typename merge_dots<list_1, list_2>::type;
};

template<Integer Row, Integer Col, Integer K_start, Integer K_end, class Array_1, class Array_2>
struct make_array_mat_mult<Row,Col,K_start,K_end,0,Array_1, Array_2>
{
    using type = zero_sd;
};

template<Integer Row, Integer Col, Integer K_start, Integer K_end, class Array_1, class Array_2>
struct make_array_mat_mult<Row,Col,K_start,K_end,1,Array_1, Array_2>
{
    using elem_1 = typename Array_1 :: template get_element<Row, K_start>::type;
    using elem_2 = typename Array_2 :: template get_element<K_start, Col>::type;
    using type   = expr_dot_scalar_data<list::list<elem_1>, list::list<elem_2>>;
};

template<Integer Row, Integer Col, Integer K_start, Integer K_end, class Array_1, class Array_2>
struct make_array_mat_mult<Row,Col,K_start,K_end,2,Array_1, Array_2>
{
    using elem_1_1      = typename Array_1 :: template get_element<Row, K_start>::type;
    using elem_2_1      = typename Array_2 :: template get_element<K_start, Col>::type;

    using elem_1_2      = typename Array_1 :: template get_element<Row, K_start+1>::type;
    using elem_2_2      = typename Array_2 :: template get_element<K_start+1, Col>::type;

    using type          = expr_dot_scalar_data<list::list<elem_1_1,elem_1_2>, list::list<elem_2_1,elem_2_2>>;
};

template<Integer Row, Integer Col, Integer K_start, Integer K_end, class Array_1, class Array_2>
struct make_array_mat_mult<Row,Col,K_start,K_end,3,Array_1, Array_2>
{
    using elem_1_1      = typename Array_1 :: template get_element<Row, K_start>::type;
    using elem_2_1      = typename Array_2 :: template get_element<K_start, Col>::type;

    using elem_1_2      = typename Array_1 :: template get_element<Row, K_start + 1>::type;
    using elem_2_2      = typename Array_2 :: template get_element<K_start + 1, Col>::type;

    using elem_1_3      = typename Array_1 :: template get_element<Row, K_start+2>::type;
    using elem_2_3      = typename Array_2 :: template get_element<K_start+2, Col>::type;

    using type          = expr_dot_scalar_data<list::list<elem_1_1,elem_1_2, elem_1_3>, 
                                   list::list<elem_2_1,elem_2_2, elem_2_3>>;
};

template<Integer Row, Integer Col, Integer K_start, Integer K_end, class Array_1, class Array_2>
struct make_array_mat_mult<Row,Col,K_start,K_end,4,Array_1, Array_2>
{
    using elem_1_1      = typename Array_1 :: template get_element<Row, K_start>::type;
    using elem_2_1      = typename Array_2 :: template get_element<K_start, Col>::type;

    using elem_1_2      = typename Array_1 :: template get_element<Row, K_start + 1>::type;
    using elem_2_2      = typename Array_2 :: template get_element<K_start + 1, Col>::type;

    using elem_1_3      = typename Array_1 :: template get_element<Row, K_start + 2>::type;
    using elem_2_3      = typename Array_2 :: template get_element<K_start + 2, Col>::type;

    using elem_1_4      = typename Array_1 :: template get_element<Row, K_start+3>::type;
    using elem_2_4      = typename Array_2 :: template get_element<K_start+3, Col>::type;

    using type          = expr_dot_scalar_data<list::list<elem_1_1,elem_1_2, elem_1_3, elem_1_4>, 
                                   list::list<elem_2_1,elem_2_2, elem_2_3, elem_2_4>>;
};

//----------------------------------------------------------------------------------
//                              struct merge_dots
//----------------------------------------------------------------------------------
template<class List_1, class List_2>
struct merge_dots
{
    static_assert(md::dependent_false<List_1>::value,
                  "this type should not be instantiated");
};

template<class ... Arg_11, class ... Arg_12, class ... Arg_21, class ... Arg_22>
struct merge_dots<expr_dot_scalar_data<list::list<Arg_11...>, list::list<Arg_12...>>,
                  expr_dot_scalar_data<list::list<Arg_21...>, list::list<Arg_22...>>>
{
    using type = expr_dot_scalar_data<list::list<Arg_11..., Arg_21...>,
                          list::list<Arg_12..., Arg_22...>>;
};

}}}
