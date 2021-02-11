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
#include "mkgen/details/expressions/expr_plus.h"

namespace matcl { namespace mkgen { namespace details
{

//----------------------------------------------------------------------------------
//                              forward declarations
//----------------------------------------------------------------------------------

template<class Array, Integer Row, Integer Col>
struct mat_scal_plus_array_get_elem;

template<class Array, Integer Row, Integer Col>
struct mat_plus_array_get_elem;

template<class Array, Integer Row, Integer Col>
struct mat_minus_array_get_elem;

template<class Array, Integer Row, Integer Col>
struct mat_scal_minus_array_get_elem;

template<class Array, Integer Row, Integer Col>
struct scal_mat_minus_array_get_elem;

template<class Array, Integer Row, Integer Col>
struct mat_uminus_array_get_elem;

//----------------------------------------------------------------------------------
//                              matrix arrays
//----------------------------------------------------------------------------------
template<Integer M, Integer N, class Array1, class Array2>
struct mat_plus_array : public matrix_array<mat_plus_array<M, N, Array1, Array2>>
{
    using this_type = mat_plus_array<M, N, Array1, Array2>;

    template<Integer Row, Integer Col>
    using get_element_impl  = mat_plus_array_get_elem<this_type, Row, Col>;
};

template<Integer M, Integer N, class Array1, class Array2>
struct mat_scal_plus_array : public matrix_array<mat_scal_plus_array<M, N, Array1, Array2>>
{
    using this_type = mat_scal_plus_array<M, N, Array1, Array2>;

    template<Integer Row, Integer Col>
    using get_element_impl  = mat_scal_plus_array_get_elem<this_type, Row, Col>;
};

template<Integer M, Integer N, class Array1, class Array2>
struct mat_minus_array : public matrix_array<mat_minus_array<M, N, Array1, Array2>>
{
    using this_type = mat_minus_array<M, N, Array1, Array2>;

    template<Integer Row, Integer Col>
    using get_element_impl  = mat_minus_array_get_elem<this_type, Row, Col>;
};

template<Integer M, Integer N, class Array1, class Array2>
struct mat_scal_minus_array : public matrix_array<mat_scal_minus_array<M, N, Array1, Array2>>
{
    using this_type = mat_scal_minus_array<M, N, Array1, Array2>;

    template<Integer Row, Integer Col>
    using get_element_impl  = mat_scal_minus_array_get_elem<this_type, Row, Col>;
};

template<Integer M, Integer N, class Array1, class Array2>
struct scal_mat_minus_array : public matrix_array<scal_mat_minus_array<M, N, Array1, Array2>>
{
    using this_type = scal_mat_minus_array<M, N, Array1, Array2>;

    template<Integer Row, Integer Col>
    using get_element_impl  = scal_mat_minus_array_get_elem<this_type, Row, Col>;
};

template<Integer M, Integer N, class Array>
struct mat_uminus_array : public matrix_array<mat_uminus_array<M, N, Array>>
{
    using this_type = mat_uminus_array<M, N, Array>;

    template<Integer Row, Integer Col>
    using get_element_impl  = mat_uminus_array_get_elem<this_type, Row, Col>;
};

//----------------------------------------------------------------------------------
//                              get_elem impl
//----------------------------------------------------------------------------------
template<class Array, Integer Row, Integer Col>
struct mat_scal_plus_array_get_elem
{};

template<Integer M, Integer N, class Array1, class Array2, Integer Row, Integer Col>
struct mat_scal_plus_array_get_elem<mkd::mat_scal_plus_array<M, N, Array1, Array2>, Row, Col>
{
    using elem_1    = typename Array1 :: template get_element<Row, Col>::type;
    using elem_2    = Array2;

    using type      = typename mkd::make_plus_root<elem_1,elem_2>::type;
};

template<class Array, Integer Row, Integer Col>
struct mat_plus_array_get_elem
{};

template<Integer M, Integer N, class Array1, class Array2, Integer Row, Integer Col>
struct mat_plus_array_get_elem<mkd::mat_plus_array<M,N,Array1,Array2>, Row, Col>
{
    using elem_1    = typename Array1 :: template get_element<Row, Col>::type;
    using elem_2    = typename Array2 :: template get_element<Row, Col>::type;

    using type      = typename mkd::make_plus_root<elem_1,elem_2>::type;
};

template<class Array, Integer Row, Integer Col>
struct mat_minus_array_get_elem
{};

template<Integer M, Integer N, class Array1, class Array2, Integer Row, Integer Col>
struct mat_minus_array_get_elem<mkd::mat_minus_array<M,N,Array1,Array2>, Row, Col>
{
    using elem_1    = typename Array1 :: template get_element<Row, Col>::type;
    using elem_2    = typename Array2 :: template get_element<Row, Col>::type;
    
    using type      = typename mkd::make_minus_root<elem_1, elem_2>::type;
};

template<class Array, Integer Row, Integer Col>
struct mat_scal_minus_array_get_elem
{};

template<Integer M, Integer N, class Array1, class Array2, Integer Row, Integer Col>
struct mat_scal_minus_array_get_elem<mkd::mat_scal_minus_array<M, N, Array1, Array2>, Row, Col>
{
    using elem_1    = typename Array1 :: template get_element<Row, Col>::type;
    using elem_2    = Array2;
    
    using type      = typename mkd::make_minus_root<elem_1,elem_2>::type;
};

template<class Array, Integer Row, Integer Col>
struct scal_mat_minus_array_get_elem
{};

template<Integer M, Integer N, class Array1, class Array2, Integer Row, Integer Col>
struct scal_mat_minus_array_get_elem<mkd::scal_mat_minus_array<M,N,Array1, Array2>, Row, Col>
{
    using elem_1    = typename Array1 :: template get_element<Row, Col>::type;
    using elem_2    = Array2;

    using type      = typename mkd::make_minus_root<elem_2,elem_1>::type;
};

template<class Array, Integer Row, Integer Col>
struct mat_uminus_array_get_elem
{};

template<Integer M, Integer N, class Array, Integer Row, Integer Col>
struct mat_uminus_array_get_elem<mkd::mat_uminus_array<M,N,Array>, Row, Col>
{
    using elem      = typename Array::template get_element<Col, Row>::type;
    using type      = typename make_uminus_root<elem>::type;
};

//----------------------------------------------------------------------------------
//                              mat_plus_impl
//----------------------------------------------------------------------------------
template<class M1, class M2>
struct mat_plus_impl
{
    static_assert(md::dependent_false<M1>::value, "M1, M2 must be ct_matrix or ct_scalar");
};

template<Integer M1, Integer N1, Mat_array Array1, DPS Deps1, 
         Integer M2, Integer N2, Mat_array Array2, DPS Deps2>
struct mat_plus_impl<ct_matrix<M1, N1, Array1, Deps1>, ct_matrix<M2, N2, Array2, Deps2>>
{
    static_assert(M1 == M2, "conconformant matrix sizes");
    static_assert(N1 == N2, "conconformant matrix sizes");

    static const Integer M1_M2  = M1;
    static const Integer N1_N2  = N1;

    using array_type    = mkd::mat_plus_array<M1_M2, N1_N2, Array1, Array2>;
    using deps          = typename link_deps<Deps1, Deps2>::type;
    using type          = ct_matrix<M1_M2, N1_N2, array_type,deps>;
};

template<Integer M1, Integer N1, Mat_array Array1, DPS Deps1, 
        Scal_data Array2, DPS Deps2>
struct mat_plus_impl<ct_matrix<M1, N1, Array1, Deps1>, ct_scalar<Array2,Deps2>>
{
    using array_type    = mkd::mat_scal_plus_array<M1, N1, Array1, Array2>;
    using deps          = typename link_deps<Deps1, Deps2>::type;
    using type          = ct_matrix<M1, N1, array_type,deps>;
};

template<Integer M1, Integer N1, Mat_array Array1, DPS Deps1, 
        Scal_data Array2, DPS Deps2>
struct mat_plus_impl<ct_scalar<Array2,Deps2>, ct_matrix<M1, N1, Array1, Deps1>>
{
    using array_type    = mkd::mat_scal_plus_array<M1, N1, Array1, Array2>;
    using deps          = typename link_deps<Deps1, Deps2>::type;
    using type          = ct_matrix<M1, N1, array_type,deps>;
};

template<Scal_data Array1, DPS Deps1, Scal_data Array2, DPS Deps2>
struct mat_plus_impl<ct_scalar<Array1, Deps1>, ct_scalar<Array2, Deps2>>
{
    using scal_1    = Array1;
    using scal_2    = Array2;

    using plus_type = typename mkd::make_plus_root<scal_1, scal_2>::type;
    using deps      = typename link_deps<Deps1, Deps2>::type;
    
    using type      = ct_scalar<plus_type, deps>;
};

//----------------------------------------------------------------------------------
//                              mat_minus_impl
//----------------------------------------------------------------------------------
template<class M1, class M2>
struct mat_minus_impl
{
    static_assert(md::dependent_false<M1>::value, "M1, M2 must be ct_matrix or ct_scalar");
};

template<Integer M1, Integer N1, Mat_array Array1, DPS Deps1, 
         Integer M2, Integer N2, Mat_array Array2, DPS Deps2>
struct mat_minus_impl<ct_matrix<M1, N1, Array1, Deps1>, ct_matrix<M2, N2, Array2, Deps2>>
{
    static_assert(M1 == M2, "conconformant matrix sizes");
    static_assert(N1 == N2, "conconformant matrix sizes");

    static const Integer M1_M2  = M1;
    static const Integer N1_N2  = N1;

    using array_type    = mkd::mat_minus_array<M1_M2, N1_N2, Array1, Array2>;
    using deps          = typename link_deps<Deps1, Deps2>::type;
    using type          = ct_matrix<M1_M2, N1_N2, array_type,deps>;
};

template<Integer M1, Integer N1, Mat_array Array1, DPS Deps1, 
        Scal_data Array2, DPS Deps2>
struct mat_minus_impl<ct_matrix<M1, N1, Array1, Deps1>, ct_scalar<Array2, Deps2>>
{
    using array_type    = mkd::mat_scal_minus_array<M1, N1, Array1, Array2>;
    using deps          = typename link_deps<Deps1, Deps2>::type;
    using type          = ct_matrix<M1, N1, array_type,deps>;
};

template<Integer M1, Integer N1, Mat_array Array1, DPS Deps1, 
        Scal_data Array2, DPS Deps2>
struct mat_minus_impl<ct_scalar<Array2, Deps2>, ct_matrix<M1, N1, Array1, Deps1>>
{
    using array_type    = mkd::scal_mat_minus_array<M1, N1, Array1, Array2>;
    using deps          = typename link_deps<Deps1, Deps2>::type;
    using type          = ct_matrix<M1, N1, array_type,deps>;
};

template<Scal_data Array1, DPS Deps1, Scal_data Array2, DPS Deps2>
struct mat_minus_impl<ct_scalar<Array1,Deps1>, ct_scalar<Array2,Deps2>>
{
    using scal_1        = Array1;
    using scal_2        = Array2;

    using minus_type    = typename mkd::make_minus_root<scal_1,scal_2>::type;
    using deps          = typename link_deps<Deps1, Deps2>::type;
    
    using type          = ct_scalar<minus_type, deps>;
};

//----------------------------------------------------------------------------------
//                              mat_uminus_impl
//----------------------------------------------------------------------------------
template<Mat_or_scalar M1>
struct mat_uminus_impl
{
    static_assert(md::dependent_false<M1>::value, "M1 must be ct_matrix or ct_scalar");
};

template<Integer M, Integer N, Mat_array Array, DPS Deps1>
struct mat_uminus_impl<ct_matrix<M, N, Array, Deps1>>
{
    using array_type    = mkd::mat_uminus_array<M, N, Array>;
    using type          = ct_matrix<N, M, array_type,Deps1>;
};

template<Scal_data Array, DPS Deps>
struct mat_uminus_impl<ct_scalar<Array, Deps>>
{
    using scal_1        = Array;
    using uminus_type   = typename make_uminus_root<scal_1>::type;
    using type          = ct_scalar<uminus_type, Deps>;
};

}}}
