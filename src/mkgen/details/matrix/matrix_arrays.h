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
#include "mkgen/details/utils/mpl.h"
#include "mkgen/details/matrix/element.h"
#include "mkgen/details/matrix/scalar_data.h"
#include "mkgen/matrix/scalar.h"

namespace matcl { namespace mkgen { namespace details
{

//-----------------------------------------------------------------------
//                      forward declarations
//-----------------------------------------------------------------------
 template<class Array, Integer Row, Integer Col>
struct error_access_to_empty_array;

template<class Data, Integer Row, Integer Col>
struct matrix_array_get_elem;

//-----------------------------------------------------------------------
//                      matrix_array impl
//-----------------------------------------------------------------------
template<class Data, Integer Row, Integer Col>
struct matrix_array_get_elem
{
    using type  = typename Data :: template get_element_impl<Row, Col>::type;

    // type must be scalar_data
    static_assert(Scal_data<type>, "Data must be a valid scalar_data");
};

//-----------------------------------------------------------------------
//                      arrays
//-----------------------------------------------------------------------

template<class Array_t, Integer Offset, Integer Step>
struct sub_array_1 : public matrix_array<sub_array_1<Array_t, Offset, Step>>
{
    using this_type = sub_array_1<Array_t, Offset, Step>;

    template<Integer Row, Integer Col>
    using get_element_impl  = sub_array_1_get_elem<this_type, Row, Col>;
};

template<class Array_t, Integer Offset1, Integer Offset2, Integer Step1, Integer Step2>
struct sub_array_2 : public matrix_array<sub_array_2<Array_t, Offset1, Offset2, Step1, Step2>>
{
    using this_type = sub_array_2<Array_t, Offset1, Offset2, Step1, Step2>;

    template<Integer Row, Integer Col>
    using get_element_impl  = sub_array_2_get_elem<this_type, Row, Col>;
};

// value_mat array
template<Tag_matrix_data Tag, Value Val_t, Integer Row, Integer Col>
struct ma_value_get_elem_scalar
        : mk::scal_data_value_tag<ma_value_get_elem_scalar<Tag, Val_t, Row, Col>>
{
    template<class Val_loc>
    static Val_loc value()
    {
        return Val_loc(Tag::value<Val_t, Row, Col>());
    };
};

template<Tag_matrix_data Tag, Value Val_t, Integer Row, Integer Col>
struct matrix_array_value_get_elem
{
    using type      = scal_data_value<ma_value_get_elem_scalar
                            <Tag, Val_t, Row, Col>, Val_t>;

    template<class Val_loc>
    static Val_loc value()
    {
        return Val_loc(Tag::value<Val_t, Row, Col>());
    };
};

template<Tag_matrix_data Tag, Value Val_t>                                     
struct matrix_array_value : public matrix_array<matrix_array_value<Tag, Val_t>>
{
    template<Integer Row, Integer Col>
    using get_element_impl  = matrix_array_value_get_elem<Tag, Val_t, Row, Col>;
};

// cont_value_mat array
template<Tag_matrix_cdata Tag, Value Val_t, Integer Row, Integer Col>
struct ma_const_value_get_elem_scalar
    : mk::scal_data_const_value_tag<ma_const_value_get_elem_scalar<Tag, Val_t, Row, Col>>
{
    template<class Val_loc>
    static constexpr Val_loc value()
    {
        return Val_loc(Tag::value<Val_t, Row, Col>());
    };
};

template<Tag_matrix_cdata Tag, Value Val_t, Integer Row, Integer Col>
struct matrix_array_const_value_get_elem
{
    using type      = scal_data_const_value<ma_const_value_get_elem_scalar
                            <Tag, Val_t, Row, Col>, Val_t>;

    template<class Val_loc>
    static constexpr Val_loc value()
    {
        return Val_loc(Tag::value<Val_t, Row, Col>());
    };
};

template<Tag_matrix_cdata Tag, Value Val_t>                                     
struct matrix_array_const_value : public matrix_array<matrix_array_const_value<Tag, Val_t>>
{
    template<Integer Row, Integer Col>
    using get_element_impl  = matrix_array_const_value_get_elem<Tag, Val_t, Row, Col>;
};

template<class Tag>                                     
struct gen_array : public mkd::matrix_array<gen_array<Tag>>
{
    using this_type     = gen_array<Tag>;

    template<Integer Row, Integer Col>
    using get_element_impl  = mkd::lazy_type<mkd::element<Tag, Row, Col>>;
};

template<class Tag> 
struct output_array : public matrix_array<output_array<Tag>>
{
    using this_type     = output_array<Tag>;

    template<Integer Row, Integer Col>
    using get_element_impl  = mkd::lazy_type<mkd::element<Tag, Row, Col>>;
};

template<class Tag, Integer Mat_Rows, Integer Mat_Cols>
struct temp_output_array : public matrix_array<temp_output_array<Tag, Mat_Rows, Mat_Cols>>
{
    using this_type = temp_output_array<Tag, Mat_Rows, Mat_Cols>;

    template<Integer Row, Integer Col>
    using get_element_impl  = mkd::lazy_type<get_temporary<Tag, Mat_Rows, Mat_Cols, Row, Col>>;
};

template<class Tag, class... Assign_List>
struct virtual_array : public matrix_array<virtual_array<Tag, Assign_List ...>>
{
    using this_type     = virtual_array<Tag, Assign_List...>;

    template<Integer Row, Integer Col>
    using get_element_impl  = mkd::get_virtual_array_assignment<Row, Col, Assign_List ...>;
};

template<class Tag, Integer M, Integer N, class Array>
struct mat_ufunc_array : public matrix_array<mat_ufunc_array<Tag, M, N, Array>>
{
    using this_type = mat_ufunc_array<Tag, M, N, Array>;

    template<Integer Row, Integer Col>
    using get_element_impl  = mat_ufunc_array_get_elem<this_type, Row, Col>;
};

template<class Tag, Integer Mat_Rows, Integer Mat_Cols, bool Force>
struct mat_temp_array : public matrix_array<mat_temp_array<Tag, Mat_Rows, Mat_Cols, Force>>
{
    using this_type = mat_temp_array<Tag, Mat_Rows, Mat_Cols, Force>;

    template<Integer Row, Integer Col>
    using get_element_impl  = mkd::lazy_type<get_temporary<Tag, Mat_Rows, Mat_Cols, Row, Col>>;
};

template<class Ret_Tag>
struct empty_array : public matrix_array<empty_array<Ret_Tag>>
{
    using this_type = empty_array<Ret_Tag>;

    template<Integer Row, Integer Col>
    using get_element_impl  = error_access_to_empty_array<this_type, Row, Col>;
};

template<Integer M, Integer N, class Array1, class Array2>
struct mat_assign_array : public matrix_array<mat_assign_array<M, N, Array1, Array2>>
{
    using this_type = mat_assign_array<M, N, Array1, Array2>;

    template<Integer Row, Integer Col>
    using get_element_impl  = mat_assign_array_get_elem<this_type, Row, Col>;
};

template<Integer M,Integer N,class Array1,class Colon, Integer M2, Integer N2, class Array2> 
struct mat_assign_array_colon : public matrix_array<
                                mat_assign_array_colon<M, N, Array1, Colon, M2, N2, Array2>>
{
    using this_type = mat_assign_array_colon<M, N, Array1, Colon, M2, N2, Array2>;

    template<Integer Row, Integer Col>
    using get_element_impl  = mat_assign_array_colon_get_elem<this_type, Row, Col>;
};

template<Integer M, Integer N, class Array1, class Array2>
struct mat_scal_assign_array : public matrix_array<
                                    mat_scal_assign_array<M, N, Array1, Array2>>
{
    using this_type = mat_scal_assign_array<M, N, Array1, Array2>;

    template<Integer Row, Integer Col>
    using get_element_impl  = mat_scal_assign_array_get_elem<this_type, Row, Col>;
};

template<class Tag, Integer M, Integer N, class Array1, class Array2>
struct mat_bfunc_array : public matrix_array<mat_bfunc_array<Tag, M, N, Array1, Array2>>
{
    using this_type  = mat_bfunc_array<Tag, M, N, Array1, Array2>;

    template<Integer Row, Integer Col>
    using get_element_impl  = mat_bfunc_array_get_elem<this_type, Row, Col>;
};

template<class Tag, Integer M, Integer N, class Array1, class Array2>
struct mat_scal_bfunc_array : public matrix_array<mat_scal_bfunc_array<Tag, M, N,Array1, Array2>>
{
    using this_type  = mat_scal_bfunc_array<Tag, M, N, Array1, Array2>;

    template<Integer Row, Integer Col>
    using get_element_impl  = mat_scal_bfunc_array_get_elem<this_type, Row, Col>;
};

template<class Tag, Integer M, Integer N, class Array1, class Array2>
struct scal_mat_bfunc_array : public matrix_array<scal_mat_bfunc_array<Tag, M, N,Array1, Array2>>
{
    using this_type  = scal_mat_bfunc_array<Tag, M, N, Array1, Array2>;

    template<Integer Row, Integer Col>
    using get_element_impl  = scal_mat_bfunc_array_get_elem<this_type, Row, Col>;
};

template<class Array, Integer Row, Integer Col>
struct error_access_to_empty_array
{
    static_assert(md::dependent_false<Array>::value, 
                "access to elements of empty_array is not allowed; this array should not be used explicitly");
};

}}}
