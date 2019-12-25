#pragma once

#include "mkgen/TODO/matrix/ct_matrix.h"
#include "mkgen/TODO/expression/ct_matrix_expr.inl"
#include "mkgen/TODO/utils/utils.h"
#include "mkgen/matrix/dependency.h"
#include "mkgen/matrix/scalar.h"
#include "mkgen/details/matrix/matrix_arrays.h"
#include "mkgen/details/matrix/element.h"

namespace matcl { namespace mkgen
{

template<class Array, Integer Row, Integer Col>
struct mat_bfunc_array_get_elem
{};

template<class Tag, Integer M, Integer N, class Array1, class Array2,
        Integer Row, Integer Col>
struct mat_bfunc_array_get_elem<mkd::mat_bfunc_array<Tag,M,N,Array1,Array2>, Row, Col>
{
    using elem_1    = typename Array1 :: template get_element<Row, Col>::type;
    using elem_2    = typename Array2 :: template get_element<Row, Col>::type;
    using new_item  = typename expr_bfunc<Tag,elem_1,elem_2>::type;
    using type      = new_item;
};

template<class Array, Integer Row, Integer Col>
struct mat_scal_bfunc_array_get_elem
{};

template<class Tag,Integer M, Integer N, class Array1, Scal_data Array2, DPS Deps2,
        Integer Row, Integer Col>
struct mat_scal_bfunc_array_get_elem<mkd::mat_scal_bfunc_array<Tag,M,N,Array1, 
                    ct_scalar<Array2, Deps2>>, Row, Col>
{
    using elem_1    = typename Array1 :: template get_element<Row, Col>::type;
    using elem_2    = typename ct_scalar<Array2,Deps2>;
    using new_item  = typename expr_bfunc<Tag,elem_1,elem_2>::type;
    using type      = new_item;
};

template<class Array, Integer Row, Integer Col>
struct scal_mat_bfunc_array_get_elem
{};

template<class Tag,Integer M, Integer N, class Array1, Scal_data Array2, DPS Deps2,
        Integer Row, Integer Col>
struct scal_mat_bfunc_array_get_elem<mkd::scal_mat_bfunc_array<Tag,M,N,Array1, 
            ct_scalar<Array2, Deps2>>, Row, Col>
{
    using elem_1    = typename Array1 :: template get_element<Row, Col>::type;
    using elem_2    = typename ct_scalar<Array2,Deps2>;
    using new_item  = typename expr_bfunc<Tag,elem_2,elem_1>::type;
    using type      = new_item;
};

//----------------------------------------------------------------------------------
//                              get_elem
//----------------------------------------------------------------------------------
template<Integer M, Integer N, Mat_array Array, DPS Deps1, Integer Row, Integer Col>
struct get_elem<ct_matrix<M,N,Array,Deps1>,Row,Col>
{
    static_assert(Row >= 1 && Col >= 1 && Row <= M && Col <= N, "invalid element");
    using type = typename Array :: template get_element<Row,Col>::type;
};

template<class Elem, Integer Step>
struct make_element_step
{
    using type = mkd::element_step<Elem,Step>;
};
template<class Elem>
struct make_element_step<Elem,1>
{
    using type = Elem;
};

template<class Array, Integer Row, Integer Col>
struct sub_array_2_get_elem
{};

template<class Array_t, Integer Offset1, Integer Offset2, Integer Step1, Integer Step2,
        Integer Row, Integer Col>
struct sub_array_2_get_elem<mkd::sub_array_2<Array_t,Offset1,Offset2,Step1,Step2>, Row, Col>
{
    static const Integer row2 = Offset1 + (Row-1)*Step1 + 1;
    static const Integer col2 = Offset2 + (Col-1)*Step2 + 1;

    using type_1    = typename Array_t:: template get_element<row2,col2>::type;
    using type      = typename make_element_step<type_1,Step1>::type;
};

// get_elem requires scalar_data<> return, not scalar
// TODO: remove this function
template<class T>
struct correct_scalar_get_elem
{
    using type = T;
};

template<Scal_data Data, DPS Deps>
struct correct_scalar_get_elem<ct_scalar<Data, Deps>>
{
    using type = Data;
};

template<class Array, Integer Row, Integer Col>
struct mat_ufunc_array_get_elem
{};

template<class Tag, Integer M, Integer N, class Array, 
        Integer Row, Integer Col>
struct mat_ufunc_array_get_elem<mkd::mat_ufunc_array<Tag, M, N, Array>, Row, Col>
{
    using elem      = typename Array :: template get_element<Row, Col>::type;
    using new_item  = typename make_expr_ufunc<Tag,elem>::type;
    using type      = typename correct_scalar_get_elem<new_item>::type;
};

template<class Array, Integer Row, Integer Col>
struct sub_array_1_get_elem
{};

template<class Array_t, Integer Offset, Integer Step, Integer Row, Integer Col>
struct sub_array_1_get_elem<mkd::sub_array_1<Array_t, Offset, Step>, Row, Col>
{
    static_assert(Col == 1, "this submatrix has only one column");
    static const Integer pos = Offset + (Row -1) * Step + 1;

    using type_1    = typename Array_t :: template get_element<pos,1>::type;
    using type      = typename make_element_step<type_1,Step>::type;
};

}}

//TODO
#if 0
template<class Tag, class Array, class Deps, Integer Row, Integer Col>
struct get_array_elem<details::scalar_ufunc_array<Tag,Array,Deps>, Row, Col>
{
    using elem      = typename ct_scalar<Array, Deps>;
    using new_item  = expr_ufunc<Tag,elem>;
    using type      = new_item;
};

template<class Tag,Scal_data Array1, DPS Deps1, Scal_data Array2, DPS Deps2, Integer Row, Integer Col>
struct get_array_elem<details::scalar_bfunc_array<Tag,ct_scalar<Array1,Deps1>,ct_scalar<Array2,Deps2>>, Row, Col>
{
    using elem_1    = typename ct_scalar<Array1,Deps1>;
    using elem_2    = typename ct_scalar<Array2,Deps1>;
    using new_item  = typename expr_bfunc<Tag,elem_1,elem_2>::type;
    using type      = new_item;
};

// get a type that represents an element at row row and column col in a matrix with
// leading dimension LD. if_expr user defines a new operator, which creates a matrix with 
// symolic elements of new type stored in Array Array, then this class must be 
// specialized to define how to access to given element.
//TODO: check if return satisfy scalar_data requirement
template<class Array, Integer Row, Integer Col>
struct get_array_elem {};


//technical arrays
                                                        struct call_array_type;
#endif