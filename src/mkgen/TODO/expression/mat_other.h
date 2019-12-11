#pragma once

#include "mkgen/TODO/matrix/ct_matrix.h"
#include "mkgen/TODO/expression/ct_matrix_expr.inl"
#include "mkgen/TODO/utils/utils.h"
#include "mkgen/TODO/evaler/dependency.h"
#include "mkgen/matrix/scalar.h"

namespace matcl { namespace mkgen
{

//----------------------------------------------------------------------------------
//                              mat_trans
//----------------------------------------------------------------------------------
template<Integer M, Integer N, class Array>
struct mat_trans_array{};

template<Integer M, Integer N, class Array, Integer Row, Integer Col>
struct get_array_elem<mat_trans_array<M,N,Array>, Row, Col>
{
    using elem  = typename get_array_elem<Array, Col, Row>::type;
    using type  = elem;
};

template<Integer M, Integer N, class Array,class Deps1>
struct mat_trans<ct_matrix<M,N,Array,Deps1>>
{
    using array_type    = mat_trans_array<M, N, Array>;
    using type          = ct_matrix<N, M, array_type,Deps1>;
};

template<class Array, class Deps>
struct mat_trans<ct_scalar<Array,Deps>>
{
    using type          = ct_scalar<Array,Deps>;
};

//----------------------------------------------------------------------------------
//                              mat_ctrans
//----------------------------------------------------------------------------------
template<Integer M, Integer N, class Array>
struct mat_ctrans_array{};

template<Integer M, Integer N, class Array, Integer Row, Integer Col>
struct get_array_elem<mat_ctrans_array<M,N,Array>, Row, Col>
{
    using elem      = typename get_array_elem<Array, Col, Row>::type;
    using new_item  = typename expr_ctrans<elem>::type; 
    using type      = new_item;                       
};

template<class Array, class Deps, Integer Row, Integer Col>
struct get_array_elem<details::scalar_ctrans_array<Array, Deps>, Row, Col>
{
    //TODO:
    using elem      = ct_scalar<Array,Deps>;
    using new_item  = typename expr_ctrans<elem>::type;
    using type      = new_item;
};

template<Integer M, Integer N, class Array,class Deps1>
struct mat_ctrans<ct_matrix<M,N,Array,Deps1>>
{
    using array_type    = mat_ctrans_array<M, N, Array>;
    using type          = ct_matrix<N, M, array_type,Deps1>;
};

template<class Array, class Deps>
struct mat_ctrans<ct_scalar<Array,Deps>>
{
    using array_type    = details::scalar_ctrans_array<Array,Deps>;
    using type          = ct_scalar<details::scalar_data<array_type>,Deps>;
};

//----------------------------------------------------------------------------------
//                              func_unary
//----------------------------------------------------------------------------------
template<class Tag, Integer M, Integer N, class Array>
struct mat_ufunc_array{};

template<class Tag, Integer M, Integer N, class Array, Integer Row, Integer Col>
struct get_array_elem<mat_ufunc_array<Tag,M,N,Array>, Row, Col>
{
    using elem      = typename get_array_elem<Array, Row, Col>::type;
    using new_item  = typename expr_ufunc<Tag,elem>::type;
    using type      = new_item;
};

template<class Tag, class Array, class Deps, Integer Row, Integer Col>
struct get_array_elem<details::scalar_ufunc_array<Tag,Array,Deps>, Row, Col>
{
    //TODO
    using elem      = typename ct_scalar<details::scalar_data<Array>,Deps>;
    using new_item  = typename expr_ufunc<Tag,elem>::type;
    using type      = new_item;
};

template<class Tag, Integer M, Integer N, class Array, class Deps1>
struct func_unary<Tag, ct_matrix<M,N,Array,Deps1>>
{
    using array_type    = mat_ufunc_array<Tag, M, N, Array>;
    using type          = ct_matrix<M, N, array_type,Deps1>;
};

template<class Tag, class Array, class Deps>
struct func_unary<Tag, ct_scalar<Array,Deps>>
{
    //TODO
     using array_type   = details::scalar_ufunc_array<Tag,Array,Deps>;
     using type         = ct_scalar<details::scalar_data<array_type>,Deps>;
};

//----------------------------------------------------------------------------------
//                              func_bin
//----------------------------------------------------------------------------------
template<class Tag, Integer M, Integer N, class Array1, class Array2>
struct mat_bfunc_array{};

template<class Tag, Integer M, Integer N, class Array1, class Array2>
struct mat_scal_bfunc_array{};

template<class Tag, Integer M, Integer N, class Array1, class Array2>
struct scal_mat_bfunc_array{};

template<class Tag, Integer M, Integer N, class Array1, class Array2,
        Integer Row, Integer Col>
struct get_array_elem<mat_bfunc_array<Tag,M,N,Array1,Array2>, Row, Col>
{
    using elem_1    = typename get_array_elem<Array1, Row, Col>::type;
    using elem_2    = typename get_array_elem<Array2, Row, Col>::type;
    using new_item  = typename expr_bfunc<Tag,elem_1,elem_2>::type;
    using type      = new_item;
};

template<class Tag,Integer M, Integer N, class Array1, class Array2, class Deps2,
        Integer Row, Integer Col>
struct get_array_elem<mat_scal_bfunc_array<Tag,M,N,Array1, ct_scalar<Array2,Deps2>>, Row, Col>
{
    using elem_1    = typename get_array_elem<Array1, Row, Col>::type;
    using elem_2    = typename ct_scalar<Array2,Deps2>;
    using new_item  = typename expr_bfunc<Tag,elem_1,elem_2>::type;
    using type      = new_item;
};

template<class Tag,Integer M, Integer N, class Array1, class Array2, class Deps2,
        Integer Row, Integer Col>
struct get_array_elem<scal_mat_bfunc_array<Tag,M,N,Array1, ct_scalar<Array2,Deps2>>, Row, Col>
{
    using elem_1    = typename get_array_elem<Array1, Row, Col>::type;
    using elem_2    = typename ct_scalar<Array2,Deps2>;
    using new_item  = typename expr_bfunc<Tag,elem_2,elem_1>::type;
    using type      = new_item;
};

template<class Tag,class Array1, class Deps1, class Array2, class Deps2, Integer Row, Integer Col>
struct get_array_elem<details::scalar_bfunc_array<Tag,ct_scalar<Array1,Deps1>,ct_scalar<Array2,Deps2>>, Row, Col>
{
    using elem_1    = typename ct_scalar<Array1,Deps1>;
    using elem_2    = typename ct_scalar<Array2,Deps1>;
    using new_item  = typename expr_bfunc<Tag,elem_1,elem_2>::type;
    using type      = new_item;
};

template<class Tag,Integer M1_M2, Integer N1_N2, class Array1, class Deps1, class Array2, class Deps2>
struct func_bin<Tag, ct_matrix<M1_M2,N1_N2,Array1,Deps1>, ct_matrix<M1_M2,N1_N2,Array2,Deps2>>
{
    using array_type    = mat_bfunc_array<Tag,M1_M2, N1_N2, Array1, Array2>;
    using deps          = typename link_deps<Deps1, Deps2>::type;
    using type          = ct_matrix<M1_M2, N1_N2, array_type,deps>;
};

template<class Tag,Integer M1, Integer N1, class Array1, class Deps1, class Array2, class Deps2>
struct func_bin<Tag, ct_matrix<M1,N1,Array1,Deps1>, ct_scalar<Array2,Deps2>>
{
    using array_type    = mat_scal_bfunc_array<Tag,M1,N1,Array1,ct_scalar<Array2,Deps2>>;
    using deps          = typename link_deps<Deps1, Deps2>::type;
    using type          = ct_matrix<M1, N1, array_type,deps>;
};

template<class Tag,Integer M1, Integer N1, class Array1, class Deps1, class Array2,class Deps2>
struct func_bin<Tag, ct_scalar<Array2, Deps2>, ct_matrix<M1,N1,Array1,Deps1>>
{
    using array_type    = scal_mat_bfunc_array<Tag,M1, N1, Array1, ct_scalar<Array2,Deps2>>;
    using deps          = typename link_deps<Deps1, Deps2>::type;
    using type          = ct_matrix<M1, N1, array_type,deps>;
};

template<class Tag,class Array1, class Deps1, class Array2, class Deps2>
struct func_bin<Tag, ct_scalar<Array1,Deps1>, ct_scalar<Array2,Deps2>>
{
    using array_type    = details::scalar_bfunc_array<Tag,ct_scalar<Array1,Deps1>,ct_scalar<Array2,Deps2>>;
    using deps          = typename link_deps<Deps1, Deps2>::type;
    using type          = ct_scalar<details::scalar_data<array_type>,deps>;
};

template<class Tag,Integer M1, Integer N1, Integer M2, Integer N2, class Array1, class Deps1,
        class Array2, class Deps2>
struct func_bin<Tag,ct_matrix<M1,N1,Array1,Deps1>, ct_matrix<M2,N2,Array2,Deps2>>
{
    static_assert(M1 == M2 && N1 == N2, "invalid matrix operation, check size of matrices");
};

//----------------------------------------------------------------------------------
//                              get_elem
//----------------------------------------------------------------------------------
template<Integer M, Integer N, class Array, class Deps1, Integer Row, Integer Col>
struct get_elem<ct_matrix<M,N,Array,Deps1>,Row,Col>
{
    static_assert(Row >= 1 && Col >= 1 && Row <= M && Col <= N, "invalid element");
    using type = typename get_array_elem<Array,Row,Col>::type;
};

//----------------------------------------------------------------------------------
//                              get_array_elem
//----------------------------------------------------------------------------------
template<class Tag, Integer Row, Integer Col>
struct get_array_elem<array<Tag>,Row,Col> 
{
    using type = element<Tag,Row,Col>;
};
template<class Tag, Integer Row, Integer Col>
struct get_array_elem<output_array<Tag>,Row,Col> 
{
    using type = element<Tag,Row,Col>;
};
template<class Tag, Integer Mat_Rows, Integer Mat_Cols, Integer Row, Integer Col>
struct get_array_elem<temp_output_array<Tag,Mat_Rows,Mat_Cols>,Row,Col> 
{
    using type = get_temporary<Tag,Mat_Rows,Mat_Cols, Row,Col>;
};

template<class Tag, Integer Row, Integer Col>
struct get_array_elem<const_array<Tag>,Row,Col> 
{
    using type = decltype(Tag::get_elem<Row,Col>());
};

template<class Elem, Integer Step>
struct make_element_step
{
    using type = element_step<Elem,Step>;
};
template<class Elem>
struct make_element_step<Elem,1>
{
    using type = Elem;
};

template<class Array_t, Integer Offset, Integer Step, Integer Row, Integer Col>
struct get_array_elem<sub_array_1<Array_t,Offset,Step>,Row,Col> 
{
    static_assert(Col == 1, "this submatrix has only one column");
    static const Integer pos = Offset + (Row -1) * Step + 1;

    using type_1    = typename get_array_elem<Array_t,pos,1>::type;
    using type      = typename make_element_step<type_1,Step>::type;
};
template<class Array_t, Integer Offset1, Integer Offset2, Integer Step1, Integer Step2,
        Integer Row, Integer Col>
struct get_array_elem<sub_array_2<Array_t,Offset1,Offset2,Step1,Step2>,Row,Col> 
{
    static const Integer row2 = Offset1 + (Row-1)*Step1 + 1;
    static const Integer col2 = Offset2 + (Col-1)*Step2 + 1;

    using type_1    = typename get_array_elem<Array_t,row2,col2>::type;
    using type      = typename make_element_step<type_1,Step1>::type;
};

}}