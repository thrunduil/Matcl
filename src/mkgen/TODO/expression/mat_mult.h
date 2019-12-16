#pragma once

#include "mkgen/TODO/matrix/ct_matrix.h"
#include "mkgen/TODO/expression/ct_matrix_expr.inl"
#include "mkgen/TODO/utils/utils.h"
#include "mkgen/matrix/dependency.h"
#include "mkgen/mkgen_fwd.h"

namespace matcl { namespace mkgen
{

template<class T, bool Is_scalar_data = mkd::is_valid_scalar_data<T>::value>
struct make_scalar_data;

template<class T>
struct make_scalar_data<T, true>
{
    using type = T;
};

template<class T>
struct make_scalar_data<T, false>
{
    using type = mkd::scalar_expr_data<T>;
};

template<class Expr_Type, class Deps>
struct make_scalar
{
    using scalar_data   = typename make_scalar_data<Expr_Type>::type;
    using type = ct_scalar<scalar_data,Deps>;
};

template<class Expr_Type, class Deps1, class Deps2>
struct make_scalar<ct_scalar<Expr_Type, Deps1>, Deps2>
{
    using type = ct_scalar<Expr_Type,Deps2>;
};

//----------------------------------------------------------------------------------
//                              mat_mult_array
//----------------------------------------------------------------------------------
template<class List_1, class List_2>
struct merge_dots
{
    static_assert(md::dependent_false<List_1>::value,
                  "this type should not be instantiated");
};
template<class ... Arg_11, class ... Arg_12, class ... Arg_21, class ... Arg_22>
struct merge_dots<expr_dot<list::list<Arg_11...>, list::list<Arg_12...>>,
                  expr_dot<list::list<Arg_21...>, list::list<Arg_22...>>>
{
    using type = expr_dot<list::list<Arg_11..., Arg_21...>,
                          list::list<Arg_12..., Arg_22...>>;
};

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
    using type = zero;
};
template<Integer Row, Integer Col, Integer K_start, Integer K_end, class Array_1, class Array_2>
struct make_array_mat_mult<Row,Col,K_start,K_end,1,Array_1, Array_2>
{
    using elem_1 = typename Array_1 :: template get_element<Row, K_start>::type;
    using elem_2 = typename Array_2 :: template get_element<K_start, Col>::type;
    //using type = typename make_mult<elem_1,elem_2>::type;
    using type   = expr_dot<list::list<elem_1>, list::list<elem_2>>;
};
template<Integer Row, Integer Col, Integer K_start, Integer K_end, class Array_1, class Array_2>
struct make_array_mat_mult<Row,Col,K_start,K_end,2,Array_1, Array_2>
{
    using elem_1_1      = typename Array_1 :: template get_element<Row, K_start>::type;
    using elem_2_1      = typename Array_2 :: template get_element<K_start, Col>::type;
    //using new_item_1  = typename make_mult<elem_1_1,elem_2_1>::type;

    using elem_1_2      = typename Array_1 :: template get_element<Row, K_start+1>::type;
    using elem_2_2      = typename Array_2 :: template get_element<K_start+1, Col>::type;
    //using new_item_2  = typename make_mult<elem_1_2,elem_2_2>::type;

    //using type        = typename make_plus<new_item_1, new_item_2>::type;
    using type          = expr_dot<list::list<elem_1_1,elem_1_2>, list::list<elem_2_1,elem_2_2>>;
};
template<Integer Row, Integer Col, Integer K_start, Integer K_end, class Array_1, class Array_2>
struct make_array_mat_mult<Row,Col,K_start,K_end,3,Array_1, Array_2>
{
    using elem_1_1      = typename Array_1 :: template get_element<Row, K_start>::type;
    using elem_2_1      = typename Array_2 :: template get_element<K_start, Col>::type;
    //using new_item_1  = typename make_mult<elem_1_1, elem_2_1>::type;

    using elem_1_2      = typename Array_1 :: template get_element<Row, K_start + 1>::type;
    using elem_2_2      = typename Array_2 :: template get_element<K_start + 1, Col>::type;
    //using new_item_2  = typename make_mult<elem_1_2, elem_2_2>::type;

    using elem_1_3      = typename Array_1 :: template get_element<Row, K_start+2>::type;
    using elem_2_3      = typename Array_2 :: template get_element<K_start+2, Col>::type;
    //using new_item_3  = typename make_mult<elem_1_3,elem_2_3>::type;

    //using plus_23     = typename make_plus<new_item_2, new_item_3>::type;
    //using type        = typename make_plus<new_item_1, plus_23>::type;

    using type          = expr_dot<list::list<elem_1_1,elem_1_2, elem_1_3>, 
                                   list::list<elem_2_1,elem_2_2, elem_2_3>>;
};
template<Integer Row, Integer Col, Integer K_start, Integer K_end, class Array_1, class Array_2>
struct make_array_mat_mult<Row,Col,K_start,K_end,4,Array_1, Array_2>
{
    using elem_1_1      = typename Array_1 :: template get_element<Row, K_start>::type;
    using elem_2_1      = typename Array_2 :: template get_element<K_start, Col>::type;
    //using new_item_1  = typename make_mult<elem_1_1, elem_2_1>::type;

    using elem_1_2      = typename Array_1 :: template get_element<Row, K_start + 1>::type;
    using elem_2_2      = typename Array_2 :: template get_element<K_start + 1, Col>::type;
    //using new_item_2  = typename make_mult<elem_1_2, elem_2_2>::type;

    using elem_1_3      = typename Array_1 :: template get_element<Row, K_start + 2>::type;
    using elem_2_3      = typename Array_2 :: template get_element<K_start + 2, Col>::type;
    //using new_item_3  = typename make_mult<elem_1_3, elem_2_3>::type;

    using elem_1_4      = typename Array_1 :: template get_element<Row, K_start+3>::type;
    using elem_2_4      = typename Array_2 :: template get_element<K_start+3, Col>::type;
    //using new_item_4  = typename make_mult<elem_1_4,elem_2_4>::type;

    //using plus_12     = typename make_plus<new_item_1, new_item_2>::type;
    //using plus_34     = typename make_plus<new_item_3, new_item_4>::type;
    //using type        = typename make_plus<plus_12, plus_34>::type;

    using type          = expr_dot<list::list<elem_1_1,elem_1_2, elem_1_3, elem_1_4>, 
                                   list::list<elem_2_1,elem_2_2, elem_2_3, elem_2_4>>;
};

template<class Array, Integer Row, Integer Col>
struct mat_mult_array_get_elem
{};

template<Integer K, class Array1, class Array2, Integer Row, Integer Col>
struct mat_mult_array_get_elem<mkd::mat_mult_array<K, Array1, Array2>, Row, Col>
{
    using type0 = typename make_array_mat_mult<Row,Col,1,K,K,Array1,Array2>::type;
    using type  = typename correct_scalar_get_elem<type0>::type;
};

template<class Array1, class Array2, Integer Row, Integer Col>
struct mat_mult_array_get_elem<mkd::mat_mult_array<0, Array1, Array2>, Row, Col>
{
    using type  = typename zero :: data_type;
};

template<class Array1, class Array2, Integer Row, Integer Col>
struct mat_mult_array_get_elem<mkd::mat_mult_array<1, Array1, Array2>, Row, Col>
{
    using elem_1 = typename Array1 :: template get_element<Row, 1>::type;
    using elem_2 = typename Array2 :: template get_element<1, Col>::type;
    using elem_s = element_step<elem_2,0>;

    using type0  = typename make_mult<elem_1,elem_s>::type;
    using type   = typename correct_scalar_get_elem<type0>::type;
};

template<class Array, Integer Row, Integer Col>
struct mult_rows_array_get_elem
{};

template<class Array1, class Array2, Integer Row, Integer Col>
struct mult_rows_array_get_elem<mkd::mult_rows_array<Array1,Array2>, Row, Col>
{
    using elem_1    = typename Array1 :: template get_element<Row, Col>::type;
    using elem_2    = typename Array2 :: template get_element<Row, 1>::type;
    using type0     = typename make_mult<elem_1,elem_2>::type;
    using type      = typename correct_scalar_get_elem<type0>::type;
};

template<class Array, Integer Row, Integer Col>
struct mat_scal_mult_array_get_elem
{};

template<class Array1, class Array2, class Deps2, Integer Row, Integer Col>
struct mat_scal_mult_array_get_elem<mkd::mat_scal_mult_array<Array1, 
                            ct_scalar<Array2, Deps2>>, Row, Col>
{
    using elem_1    = typename Array1 :: template get_element<Row,Col>::type;
    using elem_2    = typename ct_scalar<Array2,Deps2>;
    using new_item  = typename make_mult<elem_1,elem_2>::type;
    using type0     = new_item;
    using type      = typename correct_scalar_get_elem<type0>::type;
};

//----------------------------------------------------------------------------------
//                              mat_mult
//----------------------------------------------------------------------------------
template<Integer M1, Integer N1_M2, class Array1, class Deps1,
        Integer N2, class Array2, class Deps2>
struct mat_mult<ct_matrix<M1,N1_M2,Array1,Deps1>, ct_matrix<N1_M2,N2,Array2,Deps2>>
{
    using array_type    = mkd::mat_mult_array<N1_M2, Array1, Array2>;
    using deps          = typename link_deps<Deps1, Deps2>::type;
    using type          = ct_matrix<M1, N2, array_type,deps>;
};

template<Integer M1, Integer N1, class Array1, class Deps1, class Array2,class Deps2>
struct mat_mult<ct_matrix<M1,N1,Array1,Deps1>, ct_scalar<Array2,Deps2>>
{
    using array_type    = mkd::mat_scal_mult_array<Array1, ct_scalar<Array2, Deps2>>;
    using deps          = typename link_deps<Deps1, Deps2>::type;
    using type          = ct_matrix<M1, N1, array_type,deps>;
};

template<class Array1, class Deps1, Integer M2, Integer N2, class Array2, class Deps2>
struct mat_mult<ct_scalar<Array1,Deps1>, ct_matrix<M2,N2,Array2,Deps2>>
{
    using array_type    = mkd::mat_scal_mult_array<Array2, ct_scalar<Array1, Deps1>>;
    using deps          = typename link_deps<Deps1, Deps2>::type;
    using type          = ct_matrix<M2, N2, array_type,deps>;
};

template<class Array1, class Deps1, class Array2, class Deps2>
struct mat_mult<ct_scalar<Array1,Deps1>, ct_scalar<Array2,Deps2>>
{
    using scal_1    = ct_scalar<Array1, Deps1>;
    using scal_2    = ct_scalar<Array2, Deps2>;

    using mult_type = typename make_mult<scal_1, scal_2>::type;
    using deps      = typename link_deps<Deps1, Deps2>::type;
    using type      = typename make_scalar<mult_type, deps>::type;
};

template<Integer M1, Integer N1, class Array1, class Deps1,
         Integer M2, Integer N2, class Array2, class Deps2>
struct mat_mult<ct_matrix<M1,N1,Array1,Deps1>, ct_matrix<M2,N2,Array2,Deps2>>
{
    static_assert(N1 == M2, "invalid matrix product, check size of matrices");
};

//----------------------------------------------------------------------------------
//                              make_mult_rows
//----------------------------------------------------------------------------------

template<Integer M1_M2, Integer N1, class Array1, class Deps1, class Array2, class Deps2>
struct make_mult_rows<ct_matrix<M1_M2,N1,Array1,Deps1>, ct_matrix<M1_M2,1,Array2,Deps2>>
{
    using array_type    = mkd::mult_rows_array<Array1, Array2>;
    using deps          = typename link_deps<Deps1, Deps2>::type;
    using type          = ct_matrix<M1_M2, N1, array_type,deps>;
};

template<Integer M1, Integer N1, class Array1, class Deps1,
         Integer M2, Integer N2, class Array2, class Deps2>
struct make_mult_rows<ct_matrix<M1,N1,Array1,Deps1>, ct_matrix<M2,N2,Array2,Deps2>>
{
    static_assert(M1 == M2 && N2 == 1, "invalid mult_rows product, check size of matrices");
};

//----------------------------------------------------------------------------------
//                              make_mult_cols
//----------------------------------------------------------------------------------
template<Integer M1, Integer N1_M2, class Array1, class Deps1, class Array2, class Deps2>
struct make_mult_cols<ct_matrix<M1,N1_M2,Array1,Deps1>, ct_matrix<N1_M2,1,Array2,Deps2>>
{
    using array_type    = mkd::mult_cols_array<Array1, Array2>;
    using deps          = typename link_deps<Deps1, Deps2>::type;
    using type          = ct_matrix<M1, N1_M2, array_type,deps>;
};

template<Integer M1, Integer N1, class Array1, class Deps1,
         Integer M2, Integer N2, class Array2, class Deps2>
struct make_mult_cols<ct_matrix<M1,N1,Array1,Deps1>, ct_matrix<M2,N2,Array2,Deps2>>
{
    static_assert(M2 == N1 && N2 == 1, "invalid mult_cols product, check size of matrices");
};

template<class Array, Integer Row, Integer Col>
struct mult_cols_array_get_elem
{};

template<class Array1, class Array2, Integer Row, Integer Col>
struct mult_cols_array_get_elem<mkd::mult_cols_array<Array1,Array2>, Row, Col>
{
    using elem_1    = typename Array1 :: template get_element<Row, Col>::type;
    using elem_2    = typename Array2 :: template get_element<1, Col>::type;
    using type0     = typename make_mult<elem_1,elem_2>::type;
    using type      = typename correct_scalar_get_elem<type0>::type;
};

//----------------------------------------------------------------------------------
//                              make_mult
//----------------------------------------------------------------------------------
template<Integer M1_M2, Integer N1_N2, class Array1, class Deps1,
         class Array2, class Deps2>
struct make_mult_mat<ct_matrix<M1_M2,N1_N2,Array1,Deps1>, ct_matrix<M1_M2,N1_N2,Array2,Deps2>>
{
    using array_type    = mkd::mult_array<Array1, Array2>;
    using deps          = typename link_deps<Deps1, Deps2>::type;
    using type          = ct_matrix<M1_M2, N1_N2, array_type,deps>;
};

template<Integer M1, Integer N1, class Array1, class Deps1,
         Integer M2, Integer N2, class Array2, class Deps2>
struct make_mult_mat<ct_matrix<M1,N1,Array1,Deps1>, ct_matrix<M2,N2,Array2,Deps2>>
{
    static_assert(M1 == M2 && N1 == N2, "invalid mult product, check size of matrices");
};

template<class Array, Integer Row, Integer Col>
struct mult_array_get_elem
{};

template<class Array1, class Array2, Integer Row, Integer Col>
struct mult_array_get_elem<mkd::mult_array<Array1,Array2>, Row, Col>
{
    using elem_1    = typename Array1 :: template get_element<Row, Col>::type;
    using elem_2    = typename Array2 :: template get_elemend<Row, Col>::type;
    using type0     = typename make_mult<elem_1,elem_2>::type;
    using type      = typename correct_scalar_get_elem<type0>::type;
};

//----------------------------------------------------------------------------------
//                              div
//----------------------------------------------------------------------------------

template<Integer M1_M2, Integer N1_N2, class Array1, class Deps1, class Array2, class Deps2>
struct mat_div<ct_matrix<M1_M2,N1_N2,Array1,Deps1>, ct_matrix<M1_M2,N1_N2,Array2,Deps2>>
{
    using array_type    = mkd::div_array<Array1, Array2>;
    using deps          = typename link_deps<Deps1, Deps2>::type;
    using type          = ct_matrix<M1_M2, N1_N2, array_type,deps>;
};

template<Integer M1, Integer N1, class Array1, class Deps1, class Array2,class Deps2>
struct mat_div<ct_matrix<M1,N1,Array1,Deps1>, ct_scalar<Array2,Deps2>>
{
    using array_type    = mkd::div_array_mat_scal<Array1, ct_scalar<Array2, Deps2>>;
    using deps          = typename link_deps<Deps1, Deps2>::type;
    using type          = ct_matrix<M1, N1, array_type,deps>;
};

template<class Array1, class Deps1, Integer M2, Integer N2, class Array2, class Deps2>
struct mat_div<ct_scalar<Array1,Deps1>, ct_matrix<M2,N2,Array2,Deps2>>
{
    using array_type    = mkd::div_array_scal_mat<Array2, ct_scalar<Array1, Deps1>>;
    using deps          = typename link_deps<Deps1, Deps2>::type;
    using type          = ct_matrix<M2, N2, array_type,deps>;
};

template<class Array1, class Deps1, class Array2, class Deps2>
struct mat_div<ct_scalar<Array1,Deps1>, ct_scalar<Array2,Deps2>>
{
    //TODO
    using div_type      = typename make_div<ct_scalar<Array1,Deps1>, ct_scalar<Array2,Deps2>>::type;
    using deps          = typename link_deps<Deps1, Deps2>::type;
    using type          = typename make_scalar<div_type,deps>::type;
};

template<Integer M1, Integer N1, class Array1, class Deps1,
         Integer M2, Integer N2, class Array2, class Deps2>
struct mat_div<ct_matrix<M1,N1,Array1,Deps1>, ct_matrix<M2,N2,Array2,Deps2>>
{
    static_assert(M1 == M2 && N1 == N2, "invalid mult product, check size of matrices");
};

template<class Array, Integer Row, Integer Col>
struct div_array_get_elem
{};

template<class Array1, class Array2, Integer Row, Integer Col>
struct div_array_get_elem<mkd::div_array<Array1,Array2>, Row, Col>
{
    using elem_1    = typename Array1 :: template get_element<Row, Col>::type;
    using elem_2    = typename Array2 :: template get_element<Row, Col>::type;
    using type0     = typename make_div<elem_1,elem_2>::type;
    using type      = typename correct_scalar_get_elem<type0>::type;
};

template<class Array, Integer Row, Integer Col>
struct div_array_mat_scal_get_elem
{};

template<class Array1, class Scal2, Integer Row, Integer Col>
struct div_array_mat_scal_get_elem<mkd::div_array_mat_scal<Array1, Scal2>, Row, Col>
{
    using elem_1    = typename Array1 :: template get_element<Row, Col>::type;
    using elem_2    = Scal2;
    using type0     = typename make_div<elem_1,elem_2>::type;
    using type      = typename correct_scalar_get_elem<type0>::type;
};

template<class Array, Integer Row, Integer Col>
struct div_array_scal_mat_get_elem
{};

template<class Array2, class Scal1, Integer Row, Integer Col>
struct div_array_scal_mat_get_elem<mkd::div_array_scal_mat<Array2, Scal1>, Row, Col>
{
    using elem_1    = Scal1;
    using elem_2    = typename Array2 :: template get_element<Row, Col>::type;    
    using type0     = typename make_div<elem_1,elem_2>::type;
    using type      = typename correct_scalar_get_elem<type0>::type;
};


}}