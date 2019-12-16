#pragma once

#include "mkgen/TODO/matrix/ct_matrix.h"
#include "mkgen/TODO/expression/ct_matrix_expr.inl"
#include "mkgen/TODO/utils/utils.h"
#include "mkgen/matrix/dependency.h"

namespace matcl { namespace mkgen
{

//----------------------------------------------------------------------------------
//                              mat_plus
//----------------------------------------------------------------------------------

template<Integer M1_M2, Integer N1_N2, class Array1, class Deps1, class Array2, class Deps2>
struct mat_plus<ct_matrix<M1_M2,N1_N2,Array1,Deps1>, ct_matrix<M1_M2,N1_N2,Array2,Deps2>>
{
    using array_type    = mkd::mat_plus_array<M1_M2, N1_N2, Array1, Array2>;
    using deps          = typename link_deps<Deps1, Deps2>::type;
    using type          = ct_matrix<M1_M2, N1_N2, array_type,deps>;
};

template<Integer M1, Integer N1, class Array1, class Deps1, class Array2,class Deps2>
struct mat_plus<ct_matrix<M1,N1,Array1,Deps1>, ct_scalar<Array2,Deps2>>
{
    using array_type    = mkd::mat_scal_plus_array<M1,N1,Array1, ct_scalar<Array2,Deps2>>;
    using deps          = typename link_deps<Deps1, Deps2>::type;
    using type          = ct_matrix<M1, N1, array_type,deps>;
};

template<Integer M1, Integer N1, class Array1, class Deps1, class Array2, class Deps2>
struct mat_plus<ct_scalar<Array2,Deps2>, ct_matrix<M1,N1,Array1,Deps1>>
{
    using array_type    = mkd::mat_scal_plus_array<M1,N1,Array1, ct_scalar<Array2,Deps2>>;
    using deps          = typename link_deps<Deps1, Deps2>::type;
    using type          = ct_matrix<M1, N1, array_type,deps>;
};

template<class Array1, class Deps1, class Array2, class Deps2>
struct mat_plus<ct_scalar<Array1, Deps1>, ct_scalar<Array2, Deps2>>
{
    //TODO:
    using scal_1    = ct_scalar<Array1, Deps1>;
    using scal_2    = ct_scalar<Array2, Deps2>;
    using plus_type = typename make_plus<scal_1, scal_2>::type;
    using deps      = typename link_deps<Deps1, Deps2>::type;
    using type      = typename make_scalar<plus_type, deps>::type;
};

template<Integer M1, Integer N1, Integer M2, Integer N2, class Array1, class Deps1,
        class Array2, class Deps2>
struct mat_plus<ct_matrix<M1,N1,Array1,Deps1>, ct_matrix<M2,N2,Array2,Deps2>>
{
    static_assert(M1 == M2 && N1 == N2, "invalid matrix operation, check size of matrices");
};

//----------------------------------------------------------------------------------
//                              mat_minus
//----------------------------------------------------------------------------------
template<Integer M1_M2, Integer N1_N2, class Array1, class Deps1, class Array2, class Deps2>
struct mat_minus<ct_matrix<M1_M2,N1_N2,Array1,Deps1>, ct_matrix<M1_M2,N1_N2,Array2,Deps2>>
{
    using array_type    = mkd::mat_minus_array<M1_M2, N1_N2, Array1, Array2>;
    using deps          = typename link_deps<Deps1, Deps2>::type;
    using type          = ct_matrix<M1_M2, N1_N2, array_type,deps>;
};

template<Integer M1, Integer N1, class Array1, class Deps1, class Array2,class Deps2>
struct mat_minus<ct_matrix<M1,N1,Array1,Deps1>, ct_scalar<Array2,Deps2>>
{
    using array_type    = mkd::mat_scal_minus_array<M1,N1,Array1,ct_scalar<Array2,Deps2>>;
    using deps          = typename link_deps<Deps1, Deps2>::type;
    using type          = ct_matrix<M1, N1, array_type,deps>;
};

template<Integer M1, Integer N1, class Array1, class Deps1, class Array2, class Deps2>
struct mat_minus<ct_scalar<Array2,Deps2>, ct_matrix<M1,N1,Array1,Deps1>>
{
    using array_type    = mkd::scal_mat_minus_array<M1, N1, Array1, ct_scalar<Array2,Deps2>>;
    using deps          = typename link_deps<Deps1, Deps2>::type;
    using type          = ct_matrix<M1, N1, array_type,deps>;
};

template<class Array1, class Deps1, class Array2, class Deps2>
struct mat_minus<ct_scalar<Array1,Deps1>, ct_scalar<Array2,Deps2>>
{
    using scal_1        = ct_scalar<Array1, Deps1>;
    using scal_2        = ct_scalar<Array2, Deps2>;
    using minus_type    = typename make_minus<scal_1,scal_2>::type;
    using deps          = typename link_deps<Deps1, Deps2>::type;
    using type          = typename make_scalar<minus_type,deps>::type;
};

template<Integer M1, Integer N1, Integer M2, Integer N2, class Array1, class Deps1,
        class Array2, class Deps2>
struct mat_minus<ct_matrix<M1,N1,Array1,Deps1>, ct_matrix<M2,N2,Array2,Deps2>>
{
    static_assert(M1 == M2 && N1 == N2, "invalid matrix operation, check size of matrices");
};

//----------------------------------------------------------------------------------
//                              unary_minus
//----------------------------------------------------------------------------------
template<Integer M, Integer N, class Array,class Deps1>
struct unary_minus<ct_matrix<M,N,Array,Deps1>>
{
    using array_type    = mkd::mat_uminus_array<M, N, Array>;
    using type          = ct_matrix<N, M, array_type,Deps1>;
};

template<class Array, class Deps>
struct unary_minus<ct_scalar<Array,Deps>>
{
    using scal_1        = ct_scalar<Array, Deps>;
    using uminus_type   = typename make_uminus<scal_1>::type;
    using type          = typename make_scalar<uminus_type,Deps>::type;
};

template<class Array, Integer Row, Integer Col>
struct mat_uminus_array_get_elem
{};

template<Integer M, Integer N, class Array, Integer Row, Integer Col>
struct mat_uminus_array_get_elem<mkd::mat_uminus_array<M,N,Array>, Row, Col>
{
    using elem      = typename Array::template get_element<Col, Row>::type;
    using new_item  = typename make_uminus<elem>::type;
    using type      = new_item;
};

template<class Array, Integer Row, Integer Col>
struct scal_mat_minus_array_get_elem
{};

template<Integer M, Integer N, class Array1, class Array2, class Deps2, Integer Row, Integer Col>
struct scal_mat_minus_array_get_elem<mkd::scal_mat_minus_array<M,N,Array1, ct_scalar<Array2,Deps2>>, Row, Col>
{
    using elem_1    = typename Array1 :: template get_element<Row, Col>::type;
    using elem_2    = ct_scalar<Array2,Deps2>;
    using new_item  = typename make_minus<elem_2,elem_1>::type;
    using type      = new_item;
};

}}