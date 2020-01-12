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

#include "mkgen/expression/expressions.h"

namespace matcl { namespace mkgen { namespace details
{

//----------------------------------------------------------------------------------
//                              forward declarations
//----------------------------------------------------------------------------------

template<class Array, Integer Row, Integer Col>
struct mat_ufunc_array_get_elem;

template<class Tag, class Scal1>
struct make_ufunc_root;

//----------------------------------------------------------------------------------
//                              matrix arrays
//----------------------------------------------------------------------------------
template<class Tag, Integer M, Integer N, class Array>
struct mat_ufunc_array : public matrix_array<mat_ufunc_array<Tag, M, N, Array>>
{
    using this_type = mat_ufunc_array<Tag, M, N, Array>;

    template<Integer Row, Integer Col>
    using get_element_impl  = mat_ufunc_array_get_elem<this_type, Row, Col>;
};

//----------------------------------------------------------------------------------
//                              get_elem impl
//----------------------------------------------------------------------------------

template<class Array, Integer Row, Integer Col>
struct mat_ufunc_array_get_elem
{};

template<class Tag, Integer M, Integer N, class Array, 
        Integer Row, Integer Col>
struct mat_ufunc_array_get_elem<mkd::mat_ufunc_array<Tag, M, N, Array>, Row, Col>
{
    using elem      = typename Array :: template get_element<Row, Col>::type;
    using type      = typename make_ufunc_root<Tag, elem>::type;
};

//----------------------------------------------------------------------------------
//                              expr
//----------------------------------------------------------------------------------

template<class Tag, class Elem>
struct expr_ufunc : public mkd::scalar_data<expr_ufunc<Tag, Elem>>
{
    using this_type = expr_ufunc<Tag, Elem>;

    //TODO: add checks

    template<class Subs_Context>
    static void print(std::ostream& os, int prior)
    {
        Tag::print(os,details::prior_start);

        os << "(";
        elem::print<Subs_Context>(os,details::prior_start);
    };

    template<class Val, class Local_Storage>
    inline_expr
    static Val eval(const Local_Storage& ls)
    {
        Val v1  = Elem::eval<Val>(ls);
        Val tmp = Tag::eval<Val>(ls,v1);
        return tmp;
    };

    template<class Visitor>
    static void accept(Visitor& vis)
    {
        Elem::accept<Visitor>(vis);
    };

    template<class Void>
    using simplify      = this_type;

    static constexpr bool is_simplified()   { return true; };
};

//----------------------------------------------------------------------------------
//                              impl scalar
//----------------------------------------------------------------------------------

// apply unary function Tag to a scalar value 
template<class Tag, class S2>
struct make_ufunc_scal
{
    static_assert(md::dependent_false<S1>::value, 
                "this type should not be instantiated");
};

template<class UTag, Integer N1, Integer D1>
struct make_ufunc_scal<UTag, scal_data_rational<N1, D1>>
{
    //TODO
    /*
    using op    = mkd::rational_plus<N1,D1,N2,D2>;
    using type  = scal_data_rational<op::nominator,op::denominator>;
    */
};

template<class UTag, class Tag1, class Val1>
struct make_ufunc_scal<UTag, mkd::scal_data_const_value<Tag1,Val1>>
{
    //TODO
    /*
    using val   = decltype(std::declval<Val1>() + std::declval<Val2>());
    using tag   = scal_data_const_value_tag_plus<Tag1, Tag2>;

    using type  = scal_data_const_value<tag, val>;
    */
};

template<class UTag, class Tag1, class Val1>
struct make_ufunc_scal<UTag, mkd::scal_data_value<Tag1,Val1>>
{
    //TODO
    /*
    using val   = decltype(std::declval<Val1>() + std::declval<Val2>());
    using tag   = scal_data_value_tag_plus<Tag1, Tag2>;

    using type  = scal_data_value<tag, val>;
    */
};

//----------------------------------------------------------------------------------
//                              impl
//----------------------------------------------------------------------------------

template<class Tag, class T1,
        bool Is_Scal_1 = is_value_scalar_data<T1>::value>
struct make_ufunc_impl
{
    using type  = expr_ufunc<Tag, T1>;
};

template<class Tag, class S1>
struct make_ufunc_impl<Tag, S1, true>
{
    static_assert(is_value_scalar_data<S1>::value == true,
                  "invalid arguments");

    using type  = typename make_ufunc_scal<Tag, S1>::type;
};

//----------------------------------------------------------------------------------
//                              impl roots
//----------------------------------------------------------------------------------

// representation of unary function Tag(Scal1), where Scal1 is scalar_data, return
// scalar_data type
template<class Tag, class Scal1>
struct make_ufunc_root
{
    static const bool is_sd1    = mkd::is_valid_scalar_data<Scal1>::value;
    static_assert(is_sd1 == true, "scalar_data required");

    using S1s       = typename Scal1 :: template simplify<void>; 

    static const bool modif = (std::is_same<S1s, Scal1>::value == false);

    using type0  = typename mkd::static_if
                        <   modif == true,
                            make_ufunc_root<Tag, S1s>,
                            make_ufunc_impl<Tag, S1s>
                        >::type;
    using type  = typename type0 :: type;
    
    static const bool is_sdret  = mkd::is_valid_scalar_data<type>::value;
    static_assert(is_sdret == true, "type should be scalar_data");
};

//----------------------------------------------------------------------------------
//                              func_unary_impl
//----------------------------------------------------------------------------------
template<class Tag, Mat_or_scalar M1>
struct func_unary_impl
{
    static_assert(md::dependent_false<M1>::value, "M1 must be ct_matrix or ct_scalar");
};

template<class Tag, Integer M, Integer N, Mat_array Array, DPS Deps1>
struct func_unary_impl<Tag, ct_matrix<M, N, Array, Deps1>>
{
    using array_type    = mkd::mat_ufunc_array<Tag, M, N, Array>;
    using type          = ct_matrix<M, N, array_type, Deps1>;
};

template<class Tag, Scal_data Array, DPS Deps>
struct func_unary_impl<Tag, ct_scalar<Array, Deps>>
{
    using ufunc_type    = typename make_ufunc_root<Tag, Array>::type;
    using type          = ct_scalar<ufunc_type, Deps>;
};

}}}

#if 0
/*

//----------------------------------------------------------------------------------
//                              func_bin
//----------------------------------------------------------------------------------

template<class Tag,Integer M1_M2, Integer N1_N2, class Array1, class Deps1, class Array2, class Deps2>
struct func_bin<Tag, ct_matrix<M1_M2,N1_N2,Array1,Deps1>, ct_matrix<M1_M2,N1_N2,Array2,Deps2>>
{
    using array_type    = mkd::mat_bfunc_array<Tag,M1_M2, N1_N2, Array1, Array2>;
    using deps          = typename link_deps<Deps1, Deps2>::type;
    using type          = ct_matrix<M1_M2, N1_N2, array_type,deps>;
};

template<class Tag,Integer M1, Integer N1, class Array1, class Deps1, class Array2, class Deps2>
struct func_bin<Tag, ct_matrix<M1,N1,Array1,Deps1>, ct_scalar<Array2,Deps2>>
{
    using array_type    = mkd::mat_scal_bfunc_array<Tag,M1,N1,Array1,ct_scalar<Array2,Deps2>>;
    using deps          = typename link_deps<Deps1, Deps2>::type;
    using type          = ct_matrix<M1, N1, array_type,deps>;
};

template<class Tag,Integer M1, Integer N1, class Array1, class Deps1, class Array2,class Deps2>
struct func_bin<Tag, ct_scalar<Array2, Deps2>, ct_matrix<M1,N1,Array1,Deps1>>
{
    using array_type    = mkd::scal_mat_bfunc_array<Tag,M1, N1, Array1, ct_scalar<Array2,Deps2>>;
    using deps          = typename link_deps<Deps1, Deps2>::type;
    using type          = ct_matrix<M1, N1, array_type,deps>;
};

template<class Tag,class Array1, class Deps1, class Array2, class Deps2>
struct func_bin<Tag, ct_scalar<Array1,Deps1>, ct_scalar<Array2,Deps2>>
{
    using array_type    = details::scalar_bfunc_array<Tag,ct_scalar<Array1,Deps1>,ct_scalar<Array2,Deps2>>;
    using deps          = typename link_deps<Deps1, Deps2>::type;
    using type          = ct_scalar<array_type, deps>;
};

template<class Tag,Integer M1, Integer N1, Integer M2, Integer N2, class Array1, class Deps1,
        class Array2, class Deps2>
struct func_bin<Tag,ct_matrix<M1,N1,Array1,Deps1>, ct_matrix<M2,N2,Array2,Deps2>>
{
    static_assert(M1 == M2 && N1 == N2, "invalid matrix operation, check size of matrices");
};
*/
#endif
