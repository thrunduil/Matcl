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

#include "mkgen/matrix/scalar.h"
#include "mkgen/details/expressions/expr_mult_scalar_data.h"
#include "mkgen/details/utils/rational.h"
#include "mkgen/details/matrix/scal_data_value_tags.h"

namespace matcl { namespace mkgen { namespace details
{

namespace mk = matcl::mkgen;

//----------------------------------------------------------------------------------
//                              forward declarations
//----------------------------------------------------------------------------------
template<bool Is_simpl, class S, class ... T>
struct make_mult_normalize_impl;

template<class S, class T>
struct make_mult_normalize_plus;

template<class S, class T>
struct make_mult_plus_impl;

template<class S, class T, bool Can_accumulate>
struct make_mult_plus_impl_acc;

// defined in expr_plus_scalar_data
template<class List_mult, class List_plus, class S, class ... PI>
struct mult_plus_items;

// defined in expr_plus_scalar_data
template<class Plus_expr>
struct can_accumulate_scalling;

template<bool Is_simpl, class S, class ... T>
struct make_expr_plus_sd;

template<bool Is_simpl, class S, class List>
struct make_expr_plus_list;

template<class S, class List>
struct make_mult_plus_add_items;

template<class Scal1, class Scal2>
struct make_plus_root;

//----------------------------------------------------------------------------------
//                              make_mult_scal
//----------------------------------------------------------------------------------
// multiply two value scalars
template<class S1, class S2>
struct make_mult_scal
{
    static_assert(md::dependent_false<S1>::value, 
                "this type should not be instantiated");
};

template<Integer N1, Integer D1, Integer N2, Integer D2>
struct make_mult_scal<scal_data_rational<N1, D1>, scal_data_rational<N2, D2>>
{
    using op    = mkd::rational_mult<N1,D1,N2,D2>;
    using type  = scal_data_rational<op::nominator,op::denominator>;
};

template<class Tag1, class Val1, class Tag2, class Val2>
struct make_mult_scal<mkd::scal_data_const_value<Tag1,Val1>, 
                      mkd::scal_data_const_value<Tag2, Val2>>
{
    using val   = decltype(std::declval<Val1>() * std::declval<Val2>());
    using tag   = scal_data_const_value_tag_mult<Tag1, Tag2>;

    using type  = scal_data_const_value<tag, val>;
};

template<class Tag1, class Val1, class Tag2, class Val2>
struct make_mult_scal<mkd::scal_data_value<Tag1,Val1>, 
                      mkd::scal_data_value<Tag2, Val2>>
{
    using val   = decltype(std::declval<Val1>() * std::declval<Val2>());
    using tag   = scal_data_value_tag_mult<Tag1, Tag2>;

    using type  = scal_data_value<tag, val>;
};

template<class Tag1, class Val1, Integer N2, Integer D2>
struct make_mult_scal<mkd::scal_data_const_value<Tag1,Val1>, 
                      mkd::scal_data_rational<N2, D2>>
{
    using tag1  = mkd::scal_data_const_value<Tag1,Val1>;
    using tag2  = mkd::scal_data_const_value
                    <scal_data_const_value_tag_rational<N2, D2>, double>;

    using type  = typename make_mult_scal<tag1, tag2> :: type;
};

template<class Tag1, class Val1, class Tag2, class Val2>
struct make_mult_scal<mkd::scal_data_const_value<Tag1,Val1>, 
                      mkd::scal_data_value<Tag2,Val2>>
{
    using tag1  = mkd::scal_data_value<scal_data_value_tag_const<Tag1, Val1>, Val1>;
    using tag2  = mkd::scal_data_value<Tag2,Val2>;

    using type  = typename make_mult_scal<tag1, tag2> :: type;
};

template<Integer N1, Integer D1, class Tag2, class Val2>
struct make_mult_scal<mkd::scal_data_rational<N1, D1>, 
                      mkd::scal_data_const_value<Tag2,Val2>>
{
    using tag1  = mkd::scal_data_const_value
                        <scal_data_const_value_tag_rational<N1, D1>, double>;
    using tag2  = mkd::scal_data_const_value<Tag2,Val2>;
    
    using type  = typename make_mult_scal<tag1, tag2> :: type;
};

template<Integer N1, Integer D1, class Tag2, class Val2>
struct make_mult_scal<mkd::scal_data_rational<N1, D1>, 
                      mkd::scal_data_value<Tag2,Val2>>
{
    using tag1  = mkd::scal_data_value
                        <scal_data_value_tag_rational<N1, D1>, double>;
    using tag2  = mkd::scal_data_value<Tag2,Val2>;
    
    using type  = typename make_mult_scal<tag1, tag2> :: type;
};

template<class Tag1, class Val1, Integer N2, Integer D2>
struct make_mult_scal<mkd::scal_data_value<Tag1,Val1>, 
                      mkd::scal_data_rational<N2, D2>>
{
    using tag1  = mkd::scal_data_value<Tag1,Val1>;
    using tag2  = mkd::scal_data_value
                    <scal_data_value_tag_rational<N2, D2>, double>;

    using type  = typename make_mult_scal<tag1, tag2> :: type;
};

template<class Tag1, class Val1, class Tag2, class Val2>
struct make_mult_scal<mkd::scal_data_value<Tag1,Val1>, 
                      mkd::scal_data_const_value<Tag2,Val2>>
{
    using tag1  = mkd::scal_data_value<Tag1, Val1>;
    using tag2  = mkd::scal_data_value<scal_data_value_tag_const<Tag2, Val2>, Val2>;
    
    using type  = typename make_mult_scal<tag1, tag2> :: type;
};

//----------------------------------------------------------------------------------
//                              make_mult_impl
//----------------------------------------------------------------------------------
template<class T1, class T2, 
        bool Is_Scal_1 = is_value_scalar_data<T1>::value, 
        bool Is_Scal_2 = is_value_scalar_data<T2>::value>
struct make_mult_impl
{
    static_assert(is_mult_expr<T1>::value == false && is_mult_expr<T2>::value == false
                  && is_value_scalar_data<T1>::value == false 
                  && is_value_scalar_data<T2>::value == false,
                  "this case should be already processed");

    using T1s       = typename T1 :: template simplify<void>; 
    using T2s       = typename T2 :: template simplify<void>;

    static const bool modif = (std::is_same<T1s, T1>::value == false)
                            || (std::is_same<T2s, T2>::value == false);

    using type0  = typename mkd::static_if
                        <   modif == true,
                            make_mult_impl<T1s, T2s>,
                            make_expr_mult_sd<false, one_sd, mult_item<T1s, 1>, mult_item<T2s, 1>>
                        >::type;
    using type  = typename type0 :: type;
};

template<bool F1, class S1, class ... T1, 
         bool F2, class S2, class ... T2>
struct make_mult_impl<expr_mult_sd<F1, S1, T1...>, 
                      expr_mult_sd<F2, S2, T2...>, false, false>
{
    using S     = typename make_mult_scal<S1, S2>::type;
    using type  = typename make_expr_mult_sd<false, S, T1..., T2...> :: type;
};

template<bool F1, class S1, class ... T1, class T2>
struct make_mult_impl<expr_mult_sd<F1, S1, T1...>, T2, false, false>
{
    static_assert(is_mult_expr<T2>::value == false
                  && is_value_scalar_data<T2>::value == false,
                  "this case should be already processed");

    using T2s       = typename T2 :: template simplify<void>;

    static const bool modif = (std::is_same<T2s, T2>::value == false);

    using type0  = typename mkd::static_if
                        <   modif == true,
                            make_mult_impl<expr_mult_sd<F1, S1, T1...>, T2s>,
                            make_expr_mult_sd<false, S1, T1 ..., mult_item<T2s, 1>>
                        >::type;
    using type  = typename type0 :: type;
};

template<class T1, bool F2, class S2, class ... T2>
struct make_mult_impl<T1, expr_mult_sd<F2, S2, T2...>, false, false>
{
    static_assert(is_mult_expr<T1>::value == false
                  && is_value_scalar_data<T1>::value == false,
                  "this case should be already processed");

    using T1s       = typename T1 :: template simplify<void>;

    static const bool modif = (std::is_same<T1s, T1>::value == false);

    using type0  = typename mkd::static_if
                        <   modif == true,
                            make_mult_impl<T1s, expr_mult_sd<F2, S2, T2 ...>>,
                            make_expr_mult_sd<false, S2, mult_item<T1s, 1>, T2 ...>
                        >::type;
    using type  = typename type0 :: type;
};

template<bool F1, class S1, class ... T1, class S2>
struct make_mult_impl<expr_mult_sd<F1, S1, T1...>, S2, false, true>
{
    static_assert(is_value_scalar_data<S2>::value == true,
                  "invalid arguments");

    using S     = typename make_mult_scal<S1, S2>::type;
    using type  = typename make_mult_normalize_impl<F1, S, T1...>::type;
};

template<class T1, class S2>
struct make_mult_impl<T1, S2, false, true>
{
    static_assert(is_mult_expr<T1>::value == false
                  && is_value_scalar_data<T1>::value == false 
                  && is_value_scalar_data<S2>::value == true,
                  "this case should be already processed");

    using T1s       = typename T1 :: template simplify<void>;

    static const bool modif = (std::is_same<T1s, T1>::value == false);

    using type0  = typename mkd::static_if
                        <   modif == true,
                            make_mult_impl<T1s, S2>,
                            make_mult_normalize_impl<true, S2, mult_item<T1s, 1>>
                        >::type;
    using type  = typename type0 :: type;
};

template<class S1, class T2>
struct make_mult_impl<S1, T2, true, false>
{
    using type  = typename make_mult_impl<T2, S1, false, true>::type;
};

template<class S1, class S2>
struct make_mult_impl<S1, S2, true, true>
{
    static_assert(is_value_scalar_data<S1>::value == true 
                  && is_value_scalar_data<S2>::value == true,
                  "invalid arguments");

    using type  = typename make_mult_scal<S1, S2>::type;
};

//----------------------------------------------------------------------------------
//                              make_mult_normalize_impl
//----------------------------------------------------------------------------------
// make S * T1 * ... * Tk, where S is a value scalar_data and Ti are nonvalue scalar data
// if S == 0, then return 0, if S == 1 and k = 1, then return T1
// for all types Ti, Ti::simplify<> should already be called; can return plus expression
// If Is_simpl == true, then expr_mult<S, T...> was simplified previously
template<bool Is_simpl, class S, class ... T>
struct make_mult_normalize_impl
{
    static const bool is_zero   = mkd::is_scalar_data_zero<S>::value;

    using type0  = typename mkd::static_if
                        <   is_zero == true,
                            mkd::lazy_type<zero_sd>,
                            make_expr_mult_sd<Is_simpl, S, T...>
                        >::type;
    using type  = typename type0 :: type;
};

template<bool Is_simpl, class S, class It>
struct make_mult_normalize_impl<Is_simpl, S, It>
{
    static const bool is_zero   = mkd::is_scalar_data_zero<S>::value;
    static const bool is_one    = mkd::is_scalar_data_one<S>::value;
    static const Integer K      = It::exponent;
    using T                     = typename It::base; 

    using type0  = typename mkd::static_if
                        <   is_zero == true,
                            mkd::lazy_type<zero_sd>,
                            typename mkd::static_if
                                <   is_one == true && K == 1,
                                    mkd::lazy_type<T>,
                                    make_mult_normalize_plus<S, It>
                                >::type
                        >::type;

    using type  = typename type0 :: type;
};

//----------------------------------------------------------------------------------
//                              make_mult_normalize_plus
//----------------------------------------------------------------------------------
// make S * (T1)^k, where S is a value scalar_data and T1 is plus expr
// convert this to plus_item, when k = 1 and T1 can accumulate constant
// (i.e. all scalling of plus_items are nontrivial)
// T1::simplify<> should be already called; may return plus expr
template<class S, class It>
struct make_mult_normalize_plus
{
    static const bool is_zero   = mkd::is_scalar_data_zero<S>::value;
    static const bool is_one    = mkd::is_scalar_data_one<S>::value;

    static_assert(is_zero == false && is_one == false, "this case should already be processed");

    static const Integer K      = It::exponent;
    using T                     = typename It::base; 

    static const bool is_add    = mkd::is_plus_expr<T>::value;

    using type0  = typename mkd::static_if
                        <   is_add == false || K != 1,
                            make_expr_mult_sd<true, S, It>,
                            make_mult_plus_impl<S, T>
                        >::type;

    using type  = typename type0 :: type;
};

//----------------------------------------------------------------------------------
//                              make_mult_plus_impl
//----------------------------------------------------------------------------------
// make S * T1, where S is a value scalar_data and T1 is plus expr
// T1::simplify must already be called; may return plus expression
template<class S, class T>
struct make_mult_plus_impl
{
    static_assert(md::dependent_false<S>::value, "plus_expr required");
};

// expr_plus must have SP != 0, or SP == 0 and sizeof(TP) > 1
// expr_plus::simplify must already be called
template<class S, bool F, class SP, class ... TP>
struct make_mult_plus_impl<S, expr_plus_sd<F, SP, TP...>>
{
    using plus_expr             = expr_plus_sd<F, SP, TP...>;
    static const bool can_acc   = can_accumulate_scalling<plus_expr>::value;

    using type  = typename make_mult_plus_impl_acc<S, plus_expr, can_acc>::type;
};

// Plus_expr::simplify<> must already be called;
// can return plus expression
template<class S, class Plus_expr, bool Can_accumulate>
struct make_mult_plus_impl_acc
{
    using type = typename make_expr_mult_sd<true, S, mult_item<Plus_expr, 1>> :: type;
};

// expr_plus_sd<...>::simplify<> must already be called
template<class S, bool F, class SP, class ... TP>
struct make_mult_plus_impl_acc<S, expr_plus_sd<F, SP, TP...>, true>
{
    static_assert(F == true, "simplify was not called");

    // expr_plus must have SP != 0, or SP == 0 and sizeof(TP) > 1, also S != 0
    // therefore resulting expression remains expr_plus_sd and no further simplifications
    // can occur (ignoring possible underflows)
    using plus_items    = list::list<>;
    using mult_items    = list::list<>;

    using expanded      = typename mult_plus_items<mult_items, plus_items, S, TP ...>::type;

    using mult_items_ret= typename list::elem_at_pos<expanded, 0>::type;
    using plus_items_ret= typename list::elem_at_pos<expanded, 1>::type;    

    using S_ret         = typename make_mult_scal<S, SP>::type;
    using type1         = typename make_expr_plus_list<F, S_ret, mult_items_ret> ::type;

    using type          = typename make_mult_plus_add_items<type1, plus_items_ret>::type;
};

//----------------------------------------------------------------------------------
//                              make_mult_plus_add_items
//----------------------------------------------------------------------------------
template<class Type, class List>
struct make_mult_plus_add_items
{
    static_assert(md::dependent_false<Type>::value, "list::list<...> expected");
};

template<class Type, class It1, class ... It>
struct make_mult_plus_add_items<Type, list::list<It1, It...>>
{
    using type1 = typename make_plus_root<Type, It1>::type;
    using type  = typename make_mult_plus_add_items<type1, list::list<It...>>::type;
};

template<class Type>
struct make_mult_plus_add_items<Type, list::list<>>
{
    using type = Type;
};

//----------------------------------------------------------------------------------
//                              make_inv_scal
//----------------------------------------------------------------------------------
// inverse of a scalar
template<class S1>
struct make_inv_scal
{
    static_assert(md::dependent_false<S1>::value, 
                "this type should not be instantiated");
};

template<Integer N1, Integer D1>
struct make_inv_scal<scal_data_rational<N1, D1>>
{
    static_assert(N1 != 0, "inversion of zero");
    static const Integer scal   = N1 > 0 ? 1 : -1;
    using type  = scal_data_rational<D1 * scal, N1 * scal>;
};

template<class Tag1, class Val1>
struct make_inv_scal<mkd::scal_data_const_value<Tag1,Val1>>
{
    using tag   = scal_data_const_value_tag_inv<Tag1>;
    using type  = scal_data_const_value<tag, Val1>;
};

template<class Tag1, class Val1>
struct make_inv_scal<mkd::scal_data_value<Tag1,Val1>>
{
    using tag   = scal_data_value_tag_inv<Tag1>;
    using type  = scal_data_value<tag, Val1>;
};

//----------------------------------------------------------------------------------
//                              make_inv_impl
//----------------------------------------------------------------------------------
template<class T1, bool Is_Scal_1 = is_value_scalar_data<T1>::value>
struct make_inv_impl
{
    static_assert(md::dependent_false<S1>::value, 
                "this type should not be instantiated");
};

template<class T1>
struct make_inv_impl<T1, false>
{
    using T1s       = typename T1 :: template simplify<void>;

    static const bool modif = (std::is_same<T1s, T1>::value == false);

    using type0  = typename mkd::static_if
                        <   modif == true,
                            make_inv_impl<T1s>,
                            make_expr_mult_sd<false, one_sd, mult_item<T1s, -1>>
                        >::type;
    using type  = typename type0 :: type;
};

template<bool F1, class S1, class ... T1>
struct make_inv_impl<expr_mult_sd<F1, S1, T1...>, false>
{
    using S1_inv    = typename make_inv_impl<S1, true>::type;
    using type      = typename make_expr_mult_sd<false, S1_inv, 
                                typename make_inv_impl<T1, false>::type ...> :: type;
};

template<class T>
struct make_inv_impl<T, true>
{
    using type = typename make_inv_scal<T>::type;
};

//----------------------------------------------------------------------------------
//                              make_div_impl
//----------------------------------------------------------------------------------
template<class T1, class T2>
struct make_div_impl
{
    static const bool is_val_2  = is_value_scalar_data<T2>::value;

    using T2_inv    = typename make_inv_impl<T2, is_val_2>::type;
    using type      = typename make_mult_impl<T1, T2_inv>::type;
};

//----------------------------------------------------------------------------------
//                              roots
//----------------------------------------------------------------------------------

// only these templates should be used 

// representation of  Scal1 x Scal2, where Scal1, Scal2 are scalar_data, return
// scalar_data type
template<class Scal1, class Scal2>
struct make_mult_root
{
    static const bool is_sd1    = mkd::is_valid_scalar_data<Scal1>::value;
    static const bool is_sd2    = mkd::is_valid_scalar_data<Scal2>::value;

    static_assert(is_sd1 == true && is_sd2 == true, "scalar_data required");

    static const bool is_val_1  = is_value_scalar_data<Scal1>::value;
    static const bool is_val_2  = is_value_scalar_data<Scal2>::value;

    using type                  = typename make_mult_impl<Scal1, Scal2, is_val_1, is_val_2>::type;

    static const bool is_sdret  = mkd::is_valid_scalar_data<type>::value;
    static_assert(is_sdret == true, "type should be scalar_data");
};

// representation of  Scal1 / Scal2, where Scal1, Scal2 are scalar_data, return
// scalar_data type
template<class Scal1, class Scal2>
struct make_div_root
{
    static const bool is_sd1    = mkd::is_valid_scalar_data<Scal1>::value;
    static const bool is_sd2    = mkd::is_valid_scalar_data<Scal2>::value;

    static_assert(is_sd1 == true && is_sd2 == true, "scalar_data required");

    using type                  = typename make_div_impl<Scal1, Scal2>::type;

    static const bool is_sdret  = mkd::is_valid_scalar_data<type>::value;
    static_assert(is_sdret == true, "type should be scalar_data");
};

}}}
