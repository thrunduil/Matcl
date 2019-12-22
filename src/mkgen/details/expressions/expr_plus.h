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
#include "mkgen/details/expressions/expr_plus_scalar_data.h"

namespace matcl { namespace mkgen { namespace details
{

//----------------------------------------------------------------------------------
//                              forward declarations
//----------------------------------------------------------------------------------
template<bool F, class S, class ... T>
struct make_plus_normalize_impl;

//----------------------------------------------------------------------------------
//                              make_plus_scal
//----------------------------------------------------------------------------------
// add two value scalars
template<class S1, class S2>
struct make_plus_scal
{
    static_assert(md::dependent_false<S1>::value, 
                "this type should not be instantiated");
};

template<Integer N1, Integer D1, Integer N2, Integer D2>
struct make_plus_scal<scal_data_rational<N1, D1>, scal_data_rational<N2, D2>>
{
    using op    = mkd::rational_plus<N1,D1,N2,D2>;
    using type  = scal_data_rational<op::nominator,op::denominator>;
};

template<class Tag1, class Val1, class Tag2, class Val2>
struct make_plus_scal<mkd::scal_data_const_value<Tag1,Val1>, 
                      mkd::scal_data_const_value<Tag2, Val2>>
{
    using val   = decltype(std::declval<Val1>() + std::declval<Val2>());
    using tag   = scal_data_const_value_tag_plus<Tag1, Tag2>;

    using type  = scal_data_const_value<tag, val>;
};

template<class Tag1, class Val1, class Tag2, class Val2>
struct make_plus_scal<mkd::scal_data_value<Tag1,Val1>, 
                      mkd::scal_data_value<Tag2, Val2>>
{
    using val   = decltype(std::declval<Val1>() + std::declval<Val2>());
    using tag   = scal_data_value_tag_plus<Tag1, Tag2>;

    using type  = scal_data_value<tag, val>;
};

template<class Tag1, class Val1, Integer N2, Integer D2>
struct make_plus_scal<mkd::scal_data_const_value<Tag1,Val1>, 
                      mkd::scal_data_rational<N2, D2>>
{
    using tag1  = mkd::scal_data_const_value<Tag1,Val1>;
    using tag2  = mkd::scal_data_const_value
                    <scal_data_const_value_tag_rational<N2, D2>, double>;

    using type  = typename make_plus_scal<tag1, tag2> :: type;
};

template<class Tag1, class Val1, class Tag2, class Val2>
struct make_plus_scal<mkd::scal_data_const_value<Tag1,Val1>, 
                      mkd::scal_data_value<Tag2,Val2>>
{
    using tag1  = mkd::scal_data_value<scal_data_value_tag_const<Tag1, Val1>, Val1>;
    using tag2  = mkd::scal_data_value<Tag2,Val2>;

    using type  = typename make_plus_scal<tag1, tag2> :: type;
};

template<Integer N1, Integer D1, class Tag2, class Val2>
struct make_plus_scal<mkd::scal_data_rational<N1, D1>, 
                      mkd::scal_data_const_value<Tag2,Val2>>
{
    using tag1  = mkd::scal_data_const_value
                        <scal_data_const_value_tag_rational<N1, D1>, double>;
    using tag2  = mkd::scal_data_const_value<Tag2,Val2>;
    
    using type  = typename make_plus_scal<tag1, tag2> :: type;
};

template<Integer N1, Integer D1, class Tag2, class Val2>
struct make_plus_scal<mkd::scal_data_rational<N1, D1>, 
                      mkd::scal_data_value<Tag2,Val2>>
{
    using tag1  = mkd::scal_data_value
                        <scal_data_value_tag_rational<N1, D1>, double>;
    using tag2  = mkd::scal_data_value<Tag2,Val2>;
    
    using type  = typename make_plus_scal<tag1, tag2> :: type;
};

template<class Tag1, class Val1, Integer N2, Integer D2>
struct make_plus_scal<mkd::scal_data_value<Tag1,Val1>, 
                      mkd::scal_data_rational<N2, D2>>
{
    using tag1  = mkd::scal_data_value<Tag1,Val1>;
    using tag2  = mkd::scal_data_value
                    <scal_data_value_tag_rational<N2, D2>, double>;

    using type  = typename make_plus_scal<tag1, tag2> :: type;
};

template<class Tag1, class Val1, class Tag2, class Val2>
struct make_plus_scal<mkd::scal_data_value<Tag1,Val1>, 
                      mkd::scal_data_const_value<Tag2,Val2>>
{
    using tag1  = mkd::scal_data_value<Tag1, Val1>;
    using tag2  = mkd::scal_data_value<scal_data_value_tag_const<Tag2, Val2>, Val2>;
    
    using type  = typename make_plus_scal<tag1, tag2> :: type;
};

//----------------------------------------------------------------------------------
//                              make_plus_impl
//----------------------------------------------------------------------------------
template<class T1, class T2,
        bool Is_Scal_1 = is_value_scalar_data<T1>::value, 
        bool Is_Scal_2 = is_value_scalar_data<T2>::value>
struct make_plus_impl
{
    static_assert(is_plus_expr<T1>::value == false && is_plus_expr<T2>::value == false
                  && is_value_scalar_data<T1>::value == false 
                  && is_value_scalar_data<T2>::value == false,
                  "this case should be already processed");

    using T1s       = typename T1 :: template simplify<void>; 
    using T2s       = typename T2 :: template simplify<void>;

    static const bool modif = (std::is_same<T1s, T1>::value == false)
                            || (std::is_same<T2s, T2>::value == false);

    using type0  = typename mkd::static_if
                        <   modif == true,
                            make_plus_impl<T1s, T2s>,
                            make_expr_plus_sd<false, zero_sd, 
                                typename make_plus_item<T1s> :: type, 
                                typename make_plus_item<T2s> :: type>
                        >::type;
    using type  = typename type0 :: type;
};

template<bool F1, class S1, class ...T1, 
         bool F2, class S2, class ...T2>
struct make_plus_impl<expr_plus_sd<F1, S1, T1...>, expr_plus_sd<F2, S2, T2...>, false, false>
{
    using S     = typename make_plus_scal<S1, S2>::type;
    using type  = typename make_expr_plus_sd<false, S, T1..., T2 ...> :: type;
};

template<bool F1, class S1, class ...T1, class T2>
struct make_plus_impl<expr_plus_sd<F1, S1,T1...>, T2, false, false>
{
    static_assert(is_plus_expr<T2>::value == false
                  && is_value_scalar_data<T2>::value == false,
                  "this case should be already processed");

    using T2s       = typename T2 :: template simplify<void>;

    static const bool modif = (std::is_same<T2s, T2>::value == false);

    using type0  = typename mkd::static_if
                        <   modif == true,
                            make_plus_impl<expr_plus_sd<F1, S1, T1...>, T2s>,
                            make_expr_plus_sd<false, S1, T1..., 
                                typename make_plus_item<T2s> :: type>
                        >::type;
    using type  = typename type0 :: type;
};

template<class T1, bool F2, class S2, class ...T2>
struct make_plus_impl<T1, expr_plus_sd<F2, S2, T2...>, false, false>
{
    static_assert(is_plus_expr<T1>::value == false
                  && is_value_scalar_data<T1>::value == false,
                  "this case should be already processed");

    using T1s       = typename T1 :: template simplify<void>;

    static const bool modif = (std::is_same<T1s, T1>::value == false);

    using type0  = typename mkd::static_if
                        <   modif == true,
                            make_plus_impl<T1s, expr_plus_sd<F2, S2, T2 ...>>,
                            make_expr_plus_sd<false, S2, typename make_plus_item<T1> :: type, T2 ...>
                        >::type;
    using type  = typename type0 :: type;
};

template<bool F1, class S1, class ...T1, class S2>
struct make_plus_impl<expr_plus_sd<F1, S1, T1...>, S2, false, true>
{
    static_assert(is_value_scalar_data<S2>::value == true,
                  "invalid arguments");

    using S     = typename make_plus_scal<S1, S2>::type;
    using type  = typename make_plus_normalize_impl<F1, S, T1...>::type;
};

template<class T1, class S2>
struct make_plus_impl<T1, S2, false, true>
{
    static_assert(is_plus_expr<T1>::value == false
                  && is_value_scalar_data<T1>::value == false 
                  && is_value_scalar_data<S2>::value == true,
                  "this case should be already processed");

    using T1s       = typename T1 :: template simplify<void>;

    static const bool modif = (std::is_same<T1s, T1>::value == false);

    using type0  = typename mkd::static_if
                        <   modif == true,
                            make_plus_impl<T1s, S2>,
                            make_plus_normalize_impl<true, S2, typename make_plus_item<T1s>::type>
                        >::type;
    using type  = typename type0 :: type;
};

template<class S1, class T2>
struct make_plus_impl<S1, T2, true, false>
{
    using type  = typename make_plus_impl<T2, S1, false, true>::type; 
};

template<class S1, class S2>
struct make_plus_impl<S1, S2, true, true>
{
    static_assert(is_value_scalar_data<S1>::value == true 
                  && is_value_scalar_data<S2>::value == true,
                  "invalid arguments");

    using type  = typename make_plus_scal<S1, S2>::type;
};

//----------------------------------------------------------------------------------
//                              make_minus_impl
//----------------------------------------------------------------------------------

template<class T1, class T2>
struct make_minus_impl
{
    using MT2       = typename make_mult_root<mone_sd,T2>::type;
    using type      = typename make_plus_impl<T1,MT2>::type;
};

//----------------------------------------------------------------------------------
//                              make_uminus_impl
//----------------------------------------------------------------------------------
template<class T1>
struct make_uminus_impl
{
    using type      = typename make_mult_root<mone_sd, T1> :: type;
};

//----------------------------------------------------------------------------------
//                              make_plus_normalize_impl
//----------------------------------------------------------------------------------
// make S + T1 + ... + Tk, where S is a value scalar_data and Ti are nonvalue scalar data
// if S == 0 and k = 1, then return T1; 
// for all types Ti, Ti::simplify<> should already be called
// If Is_simpl == true, then expr_plus<S, T...> was simplified previously
template<bool Is_simpl, class S, class ... T>
struct make_plus_normalize_impl
{
    using type  = typename make_expr_plus_sd<Is_simpl, S, T...> :: type;
};

template<bool Is_simpl, class S, class T>
struct make_plus_normalize_impl<Is_simpl, S, T>
{
    static const bool is_zero   = mkd::is_scalar_data_zero<S>::value;

    using type0  = typename mkd::static_if
                        <   is_zero == true,
                            plus_item_to_mult<T>,
                            make_expr_plus_sd<true, S, T>
                        >::type;

    using type  = typename type0 :: type;
};

//----------------------------------------------------------------------------------
//                              roots
//----------------------------------------------------------------------------------

// only these templates should be used 

// representation of  Scal1 + Scal2, where Scal1, Scal2 are scalar_data, return
// scalar_data type
template<class Scal1, class Scal2>
struct make_plus_root
{
    static const bool is_sd1    = mkd::is_valid_scalar_data<Scal1>::value;
    static const bool is_sd2    = mkd::is_valid_scalar_data<Scal2>::value;

    static_assert(is_sd1 == true && is_sd2 == true, "scalar_data required");

    using type      = typename make_plus_impl<Scal1, Scal2>::type;

    static const bool is_sdret  = mkd::is_valid_scalar_data<type>::value;
    static_assert(is_sdret == true, "type should be scalar_data");
};

// representation of  Scal1 - Scal2, where Scal1, Scal2 are scalar_data, return
// scalar_data type
template<class Scal1, class Scal2>
struct make_minus_root
{
    static const bool is_sd1    = mkd::is_valid_scalar_data<Scal1>::value;
    static const bool is_sd2    = mkd::is_valid_scalar_data<Scal2>::value;

    static_assert(is_sd1 == true && is_sd2 == true, "scalar_data required");

    using type      = typename make_minus_impl<Scal1, Scal2>::type;

    static const bool is_sdret  = mkd::is_valid_scalar_data<type>::value;
    static_assert(is_sdret == true, "type should be scalar_data");
};

// representation of  -Scal1, where Scal1 is scalar_data, return
// scalar_data type
template<class Scal1>
struct make_uminus_root
{
    static const bool is_sd1    = mkd::is_valid_scalar_data<Scal1>::value;

    static_assert(is_sd1 == true, "scalar_data required");

    using type      = typename make_uminus_impl<Scal1>::type;

    static const bool is_sdret  = mkd::is_valid_scalar_data<type>::value;
    static_assert(is_sdret == true, "type should be scalar_data");
};

}}};