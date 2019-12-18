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
template<class S, class ... T>
struct make_mult_normalize_impl;

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
    using tag   = scal_data_value_tag_mult<Tag1, Tag2>;

    using type  = scal_data_const_value<tag, val>;
};

template<class Tag1, class Val1, Integer N2, Integer D2>
struct make_mult_scal<mkd::scal_data_const_value<Tag1,Val1>, 
                      mkd::scal_data_rational<N2, D2>>
{
    using val   = decltype(std::declval<Val1>() * std::declval<double>());
    using tag2  = scal_data_value_tag_rational<N2, D2>;
    using tag   = scal_data_value_tag_mult<Tag1, tag2>;

    using type  = scal_data_const_value<tag, val>;
};

template<class Tag1, class Val1, Integer N2, Integer D2>
struct make_mult_scal<mkd::scal_data_rational<N2, D2>, 
                      mkd::scal_data_const_value<Tag1,Val1>>
{
    using S1    = mkd::scal_data_const_value<Tag1,Val1>;
    using S2    = mkd::scal_data_rational<N2, D2>;
    using type  = typename make_mult_scal<S1, S2>::type;
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

    using it1       = mult_item<T1, 1>;
    using it2       = mult_item<T2, 1>;
    using type      = expr_mult_scalar_data<one_sd, it1, it2>;
};

template<class S1, class ... T1, class S2, class ... T2>
struct make_mult_impl<expr_mult_scalar_data<S1, T1...>, 
                      expr_mult_scalar_data<S2, T2...>, false, false>
{
    using S     = typename make_mult_scal<S1, S2>::type;
    using type  = expr_mult_scalar_data<S, T1..., T2...>;
};

template<class S1, class ... T1, class T2>
struct make_mult_impl<expr_mult_scalar_data<S1, T1...>, T2, false, false>
{
    static_assert(is_mult_expr<T2>::value == false
                  && is_value_scalar_data<T2>::value == false,
                  "this case should be already processed");

    using it2       = mult_item<T2, 1>;
    using type = expr_mult_scalar_data<S1, T1..., it2>;
};

template<class T1, class S2, class ... T2>
struct make_mult_impl<T1, expr_mult_scalar_data<S2, T2...>, false, false>
{
    static_assert(is_mult_expr<T1>::value == false
                  && is_value_scalar_data<T1>::value == false,
                  "this case should be already processed");

    using it1       = mult_item<T1, 1>;
    using type = expr_mult_scalar_data<S2, T1, T2...>;
};

template<class S1, class ... T1, class S2>
struct make_mult_impl<expr_mult_scalar_data<S1, T1...>, S2, false, true>
{
    static_assert(is_value_scalar_data<S2>::value == true,
                  "invalid arguments");

    using S     = typename make_mult_scal<S1, S2>::type;
    using type  = typename make_mult_normalize_impl<S, T1...>::type;
};

template<class T1, class S2>
struct make_mult_impl<T1, S2, false, true>
{
    static_assert(is_mult_expr<T1>::value == false
                  && is_value_scalar_data<T1>::value == false 
                  && is_value_scalar_data<S2>::value == true,
                  "this case should be already processed");

    using it1   = mult_item<T1, 1>;
    using type  = typename make_mult_normalize_impl<S2, it1>::type;
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
template<class S, class ... T>
struct make_mult_normalize_impl
{
    static const bool is_zero   = mkd::is_scalar_data_zero<S>::value;

    using type  = typename mkd::static_if
                        <   is_zero == true,
                            zero_sd,
                            expr_mult_scalar_data<S, T...>
                        >::type;
};

template<class S, class It>
struct make_mult_normalize_impl<S, It>
{
    static const bool is_zero   = mkd::is_scalar_data_zero<S>::value;
    static const bool is_one    = mkd::is_scalar_data_one<S>::value;
    static const Integer K      = It::exponent;
    using T                     = typename It::base; 

    using type  = typename mkd::static_if
                        <   is_zero == true,
                            zero_sd,
                            typename mkd::static_if
                                <   is_one == true && K == 1,
                                    T,
                                    expr_mult_scalar_data<S, It>
                                >::type
                        >::type;
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
    using tag   = scal_data_value_tag_inv<Tag1>;
    using type  = scal_data_const_value<tag, Val1>;
};

//----------------------------------------------------------------------------------
//                              make_inv_impl
//----------------------------------------------------------------------------------
template<class T1, bool Is_Scal_1>
struct make_inv_impl
{
    static_assert(md::dependent_false<S1>::value, 
                "this type should not be instantiated");
};

template<class T1>
struct make_inv_impl<T1, false>
{
    using type      = expr_mult_scalar_data<one_sd, mult_item<T1, -1>>;
};

template<class S1, class ... T1>
struct make_inv_impl<expr_mult_scalar_data<S1, T1...>, false>
{
    using S1_inv    = typename make_inv_impl<S1, true>::type;
    using type      = expr_mult_scalar_data<S1_inv, 
                        typename make_inv_impl<T1, false>::type ...>;
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
