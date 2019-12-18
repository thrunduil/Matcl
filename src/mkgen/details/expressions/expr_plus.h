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
//                              make_plus_scal
//----------------------------------------------------------------------------------
// multiply two value scalars
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
    using tag   = scal_data_value_tag_plus<Tag1, Tag2>;

    using type  = scal_data_const_value<tag, val>;
};

template<class Tag1, class Val1, Integer N2, Integer D2>
struct make_plus_scal<mkd::scal_data_const_value<Tag1,Val1>, 
                      mkd::scal_data_rational<N2, D2>>
{
    using val   = decltype(std::declval<Val1>() + std::declval<double>());
    using tag2  = scal_data_value_tag_rational<N2, D2>;
    using tag   = scal_data_value_tag_plus<Tag1, tag2>;

    using type  = scal_data_const_value<tag, val>;
};

template<class Tag1, class Val1, Integer N2, Integer D2>
struct make_plus_scal<mkd::scal_data_rational<N2, D2>, 
                      mkd::scal_data_const_value<Tag1,Val1>>
{
    using S1    = mkd::scal_data_const_value<Tag1,Val1>;
    using S2    = mkd::scal_data_rational<N2, D2>;
    using type  = typename make_plus_scal<S1, S2>::type;
};


//TODO:
//----------------------------------------------------------------------------------
//                              make_plus
//----------------------------------------------------------------------------------
template<class T1, class T2,
        bool Is_Scal_1 = is_value_scalar_data<T1>::value, 
        bool Is_Scal_2 = is_value_scalar_data<T2>::value>
struct make_plus
{
    static_assert(is_mult_expr<T1>::value == false && is_mult_expr<T2>::value == false
                  && is_value_scalar_data<T1>::value == false 
                  && is_value_scalar_data<T2>::value == false,
                  "this case should be already processed");

    using type = expr_plus_scalar_data<T1,T2>;
};

template<class S1, class ...T1, class S2, class ...T2>
struct make_plus<expr_mult_scalar_data<S1,T1...>, expr_mult_scalar_data<S2,T2...>, false, false>
{
    using ex1   = expr_mult_scalar_data<S1,T1...>;
    using ex2   = expr_mult_scalar_data<S2,T2...>;
    using type  = expr_plus_scalar_data<ex1,ex2>;
};

template<class S1, class ...T1, class T2>
struct make_plus<expr_mult_scalar_data<S1,T1...>, T2, false, false>
{
    using ex1   = expr_mult_scalar_data<S1,T1...>;
    using type  = expr_plus_scalar_data<ex1, T2>;
};

template<class T1, class S2, class ...T2>
struct make_plus<T1, expr_mult_scalar_data<S2,T2...>, false, false>
{
    using ex2   = expr_mult_scalar_data<S2,T2...>;
    using type  = expr_plus_scalar_data<T1, ex2>;
};

template<class S1, class ...T1, class S2>
struct make_plus<expr_mult_scalar_data<S1,T1...>, S2, false, true>
{
    using ex1   = expr_mult_scalar_data<S1,T1...>;

    static const bool is_zero   = mkd::is_scalar_data_zero<S2>::value;

    using type  = typename mkd::static_if
                        <   is_zero == true,
                            ex1,
                            expr_plus_scalar_data<ex1, S2>
                        >::type;
};

template<class T1, class S2>
struct make_plus<T1, S2, false, true>
{
    static const bool is_zero   = mkd::is_scalar_data_zero<S2>::value;

    using type  = typename mkd::static_if
                        <   is_zero == true,
                            T1,
                            expr_plus_scalar_data<T1, S2>
                        >::type;
};

template<class S1, class T2>
struct make_plus<S1, T2, true, false>
{
    using type  = typename make_plus<T2, S1, false, true>::type; 
};

template<class S1, class S2>
struct make_plus<S1, S2, true, true>
{
    static_assert(is_value_scalar_data<S1>::value == true 
                  && is_value_scalar_data<S2>::value == true,
                  "invalid arguments");

    using type  = typename make_plus_scal<S1, S2>::type;
};

//TODO
#if 0
template<class Scal, class ... Elems>
struct normalize_mult_scal
{
    using type = expr_mult_scalar_data<Scal,Elems...>;
};
template<class Elem>
struct normalize_mult_scal<one_sd,Elem>
{
    using type = Elem;
};

template<class T2> struct make_plus<zero_sd, T2>   { using type = T2 ; };
template<class T1> struct make_plus<T1,zero_sd>    { using type = T1; };
template<>         struct make_plus<zero_sd,zero_sd>  { using type = zero_sd; };


template<class ...T1, class T2>
struct make_plus<expr_mult_scalar_data<mone_sd,T1...>,T2>
{
    using ex1   = typename normalize_mult_scal<one_sd,T1...>::type;
    using type  = expr_minus<T2, ex1>;
};

template<class T1, class ...T2>
struct make_plus<T1, expr_mult_scalar_data<mone_sd,T2...>>
{
    using ex2   = typename normalize_mult_scal<one_sd,T2...>::type;
    using type  = expr_minus<T1, ex2>;
};

template<class ...T1, class S2, class ...T2>
struct make_plus<expr_mult_scalar_data<mone_sd,T1...>, expr_mult_scalar_data<S2,T2...>>
{
    using ex1   = typename normalize_mult_scal<one_sd,T1...>::type;
    using ex2   = expr_mult_scalar_data<S2,T2...>;
    using type  = expr_minus<ex2,ex1>;
};

template<class S1, class ...T1, class ...T2>
struct make_plus<expr_mult_scalar_data<S1,T1...>, expr_mult_scalar_data<mone_sd,T2...>>
{
    using ex1   = expr_mult_scalar_data<S1,T1...>;
    using ex2   = typename normalize_mult_scal<one_sd,T2...>::type;
    using type  = expr_minus<ex1,ex2>;
};

template<class ...T1, class ...T2>
struct make_plus<expr_mult_scalar_data<mone_sd,T1...>, expr_mult_scalar_data<mone_sd,T2...>>
{
    using ex1   = typename normalize_mult_scal<one_sd,T1...>::type;
    using ex2   = typename normalize_mult_scal<one_sd,T2...>::type;
    using ex    = expr_plus_scalar_data<ex1,ex2>;
    using type  = typename make_mult<mone_sd,ex>::type;
};

//S1*T11*T12 + S2*T21*T22
template<class S1, class T11, class T12, class S2, class T21, class T22>
struct make_plus<expr_mult_scalar_data<S1,T11,T12>, expr_mult_scalar_data<S2,T21,T22>>
{
    using ex1   = expr_mult_scalar_data<S1,T11,T12>;
    using ex2   = expr_mult_scalar_data<S2,T21,T22>;
    using type  = expr_plus_scalar_data<ex1,ex2>;
};

template<class T11, class T12, class S2, class T21, class T22>
struct make_plus<expr_mult_scalar_data<mone_sd,T11,T12>, expr_mult_scalar_data<S2,T21,T22>>
{
    using ex1   = expr_mult_scalar_data<S2,T21,T22>;
    using ex2   = expr_mult_scalar_data<one_sd,T11,T12>;    
    using type  = expr_minus<ex1,ex2>;
};

template<class S1, class T11, class T12, class T21, class T22>
struct make_plus<expr_mult_scalar_data<S1,T11,T12>, expr_mult_scalar_data<mone_sd,T21,T22>>
{
    using ex1   = expr_mult_scalar_data<S1,T11,T12>;
    using ex2   = expr_mult_scalar_data<one_sd,T21,T22>;
    using type  = expr_minus<ex1,ex2>;
};

template<class T11, class T12, class T21, class T22>
struct make_plus<expr_mult_scalar_data<mone_sd,T11,T12>, expr_mult_scalar_data<mone_sd,T21,T22>>
{
    using ex1   = expr_mult_scalar_data<one_sd,T11,T12>;
    using ex2   = expr_mult_scalar_data<one_sd,T21,T22>;
    using ex    = expr_plus_scalar_data<ex1,ex2>;
    using type  = typename make_mult<mone_sd,ex>::type;
};

//S1*T*T12 + S2*T*T22
template<class S1, class T, class T12, class S2, class T22>
struct make_plus<expr_mult_scalar_data<S1,T,T12>, expr_mult_scalar_data<S2,T,T22>>
{
    using ex1   = typename make_mult<S1,T12>::type;
    using ex2   = typename make_mult<S2,T22>::type;
    using ex    = typename make_plus<ex1,ex2>::type;
    using type  = typename make_mult<T,ex>::type;
};

template<class T, class T12, class S2, class T22>
struct make_plus<expr_mult_scalar_data<mone_sd,T,T12>, expr_mult_scalar_data<S2,T,T22>>
{
    using ex1   = typename make_mult<S2,T22>::type;
    using ex2   = T12;
    using ex    = typename make_minus<ex1,ex2>::type;
    using type  = typename make_mult<T,ex>::type;
};

template<class S1, class T, class T12, class T22>
struct make_plus<expr_mult_scalar_data<S1,T,T12>, expr_mult_scalar_data<mone_sd,T,T22>>
{
    using ex1   = typename make_mult<S1,T12>::type;
    using ex2   = T22;
    using ex    = typename make_minus<ex1,ex2>::type;
    using type  = typename make_mult<T,ex>::type;
};

template<class T, class T12, class T22>
struct make_plus<expr_mult_scalar_data<mone_sd,T,T12>, expr_mult_scalar_data<mone_sd,T,T22>>
{
    using ex3   = typename make_plus<T12,T22>::type;
    using ex    = typename make_mult<T,ex3>::type;
    using type  = typename make_mult<mone_sd,ex>::type;
};

//S1*T*T12 + S2*T21*T
template<class S1, class T, class T12, class S2, class T21>
struct make_plus<expr_mult_scalar_data<S1,T,T12>, expr_mult_scalar_data<S2,T21,T>>
{
    using ex1   = typename make_mult<S1,T12>::type;
    using ex2   = typename make_mult<S2,T21>::type;
    using ex    = typename make_plus<ex1,ex2>::type;
    using type  = typename make_mult<T,ex>::type;
};

template<class T, class T12, class S2, class T21>
struct make_plus<expr_mult_scalar_data<mone_sd,T,T12>, expr_mult_scalar_data<S2,T21,T>>
{
    using ex1   = typename make_mult<S2,T21>::type;
    using ex2   = T12;
    using ex    = typename make_minus<ex1,ex2>::type;
    using type  = typename make_mult<T,ex>::type;
};

template<class S1, class T, class T12, class T21>
struct make_plus<expr_mult_scalar_data<S1,T,T12>, expr_mult_scalar_data<mone_sd,T21,T>>
{
    using ex1   = typename make_mult<S1,T12>::type;
    using ex2   = T21;
    using ex    = typename make_minus<ex1,ex2>::type;
    using type  = typename make_mult<T,ex>::type;
};

template<class T, class T12, class T21>
struct make_plus<expr_mult_scalar_data<mone_sd,T,T12>, expr_mult_scalar_data<mone_sd,T21,T>>
{
    using ex3   = typename make_plus<T12,T21>::type;
    using ex    = typename make_mult<T,ex3>::type;
    using type  = typename make_mult<mone_sd,ex>::type;
};

//S1*T11*T + S2*T21*T
template<class S1, class T, class T11, class S2, class T21>
struct make_plus<expr_mult_scalar_data<S1,T11,T>, expr_mult_scalar_data<S2,T21,T>>
{
    using ex1   = typename make_mult<S1,T11>::type;
    using ex2   = typename make_mult<S2,T21>::type;
    using ex    = typename make_plus<ex1,ex2>::type;
    using type  = typename make_mult<T,ex>::type;
};

template<class T, class T11, class S2, class T21>
struct make_plus<expr_mult_scalar_data<mone_sd,T11,T>, expr_mult_scalar_data<S2,T21,T>>
{
    using ex1   = typename make_mult<S2,T21>::type;
    using ex2   = T11;
    using ex    = typename make_minus<ex1,ex2>::type;
    using type  = typename make_mult<T,ex>::type;
};

template<class S1, class T, class T11, class T21>
struct make_plus<expr_mult_scalar_data<S1,T11,T>, expr_mult_scalar_data<mone_sd,T21,T>>
{
    using ex1   = typename make_mult<S1,T11>::type;
    using ex2   = T21;
    using ex    = typename make_minus<ex1,ex2>::type;
    using type  = typename make_mult<T,ex>::type;
};

template<class T, class T11, class T21>
struct make_plus<expr_mult_scalar_data<mone,T11,T>, expr_mult_scalar_data<mone,T21,T>>
{
    using ex3   = typename make_plus<T11,T21>::type;
    using ex    = typename make_mult<T,ex3>::type;
    using type  = typename make_mult<mone_sd,ex>::type;
};

//S1*T11*T + S2*T*T22
template<class S1, class T, class T11, class S2, class T22>
struct make_plus<expr_mult_scalar_data<S1,T11,T>, expr_mult_scalar_data<S2,T,T22>>
{
    using ex1   = typename make_mult<S1,T11>::type;
    using ex2   = typename make_mult<S2,T22>::type;
    using ex    = typename make_plus<ex1,ex2>::type;
    using type  = typename make_mult<T,ex>::type;
};

template<class T, class T11, class S2, class T22>
struct make_plus<expr_mult_scalar_data<mone_sd,T11,T>, expr_mult_scalar_data<S2,T,T22>>
{
    using ex1   = typename make_mult<S2,T22>::type;
    using ex2   = T11;
    using ex    = typename make_minus<ex1,ex2>::type;
    using type  = typename make_mult<T,ex>::type;
};

template<class S1, class T, class T11, class T22>
struct make_plus<expr_mult_scalar_data<S1,T11,T>, expr_mult_scalar_data<mone_sd,T,T22>>
{
    using ex1   = typename make_mult<S1,T11>::type;
    using ex2   = T22;
    using ex    = typename make_minus<ex1,ex2>::type;
    using type  = typename make_mult<T,ex>::type;
};

template<class T, class T11, class T22>
struct make_plus<expr_mult_scalar_data<mone_sd,T11,T>, expr_mult_scalar_data<mone_sd,T,T22>>
{
    using ex3   = typename make_plus<T11,T22>::type;
    using ex    = typename make_mult<T,ex3>::type;
    using type  = typename make_mult<mone_sd,ex>::type;
};

template<class S1, class S2, class T11, class T21, bool Neg>
struct make_plus_2_impl
{
    using ex1   = typename normalize_mult_scal<S1,T11>::type;
    using ex2   = typename normalize_mult_scal<S2,T21>::type;
    using type  = expr_plus_scalar_data<ex1,ex2>;
};
template<class S1, class S2, class T11, class T21>
struct make_plus_2_impl<S1,S2,T11,T21,true>
{
    using ex    = typename make_minus<T11,T21>::type;
    using type  = typename make_mult<S1,ex>::type;
};
#endif

//----------------------------------------------------------------------------------
//                              make_minus
//----------------------------------------------------------------------------------

template<class T1, class T2>
struct make_minus
{
    using MT2       = typename make_mult_root<mone_sd,T2>::type;
    using type      = typename make_plus<T1,MT2>::type;
};

//----------------------------------------------------------------------------------
//                              make_minus
//----------------------------------------------------------------------------------
template<class T1>
struct make_uminus
{
    using type      = typename make_mult_root<mone_sd, T1> :: type;
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

    using type0     = typename make_plus<Scal1, Scal2>::type;

    // TODO: remove this call
    using type      = typename correct_scalar_get_elem<type0>::type;
};

// representation of  Scal1 - Scal2, where Scal1, Scal2 are scalar_data, return
// scalar_data type
template<class Scal1, class Scal2>
struct make_minus_root
{
    static const bool is_sd1    = mkd::is_valid_scalar_data<Scal1>::value;
    static const bool is_sd2    = mkd::is_valid_scalar_data<Scal2>::value;

    static_assert(is_sd1 == true && is_sd2 == true, "scalar_data required");

    using type0     = typename make_minus<Scal1, Scal2>::type;

    // TODO: remove this call
    using type      = typename correct_scalar_get_elem<type0>::type;
};

// representation of  -Scal1, where Scal1 is scalar_data, return
// scalar_data type
template<class Scal1>
struct make_uminus_root
{
    static const bool is_sd1    = mkd::is_valid_scalar_data<Scal1>::value;

    static_assert(is_sd1 == true, "scalar_data required");

    using type0     = typename make_uminus<Scal1>::type;

    // TODO: remove this call
    using type      = typename correct_scalar_get_elem<type0>::type;
};

}}};

#if 0
#include <type_traits>
#include "mkgen/matrix/scalar.h"
#include "mkgen/TODO/matrix/rational.h"
#include "mkgen/TODO/utils/utils.h"

namespace matcl { namespace mkgen
{

template<class T1, class T2>
struct make_minus;

//----------------------------------------------------------------------------------
//                              expr_plus_scalar_data
//----------------------------------------------------------------------------------


template<class T1, class T2>
struct expr_minus : public mkd::scalar_data<expr_minus<T1, T2>>
{
    template<class Subs_Context>
    static void print(std::ostream& os, int prior)
    {
        if (prior > details::prior_plus)
            os << "(";

        T1::print<Subs_Context>(os, details::prior_plus);

        os << "-";

        T2::print<Subs_Context>(os, details::prior_plus);

        if (prior > details::prior_plus)
            os << ")";
    };
    template<class Ret, class Local_Storage>
    inline_expr
    static Ret eval(const Local_Storage& ls)
    {
        Ret tmp1 = T1::eval<Ret>(ls);
        Ret tmp2 = T2::eval<Ret>(ls);
        return tmp1 - tmp2;
    };

    template<class Loop_Storage, class Ret, class Local_Storage>
    inline_expr
    static void eval_loop(Ret& ret, Integer offset, const Local_Storage& cont)
    {
        Ret val1;
        Ret val2;

        T1::eval_loop<Loop_Storage,Ret>(val1, offset,cont);
        T2::eval_loop<Loop_Storage,Ret>(val2, offset,cont);

        ret = val1 - val2;
    };

    template<Integer Step, class Arr_List>
    using get_arrays    = typename expr_minus_arrays<Step, Arr_List, T1, T2>::type;

    template<class Visitor>
    static void accept(Visitor& vis)
    {
        T1::accept<Visitor>(vis);
        T2::accept<Visitor>(vis);
        vis.visit_minus();
    };
};

template<class S1, class S2>
struct is_div_mone
{
    static_assert(md::dependent_false<S1>::value, 
                "this type should not be instantiated");
};

template<Integer M1, Integer D1, Integer M2, Integer D2, class Deps1, class Deps2>
struct is_div_mone<ct_scalar<mkd::scal_data_rational<M1,D1>, Deps1>,
                   ct_scalar<mkd::scal_data_rational<M2,D2>, Deps2>>
{
    static const bool value = (D1 == D2) && (M1 == -M2);
};

//S1*T11 + S2*T21
template<class S1, class T11, class S2, class T21>
struct make_plus<expr_mult_scalar_data<S1,T11>, expr_mult_scalar_data<S2,T21>>
    :make_plus_2_impl<S1,S2,T11,T21,is_div_mone<S1,S2>::value>
{};
template<class T11, class S2, class T21>
struct make_plus<expr_mult_scalar_data<mone_sd,T11>, expr_mult_scalar_data<S2,T21>>
{
    static_assert(mkd::is_scalar_data_one<S2>::value == false,"invalid_rep");
    using ex1   = expr_mult_scalar_data<S2,T21>;
    using ex2   = T11;    
    using type  = expr_minus<ex1,ex2>;
};
template<class S1, class T11, class T21>
struct make_plus<expr_mult_scalar_data<S1,T11>, expr_mult_scalar_data<mone_sd,T21>>
{
    static_assert(mkd::is_scalar_data_one<S1>::value == false,"invalid_rep");
    using ex1   = expr_mult_scalar_data<S1,T11>;
    using ex2   = T21;
    using type  = expr_minus<ex1,ex2>;
};
template<class T11, class T21>
struct make_plus<expr_mult_scalar_data<mone_sd,T11>, expr_mult_scalar_data<mone_sd,T21>>
{
    using ex1   = T11;
    using ex2   = T21;
    using ex    = expr_plus_scalar_data<ex1,ex2>;
    using type  = typename make_mult<mone_sd,ex>::type;
};

//S1*T + S2*T
template<class S1, class T, class S2>
struct make_plus<expr_mult_scalar_data<S1,T>, expr_mult_scalar_data<S2,T>>
{
    using ex    = typename make_plus<S1,S2>::type;
    using type  = typename make_mult<ex,T>::type;
};

//S*T1 + S*T2
template<class S, class T1, class T2>
struct make_plus<expr_mult_scalar_data<S,T1>, expr_mult_scalar_data<S,T2>>
{
    using ex    = typename make_plus<T1,T2>::type;
    using type  = typename make_mult<S,ex>::type;
};

}}

#endif