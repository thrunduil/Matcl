/*
 *  This file is a part of Matrix Computation Library (MATCL)
 *
 *  Copyright (c) Pawe³ Kowal 2019 - 2021
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

namespace matcl { namespace mkgen { namespace details
{

//----------------------------------------------------------------------------------
//                              forward declarations
//----------------------------------------------------------------------------------

template<Integer Step, class Arr_List, class ...Mult_Items>
struct expr_mult_arrays;

template<class S, class ... T>
struct can_simplify_mult;

template<Scal_data T, Integer K>
struct mult_item;

template<class Expr>
struct simplify_expr_mult;

//----------------------------------------------------------------------------------
//                              checks
//----------------------------------------------------------------------------------

template<class T>
struct check_mult_item
{
    static_assert(md::dependent_false<T>::value, "mult_item required");
};

template<Scal_data T, Integer K>
struct check_mult_item<mult_item<T, K>>
{    
    using type = void;
};

//----------------------------------------------------------------------------------
//                              mult_item
//----------------------------------------------------------------------------------
template<Scal_data T, Integer K>
struct mult_item
{
    static const bool is_vs = mkd::is_value_scalar_data<T>::value;
    static_assert(is_vs == false, "value scalar unexpected");

    static_assert(K != 0, "invalid exponent");

    static const bool is_mult   = mkd::is_mult_expr<T>::value;
    static_assert(is_mult == false, "mult expression not allowed");

    static const bool is_simpl  = T::is_simplified();
    static_assert(is_simpl == true, "subexpression not simplified");

    static const Integer exponent   = K;
    using base                      = T;

    template<Integer Step, class Arr_List>
    using get_arrays                = typename T::template get_arrays<Step, Arr_List>;

    template<class Subs_Context>
    static void print(std::ostream& os, int prior)
    {
        if constexpr(K == 1)
        {
            T::print<Subs_Context>(os, prior);
        }
        else
        {
            if (prior > details::prior_pow)
                os << "(";

            T::print<Subs_Context>(os, details::prior_pow);

            os << "^" << K;

            if (prior > details::prior_pow)
                os << ")";
        };
    };

    template<class Ret, class Local_Storage>
    inline_expr
    static Ret eval(const Local_Storage& ls)
    {
        //TODO: powers
        return T::eval<Ret>(ls);
    };

    template<class Loop_Storage, class Ret, class Local_Storage>
    inline_expr
    static void eval_loop(Ret& ret, Integer offset, const Local_Storage& cont)
    {
        //TODO: powers
        Ret val1;
        T::eval_loop<Loop_Storage,Ret>(val1,offset,cont);

        ret = val1;
    };

    template<class Visitor>
    static void accept(Visitor& vis)
    {
        //TODO
        T::accept<Visitor>(vis);
        //vis.visit_mult();
    };
};

//----------------------------------------------------------------------------------
//                              expr_mult_tail
//----------------------------------------------------------------------------------
template<class T1, class ... T>
struct expr_mult_tail
{
    static const Integer size   = 1 + sizeof...(T);

    template<class Subs_Context>
    static void print(std::ostream& os, int prior)
    {
        if (prior > details::prior_mult)
            os << "(";

        T1::print<Subs_Context>(os, details::prior_mult);                

        if constexpr(size > 1)
        {
            using tail  = expr_mult_tail<T...>;

            os << "*";
            tail::print<Subs_Context>(os, details::prior_mult);
        };

        if (prior > details::prior_mult)
            os << ")";
    };

    template<class Ret, class Local_Storage>
    inline_expr
    static Ret eval(const Local_Storage& ls)
    {
        //TODO: optimize this
        if constexpr(size > 1)
        {
            using tail  = expr_mult_tail<T...>;
            return T1::eval<Ret>(ls) * tail::eval<Ret>(ls);
        }
        else
        {
            return T1::eval<Ret>(ls);
        };
    };

    template<class Loop_Storage, class Ret, class Local_Storage>
    inline_expr
    static void eval_loop(Ret& ret, Integer offset, const Local_Storage& cont)
    {
        //TODO: optimize this
        if constexpr(size > 1)
        {
            using tail  = expr_mult_tail<T...>;

            Ret val1;
            Ret val2;

            T1::eval_loop<Loop_Storage,Ret>(val1, offset,cont);            
            tail::eval_loop<Loop_Storage,Ret>(val2, offset,cont);

            ret = val1 * val2;
        }
        else
        {
            T1::eval_loop<Loop_Storage,Ret>(ret, offset, cont);
        }
    };

    template<class Visitor>
    static void accept(Visitor& vis)
    {
        T1::accept<Visitor>(vis);
        
        if constexpr(size > 1)
        {
            using tail  = expr_mult_tail<T...>;        
            tail::accept<Visitor>(vis);
        };

        vis.visit_mult();
    };
};

//----------------------------------------------------------------------------------
//                              expr_mult_sd
//----------------------------------------------------------------------------------

// representation of mult expression as S * T1^i1 * T2^i2 * ...
// where S is a value scalar, and T1, ... are scalar_data and ik is integer
// all items T should already be simplified, but this expression can be simplified
template<bool Flag, class S, class ... T>
struct expr_mult_sd : public mkd::scalar_data<expr_mult_sd<Flag, S, T...>>
{    
    static const bool is_vs = mkd::is_value_scalar_data<S>::value;
    static_assert(is_vs == true, "value scalar reqired");

    static_assert(sizeof...(T) > 0, "invalid rep");

    // checks    
    using ch1               = list::list<typename check_mult_item<T>::type ...>;

    static const bool simp  = can_simplify_mult<S, T...>::value;
    static_assert(simp == false, "invalid rep");

    using this_type         = expr_mult_sd<Flag, S, T...>;
    using scalling          = S;
    using base              = expr_mult_tail<T...>;

    template<class Subs_Context>
    static void print(std::ostream& os, int prior)
    {
        if (prior > details::prior_mult)
            os << "(";

        if constexpr(mkd::is_scalar_data_one<S>::value == false)
        {
            if constexpr(mkd::is_scalar_data_mone<S>::value == true)
            {
                os << "-";
            }
            else
            {
                S::print<Subs_Context>(os, details::prior_mult);
                os << "*";
            };
        };        

        base::print<Subs_Context>(os, details::prior_mult);

        if (prior > details::prior_mult)
            os << ")";
    };

    template<class Ret, class Local_Storage>
    inline_expr
    static Ret eval(const Local_Storage& ls)
    {
        //TODO: optimize this

        if constexpr(mkd::is_scalar_data_one<S>::value == true)
            return base::eval<Ret>(ls);
        else
            return S::eval<Ret>(ls) * base::eval<Ret>(ls);
    };

    template<class Loop_Storage, class Ret, class Local_Storage>
    inline_expr
    static void eval_loop(Ret& ret, Integer offset, const Local_Storage& cont)
    {
        //TODO: optimize this

        if constexpr(mkd::is_scalar_data_one<S>::value == true)
        {            
            return base::eval_loop<Loop_Storage, Ret>(ret, offset,cont);
        }
        else
        {
            Ret val1;
            Ret val2;

            S::eval_loop<Loop_Storage,Ret>(val1, offset,cont);            
            base::eval_loop<Loop_Storage,Ret>(val2, offset,cont);

            ret = val1 * val2;
        };
    };

    template<Integer Step, class Arr_List>
    using get_arrays    = typename expr_mult_arrays<Step, Arr_List, T...> :: type;

    template<class Void>
    using simplify      = typename simplify_expr_mult<this_type>::type;;

    static constexpr bool is_simplified()   { return Flag; };

    template<class Visitor>
    static void accept(Visitor& vis)
    {
        if constexpr(mkd::is_scalar_data_one<S>::value == true)
        {
            base::accept<Visitor>(vis);
        }
        else
        {
            S::accept<Visitor>(vis);
            base::accept<Visitor>(vis);
            vis.visit_mult();
        };
    };
};

//----------------------------------------------------------------------------------
//                              make_expr_mult_sd
//----------------------------------------------------------------------------------
// it is assumed that T::simplify<void> was called
// if Is_simpl == true, then also expr_mult was simplified previously
template<bool Is_simpl, class S, class ... T>
struct make_expr_mult_sd
{
    using type = expr_mult_sd<false, S, T...>;
};

template<class S, class ... T>
struct make_expr_mult_sd<true, S, T...>
{
    //TODO
    using type = expr_mult_sd<true, S, T...>;
};

//----------------------------------------------------------------------------------
//                              expr_mult_arrays
//----------------------------------------------------------------------------------
template<Integer Step, class Arr_List, class ... T>
struct expr_mult_arrays
{
    static_assert(md::dependent_false<Arr_List>::value, 
                "this type should not be instantiated");
};

template<Integer Step, class T1, class ... T, class Arr_List>
struct expr_mult_arrays<Step, Arr_List, T1, T...>
{
    using arr_1 = typename expr_mult_arrays<Step, Arr_List, T...>::type;
    using type  = typename T1::template get_arrays<Step, arr_1>;
};

template<Integer Step, class Arr_List>
struct expr_mult_arrays<Step, Arr_List>
{
    using type = Arr_List;
};

//----------------------------------------------------------------------------------
//                              can_simplify_mult
//----------------------------------------------------------------------------------
template<class S, class ... T>
struct can_simplify_mult
{
    // 0 * T1 * ... => 0 
    static const bool sc_zero   = mkd::is_scalar_data_zero<S>::value;

    static const bool value     = (sc_zero == true);
};

template<class S, class T>
struct can_simplify_mult<S, T>
{
    // 0 * T1 * ... => 0
    static const bool sc_zero   = mkd::is_scalar_data_zero<S>::value;

    // 1 * (T^1) => T1
    static const bool is_one    = mkd::is_scalar_data_one<S>::value;
    static const bool pow_one   = (T::exponent == 1);

    // T::exponent != 0 is already checked

    static const bool value     = (sc_zero == true) 
                                || (is_one == true && pow_one == true);
};

//----------------------------------------------------------------------------------
//                              simplify_expr_mult
//----------------------------------------------------------------------------------
template<class Expr>
struct simplify_expr_mult
{
    static_assert(md::dependent_false<Expr>::value, "expr_mult required");
};

template<class S, class ... T>
struct simplify_expr_mult<expr_mult_sd<false, S, T...>>
{
    //TODO: implement
    using type = expr_mult_sd<true, S, T...>;
};

template<class S, class ... T>
struct simplify_expr_mult<expr_mult_sd<true, S, T...>>
{
    using type = expr_mult_sd<true, S, T...>;
};

}}}
