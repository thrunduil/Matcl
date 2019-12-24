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

namespace matcl { namespace mkgen { namespace details
{

//----------------------------------------------------------------------------------
//                              forward declarations
//----------------------------------------------------------------------------------

template<Integer Step, class Arr_List, class ... Items>
struct expr_plus_arrays;

template<class S, class ... T>
struct can_simplify_plus;

template<Scal_data S, class T>
struct plus_item;

template<class M>
struct plus_item_mult;

template<class Expr>
struct simplify_expr_plus;

//----------------------------------------------------------------------------------
//                              checks
//----------------------------------------------------------------------------------
template<class T>
struct check_plus_item
{
    static_assert(md::dependent_false<T>::value, "plus_item required");
};

template<Scal_data S, class T>
struct check_plus_item<plus_item<S, T>>
{
    using type = void;
};

template<class M>
struct check_plus_item<plus_item_mult<M>>
{
    using type = void;
};

//----------------------------------------------------------------------------------
//                              plus_item
//----------------------------------------------------------------------------------
template<Scal_data S, class T>
struct plus_item
{
    static const bool is_vs1 = mkd::is_value_scalar_data<S>::value;
    static_assert(is_vs1 == true, "value scalar required");

    static const bool is_vs2 = mkd::is_value_scalar_data<T>::value;
    static_assert(is_vs2 == false, "value scalar unexpected");

    static const bool is_0  = mkd::is_scalar_data_zero<S>::value;
    static_assert(is_0 == false, "invalid scalling");

    static const bool is_add    = mkd::is_plus_expr<T>::value;
    static_assert(is_add == false, "plus expression not allowed");

    static const bool is_mult    = mkd::is_mult_expr<T>::value;
    static_assert(is_mult == false, "mult expression not allowed");

    static const bool is_simpl  = T::is_simplified();
    static_assert(is_simpl == true, "subexpression not simplified");

    using scalling      = S;
    using base          = T;

    template<Integer Step, class Arr_List>
    using get_arrays    = typename T::template get_arrays<Step, Arr_List>;

    template<class Subs_Context>
    static void print(std::ostream& os, int prior)
    {
        if constexpr(mkd::is_scalar_data_one<S>::value == true)
        {
            base::print<Subs_Context>(os, prior);
        }
        else
        {
            if (prior > details::prior_mult)
                os << "(";

            scalling::print<Subs_Context>(os, prior);

            os << "*";

            base::print<Subs_Context>(os, details::prior_mult);

            if (prior > details::prior_mult)
                os << ")";
        };
    };

    template<class Ret, class Local_Storage>
    inline_expr
    static Ret eval(const Local_Storage& ls)
    {
        //TODO: optimize

        if constexpr(mkd::is_scalar_data_one<S>::value == true)
            return base::eval<Ret>(ls);
        else
            return scalling::eval<Ret>(ls) * base::eval<Ret>(ls);
    };

    template<class Loop_Storage, class Ret, class Local_Storage>
    inline_expr
    static void eval_loop(Ret& ret, Integer offset, const Local_Storage& cont)
    {
        //TODO: optimize

        if constexpr(mkd::is_scalar_data_one<S>::value == true)
        {
            base::eval_loop<Loop_Storage, Ret>(ret, offset, cont);
        }
        else
        {
            Ret val1;
            base::eval_loop<Loop_Storage, Ret>(ret, offset, cont);

            ret = scalling::eval<Ret>(ls) * val1;
        };
    };

    template<class Visitor>
    static void accept(Visitor& vis)
    {
        scalling::accept<Visitor>(vis);
        base::accept<Visitor>(vis);
        
        vis.visit_mult();
    };
};

//----------------------------------------------------------------------------------
//                              plus_item_mult
//----------------------------------------------------------------------------------
// M can take the form S * T, where T is plus expression, when constant S cannot be
// accumulated in T
template<class M>
struct plus_item_mult
{
    static const bool is_mult    = mkd::is_mult_expr<M>::value;
    static_assert(is_mult == true, "mult expression required");

    static const bool is_simpl  = M::is_simplified();
    static_assert(is_simpl == true, "subexpression not simplified");

    using scalling      = typename M :: scalling;
    using base          = typename M :: base;

    template<Integer Step, class Arr_List>
    using get_arrays    = typename M::template get_arrays<Step, Arr_List>;

    template<class Subs_Context>
    static void print(std::ostream& os, int prior)
    {
        M::print<Subs_Context>(os, prior);
    };

    template<class Ret, class Local_Storage>
    inline_expr
    static Ret eval(const Local_Storage& ls)
    {
        //TODO: optimize
        return M::eval<Ret>(ls);
    };

    template<class Loop_Storage, class Ret, class Local_Storage>
    inline_expr
    static void eval_loop(Ret& ret, Integer offset, const Local_Storage& cont)
    {
        //TODO: optimize
        M::eval_loop<Loop_Storage, Ret>(ret, offset, cont);
    };

    template<class Visitor>
    static void accept(Visitor& vis)
    {
        M::accept<Visitor>(vis);
    };
};

//----------------------------------------------------------------------------------
//                              make_plus_item
//----------------------------------------------------------------------------------
template<class T>
struct make_plus_item
{
    using type = plus_item<one_sd, T>;
};

template<bool F, class S, class ... T>
struct make_plus_item<expr_mult_sd<F, S, T...>>
{
    using mult  = expr_mult_sd<F, S, T...>;
    using type  = plus_item_mult<mult>;
};

//----------------------------------------------------------------------------------
//                              mult_plus_item
//----------------------------------------------------------------------------------
// it is assumed, that PI is simplified; returned type is also simplified
template<class S, class PI>
struct mult_plus_item
{
    static_assert(md::dependent_false<S>::value, "plus_item required");
};

template<class S, Scal_data SP, class TP>
struct mult_plus_item<S, plus_item<SP, TP>>
{
    static_assert(mkd::is_plus_expr<TP>::value == false, "unexpected plus expression");

    using S_ret = typename make_mult_scal<S, SP>::type;
    using type  = plus_item<S_ret, TP>;

    static const bool is_plus   = false;
};

template<class S, bool FM, class SM, class ... TM>
struct mult_plus_item<S, plus_item_mult<expr_mult_sd<FM, SM, TM...>>>
{
    using mult  = expr_mult_sd<FM, SM, TM...>;
    using type0 = typename make_mult_root<mult, S>::type;

    // multiplication by scalar returns mult expression with different scalling,
    // or when scalling is 1, and number of subexpressions is 1, then subexpression
    // can be returned. However since all subexpressions of mult expression are 
    // simplified, then type0 is also simplified

    static const bool is_plus   = mkd::is_plus_expr<type0>::value;

    using type   = typename mkd::static_if
                        <   is_plus == true,
                            mkd::lazy_type<type0>,
                            make_plus_item<type0>
                        > ::type :: type;
};

//----------------------------------------------------------------------------------
//                              mult_plus_items
//----------------------------------------------------------------------------------
// form S * PI, where PI is a plus item; if result is plus expression, then add this
// result to Plus_items list; otherwise create plus_item and push it to Mult_items
// return list<Mult_items, Plus_items>
template<class Mult_items, class Plus_items, class S, class ... PI>
struct mult_plus_items
{
    static_assert(md::dependent_false<S>::value, "this type should not be instantiated");
};

template<class Mult_items, class Plus_items, class S, class PI1, class ... PI>
struct mult_plus_items<Mult_items, Plus_items, S, PI1, PI ...>
{
    using gen   = mult_plus_item<S, PI1>;
    using It    = typename gen :: type;

    static const bool is_plus   = gen::is_plus;

    // push It to Mult_items if is_plus == false
    using Mult_items2   = typename mkd::static_if
                        <   is_plus == true,
                            mkd::lazy_type<Mult_items>,
                            list::push_back<Mult_items, It>
                        > ::type :: type;

    // push It to Plus_items if is_plus == false
    using Plus_items2   = typename mkd::static_if
                        <   is_plus == false,
                            mkd::lazy_type<Plus_items>,
                            list::push_back<Plus_items, It>
                        > ::type :: type;

    // process other items
    using type          = typename mult_plus_items<Mult_items2, Plus_items2, S, PI ...>
                                        :: type;
};

template<class Mult_items, class Plus_items, class S>
struct mult_plus_items<Mult_items, Plus_items, S>
{
    using type = list::list<Mult_items, Plus_items>;
};

//----------------------------------------------------------------------------------
//                              can_accumulate_scalling
//----------------------------------------------------------------------------------
// PI must be plus_item
template<class PI>
struct plus_item_has_trivial_scalling
{
    using scal  = typename PI :: scalling;

    static const bool value = mkd::is_scalar_data_one<scal>::value
                            || mkd::is_scalar_data_mone<scal>::value;
};

// return true if all plus_items has nontrivial scalling (i.e. different from 1 or -1)
template<class Plus_expr>
struct can_accumulate_scalling
{
    static_assert(md::dependent_false<S>::value, "expr_plus_sd required");
};

template<bool Flag, class S, class ... T>
struct can_accumulate_scalling<expr_plus_sd<Flag, S, T...>>
{
    static const bool value = (plus_item_has_trivial_scalling<T>::value || ...) == false;
};

//----------------------------------------------------------------------------------
//                              plus_item_to_mult
//----------------------------------------------------------------------------------
// PI must be simplified
template<class PI>
struct plus_item_to_mult
{
    static_assert(md::dependent_false<PI>::value, "plus_item required");
};

template<class M>
struct plus_item_to_mult<plus_item_mult<M>>
{
    using type = M;
};

template<Scal_data S, class T>
struct plus_item_to_mult<plus_item<S, T>>
{
    static const bool is_one    = mkd::is_scalar_data_one<S>::value;

    using type0  = typename mkd::static_if
                        <   is_one == true,
                            mkd::lazy_type<T>,
                            make_expr_mult_sd<true, S, mult_item<T, 1>>
                        >::type;
    using type  = typename type0 :: type;
};

//----------------------------------------------------------------------------------
//                              expr_plus_tail
//----------------------------------------------------------------------------------
template<class T1, class ... T>
struct expr_plus_tail
{
    static const Integer size   = 1 + sizeof...(T);

    template<class Subs_Context>
    static void print(std::ostream& os, int prior)
    {
        if (prior > details::prior_plus)
            os << "(";

        T1::print<Subs_Context>(os, details::prior_plus);                

        if constexpr(size > 1)
        {
            using tail  = expr_plus_tail<T...>;

            os << "+";
            tail::print<Subs_Context>(os, details::prior_plus);
        };

        if (prior > details::prior_plus)
            os << ")";
    };

    template<class Ret, class Local_Storage>
    inline_expr
    static Ret eval(const Local_Storage& ls)
    {
        //TODO: optimize this
        if constexpr(size > 1)
        {
            using tail  = expr_plus_tail<T...>;
            return T1::eval<Ret>(ls) + tail::eval<Ret>(ls);
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
            using tail  = expr_plus_tail<T...>;

            Ret val1;
            Ret val2;

            T1::eval_loop<Loop_Storage,Ret>(val1, offset,cont);            
            tail::eval_loop<Loop_Storage,Ret>(val2, offset,cont);

            ret = val1 + val2;
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
            using tail  = expr_plus_tail<T...>;        
            tail::accept<Visitor>(vis);
        };

        vis.visit_plus();
    };
};

//----------------------------------------------------------------------------------
//                              expr_plus_sd
//----------------------------------------------------------------------------------

template<bool Flag, class S, class ... T>
struct expr_plus_sd : public mkd::scalar_data<expr_plus_sd<Flag, S, T ...>>
{
    static const bool is_vs = mkd::is_value_scalar_data<S>::value;
    static_assert(is_vs == true, "value scalar reqired");

    static_assert(sizeof...(T) > 0, "invalid rep size");

    // checks    
    using ch1               = list::list<typename check_plus_item<T>::type ...>;

    static const bool simp  = can_simplify_plus<S, T...>::value;
    static_assert(simp == false, "invalid rep simpl");

    using this_type         = expr_plus_sd<Flag, S, T ...>;

    template<class Subs_Context>
    static void print(std::ostream& os, int prior)
    {
        if (prior > details::prior_plus)
            os << "(";

        if constexpr(mkd::is_scalar_data_zero<S>::value == false)
        {
            S::print<Subs_Context>(os, details::prior_plus);
            os << "+";
        };

        using tail  = expr_plus_tail<T...>;
        tail::print<Subs_Context>(os, details::prior_plus);

        if (prior > details::prior_plus)
            os << ")";
    };

    //TODO: optimizations: 
    //  1. remove scalar mult, for example (-x1*x2) + (y1*y2)   => y1*y2 - x1*x2
    //                                     (-x1*x2) + (-y1*y2)  => -(y1*y2 + x1*x2)
    //                                     (-a*x1)  + (a*y1)    => a * (y1 - x1)
    //  2. collect common terms,           (x*x1*x2) + (y1*y2*x)=> (x1*x2 + y1*y2) * x
    //  3. use fma or minus when possible
    //  4. split into series of independent evaluations, i.e. x1 + x2 + x3 + x4 + ...
    //      must be evaluated as (x1 + x2) + (x3 + x4) + ...
    template<class Ret, class Local_Storage>
    inline_expr
    static Ret eval(const Local_Storage& ls)
    {
        //TODO: optimize this
        using tail  = expr_plus_tail<T...>;

        if constexpr(mkd::is_scalar_data_zero<S>::value == true)
            return tail::eval<Ret>(ls);
        else
            return S::eval<Ret>(ls) + tail::eval<Ret>(ls);
    }

    template<class Loop_Storage, class Ret, class Local_Storage>
    inline_expr
    static void eval_loop(Ret& ret, Integer offset, const Local_Storage& cont)
    {
        //TODO: optimize this
        using tail  = expr_plus_tail<T...>;

        if constexpr(mkd::is_scalar_data_zero<S>::value == true)
        {            
            return tail::eval_loop<Loop_Storage, Ret>(ret, offset, cont);
        }
        else
        {
            Ret val1;
            Ret val2;

            S::eval_loop<Loop_Storage,Ret>(val1, offset,cont);            
            tail::eval_loop<Loop_Storage,Ret>(val2, offset,cont);

            ret = val1 + val2;
        };
    };

    template<Integer Step, class Arr_List>
    using get_arrays    = typename expr_plus_arrays<Step, Arr_List, T...> :: type;

    template<class Void>
    using simplify      = typename simplify_expr_plus<this_type>::type;

    static constexpr bool is_simplified()   { return Flag; };

    template<class Visitor>
    static void accept(Visitor& vis)
    {
        using tail  = expr_plus_tail<T...>;

        if constexpr(mkd::is_scalar_data_zero<S>::value == true)
        {
            tail::accept<Visitor>(vis);
        }
        else
        {
            S::accept<Visitor>(vis);
            tail::accept<Visitor>(vis);
            vis.visit_plus();
        };
    };
};

//----------------------------------------------------------------------------------
//                              make_expr_plus_sd
//----------------------------------------------------------------------------------
// it is assumed that T::simplify<void> was called
// // if Is_simpl == true, then also expr_plus was simplified previously
template<bool Is_simpl, class S, class ... T>
struct make_expr_plus_sd
{
    static_assert(sizeof...(T) > 0, "invalid list of items");
    using type = expr_plus_sd<false, S, T...>;
};

template<class S, class ... T>
struct make_expr_plus_sd<true, S, T...>
{
    static_assert(sizeof...(T) > 0, "invalid list of items");

    //TODO:
    using type = expr_plus_sd<true, S, T...>;
};

// it is assumed that T::simplify<void> was called for every element in the list List
// List contains valid plus_items
// If Is_simpl == true, then expr_plus<S, T...> was simplified previously
template<bool Is_simpl, class S, class List>
struct make_expr_plus_list
{
    static_assert(md::dependent_false<S>::value, "list::list<...> required");
};

template<bool Is_simpl, class S, class ... It>
struct make_expr_plus_list<Is_simpl, S, list::list<It...>>
{
    static_assert(sizeof...(It) > 1, "this case should already be processed");
    using type   = typename make_expr_plus_sd<Is_simpl, S, It ...> :: type;
};

template<bool Is_simpl, class S, class It>
struct make_expr_plus_list<Is_simpl, S, list::list<It>>
{
    // expression 0 + It can be simplified further

    static const bool is_zero   = mkd::is_scalar_data_zero<S>::value;

    using type   = typename mkd::static_if
                        <   is_zero == true,
                            plus_item_to_mult<It>,
                            make_expr_plus_sd<true, S, It>
                        > ::type :: type;
};

template<bool Is_simpl, class S>
struct make_expr_plus_list<Is_simpl, S, list::list<>>
{
    using type  = S;
};

//----------------------------------------------------------------------------------
//                              expr_plus_arrays
//----------------------------------------------------------------------------------
template<Integer Step, class Arr_List, class ... T>
struct expr_plus_arrays
{
    static_assert(md::dependent_false<Arr_List>::value, 
                "this type should not be instantiated");
};

template<Integer Step, class T1, class ... T, class Arr_List>
struct expr_plus_arrays<Step, Arr_List, T1, T...>
{
    using arr_1 = typename expr_plus_arrays<Step, Arr_List, T...>::type;
    using type  = typename T1::template get_arrays<Step, arr_1>;
};

template<Integer Step, class Arr_List>
struct expr_plus_arrays<Step, Arr_List>
{
    using type = Arr_List;
};

//----------------------------------------------------------------------------------
//                              can_simplify_plus
//----------------------------------------------------------------------------------
template<class S, class ... T>
struct can_simplify_plus
{
    static const bool value     = false;
};

template<class S, class T>
struct can_simplify_plus<S, T>
{
    // 0  + T1  => T1
    static const bool sc_zero   = mkd::is_scalar_data_zero<S>::value;

    static const bool value     = (sc_zero == true);
};

//----------------------------------------------------------------------------------
//                              simplify_expr_plus
//----------------------------------------------------------------------------------
template<class Expr>
struct simplify_expr_plus
{
    static_assert(md::dependent_false<Expr>::value, "expr_plus required");
};

template<class S, class ... T>
struct simplify_expr_plus<expr_plus_sd<false, S, T...>>
{
    //TODO: implement
    using type = expr_plus_sd<true, S, T...>;
};

template<class S, class ... T>
struct simplify_expr_plus<expr_plus_sd<true, S, T...>>
{
    using type = expr_plus_sd<true, S, T...>;
};

}}}
