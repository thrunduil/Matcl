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

//TODO
namespace matcl { namespace mkgen { namespace details
{

//----------------------------------------------------------------------------------
//                              forward declarations
//----------------------------------------------------------------------------------

template<Integer Step, class Arr_List, class T1, class T2>
struct expr_plus_arrays;

template<Integer Step, class Arr_List, class T1, class T2>
struct expr_minus_arrays;

template<class T>
struct make_flat_plus;

//TODO: remove?
template<class ... T>
struct plus_list;

template<class T1, class T2>
struct merge_plus_list;

template<class T>
struct plus_list_size;

template<class Plus_List, Integer N>
struct eval_loop_plus_list;

template<class Plus_List, Integer N>
struct eval_plus_list;

//----------------------------------------------------------------------------------
//                              expr_plus_scalar_data
//----------------------------------------------------------------------------------

//TODO: add checks
template<class T1, class T2>
struct expr_plus_scalar_data : public mkd::scalar_data<expr_plus_scalar_data<T1, T2>>
{
    template<class Subs_Context>
    static void print(std::ostream& os, int prior)
    {
        if (prior > details::prior_plus)
            os << "(";

        T1::print<Subs_Context>(os, details::prior_plus);

        os << "+";

        T2::print<Subs_Context>(os, details::prior_plus);

        if (prior > details::prior_plus)
            os << ")";
    };

    template<class Ret, class Local_Storage>
    inline_expr
    static Ret eval(const Local_Storage& ls)
    {
        using flat_plus         = typename make_flat_plus<expr_plus_scalar_data>::type;
        static const Integer N  = plus_list_size<flat_plus>::value;
        return eval_plus_list<flat_plus,N>::eval<Ret>(ls);
    };

    template<class Loop_Storage, class Ret, class Local_Storage>
    inline_expr
    static void eval_loop(Ret& ret, Integer offset, const Local_Storage& cont)
    {
        using flat_plus         = typename make_flat_plus<expr_plus_scalar_data>::type;
        static const Integer N  = plus_list_size<flat_plus>::value;
        return eval_loop_plus_list<flat_plus,N>::eval<Loop_Storage,Ret>(ret,offset,cont);
    };

    template<Integer Step, class Arr_List>
    using get_arrays    = typename expr_plus_arrays<Step, Arr_List, T1, T2> :: type;

    template<class Visitor>
    static void accept(Visitor& vis)
    {
        T1::accept<Visitor>(vis);
        T2::accept<Visitor>(vis);
        vis.visit_plus();
    };
};

//----------------------------------------------------------------------------------
//                              expr_plus_arrays
//----------------------------------------------------------------------------------
template<Integer Step, class Arr_List, class T1, class T2>
struct expr_plus_arrays
{
    using arr_1 = typename T1::template get_arrays<Step, Arr_List>;
    using type  = typename T2::template get_arrays<Step, arr_1>;
};

template<Integer Step, class Arr_List, class T1, class T2>
struct expr_minus_arrays
{
    using arr_1 = typename T1::template get_arrays<Step, Arr_List>;
    using type  = typename T2::template get_arrays<Step, arr_1>;
};

//----------------------------------------------------------------------------------
//                              make_flat_plus
//----------------------------------------------------------------------------------
template<class T>
struct make_flat_plus
{
    using type  = plus_list<T>;
};

template<class T1, class T2>
struct make_flat_plus<expr_plus_scalar_data<T1,T2>>
{
    using L1    = typename make_flat_plus<T1>::type;
    using L2    = typename make_flat_plus<T2>::type;
    using type  = typename merge_plus_list<L1,L2>::type;
};

//----------------------------------------------------------------------------------
//                              plus_list
//----------------------------------------------------------------------------------
template<class ... T>
struct plus_list
{};

template<class T1, class T2>
struct merge_plus_list
{};

template<class ... T1, class ... T2>
struct merge_plus_list<plus_list<T1...>,plus_list<T2...>>
{
    using type  = plus_list<T1..., T2...>;
};

template<class T>
struct plus_list_size
{
    static_assert(md::dependent_false<T>::value, 
                "this type should not be instantiated");
};
template<class ...T>
struct plus_list_size<plus_list<T...>>
{
    static const Integer value  = sizeof...(T);
};

//----------------------------------------------------------------------------------
//                              eval_loop_plus_list
//----------------------------------------------------------------------------------

template<class Plus_List, Integer N>
struct eval_loop_plus_list
{
    static_assert(md::dependent_false<Plus_List>::value, 
                "this type should not be instantiated");
};

template<class ... T, Integer N>
struct eval_loop_plus_list<plus_list<T...>,N>
{
    using first_half        = eval_loop_plus_list<plus_list<T...>, N / 2>;
    using rem_args_1        = typename first_half::remaining_args;
    
    using second_half       = eval_loop_plus_list<rem_args_1, N - N / 2>;
    using rem_args_2        = typename second_half::remaining_args;
    
    using remaining_args    = rem_args_2;

    template<class Loop_Storage, class Ret, class Local_Storage>
    inline_expr_split
    static void eval(Ret& ret, Integer offset, const Local_Storage& cont)
    {
        Ret val1;
        Ret val2;

        first_half::eval<Loop_Storage, Ret>(val1, offset,cont);
        second_half::eval<Loop_Storage, Ret>(val2, offset,cont);

        ret = val1 + val2;
    };
};

template<class ... T>
struct eval_loop_plus_list<plus_list<T...>,0>
{
    using remaining_args    = plus_list<T...>;

    template<class Loop_Storage, class Ret, class Local_Storage>
    static void eval(Ret& ret, Integer offset, const Local_Storage& cont)
    {};
};

template<class T1, class ... T>
struct eval_loop_plus_list<plus_list<T1, T...>,1>
{
    using remaining_args    = plus_list<T...>;

    template<class Loop_Storage, class Ret, class Local_Storage>
    inline_expr
    static void eval(Ret& ret, Integer offset, const Local_Storage& cont)
    {
        T1::eval_loop<Loop_Storage,Ret>(ret, offset, cont);
    };
};

template<class T1, class T2, class ... T>
struct eval_loop_plus_list<plus_list<T1, T2, T...>, 2>
{
    using remaining_args    = plus_list<T...>;

    template<class Loop_Storage, class Ret, class Local_Storage>
    inline_expr
    static void eval(Ret& ret, Integer offset, const Local_Storage& cont)
    {
        Ret v1;
        Ret v2;
        T1::eval_loop<Loop_Storage,Ret>(v1, offset, cont);
        T2::eval_loop<Loop_Storage,Ret>(v2, offset, cont);
        ret = v1 + v2;
    };
};

//----------------------------------------------------------------------------------
//                              eval_plus_list
//----------------------------------------------------------------------------------
template<class Plus_List, Integer N>
struct eval_plus_list
{
    static_assert(md::dependent_false<Plus_List>::value, 
                "this type should not be instantiated");
};

template<class ... T, Integer N>
struct eval_plus_list<plus_list<T...>,N>
{
    using first_half        = eval_plus_list<plus_list<T...>, N / 2>;
    using rem_args_1        = typename first_half::remaining_args;
    
    using second_half       = eval_plus_list<rem_args_1, N - N / 2>;
    using rem_args_2        = typename second_half::remaining_args;
    
    using remaining_args    = rem_args_2;

    template<class Ret, class Local_Storage>
    inline_expr_split
    static Ret eval(const Local_Storage& ls)
    {
        Ret tmp1 = first_half::eval<Ret>(ls);
        Ret tmp2 = second_half::eval<Ret>(ls);
        return tmp1 + tmp2;
    };
};
template<class ... T>
struct eval_plus_list<plus_list<T...>,0>
{    
    using remaining_args    = plus_list<T...>;

    template<class Ret, class Local_Storage>
    inline_lev_1
    static Ret eval(const Local_Storage& ls)
    {
        return Ret(0);
    };
};
template<class T1, class ... T>
struct eval_plus_list<plus_list<T1, T...>,1>
{    
    using remaining_args    = plus_list<T...>;

    template<class Ret, class Local_Storage>
    inline_expr
    static Ret eval(const Local_Storage& ls)
    {
        return T1::eval<Ret>(ls);
    };
};

template<class T1, class T2, class ... T>
struct eval_plus_list<plus_list<T1, T2, T...>,2>
{    
    using remaining_args    = plus_list<T...>;

    template<class Ret, class Local_Storage>
    inline_expr
    static Ret eval(const Local_Storage& ls)
    {
        Ret v1  = T1::eval<Ret>(ls);
        Ret v2  = T2::eval<Ret>(ls);
        return v1 + v2;
    };
};
template<class S1, class T1, class T2, class ... T>
struct eval_plus_list<plus_list<expr_mult_scalar_data<S1,T1>, T2, T...>,2>
{    
    using remaining_args    = plus_list<T...>;

    template<class Ret, class Local_Storage>
    inline_expr
    static Ret eval(const Local_Storage& ls)
    {
        Ret v1  = S1::eval<Ret>(ls);
        Ret v2  = T1::eval<Ret>(ls);
        Ret v3  = T2::eval<Ret>(ls);
        return fma_f(v1,v2, v3);
    };
};
template<class T1, class S2, class T2, class ... T>
struct eval_plus_list<plus_list<T1, expr_mult_scalar_data<S2,T2>, T...>,2>
{    
    using remaining_args    = plus_list<T...>;

    template<class Ret, class Local_Storage>
    inline_expr
    static Ret eval(const Local_Storage& ls)
    {
        Ret v1  = T1::eval<Ret>(ls);
        Ret v2  = S2::eval<Ret>(ls);
        Ret v3  = T2::eval<Ret>(ls);
        return fma_f(v2,v3,v1);
    };
};
template<class S1, class T1, class S2, class T2, class ... T>
struct eval_plus_list<plus_list<expr_mult_scalar_data<S1,T1>, expr_mult_scalar_data<S2,T2>, T...>,2>
{    
    using remaining_args    = plus_list<T...>;

    template<class Ret, class Local_Storage>
    inline_expr
    static Ret eval(const Local_Storage& ls)
    {
        Ret v1  = S1::eval<Ret>(ls);
        Ret v2  = T1::eval<Ret>(ls);
        Ret v3  = S2::eval<Ret>(ls);
        Ret v4  = T2::eval<Ret>(ls);
        return fma_f(v1,v2, v3*v4);
    };
};

}}}
