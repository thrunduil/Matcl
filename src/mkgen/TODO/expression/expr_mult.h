#pragma once

#include <type_traits>
#include "mkgen/matrix/scalar.h"
#include "mkgen/TODO/expression/expr_plus.h"
#include "mkgen/TODO/evaler/dot_evaler.h"

namespace matcl { namespace mkgen
{

template<Integer Step, class Arr_List, class ...Mult_Items>
struct expr_mult_arrays;

//----------------------------------------------------------------------------------
//                              expr_mult
//----------------------------------------------------------------------------------
template<class ... T>
struct expr_mult
{
    static_assert(details::dependent_false_var<T...>::value, 
                "this type should not be instantiated");
};

template<class S, class T2, class T3, class ... T>
struct expr_mult<S,T2,T3,T...> : public mkd::scalar_data<expr_mult<S,T2,T3,T...>>
{
    template<class Subs_Context>
    static void print(std::ostream& os, int prior)
    {
        if (prior > details::prior_mult)
            os << "(";

        if constexpr(std::is_same<S,one>::value == false)
        {
            if constexpr(std::is_same<S,mone>::value == true)
            {
                os << "-";
            }
            else
            {
                S::print<Subs_Context>(os, details::prior_mult);
                os << "*";
            };
        };

        using expr = expr_mult<T2,T3,T...>;
        expr::print<Subs_Context>(os, details::prior_mult);

        if (prior > details::prior_mult)
            os << ")";
    };

    template<class Ret, class Local_Storage>
    inline_expr
    static Ret eval(const Local_Storage& ls)
    {
        using mult_type = expr_mult<T2,T3,T...>;

        if constexpr(std::is_same<S,one>::value == true)
            return mult_type::eval<Ret>(ls);
        else
            return S::eval<Ret>(ls) * mult_type::eval<Ret>(ls);
    };

    template<class Loop_Storage, class Ret, class Local_Storage>
    inline_expr
    static void eval_loop(Ret& ret, Integer offset, const Local_Storage& cont)
    {
        if constexpr(std::is_same<S,one>::value == true)
        {
            using mult_type = expr_mult<T2,T3,T...>;
            return mult_type::eval_loop<Loop_Storage, Ret>(ret, offset,cont);
        }
        else
        {
            Ret val1;
            Ret val2;

            S::eval_loop<Loop_Storage,Ret>(val1, offset,cont);

            using mult_type = expr_mult<T2,T3,T...>;
            
            mult_type::eval_loop<Loop_Storage,Ret>(val2, offset,cont);

            ret = val1 * val2;
        };
    };

    template<Integer Step, class Arr_List>
    using get_arrays    = typename expr_mult_arrays<Step, Arr_List, S,T2,T3,T...> :: type;

    template<class Visitor>
    static void accept(Visitor& vis)
    {
        if constexpr(std::is_same<S,one>::value == true)
        {
            expr_mult<T2,T3,T...>::accept<Visitor>(vis);
        }
        else
        {
            S::accept<Visitor>(vis);
            expr_mult<T2,T3,T...>::accept<Visitor>(vis);
            vis.visit_mult();
        };
    };
};

template<class S, class T2>
struct expr_mult<S,T2> : public mkd::scalar_data<expr_mult<S,T2>>
{
    static_assert(std::is_same<S,one>::value == false,"invalid rep");

    template<class Subs_Context>
    static void print(std::ostream& os, int prior)
    {
        if (prior > details::prior_mult)
            os << "(";

        if constexpr(std::is_same<S,mone>::value == true)
        {
            os << "-";
        }
        else
        {
            S::print<Subs_Context>(os, details::prior_mult);
            os << "*";
        };

        T2::print<Subs_Context>(os, details::prior_mult);

        if (prior > details::prior_mult)
            os << ")";
    };
    
    template<class Ret, class Local_Storage>
    inline_expr
    static Ret eval(const Local_Storage& ls)
    {
        return S::eval<Ret>(ls) * T2::eval<Ret>(ls);
    };

    template<class Loop_Storage, class Ret, class Local_Storage>
    inline_expr
    static void eval_loop(Ret& ret, Integer offset, const Local_Storage& cont)
    {
        Ret val1;
        Ret val2;
        S::eval_loop<Loop_Storage,Ret>(val1,offset,cont);
        T2::eval_loop<Loop_Storage,Ret>(val2,offset,cont);

        ret = val1 * val2;
    };

    template<Integer Step, class Arr_List>
    using get_arrays    = typename expr_mult_arrays<Step, Arr_List, S, T2>::type;

    template<class Visitor>
    static void accept(Visitor& vis)
    {
        S::accept<Visitor>(vis);
        T2::accept<Visitor>(vis);
        vis.visit_mult();            
    };
};

//----------------------------------------------------------------------------------
//                              make_mult
//----------------------------------------------------------------------------------
template<class S1, class S2>
struct make_mult_scal
{
    static_assert(md::dependent_false<S1>::value, 
                "this type should not be instantiated");
};

template<Integer N1, Integer D1, Integer N2, Integer D2>
struct make_mult_scal<rational_scalar<N1,D1>, rational_scalar<N2,D2>>
{
    using op    = rational_mult<N1,D1,N2,D2>;
    using type  = rational_scalar<op::nominator,op::denominator>;
};

template<class T1, class T2, bool Is_Scal_1, bool Is_Scal_2>
struct make_mult
{
    static_assert(is_mult_expr<T1>::value == false && is_mult_expr<T2>::value == false
                  && is_scalar_expr<T1>::value == false && is_scalar_expr<T2>::value == false,
                  "this case should be already removed");

    using type = expr_mult<one,T1,T2>;
};

template<class T1, class T2>
struct make_mult<T1,T2,true,false>
{
    static_assert(is_mult_expr<T1>::value == false && is_mult_expr<T2>::value == false,
                  "this case should be already removed");

    static_assert(std::is_same<T1,one>::value == false,"invalid_rep");
    using type = expr_mult<T1, T2>;
};

template<class T2>
struct make_mult<zero,T2,true,false>
{
    static_assert(is_mult_expr<T2>::value == false, "this case should be already removed");

    using type = zero;
};

template<class T2>
struct make_mult<one,T2,true,false>
{
    static_assert(is_mult_expr<T2>::value == false, "this case should be already removed");

    using type = T2;
};

template<class T1, class T2>
struct make_mult<T1, T2, false, true>
{
    static_assert(is_mult_expr<T1>::value == false && is_mult_expr<T2>::value == false,
                  "this case should be already removed");

    static_assert(std::is_same<T2,one>::value == false,"invalid_rep");
    using type = expr_mult<T2, T1>;
};

template<class T1>
struct make_mult<T1, zero, false, true>
{
    static_assert(is_mult_expr<T1>::value == false, "this case should be already removed");

    using type = zero;
};

template<class T1>
struct make_mult<T1, one, false, true>
{
    static_assert(is_mult_expr<T1>::value == false, "this case should be already removed");

    using type = T1;
};

template<Integer N1, Integer D1, Integer N2, Integer D2>
struct make_mult<rational_scalar<N1,D1>, rational_scalar<N2,D2>, true, true>
{
    using op    = rational_mult<N1,D1,N2,D2>;
    using type  = rational_scalar<op::nominator,op::denominator>;
};

template<class ...T>
struct check_mult_rest
{};
template<class T1, class ...T>
struct check_mult_rest<T1,T...>
{
    static const bool value = is_scalar_expr<T1>::value == false
                            && check_mult_rest<T...>::value;
};
template<>
struct check_mult_rest<>
{
    static const bool value = true;
};

template<class S, class ...T>
struct check_mult
{
    static const bool value = is_scalar_expr<S>::value
                            && check_mult_rest<T...>::value;
};

template<class S1, class ... T1>
struct make_mult<expr_mult<S1,T1...>, one,false,true>
{
    using type = expr_mult<S1,T1...>;
};
template<class S1, class ... T1>
struct make_mult<expr_mult<S1,T1...>, zero,false,true>
{
    using type = zero;
};
template<class S1, class ... T1, class T2>
struct make_mult<expr_mult<S1, T1...>, T2, false,true>
{
    using S     = typename make_mult_scal<S1,T2>::type;
    using type  = typename normalize_mult_scal<S, T1...>::type;
};
template<class S1, class ... T1, class T2>
struct make_mult<expr_mult<S1, T1...>, T2, false, false>
{
    using type = expr_mult<S1, T1...,T2>;
};

template<class S2, class ... T2>
struct make_mult<zero, expr_mult<S2, T2...>,true,false>
{
    using type = zero;
};
template<class S2, class ... T2>
struct make_mult<one, expr_mult<S2, T2...>,true,false>
{
    using type = expr_mult<S2, T2...>;
};
template<class T1, class S2, class ... T2>
struct make_mult<T1, expr_mult<S2, T2...>,true,false>
{
    using S     = typename make_mult_scal<T1,S2>::type;
    using type  = typename normalize_mult_scal<S, T2...>::type;
};
template<class T1, class S2, class ... T2>
struct make_mult<T1, expr_mult<S2, T2...>,false,false>
{
    using type = expr_mult<S2, T1, T2...>;
};

template<class S1, class ... T1, class S2, class ... T2>
struct make_mult<expr_mult<S1, T1...>, expr_mult<S2, T2...>,false,false>
{
    using S     = typename make_mult_scal<S1, S2>::type;
    using type  = expr_mult<S, T1..., T2...>;
};

//----------------------------------------------------------------------------------
//                              make_div
//----------------------------------------------------------------------------------
template<class T1>
struct make_inv
{};

template<>
struct make_inv<one>
{
    using type = one;
};
template<>
struct make_inv<mone>
{
    using type = mone;
};
template<>
struct make_inv<two>
{
    using type = half;
};
template<>
struct make_inv<half>
{
    using type = two;
};

template<class T1, class T2>
struct make_div
{
    using T2_inv    = typename make_inv<T2>::type;
    using type      = typename make_mult<T1,T2_inv>::type;
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
template<Integer Step, class S, class ... T, class Arr_List>
struct expr_mult_arrays<Step, Arr_List,S,T...>
{
    using arr_1 = typename expr_mult_arrays<Step, Arr_List,T...>::type;
    using type  = typename S::template get_arrays<Step, arr_1>;
};
template<Integer Step, class Arr_List>
struct expr_mult_arrays<Step, Arr_List>
{
    using type = Arr_List;
};

//----------------------------------------------------------------------------------
//                              expr_dot
//----------------------------------------------------------------------------------
template<Integer Step, class Arr_List, class List_1, class List_2>
struct expr_dot_arrays
{
    static_assert(md::dependent_false<List_1>::value,
                  "this type should not be instantiated");
};

template<Integer Step, class Arr_List, class ... Elems>
struct expr_dot_arrays_list
{
    using type = Arr_List;
};
template<Integer Step, class Arr_List, class Elem, class ... Elems>
struct expr_dot_arrays_list<Step,Arr_List, Elem, Elems...>
{
    using arr_1     = typename Elem::template get_arrays<Step,Arr_List>;
    using type      = typename expr_dot_arrays_list<Step, arr_1, Elems...>::type;
};

template<Integer Step, class Arr_List, class ... Elems_1, class ... Elems_2>
struct expr_dot_arrays<Step, Arr_List, list::list<Elems_1...>, list::list<Elems_2...>>
{
    using arr       = typename expr_dot_arrays_list<Step, Arr_List, Elems_1...>::type;
    using type      = typename expr_dot_arrays_list<0, arr, Elems_2...>::type;
};

template<class List_1, class List_2>
struct expand_dot
{
    static_assert(md::dependent_false<List_1>::value,
                  "this type should not be instantiated");
};
template<>
struct expand_dot<list::list<>,list::list<>>
{
    using type      = zero;
};
template<class Elem_1, class ... Elems_1, class Elem_2, class ... Elems_2>
struct expand_dot<list::list<Elem_1,Elems_1...>,list::list<Elem_2, Elems_2...>>
{
    using sum_2     = typename expand_dot<list::list<Elems_1...>,list::list<Elems_2...>>::type;
    using mult      = typename make_mult<Elem_1,Elem_2>::type;
    using type      = typename make_plus<mult,sum_2>::type;
};

template<class List_1, class List_2>
struct expr_dot : public mkd::scalar_data<expr_dot<List_1, List_2>>
{
    template<Integer Step, class Arr_List>
    using get_arrays    = typename expr_dot_arrays<Step, Arr_List, List_1, List_2> :: type;

    template<class Subs_Context>
    static void print(std::ostream& os, int prior)
    {
        using expanded  = typename expand_dot<List_1, List_2>::type;
        expanded::print<Subs_Context>(os,prior);
    };

    template<class Ret, class Local_Storage>
    inline_expr
    static Ret eval(const Local_Storage& ls)
    {
        using expanded  = typename expand_dot<List_1, List_2>::type;
        return expanded::eval<Ret>(ls);
    };

    template<class Visitor>
    static void accept(Visitor& vis)
    {
        using expanded  = typename expand_dot<List_1, List_2>::type;
        return expanded::accept(vis);
    };

    template<class Loop_Storage, class Ret, class Local_Storage>
    inline_expr
    static void eval_loop(Ret& ret, Integer offset, const Local_Storage& cont)
    {
        using value_type                = typename Local_Storage::value_type;
        static const Integer simd_size  = sizeof(Ret) / sizeof(value_type);
        using dot_evaler                = typename make_dot_evaler<List_1, List_2, simd_size>::type;

        dot_evaler::eval<Loop_Storage>(ret,offset,cont);
    };
};

}}
