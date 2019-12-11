#pragma once

#include <type_traits>
#include "mkgen/matrix/scalar.h"
#include "mkgen/TODO/matrix/rational.h"

namespace matcl { namespace mkgen
{

template<class T1, class T2, bool Is_Scal_1 = is_scalar_expr<T1>::value,
                             bool Is_Scal_2 = is_scalar_expr<T2>::value>
struct make_mult;

template<Integer Step, class Arr_List, class T1, class T2>
struct expr_plus_arrays;

template<Integer Step, class Arr_List, class T1, class T2>
struct expr_minus_arrays;

template<class T1, class T2>
struct make_minus;

//----------------------------------------------------------------------------------
//                              expr_plus
//----------------------------------------------------------------------------------
template<class ... T>
struct plus_list
{};

template<class T>
struct plus_list_size
{
    static_assert(details::dependent_false<T>::value, 
                "this type should not be instantiated");
};
template<class ...T>
struct plus_list_size<plus_list<T...>>
{
    static const Integer value  = sizeof...(T);
};

template<class Plus_List, Integer N>
struct eval_plus_list
{
    static_assert(details::dependent_false<Plus_List>::value, 
                "this type should not be instantiated");
};
template<class Plus_List, Integer N>
struct eval_loop_plus_list
{
    static_assert(details::dependent_false<Plus_List>::value, 
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
struct eval_plus_list<plus_list<expr_mult<S1,T1>, T2, T...>,2>
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
struct eval_plus_list<plus_list<T1, expr_mult<S2,T2>, T...>,2>
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
struct eval_plus_list<plus_list<expr_mult<S1,T1>, expr_mult<S2,T2>, T...>,2>
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

template<class T1, class T2>
struct merge_plus_list
{};
template<class ... T1, class ... T2>
struct merge_plus_list<plus_list<T1...>,plus_list<T2...>>
{
    using type  = plus_list<T1..., T2...>;
};

template<class T>
struct make_flat_plus
{
    using type  = plus_list<T>;
};
template<class T1, class T2>
struct make_flat_plus<expr_plus<T1,T2>>
{
    using L1    = typename make_flat_plus<T1>::type;
    using L2    = typename make_flat_plus<T2>::type;
    using type  = typename merge_plus_list<L1,L2>::type;
};

template<class T1, class T2>
struct expr_plus
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
        using flat_plus         = typename make_flat_plus<expr_plus>::type;
        static const Integer N  = plus_list_size<flat_plus>::value;
        return eval_plus_list<flat_plus,N>::eval<Ret>(ls);
    };

    template<class Loop_Storage, class Ret, class Local_Storage>
    inline_expr
    static void eval_loop(Ret& ret, Integer offset, const Local_Storage& cont)
    {
        using flat_plus         = typename make_flat_plus<expr_plus>::type;
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
template<class T1, class T2>
struct expr_minus
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

//----------------------------------------------------------------------------------
//                              make_plus
//----------------------------------------------------------------------------------
template<class Scal, class ... Elems>
struct normalize_mult_scal
{
    using type = expr_mult<Scal,Elems...>;
};
template<class Elem>
struct normalize_mult_scal<one,Elem>
{
    using type = Elem;
};

template<class T1, class T2>
struct make_plus
{
    static_assert(is_mult_expr<T1>::value == false && is_mult_expr<T2>::value == false,
                  "this case should be already removed");

    using type = expr_plus<T1,T2>;
};

template<class T2> struct make_plus<zero, T2>   { using type = T2 ; };
template<class T1> struct make_plus<T1,zero>    { using type = T1; };
template<>         struct make_plus<zero,zero>  { using type = zero; };

template<Integer N1, Integer M1, class D1, Integer N2, Integer M2, class D2>
struct make_plus<ct_scalar<mkd::scalar_data<mkd::scal_data_rational<N1,M1>>,D1>,
                 ct_scalar<mkd::scalar_data<mkd::scal_data_rational<N2,M2>>,D2>>
{
    using rat_op    = rational_plus<N1, M1, N2, M2>;
    using data_type = mkd::scal_data_rational<rat_op::nominator, rat_op::denominator>;
    using type      = ct_scalar<mkd::scalar_data<data_type>,empty_deps>;
};

template<class S1, class ...T1>
struct make_plus<expr_mult<S1,T1...>,zero>
{
    using type  = expr_mult<S1,T1...>;
};
template<class ...T1>
struct make_plus<expr_mult<mone,T1...>,zero>
{
    using type  = expr_mult<mone,T1...>;
};

template<class S1, class ...T1, class T2>
struct make_plus<expr_mult<S1,T1...>,T2>
{
    using ex1   = expr_mult<S1,T1...>;
    using type  = expr_plus<ex1, T2>;
};
template<class ...T1, class T2>
struct make_plus<expr_mult<mone,T1...>,T2>
{
    using ex1   = typename normalize_mult_scal<one,T1...>::type;
    using type  = expr_minus<T2, ex1>;
};

template<class S2, class ...T2>
struct make_plus<zero, expr_mult<S2,T2...>>
{
    using type  = expr_mult<S2,T2...>;
};
template<class ...T2>
struct make_plus<zero, expr_mult<mone,T2...>>
{
    using type  = expr_mult<mone,T2...>;
};

template<class T1, class S2, class ...T2>
struct make_plus<T1, expr_mult<S2,T2...>>
{
    using ex2   = expr_mult<S2,T2...>;
    using type  = expr_plus<T1, ex2>;
};
template<class T1, class ...T2>
struct make_plus<T1, expr_mult<mone,T2...>>
{
    using ex2   = typename normalize_mult_scal<one,T2...>::type;
    using type  = expr_minus<T1, ex2>;
};

template<class S1, class ...T1, class S2, class ...T2>
struct make_plus<expr_mult<S1,T1...>, expr_mult<S2,T2...>>
{
    using ex1   = expr_mult<S1,T1...>;
    using ex2   = expr_mult<S2,T2...>;
    using type  = expr_plus<ex1,ex2>;
};
template<class ...T1, class S2, class ...T2>
struct make_plus<expr_mult<mone,T1...>, expr_mult<S2,T2...>>
{
    using ex1   = typename normalize_mult_scal<one,T1...>::type;
    using ex2   = expr_mult<S2,T2...>;
    using type  = expr_minus<ex2,ex1>;
};
template<class S1, class ...T1, class ...T2>
struct make_plus<expr_mult<S1,T1...>, expr_mult<mone,T2...>>
{
    using ex1   = expr_mult<S1,T1...>;
    using ex2   = typename normalize_mult_scal<one,T2...>::type;
    using type  = expr_minus<ex1,ex2>;
};
template<class ...T1, class ...T2>
struct make_plus<expr_mult<mone,T1...>, expr_mult<mone,T2...>>
{
    using ex1   = typename normalize_mult_scal<one,T1...>::type;
    using ex2   = typename normalize_mult_scal<one,T2...>::type;
    using ex    = expr_plus<ex1,ex2>;
    using type  = typename make_mult<mone,ex>::type;
};

//S1*T11*T12 + S2*T21*T22
template<class S1, class T11, class T12, class S2, class T21, class T22>
struct make_plus<expr_mult<S1,T11,T12>, expr_mult<S2,T21,T22>>
{
    using ex1   = expr_mult<S1,T11,T12>;
    using ex2   = expr_mult<S2,T21,T22>;
    using type  = expr_plus<ex1,ex2>;
};
template<class T11, class T12, class S2, class T21, class T22>
struct make_plus<expr_mult<mone,T11,T12>, expr_mult<S2,T21,T22>>
{
    using ex1   = expr_mult<S2,T21,T22>;
    using ex2   = expr_mult<one,T11,T12>;    
    using type  = expr_minus<ex1,ex2>;
};
template<class S1, class T11, class T12, class T21, class T22>
struct make_plus<expr_mult<S1,T11,T12>, expr_mult<mone,T21,T22>>
{
    using ex1   = expr_mult<S1,T11,T12>;
    using ex2   = expr_mult<one,T21,T22>;
    using type  = expr_minus<ex1,ex2>;
};
template<class T11, class T12, class T21, class T22>
struct make_plus<expr_mult<mone,T11,T12>, expr_mult<mone,T21,T22>>
{
    using ex1   = expr_mult<one,T11,T12>;
    using ex2   = expr_mult<one,T21,T22>;
    using ex    = expr_plus<ex1,ex2>;
    using type  = typename make_mult<mone,ex>::type;
};
//S1*T*T12 + S2*T*T22
template<class S1, class T, class T12, class S2, class T22>
struct make_plus<expr_mult<S1,T,T12>, expr_mult<S2,T,T22>>
{
    using ex1   = typename make_mult<S1,T12>::type;
    using ex2   = typename make_mult<S2,T22>::type;
    using ex    = typename make_plus<ex1,ex2>::type;
    using type  = typename make_mult<T,ex>::type;
};
template<class T, class T12, class S2, class T22>
struct make_plus<expr_mult<mone,T,T12>, expr_mult<S2,T,T22>>
{
    using ex1   = typename make_mult<S2,T22>::type;
    using ex2   = T12;
    using ex    = typename make_minus<ex1,ex2>::type;
    using type  = typename make_mult<T,ex>::type;
};
template<class S1, class T, class T12, class T22>
struct make_plus<expr_mult<S1,T,T12>, expr_mult<mone,T,T22>>
{
    using ex1   = typename make_mult<S1,T12>::type;
    using ex2   = T22;
    using ex    = typename make_minus<ex1,ex2>::type;
    using type  = typename make_mult<T,ex>::type;
};
template<class T, class T12, class T22>
struct make_plus<expr_mult<mone,T,T12>, expr_mult<mone,T,T22>>
{
    using ex3   = typename make_plus<T12,T22>::type;
    using ex    = typename make_mult<T,ex3>::type;
    using type  = typename make_mult<mone,ex>::type;
};

//S1*T*T12 + S2*T21*T
template<class S1, class T, class T12, class S2, class T21>
struct make_plus<expr_mult<S1,T,T12>, expr_mult<S2,T21,T>>
{
    using ex1   = typename make_mult<S1,T12>::type;
    using ex2   = typename make_mult<S2,T21>::type;
    using ex    = typename make_plus<ex1,ex2>::type;
    using type  = typename make_mult<T,ex>::type;
};
template<class T, class T12, class S2, class T21>
struct make_plus<expr_mult<mone,T,T12>, expr_mult<S2,T21,T>>
{
    using ex1   = typename make_mult<S2,T21>::type;
    using ex2   = T12;
    using ex    = typename make_minus<ex1,ex2>::type;
    using type  = typename make_mult<T,ex>::type;
};
template<class S1, class T, class T12, class T21>
struct make_plus<expr_mult<S1,T,T12>, expr_mult<mone,T21,T>>
{
    using ex1   = typename make_mult<S1,T12>::type;
    using ex2   = T21;
    using ex    = typename make_minus<ex1,ex2>::type;
    using type  = typename make_mult<T,ex>::type;
};
template<class T, class T12, class T21>
struct make_plus<expr_mult<mone,T,T12>, expr_mult<mone,T21,T>>
{
    using ex3   = typename make_plus<T12,T21>::type;
    using ex    = typename make_mult<T,ex3>::type;
    using type  = typename make_mult<mone,ex>::type;
};
//S1*T11*T + S2*T21*T
template<class S1, class T, class T11, class S2, class T21>
struct make_plus<expr_mult<S1,T11,T>, expr_mult<S2,T21,T>>
{
    using ex1   = typename make_mult<S1,T11>::type;
    using ex2   = typename make_mult<S2,T21>::type;
    using ex    = typename make_plus<ex1,ex2>::type;
    using type  = typename make_mult<T,ex>::type;
};
template<class T, class T11, class S2, class T21>
struct make_plus<expr_mult<mone,T11,T>, expr_mult<S2,T21,T>>
{
    using ex1   = typename make_mult<S2,T21>::type;
    using ex2   = T11;
    using ex    = typename make_minus<ex1,ex2>::type;
    using type  = typename make_mult<T,ex>::type;
};
template<class S1, class T, class T11, class T21>
struct make_plus<expr_mult<S1,T11,T>, expr_mult<mone,T21,T>>
{
    using ex1   = typename make_mult<S1,T11>::type;
    using ex2   = T21;
    using ex    = typename make_minus<ex1,ex2>::type;
    using type  = typename make_mult<T,ex>::type;
};
template<class T, class T11, class T21>
struct make_plus<expr_mult<mone,T11,T>, expr_mult<mone,T21,T>>
{
    using ex3   = typename make_plus<T11,T21>::type;
    using ex    = typename make_mult<T,ex3>::type;
    using type  = typename make_mult<mone,ex>::type;
};
//S1*T11*T + S2*T*T22
template<class S1, class T, class T11, class S2, class T22>
struct make_plus<expr_mult<S1,T11,T>, expr_mult<S2,T,T22>>
{
    using ex1   = typename make_mult<S1,T11>::type;
    using ex2   = typename make_mult<S2,T22>::type;
    using ex    = typename make_plus<ex1,ex2>::type;
    using type  = typename make_mult<T,ex>::type;
};
template<class T, class T11, class S2, class T22>
struct make_plus<expr_mult<mone,T11,T>, expr_mult<S2,T,T22>>
{
    using ex1   = typename make_mult<S2,T22>::type;
    using ex2   = T11;
    using ex    = typename make_minus<ex1,ex2>::type;
    using type  = typename make_mult<T,ex>::type;
};
template<class S1, class T, class T11, class T22>
struct make_plus<expr_mult<S1,T11,T>, expr_mult<mone,T,T22>>
{
    using ex1   = typename make_mult<S1,T11>::type;
    using ex2   = T22;
    using ex    = typename make_minus<ex1,ex2>::type;
    using type  = typename make_mult<T,ex>::type;
};
template<class T, class T11, class T22>
struct make_plus<expr_mult<mone,T11,T>, expr_mult<mone,T,T22>>
{
    using ex3   = typename make_plus<T11,T22>::type;
    using ex    = typename make_mult<T,ex3>::type;
    using type  = typename make_mult<mone,ex>::type;
};

template<class S1, class S2, class T11, class T21, bool Neg>
struct make_plus_2_impl
{
    using ex1   = typename normalize_mult_scal<S1,T11>::type;
    using ex2   = typename normalize_mult_scal<S2,T21>::type;
    using type  = expr_plus<ex1,ex2>;
};
template<class S1, class S2, class T11, class T21>
struct make_plus_2_impl<S1,S2,T11,T21,true>
{
    using ex    = typename make_minus<T11,T21>::type;
    using type  = typename make_mult<S1,ex>::type;
};

template<class S1, class S2>
struct is_div_mone
{
    static_assert(details::dependent_false<S1>::value, 
                "this type should not be instantiated");
};

template<Integer M1, Integer D1, Integer M2, Integer D2, class Deps1, class Deps2>
struct is_div_mone<ct_scalar<mkd::scalar_data<mkd::scal_data_rational<M1,D1>>, Deps1>,
                   ct_scalar<mkd::scalar_data<mkd::scal_data_rational<M2,D2>>, Deps2>>
{
    static const bool value = (D1 == D2) && (M1 == -M2);
};

//S1*T11 + S2*T21
template<class S1, class T11, class S2, class T21>
struct make_plus<expr_mult<S1,T11>, expr_mult<S2,T21>>
    :make_plus_2_impl<S1,S2,T11,T21,is_div_mone<S1,S2>::value>
{};
template<class T11, class S2, class T21>
struct make_plus<expr_mult<mone,T11>, expr_mult<S2,T21>>
{
    static_assert(std::is_same<S2,one>::value == false,"invalid_rep");
    using ex1   = expr_mult<S2,T21>;
    using ex2   = T11;    
    using type  = expr_minus<ex1,ex2>;
};
template<class S1, class T11, class T21>
struct make_plus<expr_mult<S1,T11>, expr_mult<mone,T21>>
{
    static_assert(std::is_same<S1,one>::value == false,"invalid_rep");
    using ex1   = expr_mult<S1,T11>;
    using ex2   = T21;
    using type  = expr_minus<ex1,ex2>;
};
template<class T11, class T21>
struct make_plus<expr_mult<mone,T11>, expr_mult<mone,T21>>
{
    using ex1   = T11;
    using ex2   = T21;
    using ex    = expr_plus<ex1,ex2>;
    using type  = typename make_mult<mone,ex>::type;
};

//S1*T + S2*T
template<class S1, class T, class S2>
struct make_plus<expr_mult<S1,T>, expr_mult<S2,T>>
{
    using ex    = typename make_plus<S1,S2>::type;
    using type  = typename make_mult<ex,T>::type;
};

//S*T1 + S*T2
template<class S, class T1, class T2>
struct make_plus<expr_mult<S,T1>, expr_mult<S,T2>>
{
    using ex    = typename make_plus<T1,T2>::type;
    using type  = typename make_mult<S,ex>::type;
};

//----------------------------------------------------------------------------------
//                              make_minus
//----------------------------------------------------------------------------------

template<class T1, class T2>
struct make_minus
{
    using MT2   = typename make_mult<mone,T2>::type;
    using type  = typename make_plus<T1,MT2>::type;
};

template<class T1>
struct make_uminus
{
    using type = typename make_mult<mone,T1>::type;
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

}}
