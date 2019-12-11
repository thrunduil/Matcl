#pragma once

#include "mkgen/mkgen_fwd.h"
#include "mkgen/matrix/dependency.h"

namespace matcl { namespace mkgen { namespace details
{

//allow use always false conditions in static_assert 
//template<class T>
//struct dependent_false          { static const bool value = false; };

template<class... T>
struct dependent_false_var      { static const bool value = false; };

template<class T, T Val>
struct dependent_value_false    { static const bool value = false; };

}}}

namespace matcl { namespace mkgen
{

template<class T> struct is_plus_expr                   {static const bool value = false; };
template<class T1, class T2> 
                  struct is_plus_expr<expr_plus<T1,T2>> {static const bool value = true; };

template<class T> struct is_mult_expr                   {static const bool value = false; };
template<class... T> 
                  struct is_mult_expr<expr_mult<T...>>  {static const bool value = true; };


template<class T> struct is_scalar_expr                 
{ 
    static const bool value = is_value_scalar<T>::value;
};
template<class T> struct is_computation                 {static const bool value = false; };
template<class Tag, class T, class A> 
                  struct is_computation<computation<Tag, T,A>> 
                                                        {static const bool value = true; };

template <bool Cond, class Expr>
struct case_t 
{
    static const bool value = Cond;
    using type = Expr;
};

template<class ... Cases>
struct static_switch
{
};

template<bool Value, class Case, class ... Cases>
struct find_case
{};

template<class Case_1, class Case_2, class ... Cases>
struct find_case<false, Case_1, Case_2, Cases ...>
{
    using type = typename find_case<Case_2::value, Case_2, Cases...>::type;
};

template<class Case>
struct find_case<false, Case>
{
    static_assert(details::dependent_false<Case>::value, "all conditions are false");
};

template<class Case, class ... Cases>
struct find_case<true, Case, Cases...>
{
    using type = Case;
};

template<class Case, class ... Cases>
struct static_switch<Case, Cases ...> : find_case<Case::value, Case, Cases...>::type
{};

template<class T>
struct lazy_type
{
    using type = T;
};

template<bool Cond, class If_Expr_Type, class Else_Type>
struct static_if
{
    using type = If_Expr_Type;
};

template<class If_Expr_Type, class Else_Type>
struct static_if<false,If_Expr_Type,Else_Type>
{
    using type = Else_Type;
};

//-------------------------------------------------------------------------
//                      list
//-------------------------------------------------------------------------
template<class ... Elems>
struct list{};

using nil = list<>;

template<class Elem, class ... Elems>
struct is_member_tuple
{};

template<class Elem> 
struct is_member_tuple<Elem>
{ 
    static const bool value = false;
};

template<class Elem, class ... E> 
struct is_member_tuple<Elem,Elem,E...> 
{ 
    static const bool value = true;
};

template<class Elem, class E1, class ... E>
struct is_member_tuple<Elem, E1, E...>
{
    static const bool value = is_member_tuple<Elem, E...>::value;
};

template<class Elem, class List>
struct is_member_list
{
    static_assert(details::dependent_false<List>::value, 
                "this type should not be instantiated");
};

template<class Elem, class... Args>
struct is_member_list<Elem, list<Args...>>
    : is_member_tuple<Elem, Args...>
{};

template<class Elem, class Container>
struct is_member
{
    static_assert(details::dependent_false<Elem>::value, 
                "this type should not be instantiated");
};

template<class Elem, template<class ...Args> class Container, class ... Elems>
struct is_member<Elem, Container<Elems...>>
    : is_member_tuple<Elem, Elems...>
{};

template<class T>
struct get_length_container
{
    static_assert(details::dependent_false<T>::value,
                  "this type should not be instantiated");
};

template<template<class ...Args> class Container, class ... Elems>
struct get_length_container<Container<Elems...>>
{
    static const Integer value = sizeof...(Elems);
};

template<class List, class Elem>
struct push_back
{
    static_assert(details::dependent_false<List>::value,
                  "this type should not be instantiated");
};

template<class... Elems, class Elem>
struct push_back<list<Elems...>,Elem>
{
    using type = list<Elems...,Elem>;
};

template<class List, class Elem>
struct push_front{};

template<class... Elems, class Elem>
struct push_front<list<Elems...>,Elem>
{
    using type = list<Elem, Elems...>;
};

template<class List>
struct list_size{};

template<class ... Args>
struct list_size<list<Args...>>
{
    static const Integer value = sizeof... (Args);
};

//push elemented at the end of the list if is not in the list
template<class List, class Elem>
struct insert_back 
{
    using lazy_type = 
        typename static_if
        <
            is_member<List,Elem>::value,
            lazy_type<List>,
            push_back<List,Elem>
        >:: type;

    using type = typename lazy_type::type;
};

template<class List, Integer Pos>
struct get_elem_at_pos
{
    static_assert(details::dependent_false<List>::value,
                  "this type should not be instantiated");
};

template<class Elem, class ... Elems, Integer Pos>
struct get_elem_at_pos<list<Elem,Elems...>,Pos>
{
    using type = typename get_elem_at_pos<list<Elems...>,Pos-1>::type;
};

template<class Elem, class ... Elems>
struct get_elem_at_pos<list<Elem,Elems...>,0>
{
    using type = Elem;
};

template<Integer Pos>
struct get_elem_at_pos<list<>,Pos>
{
    static_assert(details::dependent_value_false<Integer, Pos>::value,
                  "invalid pos");
};

template<class List, class Elem, Integer Pos = 0>
struct get_elem_pos
{
    static_assert(details::dependent_false<List>::value,
                  "this type should not be instantiated");
};

template<class Elem, class... T, Integer Pos>
struct get_elem_pos<list<Elem, T...>, Elem, Pos> 
{
    static const Integer value = Pos;
};

template<class Elem1, class... T, class Elem2, Integer Pos>
struct get_elem_pos<list<Elem1, T...>, Elem2, Pos> 
{
    static const Integer value = get_elem_pos<list<T...>,Elem2, Pos+1>::value;
};

template<class Elem, Integer Pos>
struct get_elem_pos<list<>, Elem, Pos> 
{
    static_assert(details::dependent_false<Elem>::value,
                  "element not found in list");
};

template<class Elem, class List, bool Is_Member>
struct insert_new_elem
{
    using type = List;
};

template<class Elem, class ... Elems>
struct insert_new_elem<Elem,list<Elems...>,false>
{
    using type = list<Elem,Elems...>;
};

template<class List>
struct get_unique_list
{};

template<class Elem, class ...Elems>
struct get_unique_list<list<Elem, Elems...>>
{    
    using type_1                = typename get_unique_list<list<Elems...>>::type;
    static const bool is_member = is_member_list<Elem, type_1>::value;
    using type                  = typename insert_new_elem<Elem, type_1,is_member>::type;    
};

template<>
struct get_unique_list<list<>>
{
    using type                  = list<>;
};

template<class Elem>
struct get_unique_list<list<Elem>>
{
    using type                  = list<Elem>;
};

//----------------------------------------------------------------------------------
//                              stack_array
//----------------------------------------------------------------------------------
template<class Val, Integer Size>
struct stack_array
{
    alignas(VEC_ALIGN) 
    Val             data[Size];

    template<Integer Pos>
    inline_lev_1
    Val&            elem()          { return data[Pos]; };

    template<Integer Pos>
    inline_lev_1
    const Val&      elem() const    { return data[Pos]; };

    inline_lev_1
    Val*            get_array()     { return data; };
};

template<class Val>
struct stack_array<Val,0>
{
    inline_lev_1
    Val*            get_array()     { return nullptr; };
};

}}