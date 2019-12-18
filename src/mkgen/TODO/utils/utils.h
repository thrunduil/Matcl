#pragma once

#include "mkgen/mkgen_fwd.h"
#include "mkgen/matrix/dependency.h"
#include "mkgen/details/utils/mpl.h"

namespace matcl { namespace mkgen
{


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
    static_assert(md::dependent_false<Case>::value, "all conditions are false");
};

template<class Case, class ... Cases>
struct find_case<true, Case, Cases...>
{
    using type = Case;
};

template<class Case, class ... Cases>
struct static_switch<Case, Cases ...> : find_case<Case::value, Case, Cases...>::type
{};


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