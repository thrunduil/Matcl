#pragma once

#include "mkgen/TODO/matrix/ct_matrix.h"
#include "mkgen/TODO/expression/ct_matrix_expr.inl"
#include "mkgen/TODO/utils/utils.h"
#include "mkgen/TODO/evaler/dependency.h"

namespace matcl { namespace mkgen
{

//------------------------------------------------------------------------------
//                      for
//------------------------------------------------------------------------------
template<class It, Integer Val>
struct Item{};

template<class It, Integer Pos, template<class Subject, class Context> class Expr,
        class Subject, class Context>
struct eval_for
{
    using new_context_item  = Item<It,Pos>;
    using new_context       = typename mkgen::push_back<Context,new_context_item>::type;
    using new_subject       = typename Expr<Subject, new_context>::type;
    using type              = new_subject;

    static_assert(is_computation<Subject>::value, "subject must have computation type");
    static_assert(is_computation<new_subject>::value, "result must have computation type");
};

template<class It, Integer First, Integer Last, Integer Length, 
    template<class Subject, class Context> class Expr, class Subject, class Context>
struct for_impl
{
    static const Integer last2      = First + Length / 2 - 1;

    using new_subject   = typename for_impl<It, First, last2, Length / 2, Expr, Subject, Context>::type;
    using type          = typename for_impl<It, last2 + 1, Last, Length - Length / 2, Expr, new_subject, Context>::type;
};

template<class It, Integer First,
    template<class Subject, class Context> class Expr, class Subject, class Context>
struct for_impl<It, First, First, 1, Expr, Subject, Context>
{
    using type  = typename eval_for<It, First, Expr, Subject, Context>::type;
};
template<class It, Integer First, Integer Last,
    template<class Subject, class Context> class Expr, class Subject, class Context>
struct for_impl<It, First, Last, 2, Expr, Subject, Context>
{
    using type_1 = typename eval_for<It, First, Expr, Subject, Context>::type;
    using type_2 = typename eval_for<It, First + 1, Expr, type_1, Context>::type;
    using type = type_2;
};
template<class It, Integer First, Integer Last,
    template<class Subject, class Context> class Expr, class Subject, class Context>
struct for_impl<It, First, Last, 3, Expr, Subject, Context>
{
    using type_1    = typename eval_for<It, First, Expr, Subject, Context>::type;
    using type_2    = typename eval_for<It, First + 1, Expr, type_1, Context>::type;
    using type_3    = typename eval_for<It, First + 2, Expr, type_2, Context>::type;
    using type      = type_3;
};
template<class It, Integer First, Integer Last,
    template<class Subject, class Context> class Expr, class Subject, class Context>
struct for_impl<It, First, Last, 4, Expr, Subject, Context>
{
    using type_1    = typename eval_for<It, First, Expr, Subject, Context>::type;
    using type_2    = typename eval_for<It, First + 1, Expr, type_1, Context>::type;
    using type_3    = typename eval_for<It, First + 2, Expr, type_2, Context>::type;
    using type_4    = typename eval_for<It, First + 3, Expr, type_3, Context>::type;
    using type      = type_4;
};

template<class It, Integer First, Integer Last, template<class Subject, class Context> class Expr,
        class Subject, class Context = list<>>
struct for_expr
{
    using type      = typename for_impl<It,First,Last,Last-First+1,Expr,Subject,Context>::type;
};

//------------------------------------------------------------------------------
//                      get_from_context
//------------------------------------------------------------------------------
template<class Context, class It>
struct get_from_context
{
    static_assert(details::dependent_false<It>::value, "this type should not be instantiated");
};
template<class It1, Integer Val, class ... Items, class It>
struct get_from_context<list<Item<It1,Val>, Items...>, It>
{
    static const Integer value = Val;
};
template<class Item, class ... Items, class It>
struct get_from_context<list<Item, Items...>, It>
{
    using context_next          = list<Items...>;
    static const Integer value  = get_from_context<context_next,It>::value;
};
template<class It>
struct get_from_context<list<>, It>
{
    static_assert(details::dependent_false<It>::value, "element not found in Context");
};

//------------------------------------------------------------------------------
//                              if_expr
//------------------------------------------------------------------------------
template<bool Cond, class True_Expr, class False_Expr>
struct if_expr
{
    using type  = True_Expr ;
};
template<class True_Expr, class False_Expr>
struct if_expr<false,True_Expr,False_Expr>
{
    using type  = False_Expr;
};

}}