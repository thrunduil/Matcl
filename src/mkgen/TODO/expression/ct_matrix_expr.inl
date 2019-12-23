#pragma once

#include <type_traits>
#include "mkgen/matrix/scalar.h"

namespace matcl { namespace mkgen
{

//----------------------------------------------------------------------------------
//                              expr_ufunc
//----------------------------------------------------------------------------------

template<class Tag, class Elem>
struct expr_ufunc : public mkd::scalar_data<expr_ufunc<Tag, Elem>>
{
    template<class Subs_Context>
    static void print(std::ostream& os, int prior)
    {
        Tag::print(os,details::prior_start);

        os << "(";
        elem::print<Subs_Context>(os,details::prior_start);
    };

    template<class Val, class Local_Storage>
    inline_expr
    static Val eval(const Local_Storage& ls)
    {
        Val v1  = elem::eval<Val>(ls);
        Val tmp = Tag::eval<Val>(ls,v1);
        return tmp;
    };
};

// TODO:
template<class Tag, class Elem>
struct make_expr_ufunc
{
    using type = expr_ufunc<Tag, Elem>;
};

//----------------------------------------------------------------------------------
//                              expr_bfunc
//----------------------------------------------------------------------------------
template<class Tag, class Elem1, class Elem2>
struct expr_bfunc
{
    template<class Subs_Context>
    static void print(std::ostream& os, int prior)
    {
        Tag::print(os,details::prior_start);

        os << "(";
        elem1::print<Subs_Context>(os,details::prior_start);
        os << ","
        elem2::print<Subs_Context>(os,details::prior_start);
        os << ")";
    };
    template<class Val, class Local_Storage>
    inline_expr
    static Val eval(const Local_Storage& ls)
    {
        Val v1  = elem1::eval<Val>(ls);
        Val v2  = elem2::eval<Val>(ls);
        Val tmp = Tag::eval<Val>(ls, v1, v2);
        return tmp;
    };
};

}}