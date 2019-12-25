#pragma once

#include <type_traits>
#include "mkgen/matrix/scalar.h"

namespace matcl { namespace mkgen
{

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