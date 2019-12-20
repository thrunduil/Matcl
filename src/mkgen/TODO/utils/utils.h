#pragma once

#include "mkgen/mkgen_fwd.h"
#include "mkgen/matrix/dependency.h"
#include "mkgen/details/utils/mpl.h"

namespace matcl { namespace mkgen
{

template<class T> struct is_computation                 {static const bool value = false; };
template<class Tag, class T, class A> 
                  struct is_computation<computation<Tag, T,A>> 
                                                        {static const bool value = true; };

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