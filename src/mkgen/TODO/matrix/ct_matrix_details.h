#pragma once

#include <iosfwd>

#include "mkgen/mkgen_fwd.h"
#include "mkgen/TODO/utils/utils.h"
#include "mkgen/TODO/expression/colon.h"

namespace matcl { namespace mkgen
{

template<class Array_t, Integer Offset, Integer Step>
struct sub_array_1 {};

template<class Array_t, Integer Offset1, Integer Offset2, Integer Step1, Integer step2>
struct sub_array_2 {};

}};

namespace matcl { namespace mkgen { namespace details
{

//--------------------------------------------------------
//              print priorities
//--------------------------------------------------------
static const int prior_start    = 0;
static const int prior_assign   = 10;
static const int prior_plus     = 20;
static const int prior_mult     = 30;

//------------------------------------------------------------------------------
//                      submatrix_maker
//------------------------------------------------------------------------------
template<class Mat, class Colon_1>
struct submatrix_maker_1
{
    static_assert(dependent_false<Mat>::value, "class T must be ct_matrix");
};
template<Integer M, Integer N, class Array_t, class Deps, class Colon_1>
struct submatrix_maker_1<ct_matrix<M,N,Array_t,Deps>,Colon_1>
{
    static const Integer size   = get_size_colon<Colon_1, M>::value;
    static const Integer offset = get_offset_colon<Colon_1>::value;
    static const Integer step   = get_step_colon<Colon_1>::value;

    using new_array = sub_array_1<Array_t,offset,step>;
    using type      = ct_matrix<size,1,new_array,Deps>;
};

template<class Mat, class Colon_1, class Colon_2>
struct submatrix_maker_2
{
    static_assert(dependent_false<Mat>::value, "class T must be ct_matrix");
};

template<Integer M, Integer N, class Array_t, class Deps, class Colon_1, class Colon_2>
struct submatrix_maker_2<ct_matrix<M,N,Array_t,Deps>,Colon_1,Colon_2>
{
    static const Integer size1      = get_size_colon<Colon_1, M>::value;
    static const Integer size2      = get_size_colon<Colon_2, N>::value;
    static const Integer offset1    = get_offset_colon<Colon_1>::value;
    static const Integer offset2    = get_offset_colon<Colon_2>::value;
    static const Integer step1      = get_step_colon<Colon_1>::value;
    static const Integer step2      = get_step_colon<Colon_2>::value;

    using new_array = sub_array_2<Array_t,offset1, offset2,step1,step2>;
    using type      = ct_matrix<size1,size2,new_array,Deps>;
};

}}}
