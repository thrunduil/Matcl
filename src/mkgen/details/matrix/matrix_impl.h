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

#include "mkgen/matrix/matrix.h"
#include "mkgen/details/matrix/matrix_arrays.h"
#include "mkgen/details/matrix/colon_func.h"

namespace matcl { namespace mkgen { namespace details
{

//------------------------------------------------------------------------------
//                      submatrix_maker_1
//------------------------------------------------------------------------------
// implements ct_matrix<>::sub(Colon_1)
template<class Mat, class Colon_1>
struct submatrix_maker_1
{
    static_assert(dependent_false<Mat>::value, "class T must be ct_matrix");
};

template<Integer M, Integer N, class Array_t, class Deps, class Colon_1>
struct submatrix_maker_1<ct_matrix<M, N, Array_t, Deps>, Colon_1>
{
    using dum                   = typename colon_func::check_colon<Colon_1, M * N>::type;

    static const Integer size   = colon_func::size<Colon_1, M>::value;
    static const Integer offset = colon_func::offset<Colon_1>::value;
    static const Integer step   = colon_func::step<Colon_1>::value;

    using new_array = sub_array_1<Array_t, offset, step>;
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
    using dum1                      = typename colon_func::check_colon<Colon_1, M>::type;
    using dum2                      = typename colon_func::check_colon<Colon_2, N>::type;

    static const Integer size1      = colon_func::size<Colon_1, M>::value;
    static const Integer size2      = colon_func::size<Colon_2, N>::value;
    static const Integer offset1    = colon_func::offset<Colon_1>::value;
    static const Integer offset2    = colon_func::offset<Colon_2>::value;
    static const Integer step1      = colon_func::step<Colon_1>::value;
    static const Integer step2      = colon_func::step<Colon_2>::value;

    using new_array = mkd::sub_array_2<Array_t, offset1, offset2, step1, step2>;
    using type      = ct_matrix<size1,size2,new_array,Deps>;
};

}}}
