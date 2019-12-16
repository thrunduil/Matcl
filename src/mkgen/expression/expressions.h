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

#include <iosfwd>

#include "mkgen/matrix/matrix.h"

namespace matcl { namespace mkgen { namespace details
{

template<class M1, class M2>
struct enable_expr_2 :
    public md::enable_if
            <	(is_scalar<M1>::value || is_matrix<M1>::value)
                && (is_scalar<M2>::value || is_matrix<M2>::value),
                const void*
            >
{};

}}}

namespace matcl { namespace mkgen
{

template<class Matrix_1, class Matrix_2,
        class Enable = typename mkd::enable_expr_2<Matrix_1, Matrix_2>::type>
constexpr auto operator+(const Matrix_1&, const Matrix_2&) 
{ 
    using ret_type  = typename mat_plus<Matrix_1, Matrix_2>::type;
    return std::declval<ret_type>();
};

template<class Matrix_1, class Matrix_2,
        class Enable = typename mkd::enable_expr_2<Matrix_1, Matrix_2>::type>
constexpr auto operator-(Matrix_1, Matrix_2) 
{ 
    using ret_type  = typename mat_minus<Matrix_1, Matrix_2>::type;
    return std::declval<ret_type>();
};

template<class Matrix_1, class Matrix_2,
        class Enable = typename mkd::enable_expr_2<Matrix_1, Matrix_2>::type>
constexpr auto operator*(Matrix_1, Matrix_2) 
{ 
    using ret_type  = typename mat_mult<Matrix_1, Matrix_2>::type;
    return std::declval<ret_type>();
};

template<class Matrix_1, class Matrix_2,
        class Enable = typename mkd::enable_expr_2<Matrix_1, Matrix_2>::type>
constexpr auto operator/(Matrix_1, Matrix_2) 
{ 
    using ret_type  = typename mat_div<Matrix_1, Matrix_2>::type;
    return std::declval<ret_type>();
};

#if 0
//------------------------------------------------------------------------------
//                      mult elem by elem
//------------------------------------------------------------------------------
template<Integer M1, Integer N1, class Array1, class Deps1,
         Integer M2, Integer N2, class Array2, class Deps2>
auto        mult_rows(ct_matrix<M1,N1,Array1,Deps1>,ct_matrix<M2,N2,Array2,Deps2>)
                                        -> typename make_mult_rows<ct_matrix<M1,N1,Array1,Deps1>,
                                                    ct_matrix<M2,N2,Array2,Deps2>>::type;

template<Integer M1, Integer N1, class Array1, class Deps1,
         Integer M2, Integer N2, class Array2, class Deps2>
auto        mult_cols(ct_matrix<M1,N1,Array1,Deps1>,ct_matrix<M2,N2,Array2,Deps2>)
                                        -> typename make_mult_cols<ct_matrix<M1,N1,Array1,Deps1>,
                                                    ct_matrix<M2,N2,Array2,Deps2>>::type;

template<Integer M1, Integer N1, class Array1, class Deps1,
         Integer M2, Integer N2, class Array2, class Deps2>
auto        mult(ct_matrix<M1,N1,Array1,Deps1>,ct_matrix<M2,N2,Array2,Deps2>)
                                        -> typename make_mult_mat<ct_matrix<M1,N1,Array1,Deps1>,
                                                    ct_matrix<M2,N2,Array2,Deps2>>::type;

//------------------------------------------------------------------------------
//                      call
//------------------------------------------------------------------------------
template<class Tag, template<class Arg> class Func,
         Integer M1, Integer N1, class Array1, class Deps1>
auto        call_inline(ct_matrix<M1, N1, Array1, Deps1>)
                                        -> typename make_call_inline<Tag, Func, 
                                                ct_matrix<M1, N1, Array1, Deps1>> :: type;

template<class Tag, class Func,
         Integer M1, Integer N1, class Array1, class Deps1>
auto        call_external(ct_matrix<M1, N1, Array1, Deps1>)
                                        -> typename make_call_external<Tag, Func, 
                                                ct_matrix<M1, N1, Array1, Deps1>> :: type;

#endif
}}