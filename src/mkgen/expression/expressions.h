/*
 *  This file is a part of Matrix Computation Library (MATCL)
 *
 *  Copyright (c) Pawe³ Kowal 2019 - 2021
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
#include "mkgen/details/utils/mpl.h"

namespace matcl { namespace mkgen
{

// A + B, where A, B is a ct_matrix or ct_scalar
template<Mat_or_scalar Matrix_1, Mat_or_scalar Matrix_2>
constexpr auto operator+(const Matrix_1& A, const Matrix_2& B) 
{ 
    (void)A; 
    (void)B;
    using ret_type  = typename mkd::mat_plus_impl<Matrix_1, Matrix_2>::type;
    return std::declval<ret_type>();
};

// A - B, where A, B is a ct_matrix or ct_scalar
template<Mat_or_scalar Matrix_1, Mat_or_scalar Matrix_2>
constexpr auto operator-(const Matrix_1& A, const Matrix_2& B) 
{
    (void)A; 
    (void)B;

    using ret_type  = typename mkd::mat_minus_impl<Matrix_1, Matrix_2>::type;
    return std::declval<ret_type>();
};

// matrix multiplication A * B, where A, B is a ct_matrix or ct_scalar
template<Mat_or_scalar Matrix_1, Mat_or_scalar Matrix_2>
constexpr auto operator*(const Matrix_1& A, const Matrix_2& B) 
{ 
    (void)A; 
    (void)B;
    
    using ret_type  = typename mkd::mat_mult_impl<Matrix_1, Matrix_2>::type;
    return std::declval<ret_type>();
};

// element-by-element division A / B, where A, B is a ct_matrix or ct_scalar
template<Mat_or_scalar Matrix_1, Mat_or_scalar Matrix_2>
constexpr auto operator/(const Matrix_1& A, const Matrix_2& B) 
{
    (void)A; 
    (void)B;

    using ret_type  = typename mkd::mat_div_impl<Matrix_1, Matrix_2>::type;
    return std::declval<ret_type>();
};

// multiply k-th element in given row of a matrix A by k-th element of a matrix B,
// where B is rows x 1 matrix, and A is rows x cols
// in Matlab's notation A .* (B * J), J = ones(1, cols)
template<Mat_or_scalar Matrix_1, Mat_or_scalar Matrix_2>
constexpr auto mult_rows(const Matrix_1& A, const Matrix_2& B) 
{
    (void)A; 
    (void)B;

    using ret_type  = typename mkd::mult_rows_impl<Matrix_1, Matrix_2>::type;
    return std::declval<ret_type>();
};

// multiply k-th element in given column of a matrix A by k-th element of a matrix B,
// where B is cols x 1 matrix, and A is rows x cols
// in Matlab's notation A .* (J * B'), J = ones(1, rows), where A has size rows x cols
template<Mat_or_scalar Matrix_1, Mat_or_scalar Matrix_2>
constexpr auto mult_cols(const Matrix_1& A, const Matrix_2& B) 
{
    (void)A; 
    (void)B;

    using ret_type  = typename mkd::mult_cols_impl<Matrix_1, Matrix_2>::type;
    return std::declval<ret_type>();
};

// element-by-element multiplication
template<Mat_or_scalar Matrix_1, Mat_or_scalar Matrix_2>
constexpr auto mul(const Matrix_1& A, const Matrix_2& B) 
{
    (void)A; 
    (void)B;

    using ret_type  = typename mkd::mul_impl<Matrix_1, Matrix_2>::type;
    return std::declval<ret_type>();
};

// -A, where A is a ct_matrix or ct_scalar
template<Mat_or_scalar Matrix_1>
constexpr auto operator-(const Matrix_1& A) 
{
    (void)A; 

    using ret_type  = typename mkd::mat_uminus_impl<Matrix_1>::type;
    return std::declval<ret_type>();
};

// unary function with Tag tag applied to matrix or scalar A
template<class Tag>
struct func_unary
{
    template<Mat_or_scalar Matrix_1>
    static constexpr auto eval(const Matrix_1& A)
    {
        (void)A;
        using ret_type  = typename mkd::func_unary_impl<Tag, Matrix_1>::type;

        return std::declval<ret_type>();
    };
};

// binary element-wise function with tag tag applied to matrices or scalars
// A and B
template<class Tag>
struct func_binary
{
    template<Mat_or_scalar Matrix_1, Mat_or_scalar Matrix_2>
    constexpr auto operator()(const Matrix_1& A, const Matrix_1& B)
    {
        (void)A;
        (void)B;
        using ret_type  = typename mkd::func_binary_impl<Tag, Matrix_1, Matrix_1>::type;

        return std::declval<ret_type>();
    };
};

}}

#include "mkgen/details/expressions/mat_mult_impl.h"
#include "mkgen/details/expressions/mat_plus_impl.h"
#include "mkgen/details/expressions/func_impl.h"