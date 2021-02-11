/*
 *  This file is a part of Matrix Computation Library (MATCL)
 *
 *  Copyright (c) Pawe³ Kowal 2018 - 2021
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

#include "bin_function.h"

namespace matcl { namespace test
{

Real bin_function::eval(const Matrix& mat1,const Matrix& mat2, int code)
{
    m_is_error = false;
    m_error.clear();

    return eval_mat(mat1,mat2,code);
};
Real bin_function::eval_scal(const Scalar& s1, const Scalar& s2, int code)
{
    m_is_error = false;
    m_error.clear();
    return eval_scalar(s1,s2,code);
};

Real bin_function::check_value_type(const Matrix& mat1,const Matrix& mat2,const Matrix& out, 
                                    bool float_return)
{
    value_code vt1  = matrix_traits::unify_value_types(mat1.get_value_code(),mat2.get_value_code());
    vt1             = matrix_traits::real_value_type(vt1);

    if (float_return == true)
        vt1         = matrix_traits::unify_value_types(vt1, value_code::v_float);

    value_code vt2  = matrix_traits::real_value_type(out.get_value_code());

    if (vt1 != vt2)
        return 1;
    else
        return 0;
};

Real bin_function::check_value_type(const Matrix& mat1,const Matrix& mat2)
{
    bool is_int1 = (mat1.get_value_code() == value_code::v_integer);
    bool is_int2 = (mat2.get_value_code() == value_code::v_integer);

    if ((is_int1 == true && is_int2 == false) || (is_int2 == true && is_int1 == false))
        return 1;
    else
        return 0;
};

};};