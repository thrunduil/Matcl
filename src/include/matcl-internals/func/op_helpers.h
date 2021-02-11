/*
 *  This file is a part of Matrix Computation Library (MATCL)
 *
 *  Copyright (c) Pawe³ Kowal 2017 - 2021
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

#include "matcl-matrep/lib_functions/manip.h"

#pragma warning(push)
#pragma warning(disable: 4127)  // conditional expression is constant

namespace matcl { namespace details
{

template<class Val, bool Convert_int = true>
void to_real_mat(matcl::Matrix& ret, const Matrix& A, const Val& val, 
                 matcl::dynamic::function_name obj_fun)
{
    matcl::value_code vt = A.get_value_code();

    if (std::is_same<Val,Object>::value == true || vt == value_code::v_object)
    {
        using ti_type_1 = matcl::ti::ti_object;
        using ti_type_2 = typename ti::get_ti_type<Val>::type;

        ti_type_1 t1    = A.get_type();
        ti_type_2 t2    = matcl::ti::get_ti(val);

        ti_type_1 ti_r  = matcl::ti::get_return_ti<ti_type_1>(obj_fun,t1,t2);

        ret = convert_object(A, ti_r);
        return;
    }

    if (vt == value_code::v_integer)
    {
        if (Convert_int == false && std::is_same<Val,Integer>::value == true)
        {
            ret = A;
            return;
        };

        matcl::struct_code st   = A.get_struct_code();
        matcl::mat_code mt      = matrix_traits::get_matrix_type(value_code::v_real,st);
        ret                     = convert(A,mt);
        return;
    }
    else if ((std::is_same<Val,Float>::value || std::is_same<Val,Float_complex>::value)
                && (vt == value_code::v_float || vt == value_code::v_float_complex) )
    {
        ret = A;
        return;
    }
    else if (vt == value_code::v_real || vt == value_code::v_complex)
    {
        ret = A;
        return;
    }
    else if (vt == value_code::v_float)
    {
        //convert to real matrix
        matcl::struct_code st   = A.get_struct_code();
        matcl::mat_code mt      = matrix_traits::get_matrix_type(value_code::v_real,st);
        ret                     = convert(A,mt);
        return;
    }
    else
    {
        //convert to complex matrix
        matcl::struct_code st   = A.get_struct_code();
        matcl::mat_code mt      = matrix_traits::get_matrix_type(value_code::v_complex,st);
        ret                     = convert(A,mt);
        return;
    };
};

template<bool is_obj>
struct select_ti_type
{
    using type          = ti::ti_empty;
    using value_type    = Integer;
};

template<>
struct select_ti_type<true>
{
    using type          = ti::ti_object;
    using value_type    = Object;
};

template<bool is_int, bool is_obj>
struct correct_int_val
{
    static Matrix eval(ti::ti_empty, const Matrix& mat)
    {
        if (mat.get_value_code() == value_code::v_integer)
        {
            matcl::value_code vt    = value_code::v_real;
            matcl::struct_code st   = mat.get_struct_code();
            matcl::mat_code mt      = matrix_traits::get_matrix_type(vt,st);

            return convert(mat,mt);
        };
        return Matrix(mat);
    };
};

template<>
struct correct_int_val<true,false>
{
    static const Matrix& eval(ti::ti_empty , const Matrix& mat)
    {
        return mat;
    };

    static Matrix eval(ti::ti_empty , Matrix&& mat)
    {
        return Matrix(std::move(mat));
    };
};

template<bool is_int>
struct correct_int_val<is_int,true>
{
    static Matrix eval(ti::ti_object ret_ti, const Matrix& mat)
    {
        return convert_object(mat,ret_ti);
    };
};

}}

#pragma warning(pop)