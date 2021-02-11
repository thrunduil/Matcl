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

#include "matcl-linalg/utils/utils.h"
#include "matcl-linalg/general/linalg_exception.h"
#include "matcl-matrep/matcl_matrep.h"

namespace matcl { namespace details
{

// not available if vt is v_object
Matrix make_nan_matrix(Integer rows, Integer cols, value_code vt);

template<class Val>
Matrix make_nan_matrix(Integer rows, Integer cols);

template<class Val>
inline Matrix make_nan_matrix_not_object(Integer rows, Integer cols, const std::string&)
{
    return make_nan_matrix<Val>(rows,cols);
};

template<>
inline Matrix make_nan_matrix_not_object<Object>(Integer, Integer, const std::string& func_name)
{
    throw error::object_value_type_not_allowed(func_name);
};

inline Matrix convert_value(const Matrix& A, value_code vc)
{
    value_code vc_A = A.get_value_code();
    value_code vc_r = matrix_traits::real_value_type(vc);

    if (vc_A == vc || vc_A == vc_r)
        return A;

    value_code vc_ret   = matrix_traits::unify_value_types(vc_A, vc_r);
    struct_code sc      = A.get_struct_code();
    mat_code mc         = matrix_traits::get_matrix_type(vc_ret, sc);

    return matcl::convert(A, mc);
};

void check_tridiag(const Matrix& D, const Matrix& E);

};}
