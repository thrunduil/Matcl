/*
 *  This file is a part of Matrix Computation Library (MATCL)
 *
 *  Copyright (c) Pawe³ Kowal 2017 - 2018
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

#include "matcl-core/matrix/scalar_types.h"
#include "matcl-internals/container/mat_d.h"
#include "matcl-internals/container/mat_s.h"
#include "matcl-internals/container/mat_b.h"

namespace matcl { namespace raw { namespace details
{

namespace md = matcl::details;

template<class MP>
struct find_helper
{	
    using val_type  = typename MP::value_type;
    using int_mat   = Matrix<Integer,struct_dense>;
    using val_mat   = Matrix<val_type,struct_dense>;

    static void     eval_find(matcl::Matrix& i, const MP& m);
    static void     eval_find(matcl::Matrix& i, const MP& m, const test_function& t);

    static void     eval_find_2(matcl::Matrix& i, matcl::Matrix& j, const MP& m);
    static void     eval_find_2(matcl::Matrix& i, matcl::Matrix& j, const MP& m, 
                                const test_function& t);

    static void     eval_find_3(matcl::Matrix& i, matcl::Matrix& j, matcl::Matrix& x, 
                                const MP& m);
    static void     eval_find_3(matcl::Matrix& i, matcl::Matrix& j, matcl::Matrix& x, 
                                const MP& m, const test_function& t);
};

};};};