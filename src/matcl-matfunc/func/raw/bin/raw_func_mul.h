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

#include "matcl-matfunc/func/raw/bin/raw_func_helpers.h"

namespace matcl { namespace raw { namespace details
{

template<class M1,class M2>
struct mul_helper_mat_mat_inpl
{
    using ret_type_mul = typename ret_type_constructor<M1,M2,1,1,1,false,false>::type;
    static void eval(matcl::Matrix& ret, const M1& A, const M2& B);
};

}}}
