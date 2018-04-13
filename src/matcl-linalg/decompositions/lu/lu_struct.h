/*
 *  This file is a part of Matrix Computation Library (MATCL)
 *
 *  Copyright (c) Pawe³ Kowal 2018
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

#include "matcl-linalg/decompositions/lu.h"
#include "matcl-matrep/matrix/permvec.h"

#include "matcl-internals/container/mat_s.h"
#include "matcl-internals/container/mat_d.h"
#include "matcl-internals/container/mat_b.h"

namespace matcl { namespace details
{

template<class value_type, class struct_type>
struct lu_diag
{
    using VR    = typename details::real_type<value_type>::type;
    using Mat   = raw::Matrix<value_type,struct_type>;
    static void eval(lu_return_type& ret, const Mat& A, opt::linsolve::pivot_type piv, VR TOL);
};

void lu_tril_nonunit(lu_return_type& ret, const matcl::Matrix& A);
void lu_triu(lu_return_type& ret, const matcl::Matrix& A);

};};