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

#include "matcl-specfunc/general/config.h"
#include "matcl-matrep/matrix/matrix.h"

namespace matcl
{

namespace md = matcl::details;

// return Riemann zeta function for every elements in the matrix A;
// zeta function for scalar x is defined as
//     zeta(s) = sum_{k=1}^infty 1/k^s
// //not available for complex values
MATCL_SF_EXPORT Matrix          zeta(const Matrix& A);

template<class S1, class Enable = typename md::enable_if_scalar<S1,void>::type>
typename md::unify_types_promote<S1,Float>::type
                                zeta(const S1& x);

};

#include "matcl-specfunc/details/zeta.inl"
