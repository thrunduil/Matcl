/*
 *  This file is a part of Matrix Computation Library (MATCL)
 *
 *  Copyright (c) Pawe� Kowal 2017 - 2021
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

// returns the exponential integral:
//     expint(x,n) = int_1^inf exp(-x*t)/t^n dt
// for n >= 0;
// not available for complex values
template<class S1, class Enable = typename md::enable_if_scalar<S1,void>::type>
typename md::unify_types_promote<S1,Float>::type
                                expint(const S1& x, Integer n);

// returns the exponential integral:
//     expint(x) = int_{-x}^inf exp(-t)/t dt
// not available for complex values
template<class S1, class Enable = typename md::enable_if_scalar<S1,void>::type>
typename md::unify_types_promote<S1,Float>::type
                                expint(const S1& x);

// version for a matrix argument, evaluate expint for each element in A
MATCL_SF_EXPORT Matrix          expint(const Matrix& A);

};

#include "matcl-specfunc/details/exp_integrals.inl"