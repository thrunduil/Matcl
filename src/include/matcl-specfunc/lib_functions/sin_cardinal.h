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

#include "matcl-specfunc/general/config.h"
#include "matcl-matrep/matrix/matrix.h"

namespace matcl
{

namespace md = matcl::details;

// return the sinus cardinal function:
//     sinc(x) = sin(x) / x
template<class S1, class Enable = typename md::enable_if_scalar<S1,void>::type>
typename md::unify_types_promote<S1,Float>::type
                                sinc(const S1& x);

// version for a matrix argument, evaluate sinc for each element in A
MATCL_SF_EXPORT Matrix          sinc(const Matrix& A);

// return the hyperbolic sinus cardinal function:
//     sinhc(x) = sinh(x) / x
template<class S1, class Enable = typename md::enable_if_scalar<S1,void>::type>
typename md::unify_types_promote<S1,Float>::type
                                sinhc(const S1& x);

// version for a matrix argument, evaluate sinhc for each element in A
MATCL_SF_EXPORT Matrix          sinhc(const Matrix& A);

};

#include "matcl-specfunc/details/sin_cardinal.inl"
