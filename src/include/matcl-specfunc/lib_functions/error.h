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

// error function of the elements of given matrix m;
// erf function for scalar x is defined as
//     erf(x) = 2/sqrt(pi) * int_0^x exp(-t^2) dt
MATCL_SF_EXPORT Matrix      erf(const Matrix &m);

template<class S1, class Enable = typename md::enable_if_scalar<S1,void>::type>
typename md::unify_types_promote<S1,Float>::type
                            erf(const S1& n);

// complementary error function of the elements of given matrix m;
// erfc function for scalars is defined as
//     erfc(x) = 1 - erf(x)
MATCL_SF_EXPORT Matrix      erfc(const Matrix &m);

template<class S1, class Enable = typename md::enable_if_scalar<S1,void>::type>
typename md::unify_types_promote<S1,Float>::type
                            erfc(const S1& n);

// inverse of error function of the elements of given matrix m;
// erf_inv function for scalar x is defined as
//     erf_inv(x) = y, where y satisfies erf(y) = x
// function is not defined for complex values
MATCL_SF_EXPORT Matrix      erf_inv(const Matrix &m);

template<class S1, class Enable = typename md::enable_if_scalar<S1,void>::type>
typename md::unify_types_promote<S1,Float>::type
                            erf_inv(const S1& n);

// inverse of complementary error function of the elements of given matrix m;
// erfc function for scalars is defined as
//     erfc_inv(x) = y, where y satisfies erfc(y) = x
// function is not defined for complex values
MATCL_SF_EXPORT Matrix      erfc_inv(const Matrix &m);

template<class S1, class Enable = typename md::enable_if_scalar<S1,void>::type>
typename md::unify_types_promote<S1,Float>::type
                            erfc_inv(const S1& n);
};

#include "matcl-specfunc/details/error.inl"