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

// evaluate the Airy function Ai which is the first solution to the differential
// equation:
//         d^2 w/dz^2 = zw
// not available for complex values
MATCL_SF_EXPORT Matrix  airy_ai(const Matrix &m);

template<class S1, class Enable = typename md::enable_if_scalar<S1,void>::type>
typename md::unify_types_promote<S1,Float>::type
                        airy_ai(const S1& x);

// evaluate the Airy function Bi which is the second solution to the differential
// equation:
//         d^2 w/dz^2 = zw
// not available for complex values
MATCL_SF_EXPORT Matrix  airy_bi(const Matrix &m);

template<class S1, class Enable = typename md::enable_if_scalar<S1,void>::type>
typename md::unify_types_promote<S1,Float>::type
                        airy_bi(const S1& x);

// evaluate differential of the Airy function of the first kind
// not available for complex values
MATCL_SF_EXPORT Matrix  airy_ai_dif(const Matrix &m);

template<class S1, class Enable = typename md::enable_if_scalar<S1,void>::type>
typename md::unify_types_promote<S1,Float>::type
                        airy_ai_dif(const S1& x);

// evaluate differential of the Airy function of the second kind
// not available for complex values
MATCL_SF_EXPORT Matrix  airy_bi_dif(const Matrix &m);

template<class S1, class Enable = typename md::enable_if_scalar<S1,void>::type>
typename md::unify_types_promote<S1,Float>::type
                        airy_bi_dif(const S1& x);

// the Airy Ai and Bi functions have an infinite number of zeros on the negative
// real axis; return m-th zero of the Airy Ai function (m >= 1)
template<class S1, class Enable = typename md::enable_if_scalar<S1,void>::type>
typename md::unify_types_promote<S1,Float>::type
                        airy_ai_zero(Integer m);
Real                    airy_ai_zero(Integer m);
Float                   fairy_ai_zero(Integer m);

// version of airy_ai_zero; for each element in the matrix mat evaluate airy_ai_zero
// at this argument; mat is converted to integer dense matrix
MATCL_SF_EXPORT Matrix  airy_ai_zero(const Matrix &mat);

// the Airy Ai and Bi functions have an infinite number of zeros on the negative
// real axis; return m-th zero of the Airy Bi function (m >= 1)
template<class S1, class Enable = typename md::enable_if_scalar<S1,void>::type>
typename md::unify_types_promote<S1,Float>::type
                        airy_bi_zero(Integer m); 
Real                    airy_bi_zero(Integer m);
Float                   fairy_bi_zero(Integer m);

// version of airy_ai_zero; for each element in the matrix mat evaluate airy_bi_zero
// at this argument; mat is converted to integer dense matrix
MATCL_SF_EXPORT Matrix  airy_bi_zero(const Matrix &m);

};

#include "matcl-specfunc/details/airy.inl"
