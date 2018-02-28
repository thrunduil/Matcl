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

// return beta function element-by-element;
// beta function for scalars a, b is defined as
//     beta(a,b) = gamma(a) * gamma(b) / gamma(a+b), or
//     beta(a,b) = int_0^1 t^(a-1) (1-t)^(b-1) dt
MATCL_SF_EXPORT Matrix      beta(const Matrix& A, const Matrix& B);
MATCL_SF_EXPORT Matrix      beta(Matrix&& A, const Matrix& B);
MATCL_SF_EXPORT Matrix      beta(const Matrix& A, Matrix&& B);
MATCL_SF_EXPORT Matrix      beta(Matrix&& A, Matrix&& B);

template<class S1, class S2, class Enable = typename md::enable_if_scalar2<S1,S2,void>::type>
typename md::unify_types2_promote<S1,S2,Float>::type
                            beta(const S1& n,const S2& x);

// returns the incomplete beta function of a, b and x:
//     ibeta(a,b,x) = int_0^x t^(a-1)*(1-t)^(b-1)dt
template<class S1, class S2, class S3, 
        class Enable = typename md::enable_if_scalar3<S1,S2,S3,void>::type>
typename md::unify_types3_promote<S1,S2,S3,Float>::type
                            ibeta(const S1& a, const S2& b, const S3& x);

// returns the complement of incomplete beta function of a, b and x:
//     ibetac(a,b,x) = 1 - ibeta(a,b,x)
template<class S1, class S2, class S3, 
        class Enable = typename md::enable_if_scalar3<S1,S2,S3,void>::type>
typename md::unify_types3_promote<S1,S2,S3,Float>::type
                            ibetac(const S1& a, const S2& b, const S3& x);

// return the normalised incomplete beta function of a, b and x, given by:
//     ibeta_norm(a,b,x) = 1/beta(a,b) * ibeta(a,b,x)
template<class S1, class S2, class S3, 
        class Enable = typename md::enable_if_scalar3<S1,S2,S3,void>::type>
typename md::unify_types3_promote<S1,S2,S3,Float>::type
                            ibeta_norm(const S1& a, const S2& b, const S3& x);

// return the normalised complement of the incomplete beta function of a, b and x:
//     ibetac_norm(a,b,x) = 1 - ibeta_norm(a,b,x)
template<class S1, class S2, class S3, 
        class Enable = typename md::enable_if_scalar3<S1,S2,S3,void>::type>
typename md::unify_types3_promote<S1,S2,S3,Float>::type
                            ibetac_norm(const S1& a, const S2& b, const S3& x);

// return inverse of the normalised incomplete beta function of a, b and x, satisfying:
//     ibeta_inv(a,b,p) = x <=> ibeta_norm(a,b,x) = p
template<class S1, class S2, class S3, 
        class Enable = typename md::enable_if_scalar3<S1,S2,S3,void>::type>
typename md::unify_types3_promote<S1,S2,S3,Float>::type
                            ibeta_inv(const S1& a, const S2& b, const S3& p);

// return inverse of the normalised complement of incomplete beta function of a, b and x, 
// satisfying:
//     ibetac_inv(a,b,p) = x <=> ibetac_norm(a,b,x) = p
template<class S1, class S2, class S3, 
        class Enable = typename md::enable_if_scalar3<S1,S2,S3,void>::type>
typename md::unify_types3_promote<S1,S2,S3,Float>::type
                            ibetac_inv(const S1& a, const S2& b, const S3& q);

// return inverse of the normalised incomplete beta function of a, b and x,
// with respect to parameter a satisfying:
//     ibeta_inva(b,x,p) = a <=> ibeta_norm(a,b,x) = p
template<class S1, class S2, class S3, 
        class Enable = typename md::enable_if_scalar3<S1,S2,S3,void>::type>
typename md::unify_types3_promote<S1,S2,S3,Float>::type
                            ibeta_inva(const S1& b, const S2& x, const S3& p);

// return inverse of the normalised complement of incomplete beta function 
// of a, b and x, with respect to parameter a satisfying:
//     ibetac_inva(b,x,p) = a <=> ibetac_norm(a,b,x) = p
template<class S1, class S2, class S3, 
        class Enable = typename md::enable_if_scalar3<S1,S2,S3,void>::type>
typename md::unify_types3_promote<S1,S2,S3,Float>::type
                            ibetac_inva(const S1& b, const S2& x, const S3& q);

// return inverse of the normalised incomplete beta function of a, b and x,
// with respect to parameter b satisfying:
//     ibeta_invb(a,x,p) = b <=> ibeta_norm(a,b,x) = p
template<class S1, class S2, class S3, 
        class Enable = typename md::enable_if_scalar3<S1,S2,S3,void>::type>
typename md::unify_types3_promote<S1,S2,S3,Float>::type
                            ibeta_invb(const S1& a, const S2& x, const S3& p);

// return inverse of the normalised complement of incomplete beta function of a, b and x,
// with respect to parameter b satisfying:
//     ibetac_invb(a,x,p) = b <=> ibetac_norm(a,b,x) = p
template<class S1, class S2, class S3, 
        class Enable = typename md::enable_if_scalar3<S1,S2,S3,void>::type>
typename md::unify_types3_promote<S1,S2,S3,Float>::type
                            ibetac_invb(const S1& a, const S2& x, const S3& q);

// return partial derivative of inverse of the normalised incomplete beta function
// of a, b and x, with respect to x, satisfying:
//     ibeta_dif(a,b,x) = d/dx[ibeta_norm(a,b,x)] = (1-x)^(b-1) * x^(a-1) / beta(a,b)
template<class S1, class S2, class S3, 
        class Enable = typename md::enable_if_scalar3<S1,S2,S3,void>::type>
typename md::unify_types3_promote<S1,S2,S3,Float>::type
                            ibeta_dif(const S1& a, const S2& b, const S3& x);
};

#include "matcl-specfunc/details/beta.inl"