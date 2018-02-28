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

//-------------------------------------------------------------------------
//              ELLIPTIC INTEGRALS - CARLSON FORM
//-------------------------------------------------------------------------

// Carlson's elliptic integral of the first kind:
//     ellint_rf(x,y,z) = 1/2 int_0^inf [(t+x)(t+y)(t+z)]^(-1/2)dt
// not available for complex values
template<class S1, class S2, class S3, 
        class Enable = typename md::enable_if_scalar3<S1,S2,S3,void>::type>
typename md::unify_types3_promote<S1,S2,S3,Float>::type
                ellint_rf(const S1& x, const S2& y, const S3& z);

// Carlson's elliptic integral of the second kind:
//     ellint_rd(x,y,z) = 3/2 int_0^inf [(t+x)(t+y)]^(-1/2)(t+z)^(-3/2) dt
// not available for complex values
template<class S1, class S2, class S3, 
        class Enable = typename md::enable_if_scalar3<S1,S2,S3,void>::type>
typename md::unify_types3_promote<S1,S2,S3,Float>::type
                ellint_rd(const S1& x, const S2& y, const S3& z);

// Carlson's elliptic integral of the third kind:
//     ellint_rj(x,y,z,p) = 3/2 int_0^inf 1/(t+p)*[(t+x)(t+y)(t+z)]^(-1/2) dt
// not available for complex values
template<class S1, class S2, class S3, class S4, 
        class Enable = typename md::enable_if_scalar4<S1,S2,S3,S4,void>::type>
typename md::unify_types4_promote<S1,S2,S3,S4,Float>::type
                ellint_rj(const S1& x, const S2& y, const S3& z, const S4& p);

// Carlson's degenerate elliptic integral:
//     ellint_rc(x,y) = 1/2 int_0^inf (t+x)^(-1/2)(t+y)^(-1) dt
// not available for complex values
template<class S1, class S2, class Enable = typename md::enable_if_scalar2<S1,S2,void>::type>
typename md::unify_types2_promote<S1,S2,Float>::type
                ellint_rc(const S1& x, const S2& y);

// Carlson's symmetric integral:
// ellint_rg(x,y,z)= 1/(4pi) int_0^(2pi) int_0^pi sqrt[R(x,y,z,t,p)] sin(t) dt dp
// R(x,y,z,t,p)    = x sin^2(t) cos^2(p) + y sin^2(t) sin^2(p) + z cos^2(t)
// not available for complex values
template<class S1, class S2, class S3, 
        class Enable = typename md::enable_if_scalar3<S1,S2,S3,void>::type>
typename md::unify_types3_promote<S1,S2,S3,Float>::type
                ellint_rg(const S1& x, const S2& y, const S3& z);

//-------------------------------------------------------------------------
//              ELLIPTIC INTEGRALS - LEGENDRE FORM
//-------------------------------------------------------------------------
// incomplete elliptic integral of the first kind:
//     ellint_1(k,phi) = int_0^phi 1/sqrt(1-k^2*sin^2(t)) dt
// not available for complex values
template<class S1, class S2, class Enable = typename md::enable_if_scalar2<S1,S2,void>::type>
typename md::unify_types2_promote<S1,S2,Float>::type
                ellint_1(const S1& k, const S2& phi);

// complete elliptic integral of the first kind:
//     ellint_1(k) = ellint_1(k,pi/2)
// not available for complex values
template<class S1, class Enable = typename md::enable_if_scalar<S1,void>::type>
typename md::unify_types_promote<S1,Float>::type
                ellint_1(const S1& k);

// incomplete elliptic integral of the second kind:
//     ellint_2(k,phi) = int_0^phi sqrt(1-k^2*sin^2(t)) dt
// not available for complex values
template<class S1, class S2, class Enable = typename md::enable_if_scalar2<S1,S2,void>::type>
typename md::unify_types2_promote<S1,S2,Float>::type
                ellint_2(const S1& k, const S2& phi);

// complete elliptic integral of the second kind:
//     ellint_2(k) = ellint_2(k,pi/2)
// not available for complex values
template<class S1, class Enable = typename md::enable_if_scalar<S1,void>::type>
typename md::unify_types_promote<S1,Float>::type
                ellint_2(const S1& k);

// incomplete elliptic integral of the third kind:
//     ellint_3(k,n,phi) = int_0^phi 1/sqrt(1-k^2*sin^2(t))/(1-n*sin^2(t)) dt
// not available for complex values
template<class S1, class S2, class S3, 
        class Enable = typename md::enable_if_scalar3<S1,S2,S3,void>::type>
typename md::unify_types3_promote<S1,S2,S3,Float>::type
                ellint_3(const S1& k, const S2& n, const S3& phi);

// complete elliptic integral of the third kind:
//     ellint_3(k,n) = ellint_3(k,n,pi/2)
// not available for complex values
template<class S1, class S2, class Enable = typename md::enable_if_scalar2<S1,S2,void>::type>
typename md::unify_types2_promote<S1,S2,Float>::type
                ellint_3(const S1& k, const S2& n);

// incomplete elliptic integral D:
//     ellint_d(k,phi) = int_0^phi sin^2(t)/sqrt(1-k^2*sin^2(t)) dt
//                     = 1/k^2 * [ellint_1(k,phi) - ellint_2(k,phi)]
// not available for complex values
template<class S1, class S2, class Enable = typename md::enable_if_scalar2<S1,S2,void>::type>
typename md::unify_types2_promote<S1,S2,Float>::type
                ellint_d(const S1& k, const S2& phi);

// complete elliptic integral D:
//     ellint_d(k) = ellint_d(k,pi/2)
// not available for complex values
template<class S1, class Enable = typename md::enable_if_scalar<S1,void>::type>
typename md::unify_types_promote<S1,Float>::type
                ellint_d(const S1& k);


//-------------------------------------------------------------------------
//                  JACOBI ELLIPTIC FUNCTIONS
//-------------------------------------------------------------------------
// Jacobi elliptic function cn defined as the inverse of incomplete elliptic
// integral of the first kind:
//     jacobi_sn(k,u) = sin(phi), where u = ellint_1(k,phi)
// not available for complex values
template<class S1, class S2, class Enable = typename md::enable_if_scalar2<S1,S2,void>::type>
typename md::unify_types2_promote<S1,S2,Float>::type
                jacobi_sn(const S1& k, const S2& u);

// Jacobi elliptic function cn defined as the inverse of incomplete elliptic
// integral of the first kind:
//     jacobi_cn(k,u) = cos(phi), where u = ellint_1(k,phi)
// not available for complex values
template<class S1, class S2, class Enable = typename md::enable_if_scalar2<S1,S2,void>::type>
typename md::unify_types2_promote<S1,S2,Float>::type
                jacobi_cn(const S1& k, const S2& u);

// Jacobi elliptic function cn defined as the inverse of incomplete elliptic
// integral of the first kind:
//     jacobi_dn(k,u) = sqrt[1-k^2*sin^2(phi)], where u = ellint_1(k,phi)
// not available for complex values
template<class S1, class S2, class Enable = typename md::enable_if_scalar2<S1,S2,void>::type>
typename md::unify_types2_promote<S1,S2,Float>::type
                jacobi_dn(const S1& k, const S2& u);

// return Jacobi elliptic functions jacobi_sn, jacobi_cn, and jacobi_dn; 
// this version is faster than calculating functions sn, cn or dn separately;
// not available for complex values
template<class S1, class S2, class Ret, 
        class Enable = typename md::enable_if_scalar3<S1,S2,Ret,void>::type>
void            jacobi_elliptic(const S1& k, const S2& u, Ret& sn, Ret& cn, Ret& dn);

};

#include "matcl-specfunc/details/elliptic.inl"