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

//-------------------------------------------------------------------------
//                  ORDINARY BESSEL FUNCTIONS
//-------------------------------------------------------------------------

// Bessel functions are solution to Bessel's ordinary differential equation:
// 
//     z^2 d^2u/dz^2 + z du/dz + (z^2 - v^2)u = 0
// 
// where v is the order of the equation. There are two linearly independent 
// solutions known as a Bessel function of the first kind:
//     J_v(z) = (1/2z)^v sum_{k=0}^infty (-1/4z^2)^k / k! / gamma(v+k+1)
// and a Bessel function of the second kind (or a Neumann function):
//     Y_v(z) = ( J_v(z)cos(v*pi) - J_{-v}(z) ) / sin(v*pi)

namespace md = matcl::details;

// return the result of the Bessel functions of the first kind:
//     cyl_bessel_j(v,x) = J_v(x)
// not available for complex values
template<class S1, class S2, class Enable = typename md::enable_if_scalar2<S1,S2,void>::type>
typename md::unify_types2_promote<S1,S2,Float>::type
                    cyl_bessel_j(const S1& v, const S2& x);

// return the result of the Bessel functions of the second kind:
//     cyl_neumann(v,x) = Y_v(x)
// not available for complex values
template<class S1, class S2, class Enable = typename md::enable_if_scalar2<S1,S2,void>::type>
typename md::unify_types2_promote<S1,S2,Float>::type
                    cyl_neumann(const S1& v, const S2& x);

// for every real order ? cylindrical Bessel and Neumann functions 
// have an infinite number of zeros on the positive real axis; this
// function returns m-th zero (m >= 1) on positive real axis of the
// Bessel function of the first kind;
// not available for complex values
template<class S1, class Enable = typename md::enable_if_scalar<S1,void>::type>
typename md::unify_types_promote<S1,Float>::type
                    cyl_bessel_j_zero(const S1& v, Integer m);

// return m-th zero (m >= 1) on positive real axis of the Bessel function
// of the second kind;
// not available for complex values
template<class S1, class Enable = typename md::enable_if_scalar<S1,void>::type>
typename md::unify_types_promote<S1,Float>::type
                    cyl_neumann_zero(const S1& v, Integer m);

// return the first derivative with respect to x of the Bessel function
// of the first kind;
// not available for complex values
template<class S1, class S2, class Enable = typename md::enable_if_scalar2<S1,S2,void>::type>
typename md::unify_types2_promote<S1,S2,Float>::type
                    cyl_bessel_j_dif(const S1& v, const S2& x);

// return the first derivative with respect to x of the Bessel function
// of the second kind;
// not available for complex values
template<class S1, class S2, class Enable = typename md::enable_if_scalar2<S1,S2,void>::type>
typename md::unify_types2_promote<S1,S2,Float>::type
                    cyl_neumann_dif(const S1& v, const S2& x);

//-------------------------------------------------------------------------
//                  MODIFIED BESSEL FUNCTIONS
//-------------------------------------------------------------------------

// Modified Bessel functions are solution to Bessel's ordinary differential
// equation:
// 
//     z^2 d^2u/dz^2 + z du/dz + (z^2 + v^2)u = 0
// 
// where v is the order of the equation. There are two linearly independent 
// solutions known as the modified Bessel function of the first kind:
//     I_v(z) = (1/2z)^v sum_{k=0}^infty (1/4z^2)^k / k! / gamma(v+k+1)
// and the modified Bessel function of the second kind:
//     K_v(z) = pi/2 * ( I_{-v}(z) - J_{v}(z) ) / sin(v*pi)

// return the result of the modified Bessel functions of the first kind:
//     cyl_bessel_i(v,x) = I_v(x)
// not available for complex values
template<class S1, class S2, class Enable = typename md::enable_if_scalar2<S1,S2,void>::type>
typename md::unify_types2_promote<S1,S2,Float>::type
                    cyl_bessel_i(const S1& v, const S2& x);

// return the result of the modified Bessel functions of the second kind:
//     cyl_bessel_k(v,x) = K_v(x)
// not available for complex values
template<class S1, class S2, class Enable = typename md::enable_if_scalar2<S1,S2,void>::type>
typename md::unify_types2_promote<S1,S2,Float>::type
                    cyl_bessel_k(const S1& v, const S2& x);

// return the first derivative with respect to x of the modified Bessel function
// function of the first kind;
// not available for complex values
template<class S1, class S2, class Enable = typename md::enable_if_scalar2<S1,S2,void>::type>
typename md::unify_types2_promote<S1,S2,Float>::type
                    cyl_bessel_i_dif(const S1& v, const S2& x);

// return the first derivative with respect to x of the modified Bessel function
// function of the second kind;
// not available for complex values
template<class S1, class S2, class Enable = typename md::enable_if_scalar2<S1,S2,void>::type>
typename md::unify_types2_promote<S1,S2,Float>::type
                    cyl_bessel_k_dif(const S1& v, const S2& x);

//-------------------------------------------------------------------------
//                  SPHERICAL BESSEL FUNCTIONS
//-------------------------------------------------------------------------
// When solving the Helmholtz equation in spherical coordinates by separation
// of variables, the radial equation has the form:
// 
//     z^2 d^2u/dz^2 + 2z du/dz + (z^2 - n(n+1))u = 0
// 
// where n is the order of the equation. There are two linearly independent 
// solutions called as the spherical Bessel function of the first kind:
//     j_n(z) = sqrt(pi/(2z)) J_{n+1/2}(z)
// and the spherical Bessel function of the second kind (also as spherical 
// Neumann function):
//     y_n(z) = sqrt(pi/(2z)) Y_{n+1/2}(z)

// return the result of the spherical Bessel functions of the first kind:
//     sph_bessel(n,x) = j_n(x)
// not available for complex values
template<class S1, class Enable = typename md::enable_if_scalar<S1,void>::type>
typename md::unify_types_promote<S1,Float>::type
                    sph_bessel(Integer n, const S1& x);

// return the result of the spherical Bessel functions of the second kind:
//     sph_neumann(n,x) = y_n(x)
// not available for complex values
template<class S1, class Enable = typename md::enable_if_scalar<S1,void>::type>
typename md::unify_types_promote<S1,Float>::type
                    sph_neumann(Integer n, const S1& x);

// return the first derivative with respect to x of the spherical Bessel function
// function of the first kind;
// not available for complex values
template<class S1, class Enable = typename md::enable_if_scalar<S1,void>::type>
typename md::unify_types_promote<S1,Float>::type
                    sph_bessel_dif(Integer n, const S1& x);

// return the first derivative with respect to x of the spherical Bessel function
// function of the second kind;
// not available for complex values
template<class S1, class Enable = typename md::enable_if_scalar<S1,void>::type>
typename md::unify_types_promote<S1,Float>::type
                    sph_neumann_dif(Integer n, const S1& x);

};

#include "matcl-specfunc/details/bessel.inl"