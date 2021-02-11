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

// return the result of the Hankel functions of the first kind:
//     cyc_hankel_1(v,x) = J_v(x) + i * Y_v(x) := H_v^1(x)
// where J_v is the Bessel function of the first kind, and 
// Y_v(x) is the Bessel function of the second kind;
// not available for complex values
template<class S1, class S2, class Enable = typename md::enable_if_scalar2<S1,S2,void>::type>
typename md::unify_types2_promote<S1,S2,Float_complex>::type
                cyl_hankel_1(const S1& v, const S2& x);

// return the result of the Hankel functions of the second kind:
//     cyc_hankel_1(v,x) = J_v(x) - i * Y_v(x) := H_v^1(x)
// where J_v is the Bessel function of the first kind, and 
// Y_v(x) is the Bessel function of the second kind;
// not available for complex values
template<class S1, class S2, class Enable = typename md::enable_if_scalar2<S1,S2,void>::type>
typename md::unify_types2_promote<S1,S2,Float_complex>::type
                cyl_hankel_2(const S1& v, const S2& x);

// return the result of the spherical Hankel functions of the first kind:
//     sph_hankel_1(v,x) = sqrt(pi/2)/sqrt(x) * H_{v+1/2}^1(x)
// where H_v^i is the Hankel function of i-th kind;
// not available for complex values
template<class S1, class S2, class Enable = typename md::enable_if_scalar2<S1,S2,void>::type>
typename md::unify_types2_promote<S1,S2,Float_complex>::type
                sph_hankel_1(const S1& v, const S2& x);

// return the result of the spherical Hankel functions of the first kind:
//     sph_hankel_1(v,x) = sqrt(pi/2)/sqrt(x) * H_{v+1/2}^2(x)
// where H_v^i is the Hankel function of i-th kind;
// not available for complex values
template<class S1, class S2, class Enable = typename md::enable_if_scalar2<S1,S2,void>::type>
typename md::unify_types2_promote<S1,S2,Float_complex>::type
                sph_hankel_2(const S1& v, const S2& x);

};

#include "matcl-specfunc/details/hankel.inl"