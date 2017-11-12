/*
 *  This file is a part of Matrix Computation Library (MATCL)
 *
 *  Copyright (c) Pawe³ Kowal 2017
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

#include "matcl-core/float/fma_dekker.h"

//simd versions are defined in "matcl-core/details/float/fma_dekker_simd.inl"

namespace matcl { namespace details
{

// represent z = x * y  as z = xy + err (exactly)
// requires: abs(x) <= 1e995, abs(y) <= 1e995, abs(x*y) <= 1e1021
// ensures: if x*y == 0 || abs(x*y) >= 1e-969, then x*y = xy + err
// modified version of Dekker function (author: Sylvie Boldo, available at
//  http://proval.lri.fr/gallery/Dekker.en.html)
template<class T>
inline
void twofold_mult_dekker_double(const T& x, const T& y, T& xy, T& err)
{
    xy          = x * y;

    double C0   = (double)(0x8000001); //2^27 + 1
    T C         = T(C0);

    T px        = x * C;
    T qx        = x - px;
    T hx        = px + qx;
    T tx        = x - hx;

    T py        = y * C;
    T qy        = y - py;
    T hy        = py + qy;
    T ty        = y - hy;

    err         = (hx * hy) - xy;
    
    err         += hx*ty;
    err         += hy*tx;
    err         += tx*ty;
};

}}

namespace matcl
{

//-----------------------------------------------------------------------
//                   COMPATIBILITY FUNCTIONS
//-----------------------------------------------------------------------
force_inline
void matcl::twofold_mult_dekker(double x, double y, double& xy, double& err)
{
    return details::twofold_mult_dekker_double(x, y, xy, err);
};

force_inline
float matcl::fma_dekker(float x, float y, float z)
{
    //for floats we can use double arithmetics
    double xd   = x;
    double yd   = y;
    double zd   = z;

    // res is exact
    double res  = xd * yd;

    // 0.5 ulp error in double precision
    res         = res + zd;

    // almost 0.5 ulp error in single precision
    return (float)res;
};

}
