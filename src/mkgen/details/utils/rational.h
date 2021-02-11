/*
 *  This file is a part of Matrix Computation Library (MATCL)
 *
 *  Copyright (c) Pawe³ Kowal 2019 - 2021
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

#include "mkgen/mkgen_fwd.h"

namespace matcl { namespace mkgen { namespace details
{

constexpr Integer const_abs(Integer v)
{
    return (v < 0) ? -v : v;
};

//----------------------------------------------------------------------------------
//                              GDC
//----------------------------------------------------------------------------------
template<Integer A, Integer B>
struct gdc_impl
{
    static const Integer value          = gdc_impl<B, A % B>::value;
};

template<Integer A>
struct gdc_impl<A, 0>
{
    static const Integer value          = A;
};

template<Integer A>
struct gdc_impl<A,1>
{
    static const Integer value          = 1;
};

template<Integer A, Integer B>
struct gdc
{
    static const Integer A2             = const_abs(A);
    static const Integer B2             = const_abs(B);
    static const Integer value          = gdc_impl<B2, A2 % B2>::value;
};

//----------------------------------------------------------------------------------
//                              N1/D1 + N2/D2
//----------------------------------------------------------------------------------
template<Integer N1, Integer D1, Integer N2, Integer D2>
struct rational_plus
{
    static const Integer nom            = N1*D2 + N2*D1;
    static const Integer den            = D1*D2;
    static const Integer d              = gdc<nom,den>::value;

    static const Integer nominator      = nom/d;
    static const Integer denominator    = den/d;
};

//----------------------------------------------------------------------------------
//                              N1/D1 * N2/D2
//----------------------------------------------------------------------------------
template<Integer N1, Integer D1, Integer N2, Integer D2>
struct rational_mult
{
    static const Integer nom            = N1*N2;
    static const Integer den            = D1*D2;
    static const Integer d              = gdc<nom,den>::value;

    static const Integer nominator      = nom/d;
    static const Integer denominator    = den/d;
};

}}}