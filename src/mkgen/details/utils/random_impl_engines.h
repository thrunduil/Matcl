/*
 *  This file is a part of Matrix Computation Library (MATCL)
 *
 *  Copyright (c) Pawe³ Kowal 2019
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

#include <cstdint>
#include "matcl-core/details/mpl.h"

namespace matcl { namespace mkgen { namespace details
{

//----------------------------------------------------------------------------------
//                              linear congruential engine
//----------------------------------------------------------------------------------
// generate random numbers according to
//  v(i) = (a * v(i-1) + c) % m
//  v(0) = seed
template<uint64_t seed, uint64_t a, uint64_t c, uint64_t m>
struct LCE
{
    using Unsigned                  = uint64_t;
    static const Unsigned value     = seed;
    static const Unsigned maxvalue  = m-1;
};

template<uint64_t seed, uint64_t a, uint64_t c, uint64_t m>
struct random_impl::eval<LCE<seed, a, c, m>>
{
    static const uint64_t value     = (a * seed + c) % m;
    static const uint64_t new_seed  = value;

    using type  = LCE<new_seed, a, c, m>;
};

template<uint64_t seed, uint64_t a, uint64_t c>
struct random_impl::eval<LCE<seed, a, c, 0>>
{
    static const uint64_t value     = a * seed + c;
    static const uint64_t new_seed  = value;

    using type  = LCE<new_seed, a, c, 0>;
};

template<uint64_t seed, uint64_t a, uint64_t c, uint64_t m>
struct random_impl::init<LCE<seed, a, c, m>>
{
    using type  = typename eval<LCE<seed,a,c,m>>::type;
};

}}}