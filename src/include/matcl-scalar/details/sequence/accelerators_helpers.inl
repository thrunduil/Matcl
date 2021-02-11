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

#include "matcl-scalar/details/sequence/accelerators_helpers.h"

namespace matcl { namespace seq_helpers
{

template<class Float>
Float error_minus(const Float& a, const Float& b, const Float& abs_xmy)
{
    // it is assumed that a and b are exact

    // Sterbez lemma
    if ( b <= Float(2) * a && b >= Float(0.5) * a)
        return Float(0);
    else
        return abs_xmy;
}

template<class Float>
Float error_plus(const Float& x, const Float& y, const Float& abs_xpy)
{
    return error_minus(x, -y, abs_xpy);
};

}}

