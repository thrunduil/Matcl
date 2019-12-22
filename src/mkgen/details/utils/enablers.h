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

#include "mkgen/mkgen_fwd.h"
#include "mkgen/details/utils/mpl.h"
#include "mkgen/matrix/scalar.h"
#include "mkgen/matrix/matrix.h"

namespace matcl { namespace mkgen { namespace details
{

//----------------------------------------------------------------------------------
//                              enablers
//----------------------------------------------------------------------------------

// enable when M1 and M2 are matrix or scalar
template<class M1, class M2>
struct enable_matscal_2 :
    public md::enable_if
            <	(is_scalar<M1>::value || is_matrix<M1>::value)
                && (is_scalar<M2>::value || is_matrix<M2>::value),
                const void*
            >
{};

// enable when M1 is matrix or scalar
template<class M1, class M2>
struct enable_matscal_1 :
    public md::enable_if
            <	(is_scalar<M1>::value || is_matrix<M1>::value),
                const void*
            >
{};

}}}

