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

#include "mkgen/details/mkgen_fwd.h"
#include "mkgen/details/utils/has_function.h"
#include "matcl-core/details/mpl.h"

namespace matcl { namespace mkgen { namespace details
{

//-----------------------------------------------------------------------
//                      check_scalar_data
//-----------------------------------------------------------------------
// check if Data parameter supplied to ct_scalar is valid
template<class Array_t>
struct check_valid_matrix_array
{
    static const bool is_array  = std::is_base_of<mkd::matrix_array<Array_t>, Array_t>::value;

    static_assert(is_array == true, "Array_t is not matrix_array<>");

    using type  = typename Array_t::template check<void>;
};


}}}
