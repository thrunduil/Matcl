/*
 *  This file is a part of Matrix Computation Library (MATCL)
 *
 *  Copyright (c) Pawe� Kowal 2019
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

namespace mk = matcl :: mkgen;

// return true if T can be passed to ct_scalar<T, ...>;
template<class T>
struct is_valid_scalar_data
{
    static const bool value = std::is_base_of<mkd::scalar_data<T>, T>::value;
};

//-----------------------------------------------------------------------
//                      check_computation_tag
//-----------------------------------------------------------------------
// check tag supplied to compute function
template<class Tag>
struct check_computation_tag
{
    // TODO: collect requirements
    using type = void;
};

}}}
