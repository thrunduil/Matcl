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

#include <type_traits>
#include <concepts>

#include "matcl-core/matrix/scalar_types.h"
#include "mkgen/matrix/base_types.h"

//------------------------------------------------------------------------------
//                      fordward declarations
//------------------------------------------------------------------------------
namespace matcl { namespace mkgen
{

    namespace details
    {
        namespace mkd = matcl :: mkgen :: details;
    };

    namespace mk    = matcl :: mkgen;
    namespace mkd   = matcl :: mkgen :: details;
}}

namespace matcl { namespace mkgen
{

//------------------------------------------------------------------------------
//                      scalar tags
//------------------------------------------------------------------------------
// concept of Tag used to create const_value_scalar
// Tag must be derived from scal_data_const_value_tag<Tag>; 
// see scal_data_const_value_tag for other requirements
template<class Tag>
concept Tag_scalar_const_value = std::is_base_of<mk::scal_data_const_value_tag<Tag>,
                                    Tag>::value
                                && scal_data_const_value_tag_check :: template is_valid<Tag>;

// concept of Tag used to create value_scalar
// Tag must be derived from scal_data_value_tag<Tag>; see scal_data_value_tag for
// other requirements
template<class Tag>
concept Tag_scalar_value    = std::is_base_of<mk::scal_data_value_tag<Tag>, Tag>::value
                            && scal_data_value_tag_check :: template is_valid<Tag>;

// concept of Tag used to create gen_scalar
// Tag must be derived from scal_data_gen_value_tag<Tag>; see scal_data_gen_value_tag
// for other requirements
template<class Tag>
concept Tag_scalar_gen_value = std::is_base_of<mk::scal_data_gen_value_tag<Tag>, 
                                    Tag>::value
                            && mk::scal_data_gen_value_tag_check :: template is_valid<Tag>;

//------------------------------------------------------------------------------
//                      matrix tags
//------------------------------------------------------------------------------

// concept of Tag used to create const_value_mat
// Tag must be derived from matrix_data_const_value_tag<Tag>; 
// see matrix_data_const_value_tag for other requirements
template<class Tag>
concept Tag_matrix_const_data = std::is_base_of<mk::matrix_data_const_value_tag
                                    <Tag>, Tag>::value
                              && matrix_data_const_value_tag_check :: template is_valid<Tag>;

// concept of Tag used to create value_mat
// Tag must be derived from matrix_data_value_tag<Tag>; see matrix_data_value_tag
// for other requirements
template<class Tag>
concept Tag_matrix_data = std::is_base_of<mk::matrix_data_value_tag<Tag>, Tag>::value
                        && matrix_data_value_tag_check :: template is_valid<Tag>;

//------------------------------------------------------------------------------
//                      other tags
//------------------------------------------------------------------------------

//TODO
template<class Tag>
concept Tag_comp    = true;

//------------------------------------------------------------------------------
//                      matrix arrays
//------------------------------------------------------------------------------
// concept of matrix array
// Mat_array must be derived from matrix_array<Mat_array>; see matrix_array for
// other requirements
template<class Arr>
concept Mat_array   = std::is_base_of<mkd::matrix_array<Arr>, Arr>::value
                    && mkd::matrix_array_check::template is_valid<Arr>;

//------------------------------------------------------------------------------
//                      scalar arrays
//------------------------------------------------------------------------------
// concept of scalar data
// Scal_data must be derived from scalar_data<Scal_data>; see scalar_data for
// other requirements
template<class Arr>
concept Scal_data   = std::is_base_of<mkd::scalar_data<Arr>, Arr>::value
                    && mkd::scalar_data_check::template is_valid<Arr>;


//------------------------------------------------------------------------------
//                      other
//------------------------------------------------------------------------------
// concept of dependencies
// DPS must have dps<...> type
template<class Deps>
concept DPS   = mkd::dps_check<Deps>::value;

}}

