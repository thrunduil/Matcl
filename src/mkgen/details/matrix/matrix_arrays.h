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

#include "mkgen/matrix/matrix.h"

namespace matcl { namespace mkgen { namespace details
{

template<class Array_t, Integer Offset, Integer Step>
struct sub_array_1 {};

template<class Array_t, Integer Offset1, Integer Offset2, Integer Step1, Integer step2>
struct sub_array_2 {};

// const_mat array
template<class Tag>                                     
struct const_array{};

template<class Tag>                                     
struct gen_array{};

template<class Tag> 
struct output_array{};

template<class Tag, Integer Rows, Integer Cols>
struct temp_output_array{};

template<class Tag, class... Assign_List>
struct virtual_array {};

//TODO
/*
// Symbolic Array of generic elements of type element<Tag, row, col, ...> for all 
// pairs of (row, col) available in given matrix.
template<class Tag>
struct array {};

template<class Tag>
struct output_array {};
*/

}}}
