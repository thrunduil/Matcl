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

#include "mkgen/mkgen_config.h"
#include "matcl-core/matrix/scalar_types.h"

namespace matcl { namespace mkgen 
{
//------------------------------------------------------------------------
//                              arrays
//------------------------------------------------------------------------
//technical arrays
struct call_array_type;

//------------------------------------------------------------------------
//                              other
//------------------------------------------------------------------------
template<class Tag, class Mat, class Assignments = list::list<>>
struct computation;

template<Integer Pos> struct colon;
struct colon_all;

template<Integer Start, Integer End> 
struct colon2;

template<Integer Start, Integer Step, Integer End> 
struct colon3;

template<class Mat1, class Mat2>
struct mat_assign;

template<class Subject, class Scalar, class Colon>
struct comp_assign_1;

template<class Comp>
struct make_comp_result;

template<class Mat_L, class Tag, bool Force>
struct mat_temporary;

template<class T, bool With_Forced>
struct is_temporary_mat;

template<class Temp_Tag, class Ret_Tag, class Colon, Integer Rows, Integer Cols, bool Init>
struct dps_modif;

template<class Tag, class Colon, bool Init>
struct modif;

template<class Tag, class Colon, Integer Mat_Rows, Integer Mat_Cols>
struct modif2;

template <class Tag, Integer Mat_Rows, Integer Mat_Cols, Integer Row, Integer Col>
struct get_temporary;

template<class Val, class Elem>
struct loop_context_data_scalar;

}}
