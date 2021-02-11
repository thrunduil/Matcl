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

#include "matcl-dynamic/register_function.h"
#include "matcl-dynamic/object_type.h"
#include "matcl-dynamic/predefined_functions.h"
#include "matcl-dynamic/details/register_function_macro.h"

namespace matcl { namespace dynamic
{

MATCL_REGISTER_OPERATOR(eeq1, ==, bool, bool, functions::op_eeq)
MATCL_REGISTER_OPERATOR(neq1, !=, bool, bool, functions::op_neq)
MATCL_REGISTER_OPERATOR(lt1, <, bool, bool, functions::op_lt)
MATCL_REGISTER_OPERATOR(leq1, <=, bool, bool, functions::op_leq)
MATCL_REGISTER_OPERATOR(geq1, >=, bool, bool, functions::op_geq)
MATCL_REGISTER_OPERATOR(gt1, >, bool, bool, functions::op_gt)

MATCL_REGISTER_SCALAR_FUNC(op_not, !, bool, functions::op_not)
MATCL_REGISTER_SCALAR_FUNC(op_bool, (bool), bool, functions::op_bool)

};};