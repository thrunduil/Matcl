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

#include "matcl-dynamic/function_name.h"

namespace matcl { namespace dynamic { namespace special_functions
{

const function_name&   convert_numeric_promotion_1();
const function_name&   convert_numeric_promotion_2();
const function_name&   convert_numeric_equivalent_1();
const function_name&   convert_numeric_equivalent_2();
const function_name&   convert_numeric_decay_1();
const function_name&   convert_numeric_decay_2();
const function_name&   convert_numeric_int_float();
const function_name&   convert_numeric_explicit();
const function_name&   convert_numeric_cast();
const function_name&   convert_any();
const function_name&   convert_unit();
const function_name&   convert_id();

const function_name&   convert_promotion();
const function_name&   convert_equivalent();
const function_name&   convert_decay();
const function_name&   convert_explicit();
const function_name&   convert_cast();

const function_name&   unifier();
const function_name&   assign();
const function_name&   assign_id();

bool is_special(const function_name& func);

};};};

