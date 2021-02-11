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

namespace matcl { namespace dynamic
{

struct operator_eeq : register_function<operator_eeq, functions::op_eeq>
{
	static OBool eval(const ONull&, const ONull&)
	{
		return OBool(true);
	};
};

struct operator_neq : register_function<operator_neq, functions::op_neq>
{
	static OBool eval(const ONull&, const ONull&)
	{
		return OBool(false);
	};
};

struct operator_lt : register_function<operator_lt, functions::op_lt>
{
	static OBool eval(const ONull&, const ONull&)
	{
		return OBool(false);
	};
};

struct operator_leq : register_function<operator_leq, functions::op_leq>
{
	static OBool eval(const ONull&, const ONull&)
	{
		return OBool(true);
	};
};

struct operator_geq : register_function<operator_geq, functions::op_geq>
{
	static OBool eval(const ONull&, const ONull&)
	{
		return OBool(true);
	};
};

struct operator_gt : register_function<operator_gt, functions::op_gt>
{
	static OBool eval(const ONull&, const ONull&)
	{
		return OBool(false);
	};
};

struct op_bool : register_function<op_bool, functions::op_bool>
{
	static OBool eval(const ONull&)
	{
		return OBool(false);
	};
};

struct op_not : register_function<op_not, functions::op_not>
{
	static OBool eval(const ONull&)
	{
		return OBool(true);
	};
};

};};