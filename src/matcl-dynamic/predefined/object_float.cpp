/*
 *  This file is a part of Matrix Computation Library (MATCL)
 *
 *  Copyright (c) Pawe³ Kowal 2017
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

#pragma warning(push)
#pragma warning(disable: 4800) //forcing value to bool 'true' or 'false' (performance warning)

namespace matcl
{

static Float real(const Float& v)                      { return v; };
static Float imag(const Float&)                        { return 0.0f; };
static Float idiv(const Float& v1, const Float& v2)    { return v1 / v2; };

};

namespace matcl { namespace dynamic
{

MATCL_REGISTER_OPERATOR(eeq1, ==, Float, Float, functions::op_eeq)
MATCL_REGISTER_OPERATOR(neq1, !=, Float, Float, functions::op_neq)
MATCL_REGISTER_OPERATOR(lt1, <, Float, Float, functions::op_lt)
MATCL_REGISTER_OPERATOR(leq1, <=, Float, Float, functions::op_leq)
MATCL_REGISTER_OPERATOR(geq1, >=, Float, Float, functions::op_geq)
MATCL_REGISTER_OPERATOR(gt1, >, Float, Float, functions::op_gt)

MATCL_REGISTER_OPERATOR(plus1, +, Float, Float, functions::op_plus)
MATCL_REGISTER_OPERATOR(minus1, -, Float, Float, functions::op_minus)
MATCL_REGISTER_OPERATOR(mult1, *, Float, Float, functions::op_mul)
MATCL_REGISTER_OPERATOR(div1, /, Float, Float, functions::op_div)
MATCL_REGISTER_BIN_FUNC(idiv1, matcl::idiv, Float, Float, functions::idiv)

#define MATCL_REGISTER_SCALAR_FUNC_SIMP(fn)                             \
MATCL_REGISTER_SCALAR_FUNC(f_##fn, matcl::fn, Float, functions::fn)

MATCL_REGISTER_SCALAR_FUNC_SIMP(real)
MATCL_REGISTER_SCALAR_FUNC_SIMP(imag)

MATCL_REGISTER_SCALAR_FUNC(uminus, -, Float, functions::op_uminus)
MATCL_REGISTER_SCALAR_FUNC(op_not, !, Float, functions::op_not)
MATCL_REGISTER_SCALAR_FUNC(op_bool, (bool), Float, functions::op_bool)

};};

#pragma warning(pop)