/*
 *  This file is a part of Matrix Computation Library (MATCL)
 *
 *  Copyright (c) Pawe³ Kowal 2017 - 2018
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

#include "matcl-scalar/objects/object_functions.h"
#include "matcl-dynamic/register_function.h"
#include "matcl-dynamic/object_type.h"
#include "matcl-matrep/objects/omatrix.h"
#include "matcl-matrep/lib_functions/func_unary.h"
#include "matcl-matrep/lib_functions/func_binary.h"
#include "matcl-matrep/IO/matrix_io.h"
#include "matcl-dynamic/details/register_function_macro.h"
#include "matcl-core/details/IO/printer.h"

namespace matcl 
{

namespace mdd = matcl::dynamic::details;
namespace mdyf = matcl::dynamic::functions;

//---------------------------------------------------------------
//              bin functions
//---------------------------------------------------------------
#define REGISTER_BIN_FUNC(fname, ret, func, func_name)                              \
MATCL_REGISTER_BIN_FUNC(fname##1, func, Matrix, Matrix, func_name)                  \
                                                                                    \
struct fname##2 : dynamic::register_function_template_return<fname##2, func_name>   \
{                                                                                   \
	static ret eval(const OMatrix& obj1, const dynamic::Template& obj2)             \
	{                                                                               \
		return ret(func(obj1.get(), obj2.get()));                                   \
	};                                                                              \
    static dynamic::Type eval_return(int n_template, const dynamic::Type* templ,    \
                int n_arg, const dynamic::Type* args)                               \
    {                                                                               \
        (void)n_template; (void)templ; (void)n_arg;                                 \
        if (args[1] == dynamic::Any::get_static_type())                             \
            return dynamic::Type();                                                 \
        else                                                                        \
            return OMatrix::get_static_type();                                      \
    };                                                                              \
};                                                                                  \
struct fname##3 : dynamic::register_function_template_return<fname##3, func_name>   \
{                                                                                   \
	static ret eval(const dynamic::Template& obj1, const OMatrix& obj2)             \
	{                                                                               \
		return ret(func(obj1.get(), obj2.get()));                                   \
	};                                                                              \
    static dynamic::Type eval_return(int n_template, const dynamic::Type* templ,    \
                int n_arg, const dynamic::Type* args)                               \
    {                                                                               \
        (void)n_template; (void)templ; (void)n_arg;                                 \
        if (args[0] == dynamic::Any::get_static_type())                             \
            return dynamic::Type();                                                 \
        else                                                                        \
            return OMatrix::get_static_type();                                      \
    };                                                                              \
};                                                                                  \

REGISTER_BIN_FUNC(f_plus, OMatrix, operator+, mdyf::op_plus)
REGISTER_BIN_FUNC(f_minus, OMatrix, operator-, mdyf::op_minus)
REGISTER_BIN_FUNC(f_mul, OMatrix, mul, mdyf::elem_mul)
REGISTER_BIN_FUNC(f_div, OMatrix, operator/, mdyf::op_div)
REGISTER_BIN_FUNC(f_idiv, OMatrix, idiv, mdyf::idiv)
REGISTER_BIN_FUNC(f_div_0, OMatrix, div_0, mdyf::div_0)
REGISTER_BIN_FUNC(f_div_1, OMatrix, div_1, mdyf::div_1)
REGISTER_BIN_FUNC(f_atan2, OMatrix, atan2, mdyf::atan2)
REGISTER_BIN_FUNC(f_hypot, OMatrix, hypot, mdyf::hypot)
REGISTER_BIN_FUNC(f_mod, OMatrix, mod, mdyf::mod)
REGISTER_BIN_FUNC(f_rem, OMatrix, rem, mdyf::rem)
REGISTER_BIN_FUNC(f_pow, OMatrix, pow, mdyf::pow)
REGISTER_BIN_FUNC(f_pow_c, OMatrix, pow_c, mdyf::pow_c)
REGISTER_BIN_FUNC(f_elem_and, OMatrix, elem_and, mdyf::elem_and)
REGISTER_BIN_FUNC(f_elem_or, OMatrix, elem_or, mdyf::elem_or)
REGISTER_BIN_FUNC(f_elem_xor, OMatrix, elem_xor, mdyf::elem_xor)
REGISTER_BIN_FUNC(f_eeq, OMatrix, eeq, mdyf::op_eeq)
REGISTER_BIN_FUNC(f_eeq_nan, OMatrix, eeq_nan, mdyf::eeq_nan)
REGISTER_BIN_FUNC(f_neq, OMatrix, neq, mdyf::op_neq)
REGISTER_BIN_FUNC(f_neq_nan, OMatrix, neq_nan, mdyf::neq_nan)
REGISTER_BIN_FUNC(f_leq, OMatrix, leq, mdyf::op_leq)
REGISTER_BIN_FUNC(f_geq, OMatrix, geq, mdyf::op_geq)
REGISTER_BIN_FUNC(f_gt, OMatrix, gt, mdyf::op_gt)
REGISTER_BIN_FUNC(f_lt, OMatrix, lt, mdyf::op_lt)
REGISTER_BIN_FUNC(f_max, OMatrix, max, mdyf::max)
REGISTER_BIN_FUNC(f_min, OMatrix, min, mdyf::min)
REGISTER_BIN_FUNC(f_op_and, dynamic::OBool, op_and, mdyf::op_and)
REGISTER_BIN_FUNC(f_op_or, dynamic::OBool, op_or, mdyf::op_or)
REGISTER_BIN_FUNC(f_op_xor, dynamic::OBool, op_xor, mdyf::op_xor)

//not defined for a matrix:
//nextabove, nextbelow, signbit, fpclassify, ldexp, frexp, modf, inv, powm1, fdim
//copysign, nextafter, fma, dot2_ac

};
