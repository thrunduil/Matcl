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

#include "matcl-scalar/objects/object_functions.h"
#include "matcl-dynamic/register_function.h"
#include "matcl-dynamic/object_type.h"
#include "matcl-matrep/objects/omatrix.h"
#include "matcl-matrep/lib_functions/func_unary.h"
#include "matcl-matrep/lib_functions/func_binary.h"
#include "matcl-matrep/lib_functions/func_matrix.h"

#include "matcl-dynamic/details/register_function_macro.h"

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

REGISTER_BIN_FUNC(f_mmul, OMatrix, operator*, mdyf::op_mul)

};
