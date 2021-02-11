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

#include "matcl-dynamic/matcl_dynamic.h"
#include "matcl-scalar/matcl_scalar.h"

using namespace matcl;
using namespace matcl::dynamic;

// define function exp(x - y) for real values
static Real expm(Real x, Real y)
{
    return exp(x - y);
};

// define function exp(x - y) for complex values
static Complex expm(const Complex& x, const Complex& y)
{
    return exp(x - y);
};

// define function exp(x - y) for objects
object expm(const object& x, const object& y)
{
    static function_name f("expm");

    object ret;
    eval_function::eval(f, ret, x, y);
    return ret;
};

// register functions expm defined for real and complex values
MATCL_DEFINE_FUNCTION_NAME(expm)
MATCL_REGISTER_BIN_FUNC(tag1, expm, Real, Real, MATCL_FUNCTION_NAME(expm))
MATCL_REGISTER_BIN_FUNC(tag2, expm, Complex, Complex, MATCL_FUNCTION_NAME(expm))

void example_object()
{
    object x = object(2.0f);
    object y = object(Integer(3));
    object z = object(Float_complex(1.0f, 2.0f));

    object ret1 = expm(x, y);
    object ret2 = expm(y, z);

    disp(ret1);
    disp(ret2);
};