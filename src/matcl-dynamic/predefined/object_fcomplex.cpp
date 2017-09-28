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
#include "matcl-core/details/complex_details.h"
#include "matcl-dynamic/details/register_function_macro.h"

namespace matcl
{

static Float real(const Float_complex& v)  { return real(v.value); };
static Float imag(const Float_complex& v)  { return imag(v.value); };

static Float_complex idiv(const Float_complex& v1, const Float_complex& v2)
{ 
    return v1 / v2; 
};

static bool operator!(const Complex& x)
{
    return real(x) == 0.0f && imag(x) == 0.0f;
};

template<class T1, class T2>
static bool operator==(const T1& a, const T2& b)
{
    return details::eeq_c(a, b);
};

template<class T1, class T2>
static bool operator!=(const T1& a, const T2& b)
{
    return details::neq_c(a, b);
};

template<class T1, class T2>
static bool operator<(const T1& a, const T2& b)
{
    return details::lt_c(a, b);
};

template<class T1, class T2>
static bool operator>(const T1& a, const T2& b)
{
    return details::gt_c(a, b);
};

template<class T1, class T2>
static bool operator<=(const T1& a, const T2& b)
{
    return details::leq_c(a, b);
};

template<class T1, class T2>
static bool operator>=(const T1& a, const T2& b)
{
    return details::geq_c(a, b);
};

};

namespace matcl { namespace dynamic
{

MATCL_REGISTER_BIN_FUNC(eeq1, operator==, Float_complex, Float_complex, functions::op_eeq)
MATCL_REGISTER_BIN_FUNC(eeq2, operator==, Float, Float_complex, functions::op_eeq)
MATCL_REGISTER_BIN_FUNC(eeq3, operator==, Float_complex, Float, functions::op_eeq)

MATCL_REGISTER_BIN_FUNC(neq1, operator!=, Float_complex, Float_complex, functions::op_neq)
MATCL_REGISTER_BIN_FUNC(neq2, operator!=, Float_complex, Float, functions::op_neq)
MATCL_REGISTER_BIN_FUNC(neq3, operator!=, Float, Float_complex, functions::op_neq)

MATCL_REGISTER_BIN_FUNC(lt1, operator<, Float_complex, Float_complex, functions::op_lt)
MATCL_REGISTER_BIN_FUNC(lt2, operator<, Float_complex, Float, functions::op_lt)
MATCL_REGISTER_BIN_FUNC(lt3, operator<, Float, Float_complex, functions::op_lt)

MATCL_REGISTER_BIN_FUNC(leq1, operator<=, Float_complex, Float_complex, functions::op_leq)
MATCL_REGISTER_BIN_FUNC(leq2, operator<=, Float_complex, Float, functions::op_leq)
MATCL_REGISTER_BIN_FUNC(leq3, operator<=, Float, Float_complex, functions::op_leq)

MATCL_REGISTER_BIN_FUNC(geq1, operator>=, Float_complex, Float_complex, functions::op_geq)
MATCL_REGISTER_BIN_FUNC(geq2, operator>=, Float, Float_complex, functions::op_geq)
MATCL_REGISTER_BIN_FUNC(geq3, operator>=, Float_complex, Float, functions::op_geq)

MATCL_REGISTER_BIN_FUNC(gt1, operator>, Float_complex, Float_complex, functions::op_gt)
MATCL_REGISTER_BIN_FUNC(gt2, operator>, Float_complex, Float, functions::op_gt)
MATCL_REGISTER_BIN_FUNC(gt3, operator>, Float, Float_complex, functions::op_gt)

MATCL_REGISTER_BIN_FUNC(plus1, operator+, Float_complex, Float_complex, functions::op_plus)
MATCL_REGISTER_BIN_FUNC(plus2, operator+, Float_complex, Float, functions::op_plus)
MATCL_REGISTER_BIN_FUNC(plus3, operator+, Float, Float_complex, functions::op_plus)

MATCL_REGISTER_BIN_FUNC(minus1, operator-, Float_complex, Float_complex, functions::op_minus)
MATCL_REGISTER_BIN_FUNC(minus2, operator-, Float_complex, Float, functions::op_minus)
MATCL_REGISTER_BIN_FUNC(minus3, operator-, Float, Float_complex, functions::op_minus)

MATCL_REGISTER_BIN_FUNC(mult1, operator*, Float_complex, Float_complex, functions::op_mul)
MATCL_REGISTER_BIN_FUNC(mult2, operator*, Float_complex, Float, functions::op_mul)
MATCL_REGISTER_BIN_FUNC(mult3, operator*, Float, Float_complex, functions::op_mul)

MATCL_REGISTER_BIN_FUNC(div1, operator/, Float_complex, Float_complex, functions::op_div)
MATCL_REGISTER_BIN_FUNC(div2, operator/, Float_complex, Float, functions::op_div)
MATCL_REGISTER_BIN_FUNC(div3, operator/, Float, Float_complex, functions::op_div)

MATCL_REGISTER_BIN_FUNC(idiv1, idiv, Float_complex, Float_complex, functions::idiv)
MATCL_REGISTER_BIN_FUNC(idiv2, idiv, Float_complex, Float, functions::idiv)
MATCL_REGISTER_BIN_FUNC(idiv3, idiv, Float, Float_complex, functions::idiv)

#define MATCL_REGISTER_SCALAR_FUNC_SIMP(fn)                             \
MATCL_REGISTER_SCALAR_FUNC(f_##fn, matcl::fn, Float_complex, functions::fn)

MATCL_REGISTER_SCALAR_FUNC_SIMP(real)
MATCL_REGISTER_SCALAR_FUNC_SIMP(imag)

MATCL_REGISTER_SCALAR_FUNC(uminus, matcl::operator-, Float_complex, functions::op_uminus)
MATCL_REGISTER_SCALAR_FUNC(op_not, matcl::operator!, Float_complex, functions::op_not)
MATCL_REGISTER_SCALAR_FUNC(op_bool, (bool), Float_complex, functions::op_bool)

};};