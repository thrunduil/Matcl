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

#include "matcl-scalar/objects/object_functions.h"
#include "matcl-dynamic/register_function.h"
#include "matcl-dynamic/object_type.h"
#include "matcl-scalar/lib_functions/func_unary.h"
#include "matcl-scalar/lib_functions/func_binary.h"
#include "matcl-dynamic/details/register_function_macro.h"

#include "matcl-scalar/matcl_scalar.h"

namespace matcl 
{

namespace mdd = matcl::dynamic::details;
namespace mdyf = matcl::dynamic::functions;

MATCL_REGISTER_BIN_FUNC(ldexp1, ldexp, Complex, Integer, mdyf::ldexp)

#define MATCL_REGISTER_SCALAR_FUNC_SIMP(fn)                             \
MATCL_REGISTER_SCALAR_FUNC(f_##fn, matcl::fn, Complex, mdyf::fn)

MATCL_REGISTER_SCALAR_FUNC_SIMP(inv)
MATCL_REGISTER_SCALAR_FUNC_SIMP(invs)
MATCL_REGISTER_SCALAR_FUNC_SIMP(is_nan)
MATCL_REGISTER_SCALAR_FUNC_SIMP(is_inf)
MATCL_REGISTER_SCALAR_FUNC_SIMP(is_regular)
MATCL_REGISTER_SCALAR_FUNC_SIMP(is_normal)
MATCL_REGISTER_SCALAR_FUNC_SIMP(is_finite)
MATCL_REGISTER_SCALAR_FUNC_SIMP(is_int)
MATCL_REGISTER_SCALAR_FUNC_SIMP(is_real)
MATCL_REGISTER_SCALAR_FUNC_SIMP(abs)
MATCL_REGISTER_SCALAR_FUNC_SIMP(abs2)
MATCL_REGISTER_SCALAR_FUNC_SIMP(conj)
MATCL_REGISTER_SCALAR_FUNC_SIMP(arg)
MATCL_REGISTER_SCALAR_FUNC_SIMP(sign)
MATCL_REGISTER_SCALAR_FUNC_SIMP(eps)
MATCL_REGISTER_SCALAR_FUNC_SIMP(sqrt)
MATCL_REGISTER_SCALAR_FUNC_SIMP(sqrt_c)
MATCL_REGISTER_SCALAR_FUNC_SIMP(exp)
MATCL_REGISTER_SCALAR_FUNC_SIMP(expm1)
MATCL_REGISTER_SCALAR_FUNC_SIMP(expi)
MATCL_REGISTER_SCALAR_FUNC_SIMP(exp2)
MATCL_REGISTER_SCALAR_FUNC_SIMP(exp10)
MATCL_REGISTER_SCALAR_FUNC_SIMP(log)
MATCL_REGISTER_SCALAR_FUNC_SIMP(log_c)
MATCL_REGISTER_SCALAR_FUNC_SIMP(log1p)
MATCL_REGISTER_SCALAR_FUNC_SIMP(log1p_c)
MATCL_REGISTER_SCALAR_FUNC_SIMP(log2)
MATCL_REGISTER_SCALAR_FUNC_SIMP(log2_c)
MATCL_REGISTER_SCALAR_FUNC_SIMP(log10)
MATCL_REGISTER_SCALAR_FUNC_SIMP(log10_c)
MATCL_REGISTER_SCALAR_FUNC_SIMP(logb)
MATCL_REGISTER_SCALAR_FUNC_SIMP(ilogb)
MATCL_REGISTER_SCALAR_FUNC_SIMP(sin)
MATCL_REGISTER_SCALAR_FUNC_SIMP(cos)
MATCL_REGISTER_SCALAR_FUNC_SIMP(tan)
MATCL_REGISTER_SCALAR_FUNC_SIMP(cot)
MATCL_REGISTER_SCALAR_FUNC_SIMP(sec)
MATCL_REGISTER_SCALAR_FUNC_SIMP(csc)
MATCL_REGISTER_SCALAR_FUNC_SIMP(sinh)
MATCL_REGISTER_SCALAR_FUNC_SIMP(cosh)
MATCL_REGISTER_SCALAR_FUNC_SIMP(tanh)
MATCL_REGISTER_SCALAR_FUNC_SIMP(coth)
MATCL_REGISTER_SCALAR_FUNC_SIMP(sech)
MATCL_REGISTER_SCALAR_FUNC_SIMP(csch)
MATCL_REGISTER_SCALAR_FUNC_SIMP(asin)
MATCL_REGISTER_SCALAR_FUNC_SIMP(asin_c)
MATCL_REGISTER_SCALAR_FUNC_SIMP(acos)
MATCL_REGISTER_SCALAR_FUNC_SIMP(acos_c)
MATCL_REGISTER_SCALAR_FUNC_SIMP(atan)
MATCL_REGISTER_SCALAR_FUNC_SIMP(acot)
MATCL_REGISTER_SCALAR_FUNC_SIMP(asec)
MATCL_REGISTER_SCALAR_FUNC_SIMP(asec_c)
MATCL_REGISTER_SCALAR_FUNC_SIMP(acsc)
MATCL_REGISTER_SCALAR_FUNC_SIMP(acsc_c)
MATCL_REGISTER_SCALAR_FUNC_SIMP(asinh)
MATCL_REGISTER_SCALAR_FUNC_SIMP(acosh)
MATCL_REGISTER_SCALAR_FUNC_SIMP(acosh_c)
MATCL_REGISTER_SCALAR_FUNC_SIMP(atanh)
MATCL_REGISTER_SCALAR_FUNC_SIMP(atanh_c)
MATCL_REGISTER_SCALAR_FUNC_SIMP(acoth)
MATCL_REGISTER_SCALAR_FUNC_SIMP(acoth_c)
MATCL_REGISTER_SCALAR_FUNC_SIMP(asech)
MATCL_REGISTER_SCALAR_FUNC_SIMP(asech_c)
MATCL_REGISTER_SCALAR_FUNC_SIMP(acsch)
MATCL_REGISTER_SCALAR_FUNC_SIMP(floor)
MATCL_REGISTER_SCALAR_FUNC_SIMP(ceil)
MATCL_REGISTER_SCALAR_FUNC_SIMP(round)
MATCL_REGISTER_SCALAR_FUNC_SIMP(trunc)
MATCL_REGISTER_SCALAR_FUNC_SIMP(sqrt1pm1)

MATCL_REGISTER_BIN_FUNC(f_hypot_1, matcl::hypot, Complex, Complex, mdyf::hypot)
MATCL_REGISTER_BIN_FUNC(f_hypot_2, matcl::hypot, Complex, Real, mdyf::hypot)
MATCL_REGISTER_BIN_FUNC(f_hypot_3, matcl::hypot, Real, Complex, mdyf::hypot)

MATCL_REGISTER_BIN_FUNC(f_pow_1, matcl::pow, Complex, Complex, mdyf::pow)
MATCL_REGISTER_BIN_FUNC(f_pow_2, matcl::pow, Complex, Integer, mdyf::pow)
MATCL_REGISTER_BIN_FUNC(f_pow_3, matcl::pow, Complex, Real, mdyf::pow)
MATCL_REGISTER_BIN_FUNC(f_pow_4, matcl::pow, Integer, Complex, mdyf::pow)
MATCL_REGISTER_BIN_FUNC(f_pow_5, matcl::pow, Real, Complex, mdyf::pow)

MATCL_REGISTER_BIN_FUNC(f_pow_c_1, matcl::pow_c, Complex, Complex, mdyf::pow_c)
MATCL_REGISTER_BIN_FUNC(f_pow_c_2, matcl::pow_c, Complex, Integer, mdyf::pow_c)
MATCL_REGISTER_BIN_FUNC(f_pow_c_3, matcl::pow_c, Complex, Real, mdyf::pow_c)
MATCL_REGISTER_BIN_FUNC(f_pow_c_4, matcl::pow_c, Integer, Complex, mdyf::pow_c)
MATCL_REGISTER_BIN_FUNC(f_pow_c_5, matcl::pow_c, Real, Complex, mdyf::pow_c)

MATCL_REGISTER_BIN_FUNC(eeq1, eeq_nan, Complex, Complex, mdyf::eeq_nan)
MATCL_REGISTER_BIN_FUNC(eeq2, eeq_nan, Real, Complex, mdyf::eeq_nan)
MATCL_REGISTER_BIN_FUNC(eeq3, eeq_nan, Complex, Real, mdyf::eeq_nan)

MATCL_REGISTER_BIN_FUNC(neq1, neq_nan, Complex, Complex, mdyf::neq_nan)
MATCL_REGISTER_BIN_FUNC(neq2, neq_nan, Real, Complex, mdyf::neq_nan)
MATCL_REGISTER_BIN_FUNC(neq3, neq_nan, Complex, Real, mdyf::neq_nan)

MATCL_REGISTER_TEMPL_1_FUNC(f_bernoulli_b2n_c,bernoulli_b2n,Complex,Integer,mdyf::bernoulli_b2n)
MATCL_REGISTER_TEMPL_0_FUNC(f_max_bernoulli_b2n_c,max_bernoulli_b2n,Complex,mdyf::max_bernoulli_b2n)
MATCL_REGISTER_TEMPL_1_FUNC(f_factorial_c,factorial,Complex,Integer,mdyf::factorial)
MATCL_REGISTER_TEMPL_1_FUNC(f_double_factorial_c,double_factorial,Complex,Integer,mdyf::double_factorial)
MATCL_REGISTER_TEMPL_2_FUNC(f_binomial_coefficient_c,binomial_coefficient,Complex,Integer,Integer,mdyf::binomial_coefficient)

};
