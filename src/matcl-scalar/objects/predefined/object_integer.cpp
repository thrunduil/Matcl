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

#define MATCL_REGISTER_SCALAR_FUNC_SIMP(fn)                             \
MATCL_REGISTER_SCALAR_FUNC(f_##fn, matcl::fn, Integer, mdyf::fn)

MATCL_REGISTER_SCALAR_FUNC_SIMP(inv)
MATCL_REGISTER_SCALAR_FUNC_SIMP(invs)
MATCL_REGISTER_SCALAR_FUNC_SIMP(is_nan)
MATCL_REGISTER_SCALAR_FUNC_SIMP(is_inf)
MATCL_REGISTER_SCALAR_FUNC_SIMP(is_regular)
MATCL_REGISTER_SCALAR_FUNC_SIMP(is_normal)
MATCL_REGISTER_SCALAR_FUNC_SIMP(is_finite)
MATCL_REGISTER_SCALAR_FUNC_SIMP(is_int)
MATCL_REGISTER_SCALAR_FUNC_SIMP(is_real)
MATCL_REGISTER_SCALAR_FUNC_SIMP(fpclassify)
MATCL_REGISTER_SCALAR_FUNC_SIMP(abs)
MATCL_REGISTER_SCALAR_FUNC_SIMP(abs2)
MATCL_REGISTER_SCALAR_FUNC_SIMP(conj)
MATCL_REGISTER_SCALAR_FUNC_SIMP(arg)
MATCL_REGISTER_SCALAR_FUNC_SIMP(signbit)
MATCL_REGISTER_SCALAR_FUNC_SIMP(sign)
MATCL_REGISTER_SCALAR_FUNC_SIMP(isign)
MATCL_REGISTER_SCALAR_FUNC_SIMP(nextabove)
MATCL_REGISTER_SCALAR_FUNC_SIMP(nextbelow)
MATCL_REGISTER_SCALAR_FUNC_SIMP(eps)
MATCL_REGISTER_SCALAR_FUNC_SIMP(sqrt)
MATCL_REGISTER_SCALAR_FUNC_SIMP(sqrt_c)
MATCL_REGISTER_SCALAR_FUNC_SIMP(cbrt)
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
MATCL_REGISTER_SCALAR_FUNC_SIMP(ifloor)
MATCL_REGISTER_SCALAR_FUNC_SIMP(iceil)
MATCL_REGISTER_SCALAR_FUNC_SIMP(iround)
MATCL_REGISTER_SCALAR_FUNC_SIMP(itrunc)
MATCL_REGISTER_SCALAR_FUNC_SIMP(sqrt1pm1)

MATCL_REGISTER_BIN_FUNC(ldexp1, matcl::ldexp, Integer, Integer, mdyf::ldexp)
MATCL_REGISTER_BIN_FUNC(frexp1, matcl::frexp, Integer, Integer&, mdyf::frexp)
MATCL_REGISTER_BIN_FUNC(modf1, matcl::modf, Integer, Real&, mdyf::modf)

MATCL_REGISTER_BIN_FUNC(f_atan2, matcl::atan2, Integer, Integer, mdyf::atan2)
MATCL_REGISTER_BIN_FUNC(f_hypot, matcl::hypot, Integer, Integer, mdyf::hypot)
MATCL_REGISTER_BIN_FUNC(f_mod, matcl::mod, Integer, Integer, mdyf::mod)
MATCL_REGISTER_BIN_FUNC(f_rem, matcl::rem, Integer, Integer, mdyf::rem)
MATCL_REGISTER_BIN_FUNC(f_dim, matcl::fdim, Integer, Integer, mdyf::fdim)

MATCL_REGISTER_BIN_FUNC(f_copysign, matcl::copysign, Integer, Integer, mdyf::copysign)
MATCL_REGISTER_BIN_FUNC(f_nextafter, matcl::nextafter, Integer, Integer, mdyf::nextafter)
MATCL_REGISTER_BIN_FUNC(f_float_distance, matcl::float_distance, Integer, Integer, mdyf::float_distance)

MATCL_REGISTER_BIN_FUNC(f_pow, matcl::pow, Integer, Integer, mdyf::pow)
MATCL_REGISTER_BIN_FUNC(f_pow_c, matcl::pow_c, Integer, Integer, mdyf::pow_c)
MATCL_REGISTER_BIN_FUNC(f_powm1, matcl::powm1, Integer, Integer, mdyf::powm1)

MATCL_REGISTER_BIN_FUNC(eeq1, eeq_nan, Integer, Integer, mdyf::eeq_nan)
MATCL_REGISTER_BIN_FUNC(neq1, neq_nan, Integer, Integer, mdyf::neq_nan)

MATCL_REGISTER_BIN_FUNC(f_rising_factorial, matcl::rising_factorial, Integer, Integer, mdyf::rising_factorial)
MATCL_REGISTER_BIN_FUNC(f_falling_factorial, matcl::falling_factorial, Integer, Integer, mdyf::falling_factorial)

MATCL_REGISTER_TEMPL_1_FUNC(f_bernoulli_b2n_i,bernoulli_b2n,Integer,Integer,mdyf::bernoulli_b2n)
MATCL_REGISTER_TEMPL_0_FUNC(f_max_bernoulli_b2n_i,max_bernoulli_b2n,Integer,mdyf::max_bernoulli_b2n)
MATCL_REGISTER_TEMPL_1_FUNC(f_factorial_i,factorial,Integer,Integer,mdyf::factorial)
MATCL_REGISTER_TEMPL_1_FUNC(f_double_factorial_i,double_factorial,Integer,Integer,mdyf::double_factorial)
MATCL_REGISTER_TEMPL_2_FUNC(f_binomial_coefficient_i,binomial_coefficient,Integer,Integer,Integer,mdyf::binomial_coefficient)

};
