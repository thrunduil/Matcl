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

#pragma once

#include "matcl-core/details/result_of_impl.h"

// definition of templates that gives result type of a function
//     func_name(arg1, ...)
// called with args arg1,... with types Arg1, ... in the namespace matcl
// this template is defined in the namespace matcl::result_of:
//
// template <Arg1 ...>
// struct result_of_(function name)
// {
//     // return type of function call
//     using type          = ( .. implementation .. ); 
//
//     // typed object storing values of type type
//     using type_object   = object_type<type>;
// }
//
// template <Arg1 ...>
// struct has_(function name)
// {
//     //value is true if function call can be performed
//     static const bool value = ( .. implementation .. );
// }
//
// if the function (function name) with arguments of given type cannot be
// found by unqualified call in the nemaspace matcl then 
// result_of_(function name)<Args...>::type is not defined and
// has_(function name)<Args...>::value is false

// create result_of_[fun] and has_[fun] templates

MATCL_RESULT_OF_UNARY(real)
MATCL_RESULT_OF_UNARY(imag)
MATCL_RESULT_OF_UNARY(arg)
MATCL_RESULT_OF_UNARY(conj)
MATCL_RESULT_OF_UNARY(abs2)
MATCL_RESULT_OF_UNARY(abs)
MATCL_RESULT_OF_UNARY(nextabove)
MATCL_RESULT_OF_UNARY(nextbelow)
MATCL_RESULT_OF_UNARY(signbit)
MATCL_RESULT_OF_UNARY(sign)
MATCL_RESULT_OF_UNARY(isign)
MATCL_RESULT_OF_UNARY(eps)
MATCL_RESULT_OF_UNARY(sqrt)
MATCL_RESULT_OF_UNARY(sqrt_c)
MATCL_RESULT_OF_UNARY(sqrt1pm1)
MATCL_RESULT_OF_UNARY(sqrt1pm1_c)
MATCL_RESULT_OF_UNARY(cbrt)
MATCL_RESULT_OF_UNARY(exp)
MATCL_RESULT_OF_UNARY(expm1)
MATCL_RESULT_OF_UNARY(expi)
MATCL_RESULT_OF_UNARY(exp2)
MATCL_RESULT_OF_UNARY(exp10)
MATCL_RESULT_OF_UNARY(log)
MATCL_RESULT_OF_UNARY(log_c)
MATCL_RESULT_OF_UNARY(log1p)
MATCL_RESULT_OF_UNARY(log1p_c)
MATCL_RESULT_OF_UNARY(log2)
MATCL_RESULT_OF_UNARY(log2_c)
MATCL_RESULT_OF_UNARY(log10)
MATCL_RESULT_OF_UNARY(log10_c)
MATCL_RESULT_OF_UNARY(sin)
MATCL_RESULT_OF_UNARY(cos)
MATCL_RESULT_OF_UNARY(tan)
MATCL_RESULT_OF_UNARY(cot)
MATCL_RESULT_OF_UNARY(sec)
MATCL_RESULT_OF_UNARY(csc)
MATCL_RESULT_OF_UNARY(sinh)
MATCL_RESULT_OF_UNARY(cosh)
MATCL_RESULT_OF_UNARY(tanh)
MATCL_RESULT_OF_UNARY(coth)
MATCL_RESULT_OF_UNARY(sech)
MATCL_RESULT_OF_UNARY(csch)
MATCL_RESULT_OF_UNARY(asin)
MATCL_RESULT_OF_UNARY(asin_c)
MATCL_RESULT_OF_UNARY(acos)
MATCL_RESULT_OF_UNARY(acos_c)
MATCL_RESULT_OF_UNARY(atan)
MATCL_RESULT_OF_UNARY(acot)
MATCL_RESULT_OF_UNARY(asec)
MATCL_RESULT_OF_UNARY(asec_c)
MATCL_RESULT_OF_UNARY(acsc)
MATCL_RESULT_OF_UNARY(acsc_c)
MATCL_RESULT_OF_UNARY(asinh)
MATCL_RESULT_OF_UNARY(acosh)
MATCL_RESULT_OF_UNARY(acosh_c)
MATCL_RESULT_OF_UNARY(atanh)
MATCL_RESULT_OF_UNARY(atanh_c)
MATCL_RESULT_OF_UNARY(acoth)
MATCL_RESULT_OF_UNARY(acoth_c)
MATCL_RESULT_OF_UNARY(asech)
MATCL_RESULT_OF_UNARY(asech_c)
MATCL_RESULT_OF_UNARY(acsch)
MATCL_RESULT_OF_UNARY(inv)
MATCL_RESULT_OF_UNARY(invs)
MATCL_RESULT_OF_UNARY(logb)
MATCL_RESULT_OF_UNARY(ilogb)
MATCL_RESULT_OF_UNARY(floor)
MATCL_RESULT_OF_UNARY(ceil)
MATCL_RESULT_OF_UNARY(trunc)
MATCL_RESULT_OF_UNARY(round)
MATCL_RESULT_OF_UNARY(ifloor)
MATCL_RESULT_OF_UNARY(iceil)
MATCL_RESULT_OF_UNARY(itrunc)
MATCL_RESULT_OF_UNARY(iround)
MATCL_RESULT_OF_UNARY(cast_bool)

MATCL_RESULT_OF_UNARY(is_nan)
MATCL_RESULT_OF_UNARY(is_inf)
MATCL_RESULT_OF_UNARY(is_finite)
MATCL_RESULT_OF_UNARY(is_regular)
MATCL_RESULT_OF_UNARY(is_normal)
MATCL_RESULT_OF_UNARY(is_int)
MATCL_RESULT_OF_UNARY(is_real)
MATCL_RESULT_OF_UNARY(fpclassify)
MATCL_RESULT_OF_UNARY(is_true)
MATCL_RESULT_OF_UNARY(neg)

MATCL_RESULT_OF_BINARY(atan2)
MATCL_RESULT_OF_BINARY(hypot)
MATCL_RESULT_OF_BINARY(mod)
MATCL_RESULT_OF_BINARY(rem)
MATCL_RESULT_OF_BINARY(copysign)
MATCL_RESULT_OF_BINARY(nextafter)
MATCL_RESULT_OF_BINARY(ldexp)
MATCL_RESULT_OF_BINARY(scalbn)
MATCL_RESULT_OF_BINARY(frexp)
MATCL_RESULT_OF_BINARY(modf)
MATCL_RESULT_OF_BINARY(fdim)
MATCL_RESULT_OF_BINARY(pow)
MATCL_RESULT_OF_BINARY(pow_c)
MATCL_RESULT_OF_BINARY(powm1)
MATCL_RESULT_OF_BINARY(idiv)
MATCL_RESULT_OF_BINARY(op_xor)
MATCL_RESULT_OF_BINARY(elem_and)
MATCL_RESULT_OF_BINARY(elem_or)
MATCL_RESULT_OF_BINARY(elem_xor)
MATCL_RESULT_OF_BINARY(min)
MATCL_RESULT_OF_BINARY(max)
MATCL_RESULT_OF_BINARY(mul)
MATCL_RESULT_OF_BINARY(rising_factorial)
MATCL_RESULT_OF_BINARY(falling_factorial)
MATCL_RESULT_OF_BINARY(eeq_nan)
MATCL_RESULT_OF_BINARY(neq_nan)
MATCL_RESULT_OF_BINARY(div_0)
MATCL_RESULT_OF_BINARY(div_1)

MATCL_RESULT_OF_TEMPL_1_1(bernoulli_b2n)
MATCL_RESULT_OF_TEMPL_1_1(factorial)
MATCL_RESULT_OF_TEMPL_1_1(double_factorial)
MATCL_RESULT_OF_TEMPL_1_0(max_bernoulli_b2n)
MATCL_RESULT_OF_TEMPL_1_2(binomial_coefficient)

MATCL_RESULT_OF_UNARY_OP(uminus,-)
MATCL_RESULT_OF_UNARY_OP(bool, (bool))
MATCL_RESULT_OF_UNARY_OP(op_not,!)

MATCL_RESULT_OF_BINARY_OP(eeq,==)
MATCL_RESULT_OF_BINARY_OP(neq,!=)
MATCL_RESULT_OF_BINARY_OP(geq,>=)
MATCL_RESULT_OF_BINARY_OP(leq,<=)
MATCL_RESULT_OF_BINARY_OP(gt,>)
MATCL_RESULT_OF_BINARY_OP(lt,<)
MATCL_RESULT_OF_BINARY_OP(op_and,&&)
MATCL_RESULT_OF_BINARY_OP(op_or,||)

MATCL_RESULT_OF_BINARY_OP(plus,+)
MATCL_RESULT_OF_BINARY_FUN_ADD_INT(plus)

MATCL_RESULT_OF_BINARY_OP(minus,-)
MATCL_RESULT_OF_BINARY_FUN_ADD_INT(minus)

MATCL_RESULT_OF_BINARY_OP(op_mul,*)
MATCL_RESULT_OF_BINARY_FUN_ADD_INT(op_mul)

MATCL_RESULT_OF_BINARY_OP(div,/)
MATCL_RESULT_OF_BINARY_FUN_ADD_INT2(div)
