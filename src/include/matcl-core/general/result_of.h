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

matcl_RESULT_OF_UNARY(real)
matcl_RESULT_OF_UNARY(imag)
matcl_RESULT_OF_UNARY(arg)
matcl_RESULT_OF_UNARY(conj)
matcl_RESULT_OF_UNARY(abs2)
matcl_RESULT_OF_UNARY(abs)
matcl_RESULT_OF_UNARY(nextabove)
matcl_RESULT_OF_UNARY(nextbelow)
matcl_RESULT_OF_UNARY(signbit)
matcl_RESULT_OF_UNARY(sign)
matcl_RESULT_OF_UNARY(isign)
matcl_RESULT_OF_UNARY(eps)
matcl_RESULT_OF_UNARY(sqrt)
matcl_RESULT_OF_UNARY(sqrt_c)
matcl_RESULT_OF_UNARY(sqrt1pm1)
matcl_RESULT_OF_UNARY(sqrt1pm1_c)
matcl_RESULT_OF_UNARY(cbrt)
matcl_RESULT_OF_UNARY(exp)
matcl_RESULT_OF_UNARY(expm1)
matcl_RESULT_OF_UNARY(expi)
matcl_RESULT_OF_UNARY(exp2)
matcl_RESULT_OF_UNARY(exp10)
matcl_RESULT_OF_UNARY(log)
matcl_RESULT_OF_UNARY(log_c)
matcl_RESULT_OF_UNARY(log1p)
matcl_RESULT_OF_UNARY(log1p_c)
matcl_RESULT_OF_UNARY(log2)
matcl_RESULT_OF_UNARY(log2_c)
matcl_RESULT_OF_UNARY(log10)
matcl_RESULT_OF_UNARY(log10_c)
matcl_RESULT_OF_UNARY(sin)
matcl_RESULT_OF_UNARY(cos)
matcl_RESULT_OF_UNARY(tan)
matcl_RESULT_OF_UNARY(cot)
matcl_RESULT_OF_UNARY(sec)
matcl_RESULT_OF_UNARY(csc)
matcl_RESULT_OF_UNARY(sinh)
matcl_RESULT_OF_UNARY(cosh)
matcl_RESULT_OF_UNARY(tanh)
matcl_RESULT_OF_UNARY(coth)
matcl_RESULT_OF_UNARY(sech)
matcl_RESULT_OF_UNARY(csch)
matcl_RESULT_OF_UNARY(asin)
matcl_RESULT_OF_UNARY(asin_c)
matcl_RESULT_OF_UNARY(acos)
matcl_RESULT_OF_UNARY(acos_c)
matcl_RESULT_OF_UNARY(atan)
matcl_RESULT_OF_UNARY(acot)
matcl_RESULT_OF_UNARY(asec)
matcl_RESULT_OF_UNARY(asec_c)
matcl_RESULT_OF_UNARY(acsc)
matcl_RESULT_OF_UNARY(acsc_c)
matcl_RESULT_OF_UNARY(asinh)
matcl_RESULT_OF_UNARY(acosh)
matcl_RESULT_OF_UNARY(acosh_c)
matcl_RESULT_OF_UNARY(atanh)
matcl_RESULT_OF_UNARY(atanh_c)
matcl_RESULT_OF_UNARY(acoth)
matcl_RESULT_OF_UNARY(acoth_c)
matcl_RESULT_OF_UNARY(asech)
matcl_RESULT_OF_UNARY(asech_c)
matcl_RESULT_OF_UNARY(acsch)
matcl_RESULT_OF_UNARY(inv)
matcl_RESULT_OF_UNARY(invs)
matcl_RESULT_OF_UNARY(logb)
matcl_RESULT_OF_UNARY(ilogb)
matcl_RESULT_OF_UNARY(floor)
matcl_RESULT_OF_UNARY(ceil)
matcl_RESULT_OF_UNARY(trunc)
matcl_RESULT_OF_UNARY(round)
matcl_RESULT_OF_UNARY(ifloor)
matcl_RESULT_OF_UNARY(iceil)
matcl_RESULT_OF_UNARY(itrunc)
matcl_RESULT_OF_UNARY(iround)
matcl_RESULT_OF_UNARY(cast_bool)

matcl_RESULT_OF_UNARY(is_nan)
matcl_RESULT_OF_UNARY(is_inf)
matcl_RESULT_OF_UNARY(is_finite)
matcl_RESULT_OF_UNARY(is_regular)
matcl_RESULT_OF_UNARY(is_normal)
matcl_RESULT_OF_UNARY(is_int)
matcl_RESULT_OF_UNARY(is_real)
matcl_RESULT_OF_UNARY(fpclassify)
matcl_RESULT_OF_UNARY(is_true)
matcl_RESULT_OF_UNARY(neg)

matcl_RESULT_OF_BINARY(atan2)
matcl_RESULT_OF_BINARY(hypot)
matcl_RESULT_OF_BINARY(mod)
matcl_RESULT_OF_BINARY(rem)
matcl_RESULT_OF_BINARY(copysign)
matcl_RESULT_OF_BINARY(nextafter)
matcl_RESULT_OF_BINARY(ldexp)
matcl_RESULT_OF_BINARY(scalbn)
matcl_RESULT_OF_BINARY(frexp)
matcl_RESULT_OF_BINARY(modf)
matcl_RESULT_OF_BINARY(fdim)
matcl_RESULT_OF_BINARY(pow)
matcl_RESULT_OF_BINARY(pow_c)
matcl_RESULT_OF_BINARY(powm1)
matcl_RESULT_OF_BINARY(idiv)
matcl_RESULT_OF_BINARY(op_xor)
matcl_RESULT_OF_BINARY(elem_and)
matcl_RESULT_OF_BINARY(elem_or)
matcl_RESULT_OF_BINARY(elem_xor)
matcl_RESULT_OF_BINARY(min)
matcl_RESULT_OF_BINARY(max)
matcl_RESULT_OF_BINARY(mul)
matcl_RESULT_OF_BINARY(rising_factorial)
matcl_RESULT_OF_BINARY(falling_factorial)
matcl_RESULT_OF_BINARY(eeq_nan)
matcl_RESULT_OF_BINARY(neq_nan)
matcl_RESULT_OF_BINARY(div_0)
matcl_RESULT_OF_BINARY(div_1)

matcl_RESULT_OF_TEMPL_1_1(bernoulli_b2n)
matcl_RESULT_OF_TEMPL_1_1(factorial)
matcl_RESULT_OF_TEMPL_1_1(double_factorial)
matcl_RESULT_OF_TEMPL_1_0(max_bernoulli_b2n)
matcl_RESULT_OF_TEMPL_1_2(binomial_coefficient)

matcl_RESULT_OF_UNARY_OP(uminus,-)
matcl_RESULT_OF_UNARY_OP(bool, (bool))
matcl_RESULT_OF_UNARY_OP(op_not,!)

matcl_RESULT_OF_BINARY_OP(eeq,==)
matcl_RESULT_OF_BINARY_OP(neq,!=)
matcl_RESULT_OF_BINARY_OP(geq,>=)
matcl_RESULT_OF_BINARY_OP(leq,<=)
matcl_RESULT_OF_BINARY_OP(gt,>)
matcl_RESULT_OF_BINARY_OP(lt,<)
matcl_RESULT_OF_BINARY_OP(op_and,&&)
matcl_RESULT_OF_BINARY_OP(op_or,||)

matcl_RESULT_OF_BINARY_OP(plus,+)
matcl_RESULT_OF_BINARY_FUN_ADD_INT(plus)

matcl_RESULT_OF_BINARY_OP(minus,-)
matcl_RESULT_OF_BINARY_FUN_ADD_INT(minus)

matcl_RESULT_OF_BINARY_OP(op_mul,*)
matcl_RESULT_OF_BINARY_FUN_ADD_INT(op_mul)

matcl_RESULT_OF_BINARY_OP(div,/)
matcl_RESULT_OF_BINARY_FUN_ADD_INT2(div)
