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

#include "matcl-mp-obj/mp_object.h"
#include "matcl-mp/func_binary.h"
#include "matcl-scalar/lib_functions/func_forwarding.h"
#include "matcl-dynamic/details/register_function_macro.h"
#include "matcl-dynamic/matcl_function_names.h"

namespace matcl
{

namespace mdy = matcl::dynamic;

//---------------------------------------------------------------
//                  OPERATOR + 
//---------------------------------------------------------------

#define REGISTER_ADD(num, arg1, arg2)                                       \
MATCL_REGISTER_BIN_FUNC(add_##num, operator+, arg1, arg2, mdy::functions::op_plus)

REGISTER_ADD(1, mp_int, mp_int);
REGISTER_ADD(2, mp_int, Integer);
REGISTER_ADD(3, Integer, mp_int);
REGISTER_ADD(4, mp_int, Real);
REGISTER_ADD(5, Real, mp_int);
REGISTER_ADD(6, mp_int  , Complex);
REGISTER_ADD(7, Complex, mp_int);
REGISTER_ADD(8, mp_int, Float_complex);
REGISTER_ADD(9, Float_complex, mp_int);
REGISTER_ADD(10, mp_rational, mp_rational);
REGISTER_ADD(11, mp_rational, Integer);
REGISTER_ADD(12, Integer, mp_rational);
REGISTER_ADD(13, mp_rational, Real);
REGISTER_ADD(14, Real, mp_rational);
REGISTER_ADD(15, mp_rational, Complex);
REGISTER_ADD(16, Complex, mp_rational);
REGISTER_ADD(17, mp_rational, Float_complex);
REGISTER_ADD(18, Float_complex, mp_rational);
REGISTER_ADD(19, mp_rational, mp_int);
REGISTER_ADD(20, mp_int, mp_rational);
REGISTER_ADD(21, mp_float, mp_float);
REGISTER_ADD(22, mp_float, Integer);
REGISTER_ADD(23, Integer, mp_float);
REGISTER_ADD(24, mp_float, Real);
REGISTER_ADD(25, Real, mp_float);
REGISTER_ADD(26, mp_float, Complex);
REGISTER_ADD(27, Complex, mp_float);
REGISTER_ADD(28, mp_float, Float_complex);
REGISTER_ADD(29, Float_complex, mp_float);
REGISTER_ADD(30, mp_float, mp_int);
REGISTER_ADD(31, mp_int, mp_float);
REGISTER_ADD(32, mp_float, mp_rational);
REGISTER_ADD(33, mp_rational, mp_float);
REGISTER_ADD(34, mp_complex, mp_complex);
REGISTER_ADD(35, mp_complex, Integer);
REGISTER_ADD(36, Integer, mp_complex);
REGISTER_ADD(37, mp_complex, Real);
REGISTER_ADD(38, Real, mp_complex);
REGISTER_ADD(39, mp_complex, Complex);
REGISTER_ADD(40, Complex, mp_complex);
REGISTER_ADD(41, mp_complex, Float_complex);
REGISTER_ADD(42, Float_complex, mp_complex);
REGISTER_ADD(43, mp_complex, mp_int);
REGISTER_ADD(44, mp_int, mp_complex);
REGISTER_ADD(45, mp_complex, mp_float);
REGISTER_ADD(46, mp_float, mp_complex);
REGISTER_ADD(47, mp_complex, mp_rational);
REGISTER_ADD(48, mp_rational, mp_complex);

//---------------------------------------------------------------
//                  OPERATOR -
//---------------------------------------------------------------

#define REGISTER_SUB(num, arg1, arg2)                                       \
MATCL_REGISTER_BIN_FUNC(sub_##num, operator-, arg1, arg2, mdy::functions::op_minus)

REGISTER_SUB(1, mp_int, mp_int);
REGISTER_SUB(2, mp_int, Integer);
REGISTER_SUB(3, Integer, mp_int);
REGISTER_SUB(4, mp_int, Real);
REGISTER_SUB(5, Real, mp_int);
REGISTER_SUB(6, mp_int  , Complex);
REGISTER_SUB(7, Complex, mp_int);
REGISTER_SUB(8, mp_int, Float_complex);
REGISTER_SUB(9, Float_complex, mp_int);
REGISTER_SUB(10, mp_rational, mp_rational);
REGISTER_SUB(11, mp_rational, Integer);
REGISTER_SUB(12, Integer, mp_rational);
REGISTER_SUB(13, mp_rational, Real);
REGISTER_SUB(14, Real, mp_rational);
REGISTER_SUB(15, mp_rational, Complex);
REGISTER_SUB(16, Complex, mp_rational);
REGISTER_SUB(17, mp_rational, Float_complex);
REGISTER_SUB(18, Float_complex, mp_rational);
REGISTER_SUB(19, mp_rational, mp_int);
REGISTER_SUB(20, mp_int, mp_rational);
REGISTER_SUB(21, mp_float, mp_float);
REGISTER_SUB(22, mp_float, Integer);
REGISTER_SUB(23, Integer, mp_float);
REGISTER_SUB(24, mp_float, Real);
REGISTER_SUB(25, Real, mp_float);
REGISTER_SUB(26, mp_float, Complex);
REGISTER_SUB(27, Complex, mp_float);
REGISTER_SUB(28, mp_float, Float_complex);
REGISTER_SUB(29, Float_complex, mp_float);
REGISTER_SUB(30, mp_float, mp_int);
REGISTER_SUB(31, mp_int, mp_float);
REGISTER_SUB(32, mp_float, mp_rational);
REGISTER_SUB(33, mp_rational, mp_float);
REGISTER_SUB(34, mp_complex, mp_complex);
REGISTER_SUB(35, mp_complex, Integer);
REGISTER_SUB(36, Integer, mp_complex);
REGISTER_SUB(37, mp_complex, Real);
REGISTER_SUB(38, Real, mp_complex);
REGISTER_SUB(39, mp_complex, Complex);
REGISTER_SUB(40, Complex, mp_complex);
REGISTER_SUB(41, mp_complex, Float_complex);
REGISTER_SUB(42, Float_complex, mp_complex);
REGISTER_SUB(43, mp_complex, mp_int);
REGISTER_SUB(44, mp_int, mp_complex);
REGISTER_SUB(45, mp_complex, mp_float);
REGISTER_SUB(46, mp_float, mp_complex);
REGISTER_SUB(47, mp_complex, mp_rational);
REGISTER_SUB(48, mp_rational, mp_complex);

//---------------------------------------------------------------
//                  OPERATOR *
//---------------------------------------------------------------

#define REGISTER_MUL(num, arg1, arg2)                                       \
MATCL_REGISTER_BIN_FUNC(mul_##num, operator*, arg1, arg2, mdy::functions::op_mul)

REGISTER_MUL(1, mp_int, mp_int);
REGISTER_MUL(2, mp_int, Integer);
REGISTER_MUL(3, Integer, mp_int);
REGISTER_MUL(4, mp_int, Real);
REGISTER_MUL(5, Real, mp_int);
REGISTER_MUL(6, mp_int  , Complex);
REGISTER_MUL(7, Complex, mp_int);
REGISTER_MUL(8, mp_int, Float_complex);
REGISTER_MUL(9, Float_complex, mp_int);
REGISTER_MUL(10, mp_rational, mp_rational);
REGISTER_MUL(11, mp_rational, Integer);
REGISTER_MUL(12, Integer, mp_rational);
REGISTER_MUL(13, mp_rational, Real);
REGISTER_MUL(14, Real, mp_rational);
REGISTER_MUL(15, mp_rational, Complex);
REGISTER_MUL(16, Complex, mp_rational);
REGISTER_MUL(17, mp_rational, Float_complex);
REGISTER_MUL(18, Float_complex, mp_rational);
REGISTER_MUL(19, mp_rational, mp_int);
REGISTER_MUL(20, mp_int, mp_rational);
REGISTER_MUL(21, mp_float, mp_float);
REGISTER_MUL(22, mp_float, Integer);
REGISTER_MUL(23, Integer, mp_float);
REGISTER_MUL(24, mp_float, Real);
REGISTER_MUL(25, Real, mp_float);
REGISTER_MUL(26, mp_float, Complex);
REGISTER_MUL(27, Complex, mp_float);
REGISTER_MUL(28, mp_float, Float_complex);
REGISTER_MUL(29, Float_complex, mp_float);
REGISTER_MUL(30, mp_float, mp_int);
REGISTER_MUL(31, mp_int, mp_float);
REGISTER_MUL(32, mp_float, mp_rational);
REGISTER_MUL(33, mp_rational, mp_float);
REGISTER_MUL(34, mp_complex, mp_complex);
REGISTER_MUL(35, mp_complex, Integer);
REGISTER_MUL(36, Integer, mp_complex);
REGISTER_MUL(37, mp_complex, Real);
REGISTER_MUL(38, Real, mp_complex);
REGISTER_MUL(39, mp_complex, Complex);
REGISTER_MUL(40, Complex, mp_complex);
REGISTER_MUL(41, mp_complex, Float_complex);
REGISTER_MUL(42, Float_complex, mp_complex);
REGISTER_MUL(43, mp_complex, mp_int);
REGISTER_MUL(44, mp_int, mp_complex);
REGISTER_MUL(45, mp_complex, mp_float);
REGISTER_MUL(46, mp_float, mp_complex);
REGISTER_MUL(47, mp_complex, mp_rational);
REGISTER_MUL(48, mp_rational, mp_complex);

//---------------------------------------------------------------
//                  OPERATOR /
//---------------------------------------------------------------

#define REGISTER_DIV(num, arg1, arg2)                                       \
MATCL_REGISTER_BIN_FUNC(div_##num, operator/, arg1, arg2, mdy::functions::op_div)

REGISTER_DIV(1, mp_int, mp_int);
REGISTER_DIV(2, mp_int, Integer);
REGISTER_DIV(3, Integer, mp_int);
REGISTER_DIV(4, mp_int, Real);
REGISTER_DIV(5, Real, mp_int);
REGISTER_DIV(6, mp_int  , Complex);
REGISTER_DIV(7, Complex, mp_int);
REGISTER_DIV(8, mp_int, Float_complex);
REGISTER_DIV(9, Float_complex, mp_int);
REGISTER_DIV(10, mp_rational, mp_rational);
REGISTER_DIV(11, mp_rational, Integer);
REGISTER_DIV(12, Integer, mp_rational);
REGISTER_DIV(13, mp_rational, Real);
REGISTER_DIV(14, Real, mp_rational);
REGISTER_DIV(15, mp_rational, Complex);
REGISTER_DIV(16, Complex, mp_rational);
REGISTER_DIV(17, mp_rational, Float_complex);
REGISTER_DIV(18, Float_complex, mp_rational);
REGISTER_DIV(19, mp_rational, mp_int);
REGISTER_DIV(20, mp_int, mp_rational);
REGISTER_DIV(21, mp_float, mp_float);
REGISTER_DIV(22, mp_float, Integer);
REGISTER_DIV(23, Integer, mp_float);
REGISTER_DIV(24, mp_float, Real);
REGISTER_DIV(25, Real, mp_float);
REGISTER_DIV(26, mp_float, Complex);
REGISTER_DIV(27, Complex, mp_float);
REGISTER_DIV(28, mp_float, Float_complex);
REGISTER_DIV(29, Float_complex, mp_float);
REGISTER_DIV(30, mp_float, mp_int);
REGISTER_DIV(31, mp_int, mp_float);
REGISTER_DIV(32, mp_float, mp_rational);
REGISTER_DIV(33, mp_rational, mp_float);
REGISTER_DIV(34, mp_complex, mp_complex);
REGISTER_DIV(35, mp_complex, Integer);
REGISTER_DIV(36, Integer, mp_complex);
REGISTER_DIV(37, mp_complex, Real);
REGISTER_DIV(38, Real, mp_complex);
REGISTER_DIV(39, mp_complex, Complex);
REGISTER_DIV(40, Complex, mp_complex);
REGISTER_DIV(41, mp_complex, Float_complex);
REGISTER_DIV(42, Float_complex, mp_complex);
REGISTER_DIV(43, mp_complex, mp_int);
REGISTER_DIV(44, mp_int, mp_complex);
REGISTER_DIV(45, mp_complex, mp_float);
REGISTER_DIV(46, mp_float, mp_complex);
REGISTER_DIV(47, mp_complex, mp_rational);
REGISTER_DIV(48, mp_rational, mp_complex);

//---------------------------------------------------------------
//                  IDIV
//---------------------------------------------------------------

#define REGISTER_IDIV(num, arg1, arg2)                                       \
MATCL_REGISTER_BIN_FUNC(idiv_##num, idiv, arg1, arg2, mdy::functions::idiv)

REGISTER_IDIV(1, mp_int, mp_int);
REGISTER_IDIV(2, mp_int, Integer);
REGISTER_IDIV(3, Integer, mp_int);
REGISTER_IDIV(4, mp_int, Real);
REGISTER_IDIV(5, Real, mp_int);
REGISTER_IDIV(6, mp_int  , Complex);
REGISTER_IDIV(7, Complex, mp_int);
REGISTER_IDIV(8, mp_int, Float_complex);
REGISTER_IDIV(9, Float_complex, mp_int);
REGISTER_IDIV(10, mp_rational, mp_rational);
REGISTER_IDIV(11, mp_rational, Integer);
REGISTER_IDIV(12, Integer, mp_rational);
REGISTER_IDIV(13, mp_rational, Real);
REGISTER_IDIV(14, Real, mp_rational);
REGISTER_IDIV(15, mp_rational, Complex);
REGISTER_IDIV(16, Complex, mp_rational);
REGISTER_IDIV(17, mp_rational, Float_complex);
REGISTER_IDIV(18, Float_complex, mp_rational);
REGISTER_IDIV(19, mp_rational, mp_int);
REGISTER_IDIV(20, mp_int, mp_rational);
REGISTER_IDIV(21, mp_float, mp_float);
REGISTER_IDIV(22, mp_float, Integer);
REGISTER_IDIV(23, Integer, mp_float);
REGISTER_IDIV(24, mp_float, Real);
REGISTER_IDIV(25, Real, mp_float);
REGISTER_IDIV(26, mp_float, Complex);
REGISTER_IDIV(27, Complex, mp_float);
REGISTER_IDIV(28, mp_float, Float_complex);
REGISTER_IDIV(29, Float_complex, mp_float);
REGISTER_IDIV(30, mp_float, mp_int);
REGISTER_IDIV(31, mp_int, mp_float);
REGISTER_IDIV(32, mp_float, mp_rational);
REGISTER_IDIV(33, mp_rational, mp_float);
REGISTER_IDIV(34, mp_complex, mp_complex);
REGISTER_IDIV(35, mp_complex, Integer);
REGISTER_IDIV(36, Integer, mp_complex);
REGISTER_IDIV(37, mp_complex, Real);
REGISTER_IDIV(38, Real, mp_complex);
REGISTER_IDIV(39, mp_complex, Complex);
REGISTER_IDIV(40, Complex, mp_complex);
REGISTER_IDIV(41, mp_complex, Float_complex);
REGISTER_IDIV(42, Float_complex, mp_complex);
REGISTER_IDIV(43, mp_complex, mp_int);
REGISTER_IDIV(44, mp_int, mp_complex);
REGISTER_IDIV(45, mp_complex, mp_float);
REGISTER_IDIV(46, mp_float, mp_complex);
REGISTER_IDIV(47, mp_complex, mp_rational);
REGISTER_IDIV(48, mp_rational, mp_complex);

//---------------------------------------------------------------
//                  POW
//---------------------------------------------------------------

#define REGISTER_POW(num, arg1, arg2)                                       \
MATCL_REGISTER_BIN_FUNC(pow_##num, pow, arg1, arg2, mdy::functions::pow)

REGISTER_POW(1, mp_int, mp_int)
REGISTER_POW(2, mp_int, Integer)
REGISTER_POW(3, Integer, mp_int)
REGISTER_POW(4, mp_rational, Integer)
REGISTER_POW(5, mp_rational, mp_int)
REGISTER_POW(6, mp_rational, mp_rational)
REGISTER_POW(7, Integer, mp_rational)
REGISTER_POW(8, mp_int, mp_rational)
REGISTER_POW(9, mp_float, mp_float)
REGISTER_POW(10, mp_float, Integer)
REGISTER_POW(11, Integer, mp_float)
REGISTER_POW(12, mp_float, mp_int)
REGISTER_POW(13, mp_int, mp_float)
REGISTER_POW(14, mp_float, mp_rational)
REGISTER_POW(15, mp_rational, mp_float)
REGISTER_POW(16, mp_complex, mp_complex)
REGISTER_POW(17, mp_complex, Integer)
REGISTER_POW(18, Integer, mp_complex)
REGISTER_POW(19, mp_complex, mp_int)
REGISTER_POW(20, mp_int, mp_complex)
REGISTER_POW(21, mp_complex, mp_float)
REGISTER_POW(22, mp_float, mp_complex)
REGISTER_POW(23, mp_complex, mp_rational)
REGISTER_POW(24, mp_rational, mp_complex)

//---------------------------------------------------------------
//                  POW_C
//---------------------------------------------------------------

#define REGISTER_POW_C(num, arg1, arg2)                                       \
MATCL_REGISTER_BIN_FUNC(pow_c_##num, pow_c, arg1, arg2, mdy::functions::pow_c)

REGISTER_POW_C(1, mp_int, mp_int)
REGISTER_POW_C(2, mp_int, Integer)
REGISTER_POW_C(3, Integer, mp_int)
REGISTER_POW_C(4, mp_rational, Integer)
REGISTER_POW_C(5, mp_rational, mp_int)
REGISTER_POW_C(6, mp_rational, mp_rational)
REGISTER_POW_C(7, Integer, mp_rational)
REGISTER_POW_C(8, mp_int, mp_rational)
REGISTER_POW_C(9, mp_float, mp_float)
REGISTER_POW_C(10, mp_float, Integer)
REGISTER_POW_C(11, Integer, mp_float)
REGISTER_POW_C(12, mp_float, mp_int)
REGISTER_POW_C(13, mp_int, mp_float)
REGISTER_POW_C(14, mp_float, mp_rational)
REGISTER_POW_C(15, mp_rational, mp_float)
REGISTER_POW_C(16, mp_complex, mp_complex)
REGISTER_POW_C(17, mp_complex, Integer)
REGISTER_POW_C(18, Integer, mp_complex)
REGISTER_POW_C(19, mp_complex, mp_int)
REGISTER_POW_C(20, mp_int, mp_complex)
REGISTER_POW_C(21, mp_complex, mp_float)
REGISTER_POW_C(22, mp_float, mp_complex)
REGISTER_POW_C(23, mp_complex, mp_rational)
REGISTER_POW_C(24, mp_rational, mp_complex)


};