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
#include "matcl-dynamic/matcl_function_names.h"

namespace matcl
{

namespace mdy = matcl::dynamic;

//---------------------------------------------------------------
//                  EQUALITY COMPARISON
//---------------------------------------------------------------

#define REGISTER_EEQ(num, arg1, arg2)                       \
struct reg_eeq_##num                                        \
    : mdy::register_function_ptr<bool (arg1, arg2),         \
                    &operator==, mdy::functions::op_eeq>{};

REGISTER_EEQ(1, const mp_int& val1, const mp_int& val2)
REGISTER_EEQ(2, const mp_int& val1, const Integer& val2)
REGISTER_EEQ(3, const Integer& val1, const mp_int& val2)
REGISTER_EEQ(4, const mp_int& val1, const Real& val2)
REGISTER_EEQ(5, const Real& val1, const mp_int& val2)
REGISTER_EEQ(6, const mp_int& val1, const Complex& val2)
REGISTER_EEQ(7, const Complex& val1, const mp_int& val2)
REGISTER_EEQ(8, const mp_int& val1, const Float_complex& val2)
REGISTER_EEQ(9, const Float_complex& val1, const mp_int& val2)
REGISTER_EEQ(10, const mp_rational& val1, const mp_rational& val2)
REGISTER_EEQ(11, const mp_rational& val1, const Integer& val2)
REGISTER_EEQ(12, const Integer& val1, const mp_rational& val2)
REGISTER_EEQ(13, const mp_rational& val1, const Real& val2)
REGISTER_EEQ(14, const Real& val1, const mp_rational& val2)
REGISTER_EEQ(15, const mp_rational& val1, const Complex& val2)
REGISTER_EEQ(16, const Complex& val1, const mp_rational& val2)
REGISTER_EEQ(17, const mp_rational& val1, const Float_complex& val2)
REGISTER_EEQ(18, const Float_complex& val1, const mp_rational& val2)
REGISTER_EEQ(19, const mp_rational& val1, const mp_int& val2)
REGISTER_EEQ(20, const mp_int& val1, const mp_rational& val2)
REGISTER_EEQ(21, const mp_float& val1, const mp_float& val2)
REGISTER_EEQ(22, const mp_float& val1, const Integer& val2)
REGISTER_EEQ(23, const Integer& val1, const mp_float& val2)
REGISTER_EEQ(24, const mp_float& val1, const Real& val2)
REGISTER_EEQ(25, const Real& val1, const mp_float& val2)
REGISTER_EEQ(26, const mp_float& val1, const Complex& val2)
REGISTER_EEQ(27, const Complex& val1, const mp_float& val2)
REGISTER_EEQ(28, const mp_float& val1, const Float_complex& val2)
REGISTER_EEQ(29, const Float_complex& val1, const mp_float& val2)
REGISTER_EEQ(30, const mp_float& val1, const mp_int& val2)
REGISTER_EEQ(31, const mp_int& val1, const mp_float& val2)
REGISTER_EEQ(32, const mp_float& val1, const mp_rational& val2)
REGISTER_EEQ(33, const mp_rational& val1, const mp_float& val2)
REGISTER_EEQ(34, const mp_complex& val1, const mp_complex& val2)
REGISTER_EEQ(35, const mp_complex& val1, const Integer& val2)
REGISTER_EEQ(36, const Integer& val1, const mp_complex& val2)
REGISTER_EEQ(37, const mp_complex& val1, const Real& val2)
REGISTER_EEQ(38, const Real& val1, const mp_complex& val2)
REGISTER_EEQ(39, const mp_complex& val1, const Complex& val2)
REGISTER_EEQ(40, const Complex& val1, const mp_complex& val2)
REGISTER_EEQ(41, const mp_complex& val1, const Float_complex& val2)
REGISTER_EEQ(42, const Float_complex& val1, const mp_complex& val2)
REGISTER_EEQ(43, const mp_complex& val1, const mp_int& val2)
REGISTER_EEQ(44, const mp_int& val1, const mp_complex& val2)
REGISTER_EEQ(45, const mp_complex& val1, const mp_float& val2)
REGISTER_EEQ(46, const mp_float& val1, const mp_complex& val2)
REGISTER_EEQ(47, const mp_complex& val1, const mp_rational& val2)
REGISTER_EEQ(48, const mp_rational& val1, const mp_complex& val2)

//---------------------------------------------------------------
//                  EQUALITY COMPARISON NAN
//---------------------------------------------------------------

#define REGISTER_EEQ_NAN(num, arg1, arg2)                   \
struct reg_eeq_nan_##num                                    \
    : mdy::register_function_ptr<bool (arg1, arg2),         \
                    &eeq_nan, mdy::functions::eeq_nan>{};

REGISTER_EEQ_NAN(1, const mp_int& val1, const mp_int& val2)
REGISTER_EEQ_NAN(2, const mp_int& val1, const Integer& val2)
REGISTER_EEQ_NAN(3, const Integer& val1, const mp_int& val2)
REGISTER_EEQ_NAN(4, const mp_int& val1, const Real& val2)
REGISTER_EEQ_NAN(5, const Real& val1, const mp_int& val2)
REGISTER_EEQ_NAN(6, const mp_int& val1, const Complex& val2)
REGISTER_EEQ_NAN(7, const Complex& val1, const mp_int& val2)
REGISTER_EEQ_NAN(8, const mp_int& val1, const Float_complex& val2)
REGISTER_EEQ_NAN(9, const Float_complex& val1, const mp_int& val2)
REGISTER_EEQ_NAN(10, const mp_rational& val1, const mp_rational& val2)
REGISTER_EEQ_NAN(11, const mp_rational& val1, const Integer& val2)
REGISTER_EEQ_NAN(12, const Integer& val1, const mp_rational& val2)
REGISTER_EEQ_NAN(13, const mp_rational& val1, const Real& val2)
REGISTER_EEQ_NAN(14, const Real& val1, const mp_rational& val2)
REGISTER_EEQ_NAN(15, const mp_rational& val1, const Complex& val2)
REGISTER_EEQ_NAN(16, const Complex& val1, const mp_rational& val2)
REGISTER_EEQ_NAN(17, const mp_rational& val1, const Float_complex& val2)
REGISTER_EEQ_NAN(18, const Float_complex& val1, const mp_rational& val2)
REGISTER_EEQ_NAN(19, const mp_rational& val1, const mp_int& val2)
REGISTER_EEQ_NAN(20, const mp_int& val1, const mp_rational& val2)
REGISTER_EEQ_NAN(21, const mp_float& val1, const mp_float& val2)
REGISTER_EEQ_NAN(22, const mp_float& val1, const Integer& val2)
REGISTER_EEQ_NAN(23, const Integer& val1, const mp_float& val2)
REGISTER_EEQ_NAN(24, const mp_float& val1, const Real& val2)
REGISTER_EEQ_NAN(25, const Real& val1, const mp_float& val2)
REGISTER_EEQ_NAN(26, const mp_float& val1, const Complex& val2)
REGISTER_EEQ_NAN(27, const Complex& val1, const mp_float& val2)
REGISTER_EEQ_NAN(28, const mp_float& val1, const Float_complex& val2)
REGISTER_EEQ_NAN(29, const Float_complex& val1, const mp_float& val2)
REGISTER_EEQ_NAN(30, const mp_float& val1, const mp_int& val2)
REGISTER_EEQ_NAN(31, const mp_int& val1, const mp_float& val2)
REGISTER_EEQ_NAN(32, const mp_float& val1, const mp_rational& val2)
REGISTER_EEQ_NAN(33, const mp_rational& val1, const mp_float& val2)
REGISTER_EEQ_NAN(34, const mp_complex& val1, const mp_complex& val2)
REGISTER_EEQ_NAN(35, const mp_complex& val1, const Integer& val2)
REGISTER_EEQ_NAN(36, const Integer& val1, const mp_complex& val2)
REGISTER_EEQ_NAN(37, const mp_complex& val1, const Real& val2)
REGISTER_EEQ_NAN(38, const Real& val1, const mp_complex& val2)
REGISTER_EEQ_NAN(39, const mp_complex& val1, const Complex& val2)
REGISTER_EEQ_NAN(40, const Complex& val1, const mp_complex& val2)
REGISTER_EEQ_NAN(41, const mp_complex& val1, const Float_complex& val2)
REGISTER_EEQ_NAN(42, const Float_complex& val1, const mp_complex& val2)
REGISTER_EEQ_NAN(43, const mp_complex& val1, const mp_int& val2)
REGISTER_EEQ_NAN(44, const mp_int& val1, const mp_complex& val2)
REGISTER_EEQ_NAN(45, const mp_complex& val1, const mp_float& val2)
REGISTER_EEQ_NAN(46, const mp_float& val1, const mp_complex& val2)
REGISTER_EEQ_NAN(47, const mp_complex& val1, const mp_rational& val2)
REGISTER_EEQ_NAN(48, const mp_rational& val1, const mp_complex& val2)

//---------------------------------------------------------------
//                  INEQUALITY COMPARISON
//---------------------------------------------------------------

#define REGISTER_NEQ(num, arg1, arg2)                       \
struct reg_neq_##num                                        \
    : mdy::register_function_ptr<bool (arg1, arg2),         \
                    &operator!=, mdy::functions::op_neq>{};

REGISTER_NEQ(1, const mp_int& val1, const mp_int& val2)
REGISTER_NEQ(2, const mp_int& val1, const Integer& val2)
REGISTER_NEQ(3, const Integer& val1, const mp_int& val2)
REGISTER_NEQ(4, const mp_int& val1, const Real& val2)
REGISTER_NEQ(5, const Real& val1, const mp_int& val2)
REGISTER_NEQ(6, const mp_int& val1, const Complex& val2)
REGISTER_NEQ(7, const Complex& val1, const mp_int& val2)
REGISTER_NEQ(8, const mp_int& val1, const Float_complex& val2)
REGISTER_NEQ(9, const Float_complex& val1, const mp_int& val2)
REGISTER_NEQ(10, const mp_rational& val1, const mp_rational& val2)
REGISTER_NEQ(11, const mp_rational& val1, const Integer& val2)
REGISTER_NEQ(12, const Integer& val1, const mp_rational& val2)
REGISTER_NEQ(13, const mp_rational& val1, const Real& val2)
REGISTER_NEQ(14, const Real& val1, const mp_rational& val2)
REGISTER_NEQ(15, const mp_rational& val1, const Complex& val2)
REGISTER_NEQ(16, const Complex& val1, const mp_rational& val2)
REGISTER_NEQ(17, const mp_rational& val1, const Float_complex& val2)
REGISTER_NEQ(18, const Float_complex& val1, const mp_rational& val2)
REGISTER_NEQ(19, const mp_rational& val1, const mp_int& val2)
REGISTER_NEQ(20, const mp_int& val1, const mp_rational& val2)
REGISTER_NEQ(21, const mp_float& val1, const mp_float& val2)
REGISTER_NEQ(22, const mp_float& val1, const Integer& val2)
REGISTER_NEQ(23, const Integer& val1, const mp_float& val2)
REGISTER_NEQ(24, const mp_float& val1, const Real& val2)
REGISTER_NEQ(25, const Real& val1, const mp_float& val2)
REGISTER_NEQ(26, const mp_float& val1, const Complex& val2)
REGISTER_NEQ(27, const Complex& val1, const mp_float& val2)
REGISTER_NEQ(28, const mp_float& val1, const Float_complex& val2)
REGISTER_NEQ(29, const Float_complex& val1, const mp_float& val2)
REGISTER_NEQ(30, const mp_float& val1, const mp_int& val2)
REGISTER_NEQ(31, const mp_int& val1, const mp_float& val2)
REGISTER_NEQ(32, const mp_float& val1, const mp_rational& val2)
REGISTER_NEQ(33, const mp_rational& val1, const mp_float& val2)
REGISTER_NEQ(34, const mp_complex& val1, const mp_complex& val2)
REGISTER_NEQ(35, const mp_complex& val1, const Integer& val2)
REGISTER_NEQ(36, const Integer& val1, const mp_complex& val2)
REGISTER_NEQ(37, const mp_complex& val1, const Real& val2)
REGISTER_NEQ(38, const Real& val1, const mp_complex& val2)
REGISTER_NEQ(39, const mp_complex& val1, const Complex& val2)
REGISTER_NEQ(40, const Complex& val1, const mp_complex& val2)
REGISTER_NEQ(41, const mp_complex& val1, const Float_complex& val2)
REGISTER_NEQ(42, const Float_complex& val1, const mp_complex& val2)
REGISTER_NEQ(43, const mp_complex& val1, const mp_int& val2)
REGISTER_NEQ(44, const mp_int& val1, const mp_complex& val2)
REGISTER_NEQ(45, const mp_complex& val1, const mp_float& val2)
REGISTER_NEQ(46, const mp_float& val1, const mp_complex& val2)
REGISTER_NEQ(47, const mp_complex& val1, const mp_rational& val2)
REGISTER_NEQ(48, const mp_rational& val1, const mp_complex& val2)

//---------------------------------------------------------------
//                  INEQUALITY COMPARISON NAN
//---------------------------------------------------------------

#define REGISTER_NEQ_NAN(num, arg1, arg2)                   \
struct reg_neq_nan##num                                     \
    : mdy::register_function_ptr<bool (arg1, arg2),         \
                    &neq_nan, mdy::functions::neq_nan>{};

REGISTER_NEQ_NAN(1, const mp_int& val1, const mp_int& val2)
REGISTER_NEQ_NAN(2, const mp_int& val1, const Integer& val2)
REGISTER_NEQ_NAN(3, const Integer& val1, const mp_int& val2)
REGISTER_NEQ_NAN(4, const mp_int& val1, const Real& val2)
REGISTER_NEQ_NAN(5, const Real& val1, const mp_int& val2)
REGISTER_NEQ_NAN(6, const mp_int& val1, const Complex& val2)
REGISTER_NEQ_NAN(7, const Complex& val1, const mp_int& val2)
REGISTER_NEQ_NAN(8, const mp_int& val1, const Float_complex& val2)
REGISTER_NEQ_NAN(9, const Float_complex& val1, const mp_int& val2)
REGISTER_NEQ_NAN(10, const mp_rational& val1, const mp_rational& val2)
REGISTER_NEQ_NAN(11, const mp_rational& val1, const Integer& val2)
REGISTER_NEQ_NAN(12, const Integer& val1, const mp_rational& val2)
REGISTER_NEQ_NAN(13, const mp_rational& val1, const Real& val2)
REGISTER_NEQ_NAN(14, const Real& val1, const mp_rational& val2)
REGISTER_NEQ_NAN(15, const mp_rational& val1, const Complex& val2)
REGISTER_NEQ_NAN(16, const Complex& val1, const mp_rational& val2)
REGISTER_NEQ_NAN(17, const mp_rational& val1, const Float_complex& val2)
REGISTER_NEQ_NAN(18, const Float_complex& val1, const mp_rational& val2)
REGISTER_NEQ_NAN(19, const mp_rational& val1, const mp_int& val2)
REGISTER_NEQ_NAN(20, const mp_int& val1, const mp_rational& val2)
REGISTER_NEQ_NAN(21, const mp_float& val1, const mp_float& val2)
REGISTER_NEQ_NAN(22, const mp_float& val1, const Integer& val2)
REGISTER_NEQ_NAN(23, const Integer& val1, const mp_float& val2)
REGISTER_NEQ_NAN(24, const mp_float& val1, const Real& val2)
REGISTER_NEQ_NAN(25, const Real& val1, const mp_float& val2)
REGISTER_NEQ_NAN(26, const mp_float& val1, const Complex& val2)
REGISTER_NEQ_NAN(27, const Complex& val1, const mp_float& val2)
REGISTER_NEQ_NAN(28, const mp_float& val1, const Float_complex& val2)
REGISTER_NEQ_NAN(29, const Float_complex& val1, const mp_float& val2)
REGISTER_NEQ_NAN(30, const mp_float& val1, const mp_int& val2)
REGISTER_NEQ_NAN(31, const mp_int& val1, const mp_float& val2)
REGISTER_NEQ_NAN(32, const mp_float& val1, const mp_rational& val2)
REGISTER_NEQ_NAN(33, const mp_rational& val1, const mp_float& val2)
REGISTER_NEQ_NAN(34, const mp_complex& val1, const mp_complex& val2)
REGISTER_NEQ_NAN(35, const mp_complex& val1, const Integer& val2)
REGISTER_NEQ_NAN(36, const Integer& val1, const mp_complex& val2)
REGISTER_NEQ_NAN(37, const mp_complex& val1, const Real& val2)
REGISTER_NEQ_NAN(38, const Real& val1, const mp_complex& val2)
REGISTER_NEQ_NAN(39, const mp_complex& val1, const Complex& val2)
REGISTER_NEQ_NAN(40, const Complex& val1, const mp_complex& val2)
REGISTER_NEQ_NAN(41, const mp_complex& val1, const Float_complex& val2)
REGISTER_NEQ_NAN(42, const Float_complex& val1, const mp_complex& val2)
REGISTER_NEQ_NAN(43, const mp_complex& val1, const mp_int& val2)
REGISTER_NEQ_NAN(44, const mp_int& val1, const mp_complex& val2)
REGISTER_NEQ_NAN(45, const mp_complex& val1, const mp_float& val2)
REGISTER_NEQ_NAN(46, const mp_float& val1, const mp_complex& val2)
REGISTER_NEQ_NAN(47, const mp_complex& val1, const mp_rational& val2)
REGISTER_NEQ_NAN(48, const mp_rational& val1, const mp_complex& val2)

//---------------------------------------------------------------
//                  GREATER OR EQUAL COMPARISON
//---------------------------------------------------------------
#define REGISTER_GEQ(num, arg1, arg2)                       \
struct reg_geq_##num                                        \
    : mdy::register_function_ptr<bool (arg1, arg2),         \
                    &operator>=, mdy::functions::op_geq>{};

REGISTER_GEQ(1, const mp_int& val1, const mp_int& val2)
REGISTER_GEQ(2, const mp_int& val1, const Integer& val2)
REGISTER_GEQ(3, const Integer& val1, const mp_int& val2)
REGISTER_GEQ(4, const mp_int& val1, const Real& val2)
REGISTER_GEQ(5, const Real& val1, const mp_int& val2)
REGISTER_GEQ(6, const mp_int& val1, const Complex& val2)
REGISTER_GEQ(7, const Complex& val1, const mp_int& val2)
REGISTER_GEQ(8, const mp_int& val1, const Float_complex& val2)
REGISTER_GEQ(9, const Float_complex& val1, const mp_int& val2)
REGISTER_GEQ(10, const mp_rational& val1, const mp_rational& val2)
REGISTER_GEQ(11, const mp_rational& val1, const Integer& val2)
REGISTER_GEQ(12, const Integer& val1, const mp_rational& val2)
REGISTER_GEQ(13, const mp_rational& val1, const Real& val2)
REGISTER_GEQ(14, const Real& val1, const mp_rational& val2)
REGISTER_GEQ(15, const mp_rational& val1, const Complex& val2)
REGISTER_GEQ(16, const Complex& val1, const mp_rational& val2)
REGISTER_GEQ(17, const mp_rational& val1, const Float_complex& val2)
REGISTER_GEQ(18, const Float_complex& val1, const mp_rational& val2)
REGISTER_GEQ(19, const mp_rational& val1, const mp_int& val2)
REGISTER_GEQ(20, const mp_int& val1, const mp_rational& val2)
REGISTER_GEQ(21, const mp_float& val1, const mp_float& val2)
REGISTER_GEQ(22, const mp_float& val1, const Integer& val2)
REGISTER_GEQ(23, const Integer& val1, const mp_float& val2)
REGISTER_GEQ(24, const mp_float& val1, const Real& val2)
REGISTER_GEQ(25, const Real& val1, const mp_float& val2)
REGISTER_GEQ(26, const mp_float& val1, const Complex& val2)
REGISTER_GEQ(27, const Complex& val1, const mp_float& val2)
REGISTER_GEQ(28, const mp_float& val1, const Float_complex& val2)
REGISTER_GEQ(29, const Float_complex& val1, const mp_float& val2)
REGISTER_GEQ(30, const mp_float& val1, const mp_int& val2)
REGISTER_GEQ(31, const mp_int& val1, const mp_float& val2)
REGISTER_GEQ(32, const mp_float& val1, const mp_rational& val2)
REGISTER_GEQ(33, const mp_rational& val1, const mp_float& val2)
REGISTER_GEQ(34, const mp_complex& val1, const mp_complex& val2)
REGISTER_GEQ(35, const mp_complex& val1, const Integer& val2)
REGISTER_GEQ(36, const Integer& val1, const mp_complex& val2)
REGISTER_GEQ(37, const mp_complex& val1, const Real& val2)
REGISTER_GEQ(38, const Real& val1, const mp_complex& val2)
REGISTER_GEQ(39, const mp_complex& val1, const Complex& val2)
REGISTER_GEQ(40, const Complex& val1, const mp_complex& val2)
REGISTER_GEQ(41, const mp_complex& val1, const Float_complex& val2)
REGISTER_GEQ(42, const Float_complex& val1, const mp_complex& val2)
REGISTER_GEQ(43, const mp_complex& val1, const mp_int& val2)
REGISTER_GEQ(44, const mp_int& val1, const mp_complex& val2)
REGISTER_GEQ(45, const mp_complex& val1, const mp_float& val2)
REGISTER_GEQ(46, const mp_float& val1, const mp_complex& val2)
REGISTER_GEQ(47, const mp_complex& val1, const mp_rational& val2)
REGISTER_GEQ(48, const mp_rational& val1, const mp_complex& val2)

//---------------------------------------------------------------
//                  LESS OR EQUAL COMPARISON
//---------------------------------------------------------------
#define REGISTER_LEQ(num, arg1, arg2)                       \
struct reg_leq_##num                                        \
    : mdy::register_function_ptr<bool (arg1, arg2),         \
                    &operator<=, mdy::functions::op_leq>{};

REGISTER_LEQ(1, const mp_int& val1, const mp_int& val2)
REGISTER_LEQ(2, const mp_int& val1, const Integer& val2)
REGISTER_LEQ(3, const Integer& val1, const mp_int& val2)
REGISTER_LEQ(4, const mp_int& val1, const Real& val2)
REGISTER_LEQ(5, const Real& val1, const mp_int& val2)
REGISTER_LEQ(6, const mp_int& val1, const Complex& val2)
REGISTER_LEQ(7, const Complex& val1, const mp_int& val2)
REGISTER_LEQ(8, const mp_int& val1, const Float_complex& val2)
REGISTER_LEQ(9, const Float_complex& val1, const mp_int& val2)
REGISTER_LEQ(10, const mp_rational& val1, const mp_rational& val2)
REGISTER_LEQ(11, const mp_rational& val1, const Integer& val2)
REGISTER_LEQ(12, const Integer& val1, const mp_rational& val2)
REGISTER_LEQ(13, const mp_rational& val1, const Real& val2)
REGISTER_LEQ(14, const Real& val1, const mp_rational& val2)
REGISTER_LEQ(15, const mp_rational& val1, const Complex& val2)
REGISTER_LEQ(16, const Complex& val1, const mp_rational& val2)
REGISTER_LEQ(17, const mp_rational& val1, const Float_complex& val2)
REGISTER_LEQ(18, const Float_complex& val1, const mp_rational& val2)
REGISTER_LEQ(19, const mp_rational& val1, const mp_int& val2)
REGISTER_LEQ(20, const mp_int& val1, const mp_rational& val2)
REGISTER_LEQ(21, const mp_float& val1, const mp_float& val2)
REGISTER_LEQ(22, const mp_float& val1, const Integer& val2)
REGISTER_LEQ(23, const Integer& val1, const mp_float& val2)
REGISTER_LEQ(24, const mp_float& val1, const Real& val2)
REGISTER_LEQ(25, const Real& val1, const mp_float& val2)
REGISTER_LEQ(26, const mp_float& val1, const Complex& val2)
REGISTER_LEQ(27, const Complex& val1, const mp_float& val2)
REGISTER_LEQ(28, const mp_float& val1, const Float_complex& val2)
REGISTER_LEQ(29, const Float_complex& val1, const mp_float& val2)
REGISTER_LEQ(30, const mp_float& val1, const mp_int& val2)
REGISTER_LEQ(31, const mp_int& val1, const mp_float& val2)
REGISTER_LEQ(32, const mp_float& val1, const mp_rational& val2)
REGISTER_LEQ(33, const mp_rational& val1, const mp_float& val2)
REGISTER_LEQ(34, const mp_complex& val1, const mp_complex& val2)
REGISTER_LEQ(35, const mp_complex& val1, const Integer& val2)
REGISTER_LEQ(36, const Integer& val1, const mp_complex& val2)
REGISTER_LEQ(37, const mp_complex& val1, const Real& val2)
REGISTER_LEQ(38, const Real& val1, const mp_complex& val2)
REGISTER_LEQ(39, const mp_complex& val1, const Complex& val2)
REGISTER_LEQ(40, const Complex& val1, const mp_complex& val2)
REGISTER_LEQ(41, const mp_complex& val1, const Float_complex& val2)
REGISTER_LEQ(42, const Float_complex& val1, const mp_complex& val2)
REGISTER_LEQ(43, const mp_complex& val1, const mp_int& val2)
REGISTER_LEQ(44, const mp_int& val1, const mp_complex& val2)
REGISTER_LEQ(45, const mp_complex& val1, const mp_float& val2)
REGISTER_LEQ(46, const mp_float& val1, const mp_complex& val2)
REGISTER_LEQ(47, const mp_complex& val1, const mp_rational& val2)
REGISTER_LEQ(48, const mp_rational& val1, const mp_complex& val2)

//---------------------------------------------------------------
//                  GREATER THAN COMPARISON
//---------------------------------------------------------------
#define REGISTER_GT(num, arg1, arg2)                        \
struct reg_gt_##num                                         \
    : mdy::register_function_ptr<bool (arg1, arg2),         \
                    &(operator>), mdy::functions::op_gt>{};

REGISTER_GT(1, const mp_int& val1, const mp_int& val2)
REGISTER_GT(2, const mp_int& val1, const Integer& val2)
REGISTER_GT(3, const Integer& val1, const mp_int& val2)
REGISTER_GT(4, const mp_int& val1, const Real& val2)
REGISTER_GT(5, const Real& val1, const mp_int& val2)
REGISTER_GT(6, const mp_int& val1, const Complex& val2)
REGISTER_GT(7, const Complex& val1, const mp_int& val2)
REGISTER_GT(8, const mp_int& val1, const Float_complex& val2)
REGISTER_GT(9, const Float_complex& val1, const mp_int& val2)
REGISTER_GT(10, const mp_rational& val1, const mp_rational& val2)
REGISTER_GT(11, const mp_rational& val1, const Integer& val2)
REGISTER_GT(12, const Integer& val1, const mp_rational& val2)
REGISTER_GT(13, const mp_rational& val1, const Real& val2)
REGISTER_GT(14, const Real& val1, const mp_rational& val2)
REGISTER_GT(15, const mp_rational& val1, const Complex& val2)
REGISTER_GT(16, const Complex& val1, const mp_rational& val2)
REGISTER_GT(17, const mp_rational& val1, const Float_complex& val2)
REGISTER_GT(18, const Float_complex& val1, const mp_rational& val2)
REGISTER_GT(19, const mp_rational& val1, const mp_int& val2)
REGISTER_GT(20, const mp_int& val1, const mp_rational& val2)
REGISTER_GT(21, const mp_float& val1, const mp_float& val2)
REGISTER_GT(22, const mp_float& val1, const Integer& val2)
REGISTER_GT(23, const Integer& val1, const mp_float& val2)
REGISTER_GT(24, const mp_float& val1, const Real& val2)
REGISTER_GT(25, const Real& val1, const mp_float& val2)
REGISTER_GT(26, const mp_float& val1, const Complex& val2)
REGISTER_GT(27, const Complex& val1, const mp_float& val2)
REGISTER_GT(28, const mp_float& val1, const Float_complex& val2)
REGISTER_GT(29, const Float_complex& val1, const mp_float& val2)
REGISTER_GT(30, const mp_float& val1, const mp_int& val2)
REGISTER_GT(31, const mp_int& val1, const mp_float& val2)
REGISTER_GT(32, const mp_float& val1, const mp_rational& val2)
REGISTER_GT(33, const mp_rational& val1, const mp_float& val2)
REGISTER_GT(34, const mp_complex& val1, const mp_complex& val2)
REGISTER_GT(35, const mp_complex& val1, const Integer& val2)
REGISTER_GT(36, const Integer& val1, const mp_complex& val2)
REGISTER_GT(37, const mp_complex& val1, const Real& val2)
REGISTER_GT(38, const Real& val1, const mp_complex& val2)
REGISTER_GT(39, const mp_complex& val1, const Complex& val2)
REGISTER_GT(40, const Complex& val1, const mp_complex& val2)
REGISTER_GT(41, const mp_complex& val1, const Float_complex& val2)
REGISTER_GT(42, const Float_complex& val1, const mp_complex& val2)
REGISTER_GT(43, const mp_complex& val1, const mp_int& val2)
REGISTER_GT(44, const mp_int& val1, const mp_complex& val2)
REGISTER_GT(45, const mp_complex& val1, const mp_float& val2)
REGISTER_GT(46, const mp_float& val1, const mp_complex& val2)
REGISTER_GT(47, const mp_complex& val1, const mp_rational& val2)
REGISTER_GT(48, const mp_rational& val1, const mp_complex& val2)

//---------------------------------------------------------------
//                  LESS THAN COMPARISON
//---------------------------------------------------------------
#define REGISTER_LT(num, arg1, arg2)                        \
struct reg_lt_##num                                         \
    : mdy::register_function_ptr<bool (arg1, arg2),         \
                    &operator<, mdy::functions::op_lt>{};

REGISTER_LT(1, const mp_int& val1, const mp_int& val2)
REGISTER_LT(2, const mp_int& val1, const Integer& val2)
REGISTER_LT(3, const Integer& val1, const mp_int& val2)
REGISTER_LT(4, const mp_int& val1, const Real& val2)
REGISTER_LT(5, const Real& val1, const mp_int& val2)
REGISTER_LT(6, const mp_int& val1, const Complex& val2)
REGISTER_LT(7, const Complex& val1, const mp_int& val2)
REGISTER_LT(8, const mp_int& val1, const Float_complex& val2)
REGISTER_LT(9, const Float_complex& val1, const mp_int& val2)
REGISTER_LT(10, const mp_rational& val1, const mp_rational& val2)
REGISTER_LT(11, const mp_rational& val1, const Integer& val2)
REGISTER_LT(12, const Integer& val1, const mp_rational& val2)
REGISTER_LT(13, const mp_rational& val1, const Real& val2)
REGISTER_LT(14, const Real& val1, const mp_rational& val2)
REGISTER_LT(15, const mp_rational& val1, const Complex& val2)
REGISTER_LT(16, const Complex& val1, const mp_rational& val2)
REGISTER_LT(17, const mp_rational& val1, const Float_complex& val2)
REGISTER_LT(18, const Float_complex& val1, const mp_rational& val2)
REGISTER_LT(19, const mp_rational& val1, const mp_int& val2)
REGISTER_LT(20, const mp_int& val1, const mp_rational& val2)
REGISTER_LT(21, const mp_float& val1, const mp_float& val2)
REGISTER_LT(22, const mp_float& val1, const Integer& val2)
REGISTER_LT(23, const Integer& val1, const mp_float& val2)
REGISTER_LT(24, const mp_float& val1, const Real& val2)
REGISTER_LT(25, const Real& val1, const mp_float& val2)
REGISTER_LT(26, const mp_float& val1, const Complex& val2)
REGISTER_LT(27, const Complex& val1, const mp_float& val2)
REGISTER_LT(28, const mp_float& val1, const Float_complex& val2)
REGISTER_LT(29, const Float_complex& val1, const mp_float& val2)
REGISTER_LT(30, const mp_float& val1, const mp_int& val2)
REGISTER_LT(31, const mp_int& val1, const mp_float& val2)
REGISTER_LT(32, const mp_float& val1, const mp_rational& val2)
REGISTER_LT(33, const mp_rational& val1, const mp_float& val2)
REGISTER_LT(34, const mp_complex& val1, const mp_complex& val2)
REGISTER_LT(35, const mp_complex& val1, const Integer& val2)
REGISTER_LT(36, const Integer& val1, const mp_complex& val2)
REGISTER_LT(37, const mp_complex& val1, const Real& val2)
REGISTER_LT(38, const Real& val1, const mp_complex& val2)
REGISTER_LT(39, const mp_complex& val1, const Complex& val2)
REGISTER_LT(40, const Complex& val1, const mp_complex& val2)
REGISTER_LT(41, const mp_complex& val1, const Float_complex& val2)
REGISTER_LT(42, const Float_complex& val1, const mp_complex& val2)
REGISTER_LT(43, const mp_complex& val1, const mp_int& val2)
REGISTER_LT(44, const mp_int& val1, const mp_complex& val2)
REGISTER_LT(45, const mp_complex& val1, const mp_float& val2)
REGISTER_LT(46, const mp_float& val1, const mp_complex& val2)
REGISTER_LT(47, const mp_complex& val1, const mp_rational& val2)
REGISTER_LT(48, const mp_rational& val1, const mp_complex& val2)

};