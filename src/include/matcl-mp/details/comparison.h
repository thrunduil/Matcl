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

#include "matcl-core/config.h"
#include "matcl-mp/config.h"
#include "matcl-mp/details/fwd_decls.h"
#include "matcl-core/general/fwd_decls.h"
#include "matcl-core/matrix/complex_type.h"

namespace matcl { namespace mp { namespace details
{

//---------------------------------------------------------------
//                  GENERAL COMPARISON
//---------------------------------------------------------------
enum class cmp_type
{
    not_ordered,    // one of argument is nan
    less,           // a < b
    equal,          // a == b
    greater         // a > b
};

// compare two numbers
cmp_type    compare(const mp_float& a, const mp_float& b);

// compare abs(a) and abs(b)
cmp_type    compare_abs(const mp_float& a, const mp_float& b);

//---------------------------------------------------------------
//                  EQUALITY COMPARISON
//---------------------------------------------------------------
MATCL_MP_EXPORT bool    eeq(const mp_int& val1, const mp_int& val2);
MATCL_MP_EXPORT bool    eeq(const mp_int& val1, Integer val2);
MATCL_MP_EXPORT bool    eeq(Integer val1, const mp_int& val2);
MATCL_MP_EXPORT bool    eeq(const mp_int& val1, Real val2);
MATCL_MP_EXPORT bool    eeq(Real val1, const mp_int& val2);
MATCL_MP_EXPORT bool    eeq(const mp_int& val1, const Complex& val2);
MATCL_MP_EXPORT bool    eeq(const Complex& val1, const mp_int& val2);
MATCL_MP_EXPORT bool    eeq(const mp_int& val1, const Float_complex& val2);
MATCL_MP_EXPORT bool    eeq(const Float_complex& val1, const mp_int& val2);

MATCL_MP_EXPORT bool    eeq(const mp_rational& val1, const mp_rational& val2);
MATCL_MP_EXPORT bool    eeq(const mp_rational& val1, Integer val2);
MATCL_MP_EXPORT bool    eeq(Integer val1, const mp_rational& val2);
MATCL_MP_EXPORT bool    eeq(const mp_rational& val1, Real val2);
MATCL_MP_EXPORT bool    eeq(Real val1, const mp_rational& val2);
MATCL_MP_EXPORT bool    eeq(const mp_rational& val1, const Complex& val2);
MATCL_MP_EXPORT bool    eeq(const Complex& val1, const mp_rational& val2);
MATCL_MP_EXPORT bool    eeq(const mp_rational& val1, const Float_complex& val2);
MATCL_MP_EXPORT bool    eeq(const Float_complex& val1, const mp_rational& val2);
MATCL_MP_EXPORT bool    eeq(const mp_rational& val1, const mp_int& val2);
MATCL_MP_EXPORT bool    eeq(const mp_int& val1, const mp_rational& val2);

MATCL_MP_EXPORT bool    eeq(const mp_float& val1, const mp_float& val2);
MATCL_MP_EXPORT bool    eeq(const mp_float& val1, Integer val2);
MATCL_MP_EXPORT bool    eeq(Integer val1, const mp_float& val2);
MATCL_MP_EXPORT bool    eeq(const mp_float& val1, Real val2);
MATCL_MP_EXPORT bool    eeq(Real val1, const mp_float& val2);
MATCL_MP_EXPORT bool    eeq(const mp_float& val1, const Complex& val2);
MATCL_MP_EXPORT bool    eeq(const Complex& val1, const mp_float& val2);
MATCL_MP_EXPORT bool    eeq(const mp_float& val1, const Float_complex& val2);
MATCL_MP_EXPORT bool    eeq(const Float_complex& val1, const mp_float& val2);
MATCL_MP_EXPORT bool    eeq(const mp_float& val1, const mp_int& val2);
MATCL_MP_EXPORT bool    eeq(const mp_int& val1, const mp_float& val2);
MATCL_MP_EXPORT bool    eeq(const mp_float& val1, const mp_rational& val2);
MATCL_MP_EXPORT bool    eeq(const mp_rational& val1, const mp_float& val2);

MATCL_MP_EXPORT bool    eeq(const mp_complex& val1, const mp_complex& val2);
MATCL_MP_EXPORT bool    eeq(const mp_complex& val1, Integer val2);
MATCL_MP_EXPORT bool    eeq(Integer val1, const mp_complex& val2);
MATCL_MP_EXPORT bool    eeq(const mp_complex& val1, Real val2);
MATCL_MP_EXPORT bool    eeq(Real val1, const mp_complex& val2);
MATCL_MP_EXPORT bool    eeq(const mp_complex& val1, const Complex& val2);
MATCL_MP_EXPORT bool    eeq(const Complex& val1, const mp_complex& val2);
MATCL_MP_EXPORT bool    eeq(const mp_complex& val1, const Float_complex& val2);
MATCL_MP_EXPORT bool    eeq(const Float_complex& val1, const mp_complex& val2);
MATCL_MP_EXPORT bool    eeq(const mp_complex& val1, const mp_int& val2);
MATCL_MP_EXPORT bool    eeq(const mp_int& val1, const mp_complex& val2);
MATCL_MP_EXPORT bool    eeq(const mp_complex& val1, const mp_float& val2);
MATCL_MP_EXPORT bool    eeq(const mp_float& val1, const mp_complex& val2);
MATCL_MP_EXPORT bool    eeq(const mp_complex& val1, const mp_rational& val2);
MATCL_MP_EXPORT bool    eeq(const mp_rational& val1, const mp_complex& val2);

//---------------------------------------------------------------
//                  INEQUALITY COMPARISON
//---------------------------------------------------------------
MATCL_MP_EXPORT bool    neq(const mp_int& val1, const mp_int& val2);
MATCL_MP_EXPORT bool    neq(const mp_int& val1, Real val2);
MATCL_MP_EXPORT bool    neq(Real val1, const mp_int& val2);
MATCL_MP_EXPORT bool    neq(const mp_int& val1, Integer val2);
MATCL_MP_EXPORT bool    neq(const mp_int& val1, const Complex& val2);
MATCL_MP_EXPORT bool    neq(const Complex& val1, const mp_int& val2);
MATCL_MP_EXPORT bool    neq(Integer val1, const mp_int& val2);
MATCL_MP_EXPORT bool    neq(const mp_int& val1, const Float_complex& val2);
MATCL_MP_EXPORT bool    neq(const Float_complex& val1, const mp_int& val2);

MATCL_MP_EXPORT bool    neq(const mp_rational& val1, const mp_rational& val2);
MATCL_MP_EXPORT bool    neq(const mp_rational& val1, Integer val2);
MATCL_MP_EXPORT bool    neq(Integer val1, const mp_rational& val2);
MATCL_MP_EXPORT bool    neq(const mp_rational& val1, Real val2);
MATCL_MP_EXPORT bool    neq(Real val1, const mp_rational& val2);
MATCL_MP_EXPORT bool    neq(const mp_rational& val1, const Complex& val2);
MATCL_MP_EXPORT bool    neq(const Complex& val1, const mp_rational& val2);
MATCL_MP_EXPORT bool    neq(const mp_rational& val1, const Float_complex& val2);
MATCL_MP_EXPORT bool    neq(const Float_complex& val1, const mp_rational& val2);
MATCL_MP_EXPORT bool    neq(const mp_rational& val1, const mp_int& val2);
MATCL_MP_EXPORT bool    neq(const mp_int& val1, const mp_rational& val2);

MATCL_MP_EXPORT bool    neq(const mp_float& val1, const mp_float& val2);
MATCL_MP_EXPORT bool    neq(const mp_float& val1, Integer val2);
MATCL_MP_EXPORT bool    neq(Integer val1, const mp_float& val2);
MATCL_MP_EXPORT bool    neq(const mp_float& val1, Real val2);
MATCL_MP_EXPORT bool    neq(Real val1, const mp_float& val2);
MATCL_MP_EXPORT bool    neq(const mp_float& val1, const Complex& val2);
MATCL_MP_EXPORT bool    neq(const Complex& val1, const mp_float& val2);
MATCL_MP_EXPORT bool    neq(const mp_float& val1, const Float_complex& val2);
MATCL_MP_EXPORT bool    neq(const Float_complex& val1, const mp_float& val2);
MATCL_MP_EXPORT bool    neq(const mp_float& val1, const mp_int& val2);
MATCL_MP_EXPORT bool    neq(const mp_int& val1, const mp_float& val2);
MATCL_MP_EXPORT bool    neq(const mp_float& val1, const mp_rational& val2);
MATCL_MP_EXPORT bool    neq(const mp_rational& val1, const mp_float& val2);

MATCL_MP_EXPORT bool    neq(const mp_complex& val1, const mp_complex& val2);
MATCL_MP_EXPORT bool    neq(const mp_complex& val1, Integer val2);
MATCL_MP_EXPORT bool    neq(Integer val1, const mp_complex& val2);
MATCL_MP_EXPORT bool    neq(const mp_complex& val1, Real val2);
MATCL_MP_EXPORT bool    neq(Real val1, const mp_complex& val2);
MATCL_MP_EXPORT bool    neq(const mp_complex& val1, const Complex& val2);
MATCL_MP_EXPORT bool    neq(const Complex& val1, const mp_complex& val2);
MATCL_MP_EXPORT bool    neq(const mp_complex& val1, const Float_complex& val2);
MATCL_MP_EXPORT bool    neq(const Float_complex& val1, const mp_complex& val2);
MATCL_MP_EXPORT bool    neq(const mp_complex& val1, const mp_int& val2);
MATCL_MP_EXPORT bool    neq(const mp_int& val1, const mp_complex& val2);
MATCL_MP_EXPORT bool    neq(const mp_complex& val1, const mp_float& val2);
MATCL_MP_EXPORT bool    neq(const mp_float& val1, const mp_complex& val2);
MATCL_MP_EXPORT bool    neq(const mp_complex& val1, const mp_rational& val2);
MATCL_MP_EXPORT bool    neq(const mp_rational& val1, const mp_complex& val2);

//---------------------------------------------------------------
//                  GREATER OR EQUAL COMPARISON
//---------------------------------------------------------------
MATCL_MP_EXPORT bool    geq(const mp_int& val1, const mp_int& val2);
MATCL_MP_EXPORT bool    geq(const mp_int& val1, Real val2);
MATCL_MP_EXPORT bool    geq(Real val1, const mp_int& val2);
MATCL_MP_EXPORT bool    geq(const mp_int& val1, Integer val2);
MATCL_MP_EXPORT bool    geq(Integer val1, const mp_int& val2);
MATCL_MP_EXPORT bool    geq(const mp_int& val1, const Complex& val2);
MATCL_MP_EXPORT bool    geq(const Complex& val1, const mp_int& val2);
MATCL_MP_EXPORT bool    geq(const mp_int& val1, const Float_complex& val2);
MATCL_MP_EXPORT bool    geq(const Float_complex& val1, const mp_int& val2);

MATCL_MP_EXPORT bool    geq(const mp_rational& val1, const mp_rational& val2);
MATCL_MP_EXPORT bool    geq(const mp_rational& val1, Integer val2);
MATCL_MP_EXPORT bool    geq(Integer val1, const mp_rational& val2);
MATCL_MP_EXPORT bool    geq(const mp_rational& val1, Real val2);
MATCL_MP_EXPORT bool    geq(Real val1, const mp_rational& val2);
MATCL_MP_EXPORT bool    geq(const mp_rational& val1, const Complex& val2);
MATCL_MP_EXPORT bool    geq(const Complex& val1, const mp_rational& val2);
MATCL_MP_EXPORT bool    geq(const mp_rational& val1, const Float_complex& val2);
MATCL_MP_EXPORT bool    geq(const Float_complex& val1, const mp_rational& val2);
MATCL_MP_EXPORT bool    geq(const mp_rational& val1, const mp_int& val2);
MATCL_MP_EXPORT bool    geq(const mp_int& val1, const mp_rational& val2);

MATCL_MP_EXPORT bool    geq(const mp_float& val1, const mp_float& val2);
MATCL_MP_EXPORT bool    geq(const mp_float& val1, Integer val2);
MATCL_MP_EXPORT bool    geq(Integer val1, const mp_float& val2);
MATCL_MP_EXPORT bool    geq(const mp_float& val1, Real val2);
MATCL_MP_EXPORT bool    geq(Real val1, const mp_float& val2);
MATCL_MP_EXPORT bool    geq(const mp_float& val1, const Complex& val2);
MATCL_MP_EXPORT bool    geq(const Complex& val1, const mp_float& val2);
MATCL_MP_EXPORT bool    geq(const mp_float& val1, const Float_complex& val2);
MATCL_MP_EXPORT bool    geq(const Float_complex& val1, const mp_float& val2);
MATCL_MP_EXPORT bool    geq(const mp_float& val1, const mp_int& val2);
MATCL_MP_EXPORT bool    geq(const mp_int& val1, const mp_float& val2);
MATCL_MP_EXPORT bool    geq(const mp_float& val1, const mp_rational& val2);
MATCL_MP_EXPORT bool    geq(const mp_rational& val1, const mp_float& val2);

MATCL_MP_EXPORT bool    geq(const mp_complex& val1, const mp_complex& val2);
MATCL_MP_EXPORT bool    geq(const mp_complex& val1, Integer val2);
MATCL_MP_EXPORT bool    geq(Integer val1, const mp_complex& val2);
MATCL_MP_EXPORT bool    geq(const mp_complex& val1, Real val2);
MATCL_MP_EXPORT bool    geq(Real val1, const mp_complex& val2);
MATCL_MP_EXPORT bool    geq(const mp_complex& val1, const Complex& val2);
MATCL_MP_EXPORT bool    geq(const Complex& val1, const mp_complex& val2);
MATCL_MP_EXPORT bool    geq(const mp_complex& val1, const Float_complex& val2);
MATCL_MP_EXPORT bool    geq(const Float_complex& val1, const mp_complex& val2);
MATCL_MP_EXPORT bool    geq(const mp_complex& val1, const mp_int& val2);
MATCL_MP_EXPORT bool    geq(const mp_int& val1, const mp_complex& val2);
MATCL_MP_EXPORT bool    geq(const mp_complex& val1, const mp_float& val2);
MATCL_MP_EXPORT bool    geq(const mp_float& val1, const mp_complex& val2);
MATCL_MP_EXPORT bool    geq(const mp_complex& val1, const mp_rational& val2);
MATCL_MP_EXPORT bool    geq(const mp_rational& val1, const mp_complex& val2);

//---------------------------------------------------------------
//                  LESS OR EQUAL COMPARISON
//---------------------------------------------------------------
MATCL_MP_EXPORT bool    leq(const mp_int& val1, const mp_int& val2);
MATCL_MP_EXPORT bool    leq(const mp_int& val1, Real val2);
MATCL_MP_EXPORT bool    leq(Real val1, const mp_int& val2);
MATCL_MP_EXPORT bool    leq(const mp_int& val1, Integer val2);
MATCL_MP_EXPORT bool    leq(Integer val1, const mp_int& val2);
MATCL_MP_EXPORT bool    leq(const mp_int& val1, const Complex& val2);
MATCL_MP_EXPORT bool    leq(const Complex& val1, const mp_int& val2);
MATCL_MP_EXPORT bool    leq(const mp_int& val1, const Float_complex& val2);
MATCL_MP_EXPORT bool    leq(const Float_complex& val1, const mp_int& val2);

MATCL_MP_EXPORT bool    leq(const mp_rational& val1, const mp_rational& val2);
MATCL_MP_EXPORT bool    leq(const mp_rational& val1, Integer val2);
MATCL_MP_EXPORT bool    leq(Integer val1, const mp_rational& val2);
MATCL_MP_EXPORT bool    leq(const mp_rational& val1, Real val2);
MATCL_MP_EXPORT bool    leq(Real val1, const mp_rational& val2);
MATCL_MP_EXPORT bool    leq(const mp_rational& val1, const Complex& val2);
MATCL_MP_EXPORT bool    leq(const Complex& val1, const mp_rational& val2);
MATCL_MP_EXPORT bool    leq(const mp_rational& val1, const Float_complex& val2);
MATCL_MP_EXPORT bool    leq(const Float_complex& val1, const mp_rational& val2);
MATCL_MP_EXPORT bool    leq(const mp_rational& val1, const mp_int& val2);
MATCL_MP_EXPORT bool    leq(const mp_int& val1, const mp_rational& val2);

MATCL_MP_EXPORT bool    leq(const mp_float& val1, const mp_float& val2);
MATCL_MP_EXPORT bool    leq(const mp_float& val1, Integer val2);
MATCL_MP_EXPORT bool    leq(Integer val1, const mp_float& val2);
MATCL_MP_EXPORT bool    leq(const mp_float& val1, Real val2);
MATCL_MP_EXPORT bool    leq(Real val1, const mp_float& val2);
MATCL_MP_EXPORT bool    leq(const mp_float& val1, const Complex& val2);
MATCL_MP_EXPORT bool    leq(const Complex& val1, const mp_float& val2);
MATCL_MP_EXPORT bool    leq(const mp_float& val1, const Float_complex& val2);
MATCL_MP_EXPORT bool    leq(const Float_complex& val1, const mp_float& val2);
MATCL_MP_EXPORT bool    leq(const mp_float& val1, const mp_int& val2);
MATCL_MP_EXPORT bool    leq(const mp_int& val1, const mp_float& val2);
MATCL_MP_EXPORT bool    leq(const mp_float& val1, const mp_rational& val2);
MATCL_MP_EXPORT bool    leq(const mp_rational& val1, const mp_float& val2);

MATCL_MP_EXPORT bool    leq(const mp_complex& val1, const mp_complex& val2);
MATCL_MP_EXPORT bool    leq(const mp_complex& val1, Integer val2);
MATCL_MP_EXPORT bool    leq(Integer val1, const mp_complex& val2);
MATCL_MP_EXPORT bool    leq(const mp_complex& val1, Real val2);
MATCL_MP_EXPORT bool    leq(Real val1, const mp_complex& val2);
MATCL_MP_EXPORT bool    leq(const mp_complex& val1, const Complex& val2);
MATCL_MP_EXPORT bool    leq(const Complex& val1, const mp_complex& val2);
MATCL_MP_EXPORT bool    leq(const mp_complex& val1, const Float_complex& val2);
MATCL_MP_EXPORT bool    leq(const Float_complex& val1, const mp_complex& val2);
MATCL_MP_EXPORT bool    leq(const mp_complex& val1, const mp_int& val2);
MATCL_MP_EXPORT bool    leq(const mp_int& val1, const mp_complex& val2);
MATCL_MP_EXPORT bool    leq(const mp_complex& val1, const mp_float& val2);
MATCL_MP_EXPORT bool    leq(const mp_float& val1, const mp_complex& val2);
MATCL_MP_EXPORT bool    leq(const mp_complex& val1, const mp_rational& val2);
MATCL_MP_EXPORT bool    leq(const mp_rational& val1, const mp_complex& val2);

//---------------------------------------------------------------
//                  GREATER THAN COMPARISON
//---------------------------------------------------------------
MATCL_MP_EXPORT bool    gt(const mp_int& val1, const mp_int& val2);
MATCL_MP_EXPORT bool    gt(const mp_int& val1, Real val2);
MATCL_MP_EXPORT bool    gt(Real val1, const mp_int& val2);
MATCL_MP_EXPORT bool    gt(const mp_int& val1, Integer val2);
MATCL_MP_EXPORT bool    gt(Integer val1, const mp_int& val2);
MATCL_MP_EXPORT bool    gt(const mp_int& val1, const Complex& val2);
MATCL_MP_EXPORT bool    gt(const Complex& val1, const mp_int& val2);
MATCL_MP_EXPORT bool    gt(const mp_int& val1, const Float_complex& val2);
MATCL_MP_EXPORT bool    gt(const Float_complex& val1, const mp_int& val2);

MATCL_MP_EXPORT bool    gt(const mp_rational& val1, const mp_rational& val2);
MATCL_MP_EXPORT bool    gt(const mp_rational& val1, Integer val2);
MATCL_MP_EXPORT bool    gt(Integer val1, const mp_rational& val2);
MATCL_MP_EXPORT bool    gt(const mp_rational& val1, Real val2);
MATCL_MP_EXPORT bool    gt(Real val1, const mp_rational& val2);
MATCL_MP_EXPORT bool    gt(const mp_rational& val1, const Complex& val2);
MATCL_MP_EXPORT bool    gt(const Complex& val1, const mp_rational& val2);
MATCL_MP_EXPORT bool    gt(const mp_rational& val1, const Float_complex& val2);
MATCL_MP_EXPORT bool    gt(const Float_complex& val1, const mp_rational& val2);
MATCL_MP_EXPORT bool    gt(const mp_rational& val1, const mp_int& val2);
MATCL_MP_EXPORT bool    gt(const mp_int& val1, const mp_rational& val2);

MATCL_MP_EXPORT bool    gt(const mp_float& val1, const mp_float& val2);
MATCL_MP_EXPORT bool    gt(const mp_float& val1, Integer val2);
MATCL_MP_EXPORT bool    gt(Integer val1, const mp_float& val2);
MATCL_MP_EXPORT bool    gt(const mp_float& val1, Real val2);
MATCL_MP_EXPORT bool    gt(Real val1, const mp_float& val2);
MATCL_MP_EXPORT bool    gt(const mp_float& val1, const Complex& val2);
MATCL_MP_EXPORT bool    gt(const Complex& val1, const mp_float& val2);
MATCL_MP_EXPORT bool    gt(const mp_float& val1, const Float_complex& val2);
MATCL_MP_EXPORT bool    gt(const Float_complex& val1, const mp_float& val2);
MATCL_MP_EXPORT bool    gt(const mp_float& val1, const mp_int& val2);
MATCL_MP_EXPORT bool    gt(const mp_int& val1, const mp_float& val2);
MATCL_MP_EXPORT bool    gt(const mp_float& val1, const mp_rational& val2);
MATCL_MP_EXPORT bool    gt(const mp_rational& val1, const mp_float& val2);

MATCL_MP_EXPORT bool    gt(const mp_complex& val1, const mp_complex& val2);
MATCL_MP_EXPORT bool    gt(const mp_complex& val1, Integer val2);
MATCL_MP_EXPORT bool    gt(Integer val1, const mp_complex& val2);
MATCL_MP_EXPORT bool    gt(const mp_complex& val1, Real val2);
MATCL_MP_EXPORT bool    gt(Real val1, const mp_complex& val2);
MATCL_MP_EXPORT bool    gt(const mp_complex& val1, const Complex& val2);
MATCL_MP_EXPORT bool    gt(const Complex& val1, const mp_complex& val2);
MATCL_MP_EXPORT bool    gt(const mp_complex& val1, const Float_complex& val2);
MATCL_MP_EXPORT bool    gt(const Float_complex& val1, const mp_complex& val2);
MATCL_MP_EXPORT bool    gt(const mp_complex& val1, const mp_int& val2);
MATCL_MP_EXPORT bool    gt(const mp_int& val1, const mp_complex& val2);
MATCL_MP_EXPORT bool    gt(const mp_complex& val1, const mp_float& val2);
MATCL_MP_EXPORT bool    gt(const mp_float& val1, const mp_complex& val2);
MATCL_MP_EXPORT bool    gt(const mp_complex& val1, const mp_rational& val2);
MATCL_MP_EXPORT bool    gt(const mp_rational& val1, const mp_complex& val2);

//---------------------------------------------------------------
//                  LESS THAN COMPARISON
//---------------------------------------------------------------
MATCL_MP_EXPORT bool    lt(const mp_int& val1, const mp_int& val2);
MATCL_MP_EXPORT bool    lt(const mp_int& val1, Real val2);
MATCL_MP_EXPORT bool    lt(Real val1, const mp_int& val2);
MATCL_MP_EXPORT bool    lt(const mp_int& val1, Integer val2);
MATCL_MP_EXPORT bool    lt(Integer val1, const mp_int& val2);
MATCL_MP_EXPORT bool    lt(const mp_int& val1, const Complex& val2);
MATCL_MP_EXPORT bool    lt(const Complex& val1, const mp_int& val2);
MATCL_MP_EXPORT bool    lt(const mp_int& val1, const Float_complex& val2);
MATCL_MP_EXPORT bool    lt(const Float_complex& val1, const mp_int& val2);

MATCL_MP_EXPORT bool    lt(const mp_rational& val1, const mp_rational& val2);
MATCL_MP_EXPORT bool    lt(const mp_rational& val1, Integer val2);
MATCL_MP_EXPORT bool    lt(Integer val1, const mp_rational& val2);
MATCL_MP_EXPORT bool    lt(const mp_rational& val1, Real val2);
MATCL_MP_EXPORT bool    lt(Real val1, const mp_rational& val2);
MATCL_MP_EXPORT bool    lt(const mp_rational& val1, const Complex& val2);
MATCL_MP_EXPORT bool    lt(const Complex& val1, const mp_rational& val2);
MATCL_MP_EXPORT bool    lt(const mp_rational& val1, const Float_complex& val2);
MATCL_MP_EXPORT bool    lt(const Float_complex& val1, const mp_rational& val2);
MATCL_MP_EXPORT bool    lt(const mp_rational& val1, const mp_int& val2);
MATCL_MP_EXPORT bool    lt(const mp_int& val1, const mp_rational& val2);

MATCL_MP_EXPORT bool    lt(const mp_float& val1, const mp_float& val2);
MATCL_MP_EXPORT bool    lt(const mp_float& val1, Integer val2);
MATCL_MP_EXPORT bool    lt(Integer val1, const mp_float& val2);
MATCL_MP_EXPORT bool    lt(const mp_float& val1, Real val2);
MATCL_MP_EXPORT bool    lt(Real val1, const mp_float& val2);
MATCL_MP_EXPORT bool    lt(const mp_float& val1, const Complex& val2);
MATCL_MP_EXPORT bool    lt(const Complex& val1, const mp_float& val2);
MATCL_MP_EXPORT bool    lt(const mp_float& val1, const Float_complex& val2);
MATCL_MP_EXPORT bool    lt(const Float_complex& val1, const mp_float& val2);
MATCL_MP_EXPORT bool    lt(const mp_float& val1, const mp_int& val2);
MATCL_MP_EXPORT bool    lt(const mp_int& val1, const mp_float& val2);
MATCL_MP_EXPORT bool    lt(const mp_float& val1, const mp_rational& val2);
MATCL_MP_EXPORT bool    lt(const mp_rational& val1, const mp_float& val2);

MATCL_MP_EXPORT bool    lt(const mp_complex& val1, const mp_complex& val2);
MATCL_MP_EXPORT bool    lt(const mp_complex& val1, Integer val2);
MATCL_MP_EXPORT bool    lt(Integer val1, const mp_complex& val2);
MATCL_MP_EXPORT bool    lt(const mp_complex& val1, Real val2);
MATCL_MP_EXPORT bool    lt(Real val1, const mp_complex& val2);
MATCL_MP_EXPORT bool    lt(const mp_complex& val1, const Complex& val2);
MATCL_MP_EXPORT bool    lt(const Complex& val1, const mp_complex& val2);
MATCL_MP_EXPORT bool    lt(const mp_complex& val1, const Float_complex& val2);
MATCL_MP_EXPORT bool    lt(const Float_complex& val1, const mp_complex& val2);
MATCL_MP_EXPORT bool    lt(const mp_complex& val1, const mp_int& val2);
MATCL_MP_EXPORT bool    lt(const mp_int& val1, const mp_complex& val2);
MATCL_MP_EXPORT bool    lt(const mp_complex& val1, const mp_float& val2);
MATCL_MP_EXPORT bool    lt(const mp_float& val1, const mp_complex& val2);
MATCL_MP_EXPORT bool    lt(const mp_complex& val1, const mp_rational& val2);
MATCL_MP_EXPORT bool    lt(const mp_rational& val1, const mp_complex& val2);

//---------------------------------------------------------------
//                  EQUALITY COMPARISON
// nan values are equal
//---------------------------------------------------------------
MATCL_MP_EXPORT bool    eeq_nan(const mp_int& val1, const mp_int& val2);
MATCL_MP_EXPORT bool    eeq_nan(const mp_int& val1, Integer val2);
MATCL_MP_EXPORT bool    eeq_nan(Integer val1, const mp_int& val2);
MATCL_MP_EXPORT bool    eeq_nan(const mp_int& val1, Float val2);
MATCL_MP_EXPORT bool    eeq_nan(const mp_int& val1, Real val2);
MATCL_MP_EXPORT bool    eeq_nan(Float val1, const mp_int& val2);
MATCL_MP_EXPORT bool    eeq_nan(Real val1, const mp_int& val2);
MATCL_MP_EXPORT bool    eeq_nan(const mp_int& val1, const Complex& val2);
MATCL_MP_EXPORT bool    eeq_nan(const Complex& val1, const mp_int& val2);
MATCL_MP_EXPORT bool    eeq_nan(const mp_int& val1, const Float_complex& val2);
MATCL_MP_EXPORT bool    eeq_nan(const Float_complex& val1, const mp_int& val2);

MATCL_MP_EXPORT bool    eeq_nan(const mp_rational& val1, const mp_rational& val2);
MATCL_MP_EXPORT bool    eeq_nan(const mp_rational& val1, Integer val2);
MATCL_MP_EXPORT bool    eeq_nan(Integer val1, const mp_rational& val2);
MATCL_MP_EXPORT bool    eeq_nan(const mp_rational& val1, Float val2);
MATCL_MP_EXPORT bool    eeq_nan(const mp_rational& val1, Real val2);
MATCL_MP_EXPORT bool    eeq_nan(Float val1, const mp_rational& val2);
MATCL_MP_EXPORT bool    eeq_nan(Real val1, const mp_rational& val2);
MATCL_MP_EXPORT bool    eeq_nan(const mp_rational& val1, const Complex& val2);
MATCL_MP_EXPORT bool    eeq_nan(const Complex& val1, const mp_rational& val2);
MATCL_MP_EXPORT bool    eeq_nan(const mp_rational& val1, const Float_complex& val2);
MATCL_MP_EXPORT bool    eeq_nan(const Float_complex& val1, const mp_rational& val2);
MATCL_MP_EXPORT bool    eeq_nan(const mp_rational& val1, const mp_int& val2);
MATCL_MP_EXPORT bool    eeq_nan(const mp_int& val1, const mp_rational& val2);

MATCL_MP_EXPORT bool    eeq_nan(const mp_float& val1, const mp_float& val2);
MATCL_MP_EXPORT bool    eeq_nan(const mp_float& val1, Integer val2);
MATCL_MP_EXPORT bool    eeq_nan(Integer val1, const mp_float& val2);
MATCL_MP_EXPORT bool    eeq_nan(const mp_float& val1, Float val2);
MATCL_MP_EXPORT bool    eeq_nan(const mp_float& val1, Real val2);
MATCL_MP_EXPORT bool    eeq_nan(Float val1, const mp_float& val2);
MATCL_MP_EXPORT bool    eeq_nan(Real val1, const mp_float& val2);
MATCL_MP_EXPORT bool    eeq_nan(const mp_float& val1, const Complex& val2);
MATCL_MP_EXPORT bool    eeq_nan(const Complex& val1, const mp_float& val2);
MATCL_MP_EXPORT bool    eeq_nan(const mp_float& val1, const Float_complex& val2);
MATCL_MP_EXPORT bool    eeq_nan(const Float_complex& val1, const mp_float& val2);
MATCL_MP_EXPORT bool    eeq_nan(const mp_float& val1, const mp_int& val2);
MATCL_MP_EXPORT bool    eeq_nan(const mp_int& val1, const mp_float& val2);
MATCL_MP_EXPORT bool    eeq_nan(const mp_float& val1, const mp_rational& val2);
MATCL_MP_EXPORT bool    eeq_nan(const mp_rational& val1, const mp_float& val2);

MATCL_MP_EXPORT bool    eeq_nan(const mp_complex& val1, const mp_complex& val2);
MATCL_MP_EXPORT bool    eeq_nan(const mp_complex& val1, Integer val2);
MATCL_MP_EXPORT bool    eeq_nan(Integer val1, const mp_complex& val2);
MATCL_MP_EXPORT bool    eeq_nan(const mp_complex& val1, Float val2);
MATCL_MP_EXPORT bool    eeq_nan(const mp_complex& val1, Real val2);
MATCL_MP_EXPORT bool    eeq_nan(Float val1, const mp_complex& val2);
MATCL_MP_EXPORT bool    eeq_nan(Real val1, const mp_complex& val2);
MATCL_MP_EXPORT bool    eeq_nan(const mp_complex& val1, const Complex& val2);
MATCL_MP_EXPORT bool    eeq_nan(const Complex& val1, const mp_complex& val2);
MATCL_MP_EXPORT bool    eeq_nan(const mp_complex& val1, const Float_complex& val2);
MATCL_MP_EXPORT bool    eeq_nan(const Float_complex& val1, const mp_complex& val2);
MATCL_MP_EXPORT bool    eeq_nan(const mp_complex& val1, const mp_int& val2);
MATCL_MP_EXPORT bool    eeq_nan(const mp_int& val1, const mp_complex& val2);
MATCL_MP_EXPORT bool    eeq_nan(const mp_complex& val1, const mp_float& val2);
MATCL_MP_EXPORT bool    eeq_nan(const mp_float& val1, const mp_complex& val2);
MATCL_MP_EXPORT bool    eeq_nan(const mp_complex& val1, const mp_rational& val2);
MATCL_MP_EXPORT bool    eeq_nan(const mp_rational& val1, const mp_complex& val2);

//---------------------------------------------------------------
//                  INEQUALITY COMPARISON
// nan values are equal
//---------------------------------------------------------------
MATCL_MP_EXPORT bool    neq_nan(const mp_int& val1, const mp_int& val2);
MATCL_MP_EXPORT bool    neq_nan(const mp_int& val1, Float val2);
MATCL_MP_EXPORT bool    neq_nan(const mp_int& val1, Real val2);
MATCL_MP_EXPORT bool    neq_nan(Float val1, const mp_int& val2);
MATCL_MP_EXPORT bool    neq_nan(Real val1, const mp_int& val2);
MATCL_MP_EXPORT bool    neq_nan(const mp_int& val1, Integer val2);
MATCL_MP_EXPORT bool    neq_nan(const mp_int& val1, const Complex& val2);
MATCL_MP_EXPORT bool    neq_nan(const Complex& val1, const mp_int& val2);
MATCL_MP_EXPORT bool    neq_nan(Integer val1, const mp_int& val2);
MATCL_MP_EXPORT bool    neq_nan(const mp_int& val1, const Float_complex& val2);
MATCL_MP_EXPORT bool    neq_nan(const Float_complex& val1, const mp_int& val2);

MATCL_MP_EXPORT bool    neq_nan(const mp_rational& val1, const mp_rational& val2);
MATCL_MP_EXPORT bool    neq_nan(const mp_rational& val1, Integer val2);
MATCL_MP_EXPORT bool    neq_nan(Integer val1, const mp_rational& val2);
MATCL_MP_EXPORT bool    neq_nan(const mp_rational& val1, Float val2);
MATCL_MP_EXPORT bool    neq_nan(const mp_rational& val1, Real val2);
MATCL_MP_EXPORT bool    neq_nan(Float val1, const mp_rational& val2);
MATCL_MP_EXPORT bool    neq_nan(Real val1, const mp_rational& val2);
MATCL_MP_EXPORT bool    neq_nan(const mp_rational& val1, const Complex& val2);
MATCL_MP_EXPORT bool    neq_nan(const Complex& val1, const mp_rational& val2);
MATCL_MP_EXPORT bool    neq_nan(const mp_rational& val1, const Float_complex& val2);
MATCL_MP_EXPORT bool    neq_nan(const Float_complex& val1, const mp_rational& val2);
MATCL_MP_EXPORT bool    neq_nan(const mp_rational& val1, const mp_int& val2);
MATCL_MP_EXPORT bool    neq_nan(const mp_int& val1, const mp_rational& val2);

MATCL_MP_EXPORT bool    neq_nan(const mp_float& val1, const mp_float& val2);
MATCL_MP_EXPORT bool    neq_nan(const mp_float& val1, Integer val2);
MATCL_MP_EXPORT bool    neq_nan(Integer val1, const mp_float& val2);
MATCL_MP_EXPORT bool    neq_nan(const mp_float& val1, Float val2);
MATCL_MP_EXPORT bool    neq_nan(const mp_float& val1, Real val2);
MATCL_MP_EXPORT bool    neq_nan(Float val1, const mp_float& val2);
MATCL_MP_EXPORT bool    neq_nan(Real val1, const mp_float& val2);
MATCL_MP_EXPORT bool    neq_nan(const mp_float& val1, const Complex& val2);
MATCL_MP_EXPORT bool    neq_nan(const Complex& val1, const mp_float& val2);
MATCL_MP_EXPORT bool    neq_nan(const mp_float& val1, const Float_complex& val2);
MATCL_MP_EXPORT bool    neq_nan(const Float_complex& val1, const mp_float& val2);
MATCL_MP_EXPORT bool    neq_nan(const mp_float& val1, const mp_int& val2);
MATCL_MP_EXPORT bool    neq_nan(const mp_int& val1, const mp_float& val2);
MATCL_MP_EXPORT bool    neq_nan(const mp_float& val1, const mp_rational& val2);
MATCL_MP_EXPORT bool    neq_nan(const mp_rational& val1, const mp_float& val2);

MATCL_MP_EXPORT bool    neq_nan(const mp_complex& val1, const mp_complex& val2);
MATCL_MP_EXPORT bool    neq_nan(const mp_complex& val1, Integer val2);
MATCL_MP_EXPORT bool    neq_nan(Integer val1, const mp_complex& val2);
MATCL_MP_EXPORT bool    neq_nan(const mp_complex& val1, Float val2);
MATCL_MP_EXPORT bool    neq_nan(const mp_complex& val1, Real val2);
MATCL_MP_EXPORT bool    neq_nan(Float val1, const mp_complex& val2);
MATCL_MP_EXPORT bool    neq_nan(Real val1, const mp_complex& val2);
MATCL_MP_EXPORT bool    neq_nan(const mp_complex& val1, const Complex& val2);
MATCL_MP_EXPORT bool    neq_nan(const Complex& val1, const mp_complex& val2);
MATCL_MP_EXPORT bool    neq_nan(const mp_complex& val1, const Float_complex& val2);
MATCL_MP_EXPORT bool    neq_nan(const Float_complex& val1, const mp_complex& val2);
MATCL_MP_EXPORT bool    neq_nan(const mp_complex& val1, const mp_int& val2);
MATCL_MP_EXPORT bool    neq_nan(const mp_int& val1, const mp_complex& val2);
MATCL_MP_EXPORT bool    neq_nan(const mp_complex& val1, const mp_float& val2);
MATCL_MP_EXPORT bool    neq_nan(const mp_float& val1, const mp_complex& val2);
MATCL_MP_EXPORT bool    neq_nan(const mp_complex& val1, const mp_rational& val2);
MATCL_MP_EXPORT bool    neq_nan(const mp_rational& val1, const mp_complex& val2);

};}}
