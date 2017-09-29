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
#include "matcl-mp/mp_float.h"
#include "matcl-mp/details/utils.h"

namespace matcl { namespace mp { namespace details
{

//---------------------------------------------------------------
//                  operator +
//---------------------------------------------------------------
MATCL_MP_EXPORT mp_int      plus_impl(const mp_int& a, const mp_int& b, precision p);
MATCL_MP_EXPORT mp_int      plus_impl(const mp_int& a, Integer b, precision p);
MATCL_MP_EXPORT mp_int      plus_impl(Integer a, const mp_int& b, precision p);
MATCL_MP_EXPORT mp_float    plus_impl(const mp_int& a, Real b, precision p);
MATCL_MP_EXPORT mp_float    plus_impl(Real a, const mp_int& b, precision p);
MATCL_MP_EXPORT mp_complex  plus_impl(const mp_int& a, const Complex& b, precision p);
MATCL_MP_EXPORT mp_complex  plus_impl(const Complex& a, const mp_int& b, precision p);
MATCL_MP_EXPORT mp_complex  plus_impl(const mp_int& a, const Float_complex& b, precision p);
MATCL_MP_EXPORT mp_complex  plus_impl(const Float_complex& a, const mp_int& b, precision p);

MATCL_MP_EXPORT mp_rational plus_impl(const mp_rational& a, const mp_rational& b, precision p);
MATCL_MP_EXPORT mp_rational plus_impl(const mp_rational& a, Integer b, precision p);
MATCL_MP_EXPORT mp_rational plus_impl(Integer a, const mp_rational& b, precision p);
MATCL_MP_EXPORT mp_float    plus_impl(const mp_rational& a, Real b, precision p);
MATCL_MP_EXPORT mp_float    plus_impl(Real a, const mp_rational& b, precision p);
MATCL_MP_EXPORT mp_complex  plus_impl(const mp_rational& a, const Complex& b, precision p);
MATCL_MP_EXPORT mp_complex  plus_impl(const Complex& a, const mp_rational& b, precision p);
MATCL_MP_EXPORT mp_complex  plus_impl(const mp_rational& a, const Float_complex& b, precision p);
MATCL_MP_EXPORT mp_complex  plus_impl(const Float_complex& a, const mp_rational& b, precision p);
MATCL_MP_EXPORT mp_rational plus_impl(const mp_rational& a, const mp_int& b, precision p);
MATCL_MP_EXPORT mp_rational plus_impl(const mp_int& a, const mp_rational& b, precision p);

MATCL_MP_EXPORT mp_float    plus_impl(const mp_float& a, const mp_float& b, precision p);
MATCL_MP_EXPORT mp_float    plus_impl(const mp_float& a, Integer b, precision p);
MATCL_MP_EXPORT mp_float    plus_impl(Integer a, const mp_float& b, precision p);
MATCL_MP_EXPORT mp_float    plus_impl(const mp_float& a, Real b, precision p);
MATCL_MP_EXPORT mp_float    plus_impl(Real a, const mp_float& b, precision p);
MATCL_MP_EXPORT mp_complex  plus_impl(const mp_float& a, const Complex& b, precision p);
MATCL_MP_EXPORT mp_complex  plus_impl(const Complex& a, const mp_float& b, precision p);
MATCL_MP_EXPORT mp_complex  plus_impl(const mp_float& a, const Float_complex& b, precision p);
MATCL_MP_EXPORT mp_complex  plus_impl(const Float_complex& a, const mp_float& b, precision p);
MATCL_MP_EXPORT mp_float    plus_impl(const mp_float& a, const mp_int& b, precision p);
MATCL_MP_EXPORT mp_float    plus_impl(const mp_int& a, const mp_float& b, precision p);
MATCL_MP_EXPORT mp_float    plus_impl(const mp_float& a, const mp_rational& b, precision p);
MATCL_MP_EXPORT mp_float    plus_impl(const mp_rational& a, const mp_float& b, precision p);

MATCL_MP_EXPORT mp_complex  plus_impl(const mp_complex& a, const mp_complex& b, precision p);
MATCL_MP_EXPORT mp_complex  plus_impl(const mp_complex& a, Integer b, precision p);
MATCL_MP_EXPORT mp_complex  plus_impl(Integer a, const mp_complex& b, precision p);
MATCL_MP_EXPORT mp_complex  plus_impl(const mp_complex& a, Real b, precision p);
MATCL_MP_EXPORT mp_complex  plus_impl(Real a, const mp_complex& b, precision p);
MATCL_MP_EXPORT mp_complex  plus_impl(const mp_complex& a, const Complex& b, precision p);
MATCL_MP_EXPORT mp_complex  plus_impl(const Complex& a, const mp_complex& b, precision p);
MATCL_MP_EXPORT mp_complex  plus_impl(const mp_complex& a, const Float_complex& b, precision p);
MATCL_MP_EXPORT mp_complex  plus_impl(const Float_complex& a, const mp_complex& b, precision p);
MATCL_MP_EXPORT mp_complex  plus_impl(const mp_complex& a, const mp_int& b, precision p);
MATCL_MP_EXPORT mp_complex  plus_impl(const mp_int& a, const mp_complex& b, precision p);
MATCL_MP_EXPORT mp_complex  plus_impl(const mp_complex& a, const mp_float& b, precision p);
MATCL_MP_EXPORT mp_complex  plus_impl(const mp_float& a, const mp_complex& b, precision p);
MATCL_MP_EXPORT mp_complex  plus_impl(const mp_complex& a, const mp_rational& b, precision p);
MATCL_MP_EXPORT mp_complex  plus_impl(const mp_rational& a, const mp_complex& b, precision p);

//---------------------------------------------------------------
//                  operator -
//---------------------------------------------------------------
MATCL_MP_EXPORT mp_int      minus_impl(const mp_int& a, const mp_int& b, precision p);
MATCL_MP_EXPORT mp_int      minus_impl(const mp_int& a, Integer b, precision p);
MATCL_MP_EXPORT mp_int      minus_impl(Integer a, const mp_int& b, precision p);
MATCL_MP_EXPORT mp_float    minus_impl(const mp_int& a, Real b, precision p);
MATCL_MP_EXPORT mp_float    minus_impl(Real a, const mp_int& b, precision p);
MATCL_MP_EXPORT mp_complex  minus_impl(const mp_int& a, const Complex& b, precision p);
MATCL_MP_EXPORT mp_complex  minus_impl(const Complex& a, const mp_int& b, precision p);
MATCL_MP_EXPORT mp_complex  minus_impl(const mp_int& a, const Float_complex& b, precision p);
MATCL_MP_EXPORT mp_complex  minus_impl(const Float_complex& a, const mp_int& b, precision p);

MATCL_MP_EXPORT mp_rational minus_impl(const mp_rational& a, const mp_rational& b, precision p);
MATCL_MP_EXPORT mp_rational minus_impl(const mp_rational& a, Integer b, precision p);
MATCL_MP_EXPORT mp_rational minus_impl(Integer a, const mp_rational& b, precision p);
MATCL_MP_EXPORT mp_float    minus_impl(const mp_rational& a, Real b, precision p);
MATCL_MP_EXPORT mp_float    minus_impl(Real a, const mp_rational& b, precision p);
MATCL_MP_EXPORT mp_complex  minus_impl(const mp_rational& a, const Complex& b, precision p);
MATCL_MP_EXPORT mp_complex  minus_impl(const Complex& a, const mp_rational& b, precision p);
MATCL_MP_EXPORT mp_complex  minus_impl(const mp_rational& a, const Float_complex& b, precision p);
MATCL_MP_EXPORT mp_complex  minus_impl(const Float_complex& a, const mp_rational& b, precision p);
MATCL_MP_EXPORT mp_rational minus_impl(const mp_rational& a, const mp_int& b, precision p);
MATCL_MP_EXPORT mp_rational minus_impl(const mp_int& a, const mp_rational& b, precision p);

MATCL_MP_EXPORT mp_float    minus_impl(const mp_float& a, const mp_float& b, precision p);
MATCL_MP_EXPORT mp_float    minus_impl(const mp_float& a, Integer b, precision p);
MATCL_MP_EXPORT mp_float    minus_impl(Integer a, const mp_float& b, precision p);
MATCL_MP_EXPORT mp_float    minus_impl(const mp_float& a, Real b, precision p);
MATCL_MP_EXPORT mp_float    minus_impl(Real a, const mp_float& b, precision p);
MATCL_MP_EXPORT mp_complex  minus_impl(const mp_float& a, const Complex& b, precision p);
MATCL_MP_EXPORT mp_complex  minus_impl(const Complex& a, const mp_float& b, precision p);
MATCL_MP_EXPORT mp_complex  minus_impl(const mp_float& a, const Float_complex& b, precision p);
MATCL_MP_EXPORT mp_complex  minus_impl(const Float_complex& a, const mp_float& b, precision p);
MATCL_MP_EXPORT mp_float    minus_impl(const mp_float& a, const mp_int& b, precision p);
MATCL_MP_EXPORT mp_float    minus_impl(const mp_int& a, const mp_float& b, precision p);
MATCL_MP_EXPORT mp_float    minus_impl(const mp_float& a, const mp_rational& b, precision p);
MATCL_MP_EXPORT mp_float    minus_impl(const mp_rational& a, const mp_float& b, precision p);

MATCL_MP_EXPORT mp_complex  minus_impl(const mp_complex& a, const mp_complex& b, precision p);
MATCL_MP_EXPORT mp_complex  minus_impl(const mp_complex& a, Integer b, precision p);
MATCL_MP_EXPORT mp_complex  minus_impl(Integer a, const mp_complex& b, precision p);
MATCL_MP_EXPORT mp_complex  minus_impl(const mp_complex& a, Real b, precision p);
MATCL_MP_EXPORT mp_complex  minus_impl(Real a, const mp_complex& b, precision p);
MATCL_MP_EXPORT mp_complex  minus_impl(const mp_complex& a, const Complex& b, precision p);
MATCL_MP_EXPORT mp_complex  minus_impl(const Complex& a, const mp_complex& b, precision p);
MATCL_MP_EXPORT mp_complex  minus_impl(const mp_complex& a, const Float_complex& b, precision p);
MATCL_MP_EXPORT mp_complex  minus_impl(const Float_complex& a, const mp_complex& b, precision p);
MATCL_MP_EXPORT mp_complex  minus_impl(const mp_complex& a, const mp_int& b, precision p);
MATCL_MP_EXPORT mp_complex  minus_impl(const mp_int& a, const mp_complex& b, precision p);
MATCL_MP_EXPORT mp_complex  minus_impl(const mp_complex& a, const mp_float& b, precision p);
MATCL_MP_EXPORT mp_complex  minus_impl(const mp_float& a, const mp_complex& b, precision p);
MATCL_MP_EXPORT mp_complex  minus_impl(const mp_complex& a, const mp_rational& b, precision p);
MATCL_MP_EXPORT mp_complex  minus_impl(const mp_rational& a, const mp_complex& b, precision p);

//---------------------------------------------------------------
//                  operator *
//---------------------------------------------------------------
MATCL_MP_EXPORT mp_int      mul_impl(const mp_int& a, const mp_int& b, precision p);
MATCL_MP_EXPORT mp_int      mul_impl(const mp_int& a, Integer b, precision p);
MATCL_MP_EXPORT mp_int      mul_impl(Integer a, const mp_int& b, precision p);
MATCL_MP_EXPORT mp_float    mul_impl(const mp_int& a, Real b, precision p);
MATCL_MP_EXPORT mp_float    mul_impl(Real a, const mp_int& b, precision p);
MATCL_MP_EXPORT mp_complex  mul_impl(const mp_int& a, const Complex& b, precision p);
MATCL_MP_EXPORT mp_complex  mul_impl(const Complex& a, const mp_int& b, precision p);
MATCL_MP_EXPORT mp_complex  mul_impl(const mp_int& a, const Float_complex& b, precision p);
MATCL_MP_EXPORT mp_complex  mul_impl(const Float_complex& a, const mp_int& b, precision p);

MATCL_MP_EXPORT mp_rational mul_impl(const mp_rational& a, const mp_rational& b, precision p);
MATCL_MP_EXPORT mp_rational mul_impl(const mp_rational& a, Integer b, precision p);
MATCL_MP_EXPORT mp_rational mul_impl(Integer a, const mp_rational& b, precision p);
MATCL_MP_EXPORT mp_float    mul_impl(const mp_rational& a, Real b, precision p);
MATCL_MP_EXPORT mp_float    mul_impl(Real a, const mp_rational& b, precision p);
MATCL_MP_EXPORT mp_complex  mul_impl(const mp_rational& a, const Complex& b, precision p);
MATCL_MP_EXPORT mp_complex  mul_impl(const Complex& a, const mp_rational& b, precision p);
MATCL_MP_EXPORT mp_complex  mul_impl(const mp_rational& a, const Float_complex& b, precision p);
MATCL_MP_EXPORT mp_complex  mul_impl(const Float_complex& a, const mp_rational& b, precision p);
MATCL_MP_EXPORT mp_rational mul_impl(const mp_rational& a, const mp_int& b, precision p);
MATCL_MP_EXPORT mp_rational mul_impl(const mp_int& a, const mp_rational& b, precision p);

MATCL_MP_EXPORT mp_float    mul_impl(const mp_float& a, const mp_float& b, precision p);
MATCL_MP_EXPORT mp_float    mul_impl(const mp_float& a, Integer b, precision p);
MATCL_MP_EXPORT mp_float    mul_impl(Integer a, const mp_float& b, precision p);
MATCL_MP_EXPORT mp_float    mul_impl(const mp_float& a, Real b, precision p);
MATCL_MP_EXPORT mp_float    mul_impl(Real a, const mp_float& b, precision p);
MATCL_MP_EXPORT mp_complex  mul_impl(const mp_float& a, const Complex& b, precision p);
MATCL_MP_EXPORT mp_complex  mul_impl(const Complex& a, const mp_float& b, precision p);
MATCL_MP_EXPORT mp_complex  mul_impl(const mp_float& a, const Float_complex& b, precision p);
MATCL_MP_EXPORT mp_complex  mul_impl(const Float_complex& a, const mp_float& b, precision p);
MATCL_MP_EXPORT mp_float    mul_impl(const mp_float& a, const mp_int& b, precision p);
MATCL_MP_EXPORT mp_float    mul_impl(const mp_int& a, const mp_float& b, precision p);
MATCL_MP_EXPORT mp_float    mul_impl(const mp_float& a, const mp_rational& b, precision p);
MATCL_MP_EXPORT mp_float    mul_impl(const mp_rational& a, const mp_float& b, precision p);

MATCL_MP_EXPORT mp_complex  mul_impl(const mp_complex& a, const mp_complex& b, precision p);
MATCL_MP_EXPORT mp_complex  mul_impl(const mp_complex& a, Integer b, precision p);
MATCL_MP_EXPORT mp_complex  mul_impl(Integer a, const mp_complex& b, precision p);
MATCL_MP_EXPORT mp_complex  mul_impl(const mp_complex& a, Real b, precision p);
MATCL_MP_EXPORT mp_complex  mul_impl(Real a, const mp_complex& b, precision p);
MATCL_MP_EXPORT mp_complex  mul_impl(const mp_complex& a, const Complex& b, precision p);
MATCL_MP_EXPORT mp_complex  mul_impl(const Complex& a, const mp_complex& b, precision p);
MATCL_MP_EXPORT mp_complex  mul_impl(const mp_complex& a, const Float_complex& b, precision p);
MATCL_MP_EXPORT mp_complex  mul_impl(const Float_complex& a, const mp_complex& b, precision p);
MATCL_MP_EXPORT mp_complex  mul_impl(const mp_complex& a, const mp_int& b, precision p);
MATCL_MP_EXPORT mp_complex  mul_impl(const mp_int& a, const mp_complex& b, precision p);
MATCL_MP_EXPORT mp_complex  mul_impl(const mp_complex& a, const mp_float& b, precision p);
MATCL_MP_EXPORT mp_complex  mul_impl(const mp_float& a, const mp_complex& b, precision p);
MATCL_MP_EXPORT mp_complex  mul_impl(const mp_complex& a, const mp_rational& b, precision p);
MATCL_MP_EXPORT mp_complex  mul_impl(const mp_rational& a, const mp_complex& b, precision p);

//---------------------------------------------------------------
//                  operator /
//---------------------------------------------------------------
// integer division gives rational result
// if mp_rational is returned and b is zero, then result is zero
MATCL_MP_EXPORT mp_rational div_impl(const mp_int& a, const mp_int& b, precision p);
MATCL_MP_EXPORT mp_rational div_impl(const mp_int& a, Integer b, precision p);
MATCL_MP_EXPORT mp_rational div_impl(Integer a, const mp_int& b, precision p);
MATCL_MP_EXPORT mp_float    div_impl(const mp_int& a, Real b, precision p);
MATCL_MP_EXPORT mp_float    div_impl(Real a, const mp_int& b, precision p);
MATCL_MP_EXPORT mp_complex  div_impl(const mp_int& a, const Complex& b, precision p);
MATCL_MP_EXPORT mp_complex  div_impl(const Complex& a, const mp_int& b, precision p);
MATCL_MP_EXPORT mp_complex  div_impl(const mp_int& a, const Float_complex& b, precision p);
MATCL_MP_EXPORT mp_complex  div_impl(const Float_complex& a, const mp_int& b, precision p);

MATCL_MP_EXPORT mp_rational div_impl(const mp_rational& a, const mp_rational& b, precision p);
MATCL_MP_EXPORT mp_rational div_impl(const mp_rational& a, Integer b, precision p);
MATCL_MP_EXPORT mp_rational div_impl(Integer a, const mp_rational& b, precision p);
MATCL_MP_EXPORT mp_float    div_impl(const mp_rational& a, Real b, precision p);
MATCL_MP_EXPORT mp_float    div_impl(Real a, const mp_rational& b, precision p);
MATCL_MP_EXPORT mp_complex  div_impl(const mp_rational& a, const Complex& b, precision p);
MATCL_MP_EXPORT mp_complex  div_impl(const Complex& a, const mp_rational& b, precision p);
MATCL_MP_EXPORT mp_complex  div_impl(const mp_rational& a, const Float_complex& b, precision p);
MATCL_MP_EXPORT mp_complex  div_impl(const Float_complex& a, const mp_rational& b, precision p);
MATCL_MP_EXPORT mp_rational div_impl(const mp_rational& a, const mp_int& b, precision p);
MATCL_MP_EXPORT mp_rational div_impl(const mp_int& a, const mp_rational& b, precision p);

MATCL_MP_EXPORT mp_float    div_impl(const mp_float& a, const mp_float& b, precision p);
MATCL_MP_EXPORT mp_float    div_impl(const mp_float& a, Integer b, precision p);
MATCL_MP_EXPORT mp_float    div_impl(Integer a, const mp_float& b, precision p);
MATCL_MP_EXPORT mp_float    div_impl(const mp_float& a, Real b, precision p);
MATCL_MP_EXPORT mp_float    div_impl(Real a, const mp_float& b, precision p);
MATCL_MP_EXPORT mp_complex  div_impl(const mp_float& a, const Complex& b, precision p);
MATCL_MP_EXPORT mp_complex  div_impl(const Complex& a, const mp_float& b, precision p);
MATCL_MP_EXPORT mp_complex  div_impl(const mp_float& a, const Float_complex& b, precision p);
MATCL_MP_EXPORT mp_complex  div_impl(const Float_complex& a, const mp_float& b, precision p);
MATCL_MP_EXPORT mp_float    div_impl(const mp_float& a, const mp_int& b, precision p);
MATCL_MP_EXPORT mp_float    div_impl(const mp_int& a, const mp_float& b, precision p);
MATCL_MP_EXPORT mp_float    div_impl(const mp_float& a, const mp_rational& b, precision p);
MATCL_MP_EXPORT mp_float    div_impl(const mp_rational& a, const mp_float& b, precision p);

MATCL_MP_EXPORT mp_complex  div_impl(const mp_complex& a, const mp_complex& b, precision p);
MATCL_MP_EXPORT mp_complex  div_impl(const mp_complex& a, Integer b, precision p);
MATCL_MP_EXPORT mp_complex  div_impl(Integer a, const mp_complex& b, precision p);
MATCL_MP_EXPORT mp_complex  div_impl(const mp_complex& a, Real b, precision p);
MATCL_MP_EXPORT mp_complex  div_impl(Real a, const mp_complex& b, precision p);
MATCL_MP_EXPORT mp_complex  div_impl(const mp_complex& a, const Complex& b, precision p);
MATCL_MP_EXPORT mp_complex  div_impl(const Complex& a, const mp_complex& b, precision p);
MATCL_MP_EXPORT mp_complex  div_impl(const mp_complex& a, const Float_complex& b, precision p);
MATCL_MP_EXPORT mp_complex  div_impl(const Float_complex& a, const mp_complex& b, precision p);
MATCL_MP_EXPORT mp_complex  div_impl(const mp_complex& a, const mp_int& b, precision p);
MATCL_MP_EXPORT mp_complex  div_impl(const mp_int& a, const mp_complex& b, precision p);
MATCL_MP_EXPORT mp_complex  div_impl(const mp_complex& a, const mp_float& b, precision p);
MATCL_MP_EXPORT mp_complex  div_impl(const mp_float& a, const mp_complex& b, precision p);
MATCL_MP_EXPORT mp_complex  div_impl(const mp_complex& a, const mp_rational& b, precision p);
MATCL_MP_EXPORT mp_complex  div_impl(const mp_rational& a, const mp_complex& b, precision p);

MATCL_MP_EXPORT mp_complex  inv_impl(const mp_complex& a, precision p);

//---------------------------------------------------------------
//                  idiv_impl
//---------------------------------------------------------------
// integer division gives integer result
// if mp_rational or mp_int is returned and b is zero, then result
// is zero
MATCL_MP_EXPORT mp_int      idiv_impl(const mp_int& a, const mp_int& b, precision p);
MATCL_MP_EXPORT mp_int      idiv_impl(const mp_int& a, Integer b, precision p);
MATCL_MP_EXPORT mp_int      idiv_impl(Integer a, const mp_int& b, precision p);
MATCL_MP_EXPORT mp_float    idiv_impl(const mp_int& a, Real b, precision p);
MATCL_MP_EXPORT mp_float    idiv_impl(Real a, const mp_int& b, precision p);
MATCL_MP_EXPORT mp_complex  idiv_impl(const mp_int& a, const Complex& b, precision p);
MATCL_MP_EXPORT mp_complex  idiv_impl(const Complex& a, const mp_int& b, precision p);
MATCL_MP_EXPORT mp_complex  idiv_impl(const mp_int& a, const Float_complex& b, precision p);
MATCL_MP_EXPORT mp_complex  idiv_impl(const Float_complex& a, const mp_int& b, precision p);

MATCL_MP_EXPORT mp_rational idiv_impl(const mp_rational& a, const mp_rational& b, precision p);
MATCL_MP_EXPORT mp_rational idiv_impl(const mp_rational& a, Integer b, precision p);
MATCL_MP_EXPORT mp_rational idiv_impl(Integer a, const mp_rational& b, precision p);
MATCL_MP_EXPORT mp_float    idiv_impl(const mp_rational& a, Real b, precision p);
MATCL_MP_EXPORT mp_float    idiv_impl(Real a, const mp_rational& b, precision p);
MATCL_MP_EXPORT mp_complex  idiv_impl(const mp_rational& a, const Complex& b, precision p);
MATCL_MP_EXPORT mp_complex  idiv_impl(const Complex& a, const mp_rational& b, precision p);
MATCL_MP_EXPORT mp_complex  idiv_impl(const mp_rational& a, const Float_complex& b, precision p);
MATCL_MP_EXPORT mp_complex  idiv_impl(const Float_complex& a, const mp_rational& b, precision p);
MATCL_MP_EXPORT mp_rational idiv_impl(const mp_rational& a, const mp_int& b, precision p);
MATCL_MP_EXPORT mp_rational idiv_impl(const mp_int& a, const mp_rational& b, precision p);

MATCL_MP_EXPORT mp_float    idiv_impl(const mp_float& a, const mp_float& b, precision p);
MATCL_MP_EXPORT mp_float    idiv_impl(const mp_float& a, Integer b, precision p);
MATCL_MP_EXPORT mp_float    idiv_impl(Integer a, const mp_float& b, precision p);
MATCL_MP_EXPORT mp_float    idiv_impl(const mp_float& a, Real b, precision p);
MATCL_MP_EXPORT mp_float    idiv_impl(Real a, const mp_float& b, precision p);
MATCL_MP_EXPORT mp_complex  idiv_impl(const mp_float& a, const Complex& b, precision p);
MATCL_MP_EXPORT mp_complex  idiv_impl(const Complex& a, const mp_float& b, precision p);
MATCL_MP_EXPORT mp_complex  idiv_impl(const mp_float& a, const Float_complex& b, precision p);
MATCL_MP_EXPORT mp_complex  idiv_impl(const Float_complex& a, const mp_float& b, precision p);
MATCL_MP_EXPORT mp_float    idiv_impl(const mp_float& a, const mp_int& b, precision p);
MATCL_MP_EXPORT mp_float    idiv_impl(const mp_int& a, const mp_float& b, precision p);
MATCL_MP_EXPORT mp_float    idiv_impl(const mp_float& a, const mp_rational& b, precision p);
MATCL_MP_EXPORT mp_float    idiv_impl(const mp_rational& a, const mp_float& b, precision p);

MATCL_MP_EXPORT mp_complex  idiv_impl(const mp_complex& a, const mp_complex& b, precision p);
MATCL_MP_EXPORT mp_complex  idiv_impl(const mp_complex& a, Integer b, precision p);
MATCL_MP_EXPORT mp_complex  idiv_impl(Integer a, const mp_complex& b, precision p);
MATCL_MP_EXPORT mp_complex  idiv_impl(const mp_complex& a, Real b, precision p);
MATCL_MP_EXPORT mp_complex  idiv_impl(Real a, const mp_complex& b, precision p);
MATCL_MP_EXPORT mp_complex  idiv_impl(const mp_complex& a, const Complex& b, precision p);
MATCL_MP_EXPORT mp_complex  idiv_impl(const Complex& a, const mp_complex& b, precision p);
MATCL_MP_EXPORT mp_complex  idiv_impl(const mp_complex& a, const Float_complex& b, precision p);
MATCL_MP_EXPORT mp_complex  idiv_impl(const Float_complex& a, const mp_complex& b, precision p);
MATCL_MP_EXPORT mp_complex  idiv_impl(const mp_complex& a, const mp_int& b, precision p);
MATCL_MP_EXPORT mp_complex  idiv_impl(const mp_int& a, const mp_complex& b, precision p);
MATCL_MP_EXPORT mp_complex  idiv_impl(const mp_complex& a, const mp_float& b, precision p);
MATCL_MP_EXPORT mp_complex  idiv_impl(const mp_float& a, const mp_complex& b, precision p);
MATCL_MP_EXPORT mp_complex  idiv_impl(const mp_complex& a, const mp_rational& b, precision p);
MATCL_MP_EXPORT mp_complex  idiv_impl(const mp_rational& a, const mp_complex& b, precision p);

//---------------------------------------------------------------
//                  pow_impl
//---------------------------------------------------------------
MATCL_MP_EXPORT mp_float    pow_impl(const mp_int& a, const mp_int& b, precision p);
MATCL_MP_EXPORT mp_float    pow_impl(const mp_int& a, Integer b, precision p);
MATCL_MP_EXPORT mp_float    pow_impl(Integer a, const mp_int& b, precision p);
MATCL_MP_EXPORT mp_float    pow_impl(const mp_rational& a, Integer b, precision p);
MATCL_MP_EXPORT mp_float    pow_impl(const mp_rational& a, const mp_int& b, precision p);

MATCL_MP_EXPORT mp_float    pow_impl(const mp_rational& a, const mp_rational& b, precision p);
MATCL_MP_EXPORT mp_float    pow_impl(Integer a, const mp_rational& b, precision p);
MATCL_MP_EXPORT mp_float    pow_impl(const mp_int& a, const mp_rational& b, precision p);

MATCL_MP_EXPORT mp_float    pow_impl(const mp_float& a, const mp_float& b, precision p);
MATCL_MP_EXPORT mp_float    pow_impl(const mp_float& a, Integer b, precision p);
MATCL_MP_EXPORT mp_float    pow_impl(Integer a, const mp_float& b, precision p);
MATCL_MP_EXPORT mp_float    pow_impl(const mp_float& a, const mp_int& b, precision p);
MATCL_MP_EXPORT mp_float    pow_impl(const mp_int& a, const mp_float& b, precision p);
MATCL_MP_EXPORT mp_float    pow_impl(const mp_float& a, const mp_rational& b, precision p);
MATCL_MP_EXPORT mp_float    pow_impl(const mp_rational& a, const mp_float& b, precision p);

MATCL_MP_EXPORT mp_complex  pow_impl(const mp_complex& a, const mp_complex& b, precision p);
MATCL_MP_EXPORT mp_complex  pow_impl(const mp_complex& a, Integer b, precision p);
MATCL_MP_EXPORT mp_complex  pow_impl(Integer a, const mp_complex& b, precision p);
MATCL_MP_EXPORT mp_complex  pow_impl(const mp_complex& a, const mp_int& b, precision p);
MATCL_MP_EXPORT mp_complex  pow_impl(const mp_int& a, const mp_complex& b, precision p);
MATCL_MP_EXPORT mp_complex  pow_impl(const mp_complex& a, const mp_float& b, precision p);
MATCL_MP_EXPORT mp_complex  pow_impl(const mp_float& a, const mp_complex& b, precision p);
MATCL_MP_EXPORT mp_complex  pow_impl(const mp_complex& a, const mp_rational& b, precision p);
MATCL_MP_EXPORT mp_complex  pow_impl(const mp_rational& a, const mp_complex& b, precision p);

template<class T1, class T2>
struct result_of_pow
{
    using T1P   = typename promote_floats<T1>::type;
    using T2P   = typename promote_floats<T2>::type;
    using type  = decltype(pow_impl(T1P(),T2P(),precision()));
};

};}}