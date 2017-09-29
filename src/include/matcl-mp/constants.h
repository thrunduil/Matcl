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

#include "matcl-mp/mp_float.h"

namespace matcl { namespace constants
{

// constants can be reevaluated at each function call

//-----------------------------------------------------------------
//                      float constants
//-----------------------------------------------------------------
// estimate machine epsilon for the given precision;
// return smallest eps such that 1.0 + eps != 1.0
MATCL_MP_EXPORT mp_float    mp_eps(precision prec = precision(0));

// infinite value
MATCL_MP_EXPORT mp_float    mp_inf(precision prec = precision(0));

// not-a-number value
MATCL_MP_EXPORT mp_float    mp_nan(precision prec = precision(0));

// e value; e = exp(1)
MATCL_MP_EXPORT mp_float    mp_e(precision prec = precision(0));

// pi value
MATCL_MP_EXPORT mp_float    mp_pi(precision prec = precision(0));

// pi/2 value
MATCL_MP_EXPORT mp_float    mp_pi_2(precision prec = precision(0));

// pi/4 value
MATCL_MP_EXPORT mp_float    mp_pi_4(precision prec = precision(0));

// the logarithm of 2
MATCL_MP_EXPORT mp_float    mp_ln2(precision prec = precision(0));

// the logarithm of 2
MATCL_MP_EXPORT mp_float    mp_ln10(precision prec = precision(0));

// the logarithm to base 2 of e
MATCL_MP_EXPORT mp_float    mp_log2e(precision prec = precision(0));

// the logarithm to base 10 of e
MATCL_MP_EXPORT mp_float    mp_log10e(precision prec = precision(0));

//-----------------------------------------------------------------
//                      complex constants
//-----------------------------------------------------------------
// imaginary unit
MATCL_MP_EXPORT mp_complex  mp_i(precision prec = precision(0));

};};