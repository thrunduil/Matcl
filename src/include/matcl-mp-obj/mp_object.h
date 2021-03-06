/*
 *  This file is a part of Matrix Computation Library (MATCL)
 *
 *  Copyright (c) Pawe� Kowal 2017 - 2021
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

#include "matcl-mp-obj/config.h"
#include "matcl-mp/mp_int.h"
#include "matcl-mp/mp_float.h"
#include "matcl-mp/mp_complex.h"
#include "matcl-mp/mp_rational.h"
#include "matcl-mp/func_unary.h"
#include "matcl-mp/func_binary.h"
#include "matcl-dynamic/matcl_dynamic.h"
#include "matcl-mp-obj/details/object_type_traits.h"

namespace matcl
{

using MP_int        = dynamic::object_type<mp_int>;
using MP_float      = dynamic::object_type<mp_float>;
using MP_complex    = dynamic::object_type<mp_complex>;
using MP_rational   = dynamic::object_type<mp_rational>;

};
