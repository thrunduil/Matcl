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

#pragma once

#include "matcl-core/config.h"
#include <stdint.h>

namespace matcl
{

template<class T> struct complex;

// signed integer
using Integer       = int;

// single precision floating point number
using Float         = float;

// doule precision floating point number
using Real          = double;

// single precision complex number
using Float_complex = matcl::complex<Float>;

// double precision complex number
using Complex       = matcl::complex<Real>;

// 64-bit signed integer
using Integer_64    = int64_t;

};
