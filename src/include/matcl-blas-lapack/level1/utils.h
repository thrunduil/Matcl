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

#include "matcl-blas-lapack/level1/config.h"
#include "matcl-blas-lapack/level1/details/simd_details.h"

#include "matcl-simd/simd.h"

namespace matcl { namespace level1
{

//get integer value from template argument, if Val != 0, or from dynamic argument val
// if Val == 0.
template<Integer Val>
struct get_int      { static Integer eval(Integer val) { (void)val; return Val; }; };

template<>
struct get_int<0>   { static Integer eval(Integer val) { return val; }; };

}}
