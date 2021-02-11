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

namespace matcl { namespace simd { namespace details
{

#if MATCL_ARCHITECTURE_HAS_FMA

// return x * y + z; use FMA instruction
inline double fma_f(double x, double y, double z);
inline float  fma_f(float x, float y, float z);

// return x * y - z; use FMA instruction
inline double fms_f(double x, double y, double z);
inline float  fms_f(float x, float y, float z);

// return -x * y + z; use FMA instruction
inline double fnma_f(double x, double y, double z);
inline float  fnma_f(float x, float y, float z);

// return -x * y - z; use FMA instruction
inline double fnms_f(double x, double y, double z);
inline float  fnms_f(float x, float y, float z);

#endif

}}}

#include "matcl-core/details/simd/simd.inl"