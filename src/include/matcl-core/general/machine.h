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

//--------------------------------------------------------------------------
//                      Machine dependent parameters
//--------------------------------------------------------------------------
// cache line size in bytes; this is optimization parameters, invalid value
// may have negative impact on performance
#define MATCL_CACHE_LINE_SIZE   64

// set value of this macro to 1 if SSE2 instruction set is available and 0
// otherwise
#define MATCL_ARCHITECTURE_HAS_SSE2 1

// set value of this macro to 1 if SSE3 instruction set is available and 0
// otherwise
#define MATCL_ARCHITECTURE_HAS_SSE3 1

// set value of this macro to 1 if SSE4.1 instruction set is available and 0
// otherwise
#define MATCL_ARCHITECTURE_HAS_SSE41 1

// set value of this macro to 1 if AVX instruction set is available and 0
// otherwise
#define MATCL_ARCHITECTURE_HAS_AVX 1

// set value of this macro to 1 if AVX2 instruction set is available and 0
// otherwise
#define MATCL_ARCHITECTURE_HAS_AVX2 1

// set value of this macro to 1 if FMA instruction set is available and 0
// otherwise
#define MATCL_ARCHITECTURE_HAS_FMA 1

// alignment for SIMD related types; this is optimization parameters, invalid
// value may have negative impact on performance
// TODO
#define MATCL_SIMD_ALIGNMENT    32