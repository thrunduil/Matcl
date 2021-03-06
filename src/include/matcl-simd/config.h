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

#include "matcl_config.h"

// export macros
#ifdef MATCL_SIMD_EXPORTS
    #define MATCL_SIMD_EXPORT  __declspec(dllexport)
#else
    #define MATCL_SIMD_EXPORT  __declspec(dllimport)
#endif

#ifdef __unix__
    #undef  MATCL_SIMD_EXPORT
    #define MATCL_SIMD_EXPORT
#endif


// define MATCL_USE_MATCL_COMPLEX in order to use
// nondefault complex types; these types must be defined
// before including simd headers

#ifndef MATCL_USE_MATCL_COMPLEX
    #define MATCL_USE_MATCL_COMPLEX 1
#else
    #define MATCL_USE_MATCL_COMPLEX 0
#endif

// enable this macro in order to define functions generating lookup tables
//#define MATCL_SIMD_GENERATE_TABLES

// machine dependent parameters
#include "matcl-core/general/machine.h"