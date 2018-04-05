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

#include "matcl_config.h"

// export macros
#ifdef MATCL_MP_EXPORTS
    #define MATCL_MP_EXPORT  __declspec(dllexport)
#else
    #define MATCL_MP_EXPORT  __declspec(dllimport)
#endif

#ifdef __unix__
    #undef  MATCL_MP_EXPORT
    #define MATCL_MP_EXPORT
#endif

// if this macro is 1, then mp_float class contains a string member
// representing stored value; this allows to inspect mp_float values
// in debugger
#define MATCL_DEBUG_MP_FLOAT 0

// if this macro is 1, then testing mode is enabled (internal use only)
#define MATCL_TEST_MP 0
