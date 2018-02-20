/*
 *  This file is a part of Matrix Computation Library (MATCL)
 *
 *  Copyright (c) Pawe� Kowal 2017 - 2018
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

// export macros
#ifdef MATCL_MATREP_EXPORTS
    #define MATCL_MATREP_EXPORT  __declspec(dllexport)
#else
    #define MATCL_MATREP_EXPORT  __declspec(dllimport)
#endif

#ifdef __unix__
    #undef  MATCL_MATREP_EXPORT
    #define MATCL_MATREP_EXPORT
#endif

// internal debug
#define DEBUG_MEMORY 0

// version info
#define MATCL_VERSION 0203
