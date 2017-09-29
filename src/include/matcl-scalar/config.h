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

// portability issues
#ifdef _MSC_VER
    #define force_inline __forceinline
#else
    #define force_inline inline __attribute__((always_inline))
#endif

// export macros
#ifdef MATCL_SCALAR_EXPORTS
    #define MATCL_SCALAR_EXPORT  __declspec(dllexport)
#else
    #define MATCL_SCALAR_EXPORT  __declspec(dllimport)
#endif

#ifdef __unix__
    #undef  MATCL_SCALAR_EXPORT
    #define MATCL_SCALAR_EXPORT
#endif

