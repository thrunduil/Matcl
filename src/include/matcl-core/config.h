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

// keyword indicating that a symbol is not aliased in the current scope
#define MATCL_RESTRICTED        __restrict

// export macros
#ifdef MATCL_CORE_EXPORTS
    #define MATCL_CORE_EXPORT  __declspec(dllexport)
#else
    #define MATCL_CORE_EXPORT  __declspec(dllimport)
#endif

#ifdef MATCL_EXPORTS
    #define MATCL_CORE_EXTERN   __declspec(dllexport)
#else
    #define MATCL_CORE_EXTERN   __declspec(dllimport)
#endif

#ifdef __unix__
    #undef  MATCL_CORE_EXPORT
    #define MATCL_CORE_EXPORT

    #undef  MATCL_CORE_EXTERN
    #define MATCL_CORE_EXTERN
#endif

//--------------------------------------------------------------------------
//                      Configuration parameters
//--------------------------------------------------------------------------

// threading
#define MATCL_MULTITHREAD_MODE  0
#define MATCL_SINGLETHREAD_MODE 1

// enable multithreaded mode; compiler must implement "magic statics"
// otherwise it will not be save to use matcl in multithreaded environment
#define MATCL_THREADING_MODE MATCL_MULTITHREAD_MODE
//#define MATCL_THREADING_MODE MATCL_SINGLETHREAD_MODE

// when this macro is set to 1, then dllmalloc memory manager is used
// instead of the default one
#define MATCL_USE_DLMALLOC      1

// when this macro is set to 1, then memory debugging and leak detection
// is enabled
#ifdef _DEBUG
    #define MATCL_DEBUG_MEMORY  1
#else
    #define MATCL_DEBUG_MEMORY  0
#endif

// machine dependent parameters
#include "matcl-core/general/machine.h"