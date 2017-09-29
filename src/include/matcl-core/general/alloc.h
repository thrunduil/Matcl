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

#include "matcl-core/matrix/scalar_types.h"
#include "matcl-core/matrix/complex_type.h"

namespace matcl { namespace details
{

// alligned allocation is not strictly necessary in matcl but may
// increase performance (aligning to cache line is recommended in MKL)
struct MATCL_CORE_EXPORT default_allocator
{
    static void*    malloc(size_t n);
    static void*    aligned_malloc(size_t n);

    static void*    realloc(void* ptr, size_t old_size, size_t n);
    static void*    aligned_realloc(void* ptr, size_t old_size, size_t n);

    static void     free(void* ptr, size_t bytes);
    static void     aligned_free(void* ptr, size_t bytes);
};

};};
