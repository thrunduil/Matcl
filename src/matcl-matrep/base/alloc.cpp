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

#include "matcl-matrep/base/alloc.h"
#include "matcl-scalar/object.h"

namespace matcl { namespace details
{

Object* allocator<Object>::malloc(ti::ti_object ti, size_t n_elem)
{
    using alloc = default_allocator<true>;
    Object* ptr = (Object*)alloc::aligned_malloc(imult_size_t(n_elem,sizeof(Object)));
    Object* ptr_sav = ptr;
    
    for (size_t i = 0; i < n_elem; ++i, ++ptr)
        new(ptr) Object(ti);

    return ptr_sav;
};

void allocator<Object>::free(Object* ptr, size_t n_elem)
{
    for (size_t i = 0; i < n_elem; ++i)
        ptr[i].~Object();

    using alloc = default_allocator<true>;
    alloc::aligned_free(ptr, n_elem * sizeof(Object));
};

Object* allocator<Object>::realloc(ti::ti_object ti, Object* old_ptr, size_t old_size, 
                                           size_t n_elem)
{
    if (old_ptr == nullptr)
        return allocator<Object>::malloc(ti,n_elem);

    if (n_elem == 0)
    {
        allocator<Object>::free(old_ptr,old_size);
        return nullptr;
    };

    if (n_elem == old_size)
        return old_ptr;

    using alloc = default_allocator<true>;

    if (n_elem < old_size)
    {
        for (size_t i = n_elem; i < old_size; ++i)
            old_ptr[i].~Object();

        Object* new_ptr = (Object*)alloc::aligned_realloc(old_ptr,
                                        imult_size_t(old_size,sizeof(Object)),
                                        imult_size_t(n_elem,sizeof(Object)));
        return new_ptr;
    };

    Object* new_ptr = (Object*)alloc::aligned_realloc(old_ptr,
                                    imult_size_t(old_size,sizeof(Object)), 
                                    imult_size_t(n_elem,sizeof(Object)));

    for (size_t i = old_size; i < n_elem; ++i)
        new(&new_ptr[i]) Object(ti);

    return new_ptr;
};

};};
