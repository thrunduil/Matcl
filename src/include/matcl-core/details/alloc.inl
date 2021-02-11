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

#include "matcl-core/memory/alloc.h"

namespace matcl
{

force_inline 
void* matcl_new_delete::operator new(size_t size)
{
    return default_allocator<true>::simple_malloc(size);
}

force_inline 
void* matcl_new_delete::operator new[](size_t size)
{
    return default_allocator<true>::simple_malloc(size);
}

force_inline
void matcl_new_delete::operator delete(void* ptr)
{
    return default_allocator<true>::simple_free(ptr);
}

force_inline
void matcl_new_delete::operator delete[](void* ptr)
{
    return default_allocator<true>::simple_free(ptr);
}

};

template<class Ty, class ... Args>
force_inline Ty* matcl_new(Args&& ... args)
{
    using allocator = matcl::default_allocator<true>;
    Ty* ptr = (Ty*)allocator::simple_malloc(sizeof(Ty));
    new(ptr) Ty(std::forward<Args>(args)...);

    return ptr;
}

template<class Ty>
force_inline void matcl_delete(Ty* ptr)
{
    using allocator = matcl::default_allocator<true>;

    ptr->~Ty();
    allocator::simple_free(ptr);
    return;
}

