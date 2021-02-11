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

#include "matcl-core/matrix/scalar_types.h"
#include "matcl-core/matrix/complex_type.h"

namespace matcl
{

// default allocator;
// if Throw_bad_alloc = true, then error is thrown when malloc fails
// otherwise nullptr is returned
template<bool Throw_bad_alloc>
struct MATCL_CORE_EXPORT default_allocator
{
    using size_type         = std::size_t;
    using difference_type   = std::ptrdiff_t;

    // to use these it is required to track number of allocated bytes
    static void*    malloc(size_t bytes);    
    static void*    realloc(void* ptr, size_t old_size, size_t new_size);
    static void     free(void* ptr, size_t bytes);

    // these functions do not require tracking a number of allocated
    // bytes
    static void*    simple_malloc(size_t bytes);    
    static void*    simple_realloc(void* ptr, size_t new_size);
    static void     simple_free(void* ptr);

    // alligned allocation is not strictly necessary in matcl but may
    // increase performance (aligning to cache line is recommended in MKL)
    static void*    aligned_malloc(size_t n);
    static void*    aligned_realloc(void* ptr, size_t old_size, size_t n);
    static void     aligned_free(void* ptr, size_t bytes);        
};

// default allocator;
// if Throw_bad_alloc = true, then error is thrown when malloc fails
// otherwise nullptr is returned
template<bool Throw_bad_alloc, class Type = void>
struct default_allocator_simple 
    : protected default_allocator<Throw_bad_alloc>
{
    using size_type         = std::size_t;
    using difference_type   = std::ptrdiff_t;

    static Type*    malloc(size_t n_bytes);
    static void     free(void* ptr);
};

template<bool Throw_bad_alloc, class Type>
inline Type* default_allocator_simple<Throw_bad_alloc, Type>
                ::malloc(size_t n_bytes)
{
    using base_type = default_allocator<Throw_bad_alloc>;
    return (Type*)base_type::simple_malloc(n_bytes);
}

template<bool Throw_bad_alloc, class Type>
inline void default_allocator_simple<Throw_bad_alloc, Type>::free(void* ptr)
{
    using base_type = default_allocator<Throw_bad_alloc>;
    return base_type::simple_free(ptr);
}

// defines new and delete operators used for classes derived from
// matcl_new_delete
struct matcl_new_delete
{
    // new operator will call default_allocator<true>::simple_malloc
    void*   operator new(size_t size);

    // new[] operator will call default_allocator<true>::simple_malloc
    void*   operator new[](std::size_t size);

    // delete operator will call default_allocator<true>::simple_free
    void    operator delete(void* ptr);

    // delete operator will call default_allocator<true>::simple_free
    void    operator delete[](void* ptr);
};

};

// call new Ty(args); 
// function default_allocator<true>::simple_malloc will be called
template<class Ty, class ... Args>
Ty*     matcl_new(Args&& ... args);

// call delete ptr, where ptr has type Ty and was created using matcl_new
// function
template<class Ty>
void    matcl_delete(Ty* ptr);

#include "matcl-core/details/alloc.inl"
