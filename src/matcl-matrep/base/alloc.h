/*
 *  This file is a part of Matrix Computation Library (MATCL)
 *
 *  Copyright (c) Pawe³ Kowal 2017 - 2018
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
#include "matcl-matrep/objects/details/type_info_object.h"
#include "matcl-core/details/integer.h"

namespace matcl { namespace details
{

template<class T>
struct allocator
{
    using alloc = default_allocator<false>;

    static T*   malloc(ti::ti_type<T> ti, size_t n);
    static T*   realloc(ti::ti_type<T> ti, T* old, size_t old_n, size_t new_n);
    static void free(T* ptr, size_t n_elems);
};

template<>
struct allocator<Object>
{
    using T     = Object;

    static T*   malloc(ti::ti_type<T> ti, size_t n);
    static T*   realloc(ti::ti_type<T> ti, T* old, size_t old_n, size_t new_n);
    static void free(T* ptr, size_t n_elems);
};

template<class T>
struct allocator<T*>
{
    using alloc = default_allocator<false>;

    static T**  malloc(size_t n);
    static void free(T** ptr, size_t n_elems);
};

template<class T>
T* allocator<T>::malloc(ti::ti_type<T> ti, size_t n)
{    
    (void)ti;
    void* ptr = alloc::aligned_malloc(imult_size_t(n, sizeof(T)));
    return static_cast<T*>(ptr);
};

template<class T>
T* allocator<T>::realloc(ti::ti_type<T> ti, T* old, size_t old_n, size_t new_n)
{
    (void)ti;
    (void)old_n;
    void* ptr = alloc::aligned_realloc(old, imult_size_t(old_n, sizeof(T)),
                                                   imult_size_t(new_n, sizeof(T)));
    return static_cast<T*>(ptr);
};

template<class T>
void allocator<T>::free(T* ptr, size_t n_elems)
{
    alloc::aligned_free(ptr, n_elems*sizeof(T));
};

template<class T>
T** allocator<T*>::malloc(size_t n)
{
    void* ptr = alloc::malloc(imult_size_t(n, sizeof(T*)));
    return static_cast<T**>(ptr);
};

template<class T>
void allocator<T*>::free(T** ptr, size_t n_elems)
{
    alloc::free(ptr, sizeof(T*) * n_elems);
};

};};
