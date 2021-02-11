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

#include "matcl-core/config.h"
#include "matcl-core/details/workspace_details.h"

#include <type_traits>

namespace matcl
{

//helper for creation pod types of the same size as given type
template<class T>
struct pod_type
{    
    char    data[sizeof(T)];
};

// pod_workspace allows for fast creation of workspace arrays.
// Objects of this class should be used as local (nonstatic) function
// variables only. This class is not thread safe. One should avoid creating
// many workspaces by pooling then into larger one.
// workspace pointers are aligned at least to MATCL_SIMD_ALIGNMENT bytes.
// workspace is allocated on stack if required size is small; otherwise is
// allocated on heap; released workspace is added to the list of available
// workspaces for future use in order to avoid expensive heap allocations.
template<class T>
class pod_workspace : private details::workspace_base
{
    static_assert(std::is_pod<T>::value, "pod type required!");

    private:
        using base_type = details::workspace_base;

    public:
        // create empty array
        pod_workspace();
  
        // create array of given size
        explicit pod_workspace(size_t size);
        
        // create array of given size and initialize all elements
        // to initial_value
        pod_workspace(size_t size, const T& initial_value);

        // move constructor
        pod_workspace(pod_workspace&& other);

        // number of elements in workspace
        size_t      size() const                    { return details::workspace_base::size_of() / sizeof(T); }

        // return pointer to data
        T*          ptr()                           { return m_ptr; };
        const T*    ptr() const                     { return m_ptr; };

        // get or set value of given element
        T&          operator[](size_t pos)          { return m_ptr[pos]; };
        const T&    operator[](size_t pos) const    { return m_ptr[pos]; };

        // change size to new_size, greater or lower then this->size()
        void        resize(size_t new_size);
        
        // change size to new_size and set value of all elements to val
        void        resize(size_t new_size, const T& val);

    private:
        pod_workspace(const pod_workspace&)             = delete;
        pod_workspace& operator=(const pod_workspace&)  = delete;

        void        set_value(const T& val);

    private:
        T*          m_ptr;
};

};

#include "matcl-core/details/workspace.inl"
