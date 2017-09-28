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

#include "matcl-dynamic/function_name.h"
#include "matcl-core/details/hash_table/object_table.inl"
#include "matcl-core/details/hash_table/hash_equal.inl"
#include "matcl-core/general/thread.h"
#include "matcl-dynamic/function.h"
#include "numeric_conversions_table.h"

namespace matcl { namespace dynamic { namespace details
{

class null_conversion 
{
    private:
        function        m_func;
        Type            m_to;

    public:
        null_conversion(Type to);

        function        get_function() const    { return m_func; };
        std::size_t     hash_value() const      { return m_to.hash_value(); };

        static size_t   eval_hash(Type to)      { return to.hash_value(); };
        bool            equal(Type to) const    { return to == this->m_to; };
};

struct null_conversion_ptr
{
    using ptr_type              = null_conversion_ptr;
    using value_type            = null_conversion*;
    using obj_type              = null_conversion;
    using type_traits           = ptr_traits<null_conversion>;

    null_conversion*            m_ptr;

    null_conversion*            operator->()        { return m_ptr; };
    const null_conversion*      operator->() const  { return m_ptr; };

    static ptr_type             make(obj_type* p)       { return ptr_type{p}; };
    static ptr_type             make(const obj_type* p) { return ptr_type{const_cast<obj_type*>(p)}; };
    static ptr_type             make(nullptr_t)         { return ptr_type{nullptr}; };
};

//---------------------------------------------------------------
//                  null conversions
//---------------------------------------------------------------

class null_conversions_table
{
    private:        
        using string_hash   = matcl::details::obj_hasher<null_conversion>;
        using string_equal  = matcl::details::obj_equaler<null_conversion>;
        using conv_table    = matcl::details::object_table<null_conversion_ptr,
                                string_hash,string_equal, matcl::details::default_allocator>;
        using mutex         = matcl::default_spinlock_mutex;

    private:
        mutable conv_table  m_table;
        mutable mutex       m_mutex;

    public:
        function            get(Type to) const;
};

inline function null_conversions_table::get(Type to) const
{
    //table can be modified and this function can be called from multiple threads
    std::unique_lock<mutex> lock(m_mutex);
    return m_table.get(to)->get_function(); 
};

};};};

