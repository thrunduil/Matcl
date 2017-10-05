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

//---------------------------------------------------------------
//                  null conversions
//---------------------------------------------------------------

class null_conversion 
{
    private:
        function        m_func;
        Type            m_to;

        null_conversion(Type to);
        ~null_conversion();

        null_conversion(const null_conversion&) = delete;
        null_conversion& operator=(const null_conversion&) = delete;

        template<class Ptr_type, class Hasher, class Equaler, class Allocator>
        friend class matcl::details::object_table;

    public:        

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

class null_conversions_table
{
    private:        
        using string_hash   = matcl::details::obj_hasher<null_conversion>;
        using string_equal  = matcl::details::obj_equaler<null_conversion>;
        using allocator     = matcl::default_allocator_simple<true, char>;
        using conv_table    = matcl::details::object_table<null_conversion_ptr,
                                string_hash,string_equal, allocator>;
        using mutex         = matcl::default_spinlock_mutex;

    private:
        mutable conv_table  m_table;
        mutable mutex       m_mutex;

        null_conversions_table(const null_conversions_table&) = delete;
        null_conversions_table& operator=(const null_conversions_table&) = delete;

    public:
        null_conversions_table();

        void                clear_global();

        function            get(Type to) const;
};

inline null_conversions_table::null_conversions_table()
    :m_table(true)
{};

inline void null_conversions_table::clear_global()
{
    m_table.clear();
}

inline function null_conversions_table::get(Type to) const
{
    //table can be modified and this function can be called from multiple threads
    std::unique_lock<mutex> lock(m_mutex);
    return m_table.get(to)->get_function(); 
};

//---------------------------------------------------------------
//                  link conversions
//---------------------------------------------------------------
class fun_conv_link;

class link_conversion 
{    
    private:
        fun_conv_link*  m_func;
        size_t          m_hash;

    private:
        link_conversion(const std::vector<function>& conv);
        ~link_conversion();

        link_conversion(const link_conversion&) = delete;
        link_conversion& operator=(const link_conversion&) = delete;

        template<class Ptr_type, class Hasher, class Equaler, class Allocator>
        friend class matcl::details::object_table;

    public:
        function        get_function() const;
        std::size_t     hash_value() const      { return m_hash; };

        static size_t   eval_hash(const std::vector<function>& conv);
        bool            equal(const std::vector<function>& conv) const;
};

struct link_conversion_ptr
{
    using ptr_type              = link_conversion_ptr;
    using value_type            = link_conversion*;
    using obj_type              = link_conversion;
    using type_traits           = ptr_traits<link_conversion>;

    link_conversion*            m_ptr;

    link_conversion*            operator->()        { return m_ptr; };
    const link_conversion*      operator->() const  { return m_ptr; };

    static ptr_type             make(obj_type* p)       { return ptr_type{p}; };
    static ptr_type             make(const obj_type* p) { return ptr_type{const_cast<obj_type*>(p)}; };
    static ptr_type             make(nullptr_t)         { return ptr_type{nullptr}; };
};

class link_conversions_table
{
    private:        
        using string_hash   = matcl::details::obj_hasher<link_conversion>;
        using string_equal  = matcl::details::obj_equaler<link_conversion>;
        using allocator     = matcl::default_allocator_simple<true, char>;
        using conv_table    = matcl::details::object_table<link_conversion_ptr,
                                string_hash,string_equal, allocator>;
        using mutex         = matcl::default_spinlock_mutex;

    private:
        mutable conv_table  m_table;
        mutable mutex       m_mutex;

        link_conversions_table(const link_conversions_table&) = delete;
        link_conversions_table& operator=(const link_conversions_table&) = delete;

    public:
        link_conversions_table();

        void                clear_global();
        function            get(const std::vector<function>& conv) const;
};

inline link_conversions_table::link_conversions_table()
    :m_table(true)
{};

inline void link_conversions_table::clear_global()
{
    m_table.clear();
}

inline function link_conversions_table::get(const std::vector<function>& conv) const
{
    //table can be modified and this function can be called from multiple threads
    std::unique_lock<mutex> lock(m_mutex);
    return m_table.get(conv)->get_function(); 
};

//---------------------------------------------------------------
//                  conversions
//---------------------------------------------------------------

struct build_convert_info
{
    function    m_base;
    int         n_deduced;
    const Type* m_deduced;
    Type        m_deduced_ret;

    const std::vector<function>* 
                m_arg_converters;

    mutable 
    size_t      m_hash;
};

class fun_evaler_conv;

class conversion 
{    
    private:
        size_t              m_hash;
        fun_evaler_conv*    m_func;

    private:
        conversion(const build_convert_info& info);
        ~conversion();

        template<class Ptr_type, class Hasher, class Equaler, class Allocator>
        friend class matcl::details::object_table;

        conversion(const conversion&) = delete;
        conversion& operator=(const conversion&) = delete;

    public:
        function        get_function() const;
        std::size_t     hash_value() const      { return m_hash; };

        static size_t   eval_hash(const build_convert_info& info);
        bool            equal(const build_convert_info& info) const;
};

struct conversion_ptr
{
    using ptr_type              = conversion_ptr;
    using value_type            = conversion*;
    using obj_type              = conversion;
    using type_traits           = ptr_traits<conversion>;

    conversion*                 m_ptr;

    conversion*                 operator->()        { return m_ptr; };
    const conversion*           operator->() const  { return m_ptr; };

    static ptr_type             make(obj_type* p)       { return ptr_type{p}; };
    static ptr_type             make(const obj_type* p) { return ptr_type{const_cast<obj_type*>(p)}; };
    static ptr_type             make(nullptr_t)         { return ptr_type{nullptr}; };
};

class conversions_table
{
    private:        
        using string_hash   = matcl::details::obj_hasher<conversion>;
        using string_equal  = matcl::details::obj_equaler<conversion>;
        using allocator     = matcl::default_allocator_simple<true, char>;
        using conv_table    = matcl::details::object_table<conversion_ptr,
                                string_hash,string_equal, allocator>;
        using mutex         = matcl::default_spinlock_mutex;

    private:
        mutable conv_table  m_table;
        mutable mutex       m_mutex;

        conversions_table(const conversions_table&) = delete;
        conversions_table& operator=(const conversions_table&) = delete;

    public:
        conversions_table();

        void        clear_global();

        // make function with conversions; returning function first evaluates
        // converters for each input argument and then this function;
        // return type of n-th converter must be the same as n-th argument
        // type of this function; deduced_ret is the deduced return type
        // and is different from Type() if return deduction is enabled;
        // n_deduced is the number of deduced template types with deduced 
        // types given by deduced; deduced return (if not empty) and deduced
        // template arguments must be passed as first m + n_deduced arguments
        // of type OType, where m is 1 if deduced_ret is not Type() and zero
        // otherwise; deduced return must be passed first
        function        get(function base, int n_deduced, const Type deduced[],
                            Type deduced_ret, const std::vector<function>& arg_converters) const;
};

inline conversions_table::conversions_table()
    :m_table(true)
{};

};};};

