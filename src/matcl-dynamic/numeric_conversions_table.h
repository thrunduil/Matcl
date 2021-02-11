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

#include "matcl-dynamic/function_name.h"
#include "matcl-core/details/hash_table/object_table.inl"
#include "matcl-core/details/hash_table/hash_equal.inl"
#include "matcl-core/general/thread.h"
#include "matcl-core/memory/alloc.h"

namespace matcl { namespace dynamic { namespace details
{

struct conv_types
{
    Type to;
    Type from;
};

struct conv_types_data
{
    function        m_func;
    e_match_type    m_match;
};

class numeric_conversion 
{
    private:
        function        m_func;
        e_match_type    m_match;
        size_t          m_hash;

    private:
        numeric_conversion(const numeric_conversion&) = delete;
        numeric_conversion& operator=(const numeric_conversion&) = delete;

        ~numeric_conversion();

        template<class Ptr_type, class Hasher, class Equaler, class Allocator>
        friend class matcl::details::object_table;

    public:
        numeric_conversion(conv_types_data data);        

        function        get_function() const    { return m_func; };
        e_match_type    get_match() const       { return m_match; };
        std::size_t     hash_value() const      { return m_hash; };  
        static size_t   eval_hash(conv_types types);
        static size_t   eval_hash(conv_types_data types);
        bool            equal(conv_types types) const;
        bool            equal(conv_types_data types) const;
};

inline numeric_conversion::numeric_conversion(conv_types_data data)
{
    m_hash      = eval_hash(data);
    m_func      = data.m_func;
    m_match     = data.m_match;
};

inline numeric_conversion::~numeric_conversion()
{
    const_cast<evaler*>(m_func.get_evaler())->destroy();
};

inline size_t numeric_conversion::eval_hash(conv_types types)
{
    size_t seed = types.to.hash_value();
    boost::hash_combine(seed,types.from.hash_value());        
    return seed;
};

inline size_t numeric_conversion::eval_hash(conv_types_data data)
{
    Type to     = data.m_func.return_type();
    Type from   = data.m_func.argument_type(0);
    size_t seed = to.hash_value();
    boost::hash_combine(seed,from.hash_value());        
    return seed;
};

inline bool numeric_conversion::equal(conv_types types) const
{
    Type to     = m_func.return_type();
    Type from   = m_func.argument_type(0);

    if (to != types.to)
        return false;
    if (from != types.from)
        return false;

    return true;
}

inline bool numeric_conversion::equal(conv_types_data data) const
{
    Type to     = m_func.return_type();
    Type from   = m_func.argument_type(0);

    Type to2    = data.m_func.return_type();
    Type from2  = data.m_func.argument_type(0);

    return (to == to2 && from == from2);
};

//---------------------------------------------------------------
//                  numeric conversions
//---------------------------------------------------------------
template<class T>
struct ptr_traits
{
    static void check_free(const T*) {};
    static void copy(const T*) {};
};

struct numeric_conversion_ptr
{
    using ptr_type              = numeric_conversion_ptr;
    using value_type            = numeric_conversion*;
    using obj_type              = numeric_conversion;
    using type_traits           = ptr_traits<numeric_conversion>;

    numeric_conversion*         m_ptr;
    numeric_conversion*         operator->()        { return m_ptr; };
    const numeric_conversion*   operator->() const  { return m_ptr; };

    static ptr_type             make(obj_type* p)       { return ptr_type{p}; };
    static ptr_type             make(const obj_type* p) { return ptr_type{const_cast<obj_type*>(p)}; };
    static ptr_type             make(nullptr_t)         { return ptr_type{nullptr}; };
};

class numeric_conversions_table
{
    private:        
        using string_hash   = matcl::details::obj_hasher<numeric_conversion>;
        using string_equal  = matcl::details::obj_equaler<numeric_conversion>;
        using allocator     = matcl::default_allocator_simple<true, char>;
        using conv_table    = matcl::details::object_table<numeric_conversion_ptr,
                                string_hash,string_equal, allocator>;

    private:
        conv_table          m_table;

        numeric_conversions_table(const numeric_conversions_table&) = delete;
        numeric_conversions_table& operator=(const numeric_conversions_table&) = delete;

    public:
        numeric_conversions_table();

        void                    clear_global();

        void                    insert(function fun_evl, e_match_type match);
        numeric_conversion_ptr  get(Type to, Type from) const;
};

inline numeric_conversions_table::numeric_conversions_table()
    :m_table(true)
{};

inline void numeric_conversions_table::clear_global()
{
    m_table.clear();
}

inline numeric_conversion_ptr numeric_conversions_table::get(Type to, Type from) const
{ 
    return m_table.get_existing(conv_types{to,from}); 
};

inline void numeric_conversions_table::insert(function fun_evl, e_match_type match)
{
    m_table.get(conv_types_data{fun_evl, match}); 
};

};};};

