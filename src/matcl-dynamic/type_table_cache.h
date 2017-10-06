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

#include "type_table_cache_data.h"

#include "matcl-core/details/hash_table/object_table.inl"
#include "matcl-core/details/hash_table/hash_equal.inl"
#include "matcl-core/memory/alloc.h"

namespace matcl { namespace dynamic { namespace details
{

struct last_call_info
{
    int         n_args;
    Type        m_ty_1;
    Type        m_ty_2;
    function    m_func;

    last_call_info();

    last_call_info(int n_args, Type t1, Type t2, function f);
};

class last_call_cache
{
    private:
        using vec_type  = std::vector<last_call_info>;

    private:
        vec_type        m_cache;

    public:
        void            clear();
        function        get_last_function_1(size_t code, Type t1) const;
        function        get_last_function_2(size_t code, Type t1, Type t2) const;

        void            set_last_function(size_t code, int n_args, Type t1, Type t2, 
                            function f);
};

//class is not thread safe
class type_table_cache : public matcl_new_delete, global_object
{
    private:
        using alloc             = matcl::default_allocator_simple<true, char>;
        using unifier_h         = matcl::details::obj_hasher<unifier_impl>;
        using unifier_e         = matcl::details::obj_equaler<unifier_impl>;
        using unifier_ptr       = simple_ptr<unifier_impl>;
        using unifier_table     = matcl::details::object_table<unifier_ptr,unifier_h,unifier_e,alloc>;

        using overload_h        = matcl::details::obj_hasher<overload_impl>;
        using overload_e        = matcl::details::obj_equaler<overload_impl>;
        using overload_ptr      = simple_ptr<overload_impl>;
        using overload_table    = matcl::details::object_table<overload_ptr,overload_h,overload_e,alloc>;

        using toverload_h       = matcl::details::obj_hasher<toverload_impl>;
        using toverload_e       = matcl::details::obj_equaler<toverload_impl>;
        using toverload_ptr     = simple_ptr<toverload_impl>;
        using toverload_table   = matcl::details::object_table<toverload_ptr,toverload_h,toverload_e,alloc>;

        using convert_h         = matcl::details::obj_hasher<convert_impl>;
        using convert_e         = matcl::details::obj_equaler<convert_impl>;
        using convert_ptr       = simple_ptr<convert_impl>;
        using convert_table     = matcl::details::object_table<convert_ptr,convert_h,convert_e,alloc>;

        using assign_h          = matcl::details::obj_hasher<assign_impl>;
        using assign_e          = matcl::details::obj_equaler<assign_impl>;
        using assign_ptr        = simple_ptr<assign_impl>;
        using assign_table      = matcl::details::object_table<assign_ptr,assign_h,assign_e,alloc>;

    private:
        unifier_table   m_unifiers;
        overload_table  m_overloads;
        toverload_table m_template_overloads;
        convert_table   m_convert;
        assign_table    m_assign;
        last_call_cache m_last_call;

    public:
        type_table_cache();

        Type            get_unifier(Type t1, Type t2) const;
        function        get_overload_1(const function_name& func, const Type t[]);
        function        get_overload_2(const function_name& func, const Type t[]);
        function        get_overload_n(const function_name& func, const Type t[], int n);

        function        get_template_overload(const function_name& func, int n_templ, 
                            const Type templates[], int n_args, const Type args[]) const;
        function        get_converter(Type to, Type from, converter_type type) const;
        function        get_assigner(Type to, Type from) const;

        void            set_unifier(Type t1, Type t2, Type t);
        function        set_overload(const function_name& func, const Type t[], int n_args,
                                function f);
        function        set_template_overload(const function_name& func, int n_templ, 
                            const Type templates[], int n_args, const Type t[], function f);
        function        set_converter(Type to, Type from, converter_type type, function f);
        function        set_assigner(Type to, Type from, function f);

        void            clear();
        virtual void    clear_global() override;
        virtual void    close_global() override;

    private:
        type_table_cache(const type_table_cache&) = delete;
        type_table_cache& operator=(const type_table_cache&) = delete;
};

};};};
