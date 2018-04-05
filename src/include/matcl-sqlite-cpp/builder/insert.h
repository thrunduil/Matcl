/*
 *  This file is a part of Matrix Computation Library (MATCL)
 *
 *  Copyright (c) Pawe³ Kowal 2018
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

#include "matcl-sqlite-cpp/config.h"
#include "matcl-sqlite-cpp/builder/details/type_builder.h"
#include "matcl-sqlite-cpp/builder/enums.h"
#include "matcl-sqlite-cpp/builder/column.h"
#include "matcl-sqlite-cpp/builder/expression.h"
#include "matcl-sqlite-cpp/builder/details/flags.h"
#include "matcl-sqlite-cpp/builder/details/insert_impl.h"

namespace matcl { namespace sql
{

template<unsigned flags>
class insert_template : public details::choose_base<flags, details::insert_template_base>::result
{
    public:
        insert_template() {}
        insert_template(const details::insert_template_base::impl_ptr& impl)
            : details::choose_base<flags, details::insert_template_base>::result(impl) 
        {}
};

template<unsigned flags, class T>
class insert_template_on_conflict : public T
{
    private:
        using result_type = insert_template<flags & ~ details::insert_flags::can_set_conflict_resolver>;

    public:
        // select conflict behaviour for this insert statement.
        result_type on_conflict(conflict_behaviour cb)
        {
            result_type rv(T::m_impl);
            T::m_impl->on_conflict(cb);
            T::m_impl.reset();
            return rv;
        }

        insert_template_on_conflict() {}
        insert_template_on_conflict(const typename T::impl_ptr& impl)
            : T(impl) 
        {}
};

template<unsigned flags, class T>
class insert_template_columns : public T
{
    private:
        using result_type = insert_template<details::get_insert_to_str<flags 
                                & ~details::insert_flags::can_select_columns>::flags>;

    public:
        // specify columns to insert to.
        result_type columns(const column_list& cols)
        {
            result_type rv(T::m_impl);
            T::m_impl->columns(cols);
            T::m_impl.reset();
            return rv;
        }
            
        insert_template_columns() {}
        insert_template_columns(const typename T::impl_ptr& impl) 
            : T(impl) 
        {}
};

template<unsigned flags, class T>
class insert_template_into : public T
{
    private:
        using result_type = insert_template<details::get_insert_to_str<flags 
                                & ~details::insert_flags::can_select_target>::flags>;

    public:
        // specify table to insert to.
        result_type into(const table& t)
        {
            result_type rv(T::m_impl);
            T::m_impl->into(t);
            T::m_impl.reset();
            return rv;
        }
            
        insert_template_into() {}
        insert_template_into(const typename T::impl_ptr& impl)
            : T(impl) 
        {}
};

template<unsigned flags, class T>
class insert_template_values : public T
{
    private:
        using result_type = insert_template<details::get_insert_to_str<flags 
                                & ~details::insert_flags::can_define_values>::flags>;

    public:
        // specify values to insert from select statement.
        result_type values(const details::select_template_base& sel)
        {
            result_type rv(T::m_impl);
            T::m_impl->values(sel);
            T::m_impl.reset();
            return rv;
        }
            
        // specify values to insert.
        result_type values(const expression_list& el) 
        {
            result_type rv(T::m_impl);
            T::m_impl->values(el);
            T::m_impl.reset();
            return rv;
        }

        insert_template_values() {}
        insert_template_values(const typename T::impl_ptr& impl) 
            : T(impl) 
        {}
};

template<unsigned flags, class T>
class insert_template_to_str : public T
{
    public:
        std::string to_str() const
        {
            return T::m_impl->to_str();
        }

        insert_template_to_str() {}
        insert_template_to_str(const typename T::impl_ptr& impl)
            : T(impl) 
        {}
};

using insert = insert_template<details::insert_flags::default_flags>;

}}

