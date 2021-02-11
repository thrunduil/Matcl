/*
 *  This file is a part of Matrix Computation Library (MATCL)
 *
 *  Copyright (c) Pawe³ Kowal 2018 - 2021
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
#include "matcl-sqlite-cpp/builder/column_decl.h"
#include "matcl-sqlite-cpp/builder/table_constraint.h"
#include "matcl-sqlite-cpp/builder/details/type_builder.h"
#include "matcl-sqlite-cpp/builder/details/flags.h"
#include "matcl-sqlite-cpp/builder/details/create_table_impl.h"

namespace matcl { namespace sql
{

template<unsigned flags>
class create_table_template : public details::choose_base<flags, details::create_table_template_base>::result
{
    public:
        create_table_template() {}
        create_table_template(const details::create_table_template_base::impl_ptr& impl) 
            : details::choose_base<flags, details::create_table_template_base>::result(impl) 
        {}
};

template<unsigned flags, class T>
class create_table_temp : public T
{
    private:
        using result_type = create_table_template<flags & ~details::create_table_flags::can_set_temp>;

    public:
        // make table temporary.
        result_type temp(bool b=true)
        {
            result_type rv(T::m_impl);
            T::m_impl->temp(b);
            T::m_impl.reset();
            return rv;
        }

        create_table_temp() {}
        create_table_temp(const typename T::impl_ptr& impl) 
            : T(impl) 
        {}
};

template<unsigned flags, class T>
class create_table_if_not_exists : public T
{
    private:
        using result_type = create_table_template<flags & ~details::create_table_flags::can_set_if_not_exists>;

    public:
        // do nothing if table already exists.
        result_type if_not_exists(bool b=true)
        {
            result_type rv(T::m_impl);
            T::m_impl->if_not_exists(b);
            T::m_impl.reset();
            return rv;
        }
            
        create_table_if_not_exists() {}
        create_table_if_not_exists(const typename T::impl_ptr& impl)
            : T(impl) 
        {}
};

template<unsigned flags, class T>
class create_table_tab : public T
{
    private:
        using result_type = create_table_template
                            < flags
                            & ~details::create_table_flags::can_set_table
                            | details::create_table_flags::can_set_cols
                            | details::create_table_flags::can_set_select
                            >;
    public:
        // name table to create.
        result_type tab(const table& arg)
        {
            result_type rv(T::m_impl);
            T::m_impl->tab(arg);
            T::m_impl.reset();
            return rv;
        }

        create_table_tab() {}
        create_table_tab(const typename T::impl_ptr& impl) 
            : T(impl) 
        {}
};

template<unsigned flags, class T>
class create_table_cols : public T
{
    private:
        using result_type   = create_table_template
                            < flags
                            & ~details::create_table_flags::can_set_cols
                            & ~details::create_table_flags::can_set_select
                            | details::create_table_flags::can_add_constraints
                            | details::create_table_flags::can_to_str
                            >;
    public:
        // define table columns.
        result_type cols(const column_decl_list& arg)
        {
            result_type rv(T::m_impl);
            T::m_impl->cols(arg);
            T::m_impl.reset();
            return rv;
        }

        create_table_cols() {}
        create_table_cols(const typename T::impl_ptr& impl) 
            : T(impl) 
        {}
};

template<unsigned flags, class T>
class create_table_constraints : public T
{
    private:
        using result_type   = create_table_template
                            < flags
                            & ~details::create_table_flags::can_add_constraints
                            >;
    public:
        // define table columns.
        result_type constraints(const table_constraint_list& arg)
        {
            result_type rv(T::m_impl);
            T::m_impl->constraints(arg);
            T::m_impl.reset();
            return rv;
        }

        create_table_constraints() {}
        create_table_constraints(const typename T::impl_ptr& impl)
            : T(impl) 
        {}
};

template<unsigned flags, class T>
class create_table_data : public T
{
    private:
        using result_type = create_table_template
                            < flags
                            & ~details::create_table_flags::can_set_cols
                            & ~details::create_table_flags::can_set_select
                            | details::create_table_flags::can_to_str
                            >;
    public:
        // import table schema and data from select statement.
        result_type data(const details::select_template_base& arg)
        {
            result_type rv(T::m_impl);
            T::m_impl->data(arg);
            T::m_impl.reset();
            return rv;
        }

        create_table_data() {}
        create_table_data(const typename T::impl_ptr& impl) 
            : T(impl) 
        {}
};

template<unsigned flags, class T>
class create_table_to_str : public T
{
    public:
        std::string to_str() const
        {
            return T::m_impl->to_str();
        };

        create_table_to_str() {}
        create_table_to_str(const typename T::impl_ptr& impl)
            : T(impl) 
        {}
};

using create_table = create_table_template<details::create_table_flags::default_flags>;

}};

