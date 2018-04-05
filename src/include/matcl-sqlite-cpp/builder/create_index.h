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

#include <string>
#include "matcl-sqlite-cpp/config.h"
#include "matcl-sqlite-cpp/builder/details/flags.h"
#include "matcl-sqlite-cpp/builder/details/type_builder.h"
#include "matcl-sqlite-cpp/builder/table.h"
#include "matcl-sqlite-cpp/builder/column.h"
#include "matcl-sqlite-cpp/builder/details/create_index_impl.h"

namespace matcl { namespace sql
{

template<unsigned flags>
class create_index_template : public details::choose_base<flags, details::create_index_template_base>::result
{
    public:
        create_index_template() {}
        create_index_template(const details::create_index_template_base::impl_ptr& impl) 
            : details::choose_base<flags, details::create_index_template_base>::result(impl) 
        {}
};

template<unsigned flags, class T>
class create_index_unique : public T
{
    private:
        using result_type   = create_index_template<flags & ~details::create_index_flags::can_set_unique>;

    public:
        // declare index unique.
        result_type unique(bool b=true)
        {
            result_type rv(T::m_impl);
            T::m_impl->unique(b);
            T::m_impl.reset();
            return rv;
        }
            
        create_index_unique() {}
        create_index_unique(const typename T::impl_ptr& impl)
            : T(impl) 
        {}
};

template<unsigned flags, class T>
class create_index_if_not_exists : public T
{
    private:
        using result_type = create_index_template<flags & ~details::create_index_flags::can_set_if_not_exists>;

    public:
        // do nothing if index already exists.
        result_type if_not_exists(bool b=true)
        {
            result_type rv(T::m_impl);
            T::m_impl->if_not_exists(b);
            T::m_impl.reset();
            return rv;
        }
            
        create_index_if_not_exists() {}
        create_index_if_not_exists(const typename T::impl_ptr& impl) 
            : T(impl) 
        {}
};

template<unsigned flags, class T>
class create_index_name  : public T
{
    private:
        using result_type = create_index_template
                            < flags
                            & ~details::create_index_flags::can_set_name
                            | details::create_index_flags::can_set_table
                            >;
    public:
        // set index name.
        result_type name(const std::string& arg)
        {
            result_type rv(T::m_impl);
            T::m_impl->name(arg);
            T::m_impl.reset();
            return rv;
        }

        create_index_name() {}
        create_index_name(const typename T::impl_ptr& impl) 
            : T(impl) 
        {}
};

template<unsigned flags, class T>
class create_index_on : public T
{
    private:
        using result_type = create_index_template
                            < flags
                            & ~details::create_index_flags::can_set_table
                            | details::create_index_flags::can_set_cols
                            >;
    public:
        // specify table to set index on.
        result_type on(const table& arg)
        {
            result_type rv(T::m_impl);
            T::m_impl->on(arg);
            T::m_impl.reset();
            return rv;
        }

        create_index_on() {}
        create_index_on(const typename T::impl_ptr& impl)
            : T(impl) 
        {}
};

template<unsigned flags, class T>
class create_index_cols : public T
{
    private:
        using result_type = create_index_template
                            < flags
                            & ~details::create_index_flags::can_set_cols
                            | details::create_index_flags::can_to_str
                            > ;
    public:
        // specify cols to set index on.
        result_type cols(const column_list& arg)
        {
            result_type rv(T::m_impl);
            T::m_impl->cols(arg);
            T::m_impl.reset();
            return rv;
        }

        create_index_cols() {}
        create_index_cols(const typename T::impl_ptr& impl) 
            : T(impl) 
        {}
};

template<unsigned flags, class T>
class create_index_to_str : public T
{
    public:
        // get SQL.
        std::string to_str() const
        {
            return T::m_impl->to_str();
        }

        create_index_to_str() {}
        create_index_to_str(const typename T::impl_ptr& impl) 
            : T(impl) {}
};

using create_index = create_index_template<details::create_index_flags::default_flags>;

}}
