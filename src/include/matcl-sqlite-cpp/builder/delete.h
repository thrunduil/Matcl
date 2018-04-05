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
#include "matcl-sqlite-cpp/builder/details/delete_impl.h"
#include "matcl-sqlite-cpp/builder/table.h"
#include "matcl-sqlite-cpp/builder/expression.h"
#include "matcl-sqlite-cpp/builder/details/flags.h"

namespace matcl { namespace sql
{

template<unsigned flags>
class delete_template : public details::choose_base<flags, details::delete_template_base>::result
{
    public:
        delete_template() {}
        delete_template(const details::delete_template_base::impl_ptr& impl)
            : details::choose_base<flags, details::delete_template_base>::result(impl) 
        {}
};

template<unsigned flags, class T>
class delete_template_tab : public T
{
    private:
        using result_type = delete_template<flags & ~details::delete_flags::can_set_table 
                                            | details::delete_flags::can_to_str>;

    public:
        // specify table to delete from.
        result_type tab(const table& arg)
        {
            result_type rv(T::m_impl);
            T::m_impl->tab(arg);
            T::m_impl.reset();
            return rv;
        }

        delete_template_tab() {}
        delete_template_tab(const typename T::impl_ptr& impl)
            : T(impl) 
        {}
};

template<unsigned flags, class T>
class delete_template_where : public T
{
    private:
        using result_type = delete_template<flags & ~details::delete_flags::can_set_condition>;

    public:
        // specify condition for delete.
        result_type where(const expression& arg)
        {
            result_type rv(T::m_impl);
            T::m_impl->where(arg);
            T::m_impl.reset();
            return rv;
        }

        delete_template_where() {}
        delete_template_where(const typename T::impl_ptr& impl)
            : T(impl) 
        {}
};

template<unsigned flags, class T>
class delete_template_to_str : public T
{
    public:
        std::string to_str() const
        {
            return T::m_impl->to_str();
        }

        delete_template_to_str() {}
        delete_template_to_str(const typename T::impl_ptr& impl)
            : T(impl) 
        {}
};

using delete_stmt = delete_template<details::delete_flags::default_flags>;

}}

