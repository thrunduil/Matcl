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
#include "matcl-sqlite-cpp/builder/details/type_builder.h"
#include "matcl-sqlite-cpp/builder/enums.h"
#include "matcl-sqlite-cpp/builder/table.h"
#include "matcl-sqlite-cpp/builder/assignment.h"
#include "matcl-sqlite-cpp/builder/details/flags.h"
#include "matcl-sqlite-cpp/builder/details/update_impl.h"

namespace matcl { namespace sql
{

template<unsigned flags>
class update_template : public details::choose_base<flags, details::update_template_base>::result
{
    public:
        update_template() {}
        update_template(const details::update_template_base::impl_ptr& impl) 
            : details::choose_base<flags, update_template_base>::result(impl) 
        {}
};

template<unsigned flags, class T>
class update_template_on_conflict : public T
{
    private:
        using result_type = update_template<flags & ~details::update_flags::can_set_conflict_resolver>;

    public:
        //select conflict behaviour for this update statement.
        result_type on_conflict(conflict_behaviour arg)
        {
            result_type rv(T::m_impl);
            T::m_impl->on_conflict(arg);
            T::m_impl.reset();
            return rv;
        }

        update_template_on_conflict() {}
        update_template_on_conflict(const typename T::impl_ptr& impl)
            : T(impl) 
        {}
};

template<unsigned flags, class T>
class update_template_tab : public T
{
    private:
        using result_type = update_template<details::get_update_to_str<flags 
                                & ~details::update_flags::can_set_table>::flags>;

    public:
        // define table to update.
        result_type tab(const table& arg)
        {
            result_type rv(T::m_impl);
            T::m_impl->tab(arg);
            T::m_impl.reset();
            return rv;
        }

        update_template_tab() {}
        update_template_tab(const typename T::impl_ptr& impl)
            : T(impl) 
        {}
};

template<unsigned flags, class T>
class update_template_set : public T
{
    private:
        using result_type = update_template<details::get_update_to_str<flags 
                                & ~details::update_flags::can_set_assignments>::flags>;

    public:
        // define columns to update.
        result_type set(const assignment_list& arg)
        {
            result_type rv(T::m_impl);
            T::m_impl->set(arg);
            T::m_impl.reset();
            return rv;
        }

        update_template_set() {}
        update_template_set(const typename T::impl_ptr& impl)
            : T(impl) 
        {}
};

template<unsigned flags, class T>
class update_template_where : public T
{
    private:
        using result_type = update_template<flags & ~details::update_flags::can_set_condition>;

    public:
        // define update condition.
        result_type where(const expression& arg)
        {
            result_type rv(T::m_impl);
            T::m_impl->where(arg);
            T::m_impl.reset();
            return rv;
        }

        update_template_where() {}
        update_template_where(const typename T::impl_ptr& impl)
            : T(impl) 
        {}
};

template<unsigned flags, class T>
class update_template_to_str : public T
{
    public:
        std::string to_str() const
        {
            return T::m_impl->to_str();
        }

        update_template_to_str() {}
        update_template_to_str(const typename T::impl_ptr& impl)
            : T(impl) 
        {}
};

using update = update_template<details::update_flags::default_flags>;

}}
