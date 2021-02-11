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

namespace matcl { namespace sql
{

template<unsigned flags>
class update_template;

template<unsigned flags, class T>
class update_template_on_conflict;

template<unsigned flags, class T>
class update_template_tab;

template<unsigned flags, class T>
class update_template_set;

template<unsigned flags, class T>
class update_template_where;

template<unsigned flags, class T>
class update_template_to_str;

}}

namespace matcl { namespace sql { namespace details
{

class update_base_impl
{
    public:
        virtual void on_conflict(conflict_behaviour) = 0;
        virtual void tab(const table&) = 0;
        virtual void set(const assignment_list&) = 0;
        virtual void where(const expression&) = 0;
        virtual std::string to_str() const = 0;
};

class MATCL_SQLITE_EXPORT_MACRO update_template_base
{
    public:
        update_template_base();

    protected:
        using impl_ptr  = std::shared_ptr<update_base_impl>;

        update_template_base(const impl_ptr& impl);
        impl_ptr m_impl;
};

template<unsigned flags_>
struct get_update_to_str
{
    private:
        static const bool enable =
            0 == (flags_ & update_flags::can_set_table) && // table is selected
            0 == (flags_ & update_flags::can_set_assignments); // assignments are defined

    public:
        static const unsigned flags = enable ? (flags_ | update_flags::can_to_str)
                                             : (flags_ & ~update_flags::can_to_str);
};

template<class T, unsigned flags>
struct flag_to_type<update_flags::can_set_conflict_resolver, true, T, flags, update_template_base>
{    
    using result = update_template_on_conflict<flags, T>;
};

template<class T, unsigned flags>
struct flag_to_type<update_flags::can_set_table, true, T, flags, update_template_base>
{    
    using result = update_template_tab<flags, T>;
};

template<class T, unsigned flags>
struct flag_to_type<update_flags::can_set_assignments, true, T, flags, update_template_base>
{    
    using result = update_template_set<flags, T>;
};

template<class T, unsigned flags>
struct flag_to_type<update_flags::can_set_condition, true, T, flags, update_template_base>
{    
    using result = update_template_where<flags, T>;
};

template<class T, unsigned flags>
struct flag_to_type<update_flags::can_to_str, true, T, flags, update_template_base>
{    
    using result = update_template_to_str<flags, T>;
};

}}}