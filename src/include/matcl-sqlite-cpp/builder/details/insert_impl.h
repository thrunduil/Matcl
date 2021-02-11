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
#include "matcl-sqlite-cpp/builder/column.h"
#include "matcl-sqlite-cpp/builder/expression.h"
#include "matcl-sqlite-cpp/builder/details/flags.h"

namespace matcl { namespace sql
{

template<unsigned flags>
class insert_template;

template<unsigned flags, class T>
class insert_template_on_conflict;

template<unsigned flags, class T>
class insert_template_columns;

template<unsigned flags, class T>
class insert_template_into;

template<unsigned flags, class T>
class insert_template_values;

template<unsigned flags, class T>
class insert_template_to_str;

}}

namespace matcl { namespace sql { namespace details
{

class insert_base_impl
{
    public:
        virtual void on_conflict(conflict_behaviour) = 0;
        virtual void columns(const column_list&) = 0;
        virtual void into(const table&) = 0;
        virtual void values(const details::select_template_base&) = 0;
        virtual void values(const expression_list& el) = 0;
        virtual std::string to_str() const = 0;
};

class MATCL_SQLITE_EXPORT_MACRO insert_template_base
{
    public:
        insert_template_base();

    protected:
        using impl_ptr = std::shared_ptr<insert_base_impl>;

        insert_template_base(const impl_ptr& impl);

        impl_ptr m_impl;
};

template<unsigned flags_>
struct get_insert_to_str
{
    private:
        static const bool enable =
            0 == (flags_ & insert_flags::can_select_target) && // table is selected
            ((flags_ & insert_flags::can_select_columns) || // columns are not selected
            0 == (flags_ & insert_flags::can_define_values)); // values are specified

    public:
        static const unsigned flags = enable ? (flags_ | insert_flags::can_to_str) 
                                             : (flags_ & ~insert_flags::can_to_str);
};

template<class T, unsigned flags>
struct flag_to_type<insert_flags::can_set_conflict_resolver, true, T, flags, insert_template_base>
{
    using result = insert_template_on_conflict<flags, T>;
};

template<class T, unsigned flags>
struct flag_to_type<insert_flags::can_select_columns, true, T, flags, insert_template_base>
{    
    using result = insert_template_columns<flags, T>;
};

template<class T, unsigned flags>
struct flag_to_type<insert_flags::can_select_target, true, T, flags, insert_template_base>
{    
    using result = insert_template_into<flags, T>;
};

template<class T, unsigned flags>
struct flag_to_type<insert_flags::can_define_values, true, T, flags, insert_template_base>
{
    using result = insert_template_values<flags, T>;
};

template<class T, unsigned flags>
struct flag_to_type<insert_flags::can_to_str, true, T, flags, insert_template_base>
{
    using result = insert_template_to_str<flags, T>;
};

}}}

