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

#include <string>
#include <memory>
#include "matcl-sqlite-cpp/config.h"
#include "matcl-sqlite-cpp/builder/details/type_builder.h"
#include "matcl-sqlite-cpp/builder/table.h"
#include "matcl-sqlite-cpp/builder/expression.h"
#include "matcl-sqlite-cpp/builder/column.h"
#include "matcl-sqlite-cpp/builder/ord_expression.h"
#include "matcl-sqlite-cpp/builder/details/flags.h"

namespace matcl { namespace sql
{

template<unsigned flags>
class select_template;

template<unsigned flags, class T>
class select_template_from;

template<unsigned flags, class T>
class select_template_distinct;

template<unsigned flags, class T>
class select_template_columns;

template<unsigned flags, class T>
class select_template_where;

template<unsigned flags, class T>
class select_template_group_by;

template<unsigned flags, class T>
class select_template_order_by;

template<unsigned flags, class T>
class select_template_limit;

template<class T, unsigned flags>
class select_template_composer;

}};

namespace matcl { namespace sql { namespace details
{

class select_base_impl
{
    public:
        virtual std::string to_str() const = 0;

        virtual void from(const join_source& js) = 0;
        virtual void where(const expression& e) = 0;
        virtual void add_column(const expression& e, const std::string& alias) = 0;
        virtual void add_table_columns(const table& t) = 0;
        virtual void distinct(bool b) = 0;
        virtual void group_by(const expression_list& e, const expression &constraint) = 0;
        virtual void order_by(const ord_expression_list& o) = 0;
        virtual void limit(const expression& e, const expression& offset) = 0;

        virtual void compose(const std::shared_ptr<select_base_impl>& impl, compose_type c) = 0;
};

class MATCL_SQLITE_EXPORT_MACRO select_template_base
{
    public:
        select_template_base();

        std::string to_str() const;

    protected:
        using impl_ptr = std::shared_ptr<select_base_impl>;

        select_template_base(const impl_ptr& impl);

        impl_ptr m_impl;
};

template<class T, unsigned flags>
struct flag_to_type<select_flags::can_set_source, true, T, flags, select_template_base>
{    
    using result = select_template_from<flags, T>;
};

template<class T, unsigned flags>
struct flag_to_type<select_flags::can_set_distinct, true, T, flags, select_template_base>
{    
    using result = select_template_distinct<flags, T>;
};

template<class T, unsigned flags>
struct flag_to_type<select_flags::can_set_columns, true, T, flags, select_template_base>
{
    using result = select_template_columns<flags, T>;
};

template<class T, unsigned flags>
struct flag_to_type<select_flags::can_set_where, true, T, flags, select_template_base>
{    
    using result = select_template_where<flags, T>;
};

template<class T, unsigned flags>
struct flag_to_type<select_flags::can_set_grouping, true, T, flags, select_template_base>
{    
    using result = select_template_group_by<flags, T>;
};

template<class T, unsigned flags>
struct flag_to_type<select_flags::can_set_ordering, true, T, flags, select_template_base>
{
    using result = select_template_order_by<flags, T>;
};

template<class T, unsigned flags>
struct flag_to_type<select_flags::can_set_limit, true, T, flags, select_template_base>
{    
    using result = select_template_limit<flags, T>;
};

template<class T, unsigned flags>
struct flag_to_type<select_flags::can_compose, true, T, flags, select_template_base>
{
    using result = select_template_composer<T,flags>;
};

}}}
