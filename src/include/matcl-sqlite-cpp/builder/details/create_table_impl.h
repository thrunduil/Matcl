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

namespace matcl { namespace sql
{

template<unsigned flags>
class create_table_template;

template<unsigned flags, class T>
class create_table_temp;

template<unsigned flags, class T>
class create_table_if_not_exists;

template<unsigned flags, class T>
class create_table_tab;

template<unsigned flags, class T>
class create_table_cols;

template<unsigned flags, class T>
class create_table_constraints;

template<unsigned flags, class T>
class create_table_data;

template<unsigned flags, class T>
class create_table_to_str;

}}

namespace matcl { namespace sql { namespace details
{

class create_table_base_impl
{
    public:
        virtual std::string to_str() const = 0;

        virtual void temp(bool b) = 0;
        virtual void if_not_exists(bool b) = 0;
        virtual void tab(const table& t) = 0;
        virtual void cols(const column_decl_list& cdl) = 0;
        virtual void data(const details::select_template_base& sel) = 0;
        virtual void constraints(const table_constraint_list& sel) = 0;
};

class MATCL_SQLITE_EXPORT_MACRO create_table_template_base
{
    public:
        create_table_template_base();

    protected:
        using impl_ptr  = std::shared_ptr<create_table_base_impl>;

        create_table_template_base(const impl_ptr& impl);

        impl_ptr m_impl;
};

template<class T, unsigned flags>
struct flag_to_type<create_table_flags::can_set_temp, true, T, flags, create_table_template_base>
{    
    using result = create_table_temp<flags,T>;
};

template<class T, unsigned flags>
struct flag_to_type<create_table_flags::can_set_if_not_exists, true, T, flags, create_table_template_base>
{    
    using result = create_table_if_not_exists<flags,T>;
};

template<class T, unsigned flags>
struct flag_to_type<create_table_flags::can_set_table, true, T, flags, create_table_template_base>
{
    using result = create_table_tab<flags,T>;
};

template<class T, unsigned flags>
struct flag_to_type<create_table_flags::can_set_cols, true, T, flags, create_table_template_base>
{    
    using result = create_table_cols<flags,T>;
};

template<class T, unsigned flags>
struct flag_to_type<create_table_flags::can_add_constraints, true, T, flags, create_table_template_base>
{    
    using result = create_table_constraints<flags,T>;
};


template<class T, unsigned flags>
struct flag_to_type<create_table_flags::can_set_select, true, T, flags, create_table_template_base>
{
    using result = create_table_data<flags,T>;
};

template<class T, unsigned flags>
struct flag_to_type<create_table_flags::can_to_str, true, T, flags, create_table_template_base>
{
    using result = create_table_to_str<flags,T>;
};

}}}
