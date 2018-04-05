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

namespace matcl { namespace sql
{

template<unsigned flags>
class create_index_template;

template<unsigned flags, class T>
class create_index_unique;

template<unsigned flags, class T>
class create_index_if_not_exists;

template<unsigned flags, class T>
class create_index_name;

template<unsigned flags, class T>
class create_index_on;

template<unsigned flags, class T>
class create_index_cols;

template<unsigned flags, class T>
class create_index_to_str;

}};

namespace matcl { namespace sql { namespace details
{

class create_index_base_impl
{
    public:
        virtual std::string to_str() const = 0;

        virtual void unique(bool b) = 0;
        virtual void if_not_exists(bool b) = 0;
        virtual void name(const std::string&) = 0;
        virtual void on(const table&) = 0;
        virtual void cols(const column_list&) = 0;
};

class MATCL_SQLITE_EXPORT_MACRO create_index_template_base
{
    public:
        create_index_template_base();

    protected:
        using impl_ptr  = std::shared_ptr<create_index_base_impl>;

        create_index_template_base(const impl_ptr& impl);
   
        impl_ptr m_impl;
};

template<class T, unsigned flags>
struct flag_to_type<create_index_flags::can_set_unique, true, T, flags, create_index_template_base>
{    
    using result = create_index_unique<flags, T>;
};

template<class T, unsigned flags>
struct flag_to_type<create_index_flags::can_set_if_not_exists, true, T, flags, create_index_template_base>
{    
    using result = create_index_if_not_exists<flags, T>;
};

template<class T, unsigned flags>
struct flag_to_type<create_index_flags::can_set_name, true, T, flags, create_index_template_base>
{    
    using result = create_index_name<flags, T>;
};

template<class T, unsigned flags>
struct flag_to_type<create_index_flags::can_set_table, true, T, flags, create_index_template_base>
{    
    using result = create_index_on<flags, T>;
};

template<class T, unsigned flags>
struct flag_to_type<create_index_flags::can_set_cols, true, T, flags, create_index_template_base>
{    
    using result = create_index_cols<flags, T>;
};

template<class T, unsigned flags>
struct flag_to_type<create_index_flags::can_to_str, true, T, flags, create_index_template_base>
{    
    using result = create_index_to_str<flags, T>;
};

}}}