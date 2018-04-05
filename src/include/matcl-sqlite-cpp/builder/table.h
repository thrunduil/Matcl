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
#include "matcl-sqlite-cpp/builder/details/fwd_decls.h"
#include "matcl-sqlite-cpp/builder/column.h"
#include "matcl-sqlite-cpp/builder/expression.h"

namespace matcl { namespace sql
{

class MATCL_SQLITE_EXPORT_MACRO join_source
{
    public:
        // conversion to SQL.
        std::string to_str() const;

        // join other table, t: table to join, if natural = true, then natural join, 
        // jt: join type.
        join_source join(const join_source& t, bool natural = false, join_type jt = join_type::full);

        // join other table, t: table to join, c: columns to match, jt: join type
        join_source join(const join_source& t, const column_list& c, join_type jt = join_type::full);
        
        // join other table, t: table to join, e: condition to match, jt: join type
        join_source join(const join_source& t, const expression& e, join_type jt = join_type::full);

    protected:
        using impl_ptr = std::shared_ptr<details::table_impl>;
        impl_ptr m_impl;

        join_source(const impl_ptr& impl);
};

class MATCL_SQLITE_EXPORT_MACRO table : public join_source
{
    public:
        explicit table(const std::string& name, const std::string& alias="");

        // get table reference name: alias if exists, or regular name otherwise.
        std::string get_ref() const;

        // get table regular name.
        std::string get_name() const;
};

}}
