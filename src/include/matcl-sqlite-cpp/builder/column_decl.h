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
#include "matcl-sqlite-cpp/builder/details/fwd_decls.h"

#include "matcl-sqlite-cpp/builder/comma_list.h"
#include "matcl-sqlite-cpp/builder/expression.h"
#include "matcl-sqlite-cpp/builder/enums.h"

namespace matcl { namespace sql
{

class MATCL_SQLITE_EXPORT_MACRO column_decl
{
    private:
        using impl_ptr  = std::shared_ptr<details::column_decl_impl>;
        impl_ptr    m_impl;

    public:
        column_decl(const std::string& name);

        std::string     to_str() const;

        column_decl&    type(column_type ct);
        column_decl&    primary_key(bool ascending = true, 
                            conflict_behaviour cb = conflict_behaviour::default_cb,
                            bool autoincrement = false);

        column_decl&    not_null(conflict_behaviour cb = conflict_behaviour::default_cb);
        column_decl&    unique(conflict_behaviour cb = conflict_behaviour::default_cb);
        column_decl&    check(const expression& e);
        column_decl&    default_val(const expression& e);
        column_decl&    collate(collation c);

         // specify foreign key.
         // for advanced foreign key features see http://www.sqlite.org/foreignkeys.html#fk_deferred 
        column_decl&    foreign_key(const table& t, const column_list& cl = column_list());
};

using column_decl_list  = comma_list<column_decl>;

}}
