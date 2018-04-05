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
#include "matcl-sqlite-cpp/builder/column.h"
#include "matcl-sqlite-cpp/builder/enums.h"

namespace matcl { namespace sql
{

class MATCL_SQLITE_EXPORT_MACRO table_constraint
{
    public:
        table_constraint();

        std::string to_str() const;

        table_constraint& primary_key(const column_list& cols, conflict_behaviour cb);
    
        table_constraint& unique(const column_list& cols, conflict_behaviour cb);

        table_constraint& check(const expression& e);

        table_constraint& foreign_key(const column_list& local_cols, const table& foreign_table,
                                const column_list& foreign_cols);

    private:
        using impl_ptr  = std::shared_ptr<details::table_constraint_impl>;
        impl_ptr m_impl;
};

using table_constraint_list = comma_list<table_constraint>;

}}

