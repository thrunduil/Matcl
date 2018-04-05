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
#include <vector>
#include <memory>

#include "matcl-sqlite-cpp/builder/details/fwd_decls.h"
#include "matcl-sqlite-cpp/config.h"
#include "matcl-sqlite-cpp/builder/comma_list.h"

namespace matcl { namespace sql
{

// column declaration.
class MATCL_SQLITE_EXPORT_MACRO column
{
    private:
        using impl_ptr  = std::shared_ptr<details::column_impl>;
        impl_ptr    m_impl;

    public:
        explicit column(const std::string& name);

        void        set_table(const table& t);

        // get column name, with table reference, if available.
        std::string to_str() const;    

        // get bare column name.
        std::string get_name() const;

        // this doesn't really belong here, but operator= can't be declared 
        // outside class in C++
        assignment operator=(const expression& e);
};

using column_list   = comma_list<column>;

}}
