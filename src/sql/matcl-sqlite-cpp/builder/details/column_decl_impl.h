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

#include "matcl-sqlite-cpp/builder/column_decl.h"
#include "matcl-sqlite-cpp/builder/table.h"

namespace matcl { namespace sql { namespace details
{

class column_decl_impl
{
    public:
        column_decl_impl(const std::string& name);

        std::string to_str() const;

        void    type(column_type ct);
        void    primary_key(bool ascending, conflict_behaviour cb, bool autoincrement);
        void    not_null(conflict_behaviour cb);
        void    unique(conflict_behaviour cb);
        void    check(const expression& e);
        void    default_val(const expression& e);
        void    collate(collation c);
        void    foreign_key(const table& t, const column_list& cl);

    private:
        struct primary_key_info
        {
            primary_key_info(bool a, conflict_behaviour c, bool b)
                : m_ascending(a), m_cb(c), m_autoincrement(b)
            {}

            bool                m_ascending;
            conflict_behaviour  m_cb;
            bool                m_autoincrement;
        };

        struct foreign_key_info
        {
            foreign_key_info(const table& t, const column_list& c)
                : m_tab(t), m_cols(c)
            {}

            table       m_tab;
            column_list m_cols;
        };

        std::string                         m_name;
        column_type                         m_type;
        std::shared_ptr<primary_key_info>   m_pkey;
        std::shared_ptr<conflict_behaviour> m_not_null;
        std::shared_ptr<conflict_behaviour> m_unique;
        expression                          m_check;
        expression                          m_default_val;
        std::shared_ptr<collation>          m_collate;
        std::shared_ptr<foreign_key_info>   m_fkey;
};

}}}
