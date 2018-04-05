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

#include <cassert>
#include "matcl-sqlite-cpp/builder/table_constraint.h"
#include "matcl-sqlite-cpp/builder/table.h"

namespace matcl { namespace sql { namespace details
{

class table_constraint_impl
{
    public:
        table_constraint_impl();

        std::string to_str() const;

        void        primary_key(const column_list& cols, conflict_behaviour cb);
        void        unique(const column_list& cols, conflict_behaviour cb);
        void        check(const expression& e);
        void        foreign_key(const column_list& local_cols, const table& foreign_table,
                        const column_list& foreign_cols);

    private:
        struct primary_key_info
        {
            primary_key_info(const column_list& cl, conflict_behaviour cb)
                : m_cols(cl) , m_cb(cb)
            { 
                assert(!m_cols.empty()); 
            }

            column_list         m_cols;
            conflict_behaviour  m_cb;
        };

        struct foreign_key_info
        {
            foreign_key_info(const column_list& local, const table& t, const column_list& foreign)
                : m_local_cols(local), m_tab(t), m_foreign_cols(foreign)
            {}

            column_list m_local_cols;
            table       m_tab;
            column_list m_foreign_cols;
        };

        std::shared_ptr<primary_key_info>   m_pkey;
        std::vector<primary_key_info>       m_uniques;
        std::vector<foreign_key_info>       m_fkeys;
        std::vector<expression>             m_checks;
};

}}}
