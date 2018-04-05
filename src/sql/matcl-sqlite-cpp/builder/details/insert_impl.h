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

#include "matcl-sqlite-cpp/builder/insert.h"
#include "matcl-sqlite-cpp/builder/table.h"

namespace matcl { namespace sql { namespace details
{

class insert_impl : public insert_base_impl
{
    public:
        insert_impl();

        void on_conflict(conflict_behaviour cb) override;
        void columns(const column_list& cl) override;
        void into(const table& t) override;
        void values(const details::select_template_base& sel) override;
        void values(const expression_list& el) override;

        std::string to_str() const override;

    private:
        conflict_behaviour  m_cb;
        column_list         m_cl;
        table               m_tab;
        std::shared_ptr<details::select_template_base> 
                            m_sel;
        expression_list     m_vals;
};

}}}

