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

#include "matcl-sqlite-cpp/builder/create_index.h"

namespace matcl { namespace sql { namespace details
{

class create_index_impl : public create_index_base_impl
{
    public:
        create_index_impl();

        std::string to_str() const override;
        void        unique(bool b) override;
        void        if_not_exists(bool b) override;
        void        name(const std::string&) override;
        void        on(const table&) override;
        void        cols(const column_list&) override;

    private:
        bool        m_unique;
        bool        m_if_not_exists;
        std::string m_name;
        table       m_table;
        column_list m_cols;
};

}}}
