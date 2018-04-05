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

#include "matcl-sqlite-cpp/builder/update.h"

namespace matcl { namespace sql { namespace details
{

class update_impl : public update_base_impl
{
    public:
        update_impl();

        void        on_conflict(conflict_behaviour) override;
        void        tab(const table&) override;
        void        set(const assignment_list&) override;
        void        where(const expression&) override;

        std::string to_str() const override;

    private:
        conflict_behaviour  m_cb;
        table               m_tab;
        assignment_list     m_al;
        expression          m_cond;
};

}}}
