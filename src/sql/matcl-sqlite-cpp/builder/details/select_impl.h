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

#include <vector>
#include "matcl-sqlite-cpp/builder/select.h"

namespace matcl { namespace sql { namespace details
{

class select_impl : public select_base_impl
{
    public:
        select_impl();

        std::string to_str() const override;

        void    from(const join_source& js) override;
        void    where(const expression& e) override;
        void    add_column(const expression& e, const std::string& alias) override;
        void    add_table_columns(const table& t) override;
        void    distinct(bool b) override;
        void    group_by(const expression_list& e, const expression &constraint) override;
        void    order_by(const ord_expression_list& o) override;
        void    limit(const expression& e, const expression& offset) override;
        void    compose(const std::shared_ptr<select_base_impl>& impl, compose_type c) override;

    private:
        // single select-core to str
        std::string to_str_internal() const;
    
        using result_column         = std::pair<expression, std::string>;
        using result_column_list    = std::vector<result_column>;

        void    print_col(std::string& query, const result_column &rc) const;

        // FROM
        join_source         m_source;
        bool                m_source_exists;

        // WHERE
        expression          m_cond;

        // <result columns>
        result_column_list  m_cols;

        // ORDER BY
        ord_expression_list m_ordering;

        // LIMIT
        expression          m_limit;
        expression          m_offset;

        // GROUP BY
        expression_list     m_grouping;
        expression          m_having;

        // DISTINCT
        bool                m_distinct;

        // UNION (ALL) / INTERSECT / EXCEPT
        std::vector<std::pair<std::shared_ptr<select_impl>, compose_type> >
                            m_compounds;
};

}}}

