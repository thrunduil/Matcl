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

#include "matcl-sqlite-cpp/builder/table.h"

namespace matcl { namespace sql { namespace details
{

struct join_def_base
{
    virtual std::string op() const = 0;
    virtual std::string constraint() const = 0;
};

class table_impl
{
    public:
        table_impl(const std::string& name, const std::string& alias);

        std::string     to_str() const;
        std::string     get_ref() const;
        std::string     get_name() const;

        void            join(const table_impl& t, bool natural, join_type jt);
        void            join(const table_impl& t, const column_list& c, join_type jt);
        void            join(const table_impl& t, const expression& e, join_type jt);

    private:
        void            join(const table_impl& t, join_def_base* jdef); // takes ownership!

        struct join_info
        {
            std::string                     m_joinee;
            std::shared_ptr<join_def_base>  m_join;
        };

        std::vector<join_info>  m_joins;
        std::string             m_name;
        std::string             m_alias;
};

struct simple_join : public join_def_base
{
    public:
        simple_join(join_type jt);

        std::string op() const override;
        std::string constraint() const override;

    protected:
        const join_type m_jt;
};

struct natural_join : public simple_join
{
    public:
        natural_join(join_type jt);
        std::string op() const override;
};

struct column_join : public simple_join
{
    public:
        column_join(const column_list& cl, join_type jt);
        std::string constraint() const override;

    private:
        column_list m_cols;
};

struct expr_join : public simple_join
{
    public:
        expr_join(const expression& e, join_type jt);
        std::string constraint() const override;

    private:
        expression m_expr;
};

}}}
