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

namespace matcl { namespace sql { namespace details
{

class expression_base_impl
{
    public:
        virtual ~expression_base_impl() {}
        virtual std::string to_str() const = 0;
};

class column_expression : public expression_base_impl
{
    public:
        column_expression(const column& c);
        std::string     to_str() const override;

    private:
        column          m_c;
};

class numeric_expression : public expression_base_impl
{
    public:
        numeric_expression(double d);
        std::string     to_str() const override;

    private:
        double          m_num;
};

class string_expression : public expression_base_impl
{
    public:
        string_expression(const std::string& s);
        std::string     to_str() const override;

    private:
        std::string     m_str;
};

class null_expression : public expression_base_impl
{
    public:
        std::string     to_str() const override;
};

class binary_op_expression : public expression_base_impl
{
    public:
        binary_op_expression(const expression& e1, const expression& e2, 
                const std::string& op);

        std::string     to_str() const override;

    private:
        expression      m_e1;
        expression      m_e2;
        std::string     m_op;
};

class unary_op_expression : public expression_base_impl
{
    public:
        unary_op_expression(const expression& e, const std::string& op);
        std::string     to_str() const override;

    private:
        expression      m_e;
        std::string     m_op;
};

class bind_param_expression : public expression_base_impl
{
    public:
        bind_param_expression(const bind_param& b);
        std::string     to_str() const override;

    private:
        bind_param      m_bp;
};

class call_expression : public expression_base_impl
{
    public:
        call_expression(const std::string& name, const expression_list& e, 
                bool distinct = false);

        std::string     to_str() const override;

    private:
        std::string     m_name;
        expression_list m_args;
        bool            m_distinct;
};

class star_expression : public expression_base_impl
{
    public:
        std::string     to_str() const override;
};

class collate_expression : public expression_base_impl
{
    public:
        collate_expression(const expression& e, collation c);

        std::string     to_str() const override;

    private:
        expression      m_expr;
        collation       m_collation;
};

std::string collation_to_str(collation c);

}}}

