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

#include <string>
#include <memory>
#include "matcl-sqlite-cpp/config.h"
#include "matcl-sqlite-cpp/builder/details/type_builder.h"
#include "matcl-sqlite-cpp/builder/table.h"
#include "matcl-sqlite-cpp/builder/expression.h"
#include "matcl-sqlite-cpp/builder/column.h"
#include "matcl-sqlite-cpp/builder/ord_expression.h"
#include "matcl-sqlite-cpp/builder/details/flags.h"
#include "matcl-sqlite-cpp/builder/details/select_impl.h"

namespace matcl { namespace sql
{

template<unsigned flags>
class select_template : public details::choose_base<flags, details::select_template_base>::result
{
    public:
        select_template() {}
        select_template(const details::select_template_base::impl_ptr& impl)
            : details::choose_base<flags, details::select_template_base>::result(impl) 
        {}
};

template<unsigned flags, class T>
class select_template_from : public T
{
    private:
        using result_type = select_template<flags & ~details::select_flags::can_set_source>;

    public:
        // define join source to read from.
        result_type from(const join_source& arg)
        {
            result_type rv(T::m_impl);
            T::m_impl->from(arg);
            T::m_impl.reset();
            return rv;
        }

        select_template_from() {}
        select_template_from(const typename T::impl_ptr& impl)
            : T(impl) 
        {}
};

template<unsigned flags, class T>
class select_template_distinct : public T
{
    public:
        // remove duplicates?
        select_template<flags & ~details::select_flags::can_set_distinct> 
            distinct(bool b = true)
        {
            select_template<flags & ~details::select_flags::can_set_distinct> rv(T::m_impl);
            T::m_impl->distinct(b);
            T::m_impl.reset();
            return rv;
        }

        select_template_distinct() {}
        select_template_distinct(const typename T::impl_ptr& impl)
            : T(impl) 
        {}
};

template<unsigned flags, class T>
class select_template_columns : public T
{
    public:
        // define column to select.
        select_template<flags> 
            operator()(const expression& col, const std::string& alias="")
        {
            select_template<flags> rv(T::m_impl);
            T::m_impl->add_column(col, alias);
            T::m_impl.reset();
            return rv;
        };

        // define table to select its columns.
        select_template<flags> operator()(const table& t)
        {
            select_template<flags> rv(T::m_impl);
            T::m_impl->add_table_columns(t);
            T::m_impl.reset();
            return rv;
        }
            
        select_template_columns() {}
        select_template_columns(const typename T::impl_ptr& impl)
            : T(impl) 
        {}
};

template<unsigned flags, class T>
class select_template_where : public T
{
    private:
        using result_type = select_template<flags & ~details::select_flags::can_set_where>;

    public:
        // define conditions.
        result_type where(const expression& arg)
        {
            result_type rv(T::m_impl);
            T::m_impl->where(arg);
            T::m_impl.reset();
            return rv;
        }

        select_template_where() {}
        select_template_where(const typename T::impl_ptr& impl)
            : T(impl) 
        {}
};

template<unsigned flags, class T>
class select_template_group_by : public T
{
    private:
        using result_type = select_template<flags & ~details::select_flags::can_set_grouping>;

    public:
        // define groupings.
        result_type group_by(const expression_list& e, const expression &constraint = expression())
        {
            result_type rv(T::m_impl);
            T::m_impl->group_by(e, constraint);
            T::m_impl.reset();
            return rv;
        }

        // define groupings.
        result_type group_by(const expression& e, const expression &constraint = expression())
        {
            return group_by(expression_list(e), constraint);
        }
            
        select_template_group_by() {}
        select_template_group_by(const typename T::impl_ptr& impl) 
            : T(impl) 
        {}

};

template<unsigned flags, class T>
class select_template_order_by : public T
{
    private:
        using result_type = select_template<
                                flags &
                                ~details::select_flags::can_compose &
                                ~details::select_flags::core_flags &
                                ~details::select_flags::can_set_ordering
                            >;

    public:
        // define orderings.
        result_type order_by(const ord_expression_list& o)
        {
            result_type rv(T::m_impl);
            T::m_impl->order_by(o);
            T::m_impl.reset();
            return rv;
        }
            
        // define orderings. (simplified version)
        result_type order_by(const expression& o)
        {
            return order_by(ord_expression_list(o));
        }

        select_template_order_by() {}
        select_template_order_by(const typename T::impl_ptr& impl)
            : T(impl) 
        {}
};

template<unsigned flags, class T>
class select_template_limit : public T
{
    private:
        using result_type   = select_template<
                                flags &
                                ~details::select_flags::can_compose &
                                ~details::select_flags::core_flags &
                                ~details::select_flags::can_set_limit
                            >;

    public:
        // set limits on result count.
        result_type limit(const expression& e, const expression& offset = expression())
        {
            result_type rv(T::m_impl);
            T::m_impl->limit(e, offset);
            T::m_impl.reset();
            return rv;
        }
            
        select_template_limit() {}
        select_template_limit(const typename T::impl_ptr& impl)
            : T(impl) 
        {}
};

template<class T, unsigned flags>
class select_template_composer : public T
{
    private:
        using result_type = select_template<
                                flags &
                                ~details::select_flags::core_flags
                            > ;

        template<unsigned other_flags>
        result_type _compose(select_template<other_flags> other, compose_type ct)
        {
            static_assert(other_flags & details::select_flags::can_compose, "Incompatible operand!");

            result_type rv(T::m_impl);
            T::m_impl->compose(other.m_impl, ct);
            T::m_impl.reset();
            other.m_impl.reset();
            return rv;
        }

        // COMPILER: whole class have to be declared outside flag_to_type due to these two lines
        template<class S, unsigned other_flags>
        friend class select_template_composer; 

    public:
        // union of results.
        template<unsigned other_flags>
        inline result_type union_with(select_template<other_flags> other)
        {
            return _compose(other, compose_type::union_type);
        }

        // union of results, without removing duplicates
        template<unsigned other_flags>
        inline result_type union_all(select_template<other_flags> other)
        {
            return _compose(other, compose_type::union_all);
        }
        
        // intersection of results
        template<unsigned other_flags>
        inline result_type intersect(select_template<other_flags> other)
        {
            return _compose(other, compose_type::intersect);
        }

        // set difference of results
        template<unsigned other_flags>
        inline result_type except(select_template<other_flags> other)
        {
            return _compose(other, compose_type::except);
        }

        select_template_composer() {}
        select_template_composer(const typename T::impl_ptr& impl) 
            : T(impl) 
        {}
};


using select = select_template<details::select_flags::default_flags>;

}}

