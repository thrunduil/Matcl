/*
 *  This file is a part of Matrix Computation Library (MATCL)
 *
 *  Copyright (c) Pawe³ Kowal 2017 - 2021
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

#include "matcl-matrep/general/config.h"
#include "matcl-core/matrix/scalar_types.h"
#include "matcl-matrep/details/refcount.h"
#include <algorithm>
#include <cassert>

namespace matcl { namespace algorithm 
{

class mark_type_ref;

class mark_type
{
    private:
        Integer         m_mark;

    public:
        mark_type(mark_type_ref mr);

        bool            operator<(mark_type other) const        { return m_mark < other.m_mark; };
        bool            operator<=(mark_type other) const       { return m_mark <= other.m_mark; };
        bool            operator>=(mark_type other) const       { return m_mark >= other.m_mark; };
        bool            operator>(mark_type other) const        { return m_mark > other.m_mark; };
        bool            operator==(mark_type other) const       { return m_mark == other.m_mark; };
        bool            operator!=(mark_type other) const       { return m_mark != other.m_mark; };

        Integer         value() const                           { return m_mark; };

    private:
        mark_type(Integer m)    :m_mark(m) {};

        friend class scatter;
        friend class mark_type_ref;
};

class mark_type_ref
{
    private:
        Integer*        m_mark;

    public:
        mark_type_ref&  operator=(mark_type m)                  { *m_mark = m.m_mark; return *this; };

        bool            operator<(mark_type other) const        { return *m_mark < other.m_mark; };
        bool            operator<=(mark_type other) const       { return *m_mark <= other.m_mark; };
        bool            operator>=(mark_type other) const       { return *m_mark >= other.m_mark; };
        bool            operator>(mark_type other) const        { return *m_mark > other.m_mark; };
        bool            operator==(mark_type other) const       { return *m_mark == other.m_mark; };
        bool            operator!=(mark_type other) const       { return *m_mark != other.m_mark; };

        bool            operator<(mark_type_ref other) const    { return *m_mark < *other.m_mark; };
        bool            operator<=(mark_type_ref other) const   { return *m_mark <= *other.m_mark; };
        bool            operator>=(mark_type_ref other) const   { return *m_mark >= *other.m_mark; };
        bool            operator>(mark_type_ref other) const    { return *m_mark > *other.m_mark; };
        bool            operator==(mark_type_ref other) const   { return *m_mark == *other.m_mark; };
        bool            operator!=(mark_type_ref other) const   { return *m_mark != *other.m_mark; };

    private:
        mark_type_ref(Integer* m)   :m_mark(m) {};        

        friend class scatter;
        friend class mark_type;
};

inline mark_type::mark_type(mark_type_ref mr)
    :m_mark(*mr.m_mark)
{};

class MATCL_MATREP_EXPORT scatter
{
    private:
        using refcount_str = matcl::details::refcount_str<matcl::details::scatter_refcount_pool>;

    private:
        Integer         m_mark;
        Integer         m_tmp_mark;
        Integer*        m_ptr;
        Integer         m_size;
        refcount_str*   m_refcount;

    public:
        scatter(const scatter&);
        ~scatter();

        scatter&        operator=(const scatter&);

        static scatter  get(Integer size, Integer n_increase);
        
        mark_type       current_mark() const    { return mark_type(m_mark); };
        mark_type       next_mark()             { ++m_mark; return mark_type(current_mark()); };
        mark_type       create_mark(Integer m)  { m_tmp_mark = std::max(m, m_tmp_mark) ; return mark_type(m); };
        mark_type_ref   operator[](Integer i)   { assert(i < m_size); return mark_type_ref(&m_ptr[i]); };

    private:
        scatter();                

        friend class scatter_provider;
};

};};