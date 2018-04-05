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

#include <sstream>
#include "matcl-sqlite-cpp/builder/bind.h"
#include "sql/matcl-sqlite-cpp/builder/details/bind_impl.h"

namespace matcl { namespace sql
{

bind_param::bind_param()
    : m_impl(new details::default_bind_param_impl())
{}

bind_param::bind_param(unsigned i)
    : m_impl(new details::numbered_bind_param_impl(i))
{}

bind_param::bind_param(const std::string& s)
    : m_impl(new details::named_bind_param_impl(s))
{}

std::string bind_param::to_str() const
{
    return m_impl->to_str();
}

namespace details
{

std::string default_bind_param_impl::to_str() const
{
    return "?";
}

numbered_bind_param_impl::numbered_bind_param_impl(unsigned i)
    : m_nr(i)
{}

std::string numbered_bind_param_impl::to_str() const
{
    std::stringstream ss;
    ss << "?" << m_nr;
    return ss.str();
}

named_bind_param_impl::named_bind_param_impl(const std::string& name)
    : m_name(name)
{}

std::string named_bind_param_impl::to_str() const
{
    return ":" + m_name;
}

};

}}
