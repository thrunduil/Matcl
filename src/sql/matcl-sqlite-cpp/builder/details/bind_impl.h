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

#include <string>

namespace matcl { namespace sql { namespace details
{

class bind_param_base_impl
{
    public:
        virtual ~bind_param_base_impl() {};
        virtual std::string to_str() const = 0;
};

class default_bind_param_impl : public bind_param_base_impl
{
    public:
        std::string to_str() const override;
};

class numbered_bind_param_impl: public bind_param_base_impl
{
    public:
        numbered_bind_param_impl(unsigned i);
        std::string to_str() const override;

    private:
        unsigned m_nr;
};

class named_bind_param_impl : public bind_param_base_impl
{
    public:
        named_bind_param_impl(const std::string& s);
        std::string to_str() const override;

    private:
        std::string m_name;
};
    
}}}
