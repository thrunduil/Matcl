/*
 *  This file is a part of Matrix Computation Library (MATCL)
 *
 *  Copyright (c) Pawe³ Kowal 2017
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

#include "matcl-core/options/matcl_options.h"
#include <boost/lexical_cast.hpp>
#include <iostream>

namespace matcl 
{

//-----------------------------------------------------------------------------------
//                              option
//-----------------------------------------------------------------------------------
option::option(const impl_type& impl)
    :m_impl(impl)
{};

option::~option()
{};

std::string option::name() const
{
    return m_impl->name();
};
std::string option::description() const
{
    return m_impl->description();
};
bool option::has_value() const
{
    return m_impl->has_value();
};

void option::accept(option_visitor& val) const
{
    return m_impl->accept(val, *this);
};
std::string option::get_option_type() const
{
    return m_impl->get_option_type();
}
void option::help(const disp_stream_ptr& ds,const options& disp_options) const
{
    m_impl->help(ds,*this,disp_options);
};

static std::string replace_all(std::string source, const std::string& from, const std::string& to )
{
    std::string new_string;
    new_string.reserve( source.length() );  // avoids a few memory allocations

    std::string::size_type last_pos = 0;
    std::string::size_type find_pos;

    while( std::string::npos != ( find_pos = source.find( from, last_pos )))
    {
        new_string.append( source, last_pos, find_pos - last_pos );
        new_string += to;
        last_pos = find_pos + from.length();
    }

    // Care for the rest after last occurrence
    new_string += source.substr( last_pos );

    return new_string;
}

std::string option::make_name(const std::string& type_name)
{
    if (type_name[0] == 'c')
    {
        //class ***
        return replace_all(type_name.substr(6), "::", ".");
    }
    else if (type_name[0] == 's')
    {
        //struct ***
        return replace_all(type_name.substr(7),"::",".");
    }
    else
    {
        return replace_all(type_name,"::",".");
    };
};

void option::error_invalid_option_type(const std::string& get_type) const
{
    throw error::invalid_option_type(name(), get_type, this->get_option_type());
};

//-----------------------------------------------------------------------------------
//                              options
//-----------------------------------------------------------------------------------
options::options()
    :m_options(impl_type(new details::options_impl()))
{};

options::~options()
{};

options::options(const option& one_option)
    :m_options(impl_type(new details::options_impl()))
{
    this->set(one_option);
};

options::options(std::initializer_list<option> options_vec)
    :m_options(impl_type(new details::options_impl()))
{
    init(options_vec.begin(), options_vec.end());
};
        
options::options(const std::vector<option>& options_vec)
    :m_options(impl_type(new details::options_impl()))
{
    init(options_vec.begin(), options_vec.end());
};

options::options(const options& other)
    :m_options(other.m_options)
{}
options::options(options&& other)
    :m_options(std::move(other.m_options))
{};

options& options::operator=(const options& other) &
{
    m_options = other.m_options;
    return *this;
};
options& options::operator=(options&& other) &
{
    m_options = std::move(other.m_options);
    return *this;
};

template<class Iterator>
void options::init(Iterator beg, Iterator end)
{
    while(beg != end)
    {
        this->set(*beg);
        ++beg;
    };
};

void options::make_unique()
{
    if (is_unique() == true)
        return;

    m_options = m_options->copy();
}
bool options::is_unique() const
{
    return m_options.unique();
};

void options::install_notifier(const notifier_ptr& notif)
{
    m_options->install_notifier(notif);
};

static options m_default_options;

options& options::default_options()
{
    return m_default_options;
};
void options::clear()
{
    make_unique();
    m_options->clear();
};
Integer options::size() const
{
    return m_options->size();
};

bool options::has_option(const option& option_type) const
{
    return (*m_options).has_option(option_type);
}

options& options::set(const option& option_value)
{
    make_unique();
    m_options->set(option_value);
    return *this;
};
options& options::remove(const option& option_value)
{
    make_unique();
    m_options->remove(option_value);
    return *this;
};

options& options::set(const options& other)
{
    make_unique();
    m_options->set(*other.m_options);
    return *this;
};
options& options::remove(const options& other)
{
    make_unique();
    m_options->remove(*other.m_options);
    return *this;
};

void options::disp(const disp_stream_ptr& ds, const options& disp_options) const
{
    m_options->disp(ds, disp_options);
};
void options::help(const disp_stream_ptr& ds, const options& disp_options)
{
    details::options_impl::help(ds,disp_options);
};

void matcl::disp(const options& opts, const disp_stream_ptr& ds, const options& disp_opts)
{
    return opts.disp(ds, disp_opts);
};
// Display information about all available options
void matcl::help(const options& opts, const disp_stream_ptr& ds, const options& print_options)
{
    (void)opts;
    return options::help(ds, print_options);
}
void matcl::help(const option& opts, const disp_stream_ptr& ds, const options& print_options)
{
    return opts.help(ds,print_options);
};

};
