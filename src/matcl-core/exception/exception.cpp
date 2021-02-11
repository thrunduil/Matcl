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

#include "matcl-core/error/exception.h"
#include "matcl-core/memory/global_objects.h"
#include "matcl-core/error/exception_message.h"
#include "matcl-core/general/thread.h"

#include <iostream>

#pragma warning( push )
#pragma warning(disable:4702)	// unreachable code
#include <boost/lexical_cast.hpp>
#pragma warning( pop )

namespace matcl { namespace details
{

enum class warning_state : int
{
    not_set     = 0,
    enabled     = 1,
    disabled    = 2
};

using atomic_warn_state = atomic<warning_state>;

atomic_warn_state m_global_warning_state  = warning_state::not_set;

MATCL_THREAD_LOCAL 
atomic_warn_state  m_local_warning_state   = warning_state::not_set;

static warning_state get_warning_state(bool global)
{
    if (global == true)
        return m_global_warning_state.load();
    else
        return m_local_warning_state.load();
}

static void set_warning_state(warning_state st, bool global)
{
    if (global == true)
        m_global_warning_state.store(st);
    else
        m_local_warning_state.store(st);
};

static bool check_warnings_enabled()
{
    warning_state gl    = get_warning_state(true);

    if (gl != warning_state::not_set)
        return gl == warning_state::enabled;

    warning_state loc   = get_warning_state(false);

    // if warning state is not set locally or globally, then
    // assume that warnings are enabled
    return (loc == warning_state::disabled) ? false : true;
};

enable_warnings_raii::enable_warnings_raii(bool enable, bool global)
    :m_global(global)
{
    m_old_state = static_cast<int>(get_warning_state(global));

    set_warning_state(enable ? warning_state::enabled 
                      : warning_state::disabled, global);
};

enable_warnings_raii::~enable_warnings_raii()
{
    if (m_old_state == -1)
        return;

    set_warning_state(static_cast<warning_state>(m_old_state), m_global);
};

enable_warnings_raii::enable_warnings_raii(enable_warnings_raii&& other)
    :m_old_state(other.m_old_state), m_global(other.m_global)
{
    other.m_old_state = -1;
};

}};

namespace matcl { namespace error
{

details::enable_warnings_raii error::enable_warnings(bool val, bool global)
{
    return details::enable_warnings_raii(val,global);
}

bool error::check_warnings_enabled()
{
    return details::check_warnings_enabled();
};

void error::warning(const std::string& msg)
{
    matcl::error::get_global_messanger()->warning(msg);
};

assert_exception::assert_exception(const std::string& file_, int line_, 
                                   const std::string& message_)
    :file(file_),line(line_),message(message_)
{
    message = message + " in file: " + file_ + " at line " 
            + boost::lexical_cast<std::string>(line);
};

std::string assert_exception::what() const throw()
{
    return message;
};

void assert_exception::make(const char* txt, const char* description, const char* file, int line)
{
    std::string message = "matcl_assert failed: ";
    message				+= txt;
    if (description)
    {
        message			+= ": ";
        message			+= description;
    };

    std::cerr << message << "\n";
    throw assert_exception(file,line,message);
}


const char* matcl_exception::what() const throw()
{
    m_message = what(*get_global_messanger());
    return m_message.c_str();
};

};};