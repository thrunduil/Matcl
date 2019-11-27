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

#include "matcl-file/exception.h"

#pragma warning( push )
#pragma warning(disable:4702)	// unreachable code
#include <boost/lexical_cast.hpp>
#pragma warning( pop )

namespace matcl { namespace error
{

const char* matcl_file_exception::what() const throw()
{
    m_message = what(*get_global_messanger_file());
    return m_message.c_str();
};

const char* matcl_file_exception::what(exception_message& em) const
{
    (void)em;
    m_message = "invalid exception_message";
    return m_message.c_str();
};

const char* error_open_matfile::what(exception_message_file& em) const
{
    return em.error_open_matfile(file_name,msg);
};

const char* error_matfile_not_opened::what(exception_message_file& em) const
{
    return em.error_matfile_not_opened();
};

const char* error_read_matfile::what(exception_message_file& em) const
{
    return em.error_read_matfile();
};

const char* error_read_matfile_var::what(exception_message_file& em) const
{
    return em.error_read_matfile_var(var_name);
};

const char* error_write_matfile::what(exception_message_file& em) const
{
    return em.error_write_matfile();
};

const char* error_open_matclfile::what(exception_message_file& em) const
{
    return em.error_open_matclfile(file_name,msg);
};

const char* error_matclfile_locked::what(exception_message_file& em) const
{
    return em.error_matclfile_locked();
};

const char* error_create_matclfile::what(exception_message_file& em) const
{
    return em.error_create_matclfile(msg);
};

const char* error_read_matclfile::what(exception_message_file& em) const
{
    return em.error_read_matclfile(msg);
};

const char* error_read_matclfile_mat_not_exist::what(exception_message_file& em) const
{
    return em.error_read_matclfile_mat_not_exist(mat_name);
};

const char* error_write_matclfile_mat_already_exist::what(exception_message_file& em) const
{
    return em.error_write_matclfile_mat_already_exist(mat_name);
};

};};