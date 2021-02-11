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

#include "matcl-sqlite-cpp/exception.h"

namespace matcl { namespace sql
{

const char* matcl_sqlite_exception::what() const
{
    return message.c_str();
}

const char* matcl_sqlite_exception::sqlite_message() const
{
    return raw_message.c_str();
}

matcl_sqlite_open_exception::matcl_sqlite_open_exception(const std::string& db_name, 
                                const std::string& msg)
{
    (void)db_name;

	message = "unable to open database, reason: " + msg;
	raw_message = msg;
};

matcl_sqlite_query_exception::matcl_sqlite_query_exception(const std::string& msg,int code)
{
	sqlite_code = code;
	message = "unable to execute query, reason: " + msg;
	raw_message = msg;
}

matcl_sqlite_lock_exception::matcl_sqlite_lock_exception(int code)
{
	sqlite_code = code;
	message = "unable to execute query, reason: database or table is locked";
	raw_message = "database or table is locked";
};

matcl_sqlite_interrupt_exception::matcl_sqlite_interrupt_exception(int code)
{
	sqlite_code = code;
	message = "unable to execute query, reason: execution was interrupted";
	raw_message = "execution was interrupted";
};

}}