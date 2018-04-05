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

#include "config.h"
#include <string>
#include <exception>

#pragma warning(push)
#pragma warning(disable:4251)   //needs to have dll-interface to be used by clients
#pragma warning(disable:4275)   //non dll-interface class used as base for dll-interface class

namespace matcl { namespace sql
{

class MATCL_SQLITE_EXPORT_MACRO matcl_sqlite_exception : public std::exception
{
	protected:
		std::string message;
		std::string raw_message;

	public:
		matcl_sqlite_exception(const std::string& msg)		: message(msg){};
		matcl_sqlite_exception()							: message(""){};

	    virtual const char* what() const;
		const char*			sqlite_message() const;
};

class MATCL_SQLITE_EXPORT_MACRO matcl_sqlite_open_exception : public matcl_sqlite_exception
{
	public:
		matcl_sqlite_open_exception(const std::string& db_name, const std::string& msg);		
};

class MATCL_SQLITE_EXPORT_MACRO matcl_sqlite_query_exception : public matcl_sqlite_exception
{
	public:
		int sqlite_code;

	public:
		matcl_sqlite_query_exception(const std::string& msg, int code);
};

class MATCL_SQLITE_EXPORT_MACRO matcl_sqlite_lock_exception : public matcl_sqlite_exception
{
	public:
		int sqlite_code;

	public:
		matcl_sqlite_lock_exception(int code);
};

class MATCL_SQLITE_EXPORT_MACRO matcl_sqlite_interrupt_exception : public matcl_sqlite_exception
{
	public:
		int sqlite_code;

	public:
		matcl_sqlite_interrupt_exception(int code);
};

}}

#pragma warning(pop)