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
#include <memory>

struct sqlite3;
struct sqlite3_stmt;
struct sqlite3_blob;

namespace matcl { namespace sql
{

namespace details
{
	class db_handle;
	class lock_db_handle;
	class handle_pool;
	class connection_data;
	class database_data;
	class query_data;

	using connection_data_ptr   = std::shared_ptr<connection_data>;
	using database_data_ptr     = std::shared_ptr<database_data>;
	using query_data_ptr        = std::shared_ptr<query_data>;
};

class string_getter;
class blob_getter;

class connection;
class database;
class query;
class sql;

class blob_owner;
class charptr_owner;
class database_blob;

using blob_ptr          = std::shared_ptr<blob_owner>;
using char_ptr          = std::shared_ptr<charptr_owner>;
using database_blob_ptr = std::shared_ptr<database_blob>;

namespace enums
{
	enum open_mode
	{
		readonly			= 0x00000001, 
		readwrite			= 0x00000002 , 
		readwrite_create	= 0x00000002 | 0x00000004 
	};

	enum value_type
	{
		integer_type		= 1, 
		float_type			= 2,
		text_type			= 3, 
		blob_type			= 4, 
		null_type			= 5
	};

	enum thread_mode
	{
		single_thread		= 0x00000000,
		multi_thread		= 0x00008000,
		serialized			= 0x00010000
	};
};

}}