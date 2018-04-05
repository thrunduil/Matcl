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

#include "matcl_file_data.h"
#include "matcl-file/exception.h"

namespace matcl { namespace details
{

static sql::enums::open_mode sqpp_open_mode(open_mode om)
{
	switch (om)
	{
		case matcl::open_mode::readonly:
			return sql::enums::readonly;
		case matcl::open_mode::readwrite:
			return sql::enums::readwrite;
		case matcl::open_mode::readwrite_create:
			return sql::enums::readwrite_create;
		default:
			matcl_assert(0,"unknown case");
			throw;
	};
};

static sql::enums::thread_mode sqpp_thread_mode(thread_mode tm)
{
	switch (tm)
	{
		case matcl::thread_mode::single_thread:
			return sql::enums::single_thread;
		case matcl::thread_mode::multi_thread:
			return sql::enums::multi_thread;
		case matcl::thread_mode::serialized:
			return sql::enums::serialized;
		default:
			matcl_assert(0,"unknown case");
			throw;
	};
};

mmlib_file_data::mmlib_file_data(const std::string& file_name, open_mode om, thread_mode tm)
{
	sql::enums::open_mode om_q = sqpp_open_mode(om);
	sql::enums::thread_mode tm_q = sqpp_thread_mode(tm);

	try
	{
		sql::database db = sql::open_database(file_name, om_q, tm_q);
		m_connection = create_connection(db);
	}
	catch(sql::matcl_sqlite_open_exception& ex)
	{
		throw error::error_open_mmlibfile(file_name,ex.sqlite_message());
	}
	catch(sql::matcl_sqlite_lock_exception&)
	{
		throw error::error_mmlibfile_locked();
	}
	catch(std::exception& ex)
	{
		throw error::error_general(ex.what());
	}
	catch(...)
	{
		throw error::error_general("unknown exception during opening mmlib_file");
	};
};

mmlib_file_data::~mmlib_file_data()
{};

void mmlib_file_data::add_to_save_list(const Matrix& mat, const std::string& mat_name, 
						const std::string& mat_string, bool allow_replace)
{
	mat_insert_elem item;
	item.m_matrix = mat;
	item.m_name = mat_name;
	item.m_mat_string = mat_string;
	item.allow_replace = allow_replace;

	m_save_mat_list.push_back(item);
};

void mmlib_file_data::add_data_to_save_list(const void* data, size_t bytes, const std::string& mat_name, 
						const std::string& mat_string, bool allow_replace)
{
	data_insert_elem item;
	item.m_data = data;
    item.m_bytes = bytes;
	item.m_name = mat_name;
	item.m_mat_string = mat_string;
	item.allow_replace = allow_replace;

	m_save_data_list.push_back(item);
};

};};