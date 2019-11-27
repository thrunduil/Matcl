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

#include "matcl-file/matcl_file.h"
#include "matcl-sqlite-cpp/matcl_sqlite_cpp.h"

namespace matcl { namespace details
{

struct mat_insert_elem
{
	Matrix			m_matrix;
	std::string		m_name;
	std::string		m_mat_string;
	bool			allow_replace;
};

struct data_insert_elem
{
	const void*     m_data;
    size_t          m_bytes;
	std::string		m_name;
	std::string		m_mat_string;
	bool			allow_replace;
};

class matcl_file_data
{
    private:
        using connection    = sql::connection;

	public:
		using save_mat_list     = std::list<mat_insert_elem>;
        using save_data_list    = std::list<data_insert_elem>;

	public:
		matcl_file_data(const std::string& file_name, open_mode, thread_mode);
		~matcl_file_data();

		void		add_to_save_list(const Matrix& mat, const std::string& mat_name, 
									const std::string& mat_string, bool allow_replace);
		void		add_data_to_save_list(const void* data, size_t bytes, const std::string& mat_name, 
									const std::string& mat_string, bool allow_replace);

	public:
		connection		        m_connection;
		save_mat_list			m_save_mat_list;
        save_data_list			m_save_data_list;
};

};};