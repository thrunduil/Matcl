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
#pragma warning(disable:4251)	// needs to have dll-interface to be used by clients

#include "matcl-sqlite-cpp/config.h"
#include "matcl-sqlite-cpp/database.h"

#include <map>

namespace matcl { namespace sql { namespace details
{

class query_data
{
	public:
		using column_map    = std::map<std::string,int>;

	public:		
		~query_data();

	private:
		query_data(connection_data_ptr con);

		query_data(const query_data& q) = delete;
		query_data& operator=(const query_data& ) = delete;

	public:
		database			m_db;
		connection_data_ptr m_con;
		db_handle*			m_handle;

		sqlite3_stmt*		m_raw_stm;
		column_map*			m_cmap;
		int					m_column_number;
		bool				m_get_row_succeeded;
		int					m_old_result_flag;
		bool				m_use_old_result_flag;
		bool				m_need_reset;

		friend query create_query(const connection_data_ptr&,const std::string&);
		friend sql   create_sql(const connection_data_ptr&,const std::string&);
};

};}}