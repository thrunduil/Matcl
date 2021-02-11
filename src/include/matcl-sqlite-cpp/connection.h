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

#pragma once

#pragma warning(disable:4251)	// needs to have dll-interface to be used by clients

#include "config.h"
#include "typedefs.h"

#include <string>

namespace matcl { namespace sql
{

//connection object is destroyed only if there are no queriess or blobs
//associated with given object. Given connection must be unique in
//qiven thread unless database in opened in serialize thread mode.

class MATCL_SQLITE_EXPORT_MACRO connection 
{
	private:
		using data_ptr  = details::connection_data_ptr;

	public:		
		connection() {};
		connection(details::connection_data_ptr m_data);
		~connection() {};

		connection(const connection&);
		connection& operator=(const connection&);
		
		bool				empty() const;

		query				execute(const std::string& sql_str);
		sql					prepare(const std::string& sql_str);

		void				timeout_limit(int msec);

		int					number_of_changed_rows();
		int					total_number_of_changed_rows();	

		//If the row that a BLOB handle points to is modified by an UPDATE, DELETE, 
		//or by ON CONFLICT side-effects then the BLOB handle is invalidated
		database_blob_ptr	open_blob(const char* table, const char* column, int64_t row,
						        bool allow_write, const char * attached_db_name = NULL);


		//parallel table modification may lead to invalid result
		int64_t				last_insert_id() const;

		sqlite3 *			get_db();
		const sqlite3 *		get_db() const;

	private:		
		data_ptr			m_data;
};

}}