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

#include "matcl-sqlite-cpp/connection.h"
#include "matcl-sqlite-cpp/exception.h"
#include <sqlite3.h>
#include "matcl-sqlite-cpp/blob.h"
#include "matcl-sqlite-cpp/sql.h"
#include "matcl-sqlite-cpp/query.h"
#include "details/connection_details.h"

namespace matcl { namespace sql
{

namespace details
{

connection_data::connection_data(database db)
	: m_db(db)
{
	db.get_db_handle(m_handle);
};

connection_data::~connection_data()
{};

};

connection::connection(details::connection_data_ptr data) 
:m_data(data)
{
	timeout_limit(max_wait_time_msec);
};

connection::connection(const connection& other)
:m_data(other.m_data)
{};

connection& connection::operator=(const connection& other)
{
	m_data = other.m_data;
	return *this;
};

sqlite3 * connection::get_db()
{
	return m_data->m_handle->get_db();
};

void connection::timeout_limit(int msec)
{
	sqlite3_busy_timeout(get_db(), msec);
};

const sqlite3 *	connection::get_db() const
{
	return m_data->m_handle->get_db();
};

query connection::execute(const std::string& sql_str)
{
	sql ns = create_sql(m_data,sql_str);
	return ns.make();
}

sql connection::prepare(const std::string& sql_str)
{
	sql ns = create_sql(m_data,sql_str);
	return ns;
}

int64_t connection::last_insert_id() const
{
	return m_data->m_handle->last_insert_id();
}

int	connection::number_of_changed_rows()
{
	return sqlite3_changes(m_data->m_handle->get_db());
};

int	connection::total_number_of_changed_rows()
{
	return sqlite3_total_changes(m_data->m_handle->get_db());
};

bool connection::empty() const
{
	return !m_data;
};

database_blob_ptr connection::open_blob(const char* table, const char* column, int64_t row,
								bool allow_write, const char * attached_db_name)
{
	if (!attached_db_name)
		attached_db_name = "main";

	sqlite3_blob *p_blob = NULL;
	int flags = allow_write? 1 : 0;
	int rc = sqlite3_blob_open(m_data->m_handle->get_db(),attached_db_name,table,column,row,
								flags,&p_blob);  

	if (rc != SQLITE_OK)
	{
		sqlite3_blob_close(p_blob);

		switch(rc)
		{
			case SQLITE_BUSY:
			case SQLITE_LOCKED:
				throw matcl_sqlite_lock_exception(rc);
			case SQLITE_INTERRUPT:
				throw matcl_sqlite_interrupt_exception(rc);
			default:
				throw matcl_sqlite_query_exception(m_data->m_handle->last_error(),rc);
		};
	};

	return database_blob_ptr(new database_blob(m_data,p_blob));
};

}}
