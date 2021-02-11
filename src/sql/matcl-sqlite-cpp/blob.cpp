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

#include "matcl-sqlite-cpp/blob.h"
#include <sqlite3.h>
#include "matcl-sqlite-cpp/exception.h"
#include "details/connection_details.h"

namespace matcl { namespace sql
{

blob_owner::blob_owner(void* ptr, size_t bytes)
: m_ptr(ptr), m_bytes(bytes)
{};

blob_owner::blob_owner(size_t bytes)
:m_bytes(bytes),m_ptr(NULL)
{	
	m_ptr = new char[bytes];	
};

charptr_owner::charptr_owner(size_t bytes)
:m_bytes(bytes),m_ptr(NULL)
{	
	m_ptr = new char[m_bytes];	
};

database_blob::~database_blob()
{
	sqlite3_blob_close(m_blob);
};

database_blob::database_blob(details::connection_data_ptr con, sqlite3_blob* blob)
:m_connection(con),m_blob(blob)
{
	m_bytes = sqlite3_blob_bytes(m_blob);
};

int	database_blob::bytes()
{
	return m_bytes;
};

void database_blob::read(void *buffer, int n_bytes, int start_pos)
{
	int rc = sqlite3_blob_read(m_blob, buffer, n_bytes, start_pos);

    if (rc != SQLITE_OK)
		throw matcl_sqlite_exception(m_connection->m_handle->last_error());
};

void database_blob::write(const void *source, int n_bytes, int start_pos)
{
	int rc = sqlite3_blob_write(m_blob, source, n_bytes, start_pos);
	if (rc != SQLITE_OK)
		throw matcl_sqlite_exception(m_connection->m_handle->last_error());
};

}}