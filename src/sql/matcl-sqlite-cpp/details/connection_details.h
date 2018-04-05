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

#include "matcl-sqlite-cpp/config.h"
#include "matcl-sqlite-cpp/typedefs.h"
#include "matcl-sqlite-cpp/database.h"
#include "details/database_details.h"

namespace matcl { namespace sql { namespace details
{

class connection_data
{
	public:		
		~connection_data();
		database			m_db;
		lock_db_handle		m_handle;

	private:
		connection_data(database db);

		connection_data(const connection_data&) = delete;
		connection_data& operator=(const connection_data&) = delete;

		friend connection matcl::sql::create_connection(database db);
};

};}}