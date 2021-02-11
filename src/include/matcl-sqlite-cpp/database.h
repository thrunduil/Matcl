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
#include <string>
#include "typedefs.h"

namespace matcl { namespace sql
{

class MATCL_SQLITE_EXPORT_MACRO database 
{
	private:
		using data_ptr  = details::database_data_ptr;

	public:		
		database() {};
		database(details::database_data_ptr data);
		virtual ~database() {};

		database(const database& );
		database& operator=(const database& );	

		bool				empty() const;

		bool				is_connected();		
								
		void				get_db_handle(details::lock_db_handle&);
			
		void				interrupt(connection con);
		void				close_unused_connections();	
		int					number_of_connections();

		//if there are active connections these functions will throw an exception
		void				set_temp_directory(const std::string&);	//UTF8 encoding
		void				set_default_temp_directory();

	protected:
		data_ptr			m_data;
};

MATCL_SQLITE_EXPORT_MACRO 
connection create_connection(database db);

MATCL_SQLITE_EXPORT_MACRO 
database open_database(const std::string& name, 	//UTF8 encoding
						enums::open_mode om = enums::readwrite_create, 
						enums::thread_mode tm = enums::multi_thread);

}}