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

#include <string>

#include "matcl-sqlite-cpp/config.h"
#include "matcl-sqlite-cpp/typedefs.h"

namespace matcl { namespace sql { namespace details
{

class db_handle 
{
	public:
		db_handle(const std::string& name, database_data* owner, enums::open_mode om, 
            enums::thread_mode tm );
		~db_handle();

		const char*		    last_error() const;
		int64_t				last_insert_id() const;

		sqlite3 *			get_db()				{ return m_db; };
		const sqlite3 *		get_db() const			{ return m_db; };

	private:
		db_handle(const db_handle&) = delete;
		db_handle& operator=(const db_handle&) = delete;

		void				free();

	private:
		sqlite3 *			m_db;
		database_data*		m_owner;

		friend class lock_db_handle;
};

class lock_db_handle 
{
	public:
		lock_db_handle()							: m_handle(NULL) {};
		~lock_db_handle()							{ if (m_handle) m_handle->free(); };

		db_handle*			get()					{ return m_handle; };
		const db_handle*	get() const				{ return m_handle; };
		db_handle*			operator->()			{ return m_handle; };
		const db_handle*	operator->() const		{ return m_handle; };

							operator bool()	const	{ return m_handle? 1 : 0; };
		bool				operator!() const		{ return !m_handle; };

	private:
		void				set(db_handle* handle);

		lock_db_handle(const lock_db_handle&) = delete;
		lock_db_handle& operator=(const lock_db_handle&) = delete;

		db_handle*			m_handle;		

        friend database;
};

class database_data
{
	public:		
		~database_data();

		std::string			m_name;
		handle_pool*		m_open_handles;
		enums::open_mode	m_open_mode;
		enums::thread_mode	m_thread_mode;

	public:
		void				free_handle(db_handle*);		

	private:
		database_data(const std::string& name, enums::open_mode om, enums::thread_mode tm);

		database_data(const database_data& ) = delete;
		database_data& operator=(const database_data& ) = delete;	

		friend MATCL_SQLITE_EXPORT_MACRO 
        database matcl::sql::open_database(const std::string&, enums::open_mode, enums::thread_mode);
};

};};}
