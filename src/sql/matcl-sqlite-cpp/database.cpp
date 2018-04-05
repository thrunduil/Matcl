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

#include "matcl-sqlite-cpp/database.h"
#include <sqlite3.h>
#include "matcl-sqlite-cpp/exception.h"
#include <list>
#include <set>
#include "matcl-sqlite-cpp/connection.h"
#include "matcl-sqlite-cpp/exception.h"
#include "details/database_details.h"
#include "details/connection_details.h"

#include <mutex>

extern char *sqlite3_temp_directory;

namespace matcl { namespace sql
{

//=========================================================================
//						PRIVATE CLASSES
//=========================================================================
namespace details
{

class handle_pool
{
	private:
		using handle_list   = std::list<db_handle *>;
		using handle_set    = std::set<db_handle *>;
        using mutex_type    = std::mutex;
		using scoped_lock   = std::unique_lock<mutex_type>;

		std::mutex	    m_mutex;
		handle_list		m_available_handles;
		handle_set		m_used_handles;
		database_data*	m_owner;
		scoped_lock		m_lock;

	public:
		handle_pool(database_data* owner);
		~handle_pool();

		db_handle * get_handle(const std::string& name, enums::open_mode om, enums::thread_mode tm);
		void		free_handle(db_handle* db);
		void		close_unused();
		int			number_of_connections();
		void		lock();
		void		unlock();

	private:
		void		close_unused_impl();
		bool		check_close_unused();
};

class database_manager
{
	private:
		using pool_set      = std::set<handle_pool*>;
        using mutex_type    = std::mutex;
		using scoped_lock   = std::unique_lock<mutex_type>;

		pool_set		m_set;
		std::mutex	    m_mutex;

	public:
		database_manager() {};
		~database_manager();
		void	register_handle_pool(handle_pool* p);
		void	unregister_handle_pool(handle_pool* p);

		void	set_temp_directory(const std::string& dir);

	private:

		void	lock_pools();
		void	unlock_pools();
		bool	are_open_connections();

		database_manager(const database_manager&);
		database_manager& operator=(const database_manager&);
};

static database_manager m_database_manager;

handle_pool::handle_pool(database_data* owner) 
    : m_owner(owner),m_lock(m_mutex,std::defer_lock) 
{
	m_database_manager.register_handle_pool(this);
}; 

handle_pool::~handle_pool()
{
	m_database_manager.unregister_handle_pool(this);
	
    {
		using iterator  = handle_list::iterator;

		iterator it_end = m_available_handles.end();
		for (iterator it = m_available_handles.begin(); it != it_end; it++)
		{
			db_handle *p = *it;
			delete p;
		}
	}
	{
		using iterator  = handle_set::iterator;

		iterator it_end = m_used_handles.end();
		for (iterator it = m_used_handles.begin(); it != it_end; it++)
		{
			db_handle *p = *it;
			delete p;
		}
	}			
};

db_handle *handle_pool::get_handle(const std::string& name, enums::open_mode om, enums::thread_mode tm)
{
	scoped_lock lock(m_mutex);

	if (m_available_handles.size() == 0)
	{
		db_handle * db = new db_handle(name,m_owner,om,tm);
		m_available_handles.push_back(db);
	};

	db_handle * db = m_available_handles.front();
	m_available_handles.pop_front();

	m_used_handles.insert(db);
	return db;
};

void handle_pool::free_handle(db_handle* db)
{
	scoped_lock lock(m_mutex);

	using iterator  = handle_set::iterator;

	iterator pos = m_used_handles.find(db);
	if (pos == m_used_handles.end())
	{
	}
	else
	{
		m_used_handles.erase(pos);
		m_available_handles.push_back(db);

		if (check_close_unused())
			close_unused_impl();
	};
};

void handle_pool::close_unused()
{
	scoped_lock lock(m_mutex);

	close_unused_impl();
};

void handle_pool::close_unused_impl()
{
	if (m_available_handles.size() <= 1)
		return;

	using iterator  = handle_list::iterator;
	iterator it_end = m_available_handles.end();
	iterator it = m_available_handles.begin();
	
	db_handle *first = *it;
	it++;

	for (; it != it_end; it++)
	{
		db_handle *p = *it;
		delete p;
	}

	m_available_handles.clear();	
	m_available_handles.push_back(first);
};

bool handle_pool::check_close_unused()
{
	size_t n_used = m_used_handles.size();
	size_t n_free = m_available_handles.size();

	if (n_used == 0 && n_free > 1)
		return true;

	if (n_free > 2*n_used)
		return true;

	return false;
};

int	handle_pool::number_of_connections()
{
	return (int)m_used_handles.size();
};

void handle_pool::lock()
{
	m_lock.lock();
}

void handle_pool::unlock()
{
	m_lock.unlock();
};

database_manager::~database_manager() 
{
	if (sqlite3_temp_directory)
	{
		sqlite3_free(sqlite3_temp_directory);
		sqlite3_temp_directory = NULL;
	}
};

void database_manager::register_handle_pool(handle_pool* p)
{
	scoped_lock lock(m_mutex);
	m_set.insert(p);
};

void database_manager::unregister_handle_pool(handle_pool* p)
{
	scoped_lock lock(m_mutex);
	m_set.erase(p);
};

void database_manager::set_temp_directory(const std::string& dir)
{
	scoped_lock lock(m_mutex);
	lock_pools();

	try
	{
		if (are_open_connections())
			throw matcl_sqlite_exception("cannot change temp dir if there are active connections");

        if (sqlite3_temp_directory)
		{
			sqlite3_free(sqlite3_temp_directory);
			sqlite3_temp_directory = NULL;
		}
		
        if (dir.size() == 0)
			return;

        size_t size = dir.size()+1;
		char* ptr = (char*)sqlite3_malloc(static_cast<int>(size*sizeof(char)));
		if (!ptr)
			throw std::bad_alloc();

		memcpy(ptr,dir.c_str(),size);
		sqlite3_temp_directory = ptr;
	}
	catch(...)
	{
		unlock_pools();
		throw;
	};

	unlock_pools();
};

void database_manager::lock_pools()
{
	using iterator  = pool_set::iterator;
	iterator pos = m_set.begin();

    while(pos!=m_set.end())
	{
		(*pos)->lock();
		pos++;
	};
};

void database_manager::unlock_pools()
{
	using iterator  = pool_set::iterator;
	iterator pos = m_set.begin();
	while(pos!=m_set.end())
	{
		(*pos)->unlock();
		pos++;
	};
};

bool database_manager::are_open_connections()
{
	using iterator  = pool_set::iterator;
	iterator pos = m_set.begin();
	while(pos!=m_set.end())
	{
		if ((*pos)->number_of_connections() > 0)
		{
			return true;
		};
		pos++;
	};
	return false;
};

//=========================================================================
//						DB_HANDLE
//=========================================================================

db_handle::db_handle(const std::string& name, database_data* owner, enums::open_mode om, 
    enums::thread_mode tm)
: m_db(NULL),m_owner(owner)
{
	int rc = sqlite3_open_v2(name.c_str(),&m_db, om | tm, NULL);
	std::string err_msg;
	if (rc || !m_db)
	{
		err_msg = last_error();
		sqlite3_close(m_db);

		throw matcl_sqlite_open_exception(name,err_msg);
	};
};

db_handle::~db_handle()
{
	sqlite3_close(m_db);	
};

void db_handle::free()
{ 
	m_owner->free_handle(this);
};

const char*	db_handle::last_error() const
{
	const char* msg = sqlite3_errmsg(m_db);
	return msg;
};

int64_t db_handle::last_insert_id() const
{
	return sqlite3_last_insert_rowid(m_db);
};

void lock_db_handle::set(db_handle* handle)
{ 
	m_handle = handle; 
};


//=========================================================================
//						DATABASE_DATA
//=========================================================================

database_data::database_data(const std::string& name, enums::open_mode om, enums::thread_mode tm)
:m_name(name),m_open_mode(om),m_thread_mode(tm)
{
	m_open_handles = new handle_pool(this);
}

database_data::~database_data()
{
	delete m_open_handles; 
}

void database_data::free_handle(db_handle* db)
{
	m_open_handles->free_handle(db);
};

};

//=========================================================================
//						DATABASE
//=========================================================================

database::database(details::database_data_ptr data)
:m_data(data)
{}

database::database(const database& other)
:m_data(other.m_data)
{};

database& database::operator=(const database& other)
{
	m_data = other.m_data;
	return *this;
};	

void database::get_db_handle(details::lock_db_handle& handle)
{
	details::db_handle *db = NULL;	
	db = m_data->m_open_handles->get_handle(m_data->m_name, m_data->m_open_mode, m_data->m_thread_mode);
	handle.set(db);
}

bool database::is_connected()
{
	details::lock_db_handle m_handle;
	try
	{
		get_db_handle(m_handle);
		return m_handle;
	}
	catch(...)
	{
		return false;
	};
}

void database::close_unused_connections()
{
	m_data->m_open_handles->close_unused();
};

void database::interrupt(connection connection)
{
	return sqlite3_interrupt(connection.get_db());
};

void database::set_temp_directory(const std::string& dir)
{
	details::m_database_manager.set_temp_directory(dir);
};

void database::set_default_temp_directory()
{
	set_temp_directory("");
};

int	database::number_of_connections()
{
	return m_data->m_open_handles->number_of_connections();
};

bool database::empty() const
{
	return !m_data;
};

//=========================================================================
//						FREE FUNCTIONS
//=========================================================================

connection create_connection(database db)
{
	return connection(details::connection_data_ptr( new details::connection_data(db)));	
};

database open_database(const std::string& name, enums::open_mode om, enums::thread_mode tm)
{
	details::database_data_ptr dbptr;
	dbptr = details::database_data_ptr(new details::database_data(name, om, tm));
	return database(dbptr);
};

}}
