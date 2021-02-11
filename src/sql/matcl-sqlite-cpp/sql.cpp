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

#include <sqlite3.h>
#include "matcl-sqlite-cpp/sql.h"
#include "matcl-sqlite-cpp/exception.h"
#include "matcl-sqlite-cpp/blob.h"
#include "matcl-sqlite-cpp/functors.h"
#include "matcl-sqlite-cpp/connection.h"
#include "details/query_details.h"
#include "details/connection_details.h"
#include "matcl-sqlite-cpp/query.h"

namespace matcl { namespace sql
{

sql::sql()
{};

sql::sql(details::query_data_ptr data, const std::string& sql)
:m_data(data)
{
	prepare(sql);
};

sql::sql(const query& q)
:m_data(q.m_data)
{};

sql::~sql()
{};

sql::sql(const sql& other)
:m_data(other.m_data)
{
};

sql& sql::operator=(const sql& other)
{
	m_data = other.m_data;
	return *this;
};

void sql::prepare(const std::string& sql)
{
	const char* s = NULL;
	int rc = sqlite3_prepare_v2(m_data->m_handle->get_db(), sql.c_str(), (int)sql.size(), &m_data->m_raw_stm, &s);
	if (rc != SQLITE_OK)
	{
		throw_query_exception(rc);
		throw;
	}
	if (!m_data->m_raw_stm)
	{
		throw matcl_sqlite_exception("query execution failed");
	};
};

query sql::make()
{	
	reset();
	m_data->m_need_reset = true;
	query q = query(*this);
	return q;
};

void sql::reset()
{
	if (m_data->m_need_reset)
	{
		sqlite3_reset(get_stm());
		m_data->m_need_reset = false;
	};
};

void sql::clear_all_bindings()
{
	sqlite3_clear_bindings(get_stm());
};

int sql::get_free_parameter_index(const char* param_name)
{
	int index = sqlite3_bind_parameter_index(get_stm(), param_name);
	if (index == 0)
	{
		std::string msg = std::string() + "parameter " + param_name 
						+ " is not a free variable is the query";
		throw matcl_sqlite_exception(msg);
	};
	return index;
};

std::string sql::get_free_parameter_name(int index)
{
	const char* name = sqlite3_bind_parameter_name(get_stm(), index);
	if (index == 0)
	{
		std::string msg = std::string() + "parameter index is out of range or given parameters is nameless";
		throw matcl_sqlite_exception(msg);
	};
	return name;
};

int	sql::get_number_of_free_parameters()
{
	return sqlite3_bind_parameter_count(get_stm());
};

void sql::bind_double(int param_nr,double val)
{
	reset();
	int rc = sqlite3_bind_double(get_stm(), param_nr, val);
	if (rc != SQLITE_OK)
	{
		throw_query_exception(rc);
		throw;
	};
};

void sql::bind_int(int param_nr,int val)
{
	reset();
	int rc = sqlite3_bind_int(get_stm(), param_nr, val);
	if (rc != SQLITE_OK)
	{
		throw_query_exception(rc);
		throw;
	};
};

void sql::bind_int64(int param_nr,int64_t val)
{
	reset();
	int rc = sqlite3_bind_int64(get_stm(), param_nr, val);
	if (rc != SQLITE_OK)
	{
		throw_query_exception(rc);
		throw;
	};
};

void sql::bind_null(int param_nr)
{
	reset();
	int rc = sqlite3_bind_null(get_stm(), param_nr);
	if (rc != SQLITE_OK)
	{
		throw_query_exception(rc);
		throw;
	};
};

void sql::bind_zero_blob(int param_nr, int n_bytes)
{
	reset();
	int rc = sqlite3_bind_zeroblob(get_stm(), param_nr, n_bytes);
	if (rc != SQLITE_OK)
	{
		throw_query_exception(rc);
		throw;
	};
};

void sql::bind_text(int param_nr, const std::string& txt)
{
	reset();
	return bind_text(param_nr,txt.c_str(),(int)txt.size());
};

void sql::bind_text(int param_nr, const char* txt, int n_elem)
{
	reset();
	int rc = sqlite3_bind_text(get_stm(), param_nr, txt,n_elem, SQLITE_TRANSIENT);
	if (rc != SQLITE_OK)
	{
		throw_query_exception(rc);
		throw;
	};
};

void sql::bind_text(int param_nr, const char* txt)
{
	reset();
	int rc = sqlite3_bind_text(get_stm(), param_nr, txt, -1, SQLITE_TRANSIENT);
	if (rc != SQLITE_OK)
	{
		throw_query_exception(rc);
		throw;
	};
};

void sql::bind_blob(int param_nr, const void* ptr, int size)
{
	reset();
	int rc = sqlite3_bind_blob(get_stm(), param_nr, ptr, size, SQLITE_TRANSIENT);
	if (rc != SQLITE_OK)
	{
		throw_query_exception(rc);
		throw;
	};
};


void sql::bind_double(const std::string& x,double val)
{
	return bind_double(get_free_parameter_index(x.c_str()),val);
};

void sql::bind_int(const std::string& x,int val)
{
	return bind_int(get_free_parameter_index(x.c_str()),val);
};

void sql::bind_int64(const std::string& x,int64_t val)
{
	return bind_int64(get_free_parameter_index(x.c_str()),val);
};

void sql::bind_null(const std::string& x)
{
	return bind_null(get_free_parameter_index(x.c_str()));
};

void sql::bind_zero_blob(const std::string& x, int n_bytes)
{
	return bind_zero_blob(get_free_parameter_index(x.c_str()),n_bytes);
};

void sql::bind_text(const std::string& x, const std::string& txt)
{
	return bind_text(get_free_parameter_index(x.c_str()),txt);
};

void sql::bind_text(const std::string& x, const char* txt, int n_elem)
{
	return bind_text(get_free_parameter_index(x.c_str()),txt,n_elem);
};

void sql::bind_text(const std::string& x, const char* txt)
{
	return bind_text(get_free_parameter_index(x.c_str()),txt);
};

void sql::bind_blob(const std::string& x, const void* ptr, int size)
{
	return bind_blob(get_free_parameter_index(x.c_str()),ptr,size);
};

sqlite3_stmt* sql::get_stm()
{ 
	return m_data->m_raw_stm; 
};

const sqlite3_stmt*	sql::get_stm() const		
{ 
	return m_data->m_raw_stm; 
};

void sql::throw_query_exception(int code)
{
	throw_query_exception(m_data->m_handle->last_error(),code);
};

void sql::throw_query_exception(const std::string& msg,int code )
{
	switch(code)
	{
	case SQLITE_BUSY:
	case SQLITE_LOCKED:
		throw matcl_sqlite_lock_exception(code);
	case SQLITE_INTERRUPT:
		throw matcl_sqlite_interrupt_exception(code);
	default:
		throw matcl_sqlite_query_exception(msg,code);
	};
};

bool sql::empty() const
{
	return !m_data;
};

namespace details
{

sql create_sql(const connection_data_ptr& con,const std::string& sql_txt)
{
	return sql(query_data_ptr(new query_data(con)),sql_txt);
};

};

}}