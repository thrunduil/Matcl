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

#include <string.h>
#include <sqlite3.h>
#include "matcl-sqlite-cpp/query.h"
#include "matcl-sqlite-cpp/exception.h"
#include "matcl-sqlite-cpp/blob.h"
#include "matcl-sqlite-cpp/functors.h"
#include "matcl-sqlite-cpp/connection.h"
#include "details/query_details.h"
#include "details/connection_details.h"

namespace matcl { namespace sql
{

namespace details
{
query_data::query_data(connection_data_ptr con)
	:	m_con(con),m_db(con->m_db),
		m_raw_stm(NULL),m_cmap(NULL),m_get_row_succeeded(false),m_column_number(0),
		m_old_result_flag(0), m_use_old_result_flag(false), m_need_reset(false)
{
	m_handle = con->m_handle.get();
};

query_data::~query_data()
{
	if (m_raw_stm)
		sqlite3_finalize(m_raw_stm);

    delete m_cmap;
};
};

query::query()
{};

query::query(sql m_sql)
:sql(m_sql)
{
	make_column_map();
	m_data->m_old_result_flag		= step(get_stm());
	m_data->m_use_old_result_flag	= true;
};

query::~query()
{};

query::query(const query& other)
:sql(other)
{};

query& query::operator=(const query& other)
{
	m_data = other.m_data;
	return *this;
};

int query::step(sqlite3_stmt* stm)
{
	int rc = sqlite3_step(stm);
	switch (rc)
	{
		case SQLITE_DONE:
		case SQLITE_ROW:
			return rc;
		default:
			throw_query_exception(rc);
			throw;
	}
};

void query::make_column_map()
{
	if (m_data->m_cmap)
		delete m_data->m_cmap;

	m_data->m_cmap = new details::query_data::column_map();

	int i = 0;
	for(;;)
	{
		const char *p = sqlite3_column_name(get_stm(), i);
		if (!p)
			break;

        (*m_data->m_cmap)[p] = i++;
	};

	m_data->m_column_number = sqlite3_column_count(get_stm());
};

int query::number_cols() const
{
	return m_data->m_column_number;
}

//=================================================================================
//							GET VALUE
//=================================================================================
int query::get_column_index(const std::string& x) const
{
	using iterator  = details::query_data::column_map::const_iterator;
	iterator pos = m_data->m_cmap->find(x);

	if (pos == m_data->m_cmap->end())
	{
		std::string msg = "column with name " + x + "does not exist";
		throw matcl_sqlite_exception(msg);		
	};
	return pos->second;
};

void query::check_column_index_is_valid(int x) const
{
	if (m_data->m_get_row_succeeded)		
	{
		if ( x < m_data->m_column_number )
		{
			return;
		}
		else
		{
			std::string msg = "column index exceeded number of columns";
			throw matcl_sqlite_exception(msg);
		};
	}
	else
	{
		std::string msg = "sql statement does not point to a valid row";
		throw matcl_sqlite_exception(msg);
	};
};

int query::bytes_at_col(const std::string& x)
{
	int index = get_column_index(x);
	return bytes_at_col(index);
};

int	query::bytes_at_col(int x)
{
	check_column_index_is_valid(x);
	return sqlite3_column_bytes(get_stm(), x);
};

enums::value_type query::type_at_col(const std::string& x)
{
	int index = get_column_index(x);
	return type_at_col(index);
};

enums::value_type query::type_at_col(int x)
{
	check_column_index_is_valid(x);
	int type = sqlite3_column_type(get_stm(), x);

	switch(type)
	{
		case SQLITE_INTEGER:
			return enums::integer_type;
		case SQLITE_FLOAT:
			return enums::float_type;
		case SQLITE_TEXT:
			return enums::text_type;
		case SQLITE_BLOB:
			return enums::blob_type;
		case SQLITE_NULL:
			return enums::null_type;
		default:
			std::string msg = "unknown column type";
			throw matcl_sqlite_exception(msg);
	};
};

std::string query::string_at_col(const std::string& x)
{
	int index = get_column_index(x);
	return string_at_col(index);
}

std::string query::string_at_col(int x)
{	
	check_column_index_is_valid(x);

	const unsigned char *tmp = sqlite3_column_text(get_stm(), x);
	if (!tmp)
	{
		return "";
	};
	int size = bytes_at_col(x);
	return std::string((const char*)tmp,size);
}

void query::string_at_col(int x, string_getter* functor)
{
	check_column_index_is_valid(x);

	const unsigned char *tmp = sqlite3_column_text(get_stm(), x);
	if (!tmp)
	{
		const char *tmp2 = "";
		functor->eval((const unsigned char*)tmp2,0);
	}
	else
	{
		int size = bytes_at_col(x);
		functor->eval(tmp,size);
	};
};

char_ptr query::charptr_at_col(const std::string& x)
{
	int index = get_column_index(x);
	return charptr_at_col(index);
};

char_ptr query::charptr_at_col(int x)
{
	check_column_index_is_valid(x);

	const unsigned char *tmp = sqlite3_column_text(get_stm(), x);
	int size = bytes_at_col(x);

	if (size == 0)
	{
		char_ptr cp(new charptr_owner(NULL,0));
		return cp;
	};

	char* ptr = new char[size];
	char_ptr cp(new charptr_owner(ptr,size));
	memcpy(ptr,tmp,size);

	return cp;
};

blob_ptr query::blob_at_col(const std::string& x)
{
	int index = get_column_index(x);
	return blob_at_col(index);
};

blob_ptr query::blob_at_col(int x)
{
	check_column_index_is_valid(x);

	const void *tmp = sqlite3_column_blob(get_stm(), x);
	int size = bytes_at_col(x);

	void* ptr = new char[size];
	blob_ptr bl(new blob_owner(ptr,size));
	memcpy(ptr,tmp,size);

	return bl;
};

void query::blob_at_col(int x, blob_getter* functor)
{
	check_column_index_is_valid(x);

	const void *tmp = sqlite3_column_blob(get_stm(), x);
	int size = bytes_at_col(x);
	functor->eval(tmp,size);
};

double query::double_at_col(const std::string& x)
{
	int index = get_column_index(x);
	return double_at_col(index);
}

double query::double_at_col(int x)
{
	check_column_index_is_valid(x);
	return sqlite3_column_double(get_stm(), x);
}

long query::long_at_col(const std::string& x)
{
	int index = get_column_index(x);
	return long_at_col(index);
}

long query::long_at_col(int x)
{
	check_column_index_is_valid(x);
	return sqlite3_column_int(get_stm(), x);
}

unsigned long query::ulong_at_col(const std::string& x)
{
	int index = get_column_index(x);
	return ulong_at_col(index);
}

unsigned long query::ulong_at_col(int x)
{
	check_column_index_is_valid(x);
	unsigned long l = sqlite3_column_int(get_stm(), x);
	return l;
}

int64_t query::int64_at_col(const std::string& x)
{
	int index = get_column_index(x);
	return int64_at_col(index);
}

int64_t query::int64_at_col(int x)
{
	check_column_index_is_valid(x);
	return sqlite3_column_int64(get_stm(), x);
}

uint64_t query::uint64_at_col(const std::string& x)
{
	int index = get_column_index(x);
	return uint64_at_col(index);
}

uint64_t query::uint64_at_col(int x)
{
	check_column_index_is_valid(x);
	uint64_t l = sqlite3_column_int64(get_stm(), x);
	return l;
}

bool query::step()
{
	m_data->m_get_row_succeeded		= false;
	int rc							= m_data->m_use_old_result_flag ? m_data->m_old_result_flag 
									: step(get_stm());
	m_data->m_use_old_result_flag	= false;

	switch (rc)
	{
	case SQLITE_DONE:
		return false;
	case SQLITE_ROW:
		m_data->m_get_row_succeeded = true;
		return true;
	default:
		throw_query_exception(rc);
		throw;
	}
}

void query::reset()
{
	sqlite3_reset(get_stm());
};

sqlite3_stmt* query::get_stm()
{ 
	return m_data->m_raw_stm; 
};

const sqlite3_stmt*	query::get_stm() const		
{ 
	return m_data->m_raw_stm; 
};

bool query::empty() const
{
	return !m_data;
};

}}
