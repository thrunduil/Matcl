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

#include "matcl-file/matcl_file.h"
#include "matcl_file_data.h"
#include "matcl-file/exception.h"
#include "matcl-core/IO/archive.h"
#include "matcl-core/matrix/matrix_traits.h"
#include <istream>
#include "matcl-matrep/details/rangebuf.h"
#include "matcl-sqlite-cpp/builder/Select.h"
#include "matcl-sqlite-cpp/builder/insert.h"
#include "matcl-sqlite-cpp/builder/update.h"
#include "matcl-sqlite-cpp/builder/delete.h"
#include "matcl-sqlite-cpp/builder/create_table.h"
#include "matcl-sqlite-cpp/builder/create_index.h"

namespace matcl
{

static const sql::table matrix_table = sql::table("matrix_table");
static const sql::table data_table   = sql::table("data_table");

static const std::string insert_query = sql::insert().into(matrix_table).values( sql::expression_list()
        (sql::bind_param("mat_name"))
        (sql::bind_param("flags"))
        (sql::bind_param("rows"))
        (sql::bind_param("cols"))
        (sql::bind_param("ldiags"))
        (sql::bind_param("udiags"))
        (sql::bind_param("nnz"))
        (sql::bind_param("value_type"))
        (sql::bind_param("struct_type"))
        (sql::bind_param("data"))
        (sql::bind_param("mat_string"))
        (sql::bind_param("associated_data"))
    ).to_str();

using cb = sql::conflict_behaviour;

static const std::string replace_query = sql::insert().on_conflict(cb::replace).into(matrix_table)
                                                .values( sql::expression_list()
        (sql::bind_param("mat_name"))
        (sql::bind_param("flags"))
        (sql::bind_param("rows"))
        (sql::bind_param("cols"))
        (sql::bind_param("ldiags"))
        (sql::bind_param("udiags"))
        (sql::bind_param("nnz"))
        (sql::bind_param("value_type"))
        (sql::bind_param("struct_type"))
        (sql::bind_param("data"))
        (sql::bind_param("mat_string"))
        (sql::bind_param("associated_data"))
    ).to_str();

static const int col_name				= 1;
static const int col_flags				= 2;
static const int col_rows				= 3;
static const int col_cols				= 4;
static const int col_ldiags				= 5;
static const int col_udiags				= 6;
static const int col_nnz				= 7;
static const int col_value_type			= 8;
static const int col_struct_type		= 9;
static const int col_data				= 10;
static const int col_mat_string			= 11;
static const int col_associated_data	= 12;

static const std::string data_insert_query = sql::insert().into(data_table).values( sql::expression_list()
        (sql::bind_param("mat_name"))
        (sql::bind_param("data"))
        (sql::bind_param("mat_string"))
        (sql::bind_param("associated_data"))
    ).to_str();

static const std::string data_replace_query = sql::insert().on_conflict(cb::replace).into(data_table)
                                                    .values( sql::expression_list()
        (sql::bind_param("mat_name"))
        (sql::bind_param("data"))
        (sql::bind_param("mat_string"))
        (sql::bind_param("associated_data"))
    ).to_str();

static const int data_col_name				= 1;
static const int data_col_data				= 2;
static const int data_col_mat_string		= 3;
static const int data_col_associated_data   = 4;

static const std::string update_mat_string_query = 
        sql::update().on_conflict(cb::ignore).tab(sql::table("main.matrix_table"))
        .set(sql::column("mat_string") = sql::bind_param("update_value"))
        .where(sql::column("name") == sql::bind_param("update_mat_name")).to_str();

static const std::string update_data_string_query = 
        sql::update().on_conflict(cb::ignore).tab(sql::table("main.data_table"))
        .set(sql::column("mat_string") = sql::bind_param("update_value"))
        .where(sql::column("name") == sql::bind_param("update_mat_name")).to_str();


static const int val_index_update_value		= 1;
static const int val_index_update_mat_name	= 2;

blob_data::blob_data(void* ptr, size_t bytes)
: m_ptr(ptr), m_bytes(bytes)
{};
blob_data::blob_data(size_t bytes)
:m_bytes(bytes),m_ptr(NULL)
{	
	m_ptr = new char[bytes];	
};

static blob_ptr convert_to_matclblob(const sql::blob_ptr& bp)
{
    blob_ptr new_bp(new blob_data(bp->release(), bp->bytes()));
    return new_bp;
};

 mmlib_file::mmlib_file(const std::string& file_name, open_mode om, thread_mode tm)
:m_data(new details::mmlib_file_data(file_name, om, tm))
{
	sql::connection q = m_data->m_connection;
	q.timeout_limit(1000);

	if (om == open_mode::readwrite_create)
	{
		insert_main_table_if_not_exist();
        insert_data_table_if_not_exist();
	};
};
mmlib_file::~mmlib_file()
{};

void mmlib_file::timeout_limit(int msec)
{
	sql::connection q = m_data->m_connection;
	q.timeout_limit(msec);
};

void mmlib_file::insert_main_table_if_not_exist()
{
	sql::connection q = m_data->m_connection;
	try
	{
        q.execute(sql::select().from(matrix_table).to_str());
	}
	catch(sql::matcl_sqlite_lock_exception&)
	{
		throw error::error_mmlibfile_locked();
	}
	catch(sql::matcl_sqlite_exception&)
	{
		insert_main_table();
	};
};
bool mmlib_file::exist_data_table() const
{
	sql::connection q = m_data->m_connection;
	try
	{
        q.execute(sql::select().from(data_table).to_str());
	}
	catch(sql::matcl_sqlite_lock_exception&)
	{
		throw error::error_mmlibfile_locked();
	}
	catch(sql::matcl_sqlite_exception&)
	{
		return false;
	};

    return true;
};
void mmlib_file::insert_data_table_if_not_exist()
{
	sql::connection q = m_data->m_connection;
	try
	{
        q.execute(sql::select().from(data_table).to_str());
	}
	catch(sql::matcl_sqlite_lock_exception&)
	{
		throw error::error_mmlibfile_locked();
	}
	catch(sql::matcl_sqlite_exception&)
	{
		insert_data_table();
	};
};

void mmlib_file::insert_main_table()
{
    using ct = sql::column_type;

	sql::connection q = m_data->m_connection;
	try
	{
        q.execute(
            sql::create_table().if_not_exists().tab(matrix_table).cols(
                sql::column_decl_list()
                    (sql::column_decl("name").type(ct::text).unique())
                    (sql::column_decl("flags").type(ct::integer))
                    (sql::column_decl("n_rows").type(ct::integer))
                    (sql::column_decl("n_cols").type(ct::integer))
                    (sql::column_decl("ldiags").type(ct::integer))
                    (sql::column_decl("udiags").type(ct::integer))
                    (sql::column_decl("nnz").type(ct::integer))
                    (sql::column_decl("value_type").type(ct::integer))
                    (sql::column_decl("struct_type").type(ct::integer))
                    (sql::column_decl("data").type(ct::blob))
                    (sql::column_decl("mat_string").type(ct::text))
                    (sql::column_decl("associated_data").type(ct::blob))
                ).to_str());

        q.execute(
            sql::create_index().
            if_not_exists().
            name("mat_name").
            on(matrix_table).
            cols(sql::column_list()(sql::column("name"))).to_str());
	}
	catch(sql::matcl_sqlite_lock_exception&)
	{
		throw error::error_mmlibfile_locked();
	}
	catch(sql::matcl_sqlite_exception& ex)
	{
		throw error::error_create_mmlibfile(ex.sqlite_message());
	};
}

void mmlib_file::insert_data_table()
{
	sql::connection q = m_data->m_connection;
    using ct = sql::column_type;

	try
	{
        q.execute(
            sql::create_table().if_not_exists().tab(data_table).cols(
                sql::column_decl_list()
                    (sql::column_decl("name").type(ct::text).unique())
                    (sql::column_decl("data").type(ct::blob))
                    (sql::column_decl("mat_string").type(ct::text))
                    (sql::column_decl("associated_data").type(ct::blob))
                ).to_str());

        q.execute(
            sql::create_index().
            if_not_exists().
            name("data_name").
            on(data_table).
            cols(sql::column_list()(sql::column("name"))).to_str());
	}
	catch(sql::matcl_sqlite_lock_exception&)
	{
		throw error::error_mmlibfile_locked();
	}
	catch(sql::matcl_sqlite_exception& ex)
	{
		throw error::error_create_mmlibfile(ex.sqlite_message());
	};
}

bool mmlib_file::exist(const std::string& mat_name)
{
	sql::connection q = m_data->m_connection;
	try
	{
		sql::query c = q.execute(sql::select().from(matrix_table).where(sql::column("name") == mat_name).to_str());
		return (c.step())? true : false;
	}
	catch(sql::matcl_sqlite_lock_exception&)
	{
		throw error::error_mmlibfile_locked();
	}
	catch(sql::matcl_sqlite_exception& ex)
	{
		throw error::error_read_mmlibfile(ex.sqlite_message());
	};	
};

bool mmlib_file::exist_data(const std::string& mat_name)
{
	sql::connection q = m_data->m_connection;
	try
	{
        bool et = this->exist_data_table();
        if (et == false)
            return false;

		sql::query c = q.execute(sql::select().from(data_table).where(sql::column("name") == mat_name).to_str());
		return (c.step())? true : false;
	}
	catch(sql::matcl_sqlite_lock_exception&)
	{
		throw error::error_mmlibfile_locked();
	}
	catch(sql::matcl_sqlite_exception& ex)
	{
		throw error::error_read_mmlibfile(ex.sqlite_message());
	};	
};

Matrix mmlib_file::load(const std::string& mat_name)
{
	sql::blob_ptr bp;
	try
	{
		sql::connection q = m_data->m_connection;
		sql::query c = q.execute(sql::select().from(matrix_table).where(sql::column("name") == mat_name).to_str());
		if (!c.step())
		{
			std::string msg = "matrix " + mat_name + " does not exist";
			throw error::error_read_mmlibfile(msg);
		};

		bp = c.blob_at_col(col_data-1);
	}
	catch(sql::matcl_sqlite_lock_exception&)
	{
		throw error::error_mmlibfile_locked();
	}
	catch(sql::matcl_sqlite_exception& ex)
	{
		throw error::error_read_mmlibfile(ex.sqlite_message());
	};	

	const char* first = (const char*)bp->pointer();
	const char* last = first + bp->bytes();

	matcl::details::rangebuf<const char*> buf(first,last);
    #ifdef _MSC_VER
        std::istream is(&buf,std::ios::binary);
    #else
	    std::istream is(&buf);//,std::ios::binary);
    #endif
	matcl::iarchive ia(is);

	Matrix mat;
	matcl::load(ia, mat);
	//ia >> mat;
	return mat;
};

blob_ptr mmlib_file::load_data(const std::string& mat_name)
{
	sql::blob_ptr bp;
	try
	{
        bool et = this->exist_data_table();
        if (et == false)
        {
			std::string msg = "data " + mat_name + " does not exist";
			throw error::error_read_mmlibfile(msg);
        };

		sql::connection q = m_data->m_connection;
		sql::query c = q.execute(sql::select().from(data_table).where(sql::column("name") == mat_name).to_str());
		if (!c.step())
		{
			std::string msg = "data " + mat_name + " does not exist";
			throw error::error_read_mmlibfile(msg);
		};

		bp = c.blob_at_col(data_col_data-1);
	}
	catch(sql::matcl_sqlite_lock_exception&)
	{
		throw error::error_mmlibfile_locked();
	}
	catch(sql::matcl_sqlite_exception& ex)
	{
		throw error::error_read_mmlibfile(ex.sqlite_message());
	};	

	return convert_to_matclblob(bp);
};

Matrix mmlib_file::load(const std::string& mat_name,std::string& mat_string)
{
	sql::blob_ptr bp;
	try
	{
		sql::connection q = m_data->m_connection;
		sql::query c = q.execute(sql::select().from(matrix_table).where(sql::column("name") == mat_name).to_str());
		if (!c.step())
		{
			std::string msg = "matrix " + mat_name + " does not exist";
			throw error::error_read_mmlibfile(msg);
		};

		bp = c.blob_at_col(col_data - 1);
		mat_string = c.string_at_col(col_mat_string - 1);
	}
	catch(sql::matcl_sqlite_lock_exception&)
	{
		throw error::error_mmlibfile_locked();
	}
	catch(sql::matcl_sqlite_exception& ex)
	{
		throw error::error_read_mmlibfile(ex.sqlite_message());
	};	

	const char* first = (const char*)bp->pointer();
	const char* last = first + bp->bytes();

	matcl::details::rangebuf<const char*> buf(first,last);
    #ifdef _MSC_VER
        std::istream is(&buf,std::ios::binary);
    #else
	    std::istream is(&buf);//,std::ios::binary);
    #endif
	matcl::iarchive ia(is);

	Matrix mat;
	matcl::load(ia, mat);
	//ia >> mat;
	return mat;
};

blob_ptr mmlib_file::load_data(const std::string& mat_name,std::string& mat_string)
{
	sql::blob_ptr bp;
	try
	{
        bool et = this->exist_data_table();
        if (et == false)
        {
			std::string msg = "data " + mat_name + " does not exist";
			throw error::error_read_mmlibfile(msg);
        };

		sql::connection q = m_data->m_connection;
		sql::query c = q.execute(sql::select().from(data_table).where(sql::column("name") == mat_name).to_str());
		if (!c.step())
		{
			std::string msg = "data " + mat_name + " does not exist";
			throw error::error_read_mmlibfile(msg);
		};

		bp = c.blob_at_col(data_col_data - 1);
		mat_string = c.string_at_col(data_col_mat_string - 1);
	}
	catch(sql::matcl_sqlite_lock_exception&)
	{
		throw error::error_mmlibfile_locked();
	}
	catch(sql::matcl_sqlite_exception& ex)
	{
		throw error::error_read_mmlibfile(ex.sqlite_message());
	};	

	return convert_to_matclblob(bp);
};

std::string mmlib_file::load_mat_string(const std::string& mat_name)
{
	try
	{
		sql::connection q = m_data->m_connection;
		sql::query c = q.execute(sql::select().from(matrix_table).where(sql::column("name") == mat_name).to_str());
		if (!c.step())
		{
			std::string msg = "matrix " + mat_name + " does not exist";
			throw error::error_read_mmlibfile(msg);
		};

		return c.string_at_col(col_mat_string - 1);
	}
	catch(sql::matcl_sqlite_lock_exception&)
	{
		throw error::error_mmlibfile_locked();
	}
	catch(sql::matcl_sqlite_exception& ex)
	{
		throw error::error_read_mmlibfile(ex.sqlite_message());
	};	
};

std::string mmlib_file::load_data_string(const std::string& mat_name)
{
	try
	{
        bool et = this->exist_data_table();
        if (et == false)
        {
			std::string msg = "data " + mat_name + " does not exist";
			throw error::error_read_mmlibfile(msg);
        };

		sql::connection q = m_data->m_connection;
		sql::query c = q.execute(sql::select().from(data_table).where(sql::column("name") == mat_name).to_str());
		if (!c.step())
		{
			std::string msg = "data " + mat_name + " does not exist";
			throw error::error_read_mmlibfile(msg);
		};

		return c.string_at_col(data_col_mat_string - 1);
	}
	catch(sql::matcl_sqlite_lock_exception&)
	{
		throw error::error_mmlibfile_locked();
	}
	catch(sql::matcl_sqlite_exception& ex)
	{
		throw error::error_read_mmlibfile(ex.sqlite_message());
	};	
};

mmlib_file::matrix_list	mmlib_file::load_all()
{
	matrix_list ml;
	
	try
	{
		sql::connection q = m_data->m_connection;
		sql::query c = q.execute(sql::select().from(matrix_table).to_str());
		while (c.step())
		{
			sql::blob_ptr bp = c.blob_at_col(col_data - 1);

			const char* first = (const char*)bp->pointer();
			const char* last = first + bp->bytes();

			matcl::details::rangebuf<const char*> buf(first,last);
            #ifdef _MSC_VER
                std::istream is(&buf,std::ios::binary);
            #else
	            std::istream is(&buf);//,std::ios::binary);
            #endif
			matcl::iarchive ia(is);

			Matrix mat;
			matcl::load(ia, mat);
			//ia >> mat;
			ml.push_back(mat);
		};		
	}
	catch(sql::matcl_sqlite_lock_exception&)
	{
		throw error::error_mmlibfile_locked();
	}
	catch(sql::matcl_sqlite_exception& ex)
	{
		throw error::error_read_mmlibfile(ex.sqlite_message());
	};	

	return ml;
};

mmlib_file::data_list mmlib_file::load_data_all()
{
	data_list ml;
	
	try
	{
        bool et = this->exist_data_table();
        if (et == false)
        {
			return ml;
        };

		sql::connection q = m_data->m_connection;
		sql::query c = q.execute(sql::select().from(data_table).to_str());
		while (c.step())
		{
			sql::blob_ptr bp = c.blob_at_col(data_col_data - 1);
			ml.push_back(convert_to_matclblob(bp));
		};		
	}
	catch(sql::matcl_sqlite_lock_exception&)
	{
		throw error::error_mmlibfile_locked();
	}
	catch(sql::matcl_sqlite_exception& ex)
	{
		throw error::error_read_mmlibfile(ex.sqlite_message());
	};	

	return ml;
};

void mmlib_file::load_all_mat_string(mmlib_file::string_list& sl)
{
	try
	{
		sql::connection q = m_data->m_connection;
		sql::query c = q.execute(sql::select().from(matrix_table).to_str());
		while (c.step())
		{
			sl.push_back(c.string_at_col(col_mat_string - 1));
		};		
	}
	catch(sql::matcl_sqlite_lock_exception&)
	{
		throw error::error_mmlibfile_locked();
	}
	catch(sql::matcl_sqlite_exception& ex)
	{
		throw error::error_read_mmlibfile(ex.sqlite_message());
	};	
};

void mmlib_file::load_data_all_mat_string(mmlib_file::string_list& sl)
{
	try
	{
        bool et = this->exist_data_table();
        if (et == false)
        {
			return;
        };

		sql::connection q = m_data->m_connection;
		sql::query c = q.execute(sql::select().from(data_table).to_str());
		while (c.step())
		{
			sl.push_back(c.string_at_col(data_col_mat_string - 1));
		};		
	}
	catch(sql::matcl_sqlite_lock_exception&)
	{
		throw error::error_mmlibfile_locked();
	}
	catch(sql::matcl_sqlite_exception& ex)
	{
		throw error::error_read_mmlibfile(ex.sqlite_message());
	};	
};

matrix_info	mmlib_file::load_info(const std::string& mat_name)
{
	sql::blob_ptr bp;
	try
	{
		sql::connection q = m_data->m_connection;
		sql::query c = q.execute(sql::select().from(matrix_table).where(sql::column("name") == mat_name).to_str());
		if (!c.step())
		{
			throw error::error_read_mmlibfile_mat_not_exist(mat_name);
		};

		return get_mat_info(&c);
	}
	catch(sql::matcl_sqlite_lock_exception&)
	{
		throw error::error_mmlibfile_locked();
	}
	catch(sql::matcl_sqlite_exception& ex)
	{
		throw error::error_read_mmlibfile(ex.sqlite_message());
	};	
};

matrix_info	mmlib_file::get_mat_info(void* ptr)
{
	sql::query* c = static_cast<sql::query*>(ptr);

	matrix_info mi;
	mi.m_name			= c->string_at_col(col_name-1);
	mi.m_flags			= c->long_at_col(col_flags-1);
	mi.m_rows			= c->long_at_col(col_rows-1);
	mi.m_cols			= c->long_at_col(col_cols-1);
	mi.m_ldiags			= c->long_at_col(col_ldiags-1);
	mi.m_udiags			= c->long_at_col(col_udiags-1);
	mi.m_nnz			= c->long_at_col(col_nnz-1);
	long vt				= c->long_at_col(col_value_type-1);
	long st				= c->long_at_col(col_struct_type-1);

	matcl::value_code evt;
	switch ((matcl::value_code)vt)
	{
		case value_code::v_integer:
			evt = value_code::v_integer;
			break;
		case value_code::v_float:
			evt = value_code::v_float;
			break;
		case value_code::v_real:
			evt = value_code::v_real;
			break;
		case value_code::v_float_complex:
			evt = value_code::v_float_complex;
			break;
		case value_code::v_complex:
			evt = value_code::v_complex;
			break;
		case value_code::v_object:
			evt = value_code::v_object;
			break;
		default:
			throw error::unable_to_read_matrix();
	};

	matcl::struct_code est;
	switch ((matcl::struct_code)st)
	{
		case struct_code::struct_dense:
			est = struct_code::struct_dense;
			break;
		case struct_code::struct_sparse:
			est = struct_code::struct_sparse;
			break;
		case struct_code::struct_banded:
			est = struct_code::struct_banded;
			break;
		case struct_code::struct_scalar:
			est = struct_code::struct_scalar;
			break;
		default:
			throw error::unable_to_read_matrix();
	};

	mi.m_struct_type		= est;
	mi.m_matrix_type		= matrix_traits::get_matrix_type(evt,est);

	return mi;
};

void mmlib_file::load_all_info(mmlib_file::matrix_info_list& ml)
{
	try
	{
		sql::connection q = m_data->m_connection;
		sql::query c = q.execute(sql::select().from(matrix_table).to_str());
		while (c.step())
		{
			matrix_info mi = get_mat_info(&c);

			ml.push_back(mi);
		};		
	}
	catch(sql::matcl_sqlite_lock_exception&)
	{
		throw error::error_mmlibfile_locked();
	}
	catch(sql::matcl_sqlite_exception& ex)
	{
		throw error::error_read_mmlibfile(ex.sqlite_message());
	};	
};

void mmlib_file::remove(const std::string& mat_name)
{
	try
	{
		sql::connection q = m_data->m_connection;
        sql::query c = q.execute(sql::delete_stmt().tab(matrix_table).where(sql::column("name") == mat_name).to_str());
	}
	catch(sql::matcl_sqlite_lock_exception&)
	{
		throw error::error_mmlibfile_locked();
	}
	catch(sql::matcl_sqlite_exception& ex)
	{
		throw error::error_read_mmlibfile(ex.sqlite_message());
	};	
};

void mmlib_file::remove_data(const std::string& mat_name)
{
	try
	{
        bool et = this->exist_data_table();
        if (et == false)
        {
			return;
        };

		sql::connection q = m_data->m_connection;
        sql::query c = q.execute(sql::delete_stmt().tab(data_table).where(sql::column("name") == mat_name).to_str());
	}
	catch(sql::matcl_sqlite_lock_exception&)
	{
		throw error::error_mmlibfile_locked();
	}
	catch(sql::matcl_sqlite_exception& ex)
	{
		throw error::error_read_mmlibfile(ex.sqlite_message());
	};	
};

void mmlib_file::save(const Matrix& mat, const std::string& mat_name, const std::string& mat_string,
                      bool allow_replace)
{
	bool ex = exist(mat_name);
	int flags = 0;

	if (ex)
	{
		if (allow_replace == false)
		{
			throw error::error_write_mmlibfile_mat_already_exist(mat_name);
		}
	};

	try
	{
		sql::connection q = m_data->m_connection;
		sql::sql c;
		if(allow_replace)
		{
			c	= q.prepare(replace_query);
		}
		else
		{
			c	= q.prepare(insert_query);
		};

		std::ostringstream ss(std::ios::binary);
		matcl::oarchive ia(ss);
		matcl::save(ia, mat);
		//ia << mat;
		const std::string& str = ss.str();

		const void* data_ptr	= str.c_str();
		size_t data_bytes		= str.size();

		void* adata_ptr			= 0;
		size_t adata_bytes		= 0;

		c.bind_text(col_name,mat_name);
		c.bind_int(col_flags,flags);
		c.bind_int(col_rows,mat.rows());
		c.bind_int(col_cols,mat.cols());
		c.bind_int(col_ldiags,mat.structural_ldiags(true));
		c.bind_int(col_udiags,mat.structural_udiags(true));
		c.bind_int(col_nnz,mat.structural_nnz());
		c.bind_int(col_value_type,(int)mat.get_value_code());
		c.bind_int(col_struct_type,(int)mat.get_struct_code());
		c.bind_blob(col_data,data_ptr,(int)data_bytes);
		c.bind_text(col_mat_string,mat_string);
		c.bind_blob(col_associated_data,adata_ptr,(int)adata_bytes);

		c.make();
	}
	catch(sql::matcl_sqlite_lock_exception&)
	{
		throw error::error_mmlibfile_locked();
	}
	catch(sql::matcl_sqlite_exception& ex)
	{
		throw error::error_read_mmlibfile(ex.sqlite_message());
	};	
};

void mmlib_file::save_data(const void* data, size_t bytes, const std::string& mat_name, 
                           const std::string& mat_string, bool allow_replace)
{
	bool ex = exist_data(mat_name);

	if (ex)
	{
		if (allow_replace == false)
		{
			throw error::error_write_mmlibfile_mat_already_exist(mat_name);
		}
	};

	try
	{
		sql::connection q = m_data->m_connection;
		sql::sql c;
		if(allow_replace)
		{
			c	= q.prepare(data_replace_query);
		}
		else
		{
			c	= q.prepare(data_insert_query);
		};

		const void* data_ptr	= data;
		size_t data_bytes		= bytes;

		void* adata_ptr			= 0;
		size_t adata_bytes		= 0;

		c.bind_text(data_col_name,mat_name);
		c.bind_blob(data_col_data,data_ptr,(int)data_bytes);
		c.bind_text(data_col_mat_string,mat_string);
		c.bind_blob(data_col_associated_data,adata_ptr,(int)adata_bytes);

		c.make();
	}
	catch(sql::matcl_sqlite_lock_exception&)
	{
		throw error::error_mmlibfile_locked();
	}
	catch(sql::matcl_sqlite_exception& ex)
	{
		throw error::error_read_mmlibfile(ex.sqlite_message());
	};	
};

void mmlib_file::add_to_save_list(const Matrix& mat, const std::string& mat_name, 
                                  const std::string& mat_string, bool allow_replace)
{
	m_data->add_to_save_list(mat,mat_name,mat_string,allow_replace);
};

void mmlib_file::add_data_to_save_list(const void* data, size_t bytes, const std::string& mat_name, 
                                       const std::string& mat_string, bool allow_replace)
{
	m_data->add_data_to_save_list(data,bytes,mat_name,mat_string,allow_replace);
};

void mmlib_file::save_list()
{
	try
	{
		using save_mat_list     = details::mmlib_file_data::save_mat_list;
        using save_data_list    = details::mmlib_file_data::save_data_list;
		using mat_iterator      = save_mat_list::const_iterator;
        using data_iterator     = save_data_list::const_iterator;

		mat_iterator mat_pos    = m_data->m_save_mat_list.begin();
		mat_iterator mat_end    = m_data->m_save_mat_list.end();
		data_iterator data_pos  = m_data->m_save_data_list.begin();
		data_iterator data_end  = m_data->m_save_data_list.end();

		sql::connection q = m_data->m_connection;

		sql::transaction tr(q);

		while (mat_pos != mat_end)
		{			
			sql::sql c;
			if (mat_pos->allow_replace)
			{
				c = q.prepare(replace_query);
			}
			else
			{
				c = q.prepare(insert_query);
			};

			Matrix mat = mat_pos->m_matrix;

			std::ostringstream ss(std::ios::binary);
			matcl::oarchive ia(ss);
			matcl::save(ia, mat);
			//ia << mat;
			const std::string& str = ss.str();

			const void* data_ptr	= str.c_str();
			size_t data_bytes		= str.size();

			void* adata_ptr			= 0;
			size_t adata_bytes		= 0;
			int flags				= 0;

			c.bind_text(col_name,mat_pos->m_name);
			c.bind_int(col_flags,flags);
			c.bind_int(col_rows,mat.rows());
			c.bind_int(col_cols,mat.cols());
			c.bind_int(col_ldiags,mat.structural_ldiags(true));
			c.bind_int(col_udiags,mat.structural_udiags(true));
			c.bind_int(col_nnz,mat.structural_nnz());
			c.bind_int(col_value_type,(int)mat.get_value_code());
			c.bind_int(col_struct_type,(int)mat.get_struct_code());
			c.bind_blob(col_data,data_ptr,(int)data_bytes);
			c.bind_text(col_mat_string,mat_pos->m_mat_string);
			c.bind_blob(col_associated_data,adata_ptr,(int)adata_bytes);

			tr.add_task(c);

			mat_pos++;
		};

		while (data_pos != data_end)
		{			
			sql::sql c;

			if (data_pos->allow_replace)
			{
				c = q.prepare(data_replace_query);
			}
			else
			{
				c = q.prepare(data_insert_query);
			};

			const void* data_ptr	= data_pos->m_data;
			size_t data_bytes		= data_pos->m_bytes;

			void* adata_ptr			= 0;
			size_t adata_bytes		= 0;

			c.bind_text(data_col_name,data_pos->m_name);
			c.bind_blob(data_col_data,data_ptr,(int)data_bytes);
			c.bind_text(data_col_mat_string,data_pos->m_mat_string);
			c.bind_blob(data_col_associated_data,adata_ptr,(int)adata_bytes);

			tr.add_task(c);

			data_pos++;
		};

		tr.make();
	}
	catch(sql::matcl_sqlite_lock_exception&)
	{
		throw error::error_mmlibfile_locked();
	}
	catch(sql::matcl_sqlite_exception& ex)
	{
		throw error::error_read_mmlibfile(ex.sqlite_message());
	};	
}

void mmlib_file::modify_mat_string(const std::string& mat_name, const std::string& mat_string)
{
	try
	{
		sql::connection q = m_data->m_connection;
		sql::sql c = q.prepare(update_mat_string_query);

		c.bind_text(val_index_update_value,mat_string);
		c.bind_text(val_index_update_mat_name,mat_name);

		c.make();
	}
	catch(sql::matcl_sqlite_lock_exception&)
	{
		throw error::error_mmlibfile_locked();
	}
	catch(sql::matcl_sqlite_exception& ex)
	{
		throw error::error_read_mmlibfile(ex.sqlite_message());
	};	
};

void mmlib_file::modify_data_mat_string(const std::string& mat_name, const std::string& mat_string)
{
	try
	{
		sql::connection q = m_data->m_connection;
		sql::sql c = q.prepare(update_data_string_query);

		c.bind_text(val_index_update_value,mat_string);
		c.bind_text(val_index_update_mat_name,mat_name);

		c.make();
	}
	catch(sql::matcl_sqlite_lock_exception&)
	{
		throw error::error_mmlibfile_locked();
	}
	catch(sql::matcl_sqlite_exception& ex)
	{
		throw error::error_read_mmlibfile(ex.sqlite_message());
	};	
};

};
