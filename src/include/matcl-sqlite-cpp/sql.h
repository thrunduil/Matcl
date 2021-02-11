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

class MATCL_SQLITE_EXPORT_MACRO sql
{
	protected:
		using data_ptr  = details::query_data_ptr;

	public:		
		sql();
		sql(data_ptr data, const std::string& sql);
		~sql();

		sql(const sql& q);
		sql& operator=(const sql& );

		bool				empty() const;

		int					get_free_parameter_index(const char* param_name);
		std::string			get_free_parameter_name(int index);
		int					get_number_of_free_parameters();

		void				bind_blob(int param_nr, const void* ptr, int size);
		void				bind_blob(const std::string& x, const void* ptr, int size);

		void				bind_double(int param_nr,double val);
		void				bind_double(const std::string& x,double val);

		void				bind_int(int param_nr,int val);
		void				bind_int(const std::string& x,int val);

		void				bind_int64(int param_nr,int64_t val);
		void				bind_int64(const std::string& x,int64_t val);

		void				bind_null(int param_nr);
		void				bind_null(const std::string& x);

		void				bind_text(int param_nr, const char* txt, int n_elem);
		void				bind_text(const std::string& x, const char* txt, int n_elem);

		void				bind_text(int param_nr, const char* txt);
		void				bind_text(const std::string& x, const char* txt);

		void				bind_text(int param_nr, const std::string&);
		void				bind_text(const std::string& x, const std::string&);

		void				bind_zero_blob(int param_nr, int n_bytes);
		void				bind_zero_blob(const std::string& x, int n_bytes);
		
		void				clear_all_bindings();
	
		sqlite3_stmt*		get_stm();
		const sqlite3_stmt*	get_stm() const;

		query				make();

	private:
		sql(const query& q);

		void				prepare(const std::string& sql);		
		void				reset();

	protected:
		void				throw_query_exception(int code);
		void				throw_query_exception(const std::string& msg,int code );		

	protected:
		data_ptr			m_data;

		friend class query;
};

namespace details
{
	sql create_sql(const details::connection_data_ptr& con,const std::string& sql_str);
};

}}