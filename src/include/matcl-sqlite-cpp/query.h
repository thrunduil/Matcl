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
#include "sql.h"

namespace matcl { namespace sql
{


class MATCL_SQLITE_EXPORT_MACRO query : private sql
{
	using sql::m_data;

	public:		
		query();
		~query();

		query(const query& q);
		query& operator=(const query& );

		bool				empty() const;

		int					number_cols() const;

		int					bytes_at_col(const std::string& x);
		int					bytes_at_col(int x);

		enums::value_type	type_at_col(const std::string& x);
		enums::value_type	type_at_col(int x);

		std::string			string_at_col(const std::string& x);
		std::string			string_at_col(int x);
		void				string_at_col(int x, string_getter* functor);
		char_ptr			charptr_at_col(const std::string& x);
		char_ptr			charptr_at_col(int x);

		blob_ptr			blob_at_col(const std::string& x);
		blob_ptr			blob_at_col(int x);
		void				blob_at_col(int x, blob_getter* functor);

		double				double_at_col(const std::string& x);
		double				double_at_col(int x);

		long				long_at_col(const std::string& x);
		long				long_at_col(int x);

		unsigned long		ulong_at_col(const std::string& x);
		unsigned long		ulong_at_col(int x);

		int64_t				int64_at_col(const std::string& x);
		int64_t				int64_at_col(int x);

		uint64_t			uint64_at_col(const std::string& x);
		uint64_t			uint64_at_col(int x);

		bool				step();
		void				reset();
	
		sqlite3_stmt*		get_stm();
		const sqlite3_stmt*	get_stm() const;

		int					get_column_index(const std::string& x) const;
		void				check_column_index_is_valid(int x) const;		

	private:
		query(sql m_sql);	
		operator sql() const;

		void				make_column_map();
		int 				step(sqlite3_stmt*);

		friend class sql;
};

}}
