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
#pragma warning(disable:4251)	// needs to have dll-interface to be used by clients

#include "config.h"
#include "typedefs.h"

namespace matcl { namespace sql
{

class MATCL_SQLITE_EXPORT_MACRO charptr_owner
{
	private:
		char*			m_ptr;
		size_t			m_bytes;

	public:
		charptr_owner(char* ptr, size_t bytes)
            : m_ptr(ptr), m_bytes(bytes) 
        {};
		
        explicit charptr_owner(size_t bytes);

		~charptr_owner()						{ delete[] m_ptr; };

		char*			pointer()				{ return m_ptr; };
		const char*		pointer() const			{ return m_ptr; };
		size_t			size() const			{ return m_bytes; };

	private:
		charptr_owner(const charptr_owner& bl) = delete;
		charptr_owner& operator=(const charptr_owner& bl) = delete;
};

class MATCL_SQLITE_EXPORT_MACRO blob_owner
{
	private:
		void*			m_ptr;
		size_t			m_bytes;

	public:
		blob_owner(void* ptr, size_t bytes);
		explicit blob_owner(size_t bytes);

		~blob_owner()							{ delete[] m_ptr; };

		void*			pointer()				{ return m_ptr; };
		const void*		pointer() const			{ return m_ptr; };
		size_t			bytes() const			{ return m_bytes; };

        void*           release()               { void* ret = m_ptr; m_ptr = nullptr; return ret; };

	private:
		blob_owner(const blob_owner& bl) = delete;
		blob_owner& operator=(const blob_owner& bl) = delete;
};

class MATCL_SQLITE_EXPORT_MACRO database_blob
{
	private:
		details::connection_data_ptr	m_connection;
		sqlite3_blob*					m_blob;
		int								m_bytes;

	public:		
		~database_blob();

		int				bytes();
		void			read(void *buffer, int n_bytes, int start_pos);
		void			write(const void *source, int n_bytes, int start_pos);

	private:
		database_blob(details::connection_data_ptr, sqlite3_blob* );

		database_blob(const database_blob& bl) = delete;
		database_blob& operator=(const database_blob& bl) = delete;

		friend class connection;
};

}}