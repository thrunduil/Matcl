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

#include "matcl-sqlite-cpp/transaction.h"
#include "matcl-sqlite-cpp/query.h"

namespace matcl { namespace sql
{

transaction::transaction()
{};

transaction::transaction(connection con)
:m_connection(con)
{};

transaction::~transaction()
{};

transaction::transaction(const transaction& q)
:m_connection(q.m_connection),m_task_list(q.m_task_list)
{};

transaction& transaction::operator=(const transaction& q)
{
	m_connection	= q.m_connection;
	m_task_list		= q.m_task_list;
	return *this;
};

bool transaction::empty() const
{
	return m_connection.empty();
};

void transaction::add_task(sql m_sql)
{
	m_task_list.push_back(m_sql);
};

void transaction::clear_tasks()
{
	m_task_list.clear();
};

void transaction::make()
{
	if (m_task_list.empty())
		return;

    m_connection.execute("begin transaction");
	
    try
	{
		using iterator  = task_list::iterator;
		iterator pos = m_task_list.begin();
		iterator end = m_task_list.end();
		while(pos!=end)
		{
			pos->make();
			pos++;
		};		
		m_connection.execute("end transaction");
		m_task_list.clear();
	}
	catch(...)
	{
		m_task_list.clear();
		m_connection.execute("rollback");
		throw;
	};
};

}}