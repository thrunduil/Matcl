/*
 *  This file is a part of Matrix Computation Library (MATCL)
 *
 *  Copyright (c) Pawe³ Kowal 2017
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

#include "matcl-dynamic/details/object_data_pool.h"
#include <boost/pool/pool.hpp>

namespace matcl { namespace dynamic
{

MATCL_DYN_EXPORT 
details::item_counter::ATOMIC_LONG details::item_counter::n_item(long(0));

details::object_data_pool_base::object_data_pool_base(size_t size)
{
    m_pool = new pool_type(size);
};

void details::object_data_pool_base::free(void* ptr)
{
    m_pool->free(ptr);
};

void* details::object_data_pool_base::malloc()
{
    return m_pool->malloc();
};

details::object_data_pool_base::~object_data_pool_base()
{};

};};
