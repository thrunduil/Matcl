/*
 *  This file is a part of Matrix Computation Library (MATCL)
 *
 *  Copyright (c) Pawe³ Kowal 2017 - 2021
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

#include "type_table_cache.inl"

namespace matcl { namespace dynamic { namespace details
{

//--------------------------------------------------------------------------
//                         type_table_cache
//--------------------------------------------------------------------------
type_table_cache::type_table_cache()
    :m_unifiers(false), m_overloads(false), m_template_overloads(false)
    ,m_convert(false), m_assign(false)
{};

void type_table_cache::clear()
{
    m_unifiers.clear();
    m_overloads.clear();
    m_convert.clear();
    m_assign.clear();
    m_template_overloads.clear();
    m_last_call.clear();
};

void type_table_cache::clear_global()
{
    clear();
}

void type_table_cache::close_global()
{    
    delete this;
}

};};};
