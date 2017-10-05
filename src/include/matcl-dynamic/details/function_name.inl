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
#pragma once

#include "matcl-dynamic/function_name.h"

namespace matcl { namespace dynamic
{

inline size_t identifier::hash_value() const
{ 
    return m_hash; 
};

inline size_t identifier::get_unique_code() const
{
    return m_code;
};

inline function_validator::function_validator()
    : m_validator(nullptr), m_printer(nullptr) 
{};

inline function_validator::function_validator(validator v, printer p)
    : m_validator(v), m_printer(p) 
{};

inline function_validator::validator 
function_validator::get_validator() const
{ 
    return m_validator; 
}

inline function_validator::printer 
function_validator::get_printer() const
{ 
    return m_printer; 
};

};};
