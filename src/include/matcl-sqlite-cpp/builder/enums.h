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

#include "matcl-sqlite-cpp/config.h"

namespace matcl { namespace sql
{

enum class column_type
{
    no_type,
    integer,
    text,
    blob,
    numeric,
    real
};

enum class conflict_behaviour
{
    default_cb,
    rollback,
    abort,
    replace,
    fail,
    ignore
};

enum class collation
{
    no_collation,
    binary,
    nocase,
    rtrim
};

enum class ordering
{
    asc,
    desc
};

enum class join_type
{ 
    // sqlite supports only these 2
    left,
    full
};

enum class wildcard 
{ 
    all 
};

enum class compose_type
{
    union_type,
    union_all,
    intersect,
    except
};

}}
