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

#include "sql/matcl-sqlite-cpp/builder/details/conflict_impl.h"

namespace matcl { namespace sql { namespace details
{

std::string cb_to_str(conflict_behaviour cb)
{
    if(conflict_behaviour::default_cb == cb )
        return "";

    std::string sql = " ON CONFLICT ";
    
    switch( cb )
    {
        case conflict_behaviour::rollback:
            sql += "ROLLBACK";
            break;
        case conflict_behaviour::abort:
            sql += "ABORT";
            break;
        case conflict_behaviour::fail:
            sql += "FAIL";
            break;
        case conflict_behaviour::replace:
            sql += "REPLACE";
            break;
        case conflict_behaviour::ignore:
            sql += "IGNORE";
            break;
    }

    return sql;
}

}}}

