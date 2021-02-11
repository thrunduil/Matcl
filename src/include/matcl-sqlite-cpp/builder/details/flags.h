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

namespace matcl { namespace sql { namespace details
{

namespace create_index_flags
{

    enum create_index_flag
    {
        can_set_unique        = 1<<0,
        can_set_if_not_exists = 1<<1,
        can_set_name          = 1<<2,
        can_set_table         = 1<<3,
        can_set_cols          = 1<<4,
        can_to_str            = 1<<5,

        default_flags
            = can_set_unique
            | can_set_if_not_exists
            | can_set_name
    };

};

namespace create_table_flags
{
    enum create_table_flag
    {
        can_set_temp          = 1<<0,
        can_set_if_not_exists = 1<<1,
        can_set_table         = 1<<2,
        can_set_cols          = 1<<3,
        can_set_select        = 1<<4,
        can_add_constraints   = 1<<5,
        can_to_str            = 1<<6,

        default_flags
            = can_set_temp
            | can_set_if_not_exists
            | can_set_table
    };
}

namespace delete_flags
{

    enum delete_flag
    {
        can_set_table     = 1<<0,
        can_set_condition = 1<<1,
        can_to_str        = 1<<2,

        default_flags
            = can_set_table
            | can_set_condition
    };

}

namespace insert_flags
{

    enum insert_flag
    {
        can_set_conflict_resolver = 1<<0,
        can_select_columns        = 1<<1,
        can_select_target         = 1<<2,
        can_define_values         = 1<<3,
        can_to_str                = 1<<4,

        default_flags
            = can_set_conflict_resolver
            | can_select_columns
            | can_select_target
            | can_define_values
    };

}

namespace select_flags
{
    enum select_flag
    {
        can_set_source   = 1<<0,
        can_compose      = 1<<1,
        can_set_ordering = 1<<2,
        can_set_limit    = 1<<3,
        can_set_distinct = 1<<4,
        can_set_columns  = 1<<5,
        can_set_where    = 1<<6,
        can_set_grouping = 1<<7,

        core_flags = can_set_columns | can_set_source | can_set_distinct | can_set_where | can_set_grouping,

        default_flags
            = can_set_source
            | can_compose
            | can_set_ordering
            | can_set_limit
            | can_set_distinct
            | can_set_columns
            | can_set_where
            | can_set_grouping
    };
}

namespace update_flags
{
    enum update_flag
    {
        can_set_conflict_resolver = 1<<0,
        can_set_table             = 1<<1,
        can_set_assignments       = 1<<2,
        can_set_condition         = 1<<3,
        can_to_str                = 1<<4,

        default_flags 
            = can_set_conflict_resolver
            | can_set_table 
            | can_set_assignments
            | can_set_condition
    };
}

}}}
