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

#include "matcl-dynamic/config.h"
#include "matcl-core/config.h"

namespace matcl { namespace dynamic
{

namespace details
{
    enum class e_match_type : int;
    enum class converter_type : int;

    class type_impl;
    class object_data_base;
    class type_table;
    class function_table;
    class evaler;
    struct delayed_function_register;
    struct delayed_function_template_register;
    class identifier_impl;
    class converter_candidate_set;
    class overload_set;
    struct function_name_templ;
    class candidate_set;
    class candidate_type_set;
    class evaler;

    template<class T> class mark_type;
};

class object;
class Type;
class printer;
class function;
class null_type;

template<class T> class object_type;

}};
