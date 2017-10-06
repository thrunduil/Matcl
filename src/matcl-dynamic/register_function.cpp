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

#include "matcl-dynamic/register_function.h"
#include "type_table.h"
#include "predefined/special_functions.h"

namespace matcl { namespace dynamic { namespace details
{

fun_register_helper::fun_register_helper(delayed_function_register* fr_int)
{	
	details::type_table::initial_register_function(fr_int);
};

fun_templ_register_helper::fun_templ_register_helper(delayed_function_template_register* fr_int)
{	
	details::type_table::initial_register_function_template(fr_int);
};

};};};

namespace matcl { namespace dynamic
{

const function_name& convert_promotion::eval()
{
    return special_functions::convert_promotion();
}
const function_name& convert_equivalent::eval()
{
    return special_functions::convert_equivalent();
};

const function_name& convert_decay::eval()
{
    return special_functions::convert_decay();
};

const function_name& convert_explicit::eval()
{
    return special_functions::convert_explicit();
};
const function_name& convert_cast::eval()
{
    return special_functions::convert_cast();
};

const function_name& details::unifier_tag::eval()
{
    return special_functions::unifier();
};

const function_name& details::assign_tag::eval()
{
    return special_functions::assign();
};

}};
