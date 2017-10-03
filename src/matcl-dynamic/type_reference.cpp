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

#include "type_reference.h"
#include "matcl-dynamic/function.h"
#include "matcl-dynamic/details/type_object.inl"
#include "matcl-dynamic/error_handler.h"

namespace matcl { namespace dynamic { namespace details
{

reference_type::reference_type(Type base)
    :type_impl("ref " + base.to_string()), m_base(base)
{};

function reference_type::generate_function(predef_fun fun) const
{
    (void)fun;
    return function();
}
reference_type::data_type* reference_type::clone(const data_type*) const
{
    error_handler eh;
    eh.data_operation_on_reference_type();
    eh.report();
    throw;
}
reference_type::data_type* reference_type::copy(const object_data_base*) const
{
    error_handler eh;
    eh.data_operation_on_reference_type();
    eh.report();
    throw;
}
bool reference_type::is_zero(const data_type*) const
{
    error_handler eh;
    eh.data_operation_on_reference_type();
    eh.report();
    throw;
}
bool reference_type::is_one(const data_type*) const
{
    error_handler eh;
    eh.data_operation_on_reference_type();
    eh.report();
    throw;
}
void reference_type::disp(const data_type*, md::printer&, Integer, 
                            align_type, Integer) const
{
    error_handler eh;
    eh.data_operation_on_reference_type();
    eh.report();
    throw;
}
reference_type::data_type* reference_type::create() const
{
    error_handler eh;
    eh.data_operation_on_reference_type();
    eh.report();
    throw;
}
bool reference_type::has_one() const
{
    error_handler eh;
    eh.data_operation_on_reference_type();
    eh.report();
    throw;
}
reference_type::data_type* reference_type::create_one() const
{
    error_handler eh;
    eh.data_operation_on_reference_type();
    eh.report();
    throw;
}

void reference_type::save(oarchive_impl&, unsigned int, const data_type*) const
{
    error_handler eh;
    eh.data_operation_on_reference_type();
    eh.report();
    throw;
}
reference_type::data_type* reference_type::load(iarchive_impl&, unsigned int) const
{
    error_handler eh;
    eh.data_operation_on_reference_type();
    eh.report();
    throw;
}

void reference_type::save(std::ostream&, const data_type*) const
{
    error_handler eh;
    eh.data_operation_on_reference_type();
    eh.report();
    throw;
}
reference_type::data_type* reference_type::load(std::istream&) const
{
    error_handler eh;
    eh.data_operation_on_reference_type();
    eh.report();
    throw;
}

Type get_arg_type_eval<object&>::eval(const object& obj)
{
    return operations::make_reference_type(obj.get_type());
};

};};};
