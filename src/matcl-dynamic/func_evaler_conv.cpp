/*
 *  This file is a part of Matrix Computation Library (MATCL)
 *
 *  Copyright (c) Pawe� Kowal 2017 - 2021
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

#include "matcl-dynamic/details/func_evaler_conv.h"
#include "matcl-core/details/stack_array.h"
#include "matcl-dynamic/object_type.h"
#include "matcl-dynamic/details/object.inl"

namespace matcl { namespace dynamic { namespace details
{

fun_evaler_conv* fun_evaler_conv::make(const evaler* fun, int n_deduced, const Type deduced[],
                        Type deduced_ret, const std::vector<function>& convs)
{
    return new fun_evaler_conv(fun, n_deduced, deduced, deduced_ret, convs);
};

fun_evaler_conv::fun_evaler_conv(const evaler* fun, int n_deduced, const Type deduced[],
                                 Type ded_ret, const std::vector<function>& convs)
	: m_fun(fun), m_converters(convs), m_has_deduced_ret(false)
{
	this->m_ret_ti      = fun->m_ret_ti;
	this->m_args_size   = fun->m_args_size - n_deduced;
    this->m_arg_ti      = nullptr;
    this->m_deduced_ret = ded_ret;
    int has_ret         = 0;
    
    if (this->m_deduced_ret != Type())
    {
        this->m_ret_ti      = this->m_deduced_ret;
        has_ret             = 1;
        m_has_deduced_ret   = true;
        this->m_args_size   -= 1;
    }

    if (this->m_args_size != (int)convs.size())
    {
        throw error::matcl_dynamic_exception(
                "unable to create converter, number of converters and"
                " function arguments does not match");
    };

    if (this->m_args_size > 0)
        this->m_arg_ti  = new Type[this->m_args_size];
	
    for(size_t i = 0; i < m_converters.size(); ++i)
	{
		if(m_converters[i].is_null() == false)
			m_arg_ti[i] = m_converters[i].argument_type(0);
		else
			m_arg_ti[i] = fun->m_arg_ti[i];
	}

    m_deduced.resize(n_deduced + has_ret);

    if (has_ret)
        m_deduced[0].reset(object(OType(ded_ret)));

    for (int i = 0; i < n_deduced; ++i)
        m_deduced[i + has_ret].reset(object(OType(deduced[i])));
};

void fun_evaler_conv::destroy()
{
    delete[] this->m_arg_ti;
    delete this;
};

int fun_evaler_conv::number_deduced_types() const
{
    int n = (int)m_deduced.size();

    if (m_has_deduced_ret == true)
        return n - 1;
    else
        return n;
};

Type fun_evaler_conv::get_deduced_type(size_t n) const
{
    if (m_has_deduced_ret == false)
        return m_deduced[n].get_type();
    else
        return m_deduced[1 + n].get_type();
};

Type fun_evaler_conv::get_deduced_ret() const
{
    return m_deduced_ret;
};

int fun_evaler_conv::number_converters() const
{
    return (int)m_converters.size();
}

function fun_evaler_conv::get_coverter(size_t n) const
{
    return m_converters[n];
};

bool fun_evaler_conv::make_eval(const object** _args, object& return_obj) const
{
    // avoid using new or malloc
    using obj_pod           = matcl::details::pod_type<object>;
    using obj_destructor    = obj_pod::destructor_type;
    using stack_array_ptr   = matcl::details::stack_array<const object*,10>;
    using stack_array_obj   = matcl::details::stack_array<obj_pod,10>;

    Integer n_arg           = (Integer)m_converters.size();
    Integer n_deduced       = (Integer)m_deduced.size();

    Integer size_counter    = 0;

    obj_destructor d(&size_counter);
	stack_array_ptr args(n_arg + n_deduced);
    stack_array_obj args_vec(n_arg + n_deduced, &d);

    const object** buff_ptr = args.get();
    object* buff_obj        = reinterpret_cast<object*>(args_vec.get());

    (void)buff_obj;
    make_convert(_args, buff_ptr, buff_obj, n_deduced, size_counter);

	m_fun->make_eval(buff_ptr, return_obj);

    if (m_deduced_ret != Type())
    {
        //is is ensured, that return_obj has Template type
        return_obj.reset(object(m_deduced_ret, return_obj));
    };

    return true;
};

void fun_evaler_conv::make_eval(const object** _args) const
{
    //avoid using new or malloc
    using obj_pod           = matcl::details::pod_type<object>;
    using obj_destructor    = obj_pod::destructor_type;
    using stack_array_ptr   = matcl::details::stack_array<const object*,10>;
    using stack_array_obj   = matcl::details::stack_array<obj_pod,10>;

    Integer n_arg           = (Integer)m_converters.size();
    Integer n_deduced       = (Integer)m_deduced.size();

    Integer size_counter    = 0;

    obj_destructor d(&size_counter);
	stack_array_ptr args(n_arg + n_deduced);
    stack_array_obj args_vec(n_arg + n_deduced, &d);

    const object** buff_ptr = args.get();
    object* buff_obj        = reinterpret_cast<object*>(args_vec.get());

    (void)buff_obj;
    make_convert(_args, buff_ptr, buff_obj, n_deduced, size_counter);

	m_fun->make_eval(buff_ptr);

    return;
};

void fun_evaler_conv::make_convert(const object** _args, const object** buff_ptr, 
                    object* buff_obj, Integer n_deduced, Integer& size_counter) const
{
    for (int i = 0; i < n_deduced; ++i)
        buff_ptr[i]         = &m_deduced[i];

	for(size_t i = 0; i < m_converters.size(); ++i)
	{
		if (m_converters[i].is_null() == false)
		{
			object ret;
			bool new_obj = m_converters[i].get_evaler()->make_eval(&_args[i], ret);

            if (new_obj == true)
            {
                new(buff_obj + size_counter) object(std::move(ret));
                buff_ptr[i+n_deduced]   = buff_obj + size_counter;
                ++size_counter;			
            }
            else
            {
                buff_ptr[i+n_deduced]   = _args[i];
            }
		}
		else
		{
			buff_ptr[i+n_deduced]       = _args[i];
		}
	}
};

};};};
