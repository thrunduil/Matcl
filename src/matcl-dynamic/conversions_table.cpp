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

#pragma once

#include "conversions_table.h"
#include "matcl-core/error/exception_classes.h"
#include "matcl-core/memory/alloc.h"
#include "matcl-dynamic/details/func_evaler_conv.h"

namespace matcl { namespace dynamic { namespace details
{

//---------------------------------------------------------------
//                  null conversions
//---------------------------------------------------------------

class fun_conv_null : public evaler
{
    private:
        fun_conv_null(const fun_conv_null&) = delete;
        fun_conv_null& operator=(const fun_conv_null&) = delete;

    private:
        fun_conv_null(Type to);

        static 
        fun_conv_null*  make(Type to);

        friend null_conversion;

        void            destroy() override;

    public:        
        bool            make_eval(const object** _args, object& ret) const override;
        void            make_eval(const object** _args) const override;
};

fun_conv_null::fun_conv_null(Type to)
{
    m_args_size = 1;
    m_ret_ti    = to;
    m_arg_ti    = new Type();
    m_arg_ti[0] = Type();
}

void fun_conv_null::destroy()
{
    delete m_arg_ti;
    delete this;
}

fun_conv_null* fun_conv_null::make(Type to)
{
    return new fun_conv_null(to);
};

bool fun_conv_null::make_eval(const object** _args, object& ret) const
{
    object& obj_null = *const_cast<object*>(_args[0]);
    obj_null.reset(object(m_ret_ti));
    ret.reset(obj_null);  
    return false;
};

void fun_conv_null::make_eval(const object** _args) const
{
    object& obj_null = *const_cast<object*>(_args[0]);
    obj_null.reset(object(m_ret_ti));
    return;
};

null_conversion::null_conversion(Type to)
    :m_to(to), m_func(fun_conv_null::make(to))
{};

null_conversion::~null_conversion()
{
    const_cast<evaler*>(m_func.get_evaler())->destroy();
}

//---------------------------------------------------------------
//                  link conversions
//---------------------------------------------------------------
class fun_conv_link : public evaler
{
    public:	
        std::vector<function>       m_converters;

    private:
        fun_conv_link(const fun_conv_link&) = delete;
        fun_conv_link& operator=(const fun_conv_link&) = delete;

    private:
        fun_conv_link(const std::vector<function>& convs);

        static 
        fun_conv_link*  make(const std::vector<function>& convs);

        void            destroy() override;

        friend link_conversion;

    public:     
       bool             make_eval(const object** _args, object& ret) const override;
       void             make_eval(const object** _args) const;
};

fun_conv_link::fun_conv_link(const std::vector<function>& convs)
    :m_converters(convs)
{
    m_args_size = 1;
    m_ret_ti    = convs.back().return_type();
    m_arg_ti    = new Type();
    m_arg_ti[0] = convs.front().argument_type(0);
};

void fun_conv_link::destroy()
{
    delete m_arg_ti;
    delete this;
}

fun_conv_link* fun_conv_link::make(const std::vector<function>& convs)
{
    return new fun_conv_link(convs);
};

bool fun_conv_link::make_eval(const object** _args, object& ret) const
{
    size_t n = m_converters.size();

    const object* args = _args[0];

    const auto& ev = m_converters[0].get_evaler();
    ev->make_eval(&args, ret);    

    for (size_t i = 1; i < n; ++i)
    {
        object tmp          = object(std::move(ret));
        const auto& ev2     = m_converters[i].get_evaler();
        const object* tmp2  = &tmp;
        ev2->make_eval(&tmp2, ret);    
    };

    return true;
}

void fun_conv_link::make_eval(const object** _args) const
{
    //nothing to do
    (void)_args;
    return;
}

function link_conversion::get_function() const
{ 
    return function(m_func); 
};

link_conversion::link_conversion(const std::vector<function>& conv)
    : m_func(fun_conv_link::make(conv)), m_hash(eval_hash(conv))
{};

link_conversion::~link_conversion()
{
    m_func->destroy();
}

size_t link_conversion::eval_hash(const std::vector<function>& conv)
{
    size_t seed = conv.size();

    for (size_t i = 0; i < conv.size(); ++i)
    {
        boost::hash_combine(seed, conv[i].get_evaler());
        boost::hash_combine(seed, i);
    };

    return seed;
};

bool link_conversion::equal(const std::vector<function>& conv) const
{
    size_t size = conv.size();

    if (size != m_func->m_converters.size())
        return false;

    for (size_t i = 0; i < size; ++i)
    {
        if (conv[i].get_evaler() != m_func->m_converters[i].get_evaler())
            return false;
    }

    return true;
};

//---------------------------------------------------------------
//                  conversions_table
//---------------------------------------------------------------
conversion::conversion(const build_convert_info& info)
    : m_hash(info.m_hash)
{
    m_func = fun_evaler_conv::make(info.m_base.get_evaler(), info.n_deduced, info.m_deduced,
                                info.m_deduced_ret, *info.m_arg_converters);
};

conversion::~conversion()
{
    m_func->destroy();
}

function conversion::get_function() const
{
    return function(m_func);
}

size_t conversion::eval_hash(const build_convert_info& info)
{
    size_t seed = (size_t)info.m_base.get_evaler();

    boost::hash_combine(seed, info.n_deduced);

    for (int i = 0; i < info.n_deduced; ++i)
    {
        boost::hash_combine(seed, i);
        boost::hash_combine(seed, info.m_deduced[i].hash_value());
    };

    boost::hash_combine(seed, info.m_deduced_ret.hash_value());

    size_t n    = info.m_arg_converters->size();

    for (size_t i = 0; i < n; ++i)
    {
        boost::hash_combine(seed, i);
        boost::hash_combine(seed, (*info.m_arg_converters)[i].get_evaler());
    };

    info.m_hash = seed;
    return seed;
};

bool conversion::equal(const build_convert_info& info) const
{
    if (info.m_base.get_evaler() != m_func->m_fun)
        return false;

    if (info.n_deduced != m_func->number_deduced_types())
        return false;

    int n    = info.n_deduced;

    for (int i = 0; i < n; ++i)
    {
        if (info.m_deduced[i] != m_func->get_deduced_type(i))
            return false;
    };

    if (info.m_deduced_ret != m_func->get_deduced_ret())
        return false;

    int num_conv = (int)info.m_arg_converters->size();

    if (num_conv != m_func->number_converters())
        return false;

    for (int i = 0; i < num_conv; ++i)
    {
        if ((*info.m_arg_converters)[i].get_evaler() != m_func->get_coverter(i).get_evaler())
            return false;
    };

    return true;
};

function conversions_table::get(function base, int n_deduced, const Type deduced[],
                        Type deduced_ret, const std::vector<function>& arg_converters) const
{
    //table can be modified and this function can be called from multiple threads
    std::unique_lock<mutex> lock(m_mutex);

    build_convert_info info{base, n_deduced, deduced, deduced_ret, &arg_converters, 0};

    return m_table.get(info)->get_function(); 
};

void conversions_table::clear_global()
{
    m_table.clear();
};

};};};

