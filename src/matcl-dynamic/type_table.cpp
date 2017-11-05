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

#include "type_table.h"
#include "type_name.h"
#include "matcl-dynamic/function.h"
#include "matcl-dynamic/special_types.h"
#include "matcl-dynamic/details/register_function_impl.h"
#include "matcl-dynamic/predefined_functions.h"
#include "matcl-dynamic/details/type_impl.h"
#include "predefined/special_functions.h"
#include "type_reference.h"
#include "type_table_cache.inl"

#include <sstream>

namespace matcl { namespace dynamic { namespace details
{

void type_table::register_object(const char* obj, constr_type constr)
{
    std::string obj_name = get_class_name(obj);

    using value_type = type_map::value_type;
    type_data_handle h;
    h.set_constructor(constr);

    auto pos = m_type_map.insert(value_type(obj_name,h));

    if (pos.second == true)
        add_predefined_functions(obj_name);
};

struct predefined_functions_reg : delayed_function_register
{
    predef_fun  fun;
    std::string class_name;

    predefined_functions_reg(predef_fun f, const std::string& cn)
        :fun(f), class_name(cn)
    {};

    virtual function execute_registration() override
    {
        Type t = type_table::get()->get_type_data_from_class(class_name);
        return t.get_impl()->generate_function(fun);
    };

    virtual function_name get_fun_name() const override
    {
        switch (fun)
        {
            case predef_fun::is_zero:
                return functions::is_zero::eval();
            case predef_fun::is_one:
                return functions::is_one::eval();
            default:
            {
                //assertion(0,"unknown case");
                throw;
            }
        }
    }
};

void type_table::add_predefined_functions(const std::string& type_name)
{
    //is_zero
    {
        delayed_function_register* fi
            = new predefined_functions_reg(predef_fun::is_zero, type_name);

        initial_register_function(fi);
    }
    //is_one
    {
        delayed_function_register* fi
            = new predefined_functions_reg(predef_fun::is_one, type_name);

        initial_register_function(fi);
    }
};

std::string type_table::get_class_name(const char* name)
{
    bool simpl;
    std::string obj_name = class_name_parser().make(name,simpl,false);
    return obj_name;
};

//---------------------------------------------------------------------------
//                  type_data_handle
//---------------------------------------------------------------------------        
void type_data_handle::set_constructor(constr_type cons)
{
    m_constructor = cons;
};

void type_data_handle::clear_global()
{
    m_pool->clear_global();
    delete m_type.get_impl();
};

void type_data_handle::call_constructor()
{
    m_type  = m_constructor(m_pool);
}

Type type_data_handle::get_type() const
{
    return m_type;
};

Type type_data_handle::get_data(pool_type*& pool) const
{
    pool    = m_pool;
    return m_type;
};

//---------------------------------------------------------------------------
//                  type_table
//---------------------------------------------------------------------------        
type_table* type_table::g_instance  = nullptr;

void type_table::open()
{
    if (g_instance == nullptr)
        g_instance = new type_table;
}

type_table::type_table()
{};

Type type_table::get_type_data(const char* name, object_data_pool_impl*& pool)
{
    std::string cl_name = get_class_name(name);
    return get_type_data_from_class(cl_name, pool);
};

Type type_table::get_type(const char* name)
{
    std::string cl_name = get_class_name(name);
    pool_type* p;
    return get_type_data_from_class(cl_name, p);
};

Type type_table::get_type_data_from_class(const std::string& cl_name)
{
    pool_type* p;
    return get_type_data_from_class(cl_name, p);
};

Type type_table::get_type_data_from_class(const std::string& cl_name, pool_type*& pool)
{
    using iterator  = type_map::iterator;
    iterator pos    = m_type_map.find(cl_name);

    if (pos == m_type_map.end())
        return Type();

    if (pos->second.get_type() != Type())
        return pos->second.get_data(pool);

    pos->second.call_constructor();
    return pos->second.get_data(pool);
}

Type type_table::get_predefined(predefined_type t)
{
    switch(t)
    {
        case ti_unit:       return get_type(typeid(unit_type).name());
        case ti_any:        return get_type(typeid(any_type).name());
        case ti_bool:       return get_type(typeid(bool).name());
        case ti_int:        return get_type(typeid(Integer).name());        
        case ti_float:      return get_type(typeid(Float).name());
        case ti_real:       return get_type(typeid(Real).name());
        case ti_compl:      return get_type(typeid(Complex).name());
        case ti_fcompl:     return get_type(typeid(Float_complex).name());
        case ti_string:     return get_type(typeid(std::string).name());
        default:
        {
            //assertion(0,"unknown case");
            throw;
        }
    };
};

Type type_table::make_reference_type(Type t)
{
    if (Type::is_reference(t) == true)
    {
        // reference of a reference is a reference
        // just a safe guard, this case is possible only if this function 
        // is called explicitly
        return t;
    }

    std::unique_lock<matcl::default_spinlock_mutex> lock(m_mutex_build_type);

    auto pos = m_type_map_ref.find(t);

    if (pos != m_type_map_ref.end())
        return pos->second;

    Type ref_t  = Type(new reference_type(t));
    m_type_map_ref.insert(pos, type_map_ref::value_type(t, ref_t));
    return ref_t;
};

MATCL_THREAD_LOCAL static 
type_table_cache* g_inst = new type_table_cache();

inline type_table_cache* type_table::get_cache()
{
    return g_inst;
};

void type_table::free_cache()
{
    get_cache()->clear();
};

Type type_table::unify_types(Type t1, Type t2)
{
    if (t1 == t2)
        return t1;

    Type t = get_cache()->get_unifier(t1,t2);
    if (t != Type())
        return t;

    initialize();

    //unit type is not a valid unifier
    converter_candidate_set set(1);
    converter_candidate_set set2(1);
    if (t1 != predefined::type_unit())
        m_fun_tab.get_converter(t1, t2, converter_type::conv_implicit, set);
    if (t2 != predefined::type_unit())
        m_fun_tab.get_converter(t2, t1, converter_type::conv_implicit, set2);

    set.join(set2);

    candidate_type_set ts(1);
    m_fun_tab.get_registered_unifier(t1, t2, ts);

    if (ts.size() == 0 && set.size() == 0)
    {
        Type ty = predefined::type_any();
        get_cache()->set_unifier(t1,t2, ty);
        return ty;
    }

    //if match for converters is equal to match for unifiers, then prefer converter
    bool select_converter   = (ts.size() == 0) 
                            || (set.get_match().total_match() <= ts.get_match().total_match());

    if (select_converter == true)
    {
        if (set.size() == 1)
        {
            function f   = set.get_final_function(0);
            Type ty      = f.return_type();
            get_cache()->set_unifier(t1,t2, ty);
            return ty;
        }

        // set.size() > 1
        error_handler eh;
        eh.error_unify_ambiguity(t1,t2,set);
        eh.report();

        return Type();
    }
    else
    {
        if (ts.size() == 1)
        {
            Type ty = ts.get_type(0);
            get_cache()->set_unifier(t1,t2, ty);
            return ty;
        }

        //ts.size() > 1
        error_handler eh;
        eh.error_unifiers_ambiguity(t1,t2,ts);
        eh.report();

        return Type();
    };
};

function
type_table::make_overload(const function_name& func, const Type t[], int n_args)
{
    initialize();

    error_handler eh;

    candidate_set candidates(1);
    m_fun_tab.get_overload(func, n_args, t,candidates, eh);
    eh.report();

    if (candidates.size() == 1)
    {
        function buf = candidates.get_function(0).function();
        return get_cache()->set_overload(func, t, n_args, buf);
    };

    if (candidates.size() == 0)
        eh.error_function_not_found(func, n_args, t);
    else
        eh.error_function_ambiguity(func, n_args, t, candidates);

    eh.report();
    return function();
};

function type_table::get_template_overload(const function_name& func, int n_templ, 
                        const Type templates[], int n_args, const Type targs[])
{
    function f;
    
    f = get_cache()->get_template_overload(func,n_templ, templates, n_args, targs);

    if (f.is_null() == false)
        return f;

    initialize();

    error_handler eh;

    candidate_set candidates(1);
    m_fun_tab.get_template_overload(func, n_templ, templates, n_args, targs, candidates, eh);
    eh.report();

    if (candidates.size() == 1)
    {
        function buf = candidates.get_function(0).function();
        return get_cache()->set_template_overload(func, n_templ, templates, n_args, targs, buf);
    };
    if (candidates.size() == 0)
        eh.error_template_function_not_found(func, n_templ, templates, n_args, targs);
    else
        eh.error_template_function_ambiguity(func, n_templ, templates, n_args, targs, candidates);

    eh.report();
    return function();
};

function
type_table::get_converter(Type to, Type from, bool implicit)
{
    converter_type c_type = implicit ? converter_type::conv_implicit 
                                     : converter_type::conv_explicit;

    function f;
    f = get_cache()->get_converter(to, from, c_type);

    if (f.is_null() == false)
        return f;

    initialize();
    converter_candidate_set set(1);
    m_fun_tab.get_converter(to, from, c_type, set);

    if (set.size() == 1)
    {
        function buf = set.get_final_function(0);
        return get_cache()->set_converter(to,from, c_type, buf);
    };

    error_handler eh;

    if (set.size() == 0)
        eh.error_unable_to_convert(to,from,c_type);
    else
        eh.error_convert_ambiguity(to,from,c_type,set);

    eh.report();
    return function();
};

function type_table::get_cast_function(Type to, Type from)
{
    function f;
    f = get_cache()->get_converter(to,from,converter_type::conv_cast);

    if (f.is_null() == false)
        return f;

    initialize();
    converter_candidate_set set(1);
    m_fun_tab.get_converter(to, from, converter_type::conv_cast, set);

    if (set.size() == 1)
    {
        function buf = set.get_final_function(0);
        return get_cache()->set_converter(to,from, converter_type::conv_cast, buf);
    };

    error_handler eh;

    if (set.size() == 0)
        eh.error_unable_to_convert(to,from,converter_type::conv_cast);
    else
        eh.error_convert_ambiguity(to,from,converter_type::conv_cast, set);

    eh.report();
    return function();
};

function type_table::get_assigner(Type to, Type from)
{
    function f;
    f = get_cache()->get_assigner(to,from);

    if (f.is_null() == false)
        return f;

    initialize();

    error_handler eh;

    assigner_candidate_set as(1);
    m_fun_tab.get_assigner(to, from, as, eh);

    eh.report();

    if (as.size() == 1)
    {
        function buf = as.get_final_function(0);
        return get_cache()->set_assigner(to,from, buf);
    };

    if (as.size() == 0)
        eh.error_unable_to_assign(to,from);
    else
        eh.error_assign_ambiguity(to, from, as);        

    eh.report();
    return function();
}

void type_table::initialize()
{
    static bool fun_table_initialized = false;
    
    if(!fun_table_initialized)
    {
        error_handler eh;
        fun_table_initialized = process_functions(eh);
        eh.report();
    };

    fun_table_initialized = true;
};

void type_table::finish_initialization()
{
    m_fun_tab.finish_initialization();
};

struct delayed_func_vec : matcl_new_delete, global_object
{
    std::vector<delayed_function_register*> m_data;

    void clear_global() override
    {
        for(int i = (int)m_data.size() - 1; i >= 0 ; --i)
            delete m_data[i];

        m_data.clear();
    };

    void close_global() override
    {
        delete this;
    };
};

struct delayed_templ_func_vec : matcl_new_delete, global_object
{
    std::vector<delayed_function_template_register*> m_data;

    void clear_global() override
    {
        for (size_t i = 0; i < m_data.size(); ++i)
            delete m_data[i];

        m_data.clear();
    };

    void close_global() override
    {
        delete this;
    };
};

static bool                     g_is_initialized = false;
static delayed_func_vec*        g_reg_int_vec = nullptr;
static delayed_templ_func_vec*  g_reg_templ_int_vec = nullptr;

bool type_table::process_functions(error_handler& eh)
{
    std::unique_lock<matcl::default_mutex> lock(m_mutex_reg_functions);

    if (g_is_initialized == true)
        return true;

    if (g_reg_int_vec != nullptr)
    {
        for(int i = (int)g_reg_int_vec->m_data.size() - 1; i >= 0 ; --i)
        {
            function fun_evl = g_reg_int_vec->m_data[i]->execute_registration();
            register_function(g_reg_int_vec->m_data[i]->get_fun_name(), fun_evl, eh);

            delete g_reg_int_vec->m_data[i];
            g_reg_int_vec->m_data.pop_back();
        }
    };

    std::vector<Type> types;

    if (g_reg_templ_int_vec != nullptr)
    {
        for(int i = (int)g_reg_templ_int_vec->m_data.size() - 1; i >= 0 ; --i)
        {
            make_return_fptr ret;
            function fun_evl = g_reg_templ_int_vec->m_data[i]->execute_registration(types, ret);

            register_function_template(g_reg_templ_int_vec->m_data[i]->get_fun_name(), 
                                       fun_evl, types, ret, eh);

            delete g_reg_templ_int_vec->m_data[i];
            g_reg_templ_int_vec->m_data.pop_back();
        }
    };

    finish_initialization();
    g_is_initialized = true;
    return true;
};

void type_table::register_function(const function_name& fn, function fun_evl, 
                                   error_handler& eh)
{
    m_fun_tab.insert(fn, fun_evl, nullptr, eh);
};

void type_table::register_function_template(const function_name& fn, function fun_evl, 
                    const type_vec& types, make_return_fptr ret, error_handler& eh)
{
    m_fun_tab.insert_template(fn, fun_evl, types, ret, eh);
};

void type_table::initial_register_function(delayed_function_register* fun)
{
    if (!g_reg_int_vec)
        g_reg_int_vec = new delayed_func_vec();

    g_reg_int_vec->m_data.push_back(fun);
};

void type_table::initial_register_function_template(delayed_function_template_register* fun)
{
    if (!g_reg_templ_int_vec)
        g_reg_templ_int_vec = new delayed_templ_func_vec();

    g_reg_templ_int_vec->m_data.push_back(fun);
};

void type_table::clear_global()
{
    m_fun_tab.clear_global();

    for (auto& elem : m_type_map)
        elem.second.clear_global();

    m_type_map.clear();

    for (auto& elem : m_type_map_ref)
        delete elem.second.get_impl();
    
    m_type_map_ref.clear();
};

void type_table::close_global()
{
    delete this;
};

};};};