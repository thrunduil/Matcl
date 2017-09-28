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

#include "function_table.h"
#include "matcl-dynamic/type.h"
#include "matcl-dynamic/predefined_functions.h"
#include "predefined/special_functions.h"
#include "matcl-dynamic/details/type_impl.h"
#include "func_conv_link.h"
#include "matcl-core/details/stack_array.h"
#include "type_table_cache_data.h"
#include "matcl-dynamic/object_type.h"

#include <algorithm>

namespace matcl { namespace dynamic { namespace details
{

//------------------------------------------------------------
//                  converter_set
//------------------------------------------------------------
converter_set::converter_set(Integer over_num)
{
    m_overloads_vec.reserve(over_num);
};

converter_set::~converter_set()
{};

void converter_set::push_back(function fun_evl, e_match_type match, error_handler& eh)
{
    check_converter(fun_evl,eh);
    m_overloads_vec.push_back(func_match(fun_evl,match));
}

void converter_set::check_converter(function fun_evl, error_handler& eh)
{
    if (fun_evl.number_arguments() != 1)
        eh.error_invalid_converter_number_args(fun_evl);

    //redefinition of conversions between predefined types is allowed, but existing one
    //are still available
};

Integer converter_set::size() const
{
    return m_overloads_vec.size();
};

void converter_set::clear()
{
    m_overloads_vec.clear();
};

function converter_set::get_function(Integer pos, e_match_type& match) const
{
    match = m_overloads_vec[pos].second;
    return m_overloads_vec[pos].first;
};

//------------------------------------------------------------
//                  converter_candidate_set
//------------------------------------------------------------
converter_candidate_set::converter_candidate_set(Integer over_num)
{
    m_overloads_vec.reserve(over_num);
    m_match = conversion_match::make_no_match();
};

converter_candidate_set::~converter_candidate_set()
{};

conversion_match converter_candidate_set::get_match() const
{
    return m_match;
};

Integer converter_candidate_set::size() const
{
    return m_overloads_vec.size();
};

void converter_candidate_set::clear()
{
    m_overloads_vec.clear();
    m_match = conversion_match::make_no_match();
};

function converter_candidate_set::get_final_function(Integer pos) const
{
    return m_overloads_vec[pos].first;
};
function converter_candidate_set::get_main_function(Integer pos) const
{
    return m_overloads_vec[pos].second;
};

void converter_candidate_set::add(function fun_evl, conversion_match match, function main_fun_evl)
{
    if (match < get_match())
    {
        this->clear();
        this->m_match = match;
        m_overloads_vec.push_back(func2(fun_evl, main_fun_evl));
    }
    else if (match == get_match())
    {
        m_overloads_vec.push_back(func2(fun_evl, main_fun_evl));
    };
};
void converter_candidate_set::join(const converter_candidate_set& other)
{
    if (this->get_match() < other.get_match())
        return;

    if (this->get_match() > other.get_match())
    {
        this->clear();
        this->m_match = other.get_match();
        this->m_overloads_vec = other.m_overloads_vec;
        return;
    };

    //equal
    for (const auto& pos : other.m_overloads_vec)
        this->m_overloads_vec.push_back(pos);
    return;
}

void converter_candidate_set::set(function fun_evl, conversion_match match, function main_fun_evl)
{
    this->clear();
    this->m_match = match;
    m_overloads_vec.push_back(func2(fun_evl, main_fun_evl));
};

//------------------------------------------------------------
//                  function_name_templ
//------------------------------------------------------------
function_name_templ::function_name_templ(function_name name, const std::vector<Type>& targs)
    :m_name(name), m_owns_array(true)
{
    n_templates = targs.size();
    m_templates = new Type[targs.size()];

    Type* type_ptr  = const_cast<Type*>(m_templates);

    for (int i = 0; i < n_templates; ++ i)
        type_ptr[i]  = targs[i];
};

function_name_templ::function_name_templ(function_name name, int n_templ, const Type* templates)
    :m_name(name), m_owns_array(false), n_templates(n_templ), m_templates(templates)
{
}

function_name_templ::~function_name_templ()
{
    if (m_owns_array == false)
        return;

    delete[] m_templates;
};

function_name_templ::function_name_templ(function_name_templ&& other)
    :m_name(other.m_name), m_owns_array(other.m_owns_array), n_templates(other.n_templates)
    , m_templates(other.m_templates)
{
    other.m_templates = nullptr;
}

function_name_templ::function_name_templ(const function_name_templ&& other)
    :m_name(other.m_name), m_owns_array(other.m_owns_array), n_templates(other.n_templates)
    , m_templates(other.m_templates)
{
    other.m_templates = nullptr;
}

bool function_name_templ::operator<(const function_name_templ& other) const
{
    if (m_name < other.m_name)
        return true;
    if (m_name > other.m_name)
        return false;
    if (n_templates < other.n_templates)
        return true;
    if (n_templates > other.n_templates)
        return false;

    for (int i = 0; i < n_templates; ++i)
    {
        if (m_templates[i] < other.m_templates[i])
            return true;
        if (m_templates[i] > other.m_templates[i])
            return false;
    }

    return false;
}

//------------------------------------------------------------
//                  function_table
//------------------------------------------------------------
function_table::function_table()
    :m_converters(5),m_assigners(5)
{};

function_table::~function_table()
{};

void function_table::finish_initialization()
{};

void function_table::insert_conv_numeric(function fun_evl, e_match_type match)
{
    m_num_convert_table.insert(fun_evl, match);
};

void function_table::insert_unifier(function fun_evl, error_handler& eh)
{
    check_unifier(fun_evl,eh);
    Type unif_type = fun_evl.return_type();
    m_unifiers.insert(unif_type);
};

void function_table::insert_assigner(function fun_evl, error_handler& eh)
{
    check_assigner(fun_evl, eh);

    Type to     = fun_evl.argument_type(0);
    Type from   = fun_evl.argument_type(1);

    if (to == from)
        to.get_impl()->m_has_trivial_assignment = false;

    m_assigners.push_back(fun_evl, nullptr);
    return;
};

void function_table::insert_special(const function_name& func, function fun_evl, error_handler& eh)
{
    if (func == special_functions::assign())
    {
        insert_assigner(fun_evl, eh);
        return;
    };
    if (func == special_functions::assign_id())
    {
        m_id_assign = fun_evl;
        return;
    };

    if (func == special_functions::convert_id())
    {
        m_function_id = fun_evl;
        return;
    };

    if (func == special_functions::convert_unit())
    {
        m_function_convert_to_unit = fun_evl;
        return;
    };

    if (func == special_functions::convert_any())
    {
        m_function_convert_to_any = fun_evl;
        return;
    };

    if (func == special_functions::convert_promotion())
    {
        m_converters.push_back(fun_evl, e_match_type::user_promotion, eh);
        return;
    };
    if (func == special_functions::convert_equivalent())
    {
        m_converters.push_back(fun_evl, e_match_type::user_convert, eh);
        return;
    };
    if (func == special_functions::convert_decay())
    {
        m_converters.push_back(fun_evl, e_match_type::user_decay, eh);
        return;
    };
    if (func == special_functions::convert_explicit())
    {
        m_converters.push_back(fun_evl, e_match_type::user_explicit, eh);
        return;
    };
    if (func == special_functions::convert_cast())
    {
        m_converters.push_back(fun_evl, e_match_type::user_cast, eh);
        return;
    };

    if (func == special_functions::convert_numeric_promotion_1())
    {
        insert_conv_numeric(fun_evl, e_match_type::standard_promotion_1);
        return;
    };
    if (func == special_functions::convert_numeric_promotion_2())
    {
        insert_conv_numeric(fun_evl, e_match_type::standard_promotion_2);
        return;
    };
    if (func == special_functions::convert_numeric_equivalent_1())
    {
        insert_conv_numeric(fun_evl, e_match_type::standard_convert_1);
        return;
    };
    if (func == special_functions::convert_numeric_equivalent_2())
    {
        insert_conv_numeric(fun_evl, e_match_type::standard_convert_2);
        return;
    };
    if (func == special_functions::convert_numeric_decay_1())
    {
        insert_conv_numeric(fun_evl, e_match_type::standard_decay_1);
        return;
    };
    if (func == special_functions::convert_numeric_decay_2())
    {
        insert_conv_numeric(fun_evl, e_match_type::standard_decay_2);
        return;
    };
    if (func == special_functions::convert_numeric_int_float())
    {
        insert_conv_numeric(fun_evl, e_match_type::user_decay);
        return;
    };    

    if (func == special_functions::convert_numeric_explicit())
    {
        insert_conv_numeric(fun_evl, e_match_type::user_explicit);
        return;
    };
    if (func == special_functions::convert_numeric_cast())
    {
        insert_conv_numeric(fun_evl, e_match_type::user_cast);
        return;
    };

    if (func == special_functions::unifier())
    {
        insert_unifier(fun_evl, eh);
        return;
    };
};

void function_table::insert(const function_name& func, function fun_evl, make_return_fptr ret,
                            error_handler& eh)
{
    if (special_functions::is_special(func) == true)
    {
        insert_special(func, fun_evl, eh);
        return;
    };

    auto pos = m_map.find(func);

    if (pos == m_map.end())
        pos = m_map.insert(pos, function_map::value_type(func, overload_set(2)));

    if (pos->first.get_validator() != func.get_validator())
    {
        eh.error_function_defined_with_different_validator(pos->first,func);
        return;
    };

    bool valid = pos->first.validate_function(fun_evl);

    if (valid == true)
        pos->second.push_back(fun_evl, ret);
    else
        eh.error_function_constraints_not_satisfied(func, fun_evl);
};

void function_table::insert_template(const function_name& func, function fun_evl, 
                    const type_vec& types, make_return_fptr ret, error_handler& eh)
{
    using iterator_1    = function_templ_map::iterator;
    using iterator_2    = templ_map::iterator;

    if (types.size() == 0)
    {
        insert(func, fun_evl, ret, eh);
        return;
    };

    bool has_templ  = contain_template(types);

    iterator_1 pos_1;

    if (has_templ == true)
        pos_1   = m_map_templ_deduce.find(func);
    else
        pos_1   = m_map_templ.find(func);

    if (has_templ == false && pos_1 == m_map_templ.end())
        pos_1 = m_map_templ.insert(pos_1, function_templ_map::value_type(func, templ_map()));

    if (has_templ == true && pos_1 == m_map_templ_deduce.end())
        pos_1 = m_map_templ_deduce.insert(pos_1, function_templ_map::value_type(func, templ_map()));

    function_name_templ ft(func, types);            

    templ_map& tm       = pos_1->second;
    iterator_2 pos_2    = tm.find(ft);

    if (pos_2 == tm.end())
        pos_2   = tm.insert(pos_2, templ_map::value_type(std::move(ft), overload_set(2)));

    if (pos_2->first.m_name.get_validator() != func.get_validator())
    {
        eh.error_function_defined_with_different_validator(pos_2->first.m_name,func);
        return;
    };

    bool valid = pos_2->first.m_name.validate_function(fun_evl);

    if (valid == true)
        pos_2->second.push_back(fun_evl, ret);
    else
        eh.error_function_constraints_not_satisfied(func, fun_evl);
};

bool function_table::contain_template(const type_vec& types) const
{
    Type t_templ    = Template::get_static_type();

    for (const auto& pos:types)
    {
        if (pos == t_templ)
            return true;
    }

    return false;
}
void function_table::get_overload(const function_name& func, int n_args, const Type t[],
                                  candidate_set& candidates, error_handler& eh) const
{
    auto pos = m_map.find(func);
    candidates.clear();

    if (pos == m_map.end())
        return;

    find_best_match(pos->second, n_args, t, candidates, eh);
};

void function_table::get_template_overload(const function_name& func, int n_templ, 
                            const Type templates[], int n_args, const Type args[],
                            candidate_set& candidates, error_handler& eh)
{
    candidates.clear();
    overload_resolution or(this, n_args, args, candidates);

    auto pos_func   = m_map_templ.find(func);    

    if (pos_func != m_map_templ.end())
    {
        auto pos_templ  = pos_func->second.find(function_name_templ(func, n_templ, templates));

        if (pos_templ != pos_func->second.end())
            or.find_best_match(pos_templ->second);
    }

    auto pos_templ  = m_map_templ_deduce.find(func);

    if (pos_templ != m_map_templ_deduce.end())
    {
        templ_map& tm   = pos_templ->second;

        type_vec deduced;

        for (const auto& pos_1 : tm)
        {
            const function_name_templ& fn   = pos_1.first;
            const overload_set& os          = pos_1.second;

            bool match  = match_templates(fn, n_templ, templates, deduced);

            if (match == false)
                continue;

            //we have a match
            or.find_best_match_templ(os, &fn, deduced);
        };
    };

    or.find_most_specialized(this, candidates);

    make_exact(candidates, or.get_match_type(), n_args, args, eh);
};

bool function_table::match_templates(const function_name_templ& fn, int n_templ, 
                                     const Type templates[], type_vec& deduced) const
{
    if (fn.n_templates != n_templ)
        return false;

    Type t_templ = Template::get_static_type();
    deduced.clear();

    for (int i = 0; i < n_templ; ++i)
    {
        if (fn.m_templates[i] == t_templ)
            deduced.push_back(templates[i]);
        else if (fn.m_templates[i] != templates[i])
            return false;
    };

    return true;
};
void function_table::find_best_match(const overload_set& set, int n_args, const Type t[], 
                                     candidate_set& candidates, error_handler& eh) const
{
    overload_resolution or(this, n_args, t, candidates);
    or.find_best_match(set);

    or.find_most_specialized(this, candidates);
    make_exact(candidates, or.get_match_type(), n_args, t, eh);
};

void function_table::get_converter(Type to, Type from, converter_type c_type, 
                                   converter_candidate_set& candidates) const
{
    if (to == from)
    {
        candidates.add(m_function_id, conversion_match::make_exact(), m_function_id);
        return;
    };

    e_match_type match_s;
    function ev             = get_converter_standard(to, from, match_s);

    if (match_s != e_match_type::no_match)
    {
        if (match_s == e_match_type::user_cast && allow_cast(c_type) == true)
            candidates.add(ev, conversion_match::make(match_s), ev);
        else if (match_s == e_match_type::user_explicit && allow_explicit(c_type) == true)
            candidates.add(ev, conversion_match::make(match_s), ev);
        else if (match_s < e_match_type::user_explicit)
            candidates.add(ev, conversion_match::make(match_s), ev);
        
        return;
    };

    size_t size = m_converters.size();	

	for(size_t i = 0; i < size; ++i)
	{
        e_match_type l_match_1;
        e_match_type l_match_2;
        e_match_type l_match_3;		

        function ev_2 = m_converters.get_function(i, l_match_2);

        if (l_match_2 == e_match_type::user_cast && allow_cast(c_type) == false)
        {
            //casts are ignored when searching for conversion
            continue;
        }

        function ev_1 = get_converter_standard(ev_2.argument_type(0), from, l_match_1);

        if (l_match_1 >= e_match_type::user_explicit)
            continue;

        function ev_3 = get_converter_standard(to, ev_2.return_type(), l_match_3);

        if (l_match_3 >= e_match_type::user_explicit)
            continue;

        conversion_match lmatch = conversion_match::make(l_match_1, l_match_2, l_match_3);

        if (l_match_1 == e_match_type::exact && l_match_3 == e_match_type::exact)
        {            
            if (l_match_2 == e_match_type::user_explicit && allow_explicit(c_type) == false)
            {
                //this is an exact match, but this converter is explicit;
                //no valid converters can be found
                candidates.clear();
                return;
            };

            //converter from->to is registered, use this one
            candidates.set(ev_2, lmatch, ev_2);
            return;
        }

        if (lmatch <= candidates.get_match())
        {
            function evl    = link_converters(ev_1, ev_2, ev_3, l_match_1, l_match_2, l_match_3);
            candidates.add(evl, lmatch, ev_2);
        };
	}

    if (candidates.get_match().total_match() == e_match_type::user_cast 
            && allow_cast(c_type) == false)
    {
        //best converter is cast; but casts are not allowed
        candidates.clear();
    };
    if (candidates.get_match().total_match() == e_match_type::user_explicit 
            && allow_explicit(c_type) == false)
    {
        //best converter is explicit; but explicit conversions are not allowed
        candidates.clear();
    };

    return;
};

void function_table::get_unifier_type(Type t1, Type t2, Type unif, candidate_type_set& types) const
{
    conversion_match conv_1 = overload_resolution::mach_types(this, t1, unif);
    if (conv_1.total_match() == e_match_type::no_match)
        return;

    conversion_match conv_2 = overload_resolution::mach_types(this, t2, unif);
    if (conv_2.total_match() == e_match_type::no_match)
        return;

    //total conversion is worse of two
    conversion_match total  = conv_1;
    if (total < conv_2)
        total = conv_2;

    if (total <= types.get_match())
        types.add(unif, total);
};

void function_table::get_registered_unifier(Type t1, Type t2, candidate_type_set& types) const
{
    for (Type loc :m_unifiers)
        get_unifier_type(t1, t2, loc, types);

    //base classes are implicit unifiers
    //currently only Any can be a base class, but base class is not considere as 
    //valid unifier unless no other candidates found
    //get_unifier_type(t1, t2, predefined::type_any(), types);

    overload_resolution::find_most_specialized(this, types);
};

function function_table::get_converter_standard(Type to, Type from, e_match_type& match) const
{
    //exact match if
        //the same types
        //Lvalue-to-rvalue conversion       (not applicable)
        //Array-to-pointer conversion       (not applicable)
        //Function-to-pointer conversion    (not applicable)
        //Qualification conversion          (not applicable)

    //standard promotion if
        //Integral promotions
        //Floating point promotions

    //standard conversion if
        //Integral conversions
        //Floating point conversions
        //Floating-integral conversions
        //Pointer conversions
        //Pointer to member conversions     (not applicable)
        //Boolean conversions               (not allowed)
        //derived-to-base

    if (from == to)
    {
        match = e_match_type::exact;
        return m_function_id;
    };

    //conversion to nonconst reference must be exact
    //or from ref null
    if (Type::is_reference(to) == true)
    {
        if (Type::is_reference(from) == true && Type::decay(from) == Type())
        {
            match = e_match_type::nullref_conversion;
            Type to_base = Type::decay(to);
            return m_null_convert_table.get(to_base);
        }
        else
        {
            match = e_match_type::no_match;
            return function();
        };
    };

    if (Type::is_reference(from) == true)
        from    = Type::decay(from);

    //reference conversion
    if (from == to)
    {
        match = e_match_type::reference_conversion;
        return m_function_id;
    };

    if (is_numeric(from) == true && is_numeric(to) == true)
        return get_converter_numeric(to, from, match);

    //currently conversion Type -> Any is the only inheritance relation
    if (to == predefined::type_any())
    {
        match = e_match_type::any_convert;
        return m_function_convert_to_any;
    };
    if (to == Template::get_static_type())
    {
        match = e_match_type::template_convert;
        return m_function_id;
    };

    if (to == predefined::type_unit())
    {
        match = e_match_type::user_decay;
        return m_function_convert_to_unit;
    };

    match = e_match_type::no_match;
    return function();
}

function function_table::get_converter_numeric(Type to, Type from, e_match_type& match) const
{
    numeric_conversion_ptr ptr = m_num_convert_table.get(to,from);

    function f  = ptr->get_function();
    match       = ptr->get_match();
    return f;
};

bool function_table::is_numeric(Type  t) const
{
    if      (t == predefined::type_int())           return true;
    else if (t == predefined::type_float())         return true;
    else if (t == predefined::type_real())          return true;
    else if (t == predefined::type_complex())       return true;
    else if (t == predefined::type_float_complex()) return true;

    return false;
};

function function_table::link_converters(function ev_1, function ev_2, function ev_3, 
                                    e_match_type m1, e_match_type m2, e_match_type m3) const
{
    if (m1 < e_match_type::standard_conversions)
        return link_converters(ev_2,ev_3, m2, m3);
    if (m2 < e_match_type::standard_conversions)
        return link_converters(ev_1,ev_3, m1, m3);
    if (m3 < e_match_type::standard_conversions)
        return link_converters(ev_1,ev_2, m1, m2);

    std::vector<function> conv_vec(3);
    conv_vec[0] = ev_1;
    conv_vec[1] = ev_2;
    conv_vec[2] = ev_3;

    return new details::fun_conv_link(conv_vec);
};

function function_table::link_converters(function ev_1, function ev_2,
                                    e_match_type m1, e_match_type m2) const
{
    if (m1 < e_match_type::standard_conversions)
        return ev_2;
    if (m2 < e_match_type::standard_conversions)
        return ev_1;

    std::vector<function> conv_vec(2);
    conv_vec[0] = ev_1;
    conv_vec[1] = ev_2;

    return new details::fun_conv_link(conv_vec);
};

void function_table::get_assigner(Type to, Type from, assigner_candidate_set& as,
                                  error_handler& eh) const
{
    bool has_self_assign = (to != Type() && to.get_impl()->has_trivial_assign() == false);
    as.clear();

    if (from == to && has_self_assign == false)
    {
        function assigner = m_id_assign;
        as.add(assigner,conversion_match::make_exact(), assigner);
        return;
    };

    size_t size = m_assigners.size();

	for(size_t i = 0; i < size; ++i)
	{
        //assigners cannot have make_return function
		function tmp_evl = m_assigners.get_function(i).first;

        Type as_to      = tmp_evl.argument_type(0);
        Type as_from    = tmp_evl.argument_type(1);
		
        if (as_to != to)
            continue;

        conversion_match m  = overload_resolution::mach_types(this, from, as_from);

        if (m.total_match() != e_match_type::no_match)
            as.add(tmp_evl, m, tmp_evl);
	};

    converter_type conv_type = converter_type::conv_implicit;

    if (has_self_assign == false)
    {
        converter_candidate_set cs(1);        

        get_converter(to,from,conv_type,cs);

        if (cs.get_match() < as.get_match())
        {
            if (cs.size() == 0)
            {
                eh.error_unable_to_convert(to,from,conv_type);
                return;
            };
            
            if (cs.size() > 1)
            {
                eh.error_convert_ambiguity(to,from,conv_type,cs);
                return;
            };

            function conv   = cs.get_final_function(0);
            function conv_m = cs.get_main_function(0);
            function assign = link_assign_convert(m_id_assign,conv,cs.get_match().total_match());
            as.add(assign, cs.get_match(), conv_m);
            return;
        };
    };

    if (as.size() == 0)
        return;

    if (as.size() > 1)
        return;

    //exactly one match, now make correct function
    function assign         = as.get_final_function(0);
    Type as_from            = assign.argument_type(0);

    converter_candidate_set cs(1);
    get_converter(as_from,from,conv_type, cs);

    if (cs.size() == 0)
    {
        eh.error_unable_to_convert(as_from,from,conv_type);
        return;
    }
    if (cs.size() > 1)
    {
        eh.error_convert_ambiguity(as_from,from,conv_type,cs);
        return;
    };

    function conv           = cs.get_final_function(0);
    function assign_exact   = link_assign_convert(assign, conv, cs.get_match().total_match());    
    as.clear();
    as.add(assign_exact, cs.get_match(), assign);

    return;
};

function function_table::make_exact(function fun_evl, int n_deduced, const Type deduced[], 
                                    Type deduced_ret, int n_args, const Type t[], error_handler& eh) const
{
	std::vector<function> conv_vec(n_args);

    converter_type conv = converter_type::conv_implicit;

    int ded_ret = 0;
    if (deduced_ret != Type())
        ded_ret = 1;

	for(int i = 0; i < n_args; ++i)
	{
		Type formal = fun_evl.argument_type(i + n_deduced + ded_ret);

		if(formal != t[i])
		{
            converter_candidate_set set(1);
            get_converter(formal, t[i], conv, set); 

            function converter;
            if (set.size() == 1)
                converter = set.get_final_function(0);

            if (set.size() == 0)
                eh.error_unable_to_convert(formal, t[i], conv);

            if (set.size() > 1)
                eh.error_convert_ambiguity(formal, t[i], conv, set);
            
			conv_vec[i] = converter;
		}
		else
		{
			conv_vec[i] = function();
		}
	}

	function fun = fun_evl.make_converter(n_deduced, deduced, deduced_ret, conv_vec);
	return fun;
};

void function_table::make_exact(candidate_set& candidates, e_match_type match, int n_args, 
                                const Type args[], error_handler& eh) const
{
    if (candidates.size() != 1)
        return;

    const func_templ& ft    = candidates.get_function(0);
    int n_deduced           = ft.deduced().size();
    bool deduced_ret        = ft.has_deduced_return();
    function f              = ft.function();
    function f_exact;

    if (match > e_match_type::standard_conversions || n_deduced > 0 || deduced_ret == true)
    {
        Type ret    = ft.get_deduced_return();
        f_exact     = make_exact(f, n_deduced, ft.deduced().data(), ret, n_args, args, eh);

        candidates.clear();
        candidates.push_back(func_templ(f_exact, f_exact.return_type(), ft.get_templates(), 
                                        ft.deduced()));
    };
}

function function_table::link_assign_convert(function assign, function conv, e_match_type conv_match) const
{
    if (conv_match < e_match_type::standard_conversions)
        return assign;

	std::vector<function> conv_vec(2);
    conv_vec[0] = function();
    conv_vec[1] = conv;

	function fun = assign.make_converter(0, nullptr, Type(), conv_vec);
	return fun;
};

void function_table::check_unifier(function fun, error_handler& eh) const
{
    if (fun.number_arguments() != 0)
        eh.error_invalid_unifier_number_args(fun);
};

void function_table::check_assigner(function fun, error_handler& eh) const
{
    if (fun.number_arguments() != 2)
        eh.error_invalid_assigner_number_args(fun);

    Type ret = fun.return_type();
    if (ret != predefined::type_unit())
        eh.error_invalid_assigner_return_type(fun);
};

};};};