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

#include "overload_resolution.h"
#include "matcl-dynamic/function.h"
#include "matcl-core/details/stack_array.h"
#include "matcl-core/error/exception_classes.h"
#include "function_table.h"
#include "type_table_cache_data.h"
#include "matcl-dynamic/object_type.h"

#include <algorithm>
#include <cassert>

namespace matcl { namespace dynamic { namespace details
{

conversion_match conversion_match::make_exact()
{
    return conversion_match{e_match_type::exact, e_match_type::exact,
                            e_match_type::exact};
};

conversion_match conversion_match::make_no_match()
{
    return conversion_match{e_match_type::no_match, e_match_type::no_match,
                            e_match_type::no_match};
};

conversion_match conversion_match::make(e_match_type match)
{
    return conversion_match{match,e_match_type::exact,e_match_type::exact};
};

conversion_match conversion_match::make(e_match_type m1, e_match_type m2, e_match_type m3)
{
    e_match_type lev_1;
    e_match_type lev_2;
    e_match_type lev_3;    

    if (m2 < m1)
    {
        if (m3 < m2)
        {
            lev_1   = m1;
            lev_2   = m2;
            lev_3   = m3;
        }
        else if (m3 < m1)
        {
            lev_1   = m1;
            lev_2   = m3;
            lev_3   = m2;
        }
        else
        {
            lev_1   = m3;
            lev_2   = m1;
            lev_3   = m2;
        };
    }
    else
    {
        if (m3 < m1)
        {
            lev_1   = m2;
            lev_2   = m1;
            lev_3   = m3;
        }
        else if (m3 < m2)
        {
            lev_1   = m2;
            lev_2   = m3;
            lev_3   = m1;
        }
        else
        {
            lev_1   = m3;
            lev_2   = m2;
            lev_3   = m1;
        };
    };

    return conversion_match{lev_1, lev_2, lev_3};
};

bool conversion_match::operator==(const conversion_match& other) const
{
    return (m_lev_1 == other.m_lev_1) && (m_lev_2 == other.m_lev_2) 
            && (m_lev_3 == other.m_lev_3);
};

bool conversion_match::operator!=(const conversion_match& other) const
{
    return (m_lev_1 != other.m_lev_1) || (m_lev_2 != other.m_lev_2) 
            || (m_lev_3 != other.m_lev_3);
};

bool conversion_match::operator<(const conversion_match& other) const
{
    if (m_lev_1 < other.m_lev_1)
        return true;
    if (m_lev_1 > other.m_lev_1)
        return false;
    if (m_lev_2 < other.m_lev_2)
        return true;
    if (m_lev_2 > other.m_lev_2)
        return false;

    return m_lev_3 < other.m_lev_3;
};
bool conversion_match::operator>(const conversion_match& other) const
{
    return other < *this;
};
bool conversion_match::operator<=(const conversion_match& other) const
{
    return (m_lev_1 <= other.m_lev_1) && (m_lev_2 <= other.m_lev_2)
            && (m_lev_3 <= other.m_lev_3);
};
bool conversion_match::operator>=(const conversion_match& other) const
{
    return other <= *this;
};
conversion_sequence::conversion_sequence(int n_args)
    :m_seq(n_args)
{};

void conversion_sequence::clear()
{
    m_seq.clear();
}

void conversion_sequence::push(conversion_match m)
{
    m_seq.push_back(m);
};

void conversion_sequence::swap(conversion_sequence& other)
{
    m_seq.swap(other.m_seq);
};

bool conversion_sequence::operator==(const conversion_sequence& other) const
{    
    size_t size = m_seq.size();

    matcl_assert(size == other.m_seq.size(), "invalid sequences");

    for (size_t i = 0; i < size; ++i)
    {
        if (m_seq[i] != other.m_seq[i])
            return false;
    };

    return true;
};

bool conversion_sequence::operator!=(const conversion_sequence& other) const
{
    return !operator==(other);
}

conversion_sequence::compare_result
conversion_sequence::compare(const conversion_sequence& other) const
{    
    size_t size = m_seq.size();
    matcl_assert(size == other.m_seq.size(), "invalid sequences");

    //seq_1 < seq_2 if worst conversion in seq_1 is better than worst conversion
    //in seq_2 at any position excluding positions at which conversions in
    //both sequences are the same
    conversion_match    max_dif_1 = conversion_match::make_exact();
    conversion_match    max_dif_2 = conversion_match::make_exact();

    for (size_t i = 0; i < size; ++i)
    {
        if (m_seq[i] == other.m_seq[i])
            continue;

        if (max_dif_1 < m_seq[i])
            max_dif_1   = m_seq[i];

        if (max_dif_2 < other.m_seq[i])
            max_dif_2   = other.m_seq[i];
    };

    if (max_dif_1 < max_dif_2)
        return compare_result::less;
    else if (max_dif_2 < max_dif_1)
        return compare_result::greater;
    else
        return compare_result::equivalent;
};

//------------------------------------------------------------
//                  overload_set
//------------------------------------------------------------
overload_set::overload_set(Integer over_num)
{
	m_overloads_vec.reserve(over_num);
};

overload_set::~overload_set()
{};

void overload_set::push_back(const function& fun_evl, make_return_fptr ret)
{
	m_overloads_vec.push_back(fun_ret(fun_evl, ret));
}

Integer overload_set::size() const
{
    return (Integer)m_overloads_vec.size();
};

void overload_set::clear()
{
    m_overloads_vec.clear();
};

overload_set::fun_ret overload_set::get_function(Integer pos) const
{
    return m_overloads_vec[pos];	
};

//------------------------------------------------------------
//                  func_templ
//------------------------------------------------------------
int func_templ::number_templates() const
{
    return m_templates->n_templates;
}
Type func_templ::get_template(int i) const
{
    return m_templates->m_templates[i];
}

//------------------------------------------------------------
//                  candidate_set
//------------------------------------------------------------
candidate_set::candidate_set(Integer over_num)
{
	m_overloads_vec.reserve(over_num);
};

candidate_set::~candidate_set()
{};

void candidate_set::push_back(const function& fun_evl, Type ret, func_templ_ptr templates, 
                              const templ_vec& deduced)
{
	m_overloads_vec.push_back(func_templ(fun_evl, ret, templates, deduced));
}
void candidate_set::push_back(const func_templ& ft)
{
    m_overloads_vec.push_back(ft);
};

Integer candidate_set::size() const
{
    return (Integer)m_overloads_vec.size();
};

void candidate_set::clear()
{
    m_overloads_vec.clear();
};

const func_templ& candidate_set::get_function(Integer pos) const
{
    return m_overloads_vec[pos];	
};

//------------------------------------------------------------
//                  candidate_type_set
//------------------------------------------------------------
candidate_type_set::candidate_type_set(Integer over_num)
{
    m_overloads_vec.reserve(over_num);
    m_match = conversion_match::make_no_match();
};

candidate_type_set::~candidate_type_set()
{};

conversion_match candidate_type_set::get_match() const
{
    return m_match;
};

Integer candidate_type_set::size() const
{
    return (Integer)m_overloads_vec.size();
};

void candidate_type_set::clear()
{
    m_overloads_vec.clear();
    m_match = conversion_match::make_no_match();
};

Type candidate_type_set::get_type(Integer pos) const
{
    return m_overloads_vec[pos];
};

void candidate_type_set::add(Type t, conversion_match match)
{
    if (match < get_match())
    {
        this->clear();
        this->m_match = match;
        m_overloads_vec.push_back(t);
    }
    else if (match == get_match())
    {
        m_overloads_vec.push_back(t);
    };
};

//------------------------------------------------------------
//                  overload_resolusion
//------------------------------------------------------------

overload_resolution::overload_resolution(const function_table* ft, int args,
                        const Type t[], candidate_set& candidates)
    :m_match_type(e_match_type::no_match), m_ft(ft), n_args(args)
    , m_types(t), m_best(args), m_candidates(candidates)
{};

e_match_type overload_resolution::get_match_type() const
{
    return m_match_type;
}

void overload_resolution::find_best_match(const overload_set& set)
{
    find_best_match_templ(set, nullptr, type_vec());
};

void overload_resolution::find_best_match_templ(const overload_set& set, 
                 const function_name_templ* ft, const type_vec& deduced)
{
	int ov_set_size             = (int)set.size();

    conversion_sequence conv_seq(n_args);    

    using conv_cmp_res          = conversion_sequence::compare_result;
    using fun_ret               = overload_set::fun_ret;

	for(int i = 0; i < ov_set_size; ++i)
	{
		fun_ret tmp_evaler      = set.get_function(i);
        bool ded_ret            = check_deduce_return(tmp_evaler);
		e_match_type tmp_match  = mach_types_list(m_ft, (int)deduced.size(), ded_ret, 
                                    n_args, m_types, tmp_evaler.first, conv_seq);

		if (tmp_match < m_match_type)
		{
            bool error;
            Type ret        = eval_return(tmp_evaler.first, tmp_evaler.second,
                                ft, deduced, error);

            if (error == true)
                continue;

			m_match_type    = tmp_match;
            m_best.swap(conv_seq);

            m_candidates.clear();
            m_candidates.push_back(tmp_evaler.first, ret, ft, deduced);
		}
        else if (tmp_match != e_match_type::no_match && tmp_match == m_match_type)
        {
            conv_cmp_res cmp    = conv_seq.compare(m_best);
            
            if (cmp == conv_cmp_res::less)
            {
                bool error;
                Type ret        = eval_return(tmp_evaler.first, tmp_evaler.second, ft, 
                                    deduced, error);

                if (error == true)
                    continue;

			    m_match_type    = tmp_match;
                m_best.swap(conv_seq);

                m_candidates.clear();
                m_candidates.push_back(tmp_evaler.first, ret, ft, deduced);
            }
            else if (cmp == conv_cmp_res::equivalent)
            {
                bool error;
                Type ret        = eval_return(tmp_evaler.first, tmp_evaler.second, 
                                    ft, deduced, error);

                if (error == true)
                    continue;

                m_candidates.push_back(tmp_evaler.first, ret, ft, deduced);
            };
        };
	};
}

bool overload_resolution::check_deduce_return(const overload_set::fun_ret& evaler) const
{
    if (evaler.second == nullptr)
        return false;

    if (evaler.first.return_type() != Template::get_static_type())
        return false;

    return true;
};
Type overload_resolution::eval_return(const function& f, make_return_fptr ret, 
            const function_name_templ* ft, const type_vec& deduced, bool& error) const
{
    error = false;
    if (ret == nullptr)
        return Type();    
    
    int n_template  = ft ? ft->n_templates : 0;
    int n_deduced   = (int)deduced.size();
    bool deduce_ret = f.return_type() == Template::get_static_type();

    std::vector<Type>   all_types(n_template + n_args);

    //make template arguments
    for (int i = 0, j = 0; i < n_template; ++i)
    {
        Type t          = ft->m_templates[i];

        if (t == Template::get_static_type())
            t           = deduced[j++];

        all_types[i]    = t;
    };

    //make fuction arguments
    for (int i = 0; i < n_args; ++i)
    {
        Type t          = f.argument_type(i + n_deduced + 1);

        if (t == Template::get_static_type())
            t           = m_types[i];

        all_types[i + n_template]   = t;
    };

    try
    {
        Type t  = (*ret)(n_template, all_types.data(), n_args, 
                         all_types.data() + n_template);

        if (t == Type())
            error = true;

        if (deduce_ret == true)
            return t;
        else
            return Type();
    }
    catch(...)
    {
        error = true;
        return Type();
    }
};
void overload_resolution::find_most_specialized(const function_table* ft, 
                                                candidate_set& func_set)
{
    if (func_set.size() < 2)
        return;

    Integer size            = func_set.size();
    Integer best_pos        = 0;
    const func_templ* best  = &func_set.get_function(0);
    bool is_best            = true;

    for (Integer pos = 1; pos < size; ++pos)
    {
        const func_templ* loc   = &func_set.get_function(pos);

        spec_type spec          = is_more_specialized(ft, *loc, *best);

        if (spec == spec_type::first)
        {
            best_pos    = pos;
            best        = loc;
            is_best     = true;
        }
        else if (spec == spec_type::second)
        {
            //do nothing
        }
        else
        {
            //not best since is equivalent to one of other matches
            is_best     = false;
        };
    };

    if (is_best == false)
        return;

    //relation is_more_specialized is not transitive, but this is not a problem
    //most specialized match must be better than any other

    //best won at least one match and all next matches are worse (since is_best == true)
    //all others are worse or equivalent to at least one match
    //=> if there exists the most specialized match, then this match is the found (best)

    //now check if (best) is really the most specialized match
    
    for (Integer pos = 0; pos < best_pos; ++pos)
    {
        const func_templ* loc   = &func_set.get_function(pos);
        spec_type spec          = is_more_specialized(ft, *loc, *best);

        if (spec != spec_type::second)
        {
            //this match is at least equivalent to best
            is_best     = false;
            break;
        };
    };

    if (is_best == true)
    {
        func_templ fbest    = *best;

        func_set.clear();
        func_set.push_back(fbest);
    };
};

void overload_resolution::find_most_specialized(const function_table* ft, 
                                        candidate_type_set& type_set)
{
    if (type_set.size() < 2)
        return;

    Integer size        = type_set.size();
    Integer best_pos    = 0;
    Type best           = type_set.get_type(0);
    bool is_best        = true;

    for (Integer pos = 1; pos < size; ++pos)
    {
        Type loc        = type_set.get_type(pos);

        spec_type spec  = is_more_specialized(ft, loc,best);
        if (spec == spec_type::first)
        {
            best_pos    = pos;
            best        = loc;
            is_best     = true;
        }
        else if (spec == spec_type::second)
        {
            //do nothing
        }
        else
        {
            //not best since is equivalent to one of other matches
            is_best     = false;
        };
    };

    if (is_best == false)
        return;

    //relation is_more_specialized is not transitive, but this is not a problem
    //most specialized match must be better than any other

    //best won at least one match and all next matches are worse (since is_best == true)
    //all others are worse or equivalent to at least one match
    //=> if there exists the most specialized match, then this match is the found (best)

    //now check if (best) is really the most specialized match
    
    for (Integer pos = 0; pos < best_pos; ++pos)
    {
        Type loc        = type_set.get_type(pos);
        spec_type spec  = is_more_specialized(ft, loc,best);

        if (spec != spec_type::second)
        {
            //this match is at least equivalent to best
            is_best     = false;
            break;
        };
    };

    if (is_best == true)
    {
        conversion_match m = type_set.get_match();
        type_set.clear();
        type_set.add(best,m);
    };
};

spec_type overload_resolution::is_more_specialized(const function_table* ft,
                        const func_templ& func_1, const func_templ& func_2)
{
    function f1             = func_1.function();
    function f2             = func_2.function();
    Integer n_args_1        = f1.number_arguments();
    Integer n_args_2        = f2.number_arguments();
    bool ded_ret_1          = func_1.has_deduced_return();
    bool ded_ret_2          = func_2.has_deduced_return();
    Integer n_ded_1         = (Integer)func_1.deduced().size();
    Integer n_ded_2         = (Integer)func_2.deduced().size();

    int n_ded_total_1       = n_ded_1;
    int n_ded_total_2       = n_ded_2;
    if (ded_ret_1 == true)
        n_ded_total_1       += 1;
    if (ded_ret_2 == true)
        n_ded_total_2       += 1;

    Integer n_1             = n_args_1 - n_ded_total_1;
    Integer n_2             = n_args_2 - n_ded_total_2;

    matcl_assert(n_1 == n_2, "error in overload_resolution");
    (void)n_2;

    conversion_sequence conv_seq_1(n_1);
    conversion_sequence conv_seq_2(n_1);

    using pod_type          = matcl::details::pod_type<Type>;
    using stack_array       = matcl::details::stack_array<pod_type,10>;
    stack_array types(n_1);
    Type* t_ptr             = reinterpret_cast<Type*>(types.get());

    e_match_type match_1, match_2;

    for (Integer i = 0; i < n_1; ++i)
        t_ptr[i]            = f1.argument_type(i + n_ded_total_1);

    //func_2 called with arguments type of the function func_1
    match_1                 = mach_types_list(ft, n_ded_2, ded_ret_2, 
                                n_1, t_ptr, f2, conv_seq_1);

    for (Integer i = 0; i < n_1; ++i)
        t_ptr[i]            = f2.argument_type(i + n_ded_total_2);

    //func_1 called with arguments type of the function func_2
    match_2                 = mach_types_list(ft, n_ded_1, ded_ret_1, 
                                n_1, t_ptr, f1, conv_seq_2);

    if (match_1 < match_2)
    {
        //func_1 is more specialized
        return spec_type::first;
    }
    else if (match_2 < match_1)
    {
        //func_2 is more specialized
        return spec_type::second;
    }

    if (match_1 == e_match_type::no_match)
    {
        //then also match_2 == no_match
        //functions cannot be compared
        return spec_type::not_comparable;
    };

    using cmp_res   = conversion_sequence::compare_result;
    cmp_res cmp     = conv_seq_1.compare(conv_seq_2);

    if (cmp == cmp_res::less)
    {
        //func_1 is more specialized
        return spec_type::first;
    }
    if (cmp == cmp_res::greater)
    {
        //func_2 is more specialized
        return spec_type::second;
    }

    //compare number of template argument to be deduced after
    //comparing function arguments
    if (n_ded_1 < n_ded_2)
        return spec_type::first;
    else if (n_ded_2 < n_ded_1)
        return spec_type::second;

    if (n_ded_total_1 < n_ded_total_2)
        return spec_type::first;
    else if (n_ded_total_2 < n_ded_total_1)
        return spec_type::second;

    //functions are equivalent
    return spec_type::equivalent;
};

spec_type overload_resolution::is_more_specialized(const function_table* ft,
                                                  Type type_1, Type type_2)
{
    conversion_match match_1    = mach_types(ft, type_1, type_2);
    conversion_match match_2    = mach_types(ft, type_2, type_1);

    if (match_1 < match_2)
    {
        //type_1 is more specialized
        return spec_type::first;
    }
    else if (match_2 < match_1)
    {
        //type_2 is more specialized
        return spec_type::second;
    }

    if (match_1.total_match() == e_match_type::no_match)
    {
        //then also match_2 == no_match
        //types cannot be compared
        return spec_type::not_comparable;
    };

    //types are equivalent
    return spec_type::equivalent;
};

e_match_type overload_resolution::mach_types_list(const function_table* ft,
                    int n_ded, bool ded_ret, int n_args, const Type t[], 
                    const function& f, conversion_sequence& match)
{    
    if (ded_ret)
        n_ded += 1;

    if (f.number_arguments() != n_args + n_ded)
        return e_match_type::no_match;

    match.clear();	

    //first n_ded arguments must be any (by assumption how deduced arguments
    //will be passed to the function
    for (int i = 0; i < n_ded; ++i)
    {
        if (f.argument_type(i) != OType::get_static_type())
            return e_match_type::no_match;
    };

    e_match_type match_total = n_ded == 0 ? e_match_type::exact 
                                          : e_match_type::any_convert;

	for(int i = 0; i < n_args; ++i)
	{
		conversion_match tmp_match = mach_types(ft, t[i], f.argument_type(i + n_ded));

        if (tmp_match.total_match() == e_match_type::no_match)
            return e_match_type::no_match;

        match.push(tmp_match);

		if (tmp_match.total_match() > match_total)
			match_total = tmp_match.total_match();
	}

	return match_total;
};

conversion_match overload_resolution::mach_types(const function_table* ft, 
                                                Type actual, Type formal)
{
	if(actual == formal)
		return conversion_match::make_exact();

    converter_candidate_set set(1);
    ft->get_converter(formal,actual,converter_type::conv_implicit,set);
    return set.get_match();
}

};};};
