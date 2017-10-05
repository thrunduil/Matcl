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

#include "matcl-dynamic/function.h"
#include "numeric_conversions_table.h"
#include "null_conversions_table.h"
#include "overload_resolution.h"
#include "error_handler.h"
#include "matcl-dynamic/details/register_function_impl.h"

#include <vector>
#include <map>
#include <set>

namespace matcl { namespace dynamic { namespace details
{

class converter_set
{
    private:
        using func_match    = std::pair<function,e_match_type>;
        using func_vec      = std::vector<func_match>;

    private:
        func_vec            m_overloads_vec;

    public:
        converter_set(Integer over_num);
        ~converter_set();

        void				push_back(const function& fun_evl, e_match_type match, 
                                error_handler& eh);
        Integer             size() const;
        void                clear();
        function            get_function(Integer pos, e_match_type& match) const;

    private:
        void                check_converter(const function& fun_evl, error_handler& eh);
};

class converter_candidate_set
{
    private:
        //final converter x main converter
        using func2         = std::pair<function,function>;
        using func_vec      = std::vector<func2>;

    private:
        func_vec            m_overloads_vec;
        conversion_match    m_match;

    public:
        converter_candidate_set(Integer over_num);
        ~converter_candidate_set();

        void				add(const function& fun_evl, conversion_match match, 
                                const function& main_fun_evl);
        void				set(const function& fun_evl, conversion_match match, 
                                const function& main_fun_evl);
        void                join(const converter_candidate_set& other);
        conversion_match    get_match() const;
        Integer             size() const;
        void                clear();
        //final function to call
        const function&     get_final_function(Integer pos) const;
        //main function only for error reporting purpose
        const function&     get_main_function(Integer pos) const;
};

using assigner_candidate_set = converter_candidate_set;

struct function_name_templ
{
    using type_ptr      = const Type*;

    function_name       m_name;
    int                 n_templates;
    mutable type_ptr    m_templates;
    bool                m_owns_array;

    function_name_templ(function_name name, const std::vector<Type>& targs);

    // can only be used as a temporary object; do not make copy of
    // templates
    function_name_templ(function_name name, int n_templ, const Type* templates);

    ~function_name_templ();

    function_name_templ(const function_name_templ&) = delete;
    function_name_templ& operator=(const function_name_templ&) = delete;

    function_name_templ(function_name_templ&& other);
    function_name_templ(const function_name_templ&& other);

    bool operator<(const function_name_templ& other) const;
};

class function_table
{
    private:
        using type_vec      = std::vector<Type>;
        using function_map  = std::map<function_name, overload_set>;
        using templ_map     = std::map<function_name_templ, overload_set>;
        using function_templ_map
                            = std::map<function_name, templ_map>;
        using unifier_set   = std::set<Type>;
        using num_convert   = numeric_conversions_table;
        using null_convert  = null_conversions_table;        

    private:
        function_map        m_map;
        function_templ_map  m_map_templ;
        function_templ_map  m_map_templ_deduce;
        converter_set       m_converters;
        overload_set        m_assigners;
        num_convert         m_num_convert_table;
        null_convert        m_null_convert_table;
        unifier_set         m_unifiers;

        function            m_function_id;
        function            m_id_assign;
        function            m_function_convert_to_any;
        function            m_function_convert_to_unit;

    public:
        function_table();
        ~function_table();

        //func_set is cleared at entry
        void                get_assigner(Type to, Type from, assigner_candidate_set& func_set,
                                error_handler& eh) const;

        void                get_converter(Type to, Type from, converter_type c_type, 
                                converter_candidate_set& os) const;
        //func_set is cleared at entry
        void                get_overload(const function_name& func, int n_args, const Type t[],
                                candidate_set& func_set, error_handler& eh) const;
        void                get_template_overload(const function_name& func, int n_templ, 
                                const Type templates[], int n_args, const Type args[],
                                candidate_set& func_set, error_handler& eh);
        void                get_registered_unifier(Type t1, Type t2, candidate_type_set& types) const;

        void                insert(const function_name& func, const function& fun_evl, 
                                make_return_fptr ret, error_handler& eh);
        void                insert_template(const function_name& func, const function& fun_evl, 
                                const type_vec& templates, make_return_fptr ret, error_handler& eh);
        void                finish_initialization();

        function            make_exact(const function& fun_evl, int n_deduced, const Type deduced[], 
                                Type deduced_ret, int n_args, const Type t[], error_handler& er) const;

    private:        
        void                get_unifier_type(Type t1, Type t2, Type unif, candidate_type_set& types) const;
        void                find_best_match(const overload_set& overloads, int n_args, 
                                const Type t[], candidate_set& func_set, error_handler& eh) const;
        void                insert_special(const function_name& func, const function& fun_evl, 
                                error_handler& eh);
        void                insert_assigner(const function& fun_evl, error_handler& eh);
        void                insert_conv_numeric(const function& fun_evl, e_match_type match);
        void                insert_unifier(const function& fun_evl, error_handler& eh);

        void                check_unifier(const function& fun, error_handler& eh) const;
        void                check_assigner(const function& fun, error_handler& eh) const;
        bool                contain_template(const type_vec& types) const;
        bool                match_templates(const function_name_templ& fn, int n_templ, 
                                const Type templ[], type_vec& deduced) const;
        void                make_exact(candidate_set& func_set, e_match_type match, int n_args, 
                                const Type t[], error_handler& eh) const;

    private:        
        function            get_converter_standard(Type to, Type from, e_match_type& match) const;
        function            get_converter_numeric(Type to, Type from, e_match_type& match) const;
        function            link_converters(const function& ev_1, const function& ev_2, 
                                const function& ev_3, e_match_type m1, e_match_type m2, 
                                e_match_type m3) const;
        function            link_converters(const function& ev_1, const function& ev_2, 
                                e_match_type m1, e_match_type m2) const;
        function            link_assign_convert(const function& assign, const function& conv, 
                                e_match_type conv_match) const;

        bool                is_numeric(Type  t) const;
};

};};};
