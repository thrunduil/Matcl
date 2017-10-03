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

#include "matcl-core/config.h"
#include "matcl-dynamic/details/fwd_decls.h"
//#include "matcl-core/details/fwd_decls.h"
#include "matcl-dynamic/function.h"
#include "matcl-dynamic/details/register_function_impl.h"

#include <vector>

namespace matcl { namespace dynamic { namespace details
{

enum class e_match_type : int
{
    exact = 0,

    reference_conversion,   //nonconst ref to const ref

    standard_conversions,   //mark a group; no conversion of this type

    nullref_conversion,     //nonconst ref null to nonconst ref type

    standard_promotion_1,   //floating-point promotions
    standard_promotion_2,   //floating-point promotions
    standard_convert_1,     //lossless conversions changing kind (int->float->complex)
    standard_convert_2,     //lossless conversions changing kind (int->float->complex)
    standard_decay_1,       //floating point precision loss
    standard_decay_2,       //Real->Float_complex conversion
    
    user_promotion,         //user defined promotions
    user_convert,           //user defined lossless conversions changing kind
    
    any_convert,            //conversion to Any
    template_convert,       //template match

    user_decay,             //user defined conversions with data loss
    user_explicit,          //user defined explicit conversion
    user_cast,              //user defined cast
    
    no_match
};

enum class spec_type : int
{
    first,                  //first function is more specialized
    equivalent,             //both functions are equivalent
    second,                 //second function is more specialized
    not_comparable          //functions cannot be compared    
};

struct conversion_match
{
    e_match_type            m_lev_1;
    e_match_type            m_lev_2;
    e_match_type            m_lev_3;

    e_match_type            total_match() const { return m_lev_1; };

    bool                    operator<(const conversion_match& other) const;
    bool                    operator<=(const conversion_match& other) const;
    bool                    operator>(const conversion_match& other) const;
    bool                    operator>=(const conversion_match& other) const;
    bool                    operator==(const conversion_match& other) const;
    bool                    operator!=(const conversion_match& other) const;

    static conversion_match make_exact();
    static conversion_match make_no_match();
    static conversion_match make(e_match_type match);
    static conversion_match make(e_match_type match_1, e_match_type match_u, e_match_type match_2);
};

class conversion_sequence
{
    public:
        enum class compare_result
        {
            less,           //fist match is better
            equivalent,     //matches are equivalent but may be different
            greater         //second match is better
        };

    private:
        using sequence      = std::vector<conversion_match>;

    private:
        sequence            m_seq;

    public:
        conversion_sequence(int n_args);

        void                clear();
        void                push(conversion_match m);

        void                swap(conversion_sequence& other);
        compare_result      compare(const conversion_sequence& other) const;
        bool                operator==(const conversion_sequence& other) const;
        bool                operator!=(const conversion_sequence& other) const;

    private:
        conversion_sequence(const conversion_sequence&) = delete;
        conversion_sequence& operator=(const conversion_sequence&) = delete;
};

class overload_set
{
    public:
        using fun_ret       = std::pair<function,make_return_fptr>;
        using func_vec      = std::vector<fun_ret>;

    private:
        func_vec            m_overloads_vec;

    public:
        overload_set(Integer over_num);
        ~overload_set();

        void				push_back(const function& fun_evl, make_return_fptr ret);
        Integer             size() const;
        void                clear();
        fun_ret             get_function(Integer pos) const;
};

struct func_templ
{
    using type_vec          = std::vector<Type>;
    using func_templ_ptr    = const function_name_templ*;

    function                m_function;
    func_templ_ptr          m_templates;
    type_vec                m_deduced;
    Type                    m_return;

    func_templ(const function& f, Type ret, func_templ_ptr templates, 
               const std::vector<Type>& deduced)
        :m_function(f), m_deduced(deduced), m_templates(templates), m_return(ret)
    {};

    Type                    get_return() const          { return m_return; };
    const function&         function() const            { return m_function; };
    const type_vec&         deduced() const             { return m_deduced; };
    func_templ_ptr          get_templates() const       { return m_templates; }
    bool                    has_deduced_return() const  { return m_return != nullptr; }
    Type                    get_deduced_return() const  { return m_return; };
    int                     number_templates() const;
    Type                    get_template(int i) const;        
};

class candidate_set
{
    public:
        using templ_vec     = std::vector<Type>;
        using func_vec      = std::vector<func_templ>;
        using func_templ_ptr= func_templ::func_templ_ptr;

    private:
        func_vec            m_overloads_vec;

    public:
        candidate_set(Integer over_num);
        ~candidate_set();

        void				push_back(const function& fun_evl, Type ret, 
                                func_templ_ptr templates, const templ_vec& deduced);
        void				push_back(const func_templ& ft);
        Integer             size() const;
        void                clear();
        const func_templ&   get_function(Integer pos) const;
};

class candidate_type_set
{
    private:
        using type_vec      = std::vector<Type>;

    private:
        type_vec            m_overloads_vec;
        conversion_match    m_match;

    public:
        candidate_type_set(Integer over_num);
        ~candidate_type_set();

        void				add(Type t, conversion_match match);
        conversion_match    get_match() const;
        Integer             size() const;
        void                clear();
        Type                get_type(Integer pos) const;
};

class overload_resolution
{
    private:
        using type_vec          = std::vector<Type>;        

    private:
        e_match_type            m_match_type;
        const function_table*   m_ft;
        int                     n_args;
        const Type*             m_types;
        conversion_sequence     m_best;
        candidate_set&          m_candidates;

    public:
        overload_resolution(const function_table* ft, int n_args, const Type t[],
                            candidate_set& candidates);

        e_match_type        get_match_type() const;
        void                find_best_match(const overload_set& overloads);
        void                find_best_match_templ(const overload_set& overloads, 
                                const function_name_templ* ft, const type_vec& deduced);

        static void         find_most_specialized(const function_table* ft, 
                                candidate_set& func_set);
        static void         find_most_specialized(const function_table* ft, 
                                candidate_type_set& func_set);
        static conversion_match    
                            mach_types(const function_table* ft, Type actual, Type formal);

    private:
        static e_match_type mach_types_list(const function_table* ft, int n_ded, 
                                bool ded_ret, int n_args, const Type t[], const function& f, 
                                conversion_sequence& match);

        static spec_type    is_more_specialized(const function_table* ft, 
                                const func_templ& func_1, const func_templ& func_2);
        static spec_type    is_more_specialized(const function_table* ft, Type type_1,
                                Type type_2);
        Type                eval_return(const function& f, make_return_fptr ret, 
                                const function_name_templ* ft, 
                                const type_vec& deduced, bool& error) const;
        bool                check_deduce_return(const overload_set::fun_ret& evaler) const;
};

};};};
