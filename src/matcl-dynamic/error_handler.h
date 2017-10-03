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
#include "matcl-dynamic/function_name.h"
#include "overload_resolution.h"

#include <vector>
#include <sstream>

namespace matcl { namespace dynamic { namespace details
{

class error_handler
{
    private:
        std::ostringstream  m_msg;
        bool                m_errors;

    public:
        error_handler();
        void    report();

        void    error_one_not_defined(Type t);

        void    error_function_constraints_not_satisfied(function_name func, const function& f);
        void    error_function_defined_with_different_validator(function_name old, function_name newf);
        void    error_function_must_return_template(function_name func, const function& f, 
                    int n_templ, const Type templ[]);

        void    error_function_not_found(function_name func, int n_args, const Type t[]);
        void    error_function_ambiguity(function_name func, int n_args, const Type t[], 
                    const candidate_set& candidates);

        void    error_template_function_not_found(function_name func, int n_templ, const Type templ[],
                    int n_args, const Type t[]);
        void    error_template_function_ambiguity(function_name func, int n_templ, const Type templ[],
                    int n_args, const Type t[], const candidate_set& candidates);

        void    error_invalid_converter_number_args(const function& conv);
        void    error_invalid_unifier_number_args(const function& conv);
        void    error_invalid_assigner_number_args(const function& fun);
        void    error_invalid_assigner_return_type(const function& fun);

        void    error_unable_to_convert(Type to, Type from, converter_type c_type);
        void    error_convert_ambiguity(Type to, Type from, converter_type c_type, 
                    const converter_candidate_set& cs);
        void    error_unable_to_assign(Type to, Type from);
        void    error_assign_ambiguity(Type to, Type from, const converter_candidate_set& as);
        void    error_unify_ambiguity(Type t1, Type t2, const converter_candidate_set& cs);
        void    error_unifiers_ambiguity(Type t1, Type t2, const candidate_type_set& ts);

        void    data_operation_on_reference_type();

    private:
        std::string convert_name(converter_type c_type);

        void    add_error(const std::ostringstream& os);
        void    disp_candidates(std::ostringstream& os, const converter_candidate_set& as);
        void    disp_candidates(std::ostringstream& os, const candidate_type_set& as);
        void    disp_candidates(std::ostringstream& os, const candidate_set& as);        
        void    disp_function_call(std::ostringstream& os, function_name func, int n_args, 
                        const Type t[]);
        void    disp_template_function_call(std::ostringstream& os, function_name func, int n_templ,
                        const Type templ[], int n_args, const Type t[]);
        void    disp_function_declaration(std::ostringstream& os, const function& f, 
                        const std::string& name = "");
        void    disp_function_declaration(std::ostringstream& os, const func_templ& f, 
                        const std::string& name = "");
        void    disp_function_declaration(std::ostringstream& os, const function& f, 
                        Type deduced_return, int n_templ, const Type templ[], 
                        const std::string& name = "");

        void    disp_function_name(std::ostringstream& os, const function& f, Type deduced_return, 
                        const std::string& name);
        void    disp_function_name(std::ostringstream& os, const func_templ& f, const std::string& name);
        void    disp_function_name(std::ostringstream& os, const function& f, Type deduced_return, 
                        int n_templ, const Type templ[], int& n_ded, const std::string& name);
        void    disp_function_arguments(std::ostringstream& os, const function& f, int n_deduced);

        void    disp_type(std::ostringstream& os, Type t);
        void    disp_function_name(std::ostringstream& os, function_name func, int n_templ, 
                        const Type templ[]);
};

};};};
