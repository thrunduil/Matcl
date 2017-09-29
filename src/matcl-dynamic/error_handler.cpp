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

#include "error_handler.h"
#include "matcl-dynamic/details/fwd_decls.h"
#include "matcl-core/general/fwd_decls.h"
#include "matcl-dynamic/type.h"
#include "matcl-dynamic/function.h"
#include "function_table.h"
#include "type_table_cache_data.h"
#include "matcl-dynamic/exception.h"
#include "matcl-dynamic/object_type.h"

#include <vector>
#include <sstream>

namespace matcl { namespace dynamic { namespace details
{

void error_handler::error_function_constraints_not_satisfied(function_name func, function f)
{
    std::ostringstream msg;
    msg << "function does not satisfy constraints";
    msg << "; function declaration: ";
    disp_function_declaration(msg,f, func.to_string());
    msg << "\n";
    msg << "\t" << "requirement: ";
    func.disp_requirements(msg);
    add_error(msg);
};

void error_handler::error_function_defined_with_different_validator
                            (function_name old, function_name newf)
{
    std::ostringstream msg;
    msg << "function " << newf.to_string();
    msg << " is already defined with different requirements";
    msg << "; previous requirement: ";
    old.disp_requirements(msg);
    msg << "; new requirement: ";
    newf.disp_requirements(msg);

    add_error(msg);
};

void error_handler::error_function_must_return_template(function_name func, function f, 
            int n_templ, const Type templ[])
{
    std::ostringstream msg;
    msg << "function must return Template; function declaration is ";
    disp_function_declaration(msg, f, f.return_type(), n_templ, templ, func.to_string());

    add_error(msg);
};

void error_handler::error_function_not_found(function_name func, int n_args, const Type t[])
{
    std::ostringstream msg;
    msg << "unable to call ";
    disp_function_call(msg,func,n_args,t);
    msg << "; no candidate function found";
    add_error(msg);
};

void error_handler::error_template_function_not_found(function_name func, 
            int n_templ, const Type templ[], int n_args, const Type t[])
{
    std::ostringstream msg;
    msg << "unable to call ";
    disp_template_function_call(msg,func,n_templ, templ, n_args,t);
    msg << "; no candidate function found";
    add_error(msg);
};

void error_handler::error_invalid_converter_number_args(function conv)
{
    std::ostringstream msg;
    msg << "invalid converter; expecting function with one argument, function has "
        << conv.number_arguments() << " arguments";
    add_error(msg);
};

void error_handler::error_invalid_unifier_number_args(function conv)
{
    std::ostringstream msg;
    msg << "invalid unifier; expecting function with zero arguments, function has "
        << conv.number_arguments() << " arguments";
    add_error(msg);
};

void error_handler::error_invalid_assigner_number_args(function fun)
{
    std::ostringstream msg;
    msg << "invalid asigner; expecting function with zero arguments, function has "
        << fun.number_arguments() << " arguments";
    add_error(msg);
};

void error_handler::error_invalid_assigner_return_type(function fun)
{
    Type ret = fun.return_type();

    std::ostringstream msg;
    msg << "invalid asigner; expecting function returning void or Unit, return type is "
        << ret.to_string();
    add_error(msg);
};

std::string error_handler::convert_name(converter_type c_type)
{
    switch(c_type)
    {
        case converter_type::conv_implicit:
            return "implicitly convert";
        case converter_type::conv_explicit:
            return "explicitly convert";
        case converter_type::conv_cast:
            return "cast";
        default:
            return "";
    }
}
void error_handler::error_unable_to_convert(Type to, Type from, converter_type c_type)
{
    std::ostringstream msg;
    msg << "unable to " << convert_name(c_type) << " from " 
        << from.to_string() << " to " << to.to_string()
        << "; no candidate converter found";
    add_error(msg);
};

void error_handler::error_one_not_defined(Type t)
{
    std::ostringstream msg;
    msg << "one is not defined for type: " << t.to_string();
    add_error(msg);
}

void error_handler::error_unable_to_assign(Type to, Type from)
{
    std::ostringstream msg;
    msg << "unable to assign " << from.to_string() << " to " << to.to_string()
        << "; no candidate assigner found";
    add_error(msg);
};

void error_handler::error_assign_ambiguity(Type to, Type from, 
                            const converter_candidate_set& as)
{
    std::ostringstream msg;
    msg << "unable to assign " << from.to_string() << " to " << to.to_string()
        << "; ambiguity between: \n";
    disp_candidates(msg, as);
    add_error(msg);
};

void error_handler::error_unify_ambiguity(Type t1, Type t2, 
                            const converter_candidate_set& cs)
{
    std::ostringstream msg;
    msg << "unable to unify type " << t1.to_string() 
        << " and " << t2.to_string()
        << "; ambiguity between converters: \n";
    disp_candidates(msg, cs);
    add_error(msg);
};

void error_handler::error_unifiers_ambiguity(Type t1, Type t2, 
                                const candidate_type_set& ts)
{
    std::ostringstream msg;
    msg << "unable to unify type " << t1.to_string() 
        << " and " << t2.to_string()
        << "; ambiguity between unifiers: \n";
    disp_candidates(msg, ts);
    add_error(msg);
};

void error_handler::error_function_ambiguity(function_name func, 
                int n_args, const Type t[], const candidate_set& candidates)
{
    std::ostringstream msg;
    msg << "unable to call ";
    disp_function_call(msg,func,n_args,t);
    msg << "; ambiguity between function: \n";
    disp_candidates(msg, candidates);
    add_error(msg);
};

void error_handler::error_template_function_ambiguity(function_name func,
        int n_templ, const Type templ[], int n_args, const Type t[], 
        const candidate_set& candidates)
{
    std::ostringstream msg;
    msg << "unable to call ";
    disp_template_function_call(msg,func,n_templ, templ, n_args,t);
    msg << "; ambiguity between function: \n";
    disp_candidates(msg, candidates);
    add_error(msg);
};

void error_handler::error_convert_ambiguity(Type to, Type from, 
                   converter_type c_type, const converter_candidate_set& cs)
{
    std::ostringstream msg;
    msg << "unable to " << convert_name(c_type) << " from " << from.to_string() 
        << " to " << to.to_string()
        << "; ambiguity between: \n";
    disp_candidates(msg, cs);
    add_error(msg);
};

void error_handler::disp_function_call(std::ostringstream& msg, function_name func, 
                                       int n_args, const Type t[])
{
    disp_function_name(msg, func, 0, nullptr);
    msg << " with " << n_args << " arguments";
        
    if (n_args > 0)
        msg << " of types: ";

    for (int i = 0; i < n_args; ++i)
    {
        if (i > 0)
        {
            msg << ", ";
        };

        msg << t[i].to_string();
    };
};

void error_handler::disp_template_function_call(std::ostringstream& msg, 
            function_name func,  int n_templ, const Type templ[], 
            int n_args, const Type t[])
{
    disp_function_name(msg, func, n_templ, templ);
    
    msg << " with " << n_args << " arguments";
    
    if (n_args > 0)
        msg << " of types: ";

    for (int i = 0; i < n_args; ++i)
    {
        if (i > 0)
        {
            msg << ", ";
        };

        msg << t[i].to_string();
    };
};

void error_handler::disp_function_name(std::ostringstream& msg, 
                 function_name func, int n_templ, const Type templ[])
{
    msg << "function " << func.to_string();

    if (n_templ == 0)
        return;

    msg << "<";

    for (int i = 0; i < n_templ; ++i)
    {
        if (i > 0)
        {
            msg << ", ";
        };

        msg << templ[i].to_string();
    };

    msg << ">";
};

error_handler::error_handler()
    :m_errors(false)
{};

void error_handler::report()
{
    if (m_errors == true)
    {
        m_errors = false;
        std::string err = m_msg.str();
        m_msg = std::ostringstream();
        throw error::matcl_dynamic_exception(err);
    };
};

void error_handler::add_error(const std::ostringstream& os)
{
    if (m_errors == false)
    {
        m_msg << os.str();
        m_errors = true;
        return;
    };

    m_msg << "\n\n";
    m_msg << os.str();
};

void error_handler::disp_candidates(std::ostringstream& os, 
                        const converter_candidate_set& as)
{
    Integer n = as.size();

    for (Integer i = 0; i < n; ++i)
    {
        function f = as.get_main_function(i);
        os << "\t";
        disp_function_declaration(os, f);
        os << "\n";
    };
};

void error_handler::disp_candidates(std::ostringstream& os, 
                                   const candidate_set& as)
{
    Integer n = as.size();

    for (Integer i = 0; i < n; ++i)
    {
        const func_templ& f = as.get_function(i);
        os << "\t";
        disp_function_declaration(os, f);
        os << "\n";
    };
};

void error_handler::disp_candidates(std::ostringstream& os, 
                            const candidate_type_set& as)
{
    Integer n = as.size();

    for (Integer i = 0; i < n; ++i)
    {
        Type t = as.get_type(i);
        os << "\t";
        disp_type(os, t);
        os << "\n";
    };
};

void error_handler::disp_function_declaration(std::ostringstream& os, 
                                   function f, const std::string& str)
{
    disp_function_name(os, f, Type(), str);
    disp_function_arguments(os, f, 0);
};

void error_handler::disp_function_declaration(std::ostringstream& os, 
                          const func_templ& f, const std::string& str)
{
    disp_function_name(os, f, str);

    function f0     = f.function();
    int n_deduced   = (int)f.deduced().size() + (f.has_deduced_return() ? 1 : 0);

    disp_function_arguments(os, f0, n_deduced);
};

void error_handler::disp_function_declaration(std::ostringstream& os, 
            const function& f, Type deduced_ret, int n_templ, const Type templ[], 
            const std::string& name)
{
    int n_ded;
    disp_function_name(os, f, deduced_ret, n_templ, templ, n_ded, name);

    if (deduced_ret != Type())
        n_ded   += 1;

    disp_function_arguments(os, f, n_ded);
};

void error_handler::disp_function_arguments(std::ostringstream& os, 
                                      function f, int n_deduced)
{
    //do not print first n_deduced arguments; these are purely technical
    //parameters used only to pass deduced templated to rgistered function

    os << "(";

    Integer n_args = f.number_arguments();

    for (Integer i = n_deduced; i < n_args; ++i)
    {
        Type arg    = f.argument_type(i);

        if (i > n_deduced)
            os << ", ";

        os << arg.to_string();
    };

    os << ")";
};
void error_handler::disp_function_name(std::ostringstream& os, 
            function f, Type deduced_ret, const std::string& name)
{
    if (deduced_ret != Type())
    {
        Type ret    = f.return_type();
        os << "[" << ret.to_string() << " = " << deduced_ret.to_string() << "] ";
    }
    else
    {
        Type ret    = f.return_type();
        os << ret.to_string() << " ";
    };

    if (name.empty() == true)
        os << "function";
    else
        os << name;
};

void error_handler::disp_function_name(std::ostringstream& os, 
                    const func_templ& f, const std::string& name)
{
    disp_function_name(os, f.function(), f.get_deduced_return(), name);

    int num_templ   = f.number_templates();

    if (num_templ == 0)
        return;

    os << "<";

    int deduced_pos = 0;

    for (int i = 0; i < num_templ; ++i)
    {
        if (i > 0)
            os << ", ";

        Type t  = f.get_template(i);
        os      << t.to_string();

        if (t == Template::get_static_type())
        {
            os  << " = " << f.deduced()[deduced_pos].to_string();
            ++deduced_pos;
        };
    }
    
    os << ">";
};

void error_handler::disp_function_name(std::ostringstream& os, 
            function f, Type deduced_return, int n_templ, const Type templ[],
            int& n_ded, const std::string& name)
{
    n_ded = 0;

    disp_function_name(os, f, deduced_return, name);

    if (n_templ == 0)
        return;

    os << "<";    

    for (int i = 0; i < n_templ; ++i)
    {
        if (i > 0)
            os << ", ";

        Type t  = templ[i];
        os      << t.to_string();

        if (t == Template::get_static_type())
            ++n_ded;
    }
    
    os << ">";
};
void error_handler::disp_type(std::ostringstream& os, Type t)
{
    os << t.to_string();
};

};};};
