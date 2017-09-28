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
#include "matcl-core/general/fwd_decls.h"
#include "matcl-dynamic/function_name.h"
#include "matcl-dynamic/details/evaler.h"
#include "matcl-dynamic/type.h"
#include "matcl-dynamic/object.h"

#include <vector>

namespace matcl { namespace dynamic
{

// class representing a function pointer
class MATCL_DYN_EXPORT function
{
    private:
        using evaler    = details::evaler;

        evaler*     m_evaler;        

    public:
        // construct uninitialized object
        function();

        // return true if this object is uninitialized
        bool            is_null() const;

        // return number of input arguments
        int             number_arguments() const;

        // get type of narg argument (0-based)
        Type            argument_type(int narg) const;

        // get return type
        Type            return_type() const;

        // return function evaler
        const evaler*   get_evaler() const;

        // make function with conversions; returning function first evaluates
        // converters for each input argument and then this function;
        // return type of n-th converter must be the same as n-th argument
        // type of this function; deduced_ret is the deduced return type
        // and is different from Type() if return deduction is enabled;
        // n_deduced is the number of deduced template types with deduced 
        // types given by deduced; deduced return (if not empty) and deduced
        // template arguments must be passed as first m + n_deduced arguments
        // of type OType, where m is 1 if deduced_ret is not Type() and zero
        // otherwise; deduced return must be passed first
        function        make_converter(int n_deduced, const Type deduced[],
                            Type deduced_ret, const std::vector<function>& arg_converters) const;

        // evaluate this function; arguments must be valid, i.e. input 
        // types must be equal to function args types
        object          make(int n_args, const object* args[]) const;

    //internal use
    public:
        function(evaler* ev);
};

// evaluate a registered function; 
// if one of the function argument is nonconstant reference, then non contant
// null object can be passed; see also register_function class
struct eval_function
{
    // all Object types must be derived from object or be convertible to object
    template<class ... Object>
    static object eval(const function_name& func, Object&& ... args);

    static object eval(const function_name& func);
};

// evaluate a registered template function; 
struct eval_function_template
{
    private:
        std::vector<Type>   m_types;

    public:
        // list of template arguments
        eval_function_template(std::initializer_list<Type> types);

        // all Object types must be derived from object or be convertible to object
        template<class ... Object>
        object eval(const function_name& func, Object&& ... args);

        object eval(const function_name& func);
};

};};

#include "matcl-dynamic/details/function.inl"