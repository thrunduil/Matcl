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
#include "matcl-core/memory/refptr.h"
#include "matcl-dynamic/function_name.h"
#include "matcl-dynamic/details/evaler.h"
#include "matcl-dynamic/details/eval_function_impl.h"
#include "matcl-dynamic/type.h"
#include "matcl-dynamic/object.h"

#include <vector>

#pragma warning(push)
#pragma warning(disable : 4251) // needs to have dll-interface to be used by clients

namespace matcl { namespace dynamic
{

// class representing a function pointer
class MATCL_DYN_EXPORT function
{
    private:
        using evaler    = details::evaler;
        using ptr_type  = const evaler*;

    private:
        ptr_type        m_evaler;        

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
        const ptr_type& get_evaler() const;

        // evaluate this function; arguments must be valid, i.e. input 
        // types must be equal to function args types
        void            make(int n_args, const object* args[], object& ret) const;

        // version of make function, when output is not needed
        void            make(int n_args, const object* args[]) const;

    //internal use
    public:
        explicit function(evaler* ev);
};

// evaluate a registered function; 
// if one of the function argument is nonconstant reference, then non contant
// null object can be passed; see also register_function class
struct eval_function
{
    // all Object types must be derived from object or be convertible to object
    // object 'ret' cannot be initialized (otherwise memory leaks are possible)
    template<class ... Object>
    static void eval(const function_name& func, object& ret, Object&& ... args);

    static void eval(const function_name& func, object& ret);
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
        // object 'ret' cannot be initialized (otherwise memory leaks are possible)
        template<class ... Object>
        void eval(const function_name& func, object& ret, Object&& ... args);

        void eval(const function_name& func, object& ret);

        // object 'ret' cannot be initialized (otherwise memory leaks are possible)
        // types and args are arrays of size n_args; types[i] == args[i]->get_type() 
        MATCL_DYN_EXPORT
        void eval_impl(const function_name& func, object& ret, const Type* types,
                    const object ** args, int n_args);
};

};};

#pragma warning(pop)

#include "matcl-dynamic/details/function.inl"

