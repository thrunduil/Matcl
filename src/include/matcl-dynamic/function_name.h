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

#include "matcl-dynamic/config.h"
#include "matcl-dynamic/details/fwd_decls.h"
#include "matcl-core/config.h"

#include <string>

#pragma warning(push)
#pragma warning(disable:4275)	//used as base for dll-interface class
#pragma warning(disable:4251)	//needs to have dll-interface to be used by clients

namespace matcl { namespace dynamic
{

// class repesenting an identifier; identifiers with the same string
// are given by the same object; 
class MATCL_DYN_EXPORT identifier
{
    private:
        using impl_type = details::identifier_impl;

    protected:
        impl_type*      m_impl;
        size_t          m_hash;
        size_t          m_code;

    public:
        // construct an identifier from string id
        identifier(const std::string& id);

        // return hash value
        size_t          hash_value() const;

        // convert to string
        std::string     to_string() const;

        // return unique code associated to this identifier
        size_t          get_unique_code() const;

        // equality and inequality comparison
        bool            operator==(identifier other) const;
        bool            operator!=(identifier other) const;

        // less than comparison based on pointer address comparison;
        // note that result may be different in different program run
        bool            operator<(identifier other) const;
        bool            operator>(identifier other) const;
};

// check if function satisfies requirements
struct function_validator
{
    public:
        using validator = bool (*)(const function&);
        using printer   = void (*)(std::ostream& os);

    private:
        validator       m_validator;
        printer         m_printer;

    public:
        // default constructor, no constraints on function
        function_validator();

        // validator function v checks if function satisfies constraints
        // printer function p describes these constraints
        function_validator(validator v, printer p);

        // return validator function
        validator       get_validator() const;

        // return printer function
        printer         get_printer() const;

        // equality and inequality comparison based on function pointers
        // comparison
        bool            operator==(function_validator o) const;
        bool            operator!=(function_validator o) const;
};

// functions name
class MATCL_DYN_EXPORT function_name : public identifier
{
    private:
        using validator = function_validator;

    private:
        validator       m_validator;

    public:
        // constructor taking function name and function validator
        function_name(const std::string& name, validator valid = validator());

        // return function validator
        validator       get_validator() const;

        // check if function f satisfies requirements defined by the
        // validator
        bool            validate_function(const function& f) const;

        // display function requirements
        void            disp_requirements(std::ostream& os) const;                

        // return true if this is internal function
        bool            is_special() const;
};

};};

#pragma warning(pop)

#include "matcl-dynamic/details/function_name.inl"