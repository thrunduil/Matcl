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

#include "matcl-dynamic/details/fwd_decls.h"
#include "matcl-core/config.h"
#include "matcl-core/general/fwd_decls.h"
#include "matcl-dynamic/function_name.h"

namespace matcl { namespace dynamic
{

// class representing type of an object
class MATCL_DYN_EXPORT Type 
{
    private:
        using type_impl    = matcl::dynamic::details::type_impl;

        const type_impl*   m_impl;
        friend type_impl;        

    public:
        // create undefined type
        Type();                

        // return hash value
        size_t              hash_value() const              { return (size_t)m_impl; };

        // equality and inequality comparison
        bool                operator==(Type other) const    { return m_impl == other.m_impl;};
        bool                operator!=(Type other) const    { return m_impl != other.m_impl;};

        // less than comparison based on pointer address comparison;
        // note that result may be different in different program run
        bool                operator<(Type other) const     { return m_impl < other.m_impl;};
        bool                operator>(Type other) const     { return m_impl > other.m_impl;};
        bool                operator<=(Type other) const    { return m_impl <= other.m_impl;};
        bool                operator>=(Type other) const    { return m_impl >= other.m_impl;};

        // return name of this Type
        std::string         to_string() const;

        // serialization and deserialization
        void                serialize(matcl::oarchive_impl & ar, unsigned int version) const;
        void                serialize(matcl::iarchive_impl & ar, unsigned int version);   

    public:
        // return true if t is a reference type
        static bool         is_reference(Type t);

        // equivalent to std::decay<T>
        static Type         decay(Type t);

    //internal use
    public:
        Type(const type_impl* impl);

        const type_impl*    get_impl() const    { return m_impl; };
};

// save and load Type to stream
MATCL_DYN_EXPORT std::ostream&  operator<<(std::ostream& os, Type);
MATCL_DYN_EXPORT std::istream&  operator>>(std::istream& is, Type&);

namespace predefined
{
    // predefined types

    // null type, equivalent to Type()
    MATCL_DYN_EXPORT Type   type_null();
    // unit type
    MATCL_DYN_EXPORT Type   type_unit();
    // any type
    MATCL_DYN_EXPORT Type   type_any(); 
    // bool type
    MATCL_DYN_EXPORT Type   type_bool();
    // int type
    MATCL_DYN_EXPORT Type   type_int(); 
    // double precision floating point type
    MATCL_DYN_EXPORT Type   type_real();
    // single precision floating point type
    MATCL_DYN_EXPORT Type   type_float();
    // double precision complex type
    MATCL_DYN_EXPORT Type   type_complex();
    // single precision complex type
    MATCL_DYN_EXPORT Type   type_float_complex();
    // string type
    MATCL_DYN_EXPORT Type   type_string();
};

namespace operations
{

// check if value one is defined for type t
MATCL_DYN_EXPORT bool       has_one(Type t);

// check if assignment lhs = rhs is trivial
MATCL_DYN_EXPORT bool       has_trivial_assignment(Type lhs, Type rhs);

// return a most specialized type T, such that t1 and t2 are converible to T
MATCL_DYN_EXPORT Type       unify_types(Type t1, Type t2);

// make type of nonconst reference to t
MATCL_DYN_EXPORT Type       make_reference_type(Type t);

// get function overload for functions with name func based on types of 
// input arguments
MATCL_DYN_EXPORT 
function                    get_overload(const function_name& func, int n_args,
                                const Type t[]);

// get function overload for functions with name func based on types of 
// input arguments and template arguments
MATCL_DYN_EXPORT 
function                    get_template_overload(const function_name& func, int n_templ, 
                                const Type templates[], int n_args, const Type arg_types[]);

// find return type of function overload for functions with name func based
// on types of input arguments
MATCL_DYN_EXPORT Type       return_type(const function_name& func, int n_args, const Type t[]);

};

};};

