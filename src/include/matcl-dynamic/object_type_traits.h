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
#include "matcl-dynamic/config.h"
#include "matcl-core/config.h"
#include "matcl-core/matrix/complex_type.h"
#include "matcl-core/matrix/enums.h"

#include <string>

namespace matcl { namespace dynamic
{

// default implementation of functions required by object class
struct object_type_traits_default
{
    // check if given type has value the identity element of multiplication
    static const bool has_one       = false;

    // if true, then mark this type as clonable
    static const bool is_clonable   = false;

    // check if given value represents the identity element of multiplication
    template<class T>
    static bool         is_one(const T&)                    { return false; };

    // construct the identity element of multiplication; this function
    // is called only if has_one is true; mark_type argument is not
    // referred
    template<class T>
    static T            make_one(const T* mark_type)        { (void)mark_type; return T(); };

    // check if given value represents the identity element of addition
    template<class T>
    static bool         is_zero(const T& v)                 { return v == T(); };

    // construct the identity element of addition
    template<class T>
    static T            make_zero()                         { return T(); };

    // convert to string; predefined scalar types should be 
    // converted to string using functions from printer pr
    template<class T>
    static std::string  to_string(const T& t, printer& pr)  { return t.to_string(pr); };

    // display element
    // if elem_width > 0, then width is fixed to elem_width characters; 
    // if elem_width < 0, then width is not fixed, and width gives maximal width
    // if elem_width == 0, then width is not specified
    // at is expected alignment type
    // value_pos is an index of subvalue (used only if given object is a composition
    // of many subobjects; value_pos subobject should be displayed
    template<class T>
    static void         disp(const T& t, printer& pr, Integer elem_width,
                             align_type at, Integer value_pos);

    // save and load to stream
    template<class T>
    static bool         read(std::istream&, T& t);

    template<class T>
    static void         write(std::ostream&, const T& t);

    // serialization and deserialization
    template<class T>
    static void	        load_data(iarchive_impl& ar, T& ret, unsigned int version);

    template<class T>
    static void	        save_data(oarchive_impl& ar, const T& val, unsigned int version);
};

// implementation of functions required by object class; this class
// should be specialized for a given type
template<class T>
struct object_type_traits : object_type_traits_default
{};

};};

#include "matcl-dynamic/details/object_type_traits.inl"
