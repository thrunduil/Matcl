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
#include "matcl-core/config.h"
#include "matcl-dynamic/details/fwd_decls.h"
#include "matcl-core/general/fwd_decls.h"
#include "matcl-dynamic/object.h"
#include "matcl-dynamic/details/utils.h"

namespace matcl { namespace dynamic
{

// Unit type, that store exactly one value; 
// every type is converible to Unit type
class MATCL_DYN_EXPORT unit_type
{
    public:
        // default constructor
        unit_type()     {};

        // conversion from other types
        template<class T>
        unit_type(T&&)  {};

        // assignment from other types
        template<class T>
        unit_type& operator=(T&& other) & { (void)other; return *this; };
};

// save and load of unit_type to stream
MATCL_DYN_EXPORT std::ostream& operator<<(std::ostream& os, const unit_type& );
MATCL_DYN_EXPORT std::istream& operator>>(std::istream& is, unit_type& );

// null type, that store nullptr
class MATCL_DYN_EXPORT null_type
{
    public:
        // default constructor
        null_type()                     {};

        // cast to bool
        explicit operator bool() const  { return false; };
};

// comparison operators
inline bool operator==(const null_type&, const null_type&)  { return true; };
inline bool operator!=(const null_type&, const null_type&)  { return false; };
inline bool operator>=(const null_type&, const null_type&)  { return true; };
inline bool operator<=(const null_type&, const null_type&)  { return true; };
inline bool operator>(const null_type&, const null_type&)   { return false; };
inline bool operator<(const null_type&, const null_type&)   { return false; };

// save and load of unit_type to stream
MATCL_DYN_EXPORT std::ostream& operator<<(std::ostream& os, const null_type& );
MATCL_DYN_EXPORT std::istream& operator>>(std::istream& is, null_type& );

// Any type; every value can be stored in this type;
// every type is convertible to this type
class MATCL_DYN_EXPORT any_type
{
    private:
        object      m_stored;

    public:
        // construct uninitialized object
        any_type();

        // convert object to Any type
        any_type(const object& other);
        any_type(object&& other);

        // convert typed object of type any_type to Any type
        any_type(const dynamic::object_type<any_type>& other);
        any_type(dynamic::object_type<any_type>&& other);

        // convert typed object to Any type
        template<class T>
        any_type(const dynamic::object_type<T>& other);
        template<class T>
        any_type(dynamic::object_type<T>&& other);

        // convert every other value to Any type
        template<class S, class Enable = typename details::enable_if_nonobject_any<S,void*>::type>
        any_type(S&& other);

        // copy and move constructors
        any_type(const any_type& other);
        any_type(any_type&& other);

        // destructor
        ~any_type();

        // assign from object
        any_type&       operator=(const object& other) &;
        any_type&       operator=(object&& other) &;

        // assign from typed object of type any_type
        any_type&       operator=(const dynamic::object_type<any_type>& other) &;
        any_type&       operator=(dynamic::object_type<any_type>&& other) &;

        // assign from typed object
        template<class T>
        any_type&       operator=(const dynamic::object_type<T>& other) &;
        template<class T>
        any_type&       operator=(dynamic::object_type<T>&& other) &;

        // assign from any other type
        template<class S, class Enable = typename details::enable_if_nonobject_any<S,void*>::type>
        any_type&       operator=(S&& other) &;

        //standard assign and move assign
        any_type&       operator=(const any_type& other) &;
        any_type&       operator=(any_type&& other) &;

        // get stored value
        const object&   get_stored() const;

        // check if stored value is uninitialized object
        bool            is_zero() const;

        // make independent copy
        any_type        clone() const;

        // equality and inequality comparison; two objects are considered
        // equal if stored elements are the same (i.e. pointers to data are
        // equal), thus any_type(1) != any_type(1) but 
        //         object i(1);
        //         any_type(i) == any_type(i)
        bool            operator==(const any_type& other) const;
        bool            operator!=(const any_type& other) const;

        // return true if object is initialized
        explicit        operator bool() const;

        // serialization and deserialization
        void	        serialize(iarchive_impl& ar, unsigned int version);
        void	        serialize(oarchive_impl& ar, unsigned int version) const;

    //internal use
    public:
        // display helper
        void            disp(matcl::details::printer& pr, Integer elem_width, 
                            align_type at, Integer value_pos) const;

    private:
        void            remove_indirections();
};

// save and load of any_type to stream
MATCL_DYN_EXPORT std::ostream& operator<<(std::ostream& os, const any_type& );
MATCL_DYN_EXPORT std::istream& operator>>(std::istream& is, any_type& );

// this type can only be used internally
template<class T>
class object_reference
{
    public:
        // get static type of this object (i.e. a reference type to T)
        static Type         get_static_type();
};

};};

#include "matcl-dynamic/details/special_types.inl"