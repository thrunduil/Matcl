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

#include "matcl-dynamic/object.h"
#include "matcl-dynamic/predefined_type_traits.h"
#include "matcl-dynamic/special_types.h"
#include "matcl-dynamic/details/utils.h"

namespace matcl { namespace dynamic
{

// construct object_type from an object
struct from_object{};

// object storing elements of type T
template<class T>
class object_type
{
    private:
        object  m_data;

    public:
        // type of stored data
        using value_type            = T;

        // has_one is true if type T has one value
        static const bool has_one   = object_type_traits<T>::has_one;

    public:
        // construct an object storing default value of type T
        object_type();
        
        // convert object other to object of given type by calling
        // registered explicit conversion for objects
        object_type(const object& other, from_object);
        object_type(object&& other, from_object);

        // convert objects of type S to objects of type T;
        // implicit conversion is enabled if there exist implicit conversion S->T
        template<class S>
        object_type(const object_type<S>& other,
                    typename details::enable_if_different_conv<S,T,true,void*>::type = 0);

        // convert objects of type S to objects of type T;
        // explicit conversion is enabled if there exist explicit conversion S->T
        template<class S>
        explicit object_type(const object_type<S>& other,
                    typename details::enable_if_different_conv<S,T,false,void*>::type = 0);

        // convert element of type S to object of type T
        // explicit conversion is enabled if there exist explicit conversion S->T
        template<class S>
        explicit object_type(S&& other,
                    typename details::enable_if_nonobject_conv<S,T,true,void*>::type = 0);
        template<class S>
        explicit object_type(S&& other,
                    typename details::enable_if_nonobject_conv<S,T,false,void*>::type = 0);

        // assignment from object is enabled if there exist an 
        // assignment T = object
        template<class Enable = typename details::enable_if_conv<object,T,true,void*>::type>
        object_type(const object& other);
        
        template<class Enable = typename details::enable_if_conv<object,T,true,void*>::type>
        object_type(object&& other);

        // call constructor of type T with more than one arguments arg1, arg2, ...
        // and construct object of type T storing this value
        template<class Arg1, class Arg2, class ... Args,
                class Enable = typename details::enable_if_not<Arg2, from_object, void>::type>
        object_type(Arg1&& arg1, Arg2&& arg2, Args&& ... args);

        //copy constructors
        object_type(const object_type<T>& other);
        object_type(object_type<T>&& other);

        // assignment from object is enabled if there exist an 
        // assignment T = object
        template<class Enable = typename details::enable_if_assign<object,T,void*>::type>
        object_type&        operator=(const object& other) &;

        // assignment from object is enabled if there exist an 
        // assignment T = object
        template<class Enable = typename details::enable_if_assign<object,T,void*>::type>
        object_type&        operator=(object&& other) &;

        // assignments from objects of type T;
        object_type&        operator=(const object_type<T>& rhs) &;

        // assignments from objects of type S;
        // this assignment is enabled if there exist an assignment T = S
        template<class S, class Enable = typename details::enable_if_assign<S,T,void*>::type>
        object_type&        operator=(const object_type<S>& rhs) &;

        // change given object to other; assign operator defined for type T is 
        // not called
        void                reset(const object_type& other) &;
        void                reset(object_type&& other) &;

        // get stored value
        const value_type&   get() const;

        // get nonconstant reference to stored value; make this object unique
        value_type&			get_unique();

        // const acces to stored data
        const value_type*   operator->() const;

        // make independent copy; for clonable objects, i.e. objects registered
        // as a clonable type in object_type_traits a clone function defined in
        // type T is called; otherwise call copy constructor
        object_type         clone() const;

        // return true if this object is unique
        bool                is_unique() const;

        // make this object unique
        void                make_unique();

        // get runtime type of this object
        Type                get_type() const;

        // get static type of this object
        static Type         get_static_type();

        // check if this object is initialized
        bool                is_null() const;

        // check if this object is zero (i.e. an identity element of addition)
        bool                is_zero() const;

        // check if this object is one (i.e. an identity element of multiplication)
        bool                is_one() const;

        // check if given object is zero (i.e. an identity element of addition)
        static bool         is_zero(const T& arg);

        // check if given object is one (i.e. an identity element of multiplication)
        static bool         is_one(const T& arg);

        // create value one; return default value if has_one is false
        static object_type  make_one();

        // convert to boolean value; use conversion operator defined for type T
        explicit            operator bool() const;

        // convert object this object to object of type S by calling 
        // registered conversion function
        template<class S>
        object_type<S>      convert() const;

        // convert object this object to object of type S by calling 
        // registered cast function
        template<class S>
        object_type<S>      cast() const;

        // convert to object; conversion to nonconstant reference to object
        // is not allowed
        explicit            operator const object&() const;

        // serialization and deserialization
        void                serialize(oarchive_impl & ar, const unsigned int version) const;
        void                serialize(iarchive_impl & ar, const unsigned int version);

        // save to stream;
        template<class T> friend 
        std::ostream&       operator<<(std::ostream& os, const object_type<T>& );

        // load from stream;
        template<class T> friend 
        std::istream&       operator>>(std::istream& is, object_type<T>& A);

        // swap two objects
        template<class T> friend
        void                swap(object_type<T>& o1, object_type<T>& o2);
};

// object of unit type
using Unit              = object_type<unit_type>;

// object of any type
using Any               = object_type<any_type>;

// object of bool type
using ONull             = object_type<null_type>;

// object of bool type
using OBool             = object_type<bool>;

// object of integer type
using OInteger          = object_type<Integer>;

// object of double precision floating point type
using OReal             = object_type<Real>;

// object of single precision floating point type
using OFloat            = object_type<Float>;

// object of double precision complex type
using OComplex          = object_type<Complex>;

// object of single precision complex type
using OFloat_complex    = object_type<Float_complex>;

// object of string type
using OString           = object_type<std::string>;

// special type, that can be only used to mark template
// arguments while registering functions
using Template          = object_type<object>;

// special type, that represents a Type as an object
using OType             = object_type<Type>;

};};

#include "matcl-dynamic/details/object_type.inl"
