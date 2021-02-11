/*
 *  This file is a part of Matrix Computation Library (MATCL)
 *
 *  Copyright (c) Pawe³ Kowal 2017 - 2021
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
#include "matcl-dynamic/type.h"
#include "matcl-dynamic/exception.h"
#include "matcl-core/matrix/enums.h"

namespace matcl { namespace dynamic
{

// class storing elements of any type with runtime type
// information; this type is reference counted
class MATCL_DYN_EXPORT object
{
    private:
        using data_type         = details::object_data_base;

    private:
        Type                    m_type;
        data_type*              m_data;

    public:
        // construct uninitialized object
        object();        

        // convert existing object other to object of type t
        object(Type t, const object& other);
        object(Type t, object&& other);

        // construct default initialized object of type t
        explicit object(Type t);

        // construct of type t representing value 1
        static object make_one(Type t);

        // copy and move constructor; data is not copied, only
        // reference count is increased
        object(const object& other);
        object(object&& other);

        // conversion from typed object
        template<class T>
        explicit object(const object_type<T>& other);

        // conversion from temporary typed object
        template<class T>
        explicit object(object_type<T>&& other);

        // convert predefined scalar types to object type
        explicit object(bool);
        explicit object(Integer);
        explicit object(Float);
        explicit object(Real);
        explicit object(const Float_complex&);
        explicit object(const Complex& );
        explicit object(const std::string& );
        explicit object(std::string&& );
        explicit object(const char* str);

        // destructor
        ~object();        

        // call assign function (if defined) or convert other to type of given
        // object and then perform standard assignment; type of given object
        // does not change
        object&                 operator=(const object& other) &;
        object&                 operator=(object&& other) &;
        
        // change given object to other; type may change
        void                    reset(const object& other) &;
        void                    reset(object&& other) &;

        // make independent copy; a clone function is called for clonable
        // objects (i.e. registered as clonable in object_type_traits), otherwise
        // the copy constructor is called
        object                  clone() const;

        // return true if this object is unique
        bool                    is_unique() const;

        // make this object unique
        void                    make_unique();

        // get runtime type of this object
        Type                    get_type() const        { return m_type; };

        // get static type of this object
        static Type             get_static_type()       { return predefined::type_any(); };

        // check if this object is initialized
        bool                    is_null() const         { return m_type == Type(); };

        // check if this object is zero (i.e. an identity element of addition)
        bool                    is_zero() const;

        // check if this object is one (i.e. an identity element of multiplication)
        bool                    is_one() const;

        // convert to boolean value
        explicit                operator bool() const;

        // serialization and deserialization
        void                    serialize(oarchive_impl & ar, const unsigned int version) const;
        void                    serialize(iarchive_impl & ar, const unsigned int version);

        // explicitly cast object obj to object of type new_type
        friend MATCL_DYN_EXPORT 
        object                  cast(Type new_type, const object& obj);

        // explicitly convert object obj to object of type new_type
        friend MATCL_DYN_EXPORT 
        object                  convert(Type new_type, const object& obj);

        // save to stream
        friend MATCL_DYN_EXPORT 
        std::ostream&	        operator<<(std::ostream& os, const object& A);

        // load from stream
        friend MATCL_DYN_EXPORT 
        std::istream&	        operator>>(std::istream& is, object& A);

        // save data to stream; type info is not saved
        friend MATCL_DYN_EXPORT 
        std::ostream&	        save_data(std::ostream& os, const object& A);

        // load data from stream; type info is not loaded and runtime type of A 
        // must be the same as type of data
        friend MATCL_DYN_EXPORT 
        std::istream&	        load_data(std::istream& is, object& A);

        // swap two objects
        friend MATCL_DYN_EXPORT 
        void                    swap(object& o1, object& o2);

    //internal use
    public:
        void                    disp(matcl::details::printer& pr, Integer elem_width, 
                                     matcl::align_type at, Integer value_pos) const; 

        const data_type*        get_data() const        { return m_data; };		
        data_type*              get_data_unique()       { make_unique(); return m_data; };	

    protected:
        struct not_null{};
        struct null{};

        // ti and m_data are not null
        object(Type ti, data_type* m_data, not_null);

        // ti and m_data are null
        object(Type ti, null);

        template<class T>
        friend class object_type;

    private:
        void                    make_assignment(const object& other, Type ty1, Type ty2);
        void                    make_assignment(object&& other, Type ty1, Type ty2);
};

// matcl-dynamic library initializer
struct MATCL_DYN_EXPORT object_initializer
{
    object_initializer();
    ~object_initializer();
};

static object_initializer g_object_initializer;

};};
