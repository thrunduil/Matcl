/*
 *  This file is a part of Matrix Computation Library (MATCL)
 *
 *  Copyright (c) Pawe³ Kowal 2017 - 2018
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

#include "matcl-matrep/IO/serialization_helper.h"

namespace matcl { namespace details
{

template<class Base, class Derived>
struct register_serialization_helper
{
    static serialization_helper<Base>* eval(const std::string& unique_name);
};

struct MATCL_MATREP_EXPORT register_serialization_impl
{
    using ptr_type = serialization_helper_base;

    // base_class_name is only used internally, can be given by typeid
    static ptr_type*    get_registered(const std::string& unique_id, 
                            const char* base_class_name);
    static void         make_register(const std::string& unique_id, 
                            const char* base_class_name, ptr_type* ptr);
};

}};

namespace matcl
{

template<class Base_class>
template<class Derived>
serialization_helper<Base_class>* 
serialization_helper<Base_class>::get(const std::string& unique_name)
{
    static serialization_helper* sh = details::register_serialization_helper<Base_class, Derived>
                                            ::eval(unique_name);
    return sh;
};

template<class Base_class>
Base_class* serialization_helper<Base_class>::load(iarchive& ar)
{
    std::string unique_id;
    ar.get() >> unique_id;

    serialization_helper* impl = get_impl(unique_id);
    return impl->make_load(ar);
}

template<class Base_class>
Base_class* serialization_helper<Base_class>::load(std::istream& is)
{
    std::string unique_id;
    is >> unique_id;

    serialization_helper* impl = get_impl(unique_id);
    return impl->make_load(is);
};

template<class Base_class>
serialization_helper<Base_class>* 
serialization_helper<Base_class>::get_impl(const std::string& unique_id)
{
    using ptr_type  = serialization_helper_base;
    ptr_type* ptr   = details::register_serialization_impl
                        ::get_registered(unique_id, typeid(Base_class).name());


    using helper_type = serialization_helper<Base_class>;
    return static_cast<helper_type*>(ptr);
};


};

namespace matcl { namespace details
{

template<class Base_class, class Derived_class>
struct serialization_helper_impl : serialization_helper<Base_class>
{
    private:
        std::string m_unique_id;

    public:
        serialization_helper_impl(const std::string& unique_id)
            :m_unique_id(unique_id)
        {};

        virtual void save(oarchive& ar, const Base_class* ptr) const override
        {
            ar.get() << m_unique_id;
            ptr->save(ar);
        }

        virtual void save(std::ostream& os, const Base_class* ptr) const override
        {
            os << m_unique_id << " ";
            ptr->save(os);
        }

        virtual Base_class* make_load(iarchive& ar) const override
        {
            Base_class* ptr = Derived_class::load(ar);
            return ptr;
        };

        virtual Base_class* make_load(std::istream& is) const override
        {
            Base_class* ptr = Derived_class::load(is);
            return ptr;
        };
};

template<class Base, class Derived>
serialization_helper<Base>* 
    register_serialization_helper<Base,Derived>::eval(const std::string& unique_name)
{
    using helper_type   = serialization_helper_impl<Base,Derived>;
    using ptr_type      = serialization_helper_base;

    ptr_type* helper    = details::register_serialization_impl
                            ::get_registered(unique_name, typeid(Base).name());

    if (helper != nullptr)
    {
        helper_type* ptr = dynamic_cast<helper_type*>(helper);

        if (!ptr)
            throw std::runtime_error("internal error: invalid cast of serialization_helper");

        return ptr;
    };

    helper_type* ret = new serialization_helper_impl<Base,Derived>(unique_name);

    details::register_serialization_impl
            ::make_register(unique_name, typeid(Base).name(), ret);

    return ret;
};

};};
