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

#include "matcl-dynamic/object_type_traits.h"
#include "matcl-dynamic/function.h"
#include "matcl-dynamic/details/register_function_impl.h"

namespace matcl { namespace dynamic { namespace details
{

template<class T>
type_impl::type_impl(const mark_type<T>& mt)                    
    :m_name(mt.name()), m_has_trivial_assignment(true)
{};

inline type_impl::type_impl(const std::string& name)
    :m_name(name), m_has_trivial_assignment(true)
{};

template<class T, bool is_clonable>
struct clone_helper
{
    static T eval(const T& val) { return val; };
};
template<class T>
struct clone_helper<T,true>
{
    static T eval(const T& val) { return val.clone(); };
};

template<class T>
type_impl::data_type* type_object<T>::clone(const data_type* obj) const
{
    static const bool is_clonable = object_type_traits<T>::is_clonable;
    return object_data<T>::create(clone_helper<T,is_clonable>::eval(obj->get_value<T>()));
};

template<class T>
type_impl::data_type* type_object<T>::copy(const object_data_base* obj) const
{
    return object_data<T>::create(obj->get_value<T>());
};

template<class T>
bool type_object<T>::is_zero(const data_type* obj_data) const
{
    return object_type_traits<T>::is_zero(obj_data->get_value<T>());
};

template<class T>
struct make_predef_functions
{
    using object_t  = object_type<T>;

    static OBool eval_is_zero(const object_t& arg)
    {
        return OBool(arg.is_zero());
    };
    static OBool eval_is_one(const object_t& arg)
    {
        return OBool(arg.is_one());
    };

    static function is_zero()
    {
        using Fun = OBool (*)(const object_t&);
        return function_register<Fun>::process_function(&eval_is_zero);
    }
    static function is_one()
    {
        using Fun = OBool (*)(const object_t&);
        return function_register<Fun>::process_function(&eval_is_one);
    }
};

template<class T>
function type_object<T>::generate_function(predef_fun fun) const
{
    switch (fun)
    {
        case predef_fun::is_zero:
            return make_predef_functions<T>::is_zero();
        case predef_fun::is_one:
            return make_predef_functions<T>::is_one();
        default:
            //unknown function
            return function();
    }
}

template<class T>
bool type_object<T>::is_one(const data_type* obj_data) const
{
    return object_type_traits<T>::is_one(obj_data->get_value<T>());
};

template<class T>
std::string type_object<T>::to_string(const data_type* obj_data, printer& pr) const
{
    return object_type_traits<T>::to_string(obj_data->get_value<T>(),pr);
};

template<class T>  
void type_object<T>::disp(const data_type* obj_data, printer& pr, Integer elem_width, 
                        align_type at, Integer value_pos) const
{
    return object_type_traits<T>::disp(obj_data->get_value<T>(), pr, elem_width, at, value_pos);
};

template<class T>
type_impl::data_type* type_object<T>::create() const
{
    return object_data<T>::create(object_type_traits<T>::make_zero<T>());
};

#pragma warning (push)
#pragma warning (disable: 4127) // conditional expression is constant

template<class T>
type_impl::data_type* type_object<T>::create_one() const
{
    if (object_type_traits<T>::has_one == true)
        return object_data<T>::create(object_type_traits<T>::make_one((T*)nullptr));
    else 
        return nullptr;
};

#pragma warning (pop)

template<class T>
bool type_object<T>::has_one() const
{
    return object_type_traits<T>::has_one;
};

template<class T>
void type_object<T>::save(oarchive_impl& ar, unsigned int version, 
                        const data_type* obj_data) const
{
    const object_data<T>* real_obj_data = reinterpret_cast<const object_data<T>*>(obj_data);

    object_type_traits<T>::save_data(ar,real_obj_data->get(), version);
};

template<class T>
type_impl::data_type* type_object<T>::load(iarchive_impl& ar, unsigned int version) const
{
    T val;
    object_type_traits<T>::load_data(ar,val,version);
    return object_data<T>::create(std::move(val));
};

template<class T>
void type_object<T>::save(std::ostream& os, const data_type* data) const
{
    const object_data<T>* real_obj_data = static_cast<const object_data<T>*>(data);
    object_type_traits<T>::write(os, real_obj_data->get());
};

template<class T>
type_impl::data_type* type_object<T>::load(std::istream& is) const
{
    T val;
    object_type_traits<T>::read(is, val);
    return object_data<T>::create(val);
};

};};};
