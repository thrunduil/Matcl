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

#include "matcl-dynamic/details/register_object.inl"
#include "matcl-dynamic/details/object_data.inl"
#include "matcl-dynamic/predefined_functions.h"
#include "matcl-dynamic/typed_object_functions.h"

namespace matcl { namespace dynamic
{

template<class T>
force_inline
object_type<T>::object_type()
    :m_data(dynamic::details::mark_type<T>().get())
{};

template<class T>
force_inline
object_type<T>::object_type(const object& other, from_object)
    :m_data(dynamic::details::mark_type<T>().get(), other)
{};

template<class T>
force_inline
object_type<T>::object_type(object&& other, from_object)
    :m_data(dynamic::details::mark_type<T>().get(), std::move(other))
{};

template<class T>
template<class Enable>
force_inline
object_type<T>::object_type(const object& other)
    :m_data(dynamic::details::mark_type<T>().get(), other)
{};

template<class T>
template<class Enable>
force_inline
object_type<T>::object_type(object&& other)
    :m_data(dynamic::details::mark_type<T>().get(), std::move(other))
{};

template<class T>
force_inline
object_type<T>::object_type(const object_type<T>& arg)
    :m_data(arg)
{};

template<class T>
force_inline
object_type<T>::object_type(object_type<T>&& arg)
    :m_data(std::move(arg.m_data))
{};

template<class T>
template<class S>
force_inline
object_type<T>::object_type(const object_type<S>& other, 
                            typename details::enable_if_different_conv<S,T,true,void*>::type)
    :m_data(dynamic::details::mark_type<T>().get(),
            details::object_data<T>::create(T(other.get())),
            object::not_null())
{};

template<>
template<class S>
force_inline
object_type<any_type>::object_type(const object_type<S>& other, 
                            typename details::enable_if_different_conv<S,any_type,true,void*>::type)
    :m_data(dynamic::details::mark_type<any_type>().get(), object(other))
{};

template<class T>
template<class S>
force_inline
object_type<T>::object_type(const object_type<S>& other, 
                            typename details::enable_if_different_conv<S,T,false,void*>::type)
    :m_data(dynamic::details::mark_type<T>().get(),
            details::object_data<T>::create(T(other.get())),
            object::not_null())
{};

template<>
template<class S>
force_inline
object_type<any_type>::object_type(const object_type<S>& other, 
                            typename details::enable_if_different_conv<S,any_type,false,void*>::type)
    :m_data(dynamic::details::mark_type<any_type>().get(), object(other))
{};

template<class T>
template<class S>
force_inline
object_type<T>::object_type(S&& arg, typename details::enable_if_nonobject_conv<S,T,true,void*>::type)
    :m_data(dynamic::details::mark_type<T>().get(),
            details::object_data<T>::create(T(std::forward<S>(arg))),
            object::not_null())
{};

template<class T>
template<class S>
force_inline
object_type<T>::object_type(S&& arg, typename details::enable_if_nonobject_conv<S,T,false,void*>::type)
    :m_data(dynamic::details::mark_type<T>().get(),
            details::object_data<T>::create(T(std::forward<S>(arg))),
            object::not_null())
{};

template<class T>
template<class Arg1, class Arg2, class ... Args, class Enable>
force_inline
object_type<T>::object_type(Arg1&& arg1, Arg2&& arg2, Args&& ... args)
    :m_data(dynamic::details::mark_type<T>().get(),
            details::object_data<T>::create(T(std::forward<Arg1>(arg1), std::forward<Arg2>(arg2), 
                                       std::forward<Args>(args)...)),
            object::not_null())
{};

template<class T>
force_inline
object_type<T>& object_type<T>::operator=(const object_type<T>& arg) &
{
    this->get_unique() = arg.get();
    return *this;
};

template<class T>
template<class S, class Enable>
force_inline
object_type<T>& object_type<T>::operator=(const object_type<S>& arg) &
{
    this->get_unique() = arg.get();
    return *this;
};

template<>
template<class S, class Enable>
force_inline
object_type<any_type>& object_type<any_type>::operator=(const object_type<S>& arg) &
{
    this->get_unique() = arg;
    return *this;
};

template<class T>
template<class Enable>
force_inline
object_type<T>& object_type<T>::operator=(const object& arg) &
{
    this->get_unique() = arg;
    return *this;
};

template<class T>
template<class Enable>
force_inline
object_type<T>& object_type<T>::operator=(object&& arg) &
{
    this->get_unique() = std::move(arg);
    return *this;
};

template<class T>
force_inline
const typename object_type<T>::value_type& object_type<T>::get() const
{
    return static_cast<const details::object_data<T>*>(m_data.get_data())->get();
};

template<class T>
force_inline
const typename object_type<T>::value_type* object_type<T>::operator->() const
{
    return &get();
};

template<class T>
force_inline
typename object_type<T>::value_type& object_type<T>::get_unique()
{
    return static_cast<details::object_data<T>*>(m_data.get_data_unique())->get();
};

template<>
force_inline
const object& object_type<object>::get() const
{
    return m_data;
};

template<>
force_inline
object& object_type<object>::get_unique()
{
    m_data.make_unique();
    return m_data;
};

template<class T>
force_inline
void object_type<T>::reset(const object_type& other) &
{
    m_data.reset(other.m_data);
}

template<class T>
force_inline
void object_type<T>::reset(object_type&& other) &
{
    m_data.reset(std::move(other.m_data));
};

template<class T>
object_type<T> object_type<T>::clone() const
{
    return object_type<T>(m_data.clone(), from_object());
}

template<class T>
force_inline
bool object_type<T>::is_unique() const
{
    return m_data.is_unique();
}

template<class T>
force_inline
void object_type<T>::make_unique()
{
    m_data.make_unique();
}

template<class T>
force_inline
Type object_type<T>::get_type() const
{
    return m_data.get_type();
}

template<class T>
force_inline
Type object_type<T>::get_static_type()
{
    return dynamic::details::mark_type<T>().get();
};

template<class T>
force_inline
bool object_type<T>::is_null() const
{
    return m_data.is_null();
}

template<class T>
force_inline
bool object_type<T>::is_zero() const
{
    return m_data.is_zero();
}

template<class T>
force_inline
bool object_type<T>::is_one() const
{
    return m_data.is_one();
}

template<class T>
force_inline
bool object_type<T>::is_zero(const T& arg)
{
    return object_type_traits<T>::is_zero(arg);
}

template<class T>
force_inline
bool object_type<T>::is_one(const T& arg)
{
    return object_type_traits<T>::is_one(arg);
};

template<class T>
object_type<T> object_type<T>::make_one()
{
    return object_type(object_type_traits<T>::make_one((T*)nullptr));
};

template<class T>
force_inline
object_type<T>::operator bool() const
{
    return static_cast<bool>(this->get());
}

template<class T>
template<class S>
object_type<S> object_type<T>::cast() const
{
    return dynamic::cast<S>(*this);
}

template<class T>
template<class S>
object_type<S> object_type<T>::convert() const
{
    return dynamic::convert<S>(*this);
}

template<class T>
void object_type<T>::serialize(oarchive_impl & ar, const unsigned int version) const
{
    m_data.serialize(ar, version);
}

template<class T>
void object_type<T>::serialize(iarchive_impl & ar, const unsigned int version)
{
    m_data.serialize(ar, version);
};

template<class T>
std::ostream& operator<<(std::ostream& os, const object_type<T>& val)
{
    save_data(os, object(val));
    return os;
}

template<class T>
std::istream& operator>>(std::istream& is, object_type<T>& A)
{
    load_data(is, A.m_data);
    return is;
}

template<class T>
force_inline
void swap(object_type<T>& o1, object_type<T>& o2)
{
    swap(o1.m_data, o2.m_data);
};

};};
