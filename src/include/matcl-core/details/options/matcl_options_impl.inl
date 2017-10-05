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

#include "matcl-core/memory/alloc.h"

#pragma warning(push)
#pragma warning(disable : 4702) //unreachable code

namespace matcl 
{

//-----------------------------------------------------------------------------------
//                              option
//-----------------------------------------------------------------------------------
template <class T>
optional<T> option::get() const
{
    using ptr_impl = details::option_impl_type<T>;

    const ptr_impl* impl = dynamic_cast<const ptr_impl*>(m_impl.get());

    if (impl != nullptr)
        return impl->get();

    error_invalid_option_type(m_impl->pretty_type_name(typeid(T).name()));
    return optional<T>();
};

//-----------------------------------------------------------------------------------
//                              option_type
//-----------------------------------------------------------------------------------
template<class T>
option_type<T>::option_type(const std::string& name, const std::string& descr,
                            const optional<T>& val)    
    :option(impl_type(new option_impl_type(name, descr, val)))
{}

template<class T>
option_type<T>::~option_type()
{};

//-----------------------------------------------------------------------------------
//                              string_option
//-----------------------------------------------------------------------------------

template<class T>   
string_option<T>::string_option(const std::string& name)
    :base_type(name, name, optional<T>())
{};

template<class T>
string_option<T>::string_option(const std::string& name, const optional<T>& value)
    :base_type(name, name, value)
{};
template<class T>
string_option<T>::string_option(const std::string& name, const T& value)
    :base_type(name, name, optional<T>(value))
{};

template<class T>
string_option<T>::string_option(const std::string& name, const optional<T>& value, const std::string& descr)
    :base_type(name, descr, value)
{};
template<class T>
string_option<T>::string_option(const std::string& name, const T& value, const std::string& descr)
    :base_type(name, descr, optional<T>(value))
{};

//-----------------------------------------------------------------------------------
//                              option_base
//-----------------------------------------------------------------------------------
template<class T, class Derived>
option_base<T,Derived>::option_base(const optional<T>& val)
    :base_type(this->get_name(), this->get_description(), 
            (this->get_validator())? this->get_validator()(val) : val)
{};

template<class T, class Derived>
std::string option_base<T,Derived>::m_name;

template<class T, class Derived>
std::string option_base<T,Derived>::m_description;

template<class T, class Derived>
bool option_base<T,Derived>::m_initialized = false;

template<class T, class Derived>
T option_base<T,Derived>::m_default_value = T();

template<class T, class Derived>
typename option_base<T,Derived>::mutex_type option_base<T,Derived>::m_mutex;

template<class T, class Derived>
typename option_base<T,Derived>::validator_type option_base<T,Derived>::m_validator;

template<class T, class Derived>
Derived option_base<T,Derived>::m_default_opt;

template<class T, class Derived>
option_base<T,Derived>::~option_base()
{};

template<class T, class Derived>
void option_base<T,Derived>::initialize()
{
    static_assert(sizeof(Derived) == sizeof(option), "option cannot have members");

    std::unique_lock<mutex_type> lock(m_mutex);

    if (m_initialized == false)
    {
        Derived::config();
        m_name  = make_name(typeid(Derived).name());        
        m_initialized = true;
        lock.unlock();

        option_base::m_default_opt = Derived(option_base<T,Derived>::m_default_value);
        details::options_impl::register_option(option_base::m_default_opt);

    };
};

template<class T, class Derived>
std::string option_base<T,Derived>::get_name()
{
    option_base<T,Derived>::initialize();
    return option_base<T,Derived>::m_name;
};

template<class T, class Derived>
std::string option_base<T,Derived>::get_description()
{
    option_base<T,Derived>::initialize();
    return option_base<T,Derived>::m_description;
};

template<class T, class Derived>
typename option_base<T,Derived>::validator_type option_base<T,Derived>::get_validator()
{
    option_base<T,Derived>::initialize();
    return option_base<T,Derived>::m_validator;
};;  

//-----------------------------------------------------------------------------------
//                              option_impl_type
//-----------------------------------------------------------------------------------
template<class T>
details::option_impl_type<T>::option_impl_type(const std::string& name, const std::string& descr,
                                               const impl_type& val)
    :m_name(name), m_value(val), m_description(descr)
{};

template<class T>
details::option_impl_type<T>::~option_impl_type()
{};

template<class T>
void details::option_impl_type<T>::accept(option_visitor& vis, const option& opt) const
{
    vis.visit(m_value, opt);
};

template<class T>
std::string details::option_impl_type<T>::get_option_type() const
{
    return pretty_type_name(typeid(T).name());
};

//--------------------------------------------------------------------------
//                              OPTIONS
//--------------------------------------------------------------------------
template <class T>
optional<T> options::get(const option& option_type) const
{
    return (*m_options).get<T>(option_type);
};

template <class T>
T options::get_option(const option& option_type) const
{
    optional<T> opt = this->get<T>(option_type);
    if (opt)
        return opt.value();
    else
        return get_default<T>(option_type);
};

template <class T>
T options::get_default(const option& option_type)
{
    optional<T> opt = default_options().get<T>(option_type);

    if (opt)
        return opt.value();

    return get_predefined<T>(option_type);    
};

template <class T>
T options::get_predefined(const option& option_type)
{
    T opt = details::options_impl::get_predefined<T>(option_type);
    return opt;
};

//--------------------------------------------------------------------------
//                              options_impl
//--------------------------------------------------------------------------
template<class T>
optional<T> details::options_impl::get(const option& option_value) const
{
    std::unique_lock<mutex_type> lock(m_mutex_local);

    std::string name = option_value.name();

    if (m_notifier)
        m_notifier->report(name);

    auto pos = m_options.find(name);

    optional<T> ret;

    if (pos != m_options.end())
        ret = pos->second.get<T>();

    if (!ret)
        ret = option_value.get<T>();

    return ret;
};

template <class T>
static T details::options_impl::get_predefined(const option& option_value)
{
    std::unique_lock<mutex_type> lock(get_mutex_global());

    std::string name = option_value.name();

    auto pos = get_options_predefined().find(name);

    if (pos == get_options_predefined().end())
    {
        error_unregistered_option(name);
        return T();
    };

    return pos->second.get<T>().value();
};

};

#pragma warning(pop)