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

#include "matcl-core/options/optional.h"

namespace matcl
{

template<class T>
inline optional<T>::optional()
    :m_has_value(false), m_value(T())
{};
        
template<class T>
inline optional<T>::optional(const T& value)
    :m_value(value), m_has_value(true)
{};

template<class T>
inline optional<T>::optional(T&& value)
    :m_value(std::move(value)), m_has_value(true)
{};

template<class T>
inline optional<T>::optional(const optional& other)
    :m_value(other.m_value), m_has_value(other.m_has_value)
{};

template<class T>
inline optional<T>::optional(optional&& other)
    :m_value(other.m_value), m_has_value(std::move(other.m_has_value))
{};


template<class T>
inline optional<T>& optional<T>::operator=(const optional<T>& other)
{
    m_value         = other.m_value;
    m_has_value     = other.m_has_value;
    return *this;
}

template<class T>
inline optional<T>& optional<T>::operator=(optional<T>&& other)
{
    m_value         = std::move(other.m_value);
    m_has_value     = other.m_has_value;
    return *this;
}

template<class T>
inline optional<T>& optional<T>::operator=(const T& other)
{
    m_value     = other;
    m_has_value = true;
    return *this;
}

template<class T>
inline optional<T>& optional<T>::operator=(T&& other)
{
    m_value     = std::move(other);
    m_has_value = true;
    return *this;
}

template<class T>
inline T& optional<T>::value() &
{
    check_has_value();
    return m_value;
}

template<class T>
inline const T& optional<T>::value() const &
{
    check_has_value();
    return m_value;
}

template<class T>
inline T&& optional<T>::value() &&
{
    check_has_value();
    return std::move(m_value);
}

template<class T>
inline const T&& optional<T>::value() const &&
{
    check_has_value();
    return std::move(m_value);
}

template<class T>
inline const T* optional<T>::operator->() const
{
    check_has_value();
    return &m_value;
}

template<class T>
inline T* optional<T>::operator->()
{
    check_has_value();
    return &m_value;
}

template<class T>
inline const T& optional<T>::operator*() const&
{
    check_has_value();
    return m_value;
}

template<class T>
inline T& optional<T>::operator*() &
{
    check_has_value();
    return m_value;
}

template<class T>
inline const T&& optional<T>::operator*() const&&
{
    check_has_value();
    return std::move(m_value);
};

template<class T>
inline T&& optional<T>::operator*() &&
{
    check_has_value();
    return std::move(m_value);
};

template<class T>
inline T optional<T>::value_or(const T& default_value ) const&
{
    if (has_value())
        return m_value;
    else
        return default_value;
}

template<class T>
inline T optional<T>::value_or(T&& default_value ) const&
{
    T tmp(std::move(default_value));

    if (has_value())
        return m_value;
    else
        return std::move(tmp);
}

template<class T>
inline T optional<T>::value_or(const T& default_value ) &&
{
    if (has_value())
        return std::move(m_value);
    else
        return default_value;
}

template<class T>
inline T optional<T>::value_or(T&& default_value ) &&
{
    T tmp(std::move(default_value));

    if (has_value())
        return std::move(m_value);
    else
        return std::move(tmp);
}

template<class T>
inline bool optional<T>::has_value() const
{
    return m_has_value;
}

template<class T>
inline optional<T>::operator bool() const
{
    return m_has_value;
};

template<class T>
void optional<T>::reset()
{
    if (has_value() == false)
        return;

    m_value = T();
};

template<class T>
void optional<T>::check_has_value() const
{
    if (m_has_value == false)
        throw std::runtime_error("optional does not contain a value");
};

}
