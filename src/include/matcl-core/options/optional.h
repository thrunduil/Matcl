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

#include "matcl-core/config.h"

namespace matcl
{

// manages an optional contained value, i.e. a value that may or may not be
// present; to be replaced by std::optional when available;
// this is much simplified implementation assuming, that T is default 
// constructible
template<class T>
class optional
{
    public:
        // type of stored data
        using value_type    = T;

    private:
        T               m_value;
        bool            m_has_value;

    public:
        // constructs the object that does not contain a value
        optional();
        
        // constructs an optional object that contains a value
        optional(const T& value);
        optional(T&& value);

        // copy and move constructors
        optional(const optional& other);
        optional(optional&& other);

        // assignment from other optional value
        optional&       operator=(const optional& other);
        optional&       operator=(optional&& other);

        // assignment from a value
        optional&       operator=(const T& other);
        optional&       operator=(T&& other);

        // return the contained value; throw exception if this object does
        // not contain a value
        T&              value() &;
        const T&        value() const &;
        T&&             value() &&;
        const T&&       value() const &&;

        // return a pointer the contained value; throw exception if this object
        // does not contain a value
        const T*        operator->() const;
        T*              operator->();

        // return the contained value; throw exception if this object does
        // not contain a value
        const T&        operator*() const&;
        T&              operator*() &;
        const T&&       operator*() const&&;
        T&&             operator*() &&;

        // returns the contained value if this object has a value, otherwise
        // returns default_value
        T               value_or(const T& default_value ) const&;
        T               value_or(T&& default_value ) const&;
        T               value_or(const T& default_value ) &&;
        T               value_or(T&& default_value ) &&;

        // return true if this object contain a value
        bool            has_value() const;

        // if this object contais a value, destroy that value; otherwise, there
        // are no effects
        void            reset();

        // return true if this object contain a value
        explicit        operator bool() const;

    private:
        void            check_has_value() const;
};

}

#include "matcl-core/details/options/optional.inl"