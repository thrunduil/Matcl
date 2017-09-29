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
#include "matcl-mp/config.h"
#include "matcl-mp/details/fwd_decls.h"
#include "matcl-core/general/fwd_decls.h"
#include "matcl-core/matrix/enums.h"
#include "matcl-mp/mp_float.h"

namespace matcl
{

// unlimited precision integer type
class MATCL_MP_EXPORT mp_int
{
    private:
        static const 
        size_t Size     = sizeof(details::impl_mp_int);

    //internal use
    public:
        char            m_data[Size];

        friend class mp_rational;

        struct impl_tag{};

        template<class T>
        mp_int(impl_tag, const T&);

        mp_int(impl_tag);

    public:
        // default constructor; return zero value
        mp_int();

        // conversion from integer
        mp_int(Integer val);

        // standard copy and move constructors
        mp_int(const mp_int& other);
        mp_int(mp_int&& other);

        // convert string to integer value
        explicit mp_int(const std::string& s, int base = 0);

        // standard destructor
        ~mp_int();

        // assignment
        mp_int&         operator=(Integer other) &;
        mp_int&         operator=(const mp_int& other) &;
        mp_int&         operator=(mp_int&& other) &;

        mp_int&         operator=(Float other) & = delete;
        mp_int&         operator=(Real other) & = delete;

        // convert to boolean value, return true if this value is nonzero
        explicit        operator bool() const;

        // convert to integer value; if this value is too big to be 
        // represented by Integer, then zero is returned
        Integer         cast_int() const;

        // convert to Real value
        Real            cast_float() const;

        // convert to rational type
        mp_rational     cast_mp_rational() const;

        // convert to floating point types with precision p; if p is zero
        // then default precision is used; result is rounded to nearest if 
        // rounding is required
        mp_float        cast_mp_float(precision p = precision()) const;        
        mp_complex      cast_mp_complex(precision p = precision()) const;

        // return true if this value can be casted to Integer
        bool            can_cast_int() const;

        // return string representation of this value
        std::string     to_string() const;

        // serialization and deserialization
        void            serialize(oarchive_impl & ar, unsigned int version) const;
        void            serialize(iarchive_impl & ar, unsigned int version);
};

// save and load to stream
MATCL_MP_EXPORT std::istream&   operator>>(std::istream&, mp_int& val);
MATCL_MP_EXPORT std::ostream&   operator<<(std::ostream&, const mp_int& val);

};
