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

#include "matcl-core/config.h"
#include "matcl-mp/config.h"
#include "matcl-mp/details/fwd_decls.h"
#include "matcl-core/general/fwd_decls.h"
#include "matcl-mp/details/initializer.h"
#include "matcl-core/matrix/enums.h"
#include "matcl-mp/mp_int.h"
#include "matcl-mp/mp_float.h"

namespace matcl
{

// unlimited precision rational value type
class MATCL_MP_EXPORT mp_rational
{
    private:
        static const 
        size_t Size     = sizeof(details::impl_mp_rational);

    //internal use
    public:
        char            m_data[Size];

        struct impl_tag{};

        template<class T>
        mp_rational(impl_tag, const T&);

        mp_rational(impl_tag);

    public:
        // default constructor; return zero value
        mp_rational();

        // convert integer value to rational value
        mp_rational(Integer num);
        mp_rational(const mp_int& num);        

        // explicit conversions; if floating point value is not
        // finite, then mp_rational is initialized to zero and
        // integer overflow flag is set; no rounding takes place
        explicit mp_rational(Float val);
        explicit mp_rational(Real val);
        explicit mp_rational(const mp_float& val);

        // create rational number num / den; if den is zero, then
        // integer overflow flag is set and result is zero
        mp_rational(Integer num, Integer den);
        mp_rational(const mp_int& num, const mp_int& den);

        // standard copy and move constructors
        mp_rational(const mp_rational& other);
        mp_rational(mp_rational&& other);

        // convert string to rational value
        explicit mp_rational(const std::string& s, int base = 0);

        // standard destructor
        ~mp_rational();

        // assignment
        mp_rational&    operator=(Integer other) &;
        mp_rational&    operator=(const mp_int& other) &;
        mp_rational&    operator=(mp_int&& other) &;
        mp_rational&    operator=(const mp_rational& other) &;
        mp_rational&    operator=(mp_rational&& other) &;

        mp_rational&    operator=(Float other) & = delete;
        mp_rational&    operator=(Real other) & = delete;

        // convert to boolean value, return true if this value is nonzero
        explicit        operator bool() const;

        // cast to Integer; on default round to zero if necessary; if this
        // value is too big to be represented by Integer, then zero is returned
        Integer         cast_int() const;

        // convert to Real value; round to nearest if necessary
        Real            cast_float() const;

        // cast to mp_int; on default round to zero if necessary
        mp_int          cast_mp_int() const;

        // cast to floating point value with given precision prec;
        // if prec = 0, then default precision is used; round to nearest
        mp_float        cast_mp_float(precision prec = precision(0)) const;        

        // cast to complex value with given precision prec; if prec = 0,
        // then default precision is used; round to nearest if necessary
        mp_complex      cast_mp_complex(precision prec = precision(0)) const;

        // return string representation of this value
        std::string     to_string() const;

        // return numerator
        mp_int          numerator() const;

        // return denominator
        mp_int          denominator() const;

        // serialization and deserialization
        void            serialize(oarchive_impl & ar, unsigned int version) const;
        void            serialize(iarchive_impl & ar, unsigned int version);
};

// save and load to stream
MATCL_MP_EXPORT std::istream&   operator>>(std::istream&, mp_rational& val);
MATCL_MP_EXPORT std::ostream&   operator<<(std::ostream&, const mp_rational& val);

};
