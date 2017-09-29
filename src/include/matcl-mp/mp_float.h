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

#pragma warning(push)
#pragma warning(disable: 4251) //needs to have dll-interface to be used by clients

namespace matcl
{

// strongly typed size_t type representing precision 
// of mp_float variables
class MATCL_MP_EXPORT precision
{
    private:
        size_t value;

    public:
        // undefined precision
        inline precision()              : value(0){};

        // precision in bits; precision can be adjusted
        // if is too small or too large
        explicit precision(size_t val);

        // get precision in bits
        inline size_t   get() const     { return value; };

        // convert this object to size_t type
        inline operator size_t() const  { return value; };

        // return precision of long double type in bits
        static precision precision_long_double();

        // return precision of double type in bits
        static precision precision_double();

        // return precision of float type in bits
        static precision precision_float();

        // basic arithmetic operators
        template<class T>
        precision       operator+(T p) const    { return precision(value + p); };
        template<class T>
        precision       operator-(T p) const    { return precision(value - p); };
        template<class T>
        precision       operator*(T p) const    { return precision(value * p); };
        template<class T>
        precision       operator/(T p) const    { return precision(value / p); };

        template<class T>
        precision&      operator+=(T p)         { value += p; return *this; };
        template<class T>
        precision&      operator-=(T p)         { value -= p; return *this; };
        template<class T>
        precision&      operator*=(T p)         { value *= p; return *this; };
        template<class T>
        precision&      operator/=(T p)         { value /= p; return *this; };
};

enum class round_mode
{
    floor,      // round toward -INF
    ceil,       // round toward +INF
    trunc,      // round toward 0
    round       // round to nearest with halfway cases rounded away from zero
};

// unlimited precision floating point value type
class MATCL_MP_EXPORT mp_float
{
    private:
        static const 
        size_t Size     = sizeof(details::impl_mp_float);

    //internal use
    public:
        char            m_data[Size];

        #if MATCL_DEBUG_MP_FLOAT
          std::string   m_debug;
          size_t        m_prec;
        #endif

        struct impl_tag{};

        mp_float(impl_tag);

    public:
        // default constructor; return zero value with default precision
        mp_float();

        // default constructor; return zero value; prec gives precision in
        // bits; if prec = 0 then default precision is used
        explicit mp_float(precision prec);

        // conversion from integer, single and double precision floating
        // point values; prec gives precision in bits; if prec = 0 then
        // default precision is used; result is rounded to nearest if 
        // rounding is required
        mp_float(Integer val, precision prec = precision(0));        
        mp_float(Float val, precision prec = precision(0));
        mp_float(Real val, precision prec = precision(0));
        mp_float(long double val, precision prec = precision(0));
        mp_float(const mp_int& val, precision prec = precision(0));
        mp_float(const mp_rational& val, precision prec = precision(0));

        // standard copy and move constructors
        mp_float(const mp_float& other);
        mp_float(mp_float&& other);

        // standard copy and move constructors, change precision to prec; 
        // if prec = 0 then default precision is used
        mp_float(const mp_float& other, precision prec);
        mp_float(mp_float&& other, precision prec);

        // convert string to floating point value with precision prec
        // if prec = 0, then default precision is used
        explicit mp_float(const std::string& s, precision prec = precision(0), 
                            int base = 0);

        // standard destructor
        ~mp_float();

        // assignment
        mp_float&       operator=(Integer other) &;        
        mp_float&       operator=(Float other) &;
        mp_float&       operator=(Real other) &;
        mp_float&       operator=(long double other) &;
        mp_float&       operator=(const mp_float& other) &;
        mp_float&       operator=(mp_float&& other) &;

        // set value to zero and precision to prec; if prec is zero
        // then default precision is used;
        mp_float&       clear(precision prec = precision(0));

        // convert to boolean value, return true if this value is nonzero
        explicit        operator bool() const;

        // convert to Real value; round to nearest
        Real            cast_float() const;

        // convert to Integer value; on default round to zero; if this value
        // is too big to be represented by Integer, or is not finite, then
        // zero is returned
        Integer         cast_int(round_mode rm = round_mode::trunc) const;

        // cast to mp_int; on default round to zero; if this value is not
        // finite then zero is returned
        mp_int          cast_mp_int(round_mode rm = round_mode::trunc) const;

        // cast to mp_rational; result is exact; if this value is not
        // finite then zero is returned
        mp_rational     cast_mp_rational() const;

        // cast to mp_complex; result is exact
        mp_complex      cast_mp_complex() const;

        // return true if this value can be casted to Integer using rm 
        // rounding mode
        bool            can_cast_int(round_mode rm = round_mode::trunc) const;

        // return string representation of this value with n significant 
        // digits; if n = 0 then convert with the maximum available digits
        std::string     to_string(precision n = precision(0)) const;

        // return string representation of this value using format string
        // see MPFR documentation for format description
        std::string     to_string(const std::string& format) const;

        // return precision for this object
        precision       get_precision() const;

        // set precision for this object; stored value can be rounded
        mp_float&       set_precision(precision prec);

        // return default precision for newly created objects; default precision
        // it thread local
        static precision get_default_precision();

        // set default precision for newly created objects; default precision
        // it thread local
        static void     set_default_precision(precision prec);

        // return minumum precision allowed
        static precision min_allowed_precision();

        // return maximum precision allowed
        static precision max_allowed_precision();

        // convert precision in decimal digits to bits
        static precision digits10_to_bits(precision prec);

        // convert precision in bits to decimal digits
        static precision bits_to_digits10(precision digits);

        // serialization and deserialization
        void            serialize(oarchive_impl & ar, unsigned int version) const;
        void            serialize(iarchive_impl & ar, unsigned int version);

    //internal use
    public:
        void            update_debug();
};

// save and load to stream
MATCL_MP_EXPORT std::istream&   operator>>(std::istream&, mp_float& val);
MATCL_MP_EXPORT std::ostream&   operator<<(std::ostream&, const mp_float& val);

#if MATCL_DEBUG_MP_FLOAT != 1
    inline void mp_float::update_debug(){};
#endif

};

#pragma warning(pop)