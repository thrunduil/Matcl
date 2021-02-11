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
#include "matcl-mp/mp_float.h"

#include <complex>

namespace matcl
{

// unlimited precision floating point complex value type;
// real and imaginary part have the same precision
class MATCL_MP_EXPORT mp_complex
{
    private:
        mp_float        m_re;
        mp_float        m_im;

    //internal use
    public:
        struct impl_tag{};

        mp_complex(impl_tag);

    public:
        // default constructor; return zero value with default precision
        mp_complex();

        // default constructor; return zero value; prec gives precision in
        // bits; if prec = 0 then default precision is used
        explicit mp_complex(precision prec);

        // conversion from integer, single and double precision floating
        // point values; prec gives precision in bits; if prec = 0 then
        // default precision is used; result is rounded to nearest if 
        // rounding is required
        explicit mp_complex(Integer re, precision prec = precision(0));
        explicit mp_complex(Float re, precision prec = precision(0));
        explicit mp_complex(Real re, precision prec = precision(0));
        explicit mp_complex(long double re, precision prec = precision(0));

        // conversion from single and double precision floating point 
        // complex values; prec gives precision in bits; if prec = 0 then
        // default precision is used; result is rounded to nearest if 
        // rounding is required
        mp_complex(const Complex& val, precision prec = precision(0));
        mp_complex(const Float_complex& val, precision prec = precision(0));

        template<class Type>
        mp_complex(const std::complex<Type>& val, precision prec = precision(0));

        // conversion from variable precision floating point values
        mp_complex(const mp_float& re);
        mp_complex(mp_float&& re);

        // conversion from variable precision floating point values; prec gives
        // precision in bits; if prec = 0 then default precision is used
        mp_complex(const mp_float& re, precision prec);
        mp_complex(mp_float&& re, precision prec);

        // explicit conversions; prec gives precision in bits; if prec = 0 then
        // default precision is used; result is rounded to nearest if 
        // rounding is required
        explicit mp_complex(const mp_int& val, precision prec = precision(0));
        explicit mp_complex(const mp_rational& val, precision prec = precision(0));

        // construct complex value from real and imaginary part; prec gives
        // precision in bits; if prec = 0 then default precision is used; 
        // result is rounded to nearest if rounding is required
        mp_complex(Float re, Float im, precision prec = precision(0));
        mp_complex(Real re, Real im, precision prec = precision(0));
        mp_complex(long double re, long double im, precision prec = precision(0));

        // construct complex value from real and imaginary part; precision is set
        // to maximum precision of re and im
        mp_complex(const mp_float& re, const mp_float& im);
        mp_complex(const mp_float& re, mp_float&& im);
        mp_complex(mp_float&& re, const mp_float& im);
        mp_complex(mp_float&& re, mp_float&& im);

        // construct complex value from real and imaginary part; prec gives
        // precision in bits; if prec = 0 then default precision is used
        mp_complex(const mp_float& re, const mp_float& im, precision prec);
        mp_complex(const mp_float& re, mp_float&& im, precision prec);
        mp_complex(mp_float&& re, const mp_float& im, precision prec);
        mp_complex(mp_float&& re, mp_float&& im, precision prec);

        // standard copy and move constructors
        mp_complex(const mp_complex& other);
        mp_complex(mp_complex&& other);

        // standard copy and move constructors, change precision to prec; 
        // if prec = 0 then default precision is used
        mp_complex(const mp_complex& other, precision prec);
        mp_complex(mp_complex&& other, precision prec);

        // standard destructor
        ~mp_complex();

        // assignment; set imaginary part to zero
        mp_complex&     operator=(Integer other) &;
        mp_complex&     operator=(Float other) &;
        mp_complex&     operator=(Real other) &;
        mp_complex&     operator=(const mp_float& other) &;
        mp_complex&     operator=(mp_float&& other) &;

        // assignment of complex values
        mp_complex&     operator=(const Complex& other) &;
        mp_complex&     operator=(const Float_complex& other) &;
        mp_complex&     operator=(const mp_complex& other) &;
        mp_complex&     operator=(mp_complex&& other) &;

        // return real and imaginary part
        const mp_float& real() const    { return m_re; };
        const mp_float& imag() const    { return m_im; };

        // set real part; precision of imaginary part is changed to
        // match precision of real part
        void            set_real(const mp_float& val);
        void            set_real(mp_float&& val);

        // set imaginary part; precision of real part is changed to
        // match precision of imaginary part
        void            set_imag(const mp_float& val);
        void            set_imag(mp_float&& val);

        // convert to boolean value, return true if this value is nonzero
        explicit        operator bool() const;

        // convert to Complex value; round to nearest
        Complex         cast_complex() const;

        // cast to mp_int; imaginary part is ignored; on default round to zero
        mp_int          cast_mp_int(round_mode rm = round_mode::trunc) const;

        // cast to mp_rational; imaginary part is ignored; result is exact
        mp_rational     cast_mp_rational() const;

        // cast to mp_float; imaginary part is ignored; result is exact
        mp_float        cast_mp_float() const;                

        // return string representation of this value with n significant 
        // digits; if n = 0 then convert with the maximum available digits
        std::string     to_string(precision n = precision(0)) const;

        // return string representation of this value; real and imaginary
        // part is converted to string according to format string
        // see MPFR documentation for format description
        std::string     to_string(const std::string& format) const;

        // return precision for this object
        precision       get_precision() const;

        // set precision for this object
        void            set_precision(precision prec);

        // default precision can be changed using functions from mp_float class

        // serialization and deserialization
        void            serialize(oarchive_impl & ar, unsigned int version) const;
        void            serialize(iarchive_impl & ar, unsigned int version);
};

// save and load to stream
MATCL_MP_EXPORT std::istream&   operator>>(std::istream&, mp_complex& val);
MATCL_MP_EXPORT std::ostream&   operator<<(std::ostream&, const mp_complex& val);

template<class Type>
inline mp_complex::mp_complex(const std::complex<Type>& val, precision prec)
    :mp_complex(std::real(val), std::imag(val), prec)
{};

};
