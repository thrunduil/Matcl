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

#include "matcl-mp/mp_complex.h"
#include "utils/impl_types.h"
#include "matcl-core/matrix/complex_type.h"
#include "matcl-core/IO/archive.h"
#include "matcl-mp/func_unary.h"
#include "utils/utils.h"

namespace matcl
{

//-----------------------------------------------------------------------------
//                                  mp_complex
//-----------------------------------------------------------------------------
mp_complex::mp_complex(impl_tag)
    :m_re(mp_float::impl_tag{}), m_im(mp_float::impl_tag{})
{};

mp_complex::mp_complex()
{};
mp_complex::mp_complex(precision prec)
    :m_re(prec), m_im(prec)
{};

mp_complex::mp_complex(long double re, precision prec)
    :m_re(re, prec), m_im(0.0, prec)
{};
mp_complex::mp_complex(Real re, precision prec)
    :m_re(re, prec), m_im(0.0, prec)
{};
mp_complex::mp_complex(Float re, precision prec)
    :m_re(re, prec), m_im(0.0, prec)
{};
mp_complex::mp_complex(Integer re, precision prec)
    :m_re(re, prec), m_im(0.0, prec)
{};

mp_complex::mp_complex(const Complex& val, precision prec)
    :m_re(matcl::real(val), prec), m_im(matcl::imag(val), prec)
{};
mp_complex::mp_complex(const Float_complex& val, precision prec)
    :m_re(matcl::real(val), prec), m_im(matcl::imag(val), prec)
{};

mp_complex::mp_complex(const mp_float& re)
    :m_re(re), m_im(0.0, re.get_precision())
{};

mp_complex::mp_complex(mp_float&& re)
    :m_re(std::move(re)), m_im(0.0, m_re.get_precision())
{};

mp_complex::mp_complex(const mp_float& re, precision prec)
    :m_re(re, prec), m_im(0.0, prec)
{};

mp_complex::mp_complex(mp_float&& re, precision prec)
    :m_re(re, prec), m_im(0.0, prec)
{};

mp_complex::mp_complex(const mp_int& val, precision prec)
    :mp_complex(mp_float(val, prec), prec)
{};
mp_complex::mp_complex(const mp_rational& val, precision prec)
    :mp_complex(mp_float(val, prec), prec)
{};

mp_complex::mp_complex(long double re, long double im, precision prec)
    :m_re(re, prec), m_im(im, prec)
{};

mp_complex::mp_complex(Real re, Real im, precision prec)
    :m_re(re, prec), m_im(im, prec)
{};
mp_complex::mp_complex(Float re, Float im, precision prec)
    :m_re(re, prec), m_im(im, prec)
{};

mp_complex::mp_complex(const mp_float& re, const mp_float& im)
    :m_re(re, details::max_precision(re, im))
    ,m_im(im, m_re.get_precision())
{};

mp_complex::mp_complex(mp_float&& re, const mp_float& im)
    :m_re(std::move(re), details::max_precision(re, im))
    ,m_im(im, m_re.get_precision())
{};

mp_complex::mp_complex(const mp_float& re, mp_float&& im)
    :m_re(re, details::max_precision(re, im))
    ,m_im(std::move(im), m_re.get_precision())
{};

mp_complex::mp_complex(mp_float&& re, mp_float&& im)
    :m_re(std::move(re), details::max_precision(re, im))
    ,m_im(std::move(im), m_re.get_precision())
{};

mp_complex::mp_complex(const mp_float& re, const mp_float& im, precision prec)
    :m_re(re, prec), m_im(im, prec)
{};

mp_complex::mp_complex(mp_float&& re, const mp_float& im, precision prec)
    :m_re(std::move(re), prec), m_im(im, prec)
{};
mp_complex::mp_complex(const mp_float& re, mp_float&& im, precision prec)
    :m_re(re, prec), m_im(std::move(im), prec)
{};
mp_complex::mp_complex(mp_float&& re, mp_float&& im, precision prec)
    :m_re(std::move(re), prec), m_im(std::move(im), prec)
{};

mp_complex::mp_complex(const mp_complex& other)
    :m_re(other.m_re), m_im(other.m_im)
{};
mp_complex::mp_complex(mp_complex&& other)
    :m_re(std::move(other.m_re)), m_im(std::move(other.m_im))
{};

mp_complex::mp_complex(const mp_complex& other, precision prec)
    :m_re(other.m_re, prec), m_im(other.m_im, prec)
{};
mp_complex::mp_complex(mp_complex&& other, precision prec)
    :m_re(std::move(other.m_re), prec), m_im(std::move(other.m_im), prec)
{};

mp_complex::~mp_complex()
{};

mp_complex& mp_complex::operator=(Integer other) &
{
    this->m_re = mp_float(other);
    this->m_im = mp_float(0.0);
    return *this;
};
mp_complex& mp_complex::operator=(Float other) &
{
    this->m_re = mp_float(other);
    this->m_im = mp_float(0.0);
    return *this;
};
mp_complex& mp_complex::operator=(Real other) &
{
    this->m_re = mp_float(other);
    this->m_im = mp_float(0.0);
    return *this;
};
mp_complex& mp_complex::operator=(const mp_float& other) &
{
    this->m_re = mp_float(other);
    this->m_im = mp_float(0.0, other.get_precision());
    return *this;
};
mp_complex& mp_complex::operator=(mp_float&& other) &
{
    this->m_re = mp_float(std::move(other));
    this->m_im = mp_float(0.0, this->m_re.get_precision());
    return *this;
};

mp_complex& mp_complex::operator=(const Complex& other) &
{
    this->m_re = mp_float(matcl::real(other));
    this->m_im = mp_float(matcl::imag(other));
    return *this;
};
mp_complex& mp_complex::operator=(const Float_complex& other) &
{
    this->m_re = mp_float(matcl::real(other));
    this->m_im = mp_float(matcl::imag(other));
    return *this;
};
mp_complex& mp_complex::operator=(const mp_complex& other) &
{
    this->m_re = other.m_re;
    this->m_im = other.m_im;
    return *this;
};
mp_complex& mp_complex::operator=(mp_complex&& other) &
{
    this->m_re = std::move(other.m_re);
    this->m_im = std::move(other.m_im);
    return *this;
};
void mp_complex::set_real(const mp_float& val)
{
    this->m_re = val;
    this->m_im.set_precision(val.get_precision());
};
void mp_complex::set_real(mp_float&& val)
{
    this->m_re = std::move(val);
    this->m_im.set_precision(this->m_re.get_precision());
};

void mp_complex::set_imag(const mp_float& val)
{
    this->m_re.set_precision(val.get_precision());
    this->m_im = val;    
}
void mp_complex::set_imag(mp_float&& val)
{
    this->m_im = std::move(val);    
    this->m_re.set_precision(this->m_im.get_precision());    
}

mp_complex::operator bool() const
{
    return (bool)m_re || (bool)m_im;
};

Complex mp_complex::cast_complex() const
{
    return Complex(m_re.cast_float(), m_im.cast_float());
}

mp_int mp_complex::cast_mp_int(round_mode rm) const
{
    return m_re.cast_mp_int(rm);
}
mp_rational mp_complex::cast_mp_rational() const
{
    return m_re.cast_mp_rational();
}
mp_float mp_complex::cast_mp_float() const
{
    return m_re;
}

precision mp_complex::get_precision() const
{
    return m_re.get_precision();
}

void mp_complex::set_precision(precision prec)
{
    m_re.set_precision(prec);
    m_im.set_precision(prec);
};

std::string mp_complex::to_string(precision n) const
{
    size_t digits   = (n > 0) ? n : 1 + mp_float::bits_to_digits10(this->get_precision());

    std::ostringstream format;    
    
    format << "%." << digits << "RNg";

    return to_string(format.str());
};

std::string mp_complex::to_string(const std::string& format) const
{
    if (format.empty() == true)
        return std::string();

    std::ostringstream out;
    out << m_re.to_string(format);
    out << "+";
    out << m_im.to_string(format);
    out << "i";

    return out.str();
};

void mp_complex::serialize(oarchive_impl & ar, unsigned int version) const
{
    m_re.serialize(ar, version);
    m_im.serialize(ar, version);
};

void mp_complex::serialize(iarchive_impl & ar, unsigned int version)
{
    m_re.serialize(ar, version);
    m_im.serialize(ar, version);
};

//functions

std::istream& matcl::operator>>(std::istream& is, mp_complex& val)
{
    mp_float re, im;
    is >> re;
    is >> im;

    val = mp_complex(std::move(re), std::move(im));
    return is;
}
std::ostream& matcl::operator<<(std::ostream& os, const mp_complex& val)
{
    os << val.real();
    os << " ";
    os << val.imag();
    return os;
}

};
