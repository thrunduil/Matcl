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

#include "matcl-mp/mp_rational.h"
#include "utils/impl_types.h"
#include "matcl-core/matrix/complex_type.h"
#include "matcl-core/IO/archive.h"
#include "matcl-mp/func_unary.h"
#include "matcl-mp/error_flags.h"

namespace matcl
{

namespace mmd = matcl::mp::details;

//-----------------------------------------------------------------------------
//                                  mp_rational
//-----------------------------------------------------------------------------
mp_rational::mp_rational(impl_tag)
{};

template<class T>
mp_rational::mp_rational(impl_tag, const T& val)
{
    using impl_type = mmd::impl_rat;
    new (mmd::get_ptr<mp_rational>::eval(m_data)) impl_type(val);
};

mp_rational::mp_rational()
{
    using impl_type = mmd::impl_rat;
    new (mmd::get_ptr<mp_rational>::eval(m_data)) impl_type();
};

mp_rational::mp_rational(Integer val)
{
    using impl_type = mmd::impl_rat;
    new (mmd::get_ptr<mp_rational>::eval(m_data)) mmd::impl_rat(val);
};

mp_rational::mp_rational(Float val)
{
    using impl_type = mmd::impl_rat;
    if (std::isfinite(val) == true)
    {
        new (mmd::get_ptr<mp_rational>::eval(m_data)) mmd::impl_rat(val);
    }
    else
    {
        error_flags::set_integer_overflow();
        new (mmd::get_ptr<mp_rational>::eval(m_data)) mmd::impl_rat(Integer(0));
    }
};
mp_rational::mp_rational(Real val)
{
    using impl_type = mmd::impl_rat;
    if (std::isfinite(val) == true)
    {
        new (mmd::get_ptr<mp_rational>::eval(m_data)) mmd::impl_rat(val);
    }
    else
    {
        error_flags::set_integer_overflow();
        new (mmd::get_ptr<mp_rational>::eval(m_data)) mmd::impl_rat(Integer(0));
    }
};
mp_rational::mp_rational(const mp_float& val)
{
    using impl_type     = mmd::impl_rat;

    if (is_finite(val) == false)
    {
        error_flags::set_integer_overflow();
        new (mmd::get_ptr<mp_rational>::eval(m_data)) mmd::impl_rat(Integer(0));
    }
    else
    {        
        using impl_float    = mmd::impl_float;
        const impl_float& v = *mmd::get_ptr<mp_float>::eval(val.m_data);
        impl_type& iv       = *mmd::get_ptr<mp_rational>::eval(this->m_data);

        //convert to GMP floating point type; this conversion should be exact
        mpf_t tmp;
        mpf_init2(tmp, (long)val.get_precision());
        mpfr_get_f(tmp, v, MPFR_RNDN);

        //convert to rational

        // we cannot use iv.backend().data()
        mpq_init(*reinterpret_cast<mpq_t*>(&iv.backend()));
        mpq_set_f(iv.backend().data(), tmp);

        mpf_clear(tmp);
    };
};
mp_rational::mp_rational(const mp_int& num)
{
    using impl_type     = mmd::impl_rat;
    using int_type      = mmd::impl_int;
    const int_type& v1  = *mmd::get_ptr<mp_int>::eval(num.m_data);
    new (mmd::get_ptr<mp_rational>::eval(m_data)) impl_type(v1);
};

mp_rational::mp_rational(Integer num, Integer den)
{
    using impl_type = mmd::impl_rat;
    
    if (den < 0)
    {
        den = -den;
        num = -num;
    };

    if (den == 0)
    {
        error_flags::set_integer_overflow();
        new (mmd::get_ptr<mp_rational>::eval(m_data)) mmd::impl_rat(0);
    }
    else
    {
        new (mmd::get_ptr<mp_rational>::eval(m_data)) mmd::impl_rat(num,den);
    }
};

mp_rational::mp_rational(const mp_int& num, const mp_int& den)
{
    using impl_type     = mmd::impl_rat;
    using int_type      = mmd::impl_int;
    const int_type& v1  = *mmd::get_ptr<mp_int>::eval(num.m_data);
    const int_type& v2  = *mmd::get_ptr<mp_int>::eval(den.m_data);

    if (is_zero(den) == true)
    {
        error_flags::set_integer_overflow();
        new (mmd::get_ptr<mp_rational>::eval(m_data)) impl_type(0);
    }
    else
    {
        new (mmd::get_ptr<mp_rational>::eval(m_data)) impl_type(v1,v2);
    }
};

mp_rational::mp_rational(const std::string& s, int base)
{
    using impl_type = mmd::impl_rat;
    new (mmd::get_ptr<mp_rational>::eval(m_data)) mmd::gmp_impl_rat(s,base);

    using impl_type     = mmd::impl_rat;
    impl_type& v1       = *mmd::get_ptr<mp_rational>::eval(this->m_data);

    mpq_canonicalize(v1.backend().data());
};

mp_rational::mp_rational(const mp_rational& other)
{
    using impl_type         = mmd::impl_rat;
    const impl_type& val    = *mmd::get_ptr<mp_rational>::eval(other.m_data);
    new (mmd::get_ptr<mp_rational>::eval(m_data)) impl_type(val);
};

mp_rational::mp_rational(mp_rational&& other)
{
    using impl_type = mmd::impl_rat;
    impl_type&& val = std::move(*mmd::get_ptr<mp_rational>::eval(other.m_data));
    new (mmd::get_ptr<mp_rational>::eval(m_data)) impl_type(std::move(val));
};

mp_rational::~mp_rational()
{
    using impl_type = mmd::impl_rat;
    impl_type& val  = *mmd::get_ptr<mp_rational>::eval(this->m_data);
    val.~impl_type();
};

mp_rational& mp_rational::operator=(Integer other) &
{
    using impl_type         = mmd::impl_rat;
    impl_type& this_val     = *mmd::get_ptr<mp_rational>::eval(this->m_data);

    this_val.operator=(other);
    return *this;
};

mp_rational& mp_rational::operator=(const mp_int& other) &
{
    using impl_type         = mmd::impl_rat;
    using impl_int          = mmd::impl_int;
    const impl_int& val     = *mmd::get_ptr<mp_int>::eval(other.m_data);
    impl_type& this_val     = *mmd::get_ptr<mp_rational>::eval(this->m_data);

    this_val.operator=(val);
    return *this;
};

mp_rational& mp_rational::operator=(mp_int&& other) &
{
    using impl_type         = mmd::impl_rat;
    using impl_int          = mmd::impl_int;
    impl_type& this_val     = *mmd::get_ptr<mp_rational>::eval(this->m_data);
    impl_int&& val          = std::move(*mmd::get_ptr<mp_int>::eval(other.m_data));

    this_val.operator=(std::move(val));
    return *this;
};

mp_rational& mp_rational::operator=(const mp_rational& other) &
{
    using impl_type         = mmd::impl_rat;
    const impl_type& val    = *mmd::get_ptr<mp_rational>::eval(other.m_data);
    impl_type& this_val     = *mmd::get_ptr<mp_rational>::eval(this->m_data);

    this_val.operator=(val);
    return *this;
};

mp_rational& mp_rational::operator=(mp_rational&& other) &
{
    using impl_type         = mmd::impl_rat;
    impl_type&& val         = std::move(*mmd::get_ptr<mp_rational>::eval(other.m_data));
    impl_type& this_val     = *mmd::get_ptr<mp_rational>::eval(this->m_data);

    this_val.operator=(std::move(val));
    return *this;
};

mp_rational::operator bool() const
{
    using impl_type         = mmd::impl_rat;
    const impl_type& val    = *mmd::get_ptr<mp_rational>::eval(this->m_data);
    return (bool)val;
};

Real mp_rational::cast_float() const
{
    //do not use mpq_get_d due to different rounding mode
    return this->cast_mp_float(precision::precision_double()).cast_float();
}
Integer mp_rational::cast_int() const
{
    return this->cast_mp_int().cast_int();
}

mp_int mp_rational::cast_mp_int() const
{
    mp_int ret(mp_int::impl_tag{});

    using impl_type     = mmd::impl_int;
    using impl_rat      = mmd::impl_rat;
    
    const impl_rat& v   = *mmd::get_ptr<mp_rational>::eval(this->m_data);
    
    new (mmd::get_ptr<mp_int>::eval(ret.m_data)) impl_type(v);
    return ret;
}
mp_float mp_rational::cast_mp_float(precision p) const
{
    return mp_float(*this, p);
}
mp_complex mp_rational::cast_mp_complex(precision p) const
{
    return mp_complex(*this, p);
}

mp_int mp_rational::numerator() const
{
    using impl_type         = mmd::impl_rat;
    const impl_type& val    = *mmd::get_ptr<mp_rational>::eval(this->m_data);

    return mp_int(mp_int::impl_tag(), boost::multiprecision::numerator(val));
};

mp_int mp_rational::denominator() const
{
    using impl_type         = mmd::impl_rat;
    const impl_type& val    = *mmd::get_ptr<mp_rational>::eval(this->m_data);

    return mp_int(mp_int::impl_tag(), boost::multiprecision::denominator(val));
};

std::string mp_rational::to_string() const
{
    using impl_type         = mmd::impl_rat;
    const impl_type& v1     = *mmd::get_ptr<mp_rational>::eval(this->m_data);
    return v1.str();
};

void mp_rational::serialize(iarchive_impl & ar, unsigned int version)
{
    mp_int num, den;
    num.serialize(ar, version);
    den.serialize(ar, version);

    *this = mp_rational(num, den);
};

void mp_rational::serialize(oarchive_impl & ar, unsigned int version) const
{
    this->numerator().serialize(ar, version);
    this->denominator().serialize(ar, version);
};

std::istream& matcl::operator>>(std::istream& is, mp_rational& val)
{
    using impl_type         = mmd::impl_rat;
    impl_type& v1           = *mmd::get_ptr<mp_rational>::eval(val.m_data);
    is >> v1;
    return is;
}

std::ostream& matcl::operator<<(std::ostream& os, const mp_rational& val)
{
    using impl_type         = mmd::impl_rat;
    const impl_type& v1     = *mmd::get_ptr<mp_rational>::eval(val.m_data);
    os << v1;
    return os;
};

template mp_rational::mp_rational(impl_tag, const mmd::impl_rat&);
};
