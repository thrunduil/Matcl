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

#include "matcl-mp/mp_int.h"
#include "utils/impl_types.h"
#include "matcl-core/matrix/complex_type.h"
#include "matcl-core/IO/archive.h"

namespace matcl
{

namespace mmd = matcl::mp::details;

//-----------------------------------------------------------------------------
//                                  mp_int
//-----------------------------------------------------------------------------
mp_int::mp_int()
{
    using impl_type = mmd::impl_int;
    new (mmd::get_ptr<mp_int>::eval(m_data)) impl_type();
};

template<class T>
mp_int::mp_int(impl_tag, const T& val)
{
    using impl_type = mmd::impl_int;
    new (mmd::get_ptr<mp_int>::eval(m_data)) impl_type(val);
};

mp_int::mp_int(impl_tag)
{
};

mp_int::mp_int(Integer val)
{
    using impl_type = mmd::impl_int;
    new (mmd::get_ptr<mp_int>::eval(m_data)) impl_type(val);
};

mp_int::mp_int(const std::string& s, int base)
{
    using impl_type = mmd::impl_int;
    new (mmd::get_ptr<mp_int>::eval(m_data)) mmd::gmp_impl_int(s,base);
};

mp_int::mp_int(const mp_int& other)
{
    using impl_type         = mmd::impl_int;
    const impl_type& val    = *mmd::get_ptr<mp_int>::eval(other.m_data);
    new (mmd::get_ptr<mp_int>::eval(m_data)) impl_type(val);
};
mp_int::mp_int(mp_int&& other)
{
    using impl_type = mmd::impl_int;
    impl_type&& val = std::move(*mmd::get_ptr<mp_int>::eval(other.m_data));
    new (mmd::get_ptr<mp_int>::eval(m_data)) impl_type(std::move(val));
};
mp_int::~mp_int()
{
    using impl_type = mmd::impl_int;
    impl_type& val  = *mmd::get_ptr<mp_int>::eval(this->m_data);
    val.~impl_type();
};
mp_int& mp_int::operator=(Integer other) &
{
    using impl_type         = mmd::impl_int;
    impl_type& this_val     = *mmd::get_ptr<mp_int>::eval(this->m_data);

    this_val.operator=(other);
    return *this;
};

mp_int& mp_int::operator=(const mp_int& other) &
{
    using impl_type         = mmd::impl_int;
    const impl_type& val    = *mmd::get_ptr<mp_int>::eval(other.m_data);
    impl_type& this_val     = *mmd::get_ptr<mp_int>::eval(this->m_data);

    this_val.operator=(val);
    return *this;
};
mp_int& mp_int::operator=(mp_int&& other) &
{
    using impl_type         = mmd::impl_int;
    impl_type&& val         = std::move(*mmd::get_ptr<mp_int>::eval(other.m_data));
    impl_type& this_val     = *mmd::get_ptr<mp_int>::eval(this->m_data);

    this_val.operator=(std::move(val));
    return *this;
};

mp_int::operator bool() const
{
    using impl_type             = mmd::impl_int;
    const impl_type& this_val   = *mmd::get_ptr<mp_int>::eval(this->m_data);
    return (bool)this_val;
};

Integer mp_int::cast_int() const
{
    if (can_cast_int() == false)
        return 0;

    using impl_type             = mmd::impl_int;
    const impl_type& this_val   = *mmd::get_ptr<mp_int>::eval(this->m_data);
    return (Integer)this_val;
};

bool mp_int::can_cast_int() const
{
    return mpz_fits_slong_p(mmd::impl_value(*this).backend().data()) != 0;
};

Real mp_int::cast_float() const
{
    //do not use mpz_get_d due to different rounding mode
    return this->cast_mp_float(precision::precision_double()).cast_float();
}

mp_float mp_int::cast_mp_float(precision p) const
{
    return mp_float(*this, p);
}
mp_rational mp_int::cast_mp_rational() const
{
    return mp_rational(*this);
}
mp_complex mp_int::cast_mp_complex(precision p) const
{
    return mp_complex(*this, p);
}

std::string mp_int::to_string() const
{
    using impl_type     = mmd::impl_int;
    const impl_type& v1 = *mmd::get_ptr<mp_int>::eval(this->m_data);
    return v1.str();
};

void mp_int::serialize(oarchive_impl & ar, unsigned int version) const
{
    using impl_type     = mmd::impl_int;
    const impl_type& v1 = *mmd::get_ptr<mp_int>::eval(this->m_data);
    std::string str     = v1.str();

    (void)version;
    ar << str;
};

void mp_int::serialize(iarchive_impl & ar, unsigned int version)
{
    (void)version;

    std::string str;
    ar >> str;

    *this = mp_int(str);
};

std::istream& matcl::operator>>(std::istream& is, mp_int& val)
{
    using impl_type         = mmd::impl_int;
    impl_type& v1           = *mmd::get_ptr<mp_int>::eval(val.m_data);
    is >> v1;
    return is;
}
std::ostream& matcl::operator<<(std::ostream& os, const mp_int& val)
{
    using impl_type         = mmd::impl_int;
    const impl_type& v1     = *mmd::get_ptr<mp_int>::eval(val.m_data);
    os << v1;
    return os;
};

template mp_int::mp_int(impl_tag, const mmd::impl_int&);

};
