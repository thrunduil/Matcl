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

#include "matcl-mp/mp_float.h"
#include "utils/impl_types.h"
#include "matcl-core/matrix/complex_type.h"
#include "matcl-core/IO/archive.h"
#include "utils/utils.h"

namespace matcl
{

namespace mmd = matcl::mp::details;

precision::precision(size_t val)
{
    value = std::min<size_t>(val, (size_t)MPFR_PREC_MAX - (size_t)20);
};

precision precision::precision_long_double()
{
    return precision(std::numeric_limits<long double>::digits);
}
precision precision::precision_double()
{
    return precision(std::numeric_limits<double>::digits);
}
precision precision::precision_float()
{
    return precision(std::numeric_limits<float>::digits);
}

//-----------------------------------------------------------------------------
//                                  mp_float
//-----------------------------------------------------------------------------
mp_float::mp_float(impl_tag)
{
    update_debug();
};

#if MATCL_DEBUG_MP_FLOAT
    void mp_float::update_debug()
    {    
        m_debug = to_string(precision(10));
        m_prec  = this->get_precision();
    };
#endif

mp_float::mp_float()
{
    precision prec  = details::correct_prec(precision(0));

    mpfr_init2(mmd::impl_value(*this), (long)prec);
    mpfr_set_zero(mmd::impl_value(*this), 1);

    update_debug();
};

mp_float::mp_float(precision prec)
{
    prec            = details::correct_prec(prec);    

    mpfr_init2(mmd::impl_value(*this), (long)prec);
    mpfr_set_zero(mmd::impl_value(*this), 1);

    update_debug();
};

mp_float::mp_float(Integer v, precision prec)
{
    prec            = details::correct_prec(prec);    

    mpfr_init2(mmd::impl_value(*this), (long)prec);
    mpfr_set_si(mmd::impl_value(*this), v, MPFR_RNDN);

    update_debug();
};
mp_float::mp_float(Real v, precision prec)
{
    prec            = details::correct_prec(prec);    

    mpfr_init2(mmd::impl_value(*this), (long)prec);
    mpfr_set_d(mmd::impl_value(*this), v, MPFR_RNDN);

    update_debug();
};
mp_float::mp_float(long double v, precision prec)
{
    prec            = details::correct_prec(prec);    

    mpfr_init2(mmd::impl_value(*this), (long)prec);
    mpfr_set_ld(mmd::impl_value(*this), v, MPFR_RNDN);

    update_debug();
};

mp_float::mp_float(Float v, precision prec)
{
    prec            = details::correct_prec(prec);

    mpfr_init2(mmd::impl_value(*this), (long)prec);
    mpfr_set_flt(mmd::impl_value(*this), v, MPFR_RNDN);

    update_debug();
};

mp_float::mp_float(const mp_int& val, precision prec)
{
    using impl_int      = mmd::impl_int;
    const impl_int& v   = *mmd::get_ptr<mp_int>::eval(val.m_data);
    prec                = details::correct_prec(prec);

    mpfr_init2(mmd::impl_value(*this), (long)prec);
    mpfr_set_z(mmd::impl_value(*this), v.backend().data(), MPFR_RNDN);

    update_debug();
}
mp_float::mp_float(const mp_rational& val, precision prec)
{
    using impl_rat      = mmd::impl_rat;
    const impl_rat& v   = *mmd::get_ptr<mp_rational>::eval(val.m_data);
    prec                = details::correct_prec(prec);

    mpfr_init2(mmd::impl_value(*this), (long)prec);
    mpfr_set_q(mmd::impl_value(*this), v.backend().data(), MPFR_RNDN);

    update_debug();
}

mp_float::mp_float(const mp_float& other)
{
    precision prec          = details::correct_prec(other.get_precision());    

    mpfr_init2(mmd::impl_value(*this), (long)prec);
    mpfr_set(mmd::impl_value(*this), mmd::impl_value(other), MPFR_RNDN);

    update_debug();
};

mp_float::mp_float(mp_float&& other)
{
    using impl_type         = mmd::impl_float;
    impl_type& oval         = *mmd::get_ptr<mp_float>::eval(other.m_data);
    impl_type& val          = *mmd::get_ptr<mp_float>::eval(this->m_data);

    val->_mpfr_d = 0;
    mpfr_swap(val, oval);

    update_debug();
};

mp_float::mp_float(const mp_float& other, precision prec)
{
    prec                    = details::correct_prec(prec);    

    mpfr_init2(mmd::impl_value(*this), (long)prec);
    mpfr_set(mmd::impl_value(*this), mmd::impl_value(other), MPFR_RNDN);

    update_debug();
};

mp_float::mp_float(mp_float&& other, precision prec)
{
    using impl_type         = mmd::impl_float;
    impl_type& oval         = *mmd::get_ptr<mp_float>::eval(other.m_data);
    impl_type& val          = *mmd::get_ptr<mp_float>::eval(this->m_data);
    prec                    = details::correct_prec(prec);    
    precision old_prec      = other.get_precision();    

    if (old_prec == prec)
    {        
        val->_mpfr_d        = 0;
        mpfr_swap(val, oval);
    }
    else
    {
        mpfr_init2(mmd::impl_value(*this), (long)prec);
        mpfr_set(val, oval, MPFR_RNDN);
    };

    update_debug();
};

mp_float::mp_float(const std::string& s, precision prec, int base)
{
    using impl_type = mmd::impl_float;
    impl_type& val  = *mmd::get_ptr<mp_float>::eval(this->m_data);
    prec            = details::correct_prec(prec);    

    mpfr_init2(mmd::impl_value(*this), (long)prec);

    if (mpfr_set_str(val, s.c_str(), base, MPFR_RNDN) != 0)
    {
        mpfr_clear(val);
        throw std::invalid_argument ("mpfr_set_str");
    }

    update_debug();
}

mp_float::~mp_float()
{
    using impl_type = mmd::impl_float;
    impl_type& val  = *mmd::get_ptr<mp_float>::eval(this->m_data);

    if(val->_mpfr_d != 0)
        mpfr_clear(mmd::impl_value(*this));
};

mp_float& mp_float::operator=(Integer other) &
{
    using impl_type         = mmd::impl_float;
    impl_type& this_val     = *mmd::get_ptr<mp_float>::eval(this->m_data);

    mpfr_set_si(this_val, other, MPFR_RNDN);

    update_debug();
    return *this;
}
mp_float& mp_float::operator=(long double other) &
{
    using impl_type         = mmd::impl_float;
    impl_type& this_val     = *mmd::get_ptr<mp_float>::eval(this->m_data);

    mpfr_set_ld(this_val, other, MPFR_RNDN);

    update_debug();
    return *this;
}
mp_float& mp_float::operator=(Real other) &
{
    using impl_type         = mmd::impl_float;
    impl_type& this_val     = *mmd::get_ptr<mp_float>::eval(this->m_data);

    mpfr_set_d(this_val, other, MPFR_RNDN);

    update_debug();
    return *this;
}
mp_float& mp_float::operator=(Float other) &
{
    using impl_type         = mmd::impl_float;
    impl_type& this_val     = *mmd::get_ptr<mp_float>::eval(this->m_data);

    mpfr_set_flt(this_val, other, MPFR_RNDN);

    update_debug();
    return *this;
}

mp_float& mp_float::operator=(const mp_float& other) &
{
    using impl_type         = mmd::impl_float;
    const impl_type& val    = *mmd::get_ptr<mp_float>::eval(other.m_data);
    impl_type& this_val     = *mmd::get_ptr<mp_float>::eval(this->m_data);

    if (this_val != val)
    {
		precision prec_this = this->get_precision();
		precision prec_other= other.get_precision();

		if(prec_this != prec_other)
        {
			this->~mp_float();
            mpfr_init2(mmd::impl_value(*this), (long)prec_other);
		}

        mpfr_set(this_val, val, MPFR_RNDN);
    };

    update_debug();
    return *this;
};

mp_float& mp_float::operator=(mp_float&& other) &
{
    using impl_type         = mmd::impl_float;
    impl_type& val          = *mmd::get_ptr<mp_float>::eval(other.m_data);
    impl_type& this_val     = *mmd::get_ptr<mp_float>::eval(this->m_data);

    mpfr_swap(this_val, val);

    update_debug();
    return *this;
};

mp_float::operator bool() const
{
    using impl_type             = mmd::impl_float;
    const impl_type& this_val   = *mmd::get_ptr<mp_float>::eval(this->m_data);

    return  mpfr_zero_p(this_val) == 0;
};

Real mp_float::cast_float() const
{
    using impl_type             = mmd::impl_float;
    const impl_type& this_val   = *mmd::get_ptr<mp_float>::eval(this->m_data);

    return mpfr_get_d(this_val, MPFR_RNDN);
};

Integer mp_float::cast_int(round_mode rm) const
{
    if (can_cast_int(rm) == false)
        return 0;

    switch (rm)
    {
        case round_mode::floor:
            return mpfr_get_si(mmd::impl_value(*this), MPFR_RNDD);
            break;
        case round_mode::ceil:
            return mpfr_get_si(mmd::impl_value(*this), MPFR_RNDU);
            break;
        case round_mode::trunc:
            return mpfr_get_si(mmd::impl_value(*this), MPFR_RNDZ);
            break;
        case round_mode::round:
        default:
            return mpfr_get_si(mmd::impl_value(*this), MPFR_RNDNA);
            break;
    };
};


mp_int mp_float::cast_mp_int(round_mode rm) const
{
    using impl_type     = mmd::impl_int;
    using impl_float    = mmd::impl_float;
    const impl_float& v = *mmd::get_ptr<mp_float>::eval(this->m_data);

    mp_int ret(mp_int::impl_tag{});

    impl_type* iv       = mmd::get_ptr<mp_int>::eval(ret.m_data);    

    // call mpz_init
    new (iv) impl_type();

    switch (rm)
    {
        case round_mode::floor:
            mpfr_get_z(iv->backend().data(), v, MPFR_RNDD);
            break;
        case round_mode::ceil:
            mpfr_get_z(iv->backend().data(), v, MPFR_RNDU);
            break;
        case round_mode::trunc:
            mpfr_get_z(iv->backend().data(), v, MPFR_RNDZ);
            break;
        case round_mode::round:
        default:
            mpfr_get_z(iv->backend().data(), v, MPFR_RNDNA);
            break;
    };
    return ret;
}
mp_rational mp_float::cast_mp_rational() const
{
    return mp_rational(*this);
}
mp_complex mp_float::cast_mp_complex() const
{
    return mp_complex(*this);
}
bool mp_float::can_cast_int(round_mode rm) const
{
    switch (rm)
    {
        case round_mode::floor:
            return mpfr_fits_slong_p(mmd::impl_value(*this), MPFR_RNDD) != 0;
            break;
        case round_mode::ceil:
            return mpfr_fits_slong_p(mmd::impl_value(*this), MPFR_RNDU) != 0;
            break;
        case round_mode::trunc:
            return mpfr_fits_slong_p(mmd::impl_value(*this), MPFR_RNDZ) != 0;
            break;
        case round_mode::round:
        default:
            return mpfr_fits_slong_p(mmd::impl_value(*this), MPFR_RNDNA) != 0;
            break;
    };
};

// digits10 = ceil(bits*log[10](2))
precision mp_float::digits10_to_bits(precision d)
{
    const double LOG2_10 = 3.3219280948873624;

    return precision(size_t(std::ceil( double(d) * LOG2_10 )));
}

// bits = ceil(digits*log[2](10))
precision mp_float::bits_to_digits10(precision b)
{
    const double LOG10_2 = 0.30102999566398119;

    return precision(size_t(std::ceil( double(b) * LOG10_2 )));
}

precision mp_float::get_precision() const
{
    using impl_type     = mmd::impl_float;
    const impl_type& v1 = *mmd::get_ptr<mp_float>::eval(this->m_data);

    size_t prec         = mpfr_get_prec(v1);
    return precision(prec);
};
mp_float& mp_float::set_precision(precision prec)
{
    if (prec == this->get_precision())
        return *this;

    using impl_type = mmd::impl_float;
    impl_type& v1   = *mmd::get_ptr<mp_float>::eval(this->m_data);

    mpfr_prec_round(v1, (long)prec, MPFR_RNDN);
    update_debug();

    return *this;
};

mp_float& mp_float::clear(precision prec)
{
    using impl_type = mmd::impl_float;
    impl_type& v1   = *mmd::get_ptr<mp_float>::eval(this->m_data);
    prec            = details::correct_prec(prec);    

    mpfr_set_prec(v1, (long)prec);
    update_debug();

    return *this;
};

precision mp_float::get_default_precision()
{
    return precision(mpfr_get_default_prec());
};
void mp_float::set_default_precision(precision prec)
{
    mpfr_set_default_prec((long)prec);
};

precision mp_float::max_allowed_precision()
{
    return precision((size_t)-1);
};

precision mp_float::min_allowed_precision()
{
    return precision(0);
};

std::string mp_float::to_string(precision n) const
{
    std::ostringstream format;    

    if (n > 0)
    {
        format << "%." << n << "RNg";
    }
    else
    {
        format << "%." << "Re";
    }

    return to_string(format.str());
};

std::string mp_float::to_string_binary(precision n) const
{
    size_t digits           = (n > 0) ? n : this->get_precision() - 1;

    std::ostringstream format;    
    
    format << "%." << digits << "RNb";

    return to_string(format.str());
};

std::string mp_float::to_string_hex(precision n) const
{
    size_t digits           = (n > 0) ? n : (this->get_precision() - 1) / 4;

    std::ostringstream format;        
    format << "%." << digits << "RNa";

    return to_string(format.str());
};

std::string mp_float::to_string(const std::string& format) const
{
    if (format.empty() == true)
        return std::string();

    using impl_type         = mmd::impl_float;
    const impl_type& v1     = *mmd::get_ptr<mp_float>::eval(this->m_data);

    char *s = nullptr;
    std::string out;

    if(!(mpfr_asprintf(&s, format.c_str(), v1) < 0))
    {
        out = std::string(s);

        mpfr_free_str(s);
    }

    return out;
};

void mp_float::serialize(oarchive_impl & ar, unsigned int version) const
{
    precision prec      = this->get_precision();
    std::string str     = to_string();    

    (void)version;
    ar << (size_t)prec;
    ar << str;
};

void mp_float::serialize(iarchive_impl & ar, unsigned int version)
{
    (void)version;

    size_t prec;
    std::string str;
    ar >> prec;
    ar >> str;

    *this = mp_float(str, precision(prec));
};

std::istream& matcl::operator>>(std::istream& is, mp_float& val)
{
    std::string str;
    is >> str;

    using impl_type         = mmd::impl_float;
    impl_type& v1           = *mmd::get_ptr<mp_float>::eval(val.m_data);

    mpfr_set_str(v1, str.c_str(), 10, MPFR_RNDN);

    return is;
}
std::ostream& matcl::operator<<(std::ostream& os, const mp_float& val)
{
    std::ostringstream format;
    const std::ios::fmtflags flags = os.flags();

    format << ((flags & std::ios::showpos) ? "%+" : "%");

    if (os.precision() >= 0)
    {
        format << '.' << os.precision() << "R*"
               << ((flags & std::ios::floatfield) == std::ios::fixed ? 'f' :
                   (flags & std::ios::floatfield) == std::ios::scientific ? 'e' :
                   'g');
    }
    else
    {
        format << "R*e";
    }

    using impl_type         = mmd::impl_float;
    const impl_type& v1     = *mmd::get_ptr<mp_float>::eval(val.m_data);

    char *s = nullptr;

    if(!(mpfr_asprintf(&s, format.str().c_str(), MPFR_RNDN, v1) < 0))
    {
        os << std::string(s);
        mpfr_free_str(s);
    }

    return os;
};

};
