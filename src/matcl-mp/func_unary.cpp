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

#include "matcl-mp/func_unary.h"
#include "utils/impl_types.h"
#include "utils/utils.h"
#include "utils/extend_precision.h"
#include "matcl-mp/constants.h"
#include "matcl-mp/func_binary.h"
#include "matcl-core/details/scalfunc_real.h"
#include "error.h"
#include <boost/functional/hash.hpp>

#pragma warning(push)
#pragma warning(disable: 4100) //unreferenced formal parameter

#include "gmp-impl.h"

#pragma warning(pop)
namespace matcl
{

namespace mmd = matcl::mp::details;
namespace mrd = matcl::raw::details;

static mp_float mult_zero(const mp_float& r, const mp_float& arg, precision req_prec)
{
    req_prec = mmd::result_prec(req_prec, r.get_precision());

    if (is_zero(arg))
        return mp_float(arg, req_prec);

    return mul(r, arg, req_prec);
}

//-----------------------------------------------------------------------------
//                                  mp_int
//-----------------------------------------------------------------------------
bool matcl::is_zero(const mp_int& val)
{
    using impl_type     = mmd::impl_int;
    const impl_type& v1 = *mmd::get_ptr<mp_int>::eval(val.m_data);

    int cmp             = mpz_cmp_si(v1.backend().data(), 0);
    
    return cmp == 0 ? true : false;
};

bool matcl::is_regular(const mp_int& val)
{
    return is_zero(val) == false;
}
bool matcl::is_normal(const mp_int& val)
{
    return is_zero(val) == false;
}

bool matcl::is_one(const mp_int& val)
{
    return val == 1;
}

mp_int matcl::operator-(const mp_int& a)
{
    using impl_type         = mmd::impl_int;
    const impl_type& v1     = *mmd::get_ptr<mp_int>::eval(a.m_data);

    return mp_int(mp_int::impl_tag(),impl_type(-v1));
};

mp_int matcl::uminus(const mp_int& a, precision req_prec)
{
    (void)req_prec;

    return operator-(a);
};

mp_int matcl::abs(const mp_int& a, precision req_prec)
{
    (void) req_prec;

    using impl_type         = mmd::impl_int;
    const impl_type& v1     = *mmd::get_ptr<mp_int>::eval(a.m_data);

    return mp_int(mp_int::impl_tag(),impl_type(abs(v1)));
}
mp_int matcl::abs2(const mp_int& val, precision req_prec)
{
    (void)req_prec;
    return val*val;
}
mp_float matcl::arg(const mp_int& val, precision req_prec)
{
    return arg(mp_float(val), req_prec);
}

//-----------------------------------------------------------------------------
//                                  mp_float
//-----------------------------------------------------------------------------
bool matcl::is_nan(const mp_float& val)
{
    using impl_type         = mmd::impl_float;
    const impl_type& v1     = *mmd::get_ptr<mp_float>::eval(val.m_data);
    int ret                 = mpfr_nan_p(v1);
    return (ret == 0) ? false : true;
}
bool matcl::is_inf(const mp_float& val)
{
    using impl_type         = mmd::impl_float;
    const impl_type& v1     = *mmd::get_ptr<mp_float>::eval(val.m_data);
    int ret                 = mpfr_inf_p(v1);
    return (ret == 0) ? false : true;
}
bool matcl::is_finite(const mp_float& val)
{
    using impl_type         = mmd::impl_float;
    const impl_type& v1     = *mmd::get_ptr<mp_float>::eval(val.m_data);
    int ret                 = mpfr_number_p(v1);
    return (ret == 0) ? false : true;
}
bool matcl::is_zero(const mp_float& val)
{
    using impl_type         = mmd::impl_float;
    const impl_type& v1     = *mmd::get_ptr<mp_float>::eval(val.m_data);
    int ret                 = mpfr_zero_p(v1);
    return (ret == 0) ? false : true;
}
bool matcl::is_one(const mp_float& val)
{
    return val == 1;
}

bool matcl::is_int(const mp_float& val)
{
    using impl_type         = mmd::impl_float;
    const impl_type& v1     = *mmd::get_ptr<mp_float>::eval(val.m_data);
    int ret                 = mpfr_integer_p(v1);
    return (ret == 0) ? false : true;
}
bool matcl::is_regular(const mp_float& val)
{
    using impl_type         = mmd::impl_float;
    const impl_type& v1     = *mmd::get_ptr<mp_float>::eval(val.m_data);
    int ret                 = mpfr_regular_p(v1);
    return (ret == 0) ? false : true;
}
bool matcl::is_normal(const mp_float& val)
{
    //mp_float does not have subnormal values
    return is_regular(val);
}

mp_float matcl::operator-(const mp_float& a)
{
    mp_float ret(a);

    using impl_type     = mmd::impl_float;
    impl_type& v1       = *mmd::get_ptr<mp_float>::eval(ret.m_data);

    mpfr_neg(v1, v1, MPFR_RNDN);
    ret.update_debug();
    return ret;
};
mp_float matcl::uminus(const mp_float& x, precision req_prec)
{
    req_prec = mmd::result_prec(req_prec, x.get_precision());

    if (req_prec == x.get_precision())
        return operator-(x);

    mp_float ret(req_prec);

    using impl_type     = mmd::impl_float;
    impl_type& vr       = *mmd::get_ptr<mp_float>::eval(ret.m_data);
    const impl_type& vx = *mmd::get_ptr<mp_float>::eval(x.m_data);

    mpfr_neg(vr, vx, MPFR_RNDN);
    ret.update_debug();
    return ret;
};

mp_float matcl::abs(const mp_float& x, precision req_prec)
{    
    req_prec = mmd::result_prec(req_prec, x.get_precision());
    mp_float ret(0.0, req_prec);

    using impl_type         = mmd::impl_float;
    const impl_type& v1     = *mmd::get_ptr<mp_float>::eval(x.m_data);
    impl_type& vr           = *mmd::get_ptr<mp_float>::eval(ret.m_data);

    mpfr_abs(vr, v1, MPFR_RNDN);
    ret.update_debug();
    return ret;
}

mp_float matcl::abs2(const mp_float& x, precision req_prec)
{   
    return mul(x, x, req_prec);
}

mp_float matcl::arg(const mp_float& x, precision req_prec)
{
    req_prec = mmd::result_prec(req_prec, x.get_precision());

    if (is_nan(x) == true)
        return constants::mp_nan(req_prec);
 
    return signbit(x) == false ? mp_float(req_prec) : constants::mp_pi(req_prec);
}

//-----------------------------------------------------------------------------
//                                  mp_rational
//-----------------------------------------------------------------------------
mp_rational matcl::operator-(const mp_rational& a)
{
    using impl_type         = mmd::impl_rat;
    const impl_type& v1     = *mmd::get_ptr<mp_rational>::eval(a.m_data);

    return mp_rational(mp_rational::impl_tag(),impl_type(-v1));
};
mp_rational matcl::uminus(const mp_rational& a, precision req_prec)
{
    (void)req_prec;

    return operator-(a);
};

mp_rational matcl::abs(const mp_rational& a, precision req_prec)
{
    (void)req_prec;

    using impl_type         = mmd::impl_rat;
    const impl_type& v1     = *mmd::get_ptr<mp_rational>::eval(a.m_data);

    return mp_rational(mp_rational::impl_tag(),impl_type(abs(v1)));
};

mp_rational matcl::abs2(const mp_rational& val, precision req_prec)
{
    (void)req_prec;

    return val * val;
};

mp_float matcl::arg(const mp_rational& val, precision req_prec)
{
    (void)req_prec;
    return arg(mp_float(val));
}

bool matcl::is_zero(const mp_rational& val)
{
    return is_zero(val.numerator());
};

bool matcl::is_regular(const mp_rational& val)
{
    return is_zero(val) == false;
}

bool matcl::is_normal(const mp_rational& val)
{
    return is_zero(val) == false;
}

bool matcl::is_int(const mp_rational& val)
{
    return is_one(val.denominator()) == true;
};

bool matcl::is_one(const mp_rational& val)
{
    return val == 1;
};

//-----------------------------------------------------------------------------
//                                  mp_complex
//-----------------------------------------------------------------------------
mp_complex matcl::operator-(const mp_complex& a)
{
    return mp_complex(-a.real(),-a.imag());
}

mp_complex matcl::uminus(const mp_complex& a, precision req_prec)
{
    return mp_complex(uminus(a.real(), req_prec), uminus(a.imag(), req_prec));
}

mp_complex matcl::conj(const mp_complex& a)
{
    return mp_complex(a.real(), -a.imag());
}

mp_float matcl::abs(const mp_complex& val, precision req_prec)
{
    return hypot(val.real(), val.imag(), req_prec);
}
mp_float matcl::abs2(const mp_complex& x, precision req_prec)
{   
    req_prec            = mmd::result_prec(req_prec, x.get_precision());
    precision int_prec  = mmd::extend_prec_abs2(req_prec);

    mp_float res        = plus(abs2(x.real(), int_prec), abs2(x.imag(), int_prec), req_prec);
    return res;
};

mp_float matcl::arg(const mp_complex& val, precision req_prec)
{
    return atan2(val.imag(), val.real(), req_prec);
}

bool matcl::is_nan(const mp_complex& val)
{
    return is_nan(val.real()) || is_nan(val.imag());
};

bool matcl::is_inf(const mp_complex& val)
{
    return is_inf(val.real()) || is_inf(val.imag());
};

bool matcl::is_finite(const mp_complex& val)
{
    return is_finite(val.real()) && is_finite(val.imag());
};

bool matcl::is_regular(const mp_complex& val)
{
    return is_finite(val.real()) && is_finite(val.imag()) && val != 0;
}
bool matcl::is_normal(const mp_complex& val)
{
    //mp_complex does not have subnormal values
    return is_regular(val);
}

bool matcl::is_int(const mp_complex& val)
{
    return is_real(val) && is_int(val.real());
}

bool matcl::is_real(const mp_complex& val)
{
    return is_zero(val.imag());
};

bool matcl::is_zero(const mp_complex& val)
{
    return is_zero(val.real()) && is_zero(val.imag());
}

bool matcl::is_one(const mp_complex& val)
{
    return is_zero(val.imag()) && is_one(val.real());
}

//------------------------------------------------------
fp_type matcl::fpclassify(const mp_int& val)
{
    return fpclassify(mp_float(val));
}
fp_type matcl::fpclassify(const mp_float& val)
{
    if (is_normal(val))
        return fp_type::fp_normal;
    if (is_inf(val))
        return fp_type::fp_infinite;
    if (is_nan(val))
        return fp_type::fp_nan;
    if (is_zero(val))
        return fp_type::fp_zero;

    //subnormal values are not supported by mfpr
    return fp_type::fp_unknown;
};

fp_type matcl::fpclassify(const mp_rational& val)
{
    return fpclassify(mp_float(val));
}

//------------------------------------------------------
size_t matcl::hash_value(const mp_int& val)
{
    using impl_int          = mmd::impl_int;

    const impl_int& imp     = *mmd::get_ptr<mp_int>::eval(val.m_data);
    const mpz_t& imp0       = imp.backend().data();

    int size                = imp0[0]._mp_size;
    size_t seed             = (size_t)size;

    size                    = std::abs(size);
    
    for (int i = 0; i < size; ++i)
        boost::hash_combine(seed, size_t(imp0[0]._mp_d[i]));

    return seed;
}

size_t matcl::hash_value(const mp_float& val)
{
    if (is_zero(val) == true)
        return 0;
    
    using impl_float        = mmd::impl_float;

    const impl_float& v1    = *mmd::get_ptr<mp_float>::eval(val.m_data);
    
    int size                = v1[0]._mpfr_prec / BITS_PER_MP_LIMB;
    size                    += (v1[0]._mpfr_prec % BITS_PER_MP_LIMB)? 1 : 0;

    std::size_t seed        = 0;    

    for (int i = 0; i < size; ++i)
        boost::hash_combine(seed, v1[0]._mpfr_d[i]);

    boost::hash_combine(seed, v1[0]._mpfr_exp);
    boost::hash_combine(seed, v1[0]._mpfr_sign);

    return seed;
}

size_t matcl::hash_value(const mp_rational& val)
{
    size_t seed = hash_value(val.numerator());
    boost::hash_combine(seed, hash_value(val.denominator()));

    return seed;
}

size_t matcl::hash_value(const mp_complex& val)
{
    size_t seed = hash_value(val.real());
    boost::hash_combine(seed, hash_value(val.imag()));

    return seed;
}

//------------------------------------------------------
mp_float matcl::fma(const mp_float& a, const mp_float& b, const mp_float& c,
                    precision p)
{
    precision prec  = mmd::result_prec(p, a.get_precision(), b.get_precision(), 
                                           c.get_precision());

    mp_float ret(0.0, prec);    
    mpfr_fma(mmd::impl_value(ret), mmd::impl_value(a), mmd::impl_value(b),
             mmd::impl_value(c), MPFR_RNDN);
    ret.update_debug();

    return ret;
};
mp_float matcl::fms(const mp_float& a, const mp_float& b, const mp_float& c,
                    precision p)
{
    precision prec  = mmd::result_prec(p, a.get_precision(), b.get_precision(), 
                                           c.get_precision());

    mp_float ret(0.0, prec);    
    mpfr_fms(mmd::impl_value(ret), mmd::impl_value(a), mmd::impl_value(b),
             mmd::impl_value(c), MPFR_RNDN);
    ret.update_debug();

    return ret;
};

//------------------------------------------------------
MATCL_MP_EXPORT
mp_float matcl::dot2_a(const mp_float& a, const mp_float& b, const mp_float& c,
                        const mp_float& d, precision req_prec)
{
    //Kahan's algorithm

    //error is 2 ulp as shown in
    //C.-P. Jeannerod, N. Louvet, and J.-M. Muller, “Further analysis of
    //Kahan’s algorithm for the accurate computation of 2 × 2 determinants”

    size_t p_max        = std::max<size_t>(std::max<size_t>(a.get_precision(), b.get_precision()), 
                                           std::max<size_t>(c.get_precision(), d.get_precision()));

    req_prec            = mmd::result_prec(req_prec,  precision(p_max));
    precision int_prec  = mmd::extend_prec_dot2_a(req_prec);

    mp_float w = mul(c, d, int_prec);       //w = c * d
    mp_float e = fms(c, d, w, int_prec);   //e = c * d - w
    mp_float f = fma(a, b, w, int_prec);   //f = a * b + w 
    mp_float g = plus(f, e, req_prec);     //g = f + e;
    return g;
};

//------------------------------------------------------
mp_rational matcl::inv(const mp_int& x, precision)
{
    return mp_rational(1, x);
};
mp_float matcl::inv(const mp_float& x, precision req_prec)
{
    req_prec = mmd::result_prec(req_prec, x.get_precision());
    return div(1, x, req_prec);
};
mp_rational matcl::inv(const mp_rational& x, precision req_prec)
{
    return div(1, x, req_prec);
};
mp_complex matcl::inv(const mp_complex& x, precision req_prec)
{
    return mp::details::inv_impl(x, req_prec);
};

//------------------------------------------------------
mp_float matcl::nextabove(const mp_int& x)
{
    return nextabove(mp_float(x));
}
mp_float matcl::nextabove(const mp_float& x)
{
    mp_float a(x);

    using impl_type     = mmd::impl_float;
    impl_type& ai       = *mmd::get_ptr<mp_float>::eval(a.m_data);

    mpfr_nextabove(ai);
    a.update_debug();
    return a;
};
mp_float matcl::nextabove(const mp_rational& x)
{
    return nextabove(mp_float(x));
}

//------------------------------------------------------
mp_float matcl::nextbelow(const mp_int& x)
{
    return nextbelow(mp_float(x));
}
mp_float matcl::nextbelow(const mp_float& x)
{
    mp_float a(x);

    using impl_type     = mmd::impl_float;
    impl_type& ai       = *mmd::get_ptr<mp_float>::eval(a.m_data);

    mpfr_nextbelow(ai);
    a.update_debug();
    return a;
}
mp_float matcl::nextbelow(const mp_rational& x)
{
    return nextbelow(mp_float(x));
}

//------------------------------------------------------
mp_int matcl::sign(const mp_int& x, precision)
{
    if (is_zero(x) == true)
        return 0;
    else if (x < 0)
        return -1;
    else
        return 1;
};
mp_float matcl::sign(const mp_float& x, precision req_prec)
{
    req_prec = mmd::result_prec(req_prec, x.get_precision());

    if (x > 0)      return mp_float(1., req_prec);
    else if (x < 0) return mp_float(-1., req_prec);
    else            return mp_float(x, req_prec);
}
mp_rational matcl::sign(const mp_rational& x, precision)
{
    if (is_zero(x) == true)
        return 0;
    else if (x < 0)
        return -1;
    else
        return 1;
}
mp_complex matcl::sign(const mp_complex& x, precision req_prec)
{
    req_prec = mmd::result_prec(req_prec, x.get_precision());

    if (is_regular(x) == false)
    {
        if (is_zero(x) == true || is_nan(x) == true)
            return mp_complex(x, req_prec);

        bool inf_re = is_inf(real(x));
        bool inf_im = is_inf(imag(x));

        if (inf_re == true && inf_im == false)
            return mp_complex(sign(real(x), req_prec));
        if (inf_re == false && inf_im == true)
            return mp_complex(mp_float(req_prec), sign(imag(x), req_prec));
        else
            return mp_complex(constants::mp_nan(req_prec));
    };

    precision int_prec  = mmd::extend_prec_sign(req_prec);

    mp_float xa = abs(x, int_prec);             //0.5 ulp
    mp_float re = div(real(x), xa, int_prec);   //1 ulp
    mp_float im = div(imag(x), xa, int_prec);   //1 ulp

    return mp_complex(re, im, req_prec);        //1.5 ulp
}

//------------------------------------------------------
Integer matcl::isign(const mp_int& x)
{
    if (is_zero(x) == true)
        return 0;
    else if (x < 0)
        return -1;
    else
        return 1;
};
Integer matcl::isign(const mp_float& x)
{
    if (is_zero(x))
        return 0;
    else if (signbit(x) == true)
        return -1;
    else
        return 1;
}
Integer matcl::isign(const mp_rational& x)
{
    if (is_zero(x) == true)
        return 0;
    else if (x < 0)
        return -1;
    else
        return 1;
}

//------------------------------------------------------
bool matcl::signbit(const mp_int& x)
{
    return x < 0;
}
bool matcl::signbit(const mp_float& x)
{
    using impl_type     = mmd::impl_float;
    const impl_type& ai = *mmd::get_ptr<mp_float>::eval(x.m_data);

    return mpfr_signbit(ai) != 0;
}
bool matcl::signbit(const mp_rational& x)
{
    return x < 0;
}

//------------------------------------------------------
mp_rational matcl::ldexp(const mp_int& x, Integer n)
{
    if (n >= 0)
    {
        mp_int ret;
        mpz_mul_2exp(mmd::impl_value(ret).backend().data(), 
                     mmd::impl_value(x).backend().data(), (size_t)n);
        return ret;
    }
    else
    {
        return ldexp(mp_rational(x), n);
    }
}
mp_float matcl::ldexp(const mp_float& x, Integer n)
{
    mp_float ret(x);
    mpfr_mul_2si(mmd::impl_value(ret), mmd::impl_value(x), n, MPFR_RNDN);
    return ret;
};
mp_rational matcl::ldexp(const mp_rational& x, Integer n)
{    
    if (n >= 0)
    {
        mp_rational ret;
        mpq_mul_2exp(mmd::impl_value(ret).backend().data(), 
                     mmd::impl_value(x).backend().data(), (size_t)n);
        return ret;
    }
    else
    {
        mp_rational ret;
        mpq_div_2exp(mmd::impl_value(ret).backend().data(), 
                     mmd::impl_value(x).backend().data(), (size_t)(-n));

        return ret;
    }
}
mp_complex matcl::ldexp(const mp_complex& x, Integer n)
{
    return mp_complex(ldexp(real(x),n), ldexp(imag(x),n));
}

//------------------------------------------------------
mp_rational matcl::scalbn(const mp_int& x, Integer n)
{
    return ldexp(x, n);
}
mp_float matcl::scalbn(const mp_float& x, Integer n)
{
    return ldexp(x, n);
}
mp_rational matcl::scalbn(const mp_rational& x, Integer n)
{
    return ldexp(x, n);
}
mp_complex matcl::scalbn(const mp_complex& x, Integer n)
{
    return ldexp(x, n);
}

//------------------------------------------------------
mp_float matcl::frexp(const mp_int& x, Integer& exp)
{
    return frexp(mp_float(x), exp);
}
mp_float matcl::frexp(const mp_float& x, Integer& exp)
{
    mp_float ret(x);
    mpfr_exp_t exp2;
    mpfr_frexp(&exp2, mmd::impl_value(ret), mmd::impl_value(x), MPFR_RNDN);
    exp = exp2;

    ret.update_debug();

    return ret;    
}
mp_float matcl::frexp(const mp_rational& x, Integer& exp)
{
    return frexp(mp_float(x), exp);
}

//------------------------------------------------------
mp_float matcl::modf(const mp_int& x, mp_float& int_part)
{
    return modf(mp_float(x), int_part);
}
mp_float matcl::modf(const mp_float& x, mp_float& int_part)
{
    mp_float ret(x);
    mpfr_modf(mmd::impl_value(int_part), mmd::impl_value(ret),  
              mmd::impl_value(x), MPFR_RNDN);

    int_part.update_debug();
    ret.update_debug();

    return ret;    
}
mp_float matcl::modf(const mp_rational& x, mp_float& int_part)
{
    return modf(mp_float(x), int_part);
}

//------------------------------------------------------
mp_float matcl::logb(const mp_int& x)
{
    return logb(mp_float(x));
}
mp_float matcl::logb(const mp_float& x)
{
    if (is_regular(x) == false)
    {
        if (is_zero(x) == true)
            return -constants::mp_inf(x.get_precision());
        if (is_nan(x) == true)
            return constants::mp_nan(x.get_precision());
        if (is_inf(x) == true)
            return constants::mp_inf(x.get_precision());
    };

    Integer exp;
    frexp(x, exp);
    return mp_float(exp - 1, x.get_precision());
}
mp_float matcl::logb(const mp_rational& x)
{
    return logb(mp_float(x));
}
mp_float matcl::logb(const mp_complex& x)
{
    return logb(abs(x));
}

//------------------------------------------------------
Integer matcl::ilogb(const mp_int& x)
{
    return ilogb(mp_float(x));
}
Integer matcl::ilogb(const mp_float& x)
{
    if (is_regular(x) == false)
    {
        if (is_zero(x) == true)
            return FP_ILOGB0;
        if (is_nan(x) == true)
            return FP_ILOGBNAN;
        if (is_inf(x) == true)
            return INT_MAX;
    };

    Integer exp;
    frexp(x, exp);
    return exp - 1;
}
Integer matcl::ilogb(const mp_rational& x)
{
    return ilogb(mp_float(x));
}
Integer matcl::ilogb(const mp_complex& x)
{
    return ilogb(abs(x));
}

//------------------------------------------------------
mp_float matcl::eps(const mp_int& x)
{
    return eps(mp_float(x));
}
mp_float matcl::eps(const mp_float& x)
{
    if (x < 0)
        return nextabove(-x) + x;
    else
        return nextabove( x) - x;
};
mp_float matcl::eps(const mp_rational& x)
{
    return eps(mp_float(x));
}
mp_float matcl::eps(const mp_complex& x)
{
    mp_float a = abs(x);
    return nextabove(a) - a;
};

//------------------------------------------------------
mp_float matcl::sqrt(const mp_int& x, precision prec)
{
    return sqrt(mp_float(x), prec);
}
mp_float matcl::sqrt(const mp_float& x, precision req_prec)
{
    req_prec = mmd::result_prec(req_prec, x.get_precision());
    mp_float ret(0.0, req_prec);    
    mpfr_sqrt(mmd::impl_value(ret), mmd::impl_value(x), MPFR_RNDN);
    ret.update_debug();

    return ret;
};
mp_float matcl::sqrt(const mp_rational& x, precision req_prec)
{
    return sqrt(mp_float(x), req_prec);
}
mp_complex matcl::sqrt(const mp_complex& x, precision req_prec)
{
    req_prec            = mmd::result_prec(req_prec, x.get_precision());
    precision int_prec  = mmd::extend_prec_sqrt(req_prec);

    mp_float x_abs  = abs(x, int_prec);                 // 1 ulp

    if (is_regular(x_abs) == false)
    {
        x_abs.set_precision(req_prec);
        return mp_complex(x_abs, x_abs);

        //sqrt(inf - inf*i) = inf + inf * 1i in VS
        //                  = inf - inf * 1i in Matlab
    };

    // compute in safest quadrant    
    mp_float r_abs  = abs(real(x), int_prec);           //0.5 ulp(rounding)
    mp_float re     = sqrt((r_abs + x_abs) / 2, int_prec);//sqrt(1.5 ulp) = 1/2*1.5 + 0.5 ~ 1.5
    mp_float im     = div(imag(x), 2 * re, req_prec);   // 2 ulp
    re              = re.set_precision(req_prec);       // 2 ulp

	if (real(x) >= 0)
		return mp_complex(re, im);
	else if (imag(x) < 0)
		return mp_complex(-im, -re);
	else
		return mp_complex( im, re);
}

//------------------------------------------------------
mp_float matcl::cbrt(const mp_int& x, precision p)
{
    return cbrt(mp_float(x), p);
}
mp_float matcl::cbrt(const mp_float& x, precision req_prec)
{
    req_prec = mmd::result_prec(req_prec, x.get_precision());
    mp_float ret(0.0, req_prec);    
    mpfr_cbrt(mmd::impl_value(ret), mmd::impl_value(x), MPFR_RNDN);
    ret.update_debug();

    return ret;
}
mp_float matcl::cbrt(const mp_rational& x, precision p)
{
    return cbrt(mp_float(x), p);
}

//------------------------------------------------------
mp_float matcl::exp(const mp_int& x, precision prec)
{
    return exp(mp_float(x), prec);
}
mp_float matcl::exp(const mp_float& x, precision req_prec)
{
    req_prec = mmd::result_prec(req_prec, x.get_precision());
    mp_float ret(0.0, req_prec);    
    mpfr_exp(mmd::impl_value(ret), mmd::impl_value(x), MPFR_RNDN);
    ret.update_debug();

    return ret;
};
mp_float matcl::exp(const mp_rational& x, precision req_prec)
{
    return exp(mp_float(x), req_prec);
}

//------------------------------------------------------
mp_float matcl::exp2(const mp_int& x, precision prec)
{
    return exp2(mp_float(x), prec);
}
mp_float matcl::exp2(const mp_float& x, precision req_prec)
{
    req_prec = mmd::result_prec(req_prec, x.get_precision());
    mp_float ret(0.0, req_prec);    
    mpfr_exp2(mmd::impl_value(ret), mmd::impl_value(x), MPFR_RNDN);
    ret.update_debug();

    return ret;
};
mp_float matcl::exp2(const mp_rational& x, precision req_prec)
{
    return exp2(mp_float(x), req_prec);
}

//------------------------------------------------------
mp_float matcl::exp10(const mp_int& x, precision prec)
{
    return exp10(mp_float(x), prec);
}
mp_float matcl::exp10(const mp_float& x, precision req_prec)
{
    req_prec = mmd::result_prec(req_prec, x.get_precision());
    mp_float ret(0.0, req_prec);    
    mpfr_exp10(mmd::impl_value(ret), mmd::impl_value(x), MPFR_RNDN);
    ret.update_debug();

    return ret;
};
mp_float matcl::exp10(const mp_rational& x, precision req_prec)
{
    return exp10(mp_float(x), req_prec);
}

//------------------------------------------------------
mp_float matcl::expm1(const mp_int& x, precision prec)
{
    return expm1(mp_float(x), prec);
}
mp_float matcl::expm1(const mp_float& x, precision req_prec)
{
    req_prec = mmd::result_prec(req_prec, x.get_precision());
    mp_float ret(0.0, req_prec);    
    mpfr_expm1(mmd::impl_value(ret), mmd::impl_value(x), MPFR_RNDN);
    ret.update_debug();

    return ret;
};
mp_float matcl::expm1(const mp_rational& x, precision req_prec)
{
    return expm1(mp_float(x), req_prec);
}

mp_complex matcl::expi(const mp_int& x, precision p)
{
    return expi(mp_float(x), p);
};
mp_complex matcl::expi(const mp_float& x, precision req_prec)
{
    req_prec = mmd::result_prec(req_prec, x.get_precision());

    mp_float si, co;
    sin_cos(x, si, co, req_prec);

    return mp_complex(co, si);
}
mp_complex matcl::expi(const mp_rational& x, precision p)
{
    return expi(mp_float(x), p);
};
mp_complex matcl::expi(const mp_complex& x, precision p)
{
    return exp(mp_complex(-imag(x), real(x)), p);
}

//------------------------------------------------------
mp_float matcl::log(const mp_int& x, precision req_prec)
{
    return log(mp_float(x), req_prec);
}
mp_float matcl::log(const mp_float& x, precision req_prec)
{
    req_prec = mmd::result_prec(req_prec, x.get_precision());
    mp_float ret(0.0, req_prec);    
    mpfr_log(mmd::impl_value(ret), mmd::impl_value(x), MPFR_RNDN);
    ret.update_debug();

    return ret;
};
mp_float matcl::log(const mp_rational& x, precision req_prec)
{
    return log(mp_float(x), req_prec);
}

//------------------------------------------------------
mp_float matcl::log1p(const mp_int& x, precision req_prec)
{
    return log1p(mp_float(x), req_prec);
}
mp_float matcl::log1p(const mp_float& x, precision req_prec)
{
    req_prec = mmd::result_prec(req_prec, x.get_precision());
    mp_float ret(0.0, req_prec);    
    mpfr_log1p(mmd::impl_value(ret), mmd::impl_value(x), MPFR_RNDN);
    ret.update_debug();

    return ret;
};
mp_float matcl::log1p(const mp_rational& x, precision req_prec)
{
    return log1p(mp_float(x), req_prec);
}

//------------------------------------------------------
mp_float matcl::log2(const mp_int& x, precision req_prec)
{
    return log2(mp_float(x), req_prec);
}
mp_float matcl::log2(const mp_float& x, precision req_prec)
{
    req_prec = mmd::result_prec(req_prec, x.get_precision());
    mp_float ret(0.0, req_prec);    
    mpfr_log2(mmd::impl_value(ret), mmd::impl_value(x), MPFR_RNDN);
    ret.update_debug();

    return ret;
};
mp_float matcl::log2(const mp_rational& x, precision req_prec)
{
    return log2(mp_float(x), req_prec);
}
mp_complex matcl::log2(const mp_complex& x, precision req_prec)
{
    req_prec            = mmd::result_prec(req_prec, x.get_precision());
    precision int_prec  = mmd::extend_prec_log2_c(req_prec);
    return mul(log(x, int_prec), constants::mp_log2e(int_prec), req_prec); // 2 ulp
}

//------------------------------------------------------
mp_float matcl::log10(const mp_int& x, precision req_prec)
{
    return log10(mp_float(x), req_prec);
}
mp_float matcl::log10(const mp_float& x, precision req_prec)
{
    req_prec = mmd::result_prec(req_prec, x.get_precision());
    mp_float ret(0.0, req_prec);    
    mpfr_log10(mmd::impl_value(ret), mmd::impl_value(x), MPFR_RNDN);
    ret.update_debug();

    return ret;
};
mp_float matcl::log10(const mp_rational& x, precision req_prec)
{
    return log10(mp_float(x), req_prec);
}
mp_complex matcl::log10(const mp_complex& x, precision req_prec)
{
    req_prec = mmd::result_prec(req_prec, x.get_precision());
    precision int_prec  = mmd::extend_prec_log10_c(req_prec);
    return mul(log(x, int_prec), constants::mp_log10e(int_prec), req_prec); // 1.5 ulp
}

//------------------------------------------------------
void matcl::sin_cos(const mp_int& x, mp_float& sin, mp_float& cos, precision req_prec)
{
    return sin_cos(mp_float(x), sin, cos, req_prec);
}
void matcl::sin_cos(const mp_rational& x, mp_float& sin, mp_float& cos, precision req_prec)
{
    return sin_cos(mp_float(x), sin, cos, req_prec);
}
void matcl::sin_cos(const mp_float& x, mp_float& sin, mp_float& cos, precision req_prec)
{
    req_prec = mmd::result_prec(req_prec, x.get_precision());

    sin.clear(req_prec);
    cos.clear(req_prec);

    mpfr_sin_cos(mmd::impl_value(sin), mmd::impl_value(cos), 
                 mmd::impl_value(x), MPFR_RNDN);
    sin.update_debug();
    cos.update_debug();
    return;
}

//------------------------------------------------------
mp_float matcl::sin(const mp_int& x, precision req_prec)
{
    return sin(mp_float(x), req_prec);
}
mp_float matcl::sin(const mp_float& x, precision req_prec)
{
    req_prec = mmd::result_prec(req_prec, x.get_precision());

    mp_float ret(0.0, req_prec);    
    mpfr_sin(mmd::impl_value(ret), mmd::impl_value(x), MPFR_RNDN);
    ret.update_debug();
    return ret;
};
mp_float matcl::sin(const mp_rational& x, precision req_prec)
{
    return sin(mp_float(x), req_prec);
}
mp_complex matcl::sin(const mp_complex& x, precision req_prec)
{
    //sin(z) = sin(x) * cosh(y) + i cos(x) * sinh(y)
    const mp_float& xre = real(x);
    const mp_float& xim = imag(x);

    if (is_zero(xim) == true)
        return sin(xre, req_prec);

    req_prec                = mmd::result_prec(req_prec, x.get_precision());
    precision int_prec      = mmd::extend_prec_sin(req_prec);

    mp_float sin_re, cos_re;
    sin_cos(xre, sin_re, cos_re, int_prec);         // 0.5 ulp

    mp_float sinh_im, cosh_im;
    sinh_cosh(xim, sinh_im, cosh_im, int_prec);     // 0.5 ulp

    mp_float r_re       = mult_zero(cosh_im, sin_re, req_prec); // 1.5 ulp
    mp_float r_im       = mult_zero(sinh_im, cos_re, req_prec); // 1.5 ulp

    return mp_complex(r_re, r_im);
}

//------------------------------------------------------
mp_float matcl::cos(const mp_int& x, precision req_prec)
{
    return cos(mp_float(x), req_prec);
}
mp_float matcl::cos(const mp_float& x, precision req_prec)
{
    req_prec = mmd::result_prec(req_prec, x.get_precision());
    mp_float ret(0.0, req_prec);    

    mpfr_cos(mmd::impl_value(ret), mmd::impl_value(x), MPFR_RNDN);
    ret.update_debug();
    return ret;
};
mp_float matcl::cos(const mp_rational& x, precision req_prec)
{
    return cos(mp_float(x), req_prec);
}
mp_complex matcl::cos(const mp_complex& x, precision req_prec)
{
    //cos(z) = cos(x) * cosh(y) - i sin(x) * sinh(y)
    const mp_float& xre = real(x);
    const mp_float& xim = imag(x);

    req_prec                = mmd::result_prec(req_prec, x.get_precision());
    precision int_prec      = mmd::extend_prec_cos(req_prec);

    if (is_zero(xim) == true)
        return cos(xre, req_prec);

    mp_float sin_re, cos_re;
    sin_cos(xre, sin_re, cos_re, int_prec);     // 0.5 ulp

    mp_float sinh_im, cosh_im;
    sinh_cosh(xim, sinh_im, cosh_im, int_prec); // 0.5 ulp

    mp_float r_re       = mult_zero(cosh_im, cos_re, req_prec);     //1.5 ulp
    mp_float r_im       = mult_zero(sinh_im, -sin_re, req_prec);    //1.5 ulp

    return mp_complex(r_re, r_im);
}

//------------------------------------------------------
mp_float matcl::tan(const mp_int& x, precision req_prec)
{
    return tan(mp_float(x), req_prec);
}
mp_float matcl::tan(const mp_float& x, precision req_prec)
{
    req_prec = mmd::result_prec(req_prec, x.get_precision());
    mp_float ret(0.0, req_prec);    

    mpfr_tan(mmd::impl_value(ret), mmd::impl_value(x), MPFR_RNDN);
    ret.update_debug();
    return ret;
};
mp_float matcl::tan(const mp_rational& x, precision req_prec)
{
    return tan(mp_float(x), req_prec);
}
mp_complex matcl::tan(const mp_complex& x, precision req_prec)
{
    //tan(z) = -i * tanh(i*z)
    const mp_float& x_re = real(x);
    const mp_float& x_im = imag(x);

    if (is_zero(x_im) == true)
        return tan(x_re, req_prec);

    req_prec            = mmd::result_prec(req_prec, x.get_precision());

    mp_complex tan_iz   = tanh(mp_complex(-x_im, x_re), req_prec);
    const mp_float& re  = real(tan_iz);
    const mp_float& im  = imag(tan_iz);

    return mp_complex(std::move(im), -std::move(re));
}

//------------------------------------------------------
mp_float matcl::cot(const mp_int& x, precision req_prec)
{
    return cot(mp_float(x), req_prec);
}
mp_float matcl::cot(const mp_float& x, precision req_prec)
{
    req_prec = mmd::result_prec(req_prec, x.get_precision());
    mp_float ret(0.0, req_prec);    

    mpfr_cot(mmd::impl_value(ret), mmd::impl_value(x), MPFR_RNDN);
    ret.update_debug();
    return ret;
};
mp_float matcl::cot(const mp_rational& x, precision req_prec)
{
    return cot(mp_float(x), req_prec);
}
mp_complex matcl::cot(const mp_complex& x, precision req_prec)
{
    req_prec            = mmd::result_prec(req_prec, x.get_precision());
    precision int_prec  = mmd::extend_prec_cot(req_prec);
    return inv(tan(x, int_prec), req_prec);
}

//------------------------------------------------------
mp_float matcl::sec(const mp_int& x, precision req_prec)
{
    return sec(mp_float(x), req_prec);
}
mp_float matcl::sec(const mp_float& x, precision req_prec)
{
    req_prec = mmd::result_prec(req_prec, x.get_precision());
    mp_float ret(0.0, req_prec);    

    mpfr_sec(mmd::impl_value(ret), mmd::impl_value(x), MPFR_RNDN);
    ret.update_debug();
    return ret;
};
mp_float matcl::sec(const mp_rational& x, precision req_prec)
{
    return sec(mp_float(x), req_prec);
}
mp_complex matcl::sec(const mp_complex& x, precision req_prec)
{
    req_prec            = mmd::result_prec(req_prec, x.get_precision());
    precision int_prec  = mmd::extend_prec_sec(req_prec);
    return inv(cos(x, int_prec), req_prec);    //real: 4, imag: 4
}

//------------------------------------------------------
mp_float matcl::csc(const mp_int& x, precision req_prec)
{
    return csc(mp_float(x), req_prec);
}
mp_float matcl::csc(const mp_float& x, precision req_prec)
{
    req_prec = mmd::result_prec(req_prec, x.get_precision());
    mp_float ret(0.0, req_prec);    

    mpfr_csc(mmd::impl_value(ret), mmd::impl_value(x), MPFR_RNDN);
    ret.update_debug();
    return ret;
};
mp_float matcl::csc(const mp_rational& x, precision req_prec)
{
    return csc(mp_float(x), req_prec);
}
mp_complex matcl::csc(const mp_complex& x, precision req_prec)
{
    req_prec            = mmd::result_prec(req_prec, x.get_precision());
    precision int_prec  = mmd::extend_prec_csc(req_prec);

    return inv(sin(x, int_prec), req_prec);    //real: 4, imag: 4
}

//------------------------------------------------------
void matcl::sinh_cosh(const mp_int& x, mp_float& sinh, mp_float& cosh, precision req_prec)
{
    return sinh_cosh(mp_float(x), sinh, cosh, req_prec);
}
void matcl::sinh_cosh(const mp_rational& x, mp_float& sinh, mp_float& cosh, precision req_prec)
{
    return sinh_cosh(mp_float(x), sinh, cosh, req_prec);
}
void matcl::sinh_cosh(const mp_float& x, mp_float& sinh, mp_float& cosh, precision req_prec)
{
    req_prec = mmd::result_prec(req_prec, x.get_precision());

    sinh.clear(req_prec);
    cosh.clear(req_prec);

    mpfr_sinh_cosh(mmd::impl_value(sinh), mmd::impl_value(cosh), 
                 mmd::impl_value(x), MPFR_RNDN);

    sinh.update_debug();
    cosh.update_debug();
    return;
}

//------------------------------------------------------
mp_float matcl::sinh(const mp_int& x, precision req_prec)
{
    return sinh(mp_float(x), req_prec);
}
mp_float matcl::sinh(const mp_float& x, precision req_prec)
{
    req_prec = mmd::result_prec(req_prec, x.get_precision());
    mp_float ret(0.0, req_prec);    

    mpfr_sinh(mmd::impl_value(ret), mmd::impl_value(x), MPFR_RNDN);
    ret.update_debug();
    return ret;
};
mp_float matcl::sinh(const mp_rational& x, precision req_prec)
{
    return sinh(mp_float(x), req_prec);
}
mp_complex matcl::sinh(const mp_complex& x, precision req_prec)
{
    //sinh(z) = sinh(x) * cos(y) + i cosh(x) * sin(y)
    const mp_float& xre = real(x);
    const mp_float& xim = imag(x);

    if (is_zero(xim) == true)
        return sinh(xre, req_prec);

    req_prec            = mmd::result_prec(req_prec, x.get_precision());
    precision int_prec  = mmd::extend_prec_sinh(req_prec);

    mp_float sin_im, cos_im;
    sin_cos(xim, sin_im, cos_im, int_prec);

    mp_float sinh_re, cosh_re;
    sinh_cosh(xre, sinh_re, cosh_re, int_prec);

    mp_float r_re       = mult_zero(sinh_re, cos_im, req_prec);
    mp_float r_im       = mult_zero(cosh_re, sin_im, req_prec);

    return mp_complex(r_re, r_im);
}

//------------------------------------------------------
mp_float matcl::cosh(const mp_int& x, precision req_prec)
{
    return cosh(mp_float(x), req_prec);
}
mp_float matcl::cosh(const mp_float& x, precision req_prec)
{
    req_prec = mmd::result_prec(req_prec, x.get_precision());
    mp_float ret(0.0, req_prec);    

    mpfr_cosh(mmd::impl_value(ret), mmd::impl_value(x), MPFR_RNDN);
    ret.update_debug();
    return ret;
};
mp_float matcl::cosh(const mp_rational& x, precision req_prec)
{
    return cosh(mp_float(x), req_prec);
}
mp_complex matcl::cosh(const mp_complex& x, precision req_prec)
{
    //cosh(z) = cosh(x) * cos(y) + i sinh(x) * sin(y)
    const mp_float& xre = real(x);
    const mp_float& xim = imag(x);

    if (is_zero(xim) == true)
        return cosh(xre, req_prec);

    req_prec            = mmd::result_prec(req_prec, x.get_precision());
    precision int_prec  = mmd::extend_prec_cosh(req_prec);

    mp_float sin_im, cos_im;
    sin_cos(xim, sin_im, cos_im, int_prec);

    mp_float sinh_re, cosh_re;
    sinh_cosh(xre, sinh_re, cosh_re, int_prec);

    mp_float r_re       = mult_zero(cosh_re, cos_im, req_prec);
    mp_float r_im       = mult_zero(sinh_re, sin_im, req_prec);

    return mp_complex(r_re, r_im);
}

//------------------------------------------------------
mp_float matcl::tanh(const mp_int& x, precision req_prec)
{
    return tanh(mp_float(x), req_prec);
}
mp_float matcl::tanh(const mp_float& x, precision req_prec)
{
    req_prec = mmd::result_prec(req_prec, x.get_precision());
    mp_float ret(0.0, req_prec);    

    mpfr_tanh(mmd::impl_value(ret), mmd::impl_value(x), MPFR_RNDN);
    ret.update_debug();
    return ret;
}
mp_float matcl::tanh(const mp_rational& x, precision req_prec)
{
    return tanh(mp_float(x), req_prec);
}
mp_complex matcl::tanh(const mp_complex& x, precision req_prec)
{
    //tanh(z)   = sinh(z) / cosh(z) = num / den
    //num       = cosh(x) * sinh(x) * (1 + tan(y)^2) +  i*tan(y)
    //den       = 1+sinh(x)^2 * (1 + tan(y)^2)

    const mp_float& x_re    = real(x);
    const mp_float& x_im    = imag(x);

    if (is_zero(x_im) == true)
        return tanh(x_re, req_prec);

    req_prec                = mmd::result_prec(req_prec, x.get_precision());
    precision int_prec      = mmd::extend_prec_tanh(req_prec);
    mp_float tan_im         = tan(x_im, int_prec);

    mp_float sinh_re, cosh_re;
    sinh_cosh(x_re, sinh_re, cosh_re, int_prec);

    mp_float Bv         = sinh_re * (1 + tan_im * tan_im);  // 2.5 ulp
    mp_float Dv         = 1 + Bv * sinh_re;                 // 3.5 ulp

    if (is_inf(Dv) == true)
    {
        mp_float r_re       = (sinh_re < 0 ? -1.0 : 1.0);
        mp_float r_im       = tan_im * 0;
        return mp_complex(r_re, r_im, req_prec);
    }

    mp_float r_re       = div(cosh_re * Bv, Dv, req_prec);  //7.5 ulp 
    mp_float r_im       = div(tan_im , Dv, req_prec);       //4.5 ulp 

    return mp_complex(r_re, r_im);
};

//------------------------------------------------------
mp_float matcl::coth(const mp_int& x, precision req_prec)
{
    return coth(mp_float(x), req_prec);
}
mp_float matcl::coth(const mp_float& x, precision req_prec)
{
    req_prec = mmd::result_prec(req_prec, x.get_precision());
    mp_float ret(0.0, req_prec);    

    mpfr_coth(mmd::impl_value(ret), mmd::impl_value(x), MPFR_RNDN);
    ret.update_debug();
    return ret;
}
mp_float matcl::coth(const mp_rational& x, precision req_prec)
{
    return coth(mp_float(x), req_prec);
}
mp_complex matcl::coth(const mp_complex& x, precision req_prec)
{
    req_prec            = mmd::result_prec(req_prec, x.get_precision());
    precision int_prec  = mmd::extend_prec_coth(req_prec);
    return inv(tanh(x, int_prec), req_prec);
}

//------------------------------------------------------
mp_float matcl::sech(const mp_int& x, precision req_prec)
{
    return sech(mp_float(x), req_prec);
}
mp_float matcl::sech(const mp_float& x, precision req_prec)
{
    req_prec = mmd::result_prec(req_prec, x.get_precision());
    mp_float ret(0.0, req_prec);    

    mpfr_sech(mmd::impl_value(ret), mmd::impl_value(x), MPFR_RNDN);

    return ret;
}
mp_float matcl::sech(const mp_rational& x, precision req_prec)
{
    return sech(mp_float(x), req_prec);
}
mp_complex matcl::sech(const mp_complex& x, precision req_prec)
{
    req_prec            = mmd::result_prec(req_prec, x.get_precision());
    precision int_prec  = mmd::extend_prec_sech(req_prec);
    return inv(cosh(x, int_prec), req_prec);    //real: 4, imag: 4
}

//------------------------------------------------------
mp_float matcl::csch(const mp_int& x, precision req_prec)
{
    return csch(mp_float(x), req_prec);
}
mp_float matcl::csch(const mp_float& x, precision req_prec)
{
    req_prec = mmd::result_prec(req_prec, x.get_precision());
    mp_float ret(0.0, req_prec);    

    mpfr_csch(mmd::impl_value(ret), mmd::impl_value(x), MPFR_RNDN);
    ret.update_debug();
    return ret;
}
mp_float matcl::csch(const mp_rational& x, precision req_prec)
{
    return csch(mp_float(x), req_prec);
}
mp_complex matcl::csch(const mp_complex& x, precision req_prec)
{
    req_prec            = mmd::result_prec(req_prec, x.get_precision());
    precision int_prec  = mmd::extend_prec_csch(req_prec);
    return inv(sinh(x, int_prec), req_prec);    //real: 4, imag: 4
}

//------------------------------------------------------
mp_float matcl::asin(const mp_int& x, precision p)
{
    return asin(mp_float(x), p);
}
mp_float matcl::asin(const mp_float& x, precision req_prec)
{
    req_prec = mmd::result_prec(req_prec, x.get_precision());
    mp_float ret(0.0, req_prec);    

    mpfr_asin(mmd::impl_value(ret), mmd::impl_value(x), MPFR_RNDN);
    ret.update_debug();
    return ret;
};
mp_float matcl::asin(const mp_rational& x, precision p)
{
    return asin(mp_float(x), p);
};

//------------------------------------------------------
mp_float matcl::acos(const mp_int& x, precision p)
{
    return acos(mp_float(x), p);
}
mp_float matcl::acos(const mp_float& x, precision req_prec)
{
    req_prec = mmd::result_prec(req_prec, x.get_precision());
    mp_float ret(0.0, req_prec);    

    mpfr_acos(mmd::impl_value(ret), mmd::impl_value(x), MPFR_RNDN);
    ret.update_debug();
    return ret;
};
mp_float matcl::acos(const mp_rational& x, precision p)
{
    return acos(mp_float(x), p);
};

//------------------------------------------------------
mp_float matcl::atan(const mp_int& x, precision p)
{
    return atan(mp_float(x), p);
}
mp_float matcl::atan(const mp_float& x, precision req_prec)
{
    req_prec = mmd::result_prec(req_prec, x.get_precision());
    mp_float ret(0.0, req_prec);    

    mpfr_atan(mmd::impl_value(ret), mmd::impl_value(x), MPFR_RNDN);
    ret.update_debug();
    return ret;
};
mp_float matcl::atan(const mp_rational& x, precision p)
{
    return atan(mp_float(x), p);
};
mp_complex matcl::atan(const mp_complex& x, precision p)
{
    //atan(z) = -i atanh(iz)
    mp_complex res = atanh(mp_complex(-imag(x), real(x)), p);
    return mp_complex(std::move(imag(res)), -std::move(real(res)));
};

//------------------------------------------------------
mp_float matcl::asinh(const mp_int& x, precision p)
{
    return asinh(mp_float(x), p);
}
mp_float matcl::asinh(const mp_float& x, precision req_prec)
{
    req_prec = mmd::result_prec(req_prec, x.get_precision());
    mp_float ret(0.0, req_prec);    

    mpfr_asinh(mmd::impl_value(ret), mmd::impl_value(x), MPFR_RNDN);
    ret.update_debug();
    return ret;
};
mp_float matcl::asinh(const mp_rational& x, precision p)
{
    return asinh(mp_float(x), p);
};
mp_complex matcl::asinh(const mp_complex& x, precision p)
{
   // asinh(z) = i asin(-i z);
    mp_complex res = asin(mp_complex(imag(x), -real(x)), p);
    return mp_complex(-std::move(imag(res)), std::move(real(res)));
};

//------------------------------------------------------
mp_float matcl::acosh(const mp_int& x, precision p)
{
    return acosh(mp_float(x), p);
}
mp_float matcl::acosh(const mp_float& x, precision req_prec)
{
    req_prec = mmd::result_prec(req_prec, x.get_precision());
    mp_float ret(0.0, req_prec);    

    mpfr_acosh(mmd::impl_value(ret), mmd::impl_value(x), MPFR_RNDN);
    ret.update_debug();
    return ret;
};
mp_float matcl::acosh(const mp_rational& x, precision p)
{
    return acosh(mp_float(x), p);
};
mp_complex matcl::acosh(const mp_complex& x, precision p)
{
    mp_complex res = acos(x, p);

   // We use the relation acosh(z) = +-i acos(z)
   // Choosing the sign of multiplier to give real(acosh(z)) >= 0
   // as well as compatibility with C99.
  //
    if(is_nan(res.imag()) == false && signbit(res.imag()) ==  true)
        return mp_complex(-imag(res), real(res));
    else
        return mp_complex(imag(res), -real(res));
};

//------------------------------------------------------
mp_float matcl::atanh(const mp_int& x, precision p)
{
    return atanh(mp_float(x), p);
}
mp_float matcl::atanh(const mp_float& x, precision req_prec)
{
    req_prec = mmd::result_prec(req_prec, x.get_precision());
    mp_float ret(0.0, req_prec);    

    mpfr_atanh(mmd::impl_value(ret), mmd::impl_value(x), MPFR_RNDN);
    ret.update_debug();
    return ret;
};
mp_float matcl::atanh(const mp_rational& x, precision p)
{
    return atanh(mp_float(x), p);
};

//------------------------------------------------------
mp_int matcl::floor(const mp_int& x)
{
    return x;
}
mp_float matcl::floor(const mp_float& x)
{
    mp_float ret(0.0, x.get_precision());

    mpfr_floor(mmd::impl_value(ret), mmd::impl_value(x));
    ret.update_debug();
    return ret;
}
mp_float matcl::floor(const mp_rational& x)
{
    return floor(mp_float(x));
}
mp_complex matcl::floor(const mp_complex& x)
{
    return mp_complex(floor(real(x)), floor(imag(x)));
}

//------------------------------------------------------
mp_int matcl::ceil(const mp_int& x)
{
    return x;
}
mp_float matcl::ceil(const mp_float& x)
{
    mp_float ret(0.0, x.get_precision());

    mpfr_ceil(mmd::impl_value(ret), mmd::impl_value(x));
    ret.update_debug();
    return ret;
}
mp_float matcl::ceil(const mp_rational& x)
{
    return ceil(mp_float(x));
}
mp_complex matcl::ceil(const mp_complex& x)
{
    return mp_complex(ceil(real(x)), ceil(imag(x)));
}

//------------------------------------------------------
mp_int matcl::round(const mp_int& x)
{
    return x;
}
mp_float matcl::round(const mp_float& x)
{
    mp_float ret(0.0, x.get_precision());

    mpfr_rint(mmd::impl_value(ret), mmd::impl_value(x), MPFR_RNDN);
    ret.update_debug();
    return ret;
}
mp_float matcl::round(const mp_rational& x)
{
    return round(mp_float(x));
}
mp_complex matcl::round(const mp_complex& x)
{
    return mp_complex(round(real(x)), round(imag(x)));
}

//------------------------------------------------------
mp_int matcl::trunc(const mp_int& x)
{
    return x;
}
mp_float matcl::trunc(const mp_float& x)
{
    mp_float ret(0.0, x.get_precision());

    mpfr_trunc(mmd::impl_value(ret), mmd::impl_value(x));
    ret.update_debug();
    return ret;
}
mp_float matcl::trunc(const mp_rational& x)
{
    return trunc(mp_float(x));
}
mp_complex matcl::trunc(const mp_complex& x)
{
    return mp_complex(trunc(real(x)), trunc(imag(x)));
}

//------------------------------------------------------
mp_int matcl::ifloor(const mp_int& x)
{
    return x;
}
mp_int matcl::ifloor(const mp_float& x)
{
    return floor(x).cast_mp_int();
}
mp_int matcl::ifloor(const mp_rational& x)
{
    return floor(x).cast_mp_int();
}

//------------------------------------------------------
mp_int matcl::iceil(const mp_int& x)
{
    return x;
}
mp_int matcl::iceil(const mp_float& x)
{
    return ceil(x).cast_mp_int();
}
mp_int matcl::iceil(const mp_rational& x)
{
    return ceil(x).cast_mp_int();
}

//------------------------------------------------------
mp_int matcl::iround(const mp_int& x)
{
    return x;
}
mp_int matcl::iround(const mp_float& x)
{
    return round(x).cast_mp_int();
}
mp_int matcl::iround(const mp_rational& x)
{
    return round(x).cast_mp_int();
}

//------------------------------------------------------
mp_int matcl::itrunc(const mp_int& x)
{
    return x;
}
mp_int matcl::itrunc(const mp_float& x)
{
    return trunc(x).cast_mp_int();
}
mp_int matcl::itrunc(const mp_rational& x)
{
    return trunc(x).cast_mp_int();
}

//------------------------------------------------------
mp_int matcl::mpi_factorial(Integer i, precision p)
{
    (void)p;

    if (i < 0)
        return mp_int();

    mp_int ret;

    mpz_fac_ui(mmd::impl_value(ret).backend().data(), i);
    return ret;
};

mp_float matcl::mpf_factorial(Integer i, precision p)
{
    p       = details::correct_prec(p);

    if (i < 0)
        return constants::mp_nan(p);

    mp_float ret(p);
    mpfr_fac_ui(mmd::impl_value(ret), i, MPFR_RNDN);
    ret.update_debug();

    return ret;
}

template<> MATCL_MP_EXPORT
mp_int matcl::factorial<mp_int>(Integer i, precision p)
{
    return mpi_factorial(i, p);
}
template<> MATCL_MP_EXPORT
mp_float matcl::factorial<mp_float>(Integer i, precision p)
{
    return mpf_factorial(i, p);
}
template<> MATCL_MP_EXPORT
mp_complex matcl::factorial<mp_complex>(Integer i, precision p)
{
    return mpf_factorial(i, p);
}
template<> MATCL_MP_EXPORT
mp_rational matcl::factorial<mp_rational>(Integer i, precision p)
{
    return mpi_factorial(i, p);
}

//------------------------------------------------------
mp_int matcl::mpi_double_factorial(Integer i, precision p)
{
    (void)p;

    if (i < 0)
        return mp_int();

    mp_int ret;

    mpz_2fac_ui(mmd::impl_value(ret).backend().data(), i);
    return ret;
};
mp_float matcl::mpf_double_factorial(Integer i, precision p)
{
    if (i < 0)
        return constants::mp_nan(p);

    return mp_float(mpi_double_factorial(i,p),p);
};

template<> MATCL_MP_EXPORT
mp_int matcl::double_factorial<mp_int>(Integer i, precision p)
{
    return mpi_double_factorial(i, p);
}
template<> MATCL_MP_EXPORT
mp_float matcl::double_factorial<mp_float>(Integer i, precision p)
{
    return mpf_double_factorial(i, p);
}
template<> MATCL_MP_EXPORT
mp_complex matcl::double_factorial<mp_complex>(Integer i, precision p)
{
    return mpf_double_factorial(i, p);
}
template<> MATCL_MP_EXPORT
mp_rational matcl::double_factorial<mp_rational>(Integer i, precision p)
{
    return mpi_double_factorial(i, p);
}

//------------------------------------------------------
mp_int matcl::mpi_binomial_coefficient(Integer n, Integer k, precision p)
{
    (void)p;

    if (n < 0 || k < 0 || k > n)
        return mp_int();

    mp_int ret;

    mpz_bin_uiui(mmd::impl_value(ret).backend().data(), n, k);
    return ret;
};
mp_float matcl::mpf_binomial_coefficient(Integer n, Integer k, precision p)
{
    if (n < 0 || k < 0 || k > n)
        return constants::mp_nan(p);

    return mp_float(mpi_binomial_coefficient(n,k,p), p);
}

template<> MATCL_MP_EXPORT
mp_int matcl::binomial_coefficient<mp_int>(Integer n, Integer k, precision p)
{
    return mpi_binomial_coefficient(n, k, p);
}
template<> MATCL_MP_EXPORT
mp_float matcl::binomial_coefficient<mp_float>(Integer n, Integer k, precision p)
{
    return mpf_binomial_coefficient(n, k, p);
}
template<> MATCL_MP_EXPORT
mp_complex matcl::binomial_coefficient<mp_complex>(Integer n, Integer k, precision p)
{
    return mpf_binomial_coefficient(n, k, p);
}
template<> MATCL_MP_EXPORT
mp_rational matcl::binomial_coefficient<mp_rational>(Integer n, Integer k, precision p)
{
    return mpi_binomial_coefficient(n, k, p);
}

};
