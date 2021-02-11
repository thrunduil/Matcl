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

#include "matcl-mp/func_binary.h"
#include "utils/impl_types.h"
#include "utils/utils.h"
#include "matcl-mp/func_unary.h"
#include "utils/extend_precision.h"
#include "matcl-mp/constants.h"

#include <iostream>

namespace matcl { namespace mp { namespace details
{

namespace mmd = matcl::mp::details;

//-----------------------------------------------------------------------------
//                                  mp_float
//-----------------------------------------------------------------------------

template <>
mp_float mp_hypot_helper<mp_float, mp_float>::eval(const mp_float& x, const mp_float& y, 
                                                   precision req_prec)
{
    req_prec                = mmd::result_prec(req_prec, x.get_precision(), y.get_precision());
    mp_float ret(req_prec);

    using impl_type         = mmd::impl_float;
    const impl_type& vx     = *mmd::get_ptr<mp_float>::eval(x.m_data);
    const impl_type& vy     = *mmd::get_ptr<mp_float>::eval(y.m_data);
    impl_type& vr           = *mmd::get_ptr<mp_float>::eval(ret.m_data);

    mpfr_hypot(vr, vx, vy, MPFR_RNDN);
    ret.update_debug();
    return ret;
};

template <>
mp_float mp_atan2_helper<mp_float, mp_float>::eval(const mp_float& x, const mp_float& y,
                                                   precision req_prec)
{
    req_prec                = mmd::result_prec(req_prec, x.get_precision(), y.get_precision());
    mp_float ret(req_prec);

    using impl_type         = mmd::impl_float;
    const impl_type& vx     = *mmd::get_ptr<mp_float>::eval(x.m_data);
    const impl_type& vy     = *mmd::get_ptr<mp_float>::eval(y.m_data);
    impl_type& vr           = *mmd::get_ptr<mp_float>::eval(ret.m_data);

    mpfr_atan2(vr, vx, vy, MPFR_RNDN);
    ret.update_debug();
    return ret;
};

mp_float mp_nextafter_helper::eval(const mp_float& x, const mp_float& y)
{
    mp_float a(x);

    using impl_type     = mmd::impl_float;
    impl_type& ai       = *mmd::get_ptr<mp_float>::eval(a.m_data);
    const impl_type& yi = *mmd::get_ptr<mp_float>::eval(y.m_data);

    mpfr_nexttoward(ai, yi);
    a.update_debug();
    return a;
};

mp_float mp_float_distance_helper::eval(const mp_float& x, const mp_float& y, precision req_p)
{
    precision p1    = x.get_precision();
    precision p2    = y.get_precision();
    precision p     = req_p == 0 ? std::min(p1, p2) : req_p;

    if (is_nan(x) || is_nan(y))
        return constants::mp_nan(p);

    mp_float xp     = mp_float(x, p);
    mp_float yp     = mp_float(y, p);

    if (xp == yp)
        return mp_float(0.0, p);

    mp_float axp    = abs(xp);
    mp_float ayp    = abs(yp);
    mp_float min_xy = (axp < ayp) ? axp : ayp;
    mp_float dist   = abs(x - y) / eps(mp_float(min_xy, p));
    return dist;
};

mp_float mp_ulp_distance_helper::eval(const mp_float& x, const mp_float& y, precision req_p)
{
    precision p1    = x.get_precision();
    precision p2    = y.get_precision();
    precision p     = req_p == 0 ? std::min(p1, p2) : req_p;

    if (is_nan(x) || is_nan(y))
        return constants::mp_nan(p);

    if (x == y)
        return mp_float(0.0, p);

    mp_float dist   = abs(x - y) / eps(mp_float(x, p));
    return dist;
};

mp_float mp_copysign_helper::eval(const mp_float& x, const mp_float& y)
{
    mp_float ret(x.get_precision());

    using impl_type     = mmd::impl_float;
    const impl_type& xi = *mmd::get_ptr<mp_float>::eval(x.m_data);
    const impl_type& yi = *mmd::get_ptr<mp_float>::eval(y.m_data);
    impl_type& ri       = *mmd::get_ptr<mp_float>::eval(ret.m_data);

    mpfr_copysign(ri, xi, yi, MPFR_RNDN);
    ret.update_debug();
    return ret;
};

//-----------------------------------------------------------------------------
//                                  mp_complex
//-----------------------------------------------------------------------------
template <>
mp_float mp_hypot_helper<mp_complex, mp_complex>::eval(const mp_complex& x, const mp_complex& y,
                                                       precision req_prec)
{
    req_prec            = mmd::result_prec(req_prec, x.get_precision(), y.get_precision());
    precision int_prec  = mmd::extend_prec_hypot(req_prec);

    mp_float h1 = mp_hypot_helper<mp_float, mp_float>::eval(real(x),imag(x), int_prec); // 0.5 ulp
    mp_float h2 = mp_hypot_helper<mp_float, mp_float>::eval(real(y),imag(y), int_prec); // 0.5 ulp

    return mp_hypot_helper<mp_float, mp_float>::eval(h1,h2, req_prec);  // 1 ulp
};

template <>
mp_float mp_hypot_helper<mp_float, mp_complex>::eval(const mp_float& x, const mp_complex& y,
                                                     precision req_prec)
{
    req_prec    = mmd::result_prec(req_prec, x.get_precision(), y.get_precision());
    mp_float h2 = mp_hypot_helper<mp_float, mp_float>::eval(real(y),imag(y), req_prec);

    return mp_hypot_helper<mp_float, mp_float>::eval(x,h2, req_prec);
};

template <>
mp_float mp_hypot_helper<mp_complex, mp_float>::eval(const mp_complex& x, const mp_float & y,
                                                     precision req_prec)
{
    req_prec    = mmd::result_prec(req_prec, x.get_precision(), y.get_precision());
    mp_float h1 = mp_hypot_helper<mp_float, mp_float>::eval(real(x),imag(x), req_prec);

    return mp_hypot_helper<mp_float, mp_float>::eval(h1,y, req_prec);
};

//-------------------------------------------------------------------------
template<>
mp_int mp_mod_helper<mp_int, mp_int, mp_int>::eval(const mp_int& x, const mp_int& y, precision)
{
    if (is_zero(y) || is_zero(x))
        return mp_int(0);

    mp_int ret;
    mpz_mod(mmd::impl_value(ret).backend().data(), mmd::impl_value(x).backend().data(), 
            mmd::impl_value(y).backend().data());

    if (is_zero(ret))
        return ret;

    //mpz_mod does not care about signs; result must be negative if y is negative
    if (y < Integer(0))
        return ret + y;
    else
        return ret;
};
template<>
mp_rational mp_mod_helper<mp_rational, mp_rational, mp_rational>
        ::eval(const mp_rational& x, const mp_rational& y, precision)
{
    if (is_zero(y) || is_zero(x))
        return mp_int(0);

    //n is truncated
    mp_int n = (x/y).cast_mp_int();

    mp_rational ret = is_zero(n) ? x : x - n*y;

    if (is_zero(ret))
        return ret;

    //result must have the same sign as y
    if (y < Integer(0) && ret > Integer(0))
        return ret + y;
    else if (y > Integer(0) && ret < Integer(0))
        return ret + y;
    else
        return ret;
};

template<>
mp_float mp_mod_helper<mp_float, mp_float, mp_float>
        ::eval(const mp_float& x, const mp_float& y, precision req_prec)
{
    req_prec = mmd::result_prec(req_prec, x.get_precision(), y.get_precision());

    if (signbit(x) == signbit(y))
        return mp_rem_helper<mp_float, mp_float, mp_float>::eval(x,y, req_prec);

    if (is_zero(y) || is_zero(x))
        return copysign(mp_float(req_prec), y);

    if (is_inf(y) == true)
    {
        if (is_inf(x) == true)
            return constants::mp_nan(req_prec);

        //general formula will produce nan
        if (signbit(x) == signbit(y))
            return mp_float(x, req_prec);
        else
            return mp_float(y, req_prec);
    }

    //round to nearest integer
    Real lost_prec          = 2.0;
    precision int_prec      = mmd::extend_prec_dynamic(req_prec, lost_prec);

    mp_float ret(int_prec);
    mpfr_remainder(mmd::impl_value(ret), mmd::impl_value(x), mmd::impl_value(y), MPFR_RNDN);
    ret.update_debug();

    if (is_zero(ret))
        return copysign(mp_float(ret, req_prec), y);

    if (signbit(ret) == signbit(y))
        return mp_float(ret, req_prec);

    //ret and y has different signs, but |ret + y| >= 1/2|y|
    //error = 0.5 + 0.5 * 2.0 * |ret|/|y| <= 1.5
    return plus(ret, y, req_prec);
};

//-------------------------------------------------------------------------
template<>
mp_int mp_rem_helper<mp_int, mp_int, mp_int>::eval(const mp_int& x, const mp_int& y, precision)
{
    if (is_zero(y) || is_zero(x))
        return mp_int(0);

    mp_int ret;
    mpz_mod(mmd::impl_value(ret).backend().data(), mmd::impl_value(x).backend().data(), 
            mmd::impl_value(y).backend().data());

    if (is_zero(ret))
        return ret;

    //mpz_mod does not care about signs; result must be negative if x is negative
    if (x < Integer(0))
        return (y < 0)? ret + y : ret - y;
    else
        return ret;
};
template<>
mp_rational mp_rem_helper<mp_rational, mp_rational, mp_rational>
        ::eval(const mp_rational& x, const mp_rational& y, precision)
{
    if (is_zero(y) || is_zero(x))
        return mp_int(0);

    //n is truncated
    mp_int n = (x/y).cast_mp_int();

    mp_rational ret = x - n*y;
    return ret;
};
template<>
mp_float mp_rem_helper<mp_float, mp_float, mp_float>
        ::eval(const mp_float& x, const mp_float& y, precision req_prec)
{
    req_prec = mmd::result_prec(req_prec, x.get_precision(), y.get_precision());

    mp_float ret(req_prec);

    if (is_zero(y) || is_zero(x))
        return copysign(mp_float(req_prec), x);

    //our definition of rem is the same as mpfr definition of mod    
    mpfr_fmod(mmd::impl_value(ret), mmd::impl_value(x), mmd::impl_value(y), MPFR_RNDN);
    ret.update_debug();

    return copysign(ret, x);
};

//-------------------------------------------------------------------------
template<>
mp_float mp_fdim_helper<mp_float, mp_float, mp_float>
        ::eval(const mp_float& x, const mp_float& y, precision req_prec)
{
    req_prec = mmd::result_prec(req_prec, x.get_precision(), y.get_precision());

    mp_float ret(req_prec);
    mpfr_dim(mmd::impl_value(ret), mmd::impl_value(x), mmd::impl_value(y), MPFR_RNDN);
    ret.update_debug();

    return ret;
};

}};}