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

#include "matcl-mp/func_binary.h"
#include "matcl-mp/details/utils.h"
#include "matcl-mp/details/prec_utils.h"

namespace matcl { namespace mp { namespace details
{

template <class T1, class T2>
struct MATCL_MP_EXPORT mp_hypot_helper
{
    static mp_float eval(const T1& x, const T2& y, precision req_prec);
};

template <class T1, class T2>
struct MATCL_MP_EXPORT mp_atan2_helper
{
    static mp_float eval(const T1& x, const T2& y, precision req_prec);
};

struct MATCL_MP_EXPORT mp_nextafter_helper
{
    static mp_float eval(const mp_float& x, const mp_float& y);
};

struct MATCL_MP_EXPORT mp_float_distance_helper
{
    static mp_float eval(const mp_float& x, const mp_float& y, precision p);
};

struct MATCL_MP_EXPORT mp_ulp_distance_helper
{
    static mp_float eval(const mp_float& x, const mp_float& y, precision p);
};

struct MATCL_MP_EXPORT mp_copysign_helper
{
    static mp_float eval(const mp_float& x, const mp_float& y);
};

template<class Ret>
struct make_result
{
    static Ret eval_0(precision p)  { return Ret(p); };
    static Ret eval_1(precision p)  { return Ret(1, p); };
};

template<>
struct make_result<mp_rational>
{
    using Ret = mp_rational;

    static Ret eval_0(precision)    { return Ret(); };
    static Ret eval_1(precision)    { return Ret(1); };
};

template<class Ret, class T1, class T2>
struct MATCL_MP_EXPORT mp_mod_helper
{
    static Ret eval(const T1& x, const T2& y, precision p);
};

template<class Ret, class T1, class T2>
struct MATCL_MP_EXPORT mp_rem_helper
{
    static Ret eval(const T1& x, const T2& y, precision p);
};

template<class Ret, class T1, class T2>
struct MATCL_MP_EXPORT mp_fdim_helper
{
    static Ret eval(const T1& x, const T2& y, precision p);
};

}}};

namespace matcl
{

namespace md = matcl :: details;

template<class T1, class T2, class Enable>
mp_float matcl::hypot(const T1& x, const T2& y, precision req_prec)
{
    namespace mmd = matcl :: mp :: details;
    using T1P   = typename md::unify_types<T1, mp_float>::type;
    using T2P   = typename md::unify_types<T2, mp_float>::type;

    return mmd::mp_hypot_helper<T1P, T2P>::eval(T1P(x), T2P(y), req_prec);
}

template<class T1, class T2, class Enable>
mp_float matcl::atan2(const T1& x, const T2& y, precision req_prec)
{
    static_assert(!md::is_complex<T1>::value && !md::is_complex<T2>::value, 
                  "not available for complex");

    namespace mmd = matcl :: mp :: details;
    using T1P   = typename md::unify_types<T1, mp_float>::type;
    using T2P   = typename md::unify_types<T2, mp_float>::type;

    return mmd::mp_atan2_helper<T1P, T2P>::eval(T1P(x), T2P(y), req_prec);
}

template<class T1, class T2, class Enable>
bool matcl::operator==(const T1& x, const T2& y)
{
    namespace mmd = matcl :: mp :: details;
    return mmd::eeq(x, y);
}

template<class T1, class T2, class Enable>
bool matcl::eeq_nan(const T1& x, const T2& y)
{
    namespace mmd = matcl :: mp :: details;
    return mmd::eeq_nan(x, y);
}

template<class T1, class T2, class Enable>
bool matcl::operator!=(const T1& x, const T2& y)
{
    namespace mmd = matcl :: mp :: details;
    return mmd::neq(x, y);
}

template<class T1, class T2, class Enable>
bool matcl::neq_nan(const T1& x, const T2& y)
{
    namespace mmd = matcl :: mp :: details;
    return mmd::neq_nan(x, y);
}

template<class T1, class T2, class Enable>
bool matcl::operator>=(const T1& x, const T2& y)
{
    namespace mmd = matcl :: mp :: details;
    return mmd::geq(x, y);
}

template<class T1, class T2, class Enable>
bool matcl::operator<=(const T1& x, const T2& y)
{
    namespace mmd = matcl :: mp :: details;
    return mmd::leq(x, y);
}

template<class T1, class T2, class Enable>
bool matcl::operator<(const T1& x, const T2& y)
{
    namespace mmd = matcl :: mp :: details;
    return mmd::lt(x, y);
}

template<class T1, class T2, class Enable>
bool matcl::operator>(const T1& x, const T2& y)
{
    namespace mmd = matcl :: mp :: details;
    return mmd::gt(x, y);
}

template<class T1, class T2, class Enable, class Result>
Result matcl::operator+(const T1& x, const T2& y)
{
    namespace mmd = matcl :: mp :: details;
    return mmd::plus_impl(x, y, precision());
}

template<class T1, class T2, class Enable, class Result>
Result matcl::operator-(const T1& x, const T2& y)
{
    namespace mmd = matcl :: mp :: details;
    return mmd::minus_impl(x, y, precision());
}

template<class T1, class T2, class Enable, class Result>
Result matcl::operator*(const T1& x, const T2& y)
{
    namespace mmd = matcl :: mp :: details;
    return mmd::mul_impl(x, y, precision());
}

template<class T1, class T2, class Enable, class Result>
Result matcl::operator/(const T1& x, const T2& y)
{
    namespace mmd = matcl :: mp :: details;
    return mmd::div_impl(x, y, precision());
}

template<class T1, class T2, class Enable, class Result>
Result matcl::idiv(const T1& x, const T2& y, precision req_prec)
{
    namespace mmd = matcl :: mp :: details;
    return mmd::idiv_impl(x, y, req_prec);
}

template<class T1, class T2, class Enable, class Result>
Result matcl::plus(const T1& x, const T2& y, precision req_prec)
{
    namespace mmd = matcl :: mp :: details;
    return mmd::plus_impl(x, y, req_prec);
}

template<class T1, class T2, class Enable, class Result>
Result matcl::minus(const T1& x, const T2& y, precision req_prec)
{
    namespace mmd = matcl :: mp :: details;
    return mmd::minus_impl(x, y, req_prec);
}

template<class T1, class T2, class Enable, class Result>
Result matcl::mul(const T1& x, const T2& y, precision req_prec)
{
    namespace mmd = matcl :: mp :: details;
    return mmd::mul_impl(x, y, req_prec);
}

template<class T1, class T2, class Enable, class Result>
Result matcl::div(const T1& x, const T2& y, precision req_prec)
{
    namespace mmd = matcl :: mp :: details;
    return mmd::div_impl(x, y, req_prec);
}

template<class T1, class T2, class Enable, class Result>
Result matcl::div_0(const T1& x, const T2& y, precision req_prec)
{
    namespace mmd = matcl :: mp :: details;
    req_prec  = mmd::result_prec(req_prec, mmd::get_precision(x), mmd::get_precision(y));

    if (is_zero(x) && is_zero(y))
        return mmd::make_result<Result>::eval_0(req_prec);
    else
        return mmd::div_impl(x, y, req_prec);
};

template<class T1, class T2, class Enable, class Result>
Result matcl::div_1(const T1& x, const T2& y, precision req_prec)
{
    namespace mmd = matcl :: mp :: details;
    req_prec  = mmd::result_prec(req_prec, mmd::get_precision(x), mmd::get_precision(y));

    if (is_zero(x) && is_zero(y))
        return mmd::make_result<Result>::eval_1(req_prec);
    else
        return mmd::div_impl(x, y, req_prec);
};

template<class T1, class T2, class Enable, class Result>
Result matcl::pow(const T1& x, const T2& y, precision req_prec)
{
    namespace mmd = matcl :: mp :: details;
    using T1P   = typename mmd::promote_floats<T1>::type;
    using T2P   = typename mmd::promote_floats<T2>::type;
    return mmd::pow_impl(T1P(x), T2P(y), req_prec);
}

template<class T1, class T2, class Enable>
mp_float matcl::nextafter(const T1& x, const T2& y)
{
    static_assert(!md::is_complex<T1>::value && !md::is_complex<T2>::value, 
                  "not available for complex");

    namespace mmd = matcl :: mp :: details;
    mp_float xp(x);
    mp_float yp(y);

    return mmd::mp_nextafter_helper::eval(xp,yp);
};

template<class T1, class T2, class Enable>
mp_float matcl::float_distance(const T1& x, const T2& y, precision p)
{
    static_assert(!md::is_complex<T1>::value && !md::is_complex<T2>::value, 
                  "not available for complex");

    namespace mmd = matcl :: mp :: details;
    mp_float xp(x);
    mp_float yp(y);

    return mmd::mp_float_distance_helper::eval(xp, yp, p);
};

template<class T1, class T2, class Enable>
mp_float matcl::ulp_distance(const T1& x, const T2& y, precision p)
{
    static_assert(!md::is_complex<T1>::value && !md::is_complex<T2>::value, 
                  "not available for complex");

    namespace mmd = matcl :: mp :: details;
    mp_float xp(x);
    mp_float yp(y);

    return mmd::mp_ulp_distance_helper::eval(xp, yp, p);
};

template<class T1, class T2, class Enable>
mp_float matcl::copysign(const T1& x, const T2& y)
{
    namespace mmd = matcl :: mp :: details;

    static_assert(!md::is_complex<T1>::value && !md::is_complex<T2>::value, 
                  "not available for complex");

    mp_float xp(x);
    mp_float yp(y);

    return mmd::mp_copysign_helper::eval(xp,yp);
};

template<class T1, class T2, class Enable, class Ret>
Ret matcl::mod(const T1& x, const T2& y, precision req_prec)
{
    static_assert(!md::is_complex<T1>::value && !md::is_complex<T2>::value, 
                  "not available for complex");

    using Ret_r = typename md::real_type<Ret>::type;
    using T1P   = typename md::unify_types<T1, Ret_r>::type;
    using T2P   = typename md::unify_types<T2, Ret_r>::type;

    namespace mmd = matcl :: mp :: details;

    T1P xp(x);
    T2P yp(y);

    return mmd::mp_mod_helper<Ret, T1P, T2P>::eval(xp,yp, req_prec);
};

template<class T1, class T2, class Enable, class Ret>
Ret matcl::rem(const T1& x, const T2& y, precision req_prec)
{
    static_assert(!md::is_complex<T1>::value && !md::is_complex<T2>::value, 
                  "not available for complex");

    using T1P   = typename md::unify_types<T1, Ret>::type;
    using T2P   = typename md::unify_types<T2, Ret>::type;

    namespace mmd = matcl :: mp :: details;

    T1P xp(x);
    T2P yp(y);

    return mmd::mp_rem_helper<Ret, T1P, T2P>::eval(xp,yp, req_prec);
};

template<class T1, class T2, class Enable, class Ret>
Ret matcl::fdim(const T1& x, const T2& y, precision req_prec)
{
    static_assert(!md::is_complex<T1>::value && !md::is_complex<T2>::value, 
                  "not available for complex");

    using T1P   = typename md::unify_types<T1, Ret>::type;
    using T2P   = typename md::unify_types<T2, Ret>::type;

    namespace mmd = matcl :: mp :: details;

    T1P xp(x);
    T2P yp(y);

    return mmd::mp_fdim_helper<Ret, T1P, T2P>::eval(xp,yp, req_prec);
};

};
