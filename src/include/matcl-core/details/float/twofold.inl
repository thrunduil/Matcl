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

#include "matcl-core/float/twofold.h"
#include "matcl-core/details/scalfunc_real.h"
#include "matcl-core/details/float/twofold_base_functions.h"
#include "matcl-core/details/float/helpers.h"

// most of algorithms are modified versions of algorithms
// from:
//  [1]  "Twofold fast arithmetic", E. Latkin, 2014.
//  [2]. "Tight and rigourous error bounds for basic building blocks of 
//       double-word arithmetic", M. Joldes, J.M. Muller, V. Popescu:
namespace matcl
{

template<class Float_type>
force_inline
twofold<Float_type>::twofold()
    : value(0.0), error(Float_type(0))
{};

template<class Float_type>
force_inline
twofold<Float_type>::twofold(uninitialized)
{};

template<class Float_type>
force_inline
twofold<Float_type>::twofold(const Float_type& val)
    : value(val), error(Float_type(0))
{};

template<class Float_type>
force_inline
twofold<Float_type>::twofold(const Float_type& val, const Float_type& err)
    : value(val), error(err)
{}

template<class Float_type>
force_inline twofold<Float_type> 
twofold<Float_type>::normalize_fast(const Float_type& val, const Float_type& err)
{
    return twofold_sum_sorted(val, err);
}

template<class Float_type>
force_inline twofold<Float_type>
twofold<Float_type>::normalize(const Float_type& val, const Float_type& err)
{
    return twofold_sum(val, err);
}

template<class Float_type>
force_inline
const Float_type& twofold<Float_type>::sum() const
{
    // twofold is normalized and value + error = value
    return value;
}

template<class Float_type>
force_inline
bool twofold<Float_type>::is_finite() const
{
    namespace mrds  = matcl::raw::details::scal_func;

    bool v1 = mrds::finite(value);
    bool v2 = mrds::finite(error);
    return v1 && v2;
}

template<class Float_type>
force_inline twofold<Float_type> 
matcl::twofold_sum_sorted(const Float_type& a, const Float_type& b)
{
    Float_type val  = a + b;
    Float_type b2   = val - a;
    Float_type err  = b - b2;

    return twofold<Float_type>(val, err);
};

template<class Float_type>
force_inline twofold<Float_type> 
matcl::twofold_sum(const Float_type& a, const Float_type& b)
{
    Float_type val  = a + b;
    Float_type b2   = val - a;
    Float_type a2   = val - b2;
    Float_type b3   = b - b2;
    Float_type a3   = a - a2;
    Float_type err  = a3 + b3;

    return twofold<Float_type>(val, err);
};

template<class Float_type>
force_inline twofold<Float_type> 
matcl::twofold_minus_sorted(const Float_type& a, const Float_type& b)
{
    Float_type val  = a - b;
    Float_type b2   = a - val;
    Float_type err  = b2 - b;

    return twofold<Float_type>(val, err);
};

template<class Float_type>
force_inline twofold<Float_type> 
matcl::twofold_minus(const Float_type& a, const Float_type& b)
{
    Float_type val  = a - b;
    Float_type b2   = val - a;
    Float_type a2   = val - b2;
    Float_type b3   = b + b2;
    Float_type a3   = a - a2;
    Float_type err  = a3 - b3;

    return twofold<Float_type>(val, err);
};

template<class Float_type>
force_inline twofold<Float_type> 
matcl::twofold_mult(const Float_type& a, const Float_type& b)
{
    Float_type val  = a * b;
    Float_type err  = details::func_fms_a<Float_type>::eval(a, b, val);
        
    return twofold<Float_type>(val, err);
};

template<class Float_type>
force_inline twofold<Float_type> 
matcl::twofold_div(const Float_type& a, const Float_type& b)
{
    namespace mrd   = matcl::raw::details;

    Float_type val  = a / b;
    Float_type err  = - details::func_fms_a<Float_type>::eval(val, b, a);
    err             = err / b;

    return twofold<Float_type>(val, err);
};

template<class Float_type>
force_inline twofold<Float_type> 
matcl::twofold_sqrt(const Float_type& a)
{
    namespace mrd   = matcl::raw::details;
    Float_type z0   = details::func_sqrt<Float_type>::eval(a);
    Float_type t2   = -details::func_fms_a<Float_type>::eval(z0, z0, a);
    Float_type d    = Float_type(2.0) * z0;
    Float_type z1   = t2 / d;

    return twofold<Float_type>(z0, z1);
};

template<class Float_type>
force_inline twofold<Float_type> 
matcl::operator+(const twofold<Float_type>& a, const twofold<Float_type>& b)
{
    // relative error is bounded by 3u^2/(1-4u) ~ 3u^2
    // this bound need not be optimal
    // see "Tight and rigourous error bounds for basic building blocks of 
    // double-word arithmetic", M. Joldes, J.M. Muller, V. Popescu

    twofold<Float_type> z0  = twofold_sum(a.value, b.value);
    twofold<Float_type> z1  = twofold_sum(a.error, b.error);

    Float_type c            = z0.error + z1.value;
    twofold<Float_type> v   = twofold_sum_sorted(z0.value, c);

    Float_type w            = z1.error + v.error;
    twofold<Float_type> r   = twofold_sum_sorted(v.value, w);
    return r;
};

template<class Float_type>
force_inline twofold<Float_type> 
matcl::operator+(const twofold<Float_type>& a, const Float_type& b)
{
    twofold<Float_type> z0  = twofold_sum(a.value, b);
    return twofold<Float_type>::normalize_fast(z0.value, a.error + z0.error);
}

template<class Float_type>
force_inline twofold<Float_type> 
matcl::operator+(const Float_type& a, const twofold<Float_type>& b)
{
    twofold<Float_type> z0  = twofold_sum(a, b.value);
    return twofold<Float_type>::normalize_fast(z0.value, b.error + z0.error);
}

template<class Float_type>
force_inline twofold<Float_type> 
matcl::operator-(const twofold<Float_type>& a, const twofold<Float_type>& b)
{
    twofold<Float_type> z0  = twofold_minus(a.value, b.value);
    twofold<Float_type> z1  = twofold_minus(a.error, b.error);

    Float_type c            = z0.error + z1.value;
    twofold<Float_type> v   = twofold_sum_sorted(z0.value, c);

    Float_type w            = z1.error + v.error;
    twofold<Float_type> r   = twofold_sum_sorted(v.value, w);
    return r;
};

template<class Float_type>
force_inline twofold<Float_type> 
matcl::operator-(const twofold<Float_type>& a, const Float_type& b)
{
    twofold<Float_type> z0  = twofold_minus(a.value, b);
    return twofold<Float_type>::normalize_fast(z0.value, a.error + z0.error);
};

template<class Float_type>
force_inline twofold<Float_type> 
matcl::operator-(const Float_type& a, const twofold<Float_type>& b)
{
    twofold<Float_type> z0  = twofold_minus(a, b.value);
    return twofold<Float_type>::normalize_fast(z0.value, z0.error - b.error);
};

template<class Float_type>
force_inline twofold<Float_type> 
matcl::operator-(const twofold<Float_type>& a)
{
    return twofold<Float_type>(-a.value, -a.error);
};

template<class Float_type>
force_inline twofold<Float_type> 
matcl::operator*(const twofold<Float_type>& a, const twofold<Float_type>& b)
{
    // "Tight and rigourous error bounds for basic building blocks of 
    // double-word arithmetic", M. Joldes, J.M. Muller, V. Popescu:
    // relative forward error is 6 u^2 if FMA is available, 7 u^2 otherwise,
    // error can be reduced further to 5 u^2 at the cost of one multiplication
    // and one addition
    twofold<Float_type> p00 = twofold_mult(a.value, b.value);

    Float_type p01  = a.value * b.error;
    Float_type err  = fma_f(a.error, b.value, p01);
    err             = err + p00.error;
    return twofold<Float_type>::normalize_fast(p00.value, err);
}

template<class Float_type>
force_inline twofold<Float_type>
matcl::operator*(const twofold<Float_type>& a, const Float_type& b)
{
    namespace mrd   = matcl::raw::details;

    twofold<Float_type> p00 = twofold_mult(a.value, b);

    // "Tight and rigourous error bounds for basic building blocks of 
    // double-word arithmetic", M. Joldes, J.M. Muller, V. Popescu:
    // relative forward error is 2 u^2 if FMA is available, 3 u^2 otherwise, 
    Float_type val          = p00.value;
    Float_type err          = fma_f(a.error, b, p00.error);
    return twofold<Float_type>::normalize_fast(val, err);

    // more accurate version of this function, relative error is 3/2 u^2
    //twofold<Float_type> p00 = twofold_mult(a.value, b);
    //twofold<Float_type> t   = twofold_sum_sorted(p00.value, a.error * b);
    //twofold<Float_type> z   = twofold_sum_sorted(t.value, t.error + p00.error);
    //return z;
}

template<class Float_type>
force_inline twofold<Float_type> 
matcl::operator*(const Float_type& a, const twofold<Float_type>& b)
{
    return b * a;
}

template<class Float_type>
force_inline twofold<Float_type>
matcl::operator/(const twofold<Float_type>& a, const twofold<Float_type>& b)
{
    // "Tight and rigourous error bounds for basic building blocks of 
    // double-word arithmetic", M. Joldes, J.M. Muller, V. Popescu:
    // relative forward error is 15 u^2

    using twofold_t = matcl::twofold<Float_type>;
    namespace mrd   = matcl::raw::details;

    Float_type val  = a.value / b.value;
    twofold_t r     = b * val;
    Float_type pi   = a.value - r.value;
    Float_type d_l  = a.error - r.error;
    Float_type d    = pi + d_l;
    Float_type err  = d / b.value;

    return twofold<Float_type>::normalize_fast(val, err);
}

template<class Float_type>
force_inline twofold<Float_type> 
matcl::operator/(const twofold<Float_type>& a, const Float_type& b)
{
    // "Tight and rigourous error bounds for basic building blocks of 
    // double-word arithmetic", M. Joldes, J.M. Muller, V. Popescu:
    // relative forward error is 3 u^2.

    namespace mrd   = matcl::raw::details;
    using twofold_t = twofold<Float_type>;

    Float_type val  = a.value / b;
    twofold_t pi    = twofold_mult(val, b);
    Float_type dh   = (a.value - pi.value);
    Float_type dt   = (dh - pi.error);
    Float_type d    = (dt + a.error);
    Float_type err  = d / b;
    return twofold<Float_type>::normalize_fast(val, err);
}

template<class Float_type>
force_inline twofold<Float_type> 
matcl::operator/(const Float_type& a, const twofold<Float_type>& b)
{
    // "Tight and rigourous error bounds for basic building blocks of 
    // double-word arithmetic", M. Joldes, J.M. Muller, V. Popescu:
    // relative forward error is 15 u^2

    using twofold_t = matcl::twofold<Float_type>;
    namespace mrd   = matcl::raw::details;

    Float_type val  = a / b.value;
    twofold_t r     = b * val;
    Float_type pi   = a - r.value;
    Float_type d    = pi - r.error;
    Float_type err  = d / b.value;

    return twofold<Float_type>::normalize_fast(val, err);
}

template<class Float_type>
force_inline twofold<Float_type> 
matcl::sqrt(const twofold<Float_type>& a)
{
    namespace mrd   = matcl::raw::details;

    Float_type z0   = details::func_sqrt<Float_type>::eval(a.value);
    Float_type t1   = details::func_fms_a<Float_type>::eval(z0, z0, a.value);
    Float_type t2   = a.error - t1;
    Float_type d    = Float_type(2.0) * z0;
    Float_type z1   = t2 / d;

    return twofold<Float_type>::normalize_fast(z0, z1);
}

template<class Float_type>
force_inline twofold<Float_type>
matcl::abs(const twofold<Float_type>& a)
{
    Float_type v, e;
    details::func_abs<Float_type>::eval(a.value, a.error, v, e);

    return twofold<Float_type>(v,e);
}

//-----------------------------------------------------------------------
//                      ERROR RELATED FUNCTIONS
//-----------------------------------------------------------------------
template<class Float_type>
Float_type matcl::eps(const twofold<Float_type>& x)
{
    static const bool is_simd   = details::is_simd_type<Float_type>::value;
    return details::func_eps<Float_type, is_simd>::eval(x);
};

template<class Float_type>
Float_type matcl::float_distance(const twofold<Float_type>& x, const twofold<Float_type>& y)
{
    static const bool is_simd   = details::is_simd_type<Float_type>::value;
    return details::func_float_distance<Float_type, is_simd>::eval(x, y);
};

//-----------------------------------------------------------------------
//                      IO FUNCTIONS
//-----------------------------------------------------------------------
template<class Float_type>
std::ostream& matcl::operator<<(std::ostream& os, const twofold<Float_type>& x)
{
    static const bool is_simd   = details::is_simd_type<Float_type>::value;
    details::func_save<Float_type, is_simd>::eval(os, x);

    return os;
};

template<class Float_type>
std::istream& matcl::operator>>(std::istream& is, twofold<Float_type>& x)
{
    static const bool is_simd   = details::is_simd_type<Float_type>::value;
    details::func_load<Float_type, is_simd>::eval(is, x);

    return is;
};

}