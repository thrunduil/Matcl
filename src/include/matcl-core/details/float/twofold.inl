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
// from 'Twofold fast arithmetic', E. Latkin, 2014.
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
    return twofold_plus_sorted(val, err);
}

template<class Float_type>
force_inline twofold<Float_type>
twofold<Float_type>::normalize(const Float_type& val, const Float_type& err)
{
    return twofold_plus(val, err);
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
matcl::twofold_plus_sorted(const Float_type& a, const Float_type& b)
{
    Float_type val  = a + b;
    Float_type b2   = val - a;
    Float_type err  = b - b2;

    return twofold<Float_type>(val, err);
};

template<class Float_type>
force_inline twofold<Float_type> 
matcl::twofold_plus(const Float_type& a, const Float_type& b)
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
    twofold<Float_type> z0  = twofold_plus(a.value, b.value);
    twofold<Float_type> z1  = twofold_plus(a.error, b.error);

    z1          = z0.error + z1;
    z1          = z0.value + z1;

    return z1;
};

template<class Float_type>
force_inline twofold<Float_type> 
matcl::operator+(const twofold<Float_type>& a, const Float_type& b)
{
    twofold<Float_type> z0  = twofold_plus(a.value, b);
    return twofold<Float_type>::normalize_fast(z0.value, a.error + z0.error);
}

template<class Float_type>
force_inline twofold<Float_type> 
matcl::operator+(const Float_type& a, const twofold<Float_type>& b)
{
    twofold<Float_type> z0  = twofold_plus(a, b.value);
    return twofold<Float_type>::normalize_fast(z0.value, b.error + z0.error);
}

template<class Float_type>
force_inline twofold<Float_type> 
matcl::operator-(const twofold<Float_type>& a, const twofold<Float_type>& b)
{
    twofold<Float_type> z0  = twofold_minus(a.value, b.value);
    twofold<Float_type> z1  = twofold_minus(a.error, b.error);

    z1          = z0.error + z1;
    z1          = z0.value + z1;

    return z1;
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
    twofold<Float_type> p00 = twofold_mult(a.value, b.value);

    Float_type p01  = a.value * b.error;
    Float_type p10  = a.error * b.value;
    Float_type p11  = a.error * b.error;

    Float_type val  = p00.value;
    Float_type err1 = p00.error + p11;
    Float_type err2 = p01 + p10;
    Float_type err  = err1 + err2;

    return twofold<Float_type>::normalize_fast(val, err);
}

template<class Float_type>
force_inline twofold<Float_type>
matcl::operator*(const twofold<Float_type>& a, const Float_type& b)
{
    namespace mrd   = matcl::raw::details;

    twofold<Float_type> p00 = twofold_mult(a.value, b);

    Float_type p10  = a.error * b;
    Float_type val  = p00.value;
    Float_type err  = p00.error + p10;

    return twofold<Float_type>::normalize_fast(val, err);
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
    namespace mrd   = matcl::raw::details;

    Float_type z0   = a.value / b.value;
    Float_type r0   = - details::func_fms_a<Float_type>::eval(z0, b.value, a.value);
    Float_type r1   = - details::func_fms_a<Float_type>::eval(z0, b.error, a.error);
    Float_type c0   = r0 + r1;
    Float_type d0   = b.value + b.error;
    Float_type z1   = c0 / d0;

    return twofold<Float_type>::normalize_fast(z0, z1);
}

template<class Float_type>
force_inline twofold<Float_type> 
matcl::operator/(const twofold<Float_type>& a, const Float_type& b)
{
    namespace mrd   = matcl::raw::details;

    Float_type z0   = a.value / b;
    Float_type r    = - details::func_fms_a<Float_type>::eval(z0, b, a.value);
    Float_type c    = r + a.error;
    Float_type z1   = c / b;

    return twofold<Float_type>::normalize_fast(z0, z1);
}

template<class Float_type>
force_inline twofold<Float_type> 
matcl::operator/(const Float_type& a, const twofold<Float_type>& b)
{
    namespace mrd   = matcl::raw::details;

    Float_type z0   = a / b.value;
    Float_type r0   = - details::func_fms_a<Float_type>::eval(z0, b.value, a);
    Float_type r1   = - z0 * b.error;
    Float_type c0   = r0 + r1;
    Float_type d0   = b.value + b.error;
    Float_type z1   = c0 / d0;

    return twofold<Float_type>::normalize_fast(z0, z1);
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