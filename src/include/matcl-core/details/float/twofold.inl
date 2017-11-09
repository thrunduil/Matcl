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

// most of algorithms are modified versions of algorithms
// from 'Twofold fast arithmetic', E. Latkin, 2014.
namespace matcl
{

force_inline
twofold::twofold()
    : value(0.0), error(0.0)
{};

force_inline
twofold::twofold(const double& val)
    : value(val), error(0.0)
{};

force_inline
twofold::twofold(const double& val, const double& err)
    : value(val), error(err)
{
}

force_inline
twofold twofold::normalize_fast(const double& val, const double& err)
{
    return twofold_plus_sorted(val, err);
}

force_inline
twofold twofold::normalize(const double& val, const double& err)
{
    return twofold_plus(val, err);
}

force_inline
const double& twofold::sum() const
{
    // twofold is normalized and value + error = value
    return value;
}

force_inline
bool twofold::is_finite() const
{
    namespace mrds  = matcl::raw::details::scal_func;

    bool v1 = mrds::finite(value);
    bool v2 = mrds::finite(error);
    return v1 && v2;
}

force_inline
twofold matcl::twofold_plus_sorted(double a, double b)
{
    double val  = a + b;
    double b2   = val - a;
    double err  = b - b2;

    return twofold(val, err);
};

force_inline
twofold matcl::twofold_plus(double a, double b)
{
    double val  = a + b;
    double b2   = val - a;
    double a2   = val - b2;
    double b3   = b - b2;
    double a3   = a - a2;
    double err  = a3 + b3;

    return twofold(val, err);
};

force_inline
twofold matcl::twofold_minus_sorted(double a, double b)
{
    double val  = a - b;
    double b2   = a - val;
    double err  = b2 - b;

    return twofold(val, err);
};

force_inline
twofold matcl::twofold_minus(double a, double b)
{
    double val  = a - b;
    double b2   = val - a;
    double a2   = val - b2;
    double b3   = b + b2;
    double a3   = a - a2;
    double err  = a3 - b3;

    return twofold(val, err);
};

force_inline
twofold matcl::twofold_mult(double a, double b)
{
    #if MATCL_ARCHITECTURE_HAS_FMA

        namespace mrd   = matcl::raw::details;

        double val  = a * b;
        double err  = mrd::scal_func::fms_a(a, b, val);
        
    #else

        return twofold_mult_dekker(a, b);

    #endif

    return twofold(val, err);
};

force_inline
twofold matcl::twofold_div(double a, double b)
{
    namespace mrd   = matcl::raw::details;

    double val      = a / b;
    double err      = - mrd::scal_func::fms_a(val, b, a);
    err             = err / b;

    return twofold(val, err);
};

force_inline
twofold matcl::twofold_sqrt(double a)
{
    namespace mrd   = matcl::raw::details;
    double z0       = mrd::scal_func::sqrt(a);
    double t2       = -mrd::scal_func::fms_a(z0, z0, a);
    double d        = 2.0 * z0;
    double z1       = t2 / d;

    return twofold(z0, z1);
};

force_inline
twofold matcl::operator+(const twofold& a, const twofold& b)
{
    twofold z0  = twofold_plus(a.value, b.value);
    twofold z1  = twofold_plus(a.error, b.error);
    z1          = z0.error + z1;
    z1          = z0.value + z1;

    return z1;
};

force_inline
twofold matcl::operator+(const twofold& a, double b)
{
    twofold z0  = twofold_plus(a.value, b);
    return twofold::normalize_fast(z0.value, a.error + z0.error);
}

force_inline
twofold matcl::operator+(double a, const twofold& b)
{
    twofold z0  = twofold_plus(a, b.value);
    return twofold::normalize_fast(z0.value, b.error + z0.error);
}

force_inline
twofold matcl::operator-(const twofold& a, const twofold& b)
{
    twofold z0  = twofold_minus(a.value, b.value);
    twofold z1  = twofold_minus(a.error, b.error);
    z1          = z0.error + z1;
    z1          = z0.value + z1;

    return z1;
};

force_inline
twofold matcl::operator-(const twofold& a, double b)
{
    twofold z0  = twofold_minus(a.value, b);
    return twofold::normalize_fast(z0.value, a.error + z0.error);
};

force_inline
twofold matcl::operator-(double a, const twofold& b)
{
    twofold z0  = twofold_minus(a, b.value);
    return twofold::normalize_fast(z0.value, z0.error - b.error);
};

force_inline
twofold matcl::operator-(const twofold& a)
{
    return twofold(-a.value, -a.error);
};

force_inline
twofold matcl::operator*(const twofold& a, const twofold& b)
{
    twofold p00 = twofold_mult(a.value, b.value);
    double p01  = a.value * b.error;
    double p10  = a.error * b.value;
    double p11  = a.error * b.error;

    double val  = p00.value;
    double err1 = p00.error + p11;
    double err2 = p01 + p10;
    double err  = err1 + err2;

    return twofold::normalize_fast(val, err);
}

force_inline
twofold matcl::operator*(const twofold& a, double b)
{
    namespace mrd   = matcl::raw::details;

    twofold p00 = twofold_mult(a.value, b);
    double p10  = a.error * b;

    double val  = p00.value;
    double err  = p00.error + p10;

    return twofold::normalize_fast(val, err);
}

force_inline
twofold matcl::operator*(double a, const twofold& b)
{
    return b * a;
}

force_inline
twofold matcl::operator/(const twofold& a, const twofold& b)
{
    namespace mrd   = matcl::raw::details;

    double z0       = a.value / b.value;
    double r0       = - mrd::scal_func::fms_a(z0, b.value, a.value);
    double r1       = - mrd::scal_func::fms_a(z0, b.error, a.error);
    double c0       = r0 + r1;
    double d0       = b.value + b.error;
    double z1       = c0 / d0;

    return twofold::normalize_fast(z0, z1);
}

force_inline
twofold matcl::operator/(const twofold& a, double b)
{
    namespace mrd   = matcl::raw::details;

    double z0       = a.value / b;
    double r        = - mrd::scal_func::fms_a(z0, b, a.value);
    double c        = r + a.error;
    double z1       = c / b;

    return twofold::normalize_fast(z0, z1);
}

force_inline
twofold matcl::operator/(double a, const twofold& b)
{
    namespace mrd   = matcl::raw::details;

    double z0       = a / b.value;
    double r0       = - mrd::scal_func::fms_a(z0, b.value, a);
    double r1       = - z0 * b.error;
    double c0       = r0 + r1;
    double d0       = b.value + b.error;
    double z1       = c0 / d0;

    return twofold::normalize_fast(z0, z1);
}

force_inline
twofold matcl::sqrt(const twofold& a)
{
    namespace mrd   = matcl::raw::details;

    double z0       = mrd::scal_func::sqrt(a.value);
    double t1       = mrd::scal_func::fms_a(z0, z0, a.value);
    double t2       = a.error - t1;
    double d        = 2.0 * z0;
    double z1       = t2 / d;

    return twofold::normalize_fast(z0, z1);
}

force_inline
twofold matcl::abs(const twofold& a)
{
    namespace mrd   = matcl::raw::details;

    // we cannot use comparison operator due to existence of -0.0
    bool sign       = mrd::scal_func::signbit(a.value);

    if (sign == true)
        return twofold(-a.value, -a.error);
    else
        return a;
}

//--------------------------------------------------------------------
// represent z = x * y  as z = xy + err (exactly)
// requires: abs(x) <= 1e995, abs(y) <= 1e995, abs(x*y) <= 1e1021
// ensures: if x*y == 0 || abs(x*y) >= 1e-969, then x*y = xy + err
// modified version of Dekker function (author: Sylvie Boldo, available at
//  http://proval.lri.fr/gallery/Dekker.en.html)
inline
twofold matcl::twofold_mult_dekker(double x, double y)
{
    double xy   = x * y;

    double C    = (double)(0x8000001); //2^27 + 1
    double px   = x * C;
    double qx   = x - px;
    double hx   = px + qx;
    double tx   = x - hx;

    double py   = y * C;
    double qy   = y - py;
    double hy   = py + qy;
    double ty   = y - hy;

    double err  = -xy + hx * hy;
    
    err         += hx*ty;
    err         += hy*tx;
    err         += tx*ty;

    return twofold(xy, err);
};

inline
double matcl::fma_dekker(double x, double y, double z)
{
    twofold xy  = twofold_mult_dekker(x, y);
    twofold s   = z + xy;
    
    return s.value;
};

inline
float matcl::fma_dekker(float x, float y, float z)
{
    //for floats we can use double arithmetics
    double xd   = x;
    double yd   = y;
    double zd   = z;

    // res is exact
    double res  = xd * yd;

    // 0.5 ulp error in double precision
    res         = res + zd;

    // almost 0.5 ulp error in single precision
    return (float)res;
};

}