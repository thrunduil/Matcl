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

#include "matcl-mp/mp_float.h"

namespace matcl { namespace mp { namespace details
{

// error propagation for a function f(x,y) is defined as
// max{|f(x+dx, y+dx) / f(x,y)| - 1} s.t. |dx| <= a*|x|, |dy| <= a*|y|
//
// first order approximations of error propagations for:
//  x + y           : ep_plus = |x/(x + y)| * a + |y/(x + y)| * b
//                  : for x,y >= 0: ep_plus <= max(a, b)
//  x - y           : ep_minus  = |x/(x - y)| * a + |y/(x - y)| * b
//                              <= a + |y/(x - y)| * (a + b)
//  x * y           : ep_mult   = a + b
//  x / y           : ep_div    = a + b
//  log(x)          : ep_log    = 1/log(x) * a
//  log(1+x)        : ep_log1p  = |x/(log(x + 1)*(x + 1))| * a
//                    for x ~ 0 = (1 - x/2) * a
//                    for x > 0 = ep_log1p <= a
//  exp(x)          : ep_exp    = |x| * a
//  sqrt(x)         : ep_sqrt   = 1/2 * a
//  root(x,n)       : ep_root   = 1/n * a
//  re(1/(x+i*y))   : ep_inv_r  = a + 2 * b
//  im(1/(x+i*y))   : ep_inv_i  = 2 * a + b
//  re(exp(x+i*y))  : ep_exp_r  = |x|*a + |y*tan(y)| * b
//  im(exp(x+i*y))  : ep_exp_i  = |x|*a + |y*cot(y)| * b
//  hypot(x, y)     : ep_hyp    = x^2/(x^2 + y^2) * a + y^2/(x^2 + y^2) * b <= max(a, b)
//  atan2(x,y)      : ep_x      = |x / ( y * atan(x/y) * (x^2/y^2 + 1))| * a <= a
//                  : ep_y      = |x / ( y * atan(x/y) * (x^2/y^2 + 1))| * b <= b
//  re(log(x+iy))   : ep_x      = |x / ( log(x^2 + y^2) * (x^2 + y^2)^(1/2) )| * a 
//                              <= |1 / log(x^2 + y^2)| * a
//                  : ep_y      = |(2*y)/( log(x^2 + y^2) * (x^2 + y^2)^(1/2))| * b
//                              <= |2/log(x^2 + y^2)| * b
//  im(log(x+iy))   : ep_x      = |x/(log(x^2 + y^2)*(x^2 + y^2)^(1/2))| * a 
//                              <= |1/log(x^2 + y^2)| * a
//                  : ep_y      = |(x*y)/(atan2(y, x) *(x^2 + y^2))| * b
//                              <= |1/(2 * atan2(y, x))| * b
//  asin(x)         : ep_x      = |x/(asin(x)*(1 - x^2)^(1/2))| * a
//  atan(x)         : ep_x      = |x/(atan(x)*(x^2 + 1))| * a <= a
// re[(x+iy)*(u+iw)]: ep:       = |(u*x)/(u*x - w*y)| * ex + |(w*y)/(u*x - w*y)| * ey
//                              + |(u*x)/(u*x - w*y)| * eu + |(w*y)/(u*x - w*y)| * ew
// im[(x+iy)*(u+iw)]: ep:       = |(w*x)/(u*y + w*x)| * ex + |(u*y)/(u*y + w*x)| * ey
//                              + |(u*y)/(u*y + w*x)| * eu + |(w*x)/(u*y + w*x)| * ew

#if MATCL_TEST_MP == 1
    static const bool test_mode = true;
#else
    static const bool test_mode = false;
#endif

//only a safeguard; precision 1ulp is obtained even for zero extra bits
inline size_t extra_bits()              { return test_mode ? 0 : 1; };
//additional precision for constants
inline size_t extra_bits_constants()    { return 4; };

//maximum number of precision increase
inline int get_max_precision_iters()    { return 10; };

//division real/complex
inline size_t error_exact_div_rc()      { return 7; };      // 3.5 ulp error
inline size_t error_exact_div_cc()      { return 6; };      // 3 ulp error

inline size_t error_exact_abs2()        { return 2; };      // 1 ulp error
inline size_t error_exact_hypot()       { return 2; };      // 1 ulp error
inline size_t error_exact_sign()        { return 3; };      // 1.5 ulp error
inline size_t error_exact_exp()         { return 3; };      // 1.5 ulp error
inline size_t error_exact_exp_mult()    { return 3; };      // 1.5 ulp error
inline size_t error_exact_sin()         { return 3; };      // 1.5 ulp error
inline size_t error_exact_cos()         { return 3; };      // 1.5 ulp error
inline size_t error_exact_sinh()        { return 3; };      // 1.5 ulp error
inline size_t error_exact_cosh()        { return 3; };      // 1.5 ulp error
inline size_t error_exact_log2_r()      { return 3; };      // 1.5 ulp error
inline size_t error_exact_log10_r()     { return 3; };      // 1.5 ulp error
inline size_t error_exact_sqrt()        { return 4; };      // 2 ulp error
inline size_t error_exact_log2_c()      { return 4; };      // 2 ulp error
inline size_t error_exact_log10_c()     { return 4; };      // 2 ulp error
inline size_t error_exact_dot2_a()      { return 4; };      // 2 ulp error
inline size_t error_exact_cot()         { return 8; };      // 4 ulp error
inline size_t error_exact_coth()        { return 8; };      // 4 ulp error
inline size_t error_exact_sec()         { return 8; };      // 4 ulp error
inline size_t error_exact_csc()         { return 8; };      // 4 ulp error
inline size_t error_exact_sech()        { return 8; };      // 4 ulp error
inline size_t error_exact_csch()        { return 8; };      // 4 ulp error
inline size_t error_exact_atanh()       { return 10; };     // 5 ulp error
inline size_t error_exact_tanh()        { return 15; };     // 7.5 ulp error

inline size_t ilogb_u(size_t v)
{
    if (v <= 1)
        return 0;
    else
        return std::ilogb(v-1)+1;
}

inline size_t ilogb_u(Real v)
{
    if (v <= 1)
        return 0;
    else
        return (size_t)std::ceil(std::log2(v));
}

inline precision extend_prec_div_rc(precision prec)
{
    size_t base = prec + ilogb_u(error_exact_div_rc()) + extra_bits();
    return precision(base);
};

inline precision extend_prec_div_cc(precision prec)
{
    size_t base = prec + ilogb_u(error_exact_div_cc()) + extra_bits();
    return precision(base);
};

inline precision extend_prec_abs2(precision prec)
{    
    size_t base = prec + ilogb_u(error_exact_abs2()) + extra_bits();
    return precision(base);
}

inline precision extend_prec_exp(precision prec)
{
    size_t base = prec + ilogb_u(error_exact_exp()) + extra_bits();
    return precision(base);
}

inline precision extend_prec_exp_mult(precision prec, Real mult)
{
    Real lp     = mult * Real(error_exact_exp_mult());
    size_t base = prec + ilogb_u(lp) + extra_bits();
    return precision(base);
}

inline precision extend_prec_sin(precision prec)
{ 
    size_t base = prec + ilogb_u(error_exact_sin()) + extra_bits();
    return precision(base);
};

inline precision extend_prec_cos(precision prec)
{ 
    size_t base = prec + ilogb_u(error_exact_cos()) + extra_bits();
    return precision(base);
};

inline precision extend_prec_sinh(precision prec)
{ 
    size_t base = prec + ilogb_u(error_exact_sinh()) + extra_bits();
    return precision(base);
};

inline precision extend_prec_cosh(precision prec)
{ 
    size_t base = prec + ilogb_u(error_exact_cosh()) + extra_bits();
    return precision(base);
};

inline precision extend_prec_sec(precision prec)
{ 
    size_t base = prec + ilogb_u(error_exact_sec()) + extra_bits();
    return precision(base);
};

inline precision extend_prec_csc(precision prec)
{ 
    size_t base = prec + ilogb_u(error_exact_csc()) + extra_bits();
    return precision(base);
};

inline precision extend_prec_sech(precision prec)
{ 
    size_t base = prec + ilogb_u(error_exact_sech()) + extra_bits();
    return precision(base);
};

inline precision extend_prec_csch(precision prec)
{ 
    size_t base = prec + ilogb_u(error_exact_csch()) + extra_bits();
    return precision(base);
};

inline precision extend_prec_tanh(precision prec)
{ 
    size_t base = prec + ilogb_u(error_exact_tanh()) + extra_bits();
    return precision(base);
};

inline precision extend_prec_cot(precision prec)
{ 
    size_t base = prec + ilogb_u(error_exact_cot()) + extra_bits();
    return precision(base);
};

inline precision extend_prec_coth(precision prec)
{ 
    size_t base = prec + ilogb_u(error_exact_coth()) + extra_bits();
    return precision(base);
};

inline precision extend_prec_dynamic(precision prec, Real half_ulp_error)
{ 
    size_t base = prec + ilogb_u(half_ulp_error) + extra_bits();
    return precision(base);
};

inline precision extend_prec_hypot(precision prec)
{
    size_t base = prec + ilogb_u(error_exact_hypot()) + extra_bits();
    return precision(base);
}

inline precision extend_prec_dot2_a(precision prec)
{
    size_t base = prec + ilogb_u(error_exact_dot2_a()) + extra_bits();
    return precision(base);
}

inline precision extend_prec_log2_r(precision prec)
{
    size_t base = prec + ilogb_u(error_exact_log2_r()) + extra_bits();
    return precision(base);
}

inline precision extend_prec_log10_r(precision prec)
{
    size_t base = prec + ilogb_u(error_exact_log10_r()) + extra_bits();
    return precision(base);
}

inline precision extend_prec_log2_c(precision prec)
{
    size_t base = prec + ilogb_u(error_exact_log2_c()) + extra_bits();
    return precision(base);
}

inline precision extend_prec_log10_c(precision prec)
{
    size_t base = prec + ilogb_u(error_exact_log10_c()) + extra_bits();
    return precision(base);
}

inline precision extend_prec_sqrt(precision prec)
{
    size_t base = prec + ilogb_u(error_exact_sqrt()) + extra_bits();
    return precision(base);
}

inline precision extend_prec_sign(precision prec)
{
    size_t base = prec + ilogb_u(error_exact_sign()) + extra_bits();
    return precision(base);
}

inline precision extend_prec_atanh(precision prec)
{
    size_t base = prec + ilogb_u(error_exact_atanh()) + extra_bits();
    return precision(base);
}

};};};