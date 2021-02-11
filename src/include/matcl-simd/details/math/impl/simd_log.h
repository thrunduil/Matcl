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

#include "matcl-simd/details/arch/simd_impl.h"
#include "matcl-simd/details/math/simd_math_func_def.h"
#include "matcl-simd/poly/poly_eval.h"

namespace matcl { namespace simd { namespace details
{

//-----------------------------------------------------------------------
//                              DOUBLE
//-----------------------------------------------------------------------

struct MATCL_SIMD_EXPORT simd_log_approx_double_data
{
    static const double poly[7];
    static const double poly_scalar[4*2];
};

template<int Bits, class Tag>
struct simd_log_approx_double
{
    using simd_type     = simd<double, Bits, Tag>;

    force_inline
    static simd_type eval(const simd_type& s)
    {        
        return estrin<7>(s, simd_log_approx_double_data::poly);
    };
};

#if MATCL_ARCHITECTURE_HAS_SSE2
    template<>
    struct simd_log_approx_double<128, scalar_sse_tag>
    {
        using simd_type         = simd<double, 128, scalar_sse_tag>;
        using simd_type_sse     = simd<double, 128, sse_tag>;    

        force_inline
        static simd_type eval(const simd_type& s2_0)
        {        
            simd_type_sse s2    = s2_0.as_vector();
            simd_type_sse s4    = s2 * s2;

            const double* poly  = simd_log_approx_double_data::poly_scalar;
            simd_type_sse p     = simd_type_sse::load(poly + 3 * 2);

            p                   = fma_f(p, s4, simd_type_sse::load(poly + 2 * 2));
            p                   = fma_f(p, s4, simd_type_sse::load(poly + 1 * 2));
            p                   = fma_f(p, s4, simd_type_sse::load(poly + 0 * 2));

            p                   = fma_f(s2, p.extract_low(), p.extract_high()); 
            return simd_type(p);
        };
    };
#endif

template<int Bits, class Tag> 
struct simd_log<double, Bits, Tag>
{
    using simd_type     = simd<double, Bits, Tag>;

    force_inline
    static void reduction(const simd_type& x, simd_type& k0, simd_type& k, simd_type& frac)
    {
        const simd_type one     = simd_type::one();

        // sqrt(2)/2
        const simd_type min_x   = simd_type(0.70710678118654752440084436210485);

        // frac in [0.5, 1)
        simd_type frac0     = fraction(x);
        k0                  = exponent(x);

        // if (cond) { frac = frac * 2; k = k - 1 };
        simd_type cond      = lt(frac0,  min_x);

        frac                = if_add(cond, frac0, frac0);
        k                   = if_sub(cond, k0, one);  
    };

    force_inline
    static simd_type finalization(const simd_type& p, const simd_type& k, const simd_type& s,
                                  const simd_type& f)
    {
        // high part of log(2) stored with precision 35 bits
        const simd_type log2_hi = simd_type(0.69314718057285062969);

        // low part of log(2)
        const simd_type log2_lo = simd_type(-1.2905320270077143679e-11);

        simd_type res       = fma_f(s, p, k * log2_lo);
        res                 = f + res;
        res                 = fma_f(k, log2_hi, res);

        return res;
    }

    force_inline
    static simd_type eval(const simd_type& x)
    {
        const simd_type one     = simd_type::one();
        const simd_type zero    = simd_type::zero();
        
        const simd_type max_k   = simd_type(1021);

        //-----------------------------------------------------------------
        //                      reduction
        //-----------------------------------------------------------------

        simd_type k0, k, frac;

        reduction(x, k0, k, frac);

        simd_type f         = frac - one;
        simd_type s         = f / (frac + one);

        // normal return for x > 0 && |k| <= max_k (i.e. not zero, inf, NaN
        // and denormals, but also some big numbers are excluded)

        simd_type cond_norm = leq(abs(k), max_k) && gt(x, zero);
        bool normal         = all(cond_norm);

        //-----------------------------------------------------------------
        //                      approximation
        //-----------------------------------------------------------------
        // h(s)     = [log( (1 + sqrt(s))/(1 - sqrt(s)) ) / sqrt(s) - 2]/ s
        //          s in  [0, (3 - 2*sqrt(2))^2]
        // log(x)   = s * (s^2 * h(s^2) - f) + f

        simd_type s2        = s * s;
        simd_type p         = simd_log_approx_double<Bits, Tag>::eval(s2);
        p                   = fms_f(s2, p, f);

        //-----------------------------------------------------------------
        //                      finalization
        //-----------------------------------------------------------------
        simd_type res       = finalization(p, k, s, f);

        if (normal == false)
            return process_overflows(x, k0, res);
        else 
            return res;
    }

    static simd_type process_overflows(const simd_type& x, const simd_type& k, 
                                       const simd_type& res)
    {
        const simd_type k_infnan    = simd_type(1025.0);
        const simd_type k_denorm    = simd_type(-1022.0);
        const simd_type minus_inf   = simd_type(-std::numeric_limits<double>::infinity());
        const simd_type nan         = simd_type(std::numeric_limits<double>::quiet_NaN());
        const simd_type zero        = simd_type::zero();

        simd_type denorm    = eeq(k, k_denorm) && gt(x, zero);

        simd_type res2      = res;

        // for x = +-Inf, NaN return x
        res2                = if_then_else(geq(k, k_infnan), x, res2);
        
        // for x = +-0  return -Inf
        res2                = if_then_else(eeq(x, zero), minus_inf, res2);

        // for x < 0 return nan
        res2                = if_then_else(lt(x, zero), nan, res2);
        
        if (any(denorm) == true)
        {
            // eval for denormal arguments
            simd_type res_den;
            res_den         = eval_denormal(x);

            res2            = if_then_else(denorm, res_den, res2);
        }

        return res2;
    };

    static simd_type eval_denormal(const simd_type& x0)
    {
        const simd_type one     = simd_type::one();

        // 2^52
        const simd_type pow2_52 = simd_type(4503599627370496.0);

        //-----------------------------------------------------------------
        //                      reduction
        //-----------------------------------------------------------------

        // move x to normal range of exponents
        simd_type x         = x0 * pow2_52;

        simd_type k0, k, frac;
        reduction(x, k0, k, frac);

        k                   = k - simd_type(52.0);
        simd_type f         = frac - one;
        simd_type s         = f / (frac + one);

        //-----------------------------------------------------------------
        //                      approximation
        //-----------------------------------------------------------------
        simd_type s2        = s * s;
        simd_type p         = simd_log_approx_double<Bits, Tag>::eval(s2);
        p                   = fms_f(s2, p, f);

        //-----------------------------------------------------------------
        //                      finalization
        //-----------------------------------------------------------------
        simd_type res       = finalization(p, k, s, f);

        return res;
    }
};

//-----------------------------------------------------------------------
//                              FLOAT
//-----------------------------------------------------------------------

struct MATCL_SIMD_EXPORT simd_log_approx_float_data
{
    static const float poly[9];
    static const float poly_scalar[12];
};

template<int Bits, class Tag>
struct simd_log_approx_float
{
    using simd_type     = matcl::simd::simd<float, Bits, Tag>;

    force_inline
    static simd_type eval(const simd_type& s)
    {        
        return estrin<9>(s, simd_log_approx_float_data::poly);
    };
};

template<int Bits>
struct simd_log_approx_float<Bits, scalar_sse_tag>
{
    using simd_type     = simd<float, Bits, scalar_sse_tag>;
    using simd_base     = simd<float, Bits, sse_tag>;

    force_inline
    static simd_type eval(const simd_type& x0)
    {   
        simd_base x     = x0.as_vector();
        simd_base x2    = x * x;
        simd_base x4    = x2 * x2;

        simd_base p1    = simd_base::load(simd_log_approx_float_data::poly_scalar + 8);
        simd_base p2    = simd_base::load(simd_log_approx_float_data::poly_scalar + 2);

        p1              = fma_f(p1, x2, simd_base::load(simd_log_approx_float_data::poly_scalar + 6));
        p2              = fma_f(p2, x2, simd_base::load(simd_log_approx_float_data::poly_scalar + 0));
        p1              = fma_f(p1, x2, simd_base::load(simd_log_approx_float_data::poly_scalar + 4));
        simd_base p     = fma_f(p1, x4, p2);

        simd_base pl    = p;
        simd_base ph    = p.select<1, 0, 0, 0>();

        simd_base res   = fma_f(ph, x, pl);
        return simd_type(res);
    };
};

template<int Bits, class Tag> 
struct simd_log<float, Bits, Tag>
{
    using simd_type     = matcl::simd::simd<float, Bits, Tag>;

    force_inline
    static void reduction(const simd_type& x, simd_type& k0, simd_type& k, simd_type& frac)
    {
        const simd_type one     = simd_type::one();

        // 2/3
        const simd_type min_x   = simd_type(0.6666666666666666f);

        // frac in [0.5, 1)
        simd_type frac0     = fraction(x);
        k0                  = exponent(x);

        // if (cond) { frac = frac * 2; k = k - 1 };
        simd_type cond      = lt(frac0,  min_x);

        frac                = if_add(cond, frac0, frac0);
        k                   = if_sub(cond, k0, one);  
    };

    force_inline
    static simd_type finalization(const simd_type& p, const simd_type& k, const simd_type& s)
    {
        // high part of log(2)/L stored with precision 16 bits
        const simd_type log2_hi  = simd_type(0.693145751953125f);

        // low part of log(2)/L
        const simd_type log2_lo  = simd_type(1.4286068203094172321e-06f);

        simd_type s2        = s * s;
        simd_type res       = p * s2;
        res                 = fma_f(k, log2_lo, res);
        res                 = res + s;
        res                 = fma_f(k, log2_hi, res);

        return res;
    }

    force_inline
    static simd_type eval(const simd_type& x)
    {
        const simd_type one     = simd_type::one();
        const simd_type zero    = simd_type::zero();
        
        const simd_type max_k   = simd_type(125);

        //-----------------------------------------------------------------
        //                      reduction
        //-----------------------------------------------------------------

        simd_type k0, k, frac;

        reduction(x, k0, k, frac);

        // f in [-1/3, 1/3]
        simd_type s         = frac - one;

        // normal return for x > 0 && |k| <= max_k (i.e. not zero, inf, NaN
        // and denormals, but also some big numbers are excluded)

        simd_type cond_norm = leq(abs(k), max_k) && gt(x, zero);
        bool normal         = all(cond_norm);

        //-----------------------------------------------------------------
        //                      approximation
        //-----------------------------------------------------------------
        // h(s)     = (log(1+s)/s - 1)/s
        //          s in  [-1/3, 1/3]
        // log(x)   = h(x) * s^2 + s

        simd_type p         = simd_log_approx_float<Bits, Tag>::eval(s);

        //-----------------------------------------------------------------
        //                      finalization
        //-----------------------------------------------------------------
        simd_type res       = finalization(p, k, s);

        if (normal == false)
            return process_overflows(x, k0, res);
        else 
            return res;
    }

    static simd_type process_overflows(const simd_type& x, const simd_type& k, 
                                       const simd_type& res)
    {
        const simd_type k_infnan    = simd_type(129.0f);
        const simd_type k_denorm    = simd_type(-126.0);
        const simd_type minus_inf   = simd_type(-std::numeric_limits<float>::infinity());
        const simd_type nan         = simd_type(std::numeric_limits<float>::quiet_NaN());
        const simd_type zero        = simd_type::zero();

        simd_type denorm    = eeq(k, k_denorm) && gt(x, zero);

        simd_type res2      = res;

        // for x = +-Inf, NaN return x
        res2                = if_then_else(geq(k, k_infnan), x, res2);
        
        // for x = +-0  return -Inf
        res2                = if_then_else(eeq(x, zero), minus_inf, res2);

        // for x < 0 return nan
        res2                = if_then_else(lt(x, zero), nan, res2);
        
        if (any(denorm) == true)
        {
            // eval for denormal arguments
            simd_type res_den;
            res_den         = eval_denormal(x);

            res2            = if_then_else(denorm, res_den, res2);
        }

        return res2;
    };

    static simd_type eval_denormal(const simd_type& x0)
    {
        const simd_type one     = simd_type::one();

        // 2^23
        const simd_type pow2_23 = simd_type(8388608.0f);

        //-----------------------------------------------------------------
        //                      reduction
        //-----------------------------------------------------------------

        // move x to normal range of exponents
        simd_type x         = x0 * pow2_23;

        simd_type k0, k, frac;
        reduction(x, k0, k, frac);

        // f in [-1/3, 1/3]
        simd_type s         = frac - one;
        k                   = k - simd_type(23.0f);

        //-----------------------------------------------------------------
        //                      approximation
        //-----------------------------------------------------------------
        simd_type p         = simd_log_approx_float<Bits, Tag>::eval(s);

        //-----------------------------------------------------------------
        //                      finalization
        //-----------------------------------------------------------------
        simd_type res       = finalization(p, k, s);
        return res;
    }
};

}}}
