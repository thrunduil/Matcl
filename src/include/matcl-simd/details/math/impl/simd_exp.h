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

#include "matcl-simd/details/arch/simd_impl.h"
#include "matcl-simd/details/math/simd_math_func_def.h"

namespace matcl { namespace simd { namespace details
{

struct MATCL_SIMD_EXPORT exp_table_double
{
    // lookup table size
    static const int L  = 256;

    // polynomial coefficients
    static const double exp_poly[4];

    // precomputed values of exp(i * log(2) / L)
    static double lookup_table_arr[L + 1];
    
    static double* lookup_table;

    #ifdef MATCL_SIMD_GENERATE_TABLES
        // generate lookup table parameters
        static void generate_lookup_table(std::ostream& os, int num_values_in_row);
    #endif

    template<class Arg>
    force_inline
    static Arg eval(const Arg& x)
    {
        return estrin<4>(x, exp_table_double::exp_poly);
    };
};

struct MATCL_SIMD_EXPORT exp_table_float
{
    // polynomial coefficients
    static const float exp_poly[6];

    template<class Arg>
    force_inline
    static Arg eval(const Arg& x)
    {
        Arg p   =  estrin<6>(x, exp_table_float::exp_poly);
        p       = fma_f(p, x, Arg(1));
        return p;
    };
};

template<int Bits, class Simd_tag>
struct simd_exp<double, Bits, Simd_tag>
{
    using simd_type     = simd<double, Bits, Simd_tag>;

    static const int L  = exp_table_double::L;

    force_inline
    static simd_type eval(const simd_type& a)
    {
        const simd_type inv_log2 = simd_type(L * 1.442695040888963407359924681001);

        // high part of log(2)/L stored with precision 35 bits
        const simd_type log2_hi  = simd_type(0.69314718057285062969 / L);

        // low part of log(2)/L
        const simd_type log2_lo  = simd_type(-1.2905320270077143679e-11 / L);                                        

        const simd_type M       = simd_type(double(L));
        const simd_type M_inv   = simd_type(1.0 / L);
        const simd_type max_k   = simd_type(1022);
        const simd_type max_l   = simd_type(double(L/2));

        //-----------------------------------------------------------------
        //                      reduction
        //-----------------------------------------------------------------

        // this reduction scheme will produce inaccurate d = a - k*inv_log2, when a is close
        // to k*inv_log2 (worst case for 7804143460206699 x 2^?49), but then d ~ 0 and 
        // result exp(d) ~ 1 is accurate
        simd_type kl    = round(inv_log2 * a);
        simd_type x     = fnma_f(kl, log2_hi, a);   // a - k*LH; exact for |kl| < 2^18
        x               = fnma_f(kl, log2_lo, x);   // x - k*LL

        simd_type k     = round(M_inv * kl);
        simd_type l     = fnma_f(k, M, kl);         // kl - k*M        
        l               = min(l, max_l);

        bool normal     = all(leq(abs(k), max_k));
        
        simd_type mult  = simd_type::gather(exp_table_double::lookup_table, l.convert_to_int32());

        //-----------------------------------------------------------------
        //                      approximation
        //-----------------------------------------------------------------
        // remez of (exp(x) - 1)/x on [-log(2)/L/2, log(2)/L/2]
        
        simd_type p     = exp_table_double::eval(x);
        p               = p * x;

        //-----------------------------------------------------------------
        //                      finalization
        //-----------------------------------------------------------------
        if (normal == false)
            return process_overflows(a, p, k, l);

        simd_type powk  = pow2k(k);

        // res = (p + 1) * exp(M_inv * l) * 2^k        
        simd_type c     = fma_f(p, mult, mult);
        c               = c * powk;

        return c;
    };

    static simd_type process_overflows(const simd_type& a, const simd_type& p, 
                                       const simd_type& k, const simd_type& l)
    {        
        const simd_type val_inf = simd_type(std::numeric_limits<double>::infinity());
        const simd_type max_a   = simd_type(1000.0);
        const simd_type min_a   = simd_type(-1000.0);
        
        simd_type k1    = round(k * simd_type(0.5));
        simd_type k2    = k - k1;

        simd_type mult  = simd_type::gather(exp_table_double::lookup_table, l.convert_to_int32());
        simd_type c     = fma_f(p, mult, mult);

        c               = c * pow2k(k1);  
        c               = c * pow2k(k2);

        c               = if_zero_else(leq(a, min_a), c);
        c               = if_then_else(geq(a, max_a), val_inf, c);
        c               = if_nan_else(is_nan(a), c);

        return c;
    };
};

template<int Bits, class Simd_tag>
struct simd_exp<float, Bits, Simd_tag>
{
    using simd_type     = simd<float, Bits, Simd_tag>;

    force_inline
    static simd_type eval(const simd_type& a)
    {
        const simd_type one      = simd_type::one();
        const simd_type inv_log2 = simd_type(1.442695040888963407359924681001f);

        // high part of log(2) stored with precision 16 bits
        const simd_type log2_hi  = simd_type(0.693145751953125f);

        // low part of log(2)
        const simd_type log2_lo  = simd_type(1.4286068203094172321e-06f);

        const simd_type max_k   = simd_type(126.f);

        //-----------------------------------------------------------------
        //                      reduction
        //-----------------------------------------------------------------
        simd_type k     = round(inv_log2 * a);
        simd_type x     = fnma_f(k, log2_hi, a);    // a - k*LH; exact for |k| < 256
        x               = fnma_f(k, log2_lo, x);    // x - k*LL

        bool normal     = all(leq(abs(k), max_k));
        
        //-----------------------------------------------------------------
        //                      approximation
        //-----------------------------------------------------------------
        // remez of (exp(x) - 1)/x on [-log(2)/2, log(2)/2]
        
        simd_type p     = exp_table_float::eval(x);

        //-----------------------------------------------------------------
        //                      finalization
        //-----------------------------------------------------------------
        if (normal == false)
            return process_overflows(a, p, k);

        simd_type c     = p * pow2k(k);

        return c;
    };

    static simd_type process_overflows(const simd_type& a, const simd_type& p, 
                                       const simd_type& k)
    {
        const simd_type val_inf = simd_type(std::numeric_limits<float>::infinity());
        const simd_type max_a   = simd_type(120.0f);
        const simd_type min_a   = simd_type(-120.0f);

        simd_type k1    = round(k * simd_type(0.5f));
        simd_type k2    = k - k1;

        simd_type c     = p * pow2k(k1);  
        c               = c * pow2k(k2);

        c               = if_zero_else(leq(a, min_a), c);
        c               = if_then_else(geq(a, max_a), val_inf, c);
        c               = if_nan_else(is_nan(a), c);

        return c;
    };
};

}}}
