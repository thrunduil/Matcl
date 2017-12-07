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
#include "matcl-simd/details/func/simd_func_def.h"
#include "matcl-core/details/float/fma_dekker_simd.inl"

namespace matcl { namespace simd { namespace details
{

struct exp_table_double
{
    // lookup table size
    static const int L  = 256;

    // polynomial coefficients
    static const double exp_poly[4];

    // precomputed values of exp(i * log(2) / L)
    static double* lookup_table;
};

template<int Bits, class Simd_tag>
struct simd_exp<double, Bits, Simd_tag>
{
    using simd_type = simd<double, Bits, Simd_tag>;

    static const int L          = exp_table_double::L;

    force_inline
    static simd_type eval(const simd_type& a)
    {
        const simd_type inv_log2 = simd_type(L * 1.442695040888963407359924681001);
        const simd_type log2_hi  = simd_type(0.6931471805598903 / L);
        const simd_type log2_lo  = simd_type(5.497923018708371e-14 / L);                                        

        const simd_type M       = simd_type(double(L));
        const simd_type M_inv   = simd_type(1.0 / L);
        const simd_type max_k   = simd_type(1022);

        //-----------------------------------------------------------------
        //                      reduction
        //-----------------------------------------------------------------
        simd_type kl    = round(inv_log2 * a);
        simd_type x     = fnma_f(kl, log2_hi, a);   // a - k*LH
        x               = fnma_f(kl, log2_lo, x);   // x - k*LL

        simd_type k     = round(M_inv * kl);
        simd_type l     = fnma_f(k, M, kl);         // kl - k*M

        bool normal     = all(leq(abs(k), max_k));
        
        //-----------------------------------------------------------------
        //                      approximation
        //-----------------------------------------------------------------
        // remez of (exp(x) - 1)/x on [-log(2)/L/2, log(2)/L/2]
        
        simd_type p     = simd_type::broadcast(exp_table_double::exp_poly+3);
        p               = fma_f(p, x, simd_type::broadcast(exp_table_double::exp_poly+2));
        p               = fma_f(p, x, simd_type::broadcast(exp_table_double::exp_poly+1));
        p               = fma_f(p, x, simd_type::broadcast(exp_table_double::exp_poly+0));
        
        p               = p * x;

        //-----------------------------------------------------------------
        //                      finalization
        //-----------------------------------------------------------------
        if (normal == false)
            return process_overflows(a, p, k, l);

        // res = (c + 1) * exp(M_inv * l) * 2^k
        simd_type mult  = simd_type::gather(exp_table_double::lookup_table, l.convert_to_int32());
        simd_type c     = fma_f(p, mult, mult);
        
        c               = c * pow2k(k);

        return c;
    };

    static simd_type process_overflows(simd_type a, simd_type p, simd_type k, simd_type l)
    {
        const simd_type max_l   = simd_type(double(L/2));
        const simd_type val_inf = simd_type(std::numeric_limits<double>::infinity());
        const simd_type max_a   = simd_type(1000.0);
        const simd_type min_a   = simd_type(-1000.0);

        l               = min(l, max_l);
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

}}}
