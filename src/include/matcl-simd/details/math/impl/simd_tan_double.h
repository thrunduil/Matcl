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
#include "matcl-simd/poly/poly_eval.h"
#include "matcl-core/float/twofold.h"
#include "matcl-simd/details/math/impl/payne_hanek.inl"
#include "matcl-simd/details/math/impl/simd_tancot_helpers.h"
#include "matcl-simd/details/math/impl/pi2_reduction.h"

#pragma warning(push)
#pragma warning(disable: 4127)  //conditional expression is constant

namespace matcl { namespace simd { namespace details
{

//-----------------------------------------------------------------------
//                              DOUBLE
//-----------------------------------------------------------------------
struct MATCL_SIMD_EXPORT simd_tancot_table_double_data
{
    static const double poly_tan_nom[4];
    static const double poly_tan_den[4];
    static const double poly_tan_nomden[4*2];
};

template<int Bits, class Tag>
struct simd_tancot_table_double
{
    using simd_type         = ms::simd<double, Bits, Tag>;

    force_inline
    static void eval(const simd_type& x, simd_type& p_tan)
    {
        simd_type p_tan_nom = estrin<4>(x, simd_tancot_table_double_data::poly_tan_nom);
        simd_type p_tan_den = estrin<4>(x, simd_tancot_table_double_data::poly_tan_den);
        p_tan               = p_tan_nom / p_tan_den;
    };
};

#if MATCL_ARCHITECTURE_HAS_SSE2
    template<>
    struct simd_tancot_table_double<128, ms::scalar_sse_tag>
    {
        using simd_type         = ms::simd<double, 128, ms::scalar_sse_tag>;
        using simd_sse          = ms::simd<double, 128, ms::sse_tag>;    

        force_inline
        static void eval(const simd_type& x, simd_type& p_tan)
        {
            simd_sse x_sse      = x.as_vector();
            const simd_sse* poly= reinterpret_cast<const simd_sse*>(simd_tancot_table_double_data::poly_tan_nomden);
            simd_sse p_nomden   = estrin<4>(x_sse, poly);
            simd_type p_nom     = simd_type(p_nomden.first());
            simd_type p_den     = simd_type(p_nomden.extract_high());

            p_tan               = p_nom / p_den;
        };
    };
#endif

template<int Bits, class Tag>
struct simd_tancot<double, Bits, Tag>
{
    using simd_type         = ms::simd<double, Bits, Tag>;
    using simd_int          = typename ms::simd<int32_t, Bits, Tag>::simd_half;

    force_inline
    static void approximation(const simd_type& x, const simd_type& xv, const simd_type& xe, 
                            simd_type& ret_tan, simd_type& ret_cot)
    {
        const simd_type zero    = simd_type::zero();

        // h(x)     = (tan(sqrt(x))/sqrt(x) - 1)/x - 1/2
        // tan(x)   = x^3 * (h(x^2) + 1/2) + x 
        // range    : [0, (pi/4)^2]

        using twofold       = matcl::twofold<simd_type>;

        simd_type sign_inf  = ms::copysign(simd_type(std::numeric_limits<double>::infinity()), x);

        twofold x2t         = twofold_mult_f(xv, xv);

        // eval polynomials
        simd_type p_tan;
        simd_tancot_table_double<Bits, Tag>::eval(x2t.value, p_tan);

        p_tan               = p_tan + simd_type(0.5);
        
        // p_tan * xv^3 must be computed with higher precision
        twofold pt_x2       = twofold_mult_f_without_norm(p_tan, x2t);
        twofold pt_x3       = twofold_mult_f_without_norm(pt_x2, xv);

        // we also need the error term in order to compute 1.0 / pt_t
        twofold pt_t        = twofold_sum_sorted_without_norm(xv, pt_x3);

        // tan(xv + xe) ~ tan(xv) + [tan(xv)^2 + 1] * xe
        simd_type tan_e     = fma_f(pt_t.value * pt_t.value, xe, xe);

        pt_t.error          = pt_t.error + tan_e;
        ret_tan             = pt_t.value + pt_t.error;

        twofold pc_t        = twofold_inv_f_without_norm(pt_t);

        ret_cot             = pc_t.value + pc_t.error;
        ret_cot             = if_then_else(neq(x, zero), ret_cot, sign_inf);
    }

    template<int Version>
    static simd_type process_overflow(const simd_type& x, const simd_type& in_range)
    {
        const simd_type nan_v   = simd_type(std::numeric_limits<double>::quiet_NaN());

        simd_type fin   = ms::is_finite(x);
        simd_type ret   = eval_with_full_reduction<Version>(x, in_range);
        ret             = ms::if_then_else(fin, ret, nan_v);

        return ret;
    }

    template<int Version>
    force_inline
    static simd_type eval_with_full_reduction(const simd_type& x, const simd_type& in_range)
    {
        simd_type xv, xe;
        simd_int q_lo;

        if (any(in_range) == true)
            xv  = pi2_reduction<double, Bits, Tag>::reduction_CW3(x, xe, q_lo); 

        pi2_reduction<double, Bits, Tag>::reduction_full(x, in_range, xv, xe, q_lo);

        simd_type pt, pc;
        approximation(x, xv, xe, pt, pc);

        simd_type ret       = simd_tancot_reconstruction<double, Bits, Tag>
                                ::eval<Version>(q_lo, pt, pc);
        return ret;
    }

    template<int Version>
    force_inline
    static simd_type eval_tancot(const simd_type& x)
    {        
        // 2^19 * pi/2
        const simd_type max_x   = simd_type(823549.6);

        simd_type in_range  = gt(max_x, abs(x));
        bool in_range_all   = all(in_range);

        if (in_range_all == false)
            return process_overflow<Version>(x, in_range);

        simd_type xe;
        simd_int qi;

        simd_type xv        = pi2_reduction<double, Bits, Tag>
                                ::reduction_CW3(x, xe, qi);        

        simd_type pt, pc;
        approximation(x, xv, xe, pt, pc);

        simd_type ret       = simd_tancot_reconstruction<double, Bits,Tag>
                                ::eval<Version>(qi, pt, pc);
        return ret;
    };

    force_inline
    static simd_type eval_tan(const simd_type& x)
    {
        return eval_tancot<tan_tag>(x);
    }

    force_inline
    static simd_type eval_cot(const simd_type& x)
    {
        return eval_tancot<cot_tag>(x);
    }
};

}}}

#pragma warning(pop)