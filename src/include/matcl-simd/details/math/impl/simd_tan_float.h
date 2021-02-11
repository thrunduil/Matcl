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
#include "matcl-core/float/twofold.h"
#include "matcl-simd/details/math/impl/payne_hanek.inl"
#include "matcl-simd/details/math/impl/simd_tancot_helpers.h"
#include "matcl-simd/details/math/impl/pi2_reduction.h"

#pragma warning(push)
#pragma warning(disable: 4127)  //conditional expression is constant

namespace matcl { namespace simd { namespace details
{

//-----------------------------------------------------------------------
//                              FLOAT
//-----------------------------------------------------------------------

struct MATCL_SIMD_EXPORT simd_tancot_table_float_data
{
    static const float  poly_tan[7];
    static const double poly_tan_double[7];
};

template<int Bits, class Tag>
struct simd_tancot_table_float
{
    using simd_type         = ms::simd<float, Bits, Tag>;
    
    force_inline
    static void eval(const simd_type& x, simd_type& p_tan)
    {
        p_tan   = estrin<7>(x, simd_tancot_table_float_data::poly_tan);
    };

    template<class Simd_double>
    force_inline
    static void eval_double(const Simd_double& x, Simd_double& p_tan)
    {
        p_tan   = estrin<7>(x, simd_tancot_table_float_data::poly_tan_double);
    };
};

#if MATCL_ARCHITECTURE_HAS_SSE2
    template<>
    struct simd_tancot_table_float<128, ms::scalar_sse_tag>
    {
        using simd_type         = ms::simd<float, 128, ms::scalar_sse_tag>;
        using simd_double       = ms::simd<double, 128, ms::scalar_sse_tag>;
        using simd_sse          = ms::simd<float, 128, ms::sse_tag>;    
        using simd_sse_double   = ms::simd<double, 128, ms::sse_tag>;   

        force_inline
        static void eval_double(const simd_double& x, simd_double& p_tan)
        {
            simd_sse_double x_sse   = x.as_vector();
            simd_sse_double val     = estrin<7>(x_sse, simd_tancot_table_float_data::poly_tan_double);
            p_tan                   = simd_double(val.first());
        };
    };
#endif

template<int Bits, class Tag, bool Is_scalar = is_scalar_tag<Tag>::value>
struct simd_tancot_float
{
    using simd_type         = ms::simd<float, Bits, Tag>;
    using simd_int          = ms::simd<int32_t, Bits, Tag>;    

    template<int Version>
    static simd_type process_overflow(const simd_type& x)
    {
        const simd_type nan_v   = simd_type(std::numeric_limits<float>::quiet_NaN());

        simd_type fin   = ms::is_finite(x);
        simd_type ret   = eval_with_full_reduction<Version>(x);
        ret             = ms::if_then_else(fin, ret, nan_v);

        return ret;
    }

    template<int Version>
    force_inline
    static simd_type eval_with_full_reduction(const simd_type& x)
    {
        // max_q * pi/2
        const simd_type max_x   = simd_type(421657428.26631f);

        simd_type in_range      = gt(max_x, abs(x));
        bool in_range_all       = all(in_range);
        bool in_range_any       = any(in_range);

        simd_type xv, xe;
        simd_int q_lo;

        if (in_range_any == true)
            xv  = pi2_reduction<float, Bits, Tag>::reduction_CW_double(x, xe, q_lo); 

        if (in_range_all == false)
            pi2_reduction<float, Bits, Tag>::reduction_full(x, in_range, xv, xe, q_lo);

        simd_type pt, pc;
        approximation(x, xv, xe, pt, pc);

        simd_type ret       = simd_tancot_reconstruction<float, Bits,Tag>
                                ::eval<Version>(q_lo, pt, pc);
        return ret;
    }

    force_inline
    static void approximation(const simd_type& x, const simd_type& xv, const simd_type& xe, 
                            simd_type& ret_tan, simd_type& ret_cot)
    {
        const simd_type zero    = simd_type::zero();
        simd_type sign_inf      = ms::copysign(simd_type(std::numeric_limits<float>::infinity()), x);

        // h(x)     = (tan(sqrt(x))/sqrt(x) - 1)/x - 1/2
        // tan(x)   = x^3 * (h(x^2) + 1/2) + x 
        // range    : [0, (pi/4)^2]

        using twofold       = matcl::twofold<simd_type>;

        simd_type x2t       = xv * xv;
        simd_type x3t       = x2t * xv;

        // eval polynomials
        simd_type p_tan;
        simd_tancot_table_float<Bits, Tag>::eval(x2t, p_tan);

        p_tan               = p_tan + simd_type(0.5);
        
        // p_tan * xv^3 must be computed with higher precision
        twofold pt_x3       = twofold_mult_f(p_tan, x3t);

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
    force_inline
    static simd_type eval_tancot(const simd_type& x)
    {        
        // max_q * pi/2
        const simd_type max_x   = simd_type(252.898f);

        simd_type in_range  = gt(max_x, abs(x));
        bool in_range_all   = all(in_range);

        if (in_range_all == false)
            return process_overflow<Version>(x);

        simd_type xe;
        simd_int qi;

        simd_type xv        = pi2_reduction<float, Bits, Tag>::reduction_CW3(x, xe, qi);

        simd_type pt, pc;
        approximation(x, xv, xe, pt, pc);

        simd_type ret       = simd_tancot_reconstruction<float, Bits, Tag>
                                ::eval<Version>(qi, pt, pc);
        return ret;
    };

    static simd_type eval_tan(const simd_type& x)
    {
        return eval_tancot<tan_tag>(x);
    }

    static simd_type eval_cot(const simd_type& x)
    {
        return eval_tancot<cot_tag>(x);
    }
};

// for scalars compute in double
template<int Bits, class Tag>
struct simd_tancot_float<Bits, Tag, true>
{
    using simd_type         = ms::simd<float, Bits, Tag>;
    using simd_double       = ms::simd<double, Bits, Tag>;
    using simd_int          = ms::simd<int32_t, Bits, Tag>;    

    template<int Version>
    static simd_type process_overflow(const simd_double& x)
    {
        const simd_double nan_v   = simd_double(std::numeric_limits<double>::quiet_NaN());

        simd_double fin = ms::is_finite(x);
        simd_double ret = eval_with_full_reduction<Version>(x);
        ret             = ms::if_then_else(fin, ret, nan_v);

        return ret.convert_to_float();
    }

    template<int Version>
    force_inline
    static simd_double eval_with_full_reduction(const simd_double& x)
    {
        simd_double xv;
        simd_int q_lo;

        pi2_reduction_scalar<float, Bits, Tag>::reduction_full(x, xv, q_lo);

        simd_double pt, pc;
        approximation(x, xv, pt, pc);

        simd_double ret     = simd_tancot_reconstruction<double, Bits, Tag>
                                ::eval<Version>(q_lo, pt, pc);
        return ret;
    }

    force_inline
    static void approximation(const simd_double& x, const simd_double& xv, 
                                simd_double& ret_tan, simd_double& ret_cot)
    {
        const simd_double zero  = simd_double::zero();
        simd_double sign_inf    = ms::copysign(simd_double(std::numeric_limits<double>::infinity()), x);

        // h(x)     = (tan(sqrt(x))/sqrt(x) - 1)/x - 1/2
        // tan(x)   = x^3 * (h(x^2) + 1/2) + x 
        // range    : [0, (pi/4)^2]

        simd_double x2      = xv * xv;

        // eval polynomials
        simd_double p_tan;
        simd_tancot_table_float<Bits, Tag>::eval_double(x2, p_tan);

        p_tan               = p_tan + simd_double(0.5);
        
        simd_double pt_x2   = p_tan * x2;
        simd_double pt_t    = fma_f(pt_x2, xv, xv);

        ret_tan             = pt_t;
        ret_cot             = simd_double::one() / pt_t;
        ret_cot             = if_then_else(neq(x, zero), ret_cot, sign_inf);
    }

    template<int Version>
    force_inline
    static simd_type eval_tancot(const simd_type& x)
    {        
        // max_q * pi/2
        const simd_double max_x = simd_double(421657428.26631);

        simd_double xd          = x.convert_to_double();

        simd_double in_range    = gt(max_x, abs(xd));
        bool in_range_all       = all(in_range);

        if (in_range_all == false)
            return process_overflow<Version>(xd);

        simd_int qi;
        simd_double xv      = pi2_reduction_scalar<float, Bits, Tag>::reduction_CW3(xd, qi);

        simd_double pt, pc;
        approximation(xd, xv, pt, pc);

        simd_double ret     = simd_tancot_reconstruction<double, Bits, Tag>
                                ::eval<Version>(qi, pt, pc);
        return ret.convert_to_float();
    };

    static simd_type eval_tan(const simd_type& x)
    {
        return eval_tancot<tan_tag>(x);
    }

    static simd_type eval_cot(const simd_type& x)
    {
        return eval_tancot<cot_tag>(x);
    }
};

#if MATCL_ARCHITECTURE_HAS_SSE2 && MATCL_ARCHITECTURE_HAS_AVX
    // compute in double
    template<>
    struct simd_tancot_float<128, sse_tag, false>
    {
        static const int Bits   = 128;

        using Tag               = sse_tag;
        using simd_type         = ms::simd<float, 128, sse_tag>;
        using simd_double       = ms::simd<double, 256, avx_tag>;
        using simd_int          = ms::simd<int32_t, 128, sse_tag>;    

        template<int Version>
        static simd_type process_overflow(const simd_double& x, const simd_double& in_range)
        {
            const simd_double nan_v   = simd_double(std::numeric_limits<double>::quiet_NaN());

            simd_double fin = ms::is_finite(x);
            simd_double ret = eval_with_full_reduction<Version>(x, in_range);
            ret             = ms::if_then_else(fin, ret, nan_v);

            return ret.convert_to_float();
        }

        template<int Version>
        force_inline
        static simd_double eval_with_full_reduction(const simd_double& x, const simd_double& in_range)
        {
            bool in_range_any       = any(in_range);

            simd_double xv;
            simd_int q_lo;

            if (in_range_any == true)
                xv  = pi2_reduction<float, Bits, Tag>::reduction_CW3(x, q_lo); 

            pi2_reduction<float, Bits, Tag>::reduction_full(x, in_range, xv, q_lo);

            simd_double pt, pc;
            approximation(x, xv, pt, pc);

            simd_double ret     = simd_tancot_reconstruction<double, 256, avx_tag>
                                    ::eval<Version>(q_lo, pt, pc);
            return ret;
        }

        force_inline
        static void approximation(const simd_double& x, const simd_double& xv, 
                                    simd_double& ret_tan, simd_double& ret_cot)
        {
            const simd_double zero  = simd_double::zero();
            simd_double sign_inf    = ms::copysign(simd_double(std::numeric_limits<double>::infinity()), x);

            // h(x)     = (tan(sqrt(x))/sqrt(x) - 1)/x - 1/2
            // tan(x)   = x^3 * (h(x^2) + 1/2) + x 
            // range    : [0, (pi/4)^2]

            simd_double x2      = xv * xv;

            // eval polynomials
            simd_double p_tan;
            simd_tancot_table_float<Bits, Tag>::eval_double(x2, p_tan);

            p_tan               = p_tan + simd_double(0.5);
        
            simd_double pt_x2   = p_tan * x2;
            simd_double pt_t    = fma_f(pt_x2, xv, xv);

            ret_tan             = pt_t;
            ret_cot             = simd_double::one() / pt_t;
            ret_cot             = if_then_else(neq(x, zero), ret_cot, sign_inf);
        }        

        template<int Version>
        force_inline
        static simd_type eval_tancot(const simd_type& x)
        {        
            // max_q * pi/2
            const simd_double max_x = simd_double(421657428.26631);

            simd_double xd          = x.convert_to_double();

            simd_double in_range    = gt(max_x, abs(xd));
            bool in_range_all       = all(in_range);

            if (in_range_all == false)
                return process_overflow<Version>(xd, in_range);

            simd_int qi;
            simd_double xv      = pi2_reduction<float, Bits, Tag>::reduction_CW3(xd, qi);

            simd_double pt, pc;
            approximation(xd, xv, pt, pc);

            simd_double ret     = simd_tancot_reconstruction<double, 256, avx_tag>
                                    ::eval<Version>(qi, pt, pc);
            return ret.convert_to_float();
        };

        static simd_type eval_tan(const simd_type& x)
        {
            return eval_tancot<tan_tag>(x);
        }

        static simd_type eval_cot(const simd_type& x)
        {
            return eval_tancot<cot_tag>(x);
        }
    };
#endif

template<int Bits, class Tag>
struct simd_tancot<float, Bits, Tag> : public simd_tancot_float<Bits, Tag>
{};

}}}

#pragma warning(pop)