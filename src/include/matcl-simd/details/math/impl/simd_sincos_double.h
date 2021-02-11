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
#include "matcl-simd/details/math/impl/simd_sincos_helpers.h"
#include "matcl-simd/details/math/impl/pi2_reduction.h"

#pragma warning(push)
#pragma warning(disable: 4127)  //conditional expression is constant

namespace matcl { namespace simd { namespace details
{

//-----------------------------------------------------------------------
//                              DOUBLE
//-----------------------------------------------------------------------

struct MATCL_SIMD_EXPORT simd_sincos_table_double_data
{
    static const double poly_sin[6];
    static const double poly_cos[6];
    static const double poly_sincos[6*2];
};

template<int Bits, class Tag>
struct simd_sincos_table_double
{
    using simd_type         = ms::simd<double, Bits, Tag>;

    force_inline
    static void eval(const simd_type& x, simd_type& p_sin, simd_type& p_cos)
    {
        p_sin               = estrin<6>(x, simd_sincos_table_double_data::poly_sin);
        p_cos               = estrin<6>(x, simd_sincos_table_double_data::poly_cos);
    };
};

#if MATCL_ARCHITECTURE_HAS_SSE2
    template<>
    struct simd_sincos_table_double<128, ms::scalar_sse_tag>
    {
        using simd_type         = ms::simd<double, 128, ms::scalar_sse_tag>;
        using simd_sse          = ms::simd<double, 128, ms::sse_tag>;    

        force_inline
        static void eval(const simd_type& x, simd_type& p_sin, simd_type& p_cos)
        {
            simd_sse x_sse      = x.as_vector();
            const simd_sse* poly= reinterpret_cast<const simd_sse*>(simd_sincos_table_double_data::poly_sincos);
            simd_sse p_sincos   = estrin<6>(x_sse, poly);
            p_sin               = simd_type(p_sincos.first());
            p_cos               = simd_type(p_sincos.extract_high());
        };
    };
#endif

template<int Bits, class Tag>
struct simd_sincos<double, Bits, Tag>
{
    using simd_type         = ms::simd<double, Bits, Tag>;
    using simd_int          = typename ms::simd<int32_t, Bits, Tag>::simd_half;
    using simd_int64        = ms::simd<int64_t, Bits, Tag>;    
    using twofold_type      = twofold<simd_type>;

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
    static simd_type eval_sincos(const simd_type& x)
    {        
        // 2^19 * pi/2
        const simd_type max_x   = simd_type(823549.6);

        simd_type in_range  = gt(max_x, abs(x));
        bool in_range_all   = all(in_range);

        if (in_range_all == false)
            return process_overflow<Version>(x, in_range);

        //-----------------------------------------------------------------
        //                      reduction
        //-----------------------------------------------------------------
        simd_type xe;
        simd_int qi;

        simd_type xv        = pi2_reduction<double, Bits, Tag>::reduction_CW3(x, xe, qi);        

        //-----------------------------------------------------------------
        //                      approximation
        //-----------------------------------------------------------------
        simd_type pc, ps;
        approximation(xv, xe, pc, ps);

        //-----------------------------------------------------------------
        //                      reconstruction
        //-----------------------------------------------------------------       
        simd_type ret       = simd_sincos_reconstruction<double, Bits,Tag>
                                ::eval<Version>(qi, pc, ps);
        return ret;
    };

    template<int Version>
    force_inline
    static simd_type eval_with_full_reduction(const simd_type& x, const simd_type& in_range)
    {
        //-----------------------------------------------------------------
        //                      reduction
        //-----------------------------------------------------------------
        simd_type xv, xe;
        simd_int q_lo;

        if (any(in_range) == true)
            xv  = pi2_reduction<double, Bits, Tag>::reduction_CW3(x, xe, q_lo); 

        pi2_reduction<double, Bits, Tag>::reduction_full(x, in_range, xv, xe, q_lo);

        //-----------------------------------------------------------------
        //                      approximation
        //-----------------------------------------------------------------
        simd_type pc, ps;
        approximation(xv, xe, pc, ps);

        //-----------------------------------------------------------------
        //                      reconstruction
        //-----------------------------------------------------------------
        simd_type ret       = simd_sincos_reconstruction<double, Bits, Tag>
                                ::eval<Version>(q_lo, pc, ps);
        return ret;
    }

    force_inline
    static void approximation(const simd_type& xv, const simd_type& xe, simd_type& ret_cos, 
                             simd_type& ret_sin)
    {
        const simd_type one     = simd_type::one();
        const simd_type half    = simd_type(0.5);

        // h(x)     = (sin(x)/x - one)/(x*x) = -1/6 + x^2/120 + ...
        // sin(x)   = x^3 * h(x) + x 
        //          = (xv + xe)^3 * h(xv + xe) + xv + xe
        //          ~ xv^3 * h(xv) + cos(x) * xe + xv
        // range    : [-pi/4, pi/4]

        // g(x)     = (cos(x) - one + 1/2*x^2)/ x^4
        // cos(x)   = x^4 * g(x) - 1/2*x^2 + 1
        //          ~ xv^4 * g(xv) - 1/2*xv^2 + 1 - sin(x) * xe
        //          ~ xv^4 * g(xv) - 1/2 * xv^2 - xv * xe + 1
        // range    : [-pi/4, pi/4]

        simd_type x2        = xv * xv;
        simd_type x4        = x2 * x2;
        simd_type x2_half   = half * x2;

        // eval polynomials
        simd_type p_sin, p_cos;
        simd_sincos_table_double<Bits, Tag>::eval(x2, p_sin, p_cos);

        // later p_sinv will be multiplied by xv in order to form p_sinv * xv^3
        // this is just an optimization saving one multiplication
        simd_type p_sinv    = p_sin * x2;        

        // do not use twofold_mult, we want fast multiply, not fma_a
        simd_type x2e       = fms_f(xv, xv, x2);
        simd_type xve       = xv * xe;

        // error compensations
        simd_type p_sine    = fnma_f(x2_half, xe, xe);
        simd_type p_cose    = fma_f(half, x2e, xve);

        // form p_sin * x3 + cos(x) * xe + xv ~ p_sin * x3 + (1 - 1/2 x2) * xe
        simd_type ps        = fma_f(p_sinv, xv, p_sine);
        ps                  = ps + xv;

        // form p_cos * x4 - (1/2 * x2e + xv * xe), where 1/2 * x2e is error compensation
        // for the term 1/2 * x2
        simd_type pc        = fms_f(p_cos, x4, p_cose);
        pc                  = (pc - x2_half) + one;

        ret_sin             = ps;  
        ret_cos             = pc;
    }

    force_inline
    static simd_type eval_sin(const simd_type& x)
    {
        return eval_sincos<sin_tag>(x);
    }

    force_inline
    static simd_type eval_cos(const simd_type& x)
    {
        return eval_sincos<cos_tag>(x);
    }
};

}}}

#pragma warning(pop)