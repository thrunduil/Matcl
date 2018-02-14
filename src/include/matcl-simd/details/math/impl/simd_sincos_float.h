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
#include "matcl-simd/details/math/impl/simd_sincos_helpers.h"
#include "matcl-simd/details/math/impl/pi2_reduction.h"

#pragma warning(push)
#pragma warning(disable: 4127)  //conditional expression is constant

namespace matcl { namespace simd { namespace details
{

//-----------------------------------------------------------------------
//                              FLOAT
//-----------------------------------------------------------------------

struct MATCL_SIMD_EXPORT simd_sincos_table_float_data
{
    static const float poly_sin[3];
    static const float poly_cos[3];
    static const float  poly_sincos[3*4];

    static const double poly_sin_double[3];
    static const double poly_cos_double[3];
    static const double poly_sincos_double[3*2];
};

template<int Bits, class Tag>
struct simd_sincos_table_float
{
    using simd_type         = ms::simd<float, Bits, Tag>;

    force_inline
    static void eval(const simd_type& x, simd_type& p_sin, simd_type& p_cos)
    {
        p_sin               = estrin<3>(x, simd_sincos_table_float_data::poly_sin);
        p_cos               = estrin<3>(x, simd_sincos_table_float_data::poly_cos);
    };

    template<class Simd_double>
    force_inline
    static void eval_double(const Simd_double& x, Simd_double& p_sin, Simd_double& p_cos)
    {
        p_sin               = estrin<3>(x, simd_sincos_table_float_data::poly_sin_double);
        p_cos               = estrin<3>(x, simd_sincos_table_float_data::poly_cos_double);
    };
};

#if MATCL_ARCHITECTURE_HAS_SSE2
    template<>
    struct simd_sincos_table_float<128, ms::scalar_sse_tag>
    {
        using simd_type         = ms::simd<float, 128, ms::scalar_sse_tag>;
        using simd_double       = ms::simd<double, 128, ms::scalar_sse_tag>;
        using simd_sse          = ms::simd<float, 128, ms::sse_tag>;    
        using simd_sse_double   = ms::simd<double, 128, ms::sse_tag>;    

        force_inline
        static void eval(const simd_type& x, simd_type& p_sin, simd_type& p_cos)
        {
            simd_sse x_sse      = x.as_vector();
            const simd_sse* poly= reinterpret_cast<const simd_sse*>(simd_sincos_table_float_data::poly_sincos);
            simd_sse p_sincos   = estrin<3>(x_sse, poly);
            p_sin               = simd_type(p_sincos.first());
            p_cos               = simd_type(p_sincos.extract_high());
        };

        force_inline
        static void eval_double(const simd_double& x, simd_double& p_sin, simd_double& p_cos)
        {
            simd_sse_double x_sse       = x.as_vector();
            const simd_sse_double* poly = reinterpret_cast<const simd_sse_double*>
                                            (simd_sincos_table_float_data::poly_sincos_double);

            simd_sse_double p_sincos = estrin<3>(x_sse, poly);
            p_sin                    = simd_double(p_sincos.first());
            p_cos                    = simd_double(p_sincos.extract_high());
        };
    };
#endif

template<int Bits, class Tag, bool Is_scalar = is_scalar_tag<Tag>::value>
struct simd_sincos_float
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

        simd_type pc, ps;
        approximation(xv, xe, pc, ps);

        simd_type ret       = simd_sincos_reconstruction<float, Bits,Tag>
                                ::eval<Version>(q_lo, pc, ps);
        return ret;
    }

    force_inline
    static void approximation(const simd_type& xv, const simd_type& xe, simd_type& ret_cos, 
                             simd_type& ret_sin)
    {
        const simd_type one     = simd_type::one();
        const simd_type half    = simd_type(0.5f);

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
        simd_sincos_table_float<Bits, Tag>::eval(x2, p_sin, p_cos);

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

    template<int Version>
    force_inline
    static simd_type eval_sincos(const simd_type& x)
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

        simd_type pc, ps;
        approximation(xv, xe, pc, ps);

        simd_type ret       = simd_sincos_reconstruction<float, Bits, Tag>
                                ::eval<Version>(qi, pc, ps);
        return ret;
    };

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

// for scalars compute in double
template<int Bits, class Tag>
struct simd_sincos_float<Bits, Tag, true>
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

        simd_double pc, ps;
        approximation(xv, pc, ps);

        simd_double ret     = simd_sincos_reconstruction<double, Bits, Tag>
                                ::eval<Version>(q_lo, pc, ps);
        return ret;
    }

    force_inline
    static void approximation(const simd_double& xv, simd_double& ret_cos, simd_double& ret_sin)
    {
        const simd_double one   = simd_double::one();
        const simd_double half  = simd_double(0.5);

        // h(x)     = (sin(x)/x - one)/(x*x) = -1/6 + x^2/120 + ...
        // sin(x)   = x^3 * h(x) + x 
        // range    : [-pi/4, pi/4]

        // g(x)     = (cos(x) - one + 1/2*x^2)/ x^4
        // cos(x)   = x^4 * g(x) - 1/2*x^2 + 1
        // range    : [-pi/4, pi/4]

        simd_double x2      = xv * xv;
        simd_double x4      = x2 * x2;
        simd_double x2_half = half * x2;

        // eval polynomials
        simd_double p_sin, p_cos;
        simd_sincos_table_float<Bits, Tag>::eval_double(x2, p_sin, p_cos);

        // form p_sin * x3 + xv
        simd_double p_sinv  = p_sin * x2;        
        simd_double ps      = fma_f(p_sinv, xv, xv);

        // form p_cos * x4 - 1/2 * x2 + 1
        simd_double pc      = fms_f(p_cos, x4, x2_half) + one;

        ret_sin             = ps;  
        ret_cos             = pc;
    }

    template<int Version>
    force_inline
    static simd_type eval_sincos(const simd_type& x)
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

        simd_double pc, ps;
        approximation(xv, pc, ps);

        simd_double ret     = simd_sincos_reconstruction<double, Bits, Tag>
                                ::eval<Version>(qi, pc, ps);
        return ret.convert_to_float();
    };

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

#if MATCL_ARCHITECTURE_HAS_SSE2 && MATCL_ARCHITECTURE_HAS_AVX
    // compute in double
    template<>
    struct simd_sincos_float<128, sse_tag, false>
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

            simd_double pc, ps;
            approximation(xv, pc, ps);

            simd_double ret     = simd_sincos_reconstruction<double, 256, avx_tag>
                                    ::eval<Version>(q_lo, pc, ps);
            return ret;
        }

        force_inline
        static void approximation(const simd_double& xv, simd_double& ret_cos, simd_double& ret_sin)
        {
            const simd_double one   = simd_double::one();
            const simd_double half  = simd_double(0.5);

            // h(x)     = (sin(x)/x - one)/(x*x) = -1/6 + x^2/120 + ...
            // sin(x)   = x^3 * h(x) + x 
            // range    : [-pi/4, pi/4]

            // g(x)     = (cos(x) - one + 1/2*x^2)/ x^4
            // cos(x)   = x^4 * g(x) - 1/2*x^2 + 1
            // range    : [-pi/4, pi/4]

            simd_double x2      = xv * xv;
            simd_double x4      = x2 * x2;
            simd_double x2_half = half * x2;

            // eval polynomials
            simd_double p_sin, p_cos;
            simd_sincos_table_float<Bits, Tag>::eval_double(x2, p_sin, p_cos);

            // form p_sin * x3 + xv
            simd_double p_sinv  = p_sin * x2;        
            simd_double ps      = fma_f(p_sinv, xv, xv);

            // form p_cos * x4 - 1/2 * x2 + 1
            simd_double pc      = fms_f(p_cos, x4, x2_half) + one;

            ret_sin             = ps;  
            ret_cos             = pc;
        }

        template<int Version>
        force_inline
        static simd_type eval_sincos(const simd_type& x)
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

            simd_double pc, ps;
            approximation(xv, pc, ps);

            simd_double ret     = simd_sincos_reconstruction<double, 256, avx_tag>
                                    ::eval<Version>(qi, pc, ps);
            return ret.convert_to_float();
        };

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
#endif

template<int Bits, class Tag>
struct simd_sincos<float, Bits, Tag> : public simd_sincos_float<Bits, Tag>
{};

}}}

#pragma warning(pop)