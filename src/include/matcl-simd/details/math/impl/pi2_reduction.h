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

#pragma warning(push)
#pragma warning(disable: 4127)  //conditional expression is constant

namespace matcl { namespace simd { namespace details
{

template<class Val, int Bits, class Simd_tag>
struct pi2_reduction;

template<class Val, int Bits, class Simd_tag>
struct pi2_reduction_scalar;

//-----------------------------------------------------------------------
//                              DOUBLE
//-----------------------------------------------------------------------

template<int Bits, class Tag>
struct pi2_reduction<double, Bits, Tag>
{
    using simd_type         = ms::simd<double, Bits, Tag>;
    using simd_int          = typename ms::simd<int32_t, Bits, Tag>::simd_half;

    // reduction scheme valid for |x| <= 2^19 * pi/2, relative error <= 0.0804 * eps
    // for |x| <= 145897 * pi/2, relative error <= 0.0056 * eps
    force_inline
    static simd_type reduction_CW3_normalize(const simd_type& x, simd_type& xe, simd_int& qi)
    {
        using twofold       = matcl::twofold<simd_type>;

        simd_type xv_l, xe_l;
        xv_l        = reduction_CW3(x, xe_l, qi);
        twofold xr  = twofold::normalize_fast(xv_l, xe_l);

        xe          = xr.error;
        return xr.value;
    };

    // reduction scheme valid for |x| <= 2^19 * pi/2, relative error <= 0.0804 * eps
    // for |x| <= 145897 * pi/2, relative error <= 0.0056 * eps
    force_inline
    static simd_type reduction_CW3(const simd_type& x, simd_type& xe, simd_int& qi)
    {
        //---------------------------------------------------------------------------
        //                              constants
        //---------------------------------------------------------------------------

        // 2 / pi
        const simd_type two_by_pi   = simd_type(0.63661977236758134307553505349006);

        // pi/2 ~ pi2_1 + pi2_2 + pi2_3; pi2_1, pi2_2 are stored with precision 53 - 17

        // pi2_1 has 3 additional bits correct and last two bits (storing a part of pi/2)
        // are zero
        const simd_type pi2_1       = simd_type(1.57079632679233327507972717285);

        // pi2_1 has 2 additional bits correct and last two bits (storing a part of pi/2) 
        // are zero; since last two bits of pi2_1 and pi2_2 are zero, pi_1*k and pi2_2* k 
        // are exact for |k| <= 2^19 = 524,288; in this range 
        // d = |k * pi/2 - fl(k * pi/2))|/eps(k*pi/2)
        // satisfies |d| > 1.58 * 2^-21 with minimum for 204551
        // for |k| <= 145897, |d| > 1.42 * 2^-17
        const simd_type pi2_2       = simd_type(2.56334415158395578912195467147e-12);

        // total precision of pi2 is 53 + (53 - 17) * 2 + 3 + 2 = 130 = 2*53 + 21 + 3
        const simd_type pi2_3       = simd_type(1.05629990669874271124186809806e-23);        

        //---------------------------------------------------------------------------
        //                              code
        //---------------------------------------------------------------------------

        simd_type q         = round(x * two_by_pi);
        qi                  = q.convert_to_int32();
        qi                  = bitwise_and(qi, simd_int(3));

        // Cody-Waite reduction schemes

        // exact
        simd_type x_red1    = fnma_f(q, pi2_1, x);        

        // compensated summation is required to obtain exact result; if simple
        // summation is used, then total error increases by approximately 0.5 ulp
        simd_type x_red2_v  = fnma_f(q, pi2_2, x_red1);
        simd_type tmp1      = x_red1 - x_red2_v; 

        // 1 ulp error in q * pi2_3 => err(q * pi2_3) <= (q * pi2_3) * 2^-52
        // |res| > d * eps(k*pi/2) >= d * (q*pi/2) *2^-53
        // re = err(q * pi2_3) / |res| < pi2_3 / d * 4 / pi = 1.785e-17 = 0.0804 * eps
        // for |k| <= 145897, re < pi2_3 / d * 4 / pi = 1.242e-18 = 0.0056 * eps
        // this error is not compensated; in most cases this error is much lower
        simd_type x_red3_v  = fnma_f(q, pi2_3, x_red2_v);
        simd_type tmp2      = x_red2_v - x_red3_v;

        simd_type x_red2_e  = fnma_f(q, pi2_2, tmp1);
        simd_type x_red3_e  = fnma_f(q, pi2_3, tmp2);

        // xv need not be exactly rounded => normalization is not required
        simd_type xv        = x_red3_v;
        xe                  = x_red2_e + x_red3_e;

        return xv;
    };

    // reduction using Payne-Hanek algorithm
    force_inline
    static void reduction_full(const simd_type& x, const simd_type& in_range, simd_type& xv, 
                            simd_type& xe, simd_int& q_lo)
    {
        const double* ptr_x     = x.get_raw_ptr();
        const double* ptr_inr   = in_range.get_raw_ptr();

        double* ptr_xv          = xv.get_raw_ptr();
        double* ptr_xe          = xe.get_raw_ptr();
        int* ptr_qlo            = q_lo.get_raw_ptr();

        static const int vec_size   = simd_type::vector_size;

        for (int i = 0; i < vec_size; ++i)
        {
            if (ptr_inr[i] == 0)
            {
                double value, error;
                int q           = ms::reduce_pi2_ph(ptr_x[i], value, error);

                ptr_xv[i]       = value;
                ptr_xe[i]       = error;
                ptr_qlo[i]      = q;
            }
        };

        if (vec_size == 1)
        {
            xv      = simd_type(xv.first());
            xe      = simd_type(xe.first());
            q_lo    = simd_int(q_lo.first());
        };
    };
};

//-----------------------------------------------------------------------
//                              FLOAT
//-----------------------------------------------------------------------

// reduction scheme valid for |x| < 268435456 * pi/2
template<int Bits, class Tag>
struct pi2_reduction_CW_float_double
{};

template<int Bits>
struct pi2_reduction_CW_float_double<Bits, avx_tag>
{
    using Tag               = avx_tag;
    using simd_type         = ms::simd<float, Bits, Tag>;
    using simd_double       = ms::simd<double, Bits, Tag>;
    using simd_int          = ms::simd<int32_t, Bits, Tag>;    

    // reduction scheme valid for |x| < 268435456 * pi/2
    force_inline
    static simd_type eval(const simd_type& x0, simd_type& xe, simd_int& qi)
    {
        simd_double xlo     = x0.convert_low_to_double();
        simd_double xhi     = x0.convert_high_to_double();

        //---------------------------------------------------------------------------
        //                              constants
        //---------------------------------------------------------------------------

        // 2 / pi
        const simd_double two_by_pi = simd_double(0.63661977236758134307553505349006);

        // pi/2 ~ pi2_1 + pi2_2; pi2_1 is stored with precision 53 - 28
        // therefore pi2_1 * k is exact for |k| < 2^28 = 268435456;
        // for |k| < 2^28, d = |k * pi/2 - fl(k * pi/2))|/eps(k*pi/2)
        // satisfies |d| > 1.5372 * 2^-32  with minimum for 394733961

        const simd_double pi2_1     = simd_double(1.570796310901641845703125);
        const simd_double pi2_2     = simd_double(1.58932547735281960548080350e-08);

        //---------------------------------------------------------------------------
        //                              code
        //---------------------------------------------------------------------------
        simd_double q_lo    = round(xlo * two_by_pi);
        simd_double q_hi    = round(xhi * two_by_pi);
        auto qi_lo          = q_lo.convert_to_int32();
        auto qi_hi          = q_hi.convert_to_int32();
        qi                  = simd_int(qi_lo, qi_hi);
        qi                  = bitwise_and(qi, simd_int(3));

        // Cody-Waite reduction schemes

        // exact
        simd_double xlo_red1    = fnma_f(q_lo, pi2_1, xlo);
        simd_double xhi_red1    = fnma_f(q_hi, pi2_1, xhi);

        simd_double xlo_red2    = fnma_f(q_lo, pi2_2, xlo_red1);
        simd_double xhi_red2    = fnma_f(q_hi, pi2_2, xhi_red1);

        auto xv_lo              = xlo_red2.convert_to_float();
        auto xv_hi              = xhi_red2.convert_to_float();

        auto xe_lo              = (xlo_red2 - xv_lo.convert_to_double()).convert_to_float();
        auto xe_hi              = (xhi_red2 - xv_hi.convert_to_double()).convert_to_float();

        simd_type xv            = simd_type(xv_lo, xv_hi);
        xe                      = simd_type(xe_lo, xe_hi);

        return xv;
    };
};

template<>
struct pi2_reduction_CW_float_double<128, sse_tag>
{
    using Tag               = sse_tag;
    using simd_type         = ms::simd<float, 128, Tag>;
    using simd_double       = typename simd_type::simd_double_2;
    using simd_int          = ms::simd<int32_t, 128, Tag>;    

    // reduction scheme valid for |x| < 268435456 * pi/2
    force_inline
    static simd_type eval(const simd_type& x0, simd_type& xe, simd_int& qi)
    {
        simd_double x       = x0.convert_to_double();

        //---------------------------------------------------------------------------
        //                              constants
        //---------------------------------------------------------------------------

        // 2 / pi
        const simd_double two_by_pi = simd_double(0.63661977236758134307553505349006);

        // pi/2 ~ pi2_1 + pi2_2; pi2_1 is stored with precision 53 - 28
        // therefore pi2_1 * k is exact for |k| < 2^28 = 268435456;
        // for |k| < 2^28, d = |k * pi/2 - fl(k * pi/2))|/eps(k*pi/2)
        // satisfies |d| > 1.5372 * 2^-32  with minimum for 394733961

        const simd_double pi2_1     = simd_double(1.570796310901641845703125);
        const simd_double pi2_2     = simd_double(1.58932547735281960548080350e-08);

        //---------------------------------------------------------------------------
        //                              code
        //---------------------------------------------------------------------------
        simd_double q       = round(x * two_by_pi);
        qi                  = q.convert_to_int32();
        qi                  = bitwise_and(qi, simd_int(3));

        // Cody-Waite reduction schemes

        // exact
        simd_double x_red1  = fnma_f(q, pi2_1, x);
        simd_double x_red2  = fnma_f(q, pi2_2, x_red1);

        simd_type xv        = x_red2.convert_to_float();
        xe                  = (x_red2 - xv.convert_to_double()).convert_to_float();

        return xv;
    };
};

template<int Bits, class Tag>
struct pi2_reduction<float, Bits, Tag>
{
    using simd_type         = ms::simd<float, Bits, Tag>;
    using simd_int          = ms::simd<int32_t, Bits, Tag>; 

    // reduction scheme valid for |x| < 161 * pi/2
    force_inline
    static simd_type reduction_CW3(const simd_type& x, simd_type& xe, simd_int& qi)
    {
        //---------------------------------------------------------------------------
        //                              constants
        //---------------------------------------------------------------------------

        // 2 / pi
        const simd_type two_by_pi   = simd_type(0.63661977236758134307553505349006f);

        // pi/2 ~ pi2_1 + pi2_2 + pi2_3; pi2_1, pi2_2 are stored with precision 24 - 8
        // therefore pi2_1 * k and pi2_2 * k is exact for |k| <= 2^8 = 256;
        // for |k| < 161, d = |k * pi/2 - fl(k * pi/2))|/eps(k*pi/2)
        // satisfies |d| > 1.53 * 2^-8 with minimum for 137

        const simd_type pi2_1       = simd_type(1.57080078125f);
        const simd_type pi2_2       = simd_type(-4.454399459064006805419921875e-06f);
        const simd_type pi2_3       = simd_type(-5.5644315544167710640977020375430583953857421875e-11f);

        //---------------------------------------------------------------------------
        //                              code
        //---------------------------------------------------------------------------

        simd_type q         = round(x * two_by_pi);
        qi                  = q.convert_to_int32();
        qi                  = bitwise_and(qi, simd_int(3));

        // Cody-Waite reduction schemes

        // exact
        simd_type x_red1    = fnma_f(q, pi2_1, x);        

        // compensated summation is required to obtain exact result; if simple
        // summation is used, then total error increases by approximately 0.5 ulp
        simd_type x_red2_v  = fnma_f(q, pi2_2, x_red1);
        simd_type tmp1      = x_red1 - x_red2_v; 

        // 1 ulp error in q * pi2_3 => err(q * pi2_3) <= (q * pi2_3) * 2^-23
        // |res| > d * eps(k*pi/2) >= d * (q*pi/2) *2^-24
        // err(q * pi2_3) / |res| < pi2_3 / d * 4 / pi = 1.1854397e-8 = 0.09944 * eps
        // this error is not compensated; in most cases this error is much lower
        simd_type x_red3_v  = fnma_f(q, pi2_3, x_red2_v);
        simd_type tmp2      = x_red2_v - x_red3_v;

        simd_type x_red2_e  = fnma_f(q, pi2_2, tmp1);
        simd_type x_red3_e  = fnma_f(q, pi2_3, tmp2);

        // xv need not be exactly rounded => normalization is not required
        simd_type xv        = x_red3_v;
        xe                  = x_red2_e + x_red3_e;

        return xv;
    };

    // reduction scheme valid for |x| < 268435456 * pi/2
    force_inline
    static simd_type reduction_CW_double(const simd_type& x, simd_type& xe, simd_int& qi)
    {
        return pi2_reduction_CW_float_double<Bits, Tag>::eval(x, xe, qi);
    }    

    // reduction using Payne-Hanek algorithm
    force_inline
    static void reduction_full(const simd_type& x, const simd_type& in_range, simd_type& xv, 
                            simd_type& xe, simd_int& q_lo)
    {
        // 30 bits of accuracy is enough
        static const 
        int max_bits            = 30;

        const float* ptr_x      = x.get_raw_ptr();
        const float* ptr_inr    = in_range.get_raw_ptr();

        float* ptr_xv           = xv.get_raw_ptr();
        float* ptr_xe           = xe.get_raw_ptr();
        int* ptr_qlo            = q_lo.get_raw_ptr();

        static const int vec_size   = simd_type::vector_size;

        for (int i = 0; i < vec_size; ++i)
        {
            if (ptr_inr[i] == 0)
            {
                float value, error;
                int q           = ms::reduce_pi2_ph_float(ptr_x[i], max_bits, value, error);

                ptr_xv[i]       = value;
                ptr_xe[i]       = error;
                ptr_qlo[i]      = q;
            }
        };
    };
};

#if MATCL_ARCHITECTURE_HAS_SSE2 && MATCL_ARCHITECTURE_HAS_AVX

    template<int Bits>
    struct pi2_reduction<float, Bits, sse_tag>
    {
        using Tag               = sse_tag;
        using simd_type         = ms::simd<float, 128, sse_tag>;
        using simd_double       = ms::simd<double, 256, avx_tag>;
        using simd_int          = ms::simd<int32_t, 128, sse_tag>;    

        // reduction using Payne-Hanek algorithm
        force_inline
        static void reduction_full(const simd_double& x, const simd_double& in_range, 
                                    simd_double& xv, simd_int& q_lo)
        {
            // 30 bits of accuracy is enough
            static const 
            int max_bits            = 30;

            const double* ptr_x     = x.get_raw_ptr();
            const double* ptr_inr   = in_range.get_raw_ptr();

            double* ptr_xv          = xv.get_raw_ptr();
            int* ptr_qlo            = q_lo.get_raw_ptr();

            static const int vec_size   = simd_type::vector_size;

            for (int i = 0; i < vec_size; ++i)
            {
                if (ptr_inr[i] == 0)
                {
                    double value;
                    int q           = ms::reduce_pi2_ph_float(ptr_x[i], max_bits, value);

                    ptr_xv[i]       = value;
                    ptr_qlo[i]      = q;
                }
            };
        };

        // reduction scheme valid for |x| < 268435456 * pi/2
        force_inline
        static simd_double reduction_CW3(const simd_double& x, simd_int& qi)
        {
            //---------------------------------------------------------------------------
            //                              constants
            //---------------------------------------------------------------------------

            // 2 / pi
            const simd_double two_by_pi = simd_double(0.63661977236758134307553505349006);

            // pi/2 ~ pi2_1 + pi2_2; pi2_1 is stored with precision 53 - 28
            // therefore pi2_1 * k is exact for |k| < 2^28 = 268435456;
            // for |k| < 2^28, d = |k * pi/2 - fl(k * pi/2))|/eps(k*pi/2)
            // satisfies |d| > 1.5372 * 2^-32  with minimum for 394733961

            const simd_double pi2_1     = simd_double(1.570796310901641845703125);
            const simd_double pi2_2     = simd_double(1.5893254773528196054808035e-08);

            //---------------------------------------------------------------------------
            //                              code
            //---------------------------------------------------------------------------
            simd_double q       = round(x * two_by_pi);
            qi                  = q.convert_to_int32();
            qi                  = bitwise_and(qi, simd_int(3));

            // Cody-Waite reduction schemes

            // exact
            simd_double x_red1  = fnma_f(q, pi2_1, x);
            simd_double x_red2  = fnma_f(q, pi2_2, x_red1);;

            return x_red2;
        };
    };

#endif

template<int Bits, class Tag>
struct pi2_reduction_scalar<float, Bits, Tag>
{
    using simd_double       = ms::simd<double, Bits, Tag>;
    using simd_int          = ms::simd<int32_t, Bits, Tag>;    

    // reduction scheme valid for |x| < 268435456 * pi/2
    force_inline
    static simd_double reduction_CW3(const simd_double& x, simd_int& qi)
    {
        //---------------------------------------------------------------------------
        //                              constants
        //---------------------------------------------------------------------------

        // 2 / pi
        const simd_double two_by_pi = simd_double(0.63661977236758134307553505349006);

        // pi/2 ~ pi2_1 + pi2_2; pi2_1 is stored with precision 53 - 28
        // therefore pi2_1 * k is exact for |k| < 2^28 = 268435456;
        // for |k| < 2^28, d = |k * pi/2 - fl(k * pi/2))|/eps(k*pi/2)
        // satisfies |d| > 1.5372 * 2^-32  with minimum for 394733961

        const simd_double pi2_1     = simd_double(1.570796310901641845703125);
        const simd_double pi2_2     = simd_double(1.5893254773528196054808035e-08);

        //---------------------------------------------------------------------------
        //                              code
        //---------------------------------------------------------------------------
        simd_double q       = round(x * two_by_pi);
        qi                  = q.convert_to_int32();
        qi                  = bitwise_and(qi, simd_int(3));

        // Cody-Waite reduction schemes

        // exact
        simd_double x_red1  = fnma_f(q, pi2_1, x);
        simd_double x_red2  = fnma_f(q, pi2_2, x_red1);;

        return x_red2;
    };

    force_inline
    static void reduction_full(const simd_double& x, simd_double& xv, simd_int& q_lo)
    {
        // 30 bits of accuracy is enough
        static const 
        int max_bits    = 30;

        double value;
        int q           = ms::reduce_pi2_ph_float(x.first(), max_bits, value);

        xv              = simd_double(value);
        q_lo            = simd_int(q);
    };
};

}}}

#pragma warning(pop)