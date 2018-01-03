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

#include "matcl-simd/details/math/impl/payne_hanek.inl"

namespace matcl { namespace simd { namespace details
{

// 1376 bits of 2/pi stored in 32-bit integers, each with 32 bits
const uint32_t two_by_pi_table::table[] = 
{
    0xa2f9836e, 0x4e441529, 0xfc2757d1, 0xf534ddc0, 0xdb629599, 0x3c439041, 0xfe5163ab, 0xdebbc561, 
    0xb7246e3a, 0x424dd2e0, 0x06492eea, 0x09d1921c, 0xfe1deb1c, 0xb129a73e, 0xe88235f5, 0x2ebb4484, 
    0xe99c7026, 0xb45f7e41, 0x3991d639, 0x835339f4, 0x9c845f8b, 0xbdf9283b, 0x1ff897ff, 0xde05980f, 
    0xef2f118b, 0x5a0a6d1f, 0x6d367ecf, 0x27cb09b7, 0x4f463f66, 0x9e5fea2d, 0x7527bac7, 0xebe5f17b, 
    0x3d0739f7, 0x8a5292ea, 0x6bfb5fb1, 0x1f8d5d08, 0x56033046, 0xfc7b6bab, 0xf0cfbc20, 0x9af4361d, 
    0xa9e39161, 0x5ee61b08, 0x6599855f, 
};

static const int max_blocks         = 10;
static const int max_blocks_float   = 5;

template<int Blocks, int Req_m>
struct reduce_pi2_ph_impl
{
    static const uint64_t pi2       = 0xC90FDAA22168C235;

    static int eval(double x, double& value, double& error)
    {
        payne_hanek_double<Blocks> ph;

        bool signed_x;
        uint64_t sign_x;
        int exp_x;

        ph.decompose(x, signed_x, sign_x, exp_x);

        bool frac_sign;
        uint64_t frac, res_frac;
        int64_t frac_exp, res_exp;
        int quadrant;

        quadrant        = ph.eval<Req_m>(signed_x, sign_x, exp_x, 
                                details::two_by_pi_table::table, -1, frac_sign, frac, frac_exp);

        int num_bits    = ph.number_correct_bits();

        if (num_bits < 64 && Blocks < max_blocks)
            return reduce_pi2_ph_impl<std::min(Blocks + 2, max_blocks), Req_m>::eval(x, value, error);

        ph.mult(frac, frac_exp, pi2, -63, res_frac, res_exp);
        ph.result_as_double(frac_sign, res_frac, res_exp, value, error);
        return quadrant;
    }
};

template<int Req_m>
struct reduce_pi2_ph_impl<max_blocks, Req_m>
{
    static const uint64_t pi2       = 0xC90FDAA22168C235;
    static const int Blocks         = max_blocks;
   
    static int eval(double x, double& value, double& error)
    {
        payne_hanek_double<Blocks> ph;

        bool signed_x;
        uint64_t sign_x;
        int exp_x;

        ph.decompose(x, signed_x, sign_x, exp_x);

        bool frac_sign;
        uint64_t frac, res_frac;
        int64_t frac_exp, res_exp;
        int quadrant;

        quadrant        = ph.eval<Req_m>(signed_x, sign_x, exp_x, 
                            details::two_by_pi_table::table, -1, frac_sign, frac, frac_exp);

        ph.mult(frac, frac_exp, pi2, -63, res_frac, res_exp);
        ph.result_as_double(frac_sign, res_frac, res_exp, value, error);
        return quadrant;
    }
};

template<int Blocks, int Req_m>
struct reduce_pi2_ph_float_impl
{
    static const uint64_t pi2       = 0xC90FDAA22168C235;

    static int eval(float x, int max_bits, float& value, float& error)
    {
        payne_hanek_float<Blocks> ph;

        bool signed_x;
        uint32_t sign_x;
        int exp_x;

        ph.decompose(x, signed_x, sign_x, exp_x);

        bool frac_sign;
        uint64_t frac, res_frac;
        int64_t frac_exp, res_exp;
        int quadrant;

        quadrant        = ph.eval<Req_m>(signed_x, sign_x, exp_x, details::two_by_pi_table::table, -1, 
                                        frac_sign, frac, frac_exp);
        
        int num_bits    = ph.number_correct_bits();

        if (num_bits < max_bits && Blocks < max_blocks_float)
        {
            return reduce_pi2_ph_float_impl<std::min(Blocks + 1, max_blocks_float), Req_m>
                        ::eval(x, max_bits, value, error);
        };

        ph.mult(frac, frac_exp, pi2, -63, res_frac, res_exp);
        ph.result_as_float(frac_sign, res_frac, res_exp, value, error);

        return quadrant;
    }

    static int eval_double(double x, int max_bits, double& value)
    {
        payne_hanek_float<Blocks> ph;

        bool signed_x;
        uint32_t sign_x;
        int exp_x;

        ph.decompose(x, signed_x, sign_x, exp_x);

        bool frac_sign;
        uint64_t frac, res_frac;
        int64_t frac_exp, res_exp;
        int quadrant;
        
        quadrant        = ph.eval<Req_m>(signed_x, sign_x, exp_x, details::two_by_pi_table::table, -1, 
                                        frac_sign, frac, frac_exp);

        int num_bits    = ph.number_correct_bits();

        if (num_bits < max_bits && Blocks < max_blocks_float)
        {
            return reduce_pi2_ph_float_impl<std::min(Blocks + 1, max_blocks_float), Req_m>
                        ::eval_double(x, max_bits, value);
        };

        ph.mult(frac, frac_exp, pi2, -63, res_frac, res_exp);
        ph.result_as_double(frac_sign, res_frac, res_exp, value);

        return quadrant;        
    }
};

template<int Req_m>
struct reduce_pi2_ph_float_impl<max_blocks_float, Req_m>
{
    static const uint64_t pi2       = 0xC90FDAA22168C235;
    static const int Blocks         = max_blocks_float;

    static int eval(float x, int max_bits, float& value, float& error)
    {
        (void)max_bits;

        payne_hanek_float<Blocks> ph;

        bool signed_x;
        uint32_t sign_x;
        int exp_x;

        ph.decompose(x, signed_x, sign_x, exp_x);

        bool frac_sign;
        uint64_t frac, res_frac;
        int64_t frac_exp, res_exp;
        int quadrant;

        quadrant        = ph.eval<Req_m>(signed_x, sign_x, exp_x, details::two_by_pi_table::table, -1, 
                                        frac_sign, frac, frac_exp);
        
        ph.mult(frac, frac_exp, pi2, -63, res_frac, res_exp);
        ph.result_as_float(frac_sign, res_frac, res_exp, value, error);
        return quadrant;
    }

    static int eval_double(double x, int max_bits, double& value)
    {
        (void)max_bits;

        payne_hanek_float<Blocks> ph;

        bool signed_x;
        uint32_t sign_x;
        int exp_x;

        ph.decompose(x, signed_x, sign_x, exp_x);

        bool frac_sign;
        uint64_t frac, res_frac;
        int64_t frac_exp, res_exp;
        int quadrant;
        
        quadrant        = ph.eval<Req_m>(signed_x, sign_x, exp_x, details::two_by_pi_table::table, -1, 
                                        frac_sign, frac, frac_exp);

        ph.mult(frac, frac_exp, pi2, -63, res_frac, res_exp);
        ph.result_as_double(frac_sign, res_frac, res_exp, value);

        return quadrant;        
    }
};

}}}

namespace matcl
{

int simd::reduce_pi2_ph(double x, double& value, double& error)
{
    return details::reduce_pi2_ph_impl<5, 2>::eval(x, value, error);
};

int simd::reduce_pi2_ph_float(float x, int max_bits, float& value, float& error)
{
    return details::reduce_pi2_ph_float_impl<3, 2>::eval(x, max_bits, value, error);
};

int simd::reduce_pi2_ph_float(double x, int max_bits, double& value)
{
    return details::reduce_pi2_ph_float_impl<3, 2>::eval_double(x, max_bits, value);
};

}
