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

#include <stdint.h>
#include "matcl-simd/config.h"

namespace matcl { namespace simd { namespace details
{

struct MATCL_SIMD_EXPORT two_by_pi_table
{
    static const uint32_t table[43];
};

}}};

namespace matcl { namespace simd
{

// decompose a product x * C, where C is a constant as
//      x * C = k + frac                                                        (1)
// where k is an integer, |frac| <= 1/2, using the Payne-Hanek algorithm;
// the integral part is returned modulo 2^req_m, i.e. ip = k mod (2 ^ req_m), 
// where 0 <= ip < 2 ^ req_m; 
// 
// The constant C is represented as C = sig_c * 2 ^ exp_c, where 1 <= sig_c < 2;
// bits of the significant sig_c are stored in 32-bit blocks;
// at least Num_blocks + floor((1024 + sig_c - 51)/32) blocks are required
//     
// the fractional part is represented as as: frac_app = value + error + rem
// where frac_app is the computed fractional part satisfying 
//      |frac_app - frac| < 2^(-p)                                              (2)
//      p  >= Num_blocks * 32 - std::min(exp_x + exp_c + 2,  req_m + 53 + 31)   (3)
//      exp_x is the exponent of x, i.e. |x| = sig_x * 2 ^ exp_x, 1 <= sig_x < 2
//
//      frac_app = frac_app_tr + rem, where frac_app_tr stores 64 bits of frac_app
//      value    = fl(frac_app_tr) is 53-bit floating point value (rounded to nearest),
//      error    = frac_app_tr - value
//
// Template arguments:
//  Num_blocks  - number of blocks storing the constant C used to calculate the fractional
//                part frac, the higher value the lower absolute error given by (2-3)
template<int Num_blocks>
class payne_hanek_double
{
    private:
        // number of full 32-bit blocks storing the fractional part; 
        // minimal value = 2, this gives at least 2*32 + 1 bits of the result
        //          (if available)
        // maximal value = 4, larger values will not increase precision of the
        //          result, maximum precision is 2 * 53 bits
        static const int max_full_blocks = 2;

        static const int number_bits    = 32;

        static_assert(Num_blocks >= max_full_blocks, "Num_blocks too small");

    private:
        int             n_bits;
        int             m_first_nz;
        
    public:
        // represent a double precision floating point number x as x = sign * frac_x * 2 ^ exp_x,
        // where sign = -1 if signed_x = true and 1 otherwise, 1 <= frac_x < 2, on return ifrac_x 
        // contains frac_x represented as 64-bit integer (i.e. ifrac_x = uint64_t(frac_x * 2^52))
        static void     decompose(const double& x, bool& signed_x, uint64_t& ifrac_x, int& exp_x);

        // compute the decomposition (1), x - a scalar, such that x * C >= 1; contant - array 
        // of 32-bit blocks storing the significatant of C; const_exponent - exponent of the
        // constant C (i.e. exp_c); the fractional part is represented as 
        // sign * frac * 2^frac_exponent, where sign = -1 if frac_signed = true and 1 otherwise,
        // the most significal bit of frac is 1 unless frac == 0;
        // returned value - ip = k mod (2 ^ Req_m), where  0 <= ip < 2 ^ Req_m; Req_m - integer,
        // Req_m >= 0; the scalar x is represented as x = sign_x * frac_x * 2 ^ exp_x, where 
        // sign_x = -1 if is_signed = true and sign_x = 1 otherwise; 1 <= frac_x < 2, 
        // ifrac_x = uint64_t(frac_x * 2^52); fuction decompose returns such representation
        template<int Req_m>
        int             eval(bool is_signed, uint64_t ifrac_x, int exp_x, const uint32_t* constant, 
                            int const_exponent, bool& frac_signed, uint64_t& frac, int64_t& frac_exponent);

        // multiply two floating point numbers represented as xi = significant_i * 2^exp_i, 
        // i = 1,2, significant_1 is 64-bit unsigned integer, and return y = x1 * x2, 
        // y = significant_res * 2 ^ exp_res (128-bit result is truncated to 64 bits)
        static void     mult(uint64_t significant_1, int64_t exp_1, uint64_t significant_2, 
                            int64_t exp_2, uint64_t& significant_res, int64_t& exp_res);

        // convert a number represented as val = sign * frac * 2^frac_exp, where frac is 64-bit
        // integer to floating point representation val = value + error, where value = fl(val),
        // error = val - value
        static void     result_as_double(bool frac_sign, uint64_t frac, int64_t frac_exp,
                            double& value, double& error);

        // return the maximum number of correct bits p of frac stored in value and error arguments
        // of the eval function; actual number of correct bits is equal to min(p, 64);
        // eval function must be called first
        int             number_correct_bits() const;

    private:
        void            mult_scal(uint64_t xlo, uint64_t xhi, const uint32_t* blocks, uint32_t* xC_prod);

        void            extract_result(int frac_bit_0, uint32_t* xC_prod, int& integral_part, 
                            bool& frac_sign, uint64_t& frac, int64_t& frac_exp);
        
        static double   pow2ki_sign(int64_t exp, int64_t sign);
};

// decompose a product x * C, where C is a constant as
//      x * C = k + frac                                                        (1)
// where k is an integer, |frac| <= 1/2, using the Payne-Hanek algorithm;
// the integral part is returned modulo 2^req_m, i.e. ip = k mod (2 ^ req_m), 
// where 0 <= ip < 2 ^ req_m; 
// 
// The constant C is represented as C = sig_c * 2 ^ exp_c, where 1 <= sig_c < 2;
// bits of the significant sig_c are stored in 32-bit blocks;
// at least Num_blocks + floor((128 + sig_c - 22)/32) blocks are required
//     
// the fractional part is represented as as: frac_app = value + error + rem
// where frac_app is the computed fractional part satisfying 
//      |frac_app - frac| < 2^(-p)                                              (2)
//      p  >= Num_blocks * 32 - std::min(exp_x + exp_c + 2,  req_m + 24 + 31)   (3)
//      exp_x is the exponent of x, i.e. |x| = sig_x * 2 ^ exp_x, 1 <= sig_x < 2
//
//      frac_app = frac_app_tr + rem, where frac_app_tr stores 64 bits of frac_app
//      value    = fl(frac_app_tr) is 24-bit floating point value (rounded to nearest),
//      error    = frac_app_tr - value
//
// Template arguments:
//  Num_blocks  - number of blocks storing the constant C used to calculate the fractional
//                part frac, the higher value the lower absolute error given by (2-3)
template<int Num_blocks>
class payne_hanek_float
{
    private:
        // number of full 32-bit blocks storing the fractional part; 
        // minimal value = 2, this gives at least 2*32 + 1 bits of the result
        //          (if available)
        static const int max_full_blocks = 2;

        static const int number_bits    = 32;

        static_assert(Num_blocks >= max_full_blocks, "Num_blocks too small");

    private:
        int             n_bits;
        int             m_first_nz;
        
    public:
        // represent a single precision floating point number x as x = sign * frac_x * 2 ^ exp_x
        // where sign = -1 if signed_x = true and 1 otherwise; 1 <= frac_x < 2; on return ifrac_x 
        // contains frac_x represented as 32-bit integer (i.e. ifrac_x = uint32_t(frac_x * 2^23))
        static void     decompose(const float& x, bool& signed_x, uint32_t& ifrac_x, int& exp_x);

        // represent a single precision floating point number x represented by the double type
        // as x = sign * frac_x * 2 ^ exp_x, where sign = -1 if signed_x = true and 1 otherwise;
        // 1 <= frac_x < 2; on return ifrac_x  contains frac_x represented as 32-bit integer 
        // (i.e. ifrac_x = uint32_t(frac_x * 2^23))
        static void     decompose(const double& x, bool& signed_x, uint32_t& ifrac_x, int& exp_x);        

        // compute the decomposition (1), x - a scalar, such that x * C >= 1; contant - array 
        // of 32-bit blocks storing the significatant of C; const_exponent - exponent of the
        // constant C (i.e. exp_c); the fractional part is represented as 
        // sign * frac * 2^frac_exponent, where sign = -1 if frac_signed = true and 1 otherwise,
        // the most significal bit of frac is 1 unless frac == 0;
        // returned value - ip = k mod (2 ^ Req_m), where  0 <= ip < 2 ^ Req_m; Req_m - integer,
        // Req_m >= 0; the scalar x is represented as x = sign_x * frac_x * 2 ^ exp_x, where
        // sign_x = -1 if is_signed = true and sign_x = 1 otherwise; 1 <= frac_x < 2, 
        // ifrac_x = uint32_t(frac_x * 2^23); fuction decompose returns such representation
        template<int Req_m>
        int             eval(bool is_signed, uint32_t ifrac_x, int exp_x, const uint32_t* constant, 
                            int const_exponent, bool& frac_signed, uint64_t& frac, int64_t& frac_exponent);

        // multiply two floating point numbers represented as xi = significant_i * 2^exp_i, 
        // i = 1,2, significant_1 is 64-bit unsigned integer, and return y = x1 * x2, 
        // y = significant_res * 2 ^ exp_res (128-bit result is truncated to 64 bits)
        static void     mult(uint64_t significant_1, int64_t exp_1, uint64_t significant_2, 
                            int64_t exp_2, uint64_t& significant_res, int64_t& exp_res);

        // return the maximum number of correct bits p of frac stored in value and error arguments
        // of the eval function; actual number of correct bits is equal to min(p, K), where K = 64
        // when result is represented as 64-bit integer, K = 53 when result is represented as double,
        // and K = 48, when result is represented as two single precision numbers
        // eval function must be called first
        int             number_correct_bits() const;

        // convert a number represented as val = sign * frac * 2^frac_exp, where frac is 64-bit
        // integer to floating point representation: value = double(val)
        static void     result_as_double(bool frac_sign, uint64_t frac, int64_t frac_exp,
                            double& value);

        // convert a number represented as val = sign * frac * 2^frac_exp, where frac is 64-bit
        // integer to floating point representation val = value + error, where value = fl(val),
        // error = val - value
        static void     result_as_float(bool frac_sign, uint64_t frac, int64_t frac_exp,
                            float& value, float& error);

    private:
        void            mult_scal(uint32_t x, const uint32_t* blocks, uint32_t* xC_prod);

        void            extract_result(int frac_bit_0, uint32_t* xC_prod, int& integral_part, 
                            bool& frac_sign, uint64_t& frac, int64_t& frac_exp);

        static double   pow2ki_sign(int64_t exp, int64_t sign);
};

// decompose a double precision scalar x as
//      x = (k + 4 * l) * pi/2 + frac
// where k, l are integers, 0 <= k < 4, |frac| <= pi/4, using the Payne-Hanek algorithm;
// this function returns k, the fractional part frac is represented as frac = value + error, 
// where value = fl(frac), error = frac_app - value, where frac_app is the approximation of
// frac correct up to at least 64 bits
MATCL_SIMD_EXPORT
int reduce_pi2_ph(double x, double& value, double& error);

// decompose a single precision scalar x as
//      x = (k + 4 * l) * pi/2 + frac
// where k, l are integers, 0 <= k < 4, |frac| <= pi/4, using the Payne-Hanek algorithm;
// this function returns k, the fractional part frac is represented as frac = value + error, 
// where value = fl(frac), error = frac_app - value, where frac_app is the approximation of
// frac correct up to at least 48 bits
MATCL_SIMD_EXPORT
int reduce_pi2_ph_float(float x, int max_bits, float& value, float& error);

// decompose a single precision scalar x represented by the double type as
//      x = (k + 4 * l) * pi/2 + frac
// where k, l are integers, 0 <= k < 4, |frac| <= pi/4, using the Payne-Hanek algorithm;
// this function returns k, the fractional part frac is represented as frac = value, 
// where value = fl(frac_app), where frac_app is the approximation of frac correct up to at
// least 53 bits; it is required, that float(x) == x
MATCL_SIMD_EXPORT
int reduce_pi2_ph_float(double x, int max_bits, double& value);

}};
