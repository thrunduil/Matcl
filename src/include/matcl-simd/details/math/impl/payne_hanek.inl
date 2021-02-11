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

#include "matcl-core/config.h"
#include "payne_hanek.h"

#include "matcl-simd/simd_math.h"
#include <algorithm>

#pragma warning(push)
#pragma warning(disable : 4127) // conditional expression is constant

namespace matcl { namespace simd
{

namespace ms = matcl::simd;

//-----------------------------------------------------------------------
//                      DOUBLE
//-----------------------------------------------------------------------

template<int Num_blocks>
force_inline
void payne_hanek_double<Num_blocks>::decompose(const double& a, bool& signed_x, uint64_t& ai, int& exp_a)
{
    // represent x as +- sig * 2^e, where 1 <= sig < 2
    exp_a               = (int)ms::iexponent(a);
    exp_a               = exp_a - 1;

    ai                  = reinterpret_cast<const uint64_t&>(a);

    // set exponent to 0, clear sign flag and set bit 52
    ai                  = (ai & 0x000FFFFFFFFFFFFFll) | 0x10000000000000ll;
    signed_x            = (a < 0);
};

template<int Num_blocks>
force_inline
double payne_hanek_double<Num_blocks>::pow2ki_sign(int64_t k, int64_t sign)
{
    int64_t ik          = k + int64_t(1023);
    ik                  = ik << 52;
    ik                  = ik | sign;

    return reinterpret_cast<double&>(ik);
}

template<int Num_blocks>
force_inline
void payne_hanek_double<Num_blocks>::mult_scal(uint64_t xlo, uint64_t xhi, const uint32_t* blocks, uint32_t* xC_prod)
{
    // carry stores at most 33 bits
    uint64_t carry;    

    {
        uint64_t prod       = xlo * blocks[Num_blocks - 1];
        xC_prod[Num_blocks+1] = (uint32_t)prod;
        carry               = prod >> 32;
    }

    for (int k = Num_blocks - 2; k >= 0; --k)
    {
        //[32]*[32] + [21]*[32] + [33] <= [65 bits]

        uint64_t prod1      = xlo * blocks[k];                  // 64 bits
        uint64_t prod2      = xhi * blocks[k+1];                // 53 bits
        uint32_t prod1_lo   = (uint32_t)prod1;                  // 32 bits
        uint64_t prod1_hi   = prod1 >> number_bits;             // 32 bits         
        uint64_t prod       = prod1_lo + prod2 + carry;         // 54 bits
        uint64_t prod_hi    = (prod >> number_bits) + prod1_hi; // 33 bits

        xC_prod[k+2]        = (uint32_t)prod;

        // at most 33 bits
        carry               = prod_hi;
    }
    {
        uint64_t prod       = xhi * blocks[0] + carry;          // 54 bits
        uint64_t prod_hi    = prod >> number_bits;              // 22 bits
        xC_prod[1]          = (uint32_t)prod;
        xC_prod[0]          = (uint32_t)prod_hi;
    }
};

template<int Num_blocks>
force_inline
void payne_hanek_double<Num_blocks>::extract_result(int frac_bit_0, uint32_t* xC_prod, 
                        int& integral_part, bool& frac_sign, uint64_t& frac, int64_t& frac_exp)
{
    static const uint32_t one_shift_32_m1   = (uint32_t)(- 1); 

    //-------------------------------------------------------------------
    //                  COMPUTE CONSTANTS
    //-------------------------------------------------------------------    
    int frac_bit            = std::max(frac_bit_0, 0);
    int byte_pos_frac       = frac_bit / number_bits;
    int rem_bits            = frac_bit % number_bits;

    //-------------------------------------------------------------------
    //                  EXTRACT INTAGRAL PART
    //-------------------------------------------------------------------    

    {
        uint32_t int_part_lo    = (rem_bits == 0)? 0 : (xC_prod[byte_pos_frac] >> (number_bits - rem_bits));
        uint32_t int_part_hi    = xC_prod[byte_pos_frac-1] << rem_bits;
        uint32_t int_part       = int_part_lo | int_part_hi;
        integral_part           = int_part;
    };    
    
    uint32_t leading_bit_mask   = uint32_t(1) << (31 - rem_bits);

    // check if fractional part >= 0.5, then negate relevant bits, set sign to -1 
    // and increase integral_parta; void costly if statement
    uint32_t need_negate        = (xC_prod[byte_pos_frac] & leading_bit_mask);
    integral_part               += (need_negate == 0) ? 0 : 1;    

    uint32_t neg_mask           = (need_negate == 0) ? 0 : uint32_t(-1);

    for (int i = byte_pos_frac; i < Num_blocks + 2; ++i)
        xC_prod[i]              = xC_prod[i] ^ neg_mask;
    
    // clear bits of the integral part    
    uint32_t frac_part_mask     = uint32_t(-1) >> rem_bits;
    xC_prod[byte_pos_frac]      = xC_prod[byte_pos_frac] & frac_part_mask;    

    //-------------------------------------------------------------------
    //                  EXTRACT FRACTIONAL PART
    //-------------------------------------------------------------------

    // find first nonzero block
    m_first_nz                  = byte_pos_frac;

    for (; m_first_nz < Num_blocks + 2; ++m_first_nz)
    {
        if (xC_prod[m_first_nz] != 0)
            break;
    };

    // we need at most 5 32-bit blocks 
    // => this gives at least 4*32 + 1 bits = 53 + 76
    // block with this number exists (xC_prod has additional zero blocks at the end)
    int last_nz                 = m_first_nz + max_full_blocks + 1;

    //number of bits in the most significant block
    n_bits                      = number_bits - number_leading_zeros(xC_prod[m_first_nz]);

    uint64_t part_1             = (uint64_t(xC_prod[m_first_nz + 1]) << number_bits) 
                                | uint64_t(xC_prod[m_first_nz + 2]);
    uint64_t part_0             = xC_prod[m_first_nz];
    frac                        = (part_0 << (64 - n_bits)) | (part_1 >> n_bits);
    frac_sign                   = need_negate != 0;

    // compute exponent of the least significant bit
    //int p                     = Num_blocks * number_bits - (2 * 53 + frac_bit_0 - 64);
    //int64_t frac_exp          = (Num_blocks + 2 - last_nz) * number_bits - (p + 2 * 53);

    int64_t frac_exponent       = frac_bit_0 - last_nz * number_bits;
    frac_exp                    = frac_exponent + n_bits;
};

template<int Num_blocks>
force_inline
void payne_hanek_double<Num_blocks>::result_as_double(bool frac_sign, uint64_t frac, int64_t frac_exp,
                    double& frac_hi, double& frac_lo)
{
    int64_t sign                = frac_sign ? int64_t(1ull << 63) : 0;    
    uint32_t part_hi            = frac >> 32;
    uint32_t part_lo            = (uint32_t)frac;

    double pow_lo               = pow2ki_sign(frac_exp, sign);
    double pow_hi               = pow_lo * 4294967296.0; // 2^32

    double res_lo               = double(part_lo) * pow_lo;
    double res_hi               = double(part_hi) * pow_hi;

    frac_hi                     = res_hi + res_lo;
    double tmp                  = frac_hi - res_hi;
    frac_lo                     = res_lo - tmp;
};

template<int Num_blocks>
force_inline
int payne_hanek_double<Num_blocks>::number_correct_bits() const
{
    int max_correct_digits      = (Num_blocks + 1 - m_first_nz) * number_bits + n_bits - 53;
    return max_correct_digits;
};

template<int Num_blocks>
force_inline
void payne_hanek_double<Num_blocks>::mult(uint64_t significant_1, int64_t exp_1, uint64_t significant_2, 
                    int64_t exp_2, uint64_t& significant_res, int64_t& exp_res)
{
    significant_res     = ms::mulh(significant_1, significant_2);
    exp_res             = exp_1 + exp_2 + 64;
};

template<int Num_blocks>
template<int Req_m>
force_inline
int payne_hanek_double<Num_blocks>::eval(bool is_signed, uint64_t ai, int exp_a, 
                        const uint32_t* constant, int const_exponent, bool& frac_signed, 
                        uint64_t& frac, int64_t& frac_exponent)
{
    static const uint32_t one_shift_32_m1   = (uint32_t)(- 1);

    //-------------------------------------------------------------------
    //                  PREPARE INPUTS
    //-------------------------------------------------------------------
    // exponent of a * C
    int exp             = exp_a + const_exponent;

    //-------------------------------------------------------------------
    //                  COMPUTE CONSTANTS
    //-------------------------------------------------------------------

    // the constant C is represented as C= {1. a[-1], a[-2], ... } and split as
    // C = Left x 2^(53-e+m) + Med x 2^(-53-e-1-p) + Right * 2^(-53-e-1-p)
    //      Left = {1, a[-1], ..., a[53-e+m]} or Left = 0 if 53-e+m >= 0
    //      Med  = {a[53-e+m-1], a[53-e], ..., a[-53-e-1-p] }
    //      Right= {0. a[-53-e-2-p], ...}
    // let x = X * 2^(e-52);   2^52 <= X < 2^53
    // then x * C = X*Left x 2^(1+m) + X*Med x 2^(-2*53-p) + X*Right * 2^(-2*53-p)
    //            = 2^(1+m) * r + X*Med x 2^(-2*53-p) + rem, rem < 2^(-53-p) 

    // Med has 2*53+p+m+1 - n_missing_digits digits
    //int req_bits      = (2 * 53 + m + 1 - n_missing_digits); //+p    
    // implied precision:
    //int p             = number_blocks * number_bits - req_bits;

    // get position of relevant bits of the constant Med
    int32_t start_pos_0 = exp - Req_m + 2 - 53;
    int32_t start_pos   = std::max(start_pos_0, 0); 

    uint32_t block_num  = (uint32_t)start_pos >> 5;
    uint32_t e0         = (uint32_t)start_pos & (number_bits - 1);

    // do not shift bits; just increase m
    //int m             = Req_m + e0;

    // number of digits implicitely set to zero (when 53-e+m-1 - j > 0)
    // Med has 2*53+p+m+1 - n_missing_digits digits
    //int n_missing_digits= std::max((int32_t)e0 - start_pos_0, 0);

    // position of the first bit in the fractional part (possibly zero)
    int frac_bit_0      = 64 + Req_m + std::min(start_pos_0,  (int32_t)e0);

    //-------------------------------------------------------------------
    //                  FORM PRODUCT A * C
    //-------------------------------------------------------------------

    // blocks store the Med part of the constant C    
    uint32_t blocks[Num_blocks];

    // copy required bits from the constant
    std::memcpy(blocks, constant + block_num, Num_blocks * sizeof(uint32_t));

    // decompose significant into 32-bit blocks: ai = xlo + xhi * 2^32
    uint64_t xlo        = ai & one_shift_32_m1;     // 32 bits
    uint64_t xhi        = ai >> number_bits;        // 21 bits

    // form product (xhi * 2^32 + xlo) * Med = xhi * Med * 2^32 + xlo * Med

    uint32_t xC_prod_arr[(Num_blocks + 2) + 1 + (max_full_blocks + 1)] = {0};
    uint32_t* xC_prod   = xC_prod_arr + 1;

    mult_scal(xlo, xhi, blocks, xC_prod);

    //-------------------------------------------------------------------
    //        DECOMPOSE PRODUCT INTO INTEGRAL AND FRACTIONAL PART
    //-------------------------------------------------------------------

    int integral_part;
    extract_result(frac_bit_0, xC_prod, integral_part, frac_signed, frac, frac_exponent);

    static const 
    int m_mod           = (1 << Req_m);

    integral_part       = integral_part & (m_mod - 1);

    if (is_signed == true)
    {
        frac_signed     = !frac_signed;
        integral_part   = (integral_part == 0) ? integral_part : m_mod - integral_part;
    }

    return integral_part;
}

//-----------------------------------------------------------------------
//                      FLOAT
//-----------------------------------------------------------------------

template<int Num_blocks>
template<int Req_m>
force_inline
int payne_hanek_float<Num_blocks>::eval(bool is_signed, uint32_t ai, int exp_a, const uint32_t* constant, 
                    int const_exponent, bool& frac_signed, uint64_t& frac, int64_t& frac_exponent)
{
    static const uint32_t one_shift_32_m1   = (uint32_t)(- 1);

    // exponent of a * C
    int exp             = exp_a + const_exponent;

    //-------------------------------------------------------------------
    //                  COMPUTE CONSTANTS
    //-------------------------------------------------------------------

    // the constant C is represented as C= {1. a[-1], a[-2], ... } as split as
    // C = Left x 2^(24-e+m) + Med x 2^(-24-e-1-p) + Right * 2^(-24-e-1-p)
    //      Left = {1, a[-1], ..., a[24-e+m]} or Left = 0 if 24-e+m >= 0
    //      Med  = {a[24-e+m-1], a[24-e], ..., a[-24-e-1-p] }
    //      Right= {0. a[-24-e-2-p], ...}
    // let x = X * 2^(e-23);   2^23 <= X < 2^24
    // then x * C = X*Left x 2^(1+m) + X*Med x 2^(-2*24-p) + X*Right * 2^(-2*24-p)
    //            = 2^(1+m) * r + X*Med x 2^(-2*24-p) + rem, rem < 2^(-24-p) 

    // Med has 2*24+p+m+1 - n_missing_digits digits
    //int req_bits      = (2 * 24 + m + 1 - n_missing_digits); //+p    
    // implied precision:
    //int p             = number_blocks * number_bits - req_bits;

    // get position of relevant bits of the constant Med
    int32_t start_pos_0 = exp - Req_m + 2 - 24;
    int32_t start_pos   = std::max(start_pos_0, 0); 

    uint32_t block_num  = (uint32_t)start_pos >> 5;
    uint32_t e0         = (uint32_t)start_pos & (number_bits - 1);

    // do not shift bits; just increase m
    //int m             = Req_m + e0;

    // number of digits implicitely set to zero (when 24-e+m-1 - j > 0)
    // Med has 2*23+p+m+1 - n_missing_digits digits
    //int n_missing_digits= std::max((int32_t)e0 - start_pos_0, 0);

    // position of the first bit in the fractional part (possibly zero)
    int frac_bit_0      = 32 + Req_m + std::min(start_pos_0,  (int32_t)e0);

    //-------------------------------------------------------------------
    //                  FORM PRODUCT A * C
    //-------------------------------------------------------------------

    // blocks store the Med part of the constant C    
    uint32_t blocks[Num_blocks];

    // copy required bits from the constant
    std::memcpy(blocks, constant + block_num, Num_blocks * sizeof(uint32_t));

    // form product ai * Med

    uint32_t xC_prod_arr[(Num_blocks + 1) + 1 + (max_full_blocks + 1)] = {0};
    uint32_t* xC_prod = xC_prod_arr + 1;

    mult_scal(ai, blocks, xC_prod);

    //-------------------------------------------------------------------
    //        DECOMPOSE PRODUCT INTO INTEGRAL AND FRACTIONAL PART
    //-------------------------------------------------------------------

    int integral_part;
    extract_result(frac_bit_0, xC_prod, integral_part, frac_signed, frac, frac_exponent);

    static const 
    int m_mod           = (1 << Req_m);

    integral_part       = integral_part & (m_mod - 1);

    if (is_signed == true)
    {
        frac_signed     = !frac_signed;
        integral_part   = (integral_part == 0) ? integral_part : m_mod - integral_part;
    }

    return integral_part;
}

template<int Num_blocks>
force_inline
void payne_hanek_float<Num_blocks>::mult_scal(uint32_t x0, const uint32_t* blocks, uint32_t* xC_prod)
{
    // carry stores at most 25 bits
    uint64_t carry          = 0;    
    uint64_t x              = uint64_t(x0);

    for (int k = Num_blocks - 1; k >= 0; --k)
    {
        //[24]*[32] + [25] <= [57 bits]

        uint64_t prod       = x * blocks[k] + carry;    // 57 bits
        xC_prod[k+1]        = (uint32_t)prod;
        carry               = (prod >> number_bits);    // 25 bits
    }

    xC_prod[0]              = (uint32_t)carry;
};

template<int Num_blocks>
force_inline
void payne_hanek_float<Num_blocks>::extract_result(int frac_bit_0, uint32_t* xC_prod, 
                        int& integral_part, bool& frac_sign, uint64_t& frac, int64_t& frac_exp)
{
    static const uint32_t one_shift_32_m1   = (uint32_t)(- 1); 

    //-------------------------------------------------------------------
    //                  COMPUTE CONSTANTS
    //-------------------------------------------------------------------    
    int frac_bit            = std::max(frac_bit_0, 0);
    int byte_pos_frac       = frac_bit / number_bits;
    int rem_bits            = frac_bit % number_bits;

    //-------------------------------------------------------------------
    //                  EXTRACT INTAGRAL PART
    //-------------------------------------------------------------------    

    {
        uint32_t int_part_lo     = (rem_bits == 0)? 0 : (xC_prod[byte_pos_frac] >> (number_bits - rem_bits));
        uint32_t int_part_hi    = xC_prod[byte_pos_frac-1] << rem_bits;
        uint32_t int_part       = int_part_lo | int_part_hi;
        integral_part           = int_part;
    };
    
    uint32_t leading_bit_mask   = uint32_t(1) << (31 - rem_bits);

    // check if fractional part >= 0.5, then negate relevant bits, set sign to -1 
    // and increase integral_parta; void costly if statement
    uint32_t need_negate        = (xC_prod[byte_pos_frac] & leading_bit_mask);
    integral_part               += (need_negate == 0) ? 0 : 1;    

    uint32_t neg_mask           = (need_negate == 0) ? 0 : uint32_t(-1);

    for (int i = byte_pos_frac; i < Num_blocks + 1; ++i)
        xC_prod[i]              = xC_prod[i] ^ neg_mask;
    
    // clear bits of the integral part    
    uint32_t frac_part_mask     = uint32_t(-1) >> rem_bits;
    xC_prod[byte_pos_frac]      = xC_prod[byte_pos_frac] & frac_part_mask;    

    //-------------------------------------------------------------------
    //                  EXTRACT FRACTIONAL PART
    //-------------------------------------------------------------------

    // find first nonzero block
    m_first_nz                  = byte_pos_frac;

    for (; m_first_nz < Num_blocks + 1; ++m_first_nz)
    {
        if (xC_prod[m_first_nz] != 0)
            break;
    };

    // block with this number exists (xC_prod has additional zero blocks at the end)
    int last_nz                 = m_first_nz + max_full_blocks + 1;

    //number of bits in the most significant block
    n_bits                      = number_bits - number_leading_zeros(xC_prod[m_first_nz]);

    uint64_t part_1             = (uint64_t(xC_prod[m_first_nz + 1]) << number_bits) 
                                | uint64_t(xC_prod[m_first_nz + 2]);
    uint64_t part_0             = xC_prod[m_first_nz];
    frac                        = (part_0 << (64 - n_bits)) | (part_1 >> n_bits);
    frac_sign                   = need_negate != 0;

    // compute exponent of the least significant bit
    //int p                     = Num_blocks * number_bits - (2 * 24 + frac_bit_0);
    //int64_t frac_exp          = (Num_blocks + 1 - last_nz) * number_bits - (p + 2 * 24);
    int64_t frac_exponent       = frac_bit_0 - last_nz * number_bits;
    frac_exp                    = frac_exponent + n_bits;    
};

template<int Num_blocks>
force_inline
int payne_hanek_float<Num_blocks>::number_correct_bits() const
{
    int max_correct_digits      = (Num_blocks - m_first_nz) * number_bits + n_bits - 24;
    return max_correct_digits;
};

template<int Num_blocks>
force_inline
void payne_hanek_float<Num_blocks>::mult(uint64_t significant_1, int64_t exp_1, uint64_t significant_2, 
                    int64_t exp_2, uint64_t& significant_res, int64_t& exp_res)
{
    significant_res     = ms::mulh(significant_1, significant_2);
    exp_res             = exp_1 + exp_2 + 64;
};

template<int Num_blocks>
force_inline
void payne_hanek_float<Num_blocks>::result_as_double(bool frac_sign, uint64_t frac, int64_t frac_exp,
                    double& res)
{
    int64_t sign                = frac_sign ? int64_t(1ull << 63) : 0;    
    double pow                  = pow2ki_sign(frac_exp, sign);
    res                         = double(frac) * pow;
};

template<int Num_blocks>
force_inline
void payne_hanek_float<Num_blocks>::result_as_float(bool frac_sign, uint64_t frac, int64_t frac_exp,
                    float& frac_hi, float& frac_lo)
{
    double res;
    result_as_double(frac_sign, frac, frac_exp, res);

    frac_hi     = float(res);
    frac_lo     = float(res - double(frac_hi));
};

template<int Num_blocks>
force_inline
double payne_hanek_float<Num_blocks>::pow2ki_sign(int64_t k, int64_t sign)
{
    int64_t ik          = k + int64_t(1023);
    ik                  = ik << 52;
    ik                  = ik | sign;

    return reinterpret_cast<double&>(ik);
}

template<int Num_blocks>
force_inline
void payne_hanek_float<Num_blocks>::decompose(const double& x, bool& signed_x, uint32_t& sign_x, int& exp_x)
{
    // represent x as +- sig * 2^e, where 1 <= sig < 2
    exp_x               = (int)ms::iexponent(x);
    exp_x               = exp_x - 1;

    uint64_t sign       = reinterpret_cast<const uint64_t&>(x);

    // set exponent to 0, clear sign flag and set bit 52
    sign                = (sign & 0x000FFFFFFFFFFFFFll) | 0x10000000000000ll;

    // remove unused bits
    sign                = sign >> (53 - 24);
    sign_x              = (uint32_t)sign;
    signed_x            = (x < 0.0);
};

template<int Num_blocks>
force_inline
void payne_hanek_float<Num_blocks>::decompose(const float& x, bool& signed_x, uint32_t& sign_x, int& exp_x)
{
    // represent x as +- sig * 2^e, where 1 <= sig < 2
    exp_x               = (int)ms::iexponent(x);
    exp_x               = exp_x - 1;

    sign_x              = reinterpret_cast<const uint32_t&>(x);

    // set exponent to 0, clear sign flag and set bit 24
    sign_x              = (sign_x & 0x00FFFFFFl) | 0x00800000l;
    signed_x            = (x < 0.f);
}

}}

#pragma warning(pop)