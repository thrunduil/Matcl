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

#include "test_simd_config.h"
#include "utils.h"
#include "matcl-scalar/lib_functions/scalar_gen.h"

namespace matcl { namespace test
{

float rand_scalar<float>::make(bool testing_values)
{
    if (testing_values)
        return make_special();
    else
        return 10.0f * frandn();
};

float rand_scalar<float>::make_special()
{
    int max_code    = 27;
    int code        = std::abs(irand()) % max_code;

    switch(code)
    {
        case 0: return 0.0f;
        case 1: return -0.0f;
        case 2: return 1.0f;
        case 3: return -1.0f;
        case 4: return 0.5f;
        case 5: return -0.5f;
        case 6: return 2.0f;
        case 7: return -2.0f;
        case 8: return -std::numeric_limits<float>::infinity();
        case 9: return std::numeric_limits<float>::infinity();
        case 10: return -std::numeric_limits<float>::quiet_NaN();
        case 11: return std::numeric_limits<float>::quiet_NaN();
        case 12: return -std::numeric_limits<float>::max();
        case 13: return std::numeric_limits<float>::max();
        case 14: return -std::numeric_limits<float>::min();
        case 15: return std::numeric_limits<float>::min();
        case 16: return 1.5f;
        case 17: return -1.5f;
        default:
            return frand() * rand_scale();
    };
};

float rand_scalar<float>::rand_scale()
{
    int max_exp     = std::numeric_limits<float>::max_exponent;
    int min_exp     = std::numeric_limits<float>::min_exponent;

    int num_pow     = max_exp - min_exp + 1;
    unsigned ir     = abs(irand());
    int pow         = min_exp + (ir % num_pow);

    float sign      = 1.0f;
    if (rand() < 0.5)
        sign        = -1.0f;

    float val       = sign * std::ldexp(1.0f, pow);
    return val;
};

int32_t rand_scalar<int32_t>::make(bool testing_values)
{
    (void)testing_values;
    return irand();
};

double rand_scalar<double>::make(bool testing_values)
{
    if (testing_values)
        return make_special();
    else
        return 10.0 * randn();
};

double rand_scalar<double>::make_special()
{
    int max_code    = 27;
    int code        = std::abs(irand()) % max_code;

    switch(code)
    {
        case 0: return 0.0;
        case 1: return -0.0;
        case 2: return 1.0;
        case 3: return -1.0;
        case 4: return 0.5;
        case 5: return -0.5;
        case 6: return 2.0;
        case 7: return -2.0;
        case 8: return -std::numeric_limits<double>::infinity();
        case 9: return std::numeric_limits<double>::infinity();
        case 10: return -std::numeric_limits<double>::quiet_NaN();
        case 11: return std::numeric_limits<double>::quiet_NaN();
        case 12: return -std::numeric_limits<double>::max();
        case 13: return std::numeric_limits<double>::max();
        case 14: return -std::numeric_limits<double>::min();
        case 15: return std::numeric_limits<double>::min();
        case 16: return 1.5;
        case 17: return -1.5;

        default:
            return rand() * rand_scale();
    };
};

double rand_scalar<double>::rand_scale()
{
    int max_exp = std::numeric_limits<double>::max_exponent;
    int min_exp = std::numeric_limits<double>::min_exponent;

    const int num_pow    = max_exp - min_exp + 1;

    // rand power
    unsigned ir     = abs(irand());
    int pow         = min_exp + (ir % num_pow);

    double sign     = 1.0;
    if (rand() < 0.5)
        sign        = -1.0;

    double val      = sign * std::ldexp(1.0, pow);
    return val;
};

int64_t rand_scalar<int64_t>::make(bool testing_values)
{
    (void)testing_values;
    return int64_t(irand()) * int64_t(irand());    
};

Complex rand_scalar<Complex>::make(bool testing_values)
{
    if (testing_values)
        return make_special();
    else
        return 10.0 * crandn();
};

Complex rand_scalar<Complex>::make_special()
{
    return Complex(rand_scalar<double>::make_special(), rand_scalar<double>::make_special());
};

Float_complex rand_scalar<Float_complex>::make(bool testing_values)
{
    if (testing_values)
        return make_special();
    else
        return 10.0f * fcrandn();
};

Float_complex rand_scalar<Float_complex>::make_special()
{
    return Float_complex(rand_scalar<float>::make_special(), rand_scalar<float>::make_special());
};

}};
