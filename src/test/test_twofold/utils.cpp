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

#include "utils.h"
#include "matcl-scalar/lib_functions/scalar_gen.h"

namespace matcl { namespace test
{

int rand_scalar_base::rand_pow(int min_exp, int max_exp)
{
    int num_pow     = max_exp - min_exp + 1;
    unsigned ir     = std::abs(irand());
    int pow         = min_exp + (ir % num_pow);
    return pow;
};

//--------------------------------------------------------------------------
//                         double
//--------------------------------------------------------------------------
double rand_scalar<double>::make(bool sec_od_sum)
{
    double val  = make_value(sec_od_sum);
    return val;
};

double rand_scalar<double>::rand()
{
    return matcl::rand();
};

double rand_scalar<double>::make_value(bool second_of_sum)
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
        case 8: return 1.5;
        case 9: return -1.5;
        case 10: return -std::numeric_limits<double>::infinity();
        case 11: return std::numeric_limits<double>::infinity();
        case 12: return -std::numeric_limits<double>::quiet_NaN();
        case 13: return std::numeric_limits<double>::quiet_NaN();
        default:
            return rand() * rand_scale(second_of_sum);
    };
};

double rand_scalar<double>::rand_scale(bool second_of_sum)
{
    int max_exp;
    int min_exp;

    if (second_of_sum == false)
    {
        max_exp     = std::numeric_limits<double>::max_exponent/3;
        min_exp     = std::numeric_limits<double>::min_exponent/3;
    }
    else
    {
        // should be in range [eps^2, 1], but we take larger range
        max_exp     = 5;
        min_exp     = -110;
    }

    int pow         = rand_pow(min_exp, max_exp);

    double sign     = 1.0;
    if (rand() < 0.5)
        sign        = -1.0;

    double val       = sign * std::ldexp(1.0, pow);
    return val;
};

//--------------------------------------------------------------------------
//                         float
//--------------------------------------------------------------------------
float rand_scalar<float>::make(bool sec_od_sum)
{
    float val  = make_value(sec_od_sum);
    return val;
};

float rand_scalar<float>::rand()
{
    return matcl::frand();
};

float rand_scalar<float>::make_value(bool second_of_sum)
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
        case 8: return 1.5f;
        case 9: return -1.5f;
        case 10: return -std::numeric_limits<float>::infinity();
        case 11: return std::numeric_limits<float>::infinity();
        case 12: return -std::numeric_limits<float>::quiet_NaN();
        case 13: return std::numeric_limits<float>::quiet_NaN();
        default:
            return frand() * rand_scale(second_of_sum);
    };
};

float rand_scalar<float>::rand_scale(bool second_of_sum)
{
    int max_exp;
    int min_exp;

    if (second_of_sum == false)
    {
        max_exp     = std::numeric_limits<float>::max_exponent/3;
        min_exp     = std::numeric_limits<float>::min_exponent/3;
    }
    else
    {
        // should be in range [eps^2, 1], but we take larger range
        max_exp     = 5;
        min_exp     = -2 * std::numeric_limits<float>::digits + 5;
    }

    int pow         = rand_pow(min_exp, max_exp);

    float sign      = 1.0f;
    if (rand() < 0.5)
        sign        = -1.0f;

    float val       = sign * std::ldexp(1.0f, pow);
    return val;
};

}};
