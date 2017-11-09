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

#include "matcl-core/matrix/scalar_types.h"
#include "matcl-core/float/twofold.h"

namespace matcl { namespace test
{

template<class T>
struct rand_scalar;

struct rand_scalar_base
{
    static int      rand_pow(int min_exp, int max_exp);
};

template<>
struct rand_scalar<twofold> : rand_scalar_base
{
    static twofold  make(bool normalized, bool second_of_sum);    
    static double   rand_mult(bool normalized);    
};

template<>
struct rand_scalar<double> : rand_scalar_base
{
    static double   make(bool second_of_sum);
    static double   make_value(bool second_of_sum);
    static double   rand_scale(bool second_of_sum);
};

template<>
struct rand_scalar<float> : rand_scalar_base
{
    static float    make(bool second_of_sum);
    static float    make_value(bool second_of_sum);
    static float    rand_scale(bool second_of_sum);
};

}};

