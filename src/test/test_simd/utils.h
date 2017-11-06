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

#include "test_simd_config.h"
#include "matcl-core/matrix/scalar_types.h"

namespace matcl { namespace test
{

template<class T>
struct rand_scalar;

template<>
struct rand_scalar<float>
{
    static float    make(bool testing_values);
    static float    make_special();
    static float    rand_scale();
};

template<>
struct rand_scalar<double>
{
    static double   make(bool testing_values);
    static double   make_special();
    static double   rand_scale();
};

template<>
struct rand_scalar<Float_complex>
{
    static Float_complex    make(bool testing_values);
    static Float_complex    make_special();
};

template<>
struct rand_scalar<Complex>
{
    static Complex  make(bool testing_values);
    static Complex  make_special();
};

}};

