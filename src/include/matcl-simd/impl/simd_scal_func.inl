/*
 *  This file is a part of Morfa Matrix Lib.
 *
 *  Copyright (c) Pawe³ Kowal 2015-2016
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

#include "matcl-simd/arch/reg_128/simd_128.h"
#include "matcl-simd/simd/simd_scalar_func.h"
#include "matcl-core/details/simd/simd.h"

namespace matcl { namespace simd
{

//TODO:
/*
force_inline
double fma(const double& x, const double& y, const double& z)
{
    using simd_type = simd<double,reg_128>;
    simd_type ret   = fma(simd_type(x),simd_type(y),simd_type(z));
    return ret.get(0);

    //return x*y+z;
};

force_inline
float fma(const float& x, const float& y, const float& z)
{
    using simd_type = simd<float,reg_128>;
    simd_type ret   = fma(simd_type(x),simd_type(y),simd_type(z));
    return ret.get(0);

    //return x*y+z;
};

force_inline
double fms(const double& x, const double& y, const double& z)
{
    using simd_type = simd<double,reg_128>;
    simd_type ret   = fms(simd_type(x),simd_type(y),simd_type(z));
    return ret.get(0);

    //return x*y+z;
};
force_inline
float fms(const float& x, const float& y, const float& z)
{
    using simd_type = simd<float,reg_128>;
    simd_type ret   = fms(simd_type(x),simd_type(y),simd_type(z));
    return ret.get(0);

    //return x*y+z;
};
*/

}}
