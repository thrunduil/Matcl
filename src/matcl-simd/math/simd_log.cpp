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

#include "matcl-simd/details/math/impl/simd_log.h"

namespace matcl { namespace simd { namespace details
{

const double simd_log_approx_double_data::poly_scalar[]  = 
{
        3.99999999987401741406920347834347647e-01L, 6.66666666666688547299310656675190687e-01L,
        2.22221871825960117598424175414698074e-01L, 2.85714288669539260936377812250363457e-01L,
        1.53009602706038610316834383087625899e-01L, 1.81841095988965909080780869966928340e-01L,
        0.0                                       , 1.49202865515212540004373762442151109e-01L,
};

const double simd_log_approx_double_data::poly[]   =  
{
    6.66666666666688547299310656675190687e-01L,
    3.99999999987401741406920347834347647e-01L,
    2.85714288669539260936377812250363457e-01L,
    2.22221871825960117598424175414698074e-01L,
    1.81841095988965909080780869966928340e-01L,
    1.53009602706038610316834383087625899e-01L,
    1.49202865515212540004373762442151109e-01L,
};

const float simd_log_approx_float_data::poly[] = 
{
    -5.00000001969598183278064207495559219e-01f,
    3.33332716385429135870565468876902043e-01f,
    -2.49998863941982779960396002029673983e-01f,
    2.00073054372748504562536157172019915e-01f,
    -1.66759465098372286204627390624499906e-01f,
    1.40545176272346789443221873134537537e-01f,
    -1.22505576737375893726817365117891765e-01f,
    1.37550000113041486713019781724756740e-01f,
    -1.25948724727867811279422965776243012e-01f,
};

const float simd_log_approx_float_data::poly_scalar[] = 
{
    -5.00000001969598183278064207495559219e-01f,
    3.33332716385429135870565468876902043e-01f,
    -2.49998863941982779960396002029673983e-01f,
    2.00073054372748504562536157172019915e-01f,

    -1.66759465098372286204627390624499906e-01f,
    1.40545176272346789443221873134537537e-01f,
    -1.22505576737375893726817365117891765e-01f,
    1.37550000113041486713019781724756740e-01f,
            
    -1.25948724727867811279422965776243012e-01f,
    0.0f,
    0.0f,
    0.0f,
};

}}}
