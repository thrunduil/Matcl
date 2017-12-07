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

#include "matcl-simd/simd_general.h"

namespace matcl { namespace simd { namespace details
{

// Generate a constant vector of 8 integers stored in memory,
// load as __m256
template <int I0, int I1, int I2, int I3, int I4, int I5, int I6, int I7>
static inline __m256i vector_8_int() 
{
    using m256_int8     = union 
                        { 
                            int     i[8];
                            __m256i ymm;
                        };
    
    static const m256_int8 val = {{I0,I1,I2,I3,I4,I5,I6,I7}};
    return val.ymm;
};

}}}
