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

#include "matcl-core/lib_functions/constants.h"
#include "matcl-core/memory/global_objects.h"
#include "matcl-core/matrix/matrix_traits.h"
#include "boost/math/constants/constants.hpp"

namespace matcl 
{

Real constants::eps(value_code vc)
{
    if (matrix_traits::is_single_precision(vc))
        return eps<Float>();
    else
        return eps<Real>();
}
Real constants::min_real(value_code vc)
{
    if (matrix_traits::is_single_precision(vc))
        return min_real<Float>();
    else
        return min_real<Real>();
}
Real constants::max_real(value_code vc)
{
    if (matrix_traits::is_single_precision(vc))
        return max_real<Float>();
    else
        return max_real<Real>();
}

};
