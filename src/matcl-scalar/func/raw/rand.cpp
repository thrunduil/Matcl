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

#include "matcl-scalar/func/raw/rand.h"
#include "matcl-scalar/extern/sfmt_wrap.h"

#pragma warning (disable:4127)

namespace matcl { namespace raw
{

namespace mr = matcl :: raw;

Real mr::rand(const matcl::rand_state& rand_ptr)
{
    return rand_ptr->genrand_real3();
}

Float mr::frand(const matcl::rand_state& rand_ptr)
{
    return static_cast<Float>(rand_ptr->genrand_real3());
}

Real mr::randn(const matcl::rand_state& rand_ptr)
{
    return rand_ptr->gen_norm();
}

Float mr::frandn(const matcl::rand_state& rand_ptr)
{
    return static_cast<Float>(rand_ptr->gen_norm());
}

};};