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

#include "matcl-scalar/config.h"

#include "sfmt/SFMT.h"
#include "dsfmt/dSFMT.h"
#include <memory>

namespace matcl { namespace details
{

const matcl::rand_state&    global_rand_state();
matcl::rand_state           get_rand_state();
void                        set_rand_state(const matcl::rand_state& st);

class rand_state
{
    private:
        sfmt_t      m_state_i;
        dsfmt_t     m_state_d;

        double      m_b;
        int         m_flag;

    public:
        //initialize state with seed = 0
        rand_state();

        //initialize state with given seed
        rand_state(Integer seed);

        // This function initializes the internal state array with a seed.
        void init_genrand(Integer seed);

        // This function generates and returns 32-bit pseudorandom number.
        Integer     genrand_int32();

        // generates a random number on [0,1)-real-interval
        double      genrand_real2();

        // generates a random number on (0,1)-real-interval
        double      genrand_real3();

        // Standard normally distributed random variables.
        double      gen_norm();
};

};};
