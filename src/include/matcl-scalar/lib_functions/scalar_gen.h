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

#include "matcl-scalar/config.h"
#include "matcl-scalar/object.h"
#include "matcl-scalar/details/scalar_gen_details.h"

#include <memory>

#pragma warning(push)
#pragma warning(disable: 4251)  //needs to have dll-interface to be used by clients of class

namespace matcl
{

//-----------------------------------------------------------------
//                  random numbers generator state
//-----------------------------------------------------------------
// type of random number generator state
// random number generator is thread local
class MATCL_SCALAR_EXPORT rand_state
{
    //internal use
    public:
        using rand_state_ptr = std::shared_ptr<details::rand_state>;

        rand_state_ptr  m_ptr;

        rand_state(const rand_state_ptr&);

        const rand_state_ptr&   operator->() const  { return m_ptr; };

    public:
        rand_state();
        ~rand_state();
};

// return independent copy of current state of random number generator
MATCL_SCALAR_EXPORT 
rand_state          get_rand_state();

// set current state of random number generator, random number generator is
// thread local, therefore this function affects only current thread
MATCL_SCALAR_EXPORT
void                set_rand_state(rand_state st);

// get state of random number generator in current thread
MATCL_SCALAR_EXPORT 
const rand_state&   global_rand_state();

// initialize state of random number generator pointed by rand_ptr with seed s;
// this affects only current thread; 
MATCL_SCALAR_EXPORT 
void                init_genrand(unsigned long s, const rand_state& rand_ptr 
                                    = global_rand_state());

// create local state of random number generator, restore previous state
// after exiting from current scope
MATCL_SCALAR_EXPORT 
details::local_rand_state_raii
                    local_rand_state(unsigned long s);  

//-----------------------------------------------------------------
//                      random numbers
//-----------------------------------------------------------------
// generate uniformly distributed Real random number on (0, 1)
MATCL_SCALAR_EXPORT 
Real                rand(const rand_state& rand_ptr = global_rand_state());

// generate uniformly distributed Float random number on (0, 1)
MATCL_SCALAR_EXPORT 
Float               frand(const rand_state& rand_ptr = global_rand_state());

// generate uniformly distributed Integer random number on [min,max], 
// where min, max is the range of all possible values of Integer scalars
MATCL_SCALAR_EXPORT 
Integer             irand(const rand_state& rand_ptr = global_rand_state());

// generate uniformly distributed Complex random number as Complex(rand(),rand())
MATCL_SCALAR_EXPORT 
Complex             crand(const rand_state& rand_ptr = global_rand_state());

// generate uniformly distributed Float_complex random number as 
// Float_complex(frand(),frand())
MATCL_SCALAR_EXPORT 
Float_complex       fcrand(const rand_state& rand_ptr = global_rand_state());

// generate normally distributed Real random number
MATCL_SCALAR_EXPORT 
Real                randn(const rand_state& rand_ptr = global_rand_state());

// generate normally distributed Float random number
MATCL_SCALAR_EXPORT 
Float               frandn(const rand_state& rand_ptr = global_rand_state());

// generate normally distributed Complex random number as Complex(randn(),randn())
MATCL_SCALAR_EXPORT 
Complex             crandn(const rand_state& rand_ptr = global_rand_state());

// generate normally distributed Float_complex random number as 
// Float_complex(frandn(),frandn())
MATCL_SCALAR_EXPORT 
Float_complex       fcrandn(const rand_state& rand_ptr = global_rand_state());

};

#pragma warning(pop)