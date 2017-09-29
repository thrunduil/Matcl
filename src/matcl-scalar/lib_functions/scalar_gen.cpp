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

#include "matcl-scalar/lib_functions/scalar_gen.h"
#include "matcl-scalar/func/raw/rand.h"
#include "matcl-scalar/extern/sfmt_wrap.h"
#include "matcl-core/details/complex_details.h"

namespace matcl
{

//-----------------------------------------------------------------
//                      rand
//-----------------------------------------------------------------
Real matcl::rand(const matcl::rand_state& rand_ptr)
{
    return rand_ptr->genrand_real3();
};
Integer matcl::irand(const matcl::rand_state& rand_ptr)
{
    return rand_ptr->genrand_int32();
};
Float matcl::frand(const matcl::rand_state& rand_ptr)
{
    return raw::frand(rand_ptr);
};
Complex matcl::crand(const matcl::rand_state& rand_ptr)
{
    return Complex(raw::rand(rand_ptr),raw::rand(rand_ptr));
};
Float_complex matcl::fcrand(const matcl::rand_state& rand_ptr)
{
    return Float_complex(raw::frand(rand_ptr),raw::frand(rand_ptr));
};

Real matcl::randn(const matcl::rand_state& rand_ptr)
{
    return raw::randn(rand_ptr);
};
Float matcl::frandn(const matcl::rand_state& rand_ptr)
{
    return raw::frandn(rand_ptr);
};
Complex matcl::crandn(const matcl::rand_state& rand_ptr)
{
    return Complex(raw::randn(rand_ptr),raw::randn(rand_ptr));
};
Float_complex matcl::fcrandn(const matcl::rand_state& rand_ptr)
{
    return Float_complex(raw::frandn(rand_ptr),raw::frandn(rand_ptr));
};

//-----------------------------------------------------------------
//                      rand state
//-----------------------------------------------------------------
rand_state::rand_state()
    :m_ptr(global_rand_state().m_ptr)
{};
rand_state::~rand_state()
{};
rand_state::rand_state(const rand_state_ptr& ptr)
    :m_ptr(ptr)
{};

void matcl::init_genrand(unsigned long s, const matcl::rand_state& rand_ptr)
{
    rand_ptr->init_genrand(s);
};

const matcl::rand_state& matcl::global_rand_state()
{
    return details::global_rand_state();
};

details::local_rand_state_raii matcl::local_rand_state(unsigned long s)
{
    return details::local_rand_state_raii(s);
};

namespace details
{
    class local_rand_state_raii_impl
    {
        private:
            matcl::rand_state  m_rand_state;

        public:
            local_rand_state_raii_impl(unsigned long s)
                :m_rand_state(get_rand_state())
            {
                matcl::init_genrand(s);
            };
            ~local_rand_state_raii_impl()
            {
                matcl::set_rand_state(m_rand_state);
            };
    };
};

details::local_rand_state_raii::local_rand_state_raii(unsigned long s)
    :m_impl(new details::local_rand_state_raii_impl(s))
{};
details::local_rand_state_raii::~local_rand_state_raii()
{};
details::local_rand_state_raii::local_rand_state_raii(details::local_rand_state_raii&& other)
    :m_impl(std::move(other.m_impl))
{};

matcl::rand_state matcl::get_rand_state()
{
    return details::get_rand_state();
};
void matcl::set_rand_state(matcl::rand_state st)
{
    details::set_rand_state(st);
};

};
