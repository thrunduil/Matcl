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
#include "matcl-core/memory/alloc.h"
#include "matcl-core/general/thread.h"
#include "matcl-core/memory/global_objects.h"

#include "sfmt_wrap.h"

namespace matcl { namespace details
{

namespace md = matcl::details;

rand_state::rand_state()
{
    sfmt_init_gen_rand(&m_state_i,0);
    dsfmt_init_gen_rand(&m_state_d,0);
    m_flag = 0;
};

rand_state::rand_state(Integer seed)
{
    sfmt_init_gen_rand(&m_state_i,seed);
    dsfmt_init_gen_rand(&m_state_d,seed);
    m_flag = 0;
};

void rand_state::init_genrand(Integer seed)
{
    sfmt_init_gen_rand(&m_state_i,seed);
    dsfmt_init_gen_rand(&m_state_d,seed);
    m_flag = 0;
};

Integer rand_state::genrand_int32()
{
    return sfmt_genrand_uint32(&m_state_i);
};

double rand_state::genrand_real2()
{
    return dsfmt_genrand_close_open(&m_state_d);
}
double rand_state::genrand_real3()
{
    return dsfmt_genrand_open_open(&m_state_d);
}

double rand_state::gen_norm()
{
    double u, v, s;

    if (m_flag) 
    {
        m_flag = 0;
        return m_b;
    }

    do 
    {
        u = 2. * genrand_real3() - 1.;
        v = 2. * genrand_real3() - 1.;
    } 
    while (((s = u * u + v * v) >= 1.) || (s == 0.));

    s   = sqrt(-2. * log(s) / s);

    m_b     = v * s;
    m_flag  = 1;

    return u * s;
}

struct rand_state_data : global_object, matcl_new_delete
{
    matcl::rand_state m_data;

    rand_state_data(const matcl::rand_state::rand_state_ptr& ptr);
    ~rand_state_data();

    void    clear_global() override;
    void    close_global() override;
};

static void global_rand_state_deleter(rand_state* ptr)
{
    using allocator = default_allocator<true>;
    allocator::aligned_free(ptr, sizeof(rand_state));
};

rand_state_data::rand_state_data(const matcl::rand_state::rand_state_ptr& ptr)
    :m_data(ptr)
{    
};

rand_state_data::~rand_state_data()
{
};

void rand_state_data::clear_global()
{
    m_data = matcl::rand_state();
};

void rand_state_data::close_global()
{
    delete this;
}

static void rand_state_deleter(rand_state* ptr)
{
    using allocator = default_allocator<true>;
    allocator::aligned_free(ptr, sizeof(rand_state));
};

MATCL_THREAD_LOCAL 
rand_state_data* global_rs = nullptr;

static void init_global_rand_state()
{
    using allocator = default_allocator<true>;
    rand_state* rs  = (rand_state*)allocator::aligned_malloc(sizeof(rand_state));
    auto ptr        = matcl::rand_state::rand_state_ptr(rs, &global_rand_state_deleter);
    ptr->init_genrand(0);

    global_rs = new rand_state_data(ptr);
};

const matcl::rand_state& md::global_rand_state()
{
    if (!global_rs)
        init_global_rand_state();

    return global_rs->m_data;
};

matcl::rand_state md::get_rand_state()
{
    if (!global_rs)
        init_global_rand_state();

    matcl::rand_state ret;

    using allocator = default_allocator<true>;
    rand_state* rs  = (rand_state*)allocator::aligned_malloc(sizeof(rand_state));
    ret.m_ptr       = matcl::rand_state::rand_state_ptr(rs, &rand_state_deleter);
    memcpy(rs,global_rs->m_data.m_ptr.get(),sizeof(rand_state));

    return ret;
};

void md::set_rand_state(const matcl::rand_state& st)
{
    if (!global_rs)
        init_global_rand_state();

    memcpy(global_rs->m_data.m_ptr.get(),st.m_ptr.get(),sizeof(rand_state));
};

};};

