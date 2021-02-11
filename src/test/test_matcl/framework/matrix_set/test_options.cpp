/*
 *  This file is a part of Matrix Computation Library (MATCL)
 *
 *  Copyright (c) Pawe³ Kowal 2018 - 2021
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

#include "test_options.h"

namespace matcl { namespace test
{

struct rand_params
{
    Integer m_seed;
    Real    m_binmat_prob;
    Integer m_num_colons;
    Integer m_num_scalar_groups;
    Integer m_num_matrix_groups;

    rand_params()
        :m_seed(0), m_binmat_prob(1.0), m_num_colons(5)
        , m_num_scalar_groups(1), m_num_matrix_groups(1)
    {};
};

static rand_params m_rand_params;

void test_options::set_seed(Integer seed)
{
    m_rand_params.m_seed = seed;
};

Integer test_options::get_seed()
{
    return m_rand_params.m_seed;
};

Integer test_options::num_scalar_groups()
{
    return m_rand_params.m_num_scalar_groups;
}

void test_options::set_num_scalar_groups(Integer val)
{
    m_rand_params.m_num_scalar_groups = val;
}

Integer test_options::num_matrix_groups()
{
    return m_rand_params.m_num_matrix_groups;
}

void test_options::set_num_matrix_groups(Integer val)
{
    m_rand_params.m_num_matrix_groups = val;
}

Real test_options::get_binmat_prob()
{
    return m_rand_params.m_binmat_prob;
};

void test_options::set_binmat_prob(Real prob)
{
    m_rand_params.m_binmat_prob = prob;
}

Integer test_options::num_colons()
{
    return m_rand_params.m_num_colons;
};
void test_options::set_num_colons(Integer val)
{
    m_rand_params.m_num_colons = val;
}

};};