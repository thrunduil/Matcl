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

#include "matrix_set_1.h"
#include "special_cases.h"
#include "test/test_matcl/framework/matrix_set/test_options.h"

namespace matcl { namespace test
{

//======================================================================
//                      mat_set_1
//======================================================================
void mat_set_1::init_matrices()
{
    m_matrices.reserve(2000);

    if (m_rand->generate_random_matrices())
    {
        for (int i = 0; i < 500; ++i)
            m_matrices.push_back(m_rand->rand_any_matrix());
    }
    else
    {
        int n_matrix = matcl::test::test_options::num_matrix_groups();

        std::vector<Integer> dims = get_dims_info();

        for (size_t i = 0; i < dims.size(); ++i)
        {
            for (size_t j = 0; j < dims.size(); ++j)
            {
                for (int k = 0; k < n_matrix; ++k)
                    add_matrices(m_matrices,dims[i], dims[j]);
            };
        };	

        special_cases::add_special_cases(m_matrices);

        int n_scalars = matcl::test::test_options::num_scalar_groups();

        for (int i = 0; i < n_scalars; ++i)
        {
            m_scalars.push_back(Scalar(m_rand->rand_scalar_int()));
            m_scalars.push_back(Scalar(m_rand->rand_scalar_real()));
            m_scalars.push_back(Scalar(m_rand->rand_scalar_float()));
            m_scalars.push_back(Scalar(m_rand->rand_scalar_compl()));
            m_scalars.push_back(Scalar(m_rand->rand_scalar_fcompl()));
        }
    };
};

mat_set_1::sparse_info mat_set_1::get_sparse_info()
{
    sparse_info out;
    out.push_back(.2);
    return out;
};

mat_set_1::band_info mat_set_1::get_band_info()
{
    band_info out;
    out.push_back(int2(0,0));
    out.push_back(int2(0,1));
    out.push_back(int2(1,0));
    out.push_back(int2(1,1));
    out.push_back(int2(3,2)); // to have ud + ld + 1 > N which requires broadening of matrices
    out.push_back(int2(160,0));
    out.push_back(int2(0,160));
    out.push_back(int2(80,90));
    out.push_back(int2(70,60));
    return out;
};

mat_set_1::dims_info mat_set_1::get_dims_info()
{
    dims_info out;
    out.push_back(0);
    out.push_back(1);
    out.push_back(5);
    out.push_back(7);
    out.push_back(15);
    out.push_back(200);
    return out;
};

//======================================================================
//                      mat_set_2
//======================================================================
void mat_set_2::init_matrices()
{
    m_matrices.reserve(20);

    int n_matrix = matcl::test::test_options::num_matrix_groups();

    for (int k = 0; k < n_matrix; ++k)
        add_matrices(m_matrices,5,5);

    int n_scalars = matcl::test::test_options::num_scalar_groups();

    for (int i = 0; i < n_scalars; ++i)
    {
        m_scalars.push_back(Scalar(m_rand->rand_scalar_int()));
        m_scalars.push_back(Scalar(m_rand->rand_scalar_real()));
        m_scalars.push_back(Scalar(m_rand->rand_scalar_float()));
        m_scalars.push_back(Scalar(m_rand->rand_scalar_compl()));
        m_scalars.push_back(Scalar(m_rand->rand_scalar_fcompl()));
    };
};

mat_set_2::sparse_info mat_set_2::get_sparse_info()
{
    sparse_info out;
    out.push_back(.2);
    return out;
};
mat_set_2::band_info mat_set_2::get_band_info()
{
    band_info out;
    out.push_back(int2(0,0));
    out.push_back(int2(0,1));
    out.push_back(int2(0,2));
    out.push_back(int2(1,0));
    out.push_back(int2(1,1));
    out.push_back(int2(2,0));
    out.push_back(int2(3,1));
    out.push_back(int2(1,3));
    return out;
};

};};
