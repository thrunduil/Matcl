/*
 *  This file is a part of Matrix Computation Library (MATCL)
 *
 *  Copyright (c) Pawe³ Kowal 2018
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

#include "matrix_set_bin_1.h"
#include "test/test_matcl/framework/matrix_set/test_options.h"

namespace matcl { namespace test
{

//======================================================================
//                      mat_set_bin_1
//======================================================================
void mat_set_bin_1::init_matrices()
{
    m_matrices.reserve(2000);

    std::vector<Integer> dims = get_dims_info();

    for (size_t i = 0; i < dims.size(); ++i)
    {
        for (size_t j = 0; j < dims.size(); ++j)
        {
            add_matrices(m_matrices,dims[i], dims[j], dims[i], dims[j]);
            // cases needed by sylvester equation
            add_matrices(m_matrices,dims[i], dims[i], dims[j], dims[j]);
            // cases needed by discrete algebraic riccati equation
            add_matrices(m_matrices,dims[i], dims[i], dims[i], dims[j]);
        };
    };

    int n_scalars = matcl::test::test_options::num_scalar_groups();

    for (int i = 0; i < n_scalars; ++i)
    {
        m_scalars.push_back(scalar_pair(m_rand->rand_scalar_int(),m_rand->rand_scalar_int()));
        m_scalars.push_back(scalar_pair(m_rand->rand_scalar_int(),m_rand->rand_scalar_real()));
        m_scalars.push_back(scalar_pair(m_rand->rand_scalar_int(),m_rand->rand_scalar_float()));
        m_scalars.push_back(scalar_pair(m_rand->rand_scalar_int(),m_rand->rand_scalar_compl()));
        m_scalars.push_back(scalar_pair(m_rand->rand_scalar_int(),m_rand->rand_scalar_fcompl()));

        m_scalars.push_back(scalar_pair(m_rand->rand_scalar_real(),m_rand->rand_scalar_int()));
        m_scalars.push_back(scalar_pair(m_rand->rand_scalar_real(),m_rand->rand_scalar_real()));
        m_scalars.push_back(scalar_pair(m_rand->rand_scalar_real(),m_rand->rand_scalar_float()));
        m_scalars.push_back(scalar_pair(m_rand->rand_scalar_real(),m_rand->rand_scalar_compl()));
        m_scalars.push_back(scalar_pair(m_rand->rand_scalar_real(),m_rand->rand_scalar_fcompl()));

        m_scalars.push_back(scalar_pair(m_rand->rand_scalar_float(),m_rand->rand_scalar_int()));
        m_scalars.push_back(scalar_pair(m_rand->rand_scalar_float(),m_rand->rand_scalar_real()));
        m_scalars.push_back(scalar_pair(m_rand->rand_scalar_float(),m_rand->rand_scalar_float()));
        m_scalars.push_back(scalar_pair(m_rand->rand_scalar_float(),m_rand->rand_scalar_compl()));
        m_scalars.push_back(scalar_pair(m_rand->rand_scalar_float(),m_rand->rand_scalar_fcompl()));

        m_scalars.push_back(scalar_pair(m_rand->rand_scalar_compl(),m_rand->rand_scalar_int()));
        m_scalars.push_back(scalar_pair(m_rand->rand_scalar_compl(),m_rand->rand_scalar_real()));
        m_scalars.push_back(scalar_pair(m_rand->rand_scalar_compl(),m_rand->rand_scalar_float()));
        m_scalars.push_back(scalar_pair(m_rand->rand_scalar_compl(),m_rand->rand_scalar_compl()));
        m_scalars.push_back(scalar_pair(m_rand->rand_scalar_compl(),m_rand->rand_scalar_fcompl()));

        m_scalars.push_back(scalar_pair(m_rand->rand_scalar_fcompl(),m_rand->rand_scalar_int()));
        m_scalars.push_back(scalar_pair(m_rand->rand_scalar_fcompl(),m_rand->rand_scalar_real()));
        m_scalars.push_back(scalar_pair(m_rand->rand_scalar_fcompl(),m_rand->rand_scalar_float()));
        m_scalars.push_back(scalar_pair(m_rand->rand_scalar_fcompl(),m_rand->rand_scalar_compl()));
        m_scalars.push_back(scalar_pair(m_rand->rand_scalar_fcompl(),m_rand->rand_scalar_fcompl()));
    };
};

mat_set_bin_1::sparse_info mat_set_bin_1::get_sparse_info()
{
    sparse_info out;
    out.push_back(.2);
    return out;
};
mat_set_bin_1::band_info mat_set_bin_1::get_band_info()
{
    band_info out;
    out.push_back(int2(0,0));
    out.push_back(int2(0,1));
    out.push_back(int2(0,2));
    out.push_back(int2(1,0));
    out.push_back(int2(2,0));
    out.push_back(int2(1,1));
    out.push_back(int2(1,2));
    out.push_back(int2(2,1));
    out.push_back(int2(2,2));
    return out;
};

mat_set_bin_1::dims_info mat_set_bin_1::get_dims_info()
{
    dims_info out;
    out.push_back(0);
    out.push_back(1);
    out.push_back(3);
    out.push_back(5);
    out.push_back(6);
    out.push_back(7);
    out.push_back(10);
    return out;
};

//======================================================================
//                      mat_set_bin_2
//======================================================================
void mat_set_bin_2::init_matrices()
{
    m_matrices.reserve(500);
    
    add_matrices(m_matrices,5,5,5,5);
    // cases needed by sylvester equation
    add_matrices(m_matrices,5,5,6,6);
    // cases needed by discrete algebraic riccati equation
    add_matrices(m_matrices,5,5,5,6);

    int n_scalars = matcl::test::test_options::num_scalar_groups();

    for (int i = 0; i < n_scalars; ++i)
    {
        m_scalars.push_back(scalar_pair(m_rand->rand_scalar_int(),m_rand->rand_scalar_int()));
        m_scalars.push_back(scalar_pair(m_rand->rand_scalar_int(),m_rand->rand_scalar_real()));
        m_scalars.push_back(scalar_pair(m_rand->rand_scalar_int(),m_rand->rand_scalar_float()));
        m_scalars.push_back(scalar_pair(m_rand->rand_scalar_int(),m_rand->rand_scalar_compl()));
        m_scalars.push_back(scalar_pair(m_rand->rand_scalar_int(),m_rand->rand_scalar_fcompl()));

        m_scalars.push_back(scalar_pair(m_rand->rand_scalar_real(),m_rand->rand_scalar_int()));
        m_scalars.push_back(scalar_pair(m_rand->rand_scalar_real(),m_rand->rand_scalar_real()));
        m_scalars.push_back(scalar_pair(m_rand->rand_scalar_real(),m_rand->rand_scalar_float()));
        m_scalars.push_back(scalar_pair(m_rand->rand_scalar_real(),m_rand->rand_scalar_compl()));
        m_scalars.push_back(scalar_pair(m_rand->rand_scalar_real(),m_rand->rand_scalar_fcompl()));

        m_scalars.push_back(scalar_pair(m_rand->rand_scalar_float(),m_rand->rand_scalar_int()));
        m_scalars.push_back(scalar_pair(m_rand->rand_scalar_float(),m_rand->rand_scalar_real()));
        m_scalars.push_back(scalar_pair(m_rand->rand_scalar_float(),m_rand->rand_scalar_float()));
        m_scalars.push_back(scalar_pair(m_rand->rand_scalar_float(),m_rand->rand_scalar_compl()));
        m_scalars.push_back(scalar_pair(m_rand->rand_scalar_float(),m_rand->rand_scalar_fcompl()));

        m_scalars.push_back(scalar_pair(m_rand->rand_scalar_compl(),m_rand->rand_scalar_int()));
        m_scalars.push_back(scalar_pair(m_rand->rand_scalar_compl(),m_rand->rand_scalar_real()));
        m_scalars.push_back(scalar_pair(m_rand->rand_scalar_compl(),m_rand->rand_scalar_float()));
        m_scalars.push_back(scalar_pair(m_rand->rand_scalar_compl(),m_rand->rand_scalar_compl()));
        m_scalars.push_back(scalar_pair(m_rand->rand_scalar_compl(),m_rand->rand_scalar_fcompl()));

        m_scalars.push_back(scalar_pair(m_rand->rand_scalar_fcompl(),m_rand->rand_scalar_int()));
        m_scalars.push_back(scalar_pair(m_rand->rand_scalar_fcompl(),m_rand->rand_scalar_real()));
        m_scalars.push_back(scalar_pair(m_rand->rand_scalar_fcompl(),m_rand->rand_scalar_float()));
        m_scalars.push_back(scalar_pair(m_rand->rand_scalar_fcompl(),m_rand->rand_scalar_compl()));
        m_scalars.push_back(scalar_pair(m_rand->rand_scalar_fcompl(),m_rand->rand_scalar_fcompl()));
    };
};

mat_set_bin_2::sparse_info mat_set_bin_2::get_sparse_info()
{
    sparse_info out;
    out.push_back(.2);
    return out;
};
mat_set_bin_2::band_info mat_set_bin_2::get_band_info()
{
    band_info out;
    out.push_back(int2(0,0));
    out.push_back(int2(0,1));
    out.push_back(int2(1,0));
    out.push_back(int2(1,1));
    return out;
};


//======================================================================
//                      mat_set_bin_mult_1
//======================================================================
void mat_set_bin_mult_1::init_matrices()
{
    m_matrices.reserve(100000);
    std::vector<Integer> dims = get_dims_info();

    for (size_t i = 0; i < dims.size(); ++i)
    {
        for (size_t j = 0; j < dims.size(); ++j)
        {
            for (size_t k = 0; k < dims.size(); ++k)
                add_matrices(m_matrices,dims[i], dims[k], dims[k], dims[j]);
        };
    };

    int n_scalars = matcl::test::test_options::num_scalar_groups();

    for (int i = 0; i < n_scalars; ++i)
    {
        m_scalars.push_back(scalar_pair(m_rand->rand_scalar_int(),m_rand->rand_scalar_int()));
        m_scalars.push_back(scalar_pair(m_rand->rand_scalar_int(),m_rand->rand_scalar_real()));
        m_scalars.push_back(scalar_pair(m_rand->rand_scalar_int(),m_rand->rand_scalar_float()));
        m_scalars.push_back(scalar_pair(m_rand->rand_scalar_int(),m_rand->rand_scalar_compl()));
        m_scalars.push_back(scalar_pair(m_rand->rand_scalar_int(),m_rand->rand_scalar_fcompl()));

        m_scalars.push_back(scalar_pair(m_rand->rand_scalar_real(),m_rand->rand_scalar_int()));
        m_scalars.push_back(scalar_pair(m_rand->rand_scalar_real(),m_rand->rand_scalar_real()));
        m_scalars.push_back(scalar_pair(m_rand->rand_scalar_real(),m_rand->rand_scalar_float()));
        m_scalars.push_back(scalar_pair(m_rand->rand_scalar_real(),m_rand->rand_scalar_compl()));
        m_scalars.push_back(scalar_pair(m_rand->rand_scalar_real(),m_rand->rand_scalar_fcompl()));

        m_scalars.push_back(scalar_pair(m_rand->rand_scalar_float(),m_rand->rand_scalar_int()));
        m_scalars.push_back(scalar_pair(m_rand->rand_scalar_float(),m_rand->rand_scalar_real()));
        m_scalars.push_back(scalar_pair(m_rand->rand_scalar_float(),m_rand->rand_scalar_float()));
        m_scalars.push_back(scalar_pair(m_rand->rand_scalar_float(),m_rand->rand_scalar_compl()));
        m_scalars.push_back(scalar_pair(m_rand->rand_scalar_float(),m_rand->rand_scalar_fcompl()));

        m_scalars.push_back(scalar_pair(m_rand->rand_scalar_compl(),m_rand->rand_scalar_int()));
        m_scalars.push_back(scalar_pair(m_rand->rand_scalar_compl(),m_rand->rand_scalar_real()));
        m_scalars.push_back(scalar_pair(m_rand->rand_scalar_compl(),m_rand->rand_scalar_float()));
        m_scalars.push_back(scalar_pair(m_rand->rand_scalar_compl(),m_rand->rand_scalar_compl()));
        m_scalars.push_back(scalar_pair(m_rand->rand_scalar_compl(),m_rand->rand_scalar_fcompl()));

        m_scalars.push_back(scalar_pair(m_rand->rand_scalar_fcompl(),m_rand->rand_scalar_int()));
        m_scalars.push_back(scalar_pair(m_rand->rand_scalar_fcompl(),m_rand->rand_scalar_real()));
        m_scalars.push_back(scalar_pair(m_rand->rand_scalar_fcompl(),m_rand->rand_scalar_float()));
        m_scalars.push_back(scalar_pair(m_rand->rand_scalar_fcompl(),m_rand->rand_scalar_compl()));
        m_scalars.push_back(scalar_pair(m_rand->rand_scalar_fcompl(),m_rand->rand_scalar_fcompl()));
    };
};

mat_set_bin_mult_1::dims_info mat_set_bin_mult_1::get_dims_info()
{
    dims_info out;
    out.push_back(0);
    out.push_back(1);
    out.push_back(3);
    out.push_back(5);
    out.push_back(7);
    out.push_back(10);
    return out;
};

mat_set_bin_mult_1::sparse_info mat_set_bin_mult_1::get_sparse_info()
{
    sparse_info out;
    out.push_back(.2);
    return out;
};
mat_set_bin_mult_1::band_info mat_set_bin_mult_1::get_band_info()
{
    band_info out;
    out.push_back(int2(0,0));
    out.push_back(int2(0,1));
    out.push_back(int2(0,2));
    out.push_back(int2(1,0));
    out.push_back(int2(2,0));
    out.push_back(int2(1,1));
    out.push_back(int2(1,2));
    out.push_back(int2(2,1));
    out.push_back(int2(2,2));
    return out;
};

//======================================================================
//                      mat_set_bin_mult_2
//======================================================================
void mat_set_bin_mult_2::init_matrices()
{
    m_matrices.reserve(2000);
    std::vector<Integer> dims = get_dims_info();

    for (size_t i = 0; i < dims.size(); ++i)
        add_matrices(m_matrices, 5, dims[i], dims[i], 5);

    int n_scalars = matcl::test::test_options::num_scalar_groups();

    for (int i = 0; i < n_scalars; ++i)
    {
        m_scalars.push_back(scalar_pair(m_rand->rand_scalar_int(),m_rand->rand_scalar_int()));
        m_scalars.push_back(scalar_pair(m_rand->rand_scalar_int(),m_rand->rand_scalar_real()));
        m_scalars.push_back(scalar_pair(m_rand->rand_scalar_int(),m_rand->rand_scalar_float()));
        m_scalars.push_back(scalar_pair(m_rand->rand_scalar_int(),m_rand->rand_scalar_compl()));
        m_scalars.push_back(scalar_pair(m_rand->rand_scalar_int(),m_rand->rand_scalar_fcompl()));

        m_scalars.push_back(scalar_pair(m_rand->rand_scalar_real(),m_rand->rand_scalar_int()));
        m_scalars.push_back(scalar_pair(m_rand->rand_scalar_real(),m_rand->rand_scalar_real()));
        m_scalars.push_back(scalar_pair(m_rand->rand_scalar_real(),m_rand->rand_scalar_float()));
        m_scalars.push_back(scalar_pair(m_rand->rand_scalar_real(),m_rand->rand_scalar_compl()));
        m_scalars.push_back(scalar_pair(m_rand->rand_scalar_real(),m_rand->rand_scalar_fcompl()));

        m_scalars.push_back(scalar_pair(m_rand->rand_scalar_float(),m_rand->rand_scalar_int()));
        m_scalars.push_back(scalar_pair(m_rand->rand_scalar_float(),m_rand->rand_scalar_real()));
        m_scalars.push_back(scalar_pair(m_rand->rand_scalar_float(),m_rand->rand_scalar_float()));
        m_scalars.push_back(scalar_pair(m_rand->rand_scalar_float(),m_rand->rand_scalar_compl()));
        m_scalars.push_back(scalar_pair(m_rand->rand_scalar_float(),m_rand->rand_scalar_fcompl()));

        m_scalars.push_back(scalar_pair(m_rand->rand_scalar_compl(),m_rand->rand_scalar_int()));
        m_scalars.push_back(scalar_pair(m_rand->rand_scalar_compl(),m_rand->rand_scalar_real()));
        m_scalars.push_back(scalar_pair(m_rand->rand_scalar_compl(),m_rand->rand_scalar_float()));
        m_scalars.push_back(scalar_pair(m_rand->rand_scalar_compl(),m_rand->rand_scalar_compl()));
        m_scalars.push_back(scalar_pair(m_rand->rand_scalar_compl(),m_rand->rand_scalar_fcompl()));

        m_scalars.push_back(scalar_pair(m_rand->rand_scalar_fcompl(),m_rand->rand_scalar_int()));
        m_scalars.push_back(scalar_pair(m_rand->rand_scalar_fcompl(),m_rand->rand_scalar_real()));
        m_scalars.push_back(scalar_pair(m_rand->rand_scalar_fcompl(),m_rand->rand_scalar_float()));
        m_scalars.push_back(scalar_pair(m_rand->rand_scalar_fcompl(),m_rand->rand_scalar_compl()));
        m_scalars.push_back(scalar_pair(m_rand->rand_scalar_fcompl(),m_rand->rand_scalar_fcompl()));
    };
};

mat_set_bin_mult_2::dims_info mat_set_bin_mult_2::get_dims_info()
{
    dims_info out;
    out.push_back(0);
    out.push_back(1);
    out.push_back(5);
    out.push_back(7);
    out.push_back(10);
    return out;
};
mat_set_bin_mult_2::sparse_info mat_set_bin_mult_2::get_sparse_info()
{
    sparse_info out;
    out.push_back(.2);
    return out;
};
mat_set_bin_mult_2::band_info mat_set_bin_mult_2::get_band_info()
{
    band_info out;
    out.push_back(int2(0,0));
    out.push_back(int2(0,1));
    out.push_back(int2(1,0));
    out.push_back(int2(1,1));
    return out;
};

//======================================================================
//                      mat_set_bin_kron_1
//======================================================================
void mat_set_bin_kron_1::init_matrices()
{
    m_matrices.reserve(20000);
    std::vector<Integer> dims = get_dims_info();

    for (size_t i = 0; i < dims.size(); ++i)
    {
        for (size_t j = 0; j < dims.size(); ++j)
        {
            for (size_t k = 0; k < dims.size(); ++k)
            {
                for (size_t l = 0; l < dims.size(); ++l)
                    add_matrices(m_matrices,dims[i], dims[j], dims[k], dims[l]);
            };
        };
    };

    int n_scalars = matcl::test::test_options::num_scalar_groups();

    for (int i = 0; i < n_scalars; ++i)
    {
        m_scalars.push_back(scalar_pair(m_rand->rand_scalar_int(),m_rand->rand_scalar_int()));
        m_scalars.push_back(scalar_pair(m_rand->rand_scalar_int(),m_rand->rand_scalar_real()));
        m_scalars.push_back(scalar_pair(m_rand->rand_scalar_int(),m_rand->rand_scalar_float()));
        m_scalars.push_back(scalar_pair(m_rand->rand_scalar_int(),m_rand->rand_scalar_compl()));
        m_scalars.push_back(scalar_pair(m_rand->rand_scalar_int(),m_rand->rand_scalar_fcompl()));

        m_scalars.push_back(scalar_pair(m_rand->rand_scalar_real(),m_rand->rand_scalar_int()));
        m_scalars.push_back(scalar_pair(m_rand->rand_scalar_real(),m_rand->rand_scalar_real()));
        m_scalars.push_back(scalar_pair(m_rand->rand_scalar_real(),m_rand->rand_scalar_float()));
        m_scalars.push_back(scalar_pair(m_rand->rand_scalar_real(),m_rand->rand_scalar_compl()));
        m_scalars.push_back(scalar_pair(m_rand->rand_scalar_real(),m_rand->rand_scalar_fcompl()));

        m_scalars.push_back(scalar_pair(m_rand->rand_scalar_float(),m_rand->rand_scalar_int()));
        m_scalars.push_back(scalar_pair(m_rand->rand_scalar_float(),m_rand->rand_scalar_real()));
        m_scalars.push_back(scalar_pair(m_rand->rand_scalar_float(),m_rand->rand_scalar_float()));
        m_scalars.push_back(scalar_pair(m_rand->rand_scalar_float(),m_rand->rand_scalar_compl()));
        m_scalars.push_back(scalar_pair(m_rand->rand_scalar_float(),m_rand->rand_scalar_fcompl()));

        m_scalars.push_back(scalar_pair(m_rand->rand_scalar_compl(),m_rand->rand_scalar_int()));
        m_scalars.push_back(scalar_pair(m_rand->rand_scalar_compl(),m_rand->rand_scalar_real()));
        m_scalars.push_back(scalar_pair(m_rand->rand_scalar_compl(),m_rand->rand_scalar_float()));
        m_scalars.push_back(scalar_pair(m_rand->rand_scalar_compl(),m_rand->rand_scalar_compl()));
        m_scalars.push_back(scalar_pair(m_rand->rand_scalar_compl(),m_rand->rand_scalar_fcompl()));

        m_scalars.push_back(scalar_pair(m_rand->rand_scalar_fcompl(),m_rand->rand_scalar_int()));
        m_scalars.push_back(scalar_pair(m_rand->rand_scalar_fcompl(),m_rand->rand_scalar_real()));
        m_scalars.push_back(scalar_pair(m_rand->rand_scalar_fcompl(),m_rand->rand_scalar_float()));
        m_scalars.push_back(scalar_pair(m_rand->rand_scalar_fcompl(),m_rand->rand_scalar_compl()));
        m_scalars.push_back(scalar_pair(m_rand->rand_scalar_fcompl(),m_rand->rand_scalar_fcompl()));
    };
};

mat_set_bin_kron_1::dims_info mat_set_bin_kron_1::get_dims_info()
{
    dims_info out;
    out.push_back(0);
    out.push_back(1);
    out.push_back(5);
    out.push_back(7);
    out.push_back(10);
    return out;
};
mat_set_bin_kron_1::sparse_info mat_set_bin_kron_1::get_sparse_info()
{
    sparse_info out;
    out.push_back(.2);
    return out;
};
mat_set_bin_kron_1::band_info mat_set_bin_kron_1::get_band_info()
{
    band_info out;
    out.push_back(int2(0,0));
    out.push_back(int2(0,1));
    out.push_back(int2(1,0));
    out.push_back(int2(1,1));
    return out;
};

//======================================================================
//                      mat_set_bin_kron_2
//======================================================================
void mat_set_bin_kron_2::init_matrices()
{
    m_matrices.reserve(2000);
    std::vector<Integer> dims = get_dims_info();

    for (size_t i = 0; i < dims.size(); ++i)
    {
        for (size_t j = 0; j < dims.size(); ++j)
            add_matrices(m_matrices, 5, dims[i], 5, dims[j]);
    };

    int n_scalars = matcl::test::test_options::num_scalar_groups();

    for (int i = 0; i < n_scalars; ++i)
    {
        m_scalars.push_back(scalar_pair(m_rand->rand_scalar_int(),m_rand->rand_scalar_int()));
        m_scalars.push_back(scalar_pair(m_rand->rand_scalar_int(),m_rand->rand_scalar_real()));
        m_scalars.push_back(scalar_pair(m_rand->rand_scalar_int(),m_rand->rand_scalar_float()));
        m_scalars.push_back(scalar_pair(m_rand->rand_scalar_int(),m_rand->rand_scalar_compl()));
        m_scalars.push_back(scalar_pair(m_rand->rand_scalar_int(),m_rand->rand_scalar_fcompl()));

        m_scalars.push_back(scalar_pair(m_rand->rand_scalar_real(),m_rand->rand_scalar_int()));
        m_scalars.push_back(scalar_pair(m_rand->rand_scalar_real(),m_rand->rand_scalar_real()));
        m_scalars.push_back(scalar_pair(m_rand->rand_scalar_real(),m_rand->rand_scalar_float()));
        m_scalars.push_back(scalar_pair(m_rand->rand_scalar_real(),m_rand->rand_scalar_compl()));
        m_scalars.push_back(scalar_pair(m_rand->rand_scalar_real(),m_rand->rand_scalar_fcompl()));

        m_scalars.push_back(scalar_pair(m_rand->rand_scalar_float(),m_rand->rand_scalar_int()));
        m_scalars.push_back(scalar_pair(m_rand->rand_scalar_float(),m_rand->rand_scalar_real()));
        m_scalars.push_back(scalar_pair(m_rand->rand_scalar_float(),m_rand->rand_scalar_float()));
        m_scalars.push_back(scalar_pair(m_rand->rand_scalar_float(),m_rand->rand_scalar_compl()));
        m_scalars.push_back(scalar_pair(m_rand->rand_scalar_float(),m_rand->rand_scalar_fcompl()));

        m_scalars.push_back(scalar_pair(m_rand->rand_scalar_compl(),m_rand->rand_scalar_int()));
        m_scalars.push_back(scalar_pair(m_rand->rand_scalar_compl(),m_rand->rand_scalar_real()));
        m_scalars.push_back(scalar_pair(m_rand->rand_scalar_compl(),m_rand->rand_scalar_float()));
        m_scalars.push_back(scalar_pair(m_rand->rand_scalar_compl(),m_rand->rand_scalar_compl()));
        m_scalars.push_back(scalar_pair(m_rand->rand_scalar_compl(),m_rand->rand_scalar_fcompl()));

        m_scalars.push_back(scalar_pair(m_rand->rand_scalar_fcompl(),m_rand->rand_scalar_int()));
        m_scalars.push_back(scalar_pair(m_rand->rand_scalar_fcompl(),m_rand->rand_scalar_real()));
        m_scalars.push_back(scalar_pair(m_rand->rand_scalar_fcompl(),m_rand->rand_scalar_float()));
        m_scalars.push_back(scalar_pair(m_rand->rand_scalar_fcompl(),m_rand->rand_scalar_compl()));
        m_scalars.push_back(scalar_pair(m_rand->rand_scalar_fcompl(),m_rand->rand_scalar_fcompl()));
    };
};

mat_set_bin_kron_2::dims_info mat_set_bin_kron_2::get_dims_info()
{
    dims_info out;
    out.push_back(0);
    out.push_back(1);
    out.push_back(5);
    out.push_back(7);
    out.push_back(10);
    return out;
};
mat_set_bin_kron_2::sparse_info mat_set_bin_kron_2::get_sparse_info()
{
    sparse_info out;
    out.push_back(.2);
    return out;
};
mat_set_bin_kron_2::band_info mat_set_bin_kron_2::get_band_info()
{
    band_info out;
    out.push_back(int2(0,0));
    out.push_back(int2(0,1));
    out.push_back(int2(1,0));
    out.push_back(int2(1,1));
    return out;
};

};};

