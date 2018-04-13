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

#pragma once

#include "test/test_matcl/framework/matrix_set/matrix_set.h"

namespace matcl { namespace test
{

class mat_set_bin_1 : public matrix_set_bin
{
    private:
        using dims_info = std::vector<Integer>;

    public:
        mat_set_bin_1(rand_matrix_ptr rand)
            :matrix_set_bin(rand, get_sparse_info(), get_band_info())
        {
            init_matrices();
        };

    private:
        void                init_matrices();
        static sparse_info  get_sparse_info();
        static band_info    get_band_info();
        static dims_info    get_dims_info();
};

class mat_set_bin_2 : public matrix_set_bin
{
    public:
        mat_set_bin_2(rand_matrix_ptr rand)
            :matrix_set_bin(rand, get_sparse_info(), get_band_info())
        {
            init_matrices();
        };

    private:
        void                init_matrices();
        static sparse_info  get_sparse_info();
        static band_info    get_band_info();
};

class mat_set_bin_mult_1 : public matrix_set_bin
{
    private:
        using dims_info = std::vector<Integer>;

    public:
        mat_set_bin_mult_1(rand_matrix_ptr rand)
            :matrix_set_bin(rand, get_sparse_info(), get_band_info())
        {
            init_matrices();
        };

    private:
        void                init_matrices();
        static sparse_info  get_sparse_info();
        static band_info    get_band_info();
        static dims_info    get_dims_info();
};

class mat_set_bin_mult_2 : public matrix_set_bin
{
    private:
        using dims_info = std::vector<Integer>;

    public:
        mat_set_bin_mult_2(rand_matrix_ptr rand)
            :matrix_set_bin(rand, get_sparse_info(), get_band_info())
        {
            init_matrices();
        };

    private:
        void                init_matrices();
        static sparse_info  get_sparse_info();
        static band_info    get_band_info();
        static dims_info    get_dims_info();
};

class mat_set_bin_kron_1 : public matrix_set_bin
{
    private:
        using dims_info = std::vector<Integer>;

    public:
        mat_set_bin_kron_1(rand_matrix_ptr rand)
            :matrix_set_bin(rand, get_sparse_info(), get_band_info())
        {
            init_matrices();
        };

    private:
        void                init_matrices();
        static sparse_info  get_sparse_info();
        static band_info    get_band_info();
        static dims_info    get_dims_info();
};

class mat_set_bin_kron_2 : public matrix_set_bin
{
    private:
        using dims_info = std::vector<Integer>;

    public:
        mat_set_bin_kron_2(rand_matrix_ptr rand)
            :matrix_set_bin(rand, get_sparse_info(), get_band_info())
        {
            init_matrices();
        };

    private:
        void                init_matrices();
        static sparse_info  get_sparse_info();
        static band_info    get_band_info();
        static dims_info    get_dims_info();
};

};};