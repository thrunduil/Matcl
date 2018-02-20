/*
 *  This file is a part of Matrix Computation Library (MATCL)
 *
 *  Copyright (c) Pawe³ Kowal 2017 - 2018
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

namespace matcl { namespace optim_params
{
    //TODO: check this
    static const double max_sparse_density_min      = .25;
    static const double max_sparse_density_max      = .6;
    static const long   min_size_to_test_sparsity   = 9;
    static const long   n_sample_to_test_sparsity   = 23;

    static const Real   mult_cost_dense_dense       = 1;
    static const Real   mult_cost_dense_sparse      = 2;
    static const Real   mult_cost_dense_banded      = 1.5;

    static const Real   mult_cost_sparse_dense      = 2;
    static const Real   mult_cost_sparse_sparse     = 3;
    static const Real   mult_cost_sparse_banded     = 2.5;

    static const Real   mult_cost_banded_dense      = 1.5;
    static const Real   mult_cost_banded_sparse     = 2.5;
    static const Real   mult_cost_banded_banded     = 2;

    static const int    block_size_trans            = 32;
};};
