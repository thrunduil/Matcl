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

#pragma once

#include "matcl-matrep/matcl_matrep.h"

namespace matcl { namespace test
{

struct test_options
{
    static void     set_seed(Integer seed);
    static Integer  get_seed();

    static Real     get_binmat_prob();
    static void     set_binmat_prob(Real prob);

    static Integer  num_colons();
    static void     set_num_colons(Integer val);

    static Integer  num_scalar_groups();
    static void     set_num_scalar_groups(Integer val);

    static Integer  num_matrix_groups();
    static void     set_num_matrix_groups(Integer val);
};

};};