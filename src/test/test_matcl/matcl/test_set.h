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

#include "test/test_matcl/framework/matrix_set/matrix_rand.h"

namespace matcl { namespace test
{

void    test_utils_st(const rand_matrix_ptr& rand);
void    test_utils_mt(const rand_matrix_ptr& rand);

void    test_matgen_st(const rand_matrix_ptr& rand);
void    test_matgen_mt(const rand_matrix_ptr& rand);

void    test_manip_st(const rand_matrix_ptr& rand);
void    test_manip_mt(const rand_matrix_ptr& rand);

void    test_vecfunc_st(const rand_matrix_ptr& rand);
void    test_vecfunc_mt(const rand_matrix_ptr& rand);

void    test_unary_st(const rand_matrix_ptr& rand);
void    test_unary_mt(const rand_matrix_ptr& rand);

void    test_matfunc_st(const rand_matrix_ptr& rand);
void    test_matfunc_mt(const rand_matrix_ptr& rand);
void    test_kron_st(const rand_matrix_ptr& rand);
void    test_kron_mt(const rand_matrix_ptr& rand);

void    test_assign_st(const rand_matrix_ptr& rand);
void    test_assign_mt(const rand_matrix_ptr& rand);

void    test_matrix_st(const rand_matrix_ptr& rand);
void    test_matrix_mt(const rand_matrix_ptr& rand);

void    test_bin_st(const rand_matrix_ptr& rand);
void    test_bin_mt(const rand_matrix_ptr& rand);

void    test_mult_st(const rand_matrix_ptr& rand);
void    test_mult_mt(const rand_matrix_ptr& rand);

void    test_io_st(const rand_matrix_ptr& rand);
void    test_io_mt(const rand_matrix_ptr& rand);

void    test_all(const rand_matrix_ptr& rand);

};};