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

#include "matcl-matrep/details/utils.h"
#include "matcl-internals/container/mat_d.h"
#include "matcl-internals/container/mat_s.h"
#include "matcl-internals/container/mat_b.h"

namespace matcl { namespace raw { namespace details
{

namespace md = matcl::details;

template<class MP>
struct vec_manip_helper
{
    using value_type        = typename MP::value_type;
    using struct_type       = typename MP::struct_type;
    using real_value_type   = typename md::real_type<value_type>::type;
    using value_type_real   = typename md::unify_types<value_type,Float>::type;

    using ret_type_sum      = Matrix<value_type,struct_dense>;
    using ret_type_nnz      = Matrix<Integer,struct_dense>;
    using ret_type_prod     = Matrix<value_type,struct_dense>;
    using ret_type_cumsum   = Matrix<value_type,struct_dense>;
    using ret_type_cumprod  = Matrix<value_type,struct_dense>;
    using ret_type_sumsq    = Matrix<value_type,struct_dense>;
    using ret_type_min      = Matrix<value_type,struct_dense>;
    using ret_type_max      = Matrix<value_type,struct_dense>;
    using ret_type_min_abs  = Matrix<real_value_type,struct_dense>;
    using ret_type_max_abs  = Matrix<real_value_type,struct_dense>;
    using ret_type_mean     = Matrix<value_type_real,struct_dense>;
    using ret_type_std      = Matrix<value_type_real,struct_dense>;
    using ret_type_any      = Matrix<Integer,struct_dense>;
    using ret_type_all      = Matrix<Integer,struct_dense>;
    using ret_type_min_2    = std::pair<ret_type_min,integer_dense>;
    using ret_type_max_2    = std::pair<ret_type_max,integer_dense>;
    using ret_type_min_abs2 = std::pair<ret_type_min_abs,integer_dense>;
    using ret_type_max_abs2 = std::pair<ret_type_max_abs,integer_dense>;

    static void	eval_nnz(matcl::Matrix& ret, const MP& m, int dim);
    static void eval_sum(matcl::Matrix& ret, const MP& m, int dim);
    static void eval_prod(matcl::Matrix& ret, const MP& m, int dim);
    static void eval_cumsum(matcl::Matrix& ret, const MP& m, int dim);
    static void eval_cumprod(matcl::Matrix& ret, const MP& m, int dim);
    static void eval_sumsq(matcl::Matrix& ret, const MP& m, int dim);
    static void eval_min(matcl::Matrix& ret, const MP& m, int dim);
    static void eval_max(matcl::Matrix& ret, const MP& m, int dim);
    static void eval_min_abs(matcl::Matrix& ret, const MP& m, int dim);
    static void eval_max_abs(matcl::Matrix& ret, const MP& m, int dim);
    static void eval_min2(matcl::Matrix& x, matcl::Matrix& i, const MP& m, int dim);
    static void eval_max2(matcl::Matrix& x, matcl::Matrix& i, const MP& m, int dim);
    static void eval_min_abs2(matcl::Matrix& x, matcl::Matrix& i, const MP& m, int dim);
    static void eval_max_abs2(matcl::Matrix& x, matcl::Matrix& i, const MP& m, int dim);
    static void eval_mean(matcl::Matrix& ret, const MP& m, int dim);
    static void eval_std(matcl::Matrix& ret, const MP& m, int dim,bool unbiased);
    static void eval_any(matcl::Matrix& ret, const MP& m, int dim);
    static void eval_all(matcl::Matrix& ret, const MP& m, int dim);
    static void eval_any(matcl::Matrix& ret, const MP& m, const test_function& t,int dim);
    static void eval_all(matcl::Matrix& ret, const MP& m, const test_function& t,int dim);


    static void	eval_nnz_vec(Integer& ret, const MP& m);
    static void eval_sum_vec(matcl::Matrix& ret, const MP& m);
    static void eval_prod_vec(matcl::Matrix& ret, const MP& m);
    static void eval_sumsq_vec(matcl::Matrix& ret, const MP& m);
    static void eval_min_vec(matcl::Matrix& ret, const MP& m);
    static void eval_max_vec(matcl::Matrix& ret, const MP& m);
    static void eval_min_abs_vec(matcl::Matrix& ret, const MP& m);
    static void eval_max_abs_vec(matcl::Matrix& ret, const MP& m);
    static void eval_mean_vec(matcl::Matrix& ret, const MP& m);
    static void eval_std_vec(matcl::Matrix& ret, const MP& m, bool unbiased);
    static void eval_any_vec(bool& ret, const MP& m);
    static void eval_all_vec(bool& ret, const MP& m);
    static void eval_any_vec(bool& ret, const MP& m, const test_function& t);
    static void eval_all_vec(bool& ret, const MP& m, const test_function& t);

    static void eval_min2_vec(matcl::Matrix& x, Integer& i, Integer& j, const MP& m);
    static void eval_max2_vec(matcl::Matrix& x, Integer& i, Integer& j, const MP& m);
    static void eval_min_abs2_vec(matcl::Matrix& x, Integer& i, Integer& j, const MP& m);
    static void eval_max_abs2_vec(matcl::Matrix& x, Integer& i, Integer& j, const MP& m);
};

};}};
