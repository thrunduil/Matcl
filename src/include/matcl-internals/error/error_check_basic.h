/*
 *  This file is a part of Matrix Computation Library (MATCL)
 *
 *  Copyright (c) Pawe³ Kowal 2017 - 2021
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

#include "matcl-matrep/general/config.h"
#include "matcl-core/matrix/scalar_types.h"
#include "matcl-core/error/exception_classes.h"
#include "matcl-core/details/integer.h"
#include "matcl-matrep/details/fwd_decls.h"

namespace matcl { namespace error
{

MATCL_MATREP_EXPORT void throw_error_col(Integer i, Integer r, Integer c);
MATCL_MATREP_EXPORT void throw_error_row(Integer i, Integer r, Integer c);
MATCL_MATREP_EXPORT void throw_error_single_index(Integer i, Integer size);
MATCL_MATREP_EXPORT void throw_error_size(Integer r, Integer c);
MATCL_MATREP_EXPORT void throw_error_alloc(Integer n);

force_inline void check_col(Integer j, Integer r, Integer c)
{
    if ((j < 1) || (j > c))
        throw_error_col(j, r, c);
}

force_inline void check_row(Integer i, Integer r, Integer c)
{
    if ((i < 1) || (i > r))
        throw_error_row(i, r, c);
}

force_inline void check_index(Integer i, Integer size)
{
    if (i < 1 || i > size)
        throw_error_single_index(i,size);
};

force_inline void check_size(Integer r, Integer c)
{
    Integer tmp;
    if ((c < 0) || (r < 0) || imult_t(r,c,tmp) == false)
        throw_error_size(r, c);
}

force_inline void check_alloc(void *ptr, Integer sz, Integer size_of)
{
    if (!ptr)
        throw_error_alloc(imult_c(sz, size_of));
}

MATCL_MATREP_EXPORT void check_diag(Integer d, Integer r, Integer c);
MATCL_MATREP_EXPORT void check_assign_1(Integer s1, Integer r2, Integer c2);
MATCL_MATREP_EXPORT void check_single_index(Integer i , Integer r, Integer c);
MATCL_MATREP_EXPORT void check_index(Integer i , Integer j, Integer r, Integer c);

//return true if p == 1:1:length
MATCL_MATREP_EXPORT bool check_permutation_vector(const matcl::Matrix& p, Integer length);
MATCL_MATREP_EXPORT void check_assign_2(Integer r1, Integer c1, Integer r2, Integer c2);
MATCL_MATREP_EXPORT void check_index_band(Integer i, Integer j, Integer r, Integer c, Integer l,Integer u);
MATCL_MATREP_EXPORT void check_horzcat(Integer r1, Integer c1, Integer r2, Integer c2);
MATCL_MATREP_EXPORT void check_vertcat(Integer r1, Integer c1, Integer r2, Integer c2);
MATCL_MATREP_EXPORT void check_scalar(Integer r,Integer c);
MATCL_MATREP_EXPORT void check_size_band(Integer r, Integer c, Integer l, Integer u);
MATCL_MATREP_EXPORT void check_resize(Integer r, Integer c);
MATCL_MATREP_EXPORT void check_reshape(Integer r1, Integer c1, Integer r2, Integer c2);
MATCL_MATREP_EXPORT void check_size_sp(Integer r, Integer c);
MATCL_MATREP_EXPORT void check_linear_index(Integer r, Integer c, Integer rows);
MATCL_MATREP_EXPORT void check_randperm_arg(Integer n);
MATCL_MATREP_EXPORT void check_bspdiags_2ndarg(Integer r, Integer c);
MATCL_MATREP_EXPORT void check_bspdiag_1starg(Integer r, Integer c);
MATCL_MATREP_EXPORT void check_diag_arg(Integer r, Integer c);
MATCL_MATREP_EXPORT void check_col_indices_sortrows(Integer size, Integer c);
MATCL_MATREP_EXPORT void check_col_indices_elem_sortrows(Integer size, Integer c);
MATCL_MATREP_EXPORT void check_row_indices_sortcols(Integer size, Integer c);
MATCL_MATREP_EXPORT void check_row_indices_elem_sortcols(Integer size, Integer c);
MATCL_MATREP_EXPORT void check_dim(Integer i, Integer d = 2);
MATCL_MATREP_EXPORT void check_eeop(Integer r1, Integer c1, Integer r2, Integer c2);
MATCL_MATREP_EXPORT void check_mul(Integer r1, Integer c1, Integer r2, Integer c2, 
                            trans_type t1, trans_type t2);

};};