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

#include "matcl-internals/error/error_check_basic.h"
#include "matcl-core/error/exception_classes.h"
#include "matcl-core/lib_functions/constants.h"
#include "matcl-internals/base/utils.h"
#include "matcl-core/details/integer.h"
#include "matcl-matrep/matrix/matrix.h"

#include <algorithm>

namespace matcl { namespace error
{

void error::check_single_index(Integer i , Integer r, Integer c)
{
    if (r == 0)
        error::check_index(i, 0);

    Integer j;
    matcl::details::pos2ind(i,r,i,j);
    check_index(i + 1, j + 1, r, c);
};

void error::check_index(Integer i , Integer j, Integer r, Integer c)
{
    if ((i < 1) || (i > r) || (j < 1) || (j > c))
        throw invalid_double_index(i, j, r, c);
};

void error::check_resize(Integer r, Integer c)
{
    if ( c < 0 || r < 0)
        throw invalid_resize(r, c);
};

void error::check_diag(Integer d, Integer r, Integer c)
{
    if (d != 0)
    {
        if (r == 0 || c == 0 || d + 1 > c || 1 - d > r)
            throw invalid_diag(d, r, c);
    };
};

void error::check_index_band(Integer i, Integer j, Integer r, Integer c, Integer l, Integer u)
{
    if ((i < 1) || (j < 1) || (i > r) || (j > c) ||
        (std::max(1, j - u) > i) || (i > std::min(r, j + l)))
    {
        throw invalid_index_band(i, j, r, c, l, u);
    };
}

void error::check_size_band(Integer r, Integer c, Integer l, Integer u)
{
    Integer diag_u = (c == 0)? 0 : 1;
    Integer diag_l = (r == 0)? 0 : 1;

    if (r < 0 || c < 0 || l < 0 || u < 0 || u + diag_u > c || l + diag_l > r)
        throw invalid_size_band(r, c, l, u);
}

void error::check_size_sp(Integer r, Integer c)
{
    if (r < 0 || c < 0)
        throw invalid_size_sp(r, c);
}

void error::check_horzcat(Integer r1, Integer c1, Integer r2, Integer c2)
{
    if (r1 != r2)
        throw invalid_horzcat(r1, c1, r2, c2);
}

void error::check_vertcat(Integer r1, Integer c1, Integer r2, Integer c2)
{
    if (c1 != c2)
        throw invalid_vertcat(r1, c1, r2, c2);
};

void error::check_assign_2(Integer r1, Integer c1, Integer r2, Integer c2)
{
    if ((r1 != r2 ) || (c1 != c2))
        throw invalid_assign_2(r1, c1, r2, c2);
}

void error::check_assign_1(Integer s, Integer r2, Integer c2)
{
    if (r2 == 0 || c2 == 0)
    {
        if (s != 0)
            throw invalid_assign_1(s, r2, c2);

        return;
    };

    Integer d, m;
    m = s%r2;
    if (m != 0)
        throw invalid_assign_1(s, r2, c2);

    d = s/r2;
    if (d != c2)
        throw invalid_assign_1(s, r2, c2);
}

void error::check_reshape(Integer r1, Integer c1, Integer r2, Integer c2)
{
    if (Real(r1)*Real(c1) != Real(r2)*Real(c2)) 
        throw invalid_reshape(r1, c1, r2, c2);
}

void error::check_dim(Integer i, Integer d)
{
    if ((i < 1) || (i > d))
        throw invalid_dim(i, d);
}

void error::check_randperm_arg(Integer n)
{
    if (n < 0) 
        throw randperm_arg_neg(n);
}

void error::check_bspdiag_1starg(Integer r, Integer c)
{
    if (r != 1 && c != 1)
    {
        throw bspdiag_1starg_not_vec(r, c);
    };
}

void error::check_bspdiags_2ndarg(Integer r, Integer c)
{
    if (r != 1 && c != 1)
        throw bspdiags_2ndarg_not_vec(r, c);
}

void error::check_diag_arg(Integer r, Integer c)
{
    if ((r != 1) && (c != 1))
        throw diag_arg_not_vec(r, c);
}

void error::check_linear_index(Integer r, Integer c, Integer rows)
{
    if ( icast_t(Real(r) + Real(c-1)*Real(rows)) == false)
        throw linear_index_too_large(r,c,rows);
}

void error::check_scalar(Integer r,Integer c)
{
    if (r != 1 || c != 1)
        throw scalar_required(r,c);
};

void error::check_row_indices_sortcols(Integer size, Integer c)
{
    if (size > c || size < 0)
        throw error::invalid_row_indices_sortcols(size,c);
};

void error::check_row_indices_elem_sortcols(Integer elem, Integer c)
{
    if (elem > c || elem < -c || elem == 0)
        throw error::invalid_row_indices_elem_sortcols(elem,c);
};

void error::check_col_indices_sortrows(Integer size, Integer c)
{
    if (size > c || size < 0)
        throw error::invalid_cols_indices_sortrows(size,c);
};

void error::check_col_indices_elem_sortrows(Integer elem, Integer c)
{
    if (elem > c || elem < -c || elem == 0)
        throw error::invalid_col_indices_elem_sortrows(elem,c);
};

bool error::check_permutation_vector(const matcl::Matrix& p, Integer length)
{
    if (p.rows() != 1 && p.cols() != 1)
        throw invalid_permvec();

    if (p.length() != length)
        throw invalid_permvec_length(p.length(), length);

    if (p.get_value_code() != matcl::value_code::v_integer)
    {
        if (p.length() != length)
            throw invalid_permvec_length(p.length(), length);
    }

    const Integer* ptr = p.get_array<Integer>();
    bool trivial = true;

    Integer_64 sum = 0;

    for(Integer i = 1; i <= length; ++i)
    {
        Integer val = ptr[i-1];
        if (val != i)
            trivial = false;

        if (val < 1 || val > length)
            throw invalid_permvec();

        sum += val;
    };

    Integer_64 length2 = length;

    if (sum != length2 * (length2 + 1) / 2)
        throw invalid_permvec();

    return trivial;
};

void error::throw_error_single_index(Integer i, Integer size)
{
    throw invalid_single_index(i,size);
}

void error::throw_error_alloc(Integer n)
{
    throw alloc(n); 
}

void error::throw_error_size(Integer r, Integer c)
{
    throw invalid_size(r, c);
}

void error::throw_error_row(Integer i, Integer r, Integer c)
{
    throw invalid_row(i, r, c);
}

void error::throw_error_col(Integer i, Integer r, Integer c)
{
    throw invalid_col(i, r, c);
}

void error::check_eeop(Integer r1, Integer c1, Integer r2, Integer c2)
{
    if (r1 != r2 || c1 != c2)
        throw invalid_eeop(r1, c1, r2, c2);
}

void error::check_mul(Integer r1, Integer c1, Integer r2, Integer c2, trans_type t1, trans_type t2)
{
    if (t1 == trans_type::no_trans)
    {
        if (t2 == trans_type::no_trans)
        {
            if (c1 != r2)
                throw invalid_mul(r1, c1, r2, c2, t1, t2);
        }
        else
        {
            if (c1 != c2)
                throw invalid_mul(r1, c1, r2, c2, t1, t2);
        };
    }
    else
    {
        if (t2 == trans_type::no_trans)
        {
            if (r1 != r2)
                throw invalid_mul(r1, c1, r2, c2, t1, t2);
        }
        else
        {
            if (r1 != c2)
                throw invalid_mul(r1, c1, r2, c2, t1, t2);
        };
    };
}

};};
