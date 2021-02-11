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

#include "matcl-matrep/algs/sparse_algs_utils.h"
#include "matcl-internals/container/mat_d.h"

namespace matcl { namespace algorithm { namespace details
{

sort_type is_sorted(const raw::integer_dense& vec)
{
    if (vec.size() == 0)
        return sorted_increasing;

    const Integer* ptr_vec = vec.make_explicit().ptr();

    Integer val_prev = ptr_vec[0];	
    bool b_increase = true;
    bool b_decrease = true;

    Integer vec_size = vec.size();
    for (Integer i = 1; i < vec_size; ++i)
    {
        Integer row = ptr_vec[i];
        if (row < val_prev)
        {
            b_increase = false;
            break;
        };
        val_prev = row;
    };

    if (b_increase)
        return sorted_increasing;

    val_prev = ptr_vec[0];	
    for (Integer i = 1; i < vec_size; ++i)
    {
        Integer row = ptr_vec[i];
        if (row > val_prev)
        {
            b_decrease = false;
            break;
        };
        val_prev = row;
    };

    if (b_decrease)
        return sorted_decreasing;

    return not_sorted;
};

Integer number_dupl(const raw::integer_dense& vec)
{
    if (vec.size() == 0)
        return 0;

    const Integer* ptr_vec = vec.ptr();
    Integer r           = vec.rows();
    Integer c           = vec.cols();
    Integer vec_ld      = vec.ld();
    Integer val_prev    = ptr_vec[0];	
    Integer out         = 0;

    Integer j = 0, i = 1;
    while (j < c)
    {
        while (i < r)
        {
            Integer row = ptr_vec[i];
            if (row == val_prev)
                ++out;

            val_prev = row;
            ++i;
        };

        ptr_vec += vec_ld;
        i = 0;
        ++j;
    };
    return out;
};

Integer number_dupl(const raw::integer_dense& vec_r, const raw::integer_dense& vec_c)
{
    if (vec_r.size() == 0)
        return 0;

    const Integer* ptr_r= vec_r.ptr();
    const Integer* ptr_c= vec_c.ptr();
    Integer r           = vec_r.rows();
    Integer c           = vec_r.cols();
    Integer vec_r_ld    = vec_r.ld();
    Integer vec_c_ld    = vec_c.ld();
    Integer val_r_prev  = ptr_r[0];	
    Integer val_c_prev  = ptr_c[0];	
    Integer out         = 0;

    Integer j = 0, i = 1;

    while (j < c)
    {
        while (i < r)
        {
            Integer row = ptr_r[i];
            Integer col = ptr_c[i];
            
            if (row == val_r_prev && col == val_c_prev)
                ++out;

            val_r_prev = row;
            val_c_prev = row;

            ++i;
        };

        ptr_r   += vec_r_ld;
        ptr_c   += vec_c_ld;
        i       = 0;
        ++j;
    };

    return out;
};

void details::init_row_selector(const raw::integer_dense& ri,Integer val, matcl::pod_workspace<Integer>& v_work_row,
                       Integer& r_min, Integer& r_max)
{
    if (ri.size() == 0)
        return;

    const Integer* ptr_ri = ri.ptr();

    r_min = ptr_ri[0];
    r_max = r_min;
    Integer ri_size = ri.size();
    for (Integer i = 1; i < ri_size; ++i)
    {
        Integer row = ptr_ri[i];
        if (row > r_max)
            r_max = row;

        if (row < r_min)
            r_min = row;
    };

    --r_max;
    --r_min;

    Integer work_size = r_max-r_min+1;
    v_work_row.resize(work_size);
    for (Integer i = 0; i < work_size;++i)
        v_work_row[i] = val;
};

void details::sort_rows_cols(Integer* ptr_r, Integer* ptr_c, Integer size)
{
    (void)ptr_r;
    (void)ptr_c;
    (void)size;

    //TODO: impl
};

void details::sort_rows_cols(Integer* ptr_r, Integer* ptr_c, Integer* indices, Integer size)
{
    (void)ptr_r;
    (void)ptr_c;
    (void)indices;
    (void)size;

    //TODO: impl
};

template<class T>
void details::sort_rows_cols_stable(Integer* ptr_r, Integer* ptr_c, T* ptr_x, Integer size)
{
    (void)ptr_r;
    (void)ptr_c;
    (void)ptr_x;
    (void)size;

    //TODO: impl
};

template
void sort_rows_cols_stable(Integer* ptr_r, Integer* ptr_c, double* ptr_x, Integer size);

template
void sort_rows_cols_stable(Integer* ptr_r, Integer* ptr_c, float* ptr_x, Integer size);

template
void sort_rows_cols_stable(Integer* ptr_r, Integer* ptr_c, Integer* ptr_x, Integer size);

template
void sort_rows_cols_stable(Integer* ptr_r, Integer* ptr_c, Complex* ptr_x, Integer size);

template
void sort_rows_cols_stable(Integer* ptr_r, Integer* ptr_c, Float_complex* ptr_x, Integer size);

template
void sort_rows_cols_stable(Integer* ptr_r, Integer* ptr_c, Object* ptr_x, Integer size);

};};};
