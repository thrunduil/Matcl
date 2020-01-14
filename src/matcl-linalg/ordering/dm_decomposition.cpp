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

#include "matcl-matrep/matcl_matrep.h"
#include "matcl-linalg/graph/dm_decomposition.h"
#include "matcl-core/utils/workspace.h"
#include "dmperm_pothen.h"
#include "matcl-linalg/general/linalg_exception.h"

namespace matcl { namespace details
{

template<class V, class S>
struct dmperm_impl
{};

template<class V>
struct dmperm_impl<V,struct_banded>
{
    using S     = struct_banded;
    using Mat   = raw::Matrix<V,S>;

    static void eval(const Mat& A, dm_decomp& ret)
    {
        using Mat_S = raw::Matrix<V,struct_sparse>;
        Mat_S As    = raw::converter<Mat_S, Mat>::eval(A);
        return dmperm_impl<V,struct_sparse>::eval(As, ret);
    };
};

template<class V>
struct dmperm_impl<V,struct_dense>
{
    using S     = struct_dense;
    using Mat   = raw::Matrix<V,S>;

    static void eval(const Mat& A, dm_decomp& ret)
    {
        using Mat_S = raw::Matrix<V,struct_sparse>;
        Mat_S As    = raw::converter<Mat_S, Mat>::eval(A);
        return dmperm_impl<V,struct_sparse>::eval(As, ret);
    };
};

template<class V>
struct dmperm_impl<V,struct_sparse>
{
    using S     = struct_sparse;
    using Mat   = raw::Matrix<V,S>;
    using ccs   = mrd::sparse_ccs<V>;
    using Mat_I = raw::Matrix<Integer,struct_dense>;

    static void fast_exit(const Mat& A, dm_decomp& ret)
    {
        Integer nrows           = A.rows();
        Integer ncols           = A.cols();

        permvec pr  = permvec::identity(nrows);
        permvec pc  = permvec::identity(ncols);

        if (nrows == 0 && ncols == 0)
        {
            Matrix ind_comp_r   = mat_col().add(nrows+1);
            Matrix ind_comp_c   = mat_col().add(ncols+1);
            ret.m_comp_under    = 0;
            ret.m_comp_square   = 0;
            ret.m_comp_over     = 0;
            ret.m_split_rows    = ind_comp_r;
            ret.m_split_cols    = ind_comp_c;
        }
        else if (nrows == 0)
        {
            Matrix ind_comp_r   = mat_col().add(nrows).add(nrows+1);
            Matrix ind_comp_c   = mat_col().add(1).add(ncols+1);
            ret.m_comp_under    = 0;
            ret.m_comp_square   = 0;
            ret.m_comp_over     = 1;
            ret.m_split_rows    = ind_comp_r;
            ret.m_split_cols    = ind_comp_c;
        }
        else if (ncols == 0)
        {
            Matrix ind_comp_r   = mat_col().add(1).add(nrows+1);
            Matrix ind_comp_c   = mat_col().add(ncols).add(ncols+1);
            ret.m_comp_under    = 1;
            ret.m_comp_square   = 0;
            ret.m_comp_over     = 0;
            ret.m_split_rows    = ind_comp_r;
            ret.m_split_cols    = ind_comp_c;
        }
        else
        {
            Matrix ind_comp_r   = mat_col().add(1).add(1).add(nrows+1);
            Matrix ind_comp_c   = mat_col().add(1).add(ncols).add(ncols+1);
            ret.m_comp_under    = 1;
            ret.m_comp_square   = 0;
            ret.m_comp_over     = 1;
            ret.m_split_rows    = ind_comp_r;
            ret.m_split_cols    = ind_comp_c;
        };
        ret.m_rows_perm     = pr;
        ret.m_cols_perm     = pc;
    };

    static void eval(const Mat& A, dm_decomp& ret)
    {
        if (A.nnz() == 0)
        {
            fast_exit(A, ret);
            return;
        };

        Matrix At(A, false);
        At          = trans(At);
        At          = convert(At, Mat::matrix_code);
        Mat raw_At  = At.get_impl_unique<Mat>();

        if (A.rows() >= A.cols())
            eval_impl(A, raw_At, false, ret);
        else
            eval_impl(raw_At, A, true, ret);
    };

    static void eval_impl(const Mat& A, const Mat& At, bool trans, dm_decomp& ret)
    {
        Integer nrows           = A.rows();
        Integer ncols           = A.cols();

        ccs rep                 = A.rep();
        const Integer* colstr   = rep.ptr_c();
        const Integer* rowidx   = rep.ptr_r();
        
        ccs rep_At              = At.rep();
        const Integer* rowstr   = rep_At.ptr_c();
        const Integer* colidx   = rep_At.ptr_r();

        Integer size_work       = 5*nrows + 5*ncols;
        Integer size_work_rpart = nrows + ncols + 1;
        Integer size_work_cpart = nrows + ncols + 1;

        using iworkspace        = matcl::pod_workspace<Integer>;
        iworkspace work         = iworkspace(size_work + size_work_rpart + size_work_cpart);
        Integer* work_ptr       = work.ptr();

        Mat_I perm_rows(ti::ti_empty(), nrows, 1);
        Mat_I perm_cols(ti::ti_empty(), ncols, 1);

        Integer* perm_rows_ptr  = perm_rows.ptr();
        Integer* perm_cols_ptr  = perm_cols.ptr();
        Integer* wpart_rows_ptr = work_ptr + size_work;
        Integer* wpart_cols_ptr = wpart_rows_ptr + size_work_rpart;

        // horizontal (underdetermined) block
        Integer horiz_rows      = 0;
        Integer horiz_cols      = 0;
        Integer horiz_conn_comp = 0;

        //square (exactly determined) block
        Integer square_size     = 0;
        Integer square_conn_comp= 0;

        //vertical (overdetermined) block
        Integer vert_rows       = 0;
        Integer vert_cols       = 0;
        Integer vert_conn_comp  = 0;

        int err = genbtf(nrows, ncols, colstr, rowidx, rowstr, colidx, work_ptr, perm_rows_ptr, 
                         perm_cols_ptr, horiz_rows, horiz_cols, horiz_conn_comp, square_size, 
                         square_conn_comp, vert_rows, vert_cols, vert_conn_comp, wpart_rows_ptr, 
                         wpart_cols_ptr);

        if (err != 0)
            throw error::error_general("internal error in maxmatch function");

        make_upper_triangular(perm_rows_ptr, perm_cols_ptr, wpart_rows_ptr, wpart_cols_ptr,
                              horiz_conn_comp, square_conn_comp, work_ptr);
        collect_zero_blocks(wpart_rows_ptr, wpart_cols_ptr, horiz_conn_comp, square_conn_comp, 
                            vert_conn_comp);

        permvec pr  = permvec::from_matrix(Matrix(perm_rows, false));
        permvec pc  = permvec::from_matrix(Matrix(perm_cols, false));

        Integer num_comp    = horiz_conn_comp + square_conn_comp + vert_conn_comp;

        Mat_I ind_part_r(ti::ti_empty(), num_comp + 1, 1);
        Mat_I ind_part_c(ti::ti_empty(), num_comp + 1, 1);

        Integer* ind_r  = ind_part_r.ptr();
        Integer* ind_c  = ind_part_c.ptr();

        for (Integer i = 0; i <= num_comp; ++ i)
        {
            ind_r[i]    = wpart_rows_ptr[i];
            ind_c[i]    = wpart_cols_ptr[i];
        };

        Matrix ind_comp_r(ind_part_r, false);
        Matrix ind_comp_c(ind_part_c, false);

        if (trans == false)
        {
            ret.m_rows_perm     = pr;
            ret.m_cols_perm     = pc;
            ret.m_comp_under    = horiz_conn_comp;
            ret.m_comp_square   = square_conn_comp;
            ret.m_comp_over     = vert_conn_comp;
            ret.m_split_rows    = ind_comp_r;
            ret.m_split_cols    = ind_comp_c;
        }
        else
        {
            ret.m_cols_perm     = pr;
            ret.m_rows_perm     = pc;
            ret.m_comp_over     = horiz_conn_comp;
            ret.m_comp_square   = square_conn_comp;
            ret.m_comp_under    = vert_conn_comp;
            ret.m_split_cols    = ind_comp_r;
            ret.m_split_rows    = ind_comp_c;
        }
    };

    static void make_upper_triangular(Integer* perm_rows_ptr, Integer* perm_cols_ptr, 
                    Integer* part_rows_ptr, Integer* part_cols_ptr, Integer horiz_conn_comp, 
                    Integer square_conn_comp, Integer* work)
    {
        Integer square_block_fr = part_rows_ptr[horiz_conn_comp];
        Integer square_block_lr = part_rows_ptr[horiz_conn_comp + square_conn_comp] - 1;

        Integer square_block_fc = part_cols_ptr[horiz_conn_comp];
        Integer square_block_lc = part_cols_ptr[horiz_conn_comp + square_conn_comp] - 1;

        // revert rows and cols in square block
        Integer pos_f           = square_block_fr - 1;
        Integer pos_l           = square_block_lr - 1;
        while(pos_f < pos_l)
        {
            std::swap(perm_rows_ptr[pos_f], perm_rows_ptr[pos_l]);
            ++pos_f;
            --pos_l;
        };

        pos_f                   = square_block_fc - 1;
        pos_l                   = square_block_lc - 1;
        while(pos_f < pos_l)
        {
            std::swap(perm_cols_ptr[pos_f], perm_cols_ptr[pos_l]);
            ++pos_f;
            --pos_l;
        };

        //revert component indices
        pos_f                   = horiz_conn_comp;
        Integer pos             = square_conn_comp - 1;
        for (Integer i = pos_f; pos >= 0; ++i, --pos)
            work[pos]           = part_rows_ptr[i+1] - part_rows_ptr[i];

        pos                     = horiz_conn_comp + 1;
        for (Integer i = 0; i < square_conn_comp; ++i, ++pos)
            part_rows_ptr[pos]  = part_rows_ptr[pos-1] + work[i];

        pos                     = horiz_conn_comp + 1;
        for (Integer i = 0; i < square_conn_comp; ++i, ++pos)
            part_cols_ptr[pos]  = part_cols_ptr[pos-1] + work[i];
    };

    static void collect_zero_blocks(Integer* part_rows_ptr, Integer* part_cols_ptr, 
                    Integer& horiz_conn_comp, Integer square_conn_comp, Integer& vert_conn_comp)
    {
        Integer pos         = 1;
        Integer i           = 1;

        Integer horiz_end   = horiz_conn_comp;
        Integer square_end  = horiz_conn_comp + square_conn_comp;
        Integer vert_end    = square_end + vert_conn_comp;
        Integer last_elem_r = part_rows_ptr[vert_end];
        Integer last_elem_c = part_cols_ptr[vert_end];

        //horizontal components
        for (i = 1; i < horiz_end; ++ i)
        {
            if (part_rows_ptr[i] == part_rows_ptr[i-1])
            {                
                --horiz_conn_comp;
            }
            else
            {
                part_rows_ptr[pos]  = part_rows_ptr[i];
                part_cols_ptr[pos]  = part_cols_ptr[i];
                ++pos;
            };
        };

        //square components
        for (; i < square_end; ++i)
        {
            part_rows_ptr[pos]  = part_rows_ptr[i];
            part_cols_ptr[pos]  = part_cols_ptr[i];
            ++pos;
        };

        //last element or first vertical block
        part_rows_ptr[pos]  = part_rows_ptr[i];
        part_cols_ptr[pos]  = part_cols_ptr[i];
        ++pos;
        ++i;

        //vertical components
        for (; i < vert_end; ++i)
        {
            if (part_cols_ptr[i] == part_cols_ptr[i-1])
            {                
                --vert_conn_comp;
            }
            else
            {
                part_rows_ptr[pos]  = part_rows_ptr[i];
                part_cols_ptr[pos]  = part_cols_ptr[i];
                ++pos;
            };
        };

        //last element
        Integer last        = horiz_conn_comp + square_conn_comp + vert_conn_comp;
        part_rows_ptr[last] = last_elem_r;
        part_cols_ptr[last] = last_elem_c;
    };
};

struct dmperm_vis : public extract_type_switch<void, dmperm_vis,true>
{
    template<class T>
    static void eval(const Matrix&, const T& mat, dm_decomp& ret)
    {
        using V = typename T::value_type;
        using S = typename T::struct_type;
        return dmperm_impl<V,S>::eval(mat, ret);
    };

    template<class T>
    static void eval_scalar(const Matrix& handle, const T& v, dm_decomp& ret)
    {
        using Mat = raw::Matrix<T,struct_dense>;
        Mat m(ti::get_ti(v), v, 1, 1);
        return eval<Mat>(handle, m, ret);
    };
};

}};

namespace matcl
{

dm_decomp::dm_decomp()
{
    m_rows_perm     = permvec::identity(1);
    m_cols_perm     = permvec::identity(1);
    m_comp_under    = 0;
    m_comp_square   = 1;
    m_comp_over     = 0;
    m_split_rows    = mat_col().add(1).add(2);
    m_split_cols    = mat_col().add(1).add(2);
};

// performe Dulmage-Mendelsohn decomposition of a matrix
dm_decomp::dm_decomp(const Matrix& A)
{
    details::dmperm_vis::make<const Matrix&>(A, *this);
};

dm_decomp::~dm_decomp()
{};

permvec dm_decomp::row_perms() const
{
    return m_rows_perm;
};

permvec dm_decomp::col_perms() const
{
    return m_cols_perm;
};

Integer dm_decomp::components_under() const
{
    return m_comp_under;
};

Integer dm_decomp::components_square() const
{
    return m_comp_square;
};

Integer dm_decomp::components_over() const
{
    return m_comp_over;
};

Integer dm_decomp::components_total() const
{
    return m_comp_under + m_comp_square + m_comp_over;
};

Matrix dm_decomp::row_split() const
{
    return m_split_rows;
};

Matrix dm_decomp::col_split() const
{
    return m_split_cols;
};

bool dm_decomp::has_under() const
{
    return m_comp_under > 0;
};

bool dm_decomp::has_square() const
{
    return m_comp_square > 0;
};

bool dm_decomp::has_over() const
{
    return m_comp_over > 0;
};

Integer dm_decomp::rows() const
{
    return m_split_rows(end).get_scalar<Integer>() - 1;
};

Integer dm_decomp::cols() const
{
    return m_split_cols(end).get_scalar<Integer>() - 1;
};

Integer dm_decomp::under_rows() const
{
    const Integer* arr  = m_split_rows.get_array<Integer>();
    return arr[m_comp_under] - arr[0];
};
Integer dm_decomp::under_cols() const
{
    const Integer* arr  = m_split_cols.get_array<Integer>();
    return arr[m_comp_under] - arr[0];
};
Integer dm_decomp::square_size() const
{
    const Integer* arr  = m_split_cols.get_array<Integer>();
    return arr[m_comp_under + m_comp_square] - arr[m_comp_under];
};
Integer dm_decomp::over_rows() const
{
    const Integer* arr  = m_split_rows.get_array<Integer>();
    return arr[m_comp_under + m_comp_square + m_comp_over] - arr[m_comp_under + m_comp_square];
};
Integer dm_decomp::over_cols() const
{
    const Integer* arr  = m_split_cols.get_array<Integer>();
    return arr[m_comp_under + m_comp_square + m_comp_over] - arr[m_comp_under + m_comp_square];
};

dm_decomp::int_tup_4 dm_decomp::get_component_under(Integer comp) const
{
    if (comp < 1 || comp > components_under())
        throw error::invalid_dm_block(comp, components_under(), 1);

    const Integer* ar   = m_split_rows.get_array<Integer>();
    const Integer* ac   = m_split_cols.get_array<Integer>();

    Integer first_row   = ar[comp - 1];
    Integer last_row    = ar[comp] - 1;
    Integer first_col   = ac[comp - 1];
    Integer last_col    = ac[comp] - 1;

    return int_tup_4(first_row, last_row, first_col, last_col);
};
dm_decomp::int_tup_4 dm_decomp::get_component_square(Integer comp) const
{
    if (comp < 1 || comp > components_under())
        throw error::invalid_dm_block(comp, components_square(), 2);

    const Integer* ar   = m_split_rows.get_array<Integer>();
    const Integer* ac   = m_split_cols.get_array<Integer>();

    Integer first_row   = ar[m_comp_under + comp - 1];
    Integer last_row    = ar[m_comp_under + comp] - 1;
    Integer first_col   = ac[m_comp_under + comp - 1];
    Integer last_col    = ac[m_comp_under + comp] - 1;

    return int_tup_4(first_row, last_row, first_col, last_col);
};

dm_decomp::int_tup_4 dm_decomp::get_component_over(Integer comp) const
{
    if (comp < 1 || comp > components_over())
        throw error::invalid_dm_block(comp, components_over(), 3);

    const Integer* ar   = m_split_rows.get_array<Integer>();
    const Integer* ac   = m_split_cols.get_array<Integer>();

    Integer first_row   = ar[m_comp_under + m_comp_square + comp - 1];
    Integer last_row    = ar[m_comp_under + m_comp_square + comp] - 1;
    Integer first_col   = ac[m_comp_under + m_comp_square + comp - 1];
    Integer last_col    = ac[m_comp_under + m_comp_square + comp] - 1;

    return int_tup_4(first_row, last_row, first_col, last_col);
};

dm_decomp::int_tup_4 dm_decomp::get_under_block() const
{
    const Integer* ar   = m_split_rows.get_array<Integer>();
    const Integer* ac   = m_split_cols.get_array<Integer>();

    Integer first_row   = 1;
    Integer last_row    = ar[m_comp_under] - 1;
    Integer first_col   = 1;
    Integer last_col    = ac[m_comp_under] - 1;

    return int_tup_4(first_row, last_row, first_col, last_col);
};

dm_decomp::int_tup_4 dm_decomp::get_square_block() const
{
    const Integer* ar   = m_split_rows.get_array<Integer>();
    const Integer* ac   = m_split_cols.get_array<Integer>();

    Integer first_row   = ar[m_comp_under];
    Integer last_row    = ar[m_comp_under + m_comp_square] - 1;
    Integer first_col   = ac[m_comp_under];
    Integer last_col    = ac[m_comp_under + m_comp_square] - 1;

    return int_tup_4(first_row, last_row, first_col, last_col);
};

dm_decomp::int_tup_4 dm_decomp::get_over_block() const
{
    const Integer* ar   = m_split_rows.get_array<Integer>();
    const Integer* ac   = m_split_cols.get_array<Integer>();

    Integer first_row   = ar[m_comp_under + m_comp_square];
    Integer last_row    = ar[m_comp_under + m_comp_square + m_comp_over] - 1;
    Integer first_col   = ac[m_comp_under + m_comp_square];
    Integer last_col    = ac[m_comp_under + m_comp_square + m_comp_over] - 1;

    return int_tup_4(first_row, last_row, first_col, last_col);
};

Integer dm_decomp::sprank() const
{
    return under_rows() + square_size() + over_cols();
};

};