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

#include "matcl-linalg/graph/matcl_graph.h"
#include "matcl-core/utils/workspace.h"
#include "dmperm_pothen.h"
#include "maxmatch_weighted.h"

namespace matcl { namespace details
{

template<class V, class S>
struct maxmatch_impl
{
    using Mat   = raw::Matrix<V,S>;

    static void eval(const Mat& A, Matrix& row_set, Matrix& col_set)
    {
        using Mat_S = raw::Matrix<V,struct_sparse>;
        Mat_S As    = raw::converter<Mat_S, Mat>::eval(A);
        return maxmatch_impl<V,struct_sparse>::eval(As, row_set, col_set);
    };
};

template<class V>
struct maxmatch_impl<V,struct_sparse>
{
    using S     = struct_sparse;
    using Mat   = raw::Matrix<V,S>;

    static void eval(const Mat& A, Matrix& row_set, Matrix& col_set)
    {
        Integer nrows           = A.rows();
        Integer ncols           = A.cols();

        mrd::sparse_ccs<V> rep  = A.rep();
        const Integer* colstr   = rep.ptr_c();
        const Integer* rowind   = rep.ptr_r();

        row_set                 = matcl::izeros(nrows,1);
        col_set                 = matcl::izeros(ncols,1);

        if (A.nnz() == 0)
            return;

        Integer*  rowset        = row_set.get_array_unique<Integer>();
        Integer*  colset        = col_set.get_array_unique<Integer>();

        using iworkspace        = matcl::pod_workspace<Integer>;
        iworkspace work         = iworkspace(4*ncols+nrows);
        Integer* work_ptr       = work.ptr();

        Integer* prevcl         = work_ptr;
        Integer* prevrw         = prevcl + ncols;
        Integer* marker         = prevrw + ncols;
        Integer* tryrow         = marker + nrows;
        Integer* nxtchp         = tryrow + ncols;

        int err = maxmatch(nrows, ncols, colstr, rowind, prevcl, prevrw, marker, 
                           tryrow, nxtchp, rowset, colset);

        if (err == 0)
            return;
 
        throw error::error_general("internal error in maxmatch function");
    }
};

template<class V, class S, bool Is_complex = md::is_complex<V>::value>
struct maxmatch_weight_impl
{
    using Mat       = raw::Matrix<V,S>;
    using VR        = typename md::real_type<V>::type;
    using ret_type  = tuple<Matrix,Matrix,Matrix,Matrix, Real, Integer>;

    static void eval(const Mat& A, ret_type& ret)
    {
        using Mat_S = raw::Matrix<VR,struct_sparse>;
        Mat_S As    = raw::converter<Mat_S, Mat>::eval(A);
        return maxmatch_weight_impl<VR,struct_sparse>::eval(As, ret);
    };
};

template<class V>
struct maxmatch_weight_impl<V,struct_sparse, false>
{
    using S         = struct_sparse;
    using Mat       = raw::Matrix<V,S>;
    using Mat_D     = raw::Matrix<V,struct_dense>;
    using VR        = typename md::real_type<V>::type;
    using ret_type  = tuple<Matrix,Matrix,Matrix,Matrix, Real, Integer>;

    static void eval(const Mat& A, ret_type& ret)
    {
        if (A.nnz() == 0)
        {
            value_code vc   = matrix_traits::value_code<V>::value;

            Matrix mr   = izeros(A.rows(), 1);
            Matrix mc   = izeros(A.cols(), 1);
            Matrix dr   = zeros(A.rows(), 1, vc);
            Matrix dc   = zeros(A.cols(), 1, vc);
            Real sum    = 0.0;
            Integer siz = 0;

            ret = ret_type(mr, mc, dr, dc, sum, siz);
            return;
        };

        const Integer* ind_c    = A.rep().ptr_c();
        const Integer* ind_r    = A.rep().ptr_r();
        const V* ptr_x          = A.rep().ptr_x();

        Integer M               = A.rows();
        Integer N               = A.cols();
        Matrix row_set          = matcl::izeros(M,1);
        Matrix col_set          = matcl::izeros(N,1);

        Mat_D dual_rows(ti::ti_empty(), M, 1);
        Mat_D dual_cols(ti::ti_empty(), N, 1);

        Integer* ptr_row        = row_set.get_array_unique<Integer>();
        Integer* ptr_col        = col_set.get_array_unique<Integer>();

        Integer iwork_size      = 3 * N + 2 * M;
        Integer rwork_size      = std::max(M, N);

        using iworkspace        = matcl::pod_workspace<Integer>;
        iworkspace iwork        = iworkspace(iwork_size);
        Integer* iwork_ptr      = iwork.ptr();

        using workspace         = matcl::pod_workspace<VR>;
        workspace rwork         = workspace(rwork_size);
        VR* rwork_ptr           = rwork.ptr();

        V* ptr_du               = dual_rows.ptr();
        V* ptr_dv               = dual_cols.ptr();
        Integer num_match;

        hungarian_match(false, M, N, ind_c, ind_r, ptr_x, ptr_row, ptr_col, num_match, ptr_du, 
                        ptr_dv, iwork_ptr, rwork_ptr);

        Real sum    = 0.0;

        for (Integer i = 0; i < M; ++i)
        {
            V val       = ptr_du[i];
            sum         += val;
        };

        for (Integer i = 0; i < N; ++i)
        {
            V val       = ptr_dv[i];
            sum         += val;
        };

        ret = ret_type(row_set, col_set, Matrix(dual_rows,false), Matrix(dual_cols,false), sum, num_match);
    }
};

struct maxmatch_vis : public extract_type_switch<void, maxmatch_vis,true>
{
    template<class T>
    static void eval(const Matrix&, const T& mat, Matrix& row_set, Matrix& col_set)
    {
        using V = typename T::value_type;
        using S = typename T::struct_type;
        return maxmatch_impl<V,S>::eval(mat, row_set, col_set);
    };

    template<class T>
    static void eval_scalar(const Matrix& handle, const T& v, Matrix& row_set, Matrix& col_set)
    {
        using Mat = raw::Matrix<T,struct_dense>;
        Mat m(ti::get_ti(v), v, 1, 1);
        return eval<Mat>(handle, m, row_set, col_set);
    };
};

struct maxmatch_weight_vis : public extract_type_switch<void, maxmatch_weight_vis,true>
{
    using ret_type = tuple<Matrix,Matrix,Matrix,Matrix, Real, Integer>;

    template<class T>
    static void eval(const Matrix&, const T& mat, ret_type& ret)
    {
        using V = typename T::value_type;
        using S = typename T::struct_type;
        return maxmatch_weight_impl<V,S>::eval(mat, ret);
    };

    template<class S>
    static void eval(const Matrix&, const raw::Matrix<Object,S>&, ret_type&)
    {
        throw error::object_value_type_not_allowed("maxmatch_weight");
    };

    static void eval_scalar(const Matrix&, const Object&, ret_type&)
    {
        throw error::object_value_type_not_allowed("maxmatch_weight");
    };

    template<class T>
    static void eval_scalar(const Matrix& handle, const T& v, ret_type& ret)
    {
        using Mat = raw::Matrix<T,struct_dense>;
        Mat m(ti::get_ti(v), v, 1, 1);
        return eval<Mat>(handle, m, ret);
    };
};

//TODO
template<class V, class S>
struct separator_perm_impl
{
    using Mat       = raw::Matrix<V,S>;
    using VR        = typename md::real_type<V>::type;
    using ret_type  = tuple<permvec, permvec>;

    static void eval(const Mat& A, ret_type& ret)
    {
        using Mat_S = raw::Matrix<VR,struct_sparse>;
        Mat_S As    = raw::converter<Mat_S, Mat>::eval(A);
        return separator_perm_impl<VR,struct_sparse>::eval(As, ret);
    };
};

template<class V>
struct separator_perm_impl<V,struct_sparse>
{
    using S         = struct_sparse;
    using Mat       = raw::Matrix<V,S>;
    using Mat_D     = raw::Matrix<V,struct_dense>;
    using VR        = typename md::real_type<V>::type;
    using ret_type  = tuple<permvec,permvec>;

    static void eval(const Mat& A, ret_type& ret)
    {
        if (A.nnz() == 0)
        {
            permvec p   = permvec::identity(A.rows());
            permvec q   = permvec::identity(A.cols());

            ret = ret_type(p, q);
            return;
        };

        //TODO
        throw;
    }
};

struct separator_perm_vis : public extract_type_switch<void, separator_perm_vis,true>
{
    using ret_type = tuple<permvec, permvec>;

    template<class T>
    static void eval(const Matrix&, const T& mat, ret_type& ret)
    {
        using V = typename T::value_type;
        using S = typename T::struct_type;
        return separator_perm_impl<V,S>::eval(mat, ret);
    };

    template<class S>
    static void eval(const Matrix&, const raw::Matrix<Object,S>&, ret_type&)
    {
        throw error::object_value_type_not_allowed("separator_perm");
    };

    static void eval_scalar(const Matrix&, const Object&, ret_type&)
    {
        throw error::object_value_type_not_allowed("separator_perm");
    };

    template<class T>
    static void eval_scalar(const Matrix& handle, const T& v, ret_type& ret)
    {
        using Mat = raw::Matrix<T,struct_dense>;
        Mat m(ti::get_ti(v), v, 1, 1);
        return eval<Mat>(handle, m, ret);
    };
};

};};

namespace matcl
{

mat_tup_2 matcl::max_match(const Matrix& mat)
{
    Matrix row_set;
    Matrix col_set;

    details::maxmatch_vis::make<const Matrix&>(mat, row_set, col_set);
    return mat_tup_2(row_set, col_set);
};

mat_tup_2 matcl::max_match_weighted(const Matrix& mat)
{
    tuple<Matrix,Matrix,Matrix,Matrix, Real, Integer> ret;
    details::maxmatch_weight_vis::make<const Matrix&>(mat, ret);
    return mat_tup_2(ret.get<1>(), ret.get<2>());
};

tuple<Matrix,Matrix,Matrix,Matrix, Real, Integer>
matcl::max_match_weighted2(const matcl::Matrix& mat)
{
    tuple<Matrix,Matrix,Matrix,Matrix, Real, Integer> ret;
    details::maxmatch_weight_vis::make<const Matrix&>(mat, ret);
    return ret;
};

Integer matcl::sprank(const Matrix& A)
{
    Matrix row_set;
    Matrix col_set;

    tie(row_set, col_set) = max_match(A);
    return nnz_vec(row_set);
};

tuple<permvec, permvec> matcl::colsmatch_to_perms(const Matrix& c, Integer n_rows)
{
    permvec pr, pc;
    tie(pr, pc) = rowsmatch_to_perms(c, n_rows);
    return tuple<permvec, permvec>(pc, pr);
};
tuple<permvec, permvec> matcl::rowsmatch_to_perms(const Matrix& r, Integer n_cols)
{
    if (r.is_vector() == false)
        throw error::vector_required(r.rows(), r.cols());

    Integer M           = r.length();
    Integer N           = n_cols;
    const Integer* ptr  = r.get_array<Integer>();

    using Mat_I = raw::Matrix<Integer, struct_dense>;

    Mat_I pr(ti::ti_empty(), M, 1);
    Mat_I pc(ti::ti_empty(), N, 1);

    Integer* perm_r     = pr.ptr();
    Integer* perm_c     = pc.ptr();

    Integer pos_r       = 0;
    Integer pos_c       = 0;

    using iworkspace    = matcl::pod_workspace<Integer>;
    iworkspace work     = iworkspace(N);
    Integer* work_c     = work.ptr();

    for (Integer i = 0; i < N; ++i)
        work_c[i]       = 0;

    for (Integer i = 0; i < M; ++i)
    {
        Integer val     = ptr[i];

        if (val > 0)
        {
            if (val > N)
                throw error::invalid_single_index(val, N);

            perm_r[pos_r++] = i+1;
            perm_c[pos_c++] = val;
            work_c[val-1]   = 1;
        };
    };

    //unmatched rows
    for (Integer i = 0; i < M; ++i)
    {
        Integer val     = ptr[i];

        if (val <= 0)
            perm_r[pos_r++] = i+1;
    };
    //unmatched cols
    for (Integer i = 0; i < N; ++i)
    {
        Integer val     = work_c[i];

        if (val == 0)
            perm_c[pos_c++] = i+1;
    };

    permvec ret_r   = permvec::from_matrix(Matrix(pr,false));
    permvec ret_c   = permvec::from_matrix(Matrix(pc,false));

    return tuple<permvec, permvec>(ret_r, ret_c);
};

//TODO
tuple<permvec, permvec> matcl::separator_perm(const Matrix& mat)
{
    tuple<permvec, permvec> ret;
    details::separator_perm_vis::make<const Matrix&>(mat, ret);
    return ret;
};

};