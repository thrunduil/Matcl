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

#include "matcl-linalg/graph/graph_manip.h"
//#include "matcl-linalg/petsc_utils/petsc_algs.h"
#include "matcl-linalg/general/linalg_exception.h"
#include "matcl-linalg/decompositions/eig_functions.h"

namespace matcl { namespace details
{

template<class V>
struct make_laplacian_impl
{
    using Mat = raw::Matrix<V,struct_sparse>;
    static void eval(const Mat& mat, Matrix& L, bool use_weights)
    {
        if (mat.nnz() == 0)
        {
            L           = Matrix(mat,false);
            return;
        };

        if (use_weights == false)
            return eval_noweights(mat, L);
        
        using VR        = typename md::unify_types<V,Float>::type;
        using Mat_ret   = raw::Matrix<VR, struct_sparse>;

        Integer N       = mat.rows();        

        using ccs       = mrd::sparse_ccs<V>;
        using ccs_ret   = mrd::sparse_ccs<VR>;

        Mat_ret L_rep(mat.get_type(), N, N, mat.nnz() + N);

        const ccs& md   = mat.rep();
        ccs_ret ret_d   = L_rep.rep();

        const Integer* A_r  = md.ptr_r();
        const Integer* A_c  = md.ptr_c();
        const V* A_x        = md.ptr_x();

        Integer* L_r        = ret_d.ptr_r();
        Integer* L_c        = ret_d.ptr_c();
        VR* L_x             = ret_d.ptr_x();

        Integer nz          = 0;

        for (Integer i = 0; i < N; ++i)
        {
            L_c[i]          = nz;

            Integer kf      = A_c[i];
            Integer kl      = A_c[i+1];
            bool has_diag   = false;
            V deg           = V(0);
            Integer pos_diag= 0;

            Integer k       = kf;
            for (; k < kl; ++k)
            {
                Integer r   = A_r[k];
                if (r < i)
                {
                    V val   = A_x[k];
                    L_r[nz] = r;
                    L_x[nz] = -val;
                    deg     = deg + val;
                    ++nz;
                }
                else if (r == i)
                {
                    L_r[nz] = r;
                    has_diag= true;
                    pos_diag= nz;
                    ++nz;
                }
                else
                {
                    break;
                }
            };

            if (has_diag == false)
            {
                L_r[nz]     = i;
                pos_diag    = nz;
                ++nz;
            }
            for (; k < kl; ++k)
            {
                V val       = A_x[k];
                Integer r   = A_r[k];
                L_r[nz]     = r;
                L_x[nz]     = -val;
                deg         = deg + val;
                ++nz;
            };

            L_x[pos_diag]   = deg;
        };

        L_c[N]              = nz;
        L                   = Matrix(L_rep,false);
        return;
    };

    static void eval_noweights(const Mat& mat, Matrix& L)
    {
        using VR0       = typename md::real_type<V>::type;
        using VR        = typename md::unify_types<VR0,Float>::type;
        using Mat_ret   = raw::Matrix<VR, struct_sparse>;

        Integer N   = mat.rows();

        Mat_ret L_rep(mat.get_type(), N, N, mat.nnz() + N);

        using ccs       = mrd::sparse_ccs<V>;
        using ccs_ret   = mrd::sparse_ccs<VR>;

        const ccs& md   = mat.rep();
        ccs_ret ret_d   = L_rep.rep();

        const Integer* A_r  = md.ptr_r();
        const Integer* A_c  = md.ptr_c();

        Integer* L_r        = ret_d.ptr_r();
        Integer* L_c        = ret_d.ptr_c();
        VR* L_x             = ret_d.ptr_x();

        Integer nz          = 0;

        for (Integer i = 0; i < N; ++i)
        {
            L_c[i]          = nz;

            Integer kf      = A_c[i];
            Integer kl      = A_c[i+1];
            bool has_diag   = false;
            Integer deg     = kl - kf;

            Integer k       = kf;
            for (; k < kl; ++k)
            {
                Integer r   = A_r[k];
                if (r < i)
                {
                    L_r[nz] = r;
                    L_x[nz] = VR(-1.0);
                    ++nz;
                }
                else if (r == i)
                {
                    L_r[nz] = r;
                    L_x[nz] = VR(deg - 1);
                    has_diag= true;
                    ++nz;
                }
                else
                {
                    break;
                }
            };

            if (has_diag == false)
            {
                L_r[nz]     = i;
                L_x[nz]     = VR(deg);
                ++nz;
            }
            for (; k < kl; ++k)
            {
                Integer r   = A_r[k];
                L_r[nz]     = r;
                L_x[nz]     = VR(-1.0);
                ++nz;
            };
        };

        L_c[N]              = nz;
        L                   = Matrix(L_rep,false);
        return;
    };
};

template<class V>
struct make_adj_impl
{
    using Mat = raw::Matrix<V,struct_sparse>;

    static void eval(const Mat& mat, Matrix& ret)
    {
        if (mat.nnz() == 0)
        {
            ret         = Matrix(mat, false);
            return;
        };

        Integer N       = mat.rows();

        Mat ret_rep(mat.get_type(), N, N, mat.nnz());

        using ccs       = mrd::sparse_ccs<V>;

        const ccs& md   = mat.rep();
        ccs ret_d       = ret_rep.rep();

        const Integer* A_r  = md.ptr_r();
        const Integer* A_c  = md.ptr_c();

        Integer* ret_r      = ret_d.ptr_r();
        Integer* ret_c      = ret_d.ptr_c();        

        //symbolic structure
        for (Integer j = 0; j <= N; ++j)
            ret_c[j]        = 0;

        for (Integer j = 0; j < N; ++j)
        {
            for (Integer k = A_c[j]; k < A_c[j+1]; ++k)
            {
                Integer r   = A_r[k];

                if (r <= j)
                    continue;                

                ++ret_c[j];
                ++ret_c[r];
            };
        };

        Integer nz  = 0;
        for (Integer j = 0; j < N; ++j)
        {
            Integer lnz     = ret_c[j];
            ret_c[j]        = nz;
            nz              += lnz;
        };
        ret_c[N]            = nz;

        //numeric phase
        const V* A_x        = md.ptr_x();
        V* ret_x            = ret_d.ptr_x();

        for (Integer j = 0; j < N; ++j)
        {
            for (Integer k = A_c[j]; k < A_c[j+1]; ++k)
            {
                Integer r   = A_r[k];

                if (r <= j)
                    continue;                

                Integer pos = ret_c[j];
                ret_x[pos]  = A_x[k];
                ret_r[pos]  = r;
                ++ret_c[j];

                pos         = ret_c[r];
                ret_x[pos]  = conj(A_x[k]);
                ret_r[pos]  = j;
                ++ret_c[r];
            };
        };        

        for (Integer j = N; j >= 1; --j)
            ret_c[j]    = ret_c[j-1];

        ret_c[0]        = 0;

        ret = Matrix(ret_rep, false);
    };
};

template<class V>
struct make_aggreg_impl
{
    using Mat = raw::Matrix<V,struct_sparse>;

    static void eval(const Mat& A, permvec& ret)
    {
        using Mat_I     = raw::Matrix<Integer, struct_dense>;
        using ccs_rep   = raw::details::sparse_ccs<V>;

        ccs_rep d           = A.rep();
        Integer off         = d.offset();
        const Integer* d_r  = d.ptr_r() + off;
        Integer nz          = d.nnz();

        if (nz != A.rows())
            throw error::invalid_aggregation_matrix(A.rows(), A.nnz());

        Mat_I p             = Mat_I(ti::ti_empty(), nz, 1);
        Integer* ptr_p      = p.ptr();

        for (Integer i = 0; i < nz; ++i)
            ptr_p[i]        = d_r[i] + 1;

        Matrix mat_p        = Matrix(p,false);
        ret = permvec::from_matrix(mat_p);
    };
};

struct visit_sparse_lap : public md::extract_type_switch<void, visit_sparse_lap,true>
{
    using base_type = md::extract_type_switch<void, visit_sparse_lap,true>;

    template<class V>
    static void eval(const matcl::Matrix&, const raw::Matrix<V,struct_sparse>& mat,
                     Matrix& L, bool use_weights)
    {
        return make_laplacian_impl<V>::eval(mat, L, use_weights);
    };

    template<class Mat, class ... Arg>
    static void eval(const matcl::Matrix& h, const Mat&, Arg&& ... args)
    {
        return base_type::make<const Matrix&>(sparse(h), std::forward<Arg>(args)...);
    };

    template<class Mat, class ... Arg>
    static void eval_scalar(const matcl::Matrix& h, const Mat&, Arg&& ... args)
    {
        return base_type::make<const Matrix&>(sparse(h), std::forward<Arg>(args)...);
    };
};

struct visit_sparse_adj : public md::extract_type_switch<void, visit_sparse_adj,true>
{
    using base_type = md::extract_type_switch<void, visit_sparse_adj,true>;

    template<class V>
    static void eval(const matcl::Matrix&, const raw::Matrix<V,struct_sparse>& mat,
                     Matrix& ret)
    {
        return make_adj_impl<V>::eval(mat, ret);
    };

    template<class Mat, class ... Arg>
    static void eval(const matcl::Matrix& h, const Mat&, Arg&& ... args)
    {
        return base_type::make<const Matrix&>(sparse(h), std::forward<Arg>(args)...);
    };

    template<class Mat, class ... Arg>
    static void eval_scalar(const matcl::Matrix& h, const Mat&, Arg&& ... args)
    {
        return base_type::make<const Matrix&>(sparse(h), std::forward<Arg>(args)...);
    };
};

struct visit_aggreg : public md::extract_type_switch<void, visit_aggreg,true>
{
    using base_type = md::extract_type_switch<void, visit_aggreg,true>;

    template<class V>
    static void eval(const matcl::Matrix&, const raw::Matrix<V,struct_sparse>& mat,
                     permvec& ret)
    {
        return make_aggreg_impl<V>::eval(mat, ret);
    };

    template<class Mat, class ... Arg>
    static void eval(const matcl::Matrix& h, const Mat&, Arg&& ... args)
    {
        return base_type::make<const Matrix&>(sparse(h), std::forward<Arg>(args)...);
    };

    template<class Mat, class ... Arg>
    static void eval_scalar(const matcl::Matrix& h, const Mat&, Arg&& ... args)
    {
        return base_type::make<const Matrix&>(sparse(h), std::forward<Arg>(args)...);
    };
};

}};
namespace matcl
{

Matrix matcl::symgraph_sum(const Matrix& A)
{
    Matrix ret = hersum(A);
    ret.add_struct(predefined_struct_type::her);
    return ret;
}

Matrix matcl::symgraph_prod(const Matrix& A, bool trans)
{
    return herprod(A, trans);
}

Matrix matcl::symgraph_bipart(const Matrix& A, bool trans)
{
    Integer M   = A.rows();
    Integer N   = A.cols();

    if (trans == false)
    {
        Matrix ret = mat_col().add(mat_row().add(spzeros(N, M)).add(ctrans(A)))
                              .add(mat_row().add(A).add(spzeros(M, N)));
        ret.add_struct(predefined_struct_type::her);
        return ret;
    }
    else
    {
        Matrix ret = mat_col().add(mat_row().add(spzeros(M, N)).add(A))
                              .add(mat_row().add(ctrans(A)).add(spzeros(N, M)));
        ret.add_struct(predefined_struct_type::her);
        return ret;
    }
};

permvec matcl::aggreg_to_perm(const Matrix& aggreg)
{
    permvec p;
    details::visit_aggreg::make<const Matrix&>(aggreg, p);
    return p;
};

tuple<permvec,permvec> matcl::biperm_to_perms(const permvec& p, Integer M, Integer N, bool trans)
{
    if (p.length() != M + N)
        throw error::invalid_permvec_length(p.length(), M+N);

    using Mat_I             = raw::Matrix<Integer, struct_dense>;

    Mat_I ir(ti::ti_empty(), M, 1);
    Mat_I ic(ti::ti_empty(), N, 1);

    const Integer* ptr_p    = p.to_array();
    Integer* ptr_r          = ir.ptr();
    Integer* ptr_c          = ic.ptr();

    // first N indices: columns, next M indices: rows
    Integer pos_r           = 0;
    Integer pos_c           = 0;

    if (trans == false)
    {
        for (Integer i = 0; i < M + N; ++i)
        {
            Integer val         = ptr_p[i];
            if (val <= N)
                ptr_c[pos_c++]  = val;
            else
                ptr_r[pos_r++]  = val - N;
        };
    }
    else
    {
        for (Integer i = 0; i < M + N; ++i)
        {
            Integer val         = ptr_p[i];
            if (val <= M)
                ptr_r[pos_r++]  = val;
            else
                ptr_c[pos_c++]  = val - M;
        };
    };

    Matrix mr(ir,false);
    Matrix mc(ic,false);

    return tuple<permvec,permvec>(permvec::from_matrix(mr), permvec::from_matrix(mc));
};

Matrix matcl::symgraph_laplacian(const Matrix& A, bool use_weights)
{
    if (A.is_square() == false)
        throw error::square_matrix_required(A.rows(), A.cols());

    Matrix L;
    details::visit_sparse_lap::make<const Matrix&>(A, L, use_weights);
    L.add_struct(predefined_struct_type::her);
    return L;
};

Matrix matcl::make_adjancency_matrix(const Matrix& A)
{
    if (A.is_square() == false)
        throw error::square_matrix_required(A.rows(), A.cols());

    Matrix ret;
    details::visit_sparse_adj::make<const Matrix&>(A, ret);
    ret.add_struct(predefined_struct_type::her);
    return ret;
};

mat_tup_2 matcl::fiedler_vectors(const Matrix& L, Integer N, bool direct)
{
    if (L.is_square() == false)
        throw error::square_matrix_required(L.rows(), L.cols());

    if (direct == true)
    {
        Matrix E, V;
        tie(E, V)   = eigsel_index2(L, 1, std::max(N,1));

        return mat_tup_2(V, E);
    }
    else
    {
        pschur_decomposition pd(L, std::max(N,1), cluster_type::SM);
        Matrix E = pd.eig();
        Matrix V = pd.U();

        return mat_tup_2(V, E);
    };
};

};