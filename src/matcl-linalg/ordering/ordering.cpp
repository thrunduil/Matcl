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

#include "matcl-linalg/graph/matcl_graph.h"
//#include "matcl-linalg/petsc_utils/petsc_algs.h"
#include "matcl-linalg/general/linalg_exception.h"
#include "matcl-core/utils/workspace.h"
#include "matcl-linalg/decompositions/lu/lu_superlu.h"
#include "matcl-linalg/decompositions/balancing.h"

extern "C"
{
    #include "extern/cholmod/colamd.h"

    typedef int shortint;

    int genmmd_(int *neqns, int *xadj, shortint *adjncy, 
	    shortint *invp, shortint *perm, int *delta, shortint *dhead, 
	    shortint *qsize, shortint *llist, shortint *marker, int *maxint, 
	    int *nofsub);
}

namespace matcl { namespace details
{

enum class mmd_type
{
    mmd,
    colamd,
    symamd
};

template<class V>
struct make_mmd_impl
{
    using Mat   = raw::Matrix<V,struct_sparse>;
    using Mat_I = raw::Matrix<Integer,struct_dense>;

    static void eval(const Mat& A, permvec& ret, mmd_type type)
    {
        switch (type)
        {
            case mmd_type::mmd:
                return eval_mmd(A, ret);
            case mmd_type::colamd:
                return eval_colamd(A, ret);
            case mmd_type::symamd:
                return eval_symamd(A, ret);
        }
    };

    static void eval_colamd(const Mat& A, permvec& ret)
    {
        Integer M               = A.rows();
        Integer N               = A.cols();
        Integer nz              = A.nnz();

        if (nz == 0)
        {
            ret = permvec::identity(N);
            return;
        }

        const Integer* ptr_c    = A.rep().ptr_c();
        const Integer* ptr_r    = A.rep().ptr_r();
        Integer off             = A.rep().offset();

        ptr_r                   += off;

        Mat_I perm(ti::ti_empty(), N, 1);

        double knobs[COLAMD_KNOBS];
        int stats[COLAMD_STATS];

        int Alen                = (int)colamd_recommended(nz, M, N);
        colamd_set_defaults(knobs);

        Integer size_work       = Alen + (N + 1);
        using iworkspace        = matcl::pod_workspace<Integer>;
        iworkspace work         = iworkspace(size_work);
        Integer* work_ptr       = work.ptr();

        Integer* ptr_A          = work_ptr;
        Integer* ptr_p          = work_ptr + Alen;

        for (Integer i = 0; i <= N; ++i)
            ptr_p[i]            = ptr_c[i] - off;
        for (Integer i = 0; i < nz; ++i)
            ptr_A[i]            = ptr_r[i];
        
        int info    = colamd(M, N, Alen, ptr_A, ptr_p, knobs, stats);
        
        if (info == 0)
        {
            throw error::error_general("colamd failed");
        }

        Integer* ptr_perm       = perm.ptr();

        for (Integer i = 0; i < N; ++i)
            ptr_perm[i]         = ptr_p[i] + 1;

        ret = permvec::from_matrix(Matrix(perm,false));
    };
    static void eval_symamd(const Mat& A, permvec& ret)
    {
        //A must be symmetric
        Integer N               = A.cols();
        Integer nz              = A.nnz();

        if (nz == 0)
        {
            ret = permvec::identity(N);
            return;
        }

        const Integer* ptr_c    = A.rep().ptr_c();
        const Integer* ptr_r    = A.rep().ptr_r();
        Integer off             = A.rep().offset();
        Mat_I sav               = Mat_I(ti::ti_empty());

        if (off != 0)
        {
            sav.resize(N+1, 1);
            Integer* ptr_c2     = sav.ptr();

            for (Integer i = 0; i <= N; ++i)
                ptr_c2[i]       = ptr_c[i] - off;

            ptr_c               = ptr_c2;
            ptr_r               += off;
        };

        Mat_I perm(ti::ti_empty(), N + 1, 1);

        double knobs[COLAMD_KNOBS];
        int stats[COLAMD_STATS];

        colamd_set_defaults(knobs);
        
        Integer* ptr_perm       = perm.ptr();

        int info    = symamd(N, (Integer*)ptr_r, (Integer*)ptr_c, ptr_perm, knobs, 
                             stats, &alloc, &release);
        
        if (info == 0)
        {
            throw error::error_general("symamd failed");
        }

        for (Integer i = 0; i < N; ++i)
            ++ptr_perm[i];

        Mat_I perm_ret  = perm.resize(N,1);

        ret = permvec::from_matrix(Matrix(perm_ret,false));
    }

    static void* alloc(size_t count, size_t size)
    {
        return calloc(count, size);
    };
    static void release(void * ptr)
    {
        free(ptr);
    };
    static void eval_mmd(const Mat& A, permvec& ret)
    {
        //A must be square
        Integer N               = A.cols();
        Integer nz              = A.nnz();

        if (nz == 0)
        {
            ret = permvec::identity(N);
            return;
        }

        const Integer* ptr_c    = A.rep().ptr_c();
        const Integer* ptr_r    = A.rep().ptr_r();
        Integer off             = A.rep().offset();
        ptr_r                   += off;        

        // delta is a parameter to allow the choice of nodes whose d
        // egree <= min-degree + DELTA
        Integer delta           = 0;
        Integer size_work       = (N + delta) * 3 + N * 2 + (N+1) + nz;
        using iworkspace        = matcl::pod_workspace<Integer>;
        iworkspace work         = iworkspace(size_work);
        Integer* work_ptr       = work.ptr();

        Integer* ptr_perm       = work_ptr;
        Integer* dhead          = work_ptr + (N + delta);
        Integer* qsize          = work_ptr + 2*(N + delta);
        Integer* llist          = work_ptr + 3*(N + delta);
        Integer* marker         = work_ptr + 3*(N + delta) + N;
        Integer* ptr_c1         = work_ptr + 3*(N + delta) + 2*N;
        Integer* ptr_r1         = work_ptr + 3*(N + delta) + 2*N + (N+1);

        /* Initialize and allocate storage for GENMMD. */
        Integer maxint          = constants::max_int();
        
        // Transform adjacency list into 1-based indexing required by GENMMD
        for (Integer i = 0; i <= N; ++i) 
            ptr_c1[i]           = ptr_c[i] - off + 1;
        for (Integer i = 0; i < nz; ++i) 
            ptr_r1[i]           = ptr_r[i] + 1;
	
        Mat_I perm(ti::ti_empty(), N + delta, 1);
        Integer* ptr_iperm      = perm.ptr();

        Integer nofsub;

        Integer info = genmmd_(&N, ptr_c1, ptr_r1, ptr_perm, ptr_iperm, &delta, dhead, 
	                        qsize, llist, marker, &maxint, &nofsub);
        
        if (info != 0)
            throw error::error_general("mmd failed");

        ret = permvec::from_matrix(Matrix(perm,false));
    }
};

template<class V, bool Is_complex = md::is_complex<V>::value>
struct make_diag_impl
{
    using VR        = typename md::real_type<V>::type;
    using Mat       = raw::Matrix<V,struct_sparse>;
    using Mat_I     = raw::Matrix<Integer,struct_dense>;
    using Mat_R     = raw::Matrix<VR,struct_dense>;
    using ret_type2 = tuple<permvec, Matrix, Matrix> ;

    static void eval_diag(const Mat& A0, permvec& p, permvec& q, Matrix& Dr, Matrix& Dc,
                          bool red_off)
    {
        Integer M               = A0.rows();
        Integer N               = A0.cols();
        Integer nz              = A0.nnz();

        if (nz == 0)
        {
            p = permvec::identity(M);
            q = permvec::identity(N);
            return;
        }
     
        bool need_A     = (red_off == true && A0.rows() == A0.cols());

        Mat A           = need_A ? A0.copy() : A0.make_unique();
        Integer* ptr_c  = A.rep().ptr_c();
        V* ptr_x        = A.rep().ptr_x();

        using workspace         = matcl::pod_workspace<V>;
        workspace work          = workspace(N);
        V* work_ptr             = work.ptr();

        for (Integer j = 0; j < N; ++j)
        {
            V max_val   = V(0.0);
            for (Integer k = ptr_c[j]; k < ptr_c[j+1]; ++k)
            {
                V val   = abs(ptr_x[k]);
                max_val = std::max(max_val, val);
            };
                        
            if (max_val == V(0.0))
            {
                work_ptr[j] = V(0.0);
                continue;
            };

            V max_log   = log(max_val);
            work_ptr[j] = max_log;

            for (Integer k = ptr_c[j]; k < ptr_c[j+1]; ++k)
            {
                ptr_x[k]= log(abs(ptr_x[k])) - max_log;
            };
        };

        Matrix row_set, col_set, u, v;
        Real f;
        Integer k;

        tie(row_set, col_set, u, v, f, k) = max_match_weighted2(Matrix(A,false));

        tie(p,q)    = rowsmatch_to_perms(row_set, A.cols());
        V* ptr_u    = u.get_array_unique<V>();
        V* ptr_v    = v.get_array_unique<V>();

        V mean_u    = V(0.0);
        V mean_v    = V(0.0);

        for (Integer i = 0; i < M; ++i)
            mean_u  += ptr_u[i];

        for (Integer i = 0; i < N; ++i)
        {
            ptr_v[i] = ptr_v[i] + work_ptr[i];
            mean_v  += ptr_v[i];
        }

        mean_u      = mean_u / M;
        mean_v      = mean_v / N;
        V adj       = (mean_u - mean_v) / 2;

        for (Integer i = 0; i < M; ++i)
            ptr_u[i]    = exp(-ptr_u[i] + adj);
        
        for (Integer i = 0; i < N; ++i)
            ptr_v[i]    = exp(-ptr_v[i] - adj);

        Dr              = u;
        Dc              = v;

        if (red_off == false || A.rows() != A.cols())
            return;

        Matrix Ac       = scale_rowscols(Matrix(A0,false), Dr, Dc);
        Ac              = Ac(p,q);

        Matrix D_off;
        tie(Ac, D_off)  = balance_offdiag(std::move(Ac));

        const VR* ptr_D         = D_off.get_array<VR>();

        const Integer* ptr_p    = p.to_array();
        const Integer* ptr_q    = q.to_array();

        for (Integer i = 0; i < N; ++i)
        {
            Integer pos = ptr_p[i] - 1;
            ptr_u[pos]  = ptr_u[pos] / ptr_D[i];
        };

        for (Integer i = 0; i < N; ++i)
        {
            Integer pos = ptr_q[i] - 1;
            ptr_v[pos]  = ptr_v[pos] * ptr_D[i];
        };
    }
};

template<class V>
struct make_diag_impl<V, true>
{
    using VR        = typename md::real_type<V>::type;
    using Mat       = raw::Matrix<V,struct_sparse>;
    using Mat_R     = raw::Matrix<VR,struct_sparse>;
    using ret_type2 = tuple<permvec, Matrix, Matrix> ;

    static void eval_diag(const Mat& A0, permvec& p, permvec& q, Matrix& Dr, Matrix& Dc, 
                          bool red_off)
    {
        Matrix Ac       = abs(Matrix(A0, false));
        const Mat_R& A  = Ac.impl<Mat_R>();
        return make_diag_impl<VR>::eval_diag(A, p, q, Dr, Dc, red_off);
    };
};

template<>
struct make_diag_impl<Object>
{
    using Mat       = raw::Matrix<Object,struct_sparse>;
    using ret_type2 = tuple<permvec, Matrix, Matrix> ;

    static void eval_diag(const Mat&, permvec&, permvec&, Matrix&, Matrix&, bool)
    {
        throw error::object_value_type_not_allowed("order_diag");
    };
};

template<>
struct make_diag_impl<Integer>
{
    using Mat       = raw::Matrix<Integer,struct_sparse>;
    using ret_type2 = tuple<permvec, Matrix, Matrix> ;

    static void eval_diag(const Mat& A0, permvec& p, permvec& q, Matrix& Dr, Matrix& Dc,
                          bool red_off)
    {
        using Mat_R = raw::Matrix<Real,struct_sparse>;
        Mat_R AR    = raw::converter<Mat_R,Mat>::eval(A0);
        return make_diag_impl<Real>::eval_diag(AR, p, q, Dr, Dc, red_off);
    };
};

struct visit_mmd : public md::extract_type_switch<void, visit_mmd,true>
{
    using base_type = md::extract_type_switch<void, visit_mmd,true>;

    template<class V>
    static void eval(const matcl::Matrix&, const raw::Matrix<V,struct_sparse>& mat,
                     permvec& ret, details::mmd_type type)
    {
        return make_mmd_impl<V>::eval(mat, ret, type);
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

struct visit_diag : public md::extract_type_switch<void, visit_diag,true>
{
    using base_type = md::extract_type_switch<void, visit_diag,true>;

    template<class V>
    static void eval(const matcl::Matrix&, const raw::Matrix<V,struct_sparse>& mat,
                     permvec& p, permvec& q, Matrix& Dr, Matrix& Dc, bool red_off)
    {
        return make_diag_impl<V>::eval_diag(mat, p, q, Dr, Dc, red_off);
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

/// nested dissection ordering of a square, symmetric matrix
permvec matcl::order_nested_dissection(const Matrix& A)
{
    if (A.is_square() == false)
        throw error::square_matrix_required(A.rows(), A.cols());

    //TODO
    throw;
    //return petsc::reorder_sym(A, petsc_ordering::nd);
};

/// one-way dissection ordering of a square, symmetric matrix
permvec matcl::order_one_way_dissection(const Matrix& A)
{
    if (A.is_square() == false)
        throw error::square_matrix_required(A.rows(), A.cols());

    //TODO
    throw;
    //return petsc::reorder_sym(A, petsc_ordering::w1d);
}

/// reverse Cuthill-McKee ordering of a square, symmetric matrix
permvec matcl::order_rcm(const Matrix& A)
{
    if (A.is_square() == false)
        throw error::square_matrix_required(A.rows(), A.cols());

    //TODO
    throw;
    //return petsc::reorder_sym(A, petsc_ordering::rcm);
}

permvec matcl::order_spectral(const Matrix& A, bool use_weights, bool direct)
{
    Matrix L    = symgraph_laplacian(A, use_weights);
    
    // when the graph is connected, then one Fielder vectors is enough
    // to determine ordering uniquely, more vectors do not hurt and they are not
    // costly
    Matrix E, V;
    tie(V,E)    = fiedler_vectors(L, 4, direct);

    Matrix x, I;
    tie(x,I)    = sortrows2(V(colon(), colon(2,4)));

    return permvec::from_matrix(I);
};

permvec matcl::order_scotch(const Matrix& A)
{
    scotch sc(A);
    return sc.ordering();
}

permvec matcl::order_metis(const Matrix& A)
{
    metis sc(A);
    return sc.ordering();
}

tuple<permvec, permvec> matcl::dmperm(const Matrix& A)
{
    dm_decomp dm(A);
    permvec p   = dm.row_perms();
    permvec q   = dm.col_perms();

    return tuple<permvec, permvec>(p,q);
};

permvec matcl::colperm(const Matrix& A)
{
    Matrix nz   = matcl::nnz(A,1);
    Matrix v, I;
    tie(v,I)    = sort2(nz, 2);

    return permvec::from_matrix(I);
};

permvec matcl::rowperm(const Matrix& A)
{
    Matrix nz   = matcl::nnz(A,2);
    Matrix v, I;
    tie(v,I)    = sort2(nz, 1);

    return permvec::from_matrix(I);
};

permvec matcl::order_mmd(const Matrix& A)
{
    if (A.is_square() == false)
        throw error::square_matrix_required(A.rows(), A.cols());

    Matrix As   = make_adjancency_matrix(A);

    permvec p;
    details::visit_mmd::make<const Matrix&>(As, p, details::mmd_type::mmd);
    return p;
};

permvec matcl::order_colamd(const Matrix& A)
{
    permvec p;
    details::visit_mmd::make<const Matrix&>(A, p, details::mmd_type::colamd);
    return p;
};

permvec matcl::order_symamd(const Matrix& A)
{
    if (A.is_square() == false)
        throw error::square_matrix_required(A.rows(), A.cols());

    Matrix As;
    if (A.get_struct().has_symher_flag())
        As  = A;
    else
        As = symgraph_sum(A);

    permvec p;
    details::visit_mmd::make<const Matrix&>(As, p, details::mmd_type::symamd);
    return p;
};

tuple<permvec,permvec, Matrix, Matrix> matcl::order_diag(const Matrix& A, bool reduce_offdiags)
{
    Matrix As(abs(A));

    permvec p, q;
    Matrix Dr, Dc;
    details::visit_diag::make<const Matrix&>(As, p, q, Dr, Dc, reduce_offdiags);

    return tuple<permvec,permvec, Matrix, Matrix>(p,q,Dr, Dc);
};

tuple<permvec,permvec, Matrix, Matrix> matcl::order_diag(Matrix&& A, bool reduce_offdiags)
{
    Matrix As = abs(std::move(A));

    permvec p, q;
    Matrix Dr, Dc;
    details::visit_diag::make<const Matrix&>(As, p, q, Dr, Dc, reduce_offdiags);

    return tuple<permvec,permvec, Matrix, Matrix>(p,q,Dr, Dc);
};

};