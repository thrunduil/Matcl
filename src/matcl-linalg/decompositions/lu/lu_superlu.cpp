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

#include "lu_superlu.h"

#include "matcl-core/utils/workspace.h"
#include "matcl-linalg/options/options_linsolve.h"
#include "matcl-linalg/decompositions/lu/lu_struct.h"
#include "matcl-linalg/special_matrices/struct_flag_linalg.h"
#include "matcl-linalg/graph/matcl_graph.h"

namespace matcl { namespace details
{

template<class V>
void superlu_wrap<V>::eval_ilu(const Mat& A, ret_type& ret, const options& opts)
{
    Integer M               = A.rows();
    Integer N               = A.cols();
    Integer K               = std::min(M, N);

    //superlu aborts on horizontally rectangular matrices
    if (M < N)
        throw error::invalid_size2(M, N, N, N);

    if (A.nnz() == 0)
    {
        using VR            = typename md::real_type<V>::type;
        value_code vc       = matrix_traits::value_code<VR>::value;
        Matrix L            = speye(M, K, vc);
        Matrix U            = spzeros(K, N, 0, vc);
        permvec p           = permvec::identity(M);
        permvec q           = permvec::identity(N);

        ret                 = ret_type(L, U, p, q);
        return;
    };

    Integer A_ldiags        = get_ld(Matrix(A,false),0);
    Integer A_udiags        = get_ud(Matrix(A,false),0);  

    namespace ol = opt::linsolve;

    bool allow_nonunit      = opts.get_option<bool>(ol::allow_nonunit_L());

    if (A_ldiags == 0 && A_udiags == 0)
        return lu_diag<V,struct_sparse>::eval(ret, A, ol::pivot_type::partial, 0.0);

    if (A_udiags == 0)
    {
        if (allow_nonunit)
            return lu_tril_nonunit(ret, Matrix(A,false));
    }
    else if (A_ldiags == 0)
    {
        return lu_triu(ret, Matrix(A,false));
    }

    superlu_options_t options;
    ilu_set_default_options(&options);
    make_options(A, options, opts, true);

    // Initialize the statistics variables
    SuperLUStat_t* stat     = new SuperLUStat_t();
    StatInit(stat);

    stat_ptr sp;
    sp.reset(stat, &deleter_stat);

    // Create matrix A in the format expected by SuperLU
    supermatrix_ptr sl_A;
    create_supermatrix(sl_A, A, M);

    // fill explicit zeros on the diagonal entries, so that the matrix is not
    // structurally singular        
    superlu_interface<V>::fill_diag(N, (NCformat*)sl_A->Store);
        
    int panel_size  = sp_ienv(1);
    int relax       = sp_ienv(2);
    
    Mat_I perm_r(ti::ti_empty(), M, 1);
    Integer* ptr_perm_r = perm_r.ptr();

    // column permutation
    Mat_I perm_c    = Mat_I(ti::ti_empty());

    do_permutation(perm_c, options, opts, Matrix(A,false), sl_A);
    Integer* ptr_perm_c     = perm_c.ptr();

    Integer isize_work      = std::max(M,N);
    using iworkspace        = matcl::pod_workspace<Integer>;
    iworkspace iwork        = iworkspace(isize_work);
    Integer* iwork_ptr      = iwork.ptr();
    Integer* etree          = iwork_ptr;

    SuperMatrix* sl_Ac      = new SuperMatrix();
    sp_preorder(&options, sl_A.get(), ptr_perm_c, etree, sl_Ac);

    supermatrix_ptr sm_Ac;
    sm_Ac.reset(sl_Ac, &deleter_permuted);

    GlobalLU_t Glu;
    SuperMatrix* L  = new SuperMatrix();
    SuperMatrix* U  = new SuperMatrix();

    int info;
    superlu_interface<V>::gsitrf(&options, sm_Ac.get(), relax, panel_size, etree, nullptr, 0, 
                                ptr_perm_c, ptr_perm_r, L, U, &Glu, stat, &info);

    supermatrix_ptr sm_L;
    supermatrix_ptr sm_U;

    sm_L.reset(L, &deleter_supernode);
    sm_U.reset(U, &deleter_comp_col);

    if (info < 0)
        throw error::error_general("invalid argument supplied to gsitrf");

    Matrix mat_L, mat_U;

    make_factors_matcl(sm_L, sm_U, mat_L, mat_U); 
    //correct_rperm(M, ptr_perm_r, iwork_ptr);

    for (Integer i = 0; i < M; ++i)
        ++ptr_perm_r[i];
    for (Integer i = 0; i < N; ++i)
        ++ptr_perm_c[i];

    permvec p   = permvec::from_matrix(Matrix(perm_r,false));
    permvec q   = permvec::from_matrix(Matrix(perm_c,false));

    p           = p.invperm();
    q           = q.invperm();

    mat_L.add_struct(predefined_struct_type::tril);
    mat_U.add_struct(predefined_struct_type::triu);

    ret         = ret_type(mat_L, mat_U, p, q);
};

template<class V>
void superlu_wrap<V>::eval_lu(const Mat& A, ret_type& ret, const options& opts)
{
    Integer M               = A.rows();
    Integer N               = A.cols();
    Integer K               = std::min(M, N);

    //superlu is broken for M < N
    if (M < N)
        throw error::invalid_size2(M, N, N, N);

    if (A.nnz() == 0)
    {
        using VR            = typename md::real_type<V>::type;
        value_code vc       = matrix_traits::value_code<VR>::value;
        Matrix L            = speye(M, K, vc);
        Matrix U            = spzeros(K, N, 0, vc);
        permvec p           = permvec::identity(M);
        permvec q           = permvec::identity(N);

        ret                 = ret_type(L, U, p, q);
        return;
    };

    Integer A_ldiags        = get_ld(Matrix(A,false),0);
    Integer A_udiags        = get_ud(Matrix(A,false),0);  

    namespace ol = opt::linsolve;

    bool allow_nonunit      = opts.get_option<bool>(ol::allow_nonunit_L());

    if (A_ldiags == 0 && A_udiags == 0)
        return lu_diag<V,struct_sparse>::eval(ret, A, ol::pivot_type::partial, 0.0);

    if (A_udiags == 0)
    {
        if (allow_nonunit)
            return lu_tril_nonunit(ret, Matrix(A,false));
    }
    else if (A_ldiags == 0)
    {
        return lu_triu(ret, Matrix(A,false));
    }
    
    superlu_options_t options;
    set_default_options(&options);
    make_options(A, options, opts, false);

    // Initialize the statistics variables
    SuperLUStat_t* stat     = new SuperLUStat_t();
    StatInit(stat);

    stat_ptr sp;
    sp.reset(stat, &deleter_stat);

    // Create matrix A in the format expected by SuperLU
    supermatrix_ptr sl_A;
    create_supermatrix(sl_A, A, M);
        
    int panel_size  = sp_ienv(1);
    int relax       = sp_ienv(2);

    // fill explicit zeros on the diagonal entries, so that the matrix is not
    // structurally singular        
    superlu_interface<V>::fill_diag(std::min(N, M), (NCformat*)sl_A->Store);

    Mat_I perm_r(ti::ti_empty(), M, 1);
    Integer* ptr_perm_r = perm_r.ptr();

    // column permutation
    Mat_I perm_c    = Mat_I(ti::ti_empty());

    do_permutation(perm_c, options, opts, Matrix(A,false), sl_A);
    Integer* ptr_perm_c     = perm_c.ptr();

    Integer isize_work      = std::max(N, M);
    using iworkspace        = matcl::pod_workspace<Integer>;
    iworkspace iwork        = iworkspace(isize_work);
    Integer* iwork_ptr      = iwork.ptr();
    Integer* etree          = iwork_ptr;

    SuperMatrix* sl_Ac      = new SuperMatrix();
    sp_preorder(&options, sl_A.get(), ptr_perm_c, etree, sl_Ac);

    supermatrix_ptr sm_Ac;
    sm_Ac.reset(sl_Ac, &deleter_permuted);

    GlobalLU_t Glu;
    SuperMatrix* L  = new SuperMatrix();
    SuperMatrix* U  = new SuperMatrix();

    int info;
    superlu_interface<V>::gstrf(&options, sm_Ac.get(), relax, panel_size, etree, nullptr, 0, 
                                ptr_perm_c, ptr_perm_r, L, U, &Glu, stat, &info);

    supermatrix_ptr sm_L;
    supermatrix_ptr sm_U;

    sm_L.reset(L, &deleter_supernode);
    sm_U.reset(U, &deleter_comp_col);

    if (info < 0)
        throw error::error_general("invalid argument supplied to gsitrf");

    Matrix mat_L, mat_U;

    make_factors_matcl(sm_L, sm_U, mat_L, mat_U); 
    //correct_rperm(M, ptr_perm_r, iwork_ptr);

    for (Integer i = 0; i < M; ++i)
        ++ptr_perm_r[i];
    for (Integer i = 0; i < N; ++i)
        ++ptr_perm_c[i];

    permvec p   = permvec::from_matrix(Matrix(perm_r,false));
    permvec q   = permvec::from_matrix(Matrix(perm_c,false));

    p           = p.invperm();
    q           = q.invperm();

    mat_L.add_struct(predefined_struct_type::tril);
    mat_U.add_struct(predefined_struct_type::triu);

    ret         = ret_type(mat_L, mat_U, p, q);
};

template<class V>
void superlu_wrap<V>::correct_rperm(Integer M, Integer* ptr_perm_r, Integer* iwork_ptr)
{
    for (Integer i = 0; i < M; ++i)
        iwork_ptr[i]    = 0;

    bool missing        = false;
    for (Integer i = 0; i < M; ++i)
    {
        Integer p       = ptr_perm_r[i];
        if (p >= M)
        {
            ptr_perm_r[i]= -1;
            missing     = true;
        }
        else if (p >= 0)
        {
            iwork_ptr[p]= 1;
        }
        else
        {
            missing     = true;
        }
    };

    if (missing == false)
        return;

    Integer pos         = 0;

    for (Integer i = 0; i < M; ++i)
    {
        Integer p       = ptr_perm_r[i];
        if (p >= 0)
            continue;

        while(iwork_ptr[pos] == 1)
            ++pos;

        ptr_perm_r[i]   = pos;
        ++pos;
    };

};

template<class V>
void superlu_wrap<V>::do_permutation(Mat_I& perm_c, const superlu_options_t& opts, const options& matcl_opts, 
                                     const Matrix& A, supermatrix_ptr& sl_A)
{
    Integer* ptr_perm_c = nullptr;
    int permc_spec      = opts.ColPerm;
    bool do_sl_perm     = false;
    Integer N           = A.cols();

    namespace ol = opt::linsolve;

    if (permc_spec == NATURAL)
    {
        Integer opt_val     = matcl_opts.get_option<Integer>(ol::lu_ordering());
        ol::lu_ordering_type ord  = (ol::lu_ordering_type)opt_val;

        permvec pc;
        switch (ord)
        {
            case ol::lu_ordering_type::natural:
            default:
            {
                do_sl_perm  = true;
                break;            
            }
            case ol::lu_ordering_type::colamd:
            {
                pc  = order_colamd(A);
                break;
            }
        };

        if (do_sl_perm == false)
        {
            pc          = pc.invperm();
            Matrix pm   = pc.to_matrix();
            perm_c.assign_to_fresh(pm.impl<Mat_I>());
            ptr_perm_c  = perm_c.ptr();

            for (Integer i = 0; i < N; ++i)
                --ptr_perm_c[i];
        };
    }
    else
    {
        do_sl_perm = true;
    };

    if (do_sl_perm == true)
    {
        perm_c.assign_to_fresh(Mat_I(ti::ti_empty(), N, 1));
        ptr_perm_c          = perm_c.ptr();

        get_perm_c(permc_spec, sl_A.get(), ptr_perm_c);
    };
};

template<class V>
void superlu_wrap<V>::make_options(const Mat& A, superlu_options_t& opts, const options& matcl_opts,
                                   bool ilu)
{    
    namespace ol = opt::linsolve;

    bool is_pos             = is_posdef(A.get_struct()) || is_semi_posdef(A.get_struct());

    if (matcl_opts.has_option(ol::tol_r()))
    {
        Real diag_thresh    = matcl_opts.get_option<Real>(ol::tol_r());
        opts.DiagPivotThresh= diag_thresh;
    };
    if (is_pos)
    {
        opts.SymmetricMode  = YES;
    }
        
    opts.ColPerm            = NATURAL;
    opts.Equil              = NO;
    opts.IterRefine         = NOREFINE;
    opts.PivotGrowth        = NO;
    opts.ConditionNumber    = NO;
    opts.RowPerm            = LargeDiag;
    opts.PrintStat          = NO;

    if (ilu == false)
        return;

    if (matcl_opts.has_option(ol::drop_tol()))
    {
        Real tol            = matcl_opts.get_option<Real>(ol::drop_tol());
        opts.ILU_DropTol    = tol;
    };
    if (matcl_opts.has_option(ol::diag_fill_tol()))
    {
        Real tol            = matcl_opts.get_option<Real>(ol::diag_fill_tol());
        opts.ILU_FillTol    = tol;
    };
    if (matcl_opts.has_option(ol::fill_factor()))
    {
        Real fill           = matcl_opts.get_option<Real>(ol::fill_factor());
        opts.ILU_FillFactor = fill;
    };    
    if (matcl_opts.has_option(ol::ilu_type()))
    {
        Integer val         = matcl_opts.get_option<Integer>(ol::ilu_type());
        ol::ilu_method met  = (ol::ilu_method)val;

        switch (met)
        {
            case ol::ilu_method::silu:
                opts.ILU_MILU = SILU;
                break;
            case ol::ilu_method::milu:
                opts.ILU_MILU = SMILU_2;
                break;
            case ol::ilu_method::milu_abs:
                opts.ILU_MILU = SMILU_3;
                break;
        }
    };    
};

template<class V>
void superlu_wrap<V>::create_supermatrix(supermatrix_ptr& sm, const Mat& A, Integer n_rows)
{
    //we need to make the matrix A square (or vertically rectangular)
    Integer M               = std::max(n_rows, A.rows());
    Integer N               = A.cols();
    Integer nz              = A.nnz();
    Integer off             = A.rep().offset();

    const V* ptr_x          = A.rep().ptr_x() + off;
    const Integer* ptr_c    = A.rep().ptr_c();
    const Integer* ptr_r    = A.rep().ptr_r() + off;        

    // copy A to superlu arrays
    Integer* ptr_c_su       = intMalloc(N+1);
    if (!ptr_c_su)
        throw error::alloc(N+1);

    Integer* ptr_r_su       = intMalloc(nz);
    if (!ptr_r_su)
    {
        superlu_free(ptr_c_su);
        throw error::alloc(nz);
    }

    V* ptr_x_su             = (V*)superlu_malloc(nz * sizeof(V));
    if (!ptr_x_su)
    {
        superlu_free(ptr_c_su);
        superlu_free(ptr_r_su);
        throw error::alloc(nz);
    }

    for (Integer i = 0; i <= N; ++i)
        ptr_c_su[i]         = ptr_c[i] - off;

    for (Integer i = 0; i < nz; ++i)
        ptr_r_su[i]         = ptr_r[i];

    for (Integer i = 0; i < nz; ++i)
        ptr_x_su[i]         = ptr_x[i];

    SuperMatrix* sl_A       = new SuperMatrix();

    superlu_interface<V>::Create_CompCol_Matrix(sl_A, M, N, nz, ptr_x_su, ptr_r_su, ptr_c_su, 
                                                SLU_NC, SLU_GE);

    sm.reset(sl_A, &deleter_comp_col);
};

template<class V>
void superlu_wrap<V>::deleter_comp_col(SuperMatrix* sm)
{
    Destroy_CompCol_Matrix(sm);
    delete sm;
};
template<class V>
void superlu_wrap<V>::deleter_supernode(SuperMatrix* sm)
{
    Destroy_SuperNode_Matrix(sm);
    delete sm;
};
template<class V>
void superlu_wrap<V>::deleter_permuted(SuperMatrix* sm)
{
    Destroy_CompCol_Permuted(sm);
    delete sm;
};
template<class V>
void superlu_wrap<V>::deleter_stat(SuperLUStat_t* stat)
{
    StatFree(stat);
    delete stat;
};

template<class V>
void superlu_wrap<V>::make_factors_matcl(const supermatrix_ptr& L, const supermatrix_ptr& U, 
                                   Matrix& mat_L, Matrix& mat_U)
{
    //dPrint_SuperNode_Matrix("L", (SuperMatrix*)SN.get());

    Integer ML              = L->nrow;
    Integer NL              = L->ncol;    
    Integer MU              = U->nrow;
    Integer NU              = U->ncol;    
    SCformat* Lstore        = (SCformat *) L->Store;
    NCformat* Ustore        = (NCformat *) U->Store;    

    const V* ptr_x_sl           = (V*)Lstore->nzval;
    const Integer* ptr_r_sl     = Lstore->rowind;
    const Integer* ptr_cr_sl    = Lstore->rowind_colptr;
    const Integer* ptr_cx_sl    = Lstore->nzval_colptr;
    const Integer* sup_to_col   = Lstore->sup_to_col;
    const V* ptr_x_su           = (V*)Ustore->nzval;
    const Integer* ptr_c_su     = Ustore->colptr;
    const Integer* ptr_r_su     = Ustore->rowind;

    Integer nz_l                = Lstore->nnz;
    Integer nz_u                = Ustore->nnz;
    Integer nsuper              = Lstore->nsuper;

    Mat rep_L(ti::ti_empty(), ML, NL, nz_l);
    Mat rep_U(ti::ti_empty(), MU, NU, nz_u);
    Integer* ptr_cl         = rep_L.rep().ptr_c();
    Integer* ptr_rl         = rep_L.rep().ptr_r();
    V* ptr_xl               = rep_L.rep().ptr_x();
    Integer* ptr_cu         = rep_U.rep().ptr_c();
    Integer* ptr_ru         = rep_U.rep().ptr_r();
    V* ptr_xu               = rep_U.rep().ptr_x();

    nz_l                    = 0;
    nz_u                    = 0;
    Integer c               = 0;    

    for (Integer ns = 0; ns <= nsuper; ++ns) 
    {
        Integer upper       = 1;
        Integer col         = sup_to_col[ns];
        Integer nsup        = sup_to_col[ns+1] - col;
        Integer cl          = std::min(NU, col + nsup);

        for (c = col; c < cl; ++c) 
        {
            ptr_cl[c]       = nz_l;
            ptr_cu[c]       = nz_u;

            for (Integer ku = ptr_c_su[c]; ku < ptr_c_su[c+1]; ++ku)
            {
                if (mrd::is_zero(ptr_x_su[ku]) == true)
                    continue;

                ptr_ru[nz_u]    = ptr_r_su[ku];
                ptr_xu[nz_u]    = ptr_x_su[ku];
                ++nz_u;
            };

            Integer d           = ptr_cx_sl[c];

	        for (Integer kl = ptr_cr_sl[col]; kl < ptr_cr_sl[col] + upper; ++kl, ++d)
            {
                Integer r       = ptr_r_sl[kl];
                V val           = ptr_x_sl[d];

                if (mrd::is_zero(val) == true)
                    continue;

                ptr_ru[nz_u]= r;
                ptr_xu[nz_u]= val;
                ++nz_u;
                continue;
            };
	        
            ptr_rl[nz_l]        = c;
            ptr_xl[nz_l]        = V(1.0);
            ++nz_l;

	        for (Integer kl = ptr_cr_sl[col] + upper; kl < ptr_cr_sl[col+1]; ++kl, ++d)
            {
                Integer r       = ptr_r_sl[kl];
                V val           = ptr_x_sl[d];

                if (mrd::is_zero(val) == true)
                    continue;

                ptr_rl[nz_l]    = r;
                ptr_xl[nz_l]    = val;

                ++nz_l;
            };

            ++upper;
        }

        for (; c < col + nsup; ++c) 
        {
            ptr_cl[c]       = nz_l;

	        Integer d       = ptr_cx_sl[c];

            ptr_rl[nz_l]    = c;
            ptr_xl[nz_l]    = V(1.0);

            ++nz_l;           

	        for (Integer kl = ptr_cr_sl[col]; kl < ptr_cr_sl[col+1]; ++kl, ++d)
            {
                Integer r       = ptr_r_sl[kl];
                V val           = ptr_x_sl[d];

                if (mrd::is_zero(val) == true)
                    continue;

                ptr_rl[nz_l]    = r;
                ptr_xl[nz_l]    = val;

                ++nz_l;
            };
        }
    }

    for (; c < NU; ++c) 
    {
        ptr_cu[c]       = nz_u;

        for (Integer ku = ptr_c_su[c]; ku < ptr_c_su[c+1]; ++ku)
        {
            if (mrd::is_zero(ptr_x_su[ku]) == true)
                continue;

            ptr_ru[nz_u]    = ptr_r_su[ku];
            ptr_xu[nz_u]    = ptr_x_su[ku];
            ++nz_u;
        };
    }

    ptr_cu[NU]  = nz_u;
    ptr_cl[NL]  = nz_l;

    //superlu is not sorted
    rep_L.rep().sort();
    rep_U.rep().sort();

    mat_L = Matrix(rep_L, true);
    mat_U = Matrix(rep_U, true);
};

template class superlu_wrap<Real>;
template class superlu_wrap<Float>;
template class superlu_wrap<Complex>;
template class superlu_wrap<Float_complex>;

}};
