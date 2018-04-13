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

#pragma once

#include "matcl-matrep/matcl_matrep.h"

#include "matcl-internals/container/mat_s.h"
#include "matcl-internals/container/mat_d.h"
#include "matcl-internals/container/mat_b.h"

extern "C"
{
    #include "extern/superlu/SuperLU_5.2.0/SRC/slu_ddefs.h"    
    #include "extern/superlu/SuperLU_5.2.0/SRC/slu_sdefs.h"    
    #include "extern/superlu/SuperLU_5.2.0/SRC/slu_zdefs.h"    
    #include "extern/superlu/SuperLU_5.2.0/SRC/slu_cdefs.h"    
};

extern "C"
{
    int dfill_diag(int n, NCformat *Astore);
    int sfill_diag(int n, NCformat *Astore);
    int zfill_diag(int n, NCformat *Astore);
    int cfill_diag(int n, NCformat *Astore);
};

namespace matcl { namespace details
{

template<class V>
class superlu_wrap
{
    private:
        using supermatrix_ptr   = std::shared_ptr<SuperMatrix>;
        using stat_ptr          = std::shared_ptr<SuperLUStat_t>;
        using ret_type          = tuple<Matrix,Matrix,permvec,permvec>;

    public:
        using Mat   = raw::Matrix<V,struct_sparse>;
        using Mat_I = raw::Matrix<Integer,struct_dense>;

        static void     eval_ilu(const Mat& A, ret_type& ret, const options& opts);
        static void     eval_lu(const Mat& A, ret_type& ret, const options& opts);

    private:
        static void     make_options(const Mat& A, superlu_options_t& opts, const options& matcl_opts,
                            bool il);
        static void     do_permutation(Mat_I& perm_c, const superlu_options_t& opts, 
                            const options& matcl_opts, const Matrix& A, supermatrix_ptr& sl_A);
        static void     create_supermatrix(supermatrix_ptr& sl_A, const Mat& A, Integer n_rows);
        static void     deleter_comp_col(SuperMatrix* sm);
        static void     deleter_supernode(SuperMatrix* sm);
        static void     deleter_permuted(SuperMatrix* sm);
        static void     deleter_stat(SuperLUStat_t* stat);
        static void     make_factors_matcl(const supermatrix_ptr& L, const supermatrix_ptr& U, 
                            Matrix& mat_L, Matrix& mat_U);
        static void     correct_rperm(Integer M, Integer* ptr_perm_r, Integer* iwork_ptr);
};

template<class V>
struct superlu_interface
{};

template<>
struct superlu_interface<Real>
{
    static void Create_CompCol_Matrix(SuperMatrix *A, int m, int n, int nnz, double *nzval, int *rowind, 
                                      int *colptr, Stype_t stype, Mtype_t mtype)
    {
        return dCreate_CompCol_Matrix(A, m, n, nnz, nzval, rowind, colptr, stype, SLU_D, mtype);
    };
    static int fill_diag(int n, NCformat *Astore)
    {
        return dfill_diag(n, Astore);
    };
    static void gsitrf(superlu_options_t *options, SuperMatrix *A, int relax, int panel_size,
	                int *etree, Real *work, int lwork, int *perm_c, int *perm_r,
	                SuperMatrix *L, SuperMatrix *U, GlobalLU_t *Glu, SuperLUStat_t *stat, int *info)
    {
        return dgsitrf(options, A, relax, panel_size, etree, work, lwork, perm_c, perm_r, 
                       L, U, Glu, stat, info);
    };
    static void gstrf(superlu_options_t *options, SuperMatrix *A, int relax, int panel_size,
	                int *etree, Real *work, int lwork, int *perm_c, int *perm_r,
	                SuperMatrix *L, SuperMatrix *U, GlobalLU_t *Glu, SuperLUStat_t *stat, int *info)
    {
        return dgstrf(options, A, relax, panel_size, etree, work, lwork, perm_c, perm_r, 
                    L, U, Glu, stat, info);
    };
};

template<>
struct superlu_interface<Float>
{
    static void Create_CompCol_Matrix(SuperMatrix *A, int m, int n, int nnz, Float *nzval, int *rowind, 
                                      int *colptr, Stype_t stype, Mtype_t mtype)
    {
        return sCreate_CompCol_Matrix(A, m, n, nnz, nzval, rowind, colptr, stype, SLU_S, mtype);
    };
    static int fill_diag(int n, NCformat *Astore)
    {
        return sfill_diag(n, Astore);
    };
    static void gsitrf(superlu_options_t *options, SuperMatrix *A, int relax, int panel_size,
	                int *etree, Float *work, int lwork, int *perm_c, int *perm_r,
	                SuperMatrix *L, SuperMatrix *U, GlobalLU_t *Glu, SuperLUStat_t *stat, int *info)
    {
        return sgsitrf(options, A, relax, panel_size, etree, work, lwork, perm_c, perm_r, 
                       L, U, Glu, stat, info);
    };
    static void gstrf(superlu_options_t *options, SuperMatrix *A, int relax, int panel_size,
	                int *etree, Float *work, int lwork, int *perm_c, int *perm_r,
	                SuperMatrix *L, SuperMatrix *U, GlobalLU_t *Glu, SuperLUStat_t *stat, int *info)
    {
        return sgstrf(options, A, relax, panel_size, etree, work, lwork, perm_c, perm_r, 
                       L, U, Glu, stat, info);
    };
};

template<>
struct superlu_interface<Complex>
{
    static void Create_CompCol_Matrix(SuperMatrix *A, int m, int n, int nnz, Complex *nzval, int *rowind, 
                                      int *colptr, Stype_t stype, Mtype_t mtype)
    {
        return zCreate_CompCol_Matrix(A, m, n, nnz, (doublecomplex*)nzval, rowind, colptr, stype, SLU_Z, mtype);
    };
    static int fill_diag(int n, NCformat *Astore)
    {
        return zfill_diag(n, Astore);
    };
    static void gsitrf(superlu_options_t *options, SuperMatrix *A, int relax, int panel_size,
	                int *etree, Complex *work, int lwork, int *perm_c, int *perm_r,
	                SuperMatrix *L, SuperMatrix *U, GlobalLU_t *Glu, SuperLUStat_t *stat, int *info)
    {
        return zgsitrf(options, A, relax, panel_size, etree, work, lwork, perm_c, perm_r, 
                       L, U, Glu, stat, info);
    };
    static void gstrf(superlu_options_t *options, SuperMatrix *A, int relax, int panel_size,
	                int *etree, Complex *work, int lwork, int *perm_c, int *perm_r,
	                SuperMatrix *L, SuperMatrix *U, GlobalLU_t *Glu, SuperLUStat_t *stat, int *info)
    {
        return zgstrf(options, A, relax, panel_size, etree, work, lwork, perm_c, perm_r, 
                       L, U, Glu, stat, info);
    };
};

template<>
struct superlu_interface<Float_complex>
{
    static void Create_CompCol_Matrix(SuperMatrix *A, int m, int n, int nnz, Float_complex *nzval, int *rowind, 
                                      int *colptr, Stype_t stype, Mtype_t mtype)
    {
        return cCreate_CompCol_Matrix(A, m, n, nnz, (::complex*)nzval, rowind, colptr, stype, SLU_C, mtype);
    };
    static int fill_diag(int n, NCformat *Astore)
    {
        return cfill_diag(n, Astore);
    };
    static void gsitrf(superlu_options_t *options, SuperMatrix *A, int relax, int panel_size,
	                int *etree, Float_complex *work, int lwork, int *perm_c, int *perm_r,
	                SuperMatrix *L, SuperMatrix *U, GlobalLU_t *Glu, SuperLUStat_t *stat, int *info)
    {
        return cgsitrf(options, A, relax, panel_size, etree, work, lwork, perm_c, perm_r, 
                       L, U, Glu, stat, info);
    };
    static void gstrf(superlu_options_t *options, SuperMatrix *A, int relax, int panel_size,
	                int *etree, Float_complex *work, int lwork, int *perm_c, int *perm_r,
	                SuperMatrix *L, SuperMatrix *U, GlobalLU_t *Glu, SuperLUStat_t *stat, int *info)
    {
        return cgstrf(options, A, relax, panel_size, etree, work, lwork, perm_c, perm_r, 
                       L, U, Glu, stat, info);
    };
};

}};

