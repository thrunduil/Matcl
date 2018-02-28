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

#include "matcl-matmult/func/raw/raw_scale.h"
#include "matcl-internals/base/instantiate.h"
#include "matcl-matrep/matrix/struct_flag_ext.h"

#include "matcl-internals/container/mat_d.h"
#include "matcl-internals/container/mat_s.h"
#include "matcl-internals/container/mat_b.h"
#include "matcl-matrep/lib_functions/func_binary.h"
#include "matcl-matrep/matrix/struct_flag_ext.h"

namespace matcl { namespace raw { namespace details
{

template<class V1, class V2, class S>
struct scaling_helper_impl
{};

template<class V1, class V2>
struct scaling_helper_impl<V1,V2,struct_dense>
{
    using MP    = raw::Matrix<V1,struct_dense>;
    using Mat_D = raw::Matrix<V2,struct_dense>;

    static void eval_rows(matcl::Matrix& ret, const MP& A, const Mat_D& Dr)
    {
        Integer M           = A.rows();
        Integer N           = A.cols();

        MP res              = A.make_unique();
        
        V1* ptr_res         = res.ptr();
        const V1* ptr_A     = A.ptr();
        const V2* Dr_ptr    = Dr.ptr();

        Integer res_ld      = res.ld();
        Integer A_ld        = A.ld();
        Integer Dr_inc      = (Dr.cols() == M)? Dr.ld() : 1;

        if (Dr_inc == 1)
        {
            for (Integer j = 0; j < N; ++j)
            {
                for (Integer i = 0; i < M; ++i)
                    ptr_res[i]  = ptr_A[i] * Dr_ptr[i];

                ptr_res         += res_ld;
                ptr_A           += A_ld;
            };
        }
        else
        {
            for (Integer j = 0; j < N; ++j)
            {
                for (Integer i = 0; i < M; ++i)
                    ptr_res[i]  = ptr_A[i] * Dr_ptr[i*Dr_inc];

                ptr_res         += res_ld;
                ptr_A           += A_ld;
            };
        };

        struct_flag sf      = res.get_struct();
        sf                  = predefined_struct_ext::get_scal_rowcol(sf);
        res.set_struct(sf);

        ret                 = matcl::Matrix(res,true);
        return;
    };

    static void eval_cols(matcl::Matrix& ret, const MP& A, const Mat_D& Dc)
    {
        Integer M           = A.rows();
        Integer N           = A.cols();

        MP res              = A.make_unique();
        
        V1* ptr_res         = res.ptr();
        const V1* ptr_A     = A.ptr();
        const V2* Dc_ptr    = Dc.ptr();

        Integer res_ld      = res.ld();
        Integer A_ld        = A.ld();
        Integer Dc_inc      = (Dc.cols() == M)? Dc.ld() : 1;

        for (Integer j = 0; j < N; ++j)
        {
            V1 dc           = Dc_ptr[j * Dc_inc];

            for (Integer i = 0; i < M; ++i)
                ptr_res[i]  = ptr_A[i] * dc;

            ptr_res         += res_ld;
            ptr_A           += A_ld;
        };

        struct_flag sf      = res.get_struct();
        sf                  = predefined_struct_ext::get_scal_rowcol(sf);
        res.set_struct(sf);

        ret                 = matcl::Matrix(res,true);
        return;
    };

    static void eval_rowscols(matcl::Matrix& ret, const MP& A, const Mat_D& Dr, const Mat_D& Dc)
    {
        Integer M           = A.rows();
        Integer N           = A.cols();

        MP res              = A.make_unique();
        
        V1* ptr_res         = res.ptr();
        const V1* ptr_A     = A.ptr();
        const V2* Dr_ptr    = Dr.ptr();
        const V2* Dc_ptr    = Dc.ptr();

        Integer res_ld      = res.ld();
        Integer A_ld        = A.ld();
        Integer Dr_inc      = (Dr.cols() == M)? Dr.ld() : 1;
        Integer Dc_inc      = (Dc.cols() == M)? Dc.ld() : 1;

        if (Dr_inc == 1)
        {
            for (Integer j = 0; j < N; ++j)
            {
                V2 dc           = Dc_ptr[j * Dc_inc];

                for (Integer i = 0; i < M; ++i)
                    ptr_res[i]  = ptr_A[i] * Dr_ptr[i] * dc;

                ptr_res         += res_ld;
                ptr_A           += A_ld;
            };
        }
        else
        {
            for (Integer j = 0; j < N; ++j)
            {
                V2 dc           = Dc_ptr[j * Dc_inc];

                for (Integer i = 0; i < M; ++i)
                    ptr_res[i]  = ptr_A[i] * Dr_ptr[i*Dr_inc] * dc;

                ptr_res         += res_ld;
                ptr_A           += A_ld;
            };
        };

        struct_flag sf      = res.get_struct();
        sf                  = predefined_struct_ext::get_scal_rowcol(sf);
        res.set_struct(sf);

        ret                 = matcl::Matrix(res,true);
        return;
    };
};

template<class V1, class V2>
struct scaling_helper_impl<V1,V2,struct_banded>
{
    using MP    = raw::Matrix<V1,struct_banded>;
    using Mat_D = raw::Matrix<V2,struct_dense>;

    static void eval_rows(matcl::Matrix& ret, const MP& A, const Mat_D& Dr)
    {
        Integer M           = A.rows();
        Integer N           = A.cols();

        MP res              = A.make_unique();       

        Integer res_ld      = res.ld();
        Integer A_ld        = A.ld();
        Integer Dr_inc      = (Dr.cols() == M)? Dr.ld() : 1;

        Integer fda         = A.first_diag();

        if (fda == A.last_diag())
        {
            V1* ptr_res     = res.rep_ptr() + res.first_elem_diag(fda);
            const V1* ptr_A = A.rep_ptr() + A.first_elem_diag(fda);
            Integer s       = A.diag_length(fda);
            Integer fr      = A.first_row_on_diag(fda);

            const V2* Dr_ptr = Dr.ptr() + fr * Dr_inc;

            for (Integer i = 0; i < s; ++i)
            {
                ptr_res[0]  = ptr_A[0] * Dr_ptr[0];

                ptr_res     += res_ld;
                ptr_A       += A_ld;
                Dr_ptr      += Dr_inc;
            };

            struct_flag sf      = res.get_struct();
            sf                  = predefined_struct_ext::get_scal_rowcol(sf);
            res.set_struct(sf);

            ret             = matcl::Matrix(res,true);
            return;
        };

        V1* ptr_res         = res.rep_ptr();
        const V1* ptr_A     = A.rep_ptr();
        const V2* Dr_ptr    = Dr.ptr();

        if (Dr_inc == 1)
        {
            for (Integer j = 0; j < N; ++j)
            {
                Integer fr      = A.first_row(j);
                Integer lr      = A.last_row(j);
                Integer pa      = A.first_elem_pos(j);
                Integer pr      = res.first_elem_pos(j);

                for (Integer i = fr; i <= lr; ++i, ++pr, ++pa)
                    ptr_res[pr] = ptr_A[pa] * Dr_ptr[i];

                ptr_res         += res_ld;
                ptr_A           += A_ld;
            };
        }
        else
        {
            for (Integer j = 0; j < N; ++j)
            {
                Integer fr      = A.first_row(j);
                Integer lr      = A.last_row(j);
                Integer pa      = A.first_elem_pos(j);
                Integer pr      = res.first_elem_pos(j);

                for (Integer i = fr; i <= lr; ++i, ++pr, ++pa)
                    ptr_res[pr] = ptr_A[pa] * Dr_ptr[i*Dr_inc];

                ptr_res         += res_ld;
                ptr_A           += A_ld;
            };
        }

        struct_flag sf      = res.get_struct();
        sf                  = predefined_struct_ext::get_scal_rowcol(sf);
        res.set_struct(sf);

        ret                 = matcl::Matrix(res,true);
        return;
    };

    static void eval_cols(matcl::Matrix& ret, const MP& A, const Mat_D& Dc)
    {
        Integer M           = A.rows();
        Integer N           = A.cols();

        MP res              = A.make_unique();       

        Integer res_ld      = res.ld();
        Integer A_ld        = A.ld();
        Integer Dc_inc      = (Dc.cols() == M)? Dc.ld() : 1;

        Integer fda         = A.first_diag();

        if (fda == A.last_diag())
        {
            V1* ptr_res     = res.rep_ptr() + res.first_elem_diag(fda);
            const V1* ptr_A = A.rep_ptr() + A.first_elem_diag(fda);
            Integer s       = A.diag_length(fda);
            Integer fc      = A.first_col_on_diag(fda);

            const V2* Dc_ptr = Dc.ptr() + fc * Dc_inc;

            for (Integer i = 0; i < s; ++i)
            {
                ptr_res[0]  = ptr_A[0] * Dc_ptr[0];

                ptr_res     += res_ld;
                ptr_A       += A_ld;
                Dc_ptr      += Dc_inc;
            };

            struct_flag sf      = res.get_struct();
            sf                  = predefined_struct_ext::get_scal_rowcol(sf);
            res.set_struct(sf);

            ret             = matcl::Matrix(res,true);
            return;
        };

        V1* ptr_res         = res.rep_ptr();
        const V1* ptr_A     = A.rep_ptr();
        const V2* Dc_ptr    = Dc.ptr();

        for (Integer j = 0; j < N; ++j)
        {
            V1 dc           = Dc_ptr[j * Dc_inc];

            Integer fr      = A.first_row(j);
            Integer lr      = A.last_row(j);
            Integer pa      = A.first_elem_pos(j);
            Integer pr      = res.first_elem_pos(j);

            for (Integer i = fr; i <= lr; ++i, ++pr, ++pa)
                ptr_res[pr] = ptr_A[pa] * dc;

            ptr_res         += res_ld;
            ptr_A           += A_ld;
        };

        struct_flag sf      = res.get_struct();
        sf                  = predefined_struct_ext::get_scal_rowcol(sf);
        res.set_struct(sf);

        ret                 = matcl::Matrix(res,true);
        return;
    };

    static void eval_rowscols(matcl::Matrix& ret, const MP& A, const Mat_D& Dr, const Mat_D& Dc)
    {
        Integer M           = A.rows();
        Integer N           = A.cols();

        MP res              = A.make_unique();       

        Integer res_ld      = res.ld();
        Integer A_ld        = A.ld();
        Integer Dr_inc      = (Dr.cols() == M)? Dr.ld() : 1;
        Integer Dc_inc      = (Dc.cols() == M)? Dc.ld() : 1;

        Integer fda         = A.first_diag();

        if (fda == A.last_diag())
        {
            V1* ptr_res     = res.rep_ptr() + res.first_elem_diag(fda);
            const V1* ptr_A = A.rep_ptr() + A.first_elem_diag(fda);
            Integer s       = A.diag_length(fda);
            Integer fr      = A.first_row_on_diag(fda);
            Integer fc      = A.first_col_on_diag(fda);

            const V2* Dr_ptr = Dr.ptr() + fr * Dr_inc;
            const V2* Dc_ptr = Dc.ptr() + fc * Dc_inc;

            for (Integer i = 0; i < s; ++i)
            {
                ptr_res[0]  = ptr_A[0] * Dr_ptr[0] * Dc_ptr[0];

                ptr_res     += res_ld;
                ptr_A       += A_ld;
                Dc_ptr      += Dc_inc;
                Dr_ptr      += Dr_inc;
            };

            struct_flag sf      = res.get_struct();
            sf                  = predefined_struct_ext::get_scal_rowcol(sf);
            res.set_struct(sf);

            ret             = matcl::Matrix(res,true);
            return;
        };

        V1* ptr_res         = res.rep_ptr();
        const V1* ptr_A     = A.rep_ptr();
        const V2* Dr_ptr    = Dr.ptr();
        const V2* Dc_ptr    = Dc.ptr();

        if (Dr_inc == 1)
        {
            for (Integer j = 0; j < N; ++j)
            {
                V2 dc           = Dc_ptr[j * Dc_inc];

                Integer fr      = A.first_row(j);
                Integer lr      = A.last_row(j);
                Integer pa      = A.first_elem_pos(j);
                Integer pr      = res.first_elem_pos(j);

                for (Integer i = fr; i <= lr; ++i, ++pr, ++pa)
                    ptr_res[pr] = ptr_A[pa] * Dr_ptr[i] * dc;

                ptr_res         += res_ld;
                ptr_A           += A_ld;
            };
        }
        else
        {
            for (Integer j = 0; j < N; ++j)
            {
                V2 dc           = Dc_ptr[j * Dc_inc];

                Integer fr      = A.first_row(j);
                Integer lr      = A.last_row(j);
                Integer pa      = A.first_elem_pos(j);
                Integer pr      = res.first_elem_pos(j);

                for (Integer i = fr; i <= lr; ++i, ++pr, ++pa)
                    ptr_res[pr] = ptr_A[pa] * Dr_ptr[i*Dr_inc] * dc;

                ptr_res         += res_ld;
                ptr_A           += A_ld;
            };
        }

        struct_flag sf      = res.get_struct();
        sf                  = predefined_struct_ext::get_scal_rowcol(sf);
        res.set_struct(sf);

        ret                 = matcl::Matrix(res,true);
        return;
    };
};

template<class V1, class V2>
struct scaling_helper_impl<V1,V2,struct_sparse>
{
    using MP    = raw::Matrix<V1,struct_sparse>;
    using Mat_D = raw::Matrix<V2,struct_dense>;

    static void eval_rows(matcl::Matrix& ret, const MP& A, const Mat_D& Dr)
    {
        Integer M           = A.rows();
        Integer N           = A.cols();

        if (A.is_unique())
            return eval_rows_inpl(ret, A, Dr);
        
        const sparse_ccs<V1>& Ad = A.rep();
        const Integer * A_c = Ad.ptr_c();
        const Integer * A_r = Ad.ptr_r();
        const V1* A_x       = Ad.ptr_x();

        MP res(A.get_type(), M, N, A.nnz());

        if (A.nnz() == 0)
        {
            ret             = matcl::Matrix(res,true);
            return;
        };

        sparse_ccs<V1>& resd = res.rep();
        Integer * d_c       = resd.ptr_c();
        Integer * d_r       = resd.ptr_r();
        V1* d_x             = resd.ptr_x();

        const V2* Dr_ptr    = Dr.ptr();

        Integer Dr_inc      = (Dr.cols() == M)? Dr.ld() : 1;
        Integer nz          = 0;

        for (Integer j = 0; j < N; ++j)
        {
            d_c[j]          = nz;

            for (Integer k = A_c[j]; k < A_c[j+1]; ++k)
            {
                Integer r   = A_r[k];
                V1 x        = Dr_ptr[r * Dr_inc] * A_x[k];

                d_x[nz]     = x;
                d_r[nz]     = r;
                ++nz;
            };
        };
        d_c[N]              = nz;

        struct_flag sf      = res.get_struct();
        sf                  = predefined_struct_ext::get_scal_rowcol(sf);
        res.set_struct(sf);

        ret                 = matcl::Matrix(res,true);
        return;
    };

    static void eval_cols(matcl::Matrix& ret, const MP& A, const Mat_D& Dc)
    {
        Integer M           = A.rows();
        Integer N           = A.cols();

        if (A.is_unique())
            return eval_cols_inpl(ret, A, Dc);
        
        const sparse_ccs<V1>& Ad = A.rep();
        const Integer * A_c = Ad.ptr_c();
        const Integer * A_r = Ad.ptr_r();
        const V1* A_x       = Ad.ptr_x();

        MP res(A.get_type(), M, N, A.nnz());

        if (A.nnz() == 0)
        {
            ret             = matcl::Matrix(res,true);
            return;
        };

        sparse_ccs<V1>& resd = res.rep();
        Integer * d_c       = resd.ptr_c();
        Integer * d_r       = resd.ptr_r();
        V1* d_x             = resd.ptr_x();

        const V2* Dc_ptr    = Dc.ptr();

        Integer Dc_inc      = (Dc.cols() == M)? Dc.ld() : 1;
        Integer nz          = 0;

        for (Integer j = 0; j < N; ++j)
        {
            d_c[j]          = nz;

            V2 dc           = Dc_ptr[j * Dc_inc];

            for (Integer k = A_c[j]; k < A_c[j+1]; ++k)
            {
                Integer r   = A_r[k];
                V1 x        = A_x[k] * dc;

                d_x[nz]     = x;
                d_r[nz]     = r;
                ++nz;
            };
        };
        d_c[N]              = nz;

        struct_flag sf      = res.get_struct();
        sf                  = predefined_struct_ext::get_scal_rowcol(sf);
        res.set_struct(sf);

        ret                 = matcl::Matrix(res,true);
        return;
    };

    static void eval_rowscols(matcl::Matrix& ret, const MP& A, const Mat_D& Dr, const Mat_D& Dc)
    {
        Integer M           = A.rows();
        Integer N           = A.cols();

        if (A.is_unique())
            return eval_rowscols_inpl(ret, A, Dr, Dc);
        
        const sparse_ccs<V1>& Ad = A.rep();
        const Integer * A_c = Ad.ptr_c();
        const Integer * A_r = Ad.ptr_r();
        const V1* A_x       = Ad.ptr_x();

        MP res(A.get_type(), M, N, A.nnz());

        if (A.nnz() == 0)
        {
            ret             = matcl::Matrix(res,true);
            return;
        };

        sparse_ccs<V1>& resd = res.rep();
        Integer * d_c       = resd.ptr_c();
        Integer * d_r       = resd.ptr_r();
        V1* d_x             = resd.ptr_x();

        const V2* Dr_ptr    = Dr.ptr();
        const V2* Dc_ptr    = Dc.ptr();

        Integer Dr_inc      = (Dr.cols() == M)? Dr.ld() : 1;
        Integer Dc_inc      = (Dc.cols() == M)? Dc.ld() : 1;
        Integer nz          = 0;

        for (Integer j = 0; j < N; ++j)
        {
            d_c[j]          = nz;

            V2 dc           = Dc_ptr[j * Dc_inc];

            for (Integer k = A_c[j]; k < A_c[j+1]; ++k)
            {
                Integer r   = A_r[k];
                V1 x        = Dr_ptr[r * Dr_inc] * A_x[k] * dc;

                d_x[nz]     = x;
                d_r[nz]     = r;
                ++nz;
            };
        };
        d_c[N]              = nz;

        struct_flag sf      = res.get_struct();
        sf                  = predefined_struct_ext::get_scal_rowcol(sf);
        res.set_struct(sf);

        ret                 = matcl::Matrix(res,true);
        return;
    };

    static void eval_rowscols_inpl(matcl::Matrix& ret, const MP& A, const Mat_D& Dr, const Mat_D& Dc)
    {
        Integer M           = A.rows();
        Integer N           = A.cols();

        //A is unique
        sparse_ccs<V1> Ad   = A.rep();
        const Integer * A_c = Ad.ptr_c();
        const Integer * A_r = Ad.ptr_r();
        V1* A_x             = Ad.ptr_x();

        const V2* Dr_ptr    = Dr.ptr();
        const V2* Dc_ptr    = Dc.ptr();

        Integer Dr_inc      = (Dr.cols() == M)? Dr.ld() : 1;
        Integer Dc_inc      = (Dc.cols() == M)? Dc.ld() : 1;

        for (Integer j = 0; j < N; ++j)
        {
            V2 dc           = Dc_ptr[j * Dc_inc];

            for (Integer k = A_c[j]; k < A_c[j+1]; ++k)
            {
                Integer r   = A_r[k];
                V1 x        = Dr_ptr[r * Dr_inc] * A_x[k] * dc;
                A_x[k]      = x;
            };
        };

        struct_flag sf      = A.get_struct();
        sf                  = predefined_struct_ext::get_scal_rowcol(sf);
        A.set_struct(sf);

        ret                 = matcl::Matrix(A,false);
        return;
    };

    static void eval_rows_inpl(matcl::Matrix& ret, const MP& A, const Mat_D& Dr)
    {
        Integer M           = A.rows();
        Integer N           = A.cols();

        //A is unique
        sparse_ccs<V1> Ad   = A.rep();
        const Integer * A_c = Ad.ptr_c();
        const Integer * A_r = Ad.ptr_r();
        V1* A_x             = Ad.ptr_x();

        const V2* Dr_ptr    = Dr.ptr();

        Integer Dr_inc      = (Dr.cols() == M)? Dr.ld() : 1;

        for (Integer j = 0; j < N; ++j)
        {
            for (Integer k = A_c[j]; k < A_c[j+1]; ++k)
            {
                Integer r   = A_r[k];
                V1 x        = Dr_ptr[r * Dr_inc] * A_x[k];
                A_x[k]      = x;
            };
        };

        struct_flag sf      = A.get_struct();
        sf                  = predefined_struct_ext::get_scal_rowcol(sf);
        A.set_struct(sf);

        ret                 = matcl::Matrix(A,false);
        return;
    };

    static void eval_cols_inpl(matcl::Matrix& ret, const MP& A, const Mat_D& Dc)
    {
        Integer M           = A.rows();
        Integer N           = A.cols();

        //A is unique
        sparse_ccs<V1> Ad   = A.rep();
        const Integer * A_c = Ad.ptr_c();
        V1* A_x             = Ad.ptr_x();

        const V2* Dc_ptr    = Dc.ptr();
        Integer Dc_inc      = (Dc.cols() == M)? Dc.ld() : 1;

        for (Integer j = 0; j < N; ++j)
        {
            V2 dc           = Dc_ptr[j * Dc_inc];

            for (Integer k = A_c[j]; k < A_c[j+1]; ++k)
            {
                V1 x        = A_x[k] * dc;
                A_x[k]      = x;
            };
        };

        struct_flag sf      = A.get_struct();
        sf                  = predefined_struct_ext::get_scal_rowcol(sf);
        A.set_struct(sf);

        ret                 = matcl::Matrix(A,false);
        return;
    };

};

template<class MT>
void scaling_helper<MT>::eval_rows(matcl::Matrix& ret, const Mat& m, const Mat_D& Dr)
{
    return scaling_helper_impl<V,V,S>::eval_rows(ret, m, Dr);
};

template<class MT>
void scaling_helper<MT>::eval_cols(matcl::Matrix& ret, const Mat& m, const Mat_D& Dc)
{
    return scaling_helper_impl<V,V,S>::eval_cols(ret, m, Dc);
};

template<class MT>
void scaling_helper<MT>::eval_rowscols(matcl::Matrix& ret, const Mat& m, const Mat_D& Dr, const Mat_D& Dc)
{
    return scaling_helper_impl<V,V,S>::eval_rowscols(ret, m, Dr, Dc);
};

template<class MT>
void scaling_real_helper<MT>::eval_rows(matcl::Matrix& ret, const Mat& m, const Mat_D& Dr)
{
    return scaling_helper_impl<V,VR,S>::eval_rows(ret, m, Dr);
};

template<class MT>
void scaling_real_helper<MT>::eval_cols(matcl::Matrix& ret, const Mat& m, const Mat_D& Dc)
{
    return scaling_helper_impl<V,VR,S>::eval_cols(ret, m, Dc);
};

template<class MT>
void scaling_real_helper<MT>::eval_rowscols(matcl::Matrix& ret, const Mat& m, const Mat_D& Dr, const Mat_D& Dc)
{
    return scaling_helper_impl<V,VR,S>::eval_rowscols(ret, m, Dr, Dc);
};

};};};

MACRO_INSTANTIATE_G_1(matcl::raw::details::scaling_helper)
MACRO_INSTANTIATE_S_1(matcl::raw::details::scaling_helper)

MACRO_INSTANTIATE_G_1(matcl::raw::details::scaling_real_helper)
MACRO_INSTANTIATE_S_1(matcl::raw::details::scaling_real_helper)
