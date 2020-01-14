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

#pragma once

#include "matcl-core/error/exception_classes.h"
#include "matcl-matrep/func/raw/raw_manip.h"
#include <vector>
#include "matcl-internals/base/optim_params.h"
#include "matcl-matrep/algs/dense_algs.h"
#include "matcl-matrep/algs/sparse_algs.h"
#include "matcl-matrep/algs/banded_algs.h"
#include "matcl-internals/algs/scatter.h"

namespace matcl { namespace details
{

namespace md = matcl::details;

template<class V>
struct get_matrix_ti
{};

template<>
struct get_matrix_ti<Integer>
{
    using V = Integer;

    static ti::ti_empty eval(const V&)                          { return ti::ti_empty(); };
    static ti::ti_empty eval(const Matrix& )                    { return ti::ti_empty(); };
    template<class S>
    static ti::ti_empty eval(const raw::Matrix<V,S>&)           { return ti::ti_empty(); };
};

template<>
struct get_matrix_ti<Real>
{
    using V = Real;
    static ti::ti_empty eval(const V&)                          { return ti::ti_empty(); };
    static ti::ti_empty eval(const Matrix& )                    { return ti::ti_empty(); };
    template<class S>
    static ti::ti_empty eval(const raw::Matrix<V,S>&)           { return ti::ti_empty(); };
};

template<>
struct get_matrix_ti<Float>
{
    using V = Float;
    static ti::ti_empty eval(const V&)                          { return ti::ti_empty(); };
    static ti::ti_empty eval(const Matrix& )                    { return ti::ti_empty(); };
    template<class S>
    static ti::ti_empty eval(const raw::Matrix<V,S>&)           { return ti::ti_empty(); };
};

template<>
struct get_matrix_ti<Complex>
{
    using V = Complex;
    static ti::ti_empty eval(const V&)                          { return ti::ti_empty(); };
    static ti::ti_empty eval(const Matrix& )                    { return ti::ti_empty(); };
    template<class S>
    static ti::ti_empty eval(const raw::Matrix<V,S>&)           { return ti::ti_empty(); };
};

template<>
struct get_matrix_ti<Float_complex>
{
    using V = Float_complex;
    static ti::ti_empty eval(const V&)                          { return ti::ti_empty(); };
    static ti::ti_empty eval(const Matrix& )                    { return ti::ti_empty(); };
    template<class S>
    static ti::ti_empty eval(const raw::Matrix<V,S>&)           { return ti::ti_empty(); };
};

template<>
struct get_matrix_ti<Object>
{
    static ti::ti_object eval(const Object& A)                  { return A.get_type(); };
    static ti::ti_object eval(const Matrix& A)                  { return A.get_type(); };
    template<class S>
    static ti::ti_object eval(const raw::Matrix<Object,S>& A)   { return A.get_type(); };
};

template<class T_ret, class T1, class T2>
struct link_ret_ti
{
    using ti_ret_t  = typename ti::get_ti_type<T_ret>::type;
    using ti_1_t    = typename ti::get_ti_type<T1>::type;
    using ti_2_t    = typename ti::get_ti_type<T2>::type;

    static ti_ret_t eval(ti_1_t t1, const T2& v2)
    { 
        return ti::unify_ti_assign<ti_ret_t,ti_1_t,ti_2_t>(t1,ti::get_ti(v2));
    }

    static ti_ret_t eval(ti_1_t t1, ti_2_t t2)
    { 
        return ti::unify_ti_assign<ti_ret_t,ti_1_t,ti_2_t>(t1,t2);
    }
};

template<class M1,class M2, class val_type,bool req_promo>
struct assign_mat_impl
{};

template<class M1,class M2,class struct_type_2>
struct assign_mat_DD
{};

template<class M1,class M2,class str_1,class str_2>
struct assign_mat_str
{};

template<class M1,class M2>
struct assign_mat_str<M1,M2,struct_dense,struct_sparse>
{
    using FullMatrix = raw::Matrix<typename M1::value_type,struct_dense>;

    static Matrix& eval(Matrix& A,const M2& B,const colon& c1)
    {
        const M1& mat_tmp = A.get_impl<M1>();

        colon_info ci;
        make_index(mat_tmp.rows(),mat_tmp.cols(),c1,ci);        

        Integer r   = B.rows();
        Integer c   = B.cols();
        Integer s   = ci.rows();

        if (ci.r_flag == 1)
            error::check_assign_1(s, r, c);
        else
            error::check_assign_1(ci.get_rim_1().size(), r, c);
        
        if (s == 0)
            return A;

        FullMatrix& mat     = A.get_impl_unique<FullMatrix>();
        using val_type_1    = typename M1::value_type;
        using val_type_2    = typename M2::value_type;

        val_type_1 Z = default_value<val_type_1>(get_matrix_ti<val_type_1>::eval(A));
        val_type_1* ptr_mat = mat.ptr();
        Integer mat_ld      = mat.ld();

        if (ci.r_flag == 0)
        {
        }
        else
        {
            if (mat.ld() == mat.rows())
            {
                for (Integer i = 0, k = ci.r_start-1; i < ci.r_size; ++i, k+= ci.r_step)
                    mrd::assign_helper(ptr_mat[k],Z);
            }
            else
            {
                for (Integer i = 0, k = ci.r_start-1; i < ci.r_size; ++i, k+= ci.r_step)
                {
                    while(k >= r)
                    {
                        k       -= r;
                        ptr_mat += mat_ld;
                    };

                    mrd::assign_helper(ptr_mat[k],Z);
                };
            };
        };

        const raw::details::sparse_ccs<val_type_2>& Bd = B.rep();

        const Integer * d_c		= Bd.ptr_c();
        const Integer * d_r		= Bd.ptr_r();
        const val_type_2 * d_x	= Bd.ptr_x();        

        if (ci.r_flag == 0)
        {
            using FullMatrix2   = raw::Matrix<val_type_2,struct_dense>;
            FullMatrix2 BF      = raw::full(B);
            const val_type_2* ptr_BF = BF.ptr();

            bool single         = ci.is_double_mat_colon() == false;

            Integer BFr         = BF.rows(); 
            Integer BFc         = BF.cols();
            Integer BFld        = BF.ld();
            Integer matr        = mat.rows();

            if (single == true)
            {   
                const mr::integer_dense& ci_ri  = ci.get_rim_1();
                const Integer* ptr_ri           = ci_ri.ptr();

                for (Integer j = 0, k = 0; j < BFc; ++j) 
                {
                    for (Integer i = 0; i < BFr; ++i, ++k)
                    {
                        Integer pos = ptr_ri[k];

                        Integer mr, mc;
                        pos2ind(pos,matr,mr,mc);

                        mrd::assign_helper(ptr_mat[mr+mc*mat_ld],val_type_1(ptr_BF[i]));
                    };
                    ptr_BF += BFld;
                };
            }
            else
            {
                const mr::integer_dense& rim    = ci.get_rim_r();
                const mr::integer_dense& cim    = ci.get_rim_c();

                const Integer* ptr_ri           = rim.ptr();
                const Integer* ptr_ci           = cim.ptr();

                for (Integer j = 0, k = 0; j < BFc; ++j) 
                {
                    for (Integer i = 0; i < BFr; ++i, ++k)
                    {
                        Integer mr  = ptr_ri[k] - 1;
                        Integer mc  = ptr_ci[k] - 1;

                        mrd::assign_helper(ptr_mat[mr+mc*mat_ld],val_type_1(ptr_BF[i]));
                    };

                    ptr_BF += BFld;
                };
            };
        }
        else
        {            
            Integer pos                 = ci.r_start - 1;            
            Integer Br                  = B.rows();
            Integer Bc                  = B.cols();

            if (mat.ld() == mat.rows())
            {
                for (Integer j = 0; j < Bc; ++j) 
                {
                    Integer row, row_l = 0;

                    for (Integer k = d_c[j]; k < d_c[j+1]; ++k)
                    {
                        row		        = d_r[k];
                        pos             += imult(row - row_l,ci.r_step);
                        row_l           = row;
                        mrd::assign_helper(ptr_mat[pos],val_type_1(d_x[k]));
                    };
                    pos                 += imult(Br - row_l,ci.r_step);                
                };
            }
            else
            {
                Integer matr            = mat.rows();

                for (Integer j = 0; j < Bc; ++j) 
                {
                    Integer row, row_l = 0;

                    for (Integer k = d_c[j]; k < d_c[j+1]; ++k)
                    {
                        row		        = d_r[k];
                        pos             += imult(row - row_l,ci.r_step);
                        row_l           = row;

                        while(pos >= matr)
                        {
                            pos         -= matr;
                            ptr_mat     += mat_ld;
                        };

                        mrd::assign_helper(ptr_mat[pos],val_type_1(d_x[k]));
                    };
                    pos                 += imult(Bc - row_l,ci.r_step);                
                };
            };
        };

        mat.get_struct().reset();
        return A;
    };

    static Matrix& eval(Matrix& A,const M2& B,const colon& c1,const colon& c2)
    {	
        const M1& mat_tmp = A.get_impl<M1>();

        colon_info ci;
        make_index(mat_tmp.rows(),mat_tmp.cols(),c1,c2,ci);

        Integer nr = ci.rows(), nc = ci.cols();

        error::check_assign_2(nr, nc, B.rows(), B.cols());

        if (nr == 0 || nc == 0)
            return A;

        FullMatrix& mat = A.get_impl_unique<FullMatrix>();
        typename M1::value_type* ptr_mat = mat.ptr();

        using value_type    = typename M2::value_type;
        using value_type_1  = typename M1::value_type;
        const raw::details::sparse_ccs<value_type>& Bd = B.rep();

        value_type_1 Z = default_value<value_type_1>(get_matrix_ti<value_type_1>::eval(A));

        const Integer * d_c		= Bd.ptr_c();
        const Integer * d_r		= Bd.ptr_r();
        const value_type * d_x	= Bd.ptr_x();

        const Integer* ptr_ri   = ci.get_ri_2_ptr();
        const Integer* ptr_ci   = ci.get_ci_2_ptr();
        Integer Bc              = B.cols();
        Integer Br              = B.rows();
        Integer mat_ld          = mat.ld();
        Integer cir             = ci.rows();

        if (ci.r_flag == 0 && ci.c_flag == 0)
        {
            using val_type_1    = typename M1::value_type;
            using scatter       = matcl::algorithm::scatter;

            workspace2<val_type_1>   x_work(get_matrix_ti<val_type_1>::eval(A),B.rows());            

            scatter sc = scatter::get(Br, Bc);

            for (Integer j = 0; j < Bc; ++j) 
            {
                Integer kf			= d_c[j];
                Integer kl			= d_c[j+1];

                Integer pos_c		= imult(ptr_ci[j]-1,mat_ld);
                ptr_mat             = mat.ptr() + pos_c;
                auto mark           = sc.next_mark();

                for (Integer k = kf; k < kl; ++k)
                {
                    Integer row		= d_r[k];
                    mrd::assign_helper(x_work[row],val_type_1(d_x[k]));
                    sc[row]		    = mark;
                };

                Integer ci_ri_size  = cir;
                for (Integer k = 0; k < ci_ri_size; ++k)
                {
                    if (sc[k] < mark)
                        mrd::assign_helper(ptr_mat[ptr_ri[k]-1],Z);
                    else
                        mrd::assign_helper(ptr_mat[ptr_ri[k]-1],x_work[k]);
                };
            };
        }
        else if (ci.r_flag == 0 && ci.c_flag == 1)
        {
            using val_type_1    = typename M1::value_type;
            using scatter       = matcl::algorithm::scatter;

            workspace2<val_type_1>   x_work(get_matrix_ti<val_type_1>::eval(A),B.rows());

            scatter sc = scatter::get(Br,Bc);

            Integer pos_c			= imult(ci.c_start-1,mat_ld);
            Integer dpos_c			= imult(ci.c_step,mat_ld);

            ptr_mat                 = mat.ptr() + pos_c;

            for (Integer j = 0; j < Bc; ++j) 
            {				
                auto mark           = sc.next_mark();

                for (Integer k = d_c[j]; k < d_c[j+1]; ++k)
                {
                    Integer row		= d_r[k];
                    mrd::assign_helper(x_work[row],val_type_1(d_x[k]));
                    sc[row]		    = mark;
                };
                Integer ci_ri_size  = cir;

                for (Integer k = 0; k < ci_ri_size; ++k)
                {
                    if (sc[k] < mark)
                        mrd::assign_helper(ptr_mat[ptr_ri[k]-1],Z);
                    else
                        mrd::assign_helper(ptr_mat[ptr_ri[k]-1],x_work[k]);
                };

                ptr_mat += dpos_c;
            };
        }
        else if (ci.r_flag == 1 && ci.c_flag == 0)
        {
            for (Integer j = 0; j < Bc; ++j) 
            {
                Integer kf			= d_c[j];
                Integer kl			= d_c[j+1];

                Integer pos_c		= imult(ptr_ci[j]-1,mat_ld);
                ptr_mat             = mat.ptr() + pos_c;

                for (Integer k = 0, pos = ci.r_start-1; k < ci.r_size; ++k, pos+=ci.r_step)
                    mrd::assign_helper(ptr_mat[pos],Z);

                if (kf == kl)
                    continue;

                for (Integer k = d_c[j]; k < d_c[j+1]; ++k)
                {
                    Integer row		= d_r[k];
                    mrd::assign_helper(ptr_mat[ci.row_index_2(row+1)-1],value_type_1(d_x[k]));
                };
            };
        }
        else if (ci.r_flag == 1 && ci.c_flag == 1)
        {
            algorithm::dense_change_entries(mat,ci,Z);

            Integer pos_c			= imult(ci.c_start-1,mat_ld);
            Integer dpos_c			= imult(ci.c_step,mat_ld);
            ptr_mat                 = mat.ptr() + pos_c;

            for (Integer j = 0; j < Bc; ++j) 
            {				
                for (Integer k = d_c[j]; k < d_c[j+1]; ++k)
                {
                    Integer row		= d_r[k];
                    mrd::assign_helper(ptr_mat[ci.row_index_2(row+1)-1],value_type_1(d_x[k]));
                };
                ptr_mat += dpos_c;
            };
        };

        mat.get_struct().reset();
        return A;
    };

    static Matrix& eval_diag(Matrix& A,const M2& B,Integer d)
    {
        using value_type    = typename M1::value_type;
        using FullMatrix    = raw::Matrix<value_type,struct_dense>;
        FullMatrix BC = mr::converter<FullMatrix,M2>::eval(B);
        algorithm::dense_change_diag<value_type>(A.get_impl_unique<M1>(),d,BC);
        return A;
    };
};

template<class M1,class M2>
struct assign_mat_str<M1,M2,struct_sparse,struct_sparse>
{
    static Matrix& eval(Matrix& A,const M2& B,const colon& c1)
    {
        M1& mat = A.get_impl_unique<M1>();

        colon_info ci;
        make_index(mat.rows(),mat.cols(),c1,ci);		

        Integer s = ci.rows();        

        if (ci.r_flag == 1)
            error::check_assign_1(s, B.rows(), B.cols());
        else
            error::check_assign_1(ci.get_rim_1().size(), B.rows(), B.cols());
        
        if (s == 0)
            return A;

        if (B.nnz() == 0)
        {
            algorithm::zero_entries_2(A, mat, ci);
            return A;
        };        

        Real density = (Real(mat.nnz()) + Real(B.nnz())) / (mat.rows()+1.) / (mat.cols()+1.);

        if (density > optim_params::max_sparse_density_max)
        {
            using value_type    = typename M1::value_type;
            using matrix_type   = raw::Matrix<value_type,struct_dense>;
            A = matcl::full(Matrix(mat,false));
            return assign_mat_str<matrix_type,M2,struct_dense,struct_sparse>::eval(A,B,c1);
        };

        const M1& mat2 = A.get_impl<M1>();

        Matrix ret;
        algorithm::change_submatrix_2(ret, mat2, ci, B);

        A = ret;
        return A;
    };

    static Matrix& eval(Matrix& A,const M2& B,const colon& c1,const colon& c2)
    {
        M1& mat = A.get_impl_unique<M1>();

        colon_info ci;
        make_index(mat.rows(),mat.cols(),c1,c2,ci);

        Integer nr = ci.rows(), nc = ci.cols();

        error::check_assign_2(nr, nc, B.rows(), B.cols());

        if (nr == 0 || nc == 0)
            return A;

        if (B.nnz() == 0)
        {
            algorithm::zero_entries(A, mat, ci);
            return A;
        };

        Real density = (Real(mat.nnz()) + Real(B.nnz())) / (mat.rows()+1.) / (mat.cols()+1.);
        if (density > optim_params::max_sparse_density_max)
        {
            using value_type    = typename M1::value_type;
            using matrix_type   = raw::Matrix<value_type,struct_dense>;
            A = matcl::full(Matrix(mat,false));
            return assign_mat_str<matrix_type,M2,struct_dense,struct_sparse>::eval(A,B,c1,c2);
        };

        const M1& mat2 = A.get_impl<M1>();

        Matrix ret;
        algorithm::change_submatrix(ret, mat2, ci, B);
        
        A = ret;
        return A;
    };

    static Matrix& eval_diag(Matrix& A,const M2& B,Integer d)
    {
        using value_type    = typename M1::value_type;
        using FullMatrix    = raw::Matrix<value_type,struct_dense>;
        FullMatrix BC = mr::converter<FullMatrix,M2>::eval(B);

        algorithm::sparse_change_diag<value_type>(A, A.get_impl_unique<M1>(), d, BC);
        return A;
    }
};

template<class M1,class M2>
struct assign_mat_str<M1,M2,struct_dense,struct_dense>
{
    static Matrix& eval(Matrix& A,const M2& B,const colon& c1)
    {
        const M1& mat_tmp   = A.get_impl<M1>();
        using val_type_1    = typename M1::value_type;
        using val_type_2    = typename M2::value_type;        

        colon_info ci;
        make_index(mat_tmp.rows(),mat_tmp.cols(),c1,ci);        

        if (ci.r_flag == 1)
            error::check_assign_1(ci.r_size, B.rows(), B.cols());
        else
            error::check_assign_1(ci.get_rim_1().size(), B.rows(), B.cols());      

        if (ci.rows() == 0)
            return A;
        
        M1& mat = A.get_impl_unique<M1>();
        val_type_1* ptr_mat = mat.ptr();

        const M2 B2             = B.make_explicit();
        const val_type_2* ptr_B = B2.ptr();
        Integer matr            = mat.rows();
        Integer mat_ld          = mat.ld();

        Integer s = ci.rows();

        if (ci.r_flag == 0)
        {
            bool single         = ci.is_double_mat_colon() == false;

            if (single == true)
            {
                const mr::integer_dense& ci_ri  = ci.get_rim_1();
                const Integer* ptr_ri           = ci_ri.ptr();

                ptr_B = B2.ptr();

                for (Integer i = 0; i < s; ++i)
                {
                    Integer pos = ptr_ri[i];

                    Integer mr, mc;
                    pos2ind(pos,matr,mr,mc);

                    mrd::assign_helper(ptr_mat[mr+mc*mat_ld],val_type_1(ptr_B[i]));
                };		
            }
            else
            {
                const mr::integer_dense& rim   = ci.get_rim_r();
                const mr::integer_dense& cim   = ci.get_rim_c();
                const Integer* ptr_ri   = rim.ptr();
                const Integer* ptr_ci   = cim.ptr();

                ptr_B = B2.ptr();

                for (Integer i = 0; i < s; ++i)
                {
                    Integer mr  = ptr_ri[i] - 1;
                    Integer mc  = ptr_ci[i] - 1;

                    mrd::assign_helper(ptr_mat[mr+mc*mat_ld],val_type_1(ptr_B[i]));
                };		
            };
        }
        else
        {            
            if (mat_ld == matr)
            {
                for (Integer i = 0, k = ci.r_start-1; i < ci.r_size; ++i, k+= ci.r_step)
                    mrd::assign_helper(ptr_mat[k],val_type_1(ptr_B[i]));
            }
            else
            {
                for (Integer i = 0, k = ci.r_start-1; i < ci.r_size; ++i, k+= ci.r_step)
                {
                    while(k >= matr)
                    {
                        k       -= matr;
                        ptr_mat += mat_ld;
                    };
                    mrd::assign_helper(ptr_mat[k],val_type_1(ptr_B[i]));
                };	
            };
        };

        mat.get_struct().reset();
        return A;
    };

    static Matrix& eval(Matrix& A,const M2& B,const colon& c1,const colon& c2)
    {
        const M1& mat_tmp   = A.get_impl<M1>();
        using val_type_1    = typename M1::value_type;
        using val_type_2    = typename M2::value_type;

        const val_type_2* ptr_B = B.ptr();
        Integer B_ld            = B.ld();

        colon_info ci;		
        make_index(mat_tmp.rows(),mat_tmp.cols(),c1,c2,ci);

        Integer nr = ci.rows(), nc = ci.cols();

        error::check_assign_2(nr, nc, B.rows(), B.cols());

        if (nr == 0 || nc == 0)
            return A;

        M1& mat = A.get_impl_unique<M1>();
        val_type_1* ptr_mat = mat.ptr();
        Integer mat_ld      = mat.ld();

        if (ci.r_flag == 0 && ci.c_flag == 0)
        {
            const mr::integer_dense& ci_ri = ci.get_rim_2();
            const mr::integer_dense& ci_ci = ci.get_cim_2();

            const Integer* ptr_ri = ci_ri.ptr();
            const Integer* ptr_ci = ci_ci.ptr();

            for (Integer j = 0; j < nc; ++j)
            {
                Integer pos = imult(ptr_ci[j]-1,mat_ld);
                ptr_mat = mat.ptr() + pos;

                for (Integer i = 0; i < nr; ++i)
                    mrd::assign_helper(ptr_mat[ptr_ri[i]-1],val_type_1(ptr_B[i]));

                ptr_B += B_ld;
            };
            mat.get_struct().reset();
            return A;
        }
        else if (ci.r_flag == 0 && ci.c_flag == 1)
        {
            Integer pos_A   = imult(ci.c_start - 1,mat.ld());
            Integer dpos_A  = imult(ci.c_step, mat.ld());

            ptr_mat = mat.ptr() + pos_A;

            const mr::integer_dense& ci_ri  = ci.get_rim_2();
            const Integer* ptr_ri           = ci_ri.ptr();

            for (Integer j = 0; j < nc; ++j)
            {
                for (Integer i = 0; i < nr; ++i)
                    mrd::assign_helper(ptr_mat[ptr_ri[i]-1],val_type_1(ptr_B[i]));

                ptr_mat += dpos_A;
                ptr_B   += B_ld;
            };

            mat.get_struct().reset();
            return A;
        }
        else if (ci.r_flag == 1 && ci.c_flag == 0)
        {
            const mr::integer_dense& ci_ci  = ci.get_cim_2();
            const Integer* ptr_ci           = ci_ci.ptr();

            for (Integer j = 0; j < nc; ++j)
            {
                Integer pos_Ac = imult(ptr_ci[j]-1,mat_ld);
                ptr_mat = mat.ptr() + pos_Ac;

                Integer pos = ci.r_start - 1;

                for (Integer i = 0; i < ci.r_size; ++i)
                {
                    mrd::assign_helper(ptr_mat[pos],val_type_1(ptr_B[i]));
                    pos += ci.r_step;
                };
                ptr_B += B_ld;
            };

            mat.get_struct().reset();
            return A;
        }
        else if (ci.r_flag == 1 && ci.c_flag == 1)
        {
            Integer pos_Ac = imult(ci.c_start - 1,mat.ld());
            Integer dpos_Ac = imult(ci.c_step, mat.ld());
            ptr_mat = mat.ptr() + pos_Ac;

            for (Integer j = 0; j < nc; ++j)
            {
                Integer pos = ci.r_start - 1;
                for (Integer i = 0; i < ci.r_size; ++i)
                {
                    mrd::assign_helper(ptr_mat[pos],val_type_1(ptr_B[i]));
                    pos += ci.r_step;
                };

                ptr_B   += B_ld;
                ptr_mat += dpos_Ac;
            };

            mat.get_struct().reset();
            return A;
        };

        matcl_assert(0,"unknown case");
        throw;
    };

    static Matrix& eval_diag(Matrix& A,const M2& B,Integer d)
    {
        using value_type    = typename M1::value_type;
        using FullMatrix    = raw::Matrix<value_type,struct_dense>;
        const FullMatrix& BC= mr::converter<FullMatrix,M2>::eval(B);
        algorithm::dense_change_diag<value_type>(A.get_impl_unique<M1>(),d,BC);
        return A;
    }
};

template<class M1,class M2>
struct assign_mat_str<M1,M2,struct_dense,struct_banded>
{
    static Matrix& eval(Matrix& A,const M2& B,const colon& c1)
    {
        const M1& mat_tmp = A.get_impl<M1>();

        colon_info ci;
        make_index(mat_tmp.rows(),mat_tmp.cols(),c1,ci);        

        Integer r   = B.rows();
        Integer s   = ci.rows();

        if (ci.r_flag == 1)
            error::check_assign_1(s, B.rows(), B.cols());
        else
            error::check_assign_1(ci.get_rim_1().size(), B.rows(), B.cols());     

        if (s == 0)
            return A;

        M1& mat = A.get_impl_unique<M1>();
        using val_type_1    = typename M1::value_type;
        using val_type_2    = typename M2::value_type;
        val_type_1* ptr_mat = mat.ptr();
        Integer matr        = mat.rows();
        Integer mat_ld      = mat.ld();

        if (ci.r_flag == 0)
        {
        }
        else
        {
            val_type_1 Z = default_value<val_type_1>(get_matrix_ti<val_type_1>::eval(A));
            if (mat_ld == matr)
            {
                for (Integer i = 0, k = ci.r_start-1; i < ci.r_size; ++i, k+= ci.r_step)
                    mrd::assign_helper(ptr_mat[k],Z);
            }
            else
            {                
                for (Integer i = 0, k = ci.r_start-1; i < ci.r_size; ++i, k+= ci.r_step)
                {
                    while(k >= matr)
                    {
                        k       -= matr;
                        ptr_mat += mat_ld;
                    };

                    mrd::assign_helper(ptr_mat[k],Z);
                };
            };
        };
        
        Integer Bc              = B.cols();
        Integer Br              = B.rows();
        Integer B_ld            = B.ld();

        if (ci.r_flag == 0)
        {
            bool single         = ci.is_double_mat_colon() == false;

            using FullMatrix    = raw::Matrix<typename M2::value_type,struct_dense>;
            FullMatrix BF       = raw::converter<FullMatrix,M2>::eval(B);
            const typename M2::value_type* ptr_BF = BF.ptr();
            Integer BFld        = BF.ld();

            if (single == true)
            {
                const Integer* ptr_ri   = ci.get_ri_1_ptr();

                for (Integer j = 0, k = 0; j < Bc; ++j)
                {
                    for (Integer i = 0; i < Br; ++i, ++k)
                    {
                        Integer pos = ptr_ri[k];

                        Integer mr, mc;
                        pos2ind(pos,matr,mr,mc);

                        mrd::assign_helper(ptr_mat[mr+mc*mat_ld],val_type_1(ptr_BF[i]));
                    };		
                    ptr_BF += BFld;
                };
            }
            else
            {
                const Integer* ptr_ri   = ci.get_ri_r_ptr();
                const Integer* ptr_ci   = ci.get_ri_c_ptr();

                for (Integer j = 0, k = 0; j < Bc; ++j)
                {
                    for (Integer i = 0; i < Br; ++i, ++k)
                    {
                        Integer mr  = ptr_ri[k] - 1;
                        Integer mc  = ptr_ci[k] - 1;

                        mrd::assign_helper(ptr_mat[mr+mc*mat_ld],val_type_1(ptr_BF[i]));
                    };		

                    ptr_BF += BFld;
                };
            };
        }
        else
        {
            if (B.first_diag() == B.last_diag())
            {
                Integer rc      = B.diag_length(B.first_diag());
                Integer pos     = ci.r_start-1;
                Integer dpos    = imult(1 + r,ci.r_step);

                const val_type_2* ptr_B = B.rep_ptr() + B.first_elem_diag(B.first_diag());

                ptr_mat         = mat.ptr();

                if (mat_ld == matr)
                {
                    for (Integer j = 0; j < rc; ++j)
                    {					
                        mrd::assign_helper(ptr_mat[pos],val_type_1(*ptr_B));
                        pos     += dpos;
                        ptr_B   += B_ld;
                    };
                }
                else
                {
                    for (Integer j = 0; j < rc; ++j)
                    {					
                        while(pos >= matr)
                        {
                            pos     -= matr;
                            ptr_mat += mat_ld;
                        };

                        mrd::assign_helper(ptr_mat[pos],val_type_1(*ptr_B));

                        pos     += dpos;
                        ptr_B   += B_ld;
                    };
                };
            }
            else
            {
                const val_type_2* ptr_B = B.rep_ptr();
                Integer pos = ci.r_start - 1;

                for (Integer j = 0; j < Bc; ++j)
                {
                    Integer row_f = B.first_row(j);
                    Integer row_l = B.last_row(j);
                    Integer row_p = B.first_elem_pos(j);

                    pos +=  imult(row_f,ci.r_step);
                    if (mat_ld == matr)
                    {
                        for (Integer i = row_f; i <= row_l; ++i, ++row_p, pos += ci.r_step)
                            mrd::assign_helper(ptr_mat[pos],val_type_1(ptr_B[row_p]));
                    }
                    else
                    {
                        for (Integer i = row_f; i <= row_l; ++i, ++row_p, pos += ci.r_step)
                        {
                            while(pos >= matr)
                            {
                                pos     -= matr;
                                ptr_mat += mat_ld;
                            };

                            mrd::assign_helper(ptr_mat[pos],val_type_1(ptr_B[row_p]));
                        };
                    };

                    ptr_B   += B_ld;			
                    pos     += imult(std::max(Br-row_l-1,0),ci.r_step);
                };
            };
        };

        mat.get_struct().reset();
        return A;
    };

    static Matrix& eval(Matrix& A,const M2& B,const colon& c1,const colon& c2)
    {
        const M1& mat_tmp = A.get_impl<M1>();

        colon_info ci;		
        make_index(mat_tmp.rows(),mat_tmp.cols(),c1,c2,ci);

        Integer nr = ci.rows(), nc = ci.cols();

        error::check_assign_2(nr, nc, B.rows(), B.cols());

        if (nr == 0 || nc == 0)
            return A;

        M1& mat             = A.get_impl_unique<M1>();	
        using val_type_1    = typename M1::value_type;
        using val_type_2    = typename M2::value_type;
        val_type_1* ptr_mat = mat.ptr();

        val_type_1 Z = default_value<val_type_1>(get_matrix_ti<val_type_1>::eval(A));

        const Integer* ptr_ri   = ci.get_ri_2_ptr();
        const Integer* ptr_ci   = ci.get_ci_2_ptr();
        Integer Bc              = B.cols();
        Integer B_ld             = B.ld();
        Integer Ar              = A.rows();
        Integer cir             = ci.rows();
        Integer mat_ld           = mat.ld();

        if (ci.r_flag == 0 && ci.c_flag == 0)
        {
            const val_type_2* ptr_B = B.rep_ptr();

            for (Integer j = 0; j < Bc; ++j)
            {				
                Integer row_f   = B.first_row(j);
                Integer row_l   = B.last_row(j);
                Integer row_p   = B.first_elem_pos(j)-row_f;

                Integer pos_A   = imult(ptr_ci[j]-1,Ar);
                ptr_mat         = mat.ptr() + pos_A;

                Integer ci_ri_size  = cir;
                for (Integer i = 0; i < ci_ri_size; ++i, ++row_p)
                {
                    Integer p   = ptr_ri[i]-1;

                    if (i < row_f || i > row_l)
                        mrd::assign_helper(ptr_mat[p],Z);
                    else
                        mrd::assign_helper(ptr_mat[p],val_type_1(ptr_B[row_p]));
                };

                ptr_B           += B_ld;
            };

            mat.get_struct().reset();
            return A;
        };

        if (ci.r_flag == 0 && ci.c_flag == 1)
        {
            Integer pos_A       = imult(ci.c_start-1,mat_ld);
            Integer dpos_A      = imult(ci.c_step,mat_ld);

            const val_type_2* ptr_B = B.rep_ptr();
            ptr_mat             = mat.ptr() + pos_A;

            for (Integer j = 0; j < Bc; ++j)
            {
                Integer row_f   = B.first_row(j);
                Integer row_l   = B.last_row(j);
                Integer row_p   = B.first_elem_pos(j) - row_f;

                Integer ci_ri_size  = cir;
                for (Integer i = 0; i < ci_ri_size; ++i, ++row_p)
                {
                    Integer p   = ptr_ri[i]-1;

                    if (i < row_f || i > row_l)
                        mrd::assign_helper(ptr_mat[p],Z);
                    else
                        mrd::assign_helper(ptr_mat[p],val_type_1(ptr_B[row_p]));
                };	

                ptr_B           += B_ld;
                ptr_mat         += dpos_A;
            };

            mat.get_struct().reset();
            return A;
        };

        if (ci.r_flag == 1 && ci.c_flag == 0)
        {
            const val_type_2* ptr_B = B.rep_ptr();

            for (Integer j = 0; j < Bc; ++j)
            {
                Integer row_f   = B.first_row(j);
                Integer row_l   = B.last_row(j);
                Integer row_p   = B.first_elem_pos(j);

                Integer pos_A   = imult(ptr_ci[j]-1,Ar);
                ptr_mat         = mat.ptr() + pos_A;

                {
                    Integer pos_r = ci.r_start - 1;
                    for (Integer i = 0; i < ci.r_size; ++i, pos_r+= ci.r_step)
                        mrd::assign_helper(ptr_mat[pos_r],Z);
                };

                Integer i       = row_f;
                Integer pos_r   = ci.r_start - 1 + imult(i,ci.r_step);

                for (; i <= row_l; ++i, ++row_p, pos_r+= ci.r_step)
                    mrd::assign_helper(ptr_mat[pos_r],val_type_1(ptr_B[row_p]));

                ptr_B           += B_ld;
            };

            mat.get_struct().reset();
            return A;
        };

        if (ci.r_flag == 1 && ci.c_flag == 1)
        {
            algorithm::dense_change_entries(mat,ci,Z);
            Integer pos_A   = imult(ci.c_start-1,mat_ld);
            Integer dpos_A  = imult(ci.c_step,mat_ld);

            const val_type_2* ptr_B = B.rep_ptr();
            ptr_mat         = mat.ptr() + pos_A;

            for (Integer j = 0; j < Bc; ++j)
            {
                Integer row_f = B.first_row(j);
                Integer row_l = B.last_row(j);
                Integer row_p = B.first_elem_pos(j);

                Integer i = row_f;
                Integer pos_r = ci.r_start - 1 + imult(i,ci.r_step);

                for (; i <= row_l; ++i, ++row_p, pos_r+= ci.r_step)
                {
                    mrd::assign_helper(ptr_mat[pos_r],val_type_1(ptr_B[row_p]));
                };		
                ptr_B       += B_ld;
                ptr_mat     += dpos_A;
            };

            mat.get_struct().reset();
            return A;
        };

        matcl_assert(0,"unknown case");
        throw;
    };

    static Matrix& eval_diag(Matrix& A,const M2& B,Integer d)
    {
        using value_type    = typename M1::value_type;
        using FullMatrix    = raw::Matrix<value_type,struct_dense>;
        FullMatrix BC       = mr::converter<FullMatrix,M2>::eval(B);
        algorithm::dense_change_diag<value_type>(A.get_impl_unique<M1>(),d,BC);
        return A;
    }
};

template<class M1,class M2>
struct assign_mat_str<M1,M2,struct_banded,struct_dense>
{
    using FullMatrix        = raw::Matrix<typename M1::value_type,struct_dense>;
    using SparseMatrix      = raw::Matrix<typename M1::value_type,struct_sparse>;

    static Matrix& eval(Matrix& A,const M2& B,const colon& c1)
    {	
        const M1 & MA = A.get_impl<M1>();

        Real density = (Real(MA.nnz()) + Real(B.size())) / (MA.rows()+1.) / (MA.cols()+1.);

        if (density > optim_params::max_sparse_density_max)
        {
            A = Matrix(raw::converter<FullMatrix,M1>::eval(MA),false);
            return assign_mat_str<FullMatrix,M2,struct_dense,struct_dense>::eval(A,B,c1);
        }
        else
        {
            A = Matrix(raw::converter<SparseMatrix,M1>::eval(MA),false);
            return assign_mat_str<SparseMatrix,M2,struct_sparse,struct_dense>::eval(A,B,c1);
        };
    };

    static Matrix& eval(Matrix& A,const M2& B,const colon& c1,const colon& c2)
    {
        const M1 & MA = A.get_impl<M1>();

        Real density = (Real(MA.nnz()) + Real(B.size())) / (MA.rows()+1.) / (MA.cols()+1.);

        if (density > optim_params::max_sparse_density_max)
        {
            A = Matrix(raw::converter<FullMatrix,M1>::eval(MA),false);
            return assign_mat_str<FullMatrix,M2,struct_dense,struct_dense>::eval(A,B,c1,c2);
        }
        else
        {
            A = Matrix(raw::converter<SparseMatrix,M1>::eval(MA),false);
            return assign_mat_str<SparseMatrix,M2,struct_sparse,struct_dense>::eval(A,B,c1,c2);
        };
    };

    static Matrix& eval_diag(Matrix& A,const M2& B,Integer d)
    {
        using value_type    = typename M1::value_type;
        using FullMatrix    = raw::Matrix<value_type,struct_dense>;
        const FullMatrix& BC= mr::converter<FullMatrix,M2>::eval(B);

        algorithm::band_change_diag<value_type>(A, A.get_impl<M1>(), d, BC);
        return A;
    }
};

template<class M1,class M2>
struct assign_mat_str<M1,M2,struct_banded,struct_sparse>
{
    using FullMatrix        = raw::Matrix<typename M1::value_type,struct_dense>;
    using SparseMatrix      = raw::Matrix<typename M1::value_type,struct_sparse>;

    static Matrix& eval(Matrix& A,const M2& B,const colon& c1)
    {	
        const M1 & MA = A.get_impl<M1>();

        Real density = (Real(MA.nnz()) + Real(B.nnz())) / (MA.rows()+1.) / (MA.cols()+1.);

        if (density > optim_params::max_sparse_density_max)
        {
            A = Matrix(raw::converter<FullMatrix,M1>::eval(MA),false);
            return assign_mat_str<FullMatrix,M2,struct_dense,struct_sparse>::eval(A,B,c1);
        }
        else
        {
            A = Matrix(raw::converter<SparseMatrix,M1>::eval(MA),false);
            return assign_mat_str<SparseMatrix,M2,struct_sparse,struct_sparse>::eval(A,B,c1);
        };
    };

    static Matrix& eval(Matrix& A,const M2& B,const colon& c1,const colon& c2)
    {
        const M1 & MA = A.get_impl<M1>();

        Real density = (Real(MA.nnz()) + Real(B.nnz())) / (MA.rows()+1.) / (MA.cols()+1.);

        if (density > optim_params::max_sparse_density_max)
        {
            A = Matrix(raw::converter<FullMatrix,M1>::eval(MA),false);
            return assign_mat_str<FullMatrix,M2,struct_dense,struct_sparse>::eval(A,B,c1,c2);
        }
        else
        {
            A = Matrix(raw::converter<SparseMatrix,M1>::eval(MA),false);
            return assign_mat_str<SparseMatrix,M2,struct_sparse,struct_sparse>::eval(A,B,c1,c2);
        };
    };

    static Matrix& eval_diag(Matrix& A,const M2& B,Integer d)
    {
        using value_type    = typename M1::value_type;
        using FullMatrix    = raw::Matrix<value_type,struct_dense>;
        FullMatrix BC       = mr::converter<FullMatrix,M2>::eval(B);

        algorithm::band_change_diag<value_type>(A, A.get_impl<M1>(), d, BC);
        return A;
    }
};

template<class M1,class M2>
struct assign_mat_str<M1,M2,struct_banded,struct_banded>
{
    using FullMatrix        = raw::Matrix<typename M1::value_type,struct_dense>;
    using SparseMatrix      = raw::Matrix<typename M1::value_type,struct_sparse>;

    static Matrix& eval(Matrix& A,const M2& B,const colon& c1)
    {	
        const M1 & MA = A.get_impl<M1>();

        Real density = (Real(MA.nnz()) + Real(B.nnz())) / (MA.rows()+1.) / (MA.cols()+1.);

        if (density > optim_params::max_sparse_density_max)
        {
            A = Matrix(raw::converter<FullMatrix,M1>::eval(MA),false);
            return assign_mat_str<FullMatrix,M2,struct_dense,struct_banded>::eval(A,B,c1);
        }
        else
        {
            A = Matrix(raw::converter<SparseMatrix,M1>::eval(MA),false);
            return assign_mat_str<SparseMatrix,M2,struct_sparse,struct_banded>::eval(A,B,c1);
        };
    };

    static Matrix& eval(Matrix& A,const M2& B,const colon& c1,const colon& c2)
    {
        const M1 & MA = A.get_impl<M1>();

        Real density = (Real(MA.nnz()) + Real(B.nnz())) / (MA.rows()+1.) / (MA.cols()+1.);

        if (density > optim_params::max_sparse_density_max)
        {
            A = Matrix(raw::converter<FullMatrix,M1>::eval(MA),false);
            return assign_mat_str<FullMatrix,M2,struct_dense,struct_banded>::eval(A,B,c1,c2);
        }
        else
        {
            A = Matrix(raw::converter<SparseMatrix,M1>::eval(MA),false);
            return assign_mat_str<SparseMatrix,M2,struct_sparse,struct_banded>::eval(A,B,c1,c2);
        };
    };

    static Matrix& eval_diag(Matrix& A,const M2& B,Integer d)
    {
        using value_type    = typename M1::value_type;
        using FullMatrix    = raw::Matrix<value_type,struct_dense>;
        FullMatrix BC       = mr::converter<FullMatrix,M2>::eval(B);

        algorithm::band_change_diag<value_type>(A, A.get_impl<M1>(), d, BC);
        return A;
    }
};

template<class M1,class M2,class val_type>
struct assign_mat_impl<M1,M2,val_type,true>
{
    using struct_type_1 = typename M1::struct_type;
    using val_type_1    = typename M1::value_type;
    using val_type_2    = typename M2::value_type;
    using new_mat_type  = raw::Matrix<val_type,struct_type_1>;

    static Matrix& eval(Matrix& A,const M2& B,const colon& c1)
    {
        ti::ti_type<val_type_1> ti_1 = get_matrix_ti<val_type_1>::eval(A);
        ti::ti_type<val_type_2> ti_2 = get_matrix_ti<val_type_2>::eval(B);
        ti::ti_type<val_type> ret_ti = link_ret_ti<val_type,val_type_1,val_type_2>::eval(ti_1,ti_2);
        new_mat_type An = raw::converter<new_mat_type,M1>::eval(A.get_impl<M1>(),ret_ti);
        A = Matrix(An,false);
        return assign_mat_impl<new_mat_type,M2,val_type,false>::eval(A,B,c1);
    };

    static Matrix& eval(Matrix& A,const M2& B,const colon& c1,const colon& c2)
    {
        ti::ti_type<val_type_1> ti_1 = get_matrix_ti<val_type_1>::eval(A);
        ti::ti_type<val_type_2> ti_2 = get_matrix_ti<val_type_2>::eval(B);
        ti::ti_type<val_type> ret_ti = link_ret_ti<val_type,val_type_1,val_type_2>::eval(ti_1,ti_2);

        new_mat_type An = raw::converter<new_mat_type,M1>::eval(A.get_impl<M1>(),ret_ti);
        A = Matrix(An,false);
        return assign_mat_impl<new_mat_type,M2,val_type,false>::eval(A,B,c1,c2);
    };

    static Matrix& eval_diag(Matrix& A,const M2& B,Integer d)
    {
        ti::ti_type<val_type_1> ti_1 = get_matrix_ti<val_type_1>::eval(A);
        ti::ti_type<val_type_2> ti_2 = get_matrix_ti<val_type_2>::eval(B);
        ti::ti_type<val_type> ret_ti = link_ret_ti<val_type,val_type_1,val_type_2>::eval(ti_1,ti_2);

        new_mat_type An = raw::converter<new_mat_type,M1>::eval(A.get_impl<M1>(),ret_ti);
        A = Matrix(An,false);
        return assign_mat_impl<new_mat_type,M2,val_type,false>::eval_diag(A,B,d);
    };
};

template<class M1,class M2,class val_type>
struct assign_mat_impl<M1,M2,val_type,false>
{
    using str_1 = typename M1::struct_type;
    using str_2 = typename M2::struct_type;

    static Matrix& eval(Matrix& A,const M2& B,const colon& c1)
    {
        return assign_mat_str<M1,M2,str_1,str_2>::eval(A,B,c1);
    };

    static Matrix& eval(Matrix& A,const M2& B,const colon& c1,const colon& c2)
    {
        return assign_mat_str<M1,M2,str_1,str_2>::eval(A,B,c1,c2);
    };

    static Matrix& eval_diag(Matrix& A,const M2& B,Integer d)
    {
        return assign_mat_str<M1,M2,str_1,str_2>::eval_diag(A,B,d);
    };
};

template<class M1, class M2, bool is_band>
struct change_submatrix_dense_help
{
    static void eval(Matrix& ret, const M1& A,const colon_info& ci, const M2& B)
    {
        return algorithm::change_submatrix_dense(ret, A, ci, B);
    };

    static void eval_2(Matrix& ret, const M1& A,const colon_info& ci, const M2& B)
    {
        return algorithm::change_submatrix_dense_2(ret, A, ci, B);
    };
};

template<class M1, class M2>
struct change_submatrix_dense_help<M1,M2,true>
{
    static void eval(Matrix& ret, const M1& A,const colon_info& ci, const M2& B)
    {
        return algorithm::change_submatrix_band(ret, A, ci, B);
    };

    static void eval_2(Matrix& ret, const M1& A,const colon_info& ci, const M2& B)
    {
        return algorithm::change_submatrix_band_2(ret, A, ci, B);
    };
};

template<class M1,class M2, class struct_type>
struct assign_mat_sp_nonsp
{
    static Matrix& eval(Matrix& A,const M2& B,const colon& c1)
    {
        const M1& mat = A.get_impl<M1>();

        colon_info ci; 
        make_index(mat.rows(),mat.cols(),c1,ci);		

        Integer s = ci.rows();

        if (ci.r_flag == 1)
            error::check_assign_1(s, B.rows(), B.cols());
        else
            error::check_assign_1(ci.get_rim_1().size(), B.rows(), B.cols());

        if (s == 0)
            return A;

        Real density = (Real(mat.nnz()) + Real(B.nnz())) / (mat.rows()+1.) / (mat.cols()+1.);

        if (density > optim_params::max_sparse_density_max)
        {
            using value_type    = typename M1::value_type;
            using matrix_type   = raw::Matrix<value_type,struct_dense>;
            A                   = matcl::full(Matrix(mat,false));

            return assign_mat_str<matrix_type,M2,struct_dense,struct_type>::eval(A,B,c1);
        };

        M1& mat2 = A.get_impl_unique<M1>();

        static const bool is_band = std::is_same<struct_type,struct_banded>::value;

        Matrix ret;
        change_submatrix_dense_help<M1,M2,is_band>::eval_2(ret, mat2, ci, B);

        A = ret;
        return A;
    };

    static Matrix& eval(Matrix& A,const M2& B,const colon& c1,const colon& c2)
    {
        const M1& mat = A.get_impl<M1>();

        colon_info ci;
        make_index(mat.rows(),mat.cols(),c1,c2,ci);

        Integer nr = ci.rows(), nc = ci.cols();

        error::check_assign_2(nr, nc, B.rows(), B.cols());

        if (nr == 0 || nc == 0)
            return A;

        Real density = (Real(mat.nnz()) + Real(B.nnz())) / (mat.rows()+1.) / (mat.cols()+1.);
        if (density > optim_params::max_sparse_density_max)
        {
            using value_type    = typename M1::value_type;
            using matrix_type   = raw::Matrix<value_type,struct_dense>;
            A                   = matcl::full(Matrix(mat,false));

            return assign_mat_str<matrix_type,M2,struct_dense,struct_type>::eval(A,B,c1,c2);
        };

        M1& mat2 = A.get_impl_unique<M1>();
        
        static const bool is_band = std::is_same<struct_type,struct_banded>::value;
        Matrix ret;
        change_submatrix_dense_help<M1,M2,is_band>::eval(ret, mat2, ci, B);

        A = ret;
        return A;
    };
};

template<class M1,class M2>
struct assign_mat_str<M1,M2,struct_sparse,struct_dense>
{
    static Matrix& eval(Matrix& A,const M2& B,const colon& c1)
    {
        return assign_mat_sp_nonsp<M1,M2,struct_dense>::eval(A,B,c1);
    }

    static Matrix& eval(Matrix& A,const M2& B,const colon& c1,const colon& c2)
    {
        return assign_mat_sp_nonsp<M1,M2,struct_dense>::eval(A,B,c1,c2);
    };

    static Matrix& eval_diag(Matrix& A,const M2& B,Integer d)
    {
        using value_type    = typename M1::value_type;
        using FullMatrix    = raw::Matrix<value_type,struct_dense>;
        const FullMatrix& BC= mr::converter<FullMatrix,M2>::eval(B);

        algorithm::sparse_change_diag<value_type>(A, A.get_impl_unique<M1>(), d, BC);
        return A;
    }
};

template<class M1,class M2>
struct assign_mat_str<M1,M2,struct_sparse,struct_banded>
{
    static Matrix& eval(Matrix& A,const M2& B,const colon& c1)
    {
        return assign_mat_sp_nonsp<M1,M2,struct_banded>::eval(A,B,c1);
    }

    static Matrix& eval(Matrix& A,const M2& B,const colon& c1,const colon& c2)
    {
        return assign_mat_sp_nonsp<M1,M2,struct_banded>::eval(A,B,c1,c2);
    };

    static Matrix& eval_diag(Matrix& A,const M2& B,Integer d)
    {
        using value_type    = typename M1::value_type;
        using FullMatrix    = raw::Matrix<value_type,struct_dense>;
        FullMatrix BC       = mr::converter<FullMatrix,M2>::eval(B);

        algorithm::sparse_change_diag<value_type>(A, A.get_impl_unique<M1>(), d, BC);
        return A;
    }
};

template<class M1,class M2>
struct assign_functor_mat
{
    using value_type_1  = typename M1::value_type;
    using value_type_2  = typename M2::value_type;
    using val_type_pr   = typename md::unify_types<value_type_1,value_type_2>::type;

    static const bool req_promo = std::is_same<value_type_1,val_type_pr>::value == false
                                || std::is_same<value_type_1,Object>::value == true;

    static Matrix& eval(Matrix& A, const M2& B,const colon& c1)
    {
        if (c1.m_flag == colon::t_all)
        {
            if (B.rows() == 0 && B.cols() == 0)
            {
                A= Matrix(raw::Matrix<value_type_1,struct_dense>(get_matrix_ti<value_type_1>::eval(A)),true);
                return A;
            }
         
            //for dense matrices we want to modify elements in the array
            if (A.get_struct_code() != struct_code::struct_dense)
            {
                const M1& mat = A.get_impl<M1>();

                if (mat.rows() == B.rows() && mat.cols() == B.cols())
                    A = Matrix(B,true);
            
                error::check_assign_1(imult_c(mat.rows(),mat.cols()), B.rows(), B.cols());

                mrd::manip_reshape_helper<M2>::eval_reshape(A, B,mat.rows(),mat.cols());
                return A;
            };
        }
        else if (B.rows() == 0 && B.cols() == 0)
        {
            colon_info ci;
            make_index(A.rows(), A.cols(), c1, ci);

            Integer nr = ci.rows();

            if (nr == 0)
            {
                return A;
            }
            else
            {
                A = delrows(vec(A),c1);
                return A;
            };
        };

        return assign_mat_impl<M1,M2,val_type_pr,req_promo>::eval(A,B,c1);
    };

    static Matrix& eval(Matrix& A,const M2& B,const colon& c1,const colon& c2)
    {
        if (c1.m_flag == colon::t_all)
        {            
            if (c2.m_flag == colon::t_all)
            {
                if (B.rows() == 0 && B.cols() == 0)
                {
                    A = Matrix(raw::Matrix<value_type_1,struct_dense>(get_matrix_ti<value_type_1>::eval(A)),true);
                    return A;
                };

                //for dense matrices we want to modify elements in the array
                if (A.get_struct_code() != struct_code::struct_dense)
                {
                    const M1& mat = A.get_impl<M1>();
                    error::check_assign_2(mat.rows(), mat.cols(), B.rows(), B.cols());

                    A = Matrix(B,true);
                    return A;
                };
            }
            else if (B.rows() == 0 && B.cols() == 0)
            {
                A = delcols(A,c2);
                return A;
            };
        }
        else if (c2.m_flag == colon::t_all && B.rows() == 0 && B.cols() == 0)
        {
            A = delrows(A,c1);
            return A;
        };

        return assign_mat_impl<M1,M2,val_type_pr,req_promo>::eval(A,B,c1,c2);
    };

    static Matrix& eval_diag(Matrix& A,const M2& B,Integer d)
    {
        return assign_mat_impl<M1,M2,val_type_pr,req_promo>::eval_diag(A,B,d);
    };
};

};};
