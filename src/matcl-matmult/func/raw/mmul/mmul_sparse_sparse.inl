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

#pragma once

namespace matcl { namespace raw { namespace details
{

//========================================================================
//                      SPARSE - SPARSE
//========================================================================

template<class Val_ret, class M1, class M2>
struct eval_mult<Val_ret,M1,M2,struct_sparse,struct_sparse> 
{ 
    using Mat_D         = raw::Matrix<Val_ret, struct_dense>;
    using ti_ret_type   = typename ti::get_ti_type<Val_ret>::type;

    static void eval(matcl::Matrix& ret, const M1& A, const M2& B, trans_type t_A, 
                          trans_type t_B)
    {
        if (t_B != trans_type::no_trans)
        {
            ti_ret_type ret_ti = ti::get_return_ti<ti_ret_type>(mdyf::op_mul::eval(),
                                                                  ti::get_ti(A), ti::get_ti(B));

            M2 Bt(ret_ti);

            if (t_B == trans_type::trans)
                mrd::manip_trans_helper<M2>::eval_trans(Bt, B);
            else
                mrd::manip_trans_helper<M2>::eval_ctrans(Bt, B);

            return eval(ret, A, Bt, t_A, trans_type::no_trans);
        };

        if (t_A == trans_type::no_trans)
            return eval_notrans(ret, A, B);
        else
            return eval_trans(ret, A, B, t_A);
    };

    static void eval_gemm(Mat_D& C, const Val_ret& alpha, const M1& A, const M2& B, 
                     trans_type t_A, trans_type t_B, const Val_ret& beta, 
                     Integer fr, Integer rows, Integer fc, Integer cols)
    {
        prepare_gemm_C<Val_ret>::eval(beta, C, fr,rows,fc,cols);

        if (mrd::is_one(alpha))
        {
            using Alpha     = alpha_one<Val_ret>;
            Alpha al;

            if (t_B == trans_type::no_trans)
            {
                if (t_A == trans_type::no_trans)
                    return eval_NN_gemm(al, C, A, B, t_A, t_B, fr, rows, fc, cols);
                else if (t_A == trans_type::trans)
                    return eval_TN_gemm<false>(alpha, C, A, B, t_A, t_B, fr, rows, fc, cols);
                else
                    return eval_TN_gemm<true>(alpha, C, A, B, t_A, t_B, fr, rows, fc, cols);
            }
            else if (t_B == trans_type::trans)
            {
                if (t_A == trans_type::no_trans)
                    return eval_NT_gemm(al, C, A, B, t_A, t_B, fr, rows, fc, cols);
                else
                    return eval_TT_gemm<false>(al, C, A, B, t_A, t_B, fr, rows, fc, cols);
            }
            else
            {
                if (t_A == trans_type::no_trans)
                    return eval_NT_gemm(al, C, A, B, t_A, t_B, fr, rows, fc, cols);
                else
                    return eval_TT_gemm<true>(al, C, A, B, t_A, t_B, fr, rows, fc, cols);
            };
        }
        else if (mrd::is_one(-alpha))
        {
            using Alpha     = alpha_mone<Val_ret>;
            Alpha al;

            if (t_B == trans_type::no_trans)
            {
                if (t_A == trans_type::no_trans)
                    return eval_NN_gemm(al, C, A, B, t_A, t_B, fr, rows, fc, cols);
                else if (t_A == trans_type::trans)
                    return eval_TN_gemm<false>(alpha, C, A, B, t_A, t_B, fr, rows, fc, cols);
                else
                    return eval_TN_gemm<true>(alpha, C, A, B, t_A, t_B, fr, rows, fc, cols);
            }
            else if (t_B == trans_type::trans)
            {
                if (t_A == trans_type::no_trans)
                    return eval_NT_gemm(al, C, A, B, t_A, t_B, fr, rows, fc, cols);
                else
                    return eval_TT_gemm<false>(al, C, A, B, t_A, t_B, fr, rows, fc, cols);
            }
            else
            {
                if (t_A == trans_type::no_trans)
                    return eval_NT_gemm(al, C, A, B, t_A, t_B, fr, rows, fc, cols);
                else
                    return eval_TT_gemm<true>(al, C, A, B, t_A, t_B, fr, rows, fc, cols);
            };
        }
        else
        {
            using Alpha     = alpha_val<Val_ret>;
            Alpha al(alpha);

            if (t_B == trans_type::no_trans)
            {
                if (t_A == trans_type::no_trans)
                    return eval_NN_gemm(al, C, A, B, t_A, t_B, fr, rows, fc, cols);
                else if (t_A == trans_type::trans)
                    return eval_TN_gemm<false>(alpha, C, A, B, t_A, t_B, fr, rows, fc, cols);
                else
                    return eval_TN_gemm<true>(alpha, C, A, B, t_A, t_B, fr, rows, fc, cols);
            }
            else if (t_B == trans_type::trans)
            {
                if (t_A == trans_type::no_trans)
                    return eval_NT_gemm(al, C, A, B, t_A, t_B, fr, rows, fc, cols);
                else
                    return eval_TT_gemm<false>(al, C, A, B, t_A, t_B, fr, rows, fc, cols);
            }
            else
            {
                if (t_A == trans_type::no_trans)
                    return eval_NT_gemm(al, C, A, B, t_A, t_B, fr, rows, fc, cols);
                else
                    return eval_TT_gemm<true>(al, C, A, B, t_A, t_B, fr, rows, fc, cols);
            };
        };
    };

    template<class Alpha>
    static void eval_NT_gemm(const Alpha& alpha, Mat_D& C, const M1& A, const M2& B, 
                     trans_type t_A, trans_type t_B, Integer fr, Integer rows, Integer fc, Integer cols)
    {
        //A * B'

        (void)t_A;
        (void)t_B;
        (void)rows;
        (void)cols;

        using V1        = typename M1::value_type;
        using V2        = typename M2::value_type;
        using V12       = typename md::unify_types<V1,V2>::type;
        using VR        = Val_ret;

        if (A.nnz() == 0 || B.nnz() == 0)
            return;

        Integer B_cols  = B.cols();
        bool conj_B     = t_B == trans_type::conj_trans;

        const sparse_ccs<V1>& Ad = A.rep();
        const sparse_ccs<V2>& Bd = B.rep();

        const Integer * A_c		= Ad.ptr_c();
        const Integer * A_r		= Ad.ptr_r();
        const V1 * A_x	        = Ad.ptr_x();

        const Integer * B_c		= Bd.ptr_c();
        const Integer * B_r		= Bd.ptr_r();
        const V2 * B_x	        = Bd.ptr_x();

        Integer C_ld            = C.ld();
        VR* ptr_C               = C.ptr() + fr + fc * C_ld;

        for (Integer j = 0; j < B_cols; ++j)
        {
            if (A_c[j] == A_c[j + 1])
                continue;

            for (Integer k = B_c[j] ; k < B_c[j + 1] ; ++k)
            {
                V2 Bx           = conj_B ? conj(B_x[k]) : B_x[k];
                Integer Br      = B_r[k];

                if (mrd::is_zero(Bx))
                    continue;

                VR* ptr_Cl      = ptr_C + Br * C_ld;

                for (Integer ka = A_c[j]; ka < A_c[j + 1]; ++ka)
                {
                    Integer i   = A_r[ka];
                    V12 tmp     = A_x[ka] * Bx;
                    ptr_Cl[i]   = alpha.eval(ptr_Cl[i], tmp);
                };
            };
        };
        
        return;
    }

    template<bool Conj_B, class Alpha>
    static void eval_TT_gemm(const Alpha& alpha, Mat_D& C, const M1& A, const M2& B, 
                     trans_type t_A, trans_type t_B, Integer fr, Integer rows, Integer fc, Integer cols)
    {
        //A' * B' = (B * A)'

        (void)t_B;
        (void)rows;
        (void)cols;

        using V1        = typename M1::value_type;
        using V2        = typename M2::value_type;
        using V12       = typename md::unify_types<V1,V2>::type;
        using VR        = Val_ret;

        Integer A_cols  = A.cols();        

        if (A.nnz() == 0 || B.nnz() == 0)
            return;

        const sparse_ccs<V1>& Ad = A.rep();
        const sparse_ccs<V2>& Bd = B.rep();

        const Integer * A_c		= Ad.ptr_c();
        const Integer * A_r		= Ad.ptr_r();
        const V1 * A_x	        = Ad.ptr_x();

        const Integer * B_c		= Bd.ptr_c();
        const Integer * B_r		= Bd.ptr_r();
        const V2 * B_x	        = Bd.ptr_x();

        Integer C_ld            = C.ld();
        VR* ptr_C               = C.ptr() + fr + fc * C_ld;

        bool conj_A             = (t_A == trans_type::conj_trans);

        for (Integer j = 0; j < A_cols; ++j)
        {
            for (Integer k = A_c[j] ; k < A_c[j + 1] ; ++k)
            {
                V1 Ax           = conj_A ? conj(A_x[k]) : A_x[k];
                Integer Ar      = A_r[k];

                if (mrd::is_zero(Ax))
                    continue;

                for (Integer kb = B_c[Ar]; kb < B_c[Ar + 1]; ++kb)
                {
                    Integer i       = B_r[kb];
                    V12 tmp         = Ax * make_conj<Conj_B>::eval(B_x[kb]);
                    ptr_C[i*C_ld]   = alpha.eval(ptr_C[i*C_ld], tmp);
                };
            };

            ptr_C               += 1;
        };
        
        return;
    };

    static void eval_trans(matcl::Matrix& ret, const M1& A, const M2& B, trans_type t_A)
    {
        using val_type_1    = typename M1::value_type;
        using val_type_2    = typename M2::value_type;

        ti_ret_type ret_ti = ti::get_return_ti<ti_ret_type>(mdyf::op_mul::eval(),
                                                                  ti::get_ti(A), ti::get_ti(B));

        Integer N = B.cols();
        Integer M = A.cols();
        Integer K = A.rows();

        if (A.nnz() == 0 || B.nnz() == 0)
        {
            using ret_type  = raw::Matrix<Val_ret,struct_sparse>;

            ret_type out(ret_ti, M, N);
            ret = matcl::Matrix(out,false);
            return;
        };

        Integer nz_est = estim_mult_nnz(A, B, t_A, trans_type::no_trans);

        using ret_type  = raw::Matrix<Val_ret,struct_sparse>;

        ret_type C(ret_ti, M, N, nz_est + M);

        sparse_ccs<Val_ret>& d           = C.rep();
        const sparse_ccs<val_type_1>& Ad = A.rep();
        const sparse_ccs<val_type_2>& Bd = B.rep();

        Integer nz		    = 0;

        Integer * d_c		= d.ptr_c();
        Integer * d_r		= d.ptr_r();
        Val_ret * d_x	    = d.ptr_x();

        const Integer * A_c		= Ad.ptr_c();
        const Integer * A_r		= Ad.ptr_r();
        const val_type_1 * A_x	= Ad.ptr_x();

        const Integer * B_c		= Bd.ptr_c();
        const Integer * B_r		= Bd.ptr_r();
        const val_type_2 * B_x	= Bd.ptr_x();

        using workspace = md::workspace2<val_type_2>;
        using scatter   = matcl::algorithm::scatter;

        scatter sc      = scatter::get(K, N);
        workspace work_x(ret_ti,K);

        Val_ret Z       = md::default_value<Val_ret>(ret_ti);

        for (Integer j = 0; j < N; ++j)
        {
            d_c[j]      = nz;
            Integer dnz = B_c[j + 1] - B_c[j];			

            if (dnz == 0)		
                continue;

            if (nz + M > d.nzmax()) 
            {
                d.add_memory( d.nzmax() + M);

                d_r		= d.ptr_r();
                d_x		= d.ptr_x();
            };

            auto mark       = sc.next_mark();

            for (Integer k = B_c[j] ; k < B_c[j + 1] ; ++k)
            {
                const val_type_2& Bx = B_x[k];                

                if (mrd::is_zero(Bx))
                    continue;

                Integer Br  = B_r[k];
                sc[Br]      = mark;
                work_x[Br]  = Bx;
            };

            if (t_A == trans_type::trans)
            {
                for (Integer i = 0; i < M; ++i)
                {
                    Integer last    = A_c[i + 1];
                    Integer first   = A_c[i];

                    if (last == first)
                        continue;

                    Val_ret dot = Z;                    

                    for (Integer ka = first; ka < last; ++ka)
                    {
                        Integer Ar = A_r[ka];
                        if (sc[Ar] == mark)
                        {
                            dot = dot + A_x[ka] * work_x[Ar];
                        };
                    };

                    if (mrd::is_zero(dot) == false)
                    {
                        d_r[nz]     = i;
                        d_x[nz]     = dot;
                        ++nz;
                    };
                };
            }
            else
            {
                for (Integer i = 0; i < M; ++i)
                {
                    Integer last    = A_c[i + 1];
                    Integer first   = A_c[i];

                    if (last == first)
                        continue;

                    Val_ret dot = Z;

                    for (Integer ka = first; ka < last; ++ka)
                    {
                        Integer Ar = A_r[ka];
                        if (sc[Ar] == mark)
                        {
                            val_type_1 a = mrd::conj_helper<val_type_1>::eval(A_x[ka]);
                            dot = dot + a * work_x[Ar];
                        };
                    };

                    if (mrd::is_zero(dot) == false)
                    {
                        d_r[nz]     = i;
                        d_x[nz]     = dot;
                        ++nz;
                    };
                };
            };
        };

        d_c[N] = nz;

        bool is_sq_A    = A.rows() == A.cols();
        bool is_sq_B    = B.rows() == B.cols();
        bool is_sq_C    = C.rows() == C.cols();
        C.set_struct(predefined_struct_ext::mult_struct(A.get_struct(),B.get_struct(),t_A,trans_type::no_trans,
                            is_real_matrix(A),is_real_matrix(B), is_sq_A, is_sq_B, is_sq_C));

        ret = matcl::Matrix(C,true);
        return;
    };

    template<bool conj_A>
    static void eval_TN_gemm(const Val_ret& alpha, Mat_D& C, const M1& A, const M2& B, 
                     trans_type t_A, trans_type t_B, Integer fr, Integer rows, Integer fc, Integer cols)
    {
        (void)t_A;
        (void)t_B;

        using V1        = typename M1::value_type;
        using V2        = typename M2::value_type;
        using V12       = typename md::unify_types<V1,V2>::type;
        using VR        = Val_ret;

        Integer M       = rows;
        Integer N       = cols;
        Integer K       = B.rows();

        if (A.nnz() == 0 || B.nnz() == 0)
            return;

        const sparse_ccs<V1>& Ad = A.rep();
        const sparse_ccs<V2>& Bd = B.rep();

        const Integer * A_c		= Ad.ptr_c();
        const Integer * A_r		= Ad.ptr_r();
        const V1 * A_x	        = Ad.ptr_x();

        const Integer * B_c		= Bd.ptr_c();
        const Integer * B_r		= Bd.ptr_r();
        const V2 * B_x	        = Bd.ptr_x();

        Integer C_ld            = C.ld();
        VR* ptr_C               = C.ptr() + fr + fc * C_ld;

        VR Z                    = md::default_value<VR>(C.get_type());

        using workspace         = md::workspace2<V2>;
        using scatter           = matcl::algorithm::scatter;

        scatter sc              = scatter::get(K, N);
        workspace work_x(B.get_type(),K);

        for (Integer j = 0; j < N; ++j)
        {
            Integer dnz         = B_c[j + 1] - B_c[j];			

            if (dnz == 0)		
            {
                ptr_C           += C_ld;
                continue;
            };

            auto mark           = sc.next_mark();

            for (Integer k = B_c[j] ; k < B_c[j + 1] ; ++k)
            {
                const V2& Bx    = B_x[k];                

                if (mrd::is_zero(Bx))
                    continue;

                Integer Br      = B_r[k];
                sc[Br]          = mark;
                work_x[Br]      = Bx;
            };

            for (Integer i = 0; i < M; ++i)
            {
                Integer last    = A_c[i + 1];
                Integer first   = A_c[i];

                if (last == first)
                    continue;

                Val_ret dot     = Z;                    

                for (Integer ka = first; ka < last; ++ka)
                {
                    Integer Ar  = A_r[ka];

                    if (sc[Ar] == mark)
                        dot = dot + make_conj<conj_A>::eval(A_x[ka]) * work_x[Ar];
                };

                if (mrd::is_zero(dot) == false)
                    ptr_C[i]    = ptr_C[i] + alpha * dot;
            };

            ptr_C   += C_ld;
        };

        return;
    };

    static void eval_notrans(matcl::Matrix& ret, const M1& A, const M2& B)
    {		
        using val_type_1    = typename M1::value_type;
        using val_type_2    = typename M2::value_type;

        ti_ret_type ret_ti = ti::get_return_ti<ti_ret_type>(mdyf::op_mul::eval(),
                                                                  ti::get_ti(A), ti::get_ti(B));

        Integer N = B.cols();
        Integer M = A.rows();

        if (A.nnz() == 0 || B.nnz() == 0)
        {
            using ret_type  = raw::Matrix<Val_ret,struct_sparse>;

            ret_type out(ret_ti, M, N);
            ret = matcl::Matrix(out,false);
            return;
        };

        Integer nz_est = estim_mult_nnz(A, B, trans_type::no_trans, trans_type::no_trans);

        using ret_type  = raw::Matrix<Val_ret,struct_sparse>;

        ret_type C(ret_ti, M, N, nz_est + M);

        sparse_ccs<Val_ret>& d           = C.rep();
        const sparse_ccs<val_type_1>& Ad = A.rep();
        const sparse_ccs<val_type_2>& Bd = B.rep();

        using workspace = md::workspace2<Val_ret>;
        using scatter   = matcl::algorithm::scatter;

        scatter sc = scatter::get(M, N);
        workspace work_x(ret_ti,M);

        Integer nz		    = 0;

        Integer * d_c		= d.ptr_c();
        Integer * d_r		= d.ptr_r();
        Val_ret * d_x	    = d.ptr_x();

        const Integer * A_c		= Ad.ptr_c();
        const Integer * A_r		= Ad.ptr_r();
        const val_type_1 * A_x	= Ad.ptr_x();

        const Integer * B_c		= Bd.ptr_c();
        const Integer * B_r		= Bd.ptr_r();
        const val_type_2 * B_x	= Bd.ptr_x();

        for (Integer j = 0; j < N; ++j)
        {
            d_c[j]      = nz;
            Integer dnz = B_c[j + 1] - B_c[j];			

            if (dnz <= 1)
            {				
                if (dnz == 0)		
                    continue;

                Integer k               = B_c[j];
                const val_type_2& Bx    = B_x[k];
                Integer Br              = B_r[k];

                if (mrd::is_zero(Bx))
                    continue;

                dnz = A_c[Br + 1] - A_c[Br];

                if (nz + dnz > d.nzmax()) 
                {
                    d.add_memory( d.nzmax() + dnz);

                    d_r		= d.ptr_r();
                    d_x		= d.ptr_x();
                };

                for (Integer ka = A_c[Br]; ka < A_c[Br + 1]; ++ka)
                {
                    d_r[nz]	= A_r[ka];
                    d_x[nz] = A_x[ka] * Bx;
                    ++nz;
                };
                continue;
            };

            if (nz + M > d.nzmax()) 
            {
                d.add_memory( d.nzmax() + M);

                d_r		= d.ptr_r();
                d_x		= d.ptr_x();
            };

            Integer nz_old  = nz;
            Integer nk      = 0;
            auto mark       = sc.next_mark();

            for (Integer k = B_c[j] ; k < B_c[j + 1] ; ++k)
            {
                const val_type_2& Bx    = B_x[k];
                Integer Br              = B_r[k];

                if (mrd::is_zero(Bx))
                    continue;

                Integer nzs = nz;

                for (Integer ka = A_c[Br]; ka < A_c[Br + 1]; ++ka)
                {
                    Integer i = A_r[ka];

                    if (sc[i] < mark)
                    {
                        sc[i]       = mark;
                        work_x[i]	= A_x[ka] * Bx;
                        d_r[nz]		= i;
                        ++nz;
                    }
                    else 
                    {
                        Val_ret tmp  = A_x[ka] * Bx;
                         work_x[i]   = mrd::plus_helper<Val_ret,Val_ret>
                                            ::eval(work_x[i], tmp);
                    };
                };

                if (nz - nzs > 0)
                    ++nk;
            }

            if (nk > 1)
                utils::sort_q(d_r+nz_old,nz - nz_old);

            for (Integer k = d_c[j]; k < nz ; ++k) 
            {
                d_x[k] = work_x[d_r[k]];
            };
        };
        
        d_c[N] = nz;

        trans_type nt = trans_type::no_trans;

        bool is_sq_A    = A.rows() == A.cols();
        bool is_sq_B    = B.rows() == B.cols();
        bool is_sq_C    = C.rows() == C.cols();
        C.set_struct(predefined_struct_ext::mult_struct(A.get_struct(),B.get_struct(),nt,nt,
                            is_real_matrix(A),is_real_matrix(B), is_sq_A, is_sq_B, is_sq_C));

        ret = matcl::Matrix(C,true);
        return;
    };

    template<class Alpha>
    static void eval_NN_gemm(const Alpha& alpha, Mat_D& C, const M1& A, const M2& B, 
                     trans_type t_A, trans_type t_B, Integer fr, Integer rows, Integer fc, Integer cols)
    {
        (void)t_A;
        (void)t_B;
        (void)rows;

        using V1        = typename M1::value_type;
        using V2        = typename M2::value_type;
        using V12       = typename md::unify_types<V1,V2>::type;
        using VR        = Val_ret;

        Integer N       = cols;        

        if (A.nnz() == 0 || B.nnz() == 0)
            return;

        const sparse_ccs<V1>& Ad = A.rep();
        const sparse_ccs<V2>& Bd = B.rep();

        const Integer * A_c		= Ad.ptr_c();
        const Integer * A_r		= Ad.ptr_r();
        const V1 * A_x	        = Ad.ptr_x();

        const Integer * B_c		= Bd.ptr_c();
        const Integer * B_r		= Bd.ptr_r();
        const V2 * B_x	        = Bd.ptr_x();

        Integer C_ld            = C.ld();
        VR* ptr_C               = C.ptr() + fr + fc * C_ld;

        for (Integer j = 0; j < N; ++j)
        {
            for (Integer k = B_c[j] ; k < B_c[j + 1] ; ++k)
            {
                const V2& Bx    = B_x[k];
                Integer Br      = B_r[k];

                if (mrd::is_zero(Bx))
                    continue;

                for (Integer ka = A_c[Br]; ka < A_c[Br + 1]; ++ka)
                {
                    Integer i   = A_r[ka];
                    V12 tmp     = A_x[ka] * Bx;
                    ptr_C[i]    = alpha.eval(ptr_C[i], tmp);
                };
            };

            ptr_C               += C_ld;
        };
        
        return;
    };
};

}}};