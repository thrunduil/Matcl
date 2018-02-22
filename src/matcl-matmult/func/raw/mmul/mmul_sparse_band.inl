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

namespace matcl { namespace raw { namespace details
{

//========================================================================
//                      SPARSE - BAND
//========================================================================
template<class Val_ret, class M1, class M2>
struct eval_mult<Val_ret,M1,M2,struct_sparse,struct_banded> 
{
    using Mat_D = raw::Matrix<Val_ret, struct_dense>;

    static void eval_diag(matcl::Matrix& ret, const M1& A, const M2& B, trans_type t_A, trans_type t_B)
    {
        using VTR           = Val_ret;
        using VT1           = typename M1::value_type;
        using VT2           = typename M2::value_type;
        using ti_ret_type   = typename ti::get_ti_type<Val_ret>::type;

        ti_ret_type ret_ti  = ti::get_return_ti<ti_ret_type>(mdyf::op_mul::eval(),
                                                                  ti::get_ti(A), ti::get_ti(B));        

        Integer M, K, N;

        if (t_B == trans_type::no_trans)
            N = B.cols();
        else
            N = B.rows();

        if (t_A == trans_type::no_trans)
        {
            M   = A.rows();
            K   = A.cols();
        }
        else
        {
            K   = A.rows();
            M   = A.cols();
        };

        if (M == 0 || K == 0 || N == 0 || A.nnz() == 0)
        {
            using val_type_ret_real = typename md::real_type<VTR>::type;
            using sparse_matrix     = Matrix<val_type_ret_real,struct_sparse>;

            sparse_matrix out(ret_ti,M,N);
            ret = matcl::Matrix(out,false);
            return;
        };

        const VT2* ptr_B    = B.rep_ptr() + B.first_elem_diag(0);
        Integer B_ld        = B.ld();

        Integer K1          = std::min(N,K);

        if (t_A == trans_type::no_trans)
        {
            using ret_type  = raw::Matrix<Val_ret,struct_sparse>;
            ret_type C(ret_ti, M, N, A.nnz());

            sparse_ccs<VTR>& d          = C.rep();
            const sparse_ccs<VT1>& Ad   = A.rep();

            Integer nz          = 0;
            Integer * d_c       = d.ptr_c();
            Integer * d_r       = d.ptr_r();
            VTR * d_x           = d.ptr_x();

            const Integer* Ad_c = Ad.ptr_c();
            const Integer* Ad_r = Ad.ptr_r();
            const VT1* Ad_x     = Ad.ptr_x();            

            bool conj_B         = (t_B == trans_type::conj_trans);

            for (Integer j = 0; j < K1; ++j, ptr_B += B_ld)
            {
                d_c[j]          = nz;
                const VT2& Bx   = (conj_B)? conj(ptr_B[0]) : ptr_B[0];

                if (mrd::is_zero(Bx))
                    continue;

                for (Integer ka = Ad_c[j]; ka < Ad_c[j+1]; ++ka)
                {
                    d_x[nz]    = Ad_x[ka] * Bx;
                    d_r[nz]    = Ad_r[ka];
                    ++nz;
                };
            }

            for (Integer j = K1; j <= N; ++j)
                d_c[j] = nz;

            d.add_memory(-1);

            bool is_sq_A    = A.rows() == A.cols();
            bool is_sq_B    = B.rows() == B.cols();
            bool is_sq_C    = C.rows() == C.cols();
            C.set_struct(predefined_struct_ext::mult_struct(A.get_struct(),B.get_struct(),t_A,t_B,
                                is_real_matrix(A),is_real_matrix(B), is_sq_A, is_sq_B, is_sq_C));

            ret = matcl::Matrix(C,true);
            return;
        };

        using ret_type      = raw::Matrix<Val_ret,struct_sparse>;
        ret_type At(ret_ti);

        if (t_A == trans_type::trans)
            mrd::manip_trans_reshaper_helper<VT1, VTR>::eval_trans(At, A, K1, A.cols(), M, N);
        else
            mrd::manip_trans_reshaper_helper<VT1, VTR>::eval_ctrans(At, A, K1, A.cols(), M, N);

        sparse_ccs<VTR>& Yd = At.rep();
        Integer * Y_c       = Yd.ptr_c();
        Integer * Y_r       = Yd.ptr_r();
        VTR* Y_x            = Yd.ptr_x();

        bool conj_B         = (t_B == trans_type::conj_trans);
        Integer j           = 0;
        Integer nz          = 0;

        for (; j < K1; ++j)
        {
            const VT2& val  = conj_B? conj(ptr_B[0]) : ptr_B[0];
            nz              = Y_c[j];

            if (mrd::is_zero(val))
            {
                ptr_B       += B_ld;
                ++j;
                break;
            };

            for (Integer i = Y_c[j]; i < Y_c[j + 1]; ++i)
                Y_x[i]      = Y_x[i] * val;

            ptr_B           += B_ld;
            nz              = Y_c[j+1];
        };

        for (; j < K1; ++j)
        {
            Integer pos_s   = Y_c[j];
            Y_c[j]          = nz;

            const VT2& val  = conj_B? conj(ptr_B[0]) : ptr_B[0];

            if (mrd::is_zero(val))
            {
                ptr_B       += B_ld;
                continue;
            };

            for (Integer i = pos_s; i < Y_c[j + 1]; ++i)
            {
                Y_r[nz]     = Y_r[i];
                Y_x[nz]     = Y_x[i] * val;
                ++nz;
            };

            ptr_B           += B_ld;
        };

        for (j = K1; j <= N; ++j)
            Y_c[j]          = nz;

        bool is_sq_A    = A.rows() == A.cols();
        bool is_sq_B    = B.rows() == B.cols();
        bool is_sq_C    = At.rows() == At.cols();
        At.set_struct(predefined_struct_ext::mult_struct(A.get_struct(),B.get_struct(),t_A,t_B,
                                is_real_matrix(A),is_real_matrix(B), is_sq_A, is_sq_B, is_sq_C));

        ret = matcl::Matrix(At,true);
        return;
    };

    static void eval_diag_gemm(const Val_ret& alpha, Mat_D& C, const M1& A, const M2& B, 
                     trans_type t_A, trans_type t_B, Integer fr, Integer rows, Integer fc, Integer cols)
    {
        using VTR           = Val_ret;
        using VT1           = typename M1::value_type;
        using VT2           = typename M2::value_type;

        Integer M           = rows;
        Integer N           = cols;

        if (M == 0 || N == 0 || A.nnz() == 0)
            return;
        
        Integer B_ld        = B.ld();
        Integer C_ld        = C.ld();

        const VT2* ptr_B    = B.rep_ptr() + B.first_elem_diag(0);
        VTR* ptr_C          = C.ptr() + fr + fc * C_ld;

        Integer K1          = std::min(B.rows(),B.cols());

        const sparse_ccs<VT1>& Ad   = A.rep();

        const Integer* Ad_c = Ad.ptr_c();
        const Integer* Ad_r = Ad.ptr_r();
        const VT1* Ad_x     = Ad.ptr_x();            

        bool conj_B         = (t_B == trans_type::conj_trans);

        if (t_A == trans_type::no_trans)
        {
            for (Integer j = 0; j < K1; ++j, ptr_B += B_ld)
            {
                const VTR& Bx   = (conj_B)? alpha * conj(ptr_B[0]) : alpha * ptr_B[0];

                if (mrd::is_zero(Bx))
                {
                    ptr_C       += C_ld;
                    continue;
                };

                for (Integer ka = Ad_c[j]; ka < Ad_c[j+1]; ++ka)
                {
                    VTR val     = Ad_x[ka] * Bx;
                    Integer r   = Ad_r[ka]; 
                    
                    ptr_C[r]    = ptr_C[r] + val;
                };

                ptr_C           += C_ld;
            }
            
            return;
        };

        bool conj_A         = (t_A == trans_type::conj_trans);

        //A.cols() == M
        if (conj_A == true)
        {
            if (conj_B == true)
            {
                for (Integer j = 0; j < M; ++j)
                {
                    for (Integer ka = Ad_c[j]; ka < Ad_c[j+1]; ++ka)
                    {
                        Integer col     = Ad_r[ka]; 

                        if (col >= K1)
                            continue;

                        VTR vd          = alpha * conj(ptr_B[col * B_ld]);
                        VTR val         = conj(Ad_x[ka]) * vd;
                    
                        ptr_C[col*C_ld] = ptr_C[col*C_ld] + val;
                    };

                    ptr_C               += 1;
                };
            }
            else
            {
                for (Integer j = 0; j < M; ++j)
                {
                    for (Integer ka = Ad_c[j]; ka < Ad_c[j+1]; ++ka)
                    {
                        Integer col     = Ad_r[ka]; 

                        if (col >= K1)
                            continue;

                        VTR vd          = alpha * ptr_B[col * B_ld];
                        VTR val         = conj(Ad_x[ka]) * vd;
                    
                        ptr_C[col*C_ld] = ptr_C[col*C_ld] + val;
                    };

                    ptr_C               += 1;
                };
            };
        }
        else
        {
            if (conj_B == true)
            {
                for (Integer j = 0; j < M; ++j)
                {
                    for (Integer ka = Ad_c[j]; ka < Ad_c[j+1]; ++ka)
                    {
                        Integer col     = Ad_r[ka]; 

                        if (col >= K1)
                            continue;

                        VTR vd          = alpha * conj(ptr_B[col * B_ld]);
                        VTR val         = Ad_x[ka] * vd;
                    
                        ptr_C[col*C_ld] = ptr_C[col*C_ld] + val;
                    };

                    ptr_C               += 1;
                };
            }
            else
            {
                for (Integer j = 0; j < M; ++j)
                {
                    for (Integer ka = Ad_c[j]; ka < Ad_c[j+1]; ++ka)
                    {
                        Integer col     = Ad_r[ka]; 

                        if (col >= K1)
                            continue;

                        VTR vd          = alpha * ptr_B[col * B_ld];
                        VTR val         = Ad_x[ka] * vd;
                    
                        ptr_C[col*C_ld] = ptr_C[col*C_ld] + val;
                    };

                    ptr_C               += 1;
                };
            };
        };
    };
    
    static void eval(matcl::Matrix& ret, const M1& A, const M2& B, trans_type t_A, trans_type t_B)
    {
        bool is_diag        = raw::is_diag(B);

        if (is_diag)
            return eval_diag(ret, A, B, t_A, t_B);

        if (t_A == trans_type::no_trans)
            return eval_notrans(ret, A, B, t_A, t_B);
        else if (t_A == trans_type::trans)
            return eval_trans<false>(ret, A, B, t_A, t_B);
        else
            return eval_trans<true>(ret, A, B, t_A, t_B);
    };

    static void eval_notrans(matcl::Matrix& ret, const M1& A, const M2& B, trans_type t_A, trans_type t_B)
    {
        using VTR           = Val_ret;
        using VT1           = typename M1::value_type;
        using VT2           = typename M2::value_type;
        using ti_ret_type   = typename ti::get_ti_type<Val_ret>::type;

        ti_ret_type ret_ti  = ti::get_return_ti<ti_ret_type>(mdyf::op_mul::eval(),
                                                                  ti::get_ti(A), ti::get_ti(B));        

        Integer M, K, N;

        if (t_B == trans_type::no_trans)
            N = B.cols();
        else
            N = B.rows();

        if (t_A == trans_type::no_trans)
        {
            M   = A.rows();
            K   = A.cols();
        }
        else
        {
            K   = A.rows();
            M   = A.cols();
        };

        if (M == 0 || K == 0 || N == 0 || A.nnz() == 0)
        {
            using val_type_ret_real = typename md::real_type<VTR>::type;
            using sparse_matrix     = Matrix<val_type_ret_real,struct_sparse>;

            sparse_matrix out(ret_ti,M,N);
            ret = matcl::Matrix(out,false);
            return;
        };

        const VT2* ptr_B    = B.rep_ptr();
        Integer B_ld        = B.ld();

        using scatter       = matcl::algorithm::scatter;

        using ret_type      = raw::Matrix<Val_ret,struct_sparse>;
        ret_type C(ret_ti, M, N, M + estim_mult_nnz(A, B, t_A, t_B));

        sparse_ccs<VTR>& d          = C.rep();
        const sparse_ccs<VT1>& Ad   = A.rep();

        Integer nz          = 0;
        Integer * d_c       = d.ptr_c();
        Integer * d_r       = d.ptr_r();
        VTR * d_x           = d.ptr_x();

        const Integer* Ad_c = Ad.ptr_c();
        const Integer* Ad_r = Ad.ptr_r();
        const VT1* Ad_x     = Ad.ptr_x();

        scatter sc          = scatter::get(M, N);
        md::workspace2<VTR>  work_x(ret_ti,M);        

        if (t_B == trans_type::no_trans)
        {
            for (Integer j = 0; j < N; ++j)
            {
                if (nz + M > d.nzmax()) 
                {
                    d.add_memory( d.nzmax() + M);

                    d_r        = d.ptr_r();
                    d_x        = d.ptr_x();
                };

                d_c[j]          = nz;
                Integer nz_old  = nz;
                Integer nk      = 0;
                auto mark       = sc.next_mark();

                Integer fr      = B.first_row(j);
                Integer fe      = B.first_elem_pos(j);
                Integer lr      = B.last_row(j);

                for (Integer k = fr, pos = fe; k <= lr; ++k, ++pos)
                {
                    const VT2& Bx = ptr_B[pos];

                    if (mrd::is_zero(Bx))
                        continue;

                    Integer nzs = nz;

                    for (Integer ka = Ad_c[k]; ka < Ad_c[k+1]; ++ka)
                    {
                        Integer i = Ad_r[ka];
                        if (sc[i] < mark)
                        {
                            sc[i]       = mark;
                            work_x[i]   = Ad_x[ka] * Bx;
                            d_r[nz]     = i;
                            ++nz;
                        }
                        else 
                        {
                            VTR tmp     = Ad_x[ka] * Bx;
                            work_x[i]   = work_x[i] + tmp;
                        };
                    }

                    if (nz - nzs > 0)
                        ++nk;
                }

                if (nk > 1)
                    utils::sort_q(d_r+nz_old,nz - nz_old);

                for (Integer k = d_c[j]; k < nz ; ++k) 
                    d_x[k] = work_x[d_r[k]] ;

                ptr_B += B_ld;
            }

            d_c[N] = nz;
        }
        else
        {
            bool conj_B = (t_B == trans_type::conj_trans);

            for (Integer j = 0; j < N; ++j)
            {
                if (nz + M > d.nzmax()) 
                {
                    d.add_memory( d.nzmax() + M);

                    d_r        = d.ptr_r();
                    d_x        = d.ptr_x();
                };

                d_c[j]          = nz;
                Integer nz_old  = nz;
                Integer nk      = 0;
                auto mark       = sc.next_mark();

                Integer fr      = B.first_col(j);
                Integer fe      = B.first_elem_pos_row(j);
                Integer lr      = B.last_col(j);

                for (Integer k = fr, pos = fe; k <= lr; ++k, pos += B_ld - 1)
                {
                    const VT2& Bx = conj_B? conj(ptr_B[pos]) : ptr_B[pos];

                    if (mrd::is_zero(Bx))
                        continue;

                    Integer nzs = nz;

                    for (Integer ka = Ad_c[k]; ka < Ad_c[k+1]; ++ka)
                    {
                        Integer i = Ad_r[ka];
                        if (sc[i] < mark)
                        {
                            sc[i]       = mark;
                            work_x[i]   = Ad_x[ka] * Bx;
                            d_r[nz]     = i;
                            ++nz;
                        }
                        else 
                        {
                            VTR tmp     = Ad_x[ka] * Bx;
                            work_x[i]   = work_x[i] + tmp;
                        };
                    }

                    if (nz - nzs > 0)
                        ++nk;
                }

                if (nk > 1)
                    utils::sort_q(d_r+nz_old,nz - nz_old);

                for (Integer k = d_c[j]; k < nz ; ++k) 
                    d_x[k] = work_x[d_r[k]] ;
            }

            d_c[N] = nz;
        };

        d.add_memory(-1);

        bool is_sq_A    = A.rows() == A.cols();
        bool is_sq_B    = B.rows() == B.cols();
        bool is_sq_C    = C.rows() == C.cols();
        C.set_struct(predefined_struct_ext::mult_struct(A.get_struct(),B.get_struct(),t_A,t_B,
                                is_real_matrix(A),is_real_matrix(B), is_sq_A, is_sq_B, is_sq_C));

        ret = matcl::Matrix(C,true);
        return;
    };

    static void eval_notrans_gemm(const Val_ret& alpha, Mat_D& C, const M1& A, const M2& B, 
                     trans_type t_A, trans_type t_B, Integer C_fr, Integer rows, Integer C_fc, Integer cols)
    {
        (void)t_A;

        using VTR           = Val_ret;
        using VT1           = typename M1::value_type;
        using VT2           = typename M2::value_type;

        Integer M           = rows;
        Integer N           = cols;

        if (M == 0 || N == 0 || A.nnz() == 0)
            return;
        
        Integer B_ld        = B.ld();
        Integer C_ld        = C.ld();
        const VT2* ptr_B    = B.rep_ptr();
        VTR* ptr_C          = C.ptr() + C_fr + C_fc * C_ld;

        const sparse_ccs<VT1>& Ad   = A.rep();

        const Integer* Ad_c = Ad.ptr_c();
        const Integer* Ad_r = Ad.ptr_r();
        const VT1* Ad_x     = Ad.ptr_x();

        if (t_B == trans_type::no_trans)
        {
            for (Integer j = 0; j < N; ++j)
            {
                Integer fr      = B.first_row(j);
                Integer lr      = B.last_row(j);
                Integer fe      = B.first_elem_pos(j);                

                for (Integer k = fr, pos = fe; k <= lr; ++k, ++pos)
                {
                    VTR Bx      = alpha * ptr_B[pos];

                    if (mrd::is_zero(Bx))
                        continue;

                    for (Integer ka = Ad_c[k]; ka < Ad_c[k+1]; ++ka)
                    {
                        Integer i   = Ad_r[ka];
                        VTR val     = Ad_x[ka] * Bx;
                        ptr_C[i]    = ptr_C[i] + val;
                    }
                }

                ptr_B   += B_ld;
                ptr_C   += C_ld;
            }
        }
        else
        {
            bool conj_B = (t_B == trans_type::conj_trans);

            for (Integer j = 0; j < N; ++j)
            {
                Integer fr      = B.first_col(j);
                Integer lr      = B.last_col(j);
                Integer fe      = B.first_elem_pos_row(j);                

                for (Integer k = fr, pos = fe; k <= lr; ++k, pos += B_ld - 1)
                {
                    VTR Bx      = conj_B? alpha * conj(ptr_B[pos]) : alpha * ptr_B[pos];

                    if (mrd::is_zero(Bx))
                        continue;

                    for (Integer ka = Ad_c[k]; ka < Ad_c[k+1]; ++ka)
                    {
                        Integer i   = Ad_r[ka];
                        VTR val     = Ad_x[ka] * Bx;
                        ptr_C[i]    = ptr_C[i] + val;
                    }
                }

                ptr_C           += C_ld;
            };
        };
        return;
    };

    template<bool Conj_A>
    static void eval_trans_gemm(const Val_ret& alpha, Mat_D& C, const M1& A, const M2& B, 
                     trans_type t_A, trans_type t_B, Integer fr, Integer rows, Integer fc, Integer cols)
    {
        if (t_B == trans_type::no_trans)
            return eval_trans_notrans_gemm<Conj_A>(alpha,C,A,B,t_A,t_B,fr,rows,fc,cols);
        else if (t_B == trans_type::trans)
            return eval_trans_trans_gemm<Conj_A,false>(alpha,C,A,B,t_A,t_B,fr,rows,fc,cols);
        else
            return eval_trans_trans_gemm<Conj_A,true>(alpha,C,A,B,t_A,t_B,fr,rows,fc,cols);
    };

    template<bool Conj_A>
    static void eval_trans(matcl::Matrix& ret, const M1& A, const M2& B, trans_type t_A, trans_type t_B)
    {
        if (t_B == trans_type::no_trans)
            return eval_trans_notrans<Conj_A>(ret,A,B,t_A,t_B);
        else if (t_B == trans_type::trans)
            return eval_trans_trans<Conj_A,false>(ret,A,B,t_A,t_B);
        else
            return eval_trans_trans<Conj_A,true>(ret,A,B,t_A,t_B);
    }

    template<bool Conj_A>
    static void eval_trans_notrans_gemm(const Val_ret& alpha, Mat_D& C, const M1& A, const M2& B, 
                    trans_type t_A, trans_type t_B, Integer C_fr, Integer rows, Integer C_fc, Integer cols)
    {
        (void)t_A;
        (void)t_B;

        using VTR           = Val_ret;
        using VT1           = typename M1::value_type;
        using VT2           = typename M2::value_type;

        VTR Z               = md::default_value<VTR>(C.get_type());

        Integer M           = rows;
        Integer N           = cols;

        if (M == 0 || N == 0 || A.nnz() == 0)
            return;

        Integer B_ld        = B.ld();
        Integer C_ld        = C.ld();
        const VT2* ptr_B    = B.rep_ptr();
        VTR* ptr_C          = C.ptr() + C_fr + C_fc * C_ld;

        const sparse_ccs<VT1>& Ad   = A.rep();
        const Integer* Ad_c = Ad.ptr_c();
        const Integer* Ad_r = Ad.ptr_r();
        const VT1* Ad_x     = Ad.ptr_x();

        for (Integer j = 0; j < N; ++j)
        {
            Integer fr      = B.first_row(j);
            Integer lr      = B.last_row(j);
            Integer fe      = B.first_elem_pos(j)-fr;            

            const VT2* ptr_B2   = ptr_B + fe;

            for (Integer i = 0; i < M; ++i)
            {
                Integer last    = Ad_c[i + 1];
                Integer first   = Ad_c[i];

                if (last == first)
                    continue;

                VTR dot         = Z;                

                for (Integer ka = first; ka < last; ++ka)
                {
                    Integer Ar  = Ad_r[ka];

                    if (Ar >= fr && Ar <= lr)
                    {
                        dot = dot + make_conj<Conj_A>::eval(Ad_x[ka]) * ptr_B2[Ar];
                    };
                };

                ptr_C[i]        = ptr_C[i] + alpha * dot;
            };

            ptr_B   += B_ld;
            ptr_C   += C_ld;
        }

        return;
    };

    template<bool Conj_A>
    static void eval_trans_notrans(matcl::Matrix& ret, const M1& A, const M2& B, trans_type t_A, trans_type t_B)
    {
        using VTR           = Val_ret;
        using VT1           = typename M1::value_type;
        using VT2           = typename M2::value_type;
        using ti_ret_type   = typename ti::get_ti_type<Val_ret>::type;

        ti_ret_type ret_ti  = ti::get_return_ti<ti_ret_type>(mdyf::op_mul::eval(),
                                                                  ti::get_ti(A), ti::get_ti(B));        

        VTR Z               = md::default_value<VTR>(ret_ti);

        Integer M, K, N;

        if (t_B == trans_type::no_trans)
            N = B.cols();
        else
            N = B.rows();

        if (t_A == trans_type::no_trans)
        {
            M   = A.rows();
            K   = A.cols();
        }
        else
        {
            K   = A.rows();
            M   = A.cols();
        };

        if (M == 0 || K == 0 || N == 0 || A.nnz() == 0)
        {
            using val_type_ret_real = typename md::real_type<VTR>::type;
            using sparse_matrix     = Matrix<val_type_ret_real,struct_sparse>;

            sparse_matrix out(ret_ti,M,N);
            ret = matcl::Matrix(out,false);
            return;
        };

        const VT2* ptr_B    = B.rep_ptr();
        Integer B_ld        = B.ld();

        using ret_type  = raw::Matrix<Val_ret,struct_sparse>;
        ret_type C(ret_ti, M, N, M + estim_mult_nnz(A, B, t_A, t_B));

        sparse_ccs<VTR>& d          = C.rep();
        const sparse_ccs<VT1>& Ad   = A.rep();

        Integer nz          = 0;
        Integer * d_c       = d.ptr_c();
        Integer * d_r       = d.ptr_r();
        VTR * d_x           = d.ptr_x();

        const Integer* Ad_c = Ad.ptr_c();
        const Integer* Ad_r = Ad.ptr_r();
        const VT1* Ad_x     = Ad.ptr_x();

        for (Integer j = 0; j < N; ++j)
        {
            if (nz + M > d.nzmax()) 
            {
                d.add_memory( d.nzmax() + M);

                d_r        = d.ptr_r();
                d_x        = d.ptr_x();
            };

            d_c[j]          = nz;

            Integer fr      = B.first_row(j);
            Integer fe      = B.first_elem_pos(j)-fr;
            Integer lr      = B.last_row(j);

            const VT2* ptr_B2   = ptr_B + fe;

            for (Integer i = 0; i < M; ++i)
            {
                Integer last    = Ad_c[i + 1];
                Integer first   = Ad_c[i];

                if (last == first)
                    continue;

                VTR dot         = Z;                

                for (Integer ka = first; ka < last; ++ka)
                {
                    Integer Ar  = Ad_r[ka];

                    if (Ar >= fr && Ar <= lr)
                    {
                        dot = dot + make_conj<Conj_A>::eval(Ad_x[ka]) * ptr_B2[Ar];
                    };
                };

                if (mrd::is_zero(dot) == false)
                {
                    d_r[nz]     = i;
                    d_x[nz]     = dot;
                    ++nz;
                };
            };

            ptr_B += B_ld;
        }

        d_c[N] = nz;

        d.add_memory(-1);

        bool is_sq_A    = A.rows() == A.cols();
        bool is_sq_B    = B.rows() == B.cols();
        bool is_sq_C    = C.rows() == C.cols();
        C.set_struct(predefined_struct_ext::mult_struct(A.get_struct(),B.get_struct(),t_A,t_B,
                            is_real_matrix(A),is_real_matrix(B), is_sq_A, is_sq_B, is_sq_C));

        ret = matcl::Matrix(C,true);
        return;
    };

    template<bool Conj_A, bool Conj_B>
    static void eval_trans_trans_gemm(const Val_ret& alpha, Mat_D& C, const M1& A, const M2& B, 
                    trans_type t_A, trans_type t_B, Integer C_fr, Integer rows, Integer C_fc, Integer cols)
    {
        (void)t_B;
        (void)t_A;

        using VTR           = Val_ret;
        using VT1           = typename M1::value_type;
        using VT2           = typename M2::value_type;

        VTR Z               = md::default_value<VTR>(C.get_type());

        Integer M           = rows;
        Integer N           = cols;

        if (rows == 0 || cols == 0 || A.nnz() == 0)
            return;

        Integer B_ld        = B.ld();
        Integer C_ld        = C.ld();

        const VT2* ptr_B    = B.rep_ptr();        
        VTR* ptr_C          = C.ptr() + C_fr + C_fc * C_ld;

        const sparse_ccs<VT1>& Ad   = A.rep();
        const Integer* Ad_c = Ad.ptr_c();
        const Integer* Ad_r = Ad.ptr_r();
        const VT1* Ad_x     = Ad.ptr_x();
        Integer B_ldm1      = B_ld - 1;

        for (Integer j = 0; j < N; ++j)
        {
            Integer fr      = B.first_col(j);
            Integer fe      = B.first_elem_pos_row(j) - fr*B_ldm1;
            Integer lr      = B.last_col(j);

            const VT2* ptr_B2   = ptr_B + fe;

            for (Integer i = 0; i < M; ++i)
            {
                Integer last    = Ad_c[i + 1];
                Integer first   = Ad_c[i];

                if (last == first)
                    continue;

                VTR dot         = Z;                

                for (Integer ka = first; ka < last; ++ka)
                {
                    Integer Ar  = Ad_r[ka];

                    if (Ar >= fr && Ar <= lr)
                    {
                        dot = dot + make_conj<Conj_A>::eval(Ad_x[ka]) * make_conj<Conj_B>::eval(ptr_B2[Ar*B_ldm1]);
                    };
                };

                ptr_C[i]        = ptr_C[i] + alpha * dot;
            };

            ptr_C   += C_ld;
        }

        return;
    };

    template<bool Conj_A, bool Conj_B>
    static void eval_trans_trans(matcl::Matrix& ret, const M1& A, const M2& B, trans_type t_A, trans_type t_B)
    {
        using VTR           = Val_ret;
        using VT1           = typename M1::value_type;
        using VT2           = typename M2::value_type;
        using ti_ret_type   = typename ti::get_ti_type<Val_ret>::type;

        ti_ret_type ret_ti  = ti::get_return_ti<ti_ret_type>(mdyf::op_mul::eval(),
                                                                  ti::get_ti(A), ti::get_ti(B));        

        VTR Z               = md::default_value<VTR>(ret_ti);

        Integer M, K, N;

        if (t_B == trans_type::no_trans)
            N = B.cols();
        else
            N = B.rows();

        if (t_A == trans_type::no_trans)
        {
            M   = A.rows();
            K   = A.cols();
        }
        else
        {
            K   = A.rows();
            M   = A.cols();
        };

        if (M == 0 || K == 0 || N == 0 || A.nnz() == 0)
        {
            using val_type_ret_real = typename md::real_type<VTR>::type;
            using sparse_matrix     = Matrix<val_type_ret_real,struct_sparse>;

            sparse_matrix out(ret_ti,M,N);
            ret = matcl::Matrix(out,false);
            return;
        };

        const VT2* ptr_B    = B.rep_ptr();
        Integer B_ld        = B.ld();

        using ret_type  = raw::Matrix<Val_ret,struct_sparse>;
        ret_type C(ret_ti, M, N, M + estim_mult_nnz(A, B, t_A, t_B));

        sparse_ccs<VTR>& d          = C.rep();
        const sparse_ccs<VT1>& Ad   = A.rep();

        Integer nz          = 0;
        Integer * d_c       = d.ptr_c();
        Integer * d_r       = d.ptr_r();
        VTR * d_x           = d.ptr_x();

        const Integer* Ad_c = Ad.ptr_c();
        const Integer* Ad_r = Ad.ptr_r();
        const VT1* Ad_x     = Ad.ptr_x();
        Integer B_ldm1      = B_ld - 1;

        for (Integer j = 0; j < N; ++j)
        {
            if (nz + M > d.nzmax()) 
            {
                d.add_memory( d.nzmax() + M);

                d_r        = d.ptr_r();
                d_x        = d.ptr_x();
            };

            d_c[j]          = nz;

            Integer fr      = B.first_col(j);
            Integer fe      = B.first_elem_pos_row(j) - fr*B_ldm1;
            Integer lr      = B.last_col(j);

            const VT2* ptr_B2   = ptr_B + fe;

            for (Integer i = 0; i < M; ++i)
            {
                Integer last    = Ad_c[i + 1];
                Integer first   = Ad_c[i];

                if (last == first)
                    continue;

                VTR dot         = Z;                

                for (Integer ka = first; ka < last; ++ka)
                {
                    Integer Ar  = Ad_r[ka];

                    if (Ar >= fr && Ar <= lr)
                    {
                        dot = dot + make_conj<Conj_A>::eval(Ad_x[ka]) * make_conj<Conj_B>::eval(ptr_B2[Ar*B_ldm1]);
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

        d_c[N] = nz;

        d.add_memory(-1);

        bool is_sq_A    = A.rows() == A.cols();
        bool is_sq_B    = B.rows() == B.cols();
        bool is_sq_C    = C.rows() == C.cols();
        C.set_struct(predefined_struct_ext::mult_struct(A.get_struct(),B.get_struct(),t_A,t_B,
                            is_real_matrix(A),is_real_matrix(B), is_sq_A, is_sq_B, is_sq_C));

        ret = matcl::Matrix(C,true);
        return;
    };

    static void eval_gemm(Mat_D& C, const Val_ret& alpha, const M1& A, const M2& B, 
                     trans_type t_A, trans_type t_B, const Val_ret& beta, 
                     Integer fr, Integer rows, Integer fc, Integer cols)
    {
        bool is_diag        = raw::is_diag(B);

        prepare_gemm_C<Val_ret>::eval(beta, C, fr,rows,fc,cols);

        if (is_diag)
            return eval_diag_gemm(alpha, C, A, B, t_A, t_B, fr, rows, fc, cols);

        if (t_A == trans_type::no_trans)
            return eval_notrans_gemm(alpha, C, A, B, t_A, t_B, fr, rows, fc, cols);
        else if (t_A == trans_type::trans)
            return eval_trans_gemm<false>(alpha, C, A, B, t_A, t_B, fr, rows, fc, cols);
        else
            return eval_trans_gemm<true>(alpha, C, A, B, t_A, t_B, fr, rows, fc, cols);
    };
};

}}};