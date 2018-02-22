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

#include "matcl-matrep/utils/workspace.h"

namespace matcl { namespace raw { namespace details
{

//========================================================================
//                      DENSE - BAND
//========================================================================

template<class V>
struct eval_dense_band_diag
{
    using M1    = raw::Matrix<V,struct_dense>;
    using M2    = raw::Matrix<V,struct_banded>;

    static void eval_diag_diag(matcl::Matrix& ret, const M1& A, const M2& B, 
                               trans_type t_A, trans_type t_B)
    {
        Integer N, M, K;

        if (t_B == trans_type::no_trans)
            N       = B.cols();
        else
            N       = B.rows();

        if (t_A == trans_type::no_trans)
        {
            M       = A.rows();
            K       = A.cols();
        }
        else
        {
            K       = A.rows();
            M       = A.cols();
        };

        Integer MN  = std::min(M,N);
        K           = std::min(MN,K);

        ti::ti_type<V> ret_ti = ti::get_return_ti<ti::ti_type<V>>(mdyf::op_mul::eval(),
                                                                  ti::get_ti(A), ti::get_ti(B));
        V Z = md::default_value<V>(ret_ti);

        M1 C(ret_ti, Z, MN, 1);

        const V* A_ptr  = A.ptr();
        const V* B_ptr  = &B(1,1);
        V* C_ptr        = C.ptr();

        Integer A_ld    = A.ld();
        Integer B_ld    = B.ld();

        if (t_A == trans_type::conj_trans)
        {
            if (t_B == trans_type::conj_trans)
            {
                for (Integer i = 0; i < K; ++i)
                {
                    C_ptr[i]    = conj(A_ptr[0]) * conj(B_ptr[0]);
                    A_ptr       += A_ld+1;
                    B_ptr       += B_ld;
                };
            }
            else
            {
                for (Integer i = 0; i < K; ++i)
                {
                    C_ptr[i]    = conj(A_ptr[0]) * B_ptr[0];
                    A_ptr       += A_ld+1;
                    B_ptr       += B_ld;
                };
            };
        }
        else
        {
            if (t_B == trans_type::conj_trans)
            {
                for (Integer i = 0; i < K; ++i)
                {
                    C_ptr[i]    = A_ptr[0] * conj(B_ptr[0]);
                    A_ptr       += A_ld+1;
                    B_ptr       += B_ld;
                };
            }
            else
            {
                for (Integer i = 0; i < K; ++i)
                {
                    C_ptr[i]    = A_ptr[0] * B_ptr[0];
                    A_ptr       += A_ld+1;
                    B_ptr       += B_ld;
                };
            };
        };

        ret = matcl::bdiags(matcl::Matrix(C,false), Integer(0), M, N);
    };

    static void eval_diag_diag_gemm(const V& alpha, M1& C, const M1& A, const M2& B, 
                trans_type t_A, trans_type t_B, Integer fr, Integer rows, Integer fc, Integer cols)
    {
        Integer MN  = std::min(rows,cols);
        Integer K   = std::min(A.rows(),A.cols());
        K           = std::min(MN,K);

        Integer A_ld    = A.ld();
        Integer B_ld    = B.ld();
        Integer C_ld    = C.ld();

        const V* A_ptr  = A.ptr();
        const V* B_ptr  = &B(1,1);
        V* C_ptr        = C.ptr() + fr + fc * C_ld;

        if (t_A == trans_type::conj_trans)
        {
            if (t_B == trans_type::conj_trans)
            {
                for (Integer i = 0; i < K; ++i)
                {
                    V tmp       = conj(A_ptr[0]) * conj(B_ptr[0]);
                    C_ptr[0]    = C_ptr[0] + alpha * tmp;
                    A_ptr       += A_ld+1;
                    B_ptr       += B_ld;
                    C_ptr       += C_ld+1;
                };
            }
            else
            {
                for (Integer i = 0; i < K; ++i)
                {
                    V tmp       = conj(A_ptr[0]) * B_ptr[0];
                    C_ptr[0]    = C_ptr[0] + alpha * tmp;
                    A_ptr       += A_ld+1;
                    B_ptr       += B_ld;
                    C_ptr       += C_ld+1;
                };
            };
        }
        else
        {
            if (t_B == trans_type::conj_trans)
            {
                for (Integer i = 0; i < K; ++i)
                {
                    V tmp       = A_ptr[0] * conj(B_ptr[0]);
                    C_ptr[0]    = C_ptr[0] + alpha * tmp;
                    A_ptr       += A_ld+1;
                    B_ptr       += B_ld;
                    C_ptr       += C_ld+1;
                };
            }
            else
            {
                for (Integer i = 0; i < K; ++i)
                {
                    V tmp       = A_ptr[0] * B_ptr[0];
                    C_ptr[0]    = C_ptr[0] + alpha * tmp;
                    A_ptr       += A_ld+1;
                    B_ptr       += B_ld;
                    C_ptr       += C_ld+1;
                };
            };
        };
    };

    static void eval1_gemm(const V& alpha, M1& C, const M1& A, const M2& B, trans_type t_A, trans_type t_B, 
                     Integer fr, Integer rows, Integer fc, Integer cols)
    {
        if (raw::is_diag(B))
            return eval_diag_diag_gemm(alpha, C, A, B, t_A, t_B, fr, rows, fc, cols);

        Integer K1  = std::min(A.rows(),A.cols());

        M1 diag(A.get_type(), K1, 1);

        Integer A_ld    = A.ld();
        Integer B_ld    = B.ld();
        Integer C_ld    = C.ld();
        Integer N       = cols;

        const V* A_ptr  = A.ptr();
        V* D_ptr        = diag.ptr();        
        const V* B_ptr  = B.rep_ptr();
        V* C_ptr        = C.ptr() + fr + fc*C_ld;

        if (t_A == trans_type::conj_trans)
        {
            for(Integer i = 0; i < K1; ++i)
            {
                D_ptr[i]    = alpha*mrd::conj_helper<V>::eval(A_ptr[0]);
                A_ptr       += A_ld + 1;
            };
        }
        else
        {
            for(Integer i = 0; i < K1; ++i)
            {
                D_ptr[i]    = alpha*A_ptr[0];
                A_ptr       += A_ld + 1;
            };
        };

        if (t_B == trans_type::no_trans)
        {
            for (Integer j = 0; j < N; ++j) 
            {	
                Integer first_row   = B.first_row(j);
                Integer last_row    = std::min(B.last_row(j),K1-1);

                Integer pos_B       = B.first_elem_pos(j);

                for (Integer l = first_row; l <= last_row; ++l, ++pos_B) 
                {
                    C_ptr[l]        = C_ptr[l] + D_ptr[l] * B_ptr[pos_B];
                };

                C_ptr       += C_ld;
                B_ptr       += B_ld;
            };
        }
        else if (t_B == trans_type::trans)
        {
            Integer Br              = B.rows();

            for (Integer j = 0; j < Br; ++j) 
            {	
                Integer first_col   = B.first_col(j);
                Integer last_col    = std::min(K1-1, B.last_col(j));

                Integer pos_B       = B.first_elem_pos_row(j);

                for (Integer l = first_col; l <= last_col; ++l) 
                {
                    C_ptr[l]        = C_ptr[l] + D_ptr[l] * B_ptr[pos_B];
                    pos_B           += B_ld - 1;
                };

                C_ptr       += C_ld;
            };
        }
        else
        {
            Integer Br              = B.rows();

            for (Integer j = 0; j < Br; ++j) 
            {	
                Integer first_col   = B.first_col(j);
                Integer last_col    = std::min(K1-1, B.last_col(j));

                Integer pos_B       = B.first_elem_pos_row(j);

                for (Integer l = first_col; l <= last_col; ++l) 
                {
                    C_ptr[l]        = C_ptr[l] + D_ptr[l] * conj(B_ptr[pos_B]);
                    pos_B           += B_ld - 1;
                };

                C_ptr       += C_ld;
            };
        };
    }

    static void eval1(matcl::Matrix& ret, const M1& A, const M2& B, trans_type t_A, trans_type t_B)
    {        
        if (raw::is_diag(B))
            return eval_diag_diag(ret, A,B, t_A, t_B);

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

        Integer K1  = std::min(M,K);

        M1 diag(A.get_type(), K1, 1);

        ti::ti_type<V> ret_ti = ti::get_return_ti<ti::ti_type<V>>(mdyf::op_mul::eval(),
                                                                  ti::get_ti(A), ti::get_ti(B));
        V Z = md::default_value<V>(ret_ti);

        Integer Bfd     = B.first_diag();
        Integer Bld     = B.last_diag();
        Integer Cfd     = (t_B == trans_type::no_trans)? Bfd : -Bld;
        Integer Cld     = (t_B == trans_type::no_trans)? Bld : -Bfd;

        M2 C(ret_ti, Z, M, N, Cfd, Cld);

        Cfd             = C.first_diag();
        Cld             = C.last_diag();

        const V* A_ptr  = A.ptr();
        V* D_ptr        = diag.ptr();
        V* C_ptr        = C.rep_ptr();
        const V* B_ptr  = B.rep_ptr();

        Integer A_ld    = A.ld();
        Integer B_ld    = B.ld();
        Integer C_ld    = C.ld();

        if (t_A == trans_type::conj_trans)
        {
            for(Integer i = 0; i < K1; ++i)
            {
                D_ptr[i]    = mrd::conj_helper<V>::eval(A_ptr[0]);
                A_ptr       += A_ld + 1;
            };
        }
        else
        {
            for(Integer i = 0; i < K1; ++i)
            {
                D_ptr[i]    = A_ptr[0];
                A_ptr       += A_ld + 1;
            };
        };

        if (t_B == trans_type::no_trans)
        {
            for (Integer j = 0; j < N; ++j) 
            {	
                Integer first_row   = B.first_row(j);
                Integer last_row    = std::min(B.last_row(j),K1-1);

                Integer pos_B       = B.first_elem_pos(j);
                Integer pos_C       = C.first_elem_pos(j);

                for (Integer l = first_row; l <= last_row; ++l, ++pos_B, ++pos_C) 
                {
                    C_ptr[pos_C]    = D_ptr[l] * B_ptr[pos_B];
                };

                C_ptr       += C_ld;
                B_ptr       += B_ld;
            };
        }
        else if (t_B == trans_type::trans)
        {
            for (Integer d = Bfd; d <= Bld; ++d)
            {
                if (C.has_diag(-d) == false)
                    continue;

                Integer pos_B   = B.first_elem_diag(d);
                Integer pos_C   = C.first_elem_diag(-d);
                Integer s1      = B.diag_length(d);
                Integer s2      = C.diag_length(-d);
                Integer s       = std::min(s1,s2);
                Integer r       = C.first_row_on_diag(-d);

                for (Integer i = 0; i < s; ++i, ++r)
                {
                    C_ptr[pos_C]    = D_ptr[r] * B_ptr[pos_B];
                    pos_B       += B_ld;
                    pos_C       += C_ld;
                };
            };
        }
        else
        {
            for (Integer d = Bfd; d <= Bld; ++d)
            {
                if (C.has_diag(-d) == false)
                    continue;

                Integer pos_B   = B.first_elem_diag(d);
                Integer pos_C   = C.first_elem_diag(-d);
                Integer s1      = B.diag_length(d);
                Integer s2      = C.diag_length(-d);
                Integer s       = std::min(s1,s2);
                Integer r       = C.first_row_on_diag(-d);

                for (Integer i = 0; i < s; ++i, ++r)
                {
                    C_ptr[pos_C]    = D_ptr[r] * conj(B_ptr[pos_B]);
                    pos_B       += B_ld;
                    pos_C       += C_ld;
                };
            };
        };

        ret = matcl::Matrix(C, true);
    };

    static void eval2(matcl::Matrix& ret, const M1& A, const M2& B, trans_type t_A, trans_type t_B)
    {
        if (raw::is_diag(A))
            return eval_diag_diag(ret,A,B, t_A, t_B);

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

        ti::ti_type<V> ret_ti = ti::get_return_ti<ti::ti_type<V>>(mdyf::op_mul::eval(),
                                                                  ti::get_ti(A), ti::get_ti(B));
        V Z = md::default_value<V>(ret_ti);

        M1 C(ret_ti, Z, M, N);

        Integer fr      = 0;
        Integer rows    = M;
        Integer fc      = 0;
        Integer cols    = N;
        V alpha         = V(1.0);

        eval2_impl(alpha, C, A, B, t_A, t_B, fr, rows, fc, cols);

        bool is_sq_A    = (A.rows() == A.cols());
        bool is_sq_B    = (B.rows() == B.cols());
        bool is_sq_C    = (C.rows() == C.cols());
        C.set_struct(predefined_struct_ext::mult_struct(A.get_struct(),B.get_struct(),t_A,t_B,
                                is_real_matrix(A),is_real_matrix(B), is_sq_A, is_sq_B, is_sq_C));

        ret = matcl::Matrix(C,true);
    };

    static void eval2_gemm(const V& alpha, M1& C, const M1& A, const M2& B, trans_type t_A, trans_type t_B, 
                     Integer fr, Integer rows, Integer fc, Integer cols)
    {
        if (raw::is_diag(A))
            return eval_diag_diag_gemm(alpha, C, A, B, t_A, t_B, fr, rows, fc, cols);

        eval2_impl(alpha, C, A, B, t_A, t_B, fr, rows, fc, cols);
    };

    static void eval2_impl(const V& alpha, M1& C, const M1& A, const M2& B, trans_type t_A, trans_type t_B,
                           Integer fr, Integer rows, Integer fc, Integer cols)
    {
        (void)cols;

        Integer K1      = std::min(B.rows(),B.cols());
        bool b_conj     = (t_B == trans_type::conj_trans);

        const V* A_ptr  = A.ptr();
        const V* B_ptr  = B.rep_ptr() + B.first_elem_diag(0);        

        Integer A_ld    = A.ld();
        Integer B_ld    = B.ld();
        Integer C_ld    = C.ld();

        Integer M       = rows;
        V* C_ptr        = C.ptr() + fr + fc * C_ld;

        if (t_A == trans_type::no_trans)
        {
            for (Integer j = 0; j < K1; ++j)
            {
                const V& tmp    = b_conj ? conj(B_ptr[0]) : B_ptr[0];

                if (mrd::is_zero(tmp))
                {
                    C_ptr   += C_ld;
                    A_ptr   += A_ld;
                    B_ptr   += B_ld;
                    continue;
                };

                for (Integer i = 0; i < M; ++i)
                {
                    V val       = A_ptr[i] * tmp;
                    C_ptr[i]    = C_ptr[i] + alpha * val;
                }

                C_ptr       += C_ld;
                A_ptr       += A_ld;
                B_ptr       += B_ld;
            };
        }
        else if (t_A == trans_type::trans)
        {
            const V* As_ptr = A.ptr();

            for (Integer j = 0; j < K1; ++j) 
            {	
                const V& tmp = b_conj ? conj(B_ptr[0]) : B_ptr[0];

                if (mrd::is_zero(tmp))
                {
                    C_ptr   += C_ld;
                    B_ptr   += B_ld;
                    continue;
                };

                A_ptr = As_ptr;
                for (Integer l = 0; l < M; ++l) 
                {
                    C_ptr[l] = C_ptr[l] + alpha * (A_ptr[j] * tmp);
                    A_ptr   += A_ld;
                };

                B_ptr += B_ld;
                C_ptr += C_ld;
            };
        }
        else
        {
            const V* As_ptr = A.ptr();

            for (Integer j = 0; j < K1; ++j) 
            {	
                const V& tmp = b_conj ? conj(B_ptr[0]) : B_ptr[0];

                if (mrd::is_zero(tmp))
                {
                    C_ptr   += C_ld;
                    B_ptr   += B_ld;
                    continue;
                };

                A_ptr = As_ptr;
                for (Integer l = 0; l < M; ++l) 
                {
                    V a         = mrd::conj_helper<V>::eval(A_ptr[j]);
                    C_ptr[l]    = C_ptr[l] + alpha * (a * tmp);
                    A_ptr       += A_ld;
                };

                B_ptr += B_ld;
                C_ptr += C_ld;
            };
        };

        return;
    };
};

template<class Val_ret>
struct eval_dense_band{};

template<class V>
struct eval_dense_band_generic
{
    using M1    = raw::Matrix<V,struct_dense>;
    using M2    = raw::Matrix<V,struct_banded>;
    
    static void eval_gemm(M1& C, const V& alpha, const M1& A, const M2& B, 
                     trans_type t_A, trans_type t_B, const V& beta, 
                     Integer fr, Integer rows, Integer fc, Integer cols)
    {        
        if (rows == 0 || cols == 0)
            return;

        prepare_gemm_C<V>::eval(beta, C, fr,rows,fc,cols);

        if (A.rows() == 0 || A.cols() == 0 || B.rows() == 0 || B.cols() == 0)
            return;

        bool is_diag_A  = raw::is_diag(A);
        bool is_diag_B  = raw::is_diag(B);

        if (is_diag_A)
            return eval_dense_band_diag<V>::eval1_gemm(alpha, C, A, B, t_A, t_B, fr, rows, fc, cols);

        if (is_diag_B)
            return eval_dense_band_diag<V>::eval2_gemm(alpha, C, A, B, t_A, t_B, fr, rows, fc, cols);

        if (mrd::is_one(alpha) == true)
        {
            using Alpha     = alpha_one<V>;
            Alpha al;

            eval_impl(al, C, A, B, t_A, t_B, fr, rows, fc, cols);
        }
        else if (mrd::is_one(-alpha) == true)
        {
            using Alpha     = alpha_mone<V>;
            Alpha al;

            eval_impl(al, C, A, B, t_A, t_B, fr, rows, fc, cols);
        }
        else
        {
            using Alpha     = alpha_val<V>;
            Alpha al(alpha);

            eval_impl(al, C, A, B, t_A, t_B, fr, rows, fc, cols);
        };

        return;
    };

    static void eval(matcl::Matrix& ret, const M1& A, const M2& B, trans_type t_A, trans_type t_B)
    {
        ti::ti_type<V> ret_ti = ti::get_return_ti<ti::ti_type<V>>(mdyf::op_mul::eval(),
                                                                  ti::get_ti(A), ti::get_ti(B));

        Integer M, K, N;

        if (t_B == trans_type::no_trans)
            N = B.cols();
        else
            N = B.rows();

        if (t_A == trans_type::no_trans)
        {
            M = A.rows();
            K = A.cols();
        }
        else
        {
            K = A.rows();
            M = A.cols();
        }

        if (M == 0 || K == 0 || N == 0)
        {
            using sparse_matrix = Matrix<V,struct_sparse>;
            sparse_matrix out(ret_ti,M, N);
            ret = matcl::Matrix(out,false);
            return;
        };		

        V Z = md::default_value<V>(ret_ti);

        if (raw::is_diag(A))
            return eval_dense_band_diag<V>::eval1(ret,A,B, t_A, t_B);

        if (raw::is_diag(B))
            return eval_dense_band_diag<V>::eval2(ret,A,B, t_A, t_B);

        M1 C(ret_ti, Z, M, N);

        Integer fr      = 0;
        Integer rows    = M;
        Integer fc      = 0;
        Integer cols    = N;

        using Alpha     = alpha_one<V>;
        Alpha al;

        eval_impl(al, C, A, B, t_A, t_B, fr, rows, fc, cols);

        bool is_full    = (fr == 0 && fc == 0 && rows == C.rows() && cols == C.cols());
        bool is_sq_A    = (A.rows() == A.cols());
        bool is_sq_B    = (B.rows() == B.cols());
        bool is_sq_C    = (C.rows() == C.cols());

        if (is_full == true)
            C.set_struct(predefined_struct_ext::mult_struct(A.get_struct(),B.get_struct(),t_A,t_B,
                                    is_real_matrix(A),is_real_matrix(B), is_sq_A, is_sq_B, is_sq_C));

        ret = matcl::Matrix(C,true);
        return;
    };

    template<class Alpha>
    static void eval_impl(const Alpha& alpha, M1& C, const M1& A, const M2& B, trans_type t_A, trans_type t_B,
                          Integer C_fr, Integer rows, Integer C_fc, Integer cols)
    {        
        const V* ptr_A;
        const V* ptr_B  = B.rep_ptr();

        Integer A_ld    = A.ld();
        Integer B_ld    = B.ld();
        Integer C_ld    = C.ld();

        V* ptr_C        = C.ptr() + C_fr + C_fc * C_ld;
        V Z             = md::default_value<V>(C.get_type());

        Integer N       = cols;
        Integer M       = rows;

        if (t_A == trans_type::no_trans)
        {
            if (t_B == trans_type::no_trans)
            {
                for (Integer j = 0; j < N; ++j) 
                {	
                    Integer fr      = B.first_row(j);
                    Integer lr      = B.last_row(j);
                    Integer fe      = B.first_elem_pos(j);
                    Integer pos_A   = imult(fr,A.ld());

                    ptr_A           = A.ptr() + pos_A;

                    for (Integer l = fr, lp = fe; l <= lr; ++l, ++lp) 
                    {
                        const V& temp = ptr_B[lp];

                        if (mrd::is_zero(temp) == false) 
                        {	
                            for (Integer k = 0; k < M; ++k) 
                            {
                                V tmp       = temp * ptr_A[k];
                                ptr_C[k]    = alpha.eval(ptr_C[k],tmp);
                            };
                        };

                        ptr_A += A_ld;
                    };

                    ptr_B += B_ld;
                    ptr_C += C_ld;
                };
            }
            else
            {
                bool conj_B         = (t_B == trans_type::conj_trans);

                for (Integer j = 0; j < N; ++j) 
                {	
                    Integer fc      = B.first_col(j);
                    Integer lc      = B.last_col(j);
                    Integer fe      = B.first_elem_pos_row(j);
                    Integer pos_A   = imult(fc,A.ld());

                    ptr_A           = A.ptr() + pos_A;

                    for (Integer l = fc, lp = fe; l <= lc; ++l) 
                    {
                        const V& temp = conj_B? conj(ptr_B[lp]) : ptr_B[lp];

                        if (mrd::is_zero(temp) == false) 
                        {	
                            for (Integer k = 0; k < M; ++k) 
                            {
                                V tmp       = temp * ptr_A[k];
                                ptr_C[k]    = alpha.eval(ptr_C[k],tmp);
                            };
                        };

                        ptr_A   += A_ld;
                        lp      += B_ld - 1;
                    };
                    
                    ptr_C += C_ld;
                };
            }
        }
        else if (t_A == trans_type::trans)
        {
            if (t_B == trans_type::no_trans)
            {
                const V* ptr_As = A.ptr();

                for (Integer j = 0; j < N; ++j) 
                {			
                    ptr_A       = ptr_As;
                    Integer fr  = B.first_row(j);
                    Integer lr  = B.last_row(j);
                    Integer fe  = B.first_elem_pos(j);

                    for (Integer l = 0; l < M; ++l) 
                    {
                        V dot = Z;		
                        for (Integer k = fr, lp = fe; k <= lr; ++k, ++lp) 
                        {
                            dot = dot + ptr_A[k] * ptr_B[lp];
                        };

                        ptr_C[l] = alpha.eval(ptr_C[l], dot);
                        ptr_A   += A_ld;
                    };

                    ptr_B += B_ld;
                    ptr_C += C_ld;
                };
            }
            else if (t_B == trans_type::trans)
            {
                const V* ptr_As = A.ptr();

                for (Integer j = 0; j < N; ++j) 
                {			
                    ptr_A       = ptr_As;
                    Integer fc  = B.first_col(j);
                    Integer lc  = B.last_col(j);
                    Integer fe  = B.first_elem_pos_row(j);

                    for (Integer l = 0; l < M; ++l) 
                    {
                        V dot = Z;		
                        for (Integer k = fc, lp = fe; k <= lc; ++k) 
                        {
                            dot = dot + ptr_A[k] * ptr_B[lp];
                            lp  += B_ld - 1;
                        };

                        ptr_C[l] = alpha.eval(ptr_C[l], dot);
                        ptr_A   += A_ld;
                    };

                    ptr_C += C_ld;
                };
            }
            else
            {
                const V* ptr_As = A.ptr();

                for (Integer j = 0; j < N; ++j) 
                {			
                    ptr_A       = ptr_As;
                    Integer fc  = B.first_col(j);
                    Integer lc  = B.last_col(j);
                    Integer fe  = B.first_elem_pos_row(j);

                    for (Integer l = 0; l < M; ++l) 
                    {
                        V dot = Z;		
                        for (Integer k = fc, lp = fe; k <= lc; ++k) 
                        {
                            dot = dot + ptr_A[k] * conj(ptr_B[lp]);
                            lp  += B_ld - 1;
                        };

                        ptr_C[l] = alpha.eval(ptr_C[l], dot);
                        ptr_A   += A_ld;
                    };

                    ptr_C += C_ld;
                };
            };
        }
        else
        {
            if (t_B == trans_type::no_trans)
            {
                const V* ptr_As = A.ptr();

                for (Integer j = 0; j < N; ++j) 
                {			
                    ptr_A       = ptr_As;
                    Integer fr  = B.first_row(j);
                    Integer lr  = B.last_row(j);
                    Integer fe  = B.first_elem_pos(j);

                    for (Integer l = 0; l < M; ++l) 
                    {
                        V dot = Z;		
                        for (Integer k = fr, lp = fe; k <= lr; ++k, ++lp) 
                        {
                            V tmp   = mrd::conj_helper<V>::eval(ptr_A[k]);
                            dot     = dot + tmp * ptr_B[lp];
                        };

                        ptr_C[l] = alpha.eval(ptr_C[l], dot);
                        ptr_A   += A_ld;
                    };

                    ptr_B += B_ld;
                    ptr_C += C_ld;
                };
            }
            else if (t_B == trans_type::trans)
            {
                const V* ptr_As = A.ptr();

                for (Integer j = 0; j < N; ++j) 
                {			
                    ptr_A       = ptr_As;
                    Integer fc  = B.first_col(j);
                    Integer lc  = B.last_col(j);
                    Integer fe  = B.first_elem_pos_row(j);

                    for (Integer l = 0; l < M; ++l) 
                    {
                        V dot = Z;		
                        for (Integer k = fc, lp = fe; k <= lc; ++k) 
                        {
                            dot = dot + conj(ptr_A[k]) * ptr_B[lp];
                            lp  += B_ld - 1;
                        };

                        ptr_C[l] = alpha.eval(ptr_C[l], dot);
                        ptr_A   += A_ld;
                    };

                    ptr_C += C_ld;
                };
            }
            else
            {
                const V* ptr_As = A.ptr();

                for (Integer j = 0; j < N; ++j) 
                {			
                    ptr_A       = ptr_As;
                    Integer fc  = B.first_col(j);
                    Integer lc  = B.last_col(j);
                    Integer fe  = B.first_elem_pos_row(j);

                    for (Integer l = 0; l < M; ++l) 
                    {
                        V dot = Z;		
                        for (Integer k = fc, lp = fe; k <= lc; ++k) 
                        {
                            dot = dot + conj(ptr_A[k]) * conj(ptr_B[lp]);
                            lp  += B_ld - 1;
                        };

                        ptr_C[l] = alpha.eval(ptr_C[l], dot);
                        ptr_A   += A_ld;
                    };

                    ptr_C += C_ld;
                };
            }
        };
    };
};

template<class V, bool Is_compl = md::is_complex<V>::value>
struct eval_dense_band_lapack_sym
{
    using M1    = raw::Matrix<V,struct_dense>;
    using M2    = raw::Matrix<V,struct_banded>;

    static void eval_sym(M1& C, const M1& A, const M2& B, trans_type t_A, trans_type t_B,
                          Integer fr, Integer rows, Integer fc, Integer cols, const V& alpha,
                          const V& beta);
};
template<class V>
struct eval_dense_band_lapack_sym<V,true>
{
    using M1    = raw::Matrix<V,struct_dense>;
    using M2    = raw::Matrix<V,struct_banded>;

    static void eval_sym(M1& C, const M1& A, const M2& B, trans_type t_A, trans_type t_B,
                          Integer fr, Integer rows, Integer fc, Integer cols, const V& alpha,
                          const V& beta);
};

template<class V>
struct eval_dense_band_lapack
{
    using M1    = raw::Matrix<V,struct_dense>;
    using M2    = raw::Matrix<V,struct_banded>;

    static void eval_gemm(M1& C, const V& alpha, const M1& A, const M2& B, 
                     trans_type t_A, trans_type t_B, const V& beta, 
                     Integer fr, Integer rows, Integer fc, Integer cols)
    {
        if (rows == 0 || cols == 0)
            return;

        if (A.rows() == 0 || B.rows() == 0 || A.cols() == 0 || B.cols() == 0) 
        {
            prepare_gemm_C<V>::eval(beta, C, fr,rows,fc,cols);
            return;
        };

        if (raw::is_diag(A))
        {
            prepare_gemm_C<V>::eval(beta, C, fr,rows,fc,cols);
            return eval_dense_band_diag<V>::eval1_gemm(alpha,C,A,B, t_A, t_B, fr, rows, fc, cols);
        }

        if (raw::is_diag(B))
        {
            prepare_gemm_C<V>::eval(beta, C, fr,rows,fc,cols);
            return eval_dense_band_diag<V>::eval2_gemm(alpha,C,A,B, t_A, t_B, fr, rows, fc, cols);
        };
        
        eval_impl(C, A, B, t_A, t_B, fr, rows, fc, cols, alpha, beta);
    };

    static void eval(matcl::Matrix& ret, const M1& A, const M2& B, trans_type t_A, trans_type t_B)
    {
        Integer M, K, N;

        if (t_B == trans_type::no_trans)
            N = B.cols();
        else 
            N = B.rows();

        if (t_A == trans_type::no_trans)
        {
            M = A.rows();
            K = A.cols();
        }
        else
        {
            K = A.rows();
            M = A.cols();
        }

        if (K == 0 || N == 0 || M == 0) 
        {
            using sparse_matrix  = Matrix<V,struct_sparse>;
            sparse_matrix out(ti::ti_empty(), M, N);
            ret = matcl::Matrix(out,false);
            return;
        };

        if (raw::is_diag(A))
            return eval_dense_band_diag<V>::eval1(ret,A,B, t_A, t_B);

        if (raw::is_diag(B))
            return eval_dense_band_diag<V>::eval2(ret,A,B, t_A, t_B);

        M1 C(ti::ti_empty(),M, N);

        Integer fr      = 0;
        Integer rows    = M;
        Integer fc      = 0;
        Integer cols    = 0;
        V alpha         = V(1.0);
        V beta          = V(0.0);
        
        eval_impl(C, A, B, t_A, t_B, fr, rows, fc, cols, alpha, beta);

        bool is_full    = (fr == 0 && fc == 0 && rows == C.rows() && cols == C.cols());
        bool is_sq_A    = (A.rows() == A.cols());
        bool is_sq_B    = (B.rows() == B.cols());
        bool is_sq_C    = (C.rows() == C.cols());

        if (is_full == true)
            C.set_struct(predefined_struct_ext::mult_struct(A.get_struct(),B.get_struct(),t_A,t_B,
                                        is_real_matrix(A),is_real_matrix(B), is_sq_A, is_sq_B, is_sq_C));

        ret = matcl::Matrix(C,true);
    };

    static void eval_impl(M1& C, const M1& A, const M2& B, trans_type t_A, trans_type t_B,
                          Integer fr, Integer rows, Integer fc, Integer cols, const V& alpha,
                          const V& beta)
    {
        Integer M   = rows;
        Integer N   = cols;

        bool is_B_sq    = B.rows() == B.cols();

        if (M == 1 && N == 1)
            return eval_dot(C,A,B,t_A, t_B, fr, rows, fc, cols, alpha, beta);

        if (B.get_struct().is_symmetric(is_B_sq, is_real_matrix(B)) 
                && md::is_float_real_scalar<V>::value == true)
        {
            return eval_sym(C,A,B,t_A, t_B, fr, rows, fc, cols, alpha, beta);
        }
        else if (B.get_struct().is_hermitian(is_B_sq, is_real_matrix(B)) 
                 && md::is_float_real_scalar<V>::value == true)
        {
            return eval_sym(C,A,B,t_A, t_B, fr, rows, fc, cols, alpha, beta);
        }

        //for hermitian complex matrices conj is required; this is probably slower than general version

        return eval_general(C,A,B,t_A, t_B, fr, rows, fc, cols, alpha, beta);
    };

    static void eval_dot(M1& C, const M1& A, const M2& B, trans_type t_A, trans_type t_B,
                            Integer fr, Integer rows, Integer fc, Integer cols, const V& alpha,const V& beta)
    {
        (void)cols;
        (void)rows;

        Integer fd      = B.first_diag();
        Integer ld      = B.last_diag();
        Integer s       = ld - fd + 1;

        Integer lda     = (t_A == trans_type::no_trans)? A.ld() : 1;
        Integer ldb     = (t_B == trans_type::no_trans)? 1 : B.ld() - 1;

        const V* A_ptr  = A.ptr();
        const V* B_ptr  = B.rep_ptr() + B.first_elem_diag(0);

        V dot;

        if (t_A != trans_type::conj_trans && t_B != trans_type::conj_trans
            || md::is_float_real_scalar<V>::value == true)
        {
            V res(lapack::dot(s, md::lap(A_ptr), lda, md::lap(B_ptr), ldb));
            dot = res;
        }
        else if (t_A == trans_type::conj_trans && t_B == trans_type::conj_trans)
        {
            V res(lapack::dot(s, md::lap(A_ptr), lda, md::lap(B_ptr), ldb));
            dot = conj(res);
        }
        else if (t_A == trans_type::conj_trans)
        {
            V res   = V(lapack::dotc(s, md::lap(A_ptr), lda, md::lap(B_ptr), ldb));
            dot     = res;
        }
        else
        {
            //t_B = conj_trans
            V res   = V(lapack::dotc(s, md::lap(B_ptr), ldb, md::lap(A_ptr), lda));
            dot     = res;
        };

        V val       = alpha * dot;
        V* ptr_C    = C.ptr() + fr + fc * C.ld();

        if (mrd::is_zero(beta))
            ptr_C[0]    = val;
        else
            ptr_C[0]    = val + beta * ptr_C[0];

        return;
    };

    static void eval_general(M1& C, const M1& A, const M2& B, trans_type t_A, trans_type t_B,
                            Integer fr, Integer rows, Integer fc, Integer cols, const V& alpha,const V& beta)
    {
        if (B.has_diag(0) == false)
        {
            if (mrd::is_one(alpha) == true)
            {
                using Alpha     = alpha_one<V>;
                Alpha al;

                eval_dense_band_generic<V>::eval_impl(al, C, A, B, t_A, t_B, fr, rows, fc, cols);
            }
            else if (mrd::is_one(-alpha) == true)
            {
                using Alpha     = alpha_mone<V>;
                Alpha al;

                eval_dense_band_generic<V>::eval_impl(al, C, A, B, t_A, t_B, fr, rows, fc, cols);
            }
            else
            {
                using Alpha     = alpha_val<V>;
                Alpha al(alpha);

                eval_dense_band_generic<V>::eval_impl(al, C, A, B, t_A, t_B, fr, rows, fc, cols);
            };

            return;
        }

        bool V_is_real = md::is_real_scalar<V>::value;

        Integer M       = rows;

        const V *A_ptr  = A.ptr();
        V *C_ptr        = C.ptr() + fr + fc * C.ld();
        Integer A_ld    = A.ld();
        Integer kl      = B.number_subdiagonals();
        Integer ku      = B.number_superdiagonals();

        if (t_B == trans_type::no_trans || t_B == trans_type::trans || md::is_float_real_scalar<V>::value)
        {
            const char* TB = (t_B == trans_type::no_trans)? "T" : "N";
            if (t_A == trans_type::no_trans)
            {
                for (Integer i = 1; i <= M; ++i, ++A_ptr, ++C_ptr)
                {
                    lapack::gbmv(TB, B.rows(), B.cols(), kl, ku, *md::lap(&alpha), md::lap(B.rep_ptr()),
                                 B.ld(), md::lap(A_ptr), A.ld(), *md::lap(&beta), md::lap(C_ptr), C.ld());
                };
            }
            else if (t_A == trans_type::trans || t_A == trans_type::conj_trans && V_is_real == true)
            {
                for (Integer i = 1; i <= M; ++i, ++C_ptr)
                {
                    lapack::gbmv(TB, B.rows(), B.cols(), kl, ku, *md::lap(&alpha), md::lap(B.rep_ptr()),
                                 B.ld(), md::lap(A_ptr), 1, *md::lap(&beta), md::lap(C_ptr), C.ld());
                    A_ptr += A_ld;
                };
            }
            else
            {
                Integer K           = A.rows();

                using VTR_pod       = matcl::pod_type<V>;
                using workspace     = matcl::pod_workspace<VTR_pod>;
                workspace WORK      = workspace(K);
                V* ptr_w            = reinterpret_cast<V*>(WORK.ptr());

                for (Integer i = 1; i <= M; ++i, ++C_ptr)
                {
                    for (Integer j = 0; j < K; ++j)
                    {
                        ptr_w[j]    = conj(A_ptr[j]);
                    };	

                    lapack::gbmv(TB, B.rows(), B.cols(), kl, ku, *md::lap(&alpha), md::lap(B.rep_ptr()), 
                                 B.ld(), md::lap(ptr_w), 1, *md::lap(&beta), md::lap(C_ptr), C.ld());

                    A_ptr += A_ld;
                };
            };
        }
        else
        {
            Integer K           = B.cols();
            bool beta_nonzero   = (mrd::is_zero(beta) == false);
            V beta_conj         = conj(beta);
            V alpha_conj        = conj(alpha);

            //take conj transpose of C = (a * A * B + b * C)
            if (t_A == trans_type::conj_trans || V_is_real == true)
            {
                for (Integer i = 1; i <= M; ++i, ++C_ptr)
                {
                    if (beta_nonzero == true)
                        lapack::lacgv(B.rows(), md::lap(C_ptr), C.ld());

                    lapack::gbmv("N", B.rows(), B.cols(), kl, ku, *md::lap(&alpha_conj), md::lap(B.rep_ptr()), 
                                 B.ld(), md::lap(A_ptr), 1, *md::lap(&beta_conj), md::lap(C_ptr), C.ld());

                    lapack::lacgv(B.rows(), md::lap(C_ptr), C.ld());
                    A_ptr += A_ld;
                };
            }
            else if (t_A == trans_type::no_trans)
            {
                using VTR_pod       = matcl::pod_type<V>;
                using workspace     = matcl::pod_workspace<VTR_pod>;
                workspace WORK      = workspace(K);
                V* ptr_w            = reinterpret_cast<V*>(WORK.ptr());

                for (Integer i = 1; i <= M; ++i, ++A_ptr, ++C_ptr)
                {
                    for (Integer j = 0; j < K; ++j)
                        ptr_w[j]    = conj(A_ptr[j*A_ld]);

                    if (beta_nonzero == true)
                        lapack::lacgv(B.rows(), md::lap(C_ptr), C.ld());

                    lapack::gbmv("N", B.rows(), B.cols(), kl, ku, *md::lap(&alpha_conj), md::lap(B.rep_ptr()),
                                 B.ld(), md::lap(ptr_w), 1, *md::lap(&beta_conj), md::lap(C_ptr), C.ld());

                    lapack::lacgv(B.rows(), md::lap(C_ptr), C.ld());
                };
            }
            else
            {
                using VTR_pod       = matcl::pod_type<V>;
                using workspace     = matcl::pod_workspace<VTR_pod>;
                workspace WORK      = workspace(K);
                V* ptr_w            = reinterpret_cast<V*>(WORK.ptr());

                for (Integer i = 1; i <= M; ++i, ++C_ptr)
                {
                    for (Integer j = 0; j < K; ++j)
                        ptr_w[j]    = conj(A_ptr[j]);

                    if (beta_nonzero == true)
                        lapack::lacgv(B.rows(), md::lap(C_ptr), C.ld());

                    lapack::gbmv("N", B.rows(), B.cols(), kl, ku, *md::lap(&alpha_conj), md::lap(B.rep_ptr()),
                                 B.ld(), md::lap(ptr_w), 1, *md::lap(&beta_conj), md::lap(C_ptr), C.ld());

                    lapack::lacgv(B.rows(), md::lap(C_ptr), C.ld());
                    A_ptr += A_ld;
                };
            }
        };
        return;
    };

    static void eval_sym(M1& C, const M1& A, const M2& B, trans_type t_A, trans_type t_B,
                        Integer fr, Integer rows, Integer fc, Integer cols, const V& alpha,const V& beta)
    {
        return eval_dense_band_lapack_sym<V>::eval_sym(C,A,B,t_A,t_B,fr,rows,fc,cols,alpha,beta);
    };
};

template<class V, bool Is_compl>
void eval_dense_band_lapack_sym<V,Is_compl>
    ::eval_sym(M1& C, const M1& A, const M2& B, trans_type t_A, trans_type t_B,
                Integer fr, Integer rows, Integer fc, Integer cols, const V& alpha, const V& beta)
{
    (void)cols;
    (void)rows;
    (void)t_B;

    bool V_is_real = md::is_real_scalar<V>::value;

    Integer kl      = -B.first_diag();
    Integer B_ld    = B.ld();

    const V* A_ptr  = A.ptr();    
    Integer A_ld    = A.ld();
    V *C_ptr        = C.ptr() + fr + fc * C.ld();

    if (t_A == trans_type::no_trans)
    {
        for (Integer i = 1; i <= A.rows(); ++i, ++C_ptr, ++A_ptr)
        {
            lapack::sbmv<V>("U",B.rows(), kl, *md::lap(&alpha), md::lap(B.rep_ptr()), B_ld, 
                            md::lap(A_ptr), A.ld(),*md::lap(&beta), md::lap(C_ptr), C.ld());
        }
    }
    else if (t_A == trans_type::trans || t_A == trans_type::conj_trans && V_is_real == true)
    {
        for (Integer i = 1; i <= A.cols(); ++i, ++C_ptr)
        {
            lapack::sbmv<V>("U",B.rows(), kl, *md::lap(&alpha), md::lap(B.rep_ptr()), B_ld, 
                            md::lap(A_ptr), 1, *md::lap(&beta), md::lap(C_ptr), C.ld());
            A_ptr += A_ld;
        }
    }
    else
    {
        Integer K           = A.rows();

        using VTR_pod       = matcl::pod_type<V>;
        using workspace     = matcl::pod_workspace<VTR_pod>;
        workspace WORK      = workspace(K);
        V* ptr_w            = reinterpret_cast<V*>(WORK.ptr());

        for (Integer i = 1; i <= A.cols(); ++i, ++C_ptr)
        {
            for (Integer j = 0; j < A.rows(); ++j)
            {
                ptr_w[j]    = conj(A_ptr[j]);
            };		

            lapack::sbmv<V>("U",B.rows(), kl, *md::lap(&alpha), md::lap(B.rep_ptr()), B_ld, 
                            md::lap(ptr_w), 1, *md::lap(&beta), md::lap(C_ptr), C.ld());

            A_ptr += A_ld;
        };
    };
};

template<class V>
void eval_dense_band_lapack_sym<V, true>
    ::eval_sym(M1& C, const M1& A, const M2& B, trans_type t_A, trans_type t_B,
                Integer fr, Integer rows, Integer fc, Integer cols, const V& alpha, const V& beta)
{
    return eval_dense_band_lapack<V>::eval_general(C,A,B, t_A, t_B,fr,rows,fc,cols,alpha,beta);
};

template<> struct eval_dense_band<Integer>
{
    using VT    = Integer;
    using Mat_D = raw::Matrix<VT,struct_dense>;
    using Mat_B = raw::Matrix<VT,struct_banded>;

    static void eval(matcl::Matrix& ret, const Mat_D& A, const Mat_B& B,
                     trans_type t_A, trans_type t_B)
    {
        return eval_dense_band_generic<VT>::eval(ret,A,B,t_A, t_B);
    };

    static void eval_gemm(Mat_D& C, const VT& alpha, const Mat_D& A, const Mat_B& B, 
                     trans_type t_A, trans_type t_B, const VT& beta, 
                     Integer fr, Integer rows, Integer fc, Integer cols)
    {
        return eval_dense_band_generic<VT>::eval_gemm(C,alpha,A,B,t_A,t_B,beta,fr,rows,fc,cols);
    };
}; 
template<> struct eval_dense_band<Object>
{
    using VT    = Object;
    using Mat_D = raw::Matrix<VT,struct_dense>;
    using Mat_B = raw::Matrix<VT,struct_banded>;

    static void eval(matcl::Matrix& ret, const Mat_D& A, const Mat_B& B,
                     trans_type t_A, trans_type t_B)
    {
        return eval_dense_band_generic<VT>::eval(ret,A,B,t_A, t_B);
    };

    static void eval_gemm(Mat_D& C, const VT& alpha, const Mat_D& A, const Mat_B& B, 
                     trans_type t_A, trans_type t_B, const VT& beta, 
                     Integer fr, Integer rows, Integer fc, Integer cols)
    {
        return eval_dense_band_generic<VT>::eval_gemm(C,alpha,A,B,t_A,t_B,beta,fr,rows,fc,cols);
    };
}; 
template<> struct eval_dense_band<Real>
{
    using VT    = Real;
    using Mat_D = raw::Matrix<VT,struct_dense>;
    using Mat_B = raw::Matrix<VT,struct_banded>;

    static void eval(matcl::Matrix& ret, const Mat_D& A, const Mat_B& B, 
                     trans_type t_A, trans_type t_B)
    {
        return eval_dense_band_lapack<Real>::eval(ret,A,B,t_A, t_B);
    };

    static void eval_gemm(Mat_D& C, const VT& alpha, const Mat_D& A, const Mat_B& B, 
                     trans_type t_A, trans_type t_B, const VT& beta, 
                     Integer fr, Integer rows, Integer fc, Integer cols)
    {
        return eval_dense_band_lapack<VT>::eval_gemm(C,alpha,A,B,t_A,t_B,beta,fr,rows,fc,cols);
    };
};
template<> struct eval_dense_band<Float>
{
    using VT    = Float;
    using Mat_D = raw::Matrix<VT,struct_dense>;
    using Mat_B = raw::Matrix<VT,struct_banded>;

    static void eval(matcl::Matrix& ret, const Mat_D& A, const Mat_B& B, 
                     trans_type t_A, trans_type t_B)
    {
        return eval_dense_band_lapack<VT>::eval(ret,A,B,t_A, t_B);
    };

    static void eval_gemm(Mat_D& C, const VT& alpha, const Mat_D& A, const Mat_B& B, 
                     trans_type t_A, trans_type t_B, const VT& beta, 
                     Integer fr, Integer rows, Integer fc, Integer cols)
    {
        return eval_dense_band_lapack<VT>::eval_gemm(C,alpha,A,B,t_A,t_B,beta,fr,rows,fc,cols);
    };
};

template<> struct eval_dense_band<Complex>
{
    using VT    = Complex;
    using Mat_D = raw::Matrix<VT,struct_dense>;
    using Mat_B = raw::Matrix<VT,struct_banded>;

    static void eval(matcl::Matrix& ret, const Mat_D& A, const Mat_B& B, 
                     trans_type t_A, trans_type t_B)
    {
        return eval_dense_band_lapack<VT>::eval(ret,A,B,t_A, t_B);
    };

    static void eval_gemm(Mat_D& C, const VT& alpha, const Mat_D& A, const Mat_B& B, 
                     trans_type t_A, trans_type t_B, const VT& beta, 
                     Integer fr, Integer rows, Integer fc, Integer cols)
    {
        return eval_dense_band_lapack<VT>::eval_gemm(C,alpha,A,B,t_A,t_B,beta,fr,rows,fc,cols);
    };
};
template<> struct eval_dense_band<Float_complex>
{
    using VT    = Float_complex;
    using Mat_D = raw::Matrix<VT,struct_dense>;
    using Mat_B = raw::Matrix<VT,struct_banded>;

    static void eval(matcl::Matrix& ret, const Mat_D& A, const Mat_B& B, 
                     trans_type t_A, trans_type t_B)
    {
        return eval_dense_band_lapack<VT>::eval(ret,A,B,t_A, t_B);
    };

    static void eval_gemm(Mat_D& C, const VT& alpha, const Mat_D& A, const Mat_B& B, 
                     trans_type t_A, trans_type t_B, const VT& beta, 
                     Integer fr, Integer rows, Integer fc, Integer cols)
    {
        return eval_dense_band_lapack<VT>::eval_gemm(C,alpha,A,B,t_A,t_B,beta,fr,rows,fc,cols);
    };
};

template<class Val_ret, class M1, class M2>
struct eval_mult<Val_ret,M1,M2,struct_dense,struct_banded> 
{ 
    using Mat_D = raw::Matrix<Val_ret, struct_dense>;

    static void eval(matcl::Matrix& ret, const M1& A, const M2& B, trans_type t_A, trans_type t_B)
    {
        using MP1           = Matrix<Val_ret,struct_dense>;
        using MP2           = Matrix<Val_ret,struct_banded>;
        return eval_dense_band<Val_ret>
            ::eval(ret,converter<MP1,M1>::eval(A),converter<MP2,M2>::eval(B), t_A, t_B);
    };

    static void eval_gemm(Mat_D& C, const Val_ret& alpha, const M1& A, const M2& B, 
                     trans_type t_A, trans_type t_B, const Val_ret& beta, 
                     Integer fr, Integer rows, Integer fc, Integer cols)
    {
        using MP1           = Matrix<Val_ret,struct_dense>;
        using MP2           = Matrix<Val_ret,struct_banded>;
        return eval_dense_band<Val_ret>
            ::eval_gemm(C,alpha,converter<MP1,M1>::eval(A),converter<MP2,M2>::eval(B), t_A, t_B,
                        beta, fr, rows, fc, cols);
    };
};

}}}