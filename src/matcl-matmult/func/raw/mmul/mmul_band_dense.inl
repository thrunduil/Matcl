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
#include "matcl-matrep/lib_functions/func_matrix.h"

namespace matcl { namespace raw { namespace details
{

//========================================================================
//                      BAND - DENSE
//========================================================================
template<class V>
struct eval_band_dense_diag
{
    using M1    = raw::Matrix<V,struct_dense>;
    using M2    = raw::Matrix<V,struct_banded>;

    static void eval_diag_diag(matcl::Matrix& ret, const M2& A, const M1& B, trans_type t_A, trans_type t_B)
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

        const V* A_ptr  = A.rep_ptr() + A.first_elem_diag(0);
        const V* B_ptr  = B.ptr();

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
                    A_ptr       += A_ld;
                    B_ptr       += B_ld+1;
                };
            }
            else
            {
                for (Integer i = 0; i < K; ++i)
                {
                    C_ptr[i]    = conj(A_ptr[0]) * B_ptr[0];
                    A_ptr       += A_ld;
                    B_ptr       += B_ld+1;
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
                    A_ptr       += A_ld;
                    B_ptr       += B_ld+1;
                };
            }
            else
            {
                for (Integer i = 0; i < K; ++i)
                {
                    C_ptr[i]    = A_ptr[0] * B_ptr[0];
                    A_ptr       += A_ld;
                    B_ptr       += B_ld+1;
                };
            };
        };

        ret = matcl::bdiags(matcl::Matrix(C,false), Integer(0), M, N);
        return;
    };

    static void eval_diag_diag_gemm(const V& alpha, M1& C, const M2& A, const M1& B, trans_type t_A, trans_type t_B, 
                                    Integer fr, Integer rows, Integer fc, Integer cols)
    {
        Integer M   = rows;
        Integer N   = cols;
        Integer K   = std::min(A.rows(), A.cols());
        Integer MN  = std::min(M,N);
        K           = std::min(MN,K);

        const V* A_ptr  = A.rep_ptr() + A.first_elem_diag(0);
        const V* B_ptr  = B.ptr();

        Integer A_ld    = A.ld();
        Integer B_ld    = B.ld();
        Integer C_ld    = C.ld();

        V* C_ptr        = C.ptr() + fr + fc * C_ld;

        if (t_A == trans_type::conj_trans)
        {
            if (t_B == trans_type::conj_trans)
            {
                for (Integer i = 0; i < K; ++i)
                {
                    V tmp       = conj(A_ptr[0]) * conj(B_ptr[0]);
                    C_ptr[0]    = C_ptr[0] + alpha * tmp;
                    A_ptr       += A_ld;
                    B_ptr       += B_ld+1;
                    C_ptr       += C_ld+1;
                };
            }
            else
            {
                for (Integer i = 0; i < K; ++i)
                {
                    V tmp       = conj(A_ptr[0]) * B_ptr[0];
                    C_ptr[0]    = C_ptr[0] + alpha * tmp;
                    A_ptr       += A_ld;
                    B_ptr       += B_ld+1;
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
                    A_ptr       += A_ld;
                    B_ptr       += B_ld+1;
                    C_ptr       += C_ld+1;
                };
            }
            else
            {
                for (Integer i = 0; i < K; ++i)
                {
                    V tmp       = A_ptr[0] * B_ptr[0];
                    C_ptr[0]    = C_ptr[0] + alpha * tmp;
                    A_ptr       += A_ld;
                    B_ptr       += B_ld+1;
                    C_ptr       += C_ld+1;
                };
            };
        };
        return;
    };

    static void eval1_impl(const V& alpha, M1& C, const M2& A, const M1& B, trans_type t_A, 
                           trans_type t_B, Integer fr, Integer rows, Integer fc, Integer cols)
    {
        (void)rows;

        Integer K1 = std::min(A.rows(),A.cols());

        M1 diag(A.get_type(), K1, 1);

        Integer A_ld    = A.ld();
        Integer B_ld    = B.ld();
        Integer C_ld    = C.ld();

        Integer N       = cols;

        const V* A_ptr  = A.rep_ptr() + A.first_elem_diag(0);
        const V* B_ptr  = B.ptr();
        V* D_ptr        = diag.ptr();
        V* C_ptr        = C.ptr() + fr + fc * C_ld;

        if (t_A == trans_type::conj_trans)
        {
            for(Integer i = 0; i < K1; ++i)
            {
                D_ptr[i]    = alpha*conj(A_ptr[0]);
                A_ptr       += A_ld;
            };
        }
        else
        {
            for(Integer i = 0; i < K1; ++i)
            {
                D_ptr[i]    = alpha*A_ptr[0];
                A_ptr       += A_ld;
            };
        };

        if (t_B == trans_type::no_trans)
        {
            for (Integer j = 0; j < N; ++j)
            {
                for (Integer i = 0; i < K1; ++i)
                {
                    C_ptr[i]    = C_ptr[i] + D_ptr[i] * B_ptr[i];
                };

                C_ptr       += C_ld;
                B_ptr       += B_ld;
            };
        }
        else if (t_B == trans_type::trans)
        {
            for (Integer j = 0; j < N; ++j)
            {
                for (Integer i = 0; i < K1; ++i)
                {
                    C_ptr[i]    = C_ptr[i] + D_ptr[i] * B_ptr[i*B_ld];
                };

                C_ptr       += C_ld;
                B_ptr       += 1;
            };
        }
        else
        {
            for (Integer j = 0; j < N; ++j)
            {
                for (Integer i = 0; i < K1; ++i)
                {
                    C_ptr[i]    = C_ptr[i] + D_ptr[i] * conj(B_ptr[i*B_ld]);
                };

                C_ptr       += C_ld;
                B_ptr       += 1;
            };
        };
    };

    static void eval1(matcl::Matrix& ret, const M2& A, const M1& B, trans_type t_A, trans_type t_B)
    {
        if (raw::is_diag(B))
            return eval_diag_diag(ret,A,B,t_A, t_B);

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

        ti::ti_type<V> ret_ti = ti::get_return_ti<ti::ti_type<V>>(mdyf::op_mul::eval(),
                                                                  ti::get_ti(A), ti::get_ti(B));
        V Z = md::default_value<V>(ret_ti);

        M1 C(ret_ti, Z, M, N);

        Integer fr      = 0;
        Integer rows    = M;
        Integer fc      = 0;
        Integer cols    = N;
        V alpha         = V(1.0);

        eval1_impl(alpha, C, A, B, t_A, t_B, fr, rows, fc, cols);

        bool is_sq_A    = (A.rows() == A.cols());
        bool is_sq_B    = (B.rows() == B.cols());
        bool is_sq_C    = (C.rows() == C.cols());
        C.set_struct(predefined_struct_ext::mult_struct(A.get_struct(),B.get_struct(),t_A,t_B,
                                        is_real_matrix(A),is_real_matrix(B), is_sq_A, is_sq_B, is_sq_C));

        ret = matcl::Matrix(C,true);
        return;
    };

    static void eval1_gemm(const V& alpha, M1& C, const M2& A, const M1& B, trans_type t_A, trans_type t_B, 
                     Integer fr, Integer rows, Integer fc, Integer cols)
    {
        if (raw::is_diag(B))
            return eval_diag_diag_gemm(alpha,C,A,B,t_A, t_B,fr,rows,fc,cols);

        eval1_impl(alpha, C, A, B, t_A, t_B, fr, rows, fc, cols);
        return;
    }

    static void eval2(matcl::Matrix& ret, const M2& A, const M1& B, trans_type t_A, trans_type t_B)
    {
        if (raw::is_diag(A))
            return eval_diag_diag(ret,A,B, t_A, t_B);

        Integer N, M, K;
        Integer fdc, ldc;

        if (t_B == trans_type::no_trans)
            N       = B.cols();
        else
            N       = B.rows();

        if (t_A == trans_type::no_trans)
        {
            M       = A.rows();
            K       = A.cols();
            fdc     = A.first_diag();
            ldc     = A.last_diag();
        }
        else
        {
            K       = A.rows();
            M       = A.cols();
            fdc     = -A.last_diag();
            ldc     = -A.first_diag();
        };        

        Integer K1 = std::min(N,K);

        M1 diag(B.get_type(), K1, 1);

        ti::ti_type<V> ret_ti = ti::get_return_ti<ti::ti_type<V>>(mdyf::op_mul::eval(),
                                                                  ti::get_ti(A), ti::get_ti(B));
        V Z = md::default_value<V>(ret_ti);

        M2 C(ret_ti, Z, M, N, fdc, ldc);

        const V* A_ptr  = A.rep_ptr();
        const V* B_ptr  = B.ptr();
        V* C_ptr        = C.rep_ptr();
        V* D_ptr        = diag.ptr();

        Integer A_ld    = A.ld();
        Integer B_ld    = B.ld();
        Integer C_ld    = C.ld();

        if (t_B == trans_type::conj_trans)
        {
            for(Integer i = 0; i < K1; ++i)
            {
                D_ptr[i]    = conj(B_ptr[0]);
                B_ptr       += B_ld + 1;
            };
        }
        else
        {
            for(Integer i = 0; i < K1; ++i)
            {
                D_ptr[i]    = B_ptr[0];
                B_ptr       += B_ld + 1;
            };
        };

        if (t_A == trans_type::no_trans)
        {
            for (Integer j = 0; j < K1; ++j)
            {
                const V& tmp    = D_ptr[j];
                if (mrd::is_zero(tmp))
                {
                    C_ptr   += C_ld;
                    A_ptr   += A_ld;
                    continue;
                };

                Integer first_row   = A.first_row(j);
                Integer last_row    = A.last_row(j);
                Integer pos_A       = A.first_elem_pos(j);
                Integer pos_C       = C.first_elem_pos(j);

                for (Integer l = first_row; l <= last_row; ++l) 
                {
                    C_ptr[pos_C]    = tmp * A_ptr[pos_A];

                    ++pos_A;
                    ++pos_C;
                };

                C_ptr       += C_ld;
                A_ptr       += A_ld;
            };
        }
        else if (t_A == trans_type::trans)
        {
            for (Integer d = A.first_diag(); d <= A.last_diag(); ++d)
            {
                Integer pos_A = A.first_elem_diag(d);
                Integer pos_C = C.first_elem_diag(-d);
                Integer fc    = C.first_col_on_diag(-d);

                Integer s     = std::min(A.diag_length(d), K1 - fc);

                for (Integer i = 0; i < s; ++i, ++fc)
                {
                    C_ptr[pos_C] = A_ptr[pos_A] * D_ptr[fc];

                    pos_A   += A_ld;
                    pos_C   += C_ld;
                };
            };
        }
        else
        {
            for (Integer d = A.first_diag(); d <= A.last_diag(); ++d)
            {
                Integer pos_A = A.first_elem_diag(d);
                Integer pos_C = C.first_elem_diag(-d);
                Integer fc    = C.first_col_on_diag(-d);

                Integer s     = std::min(A.diag_length(d), K1 - fc);

                for (Integer i = 0; i < s; ++i, ++fc)
                {
                    C_ptr[pos_C] = conj(A_ptr[pos_A]) * D_ptr[fc];

                    pos_A   += A_ld;
                    pos_C   += C_ld;
                };
            };
        };

        bool is_sq_A    = (A.rows() == A.cols());
        bool is_sq_B    = (B.rows() == B.cols());
        bool is_sq_C    = (C.rows() == C.cols());
        C.set_struct(predefined_struct_ext::mult_struct(A.get_struct(),B.get_struct(),t_A,t_B,
                        is_real_matrix(A),is_real_matrix(B), is_sq_A, is_sq_B, is_sq_C));

        ret = matcl::Matrix(C,true);
        return;
    };    
    
    static void eval2_gemm(const V& alpha, M1& C, const M2& A, const M1& B, trans_type t_A, trans_type t_B, 
                     Integer fr, Integer rows, Integer fc, Integer cols)
    {
        if (raw::is_diag(A))
            return eval_diag_diag_gemm(alpha,C,A,B,t_A, t_B,fr,rows,fc,cols);

        Integer K1 = std::min(B.rows(),B.cols());

        M1 diag(B.get_type(), K1, 1);

        Integer A_ld    = A.ld();
        Integer B_ld    = B.ld();
        Integer C_ld    = C.ld();

        const V* A_ptr  = A.rep_ptr();
        const V* B_ptr  = B.ptr();
        V* C_ptr        = C.ptr() + fr + fc * C_ld;
        V* D_ptr        = diag.ptr();

        if (t_B == trans_type::conj_trans)
        {
            for(Integer i = 0; i < K1; ++i)
            {
                D_ptr[i]    = alpha * conj(B_ptr[0]);
                B_ptr       += B_ld + 1;
            };
        }
        else
        {
            for(Integer i = 0; i < K1; ++i)
            {
                D_ptr[i]    = alpha * B_ptr[0];
                B_ptr       += B_ld + 1;
            };
        };

        if (t_A == trans_type::no_trans)
        {
            for (Integer j = 0; j < K1; ++j)
            {
                const V& tmp    = D_ptr[j];
                if (mrd::is_zero(tmp))
                {
                    C_ptr   += C_ld;
                    A_ptr   += A_ld;
                    continue;
                };

                Integer first_row   = A.first_row(j);
                Integer last_row    = A.last_row(j);
                Integer pos_A       = A.first_elem_pos(j);

                for (Integer l = first_row; l <= last_row; ++l) 
                {
                    C_ptr[l]        = C_ptr[l] + tmp * A_ptr[pos_A];
                    ++pos_A;
                };

                C_ptr       += C_ld;
                A_ptr       += A_ld;
            };
        }
        else if (t_A == trans_type::trans)
        {
            for (Integer j = 0; j < K1; ++j)
            {
                const V& tmp    = D_ptr[j];
                if (mrd::is_zero(tmp))
                {
                    C_ptr   += C_ld;
                    continue;
                };

                Integer first_col   = A.first_col(j);
                Integer last_col    = A.last_col(j);
                Integer pos_A       = A.first_elem_pos_row(j);

                for (Integer l = first_col; l <= last_col; ++l) 
                {
                    C_ptr[l]        = C_ptr[l] + tmp * A_ptr[pos_A];
                    pos_A           += A_ld - 1;
                };

                C_ptr       += C_ld;
            };
        }
        else
        {
            for (Integer j = 0; j < K1; ++j)
            {
                const V& tmp    = D_ptr[j];
                if (mrd::is_zero(tmp))
                {
                    C_ptr   += C_ld;
                    continue;
                };

                Integer first_col   = A.first_col(j);
                Integer last_col    = A.last_col(j);
                Integer pos_A       = A.first_elem_pos_row(j);

                for (Integer l = first_col; l <= last_col; ++l) 
                {
                    C_ptr[l]        = C_ptr[l] + tmp * conj(A_ptr[pos_A]);
                    pos_A           += A_ld - 1;
                };

                C_ptr       += C_ld;
            };
        };

        return;
    };
};

template<class V> 
struct eval_band_dense_generic
{
    using M1    = raw::Matrix<V,struct_dense>;
    using M2    = raw::Matrix<V,struct_banded>;

    template<class Alpha>
    static void eval_impl(const Alpha& alpha, M1& C, const M2& A, const M1& B, trans_type t_A, trans_type t_B,
                          Integer fr, Integer rows, Integer fc, Integer cols)
    {
        V Z = md::default_value<V>(C.get_type());

        const V* ptr_B  = B.ptr();        
        
        Integer A_ld    = A.ld();
        Integer B_ld    = B.ld();
        Integer C_ld    = C.ld();

        Integer M       = rows;
        Integer N       = cols;
        V* ptr_C        = C.ptr() + fr + fc*C_ld;

        if (t_A == trans_type::no_trans)
        {
            Integer K   = A.cols();

            if (t_B == trans_type::no_trans)
            {
                for (Integer j = 0; j < N; ++j) 
                {
                    const V* ptr_A = A.rep_ptr();

                    for (Integer l = 0; l < K; ++l) 
                    {
                        const V& temp = ptr_B[l];

                        if (mrd::is_zero(temp) == false) 
                        {					
                            Integer first_row   = A.first_row(l);
                            Integer last_row    = A.last_row(l);
                            Integer pos_A       = A.first_elem_pos(l);

                            for (Integer k = first_row; k <= last_row; ++k) 
                            {
                                V tmp           = ptr_A[pos_A++] *  temp;
                                ptr_C[k]        = alpha.eval(ptr_C[k],tmp);
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
                bool is_conj = (t_B == trans_type::conj_trans);

                for (Integer j = 0; j < N; ++j) 
                {
                    const V* ptr_A = A.rep_ptr();

                    for (Integer l = 0; l < K; ++l) 
                    {
                        const V& temp = is_conj ? conj(ptr_B[l*B_ld]) : ptr_B[l*B_ld];

                        if (mrd::is_zero(temp) == false) 
                        {					
                            Integer first_row   = A.first_row(l);
                            Integer last_row    = A.last_row(l);
                            Integer pos_A       = A.first_elem_pos(l);

                            for (Integer k = first_row; k <= last_row; ++k) 
                            {
                                V tmp           = ptr_A[pos_A++] * temp;
                                ptr_C[k]        = alpha.eval(ptr_C[k],tmp);
                            };
                        };

                        ptr_A += A_ld;
                    };

                    ptr_B += 1;
                    ptr_C += C_ld;
                };
            };
        }
        else if (t_A == trans_type::trans)
        {
            if (t_B == trans_type::no_trans)
            {
                for (Integer j = 0; j < N; ++j) 
                {
                    const V* ptr_A = A.rep_ptr();

                    for (Integer l = 0; l < M; ++l) 
                    {
                        V dot = Z;
        
                        Integer first_row   = A.first_row(l);
                        Integer last_row    = A.last_row(l);
                        Integer pos_A       = A.first_elem_pos(l);

                        for (Integer k = first_row; k <= last_row; ++k) 
                        {
                            V tmp   = ptr_A[pos_A++] * ptr_B[k];
                            dot     = dot + tmp;
                        };

                        ptr_C[l]    = alpha.eval(ptr_C[l],dot);
                        ptr_A       += A_ld;
                    };

                    ptr_B += B_ld;
                    ptr_C += C_ld;
                };
            }
            else
            {
                bool is_conj = (t_B == trans_type::conj_trans);

                for (Integer j = 0; j < N; ++j) 
                {
                    const V* ptr_A = A.rep_ptr();

                    for (Integer l = 0; l < M; ++l) 
                    {
                        V dot = Z;
        
                        Integer first_row   = A.first_row(l);
                        Integer last_row    = A.last_row(l);
                        Integer pos_A       = A.first_elem_pos(l);

                        if (is_conj == false)
                        {
                            for (Integer k = first_row; k <= last_row; ++k) 
                            {
                                V tmp   = ptr_A[pos_A++] * ptr_B[k*B_ld];
                                dot     = dot + tmp;
                            };
                        }
                        else
                        {
                            for (Integer k = first_row; k <= last_row; ++k) 
                            {
                                V tmp   = ptr_A[pos_A++] * conj(ptr_B[k*B_ld]);
                                dot     = dot + tmp;
                            };
                        };

                        ptr_C[l]    = alpha.eval(ptr_C[l], dot);
                        ptr_A       += A_ld;
                    };

                    ptr_B += 1;
                    ptr_C += C_ld;
                };
            };
        }
        else
        {
            if (t_B == trans_type::no_trans)
            {
                for (Integer j = 0; j < N; ++j) 
                {
                    const V* ptr_A = A.rep_ptr();

                    for (Integer l = 0; l < M; ++l) 
                    {
                        V dot = Z;
        
                        Integer first_row   = A.first_row(l);
                        Integer last_row    = A.last_row(l);
                        Integer pos_A       = A.first_elem_pos(l);

                        for (Integer k = first_row; k <= last_row; ++k) 
                        {
                            V temp          = conj(ptr_A[pos_A++]) * ptr_B[k];
                            dot             = dot + temp;
                        };

                        ptr_C[l]    = alpha.eval(ptr_C[l], dot);
                        ptr_A       += A_ld;
                    };
                    ptr_B   += B_ld;
                    ptr_C   += C_ld;
                };
            }
            else
            {
                bool is_conj = (t_B == trans_type::conj_trans);

                for (Integer j = 0; j < N; ++j) 
                {
                    const V* ptr_A = A.rep_ptr();

                    for (Integer l = 0; l < M; ++l) 
                    {
                        V dot = Z;
        
                        Integer first_row   = A.first_row(l);
                        Integer last_row    = A.last_row(l);
                        Integer pos_A       = A.first_elem_pos(l);

                        if (is_conj == false)
                        {
                            for (Integer k = first_row; k <= last_row; ++k) 
                            {
                                V temp          = conj(ptr_A[pos_A++] * ptr_B[k*B_ld]);
                                dot             = dot + temp;
                            };
                        }
                        else
                        {
                            for (Integer k = first_row; k <= last_row; ++k) 
                            {
                                V temp          = conj(ptr_A[pos_A++]) * conj(ptr_B[k*B_ld]);
                                dot             = dot + temp;
                            };
                        };

                        ptr_C[l]    = alpha.eval(ptr_C[l], dot);
                        ptr_A       += A_ld;
                    };
                    ptr_B   += 1;
                    ptr_C   += C_ld;
                };
            };
        }
    };

    static void eval(matcl::Matrix& ret, const M2& A, const M1& B, trans_type t_A, trans_type t_B)
    {
        if (raw::is_diag(A))
            return eval_band_dense_diag<V>::eval1(ret,A,B, t_A, t_B);

        if (raw::is_diag(B))
            return eval_band_dense_diag<V>::eval2(ret,A,B, t_A, t_B);

        ti::ti_type<V> ret_ti = ti::get_return_ti<ti::ti_type<V>>(mdyf::op_mul::eval(),
                                                                  ti::get_ti(A), ti::get_ti(B));

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

        V Z = md::default_value<V>(ret_ti);

        if (M == 0 || K == 0 || N == 0)
        {
            using sparse_mat = Matrix<V,struct_sparse>;
            sparse_mat out(ret_ti, M, N);
            ret = matcl::Matrix(out,false);
            return;
        };

        M1 C(ret_ti, Z, M, N);

        Integer fr      = 0;
        Integer rows    = M;
        Integer fc      = 0;
        Integer cols    = N;

        using Alpha     = alpha_one<V>;
        Alpha alpha;

        eval_impl(alpha, C, A, B, t_A, t_B, fr, rows, fc, cols);

        bool is_sq_A    = (A.rows() == A.cols());
        bool is_sq_B    = (B.rows() == B.cols());
        bool is_sq_C    = (C.rows() == C.cols());
        C.set_struct(predefined_struct_ext::mult_struct(A.get_struct(),B.get_struct(),t_A,t_B,
                            is_real_matrix(A), is_real_matrix(B), is_sq_A, is_sq_B, is_sq_C));

        ret = matcl::Matrix(C,true);
        return;
    };

    static void eval_gemm(M1& C, const V& alpha, const M2& A, const M1& B, 
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
            return eval_band_dense_diag<V>::eval1_gemm(alpha, C, A, B, t_A, t_B, fr, rows, fc, cols);

        if (is_diag_B)
            return eval_band_dense_diag<V>::eval2_gemm(alpha, C, A, B, t_A, t_B, fr, rows, fc, cols);

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
};

template<class V> 
struct eval_band_dense_lapack
{
    using M1    = raw::Matrix<V,struct_dense>;
    using M2    = raw::Matrix<V,struct_banded>;
    static const bool ISR = md::is_float_real_scalar<V>::value;
    
    static void eval_general(M1& C,const M2& A, const M1& B, trans_type t_A, trans_type t_B,
                             Integer fr, Integer rows, Integer fc, Integer cols,
                             const V& alpha, const V& beta)
    {
        if (A.has_diag(0) == false)
        {
            if (mrd::is_one(alpha) == true)
            {
                using Alpha     = alpha_one<V>;
                Alpha al;

                eval_band_dense_generic<V>::eval_impl(al, C, A, B, t_A, t_B, fr, rows, fc, cols);
            }
            else if (mrd::is_one(-alpha) == true)
            {
                using Alpha     = alpha_mone<V>;
                Alpha al;

                eval_band_dense_generic<V>::eval_impl(al, C, A, B, t_A, t_B, fr, rows, fc, cols);
            }
            else
            {
                using Alpha     = alpha_val<V>;
                Alpha al(alpha);

                eval_band_dense_generic<V>::eval_impl(al, C, A, B, t_A, t_B, fr, rows, fc, cols);
            };

            return;
        }

        Integer kl      = A.number_subdiagonals();
        Integer ku      = A.number_superdiagonals();

        Integer B_ld    = B.ld();
        Integer C_ld    = C.ld();

        const V *B_ptr  = B.ptr();
        V *C_ptr        = C.ptr() + fr + fc*C_ld;

        bool V_is_real  = md::is_real_scalar<V>::value;

        if (t_B == trans_type::no_trans)
        {
            for (Integer i = 1; i <= B.cols(); ++i, B_ptr += B_ld, C_ptr += C_ld)
            {
                lapack::gbmv(trans_code<ISR>(t_A), A.rows(), A.cols(), kl, ku, *md::lap(&alpha), 
                    md::lap(A.rep_ptr()), A.ld(), md::lap(B_ptr), 1, *md::lap(&beta), 
                    md::lap(C_ptr), 1);
            };
        }
        else if (t_B == trans_type::trans || (t_B == trans_type::conj_trans && V_is_real))
        {
            for (Integer i = 1; i <= B.rows(); ++i, B_ptr += 1, C_ptr += C_ld)
            {
                lapack::gbmv(trans_code<ISR>(t_A), A.rows(), A.cols(), kl, ku, *md::lap(&alpha), 
                    md::lap(A.rep_ptr()), A.ld(), md::lap(B_ptr), B_ld, *md::lap(&beta), 
                    md::lap(C_ptr), 1);
            };
        }
        else
        {
            using VTR_pod       = matcl::pod_type<V>;
            using workspace     = matcl::pod_workspace<VTR_pod>;
            workspace WORK      = workspace(B.cols());
            V* ptr_w            = reinterpret_cast<V*>(WORK.ptr());

            for (Integer i = 1; i <= B.rows(); ++i, B_ptr += 1, C_ptr += C_ld)
            {
                for (Integer j = 0; j < B.cols(); ++j)
                {
                    ptr_w[j]    = conj(B_ptr[j*B_ld]);
                };

                lapack::gbmv(trans_code<ISR>(t_A), A.rows(), A.cols(), kl, ku, *md::lap(&alpha), 
                        md::lap(A.rep_ptr()), A.ld(), md::lap(ptr_w), 1, *md::lap(&beta), 
                        md::lap(C_ptr), 1);
            };
        };
        return;
    };
    static void eval_her(M1& C,const M2& A, const M1& B, trans_type t_A, trans_type t_B,
                         Integer fr, Integer rows, Integer fc, Integer cols,
                         const V& alpha, const V& beta)
    {
        bool V_is_real = md::is_real_scalar<V>::value;

        if (t_A == trans_type::trans && md::is_float_real_scalar<V>::value == false
            || (t_B == trans_type::conj_trans && V_is_real == false))
        {
            return eval_general(C,A,B,t_A, t_B, fr, rows, fc, cols,alpha,beta);
        };

        Integer ld      = A.last_diag();

        Integer B_ld    = B.ld();
        Integer C_ld    = C.ld();

        const V* B_ptr  = B.ptr();
        V *C_ptr        = C.ptr() + fr + fc * C_ld;

        if (t_B == trans_type::no_trans)
        {
            for (Integer i = 1; i <= B.cols(); ++i)
            {
                lapack::hbmv("U", A.cols(), ld, *md::lap(&alpha), md::lap(A.rep_ptr()), A.ld(), 
                            md::lap(B_ptr), 1, *md::lap(&beta), md::lap(C_ptr), 1);

                B_ptr   += B_ld;
                C_ptr   += C_ld;
            };
        }
        else if (t_B == trans_type::trans || (t_B == trans_type::conj_trans && V_is_real))
        {
            for (Integer i = 1; i <= B.rows(); ++i)
            {
                lapack::hbmv("U", A.cols(), ld, *md::lap(&alpha), md::lap(A.rep_ptr()), A.ld(), 
                            md::lap(B_ptr), B_ld, *md::lap(&beta), md::lap(C_ptr), 1);

                B_ptr   += 1;
                C_ptr   += C_ld;
            };
        }
        else
        {
            //this case should already be removed
        };

        return;
    };

    static void eval_sym(M1& C,const M2& A, const M1& B, trans_type t_A, trans_type t_B,
                         Integer fr, Integer rows, Integer fc, Integer cols,
                         const V& alpha, const V& beta)
    {
        bool V_is_real = md::is_real_scalar<V>::value;

        if (t_A == trans_type::conj_trans && md::is_float_real_scalar<V>::value == false
            ||t_B == trans_type::conj_trans && V_is_real == false)
        {
            return eval_general(C,A,B,t_A, t_B,fr,rows,fc,cols,alpha,beta);
        };

        Integer ld      = A.last_diag();

        Integer B_ld    = B.ld();
        Integer C_ld    = C.ld();

        const V* B_ptr  = B.ptr();
        V *C_ptr        = C.ptr() + fr + fc * C_ld;

        if (t_B == trans_type::no_trans)
        {
            for (Integer i = 1; i <= B.cols(); ++i)
            {
                lapack::sbmv("U", A.cols(), ld, *md::lap(&alpha), md::lap(A.rep_ptr()), A.ld(), 
                             md::lap(B_ptr), 1, *md::lap(&beta), md::lap(C_ptr), 1);

                C_ptr   += C_ld;
                B_ptr   += B_ld;
            };
        }
        else if (t_B == trans_type::trans || t_B == trans_type::conj_trans && V_is_real == true)
        {
            for (Integer i = 1; i <= B.rows(); ++i)
            {
                lapack::sbmv("U", A.cols(), ld, *md::lap(&alpha), md::lap(A.rep_ptr()), A.ld(), 
                            md::lap(B_ptr), B_ld, *md::lap(&beta), md::lap(C_ptr), 1);

                B_ptr   += 1;
                C_ptr   += C_ld;
            };
        }
        else
        {
            //this case should aready be removed
        };
    };

    static void eval_impl(M1& C, const M2& A, const M1& B, trans_type t_A, trans_type t_B,
                          Integer fr, Integer rows, Integer fc, Integer cols,
                          const V& alpha, const V& beta)
    {
        bool is_A_sq    = A.rows() == A.cols();

        if (rows == 1 && cols == 1)
            return eval_dot(C,A,B,t_A, t_B, fr, rows, fc, cols, alpha, beta);

        if (A.get_struct().is_symmetric(is_A_sq, is_real_matrix(A)))
        {
            return eval_sym(C,A,B, t_A, t_B, fr, rows, fc, cols, alpha, beta);
        }
        else if (A.get_struct().is_hermitian(is_A_sq, is_real_matrix(A)))
        {
            return eval_her(C,A,B, t_A, t_B, fr, rows, fc, cols, alpha, beta); 
        };

        return eval_general(C,A,B, t_A, t_B, fr, rows, fc, cols, alpha, beta);
    };

    static void eval_dot(M1& C, const M2& A, const M1& B, trans_type t_A, trans_type t_B,
                            Integer fr, Integer rows, Integer fc, Integer cols, const V& alpha,const V& beta)
    {
        (void)cols;
        (void)rows;

        Integer fd      = A.first_diag();
        Integer ld      = A.last_diag();
        Integer s       = ld - fd + 1;
        Integer lda     = (t_A == trans_type::no_trans)? A.ld() - 1 : 1;
        Integer ldb     = (t_B == trans_type::no_trans)? 1 : B.ld();
        Integer fd_tr   = (t_A == trans_type::no_trans)? fd : ld;
        const V* A_ptr  = A.rep_ptr() + A.first_elem_diag(fd_tr);

        V dot;

        if (t_A != trans_type::conj_trans && t_B != trans_type::conj_trans
            || md::is_float_real_scalar<V>::value == true)
        {
            V res(lapack::dot(s, md::lap(A_ptr), lda, md::lap(B.ptr()), ldb));
            dot = res;
        }
        else if (t_A == trans_type::conj_trans && t_B == trans_type::conj_trans)
        {
            V res(lapack::dot(s, md::lap(A_ptr), lda, md::lap(B.ptr()), ldb));
            dot = conj(res);
        }
        else if (t_A == trans_type::conj_trans)
        {
            V res   = V(lapack::dotc(s, md::lap(A_ptr), lda, md::lap(B.ptr()), ldb));
            dot     = res;
        }
        else
        {
            //t_B = conj_trans
            V res   = V(lapack::dotc(s, md::lap(B.ptr()), ldb, md::lap(A_ptr), lda));
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

    static void eval(matcl::Matrix& ret, const M2& A, const M1& B, trans_type t_A, trans_type t_B)
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

        if (M == 0 || N == 0 || K == 0)
        {
            using sparse_mat = Matrix<V,struct_sparse>;
            sparse_mat out(ti::ti_empty(), M, N);
            ret = matcl::Matrix(out,false);
            return;
        };

        bool is_diag_A  = raw::is_diag(A);
        bool is_diag_B  = raw::is_diag(B);

        Integer fr      = 0;
        Integer rows    = M;
        Integer fc      = 0;
        Integer cols    = N;

        if (is_diag_A)
            return eval_band_dense_diag<V>::eval1(ret,A,B, t_A, t_B);

        if (is_diag_B)
            return eval_band_dense_diag<V>::eval2(ret,A,B, t_A, t_B);

        M1 C(ti::ti_empty(), rows, cols);

        V alpha         = V(1.0);
        V beta          = V(0.0);

        eval_impl(C, A, B, t_A, t_B, fr, rows, fc, cols, alpha, beta);

        bool is_sq_A    = (A.rows() == A.cols());
        bool is_sq_B    = (B.rows() == B.cols());
        bool is_sq_C    = (C.rows() == C.cols());
        C.set_struct(predefined_struct_ext::mult_struct(A.get_struct(),B.get_struct(),t_A,t_B,
                                is_real_matrix(A),is_real_matrix(B), is_sq_A, is_sq_B, is_sq_C));

        ret = matcl::Matrix(C,true);
        return;
    };

    static void eval_gemm(M1& C, const V& alpha, const M2& A, const M1& B, 
                     trans_type t_A, trans_type t_B, const V& beta, 
                     Integer fr, Integer rows, Integer fc, Integer cols)
    {
        if (rows == 0 || cols == 0)
            return;

        if (A.rows() == 0 || A.cols() == 0 || B.rows() == 0 || B.cols() == 0)
        {
            prepare_gemm_C<V>::eval(beta, C, fr,rows,fc,cols);
            return;
        }

        bool is_diag_A  = raw::is_diag(A);
        bool is_diag_B  = raw::is_diag(B);

        if (is_diag_A)
        {
            prepare_gemm_C<V>::eval(beta, C, fr,rows,fc,cols);
            return eval_band_dense_diag<V>::eval1_gemm(alpha,C,A,B, t_A, t_B, fr, rows, fc, cols);
        }

        if (is_diag_B)
        {
            prepare_gemm_C<V>::eval(beta, C, fr,rows,fc,cols);
            return eval_band_dense_diag<V>::eval2_gemm(alpha,C,A,B, t_A, t_B, fr, rows, fc, cols);
        }

        eval_impl(C, A, B, t_A, t_B, fr, rows, fc, cols, alpha, beta);
    };
};

template<class Val_ret>
struct eval_band_dense{};

template<> struct eval_band_dense<Integer>
{
    using VT        = Integer;
    using Mat_B     = raw::Matrix<VT, struct_banded>;
    using Mat_D     = raw::Matrix<VT, struct_dense>;

    static void eval(matcl::Matrix& ret, const Mat_B& A,const Mat_D& B,
                     trans_type t_A, trans_type t_B)
    {
        return eval_band_dense_generic<VT>::eval(ret,A,B, t_A, t_B);
    };
    static void eval_gemm(Mat_D& C, const VT& alpha, const Mat_B& A, const Mat_D& B, 
                     trans_type t_A, trans_type t_B, const VT& beta, 
                     Integer fr, Integer rows, Integer fc, Integer cols)
    {
        return eval_band_dense_generic<VT>::eval_gemm(C,alpha,A,B, t_A, t_B,beta,fr,rows,fc,cols);
    };
}; 

template<> struct eval_band_dense<Object>
{
    using VT        = Object;
    using Mat_B     = raw::Matrix<VT, struct_banded>;
    using Mat_D     = raw::Matrix<VT, struct_dense>;

    static void eval(matcl::Matrix& ret, const Mat_B& A,const Mat_D& B,
                     trans_type t_A, trans_type t_B)
    {
        return eval_band_dense_generic<VT>::eval(ret,A,B, t_A, t_B);
    };
    static void eval_gemm(Mat_D& C, const VT& alpha, const Mat_B& A, const Mat_D& B, 
                     trans_type t_A, trans_type t_B, const VT& beta, 
                     Integer fr, Integer rows, Integer fc, Integer cols)
    {
        return eval_band_dense_generic<VT>::eval_gemm(C,alpha,A,B, t_A, t_B,beta,fr,rows,fc,cols);
    };
}; 

template<> struct eval_band_dense<Real>
{
    using VT        = Real;
    using Mat_B     = raw::Matrix<VT, struct_banded>;
    using Mat_D     = raw::Matrix<VT, struct_dense>;

    static void eval(matcl::Matrix& ret, const Mat_B& A,const Mat_D& B,
                     trans_type t_A, trans_type t_B)
    {
        return eval_band_dense_lapack<VT>::eval(ret,A,B,t_A,t_B);
    };
    static void eval_gemm(Mat_D& C, const VT& alpha, const Mat_B& A, const Mat_D& B, 
                     trans_type t_A, trans_type t_B, const VT& beta, 
                     Integer fr, Integer rows, Integer fc, Integer cols)
    {
        return eval_band_dense_lapack<VT>::eval_gemm(C,alpha,A,B, t_A, t_B,beta,fr,rows,fc,cols);
    };
};

template<> struct eval_band_dense<Float>
{
    using VT        = Float;
    using Mat_B     = raw::Matrix<VT, struct_banded>;
    using Mat_D     = raw::Matrix<VT, struct_dense>;

    static void eval(matcl::Matrix& ret, const Mat_B& A,const Mat_D& B,
                     trans_type t_A, trans_type t_B)
    {
        return eval_band_dense_lapack<VT>::eval(ret,A,B,t_A,t_B);
    };
    static void eval_gemm(Mat_D& C, const VT& alpha, const Mat_B& A, const Mat_D& B, 
                     trans_type t_A, trans_type t_B, const VT& beta, 
                     Integer fr, Integer rows, Integer fc, Integer cols)
    {
        return eval_band_dense_lapack<VT>::eval_gemm(C,alpha,A,B, t_A, t_B,beta,fr,rows,fc,cols);
    };
};

template<> struct eval_band_dense<Complex>
{
    using VT        = Complex;
    using Mat_B     = raw::Matrix<VT, struct_banded>;
    using Mat_D     = raw::Matrix<VT, struct_dense>;

    static void eval(matcl::Matrix& ret, const Mat_B& A,const Mat_D& B,
                     trans_type t_A, trans_type t_B)
    {
        return eval_band_dense_lapack<VT>::eval(ret,A,B,t_A,t_B);
    };
    static void eval_gemm(Mat_D& C, const VT& alpha, const Mat_B& A, const Mat_D& B, 
                     trans_type t_A, trans_type t_B, const VT& beta, 
                     Integer fr, Integer rows, Integer fc, Integer cols)
    {
        return eval_band_dense_lapack<VT>::eval_gemm(C,alpha,A,B, t_A, t_B,beta,fr,rows,fc,cols);
    };
};

template<> struct eval_band_dense<Float_complex>
{
    using VT        = Float_complex;
    using Mat_B     = raw::Matrix<VT, struct_banded>;
    using Mat_D     = raw::Matrix<VT, struct_dense>;

    static void eval(matcl::Matrix& ret, const Mat_B& A,const Mat_D& B,
                     trans_type t_A, trans_type t_B)
    {
        return eval_band_dense_lapack<VT>::eval(ret,A,B,t_A,t_B);
    };
    static void eval_gemm(Mat_D& C, const VT& alpha, const Mat_B& A, const Mat_D& B, 
                     trans_type t_A, trans_type t_B, const VT& beta, 
                     Integer fr, Integer rows, Integer fc, Integer cols)
    {
        return eval_band_dense_lapack<VT>::eval_gemm(C,alpha,A,B, t_A, t_B,beta,fr,rows,fc,cols);
    };
};

template<class Val_ret, class M1, class M2>
struct eval_mult<Val_ret,M1,M2,struct_banded,struct_dense> 
{ 
    using Mat_D = raw::Matrix<Val_ret, struct_dense>;

    static void eval(matcl::Matrix& ret, const M1& A, const M2& B, trans_type t_A, trans_type t_B)
    {
        using MP1           = Matrix<Val_ret,struct_banded>;
        using MP2           = Matrix<Val_ret,struct_dense>;

        return eval_band_dense<Val_ret>::eval(ret, converter<MP1,M1>::eval(A),
                                                converter<MP2,M2>::eval(B), t_A, t_B);
    };

    static void eval_gemm(Mat_D& C, const Val_ret& alpha, const M1& A, const M2& B, 
                     trans_type t_A, trans_type t_B, const Val_ret& beta, 
                     Integer fr, Integer rows, Integer fc, Integer cols)
    {
        using MP1           = Matrix<Val_ret,struct_banded>;
        using MP2           = Matrix<Val_ret,struct_dense>;

        return eval_band_dense<Val_ret>::eval_gemm(C, alpha, converter<MP1,M1>::eval(A),
                    converter<MP2,M2>::eval(B), t_A, t_B, beta, fr, rows, fc, cols);
    };
};

template<class Val_C, class M1, class M2>
struct eval_mult_abs<Val_C, M1, M2, struct_banded>
{ 
    using Mat_C = raw::Matrix<Val_C, struct_dense>;

    static void eval_scal(matcl::Matrix& out, const M1& A, const M2& scal, trans_type t_A, const Mat_C& C)
    {
        if (C.rows() == 0 || C.cols() == 0)
        {
            out = matcl::Matrix(C, false);
            return;
        }

        if (mrd::is_zero<M2>(scal) == true)
        {
            mrd::scalfunc_real_helper<Mat_C>::eval_abs(out, C);
            return;
        };

        bool is_diag_A  = raw::is_diag(A);

        if (is_diag_A)
            return eval_diag_scal(out, A, scal, t_A, C);

        Integer M       = C.rows();
        Integer N       = C.cols();

        using V1        = typename M1::value_type;
        using V2        = M2;
        using V1R       = typename md::real_type<V1>::type;
        using VR        = typename md::real_type<Val_C>::type;
        using Mat_ret   = raw::Matrix<VR, struct_dense>;

        Mat_ret ret(C.get_type(), M, N);

        Integer A_ld    = A.ld();
        Integer C_ld    = C.ld();
        Integer ret_ld  = ret.ld();

        const Val_C* ptr_C  = C.ptr();
        VR* ptr_ret         = ret.ptr();

        if (t_A == trans_type::no_trans)
        {            
            const V1* ptr_A = A.rep_ptr();

            for (Integer j = 0; j < N; ++j) 
            {
                for (Integer i = 0; i < M; ++i)
                    ptr_ret[i]  = abs_helper<Val_C>::eval(ptr_C[i]);                

                Integer fr      = A.first_row(j);
                Integer lr      = A.last_row(j);
                Integer pos_A   = A.first_elem_pos(j);

                for (Integer l = fr; l <= lr; ++l, ++pos_A) 
                    ptr_ret[l]  = ptr_ret[l] + scal * abs_helper<V1>::eval(ptr_A[pos_A]);

                ptr_A       += A_ld;
                ptr_C       += C_ld;
                ptr_ret     += ret_ld;
            };
        }
        else
        {
            const V1* ptr_A = A.rep_ptr();

            for (Integer j = 0; j < N; ++j) 
            {
                for (Integer i = 0; i < M; ++i)
                    ptr_ret[i]  = abs_helper<Val_C>::eval(ptr_C[i]);                

                Integer fc      = A.first_col(j);
                Integer lc      = A.last_col(j);
                Integer pos_A   = A.first_elem_pos_row(j);

                for (Integer l = fc; l <= lc; ++l, pos_A += A_ld - 1) 
                    ptr_ret[l]  = ptr_ret[l] + scal * abs_helper<V1>::eval(ptr_A[pos_A]);

                ptr_C       += C_ld;
                ptr_ret     += ret_ld;
            };
        }

        out = matcl::Matrix(ret,false);
    };

    static void eval(matcl::Matrix& out, const M1& A, const M2& X, trans_type t_A, const Mat_C& C)
    {
        if (C.rows() == 0 || C.cols() == 0)
        {
            out = matcl::Matrix(C, false);
            return;
        }

        if (A.rows() == 0 || A.cols() == 0)
        {
            mrd::scalfunc_real_helper<Mat_C>::eval_abs(out, C);
            return;
        }

        bool is_diag_A  = raw::is_diag(A);

        if (is_diag_A)
            return eval_diag(out, A, X, t_A, C);

        Integer M       = C.rows();
        Integer N       = C.cols();        

        using V1        = typename M1::value_type;
        using V2        = typename M2::value_type;
        using V1R       = typename md::real_type<V1>::type;
        using VR        = typename md::real_type<Val_C>::type;
        using Mat_ret   = raw::Matrix<VR, struct_dense>;

        Mat_ret ret(C.get_type(), M, N);

        Integer A_ld    = A.ld();
        Integer X_ld    = X.ld();
        Integer C_ld    = C.ld();
        Integer ret_ld  = ret.ld();

        const V2* ptr_X     = X.ptr();        
        const Val_C* ptr_C  = C.ptr();
        VR* ptr_ret         = ret.ptr();

        if (t_A == trans_type::no_trans)
        {            
            Integer K   = A.cols();

            for (Integer j = 0; j < N; ++j) 
            {
                //make abs(C);
                for (Integer i = 0; i < M; ++i)
                    ptr_ret[i]  = abs_helper<Val_C>::eval(ptr_C[i]);

                const V1* ptr_A = A.rep_ptr();

                for (Integer l = 0; l < K; ++l) 
                {
                    const VR& temp = abs_helper<V2>::eval(ptr_X[l]);

                    if (mrd::is_zero(temp) == false) 
                    {
                        Integer first_row   = A.first_row(l);
                        Integer last_row    = A.last_row(l);
                        Integer pos_A       = A.first_elem_pos(l);

                        for (Integer k = first_row; k <= last_row; ++k) 
                        {
                            VR tmp          = abs_helper<V1>::eval(ptr_A[pos_A++]) * temp;
                            ptr_ret[k]      = ptr_ret[k] + tmp;
                        };
                    };

                    ptr_A += A_ld;
                };

                ptr_X   += X_ld;
                ptr_C   += C_ld;
                ptr_ret += ret_ld;
            };
        }
        else
        {
            for (Integer j = 0; j < N; ++j) 
            {
                const V1* ptr_A = A.rep_ptr();

                for (Integer l = 0; l < M; ++l) 
                {
                    VR dot = abs_helper<Val_C>::eval(ptr_C[l]);
        
                    Integer first_row   = A.first_row(l);
                    Integer last_row    = A.last_row(l);
                    Integer pos_A       = A.first_elem_pos(l);

                    for (Integer k = first_row; k <= last_row; ++k) 
                    {
                        VR tmp  = mrd::abs_helper<V1>::eval(ptr_A[pos_A++]) * mrd::abs_helper<V2>::eval(ptr_X[k]);
                        dot     = dot + tmp;
                    };

                    ptr_ret[l]  = dot;
                    ptr_A       += A_ld;
                };

                ptr_X   += X_ld;
                ptr_C   += C_ld;
                ptr_ret += ret_ld;
            };
        }

        out = matcl::Matrix(ret,false);
    };

    static void eval_diag_scal(matcl::Matrix& out, const M1& A, const M2& scal, trans_type t_A, const Mat_C& C)
    {
        (void)t_A;

        Integer M       = C.rows();
        Integer N       = C.cols();

        Integer K1      = std::min(A.rows(),A.cols());

        using V1        = typename M1::value_type;
        using V2        = M2;
        using V1R       = typename md::real_type<V1>::type;
        using VR        = typename md::real_type<Val_C>::type;
        using Mat_ret   = raw::Matrix<VR, struct_dense>;

        Mat_ret ret(C.get_type(), M, N);

        Integer A_ld    = A.ld();
        Integer C_ld    = C.ld();
        Integer ret_ld  = ret.ld();

        const V1* A_ptr = A.rep_ptr() + A.first_elem_diag(0);

        const Val_C* C_ptr  = C.ptr();
        VR* ret_ptr         = ret.ptr();

        for (Integer j = 0; j < K1; ++j)
        {
            for (Integer i = 0; i < M; ++i)
                ret_ptr[i]  = abs_helper<Val_C>::eval(C_ptr[i]);

            ret_ptr[j]      = ret_ptr[j] + scal * abs_helper<V1>::eval(A_ptr[0]);

            C_ptr       += C_ld;
            ret_ptr     += ret_ld;
            A_ptr       += A_ld;
        };
        for (Integer j = K1; j < N; ++j)
        {
            for (Integer i = 0; i < M; ++i)
                ret_ptr[i]  = abs_helper<Val_C>::eval(C_ptr[i]);

            C_ptr       += C_ld;
            ret_ptr     += ret_ld;
        };

        out = matcl::Matrix(ret,false);
        return;
    };

    static void eval_diag(matcl::Matrix& out, const M1& A, const M2& X, trans_type t_A, const Mat_C& C)
    {
        (void)t_A;

        Integer M       = C.rows();
        Integer N       = C.cols();

        Integer K1      = std::min(A.rows(),A.cols());

        using V1        = typename M1::value_type;
        using V2        = typename M2::value_type;
        using V1R       = typename md::real_type<V1>::type;
        using V2R       = typename md::real_type<V2>::type;
        using VR        = typename md::real_type<Val_C>::type;

        using DM        = raw::Matrix<V1R, struct_dense>;
        using Mat_ret   = raw::Matrix<VR, struct_dense>;

        Mat_ret ret(C.get_type(), M, N);


        DM diag(A.get_type(), K1, 1);

        Integer A_ld    = A.ld();
        Integer X_ld    = X.ld();
        Integer C_ld    = C.ld();
        Integer ret_ld  = ret.ld();

        const V1* A_ptr = A.rep_ptr() + A.first_elem_diag(0);
        const V2* X_ptr = X.ptr();
        V1R* D_ptr      = diag.ptr();

        const Val_C* C_ptr  = C.ptr();
        VR* ret_ptr         = ret.ptr();

        for(Integer i = 0; i < K1; ++i)
        {
            D_ptr[i]    = abs_helper<V1>::eval(A_ptr[0]);
            A_ptr       += A_ld;
        };

        for (Integer j = 0; j < N; ++j)
        {
            for (Integer i = 0; i < K1; ++i)
            {
                ret_ptr[i]  = D_ptr[i] * abs_helper<V2>::eval(X_ptr[i])
                            + abs_helper<Val_C>::eval(C_ptr[i]);
            }

            for (Integer i = K1; i < M; ++i)
                ret_ptr[i]  = abs_helper<Val_C>::eval(C_ptr[i]);

            C_ptr       += C_ld;
            X_ptr       += X_ld;
            ret_ptr     += ret_ld;
        };

        out = matcl::Matrix(ret,false);
        return;
    };
};

}}}