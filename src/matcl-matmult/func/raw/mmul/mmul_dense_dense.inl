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

#include "matcl-internals/func/raw_func_unary.h"

namespace matcl { namespace raw { namespace details
{

//========================================================================
//                      DENSE - DENSE
//========================================================================
template<class V>
struct eval_dense_diag
{
    using DM    = raw::Matrix<V,struct_dense>;

    static void eval_diag_diag(matcl::Matrix& ret, const DM& A, const DM& B, trans_type t_A, trans_type t_B)
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

        DM C(ret_ti, Z, MN, 1);

        const V* A_ptr  = A.ptr();
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
                    A_ptr       += A_ld + 1;
                    B_ptr       += B_ld + 1;
                };
            }
            else
            {
                for (Integer i = 0; i < K; ++i)
                {
                    C_ptr[i]    = conj(A_ptr[0]) * B_ptr[0];
                    A_ptr       += A_ld + 1;
                    B_ptr       += B_ld + 1;
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
                    A_ptr       += A_ld + 1;
                    B_ptr       += B_ld + 1;
                };
            }
            else
            {
                for (Integer i = 0; i < K; ++i)
                {
                    C_ptr[i]    = A_ptr[0] * B_ptr[0];
                    A_ptr       += A_ld + 1;
                    B_ptr       += B_ld + 1;
                };
            };
        };

        ret = matcl::bdiags(matcl::Matrix(C,false), Integer(0), M, N);
        return;
    };

    static void eval_diag_diag_gemm(const V& alpha, DM& C, const DM& A, const DM& B, trans_type t_A, trans_type t_B, 
                     Integer fr, Integer rows, Integer fc, Integer cols)
    {        
        Integer MN  = std::min(rows,cols);
        Integer K   = t_A == trans_type::no_trans ? A.cols() : A.rows();
        K           = std::min(MN,K);

        Integer A_ld    = A.ld();
        Integer B_ld    = B.ld();
        Integer C_ld    = C.ld();

        const V* A_ptr  = A.ptr();
        const V* B_ptr  = B.ptr();
        V* C_ptr        = C.ptr() + fr + fc * C_ld;

        if (t_A == trans_type::conj_trans)
        {
            if (t_B == trans_type::conj_trans)
            {
                for (Integer i = 0; i < K; ++i)
                {
                    V val       = conj(A_ptr[0]) * conj(B_ptr[0]);
                    C_ptr[0]    = C_ptr[0] + alpha * val;
                    A_ptr       += A_ld + 1;
                    B_ptr       += B_ld + 1;
                    C_ptr       += C_ld + 1;
                };
            }
            else
            {
                for (Integer i = 0; i < K; ++i)
                {
                    V val       = conj(A_ptr[0]) * B_ptr[0];
                    C_ptr[0]    = C_ptr[0] + alpha * val;
                    A_ptr       += A_ld + 1;
                    B_ptr       += B_ld + 1;
                    C_ptr       += C_ld + 1;
                };
            };
        }
        else
        {
            if (t_B == trans_type::conj_trans)
            {
                for (Integer i = 0; i < K; ++i)
                {
                    V val       = A_ptr[0] * conj(B_ptr[0]);
                    C_ptr[0]    = C_ptr[0] + alpha * val;
                    A_ptr       += A_ld + 1;
                    B_ptr       += B_ld + 1;
                    C_ptr       += C_ld + 1;
                };
            }
            else
            {
                for (Integer i = 0; i < K; ++i)
                {
                    V val       = A_ptr[0] * B_ptr[0];
                    C_ptr[0]    = C_ptr[0] + alpha * val;
                    A_ptr       += A_ld + 1;
                    B_ptr       += B_ld + 1;
                    C_ptr       += C_ld + 1;
                };
            };
        };
    };
    
    static void eval1_gemm(const V& alpha, DM& C, const DM& A, const DM& B, trans_type t_A, trans_type t_B, 
                     Integer fr, Integer rows, Integer fc, Integer cols)
    {
        if (raw::is_diag(B))
            return eval_diag_diag_gemm(alpha,C,A,B,t_A, t_B,fr,rows,fc,cols);

        eval1_impl(alpha, C, A, B, t_A, t_B, fr, rows, fc, cols);
    };

    static void eval1(matcl::Matrix& ret, const DM& A, const DM& B, trans_type t_A, trans_type t_B)
    {
        if (raw::is_diag(B))
            return eval_diag_diag(ret,A,B,t_A, t_B);

        Integer M, N;

        if (t_B == trans_type::no_trans)
            N = B.cols();
        else
            N = B.rows();

        if (t_A == trans_type::no_trans)
            M   = A.rows();
        else
            M   = A.cols();

        ti::ti_type<V> ret_ti = ti::get_return_ti<ti::ti_type<V>>(mdyf::op_mul::eval(),
                                                                  ti::get_ti(A), ti::get_ti(B));
        V Z = md::default_value<V>(ret_ti);

        DM C(ret_ti, Z, M, N);

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

    static void eval1_impl(const V& alpha, DM& C, const DM& A, const DM& B, trans_type t_A, trans_type t_B,
                           Integer fr, Integer rows, Integer fc, Integer cols)
    {
        (void)rows;
        Integer K1 = std::min(A.rows(),A.cols());

        DM diag(A.get_type(), K1, 1);

        Integer A_ld    = A.ld();
        Integer B_ld    = B.ld();
        Integer C_ld    = C.ld();

        const V* A_ptr  = A.ptr();
        const V* B_ptr  = B.ptr();
        V* D_ptr        = diag.ptr();
        V* C_ptr        = C.ptr() + fr + fc * C_ld;

        Integer N       = cols;

        if (t_A == trans_type::conj_trans)
        {
            for(Integer i = 0; i < K1; ++i)
            {
                D_ptr[i]    = alpha * mrd::conj_helper<V>::eval(A_ptr[0]);
                A_ptr       += A_ld + 1;
            };
        }
        else
        {
            for(Integer i = 0; i < K1; ++i)
            {
                D_ptr[i]    = alpha * A_ptr[0];
                A_ptr       += A_ld + 1;
            };
        };

        if (t_B == trans_type::no_trans)
        {
            for (Integer j = 0; j < N; ++j)
            {
                for (Integer i = 0; i < K1; ++i)
                    C_ptr[i]    = C_ptr[i] + D_ptr[i] * B_ptr[i];

                C_ptr       += C_ld;
                B_ptr       += B_ld;
            };
        }
        else if (t_B == trans_type::trans)
        {
            for (Integer j = 0; j < N; ++j)
            {
                for (Integer i = 0; i < K1; ++i)
                    C_ptr[i]    = C_ptr[i] + D_ptr[i] * B_ptr[j + i*B_ld];

                C_ptr       += C_ld;
            };
        }
        else
        {
            for (Integer j = 0; j < N; ++j)
            {
                for (Integer i = 0; i < K1; ++i)
                    C_ptr[i]    = C_ptr[i] + D_ptr[i] * conj(B_ptr[j + i*B_ld]);

                C_ptr       += C_ld;
            };
        };
    };

    static void eval2_gemm(const V& alpha, DM& C, const DM& A, const DM& B, trans_type t_A, trans_type t_B, 
                     Integer fr, Integer rows, Integer fc, Integer cols)
    {
        if (raw::is_diag(A))
            return eval_diag_diag_gemm(alpha,C,A,B, t_A, t_B, fr,rows,fc,cols);

        eval2_impl(alpha, C, A, B, t_A, t_B, fr, rows, fc, cols);
    };

    static void eval2(matcl::Matrix& ret, const DM& A, const DM& B, trans_type t_A, trans_type t_B)
    {        
        if (raw::is_diag(A))
            return eval_diag_diag(ret,A,B, t_A, t_B);        

        Integer M, N;

        if (t_B == trans_type::no_trans)
            N = B.cols();
        else
            N = B.rows();

        if (t_A == trans_type::no_trans)
            M       = A.rows();
        else
            M       = A.cols();

        ti::ti_type<V> ret_ti = ti::get_return_ti<ti::ti_type<V>>(mdyf::op_mul::eval(),
                                                                  ti::get_ti(A), ti::get_ti(B));
        V Z = md::default_value<V>(ret_ti);

        DM C(ret_ti, Z, M, N);

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
        return;
    };    

    static void eval2_impl(const V& alpha, DM& C, const DM& A, const DM& B, trans_type t_A, trans_type t_B,
                           Integer fr, Integer rows, Integer fc, Integer cols)
    {        
        (void)cols;

        Integer K1  = std::min(B.rows(),B.cols());

        DM diag(B.get_type(), K1, 1);

        Integer A_ld    = A.ld();
        Integer B_ld    = B.ld();
        Integer C_ld    = C.ld();
        Integer M       = rows;

        const V* A_ptr  = A.ptr();
        const V* B_ptr  = B.ptr();
        V* C_ptr        = C.ptr() + fr + fc * C_ld;
        V* D_ptr        = diag.ptr();
        const V* As_ptr = A_ptr;

        if (t_B == trans_type::conj_trans)
        {
            for(Integer i = 0; i < K1; ++i)
            {
                D_ptr[i]    = alpha*conj(B_ptr[0]);
                B_ptr       += B_ld + 1;
            };
        }
        else
        {
            for(Integer i = 0; i < K1; ++i)
            {
                D_ptr[i]    = alpha*B_ptr[0];
                B_ptr       += B_ld + 1;
            };
        };

        if (t_A == trans_type::no_trans)
        {
            for (Integer j = 0; j < K1; ++j)
            {
                const V& tmp = D_ptr[j];
                if (mrd::is_zero(tmp))
                {
                    C_ptr   += C_ld;
                    A_ptr   += A_ld;
                    continue;
                };

                for (Integer i = 0; i < M; ++i)
                    C_ptr[i]    = C_ptr[i] + A_ptr[i] * tmp;

                C_ptr       += C_ld;
                A_ptr       += A_ld;
            };
        }
        else if (t_A == trans_type::trans)
        {
            for (Integer j = 0; j < K1; ++j) 
            {	
                const V& tmp = D_ptr[j];

                if (mrd::is_zero(tmp))
                {
                    C_ptr   += C_ld;
                    continue;
                };

                A_ptr = As_ptr;
                for (Integer l = 0; l < M; ++l) 
                {
                    C_ptr[l] = C_ptr[l] + A_ptr[j] * tmp;
                    A_ptr   += A_ld;
                };

                C_ptr += C_ld;
            };
        }
        else
        {
            for (Integer j = 0; j < K1; ++j) 
            {	
                const V& tmp = D_ptr[j];

                if (mrd::is_zero(tmp))
                {
                    C_ptr   += C_ld;
                    continue;
                };

                A_ptr = As_ptr;
                for (Integer l = 0; l < M; ++l) 
                {
                    V a         = mrd::conj_helper<V>::eval(A_ptr[j]);
                    C_ptr[l]    = C_ptr[l] + a * tmp;
                    A_ptr       += A_ld;
                };

                C_ptr += C_ld;
            };
        };
    };    
};

template<class Val_ret>
struct eval_dense_dense{};

template<class V>
struct eval_dense_dense_generic
{
    using DM = raw::Matrix<V,struct_dense>;

    static void eval_gemm(DM& C, const V& alpha, const DM& A, const DM& B, 
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
            return eval_dense_diag<V>::eval1_gemm(alpha, C, A, B, t_A, t_B, fr, rows, fc, cols);

        if (is_diag_B)
            return eval_dense_diag<V>::eval2_gemm(alpha, C, A, B, t_A, t_B, fr, rows, fc, cols);

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
    };

    static void eval(matcl::Matrix& ret, const DM& A, const DM& B, trans_type t_A, trans_type t_B)
    {        
        Integer N, M, K;

        if (t_B == trans_type::no_trans)
            N   = B.cols();
        else
            N   = B.rows();

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

        if (M == 0 || K == 0 || N == 0)
        {
            using VR            = typename md::real_type<V>::type;
            using sparse_mat    = Matrix<VR,struct_sparse>;
            sparse_mat out(ret_ti,M, N);

            ret = matcl::Matrix(out,false);
            return;
        };

        if (raw::is_diag(A))
            return eval_dense_diag<V>::eval1(ret,A,B,t_A,t_B);

        if (raw::is_diag(B))
            return eval_dense_diag<V>::eval2(ret,A,B,t_A,t_B);

        V Z = md::default_value<V>(ret_ti);

        DM out(ret_ti,Z, M, N);

        Integer fr      = 0;
        Integer rows    = M;
        Integer fc      = 0;
        Integer cols    = N;

        using Alpha     = alpha_one<V>;
        Alpha alpha;

        eval_impl(alpha, out, A, B, t_A, t_B, fr, rows, fc, cols);

        bool is_sq_A    = (A.rows() == A.cols());
        bool is_sq_B    = (B.rows() == B.cols());
        bool is_sq_C    = (out.rows() == out.cols());
        out.set_struct(predefined_struct_ext::mult_struct(A.get_struct(),B.get_struct(),t_A,t_B,
                                is_real_matrix(A),is_real_matrix(B), is_sq_A, is_sq_B, is_sq_C));

        ret = matcl::Matrix(out,true);
        return;
    };

    template<class Alpha>
    static void eval_impl(const Alpha& alpha, DM& C, const DM& A, const DM& B, trans_type t_A, trans_type t_B,
                          Integer fr, Integer rows, Integer fc, Integer cols)
    {        
        Integer A_ld    = A.ld();
        Integer B_ld    = B.ld();
        Integer C_ld    = C.ld();

        const V* ptr_A  = A.ptr();
        const V* ptr_B  = B.ptr();
        V* ptr_C        = C.ptr() + fr + fc * C_ld;
        const V* ptr_As = ptr_A;
        V* ptr_Cs       = ptr_C;

        Integer N       = cols;
        Integer M       = rows;
        
        V Z = md::default_value<V>(C.get_type());

        if (t_A == trans_type::no_trans)
        {
            Integer K   = A.cols();

            if (t_B == trans_type::no_trans)
            {
                // C = A*B
                for (Integer j = 0; j < N; ++j) 
                {			
                    ptr_A       = ptr_As;

                    for (Integer l = 0; l < K; ++l) 
                    {
                        const V& temp = ptr_B[l];

                        if (mrd::is_zero(temp) == false) 
                        {					
                            for (Integer k = 0; k < M; ++k) 
                            {
                                V tmp       = ptr_A[k] * temp;
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
                //C = A * B'

                for (Integer l = 0; l < K; ++l)
                {
                    ptr_C   = ptr_Cs;

	                for (Integer j = 0; j < N; ++j)	
                    {
                        const V& temp = (t_B == trans_type::trans) ? ptr_B[j] : conj(ptr_B[j]);

                        if (mrd::is_zero(temp) == false) 
                        {
		                    for (Integer i = 0; i < M; ++i)
                            {
                                V tmp       = ptr_A[i] * temp;
                                ptr_C[i]    = alpha.eval(ptr_C[i],tmp);
                            };
                        };

                        ptr_C   += C_ld;
                    };

                    ptr_A   += A_ld;
                    ptr_B   += B_ld;
                };
            };
        }
        else if (t_A == trans_type::trans)
        {            
            Integer K   = A.rows();

            if (t_B == trans_type::no_trans)
            {
                //C = A'*B
                for (Integer j = 0; j < N; ++j) 
                {			
                    ptr_A = ptr_As;
                    for (Integer l = 0; l < M; ++l) 
                    {
                        V dot = Z;
                        for (Integer k = 0; k < K; ++k) 
                        {
                            V tmp   = ptr_A[k] * ptr_B[k];
                            dot     = dot + tmp;
                        };

                        ptr_C[l]    = alpha.eval(ptr_C[l],dot);
                        ptr_A       += A_ld;
                    };

                    ptr_B   += B_ld;
                    ptr_C   += C_ld;
                };
            }
            else if (t_B == trans_type::trans)
            {
                //C=  A' * B'
                for(Integer j = 0; j < N; ++j)
                {
                    ptr_A = ptr_As;

                    for (Integer i = 0; i < M; ++i)
                    {
                        V dot = Z;
                        for (Integer l = 0; l < K; ++l)
                        {
                            V tmp   = ptr_A[l] * ptr_B[j + l*B_ld];
                            dot     = dot + tmp;
                        };

                        ptr_C[i]    = alpha.eval(ptr_C[i],dot);
                        ptr_A       += A_ld;
                    };

                    ptr_C   += C_ld;
                };
            }
            else
            {
                //C=  A' * B'
                for(Integer j = 0; j < N; ++j)
                {
                    ptr_A = ptr_As;

                    for (Integer i = 0; i < M; ++i)
                    {
                        V dot = Z;
                        for (Integer l = 0; l < K; ++l)
                        {
                            V tmp   = ptr_A[l] * conj(ptr_B[j + l*B_ld]);
                            dot     = dot + tmp;
                        };

                        ptr_C[i]    = alpha.eval(ptr_C[i],dot);
                        ptr_A       += A_ld;
                    };

                    ptr_C   += C_ld;
                };
            };
        }
        else
        {
            Integer K   = A.rows();

            if (t_B == trans_type::no_trans)
            {
                for (Integer j = 0; j < N; ++j) 
                {			
                    ptr_A = ptr_As;
                    for (Integer l = 0; l < M; ++l) 
                    {
                        V dot = Z;		
                        for (Integer k = 0; k < K; ++k) 
                        {
                            V Ac    = mrd::conj_helper<V>::eval(ptr_A[k]);
                            V tmp   = Ac * ptr_B[k];
                            dot     = dot + tmp;
                        };

                        ptr_C[l]    = alpha.eval(ptr_C[l],dot);
                        ptr_A   += A_ld;
                    };

                    ptr_B += B_ld;
                    ptr_C += C_ld;
                };
            }
            else if (t_B == trans_type::trans)
            {
                //C=  A' * B'
                for(Integer j = 0; j < N; ++j)
                {
                    ptr_A = ptr_As;

                    for (Integer i = 0; i < M; ++i)
                    {
                        V dot = Z;
                        for (Integer l = 0; l < K; ++l)
                        {
                            V tmp   = conj(ptr_A[l]) * ptr_B[j + l*B_ld];
                            dot     = dot + tmp;
                        };

                        ptr_C[i]    = alpha.eval(ptr_C[i],dot);
                        ptr_A       += A_ld;
                    };

                    ptr_C   += C_ld;
                };
            }
            else
            {
                //C=  A' * B'
                for(Integer j = 0; j < N; ++j)
                {
                    ptr_A = ptr_As;

                    for (Integer i = 0; i < M; ++i)
                    {
                        V dot = Z;
                        for (Integer l = 0; l < K; ++l)
                        {
                            V tmp   = conj(ptr_A[l]) * conj(ptr_B[j + l*B_ld]);
                            dot     = dot + tmp;
                        };

                        ptr_C[i]    = alpha.eval(ptr_C[i],dot);
                        ptr_A       += A_ld;
                    };

                    ptr_C   += C_ld;
                };
            };
        };
    };
};

template<class V>
struct eval_dense_dense_lapack
{
    using DM    = raw::Matrix<V,struct_dense>;
    static const bool ISR = md::is_float_real_scalar<V>::value;

    static void eval_gemm(DM& C, const V& alpha, const DM& A, const DM& B, 
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

        Integer A_ld = raw::get_ld(A,0);
        Integer A_ud = raw::get_ud(A,0);
        Integer B_ld = raw::get_ld(B,0);
        Integer B_ud = raw::get_ud(B,0);

        if (A_ld == 0 && A_ud == 0)
        {
            prepare_gemm_C<V>::eval(beta, C, fr,rows,fc,cols);
            return eval_dense_diag<V>::eval1_gemm(alpha, C, A, B, t_A, t_B, fr, rows, fc, cols);
        }

        if (B_ld == 0 && B_ud == 0)
        {
            prepare_gemm_C<V>::eval(beta, C, fr,rows,fc,cols);
            return eval_dense_diag<V>::eval2_gemm(alpha, C, A, B, t_A, t_B, fr, rows, fc, cols);
        };

        eval_impl(C, A, B, t_A, t_B, A_ld, A_ud, B_ld, B_ud, fr, rows, fc, cols, alpha, beta);
    };

    static void eval(matcl::Matrix& ret, const DM& A, const DM& B, trans_type t_A, trans_type t_B)
    {
        Integer M, K, N;

        if (t_B == trans_type::no_trans)
            N   = B.cols();
        else
            N   = B.rows();

        if (t_A == trans_type::no_trans)
        {
            M = A.rows();
            K = A.cols();
        }
        else
        {
            K = A.rows();
            M = A.cols();
        };

        if (M == 0 || K == 0 || N == 0)
        {
            using VR            = typename md::real_type<V>::type;
            using sparse_matrix = Matrix<VR,struct_sparse>;

            sparse_matrix out(ti::ti_empty(), M, N);
            ret = matcl::Matrix(out,false);
            return;
        };

        Integer A_ld = raw::get_ld(A,0);
        Integer A_ud = raw::get_ud(A,0);
        Integer B_ld = raw::get_ld(B,0);
        Integer B_ud = raw::get_ud(B,0);

        if (A_ld == 0 && A_ud == 0)
            return eval_dense_diag<V>::eval1(ret,A,B,t_A,t_B);

        if (B_ld == 0 && B_ud == 0)
            return eval_dense_diag<V>::eval2(ret,A,B,t_A,t_B);

        DM C(ti::ti_empty(),M, N);	

        Integer fr      = 0;
        Integer rows    = M;
        Integer fc      = 0;
        Integer cols    = N;
        V alpha         = V(1.0);
        V beta          = V(0.0);

        eval_impl(C, A, B, t_A, t_B, A_ld, A_ud, B_ld, B_ud, fr, rows, fc, cols, alpha, beta);

        bool is_sq_A    = (A.rows() == A.cols());
        bool is_sq_B    = (B.rows() == B.cols());
        bool is_sq_C    = (C.rows() == C.cols());
        C.set_struct(predefined_struct_ext::mult_struct(A.get_struct(),B.get_struct(), t_A, t_B,
                                is_real_matrix(A),is_real_matrix(B), is_sq_A, is_sq_B, is_sq_C));

        ret = matcl::Matrix(C,true);
    };

    static void eval_impl(DM& C, const DM& A, const DM& B, trans_type t_A, trans_type t_B,
                          Integer A_ld, Integer A_ud, Integer B_ld, Integer B_ud, 
                          Integer fr, Integer rows, Integer fc, Integer cols, const V& alpha, const V& beta)
    {
        Integer M   = rows;
        Integer N   = cols;
        Integer K   = t_A == trans_type::no_trans ? A.cols() : A.rows();

        if ( M == 1  &&  N == 1)
            return eval_dot(C, A, B, t_A, t_B, fr, rows, fc, cols, alpha, beta);
        
        if (K == 1)
        {
            if (beta == V(1.0))
                return eval_geru(C, A, B, t_A, t_B, fr, rows, fc, cols, alpha);
            else
            {
                prepare_gemm_C<V>::eval(beta, C, fr,rows,fc,cols);
                return eval_geru(C, A, B, t_A, t_B, fr, rows, fc, cols, alpha);
            };
        }
        
        bool is_A_sq    = A.rows() == A.cols();
        bool is_B_sq    = B.rows() == B.cols();

        if (A.size() >= B.size())
        {
            if (A_ud == 0)
            {
                if (M == K) 
                    return eval_tril1(C, A, B, t_A, t_B, fr, rows, fc, cols, alpha, beta);
            }
            else if (A_ld == 0)
            {
                if (M == K) 
                    return eval_triu1(C, A, B, t_A, t_B, fr, rows, fc, cols, alpha, beta);
            }
            else if (A.get_struct().is_symmetric(is_A_sq, is_real_matrix(A)))
            {
                return eval_sym1(C,A,B, t_A, t_B, fr, rows, fc, cols, alpha, beta);
            }
            else if (A.get_struct().is_hermitian(is_A_sq, is_real_matrix(A)))
            {
                return eval_her1(C,A,B, t_A, t_B, fr, rows, fc, cols, alpha, beta);
            }

            if (B_ud == 0)
            {
                if (K == N) 
                    return eval_tril2(C,A,B,t_A, t_B, fr, rows, fc, cols, alpha, beta);
            }
            else if (B_ld == 0)
            {
                if (K == N) 
                    return eval_triu2(C,A,B,t_A, t_B, fr, rows, fc, cols, alpha, beta);
            }
            else if (B.get_struct().is_symmetric(is_B_sq, is_real_matrix(B)))
                return eval_sym2(C,A,B,t_A, t_B, fr, rows, fc, cols, alpha, beta);
            else if (B.get_struct().is_hermitian(is_B_sq, is_real_matrix(B)))
                return eval_her2(C,A,B,t_A, t_B, fr, rows, fc, cols, alpha, beta);
        }
        else
        {
            if (B_ud == 0)
            {
                if (K == N) 
                    return eval_tril2(C,A,B,t_A, t_B, fr, rows, fc, cols, alpha, beta);
            }
            else if (B_ld == 0)
            {
                if (K == N) 
                    return eval_triu2(C,A,B,t_A, t_B, fr, rows, fc, cols, alpha, beta); 
            }
            else if (B.get_struct().is_symmetric(is_B_sq, is_real_matrix(B)))
                return eval_sym2(C,A,B,t_A, t_B, fr, rows, fc, cols, alpha, beta);
            else if (B.get_struct().is_hermitian(is_B_sq, is_real_matrix(B)))
                return eval_her2(C,A,B,t_A, t_B, fr, rows, fc, cols, alpha, beta);

            if (A_ud == 0)
            {
                if (M == K) 
                    return eval_tril1(C,A,B,t_A, t_B, fr, rows, fc, cols, alpha, beta);
            }
            else if (A_ld == 0)
            {
                if (M == K) 
                    return eval_triu1(C,A,B,t_A, t_B, fr, rows, fc, cols, alpha, beta);
            }
            else if (A.get_struct().is_symmetric(is_A_sq, is_real_matrix(A)))
                return eval_sym1(C,A,B, t_A, t_B, fr, rows, fc, cols, alpha, beta);
            else if (A.get_struct().is_hermitian(is_A_sq, is_real_matrix(A)))
                return eval_her1(C,A,B, t_A, t_B, fr, rows, fc, cols, alpha, beta);
        };

        return eval_general(C,A,B, t_A, t_B, fr, rows, fc, cols, alpha, beta);
    };

    static void eval_dot(DM& C, const DM& A, const DM& B, trans_type t_A, trans_type t_B,
                             Integer fr, Integer rows, Integer fc, Integer cols, const V& alpha, const V& beta)
    {
        (void)rows;
        (void)cols;
        using diag_type = struct_flag::diag_type;

        diag_type A_ud  = (t_A == trans_type::no_trans)? A.get_struct().get_udiags()
                                                        : A.get_struct().get_ldiags();
        diag_type B_ld  = (t_B == trans_type::no_trans)? B.get_struct().get_ldiags()
                                                        : B.get_struct().get_udiags();

        Integer K       = t_A == trans_type::no_trans ? A.cols() : A.rows();

        if (A_ud == struct_flag::zero)
            K           = std::min(K,1);
        else if (A_ud != struct_flag::general)
            K           = std::min(K,2);

        if (B_ld == struct_flag::zero)
            K           = std::min(K,1);
        else if (B_ld != struct_flag::general)
            K           = std::min(K,2);

        Integer lda     = (t_A == trans_type::no_trans)? A.ld() : 1;
        Integer ldb     = (t_B == trans_type::no_trans)? 1 : B.ld();
         
        V dot;
        if (t_A != trans_type::conj_trans && t_B != trans_type::conj_trans
            || md::is_float_real_scalar<V>::value == true)
        {
            dot = V(lapack::dot(K, md::lap(A.ptr()), lda, md::lap(B.ptr()), ldb));
        }
        else if (t_A == trans_type::conj_trans && t_B == trans_type::conj_trans)
        {
            dot = V(lapack::dot(K, md::lap(A.ptr()), lda, md::lap(B.ptr()), ldb));
            dot = conj(dot);
        }
        else if (t_A == trans_type::conj_trans)
        {
            dot = V(lapack::dotc(K, md::lap(A.ptr()), lda, md::lap(B.ptr()), ldb));
        }
        else
        {
            dot = V(lapack::dotc(K, md::lap(B.ptr()), ldb, md::lap(A.ptr()), lda));
        };

        V val       = alpha * dot;
        V* ptr_C    = C.ptr() + fr + fc * C.ld();

        if (mrd::is_zero(beta))
            ptr_C[0]    = val;
        else
            ptr_C[0]    = val + beta * ptr_C[0];

        return;
    };

    static void eval_geru(DM& C, const DM& A, const DM& B, trans_type t_A, trans_type t_B,
                            Integer fr, Integer rows, Integer fc, Integer cols, const V& alpha)
    {
        // K = 1
        using diag_type = struct_flag::diag_type;

        diag_type A_ld  = (t_A == trans_type::no_trans)? A.get_struct().get_ldiags()
                                                        : A.get_struct().get_udiags();
        diag_type B_ud  = (t_B == trans_type::no_trans)? B.get_struct().get_udiags()
                                                        : B.get_struct().get_ldiags();

        Integer M1 = rows;
        Integer N1 = cols;

        if (A_ld == struct_flag::zero)
            M1          = std::min(M1,1);
        else if (A_ld != struct_flag::general)
            M1          = std::min(M1,2);

        if (B_ud == struct_flag::zero)
            N1          = std::min(N1,1);
        else if (B_ud != struct_flag::general)
            N1          = std::min(N1,2);

        Integer lda     = (t_A == trans_type::no_trans)? 1: A.ld();
        Integer ldb     = (t_B == trans_type::no_trans)? B.ld() : 1;

        Integer C_ld    = C.ld();
        V* C_ptr        = C.ptr() + fr + fc * C_ld;

        if (t_A != trans_type::conj_trans && t_B != trans_type::conj_trans
            || md::is_float_real_scalar<V>::value == true)
        {
            lapack::geru(M1, N1, *md::lap(&alpha), md::lap(A.ptr()),lda, md::lap(B.ptr()),ldb,
                        md::lap(C_ptr), C_ld);
        }
        else if (t_A != trans_type::conj_trans && t_B == trans_type::conj_trans)
        {
            lapack::gerc(M1, N1, *md::lap(&alpha), md::lap(A.ptr()),lda, md::lap(B.ptr()),ldb,
                        md::lap(C_ptr), C_ld);
        }
        else
        {
            //taking conj of gerc/geru is probably slower; not checked
            V scal1 = V(1.0);
            lapack::gemm(trans_code<ISR>(t_A), trans_code<ISR>(t_B), M1, N1, 1, *md::lap(&alpha), md::lap(A.ptr()),
                        A.ld(), md::lap(B.ptr()), B.ld(), *md::lap(&scal1), md::lap(C_ptr), C.ld());
        };
    }

    static void eval_general(DM& C, const DM& A, const DM& B, trans_type t_A, trans_type t_B,
                             Integer fr, Integer rows, Integer fc, Integer cols, const V& alpha, const V& beta)
    {
        static const bool is_real = md::is_float_real_scalar<V>::value;
        
        Integer C_ld    = C.ld();
        V* C_ptr        = C.ptr() + fr + fc * C_ld;
        Integer M       = rows;
        Integer N       = cols;

        if (N == 1)
        {
            if (t_B == trans_type::no_trans)
            {
                lapack::gemv(trans_code<ISR>(t_A), A.rows(), A.cols(), *md::lap(&alpha), md::lap(A.ptr()), A.ld(),
                            md::lap(B.ptr()), 1, *md::lap(&beta), md::lap(C_ptr), 1);
                return;
            }
            else if (t_B == trans_type::trans || t_B == trans_type::conj_trans && is_real == true)
            {
                lapack::gemv(trans_code<ISR>(t_A), A.rows(), A.cols(), *md::lap(&alpha), md::lap(A.ptr()), A.ld(),
                            md::lap(B.ptr()), B.ld(), *md::lap(&beta), md::lap(C_ptr), 1);
                return;
            };
        };
        
        if (M == 1)
        {
            if (t_B == trans_type::no_trans || t_B == trans_type::trans
                || (t_B == trans_type::conj_trans) && is_real )
            {
                const char* TB = (t_B == trans_type::no_trans)? "T" : "N";

                if (t_A == trans_type::no_trans)
                {
                    lapack::gemv(TB, B.rows(), B.cols(), *md::lap(&alpha), md::lap(B.ptr()), B.ld(),
                                md::lap(A.ptr()), A.ld(), *md::lap(&beta), md::lap(C_ptr), C_ld);
                    return;
                }
                else if (t_A == trans_type::trans || t_A == trans_type::conj_trans && is_real == true)
                {
                    lapack::gemv(TB, B.rows(), B.cols(), *md::lap(&alpha), md::lap(B.ptr()), B.ld(),
                                md::lap(A.ptr()), 1, *md::lap(&beta), md::lap(C_ptr), C_ld);
                    return;
                };
            };
        };

        Integer K = t_A == trans_type::no_trans ? A.cols() : A.rows();
        lapack::gemm(trans_code<ISR>(t_A), trans_code<ISR>(t_B), M, N, K, *md::lap(&alpha), md::lap(A.ptr()),
                        A.ld(), md::lap(B.ptr()), B.ld(), *md::lap(&beta), md::lap(C_ptr), C_ld);
        return;
    };

    static void eval_triu1(DM& C, const DM& A, const DM& B, trans_type t_A, trans_type t_B,
                           Integer fr, Integer rows, Integer fc, Integer cols, const V& alpha, const V& beta)
    {
        Integer M   = rows;
        Integer N   = cols;

        if (A.rows() != A.cols())
            return eval_general(C, A, B, t_A, t_B, fr, rows, fc, cols, alpha, beta);

        if (beta != V(0.0))
           return eval_general(C, A, B, t_A, t_B, fr, rows, fc, cols, alpha, beta);

        Integer C_ld    = C.ld();
        V* C_ptr        = C.ptr() + fr + fc * C.ld();

        const V* B_ptr  = B.ptr();
        Integer B_ld    = B.ld();

        bool is_vector  = (N == 1);

        if (t_B == trans_type::no_trans)
        {
            if (is_vector == false)
            {
                level1::copy_mat<true, V, V, 0,0,0>::eval(C_ptr, C_ld, B_ptr, B_ld, M, N);
            }
            else
            {
                //vector version requires additional scaling
                level1::axpby_test_mat<true, V, V, V, V, 0, 0, 0>
                    ::eval(C_ptr, C_ld, B_ptr, B_ld, M, N, alpha, V(0.0));
            }            
        }
        else
        {
            if (t_B == trans_type::trans)
                manip_trans_helper_mat<V>::eval_trans(C_ptr, C_ld, B);
            else
                manip_trans_helper_mat<V>::eval_ctrans(C_ptr, C_ld, B);

            if (is_vector == true)
                level1::ay_test_mat<true, V, V, 0, 0, 0>::eval(C_ptr, C_ld, M, N, alpha);
        };


        if (is_vector == true)
        {
            lapack::trmv("U", trans_code<ISR>(t_A), "N", M, md::lap(A.ptr()), A.ld(), md::lap(C_ptr), 1);
        }
        else
        {
            lapack::trmm("L", "U", trans_code<ISR>(t_A), "N", M, N, *md::lap(&alpha), md::lap(A.ptr()), A.ld(), 
                            md::lap(C_ptr), C_ld);
        }

        return;
    };
    static void eval_triu2(DM& C, const DM& A, const DM& B, trans_type t_A, trans_type t_B,
                           Integer fr, Integer rows, Integer fc, Integer cols, const V& alpha, const V& beta)
    {
        Integer M   = rows;
        Integer N   = cols;

        if (B.rows() != B.cols())
            return eval_general(C,A,B,t_A, t_B, fr, rows, fc, cols, alpha, beta);

        if (beta != V(0.0))
           return eval_general(C, A, B, t_A, t_B, fr, rows, fc, cols, alpha, beta);

        Integer C_ld    = C.ld();
        V* C_ptr        = C.ptr() + fr + fc * C.ld();

        const V* A_ptr  = A.ptr();
        Integer A_ld    = A.ld();

        bool is_mat     = (M > 1 || (t_B == trans_type::conj_trans && ISR == false));

        if (t_A == trans_type::no_trans)
        {
            if (is_mat == true)
            {
                level1::copy_mat<true, V, V, 0,0,0>::eval(C_ptr, C_ld, A_ptr, A_ld, M, N);
            }
            else
            {
                //vector version requires additional scaling
                level1::axpby_test_mat<true, V, V, V, V, 0, 0, 0>
                    ::eval(C_ptr, C_ld, A_ptr, A_ld, M, N, alpha, V(0.0));
            }   
        }
        else
        {
            if (t_A == trans_type::trans)
                manip_trans_helper_mat<V>::eval_trans(C_ptr, C_ld, A);
            else
                manip_trans_helper_mat<V>::eval_ctrans(C_ptr, C_ld, A);

            //vector version requires additional scaling
            if (is_mat == false)
                level1::ay_test_mat<true, V, V, 0, 0, 0>::eval(C_ptr, C_ld, M, N, alpha);
        };

        if (is_mat == true)
        {
            lapack::trmm("R", "U", trans_code<ISR>(t_B), "N", M, N, *md::lap(&alpha), md::lap(B.ptr()), B.ld(), 
                        md::lap(C_ptr), C_ld);
        }
        else
        {
            // t_B == trans_type::no_trans || t_B == trans_type::trans || ISR
            //x * B = (B' * x')'

            const char* TB = (t_B == trans_type::no_trans)? "T" : "N";
            lapack::trmv("U", TB, "N", N, md::lap(B.ptr()), B.ld(), md::lap(C_ptr), C_ld);
        }

        return;
    };
    static void eval_tril1(DM& C, const DM& A, const DM& B, trans_type t_A, trans_type t_B,
                           Integer fr, Integer rows, Integer fc, Integer cols, const V& alpha, const V& beta)
    {
        Integer M   = rows;
        Integer N   = cols;

        if (A.rows() != A.cols())
            return eval_general(C, A, B, t_A, t_B, fr, rows, fc, cols, alpha, beta);
        
        if (beta != V(0.0))
           return eval_general(C, A, B, t_A, t_B, fr, rows, fc, cols, alpha, beta);

        Integer C_ld    = C.ld();
        V* C_ptr        = C.ptr() + fr + fc * C.ld();

        const V* B_ptr  = B.ptr();
        Integer B_ld    = B.ld();

        bool is_vector  = (N == 1);

        if (t_B == trans_type::no_trans)
        {
            if (is_vector == false)
            {
                level1::copy_mat<true, V, V, 0,0,0>::eval(C_ptr, C_ld, B_ptr, B_ld, M, N);
            }
            else
            {
                //vector version requires additional scaling
                level1::axpby_test_mat<true, V, V, V, V, 0, 0, 0>
                    ::eval(C_ptr, C_ld, B_ptr, B_ld, M, N, alpha, V(0.0));
            };
        }
        else
        {
            if (t_B == trans_type::trans)
                manip_trans_helper_mat<V>::eval_trans(C_ptr, C_ld, B);
            else
                manip_trans_helper_mat<V>::eval_ctrans(C_ptr, C_ld, B);

            //vector version requires additional scaling
            if (is_vector == true)
                level1::ay_test_mat<true, V, V, 0, 0, 0>::eval(C_ptr, C_ld, M, N, alpha);
        };

        if (is_vector == true)
        {
            lapack::trmv("L", trans_code<ISR>(t_A), "N", M, md::lap(A.ptr()), A.ld(), md::lap(C_ptr), 1);
        }
        else
        {
            lapack::trmm("L", "L", trans_code<ISR>(t_A), "N", M, N, *md::lap(&alpha), md::lap(A.ptr()), A.ld(), 
                            md::lap(C_ptr), C_ld);
        }
    };

    static void eval_tril2(DM& C, const DM& A, const DM& B, trans_type t_A, trans_type t_B,
                           Integer fr, Integer rows, Integer fc, Integer cols, const V& alpha, const V& beta)
    {
        Integer M   = rows;
        Integer N   = cols;

        if (B.rows() != B.cols())
            return eval_general(C,A,B,t_A, t_B, fr, rows, fc, cols, alpha, beta);

        if (beta != V(0.0))
           return eval_general(C, A, B, t_A, t_B, fr, rows, fc, cols, alpha, beta);

        Integer C_ld    = C.ld();
        V* C_ptr        = C.ptr() + fr + fc * C.ld();

        const V* A_ptr  = A.ptr();
        Integer A_ld    = A.ld();

        bool is_mat     = (M > 1 || (t_B == trans_type::conj_trans && ISR == false));

        if (t_A == trans_type::no_trans)
        {
            if (is_mat == true)
            {
                level1::copy_mat<true, V, V, 0, 0, 0>::eval(C_ptr, C_ld, A_ptr, A_ld, M, N);
            }
            else
            {
                //vector version requires additional scaling
                level1::axpby_test_mat<true, V, V, V, V, 0, 0, 0>::eval(C_ptr, C_ld, A_ptr, A_ld, M, N, alpha, V(0.0));
            }
        }
        else
        {
            if (t_A == trans_type::trans)
                manip_trans_helper_mat<V>::eval_trans(C_ptr, C_ld, A);
            else
                manip_trans_helper_mat<V>::eval_ctrans(C_ptr, C_ld, A);

            //vector version requires additional scaling
            if (is_mat == false)
                level1::ay_test_mat<true, V, V, 0, 0, 0>::eval(C_ptr, C_ld, M, N, alpha);
        }

        if (is_mat == true)
        {
            lapack::trmm("R", "L", trans_code<ISR>(t_B), "N", M, N, *md::lap(&alpha), md::lap(B.ptr()), B.ld(), 
                        md::lap(C_ptr), C_ld);
        }
        else
        {
            //M = 1
            // t_B == trans_type::no_trans || t_B == trans_type::trans || ISR
            //x * B = (B' * x')'

            const char* TB = (t_B == trans_type::no_trans)? "T" : "N";
            lapack::trmv("L", TB, "N", N, md::lap(B.ptr()), B.ld(), md::lap(C_ptr), C_ld);
        }
    };

    static void eval_sym1(DM& C, const DM& A, const DM& B, trans_type t_A, trans_type t_B,
                          Integer fr, Integer rows, Integer fc, Integer cols, const V& alpha, const V& beta)
    {
        bool V_is_real = md::is_float_real_scalar<V>::value;

        if (t_A == trans_type::conj_trans && V_is_real == false
            ||t_B == trans_type::conj_trans && V_is_real == false)
        {
            return eval_general(C, A,B,t_A, t_B, fr, rows, fc, cols, alpha, beta);
        };

        Integer M   = rows;
        Integer N   = cols;

        if (N != 1 && t_B != trans_type::no_trans)
            return eval_general(C, A,B,t_A, t_B, fr, rows, fc, cols, alpha, beta);

        Integer C_ld    = C.ld();
        V* C_ptr        = C.ptr() + fr + fc * C.ld();

        if (N == 1)
        {
            if (t_B == trans_type::no_trans)
            {
                lapack::symv("U", M, *md::lap(&alpha), md::lap(A.ptr()), A.ld(), md::lap(B.ptr()), 
                    1, *md::lap(&beta), md::lap(C_ptr), 1);
            }
            else if (t_B == trans_type::trans || t_B == trans_type::conj_trans && V_is_real == true)
            {
                lapack::symv("U", M, *md::lap(&alpha), md::lap(A.ptr()), A.ld(), md::lap(B.ptr()), 
                    B.ld(), *md::lap(&beta), md::lap(C_ptr), 1);
            }
            else
            {
                //this case should already be removed
            };
        }
        else
        {
            //t_B == trans_type::no_trans
            lapack::symm("L", "U", M, N, *md::lap(&alpha), md::lap(A.ptr()), A.ld(), md::lap(B.ptr()), 
                    B.ld(), *md::lap(&beta), md::lap(C_ptr), C_ld);
        };
    };
    static void eval_sym2(DM& C, const DM& A, const DM& B, trans_type t_A, trans_type t_B,
                          Integer fr, Integer rows, Integer fc, Integer cols, const V& alpha, const V& beta)
    {
        bool V_is_real = md::is_float_real_scalar<V>::value;

        if (t_B == trans_type::conj_trans && V_is_real == false
            ||t_A == trans_type::conj_trans && V_is_real == false)
        {
            return eval_general(C, A,B,t_A, t_B, fr, rows, fc, cols, alpha, beta);
        };

        Integer M   = rows;
        Integer N   = cols;

        if (M != 1 && t_A != trans_type::no_trans)
            return eval_general(C, A,B,t_A, t_B, fr, rows, fc, cols, alpha, beta);

        Integer C_ld    = C.ld();
        V* C_ptr        = C.ptr() + fr + fc * C.ld();

        if (M == 1)
        {
            if (t_A == trans_type::no_trans)
            {
                lapack::symv("U", N, *md::lap(&alpha), md::lap(B.ptr()), B.ld(), md::lap(A.ptr()), 
                    A.ld(), *md::lap(&beta), md::lap(C_ptr), C_ld);
            }
            else if (t_A == trans_type::trans || t_A == trans_type::conj_trans && V_is_real == true)
            {
                lapack::symv("U", N, *md::lap(&alpha), md::lap(B.ptr()), B.ld(), md::lap(A.ptr()), 
                    1, *md::lap(&beta), md::lap(C_ptr), C_ld);
            }
            else
            {
                //this case should already be removed
            };
        }
        else
        {
            //t_A == trans_type::no_trans
            lapack::symm("R", "U", M, N, *md::lap(&alpha), md::lap(B.ptr()), B.ld(), md::lap(A.ptr()), 
                    A.ld(), *md::lap(&beta), md::lap(C_ptr), C_ld);
        };
    };

    static void eval_her1(DM& C, const DM& A, const DM& B, trans_type t_A, trans_type t_B,
                          Integer fr, Integer rows, Integer fc, Integer cols, const V& alpha, const V& beta)
    {
        bool V_is_real = md::is_float_real_scalar<V>::value;

        //required: op(A) = A
        if (t_A == trans_type::trans && V_is_real == false)
            return eval_general(C, A,B,t_A, t_B, fr, rows, fc, cols, alpha, beta);

        if (t_B == trans_type::conj_trans && V_is_real == false)
            return eval_general(C, A,B,t_A, t_B, fr, rows, fc, cols, alpha, beta);

        Integer M   = rows;
        Integer N   = cols;

        if (N != 1 && t_B != trans_type::no_trans)
            return eval_general(C, A,B,t_A, t_B, fr, rows, fc, cols, alpha, beta);

        Integer C_ld    = C.ld();
        V* C_ptr        = C.ptr() + fr + fc * C.ld();

        if (N == 1)
        {
            if (t_B == trans_type::no_trans)
            {
                // y := alpha*A*x + beta*y,
                lapack::hemv("U", M, *md::lap(&alpha), md::lap(A.ptr()), A.ld(), md::lap(B.ptr()), 1, 
                        *md::lap(&beta), md::lap(C_ptr), 1);
            }
            else if (t_B == trans_type::trans || t_B == trans_type::conj_trans && V_is_real == true)
            {
                lapack::hemv("U", M, *md::lap(&alpha), md::lap(A.ptr()), A.ld(), md::lap(B.ptr()), B.ld(), 
                        *md::lap(&beta), md::lap(C_ptr), 1);
            }
            else
            {
                //this case should already be removed
            };
        }
        else
        {
            //t_B == trans_type::no_trans
            lapack::hemm("L", "U", M, N, *md::lap(&alpha), md::lap(A.ptr()), A.ld(), md::lap(B.ptr()), B.ld(), 
                        *md::lap(&beta), md::lap(C_ptr), C_ld);
        }
    };

    static void eval_her2(DM& C, const DM& A, const DM& B, trans_type t_A, trans_type t_B,
                          Integer fr, Integer rows, Integer fc, Integer cols, const V& alpha, const V& beta)
    {
        bool V_is_real = md::is_float_real_scalar<V>::value;

        //required: op(B) = B
        if (t_B == trans_type::trans && V_is_real == false)
            return eval_general(C, A,B,t_A, t_B, fr, rows, fc, cols, alpha, beta);

        if (t_A != trans_type::no_trans)
            return eval_general(C, A,B,t_A, t_B, fr, rows, fc, cols, alpha, beta);

        Integer M       = rows;
        Integer N       = cols;

        Integer C_ld    = C.ld();
        V* C_ptr        = C.ptr() + fr + fc * C.ld();

        // t_A == trans_type::no_trans
        // C := alpha*B*A + beta*C, A is a hermitian matrix

        lapack::hemm("R", "U", M, N, *md::lap(&alpha), md::lap(B.ptr()), B.ld(), md::lap(A.ptr()), A.ld(), 
                    *md::lap(&beta), md::lap(C_ptr), C_ld);
    };
};

template<> struct eval_dense_dense<Integer>
{
    using VT    = Integer;
    using Mat   = raw::Matrix<VT,struct_dense>;

    static void eval(matcl::Matrix& ret, const Mat& A, const Mat& B, 
                     trans_type t_A, trans_type t_B)
    {
        return eval_dense_dense_generic<VT>::eval(ret,A,B,t_A,t_B);
    };

    static void eval_gemm(Mat& C, const VT& alpha, const Mat& A, const Mat& B, 
                     trans_type t_A, trans_type t_B, const VT& beta, 
                     Integer fr, Integer rows, Integer fc, Integer cols)
    {
        return eval_dense_dense_generic<VT>::eval_gemm(C, alpha, A, B, t_A, t_B, beta, fr, rows, fc, cols);
    };
}; 
template<> struct eval_dense_dense<Object>
{
    using VT    = Object;
    using Mat   = raw::Matrix<VT,struct_dense>;

    static void eval(matcl::Matrix& ret, const Mat& A, const Mat& B, 
                     trans_type t_A, trans_type t_B)
    {
        return eval_dense_dense_generic<VT>::eval(ret,A,B,t_A,t_B);
    };

    static void eval_gemm(Mat& C, const VT& alpha, const Mat& A, const Mat& B, 
                     trans_type t_A, trans_type t_B, const VT& beta, 
                     Integer fr, Integer rows, Integer fc, Integer cols)
    {
        return eval_dense_dense_generic<VT>::eval_gemm(C, alpha, A, B, t_A, t_B, beta, fr, rows, fc, cols);
    };
}; 
template<> struct eval_dense_dense<Real>
{
    using VT    = Real;
    using Mat   = raw::Matrix<VT,struct_dense>;

    static void eval(matcl::Matrix& ret, const Mat& A, const Mat& B, 
                     trans_type t_A, trans_type t_B)
    {
        return eval_dense_dense_lapack<VT>::eval(ret,A,B,t_A,t_B);
    };

    static void eval_gemm(Mat& C, const VT& alpha, const Mat& A, const Mat& B, 
                     trans_type t_A, trans_type t_B, const VT& beta, 
                     Integer fr, Integer rows, Integer fc, Integer cols)
    {
        return eval_dense_dense_lapack<VT>::eval_gemm(C, alpha, A, B, t_A, t_B, beta, fr, rows, fc, cols);
    };
};
template<> struct eval_dense_dense<Float>
{
    using VT    = Float;
    using Mat   = raw::Matrix<VT,struct_dense>;

    static void eval(matcl::Matrix& ret, const Mat& A, const Mat& B, 
                     trans_type t_A, trans_type t_B)
    {
        return eval_dense_dense_lapack<VT>::eval(ret,A,B,t_A,t_B);
    };

    static void eval_gemm(Mat& C, const VT& alpha, const Mat& A, const Mat& B, 
                     trans_type t_A, trans_type t_B, const VT& beta, 
                     Integer fr, Integer rows, Integer fc, Integer cols)
    {
        return eval_dense_dense_lapack<VT>::eval_gemm(C, alpha, A, B, t_A, t_B, beta, fr, rows, fc, cols);
    };

};
template<> struct eval_dense_dense<Complex>
{
    using VT    = Complex;
    using Mat   = raw::Matrix<VT,struct_dense>;

    static void eval(matcl::Matrix& ret, const Mat& A, const Mat& B, 
                     trans_type t_A, trans_type t_B)
    {
        return eval_dense_dense_lapack<VT>::eval(ret,A,B,t_A,t_B);
    };

    static void eval_gemm(Mat& C, const VT& alpha, const Mat& A, const Mat& B, 
                     trans_type t_A, trans_type t_B, const VT& beta, 
                     Integer fr, Integer rows, Integer fc, Integer cols)
    {
        return eval_dense_dense_lapack<VT>::eval_gemm(C, alpha, A, B, t_A, t_B, beta, fr, rows, fc, cols);
    };

};
template<> struct eval_dense_dense<Float_complex>
{
    using VT    = Float_complex;
    using Mat   = raw::Matrix<VT,struct_dense>;

    static void eval(matcl::Matrix& ret, const Mat& A, const Mat& B, 
                     trans_type t_A, trans_type t_B)
    {
        return eval_dense_dense_lapack<VT>::eval(ret,A,B,t_A,t_B);
    };

    static void eval_gemm(Mat& C, const VT& alpha, const Mat& A, const Mat& B, 
                     trans_type t_A, trans_type t_B, const VT& beta, 
                     Integer fr, Integer rows, Integer fc, Integer cols)
    {
        return eval_dense_dense_lapack<VT>::eval_gemm(C, alpha, A, B, t_A, t_B, beta, fr, rows, fc, cols);
    };

};

template<class Val_ret, class M1, class M2>
struct eval_mult<Val_ret,M1,M2,struct_dense,struct_dense> 
{ 
    using Mat_ret   = raw::Matrix<Val_ret,struct_dense>;

    static void eval(matcl::Matrix& ret, const M1& A, const M2& B, trans_type t_A, trans_type t_B)
    {
        Mat_ret AC  = converter<Mat_ret,M1>::eval(A);
        Mat_ret BC  = converter<Mat_ret,M2>::eval(B);
        return eval_dense_dense<Val_ret>::eval(ret,AC,BC,t_A, t_B);
    };

    static void eval_gemm(Mat_ret& C, const Val_ret& alpha, const M1& A, const M2& B, 
                     trans_type t_A, trans_type t_B, const Val_ret& beta, 
                     Integer fr, Integer rows, Integer fc, Integer cols)
    {
        Mat_ret AC  = converter<Mat_ret,M1>::eval(A);
        Mat_ret BC  = converter<Mat_ret,M2>::eval(B);
        return eval_dense_dense<Val_ret>::eval_gemm(C, alpha, AC, BC, t_A, t_B, beta, fr, rows, fc, cols);
    };
};

template<class Val_C, class M1, class M2>
struct eval_mult_abs<Val_C, M1, M2, struct_dense>
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

        using V1        = typename M1::value_type;
        using V2        = M2;
        using VR        = typename md::real_type<Val_C>::type;
        using V1R       = typename md::real_type<V1>::type;
        using Mat_ret   = raw::Matrix<VR, struct_dense>;

        Integer M       = C.rows();
        Integer N       = C.cols();

        Mat_ret ret(C.get_type(), M, N);

        const V1* ptr_A     = A.ptr();
        const Val_C* ptr_C  = C.ptr();
        VR* ptr_ret         = ret.ptr();

        Integer A_ld        = A.ld();
        Integer C_ld        = C.ld();
        Integer ret_ld      = ret.ld();

        if (t_A == trans_type::no_trans)
        {
            //C = scal * abs(A) + abs(C)
            for (Integer j = 0; j < N; ++j) 
            {
                for (Integer i = 0; i < M; ++i)
                {
                    ptr_ret[i]  = abs_helper<Val_C>::eval(ptr_C[i]) 
                                + scal * abs_helper<V1>::eval(ptr_A[i]);
                }

                ptr_A   += A_ld;
                ptr_C   += C_ld;
                ptr_ret += ret_ld;
            };
        }
        else
        {
            //C = scal * abs(A') + abs(C)

            for (Integer j = 0; j < N; ++j) 
            {
                for (Integer i = 0; i < M; ++i)
                {
                    ptr_ret[i]  = abs_helper<Val_C>::eval(ptr_C[i]) 
                                + scal * abs_helper<V1>::eval(ptr_A[i*A_ld]);
                }

                ptr_A   += 1;
                ptr_C   += C_ld;
                ptr_ret += ret_ld;
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

        using V1        = typename M1::value_type;
        using V2        = typename M2::value_type;
        using VR        = typename md::real_type<Val_C>::type;
        using V1R       = typename md::real_type<V1>::type;
        using Mat_ret   = raw::Matrix<VR, struct_dense>;

        Integer M       = C.rows();
        Integer N       = C.cols();

        Mat_ret ret(C.get_type(), M, N);

        const V1* ptr_A     = A.ptr();
        const V2* ptr_X     = X.ptr();
        const Val_C* ptr_C  = C.ptr();
        VR* ptr_ret         = ret.ptr();

        Integer A_ld        = A.ld();
        Integer X_ld        = X.ld();
        Integer C_ld        = C.ld();
        Integer ret_ld      = ret.ld();

        const V1* ptr_As    = ptr_A;

        if (t_A == trans_type::no_trans)
        {
            Integer K   = A.cols();

            // C = A*B
            for (Integer j = 0; j < N; ++j) 
            {
                //make abs(C);
                for (Integer i = 0; i < M; ++i)
                    ptr_ret[i]  = abs_helper<Val_C>::eval(ptr_C[i]);

                ptr_A       = ptr_As;

                for (Integer l = 0; l < K; ++l) 
                {
                    const VR& temp  = abs_helper<V2>::eval(ptr_X[l]);

                    if (mrd::is_zero(temp) == false) 
                    {					
                        for (Integer k = 0; k < M; ++k) 
                        {
                            VR tmp      = abs_helper<V1>::eval(ptr_A[k]) * temp;
                            ptr_ret[k]  = ptr_ret[k] + tmp;
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
            Integer K       = A.rows();

            //C = A'*B
            for (Integer j = 0; j < N; ++j) 
            {			
                ptr_A = ptr_As;
                for (Integer l = 0; l < M; ++l) 
                {
                    VR dot      = abs_helper<Val_C>::eval(ptr_C[l]);

                    for (Integer k = 0; k < K; ++k) 
                    {
                        VR tmp  = mrd::abs_helper<V1>::eval(ptr_A[k]) * mrd::abs_helper<V2>::eval(ptr_X[k]);
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
    }

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
        using DM        = raw::Matrix<V1R, struct_dense>;
        using Mat_ret   = raw::Matrix<VR, struct_dense>;

        Mat_ret ret(C.get_type(), M, N);

        Integer A_ld    = A.ld();
        Integer C_ld    = C.ld();
        Integer ret_ld  = ret.ld();

        const V1* A_ptr     = A.ptr();
        const Val_C* C_ptr  = C.ptr();
        VR* ret_ptr         = ret.ptr();

        for (Integer j = 0; j < K1; ++j)
        {
            for (Integer i = 0; i < M; ++i)
                ret_ptr[i]  = abs_helper<Val_C>::eval(C_ptr[i]);

            ret_ptr[j]      = ret_ptr[j] + scal * abs_helper<V1>::eval(A_ptr[j]);

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
        using VR        = typename md::real_type<Val_C>::type;
        using DM        = raw::Matrix<V1R, struct_dense>;
        using Mat_ret   = raw::Matrix<VR, struct_dense>;

        Mat_ret ret(C.get_type(), M, N);

        DM diag(A.get_type(), K1, 1);

        Integer A_ld    = A.ld();
        Integer X_ld    = X.ld();
        Integer C_ld    = C.ld();
        Integer ret_ld  = ret.ld();

        const V1* A_ptr     = A.ptr();
        const V2* X_ptr     = X.ptr();
        V1R* D_ptr          = diag.ptr();
        const Val_C* C_ptr  = C.ptr();
        VR* ret_ptr         = ret.ptr();

        for(Integer i = 0; i < K1; ++i)
        {
            D_ptr[i]    = abs_helper<V1>::eval(A_ptr[0]);
            A_ptr       += A_ld + 1;
        };

        for (Integer j = 0; j < N; ++j)
        {
            for (Integer i = 0; i < K1; ++i)
            {
                ret_ptr[i]  = abs_helper<Val_C>::eval(C_ptr[i]) 
                            + D_ptr[i] * abs_helper<V2>::eval(X_ptr[i]);
            }

            for (Integer i = K1; i < M; ++i)
                ret_ptr[i]  = abs_helper<Val_C>::eval(C_ptr[i]);

            C_ptr       += C_ld;
            ret_ptr     += ret_ld;
            X_ptr       += X_ld;
        };

        out = matcl::Matrix(ret,false);
        return;
    };
};

}}}
