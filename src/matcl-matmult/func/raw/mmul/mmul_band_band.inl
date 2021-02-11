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
//                      BAND - BAND
//========================================================================
template<class Val_ret, class M1, class M2>
struct eval_mult<Val_ret,M1,M2,struct_banded,struct_banded> 
{     
    using Mat_D     = raw::Matrix<Val_ret,struct_dense>;

    template<bool Conj_A, bool Conj_B, class Alpha>
    static void eval_diag_diag_gemm(const Alpha& alpha, Mat_D& C, const M1& A, const M2& B, 
                        trans_type t_A, trans_type t_B, Integer C_fr, Integer rows, Integer C_fc, Integer cols)
    {
        (void)t_A;
        (void)t_B;
        (void)rows;
        (void)cols;

        Integer rc      = std::min(A.diag_length(0),B.diag_length(0));

        using V1        = typename M1::value_type;
        using V2        = typename M2::value_type;
        using VR        = Val_ret;
        using V12       = typename md::unify_types<V1,V2>::type;

        Integer A_ld    = A.ld();
        Integer B_ld    = B.ld();
        Integer C_ld    = C.ld();

        const V1* ptr_A = A.rep_ptr() + A.first_elem_diag(0);
        const V2* ptr_B = B.rep_ptr() + B.first_elem_diag(0);
        VR * ptr_C      = C.ptr() + C_fr + C_fc * C_ld;

        for (Integer j = 0; j < rc; ++j)
        {
            V12 ab      = make_conj<Conj_A>::eval(ptr_A[0]) * make_conj<Conj_B>::eval(ptr_B[0]);

            ptr_C[0]    = alpha.eval(ptr_C[0], ab);

            ptr_A       += A_ld;
            ptr_B       += B_ld;
            ptr_C       += C_ld + 1;
        };

        return;
    };

    static void eval_diag_diag(matcl::Matrix& ret, const M1& A, const M2& B, 
                               trans_type t_A, trans_type t_B)
    {
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

        Integer rc      = std::min(A.diag_length(0),B.diag_length(0));

        using ti_ret_type   = typename ti::get_ti_type<Val_ret>::type;
        using val_type_1    = typename M1::value_type;
        using val_type_2    = typename M2::value_type;

        ti_ret_type ret_ti  = ti::get_return_ti<ti_ret_type>(mdyf::op_mul::eval(),
                                                                  ti::get_ti(A), ti::get_ti(B));
        
        Val_ret Z           = md::default_value<Val_ret>(ret_ti);

        using Mat_B         = Matrix<Val_ret,struct_banded>;

        Mat_B C(ret_ti, M, N, 0, 0);

        const val_type_1* ptr_A = A.rep_ptr() + A.first_elem_diag(0);
        const val_type_2* ptr_B = B.rep_ptr() + B.first_elem_diag(0);
        Val_ret * ptr_C         = C.rep_ptr() + C.first_elem_diag(0);

        Integer A_ld    = A.ld();
        Integer B_ld    = B.ld();

        if (t_A == trans_type::conj_trans)
        {
            if (t_B == trans_type::conj_trans)
            {
                for (Integer j = 0; j < rc; ++j)
                {
                    ptr_C[j]    = conj(ptr_A[0]) * conj(ptr_B[0]);
                    ptr_A       += A_ld;
                    ptr_B       += B_ld;
                };
            }
            else
            {
                for (Integer j = 0; j < rc; ++j)
                {
                    ptr_C[j]    = conj(ptr_A[0]) * ptr_B[0];
                    ptr_A       += A_ld;
                    ptr_B       += B_ld;
                };
            };
        }
        else
        {
            if (t_B == trans_type::conj_trans)
            {
                for (Integer j = 0; j < rc; ++j)
                {
                    ptr_C[j]    = ptr_A[0] * conj(ptr_B[0]);
                    ptr_A       += A_ld;
                    ptr_B       += B_ld;
                };
            }
            else
            {
                for (Integer j = 0; j < rc; ++j)
                {
                    ptr_C[j]    = ptr_A[0] * ptr_B[0];
                    ptr_A       += A_ld;
                    ptr_B       += B_ld;
                };
            };
        };

        if (M > 0)
        {
            for (Integer j = rc; j < C.cols(); ++j)
                ptr_C[j]    = Z;
        };

        bool is_sq_A    = A.rows() == A.cols();
        bool is_sq_B    = B.rows() == B.cols();
        bool is_sq_C    = C.rows() == C.cols();
        C.set_struct(predefined_struct_ext::mult_struct(A.get_struct(),B.get_struct(),t_A,t_B,
                                is_real_matrix(A),is_real_matrix(B), is_sq_A, is_sq_B, is_sq_C));

        ret = matcl::Matrix(C,true);
        return;
    };

    static void eval_diag_1(matcl::Matrix& ret, const M1& A, const M2& B, trans_type t_A, trans_type t_B)
    {
        if (raw::is_diag(B) == true)
            return eval_diag_diag(ret, A, B, t_A, t_B);

        Integer M, K, N, fd_C, ld_C;

        if (t_B == trans_type::no_trans)
        {
            N       = B.cols();
            fd_C    = B.first_diag();
            ld_C    = B.last_diag();
        }
        else
        {
            N       = B.rows();
            fd_C    = -B.last_diag();
            ld_C    = -B.first_diag();
        }

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
        
        Integer K1          = std::min(M,K);

        using ti_ret_type   = typename ti::get_ti_type<Val_ret>::type;
        using val_type_1    = typename M1::value_type;
        using val_type_2    = typename M2::value_type;

        ti_ret_type ret_ti  = ti::get_return_ti<ti_ret_type>(mdyf::op_mul::eval(),
                                                                  ti::get_ti(A), ti::get_ti(B));
        
        Val_ret Z           = md::default_value<Val_ret>(ret_ti);
        using Mat_B         = Matrix<Val_ret,struct_banded>;

        Mat_B C(ret_ti, Z, M, N, fd_C, ld_C);

        const val_type_1* ptr_A = A.rep_ptr() + A.first_elem_diag(0);
        const val_type_2* ptr_B = B.rep_ptr();
        Val_ret * ptr_C         = C.rep_ptr();

        Integer A_ld        = A.ld();
        Integer B_ld        = B.ld();
        Integer C_ld        = C.ld();

        bool conj_A         = (t_A == trans_type::conj_trans);

        if (t_B == trans_type::no_trans)
        {
            for (Integer j = 0; j < N; ++j)
            {
                Integer fr      = B.first_row(j);
                Integer lr      = B.last_row(j);
                Integer lr2     = std::min(lr,K1 - 1);

                Integer pos_B   = B.first_elem_pos(j);
                Integer pos_C   = C.first_elem_pos(j);

                const val_type_1* ptr_A0 = ptr_A + fr*A_ld;

                if (conj_A == true)
                {
                    for (Integer i = fr; i <= lr2; ++i, ++pos_B, ++pos_C)
                    {
                        ptr_C[pos_C] = conj(ptr_A0[0]) * ptr_B[pos_B];
                        ptr_A0      += A_ld;
                    };
                }
                else
                {
                    for (Integer i = fr; i <= lr2; ++i, ++pos_B, ++pos_C)
                    {
                        ptr_C[pos_C] = ptr_A0[0] * ptr_B[pos_B];
                        ptr_A0      += A_ld;
                    };
                };

                ptr_B           += B_ld;
                ptr_C           += C_ld;
            };
        }
        else
        {
            bool conj_B         = (t_B == trans_type::conj_trans);

            Integer Bfd         = B.first_diag();
            Integer Bld         = B.last_diag();

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

                const val_type_1* ptr_A0 = ptr_A + r*A_ld;

                if (conj_A == true && conj_B == true)
                {
                    for (Integer i = 0; i < s; ++i, ++r)
                    {
                        ptr_C[pos_C]    = conj(ptr_A0[0]) * conj(ptr_B[pos_B]);
                        pos_B       += B_ld;
                        pos_C       += C_ld;
                        ptr_A0      += A_ld;
                    };
                }
                else if (conj_A == true && conj_B == false)
                {
                    for (Integer i = 0; i < s; ++i, ++r)
                    {
                        ptr_C[pos_C]    = conj(ptr_A0[0]) * ptr_B[pos_B];
                        pos_B       += B_ld;
                        pos_C       += C_ld;
                        ptr_A0      += A_ld;
                    };
                }
                else if (conj_A == false && conj_B == true)
                {
                    for (Integer i = 0; i < s; ++i, ++r)
                    {
                        ptr_C[pos_C]    = ptr_A0[0] * conj(ptr_B[pos_B]);
                        pos_B       += B_ld;
                        pos_C       += C_ld;
                        ptr_A0      += A_ld;
                    };
                }
                else
                {
                    for (Integer i = 0; i < s; ++i, ++r)
                    {
                        ptr_C[pos_C]    = ptr_A0[0] * ptr_B[pos_B];
                        pos_B       += B_ld;
                        pos_C       += C_ld;
                        ptr_A0      += A_ld;
                    };
                };
            };
        };

        bool is_sq_A    = A.rows() == A.cols();
        bool is_sq_B    = B.rows() == B.cols();
        bool is_sq_C    = C.rows() == C.cols();
        C.set_struct(predefined_struct_ext::mult_struct(A.get_struct(),B.get_struct(),t_A,t_B,
                                is_real_matrix(A),is_real_matrix(B), is_sq_A, is_sq_B, is_sq_C));

        ret = matcl::Matrix(C,true);
        return;
    };

    static void eval_diag_2(matcl::Matrix& ret, const M1& A, const M2& B, trans_type t_A, trans_type t_B)
    {
        if (raw::is_diag(A) == true)
            return eval_diag_diag(ret, A, B, t_A, t_B);

        Integer M, K, N, fd_C, ld_C;

        if (t_B == trans_type::no_trans)
            N       = B.cols();
        else
            N       = B.rows();

        if (t_A == trans_type::no_trans)
        {
            M       = A.rows();
            K       = A.cols();
            fd_C    = A.first_diag();
            ld_C    = A.last_diag();
        }
        else
        {
            K       = A.rows();
            M       = A.cols();
            fd_C    = -A.last_diag();
            ld_C    = -A.first_diag();
        };
        
        Integer K1  = std::min(N,K);

        using ti_ret_type   = typename ti::get_ti_type<Val_ret>::type;
        using val_type_1    = typename M1::value_type;
        using val_type_2    = typename M2::value_type;

        ti_ret_type ret_ti  = ti::get_return_ti<ti_ret_type>(mdyf::op_mul::eval(),
                                                                  ti::get_ti(A), ti::get_ti(B));
        
        Val_ret Z           = md::default_value<Val_ret>(ret_ti);

        using Mat_B         = Matrix<Val_ret,struct_banded>;
        Mat_B C(ret_ti, Z, M, N, fd_C, ld_C);

        const val_type_1* ptr_A = A.rep_ptr();
        const val_type_2* ptr_B = B.rep_ptr() + B.first_elem_diag(0);
        Val_ret * ptr_C         = C.rep_ptr();

        Integer A_ld    = A.ld();
        Integer B_ld    = B.ld();
        Integer C_ld    = C.ld();

        bool conj_B     = (t_B == trans_type::conj_trans);

        if (t_A == trans_type::no_trans)
        {
            for (Integer j = 0; j < K1; ++j)
            {
                const val_type_2& v = conj_B ? conj(ptr_B[0]) : ptr_B[0];

                if (mrd::is_zero(v))
                {
                    ptr_A   += A_ld;
                    ptr_B   += B_ld;
                    ptr_C   += C_ld;
                    continue;
                };

                Integer fr      = A.first_row(j);
                Integer lr      = A.last_row(j);
                Integer pos_A   = A.first_elem_pos(j);
                Integer pos_C   = C.first_elem_pos(j);

                for (Integer i = fr; i <= lr; ++i, ++pos_A, ++pos_C)
                    ptr_C[pos_C] = ptr_A[pos_A] * v;

                ptr_A   += A_ld;
                ptr_B   += B_ld;
                ptr_C   += C_ld;
            };
        }
        else
        {
            bool conj_A     = (t_A == trans_type::conj_trans);

            for (Integer d = A.first_diag(); d <= A.last_diag(); ++d)
            {
                if (C.has_diag(-d) == false)
                    continue;

                Integer pos_A = A.first_elem_diag(d);
                Integer pos_C = C.first_elem_diag(-d);
                Integer fc    = C.first_col_on_diag(-d);

                Integer s     = std::min(A.diag_length(d), C.diag_length(-d));

                const val_type_2* ptr_B0 = ptr_B + fc * B_ld;

                if (conj_A == true && conj_B == true)
                {
                    for (Integer i = 0; i < s; ++i, ++fc)
                    {
                        ptr_C[pos_C] = conj(ptr_A[pos_A]) * conj(ptr_B0[0]);

                        pos_A   += A_ld;
                        pos_C   += C_ld;
                        ptr_B0  += B_ld;
                    };
                }
                else if (conj_A == true && conj_B == false)
                {
                    for (Integer i = 0; i < s; ++i, ++fc)
                    {
                        ptr_C[pos_C] = conj(ptr_A[pos_A]) * ptr_B0[0];

                        pos_A   += A_ld;
                        pos_C   += C_ld;
                        ptr_B0  += B_ld;
                    };
                }
                else if (conj_A == false && conj_B == true)
                {
                    for (Integer i = 0; i < s; ++i, ++fc)
                    {
                        ptr_C[pos_C] = ptr_A[pos_A] * conj(ptr_B0[0]);

                        pos_A   += A_ld;
                        pos_C   += C_ld;
                        ptr_B0  += B_ld;
                    };
                }
                else
                {
                    for (Integer i = 0; i < s; ++i, ++fc)
                    {
                        ptr_C[pos_C] = ptr_A[pos_A] * ptr_B0[0];

                        pos_A   += A_ld;
                        pos_C   += C_ld;
                        ptr_B0  += B_ld;
                    };
                };
            };
        };

        bool is_sq_A    = A.rows() == A.cols();
        bool is_sq_B    = B.rows() == B.cols();
        bool is_sq_C    = C.rows() == C.cols();
        C.set_struct(predefined_struct_ext::mult_struct(A.get_struct(),B.get_struct(),t_A,t_B,
                                is_real_matrix(A),is_real_matrix(B), is_sq_A, is_sq_B, is_sq_C));

        ret = matcl::Matrix(C,true);
        return;
    };

    template<bool Conj_A, bool Conj_B, class Alpha>
    static void eval_diag_1_gemm(const Alpha& alpha, Mat_D& C, const M1& A, const M2& B, 
                        trans_type t_A, trans_type t_B, Integer C_fr, Integer rows, Integer C_fc, Integer cols)
    {
        if (raw::is_diag(B) == true)
            return eval_diag_diag_gemm<Conj_A, Conj_B>(alpha, C, A, B, t_A, t_B, C_fr, rows, C_fc, cols);

        Integer N           = cols;
        
        Integer K1          = std::min(A.rows(),A.cols());

        using V1            = typename M1::value_type;
        using V2            = typename M2::value_type;
        using V12           = typename md::unify_types<V1,V2>::type;
        using VR            = Val_ret;
        
        Integer A_ld        = A.ld();
        Integer B_ld        = B.ld();
        Integer C_ld        = C.ld();

        const V1* ptr_A     = A.rep_ptr() + A.first_elem_diag(0);
        const V2* ptr_B     = B.rep_ptr();
        Val_ret * ptr_C     = C.ptr() + C_fr + C_fc * C_ld;

        if (t_B == trans_type::no_trans)
        {
            for (Integer j = 0; j < N; ++j)
            {
                Integer fr      = B.first_row(j);
                Integer lr      = B.last_row(j);
                Integer lr2     = std::min(lr,K1 - 1);

                Integer pos_B   = B.first_elem_pos(j);

                const V1* ptr_Al    = ptr_A + fr*A_ld;

                for (Integer i = fr; i <= lr2; ++i, ++pos_B)
                {
                    V12 tmp     = make_conj<Conj_A>::eval(ptr_Al[0]) 
                                    * make_conj<Conj_B>::eval(ptr_B[pos_B]);
                    ptr_C[i]    = alpha.eval(ptr_C[i], tmp);
                    ptr_Al      += A_ld;
                };

                ptr_B           += B_ld;
                ptr_C           += C_ld;
            };
        }
        else
        {
            for (Integer j = 0; j < N; ++j)
            {
                Integer fc      = B.first_col(j);
                Integer lc      = B.last_col(j);
                Integer lc2     = std::min(lc,K1 - 1);

                Integer pos_B   = B.first_elem_pos_row(j);

                const V1* ptr_Al    = ptr_A + fc*A_ld;

                for (Integer i = fc; i <= lc2; ++i, pos_B += B_ld-1)
                {
                    V12 tmp     = make_conj<Conj_A>::eval(ptr_Al[0]) 
                                    * make_conj<Conj_B>::eval(ptr_B[pos_B]);
                    ptr_C[i]    = alpha.eval(ptr_C[i], tmp);
                    ptr_Al      += A_ld;
                };

                ptr_C           += C_ld;
            };
        };

        return;
    };

    template<bool Conj_A, bool Conj_B, class Alpha>
    static void eval_diag_2_gemm(const Alpha& alpha, Mat_D& C, const M1& A, const M2& B, 
                        trans_type t_A, trans_type t_B, Integer C_fr, Integer rows, Integer C_fc, Integer cols)
    {
        if (raw::is_diag(A) == true)
            return eval_diag_diag_gemm<Conj_A, Conj_B>(alpha, C, A, B, t_A, t_B, C_fr, rows, C_fc, cols);
        
        Integer K1          = std::min(B.rows(),B.cols());

        using V1            = typename M1::value_type;
        using V2            = typename M2::value_type;
        using V12           = typename md::unify_types<V1,V2>::type;
        using VR            = Val_ret;
        
        Integer A_ld        = A.ld();
        Integer B_ld        = B.ld();
        Integer C_ld        = C.ld();

        const V1* ptr_A     = A.rep_ptr();
        const V2* ptr_B     = B.rep_ptr() + B.first_elem_diag(0);
        Val_ret * ptr_C     = C.ptr() + C_fr + C_fc * C_ld;

        if (t_A == trans_type::no_trans)
        {
            for (Integer j = 0; j < K1; ++j)
            {
                V2 val_B        = make_conj<Conj_B>::eval(ptr_B[0]);

                if (mrd::is_zero(val_B) == false)
                {
                    Integer fr      = A.first_row(j);
                    Integer lr      = A.last_row(j);

                    Integer pos_A   = A.first_elem_pos(j);

                    for (Integer i = fr; i <= lr; ++i, ++pos_A)
                    {
                        V12 tmp     = make_conj<Conj_A>::eval(ptr_A[pos_A]) * val_B;
                        ptr_C[i]    = alpha.eval(ptr_C[i], tmp);
                    };
                };

                ptr_A           += A_ld;
                ptr_B           += B_ld;
                ptr_C           += C_ld;
            };
        }
        else
        {
            for (Integer j = 0; j < K1; ++j)
            {
                V2 val_B        = make_conj<Conj_B>::eval(ptr_B[0]);

                if (mrd::is_zero(val_B) == false)
                {
                    Integer fc      = A.first_col(j);
                    Integer lc      = A.last_col(j);

                    Integer pos_A   = A.first_elem_pos_row(j);

                    for (Integer i = fc; i <= lc; ++i, pos_A += A_ld-1)
                    {
                        V12 tmp     = make_conj<Conj_A>::eval(ptr_A[pos_A]) * val_B;
                        ptr_C[i]    = alpha.eval(ptr_C[i], tmp);
                    };
                };

                ptr_B           += B_ld;
                ptr_C           += C_ld;
            };
        };

        return;
    };

    static void eval(matcl::Matrix& ret, const M1& A, const M2& B, trans_type t_A, trans_type t_B)
    {
        if (raw::is_diag(A) == true)
            return eval_diag_1(ret, A, B, t_A, t_B);

        if (raw::is_diag(B) == true)
            return eval_diag_2(ret, A, B, t_A, t_B);

        using ti_ret_type   = typename ti::get_ti_type<Val_ret>::type;
        ti_ret_type ret_ti  = ti::get_return_ti<ti_ret_type>(mdyf::op_mul::eval(),
                                                                  ti::get_ti(A), ti::get_ti(B));

        Integer M, K, N, fd_A, ld_A, fd_B, ld_B;

        if (t_B == trans_type::no_trans)
        {
            N       = B.cols();
            fd_B    = B.first_diag();
            ld_B    = B.last_diag();
        }
        else
        {
            N       = B.rows();
            fd_B    = -B.last_diag();
            ld_B    = -B.first_diag();
        };

        if (t_A == trans_type::no_trans)
        {
            M       = A.rows();
            K       = A.cols();
            fd_A    = A.first_diag();
            ld_A    = A.last_diag();
        }
        else
        {
            K       = A.rows();
            M       = A.cols();
            fd_A    = -A.last_diag();
            ld_A    = -A.first_diag();
        }

        Val_ret Z   = md::default_value<Val_ret>(ret_ti);

        using Mat_B = Matrix<Val_ret,struct_banded>;
        Mat_B C(ret_ti, M, N, fd_A + fd_B, ld_A + ld_B);

        using val_type_1    = typename M1::value_type;
        using val_type_2    = typename M2::value_type;

        const val_type_1* ptr_A = A.rep_ptr();
        const val_type_2* ptr_B = B.rep_ptr();
        Val_ret * ptr_C         = C.rep_ptr();

        Integer A_ld    = A.ld();
        Integer B_ld    = B.ld();
        Integer C_ld    = C.ld();

        if (t_A == trans_type::no_trans)
        {
            if (t_B == trans_type::no_trans)
            {
                for (Integer j = 0; j < N; ++j)
                {
                    Integer ic      = C.first_elem_pos(j);
                    Integer B_fr    = B.first_row(j);
                    Integer B_lr    = B.last_row(j);
                    Integer C_fr    = C.first_row(j);
                    Integer C_lr    = C.last_row(j);

                    Integer ib0     = B.first_elem_pos(j) - B.first_row(j);

                    for (Integer i = C_fr; i <= C_lr; ++i, ++ic)
                    {				
                        ptr_C[ic]   = Z;

                        Integer ks  = std::max(A.first_col(i), B_fr);
                        Integer ke  = std::min(A.last_col(i), B_lr);
                        Integer ia  = A.first_elem_pos(ks) + i - A.first_row(ks) + imult(ks, A_ld);
                        Integer ib  = ib0 + ks;

                        //dot product
                        for (Integer k = ks; k <= ke; ++k, ia += A_ld-1, ++ib)
                        {
                            Val_ret tmp = ptr_A[ia] * ptr_B[ib];
                            ptr_C[ic]   = ptr_C[ic] + tmp;
                        };
                    };

                    ptr_B   += B_ld;
                    ptr_C   += C_ld;
                };
            }
            else
            {
                bool conj_B = (t_B == trans_type::conj_trans);

                for (Integer j = 0; j < N; ++j)
                {
                    Integer ic      = C.first_elem_pos(j);
                    Integer B_fr    = B.first_col(j);
                    Integer B_lr    = B.last_col(j);
                    Integer C_fr    = C.first_row(j);
                    Integer C_lr    = C.last_row(j);

                    for (Integer i = C_fr; i <= C_lr; ++i, ++ic)
                    {				
                        ptr_C[ic]   = Z;

                        Integer ks  = std::max(A.first_col(i), B_fr);
                        Integer ke  = std::min(A.last_col(i), B_lr);
                        Integer ia  = A.first_elem_pos(ks) + i - A.first_row(ks) + imult(ks, A_ld);
                        Integer ib  = B.element_pos(j,ks);

                        if (conj_B == true)
                        {
                            for (Integer k = ks; k <= ke; ++k, ia += A_ld-1, ib += B_ld - 1)
                            {
                                Val_ret tmp = ptr_A[ia] * conj(ptr_B[ib]);
                                ptr_C[ic]   = ptr_C[ic] + tmp;
                            };
                        }
                        else
                        {
                            for (Integer k = ks; k <= ke; ++k, ia += A_ld-1, ib += B_ld - 1)
                            {
                                Val_ret tmp = ptr_A[ia] * ptr_B[ib];
                                ptr_C[ic]   = ptr_C[ic] + tmp;
                            };
                        };
                    };

                    ptr_C   += C_ld;
                };
            };
        }
        else
        {
            //t_A = conj_trans or trans
            if (t_B == trans_type::no_trans)
            {
                bool conj_A = (t_A == trans_type::conj_trans);

                for (Integer j = 0; j < N; ++j)
                {
                    Integer ic      = C.first_elem_pos(j);
                    Integer B_fr    = B.first_row(j);
                    Integer B_lr    = B.last_row(j);
                    Integer C_fr    = C.first_row(j);
                    Integer C_lr    = C.last_row(j);

                    Integer ib0     = B.first_elem_pos(j) - B.first_row(j);

                    for (Integer i = C_fr; i <= C_lr; ++i, ++ic)
                    {				
                        ptr_C[ic]   = Z;

                        Integer ks  = std::max(A.first_row(i), B_fr);
                        Integer ke  = std::min(A.last_row(i), B_lr);
                        Integer ia  = A.first_elem_pos(i) + ks - A.first_row(i) + imult(i, A_ld);
                        Integer ib  = ib0 + ks;

                        //dot product
                        if (conj_A == true)
                        {
                            for (Integer k = ks; k <= ke; ++k, ia += 1, ++ib)
                            {
                                Val_ret tmp = conj(ptr_A[ia]) * ptr_B[ib];
                                ptr_C[ic]   = ptr_C[ic] + tmp;
                            };
                        }
                        else
                        {
                            for (Integer k = ks; k <= ke; ++k, ia += 1, ++ib)
                            {
                                Val_ret tmp = ptr_A[ia] * ptr_B[ib];
                                ptr_C[ic]   = ptr_C[ic] + tmp;
                            };
                        };
                    };

                    ptr_B   += B_ld;
                    ptr_C   += C_ld;
                };
            }
            else
            {
                bool conj_A = (t_A == trans_type::conj_trans);
                bool conj_B = (t_B == trans_type::conj_trans);

                for (Integer j = 0; j < N; ++j)
                {
                    Integer ic      = C.first_elem_pos(j);
                    Integer B_fr    = B.first_col(j);
                    Integer B_lr    = B.last_col(j);
                    Integer C_fr    = C.first_row(j);
                    Integer C_lr    = C.last_row(j);

                    for (Integer i = C_fr; i <= C_lr; ++i, ++ic)
                    {				
                        ptr_C[ic]   = Z;

                        Integer ks  = std::max(A.first_row(i), B_fr);
                        Integer ke  = std::min(A.last_row(i), B_lr);
                        Integer ia  = A.first_elem_pos(i) + ks - A.first_row(i) + imult(i, A_ld);
                        Integer ib  = B.element_pos(j,ks);

                        if (conj_A == true && conj_B == true)
                        {
                            for (Integer k = ks; k <= ke; ++k, ++ia, ib += B_ld-1)
                            {
                                Val_ret tmp = conj(ptr_A[ia]) * conj(ptr_B[ib]);
                                ptr_C[ic]   = ptr_C[ic] + tmp;
                            };
                        }
                        else if (conj_A == true && conj_B == false)
                        {
                            for (Integer k = ks; k <= ke; ++k, ++ia, ib += B_ld-1)
                            {
                                Val_ret tmp = conj(ptr_A[ia]) * ptr_B[ib];
                                ptr_C[ic]   = ptr_C[ic] + tmp;
                            };
                        }
                        else if (conj_A == false && conj_B == true)
                        {
                            for (Integer k = ks; k <= ke; ++k, ++ia, ib += B_ld-1)
                            {
                                Val_ret tmp = ptr_A[ia] * conj(ptr_B[ib]);
                                ptr_C[ic]   = ptr_C[ic] + tmp;
                            };
                        }
                        else
                        {
                            for (Integer k = ks; k <= ke; ++k, ++ia, ib += B_ld-1)
                            {
                                Val_ret tmp = ptr_A[ia] * ptr_B[ib];
                                ptr_C[ic]   = ptr_C[ic] + tmp;
                            };
                        };
                    };

                    ptr_C   += C_ld;
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

    static void eval_gemm(Mat_D& C, const Val_ret& alpha, const M1& A, const M2& B, 
                     trans_type t_A, trans_type t_B, const Val_ret& beta, 
                     Integer fr, Integer rows, Integer fc, Integer cols)
    {
        prepare_gemm_C<Val_ret>::eval(beta, C, fr,rows,fc,cols);

        bool is_diag_A  = raw::is_diag(A);
        bool is_diag_B  = raw::is_diag(B);

        bool conj_A     = t_A == trans_type::conj_trans;
        bool conj_B     = t_B == trans_type::conj_trans;

        if (mrd::is_one(alpha))
        {
            using Alpha     = alpha_one<Val_ret>;
            Alpha al;

            if (is_diag_A || is_diag_B)
            {
                if (conj_A == true)
                {
                    if (conj_B == true)
                    {
                        if (is_diag_A)
                            eval_diag_1_gemm<true,true>(al, C, A, B, t_A, t_B, fr, rows, fc, cols);
                        else
                            eval_diag_2_gemm<true,true>(al, C, A, B, t_A, t_B, fr, rows, fc, cols);
                    }
                    else
                    {
                        if (is_diag_A)
                            eval_diag_1_gemm<true,false>(al, C, A, B, t_A, t_B, fr, rows, fc, cols);
                        else
                            eval_diag_2_gemm<true,false>(al, C, A, B, t_A, t_B, fr, rows, fc, cols);
                    }
                }
                else
                {
                    if (conj_B == true)
                    {
                        if (is_diag_A)
                            eval_diag_1_gemm<false,true>(al, C, A, B, t_A, t_B, fr, rows, fc, cols);
                        else
                            eval_diag_2_gemm<false,true>(al, C, A, B, t_A, t_B, fr, rows, fc, cols);
                    }
                    else
                    {
                        if (is_diag_A)
                            eval_diag_1_gemm<false,false>(al, C, A, B, t_A, t_B, fr, rows, fc, cols);
                        else
                            eval_diag_2_gemm<false,false>(al, C, A, B, t_A, t_B, fr, rows, fc, cols);
                    }
                };
            }
            else
            {
                eval_impl_gemm(al, C, A, B, t_A, t_B, fr, rows, fc, cols);
            };
        }
        else if (mrd::is_one(-alpha))
        {
            using Alpha     = alpha_mone<Val_ret>;
            Alpha al;

            if (is_diag_A || is_diag_B)
            {
                if (conj_A == true)
                {
                    if (conj_B == true)
                    {
                        if (is_diag_A)
                            eval_diag_1_gemm<true,true>(al, C, A, B, t_A, t_B, fr, rows, fc, cols);
                        else
                            eval_diag_2_gemm<true,true>(al, C, A, B, t_A, t_B, fr, rows, fc, cols);
                    }
                    else
                    {
                        if (is_diag_A)
                            eval_diag_1_gemm<true,false>(al, C, A, B, t_A, t_B, fr, rows, fc, cols);
                        else
                            eval_diag_2_gemm<true,false>(al, C, A, B, t_A, t_B, fr, rows, fc, cols);
                    }
                }
                else
                {
                    if (conj_B == true)
                    {
                        if (is_diag_A)
                            eval_diag_1_gemm<false,true>(al, C, A, B, t_A, t_B, fr, rows, fc, cols);
                        else
                            eval_diag_2_gemm<false,true>(al, C, A, B, t_A, t_B, fr, rows, fc, cols);
                    }
                    else
                    {
                        if (is_diag_A)
                            eval_diag_1_gemm<false,false>(al, C, A, B, t_A, t_B, fr, rows, fc, cols);
                        else
                            eval_diag_2_gemm<false,false>(al, C, A, B, t_A, t_B, fr, rows, fc, cols);
                    }
                };
            }
            else
            {
                eval_impl_gemm(al, C, A, B, t_A, t_B, fr, rows, fc, cols);
            };
        }
        else
        {
            using Alpha     = alpha_val<Val_ret>;
            Alpha al(alpha);

            if (is_diag_A || is_diag_B)
            {
                if (conj_A == true)
                {
                    if (conj_B == true)
                    {
                        if (is_diag_A)
                            eval_diag_1_gemm<true,true>(al, C, A, B, t_A, t_B, fr, rows, fc, cols);
                        else
                            eval_diag_2_gemm<true,true>(al, C, A, B, t_A, t_B, fr, rows, fc, cols);
                    }
                    else
                    {
                        if (is_diag_A)
                            eval_diag_1_gemm<true,false>(al, C, A, B, t_A, t_B, fr, rows, fc, cols);
                        else
                            eval_diag_2_gemm<true,false>(al, C, A, B, t_A, t_B, fr, rows, fc, cols);
                    }
                }
                else
                {
                    if (conj_B == true)
                    {
                        if (is_diag_A)
                            eval_diag_1_gemm<false,true>(al, C, A, B, t_A, t_B, fr, rows, fc, cols);
                        else
                            eval_diag_2_gemm<false,true>(al, C, A, B, t_A, t_B, fr, rows, fc, cols);
                    }
                    else
                    {
                        if (is_diag_A)
                            eval_diag_1_gemm<false,false>(al, C, A, B, t_A, t_B, fr, rows, fc, cols);
                        else
                            eval_diag_2_gemm<false,false>(al, C, A, B, t_A, t_B, fr, rows, fc, cols);
                    }
                };
            }
            else
            {
                eval_impl_gemm(al, C, A, B, t_A, t_B, fr, rows, fc, cols);
            };
        };
    };

    template<class Alpha>
    static void eval_impl_gemm(const Alpha& alpha, Mat_D& C, const M1& A, const M2& B, 
                     trans_type t_A, trans_type t_B, Integer fr, Integer rows, Integer fc, Integer cols)
    {
        (void)rows;

        Integer N   = cols;

        Val_ret Z   = md::default_value<Val_ret>(C.get_type());

        using V1    = typename M1::value_type;
        using V2    = typename M2::value_type;
        using V12   = typename md::unify_types<V1,V2>::type;

        Integer A_ld        = A.ld();
        Integer B_ld        = B.ld();
        Integer C_ld        = C.ld();
        const V1* ptr_A     = A.rep_ptr();
        const V2* ptr_B     = B.rep_ptr();
        Val_ret * ptr_C     = C.ptr() + fr + fc * C_ld;

        if (t_B == trans_type::no_trans)
        {
            if (t_A == trans_type::no_trans)
            {
                for (Integer j = 0; j < N; ++j)
                {
                    Integer B_fr    = B.first_row(j);
                    Integer B_lr    = B.last_row(j);
                    Integer B_pos   = B.first_elem_pos(j);

                    const V1* ptr_Al    = ptr_A + B_fr * A_ld;

                    for (Integer i = B_fr; i <= B_lr; ++i, ++B_pos, ptr_Al += A_ld)
                    {				
                        const V2& Bi    = ptr_B[B_pos];

                        if (mrd::is_zero(Bi))
                            continue;

                        Integer A_fr    = A.first_row(i);
                        Integer A_lr    = A.last_row(i);
                        Integer A_pos   = A.first_elem_pos(i);

                        for (Integer k = A_fr; k <= A_lr; ++k, ++A_pos)
                        {
                            V12 tmp     = ptr_Al[A_pos] * Bi;
                            ptr_C[k]    = alpha.eval(ptr_C[k], tmp);
                        };
                    };

                    ptr_B   += B_ld;
                    ptr_C   += C_ld;
                };
            }
            else
            {
                bool conj_A = (t_A == trans_type::conj_trans);

                for (Integer j = 0; j < N; ++j)
                {
                    Integer B_fr    = B.first_row(j);
                    Integer B_lr    = B.last_row(j);
                    Integer B_pos   = B.first_elem_pos(j) - B_fr;

                    Integer A_fc    = A.first_col(B_fr);
                    Integer A_lc    = A.last_col(B_lr);

                    const V1* ptr_Al    = ptr_A + A_fc * A_ld;

                    if (conj_A)
                    {
                        for (Integer i  = A_fc; i <= A_lc; ++i)
                        {
                            Integer A_fr    = A.first_row(i);
                            Integer A_lr    = A.last_row(i);
                            Integer A_pos   = A.first_elem_pos(i) - A_fr;

                            Integer s       = std::max(B_fr, A_fr);
                            Integer e       = std::min(B_lr, A_lr);
                            Integer ka      = A_pos + s;
                            Integer kb      = B_pos + s;

                            Val_ret dot     = Z;

                            for (Integer k = s; k <= e; ++k, ++ka, ++kb)
                                dot         = dot + Val_ret(conj(ptr_Al[ka]) * ptr_B[kb]);

                            ptr_C[i]        = alpha.eval(ptr_C[i], dot);
                            ptr_Al          += A_ld;
                        };
                    }
                    else
                    {
                        for (Integer i  = A_fc; i <= A_lc; ++i)
                        {                        
                            Integer A_fr    = A.first_row(i);
                            Integer A_lr    = A.last_row(i);
                            Integer A_pos   = A.first_elem_pos(i) - A_fr;

                            Integer s       = std::max(B_fr, A_fr);
                            Integer e       = std::min(B_lr, A_lr);
                            Integer ka      = A_pos + s;
                            Integer kb      = B_pos + s;

                            Val_ret dot     = Z;

                            for (Integer k = s; k <= e; ++k, ++ka, ++kb)
                                dot         = dot + Val_ret(ptr_Al[ka] * ptr_B[kb]);

                            ptr_C[i]        = alpha.eval(ptr_C[i], dot);
                            ptr_Al          += A_ld;
                        };
                    };

                    ptr_B   += B_ld;
                    ptr_C   += C_ld;
                };
            };
        }
        else
        {
            bool conj_B = (t_B == trans_type::conj_trans);

            if (t_A == trans_type::no_trans)
            {                
                for (Integer j = 0; j < N; ++j)
                {
                    Integer B_fc    = B.first_col(j);
                    Integer B_lc    = B.last_col(j);
                    Integer B_pos   = B.first_elem_pos_row(j);

                    const V1* ptr_Al    = ptr_A + B_fc * A_ld;

                    for (Integer i = B_fc; i <= B_lc; ++i, B_pos += B_ld - 1, ptr_Al += A_ld)
                    {				
                        V2 Bi       = conj_B ? conj(ptr_B[B_pos]) : ptr_B[B_pos];

                        if (mrd::is_zero(Bi))
                            continue;

                        Integer A_fr    = A.first_row(i);
                        Integer A_lr    = A.last_row(i);
                        Integer A_pos   = A.first_elem_pos(i);

                        for (Integer k = A_fr; k <= A_lr; ++k, ++A_pos)
                        {
                            V12 tmp     = ptr_Al[A_pos] * Bi;
                            ptr_C[k]    = alpha.eval(ptr_C[k], tmp);
                        };
                    };

                    ptr_C   += C_ld;
                };
            }
            else
            {
                bool conj_A = (t_A == trans_type::conj_trans);

                if (conj_A == true)
                {
                    if (conj_B == true)
                        eval_TT_gemm<true,true>(alpha, C, A, B, fr, rows, fc, cols);
                    else
                        eval_TT_gemm<true,false>(alpha, C, A, B, fr, rows, fc, cols);
                }
                else
                {
                    if (conj_B == true)
                        eval_TT_gemm<false,true>(alpha, C, A, B, fr, rows, fc, cols);
                    else
                        eval_TT_gemm<false,false>(alpha, C, A, B, fr, rows, fc, cols);
                };
            };
        };        

        return;
    };

    template<bool Conj_A, bool Conj_B, class Alpha>
    static void eval_TT_gemm(const Alpha& alpha, Mat_D& C, const M1& A, const M2& B, 
                     Integer fr, Integer rows, Integer fc, Integer cols)
    {
        (void)rows;

        Integer N   = cols;

        Val_ret Z   = md::default_value<Val_ret>(C.get_type());

        using V1    = typename M1::value_type;
        using V2    = typename M2::value_type;
        using V12   = typename md::unify_types<V1,V2>::type;

        Integer A_ld        = A.ld();
        Integer B_ld        = B.ld();
        Integer C_ld        = C.ld();
        const V1* ptr_A     = A.rep_ptr();
        const V2* ptr_B     = B.rep_ptr();
        Val_ret * ptr_C     = C.ptr() + fr + fc * C_ld;

        for (Integer j = 0; j < N; ++j)
        {
            Integer B_fc    = B.first_col(j);
            Integer B_lc    = B.last_col(j);

            Integer A_fc    = A.first_col(B_fc);
            Integer A_lc    = A.last_col(B_lc);

            const V1* ptr_Al    = ptr_A + A_fc * A_ld;

            for (Integer i  = A_fc; i <= A_lc; ++i)
            {                        
                Integer A_fr    = A.first_row(i);
                Integer A_lr    = A.last_row(i);
                Integer A_pos   = A.first_elem_pos(i) - A_fr;

                Integer s       = std::max(B_fc, A_fr);
                Integer e       = std::min(B_lc, A_lr);
                Integer ka      = A_pos + s;
                Integer kb      = B.element_pos(j,s);

                Val_ret dot     = Z;

                for (Integer k = s; k <= e; ++k, ++ka, kb += B_ld-1)
                {
                    dot         = dot + Val_ret(make_conj<Conj_A>::eval(ptr_Al[ka]) 
                                                * make_conj<Conj_B>::eval(ptr_B[kb]));
                };

                ptr_C[i]        = alpha.eval(ptr_C[i], dot);
                ptr_Al          += A_ld;
            };

            ptr_C   += C_ld;
        };
    };
};

}}}