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
//                      BAND - SPARSE
//========================================================================

//get integer value from template argument, if Val != 0, or from dynamic argument val
// if Val == 0.
template<Integer Val>
struct get_int      { static Integer eval(Integer)      { return Val; }; };
template<>
struct get_int<0>   { static Integer eval(Integer val)  { return val; }; };

template<class Val_ret, class M1, class M2>
struct eval_mult<Val_ret,M1,M2,struct_banded,struct_sparse> 
{
    using Mat_D = raw::Matrix<Val_ret, struct_dense>;

    template<bool conj_A, Integer A_ld0>
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

        if (M == 0 || K == 0 || N == 0 || B.nnz() == 0)
        {
            using val_type_ret_real = typename md::real_type<VTR>::type;
            using sparse_matrix     = Matrix<val_type_ret_real,struct_sparse>;

            sparse_matrix out(ret_ti,M,N);
            ret = matcl::Matrix(out,false);
            return;
        };

        Integer K1              = std::min(M,K);
        const VT1* ptr_A        = A.rep_ptr() + A.first_elem_diag(0);
        Integer A_ld            = get_int<A_ld0>::eval(A.ld());

        if (t_B == trans_type::no_trans)
        {
            using ret_type  = raw::Matrix<Val_ret, struct_sparse>;
            ret_type C(ret_ti, M, N, B.nnz());

            sparse_ccs<VTR>& d          = C.rep();
            const sparse_ccs<VT2>& Bd   = B.rep();

            Integer nz              = 0;

            Integer * d_c           = d.ptr_c();
            Integer * d_r           = d.ptr_r();
            VTR * d_x               = d.ptr_x();

            const Integer* Bd_c     = Bd.ptr_c();
            const Integer* Bd_r     = Bd.ptr_r();
            const VT2* Bd_x         = Bd.ptr_x();            

            for (Integer j = 0; j < N; ++j)
            {
                d_c[j]              = nz;

                for (Integer k = Bd_c[j] ; k < Bd_c[j + 1] ; ++k)
                {
                    const VT2& Bx   = Bd_x[k];
                    Integer Br      = Bd_r[k];

                    if (Br >= K1)
                        continue;
                    
                    const VT1& v    = make_conj<conj_A>::eval(ptr_A[Br*A_ld]);
                    VTR tmp         = v * Bx;
                    
                    if (mrd::is_zero(tmp))
                        continue;

                    d_r[nz]         = Br;
                    d_x[nz]         = tmp;
                    ++nz;
                }
            };
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

        using ret_type  = raw::Matrix<Val_ret, struct_sparse>;
        ret_type Bt(ret_ti);

        if (t_B == trans_type::trans)
            mrd::manip_trans_reshaper_helper<VT2, VTR>::eval_trans(Bt, B, B.rows(), K1, M, N);
        else
            mrd::manip_trans_reshaper_helper<VT2, VTR>::eval_ctrans(Bt, B, B.rows(), K1, M, N);

        sparse_ccs<VTR>& Yd = Bt.rep();
        Integer * Y_c       = Yd.ptr_c();
        Integer * Y_r       = Yd.ptr_r();
        VTR* Y_x            = Yd.ptr_x();

        Integer j           = 0;
        Integer i           = 0;
        Integer nz          = 0;

        for (; j < N; ++j)
        {
            nz              = Y_c[j];

            for (i = Y_c[j]; i < Y_c[j + 1]; ++i)
            {
                Integer r       = Y_r[i];
                const VT1& v    = make_conj<conj_A>::eval(ptr_A[r*A_ld]);
                VTR tmp         = v * Y_x[i];

                if (mrd::is_zero(tmp))
                    goto lab_remove;

                Y_x[i]      = tmp;
                ++nz;
            }

            continue;

          lab_remove:
            nz              = i;

            for (; i < Y_c[j + 1]; ++i)
            {
                Integer r       = Y_r[i];
                const VT1& v    = make_conj<conj_A>::eval(ptr_A[r*A_ld]);
                VTR tmp         = v * Y_x[i];

                if (mrd::is_zero(tmp))
                    continue;

                Y_x[nz]     = tmp;
                Y_r[nz]     = r;
                ++nz;
            }

            ++j;
            break;
        };

        for (; j < N; ++j)
        {
            Integer pos_s   = Y_c[j];
            Integer pos_end = Y_c[j + 1];
            Y_c[j]          = nz;

            for (i = pos_s; i < pos_end; ++i)
            {
                Integer r       = Y_r[i];
                const VT1& v    = make_conj<conj_A>::eval(ptr_A[r*A_ld]);
                VTR tmp         = v * Y_x[i];

                if (mrd::is_zero(tmp))
                    continue;

                Y_x[nz]     = tmp;
                Y_r[nz]     = r;
                ++nz;
            }
        };
        Y_c[N]              = nz;

        bool is_sq_A        = A.rows() == A.cols();
        bool is_sq_B        = B.rows() == B.cols();
        bool is_sq_C        = Bt.rows() == Bt.cols();

        Bt.set_struct(predefined_struct_ext::mult_struct(A.get_struct(),B.get_struct(),t_A,t_B,
                            is_real_matrix(A),is_real_matrix(B), is_sq_A, is_sq_B, is_sq_C));

        ret = matcl::Matrix(Bt,true);
        return;
    };

    template<bool conj_A, Integer A_ld0>
    static void eval_diag_gemm(const Val_ret& alpha, Mat_D& C, const M1& A, const M2& B, trans_type t_A, 
                               trans_type t_B, Integer fr, Integer rows, Integer fc, Integer cols)
    {
        (void)t_A;

        using VTR           = Val_ret;
        using VT1           = typename M1::value_type;
        using VT2           = typename M2::value_type;
        using V12           = typename md::unify_types<VT1, VT2>::type;

        Integer M           = rows;
        Integer N           = cols;

        if (M == 0 || N == 0 || B.nnz() == 0)
            return;

        Integer K1          = std::min(A.rows(), A.cols());
        const VT1* ptr_A    = A.rep_ptr() + A.first_elem_diag(0);
        Integer A_ld        = get_int<A_ld0>::eval(A.ld());

        Integer C_ld        = C.ld();
        VTR* ptr_C          = C.ptr() + fr + fc * C_ld;

        const sparse_ccs<VT2>& Bd   = B.rep();

        const Integer* Bd_c = Bd.ptr_c();
        const Integer* Bd_r = Bd.ptr_r();
        const VT2* Bd_x     = Bd.ptr_x();            

        if (t_B == trans_type::no_trans)
        {
            for (Integer j = 0; j < N; ++j)
            {
                for (Integer k = Bd_c[j] ; k < Bd_c[j + 1] ; ++k)
                {
                    const VT2& Bx   = Bd_x[k];
                    Integer Br      = Bd_r[k];

                    if (Br >= K1)
                        continue;
                    
                    const VT1& v    = make_conj<conj_A>::eval(ptr_A[Br*A_ld]);
                    V12 tmp         = v * Bx;
                    
                    if (mrd::is_zero(tmp))
                        continue;

                    ptr_C[Br]       = ptr_C[Br] + alpha * tmp;
                }

                ptr_C               += C_ld;
            };

            return;
        }
        else if (t_B == trans_type::trans)
        {
            for (Integer j = 0; j < K1; ++j)
            {
                const VT1& v        = make_conj<conj_A>::eval(ptr_A[j*A_ld]);

                if (mrd::is_zero(v))
                {
                    ptr_C           += 1;
                    continue;
                };

                for (Integer k = Bd_c[j] ; k < Bd_c[j + 1] ; ++k)
                {
                    const VT2& Bx   = Bd_x[k];
                    Integer Br      = Bd_r[k];

                    V12 tmp         = v * Bx;
                    
                    ptr_C[Br*C_ld]  = ptr_C[Br*C_ld] + alpha * tmp;
                }

                ptr_C               += 1;
            };

            return;
        }
        else
        {
            //t_B == trans_type::conj_trans
            for (Integer j = 0; j < K1; ++j)
            {
                const VT1& v        = make_conj<conj_A>::eval(ptr_A[j*A_ld]);

                if (mrd::is_zero(v))
                {
                    ptr_C           += 1;
                    continue;
                };

                for (Integer k = Bd_c[j] ; k < Bd_c[j + 1] ; ++k)
                {
                    const VT2& Bx   = Bd_x[k];
                    Integer Br      = Bd_r[k];

                    V12 tmp         = v * conj(Bx);
                    
                    ptr_C[Br*C_ld]  = ptr_C[Br*C_ld] + alpha * tmp;
                }

                ptr_C               += 1;
            };

            return;
        };
    };

    static void eval(matcl::Matrix& ret, const M1& A, const M2& B, trans_type t_A, trans_type t_B)
    {
        bool is_diag        = raw::is_diag(A);

        if (is_diag)
        {
            if (A.ld() == 1)
            {
                if (t_A == trans_type::conj_trans)
                    return eval_diag<true,1>(ret, A, B, t_A, t_B);
                else
                    return eval_diag<false,1>(ret, A, B, t_A, t_B);
            }
            else
            {
                if (t_A == trans_type::conj_trans)
                    return eval_diag<true,0>(ret, A, B, t_A, t_B);
                else
                    return eval_diag<false,0>(ret, A, B, t_A, t_B);
            };
        };

        if (t_B == trans_type::no_trans)
            return eval_notrans(ret, A, B, t_A, t_B);
        else
            return eval_trans(ret, A, B, t_A, t_B);
    };

    static void eval_gemm(Mat_D& C, const Val_ret& alpha, const M1& A, const M2& B, 
                     trans_type t_A, trans_type t_B, const Val_ret& beta, 
                     Integer fr, Integer rows, Integer fc, Integer cols)
    {
        bool is_diag        = raw::is_diag(A);

        prepare_gemm_C<Val_ret>::eval(beta, C, fr,rows,fc,cols);

        if (is_diag)
        {
            if (A.ld() == 1)
            {
                if (t_A == trans_type::conj_trans)
                    return eval_diag_gemm<true,1>(alpha, C, A, B, t_A, t_B, fr, rows, fc, cols);
                else
                    return eval_diag_gemm<false,1>(alpha, C, A, B, t_A, t_B, fr, rows, fc, cols);
            }
            else
            {
                if (t_A == trans_type::conj_trans)
                    return eval_diag_gemm<true,0>(alpha, C, A, B, t_A, t_B, fr, rows, fc, cols);
                else
                    return eval_diag_gemm<false,0>(alpha, C, A, B, t_A, t_B, fr, rows, fc, cols);
            };
        };

        if (mrd::is_one(alpha))
        {
            using Alpha     = alpha_one<Val_ret>;
            Alpha al;

            if (t_B == trans_type::no_trans)
                return eval_notrans_gemm(al, C, A, B, t_A, t_B, fr, rows, fc, cols);
            else if (t_A == trans_type::no_trans)
                return eval_notrans_trans_gemm(al, C, A, B, t_A, t_B, fr, rows, fc, cols);
            else
                return eval_trans_trans_gemm(al, C, A, B, t_A, t_B, fr, rows, fc, cols);
        }
        else if (mrd::is_one(-alpha))
        {
            using Alpha     = alpha_mone<Val_ret>;
            Alpha al;

            if (t_B == trans_type::no_trans)
                return eval_notrans_gemm(al, C, A, B, t_A, t_B, fr, rows, fc, cols);
            else if (t_A == trans_type::no_trans)
                return eval_notrans_trans_gemm(al, C, A, B, t_A, t_B, fr, rows, fc, cols);
            else
                return eval_trans_trans_gemm(al, C, A, B, t_A, t_B, fr, rows, fc, cols);
        }
        else
        {
            using Alpha     = alpha_val<Val_ret>;
            Alpha al(alpha);

            if (t_B == trans_type::no_trans)
                return eval_notrans_gemm(al, C, A, B, t_A, t_B, fr, rows, fc, cols);
            else if (t_A == trans_type::no_trans)
                return eval_notrans_trans_gemm(al, C, A, B, t_A, t_B, fr, rows, fc, cols);
            else
                return eval_trans_trans_gemm(al, C, A, B, t_A, t_B, fr, rows, fc, cols);
        };
    };

    static void eval_trans(matcl::Matrix& ret, const M1& A, const M2& B, trans_type t_A, trans_type t_B)
    {
        using VTR           = Val_ret;
        using VT1           = typename M1::value_type;
        using VT2           = typename M2::value_type;
        using ti_ret_type   = typename ti::get_ti_type<Val_ret>::type;

        ti_ret_type ret_ti  = ti::get_return_ti<ti_ret_type>(mdyf::op_mul::eval(),
                                                                  ti::get_ti(A), ti::get_ti(B));        

        M2 Bt(ret_ti);

        if (t_B == trans_type::trans)
            mrd::manip_trans_helper<M2>::eval_trans(Bt, B);
        else
            mrd::manip_trans_helper<M2>::eval_ctrans(Bt, B);

        return eval_notrans(ret, A, Bt, t_A, trans_type::no_trans);
    };

    template<class Alpha>
    static void eval_notrans_trans_gemm(const Alpha& alpha, Mat_D& C, const M1& A, const M2& B, trans_type t_A, 
                               trans_type t_B, Integer fr, Integer rows, Integer fc, Integer cols)
    {
        (void)t_A;

        using VTR           = Val_ret;
        using VT1           = typename M1::value_type;
        using VT2           = typename M2::value_type;
        using V12           = typename md::unify_types<VT1, VT2>::type;

        Integer M           = rows;
        Integer N           = cols;
        Integer B_cols      = B.cols();
        (void)B_cols;

        if (M == 0 || N == 0 || B.nnz() == 0)        
            return;

        Integer A_ld        = A.ld();
        Integer C_ld        = C.ld();        
        VTR* ptr_C          = C.ptr() + fr + fc * C_ld;

        const sparse_ccs<VT2>& Bd   = B.rep();
        const Integer* Bd_c = Bd.ptr_c();
        const Integer* Bd_r = Bd.ptr_r();
        const VT2* Bd_x     = Bd.ptr_x();

        bool conj_B         = (t_B == trans_type::conj_trans);

        Integer first_col_A = A.first_nonzero_column();
        Integer last_col_A  = A.last_nonzero_column();

        const VT1* ptr_A    = A.rep_ptr() + first_col_A * A_ld;

        for (Integer j = first_col_A; j <= last_col_A; ++j)
        {
            Integer first_row   = A.first_row(j);
            Integer last_row    = A.last_row(j);
            Integer pos0        = A.first_elem_pos(j);

            for (Integer k = Bd_c[j] ; k < Bd_c[j + 1] ; ++k)
            {
                VT2 Bx          = conj_B ? conj(Bd_x[k]) : Bd_x[k];
                Integer Br      = Bd_r[k];

                if (mrd::is_zero(Bx))
                    continue;

                VTR* ptr_Cl     = ptr_C + Br * C_ld;

                for (Integer i = first_row, pos = pos0; i <= last_row; ++i, ++pos)
                {
                    V12 tmp     = ptr_A[pos] * Bx;
                    ptr_Cl[i]   = alpha.eval(ptr_Cl[i], tmp);
                };
            };

            ptr_A               += A_ld;
        };

        return;
    };

    template<class Alpha>
    static void eval_trans_trans_gemm(const Alpha& alpha, Mat_D& C, const M1& A, const M2& B, trans_type t_A, 
                               trans_type t_B, Integer fr, Integer rows, Integer fc, Integer cols)
    {
        using VTR           = Val_ret;
        using VT1           = typename M1::value_type;
        using VT2           = typename M2::value_type;
        using V12           = typename md::unify_types<VT1, VT2>::type;

        Integer M           = rows;
        Integer N           = cols;
        Integer B_cols      = B.cols();
        (void)B_cols;

        if (M == 0 || N == 0 || B.nnz() == 0)        
            return;

        Integer A_ld        = A.ld();
        Integer C_ld        = C.ld();
        const VT1* ptr_A    = A.rep_ptr();
        VTR* ptr_C          = C.ptr() + fr + fc * C_ld;

        const sparse_ccs<VT2>& Bd   = B.rep();
        const Integer* Bd_c = Bd.ptr_c();
        const Integer* Bd_r = Bd.ptr_r();
        const VT2* Bd_x     = Bd.ptr_x();

        bool conj_B         = (t_B == trans_type::conj_trans);
        bool conj_A         = (t_A == trans_type::conj_trans);

        Integer first_row_A     = A.first_nonzero_row();
        Integer last_row_A      = A.last_nonzero_row();

        if (conj_A == false)
        {
            for (Integer j = first_row_A; j <= last_row_A; ++j)
            {
                Integer first_col   = A.first_col(j);
                Integer last_col    = A.last_col(j);
                Integer pos0        = A.first_elem_pos_row(j);

                for (Integer k = Bd_c[j] ; k < Bd_c[j + 1] ; ++k)
                {
                    VT2 Bx          = conj_B ? conj(Bd_x[k]) : Bd_x[k];
                    Integer Br      = Bd_r[k];

                    if (mrd::is_zero(Bx))
                        continue;

                    VTR* ptr_Cl     = ptr_C + Br * C_ld;
                    Integer pos     = pos0;

                    for (Integer i = first_col; i <= last_col; ++i, pos += A_ld - 1)
                    {
                        V12 tmp     = ptr_A[pos] * Bx;
                        ptr_Cl[i]   = alpha.eval(ptr_Cl[i], tmp);
                    };
                };
            };

            return;
        }
        else
        {
            for (Integer j = first_row_A; j <= last_row_A; ++j)
            {
                Integer first_col   = A.first_col(j);
                Integer last_col    = A.last_col(j);
                Integer pos0        = A.first_elem_pos_row(j);

                for (Integer k = Bd_c[j] ; k < Bd_c[j + 1] ; ++k)
                {
                    VT2 Bx          = conj_B ? conj(Bd_x[k]) : Bd_x[k];
                    Integer Br      = Bd_r[k];

                    if (mrd::is_zero(Bx))
                        continue;

                    VTR* ptr_Cl     = ptr_C + Br * C_ld;
                    Integer pos     = pos0;

                    for (Integer i = first_col; i <= last_col; ++i, pos += A_ld - 1)
                    {
                        V12 tmp     = conj(ptr_A[pos]) * Bx;
                        ptr_Cl[i]   = alpha.eval(ptr_Cl[i], tmp);
                    };
                };
            };

            return;
        };
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

        if (M == 0 || K == 0 || N == 0 || B.nnz() == 0)
        {
            using val_type_ret_real = typename md::real_type<VTR>::type;
            using sparse_matrix     = Matrix<val_type_ret_real,struct_sparse>;

            sparse_matrix out(ret_ti,M,N);
            ret = matcl::Matrix(out,false);
            return;
        };

        const VT1* ptr_A    = A.rep_ptr();
        Integer A_ld        = A.ld();

        using scatter       = matcl::algorithm::scatter;
        using ret_type      = raw::Matrix<Val_ret,struct_sparse>;

        ret_type C(ret_ti, M, N, M + estim_mult_nnz(A, B, t_A, t_B));

        sparse_ccs<VTR>& d          = C.rep();
        const sparse_ccs<VT2>& Bd   = B.rep();

        md::workspace2<VTR>  work_x(ret_ti,M);

        scatter sc          = scatter::get(M, N);
        Integer nz          = 0;

        Integer * d_c       = d.ptr_c();
        Integer * d_r       = d.ptr_r();
        VTR * d_x           = d.ptr_x();

        const Integer* Bd_c = Bd.ptr_c();
        const Integer* Bd_r = Bd.ptr_r();
        const VT2* Bd_x     = Bd.ptr_x();

        if (t_A == trans_type::no_trans)
        {
            for (Integer j = 0; j < N; ++j)
            {
                if (nz + M > d.nzmax()) 
                {
                    d.add_memory( d.nzmax() + M);

                    d_r         = d.ptr_r();
                    d_x         = d.ptr_x();
                };

                d_c[j]          = nz;
                Integer nz_old  = nz;
                Integer nk      = 0;
                auto mark       = sc.next_mark();

                for (Integer k = Bd_c[j] ; k < Bd_c[j + 1] ; ++k)
                {
                    const VT2& Bx   = Bd_x[k];
                    Integer Br      = Bd_r[k];

                    if (mrd::is_zero(Bx))
                        continue;

                    Integer nzs     = nz;

                    Integer first_row   = A.first_row(Br);
                    Integer last_row    = A.last_row(Br);
                    Integer pos         = A.first_elem_pos(Br) + imult(Br,A_ld);

                    for (Integer ka = first_row; ka <= last_row; ++ka, ++pos)
                    {
                        const VT1& v    = ptr_A[pos];

                        if (mrd::is_zero(v))
                            continue;

                        if (sc[ka] < mark)
                        {
                            sc[ka]      = mark;
                            work_x[ka]  = v * Bx;
                            d_r[nz]     = ka;
                            ++nz;
                        }
                        else 
                        {
                            VTR tmp     = v * Bx;
                            work_x[ka]  = mrd::plus_helper<VTR,VTR>::eval(work_x[ka],tmp);
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
        }
        else
        {
            bool conj_A = (t_A == trans_type::conj_trans);

            for (Integer j = 0; j < N; ++j)
            {
                if (nz + M > d.nzmax()) 
                {
                    d.add_memory( d.nzmax() + M);

                    d_r         = d.ptr_r();
                    d_x         = d.ptr_x();
                };

                d_c[j]          = nz;
                Integer nz_old  = nz;
                Integer nk      = 0;
                auto mark       = sc.next_mark();

                for (Integer k = Bd_c[j] ; k < Bd_c[j + 1] ; ++k)
                {
                    const VT2& Bx   = Bd_x[k];
                    Integer Br      = Bd_r[k];

                    if (mrd::is_zero(Bx))
                        continue;

                    Integer nzs     = nz;

                    Integer first_row   = A.first_col(Br);
                    Integer last_row    = A.last_col(Br);
                    Integer pos         = A.first_elem_pos_row(Br);

                    for (Integer ka = first_row; ka <= last_row; ++ka, pos += A_ld - 1)
                    {
                        const VT1& v    = conj_A ? conj(ptr_A[pos]) : ptr_A[pos];

                        if (mrd::is_zero(v))
                            continue;

                        if (sc[ka] < mark)
                        {
                            sc[ka]      = mark;
                            work_x[ka]  = v * Bx;
                            d_r[nz]     = ka;
                            ++nz;
                        }
                        else 
                        {
                            VTR tmp     = v * Bx;
                            work_x[ka]  = mrd::plus_helper<VTR,VTR>::eval(work_x[ka],tmp);
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

    template<class Alpha>
    static void eval_notrans_gemm(const Alpha& alpha, Mat_D& C, const M1& A, const M2& B, trans_type t_A, 
                               trans_type t_B, Integer fr, Integer rows, Integer fc, Integer cols)
    {
        (void)t_B;

        using VTR           = Val_ret;
        using VT1           = typename M1::value_type;
        using VT2           = typename M2::value_type;
        using V12           = typename md::unify_types<VT1, VT2>::type;

        Integer M           = rows;
        Integer N           = cols;

        if (M == 0 || N == 0 || B.nnz() == 0)        
            return;

        Integer A_ld        = A.ld();
        Integer C_ld        = C.ld();
        const VT1* ptr_A    = A.rep_ptr();
        VTR* ptr_C          = C.ptr() + fr + fc * C_ld;

        const sparse_ccs<VT2>& Bd   = B.rep();
        const Integer* Bd_c = Bd.ptr_c();
        const Integer* Bd_r = Bd.ptr_r();
        const VT2* Bd_x     = Bd.ptr_x();

        if (t_A == trans_type::no_trans)
        {
            for (Integer j = 0; j < N; ++j)
            {
                for (Integer k = Bd_c[j] ; k < Bd_c[j + 1] ; ++k)
                {
                    const VT2& Bx   = Bd_x[k];
                    Integer Br      = Bd_r[k];

                    if (mrd::is_zero(Bx))
                        continue;

                    Integer first_row   = A.first_row(Br);
                    Integer last_row    = A.last_row(Br);
                    Integer pos         = A.first_elem_pos(Br) + imult(Br,A_ld);

                    for (Integer ka = first_row; ka <= last_row; ++ka, ++pos)
                    {
                        const VT1& v    = ptr_A[pos];
                        V12 tmp         = v * Bx;
                        ptr_C[ka]       = alpha.eval(ptr_C[ka], tmp);
                    };
                }

                ptr_C   += C_ld;
            };

            return;
        }
        else
        {
            bool conj_A = (t_A == trans_type::conj_trans);

            for (Integer j = 0; j < N; ++j)
            {
                for (Integer k = Bd_c[j] ; k < Bd_c[j + 1] ; ++k)
                {
                    const VT2& Bx   = Bd_x[k];
                    Integer Br      = Bd_r[k];

                    if (mrd::is_zero(Bx))
                        continue;

                    Integer first_row   = A.first_col(Br);
                    Integer last_row    = A.last_col(Br);
                    Integer pos         = A.first_elem_pos_row(Br);

                    for (Integer ka = first_row; ka <= last_row; ++ka, pos += A_ld - 1)
                    {
                        const VT1& v    = conj_A ? conj(ptr_A[pos]) : ptr_A[pos];
                        V12 tmp         = v * Bx;
                        ptr_C[ka]       = alpha.eval(ptr_C[ka], tmp);
                    };
                }
                
                ptr_C   += C_ld;
            };

            return;
        };        
    };
};

}}}