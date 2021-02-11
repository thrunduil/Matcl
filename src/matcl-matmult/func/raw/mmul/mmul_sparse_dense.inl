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
//                      SPARSE - DENSE
//========================================================================

template<class Val_ret, class M1, class M2>
struct eval_mult<Val_ret,M1,M2,struct_sparse,struct_dense> 
{
    using Mat_D = raw::Matrix<Val_ret, struct_dense>;

    static M2 get_diag(const M2& X, trans_type t_A)
    {
        using VT2           = typename M2::value_type;

        VT2 Z1              = md::default_value<VT2>(X.get_type());
        Integer K           = std::min(X.rows(), X.cols());

        M2 D(X.get_type(), Z1, K, 1);

        VT2* ptr_D          = D.ptr();
        const VT2* ptr_X    = X.ptr();
        Integer X_ld        = X.ld();

        if (t_A == trans_type::conj_trans)
        {
            for (Integer i = 0; i < K; ++i)
            {
                ptr_D[i]    = conj(ptr_X[0]);
                ptr_X       += X_ld + 1;
            };
        }
        else
        {
            for (Integer i = 0; i < K; ++i)
            {
                ptr_D[i]    = ptr_X[0];
                ptr_X       += X_ld + 1;
            };
        };

        return D;
    };

    static void eval_diag(matcl::Matrix& ret, const M1& A,const M2& X, trans_type t_A, trans_type t_B)
    {
        using VTR           = Val_ret;
        using VT1           = typename M1::value_type;
        using VT2           = typename M2::value_type;
        using ti_ret_type   = ti::get_ti_type<Val_ret>::type;

        Integer M, K, N;

        if (t_B == trans_type::no_trans)
            N = X.cols();
        else
            N = X.rows();

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

        Integer K1  = std::min(N,K);

        ti_ret_type ret_ti  = ti::get_return_ti<ti_ret_type>(mdyf::op_mul::eval(),
                                                                  ti::get_ti(X), ti::get_ti(A));

        if (M == 0 || N == 0 || K == 0 || A.nnz() == 0)
        {
            using val_type_ret_real = typename md::real_type<VTR>::type;
            using sparse_matrix     = Matrix<val_type_ret_real,struct_sparse>;

            sparse_matrix out(ret_ti,M,N);
            ret = matcl::Matrix(out,false);
            return;
        };

        using SM    = raw::Matrix<VTR, struct_sparse>;

        const VT2* ptr_X    = X.ptr();
        Integer X_ld        = X.ld();

        if (t_A == trans_type::no_trans)
        {
            SM Y(ret_ti, M, N, A.nnz());

            const sparse_ccs<VT1>& d= A.rep();
            const Integer * A_c     = d.ptr_c();
            const Integer * A_r     = d.ptr_r();
            const VT1* A_x          = d.ptr_x();

            sparse_ccs<VTR>& Yd = Y.rep();
            Integer * Y_c       = Yd.ptr_c();
            Integer * Y_r       = Yd.ptr_r();
            VTR* Y_x            = Yd.ptr_x();

            Integer nz          = 0;

            bool conj_X         = (t_B == trans_type::conj_trans);

            for (Integer j = 0; j < K1; ++j)
            {
                Y_c[j]          = nz;

                const VT2& val  = conj_X? conj(ptr_X[0]) : ptr_X[0];

                if (mrd::is_zero(val))
                {
                    ptr_X       += X_ld + 1;
                    continue;
                };

                for (Integer i = A_c[j]; i < A_c[j + 1]; ++i)
                {
                    Y_r[nz]     = A_r[i];
                    Y_x[nz]     = A_x[i] * val;
                    ++nz;
                };

                ptr_X           += X_ld + 1;
            };

            for (Integer j = K1; j <= N; ++j)
                Y_c[j]          = nz;

            bool is_sq_A        = A.rows() == A.cols();
            bool is_sq_X        = X.rows() == X.cols();
            bool is_sq_Y        = Y.rows() == Y.cols();
            Y.set_struct(predefined_struct_ext::mult_struct(A.get_struct(),X.get_struct(),t_A,t_B,
                                is_real_matrix(A),is_real_matrix(X), is_sq_A, is_sq_X, is_sq_Y));

            ret = matcl::Matrix(Y,true);
            return;
        }

        SM At(ret_ti);

        if (t_A == trans_type::trans)
            mrd::manip_trans_reshaper_helper<VT1, VTR>::eval_trans(At, A, K1, A.cols(), M, N);
        else
            mrd::manip_trans_reshaper_helper<VT1, VTR>::eval_ctrans(At, A, K1, A.cols(), M, N);

        sparse_ccs<VTR>& Yd = At.rep();
        Integer * Y_c       = Yd.ptr_c();
        Integer * Y_r       = Yd.ptr_r();
        VTR* Y_x            = Yd.ptr_x();

        bool conj_X         = (t_B == trans_type::conj_trans);
        Integer j           = 0;
        Integer nz          = 0;

        //inplace modification
        for (; j < K1; ++j)
        {            
            const VT2& val  = conj_X? conj(ptr_X[0]) : ptr_X[0];
            nz              = Y_c[j];

            if (mrd::is_zero(val))
            {
                ptr_X       += X_ld + 1;                
                ++j;
                break;
            };

            for (Integer i = Y_c[j]; i < Y_c[j + 1]; ++i)
                Y_x[i]      = Y_x[i] * val;

            ptr_X           += X_ld + 1;
            nz              = Y_c[j + 1];
        };

        //need to move elements
        for (; j < K1; ++j)
        {
            Integer pos_s   = Y_c[j];
            Y_c[j]          = nz;

            const VT2& val  = conj_X? conj(ptr_X[0]) : ptr_X[0];

            if (mrd::is_zero(val))
            {
                ptr_X       += X_ld + 1;
                continue;
            };

            for (Integer i = pos_s; i < Y_c[j + 1]; ++i)
            {
                Y_r[nz]     = Y_r[i];
                Y_x[nz]     = Y_x[i] * val;
                ++nz;
            };

            ptr_X           += X_ld + 1;
        };

        for (j = K1; j <= N; ++j)
            Y_c[j]          = nz;

        bool is_sq_A        = A.rows() == A.cols();
        bool is_sq_X        = X.rows() == X.cols();
        bool is_sq_C        = At.rows() == At.cols();
        At.set_struct(predefined_struct_ext::mult_struct(A.get_struct(),X.get_struct(),t_A,t_B,
                            is_real_matrix(A),is_real_matrix(X), is_sq_A, is_sq_X, is_sq_C));

        ret = matcl::Matrix(At,true);
        return;
    };

    template<class Alpha>
    static void eval_diag_gemm(Mat_D& Y, const Alpha& alpha, const M1& A, const M2& X, 
                               trans_type t_A, trans_type t_B, Integer fr, Integer rows, Integer fc, 
                               Integer cols)
    {
        using VTR           = Val_ret;
        using VT1           = typename M1::value_type;
        using VT2           = typename M2::value_type;
        using V12           = typename md::unify_types<VT1, VT2>::type;
        using ti_ret_type   = ti::get_ti_type<Val_ret>::type;

        Integer M           = rows;
        Integer N           = cols;

        Integer K1          = std::min(X.rows(),X.cols());

        if (M == 0 || N == 0 || K1 == 0 || A.nnz() == 0)
            return;

        Integer X_ld        = X.ld();
        Integer Y_ld        = Y.ld();
        const VT2* ptr_X    = X.ptr();
        Val_ret* ptr_Y      = Y.ptr() + fr + fc * Y_ld;

        const sparse_ccs<VT1>& d= A.rep();
        const Integer * A_c     = d.ptr_c();
        const Integer * A_r     = d.ptr_r();
        const VT1* A_x          = d.ptr_x();

        if (t_A == trans_type::no_trans)
        {
            bool conj_X         = (t_B == trans_type::conj_trans);

            for (Integer j = 0; j < K1; ++j)
            {
                const VT2& val  = conj_X? conj(ptr_X[0]) : ptr_X[0];

                if (mrd::is_zero(val))
                {
                    ptr_X       += X_ld + 1;
                    continue;
                };

                for (Integer i = A_c[j]; i < A_c[j + 1]; ++i)
                {
                    Integer l   = A_r[i];
                    ptr_Y[l]    = alpha.eval(ptr_Y[l], A_x[i] * val);
                };

                ptr_X           += X_ld + 1;
                ptr_Y           += Y_ld;
            };

            return;
        }

        M2 D                = get_diag(X, t_B);
        VT2* ptr_D          = D.ptr();

        if (t_A == trans_type::trans)
        {
            if (Y_ld == 1)
            {
                for (Integer j = 0; j < M; ++j)
                {
                    for (Integer i = A_c[j]; i < A_c[j + 1]; ++i)
                    {
                        Integer l   = A_r[i];

                        if (l >= K1)
                            continue;

                        const VT2& val  = ptr_D[l];
                        V12 tmp         = A_x[i] * val;
                        ptr_Y[l]        = alpha.eval(ptr_Y[l], tmp);
                    };

                    ptr_Y           += 1;
                };
            }
            else
            {
                for (Integer j = 0; j < M; ++j)
                {
                    for (Integer i = A_c[j]; i < A_c[j + 1]; ++i)
                    {
                        Integer l   = A_r[i];

                        if (l >= K1)
                            continue;

                        const VT2& val  = ptr_D[l];
                        V12 tmp         = A_x[i] * val;
                        ptr_Y[l*Y_ld]   = alpha.eval(ptr_Y[l*Y_ld], tmp);
                    };

                    ptr_Y           += 1;
                };
            };
        }
        else
        {
            if (Y_ld == 1)
            {
                for (Integer j = 0; j < M; ++j)
                {
                    for (Integer i = A_c[j]; i < A_c[j + 1]; ++i)
                    {
                        Integer l   = A_r[i];

                        if (l >= K1)
                            continue;

                        const VT2& val  = ptr_D[l];
                        V12 tmp         = conj(A_x[i]) * val;
                        ptr_Y[l]        = alpha.eval(ptr_Y[l], tmp);
                    };

                    ptr_Y           += 1;
                };
            }
            else
            {
                for (Integer j = 0; j < M; ++j)
                {
                    for (Integer i = A_c[j]; i < A_c[j + 1]; ++i)
                    {
                        Integer l   = A_r[i];

                        if (l >= K1)
                            continue;

                        const VT2& val  = ptr_D[l];
                        V12 tmp         = conj(A_x[i]) * val;
                        ptr_Y[l*Y_ld]   = alpha.eval(ptr_Y[l*Y_ld], tmp);
                    };

                    ptr_Y           += 1;
                };
            };
        };

        return;
    };
    static void eval(matcl::Matrix& ret, const M1& A,const M2& X, trans_type t_A, trans_type t_B)
    {
        bool is_diag        = raw::is_diag(X);

        if (is_diag)
            return eval_diag(ret, A, X, t_A, t_B);

        using ti_ret_type   = ti::get_ti_type<Val_ret>::type;

        ti_ret_type ret_ti  = ti::get_return_ti<ti_ret_type>(mdyf::op_mul::eval(),
                                                                  ti::get_ti(A), ti::get_ti(X));

        Integer M, K, N;

        if (t_B == trans_type::no_trans)
            N = X.cols();
        else
            N = X.rows();

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

        Val_ret Z = md::default_value<Val_ret>(ret_ti);        

        if (M == 0 || N == 0 || K == 0 || A.nnz() == 0)
        {
            using Val_ret_real  = typename md::real_type<Val_ret>::type;
            using sparse_matrix = Matrix<Val_ret,struct_sparse>;

            sparse_matrix out(ret_ti,M,N);
            ret = matcl::Matrix(out,false);
            return;
        };

        using ret_type  = raw::Matrix<Val_ret,struct_dense>;
        ret_type Y(ret_ti, Z, M, N);

        Integer fr      = 0;
        Integer rows    = M;
        Integer fc      = 0;
        Integer cols    = N;

        using Alpha     = alpha_one<Val_ret>;
        Alpha alpha;

        if (t_A == trans_type::no_trans)
        {
            if (t_B == trans_type::no_trans)
                eval_notrans_notrans(alpha, A, X, t_A, t_B,Y,fr,rows,fc,cols);
            else
                eval_notrans_trans(alpha, A, X, t_A, t_B,Y,fr,rows,fc,cols);
        }
        else
        {
            if (t_B == trans_type::no_trans)
                eval_trans_notrans(alpha, A, X, t_A, t_B,Y,fr,rows,fc,cols);
            else if (t_B == trans_type::trans)
                eval_trans_trans<Alpha,false>(alpha, A, X, t_A, t_B,Y,fr,rows,fc,cols);
            else
                eval_trans_trans<Alpha,true>(alpha, A, X, t_A, t_B,Y,fr,rows,fc,cols);
        }

        bool is_sq_A    = A.rows() == A.cols();
        bool is_sq_X    = X.rows() == X.cols();
        bool is_sq_Y    = Y.rows() == Y.cols();
        Y.set_struct(predefined_struct_ext::mult_struct(A.get_struct(),X.get_struct(),t_A,t_B,
                            is_real_matrix(A),is_real_matrix(X), is_sq_A, is_sq_X, is_sq_Y));

        ret = matcl::Matrix(Y,true);
        return;
    };

    template<class Alpha, bool Conj_X>
    static void eval_trans_trans(const Alpha& alpha, const M1& A,const M2& X, 
                            trans_type t_A, trans_type t_B, Mat_D& Y, Integer fr, Integer rows, 
                            Integer fc, Integer cols)
    {
        (void)t_B;

        using VTR           = Val_ret;
        using VT1           = typename M1::value_type;
        using VT2           = typename M2::value_type;

        const sparse_ccs<VT1>& d= A.rep();
        const Integer * A_c     = d.ptr_c();
        const Integer * A_r     = d.ptr_r();
        const VT1* A_x          = d.ptr_x();

        Integer X_ld            = X.ld();
        Integer Y_ld            = Y.ld();

        const VT2* ptr_X        = X.ptr();
        VTR* ptr_Y              = Y.ptr() + fr + fc * Y_ld;

        Integer M               = rows;
        Integer N               = cols;

        Val_ret Z               = md::default_value<Val_ret>(Y.get_type());

        if( t_A == trans_type::trans)
        {
            if (X_ld == 1)
            {
                for (Integer j = 0; j < N; ++j)
                {
                    for (Integer k = 0; k < M; ++k)
                    {
                        VTR dot         = Z;
                        Integer last    = A_c[k + 1];

                        for (Integer i = A_c[k]; i < last; ++i)
                        {
                            dot = dot + A_x[i] * make_conj<Conj_X>::eval(ptr_X[A_r[i]]);
                        };

                        ptr_Y[k]        = ptr_Y[k] + alpha.eval(dot);
                    }

                    ptr_X += 1;
                    ptr_Y += Y_ld;
                }
            }
            else
            {
                for (Integer j = 0; j < N; ++j)
                {
                    for (Integer k = 0; k < M; ++k)
                    {
                        VTR dot         = Z;
                        Integer last    = A_c[k + 1];

                        for (Integer i = A_c[k]; i < last; ++i)
                        {
                            dot = dot + A_x[i] * make_conj<Conj_X>::eval(ptr_X[A_r[i]*X_ld]);
                        };

                        ptr_Y[k]        = ptr_Y[k] + alpha.eval(dot);
                    }

                    ptr_X += 1;
                    ptr_Y += Y_ld;
                };
            }
        }
        else 
        {
            if (X_ld == 1)
            {
                for (Integer j = 0; j < N; ++j)
                {
                    for (Integer k = 0; k < M; ++k)
                    {
                        VTR dot         = Z;
                        Integer last    = A_c[k + 1];

                        for (Integer i = A_c[k]; i < last; ++i)
                        {
                            VT1 a   = mrd::conj_helper<VT1>::eval(A_x[i]);
                            dot     = dot + a * make_conj<Conj_X>::eval(ptr_X[A_r[i]]);
                        };

                        ptr_Y[k]    = ptr_Y[k] + alpha.eval(dot);
                    }
                    ptr_X += 1;
                    ptr_Y += Y_ld;
                };
            }
            else
            {
                for (Integer j = 0; j < N; ++j)
                {
                    for (Integer k = 0; k < M; ++k)
                    {
                        VTR dot         = Z;
                        Integer last    = A_c[k + 1];

                        for (Integer i = A_c[k]; i < last; ++i)
                        {
                            VT1 a   = mrd::conj_helper<VT1>::eval(A_x[i]);
                            dot     = dot + a * make_conj<Conj_X>::eval(ptr_X[A_r[i]*X_ld]);
                        };

                        ptr_Y[k]        = ptr_Y[k] + alpha.eval(dot);
                    }

                    ptr_X += 1;
                    ptr_Y += Y_ld;
                };
            };
        }

        return;
    };

    template<class Alpha>
    static void eval_trans_notrans(const Alpha& alpha, const M1& A,const M2& X, 
                            trans_type t_A, trans_type t_B, Mat_D& Y, Integer fr, Integer rows, 
                            Integer fc, Integer cols)
    {
        (void)t_B;

        using VTR           = Val_ret;
        using VT1           = typename M1::value_type;
        using VT2           = typename M2::value_type;
        using ti_ret_type   = ti::get_ti_type<Val_ret>::type;

        const sparse_ccs<VT1>& d     = A.rep();
        const Integer * A_c     = d.ptr_c();
        const Integer * A_r     = d.ptr_r();
        const VT1* A_x          = d.ptr_x();

        Integer X_ld            = X.ld();
        Integer Y_ld            = Y.ld();

        const VT2* ptr_X        = X.ptr();
        VTR* ptr_Y              = Y.ptr() + fr + fc * Y_ld;

        Integer M       = rows;
        Integer N       = cols;

        Val_ret Z       = md::default_value<Val_ret>(Y.get_type());

        if( t_A == trans_type::trans)
        {
            for (Integer j = 0; j < N; ++j)
            {
                for (Integer k = 0; k < M; ++k)
                {
                    VTR dot         = Z;
                    Integer last    = A_c[k + 1];

                    for (Integer i = A_c[k]; i < last; ++i)
                    {
                        dot = dot + A_x[i] * ptr_X[A_r[i]];
                    };

                    ptr_Y[k] = alpha.eval(ptr_Y[k], dot);
                }

                ptr_X += X_ld;
                ptr_Y += Y_ld;
            };
        }
        else 
        {
            for (Integer j = 0; j < N; ++j)
            {
                for (Integer k = 0; k < M; ++k)
                {
                    VTR dot         = Z;
                    Integer last    = A_c[k + 1];

                    for (Integer i = A_c[k]; i < last; ++i)
                    {
                        VT1 a   = mrd::conj_helper<VT1>::eval(A_x[i]);
                        dot     = dot + a * ptr_X[A_r[i]];
                    };

                    ptr_Y[k] = alpha.eval(ptr_Y[k], dot);
                }
                ptr_X += X_ld;
                ptr_Y += Y_ld;
            };
        }

        return;
    };

    template<class Alpha>
    static void eval_notrans_notrans(const Alpha& alpha, const M1& A,const M2& X, 
                    trans_type t_A, trans_type t_B,Mat_D& Y, Integer fr, Integer rows, 
                    Integer fc, Integer cols)
    {
        (void)t_A;
        (void)t_B;
        (void)rows;

        using VTR           = Val_ret;
        using VT1           = typename M1::value_type;
        using VT2           = typename M2::value_type;
        using ti_ret_type   = ti::get_ti_type<Val_ret>::type;

        const sparse_ccs<VT1>& d     = A.rep();
        const Integer * A_c     = d.ptr_c();
        const Integer * A_r     = d.ptr_r();
        const VT1* A_x          = d.ptr_x();        

        Integer X_ld            = X.ld();
        Integer Y_ld            = Y.ld();

        VTR* ptr_Y              = Y.ptr() + fr + fc * Y_ld;
        const VT2* ptr_X        = X.ptr();        

        Integer N               = cols;
        Integer K               = A.cols();

        Real density        = Real(d.nnz()) / Real(K);
        Integer j           = 0;
        Integer Nl          = N - (N % 4);

        if (density < 2.)
        {
            //matrix is very sparse, use unblocked version
            Nl              = 0;
        };

        const Integer LDY   = Y.ld();
        const Integer LDX   = X.ld();

        Integer k0          = 0;
        Integer k1          = K;

        while(A_c[k0 + 1] == A_c[k0])
            ++k0;

        while(A_c[k1 - 1] == A_c[k1])
            --k1;

        for (; j < Nl; j += 4)
        {
            for (Integer k = k0; k < k1; ++k)
            {           
                const VT2& val0         = ptr_X[k + 0 * LDX];
                const VT2& val1         = ptr_X[k + 1 * LDX];
                const VT2& val2         = ptr_X[k + 2 * LDX];
                const VT2& val3         = ptr_X[k + 3 * LDX];

                if (mrd::is_zero(val0) || mrd::is_zero(val1)
                    ||mrd::is_zero(val2) ||mrd::is_zero(val3))
                {
                    goto lab_seq;
                };

                for (Integer i = A_c[k]; i < A_c[k + 1]; ++i)
                {
                    Integer r           = A_r[i];
                    VT1 x               = A_x[i];
                        
                    ptr_Y[r + 0*LDY]    = alpha.eval(ptr_Y[r + 0*LDY], x * val0);
                    ptr_Y[r + 1*LDY]    = alpha.eval(ptr_Y[r + 1*LDY], x * val1);
                    ptr_Y[r + 2*LDY]    = alpha.eval(ptr_Y[r + 2*LDY], x * val2);
                    ptr_Y[r + 3*LDY]    = alpha.eval(ptr_Y[r + 3*LDY], x * val3);
                };

                continue;

              lab_seq:
                    
                if (mrd::is_zero(val0) == false)
                {
                    for (Integer i = A_c[k]; i < A_c[k + 1]; ++i)
                    {
                        Integer r       = A_r[i];
                        VT1 x           = A_x[i];                        
                        ptr_Y[r+0*LDY]  = alpha.eval(ptr_Y[r+0*LDY], x * val0);
                    };
                };
                if (mrd::is_zero(val1) == false)
                {
                    for (Integer i = A_c[k]; i < A_c[k + 1]; ++i)
                    {
                        Integer r       = A_r[i];
                        VT1 x           = A_x[i];                        
                        ptr_Y[r+1*LDY]  = alpha.eval(ptr_Y[r+1*LDY], x * val1);
                    };
                };
                if (mrd::is_zero(val2) == false)
                {
                    for (Integer i = A_c[k]; i < A_c[k + 1]; ++i)
                    {
                        Integer r       = A_r[i];
                        VT1 x           = A_x[i];
                        ptr_Y[r+2*LDY]  = alpha.eval(ptr_Y[r+2*LDY], x * val2);
                    };
                };
                if (mrd::is_zero(val3) == false)
                {
                    for (Integer i = A_c[k]; i < A_c[k + 1]; ++i)
                    {
                        Integer r       = A_r[i];
                        VT1 x           = A_x[i];
                        ptr_Y[r+3*LDY]  = alpha.eval(ptr_Y[r+3*LDY], x * val3);
                    };
                };
            }

            ptr_X += LDX * 4;
            ptr_Y += LDY * 4;
        };
            
        for (; j < N; ++j)
        {
            for (Integer k = k0; k < k1; ++k)
            {
                const VT2& val      = ptr_X[k];

                if (mrd::is_zero(val))
                    continue;

                for (Integer i = A_c[k]; i < A_c[k + 1]; ++i)
                {
                    Integer r       = A_r[i];
                    VT1 x           = A_x[i];
                        
                    ptr_Y[r]        = alpha.eval(ptr_Y[r], x * val);
                };
            }

            ptr_X += X_ld;
            ptr_Y += Y_ld;
        };
    };

    template<class Alpha>
    static void eval_notrans_trans(const Alpha& alpha, const M1& A,const M2& X, 
                            trans_type t_A, trans_type t_B, Mat_D& Y, Integer fr, Integer rows, 
                            Integer fc, Integer cols)
    {
        (void)rows;
        (void)t_A;
        (void)t_B;

        using VTR           = Val_ret;
        using VT1           = typename M1::value_type;
        using VT2           = typename M2::value_type;
        using ti_ret_type   = ti::get_ti_type<Val_ret>::type;

        const sparse_ccs<VT1>& d     = A.rep();
        const Integer * A_c     = d.ptr_c();
        const Integer * A_r     = d.ptr_r();
        const VT1* A_x          = d.ptr_x();

        Integer Y_ld            = Y.ld();
        const VT2* ptr_X        = X.ptr();
        VTR* ptr_Y              = Y.ptr() + fr + fc * Y_ld;               

        Integer K               = A.cols();
        Integer N               = cols;

        Real density        = Real(d.nnz()) / Real(K);
        Integer j           = 0;
        Integer Nl          = N - (N % 4);

        if (density < 2.)
        {
            //matrix is very sparse, use unblocked version
            Nl              = 0;
        };

        const Integer LDY   = Y.ld();
        const Integer LDX   = X.ld();

        Integer k0          = 0;
        Integer k1          = K;

        bool conj_X         = (t_B == trans_type::conj_trans);

        while(A_c[k0 + 1] == A_c[k0])
            ++k0;

        while(A_c[k1 - 1] == A_c[k1])
            --k1;

        for (; j < Nl; j += 4)
        {
            for (Integer k = k0; k < k1; ++k)
            {
                const VT2& val0 = (conj_X == true)? conj(ptr_X[k * LDX + 0]): ptr_X[k * LDX + 0];
                const VT2& val1 = (conj_X == true)? conj(ptr_X[k * LDX + 1]): ptr_X[k * LDX + 1];
                const VT2& val2 = (conj_X == true)? conj(ptr_X[k * LDX + 2]): ptr_X[k * LDX + 2];
                const VT2& val3 = (conj_X == true)? conj(ptr_X[k * LDX + 3]): ptr_X[k * LDX + 3];

                if (mrd::is_zero(val0) || mrd::is_zero(val1)
                    ||mrd::is_zero(val2) ||mrd::is_zero(val3))
                {
                    goto lab_seq;
                };

                for (Integer i = A_c[k]; i < A_c[k + 1]; ++i)
                {
                    Integer r           = A_r[i];
                    VT1 x               = A_x[i];
                        
                    ptr_Y[r + 0*LDY]    = alpha.eval(ptr_Y[r + 0*LDY], x * val0);
                    ptr_Y[r + 1*LDY]    = alpha.eval(ptr_Y[r + 1*LDY], x * val1);
                    ptr_Y[r + 2*LDY]    = alpha.eval(ptr_Y[r + 2*LDY], x * val2);
                    ptr_Y[r + 3*LDY]    = alpha.eval(ptr_Y[r + 3*LDY], x * val3);
                };

                continue;

              lab_seq:
                    
                if (mrd::is_zero(val0) == false)
                {
                    for (Integer i = A_c[k]; i < A_c[k + 1]; ++i)
                    {
                        Integer r       = A_r[i];
                        VT1 x           = A_x[i];                        
                        ptr_Y[r+0*LDY]  = alpha.eval(ptr_Y[r+0*LDY], x * val0);
                    };
                };
                if (mrd::is_zero(val1) == false)
                {
                    for (Integer i = A_c[k]; i < A_c[k + 1]; ++i)
                    {
                        Integer r       = A_r[i];
                        VT1 x           = A_x[i];                        
                        ptr_Y[r+1*LDY]  = alpha.eval(ptr_Y[r+1*LDY], x * val1);
                    };
                };
                if (mrd::is_zero(val2) == false)
                {
                    for (Integer i = A_c[k]; i < A_c[k + 1]; ++i)
                    {
                        Integer r       = A_r[i];
                        VT1 x           = A_x[i];
                        ptr_Y[r+2*LDY]  = alpha.eval(ptr_Y[r+2*LDY], x * val2);
                    };
                };
                if (mrd::is_zero(val3) == false)
                {
                    for (Integer i = A_c[k]; i < A_c[k + 1]; ++i)
                    {
                        Integer r       = A_r[i];
                        VT1 x           = A_x[i];
                        ptr_Y[r+3*LDY]  = alpha.eval(ptr_Y[r+3*LDY], x * val3);
                    };
                };
            }

            ptr_X += 4;
            ptr_Y += LDY * 4;
        };
            
        for (; j < N; ++j)
        {
            for (Integer k = k0; k < k1; ++k)
            {
                const VT2& val      = (conj_X == true)? conj(ptr_X[k*LDX]) : ptr_X[k*LDX];

                if (mrd::is_zero(val))
                    continue;

                for (Integer i = A_c[k]; i < A_c[k + 1]; ++i)
                {
                    Integer r       = A_r[i];
                    VT1 x           = A_x[i];
                        
                    ptr_Y[r]        = alpha.eval(ptr_Y[r], x * val);
                };
            }

            ptr_X += 1;
            ptr_Y += Y_ld;
        };
    };

    static void eval_gemm(Mat_D& Y, const Val_ret& alpha, const M1& A, const M2& X, 
                     trans_type t_A, trans_type t_B, const Val_ret& beta, 
                     Integer fr, Integer rows, Integer fc, Integer cols)
    {
        if (rows == 0 || cols == 0)
            return;

        prepare_gemm_C<Val_ret>::eval(beta, Y, fr,rows,fc,cols);

        if (A.nnz() == 0)
            return;

        bool is_diag        = raw::is_diag(X);

        if (mrd::is_one(alpha) == true)
        {
            using Alpha     = alpha_one<Val_ret>;
            Alpha al;

            if (is_diag)
                return eval_diag_gemm(Y, al, A, X, t_A, t_B, fr, rows ,fc, cols);

            if (t_A == trans_type::no_trans)
            {
                if (t_B == trans_type::no_trans)
                    eval_notrans_notrans(al, A, X, t_A, t_B,Y,fr,rows,fc,cols);
                else
                    eval_notrans_trans(al, A, X, t_A, t_B,Y,fr,rows,fc,cols);
            }
            else
            {
                if (t_B == trans_type::no_trans)
                    eval_trans_notrans(al, A, X, t_A, t_B,Y,fr,rows,fc,cols);
                else if (t_B == trans_type::trans)
                    eval_trans_trans<Alpha,false>(al, A, X, t_A, t_B,Y,fr,rows,fc,cols);
                else
                    eval_trans_trans<Alpha,true>(al, A, X, t_A, t_B,Y,fr,rows,fc,cols);
            }
        }
        else if (mrd::is_one(-alpha) == true)
        {
            using Alpha     = alpha_mone<Val_ret>;
            Alpha al;

            if (is_diag)
                return eval_diag_gemm(Y, al, A, X, t_A, t_B, fr, rows ,fc, cols);

            if (t_A == trans_type::no_trans)
            {
                if (t_B == trans_type::no_trans)
                    eval_notrans_notrans(al, A, X, t_A, t_B,Y,fr,rows,fc,cols);
                else
                    eval_notrans_trans(al, A, X, t_A, t_B,Y,fr,rows,fc,cols);
            }
            else
            {
                if (t_B == trans_type::no_trans)
                    eval_trans_notrans(al, A, X, t_A, t_B,Y,fr,rows,fc,cols);
                else if (t_B == trans_type::trans)
                    eval_trans_trans<Alpha,false>(al, A, X, t_A, t_B,Y,fr,rows,fc,cols);
                else
                    eval_trans_trans<Alpha,true>(al, A, X, t_A, t_B,Y,fr,rows,fc,cols);
            }
        }
        else
        {
            using Alpha     = alpha_val<Val_ret>;
            Alpha al(alpha);

            if (is_diag)
                return eval_diag_gemm(Y, al, A, X, t_A, t_B, fr, rows ,fc, cols);

            if (t_A == trans_type::no_trans)
            {
                if (t_B == trans_type::no_trans)
                    eval_notrans_notrans(al, A, X, t_A, t_B,Y,fr,rows,fc,cols);
                else
                    eval_notrans_trans(al, A, X, t_A, t_B,Y,fr,rows,fc,cols);
            }
            else
            {
                if (t_B == trans_type::no_trans)
                    eval_trans_notrans(al, A, X, t_A, t_B,Y,fr,rows,fc,cols);
                else if (t_B == trans_type::trans)
                    eval_trans_trans<Alpha,false>(al, A, X, t_A, t_B,Y,fr,rows,fc,cols);
                else
                    eval_trans_trans<Alpha,true>(al, A, X, t_A, t_B,Y,fr,rows,fc,cols);
            }
        };

        return;
    };
};

template<class Val_C, class M1, class M2>
struct eval_mult_abs<Val_C, M1, M2, struct_sparse>
{ 
    using Mat_C = raw::Matrix<Val_C, struct_dense>;

    static void eval_scal(matcl::Matrix& out, const M1& A0, const M2& scal, trans_type t_A, const Mat_C& C)
    {
        if (C.rows() == 0 || C.cols() == 0)
        {
            out = matcl::Matrix(C, false);
            return;
        }
        if (mrd::is_zero<M2>(scal) == true || A0.nnz() == 0)
        {
            mrd::scalfunc_real_helper<Mat_C>::eval_abs(out, C);
            return;
        };

        using V1            = typename M1::value_type;
        using V2            = M2;
        using VR            = typename md::real_type<Val_C>::type;
        using V1R           = typename md::real_type<V1>::type;
        using Mat_ret       = raw::Matrix<VR, struct_dense>;

        Integer N           = C.cols();
        Integer M           = C.rows();

        M1 A(A0.get_type());

        if (t_A != trans_type::no_trans)
            mrd::manip_trans_helper<M1>::eval_trans(A, A0);
        else
            A.assign_to_fresh(A0);

        const sparse_ccs<V1>& d = A.rep();
        const Integer * A_c     = d.ptr_c();
        const Integer * A_r     = d.ptr_r();
        const V1* A_x           = d.ptr_x();        

        Mat_ret ret(C.get_type(), M, N);

        Integer C_ld            = C.ld();
        Integer ret_ld          = ret.ld();

        const Val_C* ptr_C      = C.ptr();
        VR* ptr_ret             = ret.ptr();

        //make abs(C);
        for (Integer i = 0; i < N; ++i)
        {
            for (Integer j = 0; j < M; ++j)
                ptr_ret[j]      = abs_helper<Val_C>::eval(ptr_C[j]);

            ptr_C               += C_ld;
            ptr_ret             += ret_ld;
        };

        ptr_C                   = C.ptr();
        ptr_ret                 = ret.ptr();

        for (Integer j = 0; j < N; ++j)
        {
            for (Integer i = A_c[j]; i < A_c[j + 1]; ++i)
            {
                Integer r       = A_r[i];
                V1R x           = abs_helper<V1>::eval(A_x[i]);
                        
                ptr_ret[r]      = ptr_ret[r] + x * scal;
            };

            ptr_ret     += ret_ld;
        };

        out = matcl::Matrix(ret, false);
        return;
    }

    static void eval(matcl::Matrix& ret, const M1& A, const M2& X, trans_type t_A, const Mat_C& C)
    {
        if (C.rows() == 0 || C.cols() == 0)
        {
            ret = matcl::Matrix(C, false);
            return;
        }

        if (A.nnz() == 0)
        {
            mrd::scalfunc_real_helper<Mat_C>::eval_abs(ret, C);
            return;
        }

        if (t_A == trans_type::no_trans)
            eval_notrans(ret, A, X, C);
        else
            eval_trans(ret, A, X, C);

        return;
    }

    static void eval_notrans(matcl::Matrix& out, const M1& A, const M2& X, const Mat_C& C)
    {
        using V1            = typename M1::value_type;
        using V2            = typename M2::value_type;
        using VR            = typename md::real_type<Val_C>::type;
        using V1R           = typename md::real_type<V1>::type;
        using V2R           = typename md::real_type<V2>::type;
        using Mat_ret       = raw::Matrix<VR, struct_dense>;

        Integer N               = C.cols();
        Integer M               = C.rows();
        Integer K               = A.cols();

        const sparse_ccs<V1>& d = A.rep();
        const Integer * A_c     = d.ptr_c();
        const Integer * A_r     = d.ptr_r();
        const V1* A_x           = d.ptr_x();        

        Mat_ret ret(C.get_type(), M, N);

        Integer X_ld            = X.ld();
        Integer C_ld            = C.ld();
        Integer ret_ld          = ret.ld();

        const Val_C* ptr_C      = C.ptr();
        const V2* ptr_X         = X.ptr();
        VR* ptr_ret             = ret.ptr();

        //make abs(Y);
        for (Integer i = 0; i < N; ++i)
        {
            for (Integer j = 0; j < M; ++j)
                ptr_ret[j]      = abs_helper<Val_C>::eval(ptr_C[j]);

            ptr_C               += C_ld;
            ptr_ret             += ret_ld;
        };

        ptr_C                   = C.ptr();
        ptr_ret                 = ret.ptr();

        Real density            = Real(d.nnz()) / Real(K);
        Integer j               = 0;
        Integer Nl              = N - (N % 4);

        if (density < 2.)
        {
            //matrix is very sparse, use unblocked version
            Nl              = 0;
        };

        const Integer LDR   = ret.ld();
        const Integer LDX   = X.ld();

        Integer k0          = 0;
        Integer k1          = K;

        while(A_c[k0 + 1] == A_c[k0])
            ++k0;

        while(A_c[k1 - 1] == A_c[k1])
            --k1;

        for (; j < Nl; j += 4)
        {
            for (Integer k = k0; k < k1; ++k)
            {           
                const V2R& val0      = abs_helper<V2>::eval(ptr_X[k + 0 * LDX]);
                const V2R& val1      = abs_helper<V2>::eval(ptr_X[k + 1 * LDX]);
                const V2R& val2      = abs_helper<V2>::eval(ptr_X[k + 2 * LDX]);
                const V2R& val3      = abs_helper<V2>::eval(ptr_X[k + 3 * LDX]);

                if (mrd::is_zero(val0) || mrd::is_zero(val1)
                    ||mrd::is_zero(val2) ||mrd::is_zero(val3))
                {
                    goto lab_seq;
                };

                for (Integer i = A_c[k]; i < A_c[k + 1]; ++i)
                {
                    Integer r           = A_r[i];
                    V1R x               = abs_helper<V1>::eval(A_x[i]);
                        
                    ptr_ret[r + 0*LDR]  = ptr_ret[r + 0*LDR] + x * val0;
                    ptr_ret[r + 1*LDR]  = ptr_ret[r + 1*LDR] + x * val1;
                    ptr_ret[r + 2*LDR]  = ptr_ret[r + 2*LDR] + x * val2;
                    ptr_ret[r + 3*LDR]  = ptr_ret[r + 3*LDR] + x * val3;
                };

                continue;

              lab_seq:
                    
                if (mrd::is_zero(val0) == false)
                {
                    for (Integer i = A_c[k]; i < A_c[k + 1]; ++i)
                    {
                        Integer r       = A_r[i];
                        V1R x           = abs_helper<V1>::eval(A_x[i]);
                        ptr_ret[r+0*LDR]= ptr_ret[r+0*LDR] + x * val0;
                    };
                };
                if (mrd::is_zero(val1) == false)
                {
                    for (Integer i = A_c[k]; i < A_c[k + 1]; ++i)
                    {
                        Integer r       = A_r[i];
                        V1R x           = abs_helper<V1>::eval(A_x[i]);
                        ptr_ret[r+1*LDR]= ptr_ret[r+1*LDR] + x * val1;
                    };
                };
                if (mrd::is_zero(val2) == false)
                {
                    for (Integer i = A_c[k]; i < A_c[k + 1]; ++i)
                    {
                        Integer r       = A_r[i];
                        V1R x           = abs_helper<V1>::eval(A_x[i]);
                        ptr_ret[r+2*LDR]= ptr_ret[r+2*LDR] + x * val2;
                    };
                };
                if (mrd::is_zero(val3) == false)
                {
                    for (Integer i = A_c[k]; i < A_c[k + 1]; ++i)
                    {
                        Integer r       = A_r[i];
                        V1R x           = abs_helper<V1>::eval(A_x[i]);
                        ptr_ret[r+3*LDR]= ptr_ret[r+3*LDR] + x * val3;
                    };
                };
            }

            ptr_X   += LDX * 4;
            ptr_ret += LDR * 4;
        };
            
        for (; j < N; ++j)
        {
            for (Integer k = k0; k < k1; ++k)
            {
                const V1R& val  = abs_helper<V2>::eval(ptr_X[k]);

                if (mrd::is_zero(val))
                    continue;

                for (Integer i = A_c[k]; i < A_c[k + 1]; ++i)
                {
                    Integer r       = A_r[i];
                    V1R x           = abs_helper<V1>::eval(A_x[i]);
                        
                    ptr_ret[r]      = ptr_ret[r] + x * val;
                };
            }

            ptr_X   += X_ld;
            ptr_ret += ret_ld;
        };

        out = matcl::Matrix(ret, false);
    }

    static void eval_trans(matcl::Matrix& out, const M1& A, const M2& X, const Mat_C& C)
    {
        using V1            = typename M1::value_type;
        using V2            = typename M2::value_type;
        using VR            = typename md::real_type<Val_C>::type;
        using V1R           = typename md::real_type<V1>::type;
        using V2R           = typename md::real_type<V2>::type;
        using Mat_ret       = raw::Matrix<VR, struct_dense>;
     
        Integer M               = C.rows();
        Integer N               = C.cols();

        Mat_ret ret(C.get_type(), M, N);

        const sparse_ccs<V1>& d = A.rep();
        const Integer * A_c     = d.ptr_c();
        const Integer * A_r     = d.ptr_r();
        const V1* A_x           = d.ptr_x();

        Integer X_ld            = X.ld();
        Integer C_ld            = C.ld();
        Integer ret_ld          = ret.ld();

        const V2* ptr_X         = X.ptr();
        const Val_C * ptr_C     = C.ptr();
        VR* ptr_ret             = ret.ptr();

        for (Integer j = 0; j < N; ++j)
        {
            for (Integer k = 0; k < M; ++k)
            {
                VR dot          = abs_helper<Val_C>::eval(ptr_C[k]);
                Integer last    = A_c[k + 1];

                for (Integer i = A_c[k]; i < last; ++i)
                    dot = dot + abs_helper<V1>::eval(A_x[i]) * abs_helper<V2>::eval(ptr_X[A_r[i]]);

                ptr_ret[k]      = dot;
            }

            ptr_X   += X_ld;
            ptr_C   += C_ld;
            ptr_ret += ret_ld;
        };

        out = matcl::Matrix(ret, false);
        return;
    }
};

}}};