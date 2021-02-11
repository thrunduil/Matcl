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
//                      DENSE - SPARSE
//========================================================================
template<class Val_ret, class M1, class M2>
struct eval_mult<Val_ret,M1,M2,struct_dense,struct_sparse> 
{
    using Mat_D = raw::Matrix<Val_ret, struct_dense>;

    static M1 get_diag(const M1& X, trans_type t_A)
    {
        using VT1           = typename M1::value_type;

        VT1 Z1              = md::default_value<VT1>(X.get_type());
        Integer K           = std::min(X.rows(), X.cols());

        M1 D(X.get_type(), Z1, K, 1);

        VT1* ptr_D          = D.ptr();
        const VT1* ptr_X    = X.ptr();
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

    template<class Alpha>
    static void eval_diag_gemm(Mat_D& Y, const Alpha& alpha, const M1& X, const M2& A, 
                               trans_type t_A, trans_type t_B, Integer fr, Integer rows, Integer fc, 
                               Integer cols)
    {
        (void)rows;

        using VTR           = Val_ret;
        using VT1           = typename M1::value_type;
        using VT2           = typename M2::value_type;
        using V12           = typename md::unify_types<VT1, VT2>::type;

        Integer N           = cols;

        M1 D                = get_diag(X, t_A);
        VT1* ptr_D          = D.ptr();
        Integer ld_Y        = Y.ld();
        VTR* ptr_Y          = Y.ptr() + fr + fc * ld_Y;

        Integer K1          = D.length();

        const sparse_ccs<VT2>& d= A.rep();
        const Integer* d_c      = d.ptr_c();
        const Integer* d_r      = d.ptr_r();
        const VT2* d_x          = d.ptr_x();

        if (t_B == trans_type::no_trans)
        {            
            for (Integer j = 0; j < N; ++j)
            {        
                for (Integer k2 = d_c[j]; k2 < d_c[j+1]; ++k2)
                {
                    Integer l       = d_r[k2];
                    const VT2& val  = d_x[k2];

                    if (l >= K1)
                        break;

                    V12 tmp         = ptr_D[l] * val;
                    ptr_Y[l]        = alpha.eval(ptr_Y[l],tmp);
                };

                ptr_Y               += ld_Y;
            };

            return;
        }
        else if (t_B == trans_type::trans)
        {
            for (Integer j = 0; j < K1; ++j)
            {        
                VT1 scal            = ptr_D[j];

                for (Integer k2 = d_c[j]; k2 < d_c[j+1]; ++k2)
                {
                    Integer l       = d_r[k2];
                    const VT2& val  = d_x[k2];

                    V12 tmp         = scal * val;
                    ptr_Y[l*ld_Y]   = alpha.eval(ptr_Y[l*ld_Y], tmp);
                };

                ptr_Y               += 1;
            };

            return;
        }
        else
        {
            for (Integer j = 0; j < K1; ++j)
            {        
                VT1 scal            = ptr_D[j];

                for (Integer k2 = d_c[j]; k2 < d_c[j+1]; ++k2)
                {
                    Integer l       = d_r[k2];
                    const VT2& val  = conj(d_x[k2]);

                    V12 tmp         = scal * val;
                    ptr_Y[l*ld_Y]   = alpha.eval(ptr_Y[l*ld_Y], tmp);
                };

                ptr_Y               += 1;
            };

            return;
        };
    };

    static void eval_diag(matcl::Matrix& ret, const M1& X, const M2& A, trans_type t_A, trans_type t_B)
    {
        using VTR           = Val_ret;
        using VT1           = typename M1::value_type;
        using VT2           = typename M2::value_type;
        using ti_ret_type   = ti::get_ti_type<Val_ret>::type;

        ti_ret_type ret_ti  = ti::get_return_ti<ti_ret_type>(mdyf::op_mul::eval(),
                                                                  ti::get_ti(X), ti::get_ti(A));

        Integer M, K, N;

        if (t_B == trans_type::no_trans)
            N = A.cols();
        else
            N = A.rows();

        if (t_A == trans_type::no_trans)
        {
            M   = X.rows();
            K   = X.cols();
        }
        else
        {
            K   = X.rows();
            M   = X.cols();
        };

        Integer K1  = std::min(M,K);

        M1 D        = get_diag(X, t_A);
        VT1* ptr_D  = D.ptr();

        const sparse_ccs<VT2>& d= A.rep();
        const Integer* d_c      = d.ptr_c();
        const Integer* d_r      = d.ptr_r();
        const VT2* d_x          = d.ptr_x();        

        using SM    = raw::Matrix<VTR, struct_sparse>;

        if (t_B == trans_type::no_trans)
        {            
            SM Y(ret_ti, M, N, A.nnz());

            sparse_ccs<VTR>& dY = Y.rep();
            Integer* Y_c        = dY.ptr_c();
            Integer* Y_r        = dY.ptr_r();
            VTR* Y_x            = dY.ptr_x();        

            Integer nz          = 0;

            for (Integer j = 0; j < N; ++j)
            {
                Y_c[j]          = nz;
        
                for (Integer k2 = d_c[j]; k2 < d_c[j+1]; ++k2)
                {
                    Integer l       = d_r[k2];
                    const VT2& val  = d_x[k2];

                    if (l >= K1)
                        break;

                    VTR tmp = ptr_D[l] * val;

                    if (mrd::is_zero(tmp))
                        continue;

                    Y_r[nz]         = l;
                    Y_x[nz]         = tmp;
                    ++nz;
                };
            };
            Y_c[N]                  = nz;

            bool is_sq_X            = X.rows() == X.cols();
            bool is_sq_A            = A.rows() == A.cols();
            bool is_sq_Y            = Y.rows() == Y.cols();
            Y.set_struct(predefined_struct_ext::mult_struct(X.get_struct(),A.get_struct(),t_A,t_B,
                                is_real_matrix(X),is_real_matrix(A), is_sq_X, is_sq_A, is_sq_Y));

            ret = matcl::Matrix(Y,true);
            return;
        };

        SM At(ret_ti);

        if (t_B == trans_type::trans)
            mrd::manip_trans_reshaper_helper<VT2, VTR>::eval_trans(At, A, A.rows(), K1, M, N);
        else
            mrd::manip_trans_reshaper_helper<VT2, VTR>::eval_ctrans(At, A, A.rows(), K1, M, N);

        sparse_ccs<VTR>& Yd = At.rep();
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
                Integer r   = Y_r[i];
                VTR tmp     = ptr_D[r] * Y_x[i];

                if (mrd::is_zero(tmp))
                    goto lab_remove;

                Y_x[i]      = tmp;
            }

            nz              = Y_c[j + 1];
            continue;

          lab_remove:
            nz              = i;

            for (; i < Y_c[j + 1]; ++i)
            {
                Integer r   = Y_r[i];
                VTR tmp     = ptr_D[r] * Y_x[i];

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
                Integer r   = Y_r[i];
                VTR tmp     = ptr_D[r] * Y_x[i];

                if (mrd::is_zero(tmp))
                    continue;

                Y_x[nz]     = tmp;
                Y_r[nz]     = r;
                ++nz;
            }
        };
        Y_c[N]              = nz;

        bool is_sq_X        = X.rows() == X.cols();
        bool is_sq_A        = A.rows() == A.cols();
        bool is_sq_C        = At.rows() == At.cols();
        At.set_struct(predefined_struct_ext::mult_struct(X.get_struct(),A.get_struct(),t_A,t_B,
                            is_real_matrix(X),is_real_matrix(A), is_sq_X, is_sq_A, is_sq_C));

        ret = matcl::Matrix(At,true);
        return;
    };    

    template<class Alpha>
    static void eval_notrans(const Alpha& alpha, const M1& X, const M2& A, trans_type t_A, trans_type t_B,
                             Mat_D& Y, Integer fr, Integer rows, Integer fc, Integer cols)
    {
        (void)t_B;

        using VTR           = Val_ret;
        using VT1           = typename M1::value_type;
        using VT2           = typename M2::value_type;
        using V12           = typename md::unify_types<VT1, VT2>::type;

        const sparse_ccs<VT2>& d= A.rep();
        const Integer* d_c      = d.ptr_c();
        const Integer* d_r      = d.ptr_r();
        const VT2* d_x          = d.ptr_x();

        Integer X_ld            = X.ld();
        Integer Y_ld            = Y.ld();
        Integer N               = cols;
        Integer M               = rows;
        VTR* ptr_Y              = Y.ptr() + fr + fc * Y_ld;

        if (t_A == trans_type::no_trans)
        {
            for (Integer j = 0; j < N; ++j)
            {
                for (Integer k = d_c[j]; k < d_c[j+1]; ++k)
                {
                    Integer l       = d_r[k];
                    const VT2& val  = d_x[k];

                    if (mrd::is_zero(val))
                        continue;

                    const VT1* ptr_X = X.ptr() + l*X_ld;

                    for (Integer i = 0; i < M; ++i)
                    {
                        V12 tmp     = ptr_X[i] * val;
                        ptr_Y[i]    = alpha.eval(ptr_Y[i],tmp);
                    };
                };

                ptr_Y += Y_ld;
            };
        }
        else
        {
            bool conj_X         = (t_A == trans_type::conj_trans);
            VTR Z               = md::default_value<VTR>(Y.get_type());

            if (conj_X == true)
            {
                for (Integer j = 0; j < N; ++j)
                {
                    if (d_c[j] == d_c[j+1])
                    {
                        ptr_Y           += Y_ld;
                        continue;
                    };

                    const VT1* ptr_X    = X.ptr();

                    for (Integer i = 0; i < M; ++i)
                    {
                        VTR dot         = Z;

                        for (Integer k = d_c[j]; k < d_c[j+1]; ++k)
                        {
                            Integer l       = d_r[k];
                            const VT2& val  = d_x[k];                        

                            dot = dot + conj(ptr_X[l]) * val;
                        };

                        ptr_Y[i]    = ptr_Y[i] + alpha.eval(dot);
                        ptr_X       += X_ld;
                    };

                    ptr_Y           += Y_ld;
                }
            }
            else
            {
                for (Integer j = 0; j < N; ++j)
                {
                    if (d_c[j] == d_c[j+1])
                    {
                        ptr_Y           += Y_ld;
                        continue;
                    };

                    const VT1* ptr_X    = X.ptr();

                    for (Integer i = 0; i < M; ++i)
                    {
                        VTR dot         = Z;

                        for (Integer k = d_c[j]; k < d_c[j+1]; ++k)
                        {
                            Integer l       = d_r[k];
                            const VT2& val  = d_x[k];                        

                            dot = dot + ptr_X[l] * val;
                        };

                        ptr_Y[i]    = ptr_Y[i] + alpha.eval(dot);
                        ptr_X       += X_ld;
                    };

                    ptr_Y           += Y_ld;
                }
            };
        };

        return;
    };
    
    template<class Alpha, bool conj_A>
    static void eval_trans(const Alpha& alpha, const M1& X, const M2& A, trans_type t_A, trans_type t_B,
                           Mat_D& Y, Integer fr, Integer rows, Integer fc, Integer cols)
    {
        if (t_A == trans_type::no_trans)
            return eval_notrans_trans<Alpha,conj_A>(alpha,X,A,t_A,t_B,Y,fr,rows,fc,cols);
        else if (t_A == trans_type::trans)
            return eval_trans_trans<Alpha,conj_A,false>(alpha,X,A,t_A,t_B,Y,fr,rows,fc,cols);
        else
            return eval_trans_trans<Alpha,conj_A,true>(alpha,X,A,t_A,t_B,Y,fr,rows,fc,cols);
    };

    template<class Alpha, bool conj_A>
    static void eval_notrans_trans(const Alpha& alpha, const M1& X, const M2& A, trans_type t_A, trans_type t_B,
                                   Mat_D& Y, Integer fr, Integer rows, Integer fc, Integer cols)
    {
        (void)cols;
        (void)t_B;
        (void)t_A;

        using VTR           = Val_ret;
        using VT1           = typename M1::value_type;
        using VT2           = typename M2::value_type;
        using V12           = typename md::unify_types<VT1, VT2>::type;
        using ti_ret_type   = ti::get_ti_type<Val_ret>::type;
               
        const sparse_ccs<VT2>& d= A.rep();
        const Integer* d_c      = d.ptr_c();
        const Integer* d_r      = d.ptr_r();
        const VT2* d_x          = d.ptr_x();

        Integer X_ld            = X.ld();
        Integer Y_ld            = Y.ld();
        Integer M               = rows;
        Integer K               = A.cols();
        const VT1* ptr_X        = X.ptr();
        VTR* ptr_Y0             = Y.ptr() + fr + fc * Y_ld;

        if (Y_ld == 1)
        {
            VTR* ptr_Y          = ptr_Y0;

            for (Integer j = 0; j < K; ++j)
            {
                for (Integer k = d_c[j]; k < d_c[j+1]; ++k)
                {
                    Integer l       = d_r[k];
                    const VT2& val  = make_conj<conj_A>::eval(d_x[k]);

                    if (mrd::is_zero(val))
                        continue;                    

                    V12 tmp         = ptr_X[0] * val;
                    ptr_Y[l]        = alpha.eval(ptr_Y[l],tmp);
                };

                ptr_X += X_ld;
            };
        }
        else
        {
            for (Integer j = 0; j < K; ++j)
            {
                for (Integer k = d_c[j]; k < d_c[j+1]; ++k)
                {
                    Integer l       = d_r[k];
                    const VT2& val  = make_conj<conj_A>::eval(d_x[k]);

                    if (mrd::is_zero(val))
                        continue;

                    VTR* ptr_Y2     = ptr_Y0 + l*Y_ld;

                    for (Integer i = 0; i < M; ++i)
                    {
                        V12 tmp     = ptr_X[i] * val;
                        ptr_Y2[i]   = alpha.eval(ptr_Y2[i],tmp);
                    };
                };

                ptr_X += X_ld;
            };
        };
        return;
    };

    template<class Alpha, bool conj_A, bool conj_X>
    static void eval_trans_trans(const Alpha& alpha, const M1& X, const M2& A, trans_type t_A, trans_type t_B,
                                 Mat_D& Y, Integer fr, Integer rows, Integer fc, Integer cols)
    {
        (void)cols;
        (void)t_B;
        (void)t_A;

        using VTR           = Val_ret;
        using VT1           = typename M1::value_type;
        using VT2           = typename M2::value_type;
        using V12           = typename md::unify_types<VT1, VT2>::type;
               
        const sparse_ccs<VT2>& d= A.rep();
        const Integer* d_c      = d.ptr_c();
        const Integer* d_r      = d.ptr_r();
        const VT2* d_x          = d.ptr_x();

        Integer X_ld            = X.ld();
        Integer Y_ld            = Y.ld();
        Integer M               = rows;
        Integer K               = A.cols();
        const VT1* ptr_X        = X.ptr();
        VTR* ptr_Y0             = Y.ptr() + fr + fc * Y_ld;

        if (Y_ld == 1)
        {            
            VTR* ptr_Y          = ptr_Y0;

            for (Integer j = 0; j < K; ++j)
            {
                for (Integer k = d_c[j]; k < d_c[j+1]; ++k)
                {
                    Integer l       = d_r[k];
                    const VT2& val  = make_conj<conj_A>::eval(d_x[k]);

                    if (mrd::is_zero(val))
                        continue;

                    V12 tmp         = make_conj<conj_X>::eval(ptr_X[0]) * val;
                    ptr_Y[l]        = alpha.eval(ptr_Y[l],tmp);
                };

                ptr_X += 1;
            };
        }
        else
        {
            for (Integer j = 0; j < K; ++j)
            {
                for (Integer k = d_c[j]; k < d_c[j+1]; ++k)
                {
                    Integer l       = d_r[k];
                    const VT2& val  = make_conj<conj_A>::eval(d_x[k]);

                    if (mrd::is_zero(val))
                        continue;

                    VTR* ptr_Y2     = ptr_Y0 + l*Y_ld;

                    for (Integer i = 0; i < M; ++i)
                    {
                        V12 tmp     = make_conj<conj_X>::eval(ptr_X[i*X_ld]) * val;
                        ptr_Y2[i]   = alpha.eval(ptr_Y2[i],tmp);
                    };
                };

                ptr_X += 1;
            };
        };

        return;
    };

    static void eval(matcl::Matrix& ret, const M1& X, const M2& A, trans_type t_A, trans_type t_B)
    {
        Integer M, K, N;

        if (t_B == trans_type::no_trans)
            N = A.cols();
        else
            N = A.rows();

        if (t_A == trans_type::no_trans)
        {
            M   = X.rows();
            K   = X.cols();
        }
        else
        {
            K   = X.rows();
            M   = X.cols();
        };           
               
        using ti_ret_type   = ti::get_ti_type<Val_ret>::type;

        if (M == 0 || K == 0 || N == 0 || A.nnz() == 0)
        {
            using Val_ret_real  = typename md::real_type<Val_ret>::type;
            using sparse_matrix = Matrix<Val_ret_real,struct_sparse>;

            ti_ret_type ret_ti  = ti::get_return_ti<ti_ret_type>(mdyf::op_mul::eval(),
                                    ti::get_ti(X), ti::get_ti(A));

            sparse_matrix out(ret_ti,M,N);
            ret = matcl::Matrix(out,false);
            return;
        };

        bool is_diag        = raw::is_diag(X);

        if (is_diag)
            return eval_diag(ret, X, A, t_A, t_B);        

        ti_ret_type ret_ti  = ti::get_return_ti<ti_ret_type>(mdyf::op_mul::eval(),
                                    ti::get_ti(X), ti::get_ti(A));

        Val_ret Z           = md::default_value<Val_ret>(ret_ti);

        using ret_type  = raw::Matrix<Val_ret,struct_dense>;
        ret_type Y(ret_ti, Z, M, N);

        Integer fr      = 0;
        Integer rows    = M;
        Integer fc      = 0;
        Integer cols    = N;

        using Alpha     = alpha_one<Val_ret>;
        Alpha alpha;

        if (t_B == trans_type::no_trans)
            eval_notrans(alpha,X,A,t_A,t_B,Y,fr,rows,fc,cols);
        else if (t_B == trans_type::trans)
            eval_trans<Alpha,false>(alpha,X,A,t_A,t_B,Y,fr,rows,fc,cols);
        else
            eval_trans<Alpha,true>(alpha,X,A,t_A,t_B,Y,fr,rows,fc,cols);

        bool is_sq_X    = (X.rows() == X.cols());
        bool is_sq_A    = (A.rows() == A.cols());
        bool is_sq_C    = (Y.rows() == Y.cols());
        Y.set_struct(predefined_struct_ext::mult_struct(X.get_struct(),A.get_struct(),t_A,t_B,
                            is_real_matrix(X),is_real_matrix(A), is_sq_X, is_sq_A, is_sq_C));

        ret = matcl::Matrix(Y,true);
        return;
    };

    static void eval_gemm(Mat_D& Y, const Val_ret& alpha, const M1& X, const M2& A, 
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
                return eval_diag_gemm(Y, al, X, A, t_A, t_B, fr, rows ,fc, cols);

            if (t_B == trans_type::no_trans)
                eval_notrans(al,X,A,t_A,t_B,Y,fr,rows,fc,cols);
            else if (t_B == trans_type::trans)
                eval_trans<Alpha,false>(al,X,A,t_A,t_B,Y,fr,rows,fc,cols);
            else
                eval_trans<Alpha,true>(al,X,A,t_A,t_B,Y,fr,rows,fc,cols);
        }
        else if (mrd::is_one(-alpha) == true)
        {
            using Alpha     = alpha_mone<Val_ret>;
            Alpha al;

            if (is_diag)
                return eval_diag_gemm(Y, al, X, A, t_A, t_B, fr, rows ,fc, cols);

            if (t_B == trans_type::no_trans)
                eval_notrans(al,X,A,t_A,t_B,Y,fr,rows,fc,cols);
            else if (t_B == trans_type::trans)
                eval_trans<Alpha,false>(al,X,A,t_A,t_B,Y,fr,rows,fc,cols);
            else
                eval_trans<Alpha,true>(al,X,A,t_A,t_B,Y,fr,rows,fc,cols);
        }
        else
        {
            using Alpha     = alpha_val<Val_ret>;
            Alpha al(alpha);

            if (is_diag)
                return eval_diag_gemm(Y, al, X, A, t_A, t_B, fr, rows ,fc, cols);

            if (t_B == trans_type::no_trans)
                eval_notrans(al,X,A,t_A,t_B,Y,fr,rows,fc,cols);
            else if (t_B == trans_type::trans)
                eval_trans<Alpha,false>(al,X,A,t_A,t_B,Y,fr,rows,fc,cols);
            else
                eval_trans<Alpha,true>(al,X,A,t_A,t_B,Y,fr,rows,fc,cols);
        };

        return;
    };
};

}}};