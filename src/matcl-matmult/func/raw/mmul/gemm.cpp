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

#include "matcl-matmult/func/raw/mmul/mmul.h"
#include "matcl-matrep/matrix/matrix.h"
#include "matcl-internals/base/instantiate.h"
#include "matcl-core/matrix/matrix_traits.h"
#include "matcl-internals/base/utils.h"
#include "matcl-scalar/details/scalfunc_helpers.h"
#include "matcl-scalar/details/matfunc_helpers.h"
#include "matcl-matrep/lib_functions/func_binary.h"
#include "matcl-matrep/lib_functions/func_matrix.h"
#include "matcl-internals/func/lapack_utils.h"
#include "matcl-blas-lapack/lapack/lapack.h"
#include "matcl-blas-lapack/level1/level1.h"

#include "matcl-internals/error/error_check_basic.h"
#include "matcl-matrep/matrix/struct_flag_ext.h"
#include "matcl-internals/func/converter.h"
#include "matcl-internals/algs/scatter.h"
#include "matcl-internals/base/sort.h"
#include "matcl-matrep/lib_functions/matrix_gen.h"
#include "matcl-matrep/matrix/struct_flag_ext.h"
#include "matcl-internals/func/raw_manip.h"

#include "mmul_utils.h"

#include "matcl-matmult/func/raw/mmul/mmul_dense_sparse.inl"
#include "matcl-matmult/func/raw/mmul/mmul_sparse_dense.inl"
#include "matcl-matmult/func/raw/mmul/mmul_band_dense.inl"
#include "matcl-matmult/func/raw/mmul/mmul_dense_band.inl"
#include "matcl-matmult/func/raw/mmul/mmul_dense_dense.inl"
#include "matcl-matmult/func/raw/mmul/mmul_sparse_band.inl"
#include "matcl-matmult/func/raw/mmul/mmul_band_sparse.inl"
#include "matcl-matmult/func/raw/mmul/mmul_sparse_sparse.inl"
#include "matcl-matmult/func/raw/mmul/mmul_band_band.inl"

namespace matcl { namespace raw { namespace details
{

namespace mr = matcl::raw;
namespace mrd = matcl::raw::details;
namespace md = matcl::details;

static struct_flag make_gemm_ret_struct(value_struct_class alpha, struct_flag A, struct_flag B, 
                                 trans_type t_A, trans_type t_B, bool is_real_A, bool is_real_B,
                                 bool is_real_alpha, value_struct_class beta, struct_flag C,
                                 bool is_real_C, bool is_square_A, bool is_square_B, bool is_square_C)
{
    struct_flag AB  = predefined_struct_ext::mult_struct(A,B,t_A,t_B,is_real_A,is_real_B, 
                                                         is_square_A, is_square_B, is_square_C);
    struct_flag aAB = predefined_struct_ext::get_scal_mult(AB, alpha);

    if (beta == value_struct_class::vc_zero)
        return aAB;

    struct_flag bC  = predefined_struct_ext::get_scal_mult(C, beta);

    bool is_real_AB = is_real_A && is_real_B && is_real_alpha;    
    struct_flag ret = predefined_struct_ext::plus_struct(aAB, bC, is_real_AB, is_real_C, is_square_C);
    return ret;
};

static struct_flag make_gemm_ret_struct(value_struct_class scal_alpha, struct_flag mat, trans_type t, 
                                 value_struct_class beta, struct_flag C, bool is_real_alpha, bool is_real_mat,
                                 bool is_real_C, bool is_square_C)
{
    struct_flag A   = md::predefined_struct::get_trans(mat,t);
    struct_flag aAB = predefined_struct_ext::get_scal_mult(A, scal_alpha);

    if (beta == value_struct_class::vc_zero)
        return aAB;

    struct_flag bC  = predefined_struct_ext::get_scal_mult(C, beta);

    bool is_real_AB = is_real_alpha && is_real_mat;    
    struct_flag ret = predefined_struct_ext::plus_struct(aAB, bC, is_real_AB, is_real_C, is_square_C);
    return ret;
};

static struct_flag make_gemm_ret_struct(value_struct_class a, value_struct_class beta, struct_flag C)
{
    if (a != value_struct_class::vc_zero)
        return struct_flag();

    if (beta == value_struct_class::vc_zero)
        return predefined_struct_type::diag;

    struct_flag bC  = predefined_struct_ext::get_scal_mult(C, beta);
    return bC;
};

static struct_flag make_gemm_ret_struct(value_struct_class beta, struct_flag C)
{
    if (beta == value_struct_class::vc_zero)
        return predefined_struct_type::diag;

    struct_flag bC  = predefined_struct_ext::get_scal_mult(C, beta);
    return bC;
};

//mult

template<class Val, class M1, class M2>
struct check_need_gemm
{
    using V1    = typename M1::value_type;
    using V2    = typename M2::value_type;
    using VT    = typename md::unify_types2<Val, V1, V2>::type;

    // if Val is not unifier of V1 and V2, then this case should already be removed
    static const bool value = std::is_same<VT, Val>::value;
};

template<class Val, class M1, class M2, class S1, class S2,
        bool need_impl_gemm = check_need_gemm<Val, M1, M2>::value>
struct eval_gemm
{
    using Mat_C = raw::Matrix<Val,struct_dense>;

    static void eval(Mat_C& C, const Val& alpha, const M1& A, const M2& B, 
                     trans_type t_A, trans_type t_B, const Val& beta, 
                     Integer fr, Integer rows, Integer fc, Integer cols)
    {
        bool is_full    = (fr == 0 && fc == 0 && rows == C.rows() && cols == C.cols());

        bool is_re_a    = md::is_float_real_scalar<Val>::value 
                                || std::is_same<Val,Integer>::value;
        bool is_re_C    = is_re_a;
        bool is_re_A    = is_real_matrix(A);
        bool is_re_B    = is_real_matrix(B);

        if (mrd::is_zero(alpha))
        {
            struct_flag ret_str;

            if (is_full == true)
            {
                value_struct_class vc_b = md::predefined_struct::get_value_type(beta,true);
                ret_str                 = make_gemm_ret_struct(vc_b, C.get_struct());
            };

            prepare_gemm_C<Val>::eval(beta, C, fr, rows, fc, cols);
            
            C.set_struct(ret_str);
        }
        else
        {
            struct_flag ret_str;
            if (is_full == true)
            {
                value_struct_class vc_a = md::predefined_struct::get_value_type(alpha,true);
                value_struct_class vc_b = md::predefined_struct::get_value_type(beta,true);
                bool is_sq_A            = (A.rows() == A.cols());
                bool is_sq_B            = (B.rows() == B.cols());
                bool is_sq_C            = (C.rows() == C.cols());
                ret_str                 = make_gemm_ret_struct(vc_a, A.get_struct(), B.get_struct(), t_A, t_B, 
                                                is_re_A, is_re_B, is_re_a, vc_b, C.get_struct(), is_re_C, 
                                                is_sq_A, is_sq_B, is_sq_C);
            };

            eval_mult<Val,M1,M2,S1,S2>::eval_gemm(C, alpha, A, B, t_A, t_B, beta, fr, rows, fc, cols);

            C.set_struct(ret_str);
        };        
    };
};

template<class Val, class M1, class M2, class S1, class S2>
struct eval_gemm<Val,M1,M2,S1,S2,false>
{
    using Mat_C = raw::Matrix<Val,struct_dense>;

    static void eval(Mat_C&, const Val&, const M1&, const M2&, 
                     trans_type, trans_type, const Val&, 
                     Integer, Integer, Integer, Integer)
    {
        //nothig to do
    };
};

//scal
template<class Val, bool is_float = md::is_float_real_scalar<Val>::value>
struct scale_inplace
{
    using Mat_D = raw::Matrix<Val,struct_dense>;
    static void eval(Mat_D& C, const Val& beta, Integer fr, Integer rows, Integer fc, Integer cols)
    {
        prepare_gemm_C<Val>::eval(beta, C, fr, rows, fc, cols);
    };

    static void eval_add(Mat_D& C, const Val& a, const Val& beta, 
                         Integer fr, Integer rows, Integer fc, Integer cols)
    {
        Integer ld_C    = C.ld();
        Val* ptr_C      = C.ptr() + fr + fc * ld_C;        

        level1::apby_test_mat<true, Val, Val, 0, 0, 0>::eval(ptr_C, ld_C, rows, cols, a, beta);
    };
};

template<class Val>
struct scale_inplace<Val,false>
{
    using Mat_D = raw::Matrix<Val,struct_dense>;
    static void eval(Mat_D& C, const Val& beta, Integer fr, Integer rows, Integer fc, Integer cols)
    {
        Integer M       = rows;
        Integer N       = cols;
        Integer ld_C    = C.ld();
        Val* ptr_C      = C.ptr() + fr + fc * ld_C;        

        if (mrd::is_zero(beta))
        {
            Val Z       = md::default_value<Val>(C.get_type());

            for (Integer j = 0; j < N; ++j)
            {
                for (Integer i = 0; i < M; ++i)
                    ptr_C[i]    = Z;

                ptr_C   += ld_C;
            };

            return;
        }
        else
        {
            for (Integer j = 0; j < N; ++j)
            {
                for (Integer i = 0; i < M; ++i)
                    ptr_C[i]    = ptr_C[i] * beta;

                ptr_C   += ld_C;
            };
        };
    };

    static void eval_add(Mat_D& C, const Val& a, const Val& beta,
                         Integer fr, Integer rows, Integer fc, Integer cols)
    {
        if (mrd::is_zero(a))
        {
            return eval(C,beta,fr, rows, fc, cols);
        }

        Integer M       = rows;
        Integer N       = cols;
        Integer ld_C    = C.ld();
        Val* ptr_C      = C.ptr() + fr + fc * ld_C;        

        if (mrd::is_zero(beta))
        {
            for (Integer j = 0; j < N; ++j)
            {
                for (Integer i = 0; i < M; ++i)
                    ptr_C[i]    = a;

                ptr_C   += ld_C;
            };

            return;
        }
        else
        {
            for (Integer j = 0; j < N; ++j)
            {
                for (Integer i = 0; i < M; ++i)
                    ptr_C[i]    = a + ptr_C[i] * beta;

                ptr_C   += ld_C;
            };
        };
    };
};


template<class Val, class Mat>
struct check_need_impl_scal_add
{
    using Val2  = typename Mat::value_type;
    using VT    = typename md::unify_types<Val,Val2>::type;

    // if Val is not unifier of Val2, then this case should already be removed
    static const bool value = std::is_same<VT, Val>::value;
};

template<class Val, class Mat, class struct_type, 
    bool need_impl = check_need_impl_scal_add<Val, Mat>::value>
struct scal_add_struct
{};

template<class Val, class Mat, class struct_type>
struct scal_add_struct<Val, Mat, struct_type, false>
{
    using Mat_D = raw::Matrix<Val,struct_dense>;
    static void eval(Mat_D&, const Val&, const Mat&, trans_type, const Val&,
                     Integer, Integer, Integer, Integer)
    {
        //nothing to do
    };
};

template<class Val, class Mat>
struct scal_add_struct<Val, Mat, struct_dense, true>
{
    using Mat_D     = raw::Matrix<Val,struct_dense>;
    using Val_A     = typename Mat::value_type;

    static void eval(Mat_D& C, const Val& alpha, const Mat& A, trans_type t, const Val& beta,
                     Integer fr, Integer rows, Integer fc, Integer cols)
    {
        if (t != trans_type::no_trans)
        {
            bool is_complex = md::is_complex<Val>::value;
            bool do_conj    = (is_complex == true) && (t == trans_type::conj_trans);

            if (do_conj)
                return eval_trans<true>(C, alpha, A, beta, fr, rows, fc, cols);
            else
                return eval_trans<false>(C, alpha, A, beta, fr, rows, fc, cols);
        }

        Integer M           = rows;
        Integer N           = cols;
        const Val_A* ptr_A  = A.ptr();

        Integer A_ld        = A.ld();
        Integer C_ld        = C.ld();

        Val* ptr_C          = C.ptr() + fr + fc * C_ld;

        level1::axpby_test_mat<true, Val, Val_A, Val, Val, 0, 0, 0>
            ::eval(ptr_C, C_ld, ptr_A, A_ld, M, N, alpha, beta);

        return;
    };

    template<bool Do_conj>
    static void eval_trans(Mat_D& C, const Val& alpha, const Mat& A, const Val& beta,
                     Integer fr, Integer rows, Integer fc, Integer cols)
    {
        Integer M           = rows;
        Integer N           = cols;
        const Val_A* ptr_A  = A.ptr();

        Integer A_ld        = A.ld();
        Integer C_ld        = C.ld();

        Val* ptr_C          = C.ptr() + fr + fc * C_ld;
        const Val_A* ptr_Al = ptr_A;

        if (mrd::is_zero(beta) == true)
        {
            if (mrd::is_one(alpha) == true)
            {
                for (Integer j = 0; j < N; ++j)
                {           
                    ptr_Al          = ptr_A;

                    for (Integer i = 0; i < M; ++i)
                    {
                        ptr_C[i]    = Val(make_conj<Do_conj>::eval(ptr_Al[0]));
                        ptr_Al      += A_ld;
                    };

                    ++ptr_A;
                    ptr_C           += C_ld;
                };
            }
            else
            {
                for (Integer j = 0; j < N; ++j)
                {           
                    ptr_Al          = ptr_A;

                    for (Integer i = 0; i < M; ++i)
                    {
                        ptr_C[i]    = alpha * make_conj<Do_conj>::eval(ptr_Al[0]);
                        ptr_Al      += A_ld;
                    };

                    ++ptr_A;
                    ptr_C           += C_ld;
                };
            };
        }
        else
        {
            if (mrd::is_one(alpha))
            {
                for (Integer j = 0; j < N; ++j)
                {           
                    ptr_Al          = ptr_A;

                    for (Integer i = 0; i < M; ++i)
                    {
                        ptr_C[i]    = make_conj<Do_conj>::eval(ptr_Al[0]) 
                                    + beta * ptr_C[i];
                        ptr_Al      += A_ld;
                    };

                    ++ptr_A;
                    ptr_C           += C_ld;
                };
            }
            else
            {
                for (Integer j = 0; j < N; ++j)
                {           
                    ptr_Al          = ptr_A;

                    for (Integer i = 0; i < M; ++i)
                    {
                        ptr_C[i]    = alpha * make_conj<Do_conj>::eval(ptr_Al[0]) 
                                    + beta * ptr_C[i];
                        ptr_Al      += A_ld;
                    };

                    ++ptr_A;
                    ptr_C           += C_ld;
                };
            };
        };

        return;
    };
};

template<class Val, class Mat>
struct scal_add_struct<Val, Mat, struct_banded, true>
{
    using Mat_D = raw::Matrix<Val,struct_dense>;
    static void eval(Mat_D& C, const Val& alpha, const Mat& A, trans_type t, const Val& beta,
                     Integer fr, Integer rows, Integer fc, Integer cols)
    {
        using Val_A = typename Mat::value_type;

        matcl_assert(t == trans_type::no_trans, "trans should already be removed");

        if (mrd::is_one(beta) == false)
            scale_inplace<Val>::eval(C,beta, fr, rows ,fc, cols);

        Integer N           = cols;

        Integer ld_C        = C.ld();
        Val* ptr_C          = C.ptr() + fr + fc * ld_C;        

        const Val_A* ptr_A  = A.rep_ptr();
        Integer ld_A        = A.ld();

        if (mrd::is_one(alpha) == true)
        {
            for (Integer j = 0; j < N; ++j)
            {
                Integer first_row   = A.first_row(j);
                Integer last_row    = A.last_row(j);
                Integer pos_A       = A.first_elem_pos(j);

                for (Integer l = first_row; l <= last_row; ++l) 
                {
                    ptr_C[l]    = ptr_C[l] + ptr_A[pos_A];

                    ++pos_A;
                };

                ptr_C       += ld_C;
                ptr_A       += ld_A;
            };
        }
        else if (mrd::is_one(-alpha) == true)
        {
            for (Integer j = 0; j < N; ++j)
            {
                Integer first_row   = A.first_row(j);
                Integer last_row    = A.last_row(j);
                Integer pos_A       = A.first_elem_pos(j);

                for (Integer l = first_row; l <= last_row; ++l) 
                {
                    ptr_C[l]    = ptr_C[l] - ptr_A[pos_A];
                    ++pos_A;
                };

                ptr_C       += ld_C;
                ptr_A       += ld_A;
            };
        }
        else
        {
            for (Integer j = 0; j < N; ++j)
            {
                Integer first_row   = A.first_row(j);
                Integer last_row    = A.last_row(j);
                Integer pos_A       = A.first_elem_pos(j);

                for (Integer l = first_row; l <= last_row; ++l) 
                {
                    ptr_C[l]    = ptr_C[l] + alpha * ptr_A[pos_A];

                    ++pos_A;
                };

                ptr_C       += ld_C;
                ptr_A       += ld_A;
            };
        };
    };
};

template<class Val, class Mat>
struct scal_add_struct<Val, Mat, struct_sparse, true>
{    
    using Mat_D = raw::Matrix<Val,struct_dense>;
    static void eval(Mat_D& C, const Val& alpha, const Mat& A, trans_type t, const Val& beta,
                     Integer fr, Integer rows, Integer fc, Integer cols)
    {
        using Val_A = typename Mat::value_type;
        
        matcl_assert(t == trans_type::no_trans, "trans should already be removed");

        Integer N   = cols;

        const sparse_ccs<Val_A>& Ad    = A.rep();

        const Integer* Ad_c = Ad.ptr_c();
        const Integer* Ad_r = Ad.ptr_r();
        const Val_A* Ad_x   = Ad.ptr_x();

        Integer ld_C        = C.ld();
        Val* ptr_C          = C.ptr() + fr + fc * ld_C;        

        if (mrd::is_one(beta) == false)
        {
            scale_inplace<Val>::eval(C, beta, fr, rows, fc, cols);
        };

        if (mrd::is_one(alpha))
        {
            for (Integer j = 0; j < N; ++j)
            {
                Integer f       = Ad_c[j];
                Integer l       = Ad_c[j+1];

                for (Integer k = f; k < l; ++k)
                {
                    Integer pos = Ad_r[k];
                    ptr_C[pos]  = Ad_x[k] + ptr_C[pos];
                };

                ptr_C           += ld_C;
            };
        }
        else if (mrd::is_one(-alpha))
        {
            for (Integer j = 0; j < N; ++j)
            {
                Integer f       = Ad_c[j];
                Integer l       = Ad_c[j+1];

                for (Integer k = f; k < l; ++k)
                {
                    Integer pos = Ad_r[k];
                    ptr_C[pos]  = ptr_C[pos] - Ad_x[k];
                };

                ptr_C           += ld_C;
            };
        }
        else
        {
            for (Integer j = 0; j < N; ++j)
            {
                Integer f       = Ad_c[j];
                Integer l       = Ad_c[j+1];

                for (Integer k = f; k < l; ++k)
                {
                    Integer pos = Ad_r[k];
                    ptr_C[pos]  = alpha * Ad_x[k] + ptr_C[pos];
                };

                ptr_C           += ld_C;
            };
        };

        return;
    };
};

template<class Val, class Mat, class struct_type>
struct scal_add_impl
{
    using Mat_D = raw::Matrix<Val,struct_dense>;
    static void eval(Mat_D& C, const Val& alpha, const Mat& A, trans_type t,
                     const Val& beta, Integer fr, Integer rows, Integer fc, Integer cols)
    {
        using val_type_in   = typename Mat::value_type;
        using ti_ret_type   = typename ti::get_ti_type<Val>::type;

        if (mrd::is_zero(alpha) || A.nnz() == 0)
            return scale_inplace<Val>::eval(C, beta, fr, rows, fc, cols);

        using A_struct      = typename Mat::struct_type;
        bool is_mat_dense   = std::is_same<A_struct, struct_dense>::value;

        if (t != trans_type::no_trans && is_mat_dense == false)
        {
            //we need to take trans now
            matcl::Matrix At(A, false);
            At  = trans(At,t);

            const Mat& A2 = At.impl<Mat>();

            return scal_add_struct<Val, Mat, struct_type>
                    ::eval(C, alpha, A2, trans_type::no_trans, beta, fr, rows, fc, cols);
        };

        return scal_add_struct<Val, Mat, struct_type>
                ::eval(C, alpha, A, t, beta, fr, rows, fc, cols);
    }
};

//----------------------------------------------------------------
//              scal
//----------------------------------------------------------------

static matcl::Matrix convert_scalar(const matcl::Matrix& alpha, const matcl::Matrix& C)
{
    value_code vc_Cr    = matrix_traits::real_value_type(C.get_value_code());
    value_code vc_a     = alpha.get_value_code();
    value_code vc       = matrix_traits::unify_value_types(vc_a, vc_Cr);

    if (vc == vc_a)
        return alpha;

    matcl::Matrix ac    = convert_value(alpha, vc);
    return ac;
};

//C = beta * C + alpha *A
template<class ret_value, class Mat>
struct scal_add 
{
    using Mat_D = raw::Matrix<ret_value,struct_dense>;
    static void eval(Mat_D& C, const ret_value& alpha, const Mat& A, trans_type t,
                     const ret_value& beta, Integer fr, Integer rows, Integer fc, Integer cols)
    {
        bool is_full = (fr == 0 && fc == 0 && rows == C.rows() && cols == C.cols());
        struct_flag ret_str;

        if (is_full == true)
        {
            bool is_re_a            = md::is_float_real_scalar<ret_value>::value 
                                    || std::is_same<ret_value,Integer>::value;
            bool is_re_C            = is_re_a;
            bool is_re_mat          = is_real_matrix(A);
            bool is_square          = C.rows() == C.cols();

            value_struct_class vc_a = md::predefined_struct::get_value_type(alpha,true);
            value_struct_class vc_b = md::predefined_struct::get_value_type(beta,true);
            ret_str                 = make_gemm_ret_struct(vc_a, A.get_struct(), t, vc_b, C.get_struct(),
                                                           is_re_a, is_re_mat, is_re_C, is_square);
        };

        scal_add_impl<ret_value,Mat,typename Mat::struct_type>
                    ::eval(C,alpha,A,t,beta, fr, rows, fc, cols);

        C.set_struct(ret_str);
    };
};

template<class Val>
struct scal_add2_impl
{
    using Mat_D = raw::Matrix<Val,struct_dense>;
    static void eval(Mat_D& C, const Val& a, const Val& beta, Integer fr, Integer rows,
                     Integer fc, Integer cols)
    {
        return scale_inplace<Val>::eval_add(C, a, beta, fr, rows, fc, cols);
    }
};

//C = beta * C + a
template<class ret_value>
struct scal_add2
{
    using Mat_D = raw::Matrix<ret_value,struct_dense>;
    static void eval(Mat_D& C, const ret_value& a, const ret_value& beta,
                     Integer fr, Integer rows, Integer fc, Integer cols)
    {
        bool is_full = (fr == 0 && fc == 0 && rows == C.rows() && cols == C.cols());
        struct_flag ret_str;

        if (is_full == true)
        {
            value_struct_class vc_a = md::predefined_struct::get_value_type(a,true);
            value_struct_class vc_b = md::predefined_struct::get_value_type(beta,true);
            ret_str                 = make_gemm_ret_struct(vc_a, vc_b, C.get_struct());
        }

        scal_add2_impl<ret_value>::eval(C,a,beta, fr, rows, fc, cols);        
        C.set_struct(ret_str);
    };
};

static void check_gemm_args(value_code A, value_code B, value_code alpha, value_code beta, 
                            value_code vC, struct_code sC)
{
    value_code vc1      = matrix_traits::unify_value_types(A, B);
    value_code vc2      = matrix_traits::unify_value_types(alpha, beta);
    value_code vc_exp   = matrix_traits::unify_value_types(vc1, vc2);

    value_code vc_conv  = matrix_traits::unify_value_types(vc_exp, vC);

    if (vc_conv != vC)
        throw error::invalid_gemm_C(vC, sC, vc_exp);

    if (sC != struct_code::struct_dense)
        throw error::invalid_gemm_C(vC, sC, vc_exp);
};

static void check_gemm_args(value_code alpha, value_code beta, 
                            value_code vC, struct_code sC)
{
    value_code vc_exp   = matrix_traits::unify_value_types(alpha, beta);
    value_code vc_conv  = matrix_traits::unify_value_types(vc_exp, vC);

    if (vc_conv != vC)
        throw error::invalid_gemm_C(vC, sC, vc_exp);

    if (sC != struct_code::struct_dense)
        throw error::invalid_gemm_C(vC, sC, vc_exp);
};

template<class T1, class T2>
void gemm_helper_mat_scal<T1, T2>::eval(const T1& mat, const T2& scal, trans_type t_A, trans_type t_B,
                    const matcl::Matrix& alpha, const matcl::Matrix& beta, matcl::Matrix& C,
                    Integer fr, Integer rows, Integer fc, Integer cols)
{
    using V1        = typename T1::value_type;
    using S1        = typename T1::struct_type;

    bool is_scal_A  = (mat.rows() == 1 && mat.cols() == 1);

    if (is_scal_A)
    {
        V1 a        = (t_A == trans_type::conj_trans)? conj(mat(1,1)) : mat(1,1);
        T2 b        = (t_B == trans_type::conj_trans)? conj(scal) : scal;        

        matcl::Matrix alpha_c       = convert_scalar(alpha, C);
        matcl::Matrix scal_alpha    = (alpha_c * a) * b;

        return gemm_helper_scal_scal::eval(scal_alpha, beta, C, fr, rows, fc, cols);
    };

    value_code vc_A = matrix_traits::value_code<V1>::value;
    value_code vc_B = matrix_traits::value_code<T2>::value;
    value_code vc_C = C.get_value_code();

    check_gemm_args(vc_A, vc_B, alpha.get_value_code(), beta.get_value_code(), vc_C, C.get_struct_code());

    if (alpha.is_scalar() == false)
        throw error::scalar_required(alpha.rows(), alpha.cols());

    if (beta.is_scalar() == false)
        throw error::scalar_required(beta.rows(), beta.cols());    

    if (t_A == trans_type::no_trans)
        error::check_eeop(mat.rows(), mat.cols(), rows, cols);
    else
        error::check_eeop(mat.cols(), mat.rows(), rows, cols);

    T2 b                        = (t_B == trans_type::conj_trans)? conj(scal) : scal;
    matcl::Matrix alpha_c       = convert_scalar(alpha, C);
    matcl::Matrix scal_alpha    = alpha_c * b;

    switch (vc_C)
    {
        case value_code::v_integer:
        {
            using VT    = Integer;
            using Mat_C = raw::Matrix<VT,struct_dense>;            
            scal_add<VT, T1>::eval(C.get_impl_unique<Mat_C>(), scal_alpha.get_scalar<VT>(),
                                mat, t_A, beta.get_scalar<VT>(), fr, rows, fc, cols);
            break;
        }
        case value_code::v_real:
        {
            using VT    = Real;
            using Mat_C = raw::Matrix<VT,struct_dense>;            
            scal_add<VT, T1>::eval(C.get_impl_unique<Mat_C>(), scal_alpha.get_scalar<VT>(),
                                mat, t_A, beta.get_scalar<VT>(), fr, rows, fc, cols);
            break;
        }
        case value_code::v_float:
        {
            using VT    = Float;
            using Mat_C = raw::Matrix<VT,struct_dense>;            
            scal_add<VT, T1>::eval(C.get_impl_unique<Mat_C>(), scal_alpha.get_scalar<VT>(),
                                mat, t_A, beta.get_scalar<VT>(), fr, rows, fc, cols);
            break;
        }
        case value_code::v_complex:
        {
            using VT    = Complex;
            using Mat_C = raw::Matrix<VT,struct_dense>;            
            scal_add<VT, T1>::eval(C.get_impl_unique<Mat_C>(), scal_alpha.get_scalar<VT>(),
                                mat, t_A, beta.get_scalar<VT>(), fr, rows, fc, cols);
            break;
        }
        case value_code::v_float_complex:
        {
            using VT    = Float_complex;
            using Mat_C = raw::Matrix<VT,struct_dense>;            
            scal_add<VT, T1>::eval(C.get_impl_unique<Mat_C>(), scal_alpha.get_scalar<VT>(),
                                mat, t_A, beta.get_scalar<VT>(), fr, rows, fc, cols);
            break;
        }
        case value_code::v_object:
        {
            using VT    = Object;
            using Mat_C = raw::Matrix<VT,struct_dense>;            
            scal_add<VT, T1>::eval(C.get_impl_unique<Mat_C>(), scal_alpha.get_scalar<VT>(),
                                mat, t_A, beta.get_scalar<VT>(), fr, rows, fc, cols);
            break;
        }
    }    
};

template<class T1, class T2>
void gemm_helper_scal_mat<T1,T2>::eval(const T1& scal, const T2& mat, trans_type t_A, trans_type t_B,
                    const matcl::Matrix& alpha, const matcl::Matrix& beta, matcl::Matrix& C,
                    Integer fr, Integer rows, Integer fc, Integer cols)
{
    using V2        = typename T2::value_type;
    using S2        = typename T2::struct_type;

    bool is_scal_B  = (mat.rows() == 1 && mat.cols() == 1);

    if (is_scal_B)
    {
        T1 a        = (t_A == trans_type::conj_trans)? conj(scal) : scal;
        V2 b        = (t_B == trans_type::conj_trans)? conj(mat(1,1)) : mat(1,1);

        matcl::Matrix alpha_c       = convert_scalar(alpha, C);
        matcl::Matrix scal_alpha    = (alpha_c * a) * b;

        return gemm_helper_scal_scal::eval(scal_alpha, beta, C, fr, rows, fc, cols);
    };

    value_code vc_A = matrix_traits::value_code<T1>::value;
    value_code vc_B = matrix_traits::value_code<V2>::value;
    value_code vc_C = C.get_value_code();

    check_gemm_args(vc_A, vc_B, alpha.get_value_code(), beta.get_value_code(), vc_C, C.get_struct_code());

    if (alpha.is_scalar() == false)
        throw error::scalar_required(alpha.rows(), alpha.cols());

    if (beta.is_scalar() == false)
        throw error::scalar_required(beta.rows(), beta.cols());

    if (t_B == trans_type::no_trans)
        error::check_eeop(mat.rows(), mat.cols(), rows, cols);
    else
        error::check_eeop(mat.cols(), mat.rows(), rows, cols);

    T1 a                        = (t_A == trans_type::conj_trans)? conj(scal) : scal;
    matcl::Matrix alpha_c       = convert_scalar(alpha, C);
    matcl::Matrix scal_alpha    = alpha_c * a;

    switch (vc_C)
    {
        case value_code::v_integer:
        {
            using VT    = Integer;
            using Mat_C = raw::Matrix<VT,struct_dense>;            
            scal_add<VT, T2>::eval(C.get_impl_unique<Mat_C>(), scal_alpha.get_scalar<VT>(),
                                mat, t_B, beta.get_scalar<VT>(), fr, rows, fc, cols);
            break;
        }
        case value_code::v_real:
        {
            using VT    = Real;
            using Mat_C = raw::Matrix<VT,struct_dense>;            
            scal_add<VT, T2>::eval(C.get_impl_unique<Mat_C>(), scal_alpha.get_scalar<VT>(),
                                mat, t_B, beta.get_scalar<VT>(), fr, rows, fc, cols);
            break;
        }
        case value_code::v_float:
        {
            using VT    = Float;
            using Mat_C = raw::Matrix<VT,struct_dense>;            
            scal_add<VT, T2>::eval(C.get_impl_unique<Mat_C>(), scal_alpha.get_scalar<VT>(),
                                mat, t_B, beta.get_scalar<VT>(), fr, rows, fc, cols);
            break;
        }
        case value_code::v_complex:
        {
            using VT    = Complex;
            using Mat_C = raw::Matrix<VT,struct_dense>;            
            scal_add<VT, T2>::eval(C.get_impl_unique<Mat_C>(), scal_alpha.get_scalar<VT>(),
                                mat, t_B, beta.get_scalar<VT>(), fr, rows, fc, cols);
            break;
        }
        case value_code::v_float_complex:
        {
            using VT    = Float_complex;
            using Mat_C = raw::Matrix<VT,struct_dense>;            
            scal_add<VT, T2>::eval(C.get_impl_unique<Mat_C>(), scal_alpha.get_scalar<VT>(),
                                mat, t_B, beta.get_scalar<VT>(), fr, rows, fc, cols);
            break;
        }
        case value_code::v_object:
        {
            using VT    = Object;
            using Mat_C = raw::Matrix<VT,struct_dense>;            
            scal_add<VT, T2>::eval(C.get_impl_unique<Mat_C>(), scal_alpha.get_scalar<VT>(),
                                mat, t_B, beta.get_scalar<VT>(), fr, rows, fc, cols);
            break;
        }
    }
};

void gemm_helper_scal_scal::eval(const matcl::Matrix& a, const matcl::Matrix& beta, matcl::Matrix& C,
                                 Integer fr, Integer rows, Integer fc, Integer cols)
{
    value_code vc_C = C.get_value_code();

    check_gemm_args(a.get_value_code(), beta.get_value_code(), vc_C, C.get_struct_code());

    if (a.is_scalar() == false)
        throw error::scalar_required(a.rows(), a.cols());

    if (beta.is_scalar() == false)
        throw error::scalar_required(beta.rows(), beta.cols());    

    switch (vc_C)
    {
        case value_code::v_integer:
        {
            using VT    = Integer;
            using Mat_C = raw::Matrix<VT,struct_dense>;            
            scal_add2<VT>::eval(C.get_impl_unique<Mat_C>(), a.get_scalar<VT>(), beta.get_scalar<VT>(),
                                fr, rows, fc, cols);
            break;
        }
        case value_code::v_real:
        {
            using VT    = Real;
            using Mat_C = raw::Matrix<VT,struct_dense>;            
            scal_add2<VT>::eval(C.get_impl_unique<Mat_C>(), a.get_scalar<VT>(), beta.get_scalar<VT>(),
                                fr, rows, fc, cols);
            break;
        }
        case value_code::v_float:
        {
            using VT    = Float;
            using Mat_C = raw::Matrix<VT,struct_dense>;            
            scal_add2<VT>::eval(C.get_impl_unique<Mat_C>(), a.get_scalar<VT>(), beta.get_scalar<VT>(),
                                fr, rows, fc, cols);
            break;
        }
        case value_code::v_complex:
        {
            using VT    = Complex;
            using Mat_C = raw::Matrix<VT,struct_dense>;            
            scal_add2<VT>::eval(C.get_impl_unique<Mat_C>(), a.get_scalar<VT>(), beta.get_scalar<VT>(),
                                fr, rows, fc, cols);
            break;
        }
        case value_code::v_float_complex:
        {
            using VT    = Float_complex;
            using Mat_C = raw::Matrix<VT,struct_dense>;            
            scal_add2<VT>::eval(C.get_impl_unique<Mat_C>(), a.get_scalar<VT>(), beta.get_scalar<VT>(),
                                fr, rows, fc, cols);
            break;
        }
        case value_code::v_object:
        {
            using VT    = Object;
            using Mat_C = raw::Matrix<VT,struct_dense>;            
            scal_add2<VT>::eval(C.get_impl_unique<Mat_C>(), a.get_scalar<VT>(), beta.get_scalar<VT>(),
                                fr, rows, fc, cols);
            break;
        }
    }    
};

template<class M1,class M2>
void gemm_helper<M1,M2>::eval(const M1& A, const M2& B, trans_type t_A, trans_type t_B,
                              const matcl::Matrix& alpha, const matcl::Matrix& beta, matcl::Matrix& C,
                              Integer fr, Integer rows, Integer fc, Integer cols)
{
    using V1            = typename M1::value_type;
    using V2            = typename M2::value_type;
    using S1            = typename M1::struct_type;
    using S2            = typename M2::struct_type;    

    bool is_scal_A      = (A.rows() == 1 && A.cols() == 1);
    bool is_scal_B      = (B.rows() == 1 && B.cols() == 1);

    if (is_scal_A && is_scal_B)
    {
        bool conj_A     = (t_A == trans_type::conj_trans);
        V1 a            = conj_A ? conj(A(1,1)) : A(1,1);

        bool conj_B     = (t_B == trans_type::conj_trans);
        V2 b            = conj_B ? conj(B(1,1)) : B(1,1);

        matcl::Matrix alpha_c       = convert_scalar(alpha, C);
        matcl::Matrix scal_alpha    = (alpha_c * a) * b;

        return gemm_helper_scal_scal::eval(scal_alpha, beta, C, fr, rows, fc, cols);
    };

    if (is_scal_A) 
    {
        V1 a            = A(1,1);
        return gemm_helper_scal_mat<V1,M2>::eval(a, B, t_A, t_B, alpha, beta, C, fr, rows, fc, cols);
    }

    if (is_scal_B) 
    {        
        V2 b            = B(1,1);
        return gemm_helper_mat_scal<M1,V2>::eval(A, b, t_A, t_B, alpha, beta, C, fr, rows, fc, cols);
    };

    value_code vc_A     = matrix_traits::value_code<V1>::value;
    value_code vc_B     = matrix_traits::value_code<V2>::value;
    value_code vc_C     = C.get_value_code();

    check_gemm_args(vc_A, vc_B, alpha.get_value_code(), beta.get_value_code(), vc_C, C.get_struct_code());

    if (alpha.is_scalar() == false)
        throw error::scalar_required(alpha.rows(), alpha.cols());

    if (beta.is_scalar() == false)
        throw error::scalar_required(beta.rows(), beta.cols());    

    error::check_mul(A.rows(), A.cols(), B.rows(), B.cols(), t_A, t_B);

    Integer ret_cols, ret_rows;

    if (t_A == trans_type::no_trans)
        ret_rows        = A.rows();
    else
        ret_rows        = A.cols();

    if (t_B == trans_type::no_trans)
        ret_cols        = B.cols();
    else
        ret_cols        = B.rows();

    error::check_eeop(ret_rows, ret_cols, rows, cols);

    //now gemm product is well formed

    if (A.get_struct().is_id())
    {
        V1 one          = md::one_value<V1>(A.get_type());
        return gemm_helper_scal_mat<V1,M2>::eval(one, B, t_A, t_B, alpha, beta, C, fr, rows, fc, cols);
    };

    if (B.get_struct().is_id())
    {
        V2 one          = md::one_value<V2>(B.get_type());
        return gemm_helper_mat_scal<M1,V2>::eval(A, one, t_A, t_B, alpha, beta, C, fr, rows, fc, cols);
    };

    switch (vc_C)
    {
        case value_code::v_integer:
        {
            using VT    = Integer;
            using Mat_C = raw::Matrix<VT,struct_dense>;            
            eval_gemm<VT, M1, M2, S1, S2>
                ::eval(C.get_impl_unique<Mat_C>(), alpha.get_scalar<VT>(),
                       A, B, t_A, t_B, beta.get_scalar<VT>(), fr, rows, fc, cols);
            break;
        }
        case value_code::v_real:
        {
            using VT    = Real;
            using Mat_C = raw::Matrix<VT,struct_dense>;            
            eval_gemm<VT, M1, M2, S1, S2>
                ::eval(C.get_impl_unique<Mat_C>(), alpha.get_scalar<VT>(),
                       A, B, t_A, t_B, beta.get_scalar<VT>(), fr, rows, fc, cols);
            break;
        }
        case value_code::v_float:
        {
            using VT    = Float;
            using Mat_C = raw::Matrix<VT,struct_dense>;            
            eval_gemm<VT, M1, M2, S1, S2>
                ::eval(C.get_impl_unique<Mat_C>(), alpha.get_scalar<VT>(),
                       A, B, t_A, t_B, beta.get_scalar<VT>(), fr, rows, fc, cols);
            break;
        }
        case value_code::v_complex:
        {
            using VT    = Complex;
            using Mat_C = raw::Matrix<VT,struct_dense>;            
            eval_gemm<VT, M1, M2, S1, S2>
                ::eval(C.get_impl_unique<Mat_C>(), alpha.get_scalar<VT>(),
                       A, B, t_A, t_B, beta.get_scalar<VT>(), fr, rows, fc, cols);
            break;
        }
        case value_code::v_float_complex:
        {
            using VT    = Float_complex;
            using Mat_C = raw::Matrix<VT,struct_dense>;            
            eval_gemm<VT, M1, M2, S1, S2>
                ::eval(C.get_impl_unique<Mat_C>(), alpha.get_scalar<VT>(),
                       A, B, t_A, t_B, beta.get_scalar<VT>(), fr, rows, fc, cols);
            break;
        }
        case value_code::v_object:
        {
            using VT    = Object;
            using Mat_C = raw::Matrix<VT,struct_dense>;            
            eval_gemm<VT, M1, M2, S1, S2>
                ::eval(C.get_impl_unique<Mat_C>(), alpha.get_scalar<VT>(),
                       A, B, t_A, t_B, beta.get_scalar<VT>(), fr, rows, fc, cols);
            break;
        }
    };    
};

}}}

MACRO_INSTANTIATE_BIN_ALL(matcl::raw::details::gemm_helper)

MACRO_INSTANTIATE_SG_2(matcl::raw::details::gemm_helper_scal_mat)
MACRO_INSTANTIATE_SST_2(matcl::raw::details::gemm_helper_scal_mat)

MACRO_INSTANTIATE_GS_2(matcl::raw::details::gemm_helper_mat_scal)
MACRO_INSTANTIATE_STS_2(matcl::raw::details::gemm_helper_mat_scal)
