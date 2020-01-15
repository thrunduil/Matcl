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

#include "matcl-matmult/func/raw/mmul/mmul.h"
#include "matcl-internals/base/instantiate.h"
#include "matcl-internals/error/error_check_basic.h"
#include "matcl-matrep/matrix/matrix.h"
#include "matcl-internals/base/utils.h"
#include "matcl-scalar/details/scalfunc_helpers.h"
#include "matcl-scalar/details/matfunc_helpers.h"
#include "matcl-matrep/matrix/struct_flag_ext.h"
#include "matcl-internals/func/lapack_utils.h"
#include "matcl-internals/func/converter.h"
#include "matcl-internals/algs/scatter.h"
#include "matcl-internals/base/sort.h"
#include "matcl-matrep/lib_functions/func_binary.h"
#include "matcl-matrep/lib_functions/func_matrix.h"
#include "matcl-matrep/lib_functions/matrix_gen.h"
#include "matcl-blas-lapack/level1/level1.h"
#include "matcl-internals/func/raw_manip.h"
#include "matcl-matrep/matrix/struct_flag_ext.h"

#include "matcl-blas-lapack/lapack/lapack.h"
#include "mmul_utils.h"
#include "matcl-matmult/func/raw/mmul/mmul_dense_dense.inl"
#include "matcl-matmult/func/raw/mmul/mmul_dense_band.inl"
#include "matcl-matmult/func/raw/mmul/mmul_band_dense.inl"
#include "matcl-matmult/func/raw/mmul/mmul_band_band.inl"
#include "matcl-matmult/func/raw/mmul/mmul_sparse_sparse.inl"
#include "matcl-matmult/func/raw/mmul/mmul_sparse_dense.inl"
#include "matcl-matmult/func/raw/mmul/mmul_dense_sparse.inl"
#include "matcl-matmult/func/raw/mmul/mmul_sparse_band.inl"
#include "matcl-matmult/func/raw/mmul/mmul_band_sparse.inl"

#pragma warning(push)
#pragma warning(disable: 4127) // conditional expression is constant

namespace matcl { namespace raw { namespace details
{

namespace gr = matcl::raw;
namespace mrd = matcl::raw::details;
namespace md = matcl::details;

//scal

template<class Val_ret, class Val_mat, class Mat, class T, class struct_type>
struct scal_struct
{};

template<class Val_ret, class Mat_type, class T>
struct scal_struct<Val_ret, Val_ret, Mat_type, T, struct_banded>
{
    using TI = typename ti::get_ti_type<Val_ret>::type;

    static void eval(matcl::Matrix& ret, const Mat_type &m,const T& a, TI ret_ti)
    {        
        using VT    = Val_ret;

        Integer M   = m.rows();
        Integer N   = m.cols();

        Mat_type res(ret_ti);

        if (m.is_unique() == true && m.get_type() == ret_ti)
        {
            res.assign_to_fresh(m);
        }
        else
        {
            Mat_type res2(ret_ti, M, N, m.first_diag(),m.last_diag());
            res.assign_to_fresh(res2);
        };

        VT* ptr_res         = res.rep_ptr();
        const VT* ptr_m     = m.rep_ptr();

        Integer res_ld      = res.ld();
        Integer m_ld        = m.ld();

        for (Integer i = 0; i < N; ++i)
        {
            Integer fr      = m.first_row(i);
            Integer lr      = m.last_row(i);
            Integer pos     = m.first_elem_pos(i);

            for (Integer k = fr; k <= lr; ++k, ++pos)
                ptr_res[pos] = ptr_m[pos] * a;

            ptr_res     += res_ld;
            ptr_m       += m_ld;
        };

        value_struct_class vt = md::predefined_struct::get_value_type(a,true);
        res.get_struct() = predefined_struct_ext::get_scal_mult(m.get_struct(), vt);
        
        ret = matcl::Matrix(res,true);
        return;
    };
};

template<class Val_ret, class Val_mat, class Mat_type, class T>
struct scal_struct<Val_ret, Val_mat, Mat_type, T, struct_banded>
{
    using TI = typename ti::get_ti_type<Val_ret>::type;

    static void eval(matcl::Matrix& ret, const Mat_type &m,const T& a, TI ret_ti)
    {        
        using VT    = Val_ret;

        Integer M   = m.rows();
        Integer N   = m.cols();

        using Mat_ret   = raw::Matrix<Val_ret, struct_banded>;
        Mat_ret res(ret_ti, M, N, m.first_diag(),m.last_diag());

        VT* ptr_res             = res.rep_ptr();
        const Val_mat* ptr_m    = m.rep_ptr();

        Integer res_ld      = res.ld();
        Integer m_ld        = m.ld();

        for (Integer i = 0; i < N; ++i)
        {
            Integer fr      = m.first_row(i);
            Integer lr      = m.last_row(i);
            Integer pos     = m.first_elem_pos(i);

            for (Integer k = fr; k <= lr; ++k, ++pos)
                ptr_res[pos] = ptr_m[pos] * a;

            ptr_res     += res_ld;
            ptr_m       += m_ld;
        };

        value_struct_class vt = md::predefined_struct::get_value_type(a,true);
        res.get_struct() = predefined_struct_ext::get_scal_mult(m.get_struct(), vt);
        
        ret = matcl::Matrix(res,true);
        return;
    };
};

template<class Val_ret, class Mat_type, class T>
struct scal_struct<Val_ret, Val_ret, Mat_type, T, struct_dense>
{
    using TI = typename ti::get_ti_type<Val_ret>::type;

    static void eval(matcl::Matrix& ret, const Mat_type &m,const T& a, TI ret_ti)
    {
        using VT    = Val_ret;

        Integer M   = m.rows();
        Integer N   = m.cols();

        Mat_type res(ret_ti);

        if (m.is_unique() == true && m.get_type() == ret_ti)
        {
            res.assign_to_fresh(m);
        }
        else
        {
            Mat_type res2(ret_ti, M, N);
            res.assign_to_fresh(res2);
        };

        const VT* ptr_m     = m.ptr();
        VT* ptr_res         = res.ptr();

        Integer m_ld        = m.ld();
        Integer res_ld      = res.ld();

        //form Y = ax
        level1::ax_test_mat<true, VT, VT, T, 0, 0, 0>
            ::eval(ptr_res, res_ld, ptr_m, m_ld, M, N, a);

        using struct_val_type   = value_struct_class;
        struct_val_type vt      = md::predefined_struct::get_value_type(a,true);
        res.get_struct()        = predefined_struct_ext::get_scal_mult(m.get_struct(), vt);

        ret = matcl::Matrix(res,true);
        return;
    };
};

template<class Val_ret, class Val_mat, class Mat_type, class T>
struct scal_struct<Val_ret, Val_mat, Mat_type, T, struct_dense>
{
    using TI = typename ti::get_ti_type<Val_ret>::type;

    static void eval(matcl::Matrix& ret, const Mat_type &m,const T& a, TI ret_ti)
    {
        using VT    = Val_ret;

        Integer M   = m.rows();
        Integer N   = m.cols();

        using Mat_ret   = raw::Matrix<Val_ret, struct_dense>;

        Mat_ret res(ret_ti, M, N);

        const Val_mat* ptr_m= m.ptr();
        VT* ptr_res         = res.ptr();

        Integer m_ld        = m.ld();
        Integer res_ld      = res.ld();

        //form Y = ax
        level1::ax_test_mat<true, VT, Val_mat, T, 0, 0, 0>
            ::eval(ptr_res, res_ld, ptr_m, m_ld, M, N, a);

        using struct_val_type   = value_struct_class;
        struct_val_type vt      = md::predefined_struct::get_value_type(a,true);
        res.get_struct()        = predefined_struct_ext::get_scal_mult(m.get_struct(), vt);

        ret = matcl::Matrix(res,true);
        return;
    };
};

template<class Val_ret, class Mat_type, class T>
struct scal_struct<Val_ret, Val_ret, Mat_type, T, struct_sparse>
{
    using TI = typename ti::get_ti_type<Val_ret>::type;

    static void eval(matcl::Matrix& ret, const Mat_type& m, const T& a, TI ret_ti)
    {
        using VT    = Val_ret;

        Integer M   = m.rows();
        Integer N   = m.cols();
        Integer s   = m.nnz();

        Mat_type res(ret_ti);
        bool inplace = false;

        if (m.is_unique() == true && m.get_type() == ret_ti)
        {
            inplace = true;
            res.assign_to_fresh(m);
        }
        else
        {
            Mat_type res2(ret_ti, M, N, s);
            res.assign_to_fresh(res2);
        };

        const sparse_ccs<VT>& Ad    = m.rep();		
        sparse_ccs<VT>& d           = res.rep();

        const Integer* Ad_c = Ad.ptr_c();
        const Integer* Ad_r = Ad.ptr_r() + Ad.offset();
        const VT* Ad_x      = Ad.ptr_x() + Ad.offset();

        Integer* d_c        = d.ptr_c();
        Integer* d_r        = d.ptr_r() + d.offset();
        VT* d_x             = d.ptr_x() + d.offset();

        //form d_x = a * Ad_x
        level1::ax_test_mat<true, VT, VT, T, 0, 1, 1>
            ::eval(d_x, 1, Ad_x, 1, s, 1, a);

        if (inplace == false)
        {
            for(Integer i = 0; i < s; ++i)
                d_r[i]  = Ad_r[i];

            Integer offset_dif  = - Ad.offset() + d.offset();

            for(Integer i = 0; i <= Ad.cols(); ++i)
                d_c[i]  = Ad_c[i] + offset_dif;
        };

        value_struct_class vt = md::predefined_struct::get_value_type(a,true);
        res.get_struct() = predefined_struct_ext::get_scal_mult(m.get_struct(), vt);

        ret = matcl::Matrix(res,true);
        return;
    };
};

template<class Val_ret, class Val_mat, class Mat_type, class T>
struct scal_struct<Val_ret, Val_mat, Mat_type, T, struct_sparse>
{
    using TI = typename ti::get_ti_type<Val_ret>::type;

    static void eval(matcl::Matrix& ret, const Mat_type& m, const T& a, TI ret_ti)
    {
        using VT    = Val_ret;

        Integer M   = m.rows();
        Integer N   = m.cols();
        Integer s   = m.nnz();

        using Mat_ret   = raw::Matrix<Val_ret, struct_sparse>;
        Mat_ret res(ret_ti, M, N, s);

        const sparse_ccs<Val_mat>& Ad   = m.rep();		
        sparse_ccs<VT>& d               = res.rep();

        const Integer* Ad_c = Ad.ptr_c();
        const Integer* Ad_r = Ad.ptr_r() + Ad.offset();
        const Val_mat* Ad_x = Ad.ptr_x() + Ad.offset();

        Integer* d_c        = d.ptr_c();
        Integer* d_r        = d.ptr_r() + d.offset();
        VT* d_x             = d.ptr_x() + d.offset();

        //form d_x = a * Ad_x
        level1::ax_test_mat<true, VT, Val_mat, T, 0, 1, 1>
            ::eval(d_x, 1, Ad_x, 1, s, 1, a);

        for(Integer i = 0; i < s; ++i)
            d_r[i]  = Ad_r[i];

        Integer offset_dif  = - Ad.offset() + d.offset();

        for(Integer i = 0; i <= Ad.cols(); ++i)
            d_c[i]  = Ad_c[i] + offset_dif;

        value_struct_class vt = md::predefined_struct::get_value_type(a,true);
        res.get_struct() = predefined_struct_ext::get_scal_mult(m.get_struct(), vt);

        ret = matcl::Matrix(res,true);
        return;
    };
};

template<class Val_ret, class Mat, class T, class struct_type>
struct scal_impl
{
    static void eval(matcl::Matrix& ret, const Mat &m, const T& a, trans_type t)
    {        
        using Val           = typename Mat::value_type;
        using ti_ret_type   = typename ti::get_ti_type<Val_ret>::type;

        ti_ret_type ret_ti  = ti::get_return_ti<ti_ret_type>(mdyf::op_mul::eval(), 
                                                            ti::get_ti(m), ti::get_ti(a));

        Integer M, N;
        if (t == trans_type::no_trans)
        {
            M   = m.rows();
            N   = m.cols();
        }
        else
        {
            M   = m.cols();
            N   = m.rows();
        };

        if (mrd::is_zero(a) || m.nnz() == 0)
        {
            using VR            = typename md::real_type<Val_ret>::type;
            using sparse_mat    = Matrix<VR,struct_sparse>;
            sparse_mat out(ret_ti,M,N);
            ret = matcl::Matrix(out,false);
            return;
        };

        if (mrd::is_one(a))
        {
            using VR            = typename md::real_type<Val_ret>::type;
            using val_ret       = typename md::unify_types<VR, Val>::type;
            using ret_mat       = Matrix<val_ret,struct_type>;

            if (t == trans_type::no_trans)
            {
                auto out        = converter<ret_mat,Mat>::eval(m,ret_ti);
                ret             = matcl::Matrix(out,true);
                return;
            }
            else if (t == trans_type::trans)
            {
                ret_mat res(ret_ti);
                mrd::manip_trans_converter_helper<struct_type, Val, val_ret>
                        ::eval_trans(res, m);

                ret             = matcl::Matrix(res,true);
                return;
            }
            else
            {
                ret_mat res(ret_ti);
                mrd::manip_trans_converter_helper<struct_type, Val, val_ret>
                        ::eval_ctrans(res, m);

                ret             = matcl::Matrix(res,true);
                return;
            }
        };

        using real_t = typename matcl::details::real_type<T>::type;
        real_t im(imag_helper<T>::eval(a));
        real_t rea(real_helper<T>::eval(a));

        if (t == trans_type::no_trans)
        {
            if (mrd::is_zero(im))
                return scal_struct<Val_ret, Val, Mat, real_t, struct_type>
                        ::eval(ret, m, rea, ret_ti);
            else
                return scal_struct<Val_ret, Val, Mat, T, struct_type>
                        ::eval(ret, m, a, ret_ti);
        }
        else if (t == trans_type::trans)
        {
            using ret_type  = raw::Matrix<Val_ret,struct_type>;
            ret_type res(ret_ti);
            mrd::manip_trans_converter_helper<struct_type, Val, Val_ret>
                    ::eval_trans(res, m);

            if (mrd::is_zero(im))
                return scal_struct<Val_ret, Val_ret, ret_type, real_t, struct_type>
                        ::eval(ret, res, rea, ret_ti);
            else
                return scal_struct<Val_ret, Val_ret, ret_type, T, struct_type>
                        ::eval(ret, res, a, ret_ti);
        }
        else
        {
            using ret_type  = raw::Matrix<Val_ret,struct_type>;
            ret_type res(ret_ti);
            mrd::manip_trans_converter_helper<struct_type, Val, Val_ret>
                    ::eval_ctrans(res, m);

            if (mrd::is_zero(im))
                return scal_struct<Val_ret, Val_ret, ret_type, real_t, struct_type>
                    ::eval(ret, res, rea, ret_ti);
            else
                return scal_struct<Val_ret, Val_ret, ret_type, T, struct_type>
                    ::eval(ret, res, a, ret_ti);
        }
    };
};

//----------------------------------------------------------------
//              scal
//----------------------------------------------------------------
#pragma warning(push)
#pragma warning(disable: 4127)  // conditional expression is constant

template<class Val_ret, class Mat, class T>
struct scal 
{
    static void eval(matcl::Matrix& ret, const Mat &m,const T& a, trans_type t)
    {
        using Val       = typename Mat::value_type;
        using Struct    = typename Mat::struct_type;

        static const bool is_dense  = std::is_same<Struct, struct_dense>::value;
        bool is_fin                 = (bool)is_finite(a);

        // if a is nan or inf and m is not dense, convert m to dense
        // in order to obtain NaN * m = dense NaN matrix
        if (is_dense == false && is_fin == false )
        {            
            using Mat_F         = raw::Matrix<Val,struct_dense>;
            const Mat_F& m_f    = raw::converter<Mat_F, Mat>::eval(m);

            return scal<Val_ret, Mat_F, T>::eval(ret, m_f, a, t);
        };

        scal_impl<Val_ret,Mat,T,typename Mat::struct_type>::eval(ret,m,a,t);

        if (is_fin == false)
            ret.set_struct(struct_flag());

        return;
    };
};

#pragma warning(pop)

template<class M1,class M2>
void mult_helper<M1,M2>::eval(matcl::Matrix& ret, const M1& A, const M2& B, trans_type t_A, trans_type t_B)
{
    using val_1         = typename M1::value_type;
    using val_2         = typename M2::value_type;
    using str_type_1    = typename M1::struct_type;
    using str_type_2    = typename M2::struct_type;    

    bool is_scal_A      = (A.rows() == 1 && A.cols() == 1);
    bool is_scal_B      = (B.rows() == 1 && B.cols() == 1);

    if (is_scal_A && is_scal_B)
    {
        bool conj_A     = (t_A == trans_type::conj_trans);
        val_1 a         = conj_A ? conj(A(1,1)) : A(1,1);

        bool conj_B     = (t_B == trans_type::conj_trans);
        val_2 b         = conj_B ? conj(B(1,1)) : B(1,1);

        ret             = a * b;
        return;
    };

    if (is_scal_A) 
    {
        using VA    = typename M1::value_type;
        bool conj_A = (t_A == trans_type::conj_trans);
        VA a        = conj_A ? conj(A(1,1)) : A(1,1);

        return scal<Val_ret,M2,VA>::eval(ret, B, a, t_B);
    }

    if (is_scal_B) 
    {        
        using VB    = typename M2::value_type;
        bool conj_B = (t_B == trans_type::conj_trans);
        VB b        = conj_B ? conj(B(1,1)) : B(1,1);

        return scal<Val_ret,M1,VB>::eval(ret, A, b, t_A);
    };

    error::check_mul(A.rows(), A.cols(), B.rows(), B.cols(), t_A, t_B);

    if (A.get_struct().is_id())
    {
        val_1 one       = md::one_value<val_1>(A.get_type());
        bool is_conj    = (t_A == trans_type::conj_trans);
        val_1 a         = is_conj? conj(one) : one;
        return scal<Val_ret,M2,val_1>::eval(ret, B, a, t_B);
    };

    if (B.get_struct().is_id())
    {
        val_2 one       = md::one_value<val_2>(B.get_type());
        bool is_conj    = (t_B == trans_type::conj_trans);
        val_2 b         = is_conj? conj(one) : one;
        return scal<Val_ret,M1,val_2>::eval(ret, A, b, t_A);
    };

    return eval_mult<Val_ret,M1,M2,str_type_1,str_type_2>::eval(ret,A,B,t_A,t_B);
};

template<class Val_C, class M1, class M2>
struct eval_mult_abs_0
{
    using V1        = typename M1::value_type;
    using V2        = typename M2::value_type;
    using S1        = typename M1::struct_type;
    using S2        = typename M2::struct_type;

    using Val_ret_0 = typename md::unify_types<V1, V2>::type;
    using Val_ret_1 = typename md::unify_types<Val_ret_0, Val_C>::type;
    using VR        = typename md::real_type<Val_ret_1>::type;

    using V1_conv   = typename md::unify_types<V1, VR>::type;
    using V2_conv   = typename md::unify_types<V2, VR>::type;

    using M1_conv   = raw::Matrix<V1_conv, S1>;
    using M2_conv   = raw::Matrix<V2_conv, struct_dense>;
    using MC_conv   = raw::Matrix<Val_ret_1, struct_dense>;

    static void eval(matcl::Matrix& ret, const M1& A, const M2& X, 
                     trans_type t_A, const matcl::Matrix& C0)
    {
        const M1_conv& Ac   = raw::converter<M1_conv,M1>::eval(A);
        const M2_conv& Xc   = raw::converter<M2_conv,M2>::eval(X);
        matcl::Matrix C     = convert(C0, MC_conv::matrix_code);
        const MC_conv& Cc   = C.get_impl<MC_conv>();

        return eval_mult_abs<Val_ret_1, M1_conv, M2_conv, S1>::eval(ret, Ac, Xc, t_A, Cc);
    }
};

template<class Val_C, class M1, class S2>
struct eval_mult_abs_0_scal
{
    using V1        = typename M1::value_type;
    using V2        = S2;
    using S1        = typename M1::struct_type;

    using Val_ret_0 = typename md::unify_types<V1, V2>::type;
    using Val_ret_1 = typename md::unify_types<Val_ret_0, Val_C>::type;
    using VR        = typename md::real_type<Val_ret_1>::type;

    using V1_conv   = typename md::unify_types<V1, VR>::type;
    using V2_conv   = typename md::unify_types<V2, VR>::type;

    using M1_conv   = raw::Matrix<V1_conv, S1>;
    using M2_conv   = V2_conv;
    using M2R       = typename md::real_type<V2_conv>::type;
    using MC_conv   = raw::Matrix<Val_ret_1, struct_dense>;

    static void eval(matcl::Matrix& ret, const M1& mat, const S2& scal, 
                     trans_type t_A, const matcl::Matrix& C0)
    {
        const M1_conv& Ac   = raw::converter<M1_conv,M1>::eval(mat);
        const M2_conv& Xc   = raw::converter<M2_conv,S2>::eval(scal);
        matcl::Matrix C     = convert(C0, MC_conv::matrix_code);
        const MC_conv& Cc   = C.get_impl<MC_conv>();

        M2R Xc_re           = abs_helper<M2_conv>::eval(Xc);

        return eval_mult_abs<Val_ret_1, M1_conv, M2R, S1>::eval_scal(ret, Ac, Xc_re, t_A, Cc);
    }
};

template<class V>
struct eval_mult_abs_scal2
{
    using VR    = typename md::real_type<V>::type;
    using Mat   = raw::Matrix<V,struct_dense>;
    using Mat_R = raw::Matrix<VR,struct_dense>;

    static void eval(matcl::Matrix& out, const VR& scal, const matcl::Matrix& C0)
    {    
        matcl::Matrix C = matcl::convert(C0, Mat::matrix_code);
        const Mat& Cc   = C.get_impl<Mat>();

        Integer M       = Cc.rows();
        Integer N       = Cc.cols();

        Mat_R ret(Cc.get_type(), M, N);

        if (M == 0 || N == 0)
        {            
            out = matcl::Matrix(ret, false);
            return;
        }        

        const V* ptr_C  = Cc.ptr();
        VR* ptr_ret     = ret.ptr();

        Integer C_ld    = Cc.ld();
        Integer ret_ld  = ret.ld();

        level1::apx_abs_test_mat<true, VR, V, 0, 0, 0>
            ::eval(ptr_ret, ret_ld, ptr_C, C_ld, M, N, scal);

        out = matcl::Matrix(ret,false);
        return;
    };
};

template<class M1,class M2>
void mult_helper<M1,M2>::eval_abs(matcl::Matrix& ret, const M1& A, const M2& X, 
                                           trans_type t_A, const matcl::Matrix& C)
{
    using val_1         = typename M1::value_type;
    using val_2         = typename M2::value_type;
    using str_type_1    = typename M1::struct_type;
    using str_type_2    = typename M2::struct_type;  

    if (A.rows() == 1 && A.cols() == 1)
        return mult_helper_scal_mat<val_1,M2>::eval_abs(ret, A(1,1), X, C);

    if (X.rows() == 1 && X.cols() == 1)
        return mult_helper_mat_scal<M1, val_2>::eval_abs(ret, A, X(1,1), t_A, C);

    error::check_mul(A.rows(), A.cols(), X.rows(), X.cols(), t_A, trans_type::no_trans);

    Integer M, N;

    if (t_A == trans_type::no_trans)
    {
        M               = A.rows();
        N               = X.cols();
    }
    else
    {
        M               = A.cols();
        N               = X.cols();
    };

    error::check_eeop(M, N, C.rows(), C.cols());  

    if (A.get_struct().is_id())
    {
        val_1 one       = md::one_value<val_1>(A.get_type());

        return mult_helper_scal_mat<val_1,M2>::eval_abs(ret, one, X, C);
    };

    switch (C.get_value_code())
    {
        case value_code::v_integer:
        {
            using V = Integer;
            return eval_mult_abs_0<V,M1,M2>::eval(ret,A,X,t_A,C);
        }    
        case value_code::v_float:
        {
            using V = Float;
            return eval_mult_abs_0<V,M1,M2>::eval(ret,A,X,t_A,C);
        }    
        case value_code::v_real:
        {
            using V = Real;
            return eval_mult_abs_0<V,M1,M2>::eval(ret,A,X,t_A,C);
        }    
        case value_code::v_complex:
        {
            using V = Complex;
            return eval_mult_abs_0<V,M1,M2>::eval(ret,A,X,t_A,C);
        }    
        case value_code::v_float_complex:
        {
            using V = Float_complex;
            return eval_mult_abs_0<V,M1,M2>::eval(ret,A,X,t_A,C);
        } 
        case value_code::v_object:
        {
            using V = Object;
            return eval_mult_abs_0<V,M1,M2>::eval(ret,A,X,t_A,C);
        }    
        default:
        {
            matcl_assert(0,"invalid value_code");
            throw error::error_general("invalid value_code");
        }
    }
};

void mult_helper_scal_scal::eval_abs(matcl::Matrix& ret, const matcl::Matrix& A, 
                                     const matcl::Matrix& B, const matcl::Matrix& C)
{
    value_code vc_A     = A.get_value_code();
    value_code vc_B     = B.get_value_code();
    value_code vc_C     = C.get_value_code();
    
    value_code vc       = matrix_traits::unify_value_types(vc_A, vc_B);
    vc                  = matrix_traits::unify_value_types(vc, vc_C);
    value_code vcr      = matrix_traits::real_value_type(vc);

    value_code vc_A2    = matrix_traits::unify_value_types(vc_A, vcr);
    value_code vc_B2    = matrix_traits::unify_value_types(vc_B, vcr);
    value_code vc_C2    = matrix_traits::unify_value_types(vc_C, vcr);

    matcl::Matrix Ac    = (vc_A != vc_A2) ? convert_value(A, vc_A2) : A;
    matcl::Matrix Bc    = (vc_B != vc_B2) ? convert_value(B, vc_B2) : B;

    matcl::Matrix scal  = abs(Ac) * abs(Bc);

    if (C.is_scalar())
    {
        matcl::Matrix Cc= (vc_C != vc_C2) ? convert_value(C, vc_C2) : C;
        ret             = scal + abs(Cc);
        return;
    };

    if (C.is_empty() == true)
    {
        ret             = convert_value(C, vcr);
        return;
    }

    switch (vc_C2)
    {
        case value_code::v_integer:
        {
            using V     = Integer;

            using VR    = md::real_type<V>::type;
            VR s        = scal.get_scalar<VR>();
            return eval_mult_abs_scal2<V>::eval(ret,s,C);
        }    
        case value_code::v_float:
        {
            using V     = Float;

            using VR    = md::real_type<V>::type;
            VR s        = scal.get_scalar<VR>();
            return eval_mult_abs_scal2<V>::eval(ret,s,C);
        }    
        case value_code::v_real:
        {
            using V     = Real;

            using VR    = md::real_type<V>::type;
            VR s        = scal.get_scalar<VR>();
            return eval_mult_abs_scal2<V>::eval(ret,s,C);
        }    
        case value_code::v_complex:
        {
            using V     = Complex;

            using VR    = md::real_type<V>::type;
            VR s        = scal.get_scalar<VR>();
            return eval_mult_abs_scal2<V>::eval(ret,s,C);
        }    
        case value_code::v_float_complex:
        {
            using V     = Float_complex;

            using VR    = md::real_type<V>::type;
            VR s        = scal.get_scalar<VR>();
            return eval_mult_abs_scal2<V>::eval(ret,s,C);
        }    
        case value_code::v_object:
        {
            using V     = Object;

            using VR    = md::real_type<V>::type;
            VR s        = scal.get_scalar<VR>();
            return eval_mult_abs_scal2<V>::eval(ret,s,C);
        }    
        default:
        {
            matcl_assert(0,"invalid value_code");
            throw error::error_general("invalid value_code");
        }
    }
};

template<class T1, class T2>
void mult_helper_mat_scal<T1,T2>::eval_abs(matcl::Matrix& ret, const T1& mat, const T2& s, 
                                  trans_type t_A, const matcl::Matrix& C)
{
    using V1        = typename T1::value_type;
    using S1        = typename T1::struct_type;

    bool is_scal_A  = (mat.rows() == 1 && mat.cols() == 1);

    if (is_scal_A)
    {
        mrd::mult_helper_scal_scal::eval_abs(ret, mat(1,1), s, C);
        return;
    };

    Integer M, N;

    if (t_A == trans_type::no_trans)
    {
        M               = mat.rows();
        N               = mat.cols();
    }
    else
    {
        M               = mat.cols();
        N               = mat.rows();
    };

    error::check_eeop(M, N, C.rows(), C.cols());

    value_code vc_C     = C.get_value_code();

    if (C.is_empty() == true)
    {
        using TS        = typename md::unify_types<V1, T2>::type;
        using TSR       = typename md::real_type<TS>::type;
        value_code vc_s = matrix_traits::value_code<TSR>::value;
        value_code vc   = matrix_traits::unify_value_types(vc_C, vc_s);
        value_code vc_r = matrix_traits::real_value_type(vc);

        ret             = convert_value(C, vc_r);
        return;
    }

    switch (C.get_value_code())
    {
        case value_code::v_integer:
        {
            using V = Integer;
            return eval_mult_abs_0_scal<V,T1,T2>::eval(ret,mat,s,t_A,C);
        }    
        case value_code::v_float:
        {
            using V = Float;
            return eval_mult_abs_0_scal<V,T1,T2>::eval(ret,mat,s,t_A,C);
        }    
        case value_code::v_real:
        {
            using V = Real;
            return eval_mult_abs_0_scal<V,T1,T2>::eval(ret,mat,s,t_A,C);
        }    
        case value_code::v_complex:
        {
            using V = Complex;
            return eval_mult_abs_0_scal<V,T1,T2>::eval(ret,mat,s,t_A,C);
        }    
        case value_code::v_float_complex:
        {
            using V = Float_complex;
            return eval_mult_abs_0_scal<V,T1,T2>::eval(ret,mat,s,t_A,C);
        } 
        case value_code::v_object:
        {
            using V = Object;
            return eval_mult_abs_0_scal<V,T1,T2>::eval(ret,mat,s,t_A,C);
        }    
        default:
        {
            matcl_assert(0,"invalid value_code");
            throw error::error_general("invalid value_code");
        }
    }
};

template<class T1, class T2>
void mult_helper_mat_scal<T1,T2>::eval(matcl::Matrix& ret, const T1& mat, const T2& s, 
                                  trans_type t_A, trans_type t_B)
{
    using V1        = typename T1::value_type;
    using S1        = typename T1::struct_type;
    using VT        = typename md::unify_types<V1,T2>::type;

    T2 b            = (t_B == trans_type::conj_trans)? conj(s) : s;

    bool is_scal_A  = (mat.rows() == 1 && mat.cols() == 1);

    if (is_scal_A)
    {
        bool conj_A = (t_A == trans_type::conj_trans);
        V1 a        = conj_A ? conj(mat(1,1)) : mat(1,1);

        ret         = a * b;
        return;
    };

    return scal<VT,T1,T2>::eval(ret, mat, b, t_A);
};

template<class T1, class T2>
void mult_helper_scal_mat<T1,T2>::eval_abs(matcl::Matrix& ret, const T1& scal, const T2& mat, 
                                  const matcl::Matrix& C)
{
    return mult_helper_mat_scal<T2,T1>::eval_abs(ret, mat, scal, trans_type::no_trans, C);
};

template<class T1, class T2>
void mult_helper_scal_mat<T1,T2>::eval(matcl::Matrix& ret, const T1& s, const T2& mat, 
                                       trans_type t_A, trans_type t_B)
{
    using V2        = typename T2::value_type;
    using S2        = typename T2::struct_type;
    using VT        = typename md::unify_types<T1,V2>::type;

    T1 a            = (t_A == trans_type::conj_trans)? conj(s) : s;

    bool is_scal_B  = (mat.rows() == 1 && mat.cols() == 1);

    if (is_scal_B)
    {
        bool conj_B = (t_B == trans_type::conj_trans);
        V2 b        = conj_B ? conj(mat(1,1)) : mat(1,1);

        ret         = a * b;
        return;
    };

    return scal<VT,T2,T1>::eval(ret, mat, a, t_B);
};


}}}

MACRO_INSTANTIATE_BIN_ALL(matcl::raw::details::mult_helper)

MACRO_INSTANTIATE_SG_2(matcl::raw::details::mult_helper_scal_mat)
MACRO_INSTANTIATE_SST_2(matcl::raw::details::mult_helper_scal_mat)

MACRO_INSTANTIATE_GS_2(matcl::raw::details::mult_helper_mat_scal)
MACRO_INSTANTIATE_STS_2(matcl::raw::details::mult_helper_mat_scal)

#pragma warning(pop)