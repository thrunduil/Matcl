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

#include "matcl-matrep/lib_functions/func_binary.h"
#include "matcl-matfunc/func/bin/op_info.h"
#include "matcl-matrep/lib_functions/matrix_gen.h"
#include "matcl-matrep/lib_functions/func_unary.h"
#include "matcl-matrep/details/extract_type2_switch.h"
#include "matcl-scalar/details/matfunc_helpers.h"
#include "matcl-internals/func/op_helpers.h"
#include "matcl-core/matrix/matrix_traits.h"
#include "matcl-matfunc/func/raw/bin/raw_func_div.h"
#include "matcl-matfunc/func/raw/bin/raw_func_mul.h"
#include "matcl-internals/func/raw_func_op_helpers.h"
#include "matcl-matfunc/func/raw/scal_mul.h"

#pragma warning( push )
#pragma warning(disable:4127)	// conditional expression is constant

namespace matcl { namespace details
{

namespace md    = matcl::details;
namespace mr    = matcl::raw;
namespace mrd   = matcl::raw::details;
namespace mdyf  = matcl::dynamic::functions;

class op_ddiv : public op_info<op_ddiv>
{
    public:
        template<class Val>
        bool is_special_case_impl(const Val& val, int pos) const				
        { 
            return pos == 2 && mrd::is_one(val); 
        };

        template<class Val>
        void eval_special_case_impl(matcl::Matrix& ret, const Val& v,const Matrix& mat) const
        { 
            return md::to_real_mat<Val>(ret,mat,v,mdyf::op_div::eval()); 
        };

        template<class Val>
        void eval_special_case_impl(matcl::Matrix& ret, const Matrix& mat, const Val& v) const
        { 
            return md::to_real_mat<Val>(ret,mat,v,mdyf::op_div::eval()); 
        };

        template<class UT, class TI1, class TI2>
        void eval_special_case_00_impl(matcl::Matrix& ret, Integer m, Integer n, TI1 t1, TI2 t2) const
        {
            using UTR       = typename md::real_type<UT>::type;
            using UTRF      = md::unify_types<UTR,Float>::type;
            using ti_type   = matcl::ti::ti_type<UTRF>;
            ti_type ti      = ti::get_return_ti<ti_type>(mdyf::op_div::eval(),t1,t2);

            ret = mrd::create_matrix<UTRF>(ti, (UTRF)constants::nan(), m,n);
        };

        template<class M1, class M2>
        static void eval_mat_mat(matcl::Matrix& ret, const M1& mat1, const M2& mat2)
        {
            return mrd::div_helper_mat_mat_inpl<M1,M2>::eval(ret, mat1, mat2);
        };

        template<class M1, class S2>
        static void eval_mat_scal(matcl::Matrix& ret, const M1& mat1, const S2& scal)
        {
            return mrd::div_helper_mat_scal_inpl<M1,S2>::eval(ret, mat1, scal);
        };

        template<class S1, class M2>
        static void eval_scal_mat(matcl::Matrix& ret, const S1& scal, const M2& mat)
        {
            return mrd::div_helper_scal_mat_inpl<S1,M2>::eval(ret, scal, mat);
        };

        template<class S1, class S2>
        static void eval_scal_scal(matcl::Matrix& ret, const S1& scal1, const S2& scal2)
        {
            ret = mrd::div_helper<S1,S2>::eval(scal1,scal2);
        };
};

class op_ddiv_0 : public op_info<op_ddiv_0>
{
    public:
        template<class Val>
        bool is_special_case_impl(const Val& val, int pos) const				
        { 
            return pos == 2 && mrd::is_one(val); 
        };

        template<class Val>
        void eval_special_case_impl(matcl::Matrix& ret, const Val& v,const Matrix& mat) const
        { 
            return md::to_real_mat<Val>(ret,mat,v,mdyf::op_div::eval()); 
        };

        template<class Val>
        void eval_special_case_impl(matcl::Matrix& ret, const Matrix& mat, const Val& v) const
        { 
            return md::to_real_mat<Val>(ret,mat,v,mdyf::op_div::eval()); 
        };

        template<class UT, class TI1, class TI2>
        void eval_special_case_00_impl(matcl::Matrix& ret, Integer m, Integer n, TI1 t1, TI2 t2) const
        {
            using UTR       = typename md::real_type<UT>::type;
            using UTRF      = md::unify_types<UTR,Float>::type;
            using ti_type   = matcl::ti::ti_type<UTRF>;
            ti_type ti      = ti::get_return_ti<ti_type>(mdyf::op_div::eval(),t1,t2);

            ret = mrd::create_matrix<UTRF>(ti, m,n);
        };

        template<class M1, class M2>
        static void eval_mat_mat(matcl::Matrix& ret, const M1& mat1, const M2& mat2)
        {
            return mrd::div_0_helper_mat_mat_inpl<M1,M2>::eval(ret, mat1, mat2);
        };

        template<class M1, class S2>
        static void eval_mat_scal(matcl::Matrix& ret, const M1& mat1, const S2& scal)
        {
            return mrd::div_0_helper_mat_scal_inpl<M1,S2>::eval(ret, mat1, scal);
        };

        template<class S1, class M2>
        static void eval_scal_mat(matcl::Matrix& ret, const S1& scal, const M2& mat)
        {
            return mrd::div_0_helper_scal_mat_inpl<S1,M2>::eval(ret, scal, mat);
        };

        template<class S1, class S2>
        static void eval_scal_scal(matcl::Matrix& ret, const S1& scal1, const S2& scal2)
        {
            ret = mrd::div_0_helper<S1,S2>::eval(scal1,scal2);
        };
};

class op_ddiv_1 : public op_info<op_ddiv_1>
{
    public:
        template<class Val>
        bool is_special_case_impl(const Val& val, int pos) const				
        { 
            return pos == 2 && mrd::is_one(val); 
        };

        template<class Val>
        void eval_special_case_impl(matcl::Matrix& ret, const Val& v,const Matrix& mat) const
        { 
            return md::to_real_mat<Val>(ret,mat,v,mdyf::op_div::eval()); 
        };

        template<class Val>
        void eval_special_case_impl(matcl::Matrix& ret, const Matrix& mat, const Val& v) const
        { 
            return md::to_real_mat<Val>(ret,mat,v,mdyf::op_div::eval()); 
        };

        template<class UT, class TI1, class TI2>
        void eval_special_case_00_impl(matcl::Matrix& ret, Integer m, Integer n, TI1 t1, TI2 t2) const
        {
            using UTR       = typename md::real_type<UT>::type;
            using UTRF      = md::unify_types<UTR,Float>::type;
            using ti_type   = matcl::ti::ti_type<UTRF>;
            ti_type ti      = ti::get_return_ti<ti_type>(mdyf::op_div::eval(),t1,t2);

            ret = mrd::create_matrix<UTRF>(ti, (UTRF)1, m,n);
        };

        template<class M1, class M2>
        static void eval_mat_mat(matcl::Matrix& ret, const M1& mat1, const M2& mat2)
        {
            return mrd::div_1_helper_mat_mat_inpl<M1,M2>::eval(ret, mat1, mat2);
        };

        template<class M1, class S2>
        static void eval_mat_scal(matcl::Matrix& ret, const M1& mat1, const S2& scal)
        {
            return mrd::div_1_helper_mat_scal_inpl<M1,S2>::eval(ret, mat1, scal);
        };

        template<class S1, class M2>
        static void eval_scal_mat(matcl::Matrix& ret, const S1& scal, const M2& mat)
        {
            return mrd::div_1_helper_scal_mat_inpl<S1,M2>::eval(ret, scal, mat);
        };

        template<class S1, class S2>
        static void eval_scal_scal(matcl::Matrix& ret, const S1& scal1, const S2& scal2)
        {
            ret = mrd::div_1_helper<S1,S2>::eval(scal1,scal2);
        };
};

struct eval_ddiv : public extract_type2_switch<void,eval_ddiv,mr::val_type_corrector_int>
{
    template<class T1, class T2>
    static void eval_mat_mat(const T1& A, const T2& B, matcl::Matrix& ret)
    {
        if (B.rows() == 1 && B.cols() == 1)
        {
            using Val_2 = typename T2::value_type;
            Val_2 val(B(1,1));

            return eval_mat_scal<T1, Val_2>(A, val, ret);
        };

        op_ddiv op;
        return md::function_op<T1,T2,op_ddiv,true>::eval(ret,A,B,op);
    };

    template<class T1, class T2>
    static void eval_scal_scal(const T1& A, const T2& B, matcl::Matrix& ret)
    {
        ret = mrd::div_helper<T1,T2>::eval(A,B);
    };

    template<class T1, class T2>
    static void eval_mat_scal(const T1& mat, const T2& scal, matcl::Matrix& ret)
    {
        using TR2       = typename md::real_type<T2>::type;
        using val_type  = typename md::unify_types<TR2,Float>::type;
        using scal_type = typename md::unify_types<T2,Float>::type;
        using V1        = typename T1::value_type;
        using S1        = typename T1::struct_type;
        using V1F       = typename md::unify_types<V1,Float>::type;
        using T1F       = mr::Matrix<V1F, S1>;

        scal_type tmp   = mrd::div_helper<val_type,T2>::eval(val_type(1.),scal);

        return mrd::mult_helper_mat_scal<T1F,scal_type>
                    ::eval(ret, mr::converter<T1F, T1>::eval(mat), tmp, trans_type::no_trans,
                           trans_type::no_trans);
    };

    template<class T1, class T2>
    static void eval_scal_mat(const T1& scal, const T2& mat, matcl::Matrix& ret)
    {
        using FullMatrix = raw::Matrix<T1,struct_dense>;
        FullMatrix m_scal(ti::get_ti(scal),scal,1,1);
        return eval_mat_mat<FullMatrix,T2>(m_scal,mat,ret);
    };
};

struct eval_ddiv_0 : public extract_type2_switch<void,eval_ddiv_0,mr::val_type_corrector_int>
{
    template<class T1, class T2>
    static void eval_mat_mat(const T1& A, const T2& B, matcl::Matrix& ret)
    {
        op_ddiv_0 op;
        return md::function_op<T1,T2,op_ddiv_0,true>::eval(ret,A,B,op);
    };

    template<class T1, class T2>
    static void eval_scal_scal(const T1& A, const T2& B, matcl::Matrix& ret)
    {
        ret = mrd::div_0_helper<T1,T2>::eval(A,B);
    };

    template<class T1, class T2>
    static void eval_mat_scal(const T1& mat, const T2& scal, matcl::Matrix& ret)
    {
        using FullMatrix = raw::Matrix<T2,struct_dense>;
        FullMatrix m_scal(ti::get_ti(scal),scal,1,1);
        return eval_mat_mat<T1,FullMatrix>(mat,m_scal,ret);
    };

    template<class T1, class T2>
    static void eval_scal_mat(const T1& scal, const T2& mat, matcl::Matrix& ret)
    {
        using FullMatrix = raw::Matrix<T1,struct_dense>;
        FullMatrix m_scal(ti::get_ti(scal),scal,1,1);
        return eval_mat_mat<FullMatrix,T2>(m_scal,mat,ret);
    };
};

struct eval_ddiv_1 : public extract_type2_switch<void,eval_ddiv_1,mr::val_type_corrector_int>
{
    template<class T1, class T2>
    static void eval_mat_mat(const T1& A, const T2& B, matcl::Matrix& ret)
    {
        op_ddiv_1 op;
        return md::function_op<T1,T2,op_ddiv_1,true>::eval(ret,A,B,op);
    };

    template<class T1, class T2>
    static void eval_scal_scal(const T1& A, const T2& B, matcl::Matrix& ret)
    {
        ret = mrd::div_1_helper<T1,T2>::eval(A,B);
    };

    template<class T1, class T2>
    static void eval_mat_scal(const T1& mat, const T2& scal, matcl::Matrix& ret)
    {
        using FullMatrix = raw::Matrix<T2,struct_dense>;
        FullMatrix m_scal(ti::get_ti(scal),scal,1,1);
        return eval_mat_mat<T1,FullMatrix>(mat,m_scal,ret);
    };

    template<class T1, class T2>
    static void eval_scal_mat(const T1& scal, const T2& mat, matcl::Matrix& ret)
    {
        using FullMatrix = raw::Matrix<T1,struct_dense>;
        FullMatrix m_scal(ti::get_ti(scal),scal,1,1);
        return eval_mat_mat<FullMatrix,T2>(m_scal,mat,ret);
    };
};

}}

matcl::Matrix matcl::div(const Matrix& A0, const Matrix& B0)
{
    //increase refcount
    Matrix A(A0);
    Matrix B(B0);

    matcl::Matrix ret;
    details::eval_ddiv::make(A,B,ret);
    return ret;
};

matcl::Matrix matcl::div(Matrix&& A0, const Matrix& B0)
{
    //increase refcount
    Matrix A(std::move(A0));
    Matrix B(B0);

    matcl::Matrix ret;
    details::eval_ddiv::make(A,B,ret);
    return ret;
};

matcl::Matrix matcl::div(Matrix&& A0, Matrix&& B0)
{
    //increase refcount
    Matrix A(std::move(A0));
    Matrix B(std::move(B0));

    matcl::Matrix ret;
    details::eval_ddiv::make(A,B,ret);
    return ret;
};

matcl::Matrix matcl::div(const Matrix& A0, Matrix&& B0)
{
    //increase refcount
    Matrix A(A0);
    Matrix B(std::move(B0));

    matcl::Matrix ret;
    details::eval_ddiv::make(A,B,ret);
    return ret;
};

//
matcl::Matrix matcl::div_0(const Matrix& A0, const Matrix& B0)
{
    //increase refcount
    Matrix A(A0);
    Matrix B(B0);

    matcl::Matrix ret;
    details::eval_ddiv_0::make(A,B,ret);
    return ret;
};

matcl::Matrix matcl::div_0(Matrix&& A0, const Matrix& B0)
{
    //increase refcount
    Matrix A(std::move(A0));
    Matrix B(B0);

    matcl::Matrix ret;
    details::eval_ddiv_0::make(A,B,ret);
    return ret;
};

matcl::Matrix matcl::div_0(Matrix&& A0, Matrix&& B0)
{
    //increase refcount
    Matrix A(std::move(A0));
    Matrix B(std::move(B0));

    matcl::Matrix ret;
    details::eval_ddiv_0::make(A,B,ret);
    return ret;
};

matcl::Matrix matcl::div_0(const Matrix& A0, Matrix&& B0)
{
    //increase refcount
    Matrix A(A0);
    Matrix B(std::move(B0));

    matcl::Matrix ret;
    details::eval_ddiv_0::make(A,B,ret);
    return ret;
};

//
matcl::Matrix matcl::div_1(const Matrix& A0, const Matrix& B0)
{
    //increase refcount
    Matrix A(A0);
    Matrix B(B0);

    matcl::Matrix ret;
    details::eval_ddiv_1::make(A,B,ret);
    return ret;
};

matcl::Matrix matcl::div_1(Matrix&& A0, const Matrix& B0)
{
    //increase refcount
    Matrix A(std::move(A0));
    Matrix B(B0);

    matcl::Matrix ret;
    details::eval_ddiv_1::make(A,B,ret);
    return ret;
};

matcl::Matrix matcl::div_1(Matrix&& A0, Matrix&& B0)
{
    //increase refcount
    Matrix A(std::move(A0));
    Matrix B(std::move(B0));

    matcl::Matrix ret;
    details::eval_ddiv_1::make(A,B,ret);
    return ret;
};

matcl::Matrix matcl::div_1(const Matrix& A0, Matrix&& B0)
{
    //increase refcount
    Matrix A(A0);
    Matrix B(std::move(B0));

    matcl::Matrix ret;
    details::eval_ddiv_1::make(A,B,ret);
    return ret;
};

#pragma warning( pop )
