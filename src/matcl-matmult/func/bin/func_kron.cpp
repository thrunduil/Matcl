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

#include "matcl-matrep/lib_functions/func_matrix.h"
#include "matcl-matrep/lib_functions/matrix_gen.h"
#include "matcl-matrep/lib_functions/func_unary.h"
#include "matcl-matrep/details/extract_type2_switch.h"
#include "matcl-scalar/details/matfunc_helpers.h"
#include "matcl-core/matrix/matrix_traits.h"
#include "matcl-matmult/func/raw/kron.h"
#include "matcl-matmult/func/raw/mmul/mmul.h"
#include "matcl-internals/func/op_helpers.h"
#include "matcl-internals/func/raw_func_op_helpers.h"

#pragma warning( push )
#pragma warning(disable:4127)	// conditional expression is constant

namespace matcl { namespace details
{

namespace md    = matcl::details;
namespace mr    = matcl::raw;
namespace mrd   = matcl::raw::details;
namespace mdyf  = matcl::dynamic::functions;

template<class M1, class M2>
struct function_kron
{
    static void eval(matcl::Matrix& ret, const M1& mat1, const M2& mat2)
    {
        using value_type_1      = typename M1::value_type;
        using value_type_2      = typename M2::value_type;
        using value_type_1_real = typename details::real_type<value_type_1>::type;
        using value_type_2_real = typename details::real_type<value_type_2>::type;

        static const bool is_int = std::is_same<value_type_1,Integer>::value 
                                   && std::is_same<value_type_2,Integer>::value;
        static const bool is_obj = std::is_same<value_type_1,Object>::value 
                                   || std::is_same<value_type_2,Object>::value;

        using ret_ti_type       = typename select_ti_type<is_obj>::type;
        using ret_ti_value      = typename select_ti_type<is_obj>::value_type;

        typename ti::get_ti_type<M1>::type ti_1 = ti::get_ti(mat1);
        typename ti::get_ti_type<M2>::type ti_2 = ti::get_ti(mat2);

        if (mat1.rows() == 1 && mat1.cols() == 1)
        {			
            value_type_1 val(mat1(1,1));

            return mrd::mult_helper_mat_scal<M2,value_type_1>
                ::eval(ret, mat2, val,trans_type::no_trans,trans_type::no_trans);
        };

        if (mat2.rows() == 1 && mat2.cols() == 1)
        {			
            value_type_2 val(mat2(1,1));

            return mrd::mult_helper_mat_scal<M1,value_type_2>
                ::eval(ret, mat1, val,trans_type::no_trans,trans_type::no_trans);
        };

        if (mat1.nnz() == 0)
        {			
            using ret_value_type    = md::unify_types<value_type_1_real,value_type_2_real>::type;
            ret_ti_type ret_ti      = ti::get_return_ti<ret_ti_type>(mdyf::op_mul::eval(), 
                                                            ti_1,ti_2);
            ret = mrd::create_matrix<ret_value_type>(ret_ti,imult_c(mat1.rows(),mat2.rows()),
                                                        imult_c(mat1.cols(),mat2.cols()));
            return;
        };

        if (mat2.nnz() == 0)
        {			
            using ret_value_type    = md::unify_types<value_type_1_real,value_type_2_real>::type;

            ret_ti_type ret_ti      = ti::get_return_ti<ret_ti_type>(mdyf::op_mul::eval(), 
                                                            ti_1,ti_2);
            ret = mrd::create_matrix<ret_value_type>(ret_ti,imult_c(mat1.rows(),mat2.rows()),
                                                    imult_c(mat1.cols(),mat2.cols()));
            return;
        };

        mrd::kron_helper<M1,M2>::eval(ret, mat1, mat2);
        return;
    };
};

struct eval_kron : public extract_type2_switch<void,eval_kron,mr::val_type_corrector_int>
{
    template<class T1, class T2>
    static void eval_mat_mat(const T1& A, const T2& B, matcl::Matrix& ret)
    {
        return function_kron<T1,T2>::eval(ret,A,B);
    };

    template<class T1, class T2>
    static void eval_scal_scal(const T1& A, const T2& B, matcl::Matrix& ret)
    {
        ret = mrd::mmul_helper<T1,T2>::eval(A,B);
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

matcl::Matrix matcl::kron(const Matrix& A0, const Matrix& B0)
{
    //increase refcount
    Matrix A(A0);
    Matrix B(B0);

    matcl::Matrix ret;
    details::eval_kron::make(A,B,ret);
    return ret;
};

matcl::Matrix matcl::kron(Matrix&& A0, const Matrix& B0)
{
    //increase refcount
    Matrix A(std::move(A0));
    Matrix B(B0);

    matcl::Matrix ret;
    details::eval_kron::make(A,B,ret);
    return ret;
};

matcl::Matrix matcl::kron(Matrix&& A0, Matrix&& B0)
{
    //increase refcount
    Matrix A(std::move(A0));
    Matrix B(std::move(B0));

    matcl::Matrix ret;
    details::eval_kron::make(A,B,ret);
    return ret;
};

matcl::Matrix matcl::kron(const Matrix& A0, Matrix&& B0)
{
    //increase refcount
    Matrix A(A0);
    Matrix B(std::move(B0));

    matcl::Matrix ret;
    details::eval_kron::make(A,B,ret);
    return ret;
};

#pragma warning( pop )