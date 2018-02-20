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
#include "matcl-matrep/func/bin/op_info.h"
#include "matcl-matrep/lib_functions/matrix_gen.h"
#include "matcl-matrep/lib_functions/func_unary.h"
#include "matcl-matrep/details/extract_type2_switch.h"
#include "matcl-scalar/details/matfunc_helpers.h"
#include "matcl-matrep/func/bin/op_helpers.h"
#include "matcl-core/matrix/matrix_traits.h"
#include "matcl-matrep/func/raw/bin/raw_func_mul.h"
#include "matcl-matrep/func/raw/bin/raw_func_op_helpers.h"
#include "matcl-matrep/func/raw/scal_mul.h"

#pragma warning( push )
#pragma warning(disable:4127)	// conditional expression is constant

namespace matcl { namespace details
{

namespace md = matcl::details;
namespace mr = matcl::raw;
namespace mrd = matcl::raw::details;

template<class M1, class M2>
struct function_op_dmult
{
    static void eval(matcl::Matrix& ret, const M1& mat1, const M2& mat2)
    {
        using value_type_1      = typename M1::value_type;
        using value_type_2      = typename M2::value_type;

        if (mat1.rows() == 1 && mat1.cols() == 1)
        {			
            value_type_1 val(mat1(1,1));

            return mrd::mult_helper_mat_scal<M2,value_type_1>
                ::eval(ret,mat2,val,trans_type::no_trans,trans_type::no_trans);
        };

        if (mat2.rows() == 1 && mat2.cols() == 1)
        {			
            value_type_2 val(mat2(1,1));

            return mrd::mult_helper_mat_scal<M1,value_type_2>
                ::eval(ret,mat1,val,trans_type::no_trans,trans_type::no_trans);
        };

        return mrd::mul_helper_mat_mat_inpl<M1,M2>::eval(ret, mat1, mat2);
    }
};

struct eval_dmult : public extract_type2_switch<void,eval_dmult,mr::val_type_corrector_int>
{
    template<class T1, class T2>
    static void eval_mat_mat(const T1& A, const T2& B, matcl::Matrix& ret)
    {
        return function_op_dmult<T1,T2>::eval(ret,A,B);
    };

    template<class T1, class T2>
    static void eval_scal_scal(const T1& A, const T2& B, matcl::Matrix& ret)
    {
        ret = mrd::elem_mul_helper<T1,T2>::eval(A,B);
    };
    
    template<class T1, class T2>
    static void eval_mat_scal(const T1& mat, const T2& scal, matcl::Matrix& ret)
    {
        return mrd::mult_helper_mat_scal<T1,T2>
            ::eval(ret, mat, scal,trans_type::no_trans,trans_type::no_trans);
    };
    
    template<class T1, class T2>
    static void eval_scal_mat(const T1& scal, const T2& mat, matcl::Matrix& ret)
    {
        return mrd::mult_helper_scal_mat<T1,T2>
            ::eval(ret, scal, mat,trans_type::no_trans,trans_type::no_trans);
    };
};

}}

matcl::Matrix matcl::mul(const Matrix& A0, const Matrix& B0)
{
    //increase refcount
    Matrix A(A0);
    Matrix B(B0);

    matcl::Matrix ret;
    details::eval_dmult::make(A,B,ret);
    return ret;
};

matcl::Matrix matcl::mul(Matrix&& A0, const Matrix& B0)
{
    //increase refcount
    Matrix A(std::move(A0));
    Matrix B(B0);

    matcl::Matrix ret;
    details::eval_dmult::make(A,B,ret);
    return ret;
};

matcl::Matrix matcl::mul(Matrix&& A0, Matrix&& B0)
{
    //increase refcount
    Matrix A(std::move(A0));
    Matrix B(std::move(B0));

    matcl::Matrix ret;
    details::eval_dmult::make(A,B,ret);
    return ret;
};

matcl::Matrix matcl::mul(const Matrix& A0, Matrix&& B0)
{
    //increase refcount
    Matrix A(std::move(A0));
    Matrix B(std::move(B0));

    matcl::Matrix ret;
    details::eval_dmult::make(A,B,ret);
    return ret;
};

#pragma warning( pop )
