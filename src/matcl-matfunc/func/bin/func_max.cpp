/*
 *  This file is a part of Matrix Computation Library (MATCL)
 *
 *  Copyright (c) Pawe� Kowal 2017 - 2021
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
#include "matcl-matfunc/func/raw/bin/raw_func_max.h"
#include "matcl-internals/func/raw_func_op_helpers.h"

#pragma warning( push )
#pragma warning(disable:4127)	// conditional expression is constant

namespace matcl { namespace details
{

namespace md = matcl::details;
namespace mr = matcl::raw;
namespace mrd = matcl::raw::details;

class op_max : public op_no_cases<op_max>
{
    public:
        template<class UT, class TI1, class TI2>
        void eval_special_case_00_impl(matcl::Matrix& ret, Integer m, Integer n, TI1 t1, TI2 t2) const
        {
            using UTR       = md::real_type<UT>::type;
            using ti_type   = matcl::ti::ti_type<UTR>;
            ti_type ti      = matcl::ti::unify_ti<ti_type>(t1,t2);

            ret = mrd::create_matrix<UTR>(ti, m,n);
        };

        template<class M1, class M2>
        static void eval_mat_mat(matcl::Matrix& ret, const M1& mat1, const M2& mat2)
        {
            return mrd::max_helper_mat_mat_inpl<M1,M2>::eval(ret, mat1, mat2);
        };

        template<class M1, class S2>
        static void eval_mat_scal(matcl::Matrix& ret, const M1& mat, const S2& scal)
        {
            if (mrd::is_leq_zero(scal))
	            return mrd::max_helper_mat_scal_inpl<M1,S2>::eval2(ret,mat,scal);
            else
	            return mrd::max_helper_mat_scal_inpl<M1,S2>::eval(ret,mat,scal);
        };

        template<class S1, class M2>
        static void eval_scal_mat(matcl::Matrix& ret, const S1& scal, const M2& mat)
        {
            if (mrd::is_leq_zero(scal))
	            return mrd::max_helper_scal_mat_inpl<S1,M2>::eval2(ret,scal,mat);
            else
	            return mrd::max_helper_scal_mat_inpl<S1,M2>::eval(ret,scal,mat);
        };

        template<class S1, class S2>
        static void eval_scal_scal(matcl::Matrix& ret, const S1& scal1, const S2& scal2)
        {
            ret = mrd::max_helper<S1,S2>::eval(scal1,scal2);
        };
};

struct eval_max : public extract_type2_switch<void,eval_max,mr::val_type_corrector_diag>
{
    template<class T1, class T2>
    static void eval_mat_mat(const T1& A, const T2& B, matcl::Matrix& ret)
    {
        op_max op;
        return md::function_op<T1,T2,op_max,false>::eval(ret,A,B,op);
    };

    template<class T1, class T2>
    static void eval_scal_scal(const T1& A, const T2& B, matcl::Matrix& ret)
    {
        ret = mrd::max_helper<T1,T2>::eval(A,B);
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

matcl::Matrix matcl::max(const Matrix& A0, const Matrix& B0)
{
    Matrix A(A0);
    Matrix B(B0);

    matcl::Matrix ret;
    details::eval_max::make(A,B,ret);
    return ret;
};

matcl::Matrix matcl::max(Matrix&& A0, const Matrix& B0)
{
    Matrix A(std::move(A0));
    Matrix B(B0);

    matcl::Matrix ret;
    details::eval_max::make(A,B,ret);
    return ret;
};

matcl::Matrix matcl::max(Matrix&& A0, Matrix&& B0)
{
    Matrix A(std::move(A0));
    Matrix B(std::move(B0));

    matcl::Matrix ret;
    details::eval_max::make(A,B,ret);
    return ret;
};

matcl::Matrix matcl::max(const Matrix& A0, Matrix&& B0)
{
    Matrix A(A0);
    Matrix B(std::move(B0));

    matcl::Matrix ret;
    details::eval_max::make(A,B,ret);
    return ret;
};

#pragma warning( pop )
