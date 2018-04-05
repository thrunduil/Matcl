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
#include "matcl-matfunc/func/raw/bin/raw_func_plus.h"

#pragma warning( push )
#pragma warning(disable:4127)	// conditional expression is constant

namespace matcl { namespace details
{

namespace md    = matcl::details;
namespace mr    = matcl::raw;
namespace mrd   = matcl::raw::details;
namespace mdyf  = matcl::dynamic::functions;

template<class M1, class M2>
struct function_op_plus
{
    static void eval(matcl::Matrix& ret, const M1& A, const M2& B)
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

        typename ti::get_ti_type<M1>::type ti_1 = ti::get_ti(A);
        typename ti::get_ti_type<M2>::type ti_2 = ti::get_ti(B);

        if (A.rows() == 1 && A.cols() == 1)
        {			
            value_type_1 val(A(1,1));

            if (mrd::is_zero(val))
            {
                ret_ti_type ret_ti  = ti::get_return_ti<ret_ti_type>(mdyf::op_plus::eval(), 
                                                        ti_1,ti_2);
                ret = correct_int_val<is_int,is_obj>::eval(ret_ti,Matrix(B,false));
                return;
            };

            value_type_1_real val_im(mrd::imag_helper<value_type_1>::eval(val));

            if (mrd::is_zero(val_im))
            {
                value_type_1_real val_re(mrd::real_helper<value_type_1>::eval(val));
                return mrd::plus_helper_mat_scal_inpl<M2,value_type_1_real>::eval(ret, B, val_re);
            }
            else
            {
                return mrd::plus_helper_mat_scal_inpl<M2,value_type_1>::eval(ret, B, val);
            };
        };

        if (B.rows() == 1 && B.cols() == 1)
        {			
            value_type_2 val(B(1,1));

            if (mrd::is_zero(val))
            {
                ret_ti_type ret_ti  = ti::get_return_ti<ret_ti_type>(mdyf::op_plus::eval(), 
                                                        ti_1,ti_2);
                ret = correct_int_val<is_int,is_obj>::eval(ret_ti,Matrix(A,false));
                return;
            };

            value_type_2_real val_im(mrd::imag_helper<value_type_2>::eval(val));

            if (mrd::is_zero(val_im))
            {
                value_type_2_real val_re(mrd::real_helper<value_type_2>::eval(val));

                return mrd::plus_helper_mat_scal_inpl<M1,value_type_2_real>::eval(ret, A, val_re);
            }
            else
            {
                return mrd::plus_helper_mat_scal_inpl<M1,value_type_2>::eval(ret, A, val);
            };
        };

        Integer m1 = A.rows();
        Integer n1 = A.cols();
        Integer m2 = B.rows();
        Integer n2 = B.cols();

        error::check_eeop(m1,n1,m2,n2);

        if (A.nnz() == 0)
        {
            ret_ti_type ret_ti  = ti::get_return_ti<ret_ti_type>(mdyf::op_plus::eval(), 
                                                        ti_1,ti_2);
            ret = correct_int_val<is_int,is_obj>::eval(ret_ti,Matrix(B,false));
            return;
        }
        else if (B.nnz() == 0)
        {
            ret_ti_type ret_ti  = ti::get_return_ti<ret_ti_type>(mdyf::op_plus::eval(), 
                                                        ti_1,ti_2);
            ret = correct_int_val<is_int,is_obj>::eval(ret_ti,Matrix(A,false));
            return;
        };

        return mrd::plus_helper_mat_mat_inpl<M1,M2>::eval(ret, A,B);
    };
};

struct eval_plus : public extract_type2_switch<void,eval_plus,mr::val_type_corrector_int>
{
    template<class T1, class T2>
    static void eval_mat_mat(const T1& A, const T2& B, matcl::Matrix& ret)
    {
        return function_op_plus<T1,T2>::eval(ret,A,B);
    };
   
    template<class T1, class T2>
    static void eval_mat_scal(const T1& A, const T2& B, matcl::Matrix& ret)
    {
        return mrd::plus_helper_mat_scal_inpl<T1,T2>::eval(ret, A, B);
    };
    
    template<class T1, class T2>
    static void eval_scal_mat(const T1& A, const T2& B, matcl::Matrix& ret)
    {
        return mrd::plus_helper_mat_scal_inpl<T2,T1>::eval(ret, B, A);
    };
    
    template<class T1, class T2>
    static void eval_scal_scal(const T1& A, const T2& B, matcl::Matrix& ret)
    {
        ret = mrd::plus_helper<T1,T2>::eval(A,B);
    };
};

}};

matcl::Matrix matcl::plus(const Matrix& A0, const Matrix& B0)
{
    Matrix A(A0);
    Matrix B(B0);

    matcl::Matrix ret;
    details::eval_plus::make(A,B,ret);
    return ret;
};

matcl::Matrix matcl::plus(Matrix&& A0, const Matrix& B0)
{
    Matrix A(std::move(A0));
    Matrix B(B0);

    matcl::Matrix ret;
    details::eval_plus::make(A,B,ret);
    return ret;
};

matcl::Matrix matcl::plus(Matrix&& A0, Matrix&& B0)
{
    Matrix A(std::move(A0));
    Matrix B(std::move(B0));

    matcl::Matrix ret;
    details::eval_plus::make(A,B,ret);
    return ret;
};

matcl::Matrix matcl::plus(const Matrix& A0, Matrix&& B0)
{
    Matrix A(A0);
    Matrix B(std::move(B0));

    matcl::Matrix ret;
    details::eval_plus::make(A,B,ret);
    return ret;
};

#pragma warning( pop )
