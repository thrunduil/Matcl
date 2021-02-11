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

#include "matcl-internals/base/utils.h"

namespace matcl { namespace details
{

namespace md = matcl::details;

template<class Derived>
class op_info
{
    public:
        template<class Val>
        bool is_special_case(const Val& val, int pos) const
        {
            return static_cast<const Derived*>(this)->is_special_case_impl(val,pos);
        };

        template<class Val>
        void eval_special_case(matcl::Matrix& ret, const Val& val, const Matrix& mat) const
        {
            return static_cast<const Derived*>(this)->eval_special_case_impl(ret, val, mat);
        };

        template<class Val>
        void eval_special_case(matcl::Matrix& ret,const Matrix& mat, const Val& val) const
        {
            return static_cast<const Derived*>(this)->eval_special_case_impl(ret, mat, val);
        };

        template<class UT, class TI1, class TI2>
        void eval_special_case_00(matcl::Matrix& ret, Integer m, Integer n, TI1 t1, TI2 t2) const
        {
            return static_cast<const Derived*>(this)
                ->template eval_special_case_00_impl<UT,TI1,TI2>(ret, m, n, t1, t2);
        };
};

template<class Derived>
class op_no_cases : public op_info<Derived>
{
    public:
        template<class Val>
        bool is_special_case_impl(const Val& , int ) const { return false; };

        template<class Val>
        void eval_special_case_impl(matcl::Matrix& ret, const Val&, const Matrix& mat) const
        { 
            ret = mat; 
        };

        template<class Val>
        void eval_special_case_impl(matcl::Matrix& ret, const Matrix& mat, const Val&) const
        { 
            ret = mat; 
        };
};

template<class M1, class M2, class op_info_type, bool iss>
struct function_op_call_ss_impl
{};

template<class M1, class M2, class op_info_type>
struct function_op_call_ss_impl<M1,M2,op_info_type,true>
{
    static void eval(matcl::Matrix& ret, const M1& mat,const M2& val, const op_info_type& op)
    {
        (void)op;
        return op.eval_scal_scal<M1,M2>(ret,mat,val);
    };
};

template<class M1, class M2, class op_info_type>
struct function_op_call_ss 
    : public function_op_call_ss_impl<M1,M2,op_info_type, 
                    md::is_scalar<M1>::value && md::is_scalar<M2>::value>
{};

template<class M1, class M2, class op_info_type, bool allow_real>
struct function_op
{
    static void eval(matcl::Matrix& ret, const M1& mat1, const M2& mat2, const op_info_type& op)
    {
        using value_type_1      = typename M1::value_type;
        using value_type_2      = typename M2::value_type;
        using value_type_1_real = typename details::real_type<value_type_1>::type;
        using value_type_2_real = typename details::real_type<value_type_2>::type;
        using UT                = typename md::unify_types<value_type_1,value_type_2>::type;

        typename ti::get_ti_type<M1>::type ti_1 = ti::get_ti(mat1);
        typename ti::get_ti_type<M2>::type ti_2 = ti::get_ti(mat2);

        value_type_1 Z1 = md::default_value<value_type_1>(ti_1);
        value_type_2 Z2 = md::default_value<value_type_2>(ti_2);

        if (mat1.rows() == 1 && mat1.cols() == 1)
        {			
            value_type_1 val(mat1(1,1));

            if (mat2.rows() == 1 && mat2.cols() == 1)
            {			
                value_type_2 val2(mat2(1,1));
                return function_op_call_ss<value_type_1,value_type_2,op_info_type>
                            ::eval(ret, val, val2, op);
            };

            value_type_1_real val_im(mrd::imag_helper<value_type_1>::eval(val));

            if (mrd::is_zero(val_im))
            {
                value_type_1_real val_re(mrd::real_helper<value_type_1>::eval(val));

                if (op.is_special_case(val_re,1))
                    return op.eval_special_case(ret,val_re,Matrix(mat2,true));

                return op.eval_scal_mat<value_type_1_real,M2>(ret,val_re,mat2);
            }
            else
            {
                if (op.is_special_case(val,1))
                    return op.eval_special_case(ret,val,Matrix(mat2,true));

                return op.eval_scal_mat<value_type_1,M2>(ret,val,mat2);
            };
        };

        if (mat2.rows() == 1 && mat2.cols() == 1)
        {			
            value_type_2 val(mat2(1,1));

            value_type_2_real val_im(mrd::imag_helper<value_type_2>::eval(val));

            if (mrd::is_zero(val_im))
            {
                value_type_2_real val_re(mrd::real_helper<value_type_2>::eval(val));

                if (op.is_special_case(val_re,2))
                    return op.eval_special_case(ret,Matrix(mat1,true),val_re);

                return op.eval_mat_scal<M1,value_type_2_real>(ret,mat1,val_re);
            }
            else
            {
                if (op.is_special_case(val,2))
                    return op.eval_special_case(ret,Matrix(mat1,true),val);

                return op.eval_mat_scal<M1,value_type_2>(ret,mat1,val);
            };
        };

        Integer m1 = mat1.rows();
        Integer n1 = mat1.cols();
        Integer m2 = mat2.rows();
        Integer n2 = mat2.cols();
        error::check_eeop(m1,n1,m2,n2);
        
        if (mat1.nnz() == 0)
        {
            if (mat2.nnz() == 0)
            {
                return op.eval_special_case_00<UT>(ret,m1,n1,ti_1,ti_2);
            }
            else
            {
                if (op.is_special_case(Z1,1))
                    return op.eval_special_case(ret,Z1,Matrix(mat2,true));
            };
        }
        else if (mat2.nnz() == 0)
        {
            if (op.is_special_case(Z2,2))
                return op.eval_special_case(ret,Matrix(mat1,true),Z2);
        };

        return op.eval_mat_mat<M1,M2>(ret,mat1,mat2);
    }
};

template<class M1, class M2, class op_info_type>
struct function_op<M1,M2,op_info_type,false>
{
    static void eval(matcl::Matrix& ret, const M1& mat1, const M2& mat2, const op_info_type& op)
    {
        using value_type_1      = typename M1::value_type;
        using value_type_2      = typename M2::value_type;
        using UT                = typename md::unify_types<value_type_1,value_type_2>::type;

        typename ti::get_ti_type<M1>::type ti_1 = ti::get_ti(mat1);
        typename ti::get_ti_type<M2>::type ti_2 = ti::get_ti(mat2);

        value_type_1 Z1 = md::default_value<value_type_1>(ti_1);
        value_type_2 Z2 = md::default_value<value_type_2>(ti_2);

        if (mat1.rows() == 1 && mat1.cols() == 1)
        {			
            value_type_1 val(mat1(1,1));

            if (mat2.rows() == 1 && mat2.cols() == 1)
            {			
                value_type_2 val2(mat2(1,1));
                return function_op_call_ss<value_type_1,value_type_2,op_info_type>
                            ::eval(ret, val, val2, op);
            };

            if (op.is_special_case(val,1))
                return op.eval_special_case(ret,val,Matrix(mat2,true));

            return op.eval_scal_mat<value_type_1,M2>(ret,val,mat2);
        };

        if (mat2.rows() == 1 && mat2.cols() == 1)
        {			
            value_type_2 val(mat2(1,1));

            if (op.is_special_case(val,2))
                return op.eval_special_case(ret,Matrix(mat1,true),val);

            return op.eval_mat_scal<M1,value_type_2>(ret,mat1,val);
        };

        Integer m1 = mat1.rows();
        Integer n1 = mat1.cols();
        Integer m2 = mat2.rows();
        Integer n2 = mat2.cols();
        error::check_eeop(m1,n1,m2,n2);
        
        if (mat1.nnz() == 0)
        {
            if (mat2.nnz() == 0)
            {
                return op.eval_special_case_00<UT>(ret,m1,n1,ti_1,ti_2);
            }
            else
            {
                if (op.is_special_case(Z1,1))
                    return op.eval_special_case(ret,Z1,Matrix(mat2,true));
            };
        }
        else if (mat2.nnz() == 0)
        {
            if (op.is_special_case(Z2,2))
                return op.eval_special_case(ret,Matrix(mat1,true),Z2);
        };

        return op.eval_mat_mat<M1,M2>(ret,mat1,mat2);
    }
};

}}
