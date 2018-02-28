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

#pragma once

#include "matcl-matrep/lib_functions/func_unary.h"
#include "matcl-matrep/func/raw/raw_func_unary.h"
#include "matcl-internals/base/instantiate.h"
#include "matcl-matrep/func/raw/eval_functor.h"
#include "matcl-scalar/details/scalfunc_helpers.h"
#include "matcl-matrep/details/extract_type_switch.h"
#include "matcl-internals/func/converter.h"
#include "matcl-internals/func/test_functor.h" 

#pragma warning(push)
#pragma warning(disable : 4127) //conditional expr is constant

namespace matcl { namespace details
{

namespace md    = matcl::details;
namespace mr    = matcl::raw;
namespace mrd   = matcl::raw::details;

template<class M1> struct in_range_asech	
{	
    static bool eval(const M1& A)
    {
        test_0_1 test_object;
        return test_range<typename M1::value_type, typename M1::struct_type,test_0_1>::eval(A, test_object);
    };	
};

template<class M1> struct in_range_acoth	
{	
    static bool eval(const M1& A)
    {
        test_inf_m1_1_inf test_object;
        return test_range<typename M1::value_type, typename M1::struct_type,test_inf_m1_1_inf>::eval(A, test_object);
    };	
};

template<class M1> struct in_range_atanh	
{	
    static bool eval(const M1& A)
    {
        test_m1_1 test_object;
        return test_range<typename M1::value_type, typename M1::struct_type,test_m1_1>::eval(A, test_object);
    };	
};

template<class M1>
struct function_op_pow2
{
    static void eval(matcl::Matrix& ret, const M1& A)
    {
        return mrd::scalar_func_helper<M1>::eval_pow2(ret,A);
    };
};

template<class M1> struct in_range_acosh	
{	
    static bool eval(const M1& A)
    {
        test_1_inf test_object;
        return test_range<typename M1::value_type, typename M1::struct_type,test_1_inf>::eval(A, test_object);
    };	
};

template<class M1> struct in_range_acsc	
{	
    static bool eval(const M1& A)
    {
        test_inf_m1_1_inf test_object;
        return test_range<typename M1::value_type, typename M1::struct_type,test_inf_m1_1_inf>::eval(A, test_object);
    };	
};

template<class M1> struct in_range_asec	
{	
    static bool eval(const M1& A)
    {
        test_inf_m1_1_inf test_object;
        return test_range<typename M1::value_type, typename M1::struct_type,test_inf_m1_1_inf>::eval(A, test_object);
    };	
};

template<class M1> struct in_range_acos	
{	
    static bool eval(const M1& A)
    {
        test_m1_1 test_object;
        return test_range<typename M1::value_type, typename M1::struct_type,test_m1_1>::eval(A, test_object);
    };	
};

template<class M1> struct in_range_asin	
{	
    static bool eval(const M1& A)
    {
        test_m1_1 test_object;
        return test_range<typename M1::value_type, typename M1::struct_type,test_m1_1>::eval(A, test_object);
    };	
};

template<class M1> struct in_range_log		
{	
    static bool eval(const M1& A)
    {
        test_geq_zero test_object;
        return test_range<typename M1::value_type, typename M1::struct_type,test_geq_zero>::eval(A, test_object);
    };	
};

template<class M1> struct in_range_log1p
{	
    static bool eval(const M1& A)
    {
        test_m1_inf test_object;
        return test_range<typename M1::value_type, typename M1::struct_type,test_m1_inf>
                ::eval(A, test_object);
    };	
};

template<class M1> struct in_range_sqrt	
{	
    static bool eval(const M1& A)
    {
        test_geq_zero test_object;
        return test_range<typename M1::value_type, typename M1::struct_type,test_geq_zero>
            ::eval(A, test_object);
    };	
};

template<class M1> struct in_range_sqrt1pm1	
{	
    static bool eval(const M1& A)
    {
        test_m1_inf test_object;
        return test_range<typename M1::value_type, typename M1::struct_type,test_m1_inf>
            ::eval(A, test_object);
    };	
};

struct eval_sgn_min : public extract_type_switch<void,eval_sgn_min,true> 
{
    template<class T>
    static void eval(const Matrix&, const T& mat, matcl::Matrix& ret)
    {
        mrd::unary_helper_impl<T>::eval_minus(ret,mat);
    };

    template<class T>
    static void eval_scalar(const Matrix& , const T& mat, matcl::Matrix& ret)
    {
        ret = raw::details::uminus_helper<T>::eval(mat);
    };
};

struct eval_inv : public extract_type_switch<void,eval_inv,true> 
{
    template<class T>
    static void eval(const Matrix&, const T& mat, matcl::Matrix& ret)
    {
        mrd::unary_helper_impl<T>::eval_inv(ret,mat);
    };

    template<class T>
    static void eval_scalar(const Matrix& , const T& mat, matcl::Matrix& ret)
    {
        ret = raw::details::invs_helper<T>::eval(mat);
    };
};

struct eval_neg : public extract_type_switch<void,eval_neg,true> 
{
    template<class T>
    static void eval(const Matrix& , const T& mat, matcl::Matrix& ret)
    {
        return mrd::unary_helper_impl<T>::eval_neg(ret, mat);
    };

    template<class T>
    static void eval_scalar(const Matrix& , const T& mat, matcl::Matrix& ret)
    {
        ret = Matrix(raw::details::op_neg_helper<T>::eval(mat));
    };
};

struct eval_is_finite : public extract_type_switch<void,eval_is_finite,true> 
{
    template<class T>
    static void eval(const Matrix& , const T& mat, matcl::Matrix& ret)
    {
        return mrd::scalfunc_isa_helper<T>::eval_is_finite(ret,mat);
    };

    template<class T>
    static void eval_scalar(const Matrix& , const T& mat, matcl::Matrix& ret)
    {
        ret = Matrix(mrd::isfinite_helper<T>::eval(mat));
    };
};

struct eval_is_nan : public extract_type_switch<void,eval_is_nan,true> 
{
    template<class T>
    static void eval(const Matrix& , const T& mat, matcl::Matrix& ret)
    {
        return mrd::scalfunc_isa_helper<T>::eval_is_nan(ret, mat);
    };

    template<class T>
    static void eval_scalar(const Matrix& , const T& mat, matcl::Matrix& ret)
    {
        ret = Matrix(mrd::isnan_helper<T>::eval(mat));
    };
};

struct eval_is_inf : public extract_type_switch<void,eval_is_inf,true> 
{
    template<class T>
    static void eval(const Matrix& , const T& mat, matcl::Matrix& ret)
    {
        return mrd::scalfunc_isa_helper<T>::eval_is_inf(ret,mat);
    };

    template<class T>
    static void eval_scalar(const Matrix& , const T& mat, matcl::Matrix& ret)
    {
        ret = Matrix(mrd::isinf_helper<T>::eval(mat));
    };
};

struct eval_is_regular : public extract_type_switch<void,eval_is_regular,true> 
{
    template<class T>
    static void eval(const Matrix& , const T& mat, matcl::Matrix& ret)
    {
        return mrd::scalfunc_isa_helper<T>::eval_is_regular(ret,mat);
    };

    template<class T>
    static void eval_scalar(const Matrix& , const T& mat, matcl::Matrix& ret)
    {
        ret = Matrix(mrd::isregular_helper<T>::eval(mat));
    };
};

struct eval_is_normal : public extract_type_switch<void,eval_is_normal,true> 
{
    template<class T>
    static void eval(const Matrix& , const T& mat, matcl::Matrix& ret)
    {
        return mrd::scalfunc_isa_helper<T>::eval_is_normal(ret,mat);
    };

    template<class T>
    static void eval_scalar(const Matrix& , const T& mat, matcl::Matrix& ret)
    {
        ret = Matrix(mrd::isnormal_helper<T>::eval(mat));
    };
};

struct eval_is_int : public extract_type_switch<void,eval_is_int,true> 
{
    template<class T>
    static void eval(const Matrix& , const T& mat, matcl::Matrix& ret)
    {
        return mrd::scalfunc_isa_helper<T>::eval_is_int(ret,mat);
    };

    template<class T>
    static void eval_scalar(const Matrix& , const T& mat, matcl::Matrix& ret)
    {
        ret = Matrix(mrd::isint_helper<T>::eval(mat));
    };
};

struct eval_is_real : public extract_type_switch<void,eval_is_real,true> 
{
    template<class T>
    static void eval(const Matrix& , const T& mat, matcl::Matrix& ret)
    {
        return mrd::scalfunc_isa_helper<T>::eval_is_real(ret,mat);
    };

    template<class T>
    static void eval_scalar(const Matrix& , const T& mat, matcl::Matrix& ret)
    {
        ret = Matrix(mrd::isreal_helper<T>::eval(mat));
    };
};

struct eval_imag : public extract_type_switch<void,eval_imag,true> 
{
    template<class T>
    static void eval(const Matrix& , const T& mat, matcl::Matrix& ret)
    {
        return mrd::scalfunc_real_helper<T>::eval_imag(ret,mat);
    };

    template<class T>
    static void eval_scalar(const Matrix& , const T& mat, matcl::Matrix& ret)
    {
        ret = mrd::imag_helper<T>::eval(mat);
    };
};

struct eval_real : public extract_type_switch<void,eval_real,true> 
{
    template<class T>
    static void eval(const Matrix& , const T& A, matcl::Matrix& ret)
    {
        if (!md::is_complex<typename T::value_type>::value &&
            !std::is_same<typename T::value_type,Object>::value)
        {
            ret = Matrix(A,false);
            return;
        };
        return mrd::scalfunc_real_helper<T>::eval_real(ret,A);
    };

    template<class T>
    static void eval_scalar(const Matrix& , const T& mat, matcl::Matrix& ret)
    {
        ret = mrd::real_helper<T>::eval(mat);
    };
};

struct eval_arg : public extract_type_switch<void,eval_arg,true> 
{
    template<class T>
    static void eval(const Matrix& , const T& mat, matcl::Matrix& ret)
    {
        return mrd::scalfunc_real_helper<T>::eval_arg(ret,mat);
    };

    template<class T>
    static void eval_scalar(const Matrix& , const T& mat, matcl::Matrix& ret)
    {
        ret = mrd::arg_helper<T>::eval(mat);
    };
};

struct eval_eps : public extract_type_switch<void,eval_eps,true> 
{
    template<class T>
    static void eval(const Matrix& , const T& mat, matcl::Matrix& ret)
    {
        return mrd::scalfunc_real_helper<T>::eval_eps(ret,mat);
    };

    template<class T>
    static void eval_scalar(const Matrix& , const T& mat, matcl::Matrix& ret)
    {
        ret = mrd::eps_helper<T>::eval(mat);
    };
};

struct eval_abs : public extract_type_switch<void,eval_abs,true> 
{
    template<class T>
    static void eval(const Matrix& , const T& mat, matcl::Matrix& ret)
    {
        return mrd::scalfunc_real_helper<T>::eval_abs(ret,mat);
    };

    template<class T>
    static void eval_scalar(const Matrix& , const T& mat, matcl::Matrix& ret)
    {
        ret = mrd::abs_helper<T>::eval(mat);
    };
};

struct eval_abs2 : public extract_type_switch<void,eval_abs2,true> 
{
    template<class T>
    static void eval(const Matrix& , const T& mat, matcl::Matrix& ret)
    {
        return mrd::scalfunc_real_helper<T>::eval_abs2(ret,mat);
    };

    template<class T>
    static void eval_scalar(const Matrix& , const T& mat, matcl::Matrix& ret)
    {
        ret = mrd::abs2_helper<T>::eval(mat);
    };
};


struct eval_is_true : public extract_type_switch<void,eval_is_true,true> 
{
    template<class T>
    static void eval(const Matrix& , const T& mat, matcl::Matrix& ret)
    {
        mrd::unary_helper_impl<T>::eval_is_true(ret,mat);
    };

    template<class T>
    static void eval_scalar(const Matrix& , const T& mat, matcl::Matrix& ret)
    {
        ret = Matrix(mrd::op_true_helper<T>::eval(mat));
    };
};

struct eval_sqrt : public extract_type_switch<void,eval_sqrt,true> 
{
    template<class T>
    static void eval(const Matrix& , const T& mat, matcl::Matrix& ret)
    {
        return mrd::scalar_func_helper<T>::eval_sqrt(ret,mat);
    };

    template<class T>
    static void eval_scalar(const Matrix& , const T& mat, matcl::Matrix& ret)
    {
        ret = mrd::sqrt_helper<T>::eval(mat);
    };
};

struct eval_sqrt_c : public extract_type_switch<void,eval_sqrt_c,true> 
{
    template<class M1>
    static void eval(const Matrix& , const M1& A, matcl::Matrix& ret)
    {
        if (md::is_complex<typename M1::value_type>::value)
            return mrd::scalar_func_helper<M1>::eval_sqrt(ret,A);

        if (in_range_sqrt<M1>::eval(A))
            return mrd::scalar_func_helper<M1>::eval_sqrt(ret,A);

        if (std::is_same<typename M1::value_type,Object>::value)
            return mrd::scalar_func_helper<M1>::eval_sqrt_c(ret,A);

        using matrix_type = raw::Matrix<Complex,typename M1::struct_type>;
        return mrd::scalar_func_helper<matrix_type>
                    ::eval_sqrt(ret,raw::converter<matrix_type,M1>::eval(A));
    };

    template<class T>
    static void eval_scalar(const Matrix& , const T& mat, matcl::Matrix& ret)
    {
        ret = mrd::sqrt_c_helper<T>::eval(mat);
    };
};

struct eval_sqrt_1pm1 : public extract_type_switch<void,eval_sqrt_1pm1,true> 
{
    template<class T>
    static void eval(const Matrix& , const T& mat, matcl::Matrix& ret)
    {
        return mrd::scalar_func_helper<T>::eval_sqrt1pm1(ret,mat);
    };

    template<class T>
    static void eval_scalar(const Matrix& , const T& mat, matcl::Matrix& ret)
    {
        ret = mrd::sqrt1pm1_helper<T>::eval(mat);
    };
};

struct eval_sqrt_1pm1_c : public extract_type_switch<void,eval_sqrt_1pm1_c,true> 
{
    template<class M1>
    static void eval(const Matrix& , const M1& A, matcl::Matrix& ret)
    {
        if (md::is_complex<typename M1::value_type>::value)
            return mrd::scalar_func_helper<M1>::eval_sqrt1pm1(ret,A);

        if (in_range_sqrt1pm1<M1>::eval(A))
            return mrd::scalar_func_helper<M1>::eval_sqrt1pm1(ret,A);

        if (std::is_same<typename M1::value_type,Object>::value)
            return mrd::scalar_func_helper<M1>::eval_sqrt1pm1_c(ret,A);

        using matrix_type = raw::Matrix<Complex,typename M1::struct_type>;
        return mrd::scalar_func_helper<matrix_type>
                    ::eval_sqrt1pm1(ret,raw::converter<matrix_type,M1>::eval(A));
    };

    template<class T>
    static void eval_scalar(const Matrix& , const T& mat, matcl::Matrix& ret)
    {
        ret = mrd::sqrt1pm1_c_helper<T>::eval(mat);
    };
};

struct eval_pow2 : public extract_type_switch<void,eval_pow2,true> 
{
    template<class T>
    static void eval(const Matrix& , const T& mat, matcl::Matrix& ret)
    {
        return function_op_pow2<T>::eval(ret,mat);
    };

    template<class T>
    static void eval_scalar(const Matrix& , const T& mat, matcl::Matrix& ret)
    {
        ret = mrd::pow2_helper<T>::eval(mat);
    };
};

struct eval_log_c : public extract_type_switch<void,eval_log_c,true> 
{
    template<class M1>
    static void eval(const Matrix& , const M1& A, matcl::Matrix& ret)
    {
        if (md::is_complex<typename M1::value_type>::value)
            return mrd::scalar_func_helper<M1>::eval_log(ret,A);

        if (in_range_log<M1>::eval(A))
            return mrd::scalar_func_helper<M1>::eval_log(ret,A);

        if (std::is_same<typename M1::value_type,Object>::value)
            return mrd::scalar_func_helper<M1>::eval_log_c(ret,A);

        using matrix_type = raw::Matrix<Complex,typename M1::struct_type>;
        return mrd::scalar_func_helper<matrix_type>
                    ::eval_log(ret,raw::converter<matrix_type,M1>::eval(A));
    };

    template<class T>
    static void eval_scalar(const Matrix& , const T& mat, matcl::Matrix& ret)
    {
        ret = mrd::log_c_helper<T>::eval(mat);
    };
};

struct eval_log_1p : public extract_type_switch<void,eval_log_1p,true> 
{
    template<class M1>
    static void eval(const Matrix& , const M1& A, matcl::Matrix& ret)
    {
        return mrd::scalar_func_helper<M1>::eval_log1p(ret,A);
    };

    template<class T>
    static void eval_scalar(const Matrix& , const T& mat, matcl::Matrix& ret)
    {
        ret = mrd::log1p_helper<T>::eval(mat);
    };
};

struct eval_log_1p_c : public extract_type_switch<void,eval_log_1p_c,true> 
{
    template<class M1>
    static void eval(const Matrix& , const M1& A, matcl::Matrix& ret)
    {
        if (md::is_complex<typename M1::value_type>::value)
            return mrd::scalar_func_helper<M1>::eval_log1p(ret,A);

        if (in_range_log1p<M1>::eval(A))
            return mrd::scalar_func_helper<M1>::eval_log1p(ret,A);

        if (std::is_same<typename M1::value_type,Object>::value)
            return mrd::scalar_func_helper<M1>::eval_log1p_c(ret,A);

        using matrix_type = raw::Matrix<Complex,typename M1::struct_type>;
        return mrd::scalar_func_helper<matrix_type>
                    ::eval_log1p(ret,raw::converter<matrix_type,M1>::eval(A));
    };

    template<class T>
    static void eval_scalar(const Matrix& , const T& mat, matcl::Matrix& ret)
    {
        ret = mrd::log1p_c_helper<T>::eval(mat);
    };
};

struct eval_exp : public extract_type_switch<void,eval_exp,true> 
{
    template<class T>
    static void eval(const Matrix& , const T& mat, matcl::Matrix& ret)
    {
        return mrd::scalar_func_helper<T>::eval_exp(ret,mat);
    };

    template<class T>
    static void eval_scalar(const Matrix& , const T& mat, matcl::Matrix& ret)
    {
        ret = mrd::exp_helper<T>::eval(mat);
    };
};

struct eval_log2 : public extract_type_switch<void,eval_log2,true> 
{
    template<class T>
    static void eval(const Matrix& , const T& mat, matcl::Matrix& ret)
    {
        return mrd::scalar_func_helper<T>::eval_log2(ret,mat);
    };

    template<class T>
    static void eval_scalar(const Matrix& , const T& mat, matcl::Matrix& ret)
    {
        ret = mrd::log2_helper<T>::eval(mat);
    };
};

struct eval_log2_c : public extract_type_switch<void,eval_log2_c,true> 
{
    template<class M1>
    static void eval(const Matrix& , const M1& A, matcl::Matrix& ret)
    {
        if (md::is_complex<typename M1::value_type>::value)
            return mrd::scalar_func_helper<M1>::eval_log2(ret,A);

        if (in_range_log<M1>::eval(A))
            return mrd::scalar_func_helper<M1>::eval_log2(ret,A);

        if (std::is_same<typename M1::value_type,Object>::value)
            return mrd::scalar_func_helper<M1>::eval_log2_c(ret,A);

        using matrix_type = raw::Matrix<Complex,typename M1::struct_type>;
        return mrd::scalar_func_helper<matrix_type>
                    ::eval_log2(ret,raw::converter<matrix_type,M1>::eval(A));
    };

    template<class T>
    static void eval_scalar(const Matrix& , const T& mat, matcl::Matrix& ret)
    {
        ret = mrd::log2_c_helper<T>::eval(mat);
    };
};

struct eval_log10 : public extract_type_switch<void,eval_log10,true> 
{
    template<class T>
    static void eval(const Matrix& , const T& mat, matcl::Matrix& ret)
    {
        return mrd::scalar_func_helper<T>::eval_log10(ret,mat);
    };

    template<class T>
    static void eval_scalar(const Matrix& , const T& mat, matcl::Matrix& ret)
    {
        ret = mrd::log10_helper<T>::eval(mat);
    };
};

struct eval_log10_c : public extract_type_switch<void,eval_log10_c,true> 
{
    template<class M1>
    static void eval(const Matrix& , const M1& A, matcl::Matrix& ret)
    {
        if (md::is_complex<typename M1::value_type>::value)
            return mrd::scalar_func_helper<M1>::eval_log10(ret,A);

        if (in_range_log<M1>::eval(A))
            return mrd::scalar_func_helper<M1>::eval_log10(ret,A);

        if (std::is_same<typename M1::value_type,Object>::value)
            return mrd::scalar_func_helper<M1>::eval_log10_c(ret,A);

        using matrix_type = raw::Matrix<Complex,typename M1::struct_type>;
        return mrd::scalar_func_helper<matrix_type>
                    ::eval_log10(ret,raw::converter<matrix_type,M1>::eval(A));
    };

    template<class T>
    static void eval_scalar(const Matrix& , const T& mat, matcl::Matrix& ret)
    {
        ret = mrd::log10_c_helper<T>::eval(mat);
    };
};

struct eval_floor : public extract_type_switch<void,eval_floor,true> 
{
    template<class M1>
    static void eval(const Matrix& , const M1& A, matcl::Matrix& ret)
    {
        if (std::is_same<typename M1::value_type,Integer>::value)
        {
            ret = Matrix(A,false);
            return;
        };
        return mrd::scalar_func_helper<M1>::eval_floor(ret,A);
    };

    template<class T>
    static void eval_scalar(const Matrix& , const T& mat, matcl::Matrix& ret)
    {
        ret = mrd::floor_helper<T>::eval(mat);
    };
};

struct eval_ceil : public extract_type_switch<void,eval_ceil,true> 
{
    template<class M1>
    static void eval(const Matrix& , const M1& A, matcl::Matrix& ret)
    {
        if (std::is_same<typename M1::value_type,Integer>::value)
        {
            ret = Matrix(A,false);
            return;
        };
        return mrd::scalar_func_helper<M1>::eval_ceil(ret,A);
    };

    template<class T>
    static void eval_scalar(const Matrix& , const T& mat, matcl::Matrix& ret)
    {
        ret = mrd::ceil_helper<T>::eval(mat);
    };
};

struct eval_round : public extract_type_switch<void,eval_round,true> 
{
    template<class M1>
    static void eval(const Matrix& , const M1& A, matcl::Matrix& ret)
    {
        if (std::is_same<typename M1::value_type,Integer>::value)
        {
            ret = Matrix(A,false);
            return;
        };
        return mrd::scalar_func_helper<M1>::eval_round(ret,A);
    };

    template<class T>
    static void eval_scalar(const Matrix& , const T& mat, matcl::Matrix& ret)
    {
        ret = mrd::round_helper<T>::eval(mat);
    };
};

struct eval_fix : public extract_type_switch<void,eval_fix,true> 
{
    template<class M1>
    static void eval(const Matrix& , const M1& A, matcl::Matrix& ret)
    {
        if (std::is_same<typename M1::value_type,Integer>::value)
        {
            ret = Matrix(A,false);
            return;
        };
        return mrd::scalar_func_helper<M1>::eval_fix(ret,A);
    };

    template<class T>
    static void eval_scalar(const Matrix& , const T& mat, matcl::Matrix& ret)
    {
        ret = mrd::fix_helper<T>::eval(mat);
    };
};

struct eval_trunc : public extract_type_switch<void,eval_trunc,true> 
{
    template<class M1>
    static void eval(const Matrix& , const M1& A, matcl::Matrix& ret)
    {
        if (std::is_same<typename M1::value_type,Integer>::value)
        {
            ret = Matrix(A,false);
            return;
        };
        return mrd::scalar_func_helper<M1>::eval_trunc(ret,A);
    };

    template<class T>
    static void eval_scalar(const Matrix& , const T& mat, matcl::Matrix& ret)
    {
        ret = mrd::trunc_helper<T>::eval(mat);
    };
};

struct eval_ifloor : public extract_type_switch<void,eval_ifloor,true> 
{
    template<class M1>
    static void eval(const Matrix& , const M1& A, matcl::Matrix& ret)
    {
        if (std::is_same<typename M1::value_type,Integer>::value)
        {
            ret = Matrix(A,false);
            return;
        };
        return mrd::scalar_func_helper<M1>::eval_ifloor(ret,A);
    };

    template<class T>
    static void eval_scalar(const Matrix& , const T& mat, matcl::Matrix& ret)
    {
        ret = mrd::ifloor_helper<T>::eval(mat);
    };
};

struct eval_iceil : public extract_type_switch<void,eval_iceil,true> 
{
    template<class M1>
    static void eval(const Matrix& , const M1& A, matcl::Matrix& ret)
    {
        if (std::is_same<typename M1::value_type,Integer>::value)
        {
            ret = Matrix(A,false);
            return;
        };
        return mrd::scalar_func_helper<M1>::eval_iceil(ret,A);
    };

    template<class T>
    static void eval_scalar(const Matrix& , const T& mat, matcl::Matrix& ret)
    {
        ret = mrd::iceil_helper<T>::eval(mat);
    };
};

struct eval_iround : public extract_type_switch<void,eval_iround,true> 
{
    template<class M1>
    static void eval(const Matrix& , const M1& A, matcl::Matrix& ret)
    {
        if (std::is_same<typename M1::value_type,Integer>::value)
        {
            ret = Matrix(A,false);
            return;
        };
        return mrd::scalar_func_helper<M1>::eval_iround(ret,A);
    };

    template<class T>
    static void eval_scalar(const Matrix& , const T& mat, matcl::Matrix& ret)
    {
        ret = mrd::iround_helper<T>::eval(mat);
    };
};

struct eval_ifix : public extract_type_switch<void,eval_ifix,true> 
{
    template<class M1>
    static void eval(const Matrix& , const M1& A, matcl::Matrix& ret)
    {
        if (std::is_same<typename M1::value_type,Integer>::value)
        {
            ret = Matrix(A,false);
            return;
        };
        return mrd::scalar_func_helper<M1>::eval_ifix(ret,A);
    };

    template<class T>
    static void eval_scalar(const Matrix& , const T& mat, matcl::Matrix& ret)
    {
        ret = mrd::ifix_helper<T>::eval(mat);
    };
};

struct eval_itrunc : public extract_type_switch<void,eval_itrunc,true> 
{
    template<class M1>
    static void eval(const Matrix& , const M1& A, matcl::Matrix& ret)
    {
        if (std::is_same<typename M1::value_type,Integer>::value)
        {
            ret = Matrix(A,false);
            return;
        };
        return mrd::scalar_func_helper<M1>::eval_itrunc(ret,A);
    };

    template<class T>
    static void eval_scalar(const Matrix& , const T& mat, matcl::Matrix& ret)
    {
        ret = mrd::itrunc_helper<T>::eval(mat);
    };
};

struct eval_isign : public extract_type_switch<void,eval_isign,true> 
{
    template<class T>
    static void eval(const Matrix& , const T& mat, matcl::Matrix& ret)
    {
        return mrd::scalar_func_helper<T>::eval_isign(ret,mat);
    };

    template<class T>
    static void eval_scalar(const Matrix& , const T& mat, matcl::Matrix& ret)
    {
        ret = mrd::isign_helper<T>::eval(mat);
    };
};

struct eval_sign : public extract_type_switch<void,eval_sign,true> 
{
    template<class T>
    static void eval(const Matrix& , const T& mat, matcl::Matrix& ret)
    {
        return mrd::scalar_func_helper<T>::eval_sign(ret,mat);
    };

    template<class T>
    static void eval_scalar(const Matrix& , const T& mat, matcl::Matrix& ret)
    {
        ret = mrd::sign_helper<T>::eval(mat);
    };
};

struct eval_sin : public extract_type_switch<void,eval_sin,true> 
{
    template<class T>
    static void eval(const Matrix& , const T& mat, matcl::Matrix& ret)
    {
        return mrd::scalar_func_helper<T>::eval_sin(ret,mat);
    };

    template<class T>
    static void eval_scalar(const Matrix& , const T& mat, matcl::Matrix& ret)
    {
        ret = mrd::sin_helper<T>::eval(mat);
    };
};

struct eval_cos : public extract_type_switch<void,eval_cos,true> 
{
    template<class T>
    static void eval(const Matrix& , const T& mat, matcl::Matrix& ret)
    {
        return mrd::scalar_func_helper<T>::eval_cos(ret,mat);
    };

    template<class T>
    static void eval_scalar(const Matrix& , const T& mat, matcl::Matrix& ret)
    {
        ret = mrd::cos_helper<T>::eval(mat);
    };
};

struct eval_tan : public extract_type_switch<void,eval_tan,true> 
{
    template<class T>
    static void eval(const Matrix& , const T& mat, matcl::Matrix& ret)
    {
        return mrd::scalar_func_helper<T>::eval_tan(ret,mat);
    };

    template<class T>
    static void eval_scalar(const Matrix& , const T& mat, matcl::Matrix& ret)
    {
        ret = mrd::tan_helper<T>::eval(mat);;
    };
};

struct eval_cot : public extract_type_switch<void,eval_cot,true> 
{
    template<class T>
    static void eval(const Matrix& , const T& mat, matcl::Matrix& ret)
    {
        return mrd::scalar_func_helper<T>::eval_cot(ret,mat);
    };

    template<class T>
    static void eval_scalar(const Matrix& , const T& mat, matcl::Matrix& ret)
    {
        ret = mrd::cot_helper<T>::eval(mat);
    };
};

struct eval_sec : public extract_type_switch<void,eval_sec,true> 
{
    template<class T>
    static void eval(const Matrix& , const T& mat, matcl::Matrix& ret)
    {
        return mrd::scalar_func_helper<T>::eval_sec(ret,mat);
    };

    template<class T>
    static void eval_scalar(const Matrix& , const T& mat, matcl::Matrix& ret)
    {
        ret = mrd::sec_helper<T>::eval(mat);
    };
};

struct eval_csc : public extract_type_switch<void,eval_csc,true> 
{
    template<class T>
    static void eval(const Matrix& , const T& mat, matcl::Matrix& ret)
    {
        return mrd::scalar_func_helper<T>::eval_csc(ret,mat);
    };

    template<class T>
    static void eval_scalar(const Matrix& , const T& mat, matcl::Matrix& ret)
    {
        ret = mrd::csc_helper<T>::eval(mat);
    };
};

struct eval_asin : public extract_type_switch<void,eval_asin,true> 
{
    template<class T>
    static void eval(const Matrix& , const T& mat, matcl::Matrix& ret)
    {
        return mrd::scalar_func_helper<T>::eval_asin(ret,mat);
    };

    template<class T>
    static void eval_scalar(const Matrix& , const T& mat, matcl::Matrix& ret)
    {
        ret = mrd::asin_helper<T>::eval(mat);
    };
};

struct eval_asin_c : public extract_type_switch<void,eval_asin_c,true> 
{
    template<class M1>
    static void eval(const Matrix& , const M1& A, matcl::Matrix& ret)
    {
        if (md::is_complex<typename M1::value_type>::value)
            return mrd::scalar_func_helper<M1>::eval_asin(ret,A);

        if (in_range_asin<M1>::eval(A))
            return mrd::scalar_func_helper<M1>::eval_asin(ret,A);

        if (std::is_same<typename M1::value_type,Object>::value)
            return mrd::scalar_func_helper<M1>::eval_asin_c(ret,A);

        using matrix_type = raw::Matrix<Complex,typename M1::struct_type>;
        return mrd::scalar_func_helper<matrix_type>
                        ::eval_asin(ret,raw::converter<matrix_type,M1>::eval(A));
    };

    template<class T>
    static void eval_scalar(const Matrix& , const T& mat, matcl::Matrix& ret)
    {
        ret = mrd::asin_c_helper<T>::eval(mat);
    };
};

struct eval_acos : public extract_type_switch<void,eval_acos,true> 
{
    template<class T>
    static void eval(const Matrix& , const T& mat, matcl::Matrix& ret)
    {
        return mrd::scalar_func_helper<T>::eval_acos(ret,mat);
    };

    template<class T>
    static void eval_scalar(const Matrix& , const T& mat, matcl::Matrix& ret)
    {
        ret = mrd::acos_helper<T>::eval(mat);
    };
};

struct eval_acos_c : public extract_type_switch<void,eval_acos_c,true> 
{
    template<class M1>
    static void eval(const Matrix& , const M1& A, matcl::Matrix& ret)
    {
        if (md::is_complex<typename M1::value_type>::value)
            return mrd::scalar_func_helper<M1>::eval_acos(ret,A);

        if (in_range_acos<M1>::eval(A))
            return mrd::scalar_func_helper<M1>::eval_acos(ret,A);

        if (std::is_same<typename M1::value_type,Object>::value)
            return mrd::scalar_func_helper<M1>::eval_acos_c(ret,A);

        using matrix_type = raw::Matrix<Complex,typename M1::struct_type>;
        return mrd::scalar_func_helper<matrix_type>
                        ::eval_acos(ret,raw::converter<matrix_type,M1>::eval(A));
    };

    template<class T>
    static void eval_scalar(const Matrix& , const T& mat, matcl::Matrix& ret)
    {
        ret = mrd::acos_c_helper<T>::eval(mat);
    };
};

struct eval_atan : public extract_type_switch<void,eval_atan,true> 
{
    template<class T>
    static void eval(const Matrix& , const T& mat, matcl::Matrix& ret)
    {
        return mrd::scalar_func_helper<T>::eval_atan(ret,mat);
    };

    template<class T>
    static void eval_scalar(const Matrix& , const T& mat, matcl::Matrix& ret)
    {
        ret = mrd::atan_helper<T>::eval(mat);
    };
};

struct eval_acot : public extract_type_switch<void,eval_acot,true> 
{
    template<class T>
    static void eval(const Matrix& , const T& mat, matcl::Matrix& ret)
    {
        return mrd::scalar_func_helper<T>::eval_acot(ret,mat);
    };

    template<class T>
    static void eval_scalar(const Matrix& , const T& mat, matcl::Matrix& ret)
    {
        ret = mrd::acot_helper<T>::eval(mat);
    };
};

struct eval_asec : public extract_type_switch<void,eval_asec,true> 
{
    template<class T>
    static void eval(const Matrix& , const T& mat, matcl::Matrix& ret)
    {
        return mrd::scalar_func_helper<T>::eval_asec(ret,mat);
    };

    template<class T>
    static void eval_scalar(const Matrix& , const T& mat, matcl::Matrix& ret)
    {
        ret = mrd::asec_helper<T>::eval(mat);
    };
};

struct eval_asec_c : public extract_type_switch<void,eval_asec_c,true> 
{
    template<class M1>
    static void eval(const Matrix& , const M1& A, matcl::Matrix& ret)
    {
        if (md::is_complex<typename M1::value_type>::value)
            return mrd::scalar_func_helper<M1>::eval_asec(ret,A);

        if (in_range_asec<M1>::eval(A))
            return mrd::scalar_func_helper<M1>::eval_asec(ret,A);

        if (std::is_same<typename M1::value_type,Object>::value)
            return mrd::scalar_func_helper<M1>::eval_asec_c(ret,A);

        using matrix_type = raw::Matrix<Complex,typename M1::struct_type>;
        return mrd::scalar_func_helper<matrix_type>
                        ::eval_asec(ret,raw::converter<matrix_type,M1>::eval(A));
    };

    template<class T>
    static void eval_scalar(const Matrix& , const T& mat, matcl::Matrix& ret)
    {
        ret = mrd::asec_c_helper<T>::eval(mat);
    };
};

struct eval_acsc : public extract_type_switch<void,eval_acsc,true> 
{
    template<class M1>
    static void eval(const Matrix& , const M1& A, matcl::Matrix& ret)
    {
        return mrd::scalar_func_helper<M1>::eval_acsc(ret,A);
    };

    template<class T>
    static void eval_scalar(const Matrix& , const T& mat, matcl::Matrix& ret)
    {
        ret = mrd::acsc_helper<T>::eval(mat);
    };
};

struct eval_acsc_c : public extract_type_switch<void,eval_acsc_c,true> 
{
    template<class M1>
    static void eval(const Matrix& , const M1& A, matcl::Matrix& ret)
    {
        if (md::is_complex<typename M1::value_type>::value)
            return mrd::scalar_func_helper<M1>::eval_acsc(ret,A);

        if (in_range_acsc<M1>::eval(A))
            return mrd::scalar_func_helper<M1>::eval_acsc(ret,A);

        if (std::is_same<typename M1::value_type,Object>::value)
            return mrd::scalar_func_helper<M1>::eval_acsc_c(ret,A);

        using matrix_type = raw::Matrix<Complex,typename M1::struct_type>;
        return mrd::scalar_func_helper<matrix_type>
                        ::eval_acsc(ret,raw::converter<matrix_type,M1>::eval(A));
    };

    template<class T>
    static void eval_scalar(const Matrix& , const T& mat, matcl::Matrix& ret)
    {
        ret = mrd::acsc_c_helper<T>::eval(mat);
    };
};

struct eval_sinh : public extract_type_switch<void,eval_sinh,true> 
{
    template<class T>
    static void eval(const Matrix& , const T& mat, matcl::Matrix& ret)
    {
        return mrd::scalar_func_helper<T>::eval_sinh(ret,mat);
    };

    template<class T>
    static void eval_scalar(const Matrix& , const T& mat, matcl::Matrix& ret)
    {
        ret = mrd::sinh_helper<T>::eval(mat);
    };
};

struct eval_cosh : public extract_type_switch<void,eval_cosh,true> 
{
    template<class T>
    static void eval(const Matrix& , const T& mat, matcl::Matrix& ret)
    {
        return mrd::scalar_func_helper<T>::eval_cosh(ret,mat);
    };

    template<class T>
    static void eval_scalar(const Matrix& , const T& mat, matcl::Matrix& ret)
    {
        ret = mrd::cosh_helper<T>::eval(mat);
    };
};

struct eval_tanh : public extract_type_switch<void,eval_tanh,true> 
{
    template<class T>
    static void eval(const Matrix& , const T& mat, matcl::Matrix& ret)
    {
        return mrd::scalar_func_helper<T>::eval_tanh(ret,mat);
    };

    template<class T>
    static void eval_scalar(const Matrix& , const T& mat, matcl::Matrix& ret)
    {
        ret = mrd::tanh_helper<T>::eval(mat);
    };
};

struct eval_coth : public extract_type_switch<void,eval_coth,true> 
{
    template<class T>
    static void eval(const Matrix& , const T& mat, matcl::Matrix& ret)
    {
        return mrd::scalar_func_helper<T>::eval_coth(ret,mat);
    };

    template<class T>
    static void eval_scalar(const Matrix& , const T& mat, matcl::Matrix& ret)
    {
        ret = mrd::coth_helper<T>::eval(mat);
    };
};

struct eval_sech : public extract_type_switch<void,eval_sech,true> 
{
    template<class T>
    static void eval(const Matrix& , const T& mat, matcl::Matrix& ret)
    {
        return mrd::scalar_func_helper<T>::eval_sech(ret,mat);
    };

    template<class T>
    static void eval_scalar(const Matrix& , const T& mat, matcl::Matrix& ret)
    {
        ret = mrd::sech_helper<T>::eval(mat);
    };
};

struct eval_csch : public extract_type_switch<void,eval_csch,true> 
{
    template<class T>
    static void eval(const Matrix& , const T& mat, matcl::Matrix& ret)
    {
        return mrd::scalar_func_helper<T>::eval_csch(ret,mat);
    };

    template<class T>
    static void eval_scalar(const Matrix& , const T& mat, matcl::Matrix& ret)
    {
        ret = mrd::csch_helper<T>::eval(mat);
    };
};

struct eval_asinh : public extract_type_switch<void,eval_asinh,true> 
{
    template<class T>
    static void eval(const Matrix& , const T& mat, matcl::Matrix& ret)
    {
        return mrd::scalar_func_helper<T>::eval_asinh(ret,mat);
    };

    template<class T>
    static void eval_scalar(const Matrix& , const T& mat, matcl::Matrix& ret)
    {
        ret = mrd::asinh_helper<T>::eval(mat);
    };
};

struct eval_acosh : public extract_type_switch<void,eval_acosh,true> 
{
    template<class T>
    static void eval(const Matrix& , const T& mat, matcl::Matrix& ret)
    {
        return mrd::scalar_func_helper<T>::eval_acosh(ret,mat);
    };

    template<class T>
    static void eval_scalar(const Matrix& , const T& mat, matcl::Matrix& ret)
    {
        ret = mrd::acosh_helper<T>::eval(mat);
    };
};

struct eval_acosh_c : public extract_type_switch<void,eval_acosh_c,true> 
{
    template<class M1>
    static void eval(const Matrix& , const M1& A, matcl::Matrix& ret)
    {
        if (md::is_complex<typename M1::value_type>::value)
            return mrd::scalar_func_helper<M1>::eval_acosh(ret,A);

        if (in_range_acosh<M1>::eval(A))
            return mrd::scalar_func_helper<M1>::eval_acosh(ret,A);

        if (std::is_same<typename M1::value_type,Object>::value)
            return mrd::scalar_func_helper<M1>::eval_acosh_c(ret,A);

        using matrix_type = raw::Matrix<Complex,typename M1::struct_type>;
        return mrd::scalar_func_helper<matrix_type>
                    ::eval_acosh(ret,raw::converter<matrix_type,M1>::eval(A));
    };

    template<class T>
    static void eval_scalar(const Matrix& , const T& mat, matcl::Matrix& ret)
    {
        ret = mrd::acosh_c_helper<T>::eval(mat);
    };
};

struct eval_atanh : public extract_type_switch<void,eval_atanh,true> 
{
    template<class T>
    static void eval(const Matrix& , const T& mat, matcl::Matrix& ret)
    {
        return mrd::scalar_func_helper<T>::eval_atanh(ret,mat);
    };

    template<class T>
    static void eval_scalar(const Matrix& , const T& mat, matcl::Matrix& ret)
    {
        ret = mrd::atanh_helper<T>::eval(mat);
    };
};

struct eval_atanh_c : public extract_type_switch<void,eval_atanh_c,true> 
{
    template<class M1>
    static void eval(const Matrix& , const M1& A, matcl::Matrix& ret)
    {
        if (md::is_complex<typename M1::value_type>::value)
            return mrd::scalar_func_helper<M1>::eval_atanh(ret,A);

        if (in_range_atanh<M1>::eval(A))
            return mrd::scalar_func_helper<M1>::eval_atanh(ret,A);

        if (std::is_same<typename M1::value_type,Object>::value)
            return mrd::scalar_func_helper<M1>::eval_atanh_c(ret,A);

        using matrix_type = raw::Matrix<Complex,typename M1::struct_type>;
        return mrd::scalar_func_helper<matrix_type>
                    ::eval_atanh(ret,raw::converter<matrix_type,M1>::eval(A));
    };

    template<class T>
    static void eval_scalar(const Matrix& , const T& mat, matcl::Matrix& ret)
    {
        ret = mrd::atanh_c_helper<T>::eval(mat);
    };
};

struct eval_acoth : public extract_type_switch<void,eval_acoth,true> 
{
    template<class T>
    static void eval(const Matrix& , const T& mat, matcl::Matrix& ret)
    {
        return mrd::scalar_func_helper<T>::eval_acoth(ret,mat);
    };

    template<class T>
    static void eval_scalar(const Matrix& , const T& mat, matcl::Matrix& ret)
    {
        ret = mrd::acoth_helper<T>::eval(mat);
    };
};

struct eval_acoth_c : public extract_type_switch<void,eval_acoth_c,true> 
{
    template<class M1>
    static void eval(const Matrix& , const M1& A, matcl::Matrix& ret)
    {
        if (md::is_complex<typename M1::value_type>::value)
            return mrd::scalar_func_helper<M1>::eval_acoth(ret,A);

        if (in_range_acoth<M1>::eval(A))
            return mrd::scalar_func_helper<M1>::eval_acoth(ret,A);

        if (std::is_same<typename M1::value_type,Object>::value)
            return mrd::scalar_func_helper<M1>::eval_acoth_c(ret,A);

        using matrix_type = raw::Matrix<Complex,typename M1::struct_type>;
        return mrd::scalar_func_helper<matrix_type>
                            ::eval_acoth(ret,raw::converter<matrix_type,M1>::eval(A));
    };

    template<class T>
    static void eval_scalar(const Matrix& , const T& mat, matcl::Matrix& ret)
    {
        ret = mrd::acoth_c_helper<T>::eval(mat);
    };
};

struct eval_asech : public extract_type_switch<void,eval_asech,true> 
{
    template<class T>
    static void eval(const Matrix& , const T& mat, matcl::Matrix& ret)
    {
        return mrd::scalar_func_helper<T>::eval_asech(ret,mat);
    };

    template<class T>
    static void eval_scalar(const Matrix& , const T& mat, matcl::Matrix& ret)
    {
        ret = mrd::asech_helper<T>::eval(mat);
    };
};

struct eval_asech_c : public extract_type_switch<void,eval_asech_c,true> 
{
    template<class M1>
    static void eval(const Matrix& , const M1& A, matcl::Matrix& ret)
    {
        if (md::is_complex<typename M1::value_type>::value)
            return mrd::scalar_func_helper<M1>::eval_asech(ret,A);

        if (in_range_asech<M1>::eval(A))
            return mrd::scalar_func_helper<M1>::eval_asech(ret,A);

        if (std::is_same<typename M1::value_type,Object>::value)
            return mrd::scalar_func_helper<M1>::eval_asech_c(ret,A);

        using matrix_type = raw::Matrix<Complex,typename M1::struct_type>;
        return mrd::scalar_func_helper<matrix_type>
                        ::eval_asech(ret,mr::converter<matrix_type,M1>::eval(A));
    };

    template<class T>
    static void eval_scalar(const Matrix& , const T& mat, matcl::Matrix& ret)
    {
        ret = mrd::asech_c_helper<T>::eval(mat);
    };
};

struct eval_acsch : public extract_type_switch<void,eval_acsch,true> 
{
    template<class T>
    static void eval(const Matrix& , const T& mat, matcl::Matrix& ret)
    {
        return mrd::scalar_func_helper<T>::eval_acsch(ret,mat);
    };

    template<class T>
    static void eval_scalar(const Matrix& , const T& mat, matcl::Matrix& ret)
    {
        ret = mrd::acsch_helper<T>::eval(mat);
    };
};

struct eval_log : public extract_type_switch<void,eval_log,true> 
{
    template<class M1>
    static void eval(const Matrix& , const M1& A, matcl::Matrix& ret)
    {
        return mrd::scalar_func_helper<M1>::eval_log(ret,A);
    };

    template<class T>
    static void eval_scalar(const Matrix& , const T& mat, matcl::Matrix& ret)
    {
        ret = mrd::log_helper<T>::eval(mat);
    };
};

struct func_cbrt
{
    template<class T>
    auto eval(const T& v) const -> decltype(mrd::cbrt_helper<T>::eval(v))
    {
        return mrd::cbrt_helper<T>::eval(v);
    }
};

struct func_isnormal
{
    template<class T>
    auto eval(const T& v) const-> decltype(mrd::is_normal_helper<T>::eval(v))
    {
        return mrd::is_normal_helper<T>::eval(v);
    }
};

struct func_expm1
{
    template<class T>
    auto eval(const T& v) const -> decltype(mrd::expm1_helper<T>::eval(v))
    {
        return mrd::expm1_helper<T>::eval(v);
    }
};

struct func_expi
{
    template<class T>
    auto eval(const T& v) const -> decltype(mrd::expi_helper<T>::eval(v))
    {
        return mrd::expi_helper<T>::eval(v);
    }
};

struct func_signbit
{
    template<class T>
    auto eval(const T& v) const -> decltype(mrd::signbit_helper<T>::eval(v))
    {
        return mrd::signbit_helper<T>::eval(v);
    }
};

struct func_logb
{
    template<class T>
    auto eval(const T& v) const -> decltype(mrd::logb_helper<T>::eval(v))
    {
        return mrd::logb_helper<T>::eval(v);
    }
};

struct func_ilogb
{
    template<class T>
    auto eval(const T& v) const -> decltype(mrd::ilogb_helper<T>::eval(v))
    {
        return mrd::ilogb_helper<T>::eval(v);
    }
};

struct func_exp2
{
    template<class T>
    auto eval(const T& v) const -> decltype(mrd::exp2_helper<T>::eval(v))
    {
        return mrd::exp2_helper<T>::eval(v);
    }
};

struct func_exp10
{
    template<class T>
    auto eval(const T& v) const -> decltype(mrd::exp10_helper<T>::eval(v))
    {
        return mrd::exp10_helper<T>::eval(v);
    }
};

struct eval_conj : public extract_type_switch<void,eval_conj,true> 
{
    template<class M1>
    static void eval(const Matrix& , const M1& A, matcl::Matrix& ret)
    {
        if (!md::is_complex<typename M1::value_type>::value &&
            !std::is_same<typename M1::value_type,Object>::value)
        {
            ret = Matrix(A,false);
            return;
        };
        return mrd::scalfunc_real_helper<M1>::eval_conj(ret,A);
    };

    template<class T>
    static void eval_scalar(const Matrix& , const T& mat, matcl::Matrix& ret)
    {
        ret = mrd::conj_helper<T>::eval(mat);
    };
};

};};

namespace matcl
{

Matrix matcl::is_finite(const Matrix& A0)
{
    Matrix A(A0);

    matcl::Matrix ret;
    details::eval_is_finite::make<const Matrix&>(A,ret);
    return ret;
};

Matrix matcl::is_finite(Matrix&& A0)
{
    Matrix A(std::move(A0));

    matcl::Matrix ret;
    details::eval_is_finite::make<const Matrix&>(A,ret);
    return ret;
};

Matrix matcl::is_nan(const Matrix& A0)
{
    Matrix A(A0);

    matcl::Matrix ret;
    details::eval_is_nan::make<const Matrix&>(A,ret);
    return ret;
};

Matrix matcl::is_nan(Matrix&& A0)
{
    Matrix A(std::move(A0));

    matcl::Matrix ret;
    details::eval_is_nan::make<const Matrix&>(A,ret);
    return ret;
};

Matrix matcl::is_inf(const Matrix& A0)
{
    Matrix A(A0);

    matcl::Matrix ret;
    details::eval_is_inf::make<const Matrix&>(A,ret);
    return ret;
};

Matrix matcl::is_inf(Matrix&& A0)
{
    Matrix A(std::move(A0));

    matcl::Matrix ret;
    details::eval_is_inf::make<const Matrix&>(A,ret);
    return ret;
};

Matrix matcl::is_regular(const Matrix& A0)
{
    Matrix A(A0);

    matcl::Matrix ret;
    details::eval_is_regular::make<const Matrix&>(A,ret);
    return ret;
};

Matrix matcl::is_regular(Matrix&& A0)
{
    Matrix A(std::move(A0));

    matcl::Matrix ret;
    details::eval_is_regular::make<const Matrix&>(A,ret);
    return ret;
};

Matrix matcl::is_normal(const Matrix& A0)
{
    Matrix A(A0);

    matcl::Matrix ret;
    details::eval_is_normal::make<const Matrix&>(A,ret);
    return ret;
};

Matrix matcl::is_normal(Matrix&& A0)
{
    Matrix A(std::move(A0));

    matcl::Matrix ret;
    details::eval_is_normal::make<const Matrix&>(A,ret);
    return ret;
};

Matrix matcl::is_int(const Matrix& A0)
{
    Matrix A(A0);

    matcl::Matrix ret;
    details::eval_is_int::make<const Matrix&>(A,ret);
    return ret;
};

Matrix matcl::is_int(Matrix&& A0)
{
    Matrix A(std::move(A0));

    matcl::Matrix ret;
    details::eval_is_int::make<const Matrix&>(A,ret);
    return ret;
};

Matrix matcl::is_real(const Matrix& A0)
{
    Matrix A(A0);

    matcl::Matrix ret;
    details::eval_is_real::make<const Matrix&>(A,ret);
    return ret;
};

Matrix matcl::is_real(Matrix&& A0)
{
    Matrix A(std::move(A0));

    matcl::Matrix ret;
    details::eval_is_real::make<const Matrix&>(A,ret);
    return ret;
};

Matrix matcl::imag(const Matrix& A0)
{
    Matrix A(A0);

    matcl::Matrix ret;
    details::eval_imag::make<const Matrix&>(A,ret);
    return ret;
};

Matrix matcl::imag(Matrix&& A0)
{
    Matrix A(std::move(A0));

    matcl::Matrix ret;
    details::eval_imag::make<const Matrix&>(A,ret);
    return ret;
};

Matrix matcl::real(const Matrix& A0)
{
    Matrix A(A0);
    
    matcl::Matrix ret;
    details::eval_real::make<const Matrix&>(A,ret);
    return ret;
};

Matrix matcl::real(Matrix&& A0)
{
    Matrix A(std::move(A0));
    
    matcl::Matrix ret;
    details::eval_real::make<const Matrix&>(A,ret);
    return ret;
};

Matrix matcl::eps(const Matrix& A0)
{
    Matrix A(A0);

    matcl::Matrix ret;
    details::eval_eps::make<const Matrix&>(A,ret);
    return ret;
};

Matrix matcl::eps(Matrix&& A0)
{
    Matrix A(std::move(A0));

    matcl::Matrix ret;
    details::eval_eps::make<const Matrix&>(A,ret);
    return ret;
};

Matrix matcl::abs(const Matrix& A0)
{
    Matrix A(A0);

    matcl::Matrix ret;
    details::eval_abs::make<const Matrix&>(A,ret);
    return ret;
};

Matrix matcl::abs(Matrix&& A0)
{
    Matrix A(std::move(A0));

    matcl::Matrix ret;
    details::eval_abs::make<const Matrix&>(A,ret);
    return ret;
};

Matrix matcl::abs2(const Matrix& A0)
{
    Matrix A(A0);

    matcl::Matrix ret;
    details::eval_abs2::make<const Matrix&>(A,ret);
    return ret;
};

Matrix matcl::abs2(Matrix&& A0)
{
    Matrix A(std::move(A0));

    matcl::Matrix ret;
    details::eval_abs2::make<const Matrix&>(A,ret);
    return ret;
};

Matrix matcl::arg(const Matrix& A0)
{
    Matrix A(A0);
    
    matcl::Matrix ret;
    details::eval_arg::make<const Matrix&>(A,ret);
    return ret;
};

Matrix matcl::arg(Matrix&& A0)
{
    Matrix A(std::move(A0));
    
    matcl::Matrix ret;
    details::eval_arg::make<const Matrix&>(A,ret);
    return ret;
};

Matrix matcl::angle(const Matrix& A0)
{
    Matrix A(A0);

    matcl::Matrix ret;
    details::eval_arg::make<const Matrix&>(A,ret);
    return ret;
};

Matrix matcl::angle(Matrix&& A0)
{
    Matrix A(std::move(A0));

    matcl::Matrix ret;
    details::eval_arg::make<const Matrix&>(A,ret);
    return ret;
};

Matrix matcl::cbrt(const Matrix& A)
{
    return eval_scalar_func(A, details::func_cbrt());
};

Matrix matcl::cbrt(Matrix&& A)
{
    return eval_scalar_func(std::move(A), details::func_cbrt());
};

Matrix matcl::expm1(const Matrix& A)
{
    return eval_scalar_func(A, details::func_expm1());
};

Matrix matcl::expm1(Matrix&& A)
{
    return eval_scalar_func(std::move(A), details::func_expm1());
};

Matrix matcl::expi(const Matrix& A)
{
    return eval_scalar_func(A, details::func_expi());
};

Matrix matcl::expi(Matrix&& A)
{
    return eval_scalar_func(std::move(A), details::func_expi());
};

Matrix matcl::signbit(const Matrix& A)
{
    return eval_scalar_func(A, details::func_signbit());
};

Matrix matcl::signbit(Matrix&& A)
{
    return eval_scalar_func(std::move(A), details::func_signbit());
};

Matrix matcl::logb(const Matrix& A)
{
    return eval_scalar_func(A, details::func_logb());
};

Matrix matcl::logb(Matrix&& A)
{
    return eval_scalar_func(std::move(A), details::func_logb());
};

Matrix matcl::ilogb(const Matrix& A)
{
    return eval_scalar_func(A, details::func_ilogb());
};

Matrix matcl::ilogb(Matrix&& A)
{
    return eval_scalar_func(std::move(A), details::func_ilogb());
};

Matrix matcl::exp2(const Matrix& A)
{
    return eval_scalar_func(A, details::func_exp2());
};

Matrix matcl::exp2(Matrix&& A)
{
    return eval_scalar_func(std::move(A), details::func_exp2());
};

Matrix matcl::exp10(const Matrix& A)
{
    return eval_scalar_func(A, details::func_exp10());
};

Matrix matcl::exp10(Matrix&& A)
{
    return eval_scalar_func(std::move(A), details::func_exp10());
};

Matrix matcl::sqrt(const Matrix &A0)
{
    Matrix A(A0);

    matcl::Matrix ret;
    details::eval_sqrt::make<const Matrix&>(A,ret);
    return ret;
};

Matrix matcl::sqrt_c(const Matrix &A0)
{
    Matrix A(A0);

    matcl::Matrix ret;
    details::eval_sqrt_c::make<const Matrix&>(A,ret);
    return ret;
};

Matrix matcl::sqrt(Matrix&& A0)
{
    Matrix A(std::move(A0));

    matcl::Matrix ret;
    details::eval_sqrt::make<const Matrix&>(A,ret);
    return ret;
};

Matrix matcl::sqrt_c(Matrix&& A0)
{
    Matrix A(std::move(A0));

    matcl::Matrix ret;
    details::eval_sqrt_c::make<const Matrix&>(A,ret);
    return ret;
};

Matrix matcl::sqrt1pm1(const Matrix &A0)
{
    Matrix A(A0);

    matcl::Matrix ret;
    details::eval_sqrt_1pm1::make<const Matrix&>(A,ret);
    return ret;
};

Matrix matcl::sqrt1pm1_c(const Matrix &A0)
{
    Matrix A(A0);

    matcl::Matrix ret;
    details::eval_sqrt_1pm1_c::make<const Matrix&>(A,ret);
    return ret;
};

Matrix matcl::sqrt1pm1(Matrix&& A0)
{
    Matrix A(std::move(A0));

    matcl::Matrix ret;
    details::eval_sqrt_1pm1::make<const Matrix&>(A,ret);
    return ret;
};

Matrix matcl::sqrt1pm1_c(Matrix&& A0)
{
    Matrix A(std::move(A0));

    matcl::Matrix ret;
    details::eval_sqrt_1pm1_c::make<const Matrix&>(A,ret);
    return ret;
};

Matrix matcl::log(const Matrix &A0)
{
    Matrix A(A0);

    matcl::Matrix ret;
    details::eval_log::make<const Matrix&>(A,ret);
    return ret;
};

Matrix matcl::log_c(const Matrix &A0)
{
    Matrix A(A0);

    matcl::Matrix ret;
    details::eval_log_c::make<const Matrix&>(A,ret);
    return ret;
};

Matrix matcl::log(Matrix&& A0)
{
    Matrix A(std::move(A0));

    matcl::Matrix ret;
    details::eval_log::make<const Matrix&>(A,ret);
    return ret;
};

Matrix matcl::log_c(Matrix&& A0)
{
    Matrix A(std::move(A0));

    matcl::Matrix ret;
    details::eval_log_c::make<const Matrix&>(A,ret);
    return ret;
};

Matrix matcl::log1p(const Matrix &A0)
{
    Matrix A(A0);

    matcl::Matrix ret;
    details::eval_log_1p::make<const Matrix&>(A,ret);
    return ret;
};

Matrix matcl::log1p_c(const Matrix &A0)
{
    Matrix A(A0);

    matcl::Matrix ret;
    details::eval_log_1p_c::make<const Matrix&>(A,ret);
    return ret;
};

Matrix matcl::log1p(Matrix&& A0)
{
    Matrix A(std::move(A0));

    matcl::Matrix ret;
    details::eval_log_1p::make<const Matrix&>(A,ret);
    return ret;
};

Matrix matcl::log1p_c(Matrix&& A0)
{
    Matrix A(std::move(A0));

    matcl::Matrix ret;
    details::eval_log_1p_c::make<const Matrix&>(A,ret);
    return ret;
};

Matrix matcl::exp(const Matrix &A0)
{
    Matrix A(A0);

    matcl::Matrix ret;
    details::eval_exp::make<const Matrix&>(A,ret);
    return ret;
};

Matrix matcl::exp(Matrix&& A0)
{
    Matrix A(std::move(A0));

    matcl::Matrix ret;
    details::eval_exp::make<const Matrix&>(A,ret);
    return ret;
};

Matrix matcl::log2(const Matrix &A0)
{
    Matrix A(A0);

    matcl::Matrix ret;
    details::eval_log2::make<const Matrix&>(A,ret);
    return ret;
};

Matrix matcl::log2_c(const Matrix &A0)
{
    Matrix A(A0);

    matcl::Matrix ret;
    details::eval_log2_c::make<const Matrix&>(A,ret);
    return ret;
};

Matrix matcl::log2(Matrix&& A0)
{
    Matrix A(std::move(A0));

    matcl::Matrix ret;
    details::eval_log2::make<const Matrix&>(A,ret);
    return ret;
};

Matrix matcl::log2_c(Matrix&& A0)
{
    Matrix A(std::move(A0));

    matcl::Matrix ret;
    details::eval_log2_c::make<const Matrix&>(A,ret);
    return ret;
};

Matrix matcl::log10(const Matrix &A0)
{
    Matrix A(A0);

    matcl::Matrix ret;
    details::eval_log10::make<const Matrix&>(A,ret);
    return ret;
};

Matrix matcl::log10_c(const Matrix &A0)
{
    Matrix A(A0);

    matcl::Matrix ret;
    details::eval_log10_c::make<const Matrix&>(A,ret);
    return ret;
};

Matrix matcl::log10(Matrix&& A0)
{
    Matrix A(std::move(A0));

    matcl::Matrix ret;
    details::eval_log10::make<const Matrix&>(A,ret);
    return ret;
};

Matrix matcl::log10_c(Matrix&& A0)
{
    Matrix A(std::move(A0));

    matcl::Matrix ret;
    details::eval_log10_c::make<const Matrix&>(A,ret);
    return ret;
};

Matrix matcl::floor(const Matrix &A0)
{
    Matrix A(A0);

    matcl::Matrix ret;
    details::eval_floor::make<const Matrix&>(A,ret);
    return ret;
};

Matrix matcl::floor(Matrix&& A0)
{
    Matrix A(std::move(A0));

    matcl::Matrix ret;
    details::eval_floor::make<const Matrix&>(A,ret);
    return ret;
};

Matrix matcl::ceil(const Matrix &A0)
{
    Matrix A(A0);

    matcl::Matrix ret;
    details::eval_ceil::make<const Matrix&>(A,ret);
    return ret;
};

Matrix matcl::ceil(Matrix&& A0)
{
    Matrix A(std::move(A0));

    matcl::Matrix ret;
    details::eval_ceil::make<const Matrix&>(A,ret);
    return ret;
};

Matrix matcl::round(const Matrix &A0)
{
    Matrix A(A0);

    matcl::Matrix ret;
    details::eval_round::make<const Matrix&>(A,ret);
    return ret;
};

Matrix matcl::round(Matrix&& A0)
{
    Matrix A(std::move(A0));

    matcl::Matrix ret;
    details::eval_round::make<const Matrix&>(A,ret);
    return ret;
};

Matrix matcl::trunc(const Matrix &A0)
{
    Matrix A(A0);

    matcl::Matrix ret;
    details::eval_trunc::make<const Matrix&>(A,ret);
    return ret;
};

Matrix matcl::trunc(Matrix&& A0)
{
    Matrix A(std::move(A0));

    matcl::Matrix ret;
    details::eval_trunc::make<const Matrix&>(A,ret);
    return ret;
};

Matrix matcl::ifloor(const Matrix &A0)
{
    Matrix A(A0);

    matcl::Matrix ret;
    details::eval_ifloor::make<const Matrix&>(A,ret);
    return ret;
};

Matrix matcl::ifloor(Matrix&& A0)
{
    Matrix A(std::move(A0));

    matcl::Matrix ret;
    details::eval_ifloor::make<const Matrix&>(A,ret);
    return ret;
};

Matrix matcl::iceil(const Matrix &A0)
{
    Matrix A(A0);

    matcl::Matrix ret;
    details::eval_iceil::make<const Matrix&>(A,ret);
    return ret;
};

Matrix matcl::iceil(Matrix&& A0)
{
    Matrix A(std::move(A0));

    matcl::Matrix ret;
    details::eval_iceil::make<const Matrix&>(A,ret);
    return ret;
};


Matrix matcl::iround(const Matrix &A0)
{
    Matrix A(A0);

    matcl::Matrix ret;
    details::eval_iround::make<const Matrix&>(A,ret);
    return ret;
};

Matrix matcl::iround(Matrix&& A0)
{
    Matrix A(std::move(A0));

    matcl::Matrix ret;
    details::eval_iround::make<const Matrix&>(A,ret);
    return ret;
};

Matrix matcl::itrunc(const Matrix &A0)
{
    Matrix A(A0);

    matcl::Matrix ret;
    details::eval_itrunc::make<const Matrix&>(A,ret);
    return ret;
};

Matrix matcl::itrunc(Matrix&& A0)
{
    Matrix A(std::move(A0));

    matcl::Matrix ret;
    details::eval_itrunc::make<const Matrix&>(A,ret);
    return ret;
};

Matrix matcl::sign(const Matrix &A0)
{
    Matrix A(A0);

    matcl::Matrix ret;
    details::eval_sign::make<const Matrix&>(A,ret);
    return ret;
};

Matrix matcl::sign(Matrix&& A0)
{
    Matrix A(std::move(A0));

    matcl::Matrix ret;
    details::eval_sign::make<const Matrix&>(A,ret);
    return ret;
};

Matrix matcl::isign(const Matrix &A0)
{
    Matrix A(A0);

    matcl::Matrix ret;
    details::eval_isign::make<const Matrix&>(A,ret);
    return ret;
};

Matrix matcl::isign(Matrix&& A0)
{
    Matrix A(std::move(A0));

    matcl::Matrix ret;
    details::eval_isign::make<const Matrix&>(A,ret);
    return ret;
};

Matrix matcl::sin(const Matrix &A0)
{
    Matrix A(A0);

    matcl::Matrix ret;
    details::eval_sin::make<const Matrix&>(A,ret);
    return ret;
};

Matrix matcl::sin(Matrix&& A0)
{
    Matrix A(std::move(A0));

    matcl::Matrix ret;
    details::eval_sin::make<const Matrix&>(A,ret);
    return ret;
};

Matrix matcl::cos(const Matrix &A0)
{
    Matrix A(A0);

    matcl::Matrix ret;
    details::eval_cos::make<const Matrix&>(A,ret);
    return ret;
};

Matrix matcl::cos(Matrix&& A0)
{
    Matrix A(std::move(A0));

    matcl::Matrix ret;
    details::eval_cos::make<const Matrix&>(A,ret);
    return ret;
};

Matrix matcl::tan(const Matrix &A0)
{
    Matrix A(A0);

    matcl::Matrix ret;
    details::eval_tan::make<const Matrix&>(A,ret);
    return ret;
};

Matrix matcl::tan(Matrix&& A0)
{
    Matrix A(std::move(A0));

    matcl::Matrix ret;
    details::eval_tan::make<const Matrix&>(A,ret);
    return ret;
};

Matrix matcl::cot(const Matrix &A0)
{
    Matrix A(A0);

    matcl::Matrix ret;
    details::eval_cot::make<const Matrix&>(A,ret);
    return ret;
};

Matrix matcl::cot(Matrix&& A0)
{
    Matrix A(std::move(A0));

    matcl::Matrix ret;
    details::eval_cot::make<const Matrix&>(A,ret);
    return ret;
};

Matrix matcl::sec(const Matrix &A0)
{
    Matrix A(A0);

    matcl::Matrix ret;
    details::eval_sec::make<const Matrix&>(A,ret);
    return ret;
};

Matrix matcl::sec(Matrix&& A0)
{
    Matrix A(std::move(A0));

    matcl::Matrix ret;
    details::eval_sec::make<const Matrix&>(A,ret);
    return ret;
};

Matrix matcl::csc(const Matrix &A0)
{
    Matrix A(A0);

    matcl::Matrix ret;
    details::eval_csc::make<const Matrix&>(A,ret);
    return ret;
};

Matrix matcl::csc(Matrix&& A0)
{
    Matrix A(std::move(A0));

    matcl::Matrix ret;
    details::eval_csc::make<const Matrix&>(A,ret);
    return ret;
};

Matrix matcl::asin(const Matrix &A0)
{
    Matrix A(A0);

    matcl::Matrix ret;
    details::eval_asin::make<const Matrix&>(A,ret);
    return ret;
};

Matrix matcl::asin_c(const Matrix &A0)
{
    Matrix A(A0);

    matcl::Matrix ret;
    details::eval_asin_c::make<const Matrix&>(A,ret);
    return ret;
};

Matrix matcl::asin(Matrix&& A0)
{
    Matrix A(std::move(A0));

    matcl::Matrix ret;
    details::eval_asin::make<const Matrix&>(A,ret);
    return ret;
};

Matrix matcl::asin_c(Matrix&& A0)
{
    Matrix A(std::move(A0));

    matcl::Matrix ret;
    details::eval_asin_c::make<const Matrix&>(A,ret);
    return ret;
};

Matrix matcl::acos(const Matrix &A0)
{
    Matrix A(A0);

    matcl::Matrix ret;
    details::eval_acos::make<const Matrix&>(A,ret);
    return ret;
};

Matrix matcl::acos_c(const Matrix &A0)
{
    Matrix A(A0);

    matcl::Matrix ret;
    details::eval_acos_c::make<const Matrix&>(A,ret);
    return ret;
};

Matrix matcl::acos(Matrix&& A0)
{
    Matrix A(std::move(A0));

    matcl::Matrix ret;
    details::eval_acos::make<const Matrix&>(A,ret);
    return ret;
};

Matrix matcl::acos_c(Matrix&& A0)
{
    Matrix A(std::move(A0));

    matcl::Matrix ret;
    details::eval_acos_c::make<const Matrix&>(A,ret);
    return ret;
};

Matrix matcl::atan(const Matrix &A0)
{
    Matrix A(A0);

    matcl::Matrix ret;
    details::eval_atan::make<const Matrix&>(A,ret);
    return ret;
};

Matrix matcl::atan(Matrix&& A0)
{
    Matrix A(std::move(A0));

    matcl::Matrix ret;
    details::eval_atan::make<const Matrix&>(A,ret);
    return ret;
};

Matrix matcl::acot(const Matrix &A0)
{
    Matrix A(A0);

    matcl::Matrix ret;
    details::eval_acot::make<const Matrix&>(A,ret);
    return ret;
};

Matrix matcl::acot(Matrix&& A0)
{
    Matrix A(std::move(A0));

    matcl::Matrix ret;
    details::eval_acot::make<const Matrix&>(A,ret);
    return ret;
};

Matrix matcl::asec(const Matrix &A0)
{
    Matrix A(A0);

    matcl::Matrix ret;
    details::eval_asec::make<const Matrix&>(A,ret);
    return ret;
};

Matrix matcl::asec_c(const Matrix &A0)
{
    Matrix A(A0);

    matcl::Matrix ret;
    details::eval_asec_c::make<const Matrix&>(A,ret);
    return ret;
};

Matrix matcl::asec(Matrix&& A0)
{
    Matrix A(std::move(A0));

    matcl::Matrix ret;
    details::eval_asec::make<const Matrix&>(A,ret);
    return ret;
};

Matrix matcl::asec_c(Matrix&& A0)
{
    Matrix A(std::move(A0));

    matcl::Matrix ret;
    details::eval_asec_c::make<const Matrix&>(A,ret);
    return ret;
};

Matrix matcl::acsc(const Matrix &A0)
{
    Matrix A(A0);

    matcl::Matrix ret;
    details::eval_acsc::make<const Matrix&>(A,ret);
    return ret;
};

Matrix matcl::acsc_c(const Matrix &A0)
{
    Matrix A(A0);

    matcl::Matrix ret;
    details::eval_acsc_c::make<const Matrix&>(A,ret);
    return ret;
};

Matrix matcl::acsc(Matrix&& A0)
{
    Matrix A(std::move(A0));

    matcl::Matrix ret;
    details::eval_acsc::make<const Matrix&>(A,ret);
    return ret;
};

Matrix matcl::acsc_c(Matrix&& A0)
{
    Matrix A(std::move(A0));

    matcl::Matrix ret;
    details::eval_acsc_c::make<const Matrix&>(A,ret);
    return ret;
};

Matrix matcl::sinh(const Matrix &A0)
{
    Matrix A(A0);

    matcl::Matrix ret;
    details::eval_sinh::make<const Matrix&>(A,ret);
    return ret;
};

Matrix matcl::sinh(Matrix&& A0)
{
    Matrix A(std::move(A0));

    matcl::Matrix ret;
    details::eval_sinh::make<const Matrix&>(A,ret);
    return ret;
};

Matrix matcl::cosh(const Matrix &A0)
{
    Matrix A(A0);

    matcl::Matrix ret;
    details::eval_cosh::make<const Matrix&>(A,ret);
    return ret;
};

Matrix matcl::cosh(Matrix&& A0)
{
    Matrix A(std::move(A0));

    matcl::Matrix ret;
    details::eval_cosh::make<const Matrix&>(A,ret);
    return ret;
};

Matrix matcl::tanh(const Matrix &A0)
{
    Matrix A(A0);

    matcl::Matrix ret;
    details::eval_tanh::make<const Matrix&>(A,ret);
    return ret;
};

Matrix matcl::tanh(Matrix&& A0)
{
    Matrix A(std::move(A0));

    matcl::Matrix ret;
    details::eval_tanh::make<const Matrix&>(A,ret);
    return ret;
};

Matrix matcl::coth(const Matrix &A0)
{
    Matrix A(A0);

    matcl::Matrix ret;
    details::eval_coth::make<const Matrix&>(A,ret);
    return ret;
};

Matrix matcl::coth(Matrix&& A0)
{
    Matrix A(std::move(A0));

    matcl::Matrix ret;
    details::eval_coth::make<const Matrix&>(A,ret);
    return ret;
};

Matrix matcl::sech(const Matrix &A0)
{
    Matrix A(A0);

    matcl::Matrix ret;
    details::eval_sech::make<const Matrix&>(A,ret);
    return ret;
};

Matrix matcl::sech(Matrix&& A0)
{
    Matrix A(std::move(A0));

    matcl::Matrix ret;
    details::eval_sech::make<const Matrix&>(A,ret);
    return ret;
};

Matrix matcl::csch(const Matrix &A0)
{
    Matrix A(A0);

    matcl::Matrix ret;
    details::eval_csch::make<const Matrix&>(A,ret);
    return ret;
};

Matrix matcl::csch(Matrix&& A0)
{
    Matrix A(std::move(A0));

    matcl::Matrix ret;
    details::eval_csch::make<const Matrix&>(A,ret);
    return ret;
};

Matrix matcl::asinh(const Matrix &A0)
{
    Matrix A(A0);

    matcl::Matrix ret;
    details::eval_asinh::make<const Matrix&>(A,ret);
    return ret;
};

Matrix matcl::asinh(Matrix&& A0)
{
    Matrix A(std::move(A0));

    matcl::Matrix ret;
    details::eval_asinh::make<const Matrix&>(A,ret);
    return ret;
};

Matrix matcl::acosh(const Matrix &A0)
{
    Matrix A(A0);

    matcl::Matrix ret;
    details::eval_acosh::make<const Matrix&>(A,ret);
    return ret;
};

Matrix matcl::acosh_c(const Matrix &A0)
{
    Matrix A(A0);

    matcl::Matrix ret;
    details::eval_acosh_c::make<const Matrix&>(A,ret);
    return ret;
};

Matrix matcl::acosh(Matrix&& A0)
{
    Matrix A(std::move(A0));

    matcl::Matrix ret;
    details::eval_acosh::make<const Matrix&>(A,ret);
    return ret;
};

Matrix matcl::acosh_c(Matrix&& A0)
{
    Matrix A(std::move(A0));

    matcl::Matrix ret;
    details::eval_acosh_c::make<const Matrix&>(A,ret);
    return ret;
};

Matrix matcl::atanh(const Matrix &A0)
{
    Matrix A(A0);

    matcl::Matrix ret;
    details::eval_atanh::make<const Matrix&>(A,ret);
    return ret;
};

Matrix matcl::atanh_c(const Matrix &A0)
{
    Matrix A(A0);

    matcl::Matrix ret;
    details::eval_atanh_c::make<const Matrix&>(A,ret);
    return ret;
};

Matrix matcl::atanh(Matrix&& A0)
{
    Matrix A(std::move(A0));

    matcl::Matrix ret;
    details::eval_atanh::make<const Matrix&>(A,ret);
    return ret;
};

Matrix matcl::atanh_c(Matrix&& A0)
{
    Matrix A(std::move(A0));

    matcl::Matrix ret;
    details::eval_atanh_c::make<const Matrix&>(A,ret);
    return ret;
};

Matrix matcl::acoth(const Matrix &A0)
{
    Matrix A(A0);

    matcl::Matrix ret;
    details::eval_acoth::make<const Matrix&>(A,ret);
    return ret;
};

Matrix matcl::acoth_c(const Matrix &A0)
{
    Matrix A(A0);

    matcl::Matrix ret;
    details::eval_acoth_c::make<const Matrix&>(A,ret);
    return ret;
};

Matrix matcl::acoth(Matrix&& A0)
{
    Matrix A(std::move(A0));

    matcl::Matrix ret;
    details::eval_acoth::make<const Matrix&>(A,ret);
    return ret;
};

Matrix matcl::acoth_c(Matrix&& A0)
{
    Matrix A(std::move(A0));

    matcl::Matrix ret;
    details::eval_acoth_c::make<const Matrix&>(A,ret);
    return ret;
};

Matrix matcl::asech(const Matrix &A0)
{
    Matrix A(A0);

    matcl::Matrix ret;
    details::eval_asech::make<const Matrix&>(A,ret);
    return ret;
};

Matrix matcl::asech_c(const Matrix &A0)
{
    Matrix A(A0);

    matcl::Matrix ret;
    details::eval_asech_c::make<const Matrix&>(A,ret);
    return ret;
};

Matrix matcl::asech(Matrix&& A0)
{
    Matrix A(std::move(A0));

    matcl::Matrix ret;
    details::eval_asech::make<const Matrix&>(A,ret);
    return ret;
};

Matrix matcl::asech_c(Matrix&& A0)
{
    Matrix A(std::move(A0));

    matcl::Matrix ret;
    details::eval_asech_c::make<const Matrix&>(A,ret);
    return ret;
};

Matrix matcl::acsch(const Matrix &A0)
{
    Matrix A(A0);

    matcl::Matrix ret;
    details::eval_acsch::make<const Matrix&>(A,ret);
    return ret;
};

Matrix matcl::acsch(Matrix&& A0)
{
    Matrix A(std::move(A0));

    matcl::Matrix ret;
    details::eval_acsch::make<const Matrix&>(A,ret);
    return ret;
};

matcl::Matrix matcl::is_true(const Matrix& A0)
{
    Matrix A(A0);

    matcl::Matrix ret;
    details::eval_is_true::make<const Matrix&>(A,ret);
    return ret;
};

Matrix matcl::is_true(Matrix&& A0)
{
    Matrix A(std::move(A0));

    matcl::Matrix ret;
    details::eval_is_true::make<const Matrix&>(A,ret);
    return ret;
};

matcl::Matrix matcl::neg(const Matrix& A0)
{
    Matrix A(A0);

    matcl::Matrix ret;
    details::eval_neg::make<const Matrix&>(A,ret);
    return ret;
};

Matrix matcl::neg(Matrix&& A0)
{
    Matrix A(std::move(A0));

    matcl::Matrix ret;
    details::eval_neg::make<const Matrix&>(A,ret);
    return ret;
};

Matrix matcl::uminus(const Matrix& A0)
{
    Matrix A(A0);

    matcl::Matrix ret;
    details::eval_sgn_min::make<const Matrix&>(A,ret);
    return ret;
};

Matrix matcl::uminus(Matrix&& A0)
{
    Matrix A(std::move(A0));

    matcl::Matrix ret;
    details::eval_sgn_min::make<const Matrix&>(A,ret);
    return ret;
};

Matrix matcl::invs(const Matrix& A0)
{
    Matrix A(A0);

    matcl::Matrix ret;
    details::eval_inv::make<const Matrix&>(A,ret);
    return ret;
};

Matrix matcl::invs(Matrix&& A0)
{
    Matrix A(std::move(A0));

    matcl::Matrix ret;
    details::eval_inv::make<const Matrix&>(A,ret);
    return ret;
};

Matrix matcl::conj(const Matrix& A0)
{
    Matrix A(A0);

    matcl::Matrix ret;
    details::eval_conj::make<const Matrix&>(A,ret);
    return ret;
};

Matrix matcl::conj(Matrix&& A0)
{
    Matrix A(std::move(A0));

    matcl::Matrix ret;
    details::eval_conj::make<const Matrix&>(A,ret);
    return ret;
};

};

#pragma warning(pop)