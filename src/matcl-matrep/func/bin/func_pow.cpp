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
#include "matcl-internals/func/test_functor.h" 
#include "matcl-matrep/func/raw/bin/raw_func_pow.h"
#include "matcl-matrep/func/raw/bin/raw_func_op_helpers.h"

#pragma warning( push )
#pragma warning(disable:4127)	// conditional expression is constant

namespace matcl { namespace details
{

namespace md = matcl::details;
namespace mr = matcl::raw;
namespace mrd = matcl::raw::details;
namespace mdyf  = matcl::dynamic::functions;

template<class To, class From>
struct convert_val
{
    static To eval(const From& v)
    {
        return To(v);
    };
};

template<>
struct convert_val<Float_complex,Integer>
{
    static Float_complex eval(const Integer& v)
    {
        return Float_complex(static_cast<Float>(v));
    };
};

template<class T, class M>
struct pow_helper_scal_mat
{    
    static void eval(matcl::Matrix& ret, const T& val, const M& mat, bool complex_ver);
};

template<class M,class T>
struct pow_helper_mat_scal
{    
    static void eval(matcl::Matrix& ret, const M& mat,const T& scal, bool complex_ver);
};

template<class M1,class M2>
struct pow_helper_mat_mat
{    
    static void eval(matcl::Matrix& ret, const M1& mat1,const M2& mat2, bool complex_ver);
};

template<class M1, class M2>
struct function_op_pow
{
    static void eval(matcl::Matrix& ret, const M1& mat1, const M2& mat2, bool complex_ver)
    {
        using V1                = typename M1::value_type;
        using V2                = typename M2::value_type;
        using V1_real           = typename details::real_type<V1>::type;
        using V2_real           = typename details::real_type<V2>::type;
        using val_real          = typename md::select_if
                                <   std::is_same<V2_real,Integer>::value, 
                                    typename md::unify_types<V1_real,Float>::type, 
                                    typename md::unify_types<V1_real,V2_real>::type
                                >::type;

        static const bool is_obj = std::is_same<V1,Object>::value 
                                   || std::is_same<V2,Object>::value;

        using ret_ti_type       = typename select_ti_type<is_obj>::type;
        using ret_ti_value      = typename select_ti_type<is_obj>::value_type;

        typename ti::get_ti_type<M1>::type ti_1 = ti::get_ti(mat1);
        typename ti::get_ti_type<M2>::type ti_2 = ti::get_ti(mat2);

        if (mat1.rows() == 1 && mat1.cols() == 1)
        {
            V1 val(mat1(1,1));

            if (mat2.rows() == 1 && mat2.cols() == 1)
            {
                V2 val2(mat2(1,1));

                if (complex_ver)
                    ret = mrd::pow_c_helper<V1,V2>::eval(val,val2);
                else
                    ret = mrd::pow_helper<V1,V2>::eval(val,val2);
                
                return;
            };

            if (mrd::is_one(val))
            {
                Integer m = mat2.rows();
                Integer n = mat2.cols();

                if (is_obj)
                {
                    ti::ti_object ret_ti;
                    
                    if (complex_ver == true)
                        ret_ti  = ti::get_return_ti<ti::ti_object>(mdyf::pow_c::eval(), ti_1,ti_2);
                    else
                        ret_ti  = ti::get_return_ti<ti::ti_object>(mdyf::pow::eval(), ti_1,ti_2);

                    ret = ones(ret_ti,m,n);
                    return;
                };

                ret_ti_type ret_ti;
                val_real one    = md::one_value<val_real>(ret_ti);
                ret             = mrd::create_matrix<val_real>(ret_ti,one,m,n);
                return;
            };

            V1_real val_im(mrd::imag_helper<V1>::eval(val));

            if (mrd::is_zero(val_im))
            {
                V1_real val_re(mrd::real_helper<V1>::eval(val));
                return details::pow_helper_scal_mat<V1_real,M2>::eval(ret, val_re, mat2, complex_ver);
            }
            else
            {
                return details::pow_helper_scal_mat<V1,M2>::eval(ret, val, mat2, complex_ver);
            };
        };

        if (mat2.rows() == 1 && mat2.cols() == 1)
        {			
            V2 val(mat2(1,1));

            if (mrd::is_one(val))
            {
                ret_ti_type ret_ti;
                
                if (complex_ver)
                    ret_ti = ti::get_return_ti<ret_ti_type>(mdyf::pow_c::eval(), ti_1,ti_2);
                else
                    ret_ti = ti::get_return_ti<ret_ti_type>(mdyf::pow::eval(), ti_1,ti_2);

                ret = correct_int_val<false,is_obj>::eval(ret_ti,Matrix(mat1,false));
                return;
            };

            if (mrd::is_zero(val))
            {
                Integer m = mat1.rows();
                Integer n = mat1.cols();

                if (is_obj)
                {
                    ti::ti_object ret_ti;

                    if (complex_ver)
                        ret_ti = ti::get_return_ti<ti::ti_object>(mdyf::pow_c::eval(), ti_1,ti_2);
                    else
                        ret_ti = ti::get_return_ti<ti::ti_object>(mdyf::pow::eval(), ti_1,ti_2);
                    
                    ret = ones(ret_ti,m,n);
                    return;
                };

                ret_ti_type ret_ti;
                val_real one    = md::one_value<val_real>(ret_ti);
                ret             = mrd::create_matrix<val_real>(ret_ti,one,m,n);
                return;
            };

            V2_real val_im(mrd::imag_helper<V2>::eval(val));

            if (mrd::is_zero(val_im))
            {
                V2_real val_re(mrd::real_helper<V2>::eval(val));
                return details::pow_helper_mat_scal<M1,V2_real>::eval(ret,mat1,val_re,complex_ver);
            }
            else
            {
                return details::pow_helper_mat_scal<M1,V2>::eval(ret,mat1,val,complex_ver);
            };
        };

        Integer m1 = mat1.rows();
        Integer n1 = mat1.cols();
        Integer m2 = mat2.rows();
        Integer n2 = mat2.cols();

        error::check_eeop(m1,n1,m2,n2);

        if (mat2.nnz() == 0)
        {
            Integer m = mat1.rows();
            Integer n = mat1.cols();

            if (is_obj)
            {
                ti::ti_object ret_ti;

                if (complex_ver)
                    ret_ti = ti::get_return_ti<ti::ti_object>(mdyf::pow_c::eval(), ti_1,ti_2);
                else
                    ret_ti = ti::get_return_ti<ti::ti_object>(mdyf::pow::eval(), ti_1,ti_2);

                ret = ones(ret_ti,m,n);
                return;
            };

            ret_ti_type ret_ti;
            val_real one    = md::one_value<val_real>(ret_ti);
            ret             = mrd::create_matrix<val_real>(ret_ti,one,m,n);
            return;
        };

        return details::pow_helper_mat_mat<M1,M2>::eval(ret, mat1,mat2,complex_ver);
    }
};

template<class T1, class T2>
struct val_type_corrector_int1_tr
{
    using val_1     = typename mr::get_value_type<T1>::type;
    using val_2     = typename mr::get_value_type<T2>::type;

    using str_1     = typename mr::get_struct_type<T1>::type;
    using str_2     = typename mr::get_struct_type<T2>::type;
    
    using val_ret_1 = typename md::select_if
                    <   std::is_same<val_2,Integer>::value, 
                        val_1, 
                        typename mr::ret_value_type_int<val_1,val_2>::type
                    >::type;

    using val_ret_2 = typename md::select_if
                    <   std::is_same<val_2,Integer>::value, 
                        Integer, 
                        typename mr::ret_value_type_int<val_2,val_1>::type
                    >::type;

    using type_1    = typename mr::ret_matrix_type<val_ret_1,str_1>::type;
    using type_2    = typename mr::ret_matrix_type<val_ret_2,str_2>::type;
};

template<class T1, class T2>
struct val_type_corrector_int1  : public mr::val_type_corrector_impl<T1,T2,
                                    typename val_type_corrector_int1_tr<T1,T2>::type_1,
                                    typename val_type_corrector_int1_tr<T1,T2>::type_2>
                                , public val_type_corrector_int1_tr<T1,T2>
{};

struct eval_pow : public extract_type2_switch<void,eval_pow,val_type_corrector_int1>
{
    template<class T1, class T2>
    static void eval_mat_mat(const T1& A, const T2& B, matcl::Matrix& ret)
    {
        return function_op_pow<T1,T2>::eval(ret, A, B, false);
    };

    template<class T1, class T2>
    static void eval_scal_scal(const T1& A, const T2& B, matcl::Matrix& ret)
    {
        ret = mrd::pow_helper<T1,T2>::eval(A,B);
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

struct eval_pow_c : public extract_type2_switch<void,eval_pow_c,val_type_corrector_int1>
{
    template<class T1, class T2>
    static void eval_mat_mat(const T1& A, const T2& B, matcl::Matrix& ret)
    {
        return function_op_pow<T1,T2>::eval(ret, A, B, true);
    };

    template<class T1, class T2>
    static void eval_scal_scal(const T1& A, const T2& B, matcl::Matrix& ret)
    {
        ret = mrd::pow_c_helper<T1,T2>::eval(A,B);
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

template<bool iso,class T, class M>
struct pow_helper_scal_mat_impl
{
    static void eval(matcl::Matrix& ret, const T& val, const M& mat, bool complex_ver)
    {
        using T2            = typename M::value_type;
        using complex_type  = typename md::select_if
                            <   std::is_same<T2,Integer>::value, 
                                typename md::unify_types<T,Float_complex>::type, 
                                typename md::unify_types2<T,T2,Float_complex>::type
                            >::type;
        
        static const bool isc = md::is_complex<T>::value 
                                || md::is_complex<M>::value;	

        if (complex_ver == false)
            return mrd::pow_helper_scal_mat_impl<T,M>::eval(ret,val,mat);

        if (isc || mrd::is_geq_zero(val))
            return mrd::pow_helper_scal_mat_impl<T,M>::eval(ret,val,mat);

        test_is_int test_object;
        if (test_range<T2,typename M::struct_type,test_is_int>::eval(mat, test_object))
            return mrd::pow_helper_scal_mat_impl<T,M>::eval(ret,val,mat);

        return mrd::pow_helper_scal_mat_impl<complex_type,M>::eval(ret,complex_type(val),mat);
    };
};

template<class T, class M>
struct pow_helper_scal_mat_impl<true,T,M>
{
    static void eval(matcl::Matrix& ret, const T& val, const M& mat, bool complex_ver)
    {
        if (complex_ver == true)
            return mrd::pow_helper_scal_mat_impl<T,M>::eval_c(ret, val, mat);
        else
            return mrd::pow_helper_scal_mat_impl<T,M>::eval_c(ret, val, mat);
    };
};

template<class T, class M>
void pow_helper_scal_mat<T,M>::eval(matcl::Matrix& ret, const T& val, const M& mat, bool complex_ver)
{
    static const bool iso = std::is_same<T,Object>::value 
                            || std::is_same<M,Object>::value;	
    return pow_helper_scal_mat_impl<iso,T,M>::eval(ret,val,mat,complex_ver);
};

template<bool iso,class M, class T>
struct pow_helper_mat_scal_impl
{
    static void eval(matcl::Matrix& ret, const M& mat,const T& val, bool complex_ver)
    {
        using T1            = typename M::value_type;
        using complex_type  = typename md::select_if
                            <   std::is_same<T,Integer>::value, 
                                typename md::unify_types<T1,Float_complex>::type, 
                                typename md::unify_types2<T1,T,Float_complex>::type
                            >::type;

        static const bool isc = md::is_complex<T>::value 
                                || md::is_complex<M>::value;

        if (complex_ver == false)
            return mrd::pow_helper_mat_scal_impl<M,T>::eval(ret, mat,val);

        if (isc || is_conv_to_integer<T>::eval(val))
            return mrd::pow_helper_mat_scal_impl<M,T>::eval(ret, mat,val);

        test_geq_zero test_object;
        if (test_range<T1,typename M::struct_type,test_geq_zero>::eval(mat, test_object))
            return mrd::pow_helper_mat_scal_impl<M,T>::eval(ret, mat,val);

        if (std::is_same<T,Integer>::value)
        {
            using matrix_type = raw::Matrix<complex_type,typename M::struct_type>;
            return mrd::pow_helper_mat_scal_impl<matrix_type,T>
                ::eval(ret,raw::converter<matrix_type,M>::eval(mat),val);
        }
        else
        {
            return mrd::pow_helper_mat_scal_impl<M,complex_type>
                ::eval(ret,mat,convert_val<complex_type,T>::eval(val));
        };
    };
};

template<class M, class T>
struct pow_helper_mat_scal_impl<true,M,T>
{
    static void eval(matcl::Matrix& ret, const M& mat,const T& val, bool complex_ver)
    {
        if (complex_ver == true)
            return mrd::pow_helper_mat_scal_impl<M,T>::eval_c(ret, mat, val);
        else
            return mrd::pow_helper_mat_scal_impl<M,T>::eval(ret, mat, val);
    };
};

template<class M,class T>
void pow_helper_mat_scal<M,T>::eval(matcl::Matrix& ret, const M& mat,const T& val, bool complex_ver)
{
    static const bool iso = std::is_same<T,Object>::value 
                            || std::is_same<M,Object>::value;	
    return pow_helper_mat_scal_impl<iso,M,T>::eval(ret, mat, val, complex_ver);
};

template<bool iso,class M1, class M2>
struct pow_helper_mat_mat_impl
{
    static void eval(matcl::Matrix& ret, const M1& mat1,const M2& mat2, bool complex_ver)
    {
        using val_1         = typename M1::value_type;
        using val_2         = typename M2::value_type;
        using complex_type  = typename md::select_if
                            <   std::is_same<val_2,Integer>::value, 
                                typename md::unify_types<val_1,Float_complex>::type, 
                                typename md::unify_types2<val_1,val_2,Float_complex>::type
                            >::type;

        static const bool isc = md::is_complex<val_1>::value || md::is_complex<val_2>::value;

        if (complex_ver == false)
            return mrd::pow_helper_mat_mat_impl<M1,M2>::eval(ret,mat1,mat2);

        if (isc || test_range2<M1,M2, mrd::test_range_pow<val_1,val_2>>::eval(mat1,mat2))
            return mrd::pow_helper_mat_mat_impl<M1,M2>::eval(ret,mat1,mat2);
        
        if (mat1.nnz() < mat2.nnz() || std::is_same<val_2,Integer>::value == true)
        {
            //dont convert integer powers to complex
            using matrix_type = raw::Matrix<complex_type,typename M1::struct_type>;
            return mrd::pow_helper_mat_mat_impl<matrix_type,M2>
                ::eval(ret,raw::converter<matrix_type,M1>::eval(mat1),mat2);
        }
        else
        {
            using matrix_type = raw::Matrix<complex_type,typename M2::struct_type>;
            return mrd::pow_helper_mat_mat_impl<M1,matrix_type>
                ::eval(ret,mat1,raw::converter<matrix_type,M2>::eval(mat2));
        };
    };
};

template<class M1, class M2>
struct pow_helper_mat_mat_impl<true,M1,M2>
{
    static void eval(matcl::Matrix& ret, const M1& mat1,const M2& mat2, bool complex_ver)
    {
        if (complex_ver == true)
            return mrd::pow_helper_mat_mat_impl<M1,M2>::eval_c(ret,mat1,mat2);
        else
            return mrd::pow_helper_mat_mat_impl<M1,M2>::eval(ret,mat1,mat2);
    };
};

template<class M1,class M2>
void pow_helper_mat_mat<M1,M2>::eval(matcl::Matrix& ret, const M1& mat1,const M2& mat2, bool complex_ver)
{
    using val_1 = typename M1::value_type;
    using val_2 = typename M2::value_type;
    static const bool iso = std::is_same<val_1,Object>::value 
                            || std::is_same<val_2,Object>::value;	

    return pow_helper_mat_mat_impl<iso,M1,M2>::eval(ret,mat1,mat2,complex_ver);
};

}}

matcl::Matrix matcl::pow(const Matrix& A0, const Matrix& B0)
{
    Matrix A(A0);
    Matrix B(B0);

    matcl::Matrix ret;
    details::eval_pow::make(A,B,ret);
    return ret;
};

matcl::Matrix matcl::pow_c(const Matrix& A0, const Matrix& B0)
{
    Matrix A(A0);
    Matrix B(B0);

    matcl::Matrix ret;
    details::eval_pow_c::make(A,B,ret);	
    return ret;
};

matcl::Matrix matcl::pow(Matrix&& A0, const Matrix& B0)
{
    Matrix A(std::move(A0));
    Matrix B(B0);

    matcl::Matrix ret;
    details::eval_pow::make(A,B,ret);
    return ret;
};

matcl::Matrix matcl::pow_c(Matrix&& A0, const Matrix& B0)
{
    Matrix A(std::move(A0));
    Matrix B(B0);

    matcl::Matrix ret;
    details::eval_pow_c::make(A,B,ret);	
    return ret;
};

matcl::Matrix matcl::pow(Matrix&& A0, Matrix&& B0)
{
    Matrix A(std::move(A0));
    Matrix B(std::move(B0));

    matcl::Matrix ret;
    details::eval_pow::make(A,B,ret);
    return ret;
};

matcl::Matrix matcl::pow_c(Matrix&& A0, Matrix&& B0)
{
    Matrix A(std::move(A0));
    Matrix B(std::move(B0));

    matcl::Matrix ret;
    details::eval_pow_c::make(A,B,ret);	
    return ret;
};

matcl::Matrix matcl::pow(const Matrix& A0, Matrix&& B0)
{
    Matrix A(A0);
    Matrix B(std::move(B0));

    matcl::Matrix ret;
    details::eval_pow::make(A,B,ret);
    return ret;
};

matcl::Matrix matcl::pow_c(const Matrix& A0, Matrix&& B0)
{
    Matrix A(A0);
    Matrix B(std::move(B0));

    matcl::Matrix ret;
    details::eval_pow_c::make(A,B,ret);	
    return ret;
};

#pragma warning( pop )
