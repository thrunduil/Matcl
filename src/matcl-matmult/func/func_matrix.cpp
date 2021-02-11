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

#include "matcl-matrep/general/config.h"
#include "matcl-matrep/lib_functions/func_matrix.h"
#include "matcl-matrep/details/extract_type_switch.h"
#include "matcl-matmult/func/raw/raw_shprod.h"
#include "matcl-matmult/func/raw/raw_scale.h"
#include "matcl-scalar/details/matfunc_helpers.h"
#include "matcl-matrep/lib_functions/matrix_gen.h"

#pragma warning(push)
#pragma warning(disable: 4127)  //conditional expression is constant

namespace matcl { namespace details
{

namespace mdyf  = matcl::dynamic::functions;

struct eval_SYMPROD : public extract_type_switch<void,eval_SYMPROD,true> 
{
    template<class T, class S>
    static void eval(const Matrix& , const T& mat, Matrix& ret, const S& trans)
    {
        return mrd::symprod_helper<T>::eval_symprod(ret, mat, trans);
    };
    template<class T, class S>
    static void eval_scalar(const Matrix&, const T& mat, Matrix& ret, const S&)
    {
        ret = mrd::mmul_helper<T,T>::eval(mat,mat);
    };
};

struct eval_SYMSUM : public extract_type_switch<void,eval_SYMSUM,true> 
{
    template<class T>
    static void eval(const Matrix& , const T& mat, Matrix& ret)
    {
        return mrd::symprod_helper<T>::eval_symsum(ret, mat);
    };
    template<class T>
    static void eval_scalar(const Matrix&, const T& mat, Matrix& ret)
    {
        ret = mrd::plus_helper<T,T>::eval(mat,mat);
    };
};

struct eval_HERPROD : public extract_type_switch<void,eval_HERPROD,true> 
{
    template<class T, class S>
    static void eval(const Matrix& , const T& mat, Matrix& ret, const S& trans)
    {
        return mrd::symprod_helper<T>::eval_herprod(ret, mat, trans);
    };
    template<class T, class S>
    static void eval_scalar(const Matrix&, const T& mat, Matrix& ret, const S&)
    {
        ret = mrd::abs2_helper<T>::eval(mat);
    };
};

struct eval_HERSUM : public extract_type_switch<void,eval_HERSUM,true> 
{
    template<class T>
    static void eval(const Matrix& , const T& mat, Matrix& ret)
    {
        return mrd::symprod_helper<T>::eval_hersum(ret, mat);
    };
    template<class T>
    static void eval_scalar(const Matrix&, const T& mat, Matrix& ret)
    {
        using VR    = typename md::real_type<T>::type;
        ret = mrd::mmul_helper<VR,VR>::eval(real(mat), VR(2.0));
    };
};

struct eval_scale_rows : public extract_type_switch<void,eval_scale_rows,true> 
{
    template<class Mat>
    static void eval_impl(Matrix& ret, const Matrix& A, const Matrix& D)
    {
        using V             = typename Mat::value_type;
        using S             = typename Mat::struct_type;
        using VR            = typename md::real_type<V>::type;
        using Mat_D         = raw::Matrix<V,struct_dense>;
        using Mat_DR        = raw::Matrix<VR,struct_dense>;

        using ti_ret_type   = typename ti::get_ti_type<V>::type;
        ti_ret_type ret_ti  = ti::get_return_ti<ti_ret_type>(mdyf::op_mul::eval(), 
                                                            A.get_type(), D.get_type());

        value_code vcd      = D.get_value_code();
        value_code vc       = matrix_traits::value_code<V>::value;
        bool is_real        = matrix_traits::is_real(vcd);

        if (std::is_same<V,Object>::value && ret_ti != A.get_type())
        {
            Matrix Ac0      = convert_object(A, ti::convert_ti_object<V>(ret_ti));
            Matrix Ac       = convert(Ac0, Mat::matrix_code);
            Matrix Dc       = convert(D, Mat_D::matrix_code);
            const Mat& Am   = Ac.get_impl<Mat>();
            const Mat_D& Dm = Dc.get_impl<Mat_D>();

            return mrd::scaling_helper<Mat>::eval_rows(ret, Am, Dm);
        };

        if (vc == vcd || is_real == false)
        {
            Matrix Ac           = convert(A, Mat::matrix_code);
            Matrix Dc           = convert(D, Mat_D::matrix_code);

            const Mat& Am       = Ac.get_impl<Mat>();
            const Mat_D& Dm     = Dc.get_impl<Mat_D>();

            return mrd::scaling_helper<Mat>::eval_rows(ret, Am, Dm);
        }
        else
        {            
            Matrix Ac           = convert(A, Mat::matrix_code);
            Matrix Dc           = convert(D, Mat_DR::matrix_code);

            const Mat& Am       = Ac.get_impl<Mat>();
            const Mat_DR& Dm    = Dc.get_impl<Mat_DR>();

            return mrd::scaling_real_helper<Mat>::eval_rows(ret, Am, Dm);
        };        
    };

    template<class T>
    static void eval(const Matrix& , const T&, Matrix& ret, const Matrix& A, const Matrix& D)
    {
        return eval_impl<T>(ret, A, D);
    };

    template<class T>
    static void eval_scalar(const Matrix&, const T&, Matrix& ret, const Matrix& A, const Matrix& D)
    {
        using Mat   = raw::Matrix<T,struct_dense>;
        return eval_impl<Mat>(ret, A, D);
    };
};

struct eval_scale_cols : public extract_type_switch<void,eval_scale_cols,true> 
{
    template<class Mat>
    static void eval_impl(Matrix& ret, const Matrix& A, const Matrix& D)
    {
        using V             = typename Mat::value_type;
        using S             = typename Mat::struct_type;
        using VR            = typename md::real_type<V>::type;
        using Mat_D         = raw::Matrix<V,struct_dense>;
        using Mat_DR        = raw::Matrix<VR,struct_dense>;

        using ti_ret_type   = typename ti::get_ti_type<V>::type;
        ti_ret_type ret_ti  = ti::get_return_ti<ti_ret_type>(mdyf::op_mul::eval(), 
                                                            A.get_type(), D.get_type());

        value_code vcd      = D.get_value_code();
        value_code vc       = matrix_traits::value_code<V>::value;
        bool is_real        = matrix_traits::is_real(vcd);

        if (std::is_same<V,Object>::value && ret_ti != A.get_type())
        {
            Matrix Ac       = convert_object(A, ti::convert_ti_object<V>(ret_ti));
            Matrix Ac2      = convert(Ac, Mat::matrix_code);
            Matrix Dc       = convert(D, Mat_D::matrix_code);

            const Mat& Am   = Ac2.get_impl<Mat>();
            const Mat_D& Dm = Dc.get_impl<Mat_D>();

            return mrd::scaling_helper<Mat>::eval_cols(ret, Am, Dm);
        };

        if (vc == vcd || is_real == false)
        {
            Matrix Ac       = convert(A, Mat::matrix_code);
            Matrix Dc       = convert(D, Mat_D::matrix_code);

            const Mat& Am       = Ac.get_impl<Mat>();
            const Mat_D& Dm     = Dc.get_impl<Mat_D>();

            return mrd::scaling_helper<Mat>::eval_cols(ret, Am, Dm);
        }
        else
        {            
            Matrix Ac       = convert(A, Mat::matrix_code);
            Matrix Dc       = convert(D, Mat_DR::matrix_code);

            const Mat& Am       = Ac.get_impl<Mat>();
            const Mat_DR& Dm    = Dc.get_impl<Mat_DR>();

            return mrd::scaling_real_helper<Mat>::eval_cols(ret, Am, Dm);
        };        
    };

    template<class T>
    static void eval(const Matrix& , const T&, Matrix& ret, const Matrix& A, const Matrix& D)
    {
        return eval_impl<T>(ret, A, D);
    };

    template<class T>
    static void eval_scalar(const Matrix&, const T&, Matrix& ret, const Matrix& A, const Matrix& D)
    {
        using Mat   = raw::Matrix<T,struct_dense>;
        return eval_impl<Mat>(ret, A, D);
    };
};

struct eval_scale_rowscols : public extract_type_switch<void,eval_scale_rowscols,true> 
{
    template<class Mat>
    static void eval_impl(Matrix& ret, const Matrix& A, const Matrix& Dr, const Matrix& Dc)
    {
        using V             = typename Mat::value_type;
        using S             = typename Mat::struct_type;
        using VR            = typename md::real_type<V>::type;
        using Mat_D         = raw::Matrix<V,struct_dense>;
        using Mat_DR        = raw::Matrix<VR,struct_dense>;

        using ti_ret_type   = typename ti::get_ti_type<V>::type;
        ti_ret_type ret_ti  = ti::get_return_ti<ti_ret_type>(mdyf::op_mul::eval(), 
                                                            A.get_type(), Dr.get_type());
        ret_ti              = ti::get_return_ti<ti_ret_type>(mdyf::op_mul::eval(), 
                                                            ret_ti, Dc.get_type());

        value_code vcdr     = Dr.get_value_code();
        value_code vcdc     = Dc.get_value_code();
        value_code vcd      = matrix_traits::unify_value_types(vcdr, vcdc);
        value_code vc       = matrix_traits::value_code<V>::value;
        bool is_real        = matrix_traits::is_real(vcd);

        if (std::is_same<V,Object>::value && ret_ti != A.get_type())
        {
            Matrix Ac       = convert_object(A, ti::convert_ti_object<V>(ret_ti));
            Matrix Ac2      = convert(Ac, Mat::matrix_code);
            Matrix Dr2      = convert(Dr, Mat_D::matrix_code);
            Matrix Dc2      = convert(Dc, Mat_D::matrix_code);

            const Mat& Am   = Ac2.get_impl<Mat>();
            const Mat_D& Drm = Dr2.get_impl<Mat_D>();
            const Mat_D& Dcm = Dc2.get_impl<Mat_D>();

            return mrd::scaling_helper<Mat>::eval_rowscols(ret, Am, Drm, Dcm);
        };

        if (vc == vcd || is_real == false)
        {
            Matrix Ac       = convert(A, Mat::matrix_code);
            Matrix Dr2      = convert(Dr, Mat_D::matrix_code);
            Matrix Dc2      = convert(Dc, Mat_D::matrix_code);

            const Mat& Am       = Ac.get_impl<Mat>();
            const Mat_D& Drm    = Dr2.get_impl<Mat_D>();
            const Mat_D& Dcm    = Dc2.get_impl<Mat_D>();

            return mrd::scaling_helper<Mat>::eval_rowscols(ret, Am, Drm, Dcm);
        }
        else
        {            
            Matrix Ac       = convert(A, Mat::matrix_code);
            Matrix Dr2      = convert(Dr, Mat_DR::matrix_code);
            Matrix Dc2      = convert(Dc, Mat_DR::matrix_code);

            const Mat& Am       = Ac.get_impl<Mat>();
            const Mat_DR& Drm   = Dr2.get_impl<Mat_DR>();
            const Mat_DR& Dcm   = Dc2.get_impl<Mat_DR>();

            return mrd::scaling_real_helper<Mat>::eval_rowscols(ret, Am, Drm, Dcm);
        };        
    };

    template<class T>
    static void eval(const Matrix& , const T&, Matrix& ret, const Matrix& A, const Matrix& Dr, const Matrix& Dc)
    {
        return eval_impl<T>(ret, A, Dr, Dc);
    };

    template<class T>
    static void eval_scalar(const Matrix&, const T&, Matrix& ret, const Matrix& A, const Matrix& Dr, const Matrix& Dc)
    {
        using Mat   = raw::Matrix<T,struct_dense>;
        return eval_impl<Mat>(ret, A, Dr, Dc);
    };
};

static bool is_vector(Integer M, Integer N)
{
    if (M == 0 && N == 0)
        return true;

    if (M == 1 || N == 1)
        return true;

    return false;
};

static Matrix make_dum(value_code vc, struct_code sc)
{
    switch (sc)
    {
        case struct_code::struct_dense:
        case struct_code::struct_scalar:
        default:
            return zeros(0,0,vc);
        case struct_code::struct_sparse:
            return spzeros(0,0,0,vc);
        case struct_code::struct_banded:
            return bzeros(0,0,0,0,vc);
    };

};
static void scale_rows_impl(Matrix& ret, const Matrix& A, const Matrix& Dr)
{
    if (is_vector(Dr.rows(), Dr.cols()) == false)
        throw error::vector_required(Dr.rows(), Dr.cols());

    if (A.rows() != Dr.length())
        throw error::invalid_size2(Dr.rows(), Dr.cols(), A.rows(), 1);

    value_code vc   = matrix_traits::unify_value_types(A.get_value_code(), Dr.get_value_code());
    Matrix dum      = make_dum(vc, A.get_struct_code());

    details::eval_scale_rows::make<const Matrix&>(dum,ret,A, Dr);
};

static void scale_cols_impl(Matrix& ret, const Matrix& A, const Matrix& Dc)
{
    if (is_vector(Dc.rows(), Dc.cols()) == false)
        throw error::vector_required(Dc.rows(), Dc.cols());

    if (A.cols() != Dc.length())
        throw error::invalid_size2(Dc.rows(), Dc.cols(), A.cols(), 1);

    value_code vc   = matrix_traits::unify_value_types(A.get_value_code(), Dc.get_value_code());
    Matrix dum      = make_dum(vc, A.get_struct_code());
    details::eval_scale_cols::make<const Matrix&>(dum, ret, A, Dc);
};

static void scale_rowscols_impl(Matrix& ret, const Matrix& A, const Matrix& Dr, const Matrix& Dc)
{
    if (is_vector(Dr.rows(), Dr.cols()) == false)
        throw error::vector_required(Dr.rows(), Dr.cols());

    if (is_vector(Dc.rows(), Dc.cols()) == false)
        throw error::vector_required(Dc.rows(), Dc.cols());

    if (A.rows() != Dr.length())
        throw error::invalid_size2(Dr.rows(), Dr.cols(), A.rows(), 1);

    if (A.cols() != Dc.length())
        throw error::invalid_size2(Dc.rows(), Dc.cols(), A.cols(), 1);

    value_code vc   = matrix_traits::unify_value_types(A.get_value_code(), Dr.get_value_code());
    vc              = matrix_traits::unify_value_types(vc, Dc.get_value_code());

    Matrix dum      = make_dum(vc, A.get_struct_code());
    details::eval_scale_rowscols::make<const Matrix&>(dum, ret, A, Dr, Dc);
};

}};

namespace matcl
{

Matrix matcl::symprod(const Matrix& A0, bool trans)
{
    Matrix A(A0);

    Matrix ret;
    details::eval_SYMPROD::make<const Matrix&>(A,ret,trans);
    return ret;
}

Matrix matcl::herprod(const Matrix& A0, bool trans)
{
    Matrix A(A0);

    Matrix ret;
    details::eval_HERPROD::make<const Matrix&>(A,ret,trans);
    return ret;
}

Matrix matcl::symprod(Matrix&& A0, bool trans)
{
    Matrix A(std::move(A0));

    Matrix ret;
    details::eval_SYMPROD::make<const Matrix&>(A,ret,trans);
    return ret;
}

Matrix matcl::herprod(Matrix&& A0, bool trans)
{
    Matrix A(std::move(A0));

    Matrix ret;
    details::eval_HERPROD::make<const Matrix&>(A,ret,trans);
    return ret;
}

Matrix matcl::symsum(const Matrix& A0)
{
    Matrix A(A0);

    Matrix ret;
    details::eval_SYMSUM::make<const Matrix&>(A,ret);
    return ret;
}

Matrix matcl::symsum(Matrix&& A0)
{
    Matrix A(std::move(A0));

    Matrix ret;
    details::eval_SYMSUM::make<const Matrix&>(A,ret);
    return ret;
}

Matrix matcl::hersum(const Matrix& A0)
{
    Matrix A(A0);

    Matrix ret;
    details::eval_HERSUM::make<const Matrix&>(A,ret);
    return ret;
}

Matrix matcl::hersum(Matrix&& A0)
{
    Matrix A(std::move(A0));

    Matrix ret;
    details::eval_HERSUM::make<const Matrix&>(A,ret);
    return ret;
}

Matrix matcl::scale_rows(const Matrix& A0, const Matrix& D)
{
    Matrix A(A0);

    Matrix ret;
    details::scale_rows_impl(ret, A, D);    
    return ret;
};

Matrix matcl::scale_rows(Matrix&& A0, const Matrix& D)
{
    Matrix A(std::move(A0));

    Matrix ret;
    details::scale_rows_impl(ret, A, D);
    return ret;
};

Matrix matcl::scale_cols(const Matrix& A0, const Matrix& D)
{
    Matrix A(A0);

    Matrix ret;
    details::scale_cols_impl(ret, A, D);
    return ret;
};

Matrix matcl::scale_cols(Matrix&& A0, const Matrix& D)
{
    Matrix A(std::move(A0));

    Matrix ret;
    details::scale_cols_impl(ret, A, D);
    return ret;
};

Matrix matcl::scale_rowscols(const Matrix& A0, const Matrix& Dr, const Matrix& Dc)
{
    Matrix A(A0);

    Matrix ret;
    details::scale_rowscols_impl(ret, A, Dr, Dc);
    return ret;
};

Matrix matcl::scale_rowscols(Matrix&& A0, const Matrix& Dr, const Matrix& Dc)
{
    Matrix A(std::move(A0));

    Matrix ret;
    details::scale_rowscols_impl(ret, A, Dr, Dc);
    return ret;
};

};

#pragma warning(pop)