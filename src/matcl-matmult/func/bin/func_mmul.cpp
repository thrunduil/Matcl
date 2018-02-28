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
#include "matcl-matrep/lib_functions/func_matrix.h"
#include "matcl-matrep/lib_functions/matrix_gen.h"
#include "matcl-matrep/lib_functions/func_unary.h"
#include "matcl-matrep/details/extract_type2_switch.h"
#include "matcl-scalar/details/matfunc_helpers.h"
#include "matcl-core/matrix/matrix_traits.h"
#include "matcl-matmult/func/raw/mmul/mmul.h"
#include "matcl-matrep/matrix/unique_matrix.h"

#pragma warning( push )
#pragma warning(disable:4127)	// conditional expression is constant

namespace matcl { namespace details
{

namespace md = matcl::details;
namespace mr = matcl::raw;
namespace mrd = matcl::raw::details;

struct mult_visitor : public details::extract_type2_switch<void,mult_visitor,mr::val_type_corrector_int,
                                                            details::ver_nonstatic>
{
    template<class T1, class T2>
    void eval_mat_mat(const T1& A, const T2& B, matcl::Matrix& ret, trans_type t_A, trans_type t_B)
    {
        mrd::mult_helper<T1,T2>::eval(ret, A, B, t_A, t_B);
        return;
    };
   
    template<class T1, class T2>
    void eval_scal_scal(const T1& A, const T2& B, matcl::Matrix& ret, trans_type t_A, trans_type t_B)
    {
        ret = mrd::mmul_helper<T1,T2>::eval( (t_A == trans_type::conj_trans)? conj(A) : A,
                                             (t_B == trans_type::conj_trans)? conj(B) : B);
    };
    
    template<class T1, class T2>
    void eval_mat_scal(const T1& mat, const T2& scal, matcl::Matrix& ret, trans_type t_A, trans_type t_B)
    {
        mrd::mult_helper_mat_scal<T1,T2>::eval(ret, mat, scal, t_A, t_B);
        return;
    };
    
    template<class T1, class T2>
    void eval_scal_mat(const T1& scal, const T2& mat, matcl::Matrix& ret, trans_type t_A, trans_type t_B)
    {
        mrd::mult_helper_scal_mat<T1,T2>::eval(ret, scal, mat, t_A, t_B);
        return;
    };
};

static Matrix convert_scalar(const Matrix& alpha, const Matrix& C)
{
    value_code vc_Cr    = matrix_traits::real_value_type(C.get_value_code());
    value_code vc_a     = alpha.get_value_code();
    value_code vc       = matrix_traits::unify_value_types(vc_a, vc_Cr);

    if (vc == vc_a)
        return alpha;

    Matrix ac           = convert_value(alpha, vc);
    return ac;
};

struct mmul_abs_visitor : public details::extract_type2_switch<void,mmul_abs_visitor,mr::val_type_corrector_int,
                                                            details::ver_nonstatic>
{
    template<class T1, class T2>
    void eval_mat_mat(const T1& A, const T2& X, trans_type t_A, const Matrix& C, matcl::Matrix& ret)
    {
        mrd::mult_helper<T1,T2>::eval_abs(ret, A, X, t_A, C);
        return;
    };

    template<class T1, class T2>
    void eval_scal_scal(const T1& A, const T2& X, trans_type t_A, const Matrix& C, matcl::Matrix& ret)
    {
        (void)t_A;
        mrd::mult_helper_scal_scal::eval_abs(ret, A, X, C);
    };

    template<class T1, class T2>
    void eval_mat_scal(const T1& mat, const T2& scal, trans_type t_A, const Matrix& C, matcl::Matrix& ret)
    {
        mrd::mult_helper_mat_scal<T1,T2>::eval_abs(ret, mat, scal, t_A, C);
    };

    template<class T1, class T2>
    void eval_scal_mat(const T1& scal, const T2& mat, trans_type t_A, const Matrix& C, matcl::Matrix& ret)
    {
        (void)t_A;
        mrd::mult_helper_scal_mat<T1,T2>::eval_abs(ret, scal, mat, C);
    };
};

struct gemm_visitor : public details::extract_type2_switch<void,gemm_visitor,mr::val_type_corrector_int,
                                                            details::ver_nonstatic>
{
    template<class T1, class T2>
    void eval_mat_mat(const T1& A, const T2& B, trans_type t_A, trans_type t_B,
                      const Matrix& alpha, const Matrix& beta, Matrix& C,
                      Integer fr, Integer rows, Integer fc, Integer cols)
    {
        mrd::gemm_helper<T1,T2>::eval(A, B, t_A, t_B, alpha, beta, C, fr, rows, fc, cols);
        return;
    };

    template<class T1, class T2>
    void eval_scal_scal(const T1& A, const T2& B, trans_type t_A, trans_type t_B,
                        const Matrix& alpha, const Matrix& beta, Matrix& C,
                        Integer fr, Integer rows, Integer fc, Integer cols)
    {
        Matrix Ac   = (t_A == trans_type::conj_trans)? conj(A) : A;
        Matrix Bc   = (t_B == trans_type::conj_trans)? conj(B) : B;

        Matrix ac   = convert_scalar(alpha, C);
        Matrix ABa  = (ac * Ac) * Bc;

        mrd::gemm_helper_scal_scal::eval(ABa, beta, C, fr, rows, fc, cols);
    };

    template<class T1, class T2>
    void eval_mat_scal(const T1& mat, const T2& scal, trans_type t_A, trans_type t_B,
                       const Matrix& alpha, const Matrix& beta, Matrix& C,
                       Integer fr, Integer rows, Integer fc, Integer cols)
    {
        mrd::gemm_helper_mat_scal<T1,T2>::eval(mat, scal, t_A, t_B, alpha, beta, C, fr, rows, fc, cols);
    };

    template<class T1, class T2>
    void eval_scal_mat(const T1& scal, const T2& mat, trans_type t_A, trans_type t_B,
                       const Matrix& alpha, const Matrix& beta, Matrix& C,
                       Integer fr, Integer rows, Integer fc, Integer cols)
    {
        mrd::gemm_helper_scal_mat<T1,T2>::eval(scal, mat, t_A, t_B, alpha, beta, C, fr, rows, fc, cols);
    };
};

static void convert_precision(Matrix& Ac, Matrix& Bc, const Matrix& A, const Matrix& B, const Matrix& C)
{
    bool is_int_A       = A.get_value_code() == value_code::v_integer;
    bool is_int_B       = B.get_value_code() == value_code::v_integer;
    bool is_int_C       = C.get_value_code() == value_code::v_integer;

    bool is_single_A    = matrix_traits::is_single_precision(A.get_value_code());
    bool is_single_B    = matrix_traits::is_single_precision(B.get_value_code());
    bool is_single_C    = matrix_traits::is_single_precision(C.get_value_code());

    bool convert        = false;

    if (is_int_A == true && is_int_B == true && is_int_C == false)
        convert         = true;
    else if (is_single_C == true)
        convert         = false;
    else if (is_single_A == false || is_single_B == false)
        convert         = false;
    else
        convert         = true;

    if (convert == false)
    {
        Ac  = A;
        Bc  = B;
        return;
    };

    Integer nz_A        = A.structural_nnz();
    Integer nz_B        = B.structural_nnz();

    if (nz_A < nz_B)
    {
        //convert_A;
        value_code vc_A2    = matrix_traits::double_precision(A.get_value_code());
        Ac                  = convert_value(A, vc_A2);
        Bc                  = B;
        return;
    }
    else
    {
        //convert_B;
        value_code vc_B2    = matrix_traits::double_precision(B.get_value_code());
        Bc                  = convert_value(B, vc_B2);
        Ac                  = A;
        return;
    }
};

}}

matcl::Matrix matcl::mmul(const Matrix& A0, const Matrix& B0, trans_type t_A, trans_type t_B)
{
    //increase refcount
    Matrix A(A0);
    Matrix B(B0);

    Matrix ret;
    details::mult_visitor().make(A,B,ret,t_A,t_B);
    return ret;
};

matcl::Matrix matcl::mmul(Matrix&& A0, const Matrix& B0, trans_type t_A, trans_type t_B)
{
    //increase refcount
    Matrix A(std::move(A0));
    Matrix B(B0);

    Matrix ret;
    details::mult_visitor().make(A,B,ret,t_A,t_B);
    return ret;
};

matcl::Matrix matcl::mmul(Matrix&& A0, Matrix&& B0, trans_type t_A, trans_type t_B)
{
    //increase refcount
    Matrix A(std::move(A0));
    Matrix B(std::move(B0));

    Matrix ret;
    details::mult_visitor().make(A,B,ret,t_A,t_B);
    return ret;
};

matcl::Matrix matcl::mmul(const Matrix& A0, Matrix&& B0, trans_type t_A, trans_type t_B)
{
    //increase refcount
    Matrix A(A0);
    Matrix B(std::move(B0));

    Matrix ret;
    details::mult_visitor().make(A,B,ret,t_A,t_B);
    return ret;
};

void matcl::gemm(const Matrix& alpha, const Matrix& A, const Matrix& B, 
                        trans_type t_A, trans_type t_B, const Matrix& beta, Matrix& C)
{    
    Integer fr      = 0;
    Integer rows    = C.rows();
    Integer fc      = 0;
    Integer cols    = C.cols();

    details::prepare_for_assign(C, A, B);

    //make sure, that A or B has proper precision
    Matrix Ac, Bc;
    convert_precision(Ac, Bc, A, B, C);

    details::gemm_visitor().make(Ac, Bc, t_A, t_B, alpha, beta, C, fr, rows, fc, cols);
};

matcl::Matrix matcl::mmul_abs(const matcl::Matrix& A0, const matcl::Matrix& X0, matcl::trans_type t_A, 
                              const matcl::Matrix& B0)
{
    matcl::Matrix A(A0);
    matcl::Matrix X(X0);
    matcl::Matrix B(B0);

    matcl::Matrix ret;
    details::mmul_abs_visitor().make(A, X, t_A, B, ret);
    return ret;
};

matcl::Matrix matcl::mmul_abs(const matcl::Matrix& A0, const matcl::Matrix& X0, trans_type t_A, 
                              matcl::Matrix&& B0)
{
    Matrix A(A0);
    Matrix X(X0);
    Matrix B(std::move(B0));

    Matrix ret;
    details::mmul_abs_visitor().make(A, X, t_A, B, ret);
    return ret;
};

void matcl::gemm(const Matrix& alpha, const Matrix& A, const Matrix& B, 
                                 trans_type t_A, trans_type t_B, const Matrix& beta, sub_matrix_1&& sub_C)
{
    Integer fr, rows, fc, cols;
    Matrix& C   = sub_C.get_view_info(fr, rows, fc, cols);

    //make sure, that A or B has proper precision
    Matrix Ac, Bc;
    convert_precision(Ac, Bc, A, B, C);

    details::gemm_visitor().make(Ac, Bc, t_A, t_B, alpha, beta, C, fr-1, rows, fc-1, cols);
}

void matcl::gemm(const Matrix& alpha, const Matrix& A, const Matrix& B, 
                                 trans_type t_A, trans_type t_B, const Matrix& beta, sub_matrix_2&& sub_C)
{
    Integer fr, rows, fc, cols;
    Matrix& C   = sub_C.get_view_info(fr, rows, fc, cols);

    //make sure, that A or B has proper precision
    Matrix Ac, Bc;
    convert_precision(Ac, Bc, A, B, C);

    details::gemm_visitor().make(Ac, Bc, t_A, t_B, alpha, beta, C, fr-1, rows, fc-1, cols);
}

void matcl::gemm(const Matrix& alpha, const Matrix& A, const Matrix& B, 
                                 trans_type t_A, trans_type t_B, const Matrix& beta, sub_matrix&& sub_C)
{
    Integer fr, rows, fc, cols;
    Matrix& C   = sub_C.get_view_info(fr, rows, fc, cols);

    details::prepare_for_assign(C, A, B);

    //make sure, that A or B has proper precision
    Matrix Ac, Bc;
    convert_precision(Ac, Bc, A, B, C);

    details::gemm_visitor().make(Ac, Bc, t_A, t_B, alpha, beta, C, fr-1, rows, fc-1, cols);
};

void matcl::gemm(const Matrix& alpha, const Matrix& A, const Matrix& B, 
                                 trans_type t_A, trans_type t_B, const Matrix& beta, unique_matrix&& C)
{
    if (C.is_matrix())
        return gemm(alpha, A, B, t_A, t_B, beta, std::move(C).get_matrix());
    else if (C.is_submatrix())
        return gemm(alpha, A, B, t_A, t_B, beta, std::move(C).get_submatrix());
    else if (C.is_submatrix_1())
        return gemm(alpha, A, B, t_A, t_B, beta, std::move(C).get_submatrix_1());
    else if (C.is_submatrix_2())
        return gemm(alpha, A, B, t_A, t_B, beta, std::move(C).get_submatrix_2());
    else
    {
        matcl_assert(0,"invalid case");
        throw error::error_general("invalid case");
    }
};

namespace matcl { namespace details
{

struct mmul_impl
{
    template<class Mat1, class Mat2>
    static Matrix eval(Mat1&& A, Mat2&& B, trans_type_ext t_A, trans_type_ext t_B)
    {
        bool conj_A;
        bool conj_B;

        trans_type tA   = details::trans_manip::convert_trans(t_A, conj_A);
        trans_type tB   = details::trans_manip::convert_trans(t_B, conj_B);

        if (conj_A == true && matrix_traits::is_float_real(A.get_value_code()) == true)
            conj_A      = false;

        if (conj_B == true && matrix_traits::is_float_real(B.get_value_code()) == true)
            conj_B      = false;

        if (conj_A == false && conj_B == false)
            return mmul(std::forward<Mat1>(A), std::forward<Mat2>(B), tA, tB);

        if (conj_A == true)
        {
            if (conj_B == true)
            {
                //conj(A) * conj(B) = conj(A*B)
                Matrix ret  = mmul(std::forward<Mat1>(A), std::forward<Mat2>(B), tA, tB);
                return conj(std::move(ret));
            }
            else
            {
                //conj(A) * B = conj(A*conj(B))

                if (A.numel() < B.numel())
                {
                    //conj(A) * B
                    Matrix ret  = mmul(conj(std::forward<Mat1>(A)), std::forward<Mat2>(B), tA, tB);
                    return ret;
                }
                else
                {
                    //conj(A*conj(B))
                    Matrix ret  = mmul(std::forward<Mat1>(A), conj(std::forward<Mat2>(B)), tA, tB);
                    return conj(std::move(ret));
                }
            };
        }
        else
        {
            //conj_B == true

            //A * conj(B) = conj(conj(A)*B)
            if (A.numel() < B.numel())
            {
                //conj(conj(A)*B)
                Matrix ret  = mmul(conj(std::forward<Mat1>(A)), std::forward<Mat2>(B), tA, tB);
                return conj(std::move(ret));
            }
            else
            {
                //A * conj(B)
                Matrix ret  = mmul(std::forward<Mat1>(A), conj(std::forward<Mat2>(B)), tA, tB);
                return ret;
            }
        };   
    };
};

}};

matcl::Matrix matcl::mmul(const Matrix& A, const Matrix& B, trans_type_ext t_A, trans_type_ext t_B)
{
    return details::mmul_impl::eval(A,B,t_A, t_B);
};

matcl::Matrix matcl::mmul(Matrix&& A, const Matrix& B, trans_type_ext t_A, trans_type_ext t_B)
{
    return details::mmul_impl::eval(std::move(A),B,t_A, t_B);
};

matcl::Matrix matcl::mmul(const Matrix& A, Matrix&& B, trans_type_ext t_A, trans_type_ext t_B)
{
    return details::mmul_impl::eval(A,std::move(B),t_A, t_B);
};

matcl::Matrix matcl::mmul(Matrix&& A, Matrix&& B, trans_type_ext t_A, trans_type_ext t_B)
{
    return details::mmul_impl::eval(std::move(A),std::move(B),t_A, t_B);
};

#pragma warning( pop )
