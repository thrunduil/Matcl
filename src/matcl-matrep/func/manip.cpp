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


#include "matcl-matrep/lib_functions/manip.h"
#include "matcl-matrep/lib_functions/matrix_utils.h"
#include "matcl-matrep/lib_functions/func_unary.h"
#include "matcl-internals/container/mat_d.h"
#include "matcl-internals/container/mat_s.h"
#include "matcl-internals/container/mat_b.h"
#include "matcl-matrep/details/extract_type_switch.h"
#include "matcl-matrep/func/raw/raw_manip.h"
#include "matcl-core/matrix/matrix_traits.h"
#include "matcl-scalar/details/scalfunc_helpers.h"
#include "matcl-matrep/func/raw/find.h"
#include "matcl-matrep/func/raw/sort2.h"
#include "matcl-matrep/matrix/matrix_concat.h"
#include "matcl-internals/base/instantiate.h"
#include "matcl-internals/func/converter.h"
#include "matcl-matrep/details/details_manip.h"
#include "matcl-matrep/details/matrix.inl"
#include "matcl-matrep/lib_functions/matrix_gen.h"

namespace matcl
{

namespace md  = matcl::details;
namespace mr  = matcl::raw;
namespace mrd = matcl::raw::details;

namespace details
{

struct eval_clone : extract_type_switch<void, eval_clone, true>
{
    template<class T>
    static void eval(const Matrix& , const T& mat, Matrix& ret)
    {
        T tmp = mat.clone();
        ret = Matrix(tmp,false);
    };

    template<class T>
    static void eval_scalar(const Matrix& handle, const T& , Matrix& ret)
    {
        ret = handle;
    };

    static void eval_scalar(const Matrix&, const Object& m, Matrix& ret)
    {
        ret = m.clone();
    };
};

template<class ret_type, class in_type, bool is_proper>
struct function_op_convert_eval
{
    static void eval(const in_type& A, Matrix& ret)
    {
        const ret_type& tmp = raw::converter<ret_type, in_type>::eval(A);
        ret = Matrix(tmp, false);
    };
};

template<class M1,enum class matcl::mat_code mt>
struct function_op_convert_impl
{
    using new_mat_type  = typename matrix_traits::mat_type_info_code<mt>::matrix_type;

    static void eval(const M1& A, Matrix& ret)
    {
        static const bool is_pr_val = true;
        return function_op_convert_eval<new_mat_type,M1,is_pr_val>::eval(A, ret);
    };
};

template<class M1,class T>
struct function_op_convert_scal
{
    static void eval(const M1& A, Matrix& ret)
    {
        using old_val_type      = typename M1::value_type;
        using old_struct_type   = typename M1::struct_type;

        static const bool is_pr_val = true;

        if (A.rows() != 1 || A.cols() != 1 || !is_pr_val)
        {
            matcl::mat_code ret_code    = matrix_traits::mat_type_info_type_2<old_val_type,old_struct_type>
                                        ::matrix_code;
            matcl::mat_code in_code     = matrix_traits::get_matrix_type(value_to_code<T>::value,
                                                                         struct_code::struct_scalar);
            throw error::unable_to_convert(ret_code,in_code);
        }
        else
        {
            ret = T(matcl::raw::converter<T,old_val_type>::eval(A(1,1)));
        };
    };
};

struct eval_convert : extract_type_switch<void, eval_convert, true>
{
    template<class M1>
    static void eval(const Matrix& , const M1& A, Matrix& ret, matcl::mat_code new_type)
    {
        switch (new_type)
        {		
            case mat_code::integer_dense:
            {
                function_op_convert_impl<M1,mat_code::integer_dense>::eval(A, ret);
                return;
            }
            case mat_code::real_dense:
            {
                function_op_convert_impl<M1,mat_code::real_dense>::eval(A, ret);
                return;
            }
            case mat_code::float_dense:
            {
                function_op_convert_impl<M1,mat_code::float_dense>::eval(A, ret);
                return;
            }
            case mat_code::complex_dense:
            {
                function_op_convert_impl<M1,mat_code::complex_dense>::eval(A, ret);
                return;
            }
            case mat_code::float_complex_dense:
            {
                function_op_convert_impl<M1,mat_code::float_complex_dense>::eval(A, ret);
                return;
            }
            case mat_code::object_dense:
            {
                function_op_convert_impl<M1,mat_code::object_dense>::eval(A, ret);
                return;
            }
            case mat_code::integer_sparse:
            {
                function_op_convert_impl<M1,mat_code::integer_sparse>::eval(A, ret);
                return;
            }
            case mat_code::real_sparse:
            {
                function_op_convert_impl<M1,mat_code::real_sparse>::eval(A, ret);
                return;
            }
            case mat_code::float_sparse:
            {
                function_op_convert_impl<M1,mat_code::float_sparse>::eval(A, ret);
                return;
            }
            case mat_code::complex_sparse:
            {
                function_op_convert_impl<M1,mat_code::complex_sparse>::eval(A, ret);
                return;
            }
            case mat_code::float_complex_sparse:
            {
                function_op_convert_impl<M1,mat_code::float_complex_sparse>::eval(A, ret);
                return;
            }
            case mat_code::object_sparse:
            {
                function_op_convert_impl<M1,mat_code::object_sparse>::eval(A,ret);
                return;
            }
            case mat_code::integer_band:
            {
                function_op_convert_impl<M1,mat_code::integer_band>::eval(A,ret);
                return;
            }
            case mat_code::real_band:
            {
                function_op_convert_impl<M1,mat_code::real_band>::eval(A,ret);
                return;
            }
            case mat_code::float_band:
            {
                function_op_convert_impl<M1,mat_code::float_band>::eval(A,ret);
                return;
            }
            case mat_code::complex_band:
            {
                function_op_convert_impl<M1,mat_code::complex_band>::eval(A,ret);
                return;
            }
            case mat_code::float_complex_band:
            {
                function_op_convert_impl<M1,mat_code::float_complex_band>::eval(A,ret);
                return;
            }
            case mat_code::object_band:
            {
                function_op_convert_impl<M1,mat_code::object_band>::eval(A,ret);
                return;
            }
            case mat_code::integer_scalar:
            {
                function_op_convert_scal<M1,Integer>::eval(A,ret);
                return;
            }
            case mat_code::real_scalar:
            {
                function_op_convert_scal<M1,Real>::eval(A,ret);
                return;
            }
            case mat_code::float_scalar:
            {
                function_op_convert_scal<M1,Float>::eval(A,ret);
                return;
            }
            case mat_code::complex_scalar:
            {
                function_op_convert_scal<M1,Complex>::eval(A,ret);
                return;
            }
            case mat_code::float_complex_scalar:
            {
                function_op_convert_scal<M1,Float_complex>::eval(A,ret);
                return;
            }
            case mat_code::object_scalar:
            {
                function_op_convert_scal<M1,Object>::eval(A,ret);
                return;
            }
            default:
            {
                throw error::unable_to_convert_invalid_code();
            }
        };
    };

    template<class M1>
    static void eval_scalar(const Matrix& handle, const M1& A, Matrix& ret, matcl::mat_code new_type)
    {
        using full_matrix = raw::Matrix<M1,struct_dense>;

        switch (new_type)
        {			
            case mat_code::integer_scalar:
            {
                ret = matcl::raw::converter<Integer,M1>::eval(A);
                return;
            }
            case mat_code::real_scalar:
            {
                ret = matcl::raw::converter<Real,M1>::eval(A);
                return;
            }
            case mat_code::float_scalar:
            {
                ret = matcl::raw::converter<Float,M1>::eval(A);
                return;
            }
            case mat_code::complex_scalar:
            {
                ret = matcl::raw::converter<Complex,M1>::eval(A);
                return;
            }
            case mat_code::float_complex_scalar:
            {
                ret = matcl::raw::converter<Float_complex,M1>::eval(A);
                return;
            }
            default:
            {
                full_matrix AF(ti::get_ti(A),A,1,1);
                eval<full_matrix>(handle, AF, ret, new_type);
                return;
            }
        };
    };
};

struct eval_convert_object : extract_type_switch<void, eval_convert_object, true>
{
    template<class M1>
    static void eval(const Matrix& , const M1& A, Matrix& ret, ti::ti_object ti)
    {
        using val_type      = typename M1::value_type;
        using struct_type   = typename M1::struct_type;
        using new_type      = matcl::raw::Matrix<Object,struct_type> ;

        Matrix th;
        ret = Matrix(matcl::raw::converter<new_type,M1>::eval(A,ti, th), false);
    };

    template<class T>
    static void eval_scalar(const Matrix&, const T& A, Matrix& ret, ti::ti_object ti)
    {
        ret = convert_to_object(ti, A);
    };
};

template<class M1, class struct_type>
struct function_op_full
{
    static void eval(const M1& A, Matrix& ret)
    {
        using matrix_type = typename raw::Matrix<typename M1::value_type,struct_dense>;
        ret = Matrix(raw::converter<matrix_type,M1>::eval(A),false);
    };
};

template<class M1>
struct function_op_full<M1,struct_dense>
{
    static void eval(const M1& A, Matrix& ret)
    {
        ret = Matrix(A,false);
    };
};

struct eval_full : extract_type_switch<void, eval_full, true>
{
    template<class T>
    static void eval(const Matrix& , const T& mat, Matrix& ret)
    {
        using struct_type   = typename T::struct_type;
        function_op_full<T,struct_type>::eval(mat,ret);
    };

    template<class T>
    static void eval_scalar(const Matrix&, const T& mat, Matrix& ret)
    {
        using full_matrix   = raw::Matrix<T,struct_dense> ;
        full_matrix tmp(ti::get_ti(mat),mat,1,1);

        ret = Matrix(tmp,false);
    };
};

struct eval_vec : extract_type_switch<void, eval_vec, true>
{
    template<class T>
    static void eval(const Matrix& , const T& mat, Matrix& ret)
    {
        if (mat.cols() == 1)
        {
            ret = Matrix(mat,false);
            return;
        };

        return mrd::manip_reshape_helper<T>::eval_vec(ret, mat);
    };

    template<class T>
    static void eval_scalar(const Matrix& handle, const T& mat, Matrix& ret)
    {
        (void)mat;
        ret = handle;
    };
};

template<class M1, bool is_sp>
struct function_op_sparse
{};

template<class M1>
struct function_op_sparse<M1,true>
{
    static void eval(const M1& A, Matrix& ret)
    {
        ret = Matrix(A,false);
    };
};

template<class M1>
struct function_op_sparse<M1,false>
{
    static void eval(const M1& A, Matrix& ret)
    {
        ret = Matrix(raw::sparse(A),false);
    };
};

struct eval_sparse : extract_type_switch<void, eval_sparse, true>
{
    template<class T>
    static void eval(const Matrix& , const T& mat, Matrix& ret)
    {
        static const bool is_sp   = is_sparse<T>::value;
        function_op_sparse<T,is_sp>::eval(mat, ret);
    };

    template<class T>
    static void eval_scalar(const Matrix& handle, const T& mat, Matrix& ret)
    {
        using FullMatrix = raw::Matrix<T,struct_dense>;
        eval<FullMatrix>(handle,FullMatrix(ti::get_ti(mat),mat,1,1), ret);
    };
};

template<class M1, bool is_b>
struct function_op_band
{};

template<class M1>
struct function_op_band<M1,true>
{
    static void eval(const M1& A, Matrix& ret)
    {
        ret = Matrix(A,false);
    };
};

template<class M1>
struct function_op_band<M1,false>
{
    static void eval(const M1& A, Matrix& ret)
    {
        ret = Matrix(raw::band(A),false);
    };
};

struct eval_band : extract_type_switch<void, eval_band, true>
{
    template<class T>
    static void eval(const Matrix& , const T& mat, Matrix& ret)
    {
        static const bool is_b   = is_band<T>::value;
        function_op_band<T,is_b>::eval(mat, ret);
    };

    template<class T>
    static void eval_scalar(const Matrix& handle, const T& mat, Matrix& ret)
    {
        using FullMatrix = raw::Matrix<T,struct_dense>;
        eval<FullMatrix>(handle,FullMatrix(ti::get_ti(mat),mat,1,1),ret);
    };
};

struct eval_trans : extract_type_switch<void, eval_trans, true>
{
    template<class T>
    static void eval(const Matrix& , const T& mat, Matrix& ret)
    {        
        T ret_tmp(mat.get_type());
        mrd::manip_trans_helper<T>::eval_trans(ret_tmp, mat);
        ret = matcl::Matrix(ret_tmp,false);
    };

    template<class T>
    static void eval_scalar(const Matrix& handle, const T& mat, Matrix& ret)
    {
        (void)mat;
        ret = handle;
    };
};

template<class M1>
struct function_op_ctrans
{
    static void eval(const M1& A, Matrix& ret)
    {
        M1 ret_tmp(A.get_type());
        mrd::manip_trans_helper<M1>::eval_ctrans(ret_tmp, A);
        ret = matcl::Matrix(ret_tmp,false);
    };

    static void eval_scal(const M1& A, Matrix& ret)
    {
        ret = A;
    };
};

template<>
struct function_op_ctrans<Object>
{
    static void eval_scal(const Object& A, Matrix& ret)
    {
        ret = raw::details::conj_helper<Object>::eval(A);
    };
};

template<>
struct function_op_ctrans<Complex>
{
    static void eval_scal(const Complex& A, Matrix& ret)
    {
        ret = Matrix(Complex(matcl::real(A),-matcl::imag(A)));
    };
};

template<>
struct function_op_ctrans<Float_complex>
{
    static void eval_scal(const Float_complex& A, Matrix& ret)
    {
        ret = Matrix(Float_complex(matcl::real(A),-matcl::imag(A)));
    };
};

struct eval_ctrans : extract_type_switch<void, eval_ctrans, true>
{
    template<class T>
    static void eval(const Matrix& , const T& mat, Matrix& ret)
    {
        function_op_ctrans<T>::eval(mat, ret);
    };

    template<class T>
    static void eval_scalar(const Matrix& handle, const T& mat, Matrix& ret)
    {
        (void)handle;
        function_op_ctrans<T>::eval_scal(mat, ret);
    };
};

struct eval_fliplr : extract_type_switch<void, eval_fliplr, true>
{
    template<class T>
    static void eval(const Matrix& , const T& A, Matrix& ret)
    {
        if (A.cols() < 2)
        {
            ret = Matrix(A,false);
            return;
        };

        return mrd::manip_reshape_helper<T>::eval_fliplr(ret, A);
    };

    template<class T>
    static void eval_scalar(const Matrix& handle, const T& mat, Matrix& ret)
    {
        (void)mat;
        ret = handle;
    };
};

struct eval_flipud : extract_type_switch<void, eval_flipud, true>
{
    template<class T>
    static void eval(const Matrix& , const T& A, Matrix& ret)
    {
        if (A.rows() < 2)
        {
            ret = Matrix(A,false);
            return;
        };
        return mrd::manip_reshape_helper<T>::eval_flipud(ret, A);
    };

    template<class T>
    static void eval_scalar(const Matrix& handle, const T& mat, Matrix& ret)
    {
        (void)mat;
        ret = handle;
    };
};

struct eval_reshape : extract_type_switch<void, eval_reshape, true>
{
    template<class T>
    static void eval(const Matrix& , const T& A, Matrix& ret, Integer m, Integer n)
    {
        if (m == A.rows() && n == A.cols())
        {
            ret = Matrix(A,false);
            return;
        }

        return mrd::manip_reshape_helper<T>::eval_reshape(ret, A, m, n);
    };

    template<class T>
    static void eval_scalar(const Matrix& handle, const T& mat, Matrix& ret, Integer m, Integer n)
    {
        (void)handle;
        error::check_reshape(1, 1, m, n);
        ret = mat;
    };
};

struct eval_repmat : extract_type_switch<void, eval_repmat, true>
{
    template<class T>
    static void eval(const Matrix& handle, const T& A, Matrix& ret, Integer m, Integer n)
    {
        if (m == 1 && n == 1)
        {
            ret = Matrix(A,false);
            return;
        };
        if (A.rows() == 1 && A.cols() == 1)
        {
            eval_scalar<typename T::value_type>(handle, A(1,1),ret, m,n);
            return;
        };
        
        return mrd::manip_reshape_helper<T>::eval_repmat(ret, A, m, n);
    };

    template<class T>
    static void eval_scalar(const Matrix& handle, const T& A, Matrix& ret, Integer m, Integer n)
    {
        (void)handle;
        if (m == 1 && n == 1)
        {
            ret = A;
            return;
        };

        if (mrd::is_zero(A))
        {
            auto tmp    = md::zero_matrix<T,struct_sparse>::eval(ti::get_ti(A),m,n);
            ret         = matcl::Matrix(tmp,false);
        }
        else
        {
            using FullMatrix    = mr::Matrix<T,struct_dense>;
            ret                 = matcl::Matrix(FullMatrix(ti::get_ti(A),A,m,n),false);
        };
    };
};

struct eval_diag : extract_type_switch<void, eval_diag, true>
{
    template<class T>
    static void eval(const Matrix& handle, const T& A, Matrix& ret, Integer d)
    {
        (void)handle;
        ret = Matrix(A.get_diag(d),true);
    };

    template<class T>
    static void eval_scalar(const Matrix& handle, const T& mat, Matrix& ret, Integer d)
    {
        using FullMatrix = raw::Matrix<T,struct_dense>;
        eval<FullMatrix>(handle,FullMatrix(ti::get_ti(mat),mat,1,1),ret,d);
    };
};

struct eval_tril : extract_type_switch<void, eval_tril, true>
{
    template<class T>
    static void eval(const Matrix& handle, const T& A, Matrix& ret, Integer d, bool rvalue)
    {
        (void)handle;
        return mrd::manip_tr_helper<T>::eval_tril(ret, A, d, rvalue);
    };

    template<class T>
    static void eval_scalar(const Matrix& handle, const T& mat, Matrix& ret, Integer d, bool)
    {
        (void)handle;

        if (d < 0)
        {
            ret = md::default_value<T>(ti::get_ti(mat));
            return;
        };
        ret = mat;
    };
};

struct eval_triu : extract_type_switch<void, eval_triu, true>
{
    template<class T>
    static void eval(const Matrix& handle, const T& A, Matrix& ret, Integer d, bool rvalue)
    {
        (void)handle;
        return mrd::manip_tr_helper<T>::eval_triu(ret, A, d,rvalue);
    };

    template<class T>
    static void eval_scalar(const Matrix& handle, const T& mat, Matrix& ret, Integer d, bool)
    {
        (void)handle;
        if (d > 0)
        {
            ret = md::default_value<T>(ti::get_ti(mat));
            return;
        };
        ret = mat;
        return;
    };
};

struct eval_select_band : extract_type_switch<void, eval_select_band, true>
{
    template<class T>
    static void eval(const Matrix& handle, const T& A, Matrix& ret, Integer fd, Integer ld)
    {
        (void)handle;
        return mrd::manip_tr_helper<T>::eval_select_band(ret, A, fd, ld);
    };

    template<class T>
    static void eval_scalar(const Matrix&, const T& mat, Matrix& ret, Integer fd, Integer ld)
    {
        if (fd <= 0 && ld >= 0)
            ret = mat;
        else
            ret = md::default_value<T>(ti::get_ti(mat));

        return;
    };
};

struct eval_is_sym : extract_type_switch<bool, eval_is_sym, true>
{
    template<class T>
    static bool eval(const Matrix& handle, const T& A, Real tol)
    {
        (void)handle;
        return raw::is_sym(A, tol, true);
    };

    template<class T>
    static bool eval_scalar(const Matrix& handle, const T& A, Real tol)
    {
        (void)handle;
        (void)A;
        (void)tol;
        return true;
    };
};

template<class T>
struct cast_scalar
{
    template<class TV>
    static T eval(const TV& val)
    {
        return (T)val;
    }
};

template<>
struct cast_scalar<Integer>
{
    template<class TV>
    static Real eval(const TV& val)
    {
        return 0.0;
    }
};

struct eval_is_her : extract_type_switch<bool, eval_is_her, true>
{
    template<class T>
    static bool eval(const Matrix& handle, const T& A, Real tol)
    {
        (void)handle;

        return raw::is_her(A, tol, true);
    };

    template<class T>
    static bool eval_scalar(const Matrix& handle, const T& A, Real tol0)
    {
        using TR    = typename md::real_type_int_real<T>::type;

        (void)handle;
        if (tol0 > 0.0)
        {
            TR tol = matcl::eps(A) * cast_scalar<TR>::eval(tol0);
            return (bool)(mrd::abs_helper<TR>::eval(mrd::imag_helper<T>::eval(A)) > tol);
        }
        else
        {
            return mrd::is_zero(mrd::imag_helper<T>::eval(A)) == true;
        };
    };
};

struct eval_get_ld : extract_type_switch<Integer, eval_get_ld, true>
{
    template<class T>
    static Integer eval(const Matrix& handle, const T& A, Integer min)
    {
        (void)handle;
        return raw::get_ld(A, min);
    };

    template<class T>
    static Integer eval_scalar(const Matrix& handle, const T& A, Integer min)
    {
        (void)handle;
        (void)A;
        (void)min;
        return 0;
    };
};

struct eval_get_ud : extract_type_switch<Integer, eval_get_ud, true>
{
    template<class T>
    static Integer eval(const Matrix& handle, const T& A, Integer min)
    {
        (void)handle;
        return raw::get_ud(A,min);
    };

    template<class T>
    static Integer eval_scalar(const Matrix& handle, const T& A, Integer min)
    {
        (void)handle;
        (void)A;
        (void)min;
        return 0;
    };
};

struct eval_nnz_total : extract_type_switch<Integer, eval_nnz_total, true>
{
    template<class T>
    static Integer eval(const Matrix& handle, const T& A)
    {
        (void)handle;
        return raw::nnz_total(A);
    };

    template<class T>
    static Integer eval_scalar(const Matrix& handle, const T& A)
    {
        (void)handle;
        return mrd::is_zero<T>(A) ? 0 : 1;
    };
};

template<class struct_type, class T>
struct drop_sparse_impl
{
    static void eval(Matrix& ret, const T& mat, Real tol)
    {
        (void)tol;
        ret = matcl::Matrix(mat, false);
    };
};

template<class T>
struct drop_sparse_impl<struct_sparse,T>
{
    static void eval(Matrix& ret, const T& mat, Real tol)
    {
        ret = matcl::Matrix(mat.drop(tol).get(), false);
    };
};

struct drop_sparse_visitor : public extract_type_switch<void, drop_sparse_visitor, true>
{
    template<class T>
    static void eval(const Matrix&, const T& mat, Matrix& ret, Real tol)        
    {
        using struct_type = typename T::struct_type;

        return drop_sparse_impl<struct_type, T>::eval(ret, mat,tol);
    };

    template<class T>
    static void eval_scalar(const Matrix& handle, const T& val, Matrix& ret, Real tol)
    {
        (void)tol;
        (void)handle;
        ret = Matrix(val,false);
    };
};

struct eval_find : public extract_type_switch<void, eval_find, true>
{
    template<class T>
    static void eval(const Matrix&, const T& A, Matrix& ret)        
    {
        return mrd::find_helper<T>::eval_find(ret, A);
    };

    template<class T>
    static void eval_scalar(const Matrix& handle, const T& mat, Matrix& ret)
    {
        using full_matrix = raw::Matrix<T,struct_dense>;
        return eval<full_matrix>(handle, full_matrix(ti::get_ti(mat),mat,1,1), ret);
    };
};

struct eval_find_t : public extract_type_switch<void, eval_find_t, true>
{
    template<class T>
    static void eval(const Matrix&, const T& A, Matrix& ret, const test_function& t)        
    {
        return mrd::find_helper<T>::eval_find(ret, A, t);
    };

    template<class T>
    static void eval_scalar(const Matrix& handle, const T& mat, Matrix& ret, const test_function& t)
    {
        using full_matrix = raw::Matrix<T,struct_dense>;
        return eval<full_matrix>(handle, full_matrix(ti::get_ti(mat),mat,1,1), ret, t);
    };
};

struct eval_find2 : public extract_type_switch<void, eval_find2, true>
{
    template<class T>
    static void eval(const Matrix&, const T& mat, Matrix& i, Matrix& j)        
    {
        return mrd::find_helper<T>::eval_find_2(i, j, mat);
    };

    template<class T>
    static void eval_scalar(const Matrix& handle, const T& mat, Matrix& ret1, Matrix& ret2)
    {
        using full_matrix = raw::Matrix<T,struct_dense>;
        return eval<full_matrix>(handle, full_matrix(ti::get_ti(mat),mat,1,1), ret1, ret2);
    };
};

struct eval_find2_t : public extract_type_switch<void, eval_find2_t, true>
{
    template<class T>
    static void eval(const Matrix&, const T& mat, Matrix& i, Matrix& j, const test_function& t)        
    {
        return mrd::find_helper<T>::eval_find_2(i, j, mat, t);
    };

    template<class T>
    static void eval_scalar(const Matrix& handle, const T& mat, Matrix& ret1, Matrix& ret2, 
                            const test_function& t)
    {
        using full_matrix = raw::Matrix<T,struct_dense>;
        return eval<full_matrix>(handle, full_matrix(ti::get_ti(mat),mat,1,1), ret1, ret2, t);
    };
};

struct eval_find3 : public extract_type_switch<void, eval_find3, true>
{
    template<class T>
    static void eval(const Matrix&, const T& mat, Matrix& i, Matrix& j, Matrix& x)        
    {
        return mrd::find_helper<T>::eval_find_3(i, j, x, mat);
    };

    template<class T>
    static void eval_scalar(const Matrix& handle, const T& mat, Matrix& ret1, Matrix& ret2, 
                            Matrix& ret3)
    {
        using full_matrix = raw::Matrix<T,struct_dense>;
        return eval<full_matrix>(handle, full_matrix(ti::get_ti(mat),mat,1,1), ret1, ret2, ret3);
    };
};

struct eval_find3_t : public extract_type_switch<void, eval_find3_t, true>
{
    template<class T>
    static void eval(const Matrix&, const T& mat, Matrix& i, Matrix& j, Matrix& x, 
                     const test_function& t)        
    {
        return mrd::find_helper<T>::eval_find_3(i, j, x, mat, t);
    };

    template<class T>
    static void eval_scalar(const Matrix& handle, const T& mat, Matrix& ret1, Matrix& ret2, 
                            Matrix& ret3, const test_function& t)
    {
        using full_matrix = raw::Matrix<T,struct_dense>;
        return eval<full_matrix>(handle, full_matrix(ti::get_ti(mat),mat,1,1), ret1, ret2, ret3, t);
    };
};

struct eval_sort : public extract_type_switch<void, eval_sort, true>
{
    template<class T>
    static void eval(const Matrix&, const T& mat, Matrix& ret, int dim, bool asceding)        
    {
        return mrd::sort_helper<T>::eval_sort(ret, mat, dim, asceding);
    };

    template<class T>
    static void eval_scalar(const Matrix& handle, const T& mat, Matrix& ret, int dim, bool asceding)
    {
        using full_matrix = raw::Matrix<T,struct_dense>;
        return eval<full_matrix>(handle,full_matrix(ti::get_ti(mat),mat,1,1), ret, dim, asceding);
    };
};

struct eval_sort2 : public extract_type_switch<void, eval_sort2, true>
{
    template<class T>
    static void eval(const Matrix&, const T& mat, Matrix& i, Matrix& x, int dim, bool asceding)        
    {
        return mrd::sort_helper<T>::eval_sort_2(i, x, mat, dim, asceding);
    };

    template<class T>
    static void eval_scalar(const Matrix& handle, const T& mat, Matrix& ret1, Matrix& ret2, 
                            int dim, bool asceding)
    {
        using full_matrix = raw::Matrix<T,struct_dense>;
        return eval<full_matrix>(handle,full_matrix(ti::get_ti(mat),mat,1,1), ret1, ret2, dim, asceding);
    };
};

struct eval_sortrows : public extract_type_switch<void, eval_sortrows, true>
{
    template<class T>
    static void eval(const Matrix&, const T& mat, Matrix& ret)
    {
        return mrd::sort_helper<T>::eval_sortrows(ret, mat);
    };

    template<class T>
    static void eval_scalar(const Matrix& handle, const T& mat, Matrix& ret)
    {
        using full_matrix = raw::Matrix<T,struct_dense>;
        return eval<full_matrix>(handle,full_matrix(ti::get_ti(mat),mat,1,1), ret);
    };
};

struct eval_sortrows2 : public extract_type_switch<void, eval_sortrows2, true>
{
    template<class T>
    static void eval(const Matrix&, const T& mat, Matrix& i, Matrix& x)
    {
        return mrd::sort_helper<T>::eval_sortrows_2(i, x, mat);
    };

    template<class T>
    static void eval_scalar(const Matrix& handle, const T& mat, Matrix& ret1, Matrix& ret2)
    {
        using full_matrix = raw::Matrix<T,struct_dense>;
        return eval<full_matrix>(handle,full_matrix(ti::get_ti(mat),mat,1,1), ret1, ret2);
    };
};

struct eval_sortcols : public extract_type_switch<void, eval_sortcols, true>
{
    template<class T>
    static void eval(const Matrix&, const T& mat, Matrix& ret)
    {
        return mrd::sort_helper<T>::eval_sortcols(ret, mat);
    };

    template<class T>
    static void eval_scalar(const Matrix& handle, const T& mat, Matrix& ret)
    {
        using full_matrix = raw::Matrix<T,struct_dense>;
        return eval<full_matrix>(handle,full_matrix(ti::get_ti(mat),mat,1,1), ret);
    };
};

struct eval_sortcols2 : public extract_type_switch<void, eval_sortcols2, true>
{
    template<class T>
    static void eval(const Matrix&, const T& mat, Matrix& i, Matrix& x)
    {
        return mrd::sort_helper<T>::eval_sortcols_2(i,x,mat);
    };

    template<class T>
    static void eval_scalar(const Matrix& handle, const T& mat, Matrix& ret1, Matrix& ret2)
    {
        using full_matrix = raw::Matrix<T,struct_dense>;
        return eval<full_matrix>(handle,full_matrix(ti::get_ti(mat),mat,1,1), ret1, ret2);
    };
};
struct eval_sortrows_dim : public extract_type_switch<void, eval_sortrows_dim, true>
{
    template<class T>
    static void eval(const Matrix&, const T& mat, Matrix& ret, const raw::integer_dense& dims)
    {
        return mrd::sort_helper<T>::eval_sortrows(ret, mat, dims);
    };

    template<class T>
    static void eval_scalar(const Matrix& handle, const T& mat, Matrix& ret, 
                            const raw::integer_dense& dims)
    {
        using full_matrix = raw::Matrix<T,struct_dense>;
        return eval<full_matrix>(handle,full_matrix(ti::get_ti(mat),mat,1,1), ret, dims);
    };
};

struct eval_sortrows2_dim : public extract_type_switch<void, eval_sortrows2_dim, true>
{
    template<class T>
    static void eval(const Matrix&, const T& mat, Matrix& i, Matrix& x, 
                     const raw::integer_dense& dims)
    {
        return mrd::sort_helper<T>::eval_sortrows_2(i,x,mat,dims);
    };

    template<class T>
    static void eval_scalar(const Matrix& handle, const T& mat, Matrix& ret1, Matrix& ret2, 
                            const raw::integer_dense& dims)
    {
        using full_matrix = raw::Matrix<T,struct_dense>;
        return eval<full_matrix>(handle,full_matrix(ti::get_ti(mat),mat,1,1), ret1, ret2, dims);
    };
};

struct eval_sortcols_dim : public extract_type_switch<void, eval_sortcols_dim, true>
{
    template<class T>
    static void eval(const Matrix&, const T& mat, Matrix& i, const raw::integer_dense& dims)
    {
        return mrd::sort_helper<T>::eval_sortcols(i,mat,dims);
    };

    template<class T>
    static void eval_scalar(const Matrix& handle, const T& mat, Matrix& ret, 
                            const raw::integer_dense& dims)
    {
        using full_matrix = raw::Matrix<T,struct_dense>;
        return eval<full_matrix>(handle,full_matrix(ti::get_ti(mat),mat,1,1), ret, dims);
    };
};

struct eval_sortcols2_dim : public extract_type_switch<void, eval_sortcols2_dim, true>
{
    template<class T>
    static void eval(const Matrix&, const T& mat, Matrix& i, Matrix& x, 
                     const raw::integer_dense& dims)
    {
        return mrd::sort_helper<T>::eval_sortcols_2(i,x,mat,dims);
    };

    template<class T>
    static void eval_scalar(const Matrix& handle, const T& mat, Matrix& ret1, Matrix& ret2, 
                            const raw::integer_dense& dims)
    {
        using full_matrix = raw::Matrix<T,struct_dense>;
        return eval<full_matrix>(handle,full_matrix(ti::get_ti(mat),mat,1,1), ret1, ret2, dims);
    };
};

struct eval_issorted_rows : public extract_type_switch<bool, eval_issorted_rows, true>
{
    template<class T>
    static bool eval(const Matrix&, const T& mat)
    {
        return mrd::sort_helper<T>::eval_issorted_rows(mat);
    };

    template<class T>
    static bool eval_scalar(const Matrix& handle, const T& mat)
    {
        using full_matrix = raw::Matrix<T,struct_dense>;
        return eval<full_matrix>(handle,full_matrix(ti::get_ti(mat),mat,1,1));
    };
};

struct eval_issorted_cols : public extract_type_switch<bool, eval_issorted_cols, true>
{
    template<class T>
    static bool eval(const Matrix&, const T& mat)
    {
        return mrd::sort_helper<T>::eval_issorted_cols(mat);
    };

    template<class T>
    static bool eval_scalar(const Matrix& handle, const T& mat)
    {
        using full_matrix = raw::Matrix<T,struct_dense>;
        return eval<full_matrix>(handle,full_matrix(ti::get_ti(mat),mat,1,1));
    };
};

struct eval_issorted : public extract_type_switch<void, eval_issorted, true>
{
    template<class T>
    static void eval(const Matrix&, const T& mat, matcl::Matrix& ret,int dim, bool asceding)
    {
        return mrd::sort_helper<T>::eval_issorted(ret,mat,dim,asceding);
    };

    template<class T>
    static void eval_scalar(const Matrix& handle, const T& mat, matcl::Matrix& ret,int dim, 
                            bool asceding)
    {
        using full_matrix = raw::Matrix<T,struct_dense>;
        return eval<full_matrix>(handle,full_matrix(ti::get_ti(mat),mat,1,1), ret, dim, asceding);
    };
};

};

Matrix matcl::delrows(const Matrix& A, const colon& c)
{
    return A.delrows(c);
}

Matrix matcl::delrows(Matrix&& A, const colon& c)
{
    return std::move(A).delrows(c);
}

Matrix matcl::delcols(const Matrix& A, const colon& c)
{
    return A.delcols(c);
}

Matrix matcl::delcols(const Matrix&& A, const colon& c)
{
    return std::move(A).delcols(c);
}

Matrix matcl::delrowscols(const Matrix& A, const colon& c1, const colon& c2)
{
    return A.delrowscols(c1,c2);
}

Matrix matcl::delrowscols(Matrix&& A, const colon& c1, const colon& c2)
{
    return std::move(A).delrowscols(c1,c2);
}

const Matrix Matrix::clone() const
{
    Matrix ret;
    details::eval_clone::make<const Matrix&>(*this,ret);
    return ret;
};

Matrix convert(const Matrix& A, matcl::mat_code new_type)
{
    Matrix ret;
    details::eval_convert::make<const Matrix&>(A, ret, new_type);
    return ret;
};

Matrix convert_value(const Matrix& A,matcl::value_code vc)
{
    struct_code sc  = A.get_struct_code();
    mat_code mc     = matrix_traits::get_matrix_type(vc, sc);

    return convert(A, mc);
};

Matrix convert_object(const Matrix& A, ti::ti_object ti)
{
    Matrix ret;
    details::eval_convert_object::make<const Matrix&>(A, ret, ti);
    return ret;
};

Matrix full(const Matrix& A)
{
    Matrix ret;
    details::eval_full::make<const Matrix&>(A, ret);
    return ret;
};

Matrix vec(const Matrix& A)
{
    Matrix ret;
    details::eval_vec::make<const Matrix&>(A, ret);
    return ret;
};

Matrix sparse(const Matrix& A)
{
    Matrix ret;
    details::eval_sparse::make<const Matrix&>(A, ret);
    return ret;
};

Matrix band(const Matrix& A)
{
    Matrix ret;
    details::eval_band::make<const Matrix&>(A, ret);
    return ret;
};

Matrix trans(const Matrix& A) 
{
    Matrix ret;
    details::eval_trans::make<const Matrix&>(A, ret);
    return ret;
};

Matrix trans(const Matrix& A, trans_type t) 
{
    if (t == trans_type::no_trans)
        return A;
    else if (t == trans_type::trans)
        return trans(A);
    else
        return ctrans(A);
};

Matrix trans(const Matrix& A, trans_type_ext t) 
{
    if (t == trans_type_ext::no_trans)
        return A;
    else if (t == trans_type_ext::trans)
        return trans(A);
    else if (t == trans_type_ext::conj_trans)
        return ctrans(A);
    else
        return conj(A);
};

Matrix ctrans(const Matrix& A)
{
    Matrix ret;
    details::eval_ctrans::make<const Matrix&>(A, ret);
    return ret;
};

Matrix fliplr(const Matrix& A)
{
    Matrix ret;
    details::eval_fliplr::make<const Matrix&>(A, ret);
    return ret;
};

Matrix flipud(const Matrix& A)
{
    Matrix ret;
    details::eval_flipud::make<const Matrix&>(A, ret);
    return ret;
};

Matrix reshape(const Matrix& A,Integer m, Integer n)
{
    Matrix ret;
    details::eval_reshape::make<const Matrix&>(A, ret, m, n);
    return ret;

};

Matrix repmat(const Matrix &A, Integer m, Integer n)
{
    Matrix ret;
    details::eval_repmat::make<const Matrix&>(A, ret, m, n);
    return ret;
};

Matrix get_diag(const Matrix& A,Integer d)
{ 
    Matrix ret;
    details::eval_diag::make<const Matrix&>(A, ret, d);
    return ret;
};

Matrix matcl::tril(const Matrix& A,Integer d)
{
    Matrix ret;
    details::eval_tril::make<const Matrix&>(A, ret, d, false);
    return ret;
};

Matrix matcl::tril(Matrix&& A,Integer d)
{
    Matrix ret;
    details::eval_tril::make<const Matrix&>(A, ret, d, true);
    return ret;
};

Matrix matcl::triu(const Matrix& A,Integer d)
{
    Matrix ret;
    details::eval_triu::make<const Matrix&>(A, ret, d, false);
    return ret;
};

Matrix matcl::triu(Matrix&& A,Integer d)
{
    Matrix ret;
    details::eval_triu::make<const Matrix&>(A, ret, d, true);
    return ret;
};

Matrix matcl::select_band(const Matrix &m, Integer fd, Integer ld)
{
    Matrix ret;
    details::eval_select_band::make<const Matrix&>(m, ret, fd, ld);
    return ret;
};

Integer get_ld(const Matrix& A, Integer min)
{
    return details::eval_get_ld::make<const Matrix&>(A, min);
};

Integer get_ud(const Matrix& A, Integer min)
{
    return details::eval_get_ud::make<const Matrix&>(A, min);
};

// number of nonzero elements
Integer matcl::nnz(const Matrix& A)
{
    return details::eval_nnz_total::make<const Matrix&>(A);
};

Matrix drop_sparse(const Matrix& A, Real tol)
{
    Matrix ret;
    details::drop_sparse_visitor::make<const Matrix&>(A, ret,tol);
    return ret;
};

bool is_tril(const Matrix& A)
{
    return get_ud(A,0) > 0 ? false : true;
};

bool is_triu(const Matrix& A)
{
    return get_ld(A,0) > 0 ? false : true;
};

bool is_diag(const Matrix& A)
{
    return (get_ld(A,0) == 0) && (get_ud(A,0) == 0);
};

bool is_sym(const Matrix& A, Real tol)
{
    return details::eval_is_sym::make<const Matrix&>(A, tol);
};

bool is_her(const Matrix& A, Real tol)
{
    return details::eval_is_her::make<const Matrix&>(A, tol);
};

bool matcl::is_real_matrix(const Matrix& A)
{
    value_code vc = A.get_value_code();
    return vc == value_code::v_integer || vc == value_code::v_real || vc == value_code::v_float;
};

bool matcl::is_real_float_matrix(const Matrix& A)
{
    value_code vc = A.get_value_code();
    return vc == value_code::v_real || vc == value_code::v_float;
};

bool matcl::is_integer_matrix(const Matrix& A)
{
    value_code vc = A.get_value_code();
    return vc == value_code::v_integer;
};

bool matcl::is_complex_matrix(const Matrix& A)
{
    value_code vc = A.get_value_code();
    return vc == value_code::v_complex || vc == value_code::v_float_complex;
};

bool matcl::is_object_matrix(const Matrix& A)
{
    value_code vc = A.get_value_code();
    return vc == value_code::v_object;
};

bool matcl::is_dense_matrix(const Matrix& A)
{
    struct_code sc = A.get_struct_code();
    return sc == struct_code::struct_dense;
};

bool matcl::is_sparse_matrix(const Matrix& A)
{
    struct_code sc = A.get_struct_code();
    return sc == struct_code::struct_sparse;
};

bool matcl::is_band_matrix(const Matrix& A)
{
    struct_code sc = A.get_struct_code();
    return sc == struct_code::struct_banded;
};

bool matcl::is_scalar_matrix(const Matrix& A)
{
    return A.is_scalar_type();
};

bool matcl::is_same_matrix(const Matrix& A, const Matrix& B)
{
    const details::matrix_base& b1 = details::matrix_data_accesser::get_base(A);
    const details::matrix_base& b2 = details::matrix_data_accesser::get_base(B);

    return details::matrix_base::is_same_matrix(b1,b2);
};

Matrix find(const Matrix& A)
{
    Matrix ret;
    details::eval_find::make<const Matrix&>(A, ret);
    return ret;
};

Matrix find(const Matrix& A,const test_function& t)
{
    Matrix ret;
    details::eval_find_t::make<const Matrix&>(A, ret, t);
    return ret;
};

mat_tup_2 find2(const Matrix& A)
{
    Matrix ret1, ret2;
    details::eval_find2::make<const Matrix&>(A, ret1, ret2);
    return mat_tup_2(ret1,ret2);
};

mat_tup_2 find2(const Matrix& A,const test_function& t)
{
    Matrix ret1, ret2;
    details::eval_find2_t::make<const Matrix&>(A, ret1, ret2, t);
    return mat_tup_2(ret1,ret2);
};

mat_tup_3 find3(const Matrix& A)
{
    Matrix ret1, ret2, ret3;
    details::eval_find3::make<const Matrix&>(A, ret1, ret2, ret3);
    return mat_tup_3(ret1,ret2, ret3);
};

mat_tup_3 find3(const Matrix& A,const test_function& t)
{
    Matrix ret1, ret2, ret3;
    details::eval_find3_t::make<const Matrix&>(A, ret1, ret2, ret3, t);
    return mat_tup_3(ret1,ret2, ret3);
};

Matrix rot90(const Matrix& A, Integer n)
{
    Matrix tmp;
    n = n % 4;

    if (n < 0) 
        n += 4;
    
    switch (n)
    {
        case 0: tmp = A; break;
        case 1: tmp = matcl::flipud(matcl::trans(A)); break;
        case 2: tmp = matcl::flipud(matcl::fliplr(A)); break;
        default: tmp = matcl::trans(matcl::flipud(A));
    }
    return tmp;
};

Matrix sort(const Matrix& A, int dim, bool asceding)
{ 
    Matrix ret;
    details::eval_sort::make<const Matrix&>(A, ret, dim, asceding);
    return ret;
};

mat_tup_2 sort2(const Matrix& A,int dim, bool asceding)
{ 
    Matrix I, X;
    details::eval_sort2::make<const Matrix&>(A, I, X, dim, asceding);
    return mat_tup_2(X, I);
};

Matrix sortrows(const Matrix& A)
{ 
    Matrix ret;
    details::eval_sortrows::make<const Matrix&>(A, ret);
    return ret;
};

mat_tup_2 sortrows2(const Matrix& A)
{ 
    Matrix I, X;
    details::eval_sortrows2::make<const Matrix&>(A, I, X);
    return mat_tup_2(X,I);
};

Matrix sortrows(const Matrix& A, const Matrix& dims0)
{ 
    using Mat_I         = raw::integer_dense;
    Matrix dims         = convert(dims0, matcl::mat_code::integer_dense);
    const Mat_I& mat_d  = dims.get_impl<raw::integer_dense>().make_explicit();

    Matrix ret;
    details::eval_sortrows_dim::make<const Matrix&>(A, ret, mat_d);
    return ret;
};

mat_tup_2 sortrows2(const Matrix& A,const Matrix& dims0)
{ 
    using Mat_I         = raw::integer_dense;
    Matrix dims         = convert(dims0, matcl::mat_code::integer_dense);
    const Mat_I& mat_d  = dims.get_impl<raw::integer_dense>().make_explicit();

    Matrix I, X;
    details::eval_sortrows2_dim::make<const Matrix&>(A, I, X, mat_d);
    return mat_tup_2(X,I);
};

Matrix sortcols(const Matrix& A)
{ 
    Matrix ret;
    details::eval_sortcols::make<const Matrix&>(A, ret);
    return ret;
};

mat_tup_2 sortcols2(const Matrix& A)
{ 
    Matrix I, X;
    details::eval_sortcols2::make<const Matrix&>(A, I, X);
    return mat_tup_2(X,I);
};

Matrix sortcols(const Matrix& A,const Matrix& dims0)
{ 
    using Mat_I         = raw::integer_dense;
    Matrix dims         = convert(dims0, matcl::mat_code::integer_dense);
    const Mat_I& mat_d  = dims.get_impl<raw::integer_dense>().make_explicit();

    Matrix ret;
    details::eval_sortcols_dim::make<const Matrix&>(A, ret, mat_d);
    return ret;
};

mat_tup_2 sortcols2(const Matrix& A,const Matrix& dims0)
{ 
    using Mat_I         = raw::integer_dense;
    Matrix dims         = convert(dims0, matcl::mat_code::integer_dense);
    const Mat_I& mat_d  = dims.get_impl<raw::integer_dense>().make_explicit();

    Matrix I, X;
    details::eval_sortcols2_dim::make<const Matrix&>(A, I, X, mat_d);
    return mat_tup_2(X,I);
};

bool issorted_rows(const Matrix& A)
{
    return details::eval_issorted_rows::make<const Matrix&>(A);
}

bool issorted_cols(const Matrix& A)
{
    return details::eval_issorted_cols::make<const Matrix&>(A);
}

Matrix issorted(const Matrix& A,int dim, bool asceding)
{
    Matrix ret;
    details::eval_issorted::make<const Matrix&>(A, ret, dim, asceding);
    return ret;
};

Matrix horzcat(const Matrix& A,const Matrix& B)
{
    return (mat_row(), A, B); 
};

Matrix horzcat(std::initializer_list<Matrix> mat_list)
{
    mat_row ret;
    for( auto it: mat_list)
        ret.add(it);

    return ret;
};

Matrix horzcat(const std::vector<Matrix>& mat_list)
{
    mat_row ret;
    for( auto it: mat_list)
        ret.add(it);

    return ret;
};

Matrix vertcat(const Matrix& A,const Matrix& B)
{
    return (mat_col(), A, B); 
};

Matrix matcl::vertcat(std::initializer_list<Matrix> mat_list)
{
    mat_col ret;
    for( auto it: mat_list)
        ret.add(it);

    return ret;
};

Matrix matcl::vertcat(const std::vector<Matrix>& mat_list)
{
    mat_col ret;
    for( auto it: mat_list)
        ret.add(it);

    return ret;
};

template<class It>
static Matrix blkdiag_impl(It beg, It end)
{
    Integer M   = 0;
    for( auto it = beg; it != end; ++it)
        M       += it->rows();

    mat_row ret;
    Integer acc_size    = 0;

    for( auto it = beg; it != end; ++it)
    {
        Integer loc_M   = it->rows();
        Integer loc_N   = it->cols();
        
        Matrix top      = spzeros(acc_size, loc_N, 0, it->get_value_code());
        Matrix cur      = *it;
        Matrix bottom   = spzeros(M - acc_size - loc_M, loc_N, 0, it->get_value_code());

        ret.add(mat_col().add(top).add(cur).add(bottom));

        acc_size        += loc_M;
    };

    return ret;
};

Matrix matcl::blkdiag(const Matrix& A, const Matrix& B)
{
    return blkdiag({A,B});
}

Matrix matcl::blkdiag(std::initializer_list<Matrix> mat_list)
{
    return blkdiag_impl(mat_list.begin(), mat_list.end());
};

Matrix matcl::blkdiag(const std::vector<Matrix>& mat_list)
{
    return blkdiag_impl(mat_list.begin(), mat_list.end());
};

template<class T, class S>
T details::convert_scalar_imp<T,S>::eval(const S& val)
{
    return matcl::raw::converter<T,S>::eval(val);
};

};

MACRO_INSTANTIATE_SS_2_F(matcl::details::convert_scalar_imp)
