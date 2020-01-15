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

#include "matcl-matrep/details/fwd_decls.h"
#include "matcl-matrep/details/matrix.inl"
#include "matcl-matrep/details/promote_type.h"

namespace matcl { namespace details
{

struct ver_static{};
struct ver_nonstatic{};


#define EXPAND_CASE_MAT(macro,arg)              \
    macro(mat_code::integer_dense,arg)          \
    macro(mat_code::float_dense,arg)            \
    macro(mat_code::real_dense,arg)             \
    macro(mat_code::float_complex_dense,arg)    \
    macro(mat_code::complex_dense,arg)          \
    macro(mat_code::object_dense,arg)           \
    macro(mat_code::integer_sparse,arg)         \
    macro(mat_code::float_sparse,arg)           \
    macro(mat_code::real_sparse,arg)            \
    macro(mat_code::float_complex_sparse,arg)   \
    macro(mat_code::complex_sparse,arg)         \
    macro(mat_code::object_sparse,arg)          \
    macro(mat_code::integer_band,arg)           \
    macro(mat_code::float_band,arg)             \
    macro(mat_code::real_band,arg)              \
    macro(mat_code::float_complex_band,arg)     \
    macro(mat_code::complex_band,arg)           \
    macro(mat_code::object_band,arg)       

#define EXPAND_CASE_MAT_2(macro,arg)            \
    macro(mat_code::integer_dense,arg)          \
    macro(mat_code::float_dense,arg)            \
    macro(mat_code::real_dense,arg)             \
    macro(mat_code::float_complex_dense,arg)    \
    macro(mat_code::complex_dense,arg)          \
    macro(mat_code::object_dense,arg)           \
    macro(mat_code::integer_sparse,arg)         \
    macro(mat_code::float_sparse,arg)           \
    macro(mat_code::real_sparse,arg)            \
    macro(mat_code::float_complex_sparse,arg)   \
    macro(mat_code::complex_sparse,arg)         \
    macro(mat_code::object_sparse,arg)          \
    macro(mat_code::integer_band,arg)           \
    macro(mat_code::float_band,arg)             \
    macro(mat_code::real_band,arg)              \
    macro(mat_code::float_complex_band,arg)     \
    macro(mat_code::complex_band,arg)           \
    macro(mat_code::object_band,arg)       

#define EXPAND_CASE_SCAL(macro,arg)             \
    macro(mat_code::integer_scalar,arg)         \
    macro(mat_code::float_scalar,arg)           \
    macro(mat_code::real_scalar,arg)            \
    macro(mat_code::float_complex_scalar,arg)   \
    macro(mat_code::complex_scalar,arg)         \
    macro(mat_code::object_scalar,arg)     

#define EXPAND_CASE_SCAL_2(macro,arg)           \
    macro(mat_code::integer_scalar,arg)         \
    macro(mat_code::float_scalar,arg)           \
    macro(mat_code::real_scalar,arg)            \
    macro(mat_code::float_complex_scalar,arg)   \
    macro(mat_code::complex_scalar,arg)         \
    macro(mat_code::object_scalar,arg)     

#define EXPAND_CASE_MAT_MAT_IMPL(code_2,code_1) \
    case (int)code_1*N + (int)code_2:           \
    return make_mat_mat<code_to_type<code_1>::type,code_to_type<code_2>::type>(A,B, std::forward<Args>(args)...);

#define EXPAND_CASE_MAT_SCAL_IMPL(code_2,code_1) \
    case (int)code_1*N + (int)code_2:            \
    return make_mat_scal<code_to_type<code_1>::type,code_to_type<code_2>::type>(A,B, std::forward<Args>(args)...);

#define EXPAND_CASE_SCAL_MAT_IMPL(code_2,code_1) \
    case (int)code_1*N + (int)code_2:            \
    return make_scal_mat<code_to_type<code_1>::type,code_to_type<code_2>::type>(A,B, std::forward<Args>(args)...);

#define EXPAND_CASE_SCAL_SCAL_IMPL(code_2,code_1) \
    case (int)code_1*N + (int)code_2:             \
    return make_scal_scal<code_to_type<code_1>::type,code_to_type<code_2>::type>(A,B, std::forward<Args>(args)...);

#define EXPAND_CASE_MAT_MAT_2(code_1,arg)       \
EXPAND_CASE_MAT_2(EXPAND_CASE_MAT_MAT_IMPL,code_1)    

#define EXPAND_CASE_MAT_SCAL_2(code_1,arg)       \
EXPAND_CASE_SCAL_2(EXPAND_CASE_MAT_SCAL_IMPL,code_1)    

#define EXPAND_CASE_SCAL_MAT_2(code_1,arg)       \
EXPAND_CASE_MAT_2(EXPAND_CASE_SCAL_MAT_IMPL,code_1)    

#define EXPAND_CASE_SCAL_SCAL_2(code_1,arg)       \
EXPAND_CASE_SCAL_2(EXPAND_CASE_SCAL_SCAL_IMPL,code_1)    


#define EXPAND_CASE_MAT_MAT             \
EXPAND_CASE_MAT(EXPAND_CASE_MAT_MAT_2,)     

#define EXPAND_CASE_MAT_SCAL             \
EXPAND_CASE_MAT(EXPAND_CASE_MAT_SCAL_2,)     

#define EXPAND_CASE_SCAL_MAT             \
EXPAND_CASE_SCAL(EXPAND_CASE_SCAL_MAT_2,)     

#define EXPAND_CASE_SCAL_SCAL             \
EXPAND_CASE_SCAL(EXPAND_CASE_SCAL_SCAL_2,)     

template<class ret, class derived, template<class T1, class T2> class val_corrector,
        class T1, class T2, class ver>
struct extract_type2_switch_impl
{
    template<class ... Args>
    static ret eval_mat_mat(const T1& A, const T2& B, Args&& ... args)
    {
        using corrector = val_corrector<T1,T2>;
        using MT1       = typename corrector::type_1;
        using MT2       = typename corrector::type_2;

        matcl::Matrix tmp1;
        matcl::Matrix tmp2;

        return derived::template eval_mat_mat<MT1,MT2>(corrector::convert_1(A, tmp1), 
                            corrector::convert_2(B, tmp2), std::forward<Args>(args)...);
    };

    template<class ... Args>
    static ret eval_mat_scal(const T1& A, const T2& B, Args&& ... args)
    {
        using corrector = val_corrector<T1,T2>;
        using MT1       = typename corrector::type_1;
        using MT2       = typename corrector::type_2;

        matcl::Matrix tmp1;
        matcl::Matrix tmp2;

        return derived::template eval_mat_scal<MT1,MT2>(corrector::convert_1(A, tmp1),
                            corrector::convert_2(B, tmp2), std::forward<Args>(args)...);
    };

    template<class ... Args>
    static ret eval_scal_mat(const T1& A, const T2& B, Args&& ... args)
    {
        using corrector = val_corrector<T1,T2>;
        using MT1       = typename corrector::type_1;
        using MT2       = typename corrector::type_2;

        matcl::Matrix tmp1;
        matcl::Matrix tmp2;

        return derived::template eval_scal_mat<MT1,MT2>(corrector::convert_1(A, tmp1),
                            corrector::convert_2(B, tmp2), std::forward<Args>(args)...);
    };

    template<class ... Args>
    static ret eval_scal_scal(const T1& A, const T2& B, Args&& ... args)
    {
        using corrector = val_corrector<T1,T2>;
        using MT1       = typename corrector::type_1;
        using MT2       = typename corrector::type_2;

        matcl::Matrix tmp1;
        matcl::Matrix tmp2;

        return derived::template eval_scal_scal<MT1,MT2>(corrector::convert_1(A, tmp1),
                            corrector::convert_2(B, tmp2), std::forward<Args>(args)...);
    };
};

template<class ret, class derived, template<class T1, class T2> class val_corrector,
        class T1, class T2>
struct extract_type2_switch_impl<ret,derived,val_corrector,T1,T2,ver_nonstatic>
{
    template<class ... Args>
    static ret eval_mat_mat(derived* d, const T1& A, const T2& B, Args&& ... args)
    {
        using corrector = val_corrector<T1,T2>;
        using MT1       = typename corrector::type_1;
        using MT2       = typename corrector::type_2;

        matcl::Matrix tmp1;
        matcl::Matrix tmp2;

        return d->template eval_mat_mat<MT1,MT2>(corrector::convert_1(A, tmp1), 
                            corrector::convert_2(B, tmp2), std::forward<Args>(args)...);
    };

    template<class ... Args>
    static ret eval_mat_scal(derived* d, const T1& A, const T2& B, Args&& ... args)
    {
        using corrector = val_corrector<T1,T2>;
        using MT1       = typename corrector::type_1;
        using MT2       = typename corrector::type_2;

        matcl::Matrix tmp1;
        matcl::Matrix tmp2;

        return d->template eval_mat_scal<MT1,MT2>(corrector::convert_1(A, tmp1), 
                            corrector::convert_2(B, tmp2), std::forward<Args>(args)...);
    };

    template<class ... Args>
    static ret eval_scal_mat(derived* d, const T1& A, const T2& B, Args&& ... args)
    {
        using corrector = val_corrector<T1,T2>;
        using MT1       = typename corrector::type_1;
        using MT2       = typename corrector::type_2;

        matcl::Matrix tmp1;
        matcl::Matrix tmp2;

        return d->template eval_scal_mat<MT1,MT2>(corrector::convert_1(A, tmp1), 
                            corrector::convert_2(B, tmp2), std::forward<Args>(args)...);
    };

    template<class ... Args>
    static ret eval_scal_scal(derived* d, const T1& A, const T2& B, Args&& ... args)
    {
        using corrector = val_corrector<T1,T2>;
        using MT1       = typename corrector::type_1;
        using MT2       = typename corrector::type_2;

        matcl::Matrix tmp1;
        matcl::Matrix tmp2;

        return d->template eval_scal_scal<MT1,MT2>(corrector::convert_1(A, tmp1),
                            corrector::convert_2(B, tmp2), std::forward<Args>(args)...);
    };
};

template<class ret, class derived,bool iso,class T1, class T2>
struct extract_type2_switch_nc_impl
{
    static ret eval_mat_mat(derived* d, Matrix& h, T1& A, const T2& B)
    {
        return d->eval_mat_mat(h,A,B);
    };

    static ret eval_mat_scal(derived* d, Matrix& h, T1& A, const T2& B)
    {
        return d->eval_mat_scal(h,A,B);
    };

    static ret eval_scal_mat(derived* d, Matrix& h, T1& A, const T2& B)
    {
        return d->eval_scal_mat(h,A,B);
    };

    static ret eval_scal_scal(derived* d, Matrix& h, T1& A, const T2& B)
    {
        return d->eval_scal_scal(h,A,B);
    };
};

template<class ret, class derived,class T1, class T2>
struct extract_type2_switch_nc_impl<ret,derived,true,T1,T2>
{
    static ret eval_mat_mat(derived* d, Matrix& h, T1& A, const T2& B)
    {
        using corrector = matcl::raw::val_type_corrector<T1,T2>;
        using MT1       = typename corrector::type_1;
        using MT2       = typename corrector::type_2;

        {
            Matrix tmp;
            h = Matrix(corrector::convert_1(A, tmp), false);
        };

        Matrix tmp_B;
        return d->template eval_mat_mat<MT1,MT2>(h, h.get_impl_unique<MT1>(), 
                                                 corrector::convert_2(B, tmp_B));
    };

    static ret eval_mat_scal(derived* d, Matrix& h, T1& A, const T2& B)
    {
        using corrector = matcl::raw::val_type_corrector<T1,T2>;
        using MT1       = typename corrector::type_1;
        using MT2       = typename corrector::type_2;

        {
            Matrix tmp;
            h = Matrix(corrector::convert_1(A, tmp), false);
        };

        Matrix tmp_B;
        return d->template  eval_mat_scal<MT1,MT2>(h, h.get_impl_unique<MT1>(), 
                                                   corrector::convert_2(B, tmp_B));
    };

    static ret eval_scal_mat(derived* d, Matrix& h, T1& A, const T2& B)
    {
        using corrector = matcl::raw::val_type_corrector<T1,T2>;
        using MT1       = typename corrector::type_1;
        using MT2       = typename corrector::type_2;

        {
            Matrix tmp;
            h = Matrix(corrector::convert_1(A, tmp), false);
        };

        Matrix tmp_B;
        return d->template eval_scal_mat<MT1,MT2>(h,h.get_scalar_unique<MT1>(), 
                                                  corrector::convert_2(B, tmp_B));
    };

    static ret eval_scal_scal(derived* d, Matrix& h, T1& A, const T2& B)
    {
        using corrector = matcl::raw::val_type_corrector<T1,T2>;
        using MT1       = typename corrector::type_1;
        using MT2       = typename corrector::type_2;

        {
            Matrix tmp;
            h = Matrix(corrector::convert_1(A ,tmp), false);
        };
        
        Matrix tmp_B;
        return d->template eval_scal_scal<MT1,MT2>(h,h.get_scalar_unique<MT1>(), 
                                                   corrector::convert_2(B, tmp_B));
    };
};

template<class ret, class derived, template<class T1, class T2> class val_corrector, 
         class ver = ver_static>
struct extract_type2_switch
{
    template<class ... Args>
    static ret make(const Matrix& A, const Matrix& B, Args&& ... args)
    {
        static const size_t N = macro_last_matrix_type_code+1;
        int code = (int)A.get_matrix_code()*N + (int)B.get_matrix_code();
        switch(code)
        {
            EXPAND_CASE_MAT_MAT
            EXPAND_CASE_MAT_SCAL
            EXPAND_CASE_SCAL_MAT
            EXPAND_CASE_SCAL_SCAL
            default:
                matcl_assert(0,"unknown case");
                throw;
        };
    };		

    template<class T1, class T2, class ... Args>
    static ret make_mat_mat(const Matrix& A, const Matrix& B, Args&& ... args)
    {
        return extract_type2_switch_impl<ret,derived,val_corrector,T1,T2,ver_static>
            ::eval_mat_mat(A.get_impl<T1>(),B.get_impl<T2>(), std::forward<Args>(args)...);
    };

    template<class T1, class T2, class ... Args>
    static ret make_mat_scal(const Matrix& A, const Matrix& B, Args&& ... args)
    {
        return extract_type2_switch_impl<ret,derived,val_corrector,T1,T2,ver_static>
            ::eval_mat_scal(A.get_impl<T1>(),B.get_scalar<T2>(), std::forward<Args>(args)...);
    };

    template<class T1, class T2, class ... Args>
    static ret make_scal_mat(const Matrix& A, const Matrix& B, Args&& ... args)
    {
        return extract_type2_switch_impl<ret,derived,val_corrector,T1,T2,ver_static>
            ::eval_scal_mat(A.get_scalar<T1>(),B.get_impl<T2>(), std::forward<Args>(args)...);
    };

    template<class T1, class T2, class ... Args>
    static ret make_scal_scal(const Matrix& A, const Matrix& B, Args&& ... args)
    {
        return extract_type2_switch_impl<ret,derived,val_corrector,T1,T2,ver_static>
            ::eval_scal_scal(A.get_scalar<T1>(),B.get_scalar<T2>(), std::forward<Args>(args)...);
    };
};

template<class ret, class derived,template<class T1, class T2> class val_corrector>
struct extract_type2_switch<ret,derived,val_corrector,ver_nonstatic>
{
    template<class ... Args>
    ret make(const Matrix& A, const Matrix& B, Args&& ... args)
    {
        static const size_t N = macro_last_matrix_type_code+1;
        int code = (int)A.get_matrix_code()*N + (int)B.get_matrix_code();
        switch(code)
        {
            EXPAND_CASE_MAT_MAT
            EXPAND_CASE_MAT_SCAL
            EXPAND_CASE_SCAL_MAT
            EXPAND_CASE_SCAL_SCAL
            default:
                matcl_assert(0,"unknown case");
                throw;
        };
    };		

    template<class T1, class T2, class ... Args>
    ret make_mat_mat(const Matrix& A, const Matrix& B, Args&& ... args)
    {
        return extract_type2_switch_impl<ret,derived,val_corrector,T1,T2,ver_nonstatic>
            ::eval_mat_mat(static_cast<derived*>(this),A.get_impl<T1>(),B.get_impl<T2>(),
                           std::forward<Args>(args)...);
    };

    template<class T1, class T2, class ... Args>
    ret make_mat_scal(const Matrix& A, const Matrix& B, Args&& ... args)
    {
        return extract_type2_switch_impl<ret,derived,val_corrector,T1,T2,ver_nonstatic>
            ::eval_mat_scal(static_cast<derived*>(this),A.get_impl<T1>(),B.get_scalar<T2>(),
                            std::forward<Args>(args)...);
    };

    template<class T1, class T2, class ... Args>
    ret make_scal_mat(const Matrix& A, const Matrix& B, Args&& ... args)
    {
        return extract_type2_switch_impl<ret,derived,val_corrector,T1,T2,ver_nonstatic>
            ::eval_scal_mat(static_cast<derived*>(this),A.get_scalar<T1>(),B.get_impl<T2>(),
                            std::forward<Args>(args)...);
    };

    template<class T1, class T2, class ... Args>
    ret make_scal_scal(const Matrix& A, const Matrix& B, Args&& ... args)
    {
        return extract_type2_switch_impl<ret,derived,val_corrector,T1,T2,ver_nonstatic>
            ::eval_scal_scal(static_cast<derived*>(this),A.get_scalar<T1>(),B.get_scalar<T2>(),
                             std::forward<Args>(args)...);
    };
};

template<class T>
struct is_object_mat
{
    static const bool value = std::is_same<T,Object>::value;
};

template<class V,class S>
struct is_object_mat<raw::Matrix<V,S>>
{
    static const bool value = std::is_same<V,Object>::value;
};

template<class ret, class derived>
struct extract_type2_switch_nc
{
    template<class ... Args>
    ret make(Matrix& A, const Matrix& B, Args&& ... args)
    {
        static const size_t N = macro_last_matrix_type_code+1;
        int code = (int)A.get_matrix_code()*N + (int)B.get_matrix_code();
        switch(code)
        {
            EXPAND_CASE_MAT_MAT
            EXPAND_CASE_MAT_SCAL
            EXPAND_CASE_SCAL_MAT
            EXPAND_CASE_SCAL_SCAL
            default:
                matcl_assert(0,"unknown case");
                throw;
        };
    };		

    template<class T1, class T2, class ... Args>
    ret make_mat_mat(Matrix& A, const Matrix& B, Args&& ... args)
    {
        static const bool iso = is_object_mat<T1>::value || is_object_mat<T2>::value;
        return extract_type2_switch_nc_impl<ret,derived,iso,T1,T2>
            ::eval_mat_mat(static_cast<derived*>(this), A, A.get_impl_unique<T1>(),
                           B.get_impl<T2>(), std::forward<Args>(args)...);
    };

    template<class T1, class T2, class ... Args>
    ret make_mat_scal(Matrix& A, const Matrix& B, Args&& ... args)
    {
        static const bool iso = is_object_mat<T1>::value || is_object_mat<T2>::value;
        return extract_type2_switch_nc_impl<ret,derived,iso,T1,T2>
            ::eval_mat_scal(static_cast<derived*>(this), A, A.get_impl_unique<T1>(),
                            B.get_scalar<T2>(), std::forward<Args>(args)...);
    };

    template<class T1, class T2, class ... Args>
    ret make_scal_mat(Matrix& A, const Matrix& B, Args&& ... args)
    {
        static const bool iso = is_object_mat<T1>::value || is_object_mat<T2>::value;
        return extract_type2_switch_nc_impl<ret,derived,iso,T1,T2>
            ::eval_scal_mat(static_cast<derived*>(this),A,A.get_scalar_unique<T1>(),B.get_impl<T2>(),
                            std::forward<Args>(args)...);
    };

    template<class T1, class T2, class ... Args>
    ret make_scal_scal(Matrix& A, const Matrix& B, Args&& ... args)
    {
        static const bool iso = is_object_mat<T1>::value || is_object_mat<T2>::value;
        return extract_type2_switch_nc_impl<ret,derived,iso,T1,T2>
            ::eval_scal_scal(static_cast<derived*>(this),A,A.get_scalar_unique<T1>(),B.get_scalar<T2>(),
                             std::forward<Args>(args)...);
    };
};

};};