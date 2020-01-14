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

#include "matcl-matrep/details/isa.h"
#include "matcl-matrep/details/mpl.h"
#include "matcl-internals/func/converter.h"
#include "matcl-internals/container/mat_d.h"
#include "matcl-internals/container/mat_s.h"
#include "matcl-internals/container/mat_b.h"

namespace matcl { namespace raw
{

namespace md = matcl::details;

template<class T> struct get_value_type         { using type = typename T::value_type; };
template<> struct get_value_type<Integer>       { using type = Integer; };
template<> struct get_value_type<Real>          { using type = Real; };
template<> struct get_value_type<Float>         { using type = Float; };
template<> struct get_value_type<Complex>       { using type = Complex; };
template<> struct get_value_type<Float_complex> { using type = Float_complex; };
template<> struct get_value_type<Object>        { using type = Object; };

template<class T> struct get_struct_type        { using type = typename T::struct_type; };
template<> struct get_struct_type<Integer>      { using type = struct_scalar; };
template<> struct get_struct_type<Real>         { using type = struct_scalar; };
template<> struct get_struct_type<Float>        { using type = struct_scalar; };
template<> struct get_struct_type<Complex>      { using type = struct_scalar; };
template<> struct get_struct_type<Float_complex>{ using type = struct_scalar; };
template<> struct get_struct_type<Object>       { using type = struct_scalar; };

template<class T> struct get_struct_type_dense        { using type = struct_dense; };
template<> struct get_struct_type_dense<Integer>      { using type = struct_scalar; };
template<> struct get_struct_type_dense<Real>         { using type = struct_scalar; };
template<> struct get_struct_type_dense<Float>        { using type = struct_scalar; };
template<> struct get_struct_type_dense<Complex>      { using type = struct_scalar; };
template<> struct get_struct_type_dense<Float_complex>{ using type = struct_scalar; };
template<> struct get_struct_type_dense<Object>       { using type = struct_scalar; };

template<class V1, class V2>
struct ret_value_type
{    
    static const bool condval = std::is_same<V1,Object>::value || std::is_same<V2,Object>::value;
    using vtype = V1;
    using type  = typename matcl::details::select_if<condval, Object, vtype>::type;
};

template<class V1, class V2>
struct ret_value_type_int
{        
    using V1R   = typename md::real_type<V1>::type;
    using V2R   = typename md::real_type<V2>::type;
    using VR    = typename md::unify_types<V1R,V2R>::type;
    using VC    = matcl::complex<VR>;

    static const bool iscompl   = md::is_complex<V1>::value && md::is_object<VR>::value == false;
    using type  = typename matcl::details::select_if<iscompl, VC, VR>::type;
};

template<class V, class S>
struct ret_matrix_type
{
    using type = matcl::raw::Matrix<V,S>;
};

template<class V>
struct ret_matrix_type<V,struct_scalar>
{
    using type = V;
};

template<class T1, class T2, bool is_obj = md::is_object<T1>::value || md::is_object<T2>::value>
struct make_int_obj_ret
{
    using type = Integer;
};

template<class T1, class T2>
struct make_int_obj_ret<T1, T2, true>
{
    using type = Object;
};

template<class T1, class T2>
struct val_type_corrector_tr
{
    using val_1     = typename get_value_type<T1>::type;
    using val_2     = typename get_value_type<T2>::type;

    using str_1     = typename get_struct_type<T1>::type;
    using str_2     = typename get_struct_type<T2>::type;
    
    using val_ret_1 = typename ret_value_type<val_1,val_2>::type;
    using val_ret_2 = typename ret_value_type<val_2,val_1>::type;

    using type_1    = typename ret_matrix_type<val_ret_1,str_1>::type;
    using type_2    = typename ret_matrix_type<val_ret_2,str_2>::type;
};

template<class T1, class T2>
struct val_type_corrector_int_tr
{
    using val_1     = typename get_value_type<T1>::type;
    using val_2     = typename get_value_type<T2>::type;

    using str_1     = typename get_struct_type<T1>::type;
    using str_2     = typename get_struct_type<T2>::type;
    
    using val_ret_1 = typename ret_value_type_int<val_1,val_2>::type;
    using val_ret_2 = typename ret_value_type_int<val_2,val_1>::type;

    using type_1    = typename ret_matrix_type<val_ret_1,str_1>::type;
    using type_2    = typename ret_matrix_type<val_ret_2,str_2>::type;
};

template<class T1, class T2>
struct val_type_corrector_diag_tr
{
    using val_1     = typename get_value_type<T1>::type;
    using val_2     = typename get_value_type<T2>::type;

    using str_1     = typename get_struct_type<T1>::type;
    using str_2     = typename get_struct_type<T2>::type;
    
    using val_ret_1 = typename md::unify_types<val_1,val_2>::type;
    using val_ret_2 = val_ret_1;

    using type_1    = typename ret_matrix_type<val_ret_1,str_1>::type;
    using type_2    = typename ret_matrix_type<val_ret_2,str_2>::type;
};

template<class T1, class T2>
struct val_type_corrector_diag_dense_tr
{
    using val_1     = typename get_value_type<T1>::type;
    using val_2     = typename get_value_type<T2>::type;

    using str_1     = typename get_struct_type_dense<T1>::type;
    using str_2     = typename get_struct_type_dense<T2>::type;
    
    using val_ret_1 = typename md::unify_types<val_1,val_2>::type;
    using val_ret_2 = val_ret_1;

    using type_1    = typename ret_matrix_type<val_ret_1,str_1>::type;
    using type_2    = typename ret_matrix_type<val_ret_2,str_2>::type;
};

template<class T1, class T2>
struct val_type_corrector_dense_tr
{
    using val_1     = typename get_value_type<T1>::type;
    using val_2     = typename get_value_type<T2>::type;

    using str_1     = typename get_struct_type_dense<T1>::type;
    using str_2     = typename get_struct_type_dense<T2>::type;
    
    using val_ret_1 = typename ret_value_type<val_1,val_2>::type;
    using val_ret_2 = typename ret_value_type<val_2,val_1>::type;

    using type_1    = typename ret_matrix_type<val_ret_1,str_1>::type;
    using type_2    = typename ret_matrix_type<val_ret_2,str_2>::type;
};

template<class T1, class T2>
struct val_type_corrector_bool_tr
{
    using val_1     = typename get_value_type<T1>::type;
    using val_2     = typename get_value_type<T2>::type;

    using str_1     = typename get_struct_type<T1>::type;
    using str_2     = typename get_struct_type<T2>::type;
    
    using val_ret   = typename make_int_obj_ret<val_1,val_2>::type;
    using val_ret_1 = val_ret;
    using val_ret_2 = val_ret;

    using type_1    = typename ret_matrix_type<val_ret_1,str_1>::type;
    using type_2    = typename ret_matrix_type<val_ret_2,str_2>::type;
};

template<class T1, class T2, class T1_ret, class T2_ret>
struct val_type_corrector_impl
{
    public:
        static T1_ret convert_1(const T1& mat)
        {
            using converter = matcl::raw::converter<T1_ret,T1>;
            return converter::eval(mat);
        };
        static T2_ret convert_2(const T2& mat)
        {
            using converter = matcl::raw::converter<T2_ret,T2>;
            return converter::eval(mat);
        };

    private:
        static T1_ret convert_1(T1&& mat);
        static T2_ret convert_2(T2&& mat);
};

template<class T1, class T2>
struct val_type_corrector_impl<T1,T2,T1,T2>
{
    public:
        static const T1& convert_1(const T1& mat)
        {
            return mat;
        };

        static const T2& convert_2(const T2& mat)
        {
            return mat;
        };

    private:
        static const T1& convert_1(T1&& mat);
        static const T2& convert_2(T2&& mat);
};

template<class T1, class T2, class T2_ret>
struct val_type_corrector_impl<T1,T2,T1,T2_ret>
{
    public:
        static const T1& convert_1(const T1& mat)
        {
            return mat;
        };

        static T2_ret convert_2(const T2& mat)
        {
            using converter = matcl::raw::converter<T2_ret,T2>;
            return converter::eval(mat);
        };

    private:
        static const T1&    convert_1(T1&& mat);
        static T2_ret       convert_2(T2&& mat);
};

template<class In, class Ret, bool Is_scal, bool Is_obj>
struct cast_bool_impl
{};

template<class In, class Ret>
struct cast_bool_impl<In, Ret, false, false>
{
    using In_val        = typename In::value_type;
    using Ret_val       = typename Ret::value_type;
    using In_str        = typename In::struct_type;
    using Ret_holder    = const_matrix<Ret>;

    static void eval(Ret_holder& ret, const In& mat)
    {
        using matrix_type = Matrix<In_val, In_str>;

        matcl::Matrix ret_m;
        mrd::unary_helper_impl<matrix_type>::eval_is_true(ret_m, mat);

        ret_m = matcl::convert(ret_m, Ret::matrix_code);
        rebind<Ret, Ret_holder>::eval(ret, ret_m.get_impl<Ret>());
    };
};

template<class In, class Ret>
struct cast_bool_impl<In, Ret, false, true>
{
    using In_val        = typename In::value_type;
    using Ret_val       = typename Ret::value_type;
    using In_str        = typename In::struct_type;
    using Ret_holder    = const_matrix<Ret>;

    static void eval(Ret_holder& ret, const In& mat)
    {
        rebind<Ret, Ret_holder>::eval(ret, md::converter<Ret>::eval(mat));
    };
};

template<class In, class Ret>
struct cast_bool_impl<In, Ret, true, false>
{
    static Ret eval(const In& scal)
    {
        return cast_bool(scal);
    };
};

template<class In, class Ret>
struct cast_bool_impl<In, Ret, true, true>
{
    static Ret eval(const In& scal)
    {
        return Object(scal);
    };
};

template<class T1, class T2, class T1_ret>
struct val_type_corrector_impl<T1,T2,T1_ret, T2>
{
    public:
        static T1_ret convert_1(const T1& mat)
        {
            using converter = matcl::raw::converter<T1_ret,T1>;
            return converter::eval(mat);
        };

        static const T2&  convert_2(const T2& mat)
        {
            return mat;
        };

    private:
        static T1_ret       convert_1(T1&& mat);
        static const T2&    convert_2(T2&& mat);
};

template<class T1, class T2, class T1_ret, class T2_ret>
struct val_type_corrector_bool_impl
{
    private:
        static const bool is_scal1   = md::is_scalar<T1>::value;
        static const bool is_scal2   = md::is_scalar<T2>::value;
        static const bool is_obj1    = md::is_object<T1_ret>::value;
        static const bool is_obj2    = md::is_object<T2_ret>::value;

    public:
        static T1_ret convert_1(const T1& mat)
        {
            return cast_bool_impl<T1,T1_ret, is_scal1, is_obj1>::eval(mat);
        };

        static T2_ret convert_2(const T2& mat)
        {
            return cast_bool_impl<T2,T2_ret, is_scal2, is_obj2>::eval(mat);
        };

    private:
        static T1_ret convert_1(T1&& mat);
        static T2_ret convert_2(T2&& mat);
};

template<class T1, class T2>
struct val_type_corrector_bool_impl<T1,T2,T1,T2>
{
    public:
        static const T1& convert_1(const T1& mat)
        {
            return mat;
        };

        static const T2& convert_2(const T2& mat)
        {
            return mat;
        };

    private:
        static const T1& convert_1(T1&& mat);
        static const T2& convert_2(T2&& mat);
};

template<class T1, class T2, class T2_ret>
struct val_type_corrector_bool_impl<T1,T2,T1,T2_ret>
{
    private:
        static const bool is_scal2   = md::is_scalar<T2>::value;
        static const bool is_obj2    = md::is_object<T2_ret>::value;

    public:
        static const T1& convert_1(const T1& mat)
        {
            return mat;
        };

        static T2_ret convert_2(const T2& mat)
        {
            return cast_bool_impl<T2,T2_ret, is_scal2, is_obj2>::eval(mat);
        };

    private:
        static const T1&    convert_1(T1&& mat);
        static T2_ret       convert_2(T2&& mat);
};

template<class T1, class T2, class T1_ret>
struct val_type_corrector_bool_impl<T1,T2,T1_ret, T2>
{
    private:
        static const bool is_scal1   = md::is_scalar<T1>::value;
        static const bool is_obj1    = md::is_object<T1_ret>::value;

    public:
        static T1_ret convert_1(const T1& mat)
        {
            return cast_bool_impl<T1,T1_ret,is_scal1, is_obj1>::eval(mat);
        };

        static const T2&  convert_2(const T2& mat)
        {
            return mat;
        };

    private:
        static T1_ret       convert_1(T1&& mat);
        static const T2&    convert_2(T2&& mat);
};


template<class T1, class T2>
struct val_type_corrector : public val_type_corrector_impl<T1,T2,
                                    typename val_type_corrector_tr<T1,T2>::type_1,
                                    typename val_type_corrector_tr<T1,T2>::type_2>
                          , public val_type_corrector_tr<T1,T2>
{};

template<class T1, class T2>
struct val_type_corrector_int : public val_type_corrector_impl<T1,T2,
                                    typename val_type_corrector_int_tr<T1,T2>::type_1,
                                    typename val_type_corrector_int_tr<T1,T2>::type_2>
                              , public val_type_corrector_int_tr<T1,T2>
{};

template<class T1, class T2>
struct val_type_corrector_diag: public val_type_corrector_impl<T1,T2,
                                    typename val_type_corrector_diag_tr<T1,T2>::type_1,
                                    typename val_type_corrector_diag_tr<T1,T2>::type_2>
                              , public val_type_corrector_diag_tr<T1,T2>
{};

template<class T1, class T2>
struct val_type_corrector_diag_dense: public val_type_corrector_impl<T1,T2,
                                        typename val_type_corrector_diag_dense_tr<T1,T2>::type_1,
                                        typename val_type_corrector_diag_dense_tr<T1,T2>::type_2>
                                    , public val_type_corrector_diag_dense_tr<T1,T2>
{};

template<class T1, class T2>
struct val_type_corrector_dense: public val_type_corrector_impl<T1,T2,
                                        typename val_type_corrector_dense_tr<T1,T2>::type_1,
                                        typename val_type_corrector_dense_tr<T1,T2>::type_2>
                               , public val_type_corrector_dense_tr<T1,T2>
{};

template<class T1, class T2>
struct val_type_corrector_bool : public val_type_corrector_bool_impl<T1,T2,
                                    typename val_type_corrector_bool_tr<T1,T2>::type_1,
                                    typename val_type_corrector_bool_tr<T1,T2>::type_2>
                              , public val_type_corrector_bool_tr<T1,T2>
{};

}};