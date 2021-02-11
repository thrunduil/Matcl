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

#include "matcl-matrep/details/mpl.h"
#include "matcl-core/error/exception_classes.h"

namespace matcl { namespace details
{

// eval functor given by class derived for each matrix representation. Matrix
// representation is unknown at compile time and is given by mat_code. 
// ret - return type of the functor
// derived - a class derived from this class with functions
//     template<class Mat_Rep, ...> static ret eval(Arg... arg)
//     template<class Mat_Rep, ...> static ret eval_scalar(Arg... arg)
// 
// functions:
//         template<class ... Arg>
//          static ret make(mat_code code, Arg&& ... args)
// 
//         this function takes a representation code, and pass to the functor
//         representation type as template argument and additional arguments args.
template<class ret, class derived>
struct extract_type_from_code
{
    template<class ... Arg>
    static ret make(matcl::mat_code code, Arg&& ... args);
};

// eval functor given by class derived for each matrix representation. Matrix
// representation is unknown at compile time and is given by object of class 
// matcl::Matrix. 
// 
// ret     : return type of the functor
// derived : a class derived from this class with functions
//             template<class Mat_Rep, ...> 
//                 static ret eval(Matrix_Type Mat, Rep_Type rep, Arg... arg)
//             template<class Mat_Rep, ...>
//                 static ret eval_scalar(Matrix_Type Mat, Rep_Type rep, Arg... arg)
// 
//         and if Mat_Rep == Matrix&&
// 
//             template<class RValue_Rep, ...> 
//                 static ret eval(RValue_Rep&& rep, Arg... arg)
// 
//         with RValue_Rep = rvalue_holder<Mat_Rep>
// 
// const_access
//         : read only or write access to matrix representation
//         = true:  Rep_Type = const Mat_Rep& 
//         = false: Rep_Type = Mat_Rep&, this makes the matrix unique.
// 
// functions:
//         template<class Matrix_Type, class ... Arg>
//         static ret make(Matrix_Type mat, Arg&& ... args)
// 
//         this function takes an object of class Matrix, extract appropriate
//         representation, and pass representation and additional arguments to the
//         functor. Matrix_Type must be given explicitly and is Matrix&&, Matrix&,
//         or const Matrix&.
//         
//         if Matrix_Type == Matrix&& and const_access == true, then if mat is 
//         unique, then non const representation is taken, and if mat is not unique
//          then const representation is taken. This allows for evaluation inplace for
//         unique temporary matrices.
template<class ret, class derived, bool const_access>
struct extract_type_switch
{
    template<class Matrix_Type, class ... Arg>
    static ret make(typename enable_lvalue<Matrix_Type>::type mat, Arg&& ... args);

    template<class Matrix_Type, class ... Arg>
    static ret make(typename enable_rvalue<Matrix_Type>::type mat, Arg&& ... args);
};

//----------------------------------------------------------------------------
//                          IMPLEMENTATION
//----------------------------------------------------------------------------
template<class ret, class derived>
template<class ... Arg>
ret extract_type_from_code<ret,derived>::make(matcl::mat_code code, Arg&& ... args)
{
    switch(code)
    {
        case mat_code::integer_dense:
        {
            using type = raw::Matrix<Integer,struct_dense>;
            return derived::eval<type>(std::forward<Arg>(args)...);
        }
        case mat_code::real_dense:
        {
            using type = raw::Matrix<Real,struct_dense>;
            return derived::eval<type>(std::forward<Arg>(args)...);
        }
        case mat_code::float_dense:
        {
            using type = raw::Matrix<Float,struct_dense>;
            return derived::eval<type>(std::forward<Arg>(args)...);
        }
        case mat_code::complex_dense:
        {
            using type = raw::Matrix<Complex,struct_dense>;
            return derived::eval<type>(std::forward<Arg>(args)...);
        }
        case mat_code::float_complex_dense:
        {
            using type = raw::Matrix<Float_complex,struct_dense>;
            return derived::eval<type>(std::forward<Arg>(args)...);
        }
        case mat_code::object_dense:
        {
            using type = raw::Matrix<Object,struct_dense>;
            return derived::eval<type>(std::forward<Arg>(args)...);
        }
        case mat_code::integer_sparse:
        {
            using type = raw::Matrix<Integer,struct_sparse>;
            return derived::eval<type>(std::forward<Arg>(args)...);
        }
        case mat_code::real_sparse:
        {
            using type = raw::Matrix<Real,struct_sparse>;
            return derived::eval<type>(std::forward<Arg>(args)...);
        }
        case mat_code::float_sparse:
        {
            using type = raw::Matrix<Float,struct_sparse>;
            return derived::eval<type>(std::forward<Arg>(args)...);
        }
        case mat_code::complex_sparse:
        {
            using type = raw::Matrix<Complex,struct_sparse>;
            return derived::eval<type>(std::forward<Arg>(args)...);
        }
        case mat_code::float_complex_sparse:
        {
            using type = raw::Matrix<Float_complex,struct_sparse>;
            return derived::eval<type>(std::forward<Arg>(args)...);
        }
        case mat_code::object_sparse:
        {
            using type = raw::Matrix<Object,struct_sparse>;
            return derived::eval<type>(std::forward<Arg>(args)...);
        }
        case mat_code::integer_band:
        {
            using type = raw::Matrix<Integer,struct_banded>;
            return derived::eval<type>(std::forward<Arg>(args)...);
        }
        case mat_code::real_band:
        {
            using type = raw::Matrix<Real,struct_banded>;
            return derived::eval<type>(std::forward<Arg>(args)...);
        }
        case mat_code::float_band:
        {
            using type = raw::Matrix<Float,struct_banded>;
            return derived::eval<type>(std::forward<Arg>(args)...);
        }
        case mat_code::complex_band:
        {
            using type = raw::Matrix<Complex,struct_banded>;
            return derived::eval<type>(std::forward<Arg>(args)...);
        }
        case mat_code::float_complex_band:
        {
            using type = raw::Matrix<Float_complex,struct_banded>;
            return derived::eval<type>(std::forward<Arg>(args)...);
        }
        case mat_code::object_band:
        {
            using type = raw::Matrix<Object,struct_banded>;
            return derived::eval<type>(std::forward<Arg>(args)...);
        }
        case mat_code::integer_scalar:
        {
            using type = Integer;
            return derived::eval_scalar<type>(std::forward<Arg>(args)...);
        }
        case mat_code::real_scalar:
        {
            using type = Real;
            return derived::eval_scalar<type>(std::forward<Arg>(args)...);
        }
        case mat_code::float_scalar:
        {
            using type = Float;
            return derived::eval_scalar<type>(std::forward<Arg>(args)...);
        }
        case mat_code::complex_scalar:
        {
            using type = Complex;
            return derived::eval_scalar<type>(std::forward<Arg>(args)...);
        }
        case mat_code::float_complex_scalar:
        {
            using type = Float_complex;
            return derived::eval_scalar<type>(std::forward<Arg>(args)...);
        }
        case mat_code::object_scalar:
        {
            using type = Object;
            return derived::eval_scalar<type>(std::forward<Arg>(args)...);
        }
        default:
            matcl_assert(0,"invalid case");
            throw;
    };
};

template<class MT, class Rep_Type, bool const_access>
struct get_rep{};

template<class MT, class Rep_Type>
struct get_rep<MT,Rep_Type, true>
{
    static const Rep_Type& eval(MT mat)
    {
        return mat.template get_impl<Rep_Type>();
    };
};

template<class MT, class Rep_Type>
struct get_rep<MT,Rep_Type, false>
{
    static Rep_Type& eval(MT mat)
    {
        return mat.template get_impl_unique<Rep_Type>();
    };
};

template<class MT, class Rep_Type, bool const_access>
struct get_scal_rep{};

template<class MT, class Rep_Type>
struct get_scal_rep<MT,Rep_Type, true>
{
    static Rep_Type eval(MT mat)
    {
        return mat.template get_scalar<Rep_Type>();
    };
};

template<class MT, class Rep_Type>
struct get_scal_rep<MT,Rep_Type, false>
{
    static Rep_Type& eval(MT mat)
    {
        return mat.template get_scalar_unique<Rep_Type>();
    };
};

template<class ret, class derived, bool const_access>
template<class Matrix_Type, class ... Arg>
static ret extract_type_switch<ret,derived,const_access>
        ::make(typename enable_lvalue<Matrix_Type>::type mat, Arg&& ... args)
{
    using MT = typename hide_type<Matrix_Type>::type;

    switch(mat.get_matrix_code())
    {
        case mat_code::integer_dense:
        {
            using type = raw::Matrix<Integer,struct_dense>;
            return derived::eval(mat, get_rep<MT,type,const_access>::eval(mat), 
                                std::forward<Arg>(args)...);
        }
        case mat_code::real_dense:
        {
            using type = raw::Matrix<Real,struct_dense>;
            return derived::eval(mat, get_rep<MT,type,const_access>::eval(mat), 
                                std::forward<Arg>(args)...);
        }
        case mat_code::float_dense:
        {
            using type = raw::Matrix<Float,struct_dense>;
            return derived::eval(mat, get_rep<MT,type,const_access>::eval(mat), 
                                std::forward<Arg>(args)...);
        }
        case mat_code::complex_dense:
        {
            using type = raw::Matrix<Complex,struct_dense>;
            return derived::eval(mat, get_rep<MT,type,const_access>::eval(mat), 
                                std::forward<Arg>(args)...);
        }
        case mat_code::float_complex_dense:
        {
            using type = raw::Matrix<Float_complex,struct_dense>;
            return derived::eval(mat, get_rep<MT,type,const_access>::eval(mat), 
                                std::forward<Arg>(args)...);
        }
        case mat_code::object_dense:
        {
            using type = raw::Matrix<Object,struct_dense>;
            return derived::eval(mat, get_rep<MT,type,const_access>::eval(mat), 
                                std::forward<Arg>(args)...);
        }
        case mat_code::integer_sparse:
        {
            using type = raw::Matrix<Integer,struct_sparse>;
            return derived::eval(mat, get_rep<MT,type,const_access>::eval(mat), 
                                std::forward<Arg>(args)...);
        }
        case mat_code::real_sparse:
        {
            using type = raw::Matrix<Real,struct_sparse>;
            return derived::eval(mat, get_rep<MT,type,const_access>::eval(mat), 
                                std::forward<Arg>(args)...);
        }
        case mat_code::float_sparse:
        {
            using type = raw::Matrix<Float,struct_sparse>;
            return derived::eval(mat, get_rep<MT,type,const_access>::eval(mat), 
                                std::forward<Arg>(args)...);
        }
        case mat_code::complex_sparse:
        {
            using type = raw::Matrix<Complex,struct_sparse>;
            return derived::eval(mat, get_rep<MT,type,const_access>::eval(mat), 
                                std::forward<Arg>(args)...);
        }
        case mat_code::float_complex_sparse:
        {
            using type = raw::Matrix<Float_complex,struct_sparse>;
            return derived::eval(mat, get_rep<MT,type,const_access>::eval(mat), 
                                std::forward<Arg>(args)...);
        }
        case mat_code::object_sparse:
        {
            using type = raw::Matrix<Object,struct_sparse>;
            return derived::eval(mat, get_rep<MT,type,const_access>::eval(mat), 
                                std::forward<Arg>(args)...);
        }
        case mat_code::integer_band:
        {
            using type = raw::Matrix<Integer,struct_banded>;
            return derived::eval(mat, get_rep<MT,type,const_access>::eval(mat), 
                                std::forward<Arg>(args)...);
        }
        case mat_code::real_band:
        {
            using type = raw::Matrix<Real,struct_banded>;
            return derived::eval(mat, get_rep<MT,type,const_access>::eval(mat), 
                                std::forward<Arg>(args)...);
        }
        case mat_code::float_band:
        {
            using type = raw::Matrix<Float,struct_banded>;
            return derived::eval(mat, get_rep<MT,type,const_access>::eval(mat), 
                                std::forward<Arg>(args)...);
        }
        case mat_code::complex_band:
        {
            using type = raw::Matrix<Complex,struct_banded>;
            return derived::eval(mat, get_rep<MT,type,const_access>::eval(mat), 
                                std::forward<Arg>(args)...);
        }
        case mat_code::float_complex_band:
        {
            using type = raw::Matrix<Float_complex,struct_banded>;
            return derived::eval(mat, get_rep<MT,type,const_access>::eval(mat), 
                                std::forward<Arg>(args)...);
        }
        case mat_code::object_band:
        {
            using type = raw::Matrix<Object,struct_banded>;
            return derived::eval(mat, get_rep<MT,type,const_access>::eval(mat), 
                                std::forward<Arg>(args)...);
        }
        case mat_code::integer_scalar:
        {
            using type = Integer;
            return derived::eval_scalar(mat,get_scal_rep<MT,type,const_access>::eval(mat), 
                                std::forward<Arg>(args)...);
        }
        case mat_code::real_scalar:
        {
            using type = Real;
            return derived::eval_scalar(mat,get_scal_rep<MT,type,const_access>::eval(mat), 
                                std::forward<Arg>(args)...);
        }
        case mat_code::float_scalar:
        {
            using type = Float;
            return derived::eval_scalar(mat,get_scal_rep<MT,type,const_access>::eval(mat), 
                                std::forward<Arg>(args)...);
        }
        case mat_code::complex_scalar:
        {
            using type = Complex;
            return derived::eval_scalar(mat,get_scal_rep<MT,type,const_access>::eval(mat), 
                                std::forward<Arg>(args)...);
        }
        case mat_code::float_complex_scalar:
        {
            using type = Float_complex;
            return derived::eval_scalar(mat,get_scal_rep<MT,type,const_access>::eval(mat), 
                                std::forward<Arg>(args)...);
        }
        case mat_code::object_scalar:
        {
            using type = Object;
            return derived::eval_scalar(mat,get_scal_rep<MT,type,const_access>::eval(mat), 
                                std::forward<Arg>(args)...);
        }
        default:
            matcl_assert(0,"invalid case");
            throw;
    };
};

template<class Rep_Type>
struct get_rep2
{
    using RValue_Rep_Type = rvalue_holder<Rep_Type>;
    static RValue_Rep_Type eval(Matrix&& mat)
    {
        //mat is temporary and unique, thus rep is also temporary and unique.
        RValue_Rep_Type tmp = mat.template move_impl<Rep_Type>();
        return std::move(tmp);
    };
};

template<class Rep_Type>
struct get_scal2
{
    static Rep_Type eval(Matrix&& mat)
    {
        //mat is temporary and unique, thus rep is also temporary and unique.
        Rep_Type tmp = mat.template get_scalar<Rep_Type>();
        return tmp;
    };
};

template<class ret, class derived, bool const_access>
struct extract_type_switch_nonunique_rvalue
{};

template<class ret, class derived>
struct extract_type_switch_nonunique_rvalue<ret, derived, true>
{
    static const bool const_access = true;

    template<class Matrix_Type, class ... Arg>
    static ret make(typename enable_rvalue<Matrix_Type>::type mat, Arg&& ... args)
    {        
        return extract_type_switch<ret,derived,const_access>
                    ::make<const Matrix_Type&>(mat, std::forward<Arg>(args)...);
    };
};

template<class ret, class derived>
struct extract_type_switch_nonunique_rvalue<ret, derived, false>
{
    static const bool const_access = false;

    template<class Matrix_Type, class ... Arg>
    static ret make(typename enable_rvalue<Matrix_Type>::type mat, Arg&& ... args)
    {
        return extract_type_switch<ret,derived,const_access>
                    ::make<Matrix_Type&>(std::move(mat), std::forward<Arg>(args)...);
    };
};

template<class ret, class derived, bool const_access>
template<class Matrix_Type, class ... Arg>
static ret extract_type_switch<ret,derived,const_access>
        ::make(typename enable_rvalue<Matrix_Type>::type mat, Arg&& ... args)
{
    if (mat.is_unique() == false)
    {
        return extract_type_switch_nonunique_rvalue<ret,derived,const_access>
            ::make<Matrix_Type&&>(std::move(mat), std::forward<Arg>(args)...);
    };

    switch(mat.get_matrix_code())
    {
        case mat_code::integer_dense:
        {
            using type = raw::Matrix<Integer,struct_dense>;
            return derived::eval(get_rep2<type>::eval(std::move(mat)), std::forward<Arg>(args)...);
        }
        case mat_code::real_dense:
        {
            using type = raw::Matrix<Real,struct_dense>;
            return derived::eval(get_rep2<type>::eval(std::move(mat)), std::forward<Arg>(args)...);
        }
        case mat_code::float_dense:
        {
            using type = raw::Matrix<Float,struct_dense>;
            return derived::eval(get_rep2<type>::eval(std::move(mat)), std::forward<Arg>(args)...);
        }
        case mat_code::complex_dense:
        {
            using type = raw::Matrix<Complex,struct_dense>;
            return derived::eval(get_rep2<type>::eval(std::move(mat)), std::forward<Arg>(args)...);
        }
        case mat_code::float_complex_dense:
        {
            using type = raw::Matrix<Float_complex,struct_dense>;
            return derived::eval(get_rep2<type>::eval(std::move(mat)), std::forward<Arg>(args)...);
        }
        case mat_code::object_dense:
        {
            using type = raw::Matrix<Object,struct_dense>;
            return derived::eval(get_rep2<type>::eval(std::move(mat)), std::forward<Arg>(args)...);
        }
        case mat_code::integer_sparse:
        {
            using type = raw::Matrix<Integer,struct_sparse>;
            return derived::eval(get_rep2<type>::eval(std::move(mat)), std::forward<Arg>(args)...);
        }
        case mat_code::real_sparse:
        {
            using type = raw::Matrix<Real,struct_sparse>;
            return derived::eval(get_rep2<type>::eval(std::move(mat)), std::forward<Arg>(args)...);
        }
        case mat_code::float_sparse:
        {
            using type = raw::Matrix<Float,struct_sparse>;
            return derived::eval(get_rep2<type>::eval(std::move(mat)), std::forward<Arg>(args)...);
        }
        case mat_code::complex_sparse:
        {
            using type = raw::Matrix<Complex,struct_sparse>;
            return derived::eval(get_rep2<type>::eval(std::move(mat)), std::forward<Arg>(args)...);
        }
        case mat_code::float_complex_sparse:
        {
            using type = raw::Matrix<Float_complex,struct_sparse>;
            return derived::eval(get_rep2<type>::eval(std::move(mat)), std::forward<Arg>(args)...);
        }
        case mat_code::object_sparse:
        {
            using type = raw::Matrix<Object,struct_sparse>;
            return derived::eval(get_rep2<type>::eval(std::move(mat)), std::forward<Arg>(args)...);
        }
        case mat_code::integer_band:
        {
            using type = raw::Matrix<Integer,struct_banded>;
            return derived::eval(get_rep2<type>::eval(std::move(mat)), std::forward<Arg>(args)...);
        }
        case mat_code::real_band:
        {
            using type = raw::Matrix<Real,struct_banded>;
            return derived::eval(get_rep2<type>::eval(std::move(mat)), std::forward<Arg>(args)...);
        }
        case mat_code::float_band:
        {
            using type = raw::Matrix<Float,struct_banded>;
            return derived::eval(get_rep2<type>::eval(std::move(mat)), std::forward<Arg>(args)...);
        }
        case mat_code::complex_band:
        {
            using type = raw::Matrix<Complex,struct_banded>;
            return derived::eval(get_rep2<type>::eval(std::move(mat)), std::forward<Arg>(args)...);
        }
        case mat_code::float_complex_band:
        {
            using type = raw::Matrix<Float_complex,struct_banded>;
            return derived::eval(get_rep2<type>::eval(std::move(mat)), std::forward<Arg>(args)...);
        }
        case mat_code::object_band:
        {
            using type = raw::Matrix<Object,struct_banded>;
            return derived::eval(get_rep2<type>::eval(std::move(mat)), std::forward<Arg>(args)...);
        }
        case mat_code::integer_scalar:
        {
            using type = Integer;
            return derived::eval_scalar(get_scal2<type>::eval(std::move(mat)), std::forward<Arg>(args)...);
        }
        case mat_code::real_scalar:
        {
            using type = Real;
            return derived::eval_scalar(get_scal2<type>::eval(std::move(mat)), std::forward<Arg>(args)...);
        }
        case mat_code::float_scalar:
        {
            using type = Float;
            return derived::eval_scalar(get_scal2<type>::eval(std::move(mat)), std::forward<Arg>(args)...);
        }
        case mat_code::complex_scalar:
        {
            using type = Complex;
            return derived::eval_scalar(get_scal2<type>::eval(std::move(mat)), std::forward<Arg>(args)...);
        }
        case mat_code::float_complex_scalar:
        {
            using type = Float_complex;
            return derived::eval_scalar(get_scal2<type>::eval(std::move(mat)), std::forward<Arg>(args)...);
        }
        case mat_code::object_scalar:
        {
            using type = Object;
            return derived::eval_scalar(get_scal2<type>::eval(std::move(mat)), std::forward<Arg>(args)...);
        }
        default:
            matcl_assert(0,"invalid case");
            throw;
    };
};

};};
