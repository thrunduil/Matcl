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

#include "matcl-matrep/objects/details/type_info_object.h"
#include "matcl-internals/container/mat_d.h"
#include "matcl-internals/container/mat_b.h"
#include "matcl-internals/container/mat_s.h"

namespace matcl { namespace ti { namespace details
{

template<class T, class TI>
struct get_ti_type_impl
{
    static TI eval(const T&)
    {
        return TI();
    };
};

template<class T>
struct get_ti_type_impl<T,ti_object>
{
    static ti_object eval(const T& val)
    {
        return val.get_type();
    };
};

template<class T>
struct convert_ti_object_impl
{
};

template<>
struct convert_ti_object_impl<Integer>
{
    static ti_object eval(ti_object t)          { return t; };
    static ti_object eval(ti_empty )            { return predefined::get_ti_int();};
};

template<>
struct convert_ti_object_impl<Float>
{
    static ti_object eval(ti_object t)          { return t; };
    static ti_object eval(ti_empty )            { return predefined::get_ti_float();};
};

template<>
struct convert_ti_object_impl<Real>
{
    static ti_object eval(ti_object t)          { return t; };
    static ti_object eval(ti_empty )            { return predefined::get_ti_real();};
};

template<>
struct convert_ti_object_impl<Complex>
{
    static ti_object eval(ti_object t)          { return t; };
    static ti_object eval(ti_empty )            { return predefined::get_ti_complex();};
};

template<>
struct convert_ti_object_impl<Float_complex>
{
    static ti_object eval(ti_object t)          { return t; };
    static ti_object eval(ti_empty )            { return predefined::get_ti_float_complex();};
};

};};};

namespace matcl { namespace ti
{

std::ostream& ti::operator<<(std::ostream& os, ti_object t)
{
    return operator<<(os,static_cast<const dynamic::Type&>(t));
};

std::istream& ti::operator>>(std::istream& is, ti_object& t)
{
    dynamic::Type tmp;
    std::istream& ret = operator>>(is,tmp);
    t = tmp;
    return ret;
};

//----------------------------------------------------------------------
//                          get_ti
//----------------------------------------------------------------------
template<class T>
typename get_ti_type<T>::type get_ti(const T& m)
{
    using ti_type = typename get_ti_type<T>::type;
    return details::get_ti_type_impl<T,ti_type>::eval(m);
};

template MATCL_MATREP_EXPORT get_ti_type<Integer>::type         get_ti<Integer>(const Integer& m);
template MATCL_MATREP_EXPORT get_ti_type<Float>::type           get_ti<matcl::Float>(const Float& m);
template MATCL_MATREP_EXPORT get_ti_type<Real>::type            get_ti<matcl::Real>(const Real& m);
template MATCL_MATREP_EXPORT get_ti_type<Float_complex>::type   get_ti<matcl::Float_complex>(const Float_complex& m);
template MATCL_MATREP_EXPORT get_ti_type<Complex>::type         get_ti<matcl::Complex>(const Complex& m);
template MATCL_MATREP_EXPORT get_ti_type<Object>::type          get_ti<Object>(const Object& m);

template MATCL_MATREP_EXPORT get_ti_type<raw::Matrix<Integer,struct_dense>>::type
    get_ti<raw::Matrix<Integer,struct_dense>>(const raw::Matrix<Integer,struct_dense>& m);
template MATCL_MATREP_EXPORT get_ti_type<raw::Matrix<Float,struct_dense>>::type
    get_ti<raw::Matrix<Float,struct_dense>>(const raw::Matrix<Float,struct_dense>& m);
template MATCL_MATREP_EXPORT get_ti_type<raw::Matrix<Real,struct_dense>>::type
    get_ti<raw::Matrix<Real,struct_dense>>(const raw::Matrix<Real,struct_dense>& m);
template MATCL_MATREP_EXPORT get_ti_type<raw::Matrix<Float_complex,struct_dense>>::type
    get_ti<raw::Matrix<Float_complex,struct_dense>>(const raw::Matrix<Float_complex,struct_dense>& m);
template MATCL_MATREP_EXPORT get_ti_type<raw::Matrix<Complex,struct_dense>>::type
    get_ti<raw::Matrix<Complex,struct_dense>>(const raw::Matrix<Complex,struct_dense>& m);
template MATCL_MATREP_EXPORT get_ti_type<raw::Matrix<Object,struct_dense>>::type
    get_ti<raw::Matrix<Object,struct_dense>>(const raw::Matrix<Object,struct_dense>& m);

template MATCL_MATREP_EXPORT get_ti_type<raw::Matrix<Integer,struct_banded>>::type
    get_ti<raw::Matrix<Integer,struct_banded>>(const raw::Matrix<Integer,struct_banded>& m);
template MATCL_MATREP_EXPORT get_ti_type<raw::Matrix<Float,struct_banded>>::type
    get_ti<raw::Matrix<Float,struct_banded>>(const raw::Matrix<Float,struct_banded>& m);
template MATCL_MATREP_EXPORT get_ti_type<raw::Matrix<Real,struct_banded>>::type
    get_ti<raw::Matrix<Real,struct_banded>>(const raw::Matrix<Real,struct_banded>& m);
template MATCL_MATREP_EXPORT get_ti_type<raw::Matrix<Float_complex,struct_banded>>::type
    get_ti<raw::Matrix<Float_complex,struct_banded>>(const raw::Matrix<Float_complex,struct_banded>& m);
template MATCL_MATREP_EXPORT get_ti_type<raw::Matrix<Complex,struct_banded>>::type
    get_ti<raw::Matrix<Complex,struct_banded>>(const raw::Matrix<Complex,struct_banded>& m);
template MATCL_MATREP_EXPORT get_ti_type<raw::Matrix<Object,struct_banded>>::type
    get_ti<raw::Matrix<Object,struct_banded>>(const raw::Matrix<Object,struct_banded>& m);

template MATCL_MATREP_EXPORT get_ti_type<raw::Matrix<Integer,struct_sparse>>::type
    get_ti<raw::Matrix<Integer,struct_sparse>>(const raw::Matrix<Integer,struct_sparse>& m);
template MATCL_MATREP_EXPORT get_ti_type<raw::Matrix<Float,struct_sparse>>::type
    get_ti<raw::Matrix<Float,struct_sparse>>(const raw::Matrix<Float,struct_sparse>& m);
template MATCL_MATREP_EXPORT get_ti_type<raw::Matrix<Real,struct_sparse>>::type
    get_ti<raw::Matrix<Real,struct_sparse>>(const raw::Matrix<Real,struct_sparse>& m);
template MATCL_MATREP_EXPORT get_ti_type<raw::Matrix<Float_complex,struct_sparse>>::type
    get_ti<raw::Matrix<Float_complex,struct_sparse>>(const raw::Matrix<Float_complex,struct_sparse>& m);
template MATCL_MATREP_EXPORT get_ti_type<raw::Matrix<Complex,struct_sparse>>::type
    get_ti<raw::Matrix<Complex,struct_sparse>>(const raw::Matrix<Complex,struct_sparse>& m);
template MATCL_MATREP_EXPORT get_ti_type<raw::Matrix<Object,struct_sparse>>::type
    get_ti<raw::Matrix<Object,struct_sparse>>(const raw::Matrix<Object,struct_sparse>& m);

//----------------------------------------------------------------------
//                          convert_ti_object
//----------------------------------------------------------------------
template<class val_type>
ti_object convert_ti_object(ti_empty t)
{
    return details::convert_ti_object_impl<val_type>::eval(t);
};

ti_object ti::ti_object_type(value_code vc)
{
    switch (vc)
    {
        case value_code::v_integer:         return ti_object_type<Integer>();
        case value_code::v_real:            return ti_object_type<Real>();
        case value_code::v_float:           return ti_object_type<Float>();
        case value_code::v_complex:         return ti_object_type<Complex>();
        case value_code::v_float_complex:   return ti_object_type<Float_complex>();

        case value_code::v_object:
            throw error::object_value_type_not_allowed("ti_object_type");

        default:
        {
            matcl_assert(0,"invalid case");
            throw error::error_general("invalid case");
        }
    };
};

template MATCL_MATREP_EXPORT ti_object matcl::ti::convert_ti_object<Integer>(ti_empty t);
template MATCL_MATREP_EXPORT ti_object matcl::ti::convert_ti_object<Float>(ti_empty t);
template MATCL_MATREP_EXPORT ti_object matcl::ti::convert_ti_object<Real>(ti_empty t);
template MATCL_MATREP_EXPORT ti_object matcl::ti::convert_ti_object<Float_complex>(ti_empty t);
template MATCL_MATREP_EXPORT ti_object matcl::ti::convert_ti_object<Complex>(ti_empty t);

//----------------------------------------------------------------------
//                          get_object_ti
//----------------------------------------------------------------------
template<class T>
typename get_ti_type<T>::type get_object_ti()
{
    using ti_type =  typename get_ti_type<T>::type;
    return ti_type();
};

template MATCL_MATREP_EXPORT get_ti_type<Integer>::type         ti::get_object_ti<Integer>();
template MATCL_MATREP_EXPORT get_ti_type<Float>::type           ti::get_object_ti<Float>();
template MATCL_MATREP_EXPORT get_ti_type<Real>::type            ti::get_object_ti<Real>();
template MATCL_MATREP_EXPORT get_ti_type<Float_complex>::type   ti::get_object_ti<Float_complex>();
template MATCL_MATREP_EXPORT get_ti_type<Complex>::type         ti::get_object_ti<Complex>();
template MATCL_MATREP_EXPORT get_ti_type<Object>::type          ti::get_object_ti<Object>();

ti_object details::unify_ti_assign_obj(ti_object t1, ti_object t2)
{
    bool ta = has_trivial_assignment(t1,t2);

    if (ta == false)
        return t1;
    else
        return unify_ti_obj(t1,t2);
};

ti_object details::unify_ti_obj(ti_object t1, ti_object t2)
{
    return dynamic::operations::unify_types(t1,t2);
};

ti_object details::get_return_ti_obj(const dynamic::function_name& fn, ti_object t1, ti_object t2)
{
    dynamic::Type t[] = {t1, t2};
    return dynamic::operations::return_type(fn, 2, t);
};

ti_object details::get_return_ti_obj(const dynamic::function_name& fn, ti_object t1)
{
    dynamic::Type t[] = {t1};
    return dynamic::operations::return_type(fn, 1, t);
};

};};
