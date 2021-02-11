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

#include "matcl-core/matrix/scalar_types.h"
#include "matcl-matrep/details/mpl.h"
#include "matcl-matrep/details/isa.h"
#include "matcl-core/error/exception_message.h"
#include "matcl-core/matrix/complex_type.h"
#include "matcl-matrep/objects/details/type_info_object.h"
#include "matcl-scalar/object.h"
#include "matcl-internals/container/mat_d.h"
#include "matcl-internals/container/mat_s.h"
#include "matcl-internals/container/mat_b.h"

namespace matcl { namespace raw
{

template<class V, class S> 
class Matrix;

template<class ret, class T>
struct converter_scalar;

namespace details
{
    //----------------------------------------------------------------------------
    //                 helpers
    //----------------------------------------------------------------------------

    template<class T> struct check_type                 {};
    template<> struct check_type<Integer>               { using type = Integer; };
    template<> struct check_type<Float>                 { using type = Float; };
    template<> struct check_type<Real>                  { using type = Real; };
    template<> struct check_type<Float_complex>         { using type = Float_complex; };
    template<> struct check_type<Complex>               { using type = Complex; };
    template<> struct check_type<Object>                { using type = Object; };

    template<class V, class S> 
    struct check_type<Matrix<V,S>>                      { using type = V; };

    template<class T> struct check_scalar_type          {};
    template<> struct check_scalar_type<Integer>        { using type = Integer; };
    template<> struct check_scalar_type<Float>          { using type = Float; };
    template<> struct check_scalar_type<Real>           { using type = Real; };
    template<> struct check_scalar_type<Float_complex>  { using type = Float_complex; };
    template<> struct check_scalar_type<Complex>        { using type = Complex; };
    template<> struct check_scalar_type<Object>         { using type = Object; };

    //----------------------------------------------------------------------------
    //                 create_object_ti
    //----------------------------------------------------------------------------

    template<class T> struct create_object_ti{};
    
    template<> struct create_object_ti<Integer>  
    {
        static ti::ti_object eval() { return ti::predefined::get_ti_int();};
    };

    template<> struct create_object_ti<Float>     
    {
        static ti::ti_object eval() { return ti::predefined::get_ti_float();};
    };

    template<> struct create_object_ti<Real>     
    {
        static ti::ti_object eval() { return ti::predefined::get_ti_real();};
    };

    template<> struct create_object_ti<Float_complex>  
    {
        static ti::ti_object eval() { return ti::predefined::get_ti_float_complex();};
    };

    template<> struct create_object_ti<Complex>  
    {
        static ti::ti_object eval() { return ti::predefined::get_ti_complex();};
    };

    template<class V, class S> 
    struct create_object_ti<raw::Matrix<V,S>> 
        : create_object_ti<V>
    {};

    //----------------------------------------------------------------------------
    //                 convert_ti_impl
    //----------------------------------------------------------------------------
    template<class T, class ti_t>
    struct convert_ti_impl
    {
        static ti_t eval(const T&)
        {
            return ti_t();
        };
    };

    template<class T>
    struct convert_ti_impl<T, ti::ti_object>
    {
        static ti::ti_object eval(const T&)
        {
            return create_object_ti<T>::eval();
        };
    };

    template<>
    struct convert_ti_impl<Object, ti::ti_object>
    {
        static ti::ti_object eval(const Object& val)
        {
            return val.get_type();
        };
    };

    template<class S>
    struct convert_ti_impl<raw::Matrix<Object,S>, ti::ti_object>
    {
        static ti::ti_object eval(const raw::Matrix<Object,S>& val)
        {
            return val.get_type();
        };
    };

    //----------------------------------------------------------------------------
    //                 converter_mat_mat_impl
    //----------------------------------------------------------------------------
    template<class ret, class T>
    struct MATCL_MATREP_EXPORT converter_mat_mat_impl
    {
        using dum1          = typename check_type<ret>::type;
        using dum2          = typename check_type<T>::type;

        static const bool is_scal_ret   = md::is_scalar<ret>::value;
        static const bool is_scal_in    = md::is_scalar<T>::value;

        static_assert(is_scal_ret == false && is_scal_in == false, "");

        using tinfo_ret     = ti::ti_type<dum1>;
        using ret_ref_ti    = ret;
        using ret_ref       = ret;
        using ret_ref_holder= const ret&;
        using ret_holder    = const ret&;

        static ret          eval(const T& val, const tinfo_ret& ti);
        
        static ret eval(const T& val)
        { 
            return eval(val, convert_ti(val));
        };

        static ret_ref_holder eval(const T& val, matcl::Matrix& tmp_holder)
        {
            tmp_holder      = matcl::Matrix(eval(val), false);
            return tmp_holder.get_impl<ret>();
        };

        static ret_holder eval(T&& val, matcl::Matrix& tmp_holder)
        {
            tmp_holder      = matcl::Matrix(eval(std::move(val)), false);
            return tmp_holder.get_impl<ret>();
        };

        static ret_ref_holder eval(const T& val, const tinfo_ret& ti, matcl::Matrix& tmp_holder)
        {
            tmp_holder      = matcl::Matrix(eval(val, ti), false);
            return tmp_holder.get_impl<ret>();
        };

        static ret_holder eval(T&& val, const tinfo_ret& ti, matcl::Matrix& tmp_holder)
        {
            tmp_holder      = matcl::Matrix(eval(std::move(val), ti), false);
            return tmp_holder.get_impl<ret>();
        };

        static tinfo_ret convert_ti(const T& val)
        { 
            return convert_ti_impl<T,tinfo_ret>::eval(val); 
        };
    };

    template<class ret, class T>
    struct MATCL_MATREP_EXPORT converter_scal_mat_impl
    {
        using dum1          = typename check_type<ret>::type;
        using dum2          = typename check_type<T>::type;

        static const bool is_scal_ret   = md::is_scalar<ret>::value;
        static const bool is_scal_in    = md::is_scalar<T>::value;

        static_assert(is_scal_ret == false && is_scal_in == true, "");

        using tinfo_ret     = ti::ti_type<dum1>;
        using ret_ref_ti    = ret;
        using ret_ref       = ret;
        using ret_ref_holder= const ret&;
        using ret_holder    = const ret&;

        static ret          eval(const T& val, const tinfo_ret& ti);

        static ret eval(const T& val)
        { 
            return eval(val,convert_ti(val));
        };

        static ret_ref_holder eval(const T& val, matcl::Matrix& tmp_holder)
        {
            tmp_holder      = matcl::Matrix(eval(val), false);
            return tmp_holder.get_impl<ret>();
        };

        static ret_holder eval(T&& val, matcl::Matrix& tmp_holder)
        {
            tmp_holder      = matcl::Matrix(eval(std::move(val)), false);
            return tmp_holder.get_impl<ret>();
        };

        static ret_ref_holder eval(const T& val, const tinfo_ret& ti, matcl::Matrix& tmp_holder)
        {
            tmp_holder      = matcl::Matrix(eval(val, ti), false);
            return tmp_holder.get_impl<ret>();
        };

        static ret_holder eval(T&& val, const tinfo_ret& ti, matcl::Matrix& tmp_holder)
        {
            tmp_holder      = matcl::Matrix(eval(std::move(val), ti), false);
            return tmp_holder.get_impl<ret>();
        };

        static tinfo_ret convert_ti(const T& val)
        { 
            return convert_ti_impl<T,tinfo_ret>::eval(val); 
        };
    };

    template<class ret, class T>
    struct MATCL_MATREP_EXPORT converter_mat_scal_impl
    {
        static_assert(md::dependent_false<ret>::value, "conversion of matrix to scalar is not allowed");
    };

    template<class ret>
    struct MATCL_MATREP_EXPORT converter_scal_obj_impl
    {
        using dum1          = typename check_type<ret>::type;        
        using T             = Object;

        static const bool is_scal_ret   = md::is_scalar<ret>::value;

        static_assert(is_scal_ret == true, "");

        using tinfo_ret     = ti::ti_type<dum1>;
        using ret_ref_ti    = ret;
        using ret_ref       = ret;
        using ret_ref_holder= ret;
        using ret_holder    = ret;

        static ret          eval(const T& val, const tinfo_ret& ti);

        static ret eval(const T& val)
        { 
            return eval(val, convert_ti(val));
        };

        static ret_ref_holder eval(const T& val, matcl::Matrix&)
        {
            return eval(val);
        };

        static ret_holder eval(T&& val, matcl::Matrix&)
        {
            return eval(std::move(val));
        };

        static ret_ref_holder eval(const T& val, const tinfo_ret& ti, matcl::Matrix&)
        {
            return eval(val, ti);
        };

        static ret_holder eval(T&& val, const tinfo_ret& ti, matcl::Matrix&)
        {
            return eval(std::move(val), ti);
        };

        static tinfo_ret convert_ti(const T& val)
        { 
            return convert_ti_impl<T,tinfo_ret>::eval(val); 
        };
    };

    //----------------------------------------------------------------------------
    //                 converter_mat_mat_test
    //----------------------------------------------------------------------------
    template<class ret, class T>
    struct converter_mat_mat_test : public converter_mat_mat_impl<ret,T>
    {};

    template<class T>
    struct converter_mat_mat_test<T, T>
    {
        using tinfo_ret = ti::ti_type<typename check_type<T>::type>;

        using ret_ref_ti    = T;
        using ret_ref       = const T&;
        using ret_ref_holder= const T&;
        using ret_holder    = const T&;

        static T            eval(const T& val,ti::ti_object ti) { return converter_mat_mat_impl<T,T>::eval(val,ti);};
        static const T&		eval(const T& val,ti::ti_empty)     { return val; }
        static T			eval(T&& val,ti::ti_empty)          { return T(std::move(val)); }
        static const T&		eval(const T& val)                  { return val; };
        static T		    eval(T&& val)                       { return T(std::move(val)); };

        static ret_ref_holder eval(const T& val, matcl::Matrix& tmp_holder)
        {
            return val;
        };

        static ret_holder eval(T&& val, matcl::Matrix& tmp_holder)
        {
            tmp_holder      = matcl::Matrix(std::move(val), false);
            return tmp_holder.get_impl<ret>();
        };

        static ret_ref_holder eval(const T& val,ti::ti_object ti, matcl::Matrix& tmp_holder)
        {
            tmp_holder  = matcl::Matrix(eval(val, ti), false);
            return tmp_holder.get_impl<T>();
        };

        static ret_holder eval(T&& val,ti::ti_object ti, matcl::Matrix& tmp_holder)
        {
            tmp_holder  = matcl::Matrix(eval(std::move(val), ti), false);
            return tmp_holder.get_impl<T>();
        };

        static ret_ref_holder eval(const T& val,ti::ti_empty, matcl::Matrix&)
        {
            return val;
        };

        static ret_holder eval(T&& val,ti::ti_empty ti, matcl::Matrix& tmp_holder)
        {
            tmp_holder      = matcl::Matrix(std::move(val), false);
            return tmp_holder.get_impl<ret>();
        };
    };

    //----------------------------------------------------------------------------
    //                 converter_scalars
    //----------------------------------------------------------------------------
    template<class ret, class T, bool Is_ret_obj = std::is_same<ret, Object>::value, 
            bool Is_in_obj = std::is_same<T, Object>::value>
    struct converter_scalars
    {};

    template<class T>
    struct converter_scalars<Object, T, true, false>
    {
        static Object eval(const T& val,ti::ti_type<Object> ti)
        {
            return matcl::convert_to_object(ti, val);
        };

        static Object eval(const T& val)
        {
            return Object(val);
        };

        static ti::ti_type<Object> convert_ti(const T& )
        {
            return ti::ti_type<Object>();
        };
    };

    template<class T>
    struct converter_scalars<T, Object, false, true> : converter_scal_obj_impl<T>
    {};

    template<class ret, class T>
    struct converter_scalars<ret, T, false, false>
    {
        static ret eval(const T& val, ti::ti_type<ret>)
        {
            return ret(val);
        };

        static ret eval(const T& val)
        {
            return ret(val);
        };

        static ti::ti_type<ret> convert_ti(const T& )
        {
            return ti::ti_type<ret>();
        };
    }; 

    template<>
    struct converter_scalars<Integer,Float, false, false>
    {
        static Integer eval(const Float& val,ti::ti_empty)
        {
            return eval(val);
        };

        static Integer eval(const Float& val)
        {
            if (val - (Integer)val != 0)
                matcl::error::get_global_messanger()->warning_precision_lost_float_to_int(val);                    

            return (Integer)val;
        };

        static ti::ti_int convert_ti(const Float& )
        {
            return ti::ti_int();
        };
    };

    template<>
    struct converter_scalars<Integer,Real, false, false>
    {
        static Integer eval(const Real& val,ti::ti_empty)
        {
            return eval(val);
        };

        static Integer eval(const Real& val)
        {
            if (val - (Integer)val != 0)
                matcl::error::get_global_messanger()->warning_precision_lost_real_to_int(val);                    

            return (Integer)val;
        };

        static ti::ti_int convert_ti(const Real& )
        {
            return ti::ti_int();
        };
    };

    template<>
    struct converter_scalars<Float,Integer, false, false>
    {
        static Float eval(const Integer& val,ti::ti_empty)
        {
            return eval(val);
        };

        static Float eval(const Integer& val)
        {
            if (val - (Float)val != 0)
                matcl::error::get_global_messanger()->warning_precision_lost_int_to_float(val);                    

            return (Float)val;
        };

        static ti::ti_float convert_ti(const Integer& )
        {
            return ti::ti_float();
        };
    };

    template<>
    struct converter_scalars<Float,Real, false, false>
    {
        static Float eval(const Real& val,ti::ti_empty)
        {
            return eval(val);
        };

        static Float eval(const Real& val)
        {
            if (val - (Float)val != 0)
                matcl::error::get_global_messanger()->warning_precision_lost_real_to_float(val);                    

            return (Float)val;
        };

        static ti::ti_float convert_ti(const Real& )
        {
            return ti::ti_float();
        };
    };

    template<>
    struct converter_scalars<Float,Complex, false, false>
    {
        static Float eval(const Complex& val,ti::ti_empty)
        {
            return eval(val);
        };

        static Float eval(const Complex& val)
        {
            if (imag(val) != 0)
                matcl::error::get_global_messanger()->warning_precision_lost_compl_to_real(val);            

            return converter_scalars<Float,Real>::eval(real(val));
        };

        static ti::ti_float convert_ti(const Complex& )
        {
            return ti::ti_float();
        };
    };

    template<>
    struct converter_scalars<Float,Float_complex, false, false>
    {
        static Float eval(const Float_complex& val,ti::ti_empty)
        {
            return eval(val);
        };

        static Float eval(const Float_complex& val)
        {
            if (imag(val) != 0)
                matcl::error::get_global_messanger()->warning_precision_lost_float_compl_to_float(val);            

            return converter_scalars<Float,Real>::eval(real(val));
        };

        static ti::ti_float convert_ti(const Float_complex& )
        {
            return ti::ti_float();
        };
    };

    template<>
    struct converter_scalars<Real,Complex, false, false>
    {
        static Real eval(const Complex& val,ti::ti_empty)
        {
            return eval(val);
        };

        static Real eval(const Complex& val)
        {
            if (imag(val) != 0)
                matcl::error::get_global_messanger()->warning_precision_lost_compl_to_real(val);            

            return real(val);
        };

        static ti::ti_real convert_ti(const Complex& )
        {
            return ti::ti_real();
        };
    };

    template<>
    struct converter_scalars<Real,Float_complex, false, false>
    {
        static Real eval(const Float_complex& val,ti::ti_empty)
        {
            return eval(val);
        };

        static Real eval(const Float_complex& val)
        {
            if (imag(val) != 0)
                matcl::error::get_global_messanger()->warning_precision_lost_float_compl_to_float(val);            

            return real(val);
        };

        static ti::ti_real convert_ti(const Float_complex& )
        {
            return ti::ti_real();
        };
    };
    
    template<>
    struct converter_scalars<Integer,Complex, false, false>
    {
        static Integer eval(const Complex& val,ti::ti_empty)
        {
            return eval(val);
        };

        static Integer eval(const Complex& val)
        {
            Real tmp = converter_scalars<Real,Complex>::eval(val);
            return converter_scalars<Integer,Real>::eval(tmp);
        };

        static ti::ti_int convert_ti(const Complex& )
        {
            return ti::ti_int();
        };
    };

    template<>
    struct converter_scalars<Integer,Float_complex, false, false>
    {
        static Integer eval(const Float_complex& val,ti::ti_empty)
        {
            return eval(val);
        };

        static Integer eval(const Float_complex& val)
        {
            Float tmp = converter_scalars<Float,Float_complex>::eval(val);
            return converter_scalars<Integer,Float>::eval(tmp);
        };

        static ti::ti_int convert_ti(const Float_complex& )
        {
            return ti::ti_int();
        };
    };

    template<>
    struct converter_scalars<Float_complex,Integer, false, false>
    {
        static Float_complex eval(const Integer& val,ti::ti_empty)
        {
            return eval(val);
        };

        static Float_complex eval(const Integer& val)
        {
            if (val - (Float)val != 0)
            {
                matcl::error::get_global_messanger()->warning_precision_lost_int_to_float(val);                    
            };
            return (Float)val;
        };

        static ti::ti_float_compl convert_ti(const Integer& )
        {
            return ti::ti_float_compl();
        };
    };

    template<>
    struct converter_scalars<Float_complex,Real, false, false>
    {
        static Float_complex eval(const Real& val,ti::ti_empty)
        {
            return eval(val);
        };

        static Float_complex eval(const Real& val)
        {
            if (val - (Float)val != 0)
            {
                matcl::error::get_global_messanger()->warning_precision_lost_real_to_float(val);                    
            };
            return (Float)val;
        };

        static ti::ti_float_compl convert_ti(const Real& )
        {
            return ti::ti_float_compl();
        };
    };

    template<>
    struct converter_scalars<Float_complex,Complex, false, false>
    {
        static Float_complex eval(const Complex& val,ti::ti_empty)
        {
            return eval(val);
        };

        static Float_complex eval(const Complex& val)
        {
            if (val.value.real() - (Float)val.value.real() != 0)
                matcl::error::get_global_messanger()->warning_precision_lost_real_to_float(val.value.real());

            if (val.value.imag() - (Float)val.value.imag() != 0)
                matcl::error::get_global_messanger()->warning_precision_lost_real_to_float(val.value.imag());

            return Float_complex((Float)val.value.real(),(Float)val.value.imag());
        };

        static ti::ti_float_compl convert_ti(const Complex& )
        {
            return ti::ti_float_compl();
        };
    };

    template<>
    struct converter_scalars<Complex,Float_complex, false, false>
    {
        static Complex eval(const Float_complex& val,ti::ti_empty)
        {
            return eval(val);
        };

        static Complex eval(const Float_complex& val)
        {
            return Complex(val.value.real(),val.value.imag());
        };

        static ti::ti_compl convert_ti(const Float_complex& )
        {
            return ti::ti_compl();
        };
    };

    template<class ret, class T>
    struct converter_scalars_test : public converter_scalars<ret, T>
    {
        using ret_ref_ti    = ret;
        using ret_ref       = ret;
    };

    template<class ret>
    struct converter_scalars_test<ret, ret>
    {
        using ret_ref_ti    = const ret&;
        using ret_ref       = const ret&;

        static const ret& eval(const ret& val,ti::ti_type<ret>)
        {
            return val;
        };

        static ret eval(ret&& val,ti::ti_type<ret>)
        {
            return val;
        };

        static const ret& eval(const ret& val)
        {
            return val;
        };

        static ret eval(ret&& val)
        {
            return val;
        };
    };

    template<>
    struct converter_scalars_test<Object,Object>
    {
        using ret_ref_ti    = Object;
        using ret_ref       = const Object&;

        static Object eval(const Object& val,ti::ti_type<Object> ti)
        {
            return Object(ti,val);
        };

        static const Object& eval(const Object& val)
        {
            return val;
        };

        static Object eval(Object&& val)
        {
            return Object(std::move(val));
        };
    };

    //----------------------------------------------------------------------------
    //                 converter_selector
    //----------------------------------------------------------------------------
    template<class ret, class T>
    struct converter_selector
    {
        static const bool is_scal_ret   = md::is_scalar<ret>::value;
        static const bool is_scal_in    = md::is_scalar<T>::value;

        using type  = typename md::select_if
                    <
                        is_scal_in == true,
                        typename md::select_if
                            <
                                is_scal_ret == true,
                                converter_scalar<ret, T>,
                                converter_scal_mat_impl<ret,T>
                            >::type,
                        typename md::select_if
                            <
                                is_scal_ret == true,
                                converter_mat_scal_impl<ret,T>,
                                converter_mat_mat_test<ret,T>
                            >::type
                    >::type;
    };
};

//----------------------------------------------------------------------------
//                 converter
//----------------------------------------------------------------------------

template<class ret, class T>
struct converter
{
    private:
        using impl_type     = typename details::converter_selector<ret, T>::type;

        using dum_type_1    = typename details::check_type<ret>::type;
        using dum_type_2    = typename details::check_type<T>::type;

    public:
        using ti_type       = ti::ti_type<dum_type_1>;

        using ret_ref_ti    = typename impl_type::ret_ref_ti;
        using ret_ref       = typename impl_type::ret_ref;
        using ret_ref_holder= typename impl_type::ret_ref_holder;
        using ret_holder    = typename impl_type::ret_holder;

        static ret_ref eval(const T& val)
        {
            return impl_type::eval(val);
        };

        static ret eval(T&& val)
        {
            return impl_type::eval(std::move(val));
        };

        static ret_ref_holder eval(const T& val, matcl::Matrix& tmp_holder)
        {
            return impl_type::eval(val, tmp_holder);
        };

        static ret_holder eval(T&& val, matcl::Matrix& tmp_holder)
        {
            return impl_type::eval(std::move(val), tmp_holder);
        };

        static ret_ref_holder eval(const T& val, ti_type ti, matcl::Matrix& tmp_holder)
        {
            return impl_type::eval(val, ti, tmp_holder);
        };

        static ret_holder eval(T&& val, ti_type ti, matcl::Matrix&)
        {
            return impl_type::eval(std::move(val), ti, tmp_holder);
        };
};

//----------------------------------------------------------------------------
//                 converter_scalar
//----------------------------------------------------------------------------

template<class ret, class T>
struct converter_scalar
{
    private:
        using impl_type     = details::converter_scalars_test<ret, T>;

        using dum_type_1    = typename details::check_scalar_type<ret>::type;
        using dum_type_2    = typename details::check_scalar_type<T>::type;

    public:
        using ti_type       = ti::ti_type<dum_type_1>;

        using ret_ref_ti    = typename impl_type::ret_ref_ti;
        using ret_ref       = typename impl_type::ret_ref;
        using ret_ref_holder= ret_ref;
        using ret_holder    = ret;
        
        static ret_ref_ti eval(const T& val, ti_type ti)
        {
            return impl_type::eval(val,ti);
        };

        static ret eval(T&& val, ti_type ti)
        {
            return impl_type::eval(std::move(val),ti);
        };

        static ret_ref eval(const T& val)
        {
            return impl_type::eval(val);
        };

        static ret eval(T&& val)
        {
            return impl_type::eval(std::move(val));
        };

        static ret_ref_holder eval(const T& val, matcl::Matrix&)
        {
            return impl_type::eval(val);
        };

        static ret_holder eval(T&& val, matcl::Matrix&)
        {
            return impl_type::eval(std::move(val));
        };

        static ret_ref_holder eval(const T& val, ti_type ti, matcl::Matrix&)
        {
            return impl_type::eval(val,ti);
        };

        static ret_holder eval(T&& val, ti_type ti, matcl::Matrix&)
        {
            return impl_type::eval(std::move(val),ti);
        };
};

};};
