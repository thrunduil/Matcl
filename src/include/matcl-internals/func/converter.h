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

namespace mr = matcl::raw;

template<class V, class S> class Matrix;

namespace details
{
    template<class T> struct check_type         {};
    template<> struct check_type<Integer>       { using type = Integer; };
    template<> struct check_type<Float>         { using type = Float; };
    template<> struct check_type<Real>          { using type = Real; };
    template<> struct check_type<Float_complex> { using type = Float_complex; };
    template<> struct check_type<Complex>       { using type = Complex; };
    template<> struct check_type<Object>        { using type = Object; };

    template<class V, class S> 
    struct check_type<Matrix<V,S>>              { using type = V; };

    template<class T> struct promote_type
    {
        using type  = typename matcl::details::lazy_select_if
                    <
                        matcl::details::is_scalar<T>::value,
                        matcl::details::promote_scalar<T>,
                        matcl::details::lazy_type<T>
                    >::type;
    };

    template<class ret,class T>
    struct converter_scalars
    {
        static ret eval(const T& val,ti::ti_type<ret>)
        {
            return val;
        };

        static ret eval(const T& val)
        {
            return val;
        };

        static ti::ti_type<ret> convert_ti(const T& )
        {
            return ti::ti_type<ret>();
        };
    };

    template<class T>
    struct converter_scalars<Object,T>
    {
        static Object eval(const T& val,ti::ti_type<Object> ti)
        {
            return Object(ti,val);
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

    template<>
    struct converter_scalars<Integer,Float>
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
    struct converter_scalars<Integer,Real>
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
    struct converter_scalars<Float,Integer>
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
    struct converter_scalars<Float,Real>
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
    struct converter_scalars<Float,Complex>
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
    struct converter_scalars<Float,Float_complex>
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
    struct converter_scalars<Real,Complex>
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
    struct converter_scalars<Real,Float_complex>
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
    struct converter_scalars<Integer,Complex>
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
    struct converter_scalars<Integer,Float_complex>
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
    struct converter_scalars<Float_complex,Integer>
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
    struct converter_scalars<Float_complex,Real>
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
    struct converter_scalars<Float_complex,Complex>
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
    struct converter_scalars<Complex,Float_complex>
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
    struct converter_scalars_test : public converter_scalars<ret,T>
    {};

    template<class ret>
    struct converter_scalars_test<ret,ret>
    {
        static const ret& eval(const ret& val,ti::ti_type<ret>)
        {
            return val;
        };

        static ret eval(const ret&& val,ti::ti_type<ret>)
        {
            return val;
        };

        static const ret& eval(const ret& val)
        {
            return val;
        };

        static ret eval(const ret&& val)
        {
            return val;
        };
    };

    template<>
    struct converter_scalars_test<Object,Object>
    {
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

    template<class T> struct create_object_ti{};
    template<> struct create_object_ti<Integer>  
    {
        static ti::ti_object eval() {return ti::predefined::get_ti_int();};
    };

    template<> struct create_object_ti<Float>     
    {
        static ti::ti_object eval() {return ti::predefined::get_ti_float();};
    };

    template<> struct create_object_ti<Real>     
    {
        static ti::ti_object eval() {return ti::predefined::get_ti_real();};
    };

    template<> struct create_object_ti<Float_complex>  
    {
        static ti::ti_object eval() {return ti::predefined::get_ti_float_complex();};
    };

    template<> struct create_object_ti<Complex>  
    {
        static ti::ti_object eval() {return ti::predefined::get_ti_complex();};
    };

    template<class V, class S> struct create_object_ti<raw::Matrix<V,S>> : create_object_ti<V>{};

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

    template<class ret, class T>
    struct MATCL_MATREP_EXPORT converter_impl
    {
        using tinfo_ret = ti::ti_type<typename check_type<ret>::type>;

        static ret          eval(const T& val, const tinfo_ret& ti);
        static ret          eval(const T& val)      { return eval(val,convert_ti(val));};
        static tinfo_ret convert_ti(const T& val)   { return convert_ti_impl<T,tinfo_ret>::eval(val); };
    };

    template<class ret, class T>
    struct converter_impl_test : public converter_impl<ret,T>
    {};

    template<class T>
    struct converter_impl_test<T,T>
    {
        using tinfo_ret = ti::ti_type<typename check_type<T>::type>;

        static T            eval(const T& val,ti::ti_object ti) { return converter_impl<T,T>::eval(val,ti);};
        static const T&		eval(const T& val,ti::ti_empty)     { return val; }
        static T			eval(const T&& val,ti::ti_empty)    { return T(std::move(val)); }
        static const T&		eval(const T& val)                  { return val; };
        static T		    eval(T&& val)                       { return T(std::move(val)); };

        static tinfo_ret convert_ti(const T& val)
        {
            return convert_ti_impl<T,tinfo_ret>::eval(val);
        };
    };

    template<class ret,class T>
    struct converter_base
    {
        using type  = typename matcl::details::select_if
                    <
                        matcl::details::is_scalar<T>::value && matcl::details::is_scalar<ret>::value
                        && !std::is_same<T,Object>::value && !std::is_same<ret,Object>::value,
                        converter_scalars_test<ret,T>,
                        converter_impl_test<ret,T>
                    >::type ;
    };

    template<class ret,class T, bool with_ti>
    struct add_ref_if_allowed
    {
        using type = ret;
    };

    template<class ret>
    struct add_ref_if_allowed<ret,ret,false>
    {
        using type = const ret&;
    };

    template<class ret>
    struct add_ref_if_allowed<ret,ret,true>
    {
        using value_type = typename details::check_type<ret>::type;

        using type  = typename matcl::details::select_if
                    <
                        std::is_same<value_type,Object>::value == true,
                        ret,
                        const ret&
                    >::type;
    };
};

template<class ret, class T>
struct converter : details::converter_base<ret,T>::type
{
    public:
        using base_type     = typename details::converter_base<ret,T>::type;
        using dum_type_1    = typename details::check_type<ret>::type;
        using dum_type_2    = typename details::check_type<T>::type;
        using ti_type       = ti::ti_type<dum_type_1>;
        using ret_ref       = typename details::add_ref_if_allowed<ret,T,false>::type;
        using ret_ref_ti    = typename details::add_ref_if_allowed<ret,T,true>::type;

        static ret_ref_ti eval(const T& val,ti_type ti)
        {
            return base_type::eval(val,ti);
        };

        static ret_ref eval(const T& val)
        {
            return base_type::eval(val);
        };

        static ret eval(const T&& val,ti_type ti)
        {
            return base_type::eval(std::move(val),ti);
        };

        static ret eval(const T&& val)
        {
            return base_type::eval(std::move(val));
        };
};

};};
