/*
 *  This file is a part of Matrix Computation Library (MATCL)
 *
 *  Copyright (c) Pawe³ Kowal 2017
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

#include "matcl-core/details/scalfunc_real.h"
#include "matcl-core/details/scalfunc_complex.h"

#include "matcl-scalar/details/utils.h"
#include "matcl-dynamic/object_type.h"
#include "matcl-dynamic/predefined_functions.h"
#include "matcl-dynamic/matcl_function_names.h"
#include "matcl-scalar/objects/object_functions.h"

namespace matcl { namespace raw { namespace details
{

namespace md    = matcl::details;
namespace mdyf  = matcl::dynamic::functions;
namespace mrd   = matcl::raw::details;

template<>
struct is_zero_helper<dynamic::object>
{
    force_inline
    static bool eval(const dynamic::object& val) 
    { 
        return val.is_zero(); 
    };
};

template<>
struct is_one_helper<dynamic::object>
{
    force_inline
    static bool eval(const dynamic::object& val) 
    { 
        return val.is_one(); 
    };
};

template<>
struct is_zero_helper<bool>
{
    force_inline
    static bool eval(bool val)
    { 
        return val == false; 
    };
};

template<>
struct is_zero_helper<std::string>
{
    force_inline
    static bool eval(const std::string& val)
    { 
        return val.empty() == true; 
    };
};

template<class T> force_inline
void assign_helper(T& lhs, const T& rhs)
{ 
    lhs = rhs; 
};

template<> force_inline
void assign_helper(dynamic::object& lhs, const dynamic::object& rhs)
{ 
    lhs = rhs; 
};

template<class T> force_inline
void reset_helper(T& lhs, const T& rhs)
{ 
    lhs = rhs; 
};

template<> force_inline
void reset_helper(dynamic::object& lhs, const dynamic::object& rhs)
{ 
    lhs.reset(rhs); 
};

template<class T> force_inline
void assign_change_helper(T& lhs, const T&, const T& rnew)
{ 
    lhs = rnew; 
};   
template<> force_inline
void assign_change_helper(dynamic::object& lhs, const dynamic::object& old, 
                          const dynamic::object& rnew)
{ 
    lhs = old; 
    lhs = rnew;
};

template<class T>
struct clone_helper
{
    force_inline
    static T eval(const T& v)							
    {	
        return v;	
    };
};
template<>
struct clone_helper<dynamic::object>
{
    static dynamic::object eval(const dynamic::object& v)							
    {	
        return v.clone();	
    };
};

template<class T>
struct imag_helper
{
    force_inline
    static T eval(const T&)							
    {	
        return T();	
    };

    static dynamic::function_name name()
    { 
        return function_name::imag::eval(); 
    };
};
template<>
struct imag_helper<Complex>
{
    force_inline
    static Real eval(const Complex& arg)		
    {	
        return matcl::imag(arg);	
    };
    static dynamic::function_name name() 
    { 
        return mdyf::imag::eval(); 
    };
};
template<>
struct imag_helper<Float_complex>
{
    force_inline
    static Float eval(const Float_complex& arg)		
    {	
        return matcl::imag(arg);	
    };
    static dynamic::function_name name()
    { 
        return mdyf::imag::eval(); 
    };
};
template<>
struct imag_helper<dynamic::object>
{
    static dynamic::object eval(const dynamic::object& arg)		
    {	
        return dynamic::imag(arg);	
    };
    static dynamic::function_name name()
    { 
        return mdyf::imag::eval(); 
    };
};

template<class T>
struct real_helper	
{
    force_inline
    static T eval(const T& arg)		
    {	
        return arg;	
    };
    static dynamic::function_name name() 
    { 
        return function_name::real::eval(); 
    };
};
template<>
struct real_helper<Complex>
{
    force_inline
    static Real eval(const Complex& arg)		
    {	
        return matcl::real(arg);	
    };
    static dynamic::function_name name() 
    { 
        return mdyf::real::eval(); 
    };
};
template<>
struct real_helper<Float_complex>
{
    force_inline
    static Float eval(const Float_complex& arg)		
    {	
        return matcl::real(arg);	
    };
    static dynamic::function_name name() 
    { 
        return mdyf::real::eval(); 
    };
};
template<>
struct real_helper<dynamic::object>
{
    static dynamic::object eval(const dynamic::object& arg)		
    {	
        return dynamic::real(arg);
    };
    static dynamic::function_name name() 
    { 
        return mdyf::real::eval(); 
    };
};

template<class T>
struct cast_bool_helper
{
    force_inline
    static bool eval(const T& val)		
    { 
        return (val == T())? false : true; 
    };
    static dynamic::function_name name() 
    { 
        return mdyf::op_bool::eval(); 
    };
};
template<>
struct cast_bool_helper<Complex>
{
    force_inline
    static bool eval(const Complex& val)	
    {
        return cast_bool_helper<Real>::eval(matcl::real(val)) 
                    || cast_bool_helper<Real>::eval(matcl::imag(val)); 
    };
    static dynamic::function_name name()
    { 
        return mdyf::op_bool::eval(); 
    };
};
template<>
struct cast_bool_helper<Float_complex>
{
    force_inline
    static bool eval(const Float_complex& val)	
    {
        return cast_bool_helper<Float>::eval(matcl::real(val)) 
                    || cast_bool_helper<Float>::eval(matcl::imag(val)); 
    };
    static dynamic::function_name name()
    { 
        return mdyf::op_bool::eval(); 
    };
};
template<>
struct cast_bool_helper<dynamic::object>
{
    force_inline
    static bool eval(const dynamic::object& val)	
    { 
        return dynamic::cast_bool(val);
    };
    static dynamic::function_name name() 
    { 
        return mdyf::op_bool::eval(); 
    };
};

template<class T>
struct op_not_helper
{
    force_inline
    static bool eval(const T& val)	
    { 
        return (val == T())? true : false; 
    };
    static dynamic::function_name name()
    { 
        return mdyf::op_not::eval(); 
    };
};
template<>
struct op_not_helper<Complex>
{
    force_inline
    static bool eval(const Complex& val)		
    {
        return op_not_helper<Real>::eval(matcl::real(val)) 
                        && op_not_helper<Real>::eval(matcl::imag(val)); 
    };
    static dynamic::function_name name()
    { 
        return mdyf::op_not::eval(); 
    };
};
template<>
struct op_not_helper<Float_complex>
{
    force_inline
    static bool eval(const Float_complex& val)		
    {
        return op_not_helper<Float>::eval(matcl::real(val)) 
                        && op_not_helper<Float>::eval(matcl::imag(val)); 
    };
    static dynamic::function_name name()
    { 
        return mdyf::op_not::eval(); 
    };
};
template<>
struct op_not_helper<dynamic::object>
{
    force_inline
    static bool eval(const dynamic::object& val)		  
    {
        return dynamic::operator!(val);
    };
    static dynamic::function_name name()
    { 
        return mdyf::op_not::eval(); 
    };
};

template<class T>
struct op_neg_helper
{
    force_inline
    static bool eval(const T& val)
    { 
        return op_not_helper<T>::eval(val);
    };
    static dynamic::function_name name() { return mdyf::op_neg::eval(); };
};
template<>
struct op_neg_helper<dynamic::object>
{
    static dynamic::object eval(const dynamic::object& val)
    {
        return dynamic::op_neg(val); 
    };
    static dynamic::function_name name() { return mdyf::op_neg::eval(); };
};

template<class T>
struct op_true_helper
{
    force_inline
    static bool eval(const T& val)
    { 
        return cast_bool_helper<T>::eval(val);
    };
    static dynamic::function_name name() { return mdyf::op_true::eval(); };
};
template<>
struct op_true_helper<dynamic::object>
{
    static dynamic::object eval(const dynamic::object& val)
    {
        return dynamic::op_true(val); 
    };
    static dynamic::function_name name() { return mdyf::op_true::eval(); };
};

template<class T>
struct isnan_helper
{
    force_inline
    static bool eval(const T& arg)	
    { 
        return scal_func::isnan(arg); 
    };
    static dynamic::function_name name()
    { 
        return mdyf::is_nan::eval(); 
    };
};
template<>
struct isnan_helper<Integer>					
{
    force_inline
    static bool eval(Integer )	
    { 
        return false; 
    };
    static dynamic::function_name name()
    { 
        return mdyf::is_nan::eval(); 
    };
};
template<>
struct isnan_helper<dynamic::object>
{
    static dynamic::object eval(const dynamic::object& arg)
    {	
        return dynamic::is_nan(arg);
    };
    static dynamic::function_name name()
    { 
        return mdyf::is_nan::eval(); 
    };
};

template<class T>
struct isinf_helper
{
    force_inline
    static bool eval(const T& arg)	
    { 
        return scal_func::isinf(arg);		
    };
    static dynamic::function_name name()
    { 
        return mdyf::is_inf::eval(); 
    };
};
template<>
struct isinf_helper<Integer>					
{
    force_inline
    static bool eval(Integer)	
    { 
        return false; 
    };
    static dynamic::function_name name()
    { 
        return mdyf::is_inf::eval(); 
    };
};
template<>
struct isinf_helper<dynamic::object>
{
    static dynamic::object eval(const dynamic::object& arg)
    {	
        return dynamic::is_inf(arg);
    };
    static dynamic::function_name name()
    { 
        return mdyf::is_inf::eval(); 
    };
};

template<class T>
struct isfinite_helper
{
    force_inline
    static bool eval(const T& arg)	
    { 
        return scal_func::finite(arg); 
    };
    static dynamic::function_name name()
    { 
        return mdyf::is_finite::eval(); 
    };
};
template<>
struct isfinite_helper<Integer>		
{
    force_inline
    static bool eval(Integer)	
    { 
        return true; 
    };
    static dynamic::function_name name()
    { 
        return mdyf::is_finite::eval(); 
    };
};
template<>
struct isfinite_helper<dynamic::object>
{
    static dynamic::object eval(const dynamic::object& arg)
    {	
        return dynamic::is_finite(arg);
    };
    static dynamic::function_name name()
    { 
        return mdyf::is_finite::eval(); 
    };
};

template<class T>
struct isregular_helper
{
    force_inline
    static bool eval(const T& arg)	
    { 
        return scal_func::isregular(arg); 
    };
    static dynamic::function_name name()
    { 
        return mdyf::is_regular::eval(); 
    };
};
template<>
struct isregular_helper<Integer>		
{
    force_inline
    static bool eval(Integer v)	
    { 
        return v != 0; 
    };
    static dynamic::function_name name()
    { 
        return mdyf::is_regular::eval(); 
    };
};
template<>
struct isregular_helper<dynamic::object>
{
    static dynamic::object eval(const dynamic::object& arg)
    {	
        return dynamic::is_regular(arg);
    };
    static dynamic::function_name name()
    { 
        return mdyf::is_regular::eval(); 
    };
};

template<class T>
struct isnormal_helper
{
    using ret_type  = bool;

    force_inline
    static ret_type eval(const T& val)
    { 
        return scal_func::isnormal(val); 
    };

    static dynamic::function_name name() { return mdyf::is_normal::eval(); };
};
template<>
struct isnormal_helper<dynamic::object>
{
    static dynamic::object eval(const dynamic::object& val)
    { 
        return dynamic::is_normal(val); 
    };

    static dynamic::function_name name() { return mdyf::is_normal::eval(); };
};

template<class T>
struct isint_helper
{
    force_inline
    static bool eval(const T& arg)	
    { 
        return scal_func::is_int(arg); 
    };
    static dynamic::function_name name()
    { 
        return mdyf::is_int::eval(); 
    };
};
template<>
struct isint_helper<Integer>		
{
    force_inline
    static bool eval(Integer)	
    { 
        return true; 
    };
    static dynamic::function_name name()
    { 
        return mdyf::is_int::eval(); 
    };
};
template<>
struct isint_helper<dynamic::object>
{
    static dynamic::object eval(const dynamic::object& arg)
    {	
        return dynamic::is_int(arg);
    };
    static dynamic::function_name name()
    { 
        return mdyf::is_int::eval(); 
    };
};

template<class T>
struct isreal_helper
{
    force_inline
    static bool eval(const T& arg)	
    { 
        return scal_func::is_real(arg); 
    };
    static dynamic::function_name name()
    { 
        return mdyf::is_real::eval(); 
    };
};
template<>
struct isreal_helper<Integer>		
{
    force_inline
    static bool eval(Integer)	
    { 
        return true; 
    };
    static dynamic::function_name name()
    { 
        return mdyf::is_real::eval(); 
    };
};
template<>
struct isreal_helper<dynamic::object>
{
    static dynamic::object eval(const dynamic::object& arg)
    {	
        return dynamic::is_real(arg);
    };
    static dynamic::function_name name()
    { 
        return mdyf::is_real::eval(); 
    };
};

template<class T>
struct iszero_helper
{
    force_inline
    static bool eval(const T& arg)	
    { 
        return mrd::is_zero(arg); 
    };
    static dynamic::function_name name()
    { 
        return mdyf::is_zero::eval(); 
    };
};

template<class T>
struct isone_helper
{
    force_inline
    static bool eval(const T& arg)	
    { 
        return mrd::is_one(arg); 
    };
    static dynamic::function_name name()
    { 
        return mdyf::is_one::eval(); 
    };
};

template<class T>
struct conj_helper
{
    force_inline
    static T eval(const T& arg)						
    {	
        return arg;	
    };
    static dynamic::function_name name()
    { 
        return mdyf::conj::eval(); 
    };
};
template<>
struct conj_helper<Complex>
{
    force_inline
    static Complex eval(const Complex& arg)		
    {	
        return scal_func::conj(arg);
    };
    static dynamic::function_name name()
    { 
        return mdyf::conj::eval(); 
    };
};
template<>
struct conj_helper<Float_complex>
{
    force_inline
    static Float_complex eval(const Float_complex& arg)		
    {	
        return scal_func::conj(arg);
    };
    static dynamic::function_name name()
    { 
        return mdyf::conj::eval(); 
    };
};
template<>
struct conj_helper<dynamic::object>
{
    static dynamic::object eval(const dynamic::object& arg)		
    {	
        return dynamic::conj(arg);
    };
    static dynamic::function_name name()
    { 
        return mdyf::conj::eval(); 
    };
};

template<class T>
struct abs_helper
{
    force_inline
    static T eval(const T& arg)
    {	
        return scal_func::abs(arg); 
    };
    static dynamic::function_name name()
    { 
        return mdyf::abs::eval(); 
    };
};
template<>
struct abs_helper<Complex>
{
    force_inline
    static Real eval(const Complex& arg)
    {  
        return scal_func::abs(arg); 
    };
    static dynamic::function_name name()
    { 
        return mdyf::abs::eval(); 
    };
};
template<>
struct abs_helper<Float_complex>
{
    force_inline
    static Float eval(const Float_complex& arg)
    {  
        return scal_func::abs(arg); 
    };
    static dynamic::function_name name() 
    { 
        return mdyf::abs::eval(); 
    };
};
template<>
struct abs_helper<Object>
{
    static dynamic::object eval(const Object& arg)
    {	
        return dynamic::abs(arg);
    };
    static dynamic::function_name name()
    { 
        return mdyf::abs::eval(); 
    };
};

template<class T>
struct eps_helper
{
    using real_type = typename md::real_type_int_real<T>::type;

    force_inline
    static real_type eval(const T& arg)	
    {	
        return scal_func::eps(arg);		
    };
    static dynamic::function_name name() { return mdyf::eps::eval(); };
};
template<>
struct eps_helper<dynamic::object>
{
    static dynamic::object eval(const dynamic::object& arg)
    {	
        return dynamic::eps(arg);
    };
    static dynamic::function_name name() { return mdyf::eps::eval(); };
};

template<class T>
struct sqrt_helper
{
    force_inline
    static Real eval(Real arg)		
    {  
        return scal_func::sqrt(arg);	
    };
    static dynamic::function_name name() { return mdyf::sqrt::eval(); };
};
template<>
struct sqrt_helper<Float>
{
    force_inline
    static Float eval(Float arg)		
    {  
        return scal_func::sqrt(arg);	
    };
    static dynamic::function_name name() { return mdyf::sqrt::eval(); };
};
template<>
struct sqrt_helper<Complex>
{
    force_inline
    static Complex eval(const Complex& arg)
    {	
        return scal_func::sqrt(arg);
    };
    static dynamic::function_name name() { return mdyf::sqrt::eval(); };
};
template<>
struct sqrt_helper<Float_complex>
{
    force_inline
    static Float_complex eval(const Float_complex& arg)
    {	
        return scal_func::sqrt(arg);
    };
    static dynamic::function_name name() { return mdyf::sqrt::eval(); };
};
template<>
struct sqrt_helper<dynamic::object>
{
    static dynamic::object eval(const dynamic::object& arg)
    {	
        return dynamic::sqrt(arg);
    };
    static dynamic::function_name name() { return mdyf::sqrt::eval(); };
};

template<class T>
struct sqrt_c_helper
{
    force_inline
    static Complex eval(Real arg)		
    {  
        return scal_func::sqrt_c(arg);	
    };
    static dynamic::function_name name() { return mdyf::sqrt_c::eval(); };
};
template<>
struct sqrt_c_helper<Float>
{
    force_inline
    static Float_complex eval(Float arg)		
    {  
        return scal_func::sqrt_c(arg);	
    };
    static dynamic::function_name name() { return mdyf::sqrt_c::eval(); };
};
template<>
struct sqrt_c_helper<Complex>
{
    force_inline
    static Complex eval(const Complex& arg)
    {	
        return scal_func::sqrt(arg);
    };
    static dynamic::function_name name() { return mdyf::sqrt_c::eval(); };
};
template<>
struct sqrt_c_helper<Float_complex>
{
    force_inline
    static Float_complex eval(const Float_complex& arg)
    {	
        return scal_func::sqrt(arg);
    };
    static dynamic::function_name name() { return mdyf::sqrt_c::eval(); };
};
template<>
struct sqrt_c_helper<dynamic::object>
{
    static dynamic::object eval(const dynamic::object& arg)
    {	
        return dynamic::sqrt_c(arg);
    };
    static dynamic::function_name name() { return mdyf::sqrt_c::eval(); };
};

template<class T>
struct sqrt1pm1_helper
{
    using ret_type = typename md::unify_types<T,Float>::type;

    force_inline
    static ret_type eval(const T& arg)		
    {  
        return scal_func::sqrt1pm1((ret_type)arg);	
    };
    static dynamic::function_name name() { return mdyf::sqrt1pm1::eval(); };
};
template<>
struct sqrt1pm1_helper<dynamic::object>
{
    using ret_type = dynamic::object;

    static ret_type eval(const dynamic::object& arg)		
    {  
        return dynamic::sqrt1pm1(arg);
    };
    static dynamic::function_name name() { return mdyf::sqrt1pm1::eval(); };
};

template<class T>
struct sqrt1pm1_c_helper
{
    force_inline
    static Complex eval(Real arg)		
    {  
        return scal_func::sqrt1pm1_c(arg);	
    };
    static dynamic::function_name name() { return mdyf::sqrt1pm1_c::eval(); };
};
template<>
struct sqrt1pm1_c_helper<Float>
{
    force_inline
    static Float_complex eval(Float arg)		
    {  
        return scal_func::sqrt1pm1_c(arg);	
    };
    static dynamic::function_name name() { return mdyf::sqrt1pm1_c::eval(); };
};
template<>
struct sqrt1pm1_c_helper<Complex>
{
    force_inline
    static Complex eval(const Complex& arg)		
    {  
        return scal_func::sqrt1pm1(arg);	
    };
    static dynamic::function_name name() { return mdyf::sqrt1pm1_c::eval(); };
};
template<>
struct sqrt1pm1_c_helper<Float_complex>
{
    force_inline
    static Float_complex eval(const Float_complex& arg)		
    {  
        return scal_func::sqrt1pm1(arg);	
    };
    static dynamic::function_name name() { return mdyf::sqrt1pm1_c::eval(); };
};
template<>
struct sqrt1pm1_c_helper<dynamic::object>
{
    using ret_type = dynamic::object;

    static ret_type eval(const dynamic::object& arg)		
    {  
        return dynamic::sqrt1pm1_c(arg);
    };
    static dynamic::function_name name() { return mdyf::sqrt1pm1_c::eval(); };
};

template<class T>
struct arg_helper
{
    using real_type = typename md::real_type_int_real<T>::type;

    force_inline
    static real_type eval(const T& arg)	
    {	
        return scal_func::arg(arg);		
    };
    static dynamic::function_name name() { return mdyf::arg::eval(); };
};
template<>
struct arg_helper<dynamic::object>
{
    static dynamic::object eval(const dynamic::object& arg)
    {	
        return dynamic::arg(arg);
    };
    static dynamic::function_name name() { return mdyf::arg::eval(); };
};

template<class T>
struct abs2_helper
{
    using real_type = typename md::real_type<T>::type;

    force_inline
    static real_type eval(const T& arg)
    {	
        return scal_func::abs2(arg); 
    };
    static dynamic::function_name name() { return mdyf::abs2::eval(); };
};
template<>
struct abs2_helper<dynamic::object>
{
    static dynamic::object eval(const dynamic::object& arg)
    {	
        return dynamic::abs2(arg);
    };
    static dynamic::function_name name() { return mdyf::abs2::eval(); };
};

template<class T>
struct exp_helper
{
    force_inline
    static Real eval(Real arg)		
    {	
        return scal_func::exp(arg);	
    };
    static dynamic::function_name name() { return mdyf::exp::eval(); };
};
template<>
struct exp_helper<Float>
{
    force_inline
    static Float eval(Float arg)		
    {	
        return scal_func::exp(arg);	
    };
    static dynamic::function_name name() { return mdyf::exp::eval(); };
};
template<>
struct exp_helper<Complex>
{
    force_inline
    static Complex eval(const Complex& arg)		
    {	
        return scal_func::exp(arg);
    };
    static dynamic::function_name name() { return mdyf::exp::eval(); };
};
template<>
struct exp_helper<Float_complex>
{
    force_inline
    static Float_complex eval(const Float_complex& arg)		
    {	
        return scal_func::exp(arg);
    };
    static dynamic::function_name name() { return mdyf::exp::eval(); };
};
template<>
struct exp_helper<dynamic::object>
{
    static dynamic::object eval(const dynamic::object& arg)	
    {	
        return dynamic::exp(arg);
    };
    static dynamic::function_name name() { return mdyf::exp::eval(); };
};

template<class T>
struct exp2_helper
{
    force_inline
    static Real eval(Real arg)			
    {	
        return scal_func::exp2(arg);	
    };
    static dynamic::function_name name() { return mdyf::exp2::eval(); };
};
template<>
struct exp2_helper<Float>
{
    force_inline
    static Float eval(Float arg)			
    {	
        return scal_func::exp2(arg);	
    };
    static dynamic::function_name name() { return mdyf::exp2::eval(); };
};
template<>
struct exp2_helper<Complex>
{
    force_inline
    static Complex eval(const Complex& arg)
    {	
        return scal_func::exp2(arg);	
    };
    static dynamic::function_name name() { return mdyf::exp2::eval(); };
};
template<>
struct exp2_helper<Float_complex>
{
    force_inline
    static Float_complex eval(const Float_complex& arg)
    {	
        return scal_func::exp2(arg);	
    };
    static dynamic::function_name name() { return mdyf::exp2::eval(); };
};
template<>
struct exp2_helper<dynamic::object>
{
    static dynamic::object eval(const dynamic::object& arg)
    {	
        return dynamic::exp2(arg);	
    };
    static dynamic::function_name name() { return mdyf::exp2::eval(); };
};

//exp10
template<class T>
struct exp10_helper
{
    using ret_type  = typename md::unify_types<T,Float>::type;

    force_inline
    static ret_type eval(const T& val)
    { 
        return scal_func::exp10(val); 
    };

    static dynamic::function_name name() { return mdyf::exp10::eval(); };
};
template<>
struct exp10_helper<dynamic::object>
{
    using T         = dynamic::object;
    using ret_type  = dynamic::object;

    static dynamic::object eval(const T& val)
    { 
        return dynamic::exp10(val); 
    };

    static dynamic::function_name name() { return mdyf::exp10::eval(); };
};
template<>
struct exp10_helper<Complex>
{
    using T         = Complex;
    using ret_type  = Complex;

    force_inline
    static ret_type eval(const T& val)
    { 
        return scal_func::exp10(val);
    };

    static dynamic::function_name name() { return mdyf::exp10::eval(); };
};
template<>
struct exp10_helper<Float_complex>
{
    using T         = Float_complex;
    using ret_type  = Float_complex;

    force_inline
    static ret_type eval(const T& val)
    { 
        return scal_func::exp10(val);
    };

    static dynamic::function_name name() { return mdyf::exp10::eval(); };
};

template<class T>
struct log_helper
{
    force_inline
    static Real eval(Real arg)
    {	
        return scal_func::log(arg);	
    };
    static dynamic::function_name name() { return mdyf::log::eval(); };
};
template<>
struct log_helper<Float>
{
    force_inline
    static Float eval(Float arg)
    {	
        return scal_func::log(arg);	
    };
    static dynamic::function_name name() { return mdyf::log::eval(); };
};
template<>
struct log_helper<Complex>
{
    force_inline
    static Complex eval(const Complex& arg)	
    {	
        return scal_func::log(arg);
    };
    static dynamic::function_name name() { return mdyf::log::eval(); };
};
template<>
struct log_helper<Float_complex>
{
    force_inline
    static Float_complex eval(const Float_complex& arg)	
    {	
        return scal_func::log(arg);
    };
    static dynamic::function_name name() { return mdyf::log::eval(); };
};
template<>
struct log_helper<dynamic::object>
{
    static dynamic::object eval(const dynamic::object& arg)
    {	
        return dynamic::log(arg);
    };
    static dynamic::function_name name() { return mdyf::log::eval(); };
};

template<class T>
struct log_c_helper
{
    force_inline
    static Complex eval(Real arg)
    {	
        return scal_func::log_c(arg);	
    };
    static dynamic::function_name name() { return mdyf::log_c::eval(); };
};
template<>
struct log_c_helper<Float>
{
    force_inline
    static Float_complex eval(Float arg)
    {	
        return scal_func::log_c(arg);	
    };
    static dynamic::function_name name() { return mdyf::log_c::eval(); };
};
template<>
struct log_c_helper<Complex>
{
    force_inline
    static Complex eval(const Complex& arg)	
    {	
        return scal_func::log(arg);
    };
    static dynamic::function_name name() { return mdyf::log_c::eval(); };
};
template<>
struct log_c_helper<Float_complex>
{
    force_inline
    static Float_complex eval(const Float_complex& arg)	
    {	
        return scal_func::log(arg);
    };
    static dynamic::function_name name() { return mdyf::log_c::eval(); };
};
template<>
struct log_c_helper<dynamic::object>
{
    static dynamic::object eval(const dynamic::object& arg)
    {	
        return dynamic::log_c(arg);
    };
    static dynamic::function_name name() { return mdyf::log_c::eval(); };
};

template<class T>
struct log1p_helper
{
    using ret_type = typename md::unify_types<T,Float>::type;

    force_inline
    static ret_type eval(const T& arg)
    {	
        return scal_func::log1p((ret_type)arg);	
    };
    static dynamic::function_name name() { return mdyf::log1p::eval(); };
};
template<>
struct log1p_helper<dynamic::object>
{
    using ret_type = dynamic::object;

    static ret_type eval(const dynamic::object& arg)
    {	
        return dynamic::log1p(arg);	
    };
    static dynamic::function_name name() { return mdyf::log1p::eval(); };
};

template<class T>
struct log1p_c_helper
{
    force_inline
    static Complex eval(Real arg)
    {	
        return scal_func::log1p_c(arg);	
    };
    static dynamic::function_name name() { return mdyf::log1p_c::eval(); };
};
template<>
struct log1p_c_helper<Float>
{
    force_inline
    static Float_complex eval(Float arg)
    {	
        return scal_func::log1p_c(arg);	
    };
    static dynamic::function_name name() { return mdyf::log1p_c::eval(); };
};
template<>
struct log1p_c_helper<Complex>
{
    force_inline
    static Complex eval(const Complex& arg)
    {	
        return scal_func::log1p(arg);	
    };
    static dynamic::function_name name() { return mdyf::log1p_c::eval(); };
};
template<>
struct log1p_c_helper<Float_complex>
{
    force_inline
    static Float_complex eval(const Float_complex& arg)
    {	
        return scal_func::log1p(arg);	
    };
    static dynamic::function_name name() { return mdyf::log1p_c::eval(); };
};
template<>
struct log1p_c_helper<dynamic::object>
{
    using ret_type = dynamic::object;

    static ret_type eval(const dynamic::object& arg)
    {	
        return dynamic::log1p_c(arg);	
    };
    static dynamic::function_name name() { return mdyf::log1p_c::eval(); };
};

template<class T>
struct log2_helper
{
    force_inline
    static Real eval(Real arg)
    {	
        return scal_func::log2(arg);	
    };
    static dynamic::function_name name() { return mdyf::log2::eval(); };
};
template<>
struct log2_helper<Float>
{
    force_inline
    static Float eval(Float arg)
    {	
        return scal_func::log2(arg);	
    };
    static dynamic::function_name name() { return mdyf::log2::eval(); };
};
template<>
struct log2_helper<Complex>
{
    force_inline
    static Complex eval(const Complex& arg)
    {	
        return scal_func::log2(arg);	
    };
    static dynamic::function_name name() { return mdyf::log2::eval(); };
};
template<>
struct log2_helper<Float_complex>
{
    force_inline
    static Float_complex eval(const Float_complex& arg)
    {	
        return scal_func::log2(arg);	
    };
    static dynamic::function_name name() { return mdyf::log2::eval(); };
};

template<>
struct log2_helper<dynamic::object>
{
    static dynamic::object eval(const dynamic::object& arg)
    {	
        return dynamic::log2(arg);
    };
    static dynamic::function_name name() { return mdyf::log2::eval(); };
};

template<class T>
struct log2_c_helper
{
    force_inline
    static Complex eval(Real arg)
    {	
        return scal_func::log2_c(arg);	
    };
    static dynamic::function_name name() { return mdyf::log2_c::eval(); };
};
template<>
struct log2_c_helper<Float>
{
    force_inline
    static Float_complex eval(Float arg)
    {	
        return scal_func::log2_c(arg);	
    };
    static dynamic::function_name name() { return mdyf::log2_c::eval(); };
};
template<>
struct log2_c_helper<Complex>
{
    force_inline
    static Complex eval(const Complex& arg)
    {	
        return scal_func::log2(arg);	
    };
    static dynamic::function_name name() { return mdyf::log2_c::eval(); };
};
template<>
struct log2_c_helper<Float_complex>
{
    force_inline
    static Float_complex eval(const Float_complex& arg)
    {	
        return scal_func::log2(arg);	
    };
    static dynamic::function_name name() { return mdyf::log2_c::eval(); };
};

template<>
struct log2_c_helper<dynamic::object>
{
    static dynamic::object eval(const dynamic::object& arg)
    {	
        return dynamic::log2_c(arg);
    };
    static dynamic::function_name name() { return mdyf::log2_c::eval(); };
};

template<class T>
struct log10_helper	
{
    force_inline
    static Real eval(Real arg)
    {	
        return scal_func::log10(arg);	
    };
    static dynamic::function_name name() { return mdyf::log10::eval(); };
};
template<>
struct log10_helper<Float>
{
    force_inline
    static Float eval(Float arg)
    {	
        return scal_func::log10(arg);	
    };
    static dynamic::function_name name() { return mdyf::log10::eval(); };
};
template<>
struct log10_helper<Complex>
{
    force_inline
    static Complex eval(const Complex& arg)
    {	
        return scal_func::log10(arg);
    };
    static dynamic::function_name name() { return mdyf::log10::eval(); };
};
template<>
struct log10_helper<Float_complex>
{
    force_inline
    static Float_complex eval(const Float_complex& arg)
    {	
        return scal_func::log10(arg);
    };
    static dynamic::function_name name() { return mdyf::log10::eval(); };
};
template<>
struct log10_helper<dynamic::object>
{
    static dynamic::object eval(const dynamic::object& arg)
    {	
        return dynamic::log10(arg);
    };
    static dynamic::function_name name() { return mdyf::log10::eval(); };
};

template<class T>
struct log10_c_helper	
{
    force_inline
    static Complex eval(Real arg)
    {	
        return scal_func::log10_c(arg);	
    };
    static dynamic::function_name name() { return mdyf::log10_c::eval(); };
};
template<>
struct log10_c_helper<Float>
{
    force_inline
    static Float_complex eval(Float arg)
    {	
        return scal_func::log10_c(arg);	
    };
    static dynamic::function_name name() { return mdyf::log10_c::eval(); };
};
template<>
struct log10_c_helper<Complex>
{
    force_inline
    static Complex eval(const Complex& arg)
    {	
        return scal_func::log10(arg);
    };
    static dynamic::function_name name() { return mdyf::log10_c::eval(); };
};
template<>
struct log10_c_helper<Float_complex>
{
    force_inline
    static Float_complex eval(const Float_complex& arg)
    {	
        return scal_func::log10(arg);
    };
    static dynamic::function_name name() { return mdyf::log10_c::eval(); };
};
template<>
struct log10_c_helper<dynamic::object>
{
    static dynamic::object eval(const dynamic::object& arg)
    {	
        return dynamic::log10_c(arg);
    };
    static dynamic::function_name name() { return mdyf::log10_c::eval(); };
};

template<class T> 
struct sign_helper	 {};

template<> 
struct sign_helper<Integer>
{
    force_inline
    static Integer eval(Integer arg)
    {	
        return scal_func::sign(arg);	
    };
    static dynamic::function_name name() { return mdyf::sign::eval(); };
};
template<> 
struct sign_helper<Real>
{
    force_inline
    static Real eval(Real arg)
    {	
        return scal_func::sign(arg);	
    };
    static dynamic::function_name name() { return mdyf::sign::eval(); };
};
template<> 
struct sign_helper<Float>
{
    force_inline
    static Float eval(Float arg)
    {	
        return scal_func::sign(arg);	
    };
    static dynamic::function_name name() { return mdyf::sign::eval(); };
};
template<> 
struct sign_helper<Complex>
{
    force_inline
    static Complex eval(const Complex& arg)
    {	
        return scal_func::sign(arg);	
    };
    static dynamic::function_name name() { return mdyf::sign::eval(); };
};
template<> 
struct sign_helper<Float_complex>
{
    force_inline
    static Float_complex eval(const Float_complex& arg)
    {	
        return scal_func::sign(arg);	
    };
    static dynamic::function_name name() { return mdyf::sign::eval(); };
};
template<>
struct sign_helper<dynamic::object>
{
    static dynamic::object eval(const dynamic::object& arg)
    {	
        return dynamic::sign(arg);
    };
    static dynamic::function_name name() { return mdyf::sign::eval(); };
};

template<class T> 
struct floor_helper	 {};

template<> 
struct floor_helper<Integer>
{
    force_inline
    static Integer eval(Integer arg)
    {	
        return arg;	
    };
    static dynamic::function_name name() { return mdyf::floor::eval(); };
};
template<> 
struct floor_helper<Real>
{
    force_inline
    static Real eval(Real arg)
    {	
        return scal_func::floor(arg);
    };
    static dynamic::function_name name() { return mdyf::floor::eval(); };
};
template<> 
struct floor_helper<Float>
{
    force_inline
    static Float eval(Float arg)
    {	
        return scal_func::floor(arg);
    };
    static dynamic::function_name name() { return mdyf::floor::eval(); };
};

template<> 
struct floor_helper<Complex>
{
    force_inline
    static Complex eval(const Complex& arg)
    {	
        return scal_func::floor(arg);	
    };
    static dynamic::function_name name() { return mdyf::floor::eval(); };
};
template<> 
struct floor_helper<Float_complex>
{
    force_inline
    static Float_complex eval(const Float_complex& arg)
    {	
        return scal_func::floor(arg);	
    };
    static dynamic::function_name name() { return mdyf::floor::eval(); };
};
template<> 
struct floor_helper<dynamic::object>
{
    static dynamic::object eval(const dynamic::object& arg)
    {	
        return dynamic::floor(arg);	
    };
    static dynamic::function_name name() { return mdyf::floor::eval(); };
};

template<class T>
struct ifloor_helper
{
    force_inline
    static Integer eval(const T& arg)	
    {	
        return scal_func::ifloor(arg);
    };
    static dynamic::function_name name() { return mdyf::ifloor::eval(); };
};
template<>
struct ifloor_helper<Complex>
{
    force_inline
    static Integer eval(const Complex& arg)		
    {	
        return scal_func::ifloor(arg);
    };
    static dynamic::function_name name() { return mdyf::ifloor::eval(); };
};
template<>
struct ifloor_helper<Float_complex>
{
    force_inline
    static Integer eval(const Float_complex& arg)		
    {	
        return scal_func::ifloor(arg);
    };
    static dynamic::function_name name() { return mdyf::ifloor::eval(); };
};
template<>
struct ifloor_helper<dynamic::object>
{
    static dynamic::object eval(const dynamic::object& arg)		
    {	
        return dynamic::ifloor(arg);
    };
    static dynamic::function_name name() { return mdyf::ifloor::eval(); };
};

template<class T> 
struct ceil_helper	 {};

template<> 
struct ceil_helper<Integer>
{
    force_inline
    static Integer eval(Integer arg)
    {	
        return arg;	
    };
    static dynamic::function_name name() { return mdyf::ceil::eval(); };
};
template<> 
struct ceil_helper<Real>
{
    force_inline
    static Real eval(Real arg)
    {	
        return scal_func::ceil(arg);
    };
    static dynamic::function_name name() { return mdyf::ceil::eval(); };
};
template<> 
struct ceil_helper<Float>
{
    force_inline
    static Float eval(Float arg)
    {	
        return scal_func::ceil(arg);
    };
    static dynamic::function_name name() { return mdyf::ceil::eval(); };
};

template<> 
struct ceil_helper<Complex>
{
    force_inline
    static Complex eval(const Complex& arg)
    {	
        return scal_func::ceil(arg);	
    };
    static dynamic::function_name name() { return mdyf::ceil::eval(); };
};
template<> 
struct ceil_helper<Float_complex>
{
    force_inline
    static Float_complex eval(const Float_complex& arg)
    {	
        return scal_func::ceil(arg);	
    };
    static dynamic::function_name name() { return mdyf::ceil::eval(); };
};
template<>
struct ceil_helper<dynamic::object>
{
    static dynamic::object eval(const dynamic::object& arg)
    {	
        return dynamic::ceil(arg);
    };
    static dynamic::function_name name() { return mdyf::ceil::eval(); };
};
template<class T> 
struct round_helper	 {};

template<> 
struct round_helper<Integer>
{
    force_inline
    static Integer eval(Integer arg)
    {	
        return arg;	
    };
    static dynamic::function_name name() { return mdyf::round::eval(); };
};
template<> 
struct round_helper<Real>
{
    force_inline
    static Real eval(Real arg)
    {	
        return scal_func::round(arg);
    };
    static dynamic::function_name name() { return mdyf::round::eval(); };
};
template<> 
struct round_helper<Float>
{
    force_inline
    static Float eval(Float arg)
    {	
        return scal_func::round(arg);	
    };
    static dynamic::function_name name() { return mdyf::round::eval(); };
};

template<> 
struct round_helper<Complex>
{
    force_inline
    static Complex eval(const Complex& arg)
    {	
        return scal_func::round(arg);	
    };
    static dynamic::function_name name() { return mdyf::round::eval(); };
};
template<> 
struct round_helper<Float_complex>
{
    force_inline
    static Float_complex eval(const Float_complex& arg)
    {	
        return scal_func::round(arg);	
    };
    static dynamic::function_name name() { return mdyf::round::eval(); };
};
template<>
struct round_helper<dynamic::object>
{
    static dynamic::object eval(const dynamic::object& arg)	
    {	
        return dynamic::round(arg);
    };
    static dynamic::function_name name() { return mdyf::round::eval(); };
};

template<class T> 
struct trunc_helper	 {};

template<> 
struct trunc_helper<Integer>
{
    force_inline
    static Integer eval(Integer arg)
    {	
        return arg;	
    };
    static dynamic::function_name name() { return mdyf::trunc::eval(); };
};
template<> 
struct trunc_helper<Real>
{
    force_inline
    static Real eval(Real arg)
    {	
        return scal_func::trunc(arg);
    };
    static dynamic::function_name name() { return mdyf::trunc::eval(); };
};
template<> 
struct trunc_helper<Float>
{
    force_inline
    static Float eval(Float arg)
    {	
        return scal_func::trunc(arg);	
    };
    static dynamic::function_name name() { return mdyf::trunc::eval(); };
};

template<> 
struct trunc_helper<Complex>
{
    force_inline
    static Complex eval(const Complex& arg)
    {	
        return scal_func::trunc(arg);	
    };
    static dynamic::function_name name() { return mdyf::trunc::eval(); };
};
template<> 
struct trunc_helper<Float_complex>
{
    force_inline
    static Float_complex eval(const Float_complex& arg)
    {	
        return scal_func::trunc(arg);	
    };
    static dynamic::function_name name() { return mdyf::trunc::eval(); };
};
template<>
struct trunc_helper<dynamic::object>
{
    static dynamic::object eval(const dynamic::object& arg)
    {	
        return dynamic::trunc(arg);
    };
    static dynamic::function_name name() { return mdyf::trunc::eval(); };
};

template<class T>
struct iceil_helper
{
    force_inline
    static Integer eval(const T& arg)
    {	
        return scal_func::iceil(arg);
    };
    static dynamic::function_name name() { return mdyf::iceil::eval(); };
};
template<>
struct iceil_helper<Complex>
{
    force_inline
    static Integer eval(const Complex& arg)	
    {	
        return scal_func::iceil(arg);
    };
    static dynamic::function_name name() { return mdyf::iceil::eval(); };
};
template<>
struct iceil_helper<Float_complex>
{
    force_inline
    static Integer eval(const Float_complex& arg)	
    {	
        return scal_func::iceil(arg);
    };
    static dynamic::function_name name() { return mdyf::iceil::eval(); };
};
template<>
struct iceil_helper<dynamic::object>
{
    static dynamic::object eval(const dynamic::object& arg)
    {	
        return dynamic::iceil(arg);
    };
    static dynamic::function_name name() { return mdyf::iceil::eval(); };
};

template<class T>
struct iround_helper
{
    force_inline
    static Integer eval(const T& arg)		
    {	
        return scal_func::iround(arg);
    };
    static dynamic::function_name name() { return mdyf::iround::eval(); };
};
template<>
struct iround_helper<Complex>
{
    force_inline
    static Integer eval(const Complex& arg)	
    {	
        return scal_func::iround(arg);
    };
    static dynamic::function_name name() { return mdyf::iround::eval(); };
};
template<>
struct iround_helper<Float_complex>
{
    force_inline
    static Integer eval(const Float_complex& arg)	
    {	
        return scal_func::iround(arg);
    };
    static dynamic::function_name name() { return mdyf::iround::eval(); };
};

template<>
struct iround_helper<dynamic::object>
{
    static dynamic::object eval(const dynamic::object& arg)
    {	
        return dynamic::iround(arg);
    };
    static dynamic::function_name name() { return mdyf::iround::eval(); };
};

template<class T>
struct itrunc_helper
{
    force_inline
    static Integer eval(const T& arg)
    {	
        return scal_func::itrunc(arg);
    };
    static dynamic::function_name name() { return mdyf::itrunc::eval(); };
};
template<>
struct itrunc_helper<Complex>
{
    force_inline
    static Integer eval(const Complex& arg)
    {	
        return scal_func::itrunc(arg);
    };
    static dynamic::function_name name() { return mdyf::itrunc::eval(); };
};
template<>
struct itrunc_helper<Float_complex>
{
    force_inline
    static Integer eval(const Float_complex& arg)
    {	
        return scal_func::itrunc(arg);
    };
    static dynamic::function_name name() { return mdyf::itrunc::eval(); };
};

template<>
struct itrunc_helper<dynamic::object>
{
    static dynamic::object eval(const dynamic::object& arg)
    {	
        return dynamic::itrunc(arg);
    };
    static dynamic::function_name name() { return mdyf::itrunc::eval(); };
};

template<class T>
struct isign_helper
{
    force_inline
    static Integer eval(const T& arg)	
    {	
        return scal_func::isign(arg);	
    };
    static dynamic::function_name name() { return mdyf::isign::eval(); };
};
template<>
struct isign_helper<Complex>
{
    force_inline
    static Integer eval(const Complex& arg)
    {	
        return scal_func::isign(arg);
    };
    static dynamic::function_name name() { return mdyf::isign::eval(); };
};
template<>
struct isign_helper<Float_complex>
{
    force_inline
    static Integer eval(const Float_complex& arg)
    {	
        return scal_func::isign(arg);
    };
    static dynamic::function_name name() { return mdyf::isign::eval(); };
};
template<>
struct isign_helper<dynamic::object>
{
    static dynamic::object eval(const dynamic::object& arg)
    {	
        return dynamic::isign(arg);
    };
    static dynamic::function_name name() { return mdyf::isign::eval(); };
};

template<class T>
struct sin_helper
{
    force_inline
    static Real eval(Real arg)
    {	
        return scal_func::sin(arg);	
    };
    static dynamic::function_name name() { return mdyf::sin::eval(); };
};
template<>
struct sin_helper<Float>
{
    force_inline
    static Float eval(Float arg)
    {	
        return scal_func::sin(arg);	
    };
    static dynamic::function_name name() { return mdyf::sin::eval(); };
};
template<>
struct sin_helper<Complex>
{
    force_inline
    static Complex eval(const Complex& arg)		
    {	
        return scal_func::sin(arg);
    };
    static dynamic::function_name name() { return mdyf::sin::eval(); };
};
template<>
struct sin_helper<Float_complex>
{
    force_inline
    static Float_complex eval(const Float_complex& arg)		
    {	
        return scal_func::sin(arg);
    };
    static dynamic::function_name name() { return mdyf::sin::eval(); };
};
template<>
struct sin_helper<dynamic::object>
{
    static dynamic::object eval(const dynamic::object& arg)
    {	
        return dynamic::sin(arg);
    };
    static dynamic::function_name name() { return mdyf::sin::eval(); };
};

template<class T>
struct tan_helper
{
    force_inline
    static Real eval(Real arg)
    {	
        return scal_func::tan(arg);	
    };
    static dynamic::function_name name() { return mdyf::tan::eval(); };
};
template<>
struct tan_helper<Float>
{
    force_inline
    static Float eval(Float arg)
    {	
        return scal_func::tan(arg);	
    };
    static dynamic::function_name name() { return mdyf::tan::eval(); };
};

template<>
struct tan_helper<Complex>
{
    force_inline
    static Complex eval(const Complex& arg)
    {	
        return scal_func::tan(arg);
    };
    static dynamic::function_name name() { return mdyf::tan::eval(); };
};
template<>
struct tan_helper<Float_complex>
{
    force_inline
    static Float_complex eval(const Float_complex& arg)
    {	
        return scal_func::tan(arg);
    };
    static dynamic::function_name name() { return mdyf::tan::eval(); };
};
template<>
struct tan_helper<dynamic::object>
{
    static dynamic::object eval(const dynamic::object& arg)	
    {	
        return dynamic::tan(arg);
    };
    static dynamic::function_name name() { return mdyf::tan::eval(); };
};

template<class T>
struct cot_helper
{
    force_inline
    static Real eval(Real arg)	
    {	
        return scal_func::cot(arg);	
    };
    static dynamic::function_name name() { return mdyf::cot::eval(); };
};
template<>
struct cot_helper<Float>
{
    force_inline
    static Float eval(Float arg)	
    {	
        return scal_func::cot(arg);	
    };
    static dynamic::function_name name() { return mdyf::cot::eval(); };
};

template<>
struct cot_helper<Complex>
{
    force_inline
    static Complex eval(const Complex& arg)
    {	
        return scal_func::cot(arg);	
    };
    static dynamic::function_name name() { return mdyf::cot::eval(); };
};
template<>
struct cot_helper<Float_complex>
{
    force_inline
    static Float_complex eval(const Float_complex& arg)
    {	
        return scal_func::cot(arg);	
    };
    static dynamic::function_name name() { return mdyf::cot::eval(); };
};
template<>
struct cot_helper<dynamic::object>
{
    static dynamic::object eval(const dynamic::object& arg)
    {	
        return dynamic::cot(arg);
    };
    static dynamic::function_name name() { return mdyf::cot::eval(); };
};

template<class T>
struct sec_helper
{
    force_inline
    static Real eval(Real arg)
    {	
        return scal_func::sec(arg);	
    };
    static dynamic::function_name name() { return mdyf::sec::eval(); };
};
template<>
struct sec_helper<Float>
{
    force_inline
    static Float eval(Float arg)
    {	
        return scal_func::sec(arg);	
    };
    static dynamic::function_name name() { return mdyf::sec::eval(); };
};

template<>
struct sec_helper<Complex>
{
    force_inline
    static Complex eval(const Complex& arg)
    {	
        return scal_func::sec(arg);	
    };
    static dynamic::function_name name() { return mdyf::sec::eval(); };
};
template<>
struct sec_helper<Float_complex>
{
    force_inline
    static Float_complex eval(const Float_complex& arg)
    {	
        return scal_func::sec(arg);	
    };
    static dynamic::function_name name() { return mdyf::sec::eval(); };
};
template<>
struct sec_helper<dynamic::object>
{
    static dynamic::object eval(const dynamic::object& arg)
    {	
        return dynamic::sec(arg);
    };
    static dynamic::function_name name() { return mdyf::sec::eval(); };
};

template<class T>
struct csc_helper
{
    force_inline
    static Real eval(Real arg)
    {	
        return scal_func::csc(arg);	
    };
    static dynamic::function_name name() { return mdyf::csc::eval(); };
};
template<>
struct csc_helper<Float>
{
    force_inline
    static Float eval(Float arg)
    {	
        return scal_func::csc(arg);	
    };
    static dynamic::function_name name() { return mdyf::csc::eval(); };
};
template<>
struct csc_helper<Complex>
{
    force_inline
    static Complex eval(const Complex& arg)
    {	
        return scal_func::csc(arg);	
    };
    static dynamic::function_name name() { return mdyf::csc::eval(); };
};
template<>
struct csc_helper<Float_complex>
{
    force_inline
    static Float_complex eval(const Float_complex& arg)
    {	
        return scal_func::csc(arg);	
    };
    static dynamic::function_name name() { return mdyf::csc::eval(); };
};
template<>
struct csc_helper<dynamic::object>
{
    static dynamic::object eval(const dynamic::object& arg)
    {	
        return dynamic::csc(arg);
    };
    static dynamic::function_name name() { return mdyf::csc::eval(); };
};

template<class T>
struct asin_helper
{
    force_inline
    static Real eval(Real arg)
    {	
        return scal_func::asin(arg);	
    };
    static dynamic::function_name name() { return mdyf::asin::eval(); };
};
template<>
struct asin_helper<Float>
{
    force_inline
    static Float eval(Float arg)
    {	
        return scal_func::asin(arg);	
    };
    static dynamic::function_name name() { return mdyf::asin::eval(); };
};
template<>
struct asin_helper<Complex>
{
    force_inline
    static Complex eval(const Complex& arg)
    {	
        return scal_func::asin(arg);	
    };
    static dynamic::function_name name() { return mdyf::asin::eval(); };
};
template<>
struct asin_helper<Float_complex>
{
    force_inline
    static Float_complex eval(const Float_complex& arg)
    {	
        return scal_func::asin(arg);	
    };
    static dynamic::function_name name() { return mdyf::asin::eval(); };
};
template<>
struct asin_helper<dynamic::object>
{
    static dynamic::object eval(const dynamic::object& arg)
    {	
        return dynamic::asin(arg);
    };
    static dynamic::function_name name() { return mdyf::asin::eval(); };
};

template<class T>
struct asin_c_helper
{
    force_inline
    static Complex eval(Real arg)
    {	
        return scal_func::asin_c(arg);	
    };
    static dynamic::function_name name() { return mdyf::asin_c::eval(); };
};
template<>
struct asin_c_helper<Float>
{
    force_inline
    static Float_complex eval(Float arg)
    {	
        return scal_func::asin_c(arg);	
    };
    static dynamic::function_name name() { return mdyf::asin_c::eval(); };
};
template<>
struct asin_c_helper<Complex>
{
    force_inline
    static Complex eval(const Complex& arg)
    {	
        return scal_func::asin(arg);	
    };
    static dynamic::function_name name() { return mdyf::asin_c::eval(); };
};
template<>
struct asin_c_helper<Float_complex>
{
    force_inline
    static Float_complex eval(const Float_complex& arg)
    {	
        return scal_func::asin(arg);	
    };
    static dynamic::function_name name() { return mdyf::asin_c::eval(); };
};
template<>
struct asin_c_helper<dynamic::object>
{
    static dynamic::object eval(const dynamic::object& arg)
    {	
        return dynamic::asin_c(arg);
    };
    static dynamic::function_name name() { return mdyf::asin_c::eval(); };
};

template<class T>
struct acos_helper
{
    force_inline
    static Real eval(Real arg)
    {	
        return scal_func::acos(arg);	
    };
    static dynamic::function_name name() { return mdyf::acos::eval(); };
};
template<>
struct acos_helper<Float>
{
    force_inline
    static Float eval(Float arg)
    {	
        return scal_func::acos(arg);	
    };
    static dynamic::function_name name() { return mdyf::acos::eval(); };
};
template<>
struct acos_helper<Complex>
{
    force_inline
    static Complex eval(const Complex& arg)	
    {	
        return scal_func::acos(arg);	
    };
    static dynamic::function_name name() { return mdyf::acos::eval(); };
};
template<>
struct acos_helper<Float_complex>
{
    force_inline
    static Float_complex eval(const Float_complex& arg)	
    {	
        return scal_func::acos(arg);	
    };
    static dynamic::function_name name() { return mdyf::acos::eval(); };
};
template<>
struct acos_helper<dynamic::object>
{
    static dynamic::object eval(const dynamic::object& arg)
    {	
        return dynamic::acos(arg);
    };
    static dynamic::function_name name() { return mdyf::acos::eval(); };
};

template<class T>
struct acos_c_helper
{
    force_inline
    static Complex eval(Real arg)
    {	
        return scal_func::acos_c(arg);	
    };
    static dynamic::function_name name() { return mdyf::acos_c::eval(); };
};
template<>
struct acos_c_helper<Float>
{
    force_inline
    static Float_complex eval(Float arg)
    {	
        return scal_func::acos_c(arg);	
    };
    static dynamic::function_name name() { return mdyf::acos_c::eval(); };
};
template<>
struct acos_c_helper<Complex>
{
    force_inline
    static Complex eval(const Complex& arg)	
    {	
        return scal_func::acos(arg);	
    };
    static dynamic::function_name name() { return mdyf::acos_c::eval(); };
};
template<>
struct acos_c_helper<Float_complex>
{
    force_inline
    static Float_complex eval(const Float_complex& arg)	
    {	
        return scal_func::acos(arg);	
    };
    static dynamic::function_name name() { return mdyf::acos_c::eval(); };
};
template<>
struct acos_c_helper<dynamic::object>
{
    static dynamic::object eval(const dynamic::object& arg)
    {	
        return dynamic::acos_c(arg);
    };
    static dynamic::function_name name() { return mdyf::acos_c::eval(); };
};

template<class T>
struct atan_helper
{
    force_inline
    static Real eval(Real arg)
    {	
        return scal_func::atan(arg);	
    };
    static dynamic::function_name name() { return mdyf::atan::eval(); };
};
template<>
struct atan_helper<Float>
{
    force_inline
    static Float eval(Float arg)
    {	
        return scal_func::atan(arg);	
    };
    static dynamic::function_name name() { return mdyf::atan::eval(); };
};

template<>
struct atan_helper<Complex>
{
    force_inline
    static Complex eval(const Complex& arg)
    {	
        return scal_func::atan(arg);	
    };
    static dynamic::function_name name() { return mdyf::atan::eval(); };
};
template<>
struct atan_helper<Float_complex>
{
    force_inline
    static Float_complex eval(const Float_complex& arg)
    {	
        return scal_func::atan(arg);	
    };
    static dynamic::function_name name() { return mdyf::atan::eval(); };
};
template<>
struct atan_helper<dynamic::object>
{
    static dynamic::object eval(const dynamic::object& arg)
    {	
        return dynamic::atan(arg);
    };
    static dynamic::function_name name() { return mdyf::atan::eval(); };
};
template<class T>
struct acot_helper
{
    force_inline
    static Real eval(Real arg)
    {	
        return scal_func::acot(arg);	
    };
    static dynamic::function_name name() { return mdyf::acot::eval(); };
};
template<>
struct acot_helper<Float>
{
    force_inline
    static Float eval(Float arg)
    {	
        return scal_func::acot(arg);	
    };
    static dynamic::function_name name() { return mdyf::acot::eval(); };
};

template<>
struct acot_helper<Complex>
{
    force_inline
    static Complex eval(const Complex& arg)	
    {	
        return scal_func::acot(arg);	
    };
    static dynamic::function_name name() { return mdyf::acot::eval(); };
};
template<>
struct acot_helper<Float_complex>
{
    force_inline
    static Float_complex eval(const Float_complex& arg)	
    {	
        return scal_func::acot(arg);	
    };
    static dynamic::function_name name() { return mdyf::acot::eval(); };
};
template<>
struct acot_helper<dynamic::object>
{
    static dynamic::object eval(const dynamic::object& arg)	
    {	
        return dynamic::acot(arg);
    };
    static dynamic::function_name name() { return mdyf::acot::eval(); };
};

template<class T>
struct asec_helper
{
    force_inline
    static Real eval(Real arg)
    {	
        return scal_func::asec(arg);	
    };
    static dynamic::function_name name() { return mdyf::asec::eval(); };
};
template<>
struct asec_helper<Float>
{
    force_inline
    static Float eval(Float arg)
    {	
        return scal_func::asec(arg);	
    };
    static dynamic::function_name name() { return mdyf::asec::eval(); };
};
template<>
struct asec_helper<Complex>
{
    force_inline
    static Complex eval(const Complex& arg)
    {	
        return scal_func::asec(arg);	
    };
    static dynamic::function_name name() { return mdyf::asec::eval(); };
};
template<>
struct asec_helper<Float_complex>
{
    force_inline
    static Float_complex eval(const Float_complex& arg)
    {	
        return scal_func::asec(arg);	
    };
    static dynamic::function_name name() { return mdyf::asec::eval(); };
};
template<>
struct asec_helper<dynamic::object>
{
    static dynamic::object eval(const dynamic::object& arg)
    {	
        return dynamic::asec(arg);
    };
    static dynamic::function_name name() { return mdyf::asec::eval(); };
};

template<class T>
struct asec_c_helper
{
    force_inline
    static Complex eval(Real arg)
    {	
        return scal_func::asec_c(arg);	
    };
    static dynamic::function_name name() { return mdyf::asec_c::eval(); };
};
template<>
struct asec_c_helper<Float>
{
    force_inline
    static Float_complex eval(Float arg)
    {	
        return scal_func::asec_c(arg);	
    };
    static dynamic::function_name name() { return mdyf::asec_c::eval(); };
};
template<>
struct asec_c_helper<Complex>
{
    force_inline
    static Complex eval(const Complex& arg)
    {	
        return scal_func::asec(arg);	
    };
    static dynamic::function_name name() { return mdyf::asec_c::eval(); };
};
template<>
struct asec_c_helper<Float_complex>
{
    force_inline
    static Float_complex eval(const Float_complex& arg)
    {	
        return scal_func::asec(arg);	
    };
    static dynamic::function_name name() { return mdyf::asec_c::eval(); };
};
template<>
struct asec_c_helper<dynamic::object>
{
    static dynamic::object eval(const dynamic::object& arg)
    {	
        return dynamic::asec_c(arg);
    };
    static dynamic::function_name name() { return mdyf::asec_c::eval(); };
};

template<class T>
struct acsc_helper
{
    force_inline
    static Real eval(Real arg)
    {	
        return scal_func::acsc(arg);	
    };
    static dynamic::function_name name() { return mdyf::acsc::eval(); };
};
template<>
struct acsc_helper<Float>
{
    force_inline
    static Float eval(Float arg)
    {	
        return scal_func::acsc(arg);	
    };
    static dynamic::function_name name() { return mdyf::acsc::eval(); };
};
template<>
struct acsc_helper<Complex>
{
    force_inline
    static Complex eval(const Complex& arg)
    {	
        return scal_func::acsc(arg);	
    };
    static dynamic::function_name name() { return mdyf::acsc::eval(); };
};
template<>
struct acsc_helper<Float_complex>
{
    force_inline
    static Float_complex eval(const Float_complex& arg)
    {	
        return scal_func::acsc(arg);	
    };
    static dynamic::function_name name() { return mdyf::acsc::eval(); };
};
template<>
struct acsc_helper<dynamic::object>
{
    static dynamic::object eval(const dynamic::object& arg)
    {	
        return dynamic::acsc(arg);
    };
    static dynamic::function_name name() { return mdyf::acsc::eval(); };
};

template<class T>
struct acsc_c_helper
{
    force_inline
    static Complex eval(Real arg)
    {	
        return scal_func::acsc_c(arg);	
    };
    static dynamic::function_name name() { return mdyf::acsc_c::eval(); };
};
template<>
struct acsc_c_helper<Float>
{
    force_inline
    static Float_complex eval(Float arg)
    {	
        return scal_func::acsc_c(arg);	
    };
    static dynamic::function_name name() { return mdyf::acsc_c::eval(); };
};
template<>
struct acsc_c_helper<Complex>
{
    force_inline
    static Complex eval(const Complex& arg)
    {	
        return scal_func::acsc(arg);	
    };
    static dynamic::function_name name() { return mdyf::acsc_c::eval(); };
};
template<>
struct acsc_c_helper<Float_complex>
{
    force_inline
    static Float_complex eval(const Float_complex& arg)
    {	
        return scal_func::acsc(arg);	
    };
    static dynamic::function_name name() { return mdyf::acsc_c::eval(); };
};
template<>
struct acsc_c_helper<dynamic::object>
{
    static dynamic::object eval(const dynamic::object& arg)
    {	
        return dynamic::acsc_c(arg);
    };
    static dynamic::function_name name() { return mdyf::acsc_c::eval(); };
};

template<class T>
struct sinh_helper
{
    force_inline
    static Real eval(Real arg)
    {	
        return scal_func::sinh(arg);
    };
    static dynamic::function_name name() { return mdyf::sinh::eval(); };
};
template<>
struct sinh_helper<Float>
{
    force_inline
    static Float eval(Float arg)
    {	
        return scal_func::sinh(arg);
    };
    static dynamic::function_name name() { return mdyf::sinh::eval(); };
};

template<>
struct sinh_helper<Complex>
{
    force_inline
    static Complex eval(const Complex& arg)
    {	
        return scal_func::sinh(arg);
    };
    static dynamic::function_name name() { return mdyf::sinh::eval(); };
};
template<>
struct sinh_helper<Float_complex>
{
    force_inline
    static Float_complex eval(const Float_complex& arg)
    {	
        return scal_func::sinh(arg);
    };
    static dynamic::function_name name() { return mdyf::sinh::eval(); };
};
template<>
struct sinh_helper<dynamic::object>
{
    static dynamic::object eval(const dynamic::object& arg)
    {	
        return dynamic::sinh(arg);
    };
    static dynamic::function_name name() { return mdyf::sinh::eval(); };
};
template<class T>
struct cosh_helper
{
    force_inline
    static Real eval(Real arg)		
    {	
        return scal_func::cosh(arg);	
    };
    static dynamic::function_name name() { return mdyf::cosh::eval(); };
};
template<>
struct cosh_helper<Float>
{
    force_inline
    static Float eval(Float arg)		
    {	
        return scal_func::cosh(arg);	
    };
    static dynamic::function_name name() { return mdyf::cosh::eval(); };
};

template<>
struct cosh_helper<Complex>
{
    force_inline
    static Complex eval(const Complex& arg)
    {	
        return scal_func::cosh(arg);
    };
    static dynamic::function_name name() { return mdyf::cosh::eval(); };
};
template<>
struct cosh_helper<Float_complex>
{
    force_inline
    static Float_complex eval(const Float_complex& arg)
    {	
        return scal_func::cosh(arg);
    };
    static dynamic::function_name name() { return mdyf::cosh::eval(); };
};
template<>
struct cosh_helper<dynamic::object>
{
    static dynamic::object eval(const dynamic::object& arg)
    {	
        return dynamic::cosh(arg);
    };
    static dynamic::function_name name() { return mdyf::cosh::eval(); };
};
template<class T>
struct tanh_helper
{
    force_inline
    static Real eval(Real arg)
    {	
        return scal_func::tanh(arg);	
    };
    static dynamic::function_name name() { return mdyf::tanh::eval(); };
};
template<>
struct tanh_helper<Float>
{
    force_inline
    static Float eval(Float arg)
    {	
        return scal_func::tanh(arg);	
    };
    static dynamic::function_name name() { return mdyf::tanh::eval(); };
};

template<>
struct tanh_helper<Complex>
{
    force_inline
    static Complex eval(const Complex& arg)
    {	
        return scal_func::tanh(arg);
    };
    static dynamic::function_name name() { return mdyf::tanh::eval(); };
};
template<>
struct tanh_helper<Float_complex>
{
    force_inline
    static Float_complex eval(const Float_complex& arg)
    {	
        return scal_func::tanh(arg);
    };
    static dynamic::function_name name() { return mdyf::tanh::eval(); };
};
template<>
struct tanh_helper<dynamic::object>
{
    static dynamic::object eval(const dynamic::object& arg)
    {	
        return dynamic::tanh(arg);
    };
    static dynamic::function_name name() { return mdyf::tanh::eval(); };
};
template<class T>
struct coth_helper
{
    force_inline
    static Real eval(Real arg)
    {	
        return scal_func::coth(arg);	
    };
    static dynamic::function_name name() { return mdyf::coth::eval(); };
};
template<>
struct coth_helper<Float>
{
    force_inline
    static Float eval(Float arg)
    {	
        return scal_func::coth(arg);	
    };
    static dynamic::function_name name() { return mdyf::coth::eval(); };
};

template<>
struct coth_helper<Complex>
{
    force_inline
    static Complex eval(const Complex& arg)
    {	
        return scal_func::coth(arg);	
    };
    static dynamic::function_name name() { return mdyf::coth::eval(); };
};
template<>
struct coth_helper<Float_complex>
{
    force_inline
    static Float_complex eval(const Float_complex& arg)
    {	
        return scal_func::coth(arg);	
    };
    static dynamic::function_name name() { return mdyf::coth::eval(); };
};
template<>
struct coth_helper<dynamic::object>
{
    static dynamic::object eval(const dynamic::object& arg)
    {	
        return dynamic::coth(arg);
    };
    static dynamic::function_name name() { return mdyf::coth::eval(); };
};
template<class T>
struct sech_helper
{
    force_inline
    static Real eval(Real arg)	
    {	
        return scal_func::sech(arg);	
    };
    static dynamic::function_name name() { return mdyf::sech::eval(); };
};
template<>
struct sech_helper<Float>
{
    force_inline
    static Float eval(Float arg)	
    {	
        return scal_func::sech(arg);	
    };
    static dynamic::function_name name() { return mdyf::sech::eval(); };
};

template<>
struct sech_helper<Complex>
{
    force_inline
    static Complex eval(const Complex& arg)
    {	
        return scal_func::sech(arg);	
    };
    static dynamic::function_name name() { return mdyf::sech::eval(); };
};
template<>
struct sech_helper<Float_complex>
{
    force_inline
    static Float_complex eval(const Float_complex& arg)
    {	
        return scal_func::sech(arg);	
    };
    static dynamic::function_name name() { return mdyf::sech::eval(); };
};
template<>
struct sech_helper<dynamic::object>
{
    static dynamic::object eval(const dynamic::object& arg)
    {	
        return dynamic::sech(arg);
    };
    static dynamic::function_name name() { return mdyf::sech::eval(); };
};
template<class T>
struct csch_helper
{
    force_inline
    static Real eval(Real arg)	
    {	
        return scal_func::csch(arg);	
    };
    static dynamic::function_name name() { return mdyf::csch::eval(); };
};
template<>
struct csch_helper<Float>
{
    force_inline
    static Float eval(Float arg)	
    {	
        return scal_func::csch(arg);	
    };
    static dynamic::function_name name() { return mdyf::csch::eval(); };
};

template<>
struct csch_helper<Complex>
{
    force_inline
    static Complex eval(const Complex& arg)	
    {	
        return scal_func::csch(arg);	
    };
    static dynamic::function_name name() { return mdyf::csch::eval(); };
};
template<>
struct csch_helper<Float_complex>
{
    force_inline
    static Float_complex eval(const Float_complex& arg)	
    {	
        return scal_func::csch(arg);	
    };
    static dynamic::function_name name() { return mdyf::csch::eval(); };
};
template<>
struct csch_helper<dynamic::object>
{
    static dynamic::object eval(const dynamic::object& arg)	
    {	
        return dynamic::csch(arg);
    };
    static dynamic::function_name name() { return mdyf::csch::eval(); };
};
template<class T>
struct asinh_helper
{
    force_inline
    static Real eval(Real arg)		
    {	
        return scal_func::asinh(arg);	
    };
    static dynamic::function_name name() { return mdyf::asinh::eval(); };
};
template<>
struct asinh_helper<Float>
{
    force_inline
    static Float eval(Float arg)		
    {	
        return scal_func::asinh(arg);	
    };
    static dynamic::function_name name() { return mdyf::asinh::eval(); };
};

template<>
struct asinh_helper<Complex>
{
    force_inline
    static Complex eval(const Complex& arg)	
    {	
        return scal_func::asinh(arg);	
    };
    static dynamic::function_name name() { return mdyf::asinh::eval(); };
};
template<>
struct asinh_helper<Float_complex>
{
    force_inline
    static Float_complex eval(const Float_complex& arg)	
    {	
        return scal_func::asinh(arg);	
    };
    static dynamic::function_name name() { return mdyf::asinh::eval(); };
};
template<>
struct asinh_helper<dynamic::object>
{
    static dynamic::object eval(const dynamic::object& arg)
    {	
        return dynamic::asinh(arg);
    };
    static dynamic::function_name name() { return mdyf::asinh::eval(); };
};

template<class T>
struct acosh_helper
{
    force_inline
    static Real eval(Real arg)	
    {	
        return scal_func::acosh(arg);	
    };
    static dynamic::function_name name() { return mdyf::acosh::eval(); };
};
template<>
struct acosh_helper<Float>
{
    force_inline
    static Float eval(Float arg)	
    {	
        return scal_func::acosh(arg);	
    };
    static dynamic::function_name name() { return mdyf::acosh::eval(); };
};
template<>
struct acosh_helper<Complex>
{
    force_inline
    static Complex eval(const Complex& arg)
    {	
        return scal_func::acosh(arg);	
    };
    static dynamic::function_name name() { return mdyf::acosh::eval(); };
};
template<>
struct acosh_helper<Float_complex>
{
    force_inline
    static Float_complex eval(const Float_complex& arg)
    {	
        return scal_func::acosh(arg);	
    };
    static dynamic::function_name name() { return mdyf::acosh::eval(); };
};
template<>
struct acosh_helper<dynamic::object>
{
    static dynamic::object eval(const dynamic::object& arg)
    {	
        return dynamic::acosh(arg);
    };
    static dynamic::function_name name() { return mdyf::acosh::eval(); };
};

template<class T>
struct acosh_c_helper
{
    force_inline
    static Complex eval(Real arg)	
    {	
        return scal_func::acosh_c(arg);	
    };
    static dynamic::function_name name() { return mdyf::acosh_c::eval(); };
};
template<>
struct acosh_c_helper<Float>
{
    force_inline
    static Float_complex eval(Float arg)	
    {	
        return scal_func::acosh_c(arg);	
    };
    static dynamic::function_name name() { return mdyf::acosh_c::eval(); };
};
template<>
struct acosh_c_helper<Complex>
{
    force_inline
    static Complex eval(const Complex& arg)
    {	
        return scal_func::acosh(arg);	
    };
    static dynamic::function_name name() { return mdyf::acosh_c::eval(); };
};
template<>
struct acosh_c_helper<Float_complex>
{
    force_inline
    static Float_complex eval(const Float_complex& arg)
    {	
        return scal_func::acosh(arg);	
    };
    static dynamic::function_name name() { return mdyf::acosh_c::eval(); };
};
template<>
struct acosh_c_helper<dynamic::object>
{
    static dynamic::object eval(const dynamic::object& arg)
    {	
        return dynamic::acosh_c(arg);
    };
    static dynamic::function_name name() { return mdyf::acosh_c::eval(); };
};

template<class T>
struct atanh_helper
{
    force_inline
    static Real eval(Real arg)
    {	
        return scal_func::atanh(arg);	
    };
    static dynamic::function_name name() { return mdyf::atanh::eval(); };
};
template<>
struct atanh_helper<Float>
{
    force_inline
    static Float eval(Float arg)
    {	
        return scal_func::atanh(arg);	
    };
    static dynamic::function_name name() { return mdyf::atanh::eval(); };
};
template<>
struct atanh_helper<Complex>
{
    force_inline
    static Complex eval(const Complex& arg)
    {	
        return scal_func::atanh(arg);	
    };
    static dynamic::function_name name() { return mdyf::atanh::eval(); };
};
template<>
struct atanh_helper<Float_complex>
{
    force_inline
    static Float_complex eval(const Float_complex& arg)
    {	
        return scal_func::atanh(arg);	
    };
    static dynamic::function_name name() { return mdyf::atanh::eval(); };
};
template<>
struct atanh_helper<dynamic::object>
{
    static dynamic::object eval(const dynamic::object& arg)	
    {	
        return dynamic::atanh(arg);
    };
    static dynamic::function_name name() { return mdyf::atanh::eval(); };
};

template<class T>
struct atanh_c_helper
{
    force_inline
    static Complex eval(Real arg)
    {	
        return scal_func::atanh_c(arg);	
    };
    static dynamic::function_name name() { return mdyf::atanh_c::eval(); };
};
template<>
struct atanh_c_helper<Float>
{
    force_inline
    static Float_complex eval(Float arg)
    {	
        return scal_func::atanh_c(arg);	
    };
    static dynamic::function_name name() { return mdyf::atanh_c::eval(); };
};
template<>
struct atanh_c_helper<Complex>
{
    force_inline
    static Complex eval(const Complex& arg)
    {	
        return scal_func::atanh(arg);	
    };
    static dynamic::function_name name() { return mdyf::atanh_c::eval(); };
};
template<>
struct atanh_c_helper<Float_complex>
{
    force_inline
    static Float_complex eval(const Float_complex& arg)
    {	
        return scal_func::atanh(arg);	
    };
    static dynamic::function_name name() { return mdyf::atanh_c::eval(); };
};
template<>
struct atanh_c_helper<dynamic::object>
{
    static dynamic::object eval(const dynamic::object& arg)	
    {	
        return dynamic::atanh_c(arg);
    };
    static dynamic::function_name name() { return mdyf::atanh_c::eval(); };
};

template<class T>
struct acoth_helper
{
    force_inline
    static Real eval(Real arg)
    {	
        return scal_func::acoth(arg);	
    };
    static dynamic::function_name name() { return mdyf::acoth::eval(); };
};
template<>
struct acoth_helper<Float>
{
    force_inline
    static Float eval(Float arg)
    {	
        return scal_func::acoth(arg);	
    };
    static dynamic::function_name name() { return mdyf::acoth::eval(); };
};
template<>
struct acoth_helper<Complex>
{
    force_inline
    static Complex eval(const Complex& arg)
    {	
        return scal_func::acoth(arg);	
    };
    static dynamic::function_name name() { return mdyf::acoth::eval(); };
};
template<>
struct acoth_helper<Float_complex>
{
    force_inline
    static Float_complex eval(const Float_complex& arg)
    {	
        return scal_func::acoth(arg);	
    };
    static dynamic::function_name name() { return mdyf::acoth::eval(); };
};
template<>
struct acoth_helper<dynamic::object>
{
    static dynamic::object eval(const dynamic::object& arg)
    {	
        return dynamic::acoth(arg);
    };
    static dynamic::function_name name() { return mdyf::acoth::eval(); };
};

template<class T>
struct acoth_c_helper
{
    force_inline
    static Complex eval(Real arg)
    {	
        return scal_func::acoth_c(arg);	
    };
    static dynamic::function_name name() { return mdyf::acoth_c::eval(); };
};
template<>
struct acoth_c_helper<Float>
{
    force_inline
    static Float_complex eval(Float arg)
    {	
        return scal_func::acoth_c(arg);	
    };
    static dynamic::function_name name() { return mdyf::acoth_c::eval(); };
};
template<>
struct acoth_c_helper<Complex>
{
    force_inline
    static Complex eval(const Complex& arg)
    {	
        return scal_func::acoth(arg);	
    };
    static dynamic::function_name name() { return mdyf::acoth_c::eval(); };
};
template<>
struct acoth_c_helper<Float_complex>
{
    force_inline
    static Float_complex eval(const Float_complex& arg)
    {	
        return scal_func::acoth(arg);	
    };
    static dynamic::function_name name() { return mdyf::acoth_c::eval(); };
};
template<>
struct acoth_c_helper<dynamic::object>
{
    static dynamic::object eval(const dynamic::object& arg)
    {	
        return dynamic::acoth_c(arg);
    };
    static dynamic::function_name name() { return mdyf::acoth_c::eval(); };
};

template<class T>
struct asech_helper
{
    force_inline
    static Real eval(Real arg)
    {	
        return scal_func::asech(arg);	
    };
    static dynamic::function_name name() { return mdyf::asech::eval(); };
};
template<>
struct asech_helper<Float>
{
    force_inline
    static Float eval(Float arg)
    {	
        return scal_func::asech(arg);	
    };
    static dynamic::function_name name() { return mdyf::asech::eval(); };
};
template<>
struct asech_helper<Complex>
{
    force_inline
    static Complex eval(const Complex& arg)	
    {	
        return scal_func::asech(arg);	
    };
    static dynamic::function_name name() { return mdyf::asech::eval(); };
};
template<>
struct asech_helper<Float_complex>
{
    force_inline
    static Float_complex eval(const Float_complex& arg)	
    {	
        return scal_func::asech(arg);	
    };
    static dynamic::function_name name() { return mdyf::asech::eval(); };
};
template<>
struct asech_helper<dynamic::object>
{
    static dynamic::object eval(const dynamic::object& arg)
    {	
        return dynamic::asech(arg);
    };
    static dynamic::function_name name() { return mdyf::asech::eval(); };
};

template<class T>
struct asech_c_helper
{
    force_inline
    static Complex eval(Real arg)
    {	
        return scal_func::asech_c(arg);	
    };
    static dynamic::function_name name() { return mdyf::asech_c::eval(); };
};
template<>
struct asech_c_helper<Float>
{
    force_inline
    static Float_complex eval(Float arg)
    {	
        return scal_func::asech_c(arg);	
    };
    static dynamic::function_name name() { return mdyf::asech_c::eval(); };
};
template<>
struct asech_c_helper<Complex>
{
    force_inline
    static Complex eval(const Complex& arg)	
    {	
        return scal_func::asech(arg);	
    };
    static dynamic::function_name name() { return mdyf::asech_c::eval(); };
};
template<>
struct asech_c_helper<Float_complex>
{
    force_inline
    static Float_complex eval(const Float_complex& arg)	
    {	
        return scal_func::asech(arg);	
    };
    static dynamic::function_name name() { return mdyf::asech_c::eval(); };
};
template<>
struct asech_c_helper<dynamic::object>
{
    static dynamic::object eval(const dynamic::object& arg)
    {	
        return dynamic::asech_c(arg);
    };
    static dynamic::function_name name() { return mdyf::asech_c::eval(); };
};

template<class T>
struct acsch_helper
{
    force_inline
    static Real eval(Real arg)
    {	
        return scal_func::acsch(arg);	
    };
    static dynamic::function_name name() { return mdyf::acsch::eval(); };
};
template<>
struct acsch_helper<Float>
{
    force_inline
    static Float eval(Float arg)
    {	
        return scal_func::acsch(arg);	
    };
    static dynamic::function_name name() { return mdyf::acsch::eval(); };
};

template<>
struct acsch_helper<Complex>
{
    force_inline
    static Complex eval(const Complex& arg)
    {	
        return scal_func::acsch(arg);	
    };
    static dynamic::function_name name() { return mdyf::acsch::eval(); };
};
template<>
struct acsch_helper<Float_complex>
{
    force_inline
    static Float_complex eval(const Float_complex& arg)
    {	
        return scal_func::acsch(arg);	
    };
    static dynamic::function_name name() { return mdyf::acsch::eval(); };
};
template<>
struct acsch_helper<dynamic::object>
{
    static dynamic::object eval(const dynamic::object& arg)
    {	
        return dynamic::acsch(arg);
    };
    static dynamic::function_name name() { return mdyf::acsch::eval(); };
};

template<class T>
struct cos_helper
{
    force_inline
    static Real eval(Real arg)	
    {	
        return scal_func::cos(arg);	
    };
    static dynamic::function_name name() { return mdyf::cos::eval(); };
};
template<>
struct cos_helper<Float>
{
    force_inline
    static Float eval(Float arg)	
    {	
        return scal_func::cos(arg);	
    };
    static dynamic::function_name name() { return mdyf::cos::eval(); };
};

template<>
struct cos_helper<Complex>
{
    force_inline
    static Complex eval(const Complex& arg)
    {	
        return scal_func::cos(arg);
    };
    static dynamic::function_name name() { return mdyf::cos::eval(); };
};
template<>
struct cos_helper<Float_complex>
{
    force_inline
    static Float_complex eval(const Float_complex& arg)
    {	
        return scal_func::cos(arg);
    };
    static dynamic::function_name name() { return mdyf::cos::eval(); };
};
template<>
struct cos_helper<dynamic::object>
{
    static dynamic::object eval(const dynamic::object& arg)	
    {	
        return dynamic::cos(arg);
    };
    static dynamic::function_name name() { return mdyf::cos::eval(); };
};

template<class T>
struct uminus_helper
{
    force_inline
    static T eval(const T& val)
    { 
        return -val; 
    };
    static dynamic::function_name name() { return mdyf::op_uminus::eval(); };
};
template<>
struct uminus_helper<Complex>
{
    force_inline
    static Complex eval(const Complex& val)
    { 
        return matcl::details::uminus_c(val); 
    };
    static dynamic::function_name name() { return mdyf::op_uminus::eval(); };
};
template<>
struct uminus_helper<Float_complex>
{
    force_inline
    static Float_complex eval(const Float_complex& val)
    { 
        return matcl::details::uminus_c(val); 
    };
    static dynamic::function_name name() { return mdyf::op_uminus::eval(); };
};
template<>
struct uminus_helper<dynamic::object>
{
    static dynamic::object eval(const dynamic::object& arg)
    { 
        return dynamic::operator-(arg); 
    };
    static dynamic::function_name name() { return mdyf::op_uminus::eval(); };
};

template<class T>
struct invs_helper
{
    force_inline
    static Real eval(Integer val)
    { 
        return 1.0 / Real(val); 
    };
    static dynamic::function_name name() { return mdyf::invs::eval(); };
};
template<>
struct invs_helper<Real>
{
    force_inline
    static Real eval(Real val)
    { 
        return 1.0/val;
    };
    static dynamic::function_name name() { return mdyf::invs::eval(); };
};
template<>
struct invs_helper<Float>
{
    force_inline
    static Float eval(Float val)
    { 
        return 1.0f/val;
    };
    static dynamic::function_name name() { return mdyf::invs::eval(); };
};
template<>
struct invs_helper<Complex>
{
    force_inline
    static Complex eval(const Complex& val)
    { 
        return matcl::details::inv_c(val); 
    };
    static dynamic::function_name name() { return mdyf::invs::eval(); };
};
template<>
struct invs_helper<Float_complex>
{
    force_inline
    static Float_complex eval(const Float_complex& val)
    { 
        return matcl::details::inv_c(val); 
    };
    static dynamic::function_name name() { return mdyf::invs::eval(); };
};
template<>
struct invs_helper<dynamic::object>
{
    static dynamic::object eval(const dynamic::object& arg)
    { 
        return dynamic::invs(arg); 
    };
    static dynamic::function_name name() { return mdyf::invs::eval(); };
};

template<class T>
struct inv_helper
{
    using Ret = typename md::unify_types<T, Float>::type;

    force_inline
    static Ret eval(const T& val)
    { 
        return invs_helper<T>::eval(val);
    };
    static dynamic::function_name name() { return mdyf::invs::eval(); };
};
template<>
struct inv_helper<dynamic::object>
{
    static dynamic::object eval(const dynamic::object& arg)
    { 
        return dynamic::inv(arg); 
    };
    static dynamic::function_name name() { return mdyf::inv::eval(); };
};


template<class T> force_inline
bool is_geq_zero(const T& val)
{ 
    return val >= T(); 
};

template<> force_inline
bool is_geq_zero<dynamic::object>(const dynamic::object& val)
{ 
    return (bool)dynamic::operator>=(val,dynamic::object(val.get_type())); 
};

template<> force_inline
bool is_geq_zero<Complex>(const Complex& val)
{ 
    return matcl::details::geq_c(val,0.); 
};

template<class T> force_inline
bool is_leq_zero(const T& val)
{ 
    return val<= T(); 
};

template<> force_inline
bool is_leq_zero<dynamic::object>(const dynamic::object& val)
{ 
    return (bool)dynamic::operator<=(val, dynamic::object(val.get_type())); 
};
template<> force_inline
bool is_leq_zero<Complex>(const Complex& val)
{ 
    return matcl::details::leq_c(val,0.); 
};

template<class T> force_inline
bool is_gt_zero(const T& val)
{ 
    return val > T(); 
};
template<> force_inline
bool is_gt_zero<dynamic::object>(const dynamic::object& val)
{ 
    return (bool)dynamic::operator>(val, dynamic::object(val.get_type())); 
};
template<> force_inline
bool is_gt_zero<Complex>(const Complex& val)
{ 
    return matcl::details::gt_c(val,0.); 
};

template<class T> force_inline
bool is_lt_zero(const T& val)
{ 
    return val < T(); 
};

template<> force_inline
bool is_lt_zero<dynamic::object>(const dynamic::object& val)
{ 
    return (bool)dynamic::operator<(val,dynamic::object(val.get_type())); 
};
template<> force_inline
bool is_lt_zero<Complex>(const Complex& val)
{ 
    return matcl::details::lt_c(val,0.); 
};

//cbrt
template<class T>
struct cbrt_helper
{
    using ret_type  = typename md::unify_types<T,Float>::type;

    force_inline
    static ret_type eval(const T& val)
    { 
        return scal_func::cbrt(val); 
    };

    static dynamic::function_name name() { return mdyf::cbrt::eval(); };
};
template<>
struct cbrt_helper<Integer>
{
    using ret_type  = Real;

    force_inline
    static ret_type eval(Integer val)
    { 
        return scal_func::cbrt((Real)val); 
    };

    static dynamic::function_name name() { return mdyf::cbrt::eval(); };
};
template<>
struct cbrt_helper<dynamic::object>
{
    using ret_type  = dynamic::object;

    static ret_type eval(const dynamic::object& val)
    { 
        return dynamic::cbrt(val); 
    };

    static dynamic::function_name name() { return mdyf::cbrt::eval(); };
};

//expm1
template<class T>
struct expm1_helper
{
    using ret_type  = typename md::unify_types<T,Float>::type;

    force_inline
    static ret_type eval(const T& val)
    { 
        return scal_func::expm1(val); 
    };

    static dynamic::function_name name() { return mdyf::expm1::eval(); };
};

template<>
struct expm1_helper<Integer>
{
    using ret_type  = Real;

    force_inline
    static ret_type eval(Integer val)
    { 
        return scal_func::expm1(Real(val)); 
    };

    static dynamic::function_name name() { return mdyf::expm1::eval(); };
};

template<>
struct expm1_helper<dynamic::object>
{
    using ret_type  = dynamic::object;

    static ret_type eval(const dynamic::object& val)
    { 
        return dynamic::expm1(val); 
    };

    static dynamic::function_name name() { return mdyf::expm1::eval(); };
};

//expm1
template<class T>
struct expi_helper
{
    using ret_type  = typename md::unify_types<T,Float_complex>::type;

    force_inline
    static ret_type eval(const T& val)
    { 
        return scal_func::expi(val); 
    };

    static dynamic::function_name name() { return mdyf::expi::eval(); };
};

template<>
struct expi_helper<Integer>
{
    using ret_type  = Complex;

    force_inline
    static ret_type eval(Integer val)
    { 
        return scal_func::expi(Real(val)); 
    };

    static dynamic::function_name name() { return mdyf::expi::eval(); };
};
template<>
struct expi_helper<dynamic::object>
{
    using ret_type  = dynamic::object;

    static ret_type eval(const dynamic::object& val)
    { 
        return dynamic::expi(val); 
    };

    static dynamic::function_name name() { return mdyf::expi::eval(); };
};

//signbit
template<class T>
struct signbit_helper
{
    using ret_type  = bool;

    force_inline
    static ret_type eval(const T& val)
    { 
        return scal_func::signbit(val); 
    };

    static dynamic::function_name name() { return mdyf::signbit::eval(); };
};
template<>
struct signbit_helper<dynamic::object>
{
    static dynamic::object eval(const dynamic::object& val)
    { 
        return dynamic::signbit(val); 
    };

    static dynamic::function_name name() { return mdyf::signbit::eval(); };
};

//logb
template<class T>
struct logb_helper
{
    using ret_type  = typename md::unify_types<T,Float>::type;

    force_inline
    static ret_type eval(const T& val)
    { 
        return scal_func::logb(val); 
    };

    static dynamic::function_name name() { return mdyf::logb::eval(); };
};
template<>
struct logb_helper<dynamic::object>
{
    using T         = dynamic::object;
    using ret_type  = T;

    static ret_type eval(const T& val)
    { 
        return dynamic::logb(val); 
    };

    static dynamic::function_name name() { return mdyf::logb::eval(); };
};

template<>
struct logb_helper<Complex>
{
    force_inline
    static Real eval(const Complex& val)
    { 
        Real a  = abs_helper<Complex>::eval(val);
        return logb_helper<Real>::eval(a);
    };

    static dynamic::function_name name() { return mdyf::logb::eval(); };
};
template<>
struct logb_helper<Float_complex>
{
    force_inline
    static Float eval(const Float_complex& val)
    { 
        Float a  = abs_helper<Float_complex>::eval(val);
        return logb_helper<Float>::eval(a);
    };

    static dynamic::function_name name() { return mdyf::logb::eval(); };
};

//ilogb
template<class T>
struct ilogb_helper
{
    using ret_type  = Integer;

    force_inline
    static ret_type eval(const T& val)
    { 
        return scal_func::ilogb(val); 
    };

    static dynamic::function_name name() { return mdyf::ilogb::eval(); };
};
template<>
struct ilogb_helper<dynamic::object>
{
    using T         = dynamic::object;
    using ret_type  = Integer;

    static dynamic::object eval(const T& val)
    { 
        return dynamic::ilogb(val); 
    };

    static dynamic::function_name name() { return mdyf::ilogb::eval(); };
};

template<>
struct ilogb_helper<Complex>
{
    force_inline
    static Integer eval(const Complex& val)
    { 
        Real a = abs_helper<Complex>::eval(val);
        return ilogb_helper<Real>::eval(a);
    };

    static dynamic::function_name name() { return mdyf::ilogb::eval(); };
};
template<>
struct ilogb_helper<Float_complex>
{
    force_inline
    static Integer eval(const Float_complex& val)
    { 
        Float a = abs_helper<Float_complex>::eval(val);
        return ilogb_helper<Float>::eval(a);
    };

    static dynamic::function_name name() { return mdyf::ilogb::eval(); };
};

template<class T>
struct fpclassify_helper
{
    force_inline
    static fp_type eval(const T& val)
    {
        return scal_func::fpclassify(val);
    };
    
    static dynamic::function_name name() 
    { 
        return mdyf::fpclassify::eval(); 
    };
};

template<>
struct fpclassify_helper<Integer>
{
    force_inline
    static fp_type eval(Integer val)
    {
        if (val == 0)
            return fp_type::fp_zero;
        else
            return fp_type::fp_normal;
    };
    static dynamic::function_name name() { return mdyf::fpclassify::eval(); };
};
template<>
struct fpclassify_helper<dynamic::object>
{
    force_inline
    static fp_type eval(const dynamic::object& val)
    {
        return dynamic::fpclassify(val);
    };
    static dynamic::function_name name() { return mdyf::fpclassify::eval(); };
};

template<class T>
struct ldexp_helper
{};

template<>
struct ldexp_helper<Real>
{
    force_inline
    static Real eval(Real x, Integer exp)
    {
        return scal_func::ldexp(x,exp);
    }
    static dynamic::function_name name() { return mdyf::ldexp::eval(); };
};

template<>
struct ldexp_helper<Integer>
{
    force_inline
    static Real eval(Integer x, Integer exp)
    {
        return ldexp_helper<Real>::eval(Real(x),exp);
    }
    static dynamic::function_name name() { return mdyf::ldexp::eval(); };
};
template<>
struct ldexp_helper<Float>
{
    force_inline
    static Float eval(Float x, Integer exp)
    {
        return scal_func::ldexp(x,exp);
    }
    static dynamic::function_name name() { return mdyf::ldexp::eval(); };
};
template<>
struct ldexp_helper<Complex>
{
    force_inline
    static Complex eval(const Complex& x, Integer exp)
    {
        return Complex(ldexp_helper<Real>::eval(real(x),exp), 
                       ldexp_helper<Real>::eval(imag(x),exp));
    }
    static dynamic::function_name name() { return mdyf::ldexp::eval(); };
};
template<>
struct ldexp_helper<Float_complex>
{
    force_inline
    static Float_complex eval(const Float_complex& x, Integer exp)
    {
        return Float_complex(ldexp_helper<Float>::eval(real(x),exp), 
                             ldexp_helper<Float>::eval(imag(x),exp));
    }
    static dynamic::function_name name() { return mdyf::ldexp::eval(); };
};

template<>
struct ldexp_helper<dynamic::object>
{
    static dynamic::object eval(const dynamic::object& x, Integer exp)
    {
        return dynamic::ldexp(x, dynamic::object(dynamic::OInteger(exp)));
    }
    static dynamic::function_name name() { return mdyf::ldexp::eval(); };
};

template<class T>
struct scalbn_helper
{};
template<>
struct scalbn_helper<Real>
{
    force_inline
    static Real eval(Real x, Integer n)
    {
        return scal_func::scalbn(x,n);
    }
    static dynamic::function_name name() { return mdyf::ldexp::eval(); };
};
template<>
struct scalbn_helper<Integer>
{
    force_inline
    static Real eval(Integer x, Integer n)
    {
        return scalbn_helper<Real>::eval(Real(x),n);
    }
    static dynamic::function_name name() { return mdyf::ldexp::eval(); };
};
template<>
struct scalbn_helper<Float>
{
    force_inline
    static Float eval(Float x, Integer n)
    {
        return scal_func::scalbn(x,n);
    }
    static dynamic::function_name name() { return mdyf::ldexp::eval(); };
};
template<>
struct scalbn_helper<Complex>
{
    force_inline
    static Complex eval(const Complex& x, Integer n)
    {
        return ldexp_helper<Complex>::eval(x, n);
    }
    static dynamic::function_name name() { return mdyf::ldexp::eval(); };
};
template<>
struct scalbn_helper<Float_complex>
{
    force_inline
    static Float_complex eval(const Float_complex& x, Integer n)
    {
        return ldexp_helper<Float_complex>::eval(x, n);
    }
    static dynamic::function_name name() { return mdyf::ldexp::eval(); };
};
template<>
struct scalbn_helper<dynamic::object>
{
    static dynamic::object eval(const dynamic::object& x, Integer n)
    {
        return dynamic::scalbn(x, dynamic::object(dynamic::OInteger(n)));
    }
    static dynamic::function_name name() { return mdyf::ldexp::eval(); };
};

template<class T>
struct frexp_helper
{};
template<>
struct frexp_helper<Real>
{
    force_inline
    static Real eval(Real x, Integer& n)
    {
        return scal_func::frexp(x, n);
    }
    static dynamic::function_name name() { return mdyf::frexp::eval(); };
};
template<>
struct frexp_helper<Integer>
{
    force_inline
    static Real eval(Integer x, Integer& n)
    {
        return scal_func::frexp(Real(x), n);
    }
    static dynamic::function_name name() { return mdyf::frexp::eval(); };
};
template<>
struct frexp_helper<Float>
{
    force_inline
    static Float eval(Float x, Integer& n)
    {
        return scal_func::frexp(x, n);
    }
    static dynamic::function_name name() { return mdyf::frexp::eval(); };
};
template<>
struct frexp_helper<Complex>
{
    force_inline
    static Real eval(const Complex& x, Integer& n)
    {
        return scal_func::frexp(x, n);
    }
    static dynamic::function_name name() { return mdyf::frexp::eval(); };
};
template<>
struct frexp_helper<Float_complex>
{
    force_inline
    static Float eval(const Float_complex& x, Integer& n)
    {
        return scal_func::frexp(x, n);
    }
    static dynamic::function_name name() { return mdyf::frexp::eval(); };
};
template<>
struct frexp_helper<dynamic::object>
{
    static dynamic::object eval(const dynamic::object& x, Integer& n)
    {
        dynamic::object exp;
        dynamic::object d = dynamic::frexp(x, exp);
        n = dynamic::OInteger(exp, dynamic::from_object()).get();
        return d;
    }
    static dynamic::function_name name() { return mdyf::frexp::eval(); };
};

template<class T>
struct modf_helper
{};
template<>
struct modf_helper<Real>
{
    force_inline
    static Real eval(Real x, Real& n)
    {
        return scal_func::modf(x, n);
    }
    static dynamic::function_name name() { return mdyf::modf::eval(); };
};
template<>
struct modf_helper<Integer>
{
    force_inline
    static Real eval(Integer x, Real& n)
    {
        return modf_helper<Real>::eval(Real(x),n);
    }
    static dynamic::function_name name() { return mdyf::modf::eval(); };
};
template<>
struct modf_helper<Float>
{
    force_inline
    static Float eval(Float x, Float& n)
    {
        return scal_func::modf(x, n);
    }
    static dynamic::function_name name() { return mdyf::modf::eval(); };
};
template<>
struct modf_helper<Complex>
{
    force_inline
    static Real eval(const Complex& x, Real& n)
    {
        return scal_func::modf(x, n);
    }
    static dynamic::function_name name() { return mdyf::modf::eval(); };
};
template<>
struct modf_helper<Float_complex>
{
    force_inline
    static Float eval(const Float_complex& x, Float& n)
    {
        return scal_func::modf(x, n);
    }
    static dynamic::function_name name() { return mdyf::modf::eval(); };
};
template<>
struct modf_helper<dynamic::object>
{
    static dynamic::object eval(const dynamic::object& x, dynamic::object& n)
    {
        dynamic::object d = dynamic::modf(x, n);
        return d;
    }
    static dynamic::function_name name() { return mdyf::modf::eval(); };
};

template<class T>
struct nextabove_helper
{
    force_inline
    static T eval(const T& x)
    {
        return scal_func::nextabove(x);
    }
    static dynamic::function_name name() { return mdyf::nextabove::eval(); };
};
template<>
struct nextabove_helper<Integer>
{
    force_inline
    static Real eval(Integer x)
    {
        return scal_func::nextabove(Real(x));
    }
    static dynamic::function_name name() { return mdyf::nextabove::eval(); };
};
template<>
struct nextabove_helper<dynamic::object>
{
    static dynamic::object eval(const dynamic::object& x)
    {
        return dynamic::nextabove(x);
    }
    static dynamic::function_name name() { return mdyf::nextabove::eval(); };
};

template<class T>
struct nextbelow_helper
{
    force_inline
    static T eval(const T& x)
    {
        return scal_func::nextbelow(x);
    }
    static dynamic::function_name name() { return mdyf::nextbelow::eval(); };
};
template<>
struct nextbelow_helper<Integer>
{
    force_inline
    static Real eval(Integer x)
    {
        return scal_func::nextbelow(Real(x));
    }
    static dynamic::function_name name() { return mdyf::nextbelow::eval(); };
};
template<>
struct nextbelow_helper<dynamic::object>
{
    static dynamic::object eval(const dynamic::object& x)
    {
        return dynamic::nextbelow(x);
    }
    static dynamic::function_name name() { return mdyf::nextbelow::eval(); };
};

template<class T>
struct fma_f_helper
{
    force_inline
    static T eval(const T& a, const T& b, const T& c)
    {
        return scal_func::fma_f(a, b, c);
    }

    static dynamic::function_name name() { return mdyf::fma_f::eval(); };
};

template<>
struct fma_f_helper<dynamic::object>
{
    static dynamic::object eval(const dynamic::object& a, const dynamic::object& b, const dynamic::object& c)
    {
        return dynamic::fma_f(a, b, c);
    }
    static dynamic::function_name name() { return mdyf::fma_f::eval(); };
};

template<class T>
struct fms_f_helper
{
    force_inline
    static T eval(const T& a, const T& b, const T& c)
    {
        return scal_func::fms_f(a, b, c);
    }

    static dynamic::function_name name() { return mdyf::fms_f::eval(); };
};

template<>
struct fms_f_helper<dynamic::object>
{
    static dynamic::object eval(const dynamic::object& a, const dynamic::object& b, const dynamic::object& c)
    {
        return dynamic::fms_f(a, b, c);
    }
    static dynamic::function_name name() { return mdyf::fms_f::eval(); };
};

//
template<class T>
struct fma_a_helper
{
    force_inline
    static T eval(const T& a, const T& b, const T& c)
    {
        return scal_func::fma_a(a, b, c);
    }

    static dynamic::function_name name() { return mdyf::fma_a::eval(); };
};

template<>
struct fma_a_helper<dynamic::object>
{
    static dynamic::object eval(const dynamic::object& a, const dynamic::object& b, const dynamic::object& c)
    {
        return dynamic::fma_a(a, b, c);
    }
    static dynamic::function_name name() { return mdyf::fma_a::eval(); };
};

template<class T>
struct fms_a_helper
{
    force_inline
    static T eval(const T& a, const T& b, const T& c)
    {
        return scal_func::fms_a(a, b, c);
    }

    static dynamic::function_name name() { return mdyf::fms_f::eval(); };
};

template<>
struct fms_a_helper<dynamic::object>
{
    static dynamic::object eval(const dynamic::object& a, const dynamic::object& b, const dynamic::object& c)
    {
        return dynamic::fms_a(a, b, c);
    }
    static dynamic::function_name name() { return mdyf::fms_a::eval(); };
};

template<class T>
struct dot2_a_helper
{
    force_inline
    static T eval(const T& a, const T& b, const T& c, const T& d)
    {
        return scal_func::dot2_a(a, b, c, d);
    }
    static dynamic::function_name name() { return mdyf::fms::dot2_a(); };
};
template<>
struct dot2_a_helper<dynamic::object>
{
    static dynamic::object eval(const dynamic::object& a, const dynamic::object& b, const dynamic::object& c, const dynamic::object& d)
    {
        return dynamic::dot2_a(a, b, c, d);
    }
    static dynamic::function_name name() { return mdyf::dot2_a::eval(); };
};

}}}
