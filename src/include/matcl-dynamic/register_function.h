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

#include "matcl-dynamic/details/register_function_impl.h"

namespace matcl { namespace dynamic
{

// register unifier type; 
// in order to do registration this type must be instantiated
// Derived         - type of class, that derives from this class
//                   Derived class must implement a function
//                     static Type eval();
//                   which marks the type Type as an unifier; this function
//                   will not be called
//
// example:
//     struct unif_real : register_unifier<unif_real>
//     {
//         static OReal eval() { return OReal();};
//     };
template<class Derived>
struct register_unifier
{
    private:
        virtual void initialize()
        { 
            using hook_type = details::make_function_type<Derived, details::unifier_tag, 
                                    decltype(&Derived::eval),&Derived::eval>;
            return hook_type::get().initialize(); 
        };
};

// implicit conversion; user defined promotions
struct MATCL_DYN_EXPORT convert_promotion   { static function_name eval(); };

// implicit conversion; user defined lossless conversions changing kind
struct MATCL_DYN_EXPORT convert_equivalent  { static function_name eval(); };

// implicit conversion; user defined conversions with data loss
struct MATCL_DYN_EXPORT convert_decay       { static function_name eval(); };

// explicit conversion
struct MATCL_DYN_EXPORT convert_explicit    { static function_name eval(); };

// explicit cast, called when the function cast is called
struct MATCL_DYN_EXPORT convert_cast        { static function_name eval(); };

// register conversion from From type to To type; 
// in order to do registration this type must be instantiated
// Derived         - type of class, that derives from this class
//                   Derived class must implement a function
//                     static To eval(const From& val);
//                   for given types LHS and RHS
// Convert_type    - type of conversion, must be one of types: convert_promotion, 
//                     convert_equivalent, convert_decay, convert_explicit,
//                     convert_cast
//
// example:
//     struct conv_int_float : register_convert<conv_int_float, numeric_convert_decay>
//     {
//     	static OFloat eval(const OInteger& val) { return Float(val.get()); };
//     };
template<class Derived, class Convert_type>
struct register_convert
{
    private:
        virtual void initialize()
        { 
            using hook_type = details::make_function_type<Derived, Convert_type, 
                                    decltype(&Derived::eval),&Derived::eval>;
            return hook_type::get().initialize(); 
        };
};

// register function that should be called, when assignment lhs = rhs
// is called; in order to do registration this type must be instantiated
// Derived         - type of class, that derives from this class
//                   Derived class must implement a function
//                     static void eval(LHS& lhs, const RHS& rhs);
//                   for given types LHS and RHS
// 
// example:
//     struct assign_compl : register_assign<assign_compl>
//     {
//     	static void eval(OComplex& to, const OReal& from)
//         { 
//             to.get_unique() = Complex(from.get()); 
//         };
//     };
template<class Derived>
struct register_assign
{
    private:
        virtual void initialize()
        { 
            using hook_type = details::make_function_type<Derived, details::assign_tag, 
                                    decltype(&Derived::eval),&Derived::eval>;
            return hook_type::get().initialize(); 
        };
};

// register function that should be called, when function with given name
// is called on object of given type; in order to do registration this type
// must be instantiated;
//
// it is possible to register a function that takes non constant reference
// as an argument; however non constant null object can bind to this 
// reference (null object is converted to an object of the same type as the
// reference type and default initialized); if this nonconstant reference
// argument is both input and output, then it must be split to two arguments -
// one being constant reference (the input) and the second being nonconst
// reference (the output)
// Derived         - type of class, that derives from this class
//                   Derived class must implement a function
//                     static [Type_ret] eval([Type_1] arg1, ...);
//                   for some return type Type_ret and one or more input types
//                   Type_1, ...
// Function_name   - functor returning function name
// 
// example:
//     struct op_minus : register_function<op_minus, functions::op_minus>
//     {
//         static OFloat eval(const OFloat& obj1, const OFloat& obj2)
//         {
//             return OFloat(obj1.get() - obj2.get());
//         };
//     };
template<class Derived, class Function_name>
struct register_function
{
    private:
        virtual void initialize()    
        { 
            using hook_type = details::make_function_type<Derived, Function_name, 
                                    decltype(&Derived::eval), &Derived::eval>;
            return hook_type::get().initialize(); 
        };
};

// register functions dependent on types; 
// restrictions on registered functions are the same as for register_function;
// this function can be called by eval_function_template; if some of template
// arguments are of type Template, then these types are matched againts actual types
// given in the call site; these matched types are then supplied to this registered
// fuction as first k arguments of type OType storing infered types;
//
// Derived         - type of class, that derives from this class
//                   Derived class must implement a function
//                     static [Type_ret] eval(OType ti, ..., [Type_1] arg1, ...);
//                   for some return type Type_ret and one or more input types
//                   Type_1, ...
// Function_name   - functor returning function name
// Templates       - variadic number of additional object types; reference and cv
//                     qualifiers are ignored; if some of argument is Template, then
//                     any type will match this template argument, otherwise an
//                     exact match is required
// 
// example:
//     struct impl_rand : register_function_template<impl_rand, functions::rand, OReal>
//     {
//         static OReal eval()
//         {
//             return OReal(randn());
//         };
//     };
template<class Derived, class Function_name, class ... Templates>
struct register_function_template
{
    private:
        virtual void initialize()    
        { 
            using hook_type = details::make_function_template_type<Derived, Function_name, 
                                decltype(&Derived::eval), &Derived::eval, nullptr, Templates ...>;
            return hook_type::get().initialize(); 
        };
};

// register functions dependent on types; return type is calculated by eval_return
// function; see also register_function_template; the first argument of registered
// function must have OType type, that stores return type (if return type is Template,
// followed by zero or more arguments of type OType storing inferred template 
// arguments and finally n ordinary function arguments
//
// Derived         - type of class, that derives from this class
//                   Derived class must implement a function
//                     static [Ret] eval(OType ti, ..., [Type_1] arg1, ...);
//                   with one or more input types Type_1, ...; return type may have 
///                  Template type; and a function
//                     static Type eval_return(int n_template, const Type templates[], 
//                                     int n_arg, const Type args[]);
//                   where n_template is number of template arguments with types
//                   given by templates, n_arg is number of supplied arguments with
//                   types arg; elements in args array are the same as types of function
//                   arguments except of arguments of type Template which  are substituted
//                   by types of supplied arguments; if eval_return returns t being a null
//                   Type then this function is removed from the overload set, otherwise
//                   return type of this function is t if Ret = Template (and Ret otherwise)
// Function_name   - functor returning function name
// Templates       - variadic number of additional object types; reference and cv
//                     qualifiers are ignored; if some of argument is Template, then
//                     any type will match this template argument, otherwise an
//                     exact match is required
// 
// example:
//     struct kron : register_function_template<kron, functions::kron>
//     {
//         static Template eval(const OType& ret_t, const Template& x, const Template& y)
//         {
//             return object(ret_t.get(), x.get() * y.get());
//         };
//         static Type eval_return(int n_t, const Type t[], int n_arg, const Type args[])
//         {
//             return operations::return_type(functions::op_mul::eval(), n_arg, args);
//         }
//     };
//
//     struct impl_rand : register_function_template<impl_rand, functions::rand, Template>
//     {
//         static Template eval(const OType& ret_t, const OType& templ_t)
//         {
//             return object(ret_t.get(), randn());
//         };
//         static Type eval_return(int n_t, const Type t[], int n_arg, const Type args[])
//         {
//             return args[0];
//         }
//     };
template<class Derived, class Function_name, class ... Templates>
struct register_function_template_return
{
    private:
        virtual void initialize()    
        { 
            using hook_type = details::make_function_template_type<Derived, Function_name, 
                                decltype(&Derived::eval), &Derived::eval, &Derived::eval_return,
                                Templates ...>;
            return hook_type::get().initialize(); 
        };
};

// register function that should be called, when function with given name
// is called on object of given type; in order to do registration this type
// must be instantiated; see also register_function class
// Func            - type of function 
// f               - pointer to function of type Func
// Function_name   - functor returning function name
// 
// example:
//     double real(const double& v) { return v; };
// 
//     template struct register_function_ptr<decltype(real),real, functions::real>;
template<class Func, Func f, class Function_name>
struct register_function_ptr : register_function<details::register_function_eval<Func,f>, 
                                                Function_name>
{};

};};