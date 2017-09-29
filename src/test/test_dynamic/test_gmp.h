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

#include "matcl-dynamic/matcl_dynamic.h"
#include "matcl-scalar/lib_functions/func_unary.h"
#include "matcl-scalar/lib_functions/func_binary.h"
#include "matcl-scalar/objects/typed_object_functions.h"

#include "scalar.h"
#include "scalar_ext.h"
#include <vector>

namespace matcl { namespace test
{

void test_gmp();
void test_gmp_bin();

namespace md = matcl::details;

class gmp_tester
{
    public:
        void    make();

    private:
        void    test_func(double dif, const std::string& test_name);
        void    test_reim();
        bool    test_reim_impl();
        double  test_reim(const Scalar& s, Integer code);
        void    test_uminus();
        bool    test_uminus_impl();
        double  test_uminus(const Scalar& s, Integer code);
        void    test_is();
        void    test_next();
        void    test_signbit();
        void    test_eps();
        bool    test_is_impl();
        bool    test_next_impl();
        bool    test_sign_impl();
        bool    test_eps_impl();
        double  test_is(const Scalar& s, Integer code);
        double  test_next(const Scalar& s, Integer code);
        double  test_sign(const Scalar& s, Integer code);
        double  test_eps(const Scalar& s, Integer code);

        double  test_mp_int();
        double  test_mp_float();
        double  test_mp_complex();
        double  test_mp_rational();

        template<class Func>
        void    test_scalar_func();

        template<class Func>
        double  test_scalar(const Scalar& s, Integer code);

        void    test_ldexp_p3();
        void    test_ldexp_m3();
        void    test_scalbn_p3();
        void    test_scalbn_m3();
        void    test_frexp();
        void    test_modf_frac();
        void    test_modf_int();
        void    test_logb();
        void    test_ilogb();

        void    test_inv();
        void    test_invs();
        void    test_sqrt();
        void    test_cbrt();
        void    test_sqrt_c();
        void    test_exp();
        void    test_expm1();
        void    test_expi();
        void    test_exp2(); 
        void    test_exp10(); 
        void    test_log();
        void    test_log1p();
        void    test_log2();
        void    test_log10();
        void    test_log_c();
        void    test_log1p_c();
        void    test_log2_c();
        void    test_log10_c();
        void    test_sin();
        void    test_cos();
        void    test_tan();
        void    test_cot();
        void    test_sec();
        void    test_csc();
        void    test_sinh();
        void    test_cosh();
        void    test_tanh();
        void    test_coth();
        void    test_sech();
        void    test_csch();
        void    test_asin();
        void    test_asin_c();
        void    test_acos();
        void    test_acos_c();
        void    test_atan();
        void    test_asinh();
        void    test_acosh();
        void    test_acosh_c();
        void    test_atanh();
        void    test_atanh_c();
        void    test_floor();
        void    test_ceil();
        void    test_round();
        void    test_trunc();
        void    test_ifloor();
        void    test_iceil();
        void    test_iround();
        void    test_itrunc();
        void    test_sign();
        void    test_isign();
        void    test_fpclassify();
        void    test_combinatorics();
        void    test_combinatorics2();
};

class gmp_tester_bin
{
    public:
        void            make();

    private:
        void            test_add();
        void            test_sub();
        void            test_mul();
        void            test_mmul();
        void            test_div();
        void            test_div2();
        void            test_div_0();
        void            test_div_1();
        void            test_idiv();
        void            test_pow();
        void            test_pow_c();

        void            test_eeq();
        void            test_neq();
        void            test_eeq_nan();
        void            test_neq_nan();
        void            test_leq();
        void            test_geq();
        void            test_lt();
        void            test_gt();

        void            test_atan2();
        void            test_hypot();
        void            test_mod(); 
        void            test_rem();
        void            test_fdim();
        void            test_nextafter();
        void            test_copysign();
        void            test_pow_spec();
        void            test_min();
        void            test_max();
        void            test_op_and();
        void            test_op_or();
        void            test_op_xor();
        void            test_elem_and();
        void            test_elem_or();
        void            test_elem_xor();

        void            test_atan2_obj();
        void            test_hypot_obj();
        void            test_copysign_obj();
        void            test_nextafter_obj();
        void            test_mod_obj();
        void            test_rem_obj();
        void            test_fdim_obj();

        void            test_add_obj();
        void            test_sub_obj();
        void            test_mul_obj();
        void            test_mmul_obj();
        void            test_div_obj();
        void            test_div2_obj();
        void            test_div_0_obj();
        void            test_div_1_obj();
        void            test_idiv_obj();
        void            test_pow_obj();
        void            test_pow_c_obj();
        void            test_min_obj();
        void            test_max_obj();
        void            test_op_and_obj();
        void            test_op_or_obj();
        void            test_op_xor_obj();
        void            test_elem_and_obj();
        void            test_elem_or_obj();
        void            test_elem_xor_obj();

        void            test_eeq_obj();
        void            test_neq_obj();
        void            test_leq_obj();
        void            test_geq_obj();
        void            test_lt_obj();
        void            test_gt_obj();
        void            test_eeq_nan_obj();
        void            test_neq_nan_obj();

        double          test_atan2(const Scalar& s1, const Scalar& s2, Integer code);
        double          test_hypot(const Scalar& s1, const Scalar& s2, Integer code);
        double          test_mod(const Scalar& s1, const Scalar& s2, Integer code);
        double          test_rem(const Scalar& s1, const Scalar& s2, Integer code);
        double          test_fdim(const Scalar& s1, const Scalar& s2, Integer code);
        double          test_nextafter(const Scalar& s1, const Scalar& s2, Integer code);
        double          test_copysign(const Scalar& s1, const Scalar& s2, Integer code);

        double          test_add(const Scalar& s1, const Scalar& s2, Integer code);
        double          test_sub(const Scalar& s1, const Scalar& s2, Integer code);
        double          test_mul(const Scalar& s1, const Scalar& s2, Integer code);
        double          test_mmul(const Scalar& s1, const Scalar& s2, Integer code);
        double          test_div(const Scalar& s1, const Scalar& s2, Integer code);
        double          test_div2(const Scalar& s1, const Scalar& s2, Integer code);
        double          test_div_0(const Scalar& s1, const Scalar& s2, Integer code);
        double          test_div_1(const Scalar& s1, const Scalar& s2, Integer code);
        double          test_idiv(const Scalar& s1, const Scalar& s2, Integer code);
        double          test_pow(const Scalar& s1, const Scalar& s2, Integer code);
        double          test_pow_c(const Scalar& s1, const Scalar& s2, Integer code);
        double          test_min(const Scalar& s1, const Scalar& s2, Integer code);
        double          test_max(const Scalar& s1, const Scalar& s2, Integer code);
        double          test_op_or(const Scalar& s1, const Scalar& s2, Integer code);
        double          test_op_xor(const Scalar& s1, const Scalar& s2, Integer code);
        double          test_op_and(const Scalar& s1, const Scalar& s2, Integer code);
        double          test_elem_or(const Scalar& s1, const Scalar& s2, Integer code);
        double          test_elem_xor(const Scalar& s1, const Scalar& s2, Integer code);
        double          test_elem_and(const Scalar& s1, const Scalar& s2, Integer code);

        double          test_eeq(const Scalar& s1, const Scalar& s2, Integer code);
        double          test_neq(const Scalar& s1, const Scalar& s2, Integer code);
        double          test_eeq_nan(const Scalar& s1, const Scalar& s2, Integer code);
        double          test_neq_nan(const Scalar& s1, const Scalar& s2, Integer code);
        double          test_leq(const Scalar& s1, const Scalar& s2, Integer code);
        double          test_geq(const Scalar& s1, const Scalar& s2, Integer code);
        double          test_lt(const Scalar& s1, const Scalar& s2, Integer code);
        double          test_gt(const Scalar& s1, const Scalar& s2, Integer code);

        double          test_atan2_obj(const Scalar_ext& s1, const Scalar_ext& s2, Integer code);
        double          test_hypot_obj(const Scalar_ext& s1, const Scalar_ext& s2, Integer code);
        double          test_copysign_obj(const Scalar_ext& s1, const Scalar_ext& s2, Integer code);
        double          test_nextafter_obj(const Scalar_ext& s1, const Scalar_ext& s2, Integer code);
        double          test_mod_obj(const Scalar_ext& s1, const Scalar_ext& s2, Integer code);
        double          test_rem_obj(const Scalar_ext& s1, const Scalar_ext& s2, Integer code);
        double          test_fdim_obj(const Scalar_ext& s1, const Scalar_ext& s2, Integer code);
        double          test_pow_obj(const Scalar_ext& s1, const Scalar_ext& s2, Integer code);
        double          test_pow_c_obj(const Scalar_ext& s1, const Scalar_ext& s2, Integer code);
        double          test_min_obj(const Scalar_ext& s1, const Scalar_ext& s2, Integer code);
        double          test_max_obj(const Scalar_ext& s1, const Scalar_ext& s2, Integer code);
        double          test_op_or_obj(const Scalar_ext& s1, const Scalar_ext& s2, Integer code);
        double          test_op_xor_obj(const Scalar_ext& s1, const Scalar_ext& s2, Integer code);
        double          test_op_and_obj(const Scalar_ext& s1, const Scalar_ext& s2, Integer code);
        double          test_elem_or_obj(const Scalar_ext& s1, const Scalar_ext& s2, Integer code);
        double          test_elem_xor_obj(const Scalar_ext& s1, const Scalar_ext& s2, Integer code);
        double          test_elem_and_obj(const Scalar_ext& s1, const Scalar_ext& s2, Integer code);

        double          test_add_obj(const Scalar_ext& s1, const Scalar_ext& s2, Integer code);
        double          test_sub_obj(const Scalar_ext& s1, const Scalar_ext& s2, Integer code);
        double          test_mul_obj(const Scalar_ext& s1, const Scalar_ext& s2, Integer code);
        double          test_mmul_obj(const Scalar_ext& s1, const Scalar_ext& s2, Integer code);
        double          test_div_obj(const Scalar_ext& s1, const Scalar_ext& s2, Integer code);
        double          test_div2_obj(const Scalar_ext& s1, const Scalar_ext& s2, Integer code);
        double          test_div_0_obj(const Scalar_ext& s1, const Scalar_ext& s2, Integer code);
        double          test_div_1_obj(const Scalar_ext& s1, const Scalar_ext& s2, Integer code);
        double          test_idiv_obj(const Scalar_ext& s1, const Scalar_ext& s2, Integer code);

        double          test_eeq_obj(const Scalar_ext& s1, const Scalar_ext& s2, Integer code);
        double          test_neq_obj(const Scalar_ext& s1, const Scalar_ext& s2, Integer code);
        double          test_eeq_nan_obj(const Scalar_ext& s1, const Scalar_ext& s2, Integer code);
        double          test_neq_nan_obj(const Scalar_ext& s1, const Scalar_ext& s2, Integer code);
        double          test_leq_obj(const Scalar_ext& s1, const Scalar_ext& s2, Integer code);
        double          test_geq_obj(const Scalar_ext& s1, const Scalar_ext& s2, Integer code);
        double          test_lt_obj(const Scalar_ext& s1, const Scalar_ext& s2, Integer code);
        double          test_gt_obj(const Scalar_ext& s1, const Scalar_ext& s2, Integer code);
};

struct Allow_complex
{
    static const bool is_complex_allowed    = true;
};

struct Exp_func : Allow_complex
{
    static std::string name() { return "exp"; };

    template<class T1, class ... Args>
    static auto eval(const T1& a1, Args&& ... args) -> decltype(exp(std::declval<T1>()))
    {
        return exp(a1, std::forward<Args>(args)...);
    };
};

struct Expm1_func : Allow_complex
{
    static std::string name() { return "expm1"; };

    template<class T1, class ... Args>
    static auto eval(const T1& a1, Args&& ... args) -> decltype(expm1(std::declval<T1>()))
    {
        return expm1(a1, std::forward<Args>(args)...);
    };
};

struct Expi_func : Allow_complex
{
    static std::string name() { return "expi"; };

    template<class T1, class ... Args>
    static auto eval(const T1& a1, Args&& ... args) -> decltype(expi(std::declval<T1>()))
    {
        return expi(a1, std::forward<Args>(args)...);
    };
};

struct Exp2_func : Allow_complex
{
    static std::string name() { return "exp2"; };

    template<class T1, class ... Args>
    static auto eval(const T1& a1, Args&& ... args) -> decltype(exp2(std::declval<T1>()))
    {
        return exp2(a1, std::forward<Args>(args)...);
    };
};

struct Exp10_func : Allow_complex
{
    static std::string name() { return "exp10"; };

    template<class T1, class ... Args>
    static auto eval(const T1& a1, Args&& ... args) -> decltype(exp10(std::declval<T1>()))
    {
        return exp10(a1, std::forward<Args>(args)...);
    };
};

struct Log_func : Allow_complex
{
    static std::string name() { return "log"; };

    template<class T1, class ... Args>
    static auto eval(const T1& a1, Args&& ... args) -> decltype(log(std::declval<T1>()))
    {
        return log(a1, std::forward<Args>(args)...);
    };
};

struct Log1p_func : Allow_complex
{
    static std::string name() { return "log1p"; };

    template<class T1, class ... Args>
    static auto eval(const T1& a1, Args&& ... args) -> decltype(log1p(std::declval<T1>()))
    {
        return log1p(a1, std::forward<Args>(args)...);
    };
};

struct Log2_func : Allow_complex
{
    static std::string name() { return "log2"; };

    template<class T1, class ... Args>
    static auto eval(const T1& a1, Args&& ... args) -> decltype(log2(std::declval<T1>()))
    {
        return log2(a1, std::forward<Args>(args)...);
    };
};
struct Log10_func : Allow_complex
{
    static std::string name() { return "log10"; };

    template<class T1, class ... Args>
    static auto eval(const T1& a1, Args&& ... args) -> decltype(log10(std::declval<T1>()))
    {
        return log10(a1, std::forward<Args>(args)...);
    };
};

struct Log_c_func : Allow_complex
{
    static std::string name() { return "log_c"; };

    template<class T1, class ... Args>
    static auto eval(const T1& a1, Args&& ... args) -> decltype(log_c(std::declval<T1>()))
    {
        return log_c(a1, std::forward<Args>(args)...);
    };
};
struct Log1p_c_func : Allow_complex
{
    static std::string name() { return "log1p_c"; };

    template<class T1, class ... Args>
    static auto eval(const T1& a1, Args&& ... args) -> decltype(log1p_c(std::declval<T1>()))
    {
        return log1p_c(a1, std::forward<Args>(args)...);
    };
};

struct Log2_c_func : Allow_complex
{
    static std::string name() { return "log2_c"; };

    template<class T1, class ... Args>
    static auto eval(const T1& a1, Args&& ... args) -> decltype(log2_c(std::declval<T1>()))
    {
        return log2_c(a1, std::forward<Args>(args)...);
    };
};
struct Log10_c_func : Allow_complex
{
    static std::string name() { return "log10_c"; };

    template<class T1, class ... Args>
    static auto eval(const T1& a1, Args&& ... args) -> decltype(log10_c(std::declval<T1>()))
    {
        return log10_c(a1, std::forward<Args>(args)...);
    };
};

struct Log_c2_func : Allow_complex
{
    static std::string name() { return "log_c"; };

    template<class T1, class ... Args>
    static auto eval(const T1& a1, Args&& ... args) -> decltype(log_c(std::declval<T1>()))
    {
        return log_c(real(a1), std::forward<Args>(args)...);
    };
};

struct Sin_func : Allow_complex
{
    static std::string name() { return "sin"; };

    template<class T1, class ... Args>
    static auto eval(const T1& a1, Args&& ... args) -> decltype(sin(std::declval<T1>()))
    {
        return sin(a1, std::forward<Args>(args)...);
    };
};

struct Cos_func : Allow_complex
{
    static std::string name() { return "cos"; };

    template<class T1, class ... Args>
    static auto eval(const T1& a1, Args&& ... args) -> decltype(cos(std::declval<T1>()))
    {
        return cos(a1, std::forward<Args>(args)...);
    };
};

struct Tan_func : Allow_complex
{
    static std::string name() { return "tan"; };

    template<class T1, class ... Args>
    static auto eval(const T1& a1, Args&& ... args) -> decltype(tan(std::declval<T1>()))
    {
        return tan(a1, std::forward<Args>(args)...);
    };
};

struct Cot_func : Allow_complex
{
    static std::string name() { return "cot"; };

    template<class T1, class ... Args>
    static auto eval(const T1& a1, Args&& ... args) -> decltype(cot(std::declval<T1>()))
    {
        return cot(a1, std::forward<Args>(args)...);
    };
};

struct Sec_func : Allow_complex
{
    static std::string name() { return "sec"; };

    template<class T1, class ... Args>
    static auto eval(const T1& a1, Args&& ... args) -> decltype(sec(std::declval<T1>()))
    {
        return sec(a1, std::forward<Args>(args)...);
    };
};

struct Csc_func : Allow_complex
{
    static std::string name() { return "csc"; };

    template<class T1, class ... Args>
    static auto eval(const T1& a1, Args&& ... args) -> decltype(csc(std::declval<T1>()))
    {
        return csc(a1, std::forward<Args>(args)...);
    };
};

struct Sinh_func : Allow_complex
{
    static std::string name() { return "sinh"; };

    template<class T1, class ... Args>
    static auto eval(const T1& a1, Args&& ... args) -> decltype(sinh(std::declval<T1>()))
    {
        return sinh(a1, std::forward<Args>(args)...);
    };
};

struct Cosh_func : Allow_complex
{
    static std::string name() { return "cosh"; };

    template<class T1, class ... Args>
    static auto eval(const T1& a1, Args&& ... args) -> decltype(cosh(std::declval<T1>()))
    {
        return cosh(a1, std::forward<Args>(args)...);
    };
};

struct Tanh_func : Allow_complex
{
    static std::string name() { return "tanh"; };

    template<class T1, class ... Args>
    static auto eval(const T1& a1, Args&& ... args) -> decltype(tanh(std::declval<T1>()))
    {
        return tanh(a1, std::forward<Args>(args)...);
    };
};

struct Coth_func : Allow_complex
{
    static std::string name() { return "coth"; };

    template<class T1, class ... Args>
    static auto eval(const T1& a1, Args&& ... args) -> decltype(coth(std::declval<T1>()))
    {
        return coth(a1, std::forward<Args>(args)...);
    };
};

struct Sech_func : Allow_complex
{
    static std::string name() { return "sech"; };

    template<class T1, class ... Args>
    static auto eval(const T1& a1, Args&& ... args) -> decltype(sech(std::declval<T1>()))
    {
        return sech(a1, std::forward<Args>(args)...);
    };
};

struct Csch_func : Allow_complex
{
    static std::string name() { return "csch"; };

    template<class T1, class ... Args>
    static auto eval(const T1& a1, Args&& ... args) -> decltype(csch(std::declval<T1>()))
    {
        return csch(a1, std::forward<Args>(args)...);
    };
};

struct Floor_func : Allow_complex
{
    static std::string name() { return "floor"; };

    template<class T1, class ... Args>
    static auto eval(const T1& a1, Args&& ... args) -> decltype(floor(std::declval<T1>()))
    {
        return floor(a1, std::forward<Args>(args)...);
    };
};

struct Ceil_func : Allow_complex
{
    static std::string name() { return "ceil"; };

    template<class T1, class ... Args>
    static auto eval(const T1& a1, Args&& ... args) -> decltype(ceil(std::declval<T1>()))
    {
        return ceil(a1, std::forward<Args>(args)...);
    };
};

struct Round_func : Allow_complex
{
    static std::string name() { return "round"; };

    template<class T1, class ... Args>
    static auto eval(const T1& a1, Args&& ... args) -> decltype(round(std::declval<T1>()))
    {
        return round(a1, std::forward<Args>(args)...);
    };
};

struct Factorial_func : Allow_complex
{
    static std::string name() { return "factorial"; };

    template<class T1, class ... Args>
    static auto eval(const T1& a1, Args&& ... args) 
        -> decltype(factorial<T1>(Integer()))
    {
        return factorial<T1>(21);
    };
    template<class T1, class ... Args>
    static auto eval(const dynamic::object_type<T1>& a1, Args&& ... args) 
        -> decltype(dynamic::factorial<dynamic::object_type<T1>>(Integer()))
    {
        return dynamic::factorial<dynamic::object_type<T1>>(21);
    };
    template<class ... Args>
    static auto eval(const dynamic::object& a1, Args&& ... args) 
        -> dynamic::object
    {
        return dynamic::factorial(a1.get_type(),21);
    };
};

struct Double_factorial_func : Allow_complex
{
    static std::string name() { return "double_factorial"; };

    template<class T1, class ... Args>
    static auto eval(const T1& a1, Args&& ... args) 
        -> decltype(double_factorial<T1>(Integer()))
    {
        return double_factorial<T1>(21);
    };
    template<class T1, class ... Args>
    static auto eval(const dynamic::object_type<T1>& a1, Args&& ... args) 
        -> decltype(dynamic::double_factorial<dynamic::object_type<T1>>(Integer()))
    {
        return dynamic::double_factorial<dynamic::object_type<T1>>(21);
    };
    template<class ... Args>
    static auto eval(const dynamic::object& a1, Args&& ... args) 
        -> dynamic::object
    {
        return dynamic::double_factorial(a1.get_type(),21);
    };
};

struct Binomial_coefficient_func : Allow_complex
{
    static std::string name() { return "binomial_coefficient"; };

    template<class T1, class ... Args>
    static auto eval(const T1& a1, Args&& ... args) 
        -> decltype(binomial_coefficient<T1>(Integer(),Integer()))
    {
        return binomial_coefficient<T1>(21,15);
    };

    template<class T1, class ... Args>
    static auto eval(const dynamic::object_type<T1>& a1, Args&& ... args) 
        -> decltype(dynamic::binomial_coefficient<const dynamic::object_type<T1>>(Integer(),Integer()))
    {
        return dynamic::binomial_coefficient<const dynamic::object_type<T1>>(21,15);
    };
    template<class ... Args>
    static auto eval(const dynamic::object& a1, Args&& ... args) 
        -> dynamic::object
    {
        return dynamic::binomial_coefficient(a1.get_type(),21,15);
    };
};
struct Trunc_func : Allow_complex
{
    static std::string name() { return "trunc"; };

    template<class T1, class ... Args>
    static auto eval(const T1& a1, Args&& ... args) -> decltype(trunc(std::declval<T1>()))
    {
        return trunc(a1, std::forward<Args>(args)...);
    };
};

struct Notcomplex_func
{
    static const bool is_complex_allowed    = false;

    template<class ... Args>
    static Real eval(const mp_complex&,  Args&& ... )   { return 0.; };
    template<class ... Args>
    static Real eval(const Complex&,  Args&& ... )      { return 0.; };
    template<class ... Args>
    static Real eval(const Float_complex&,  Args&& ... ){ return 0.; };

    template<class ... Args>
    static OReal eval(const MP_complex&,  Args&& ... )      { return OReal(0.); };
    template<class ... Args>
    static OReal eval(const OComplex&,  Args&& ... )        { return OReal(0.); };
    template<class ... Args>
    static OReal eval(const OFloat_complex&,  Args&& ... )  { return OReal(0.); };
};

template<class T>
struct is_complex_val                           { static const bool value = md::is_complex<T>::value; };

template<class T>
struct is_complex_val<dynamic::object_type<T>>  { static const bool value = md::is_complex<T>::value; };

template<class T1, class T2>
struct enable_not_complex
    :std::enable_if<is_complex_val<T1>::value == false && is_complex_val<T2>::value == false>
{};

template<class T1, class T2>
struct enable_complex_obj
    :std::enable_if<(is_complex_val<T1>::value == true || is_complex_val<T2>::value == true)
                    && (md::is_object<T1>::value == true || md::is_object<T2>::value == true)>
{};

template<class T1, class T2>
struct enable_complex_nobj
    :std::enable_if<(is_complex_val<T1>::value == true || is_complex_val<T2>::value == true)
                    && (md::is_object<T1>::value == false && md::is_object<T2>::value == false)>
{};

struct IFloor_func : Notcomplex_func
{
    using Notcomplex_func::eval;

    static std::string name() { return "ifloor"; };

    template<class T1, class ... Args>
    static auto eval(const T1& a1, Args&& ... args) -> decltype(ifloor(std::declval<T1>()))
    {
        return ifloor(a1, std::forward<Args>(args)...);
    };
};

struct ICeil_func : Notcomplex_func
{
    using Notcomplex_func::eval;

    static std::string name() { return "iceil"; };

    template<class T1, class ... Args>
    static auto eval(const T1& a1, Args&& ... args) -> decltype(iceil(std::declval<T1>()))
    {
        return iceil(a1, std::forward<Args>(args)...);
    };
};

struct IRound_func : Notcomplex_func
{
    using Notcomplex_func::eval;

    static std::string name() { return "iround"; };

    template<class T1, class ... Args>
    static auto eval(const T1& a1, Args&& ... args) -> decltype(iround(std::declval<T1>()))
    {
        return iround(a1, std::forward<Args>(args)...);
    };
};

struct ITrunc_func : Notcomplex_func
{
    using Notcomplex_func::eval;

    static std::string name() { return "itrunc"; };

    template<class T1, class ... Args>
    static auto eval(const T1& a1, Args&& ... args) -> decltype(itrunc(std::declval<T1>()))
    {
        return itrunc(a1, std::forward<Args>(args)...);
    };
};

struct Sign_func : Allow_complex
{
    static std::string name() { return "sign"; };

    template<class T1, class ... Args>
    static auto eval(const T1& a1, Args&& ... args) -> decltype(sign(std::declval<T1>()))
    {
        return sign(a1, std::forward<Args>(args)...);
    };
};

struct ISign_func : Notcomplex_func
{
    using Notcomplex_func::eval;

    static std::string name() { return "isign"; };

    template<class T1, class ... Args>
    static auto eval(const T1& a1, Args&& ... args) -> decltype(isign(std::declval<T1>()))
    {
        return isign(a1, std::forward<Args>(args)...);
    };
};

template<class T, bool Is_enum = std::is_enum<T>::value>
struct convert_enum             { using type = T; };

template<class T>
struct convert_enum<T,true>     { using type = Integer; };

struct Fpclassify_func : Notcomplex_func
{
    using Notcomplex_func::eval;

    static std::string name() { return "fpclassify"; };

    template<class T1, class ... Args>
    static auto eval(const T1& a1, Args&& ... args) 
        -> typename convert_enum<decltype(fpclassify(std::declval<T1>()))>::type
    {
        using Ret = typename convert_enum<decltype(fpclassify(std::declval<T1>()))>::type;
        return Ret(fpclassify(a1, std::forward<Args>(args)...));
    };
};

struct Uminus_func : Allow_complex
{
    static std::string name() { return "uminus"; };

    template<class T1, class ... Args>
    static auto eval(const T1& a1, Args&& ... args) -> decltype(uminus(std::declval<T1>()))
    {
        return uminus(a1, std::forward<Args>(args)...);
    };
};

struct Real_func : Allow_complex
{
    static std::string name() { return "real"; };

    template<class T1, class ... Args>
    static auto eval(const T1& a1, Args&& ...) -> decltype(real(std::declval<T1>()))
    {
        return real(a1);
    };
};

struct Imag_func : Allow_complex
{
    static std::string name() { return "imag"; };

    template<class T1, class ... Args>
    static auto eval(const T1& a1, Args&& ...) -> decltype(imag(std::declval<T1>()))
    {
        return imag(a1);
    };
};

struct Conj_func : Allow_complex
{
    static std::string name() { return "conj"; };

    template<class T1, class ... Args>
    static auto eval(const T1& a1, Args&& ... ) -> decltype(conj(std::declval<T1>()))
    {
        return conj(a1);
    };
};

struct Abs_func : Allow_complex
{
    static std::string name() { return "abs"; };

    template<class T1, class ... Args>
    static auto eval(const T1& a1, Args&& ... args) -> decltype(abs(std::declval<T1>()))
    {
        return abs(a1, std::forward<Args>(args)...);
    };
};

struct Abs2_func : Allow_complex
{
    static std::string name() { return "abs2"; };

    template<class T1, class ... Args>
    static auto eval(const T1& a1, Args&& ... args) -> decltype(abs2(std::declval<T1>()))
    {
        return abs2(a1, std::forward<Args>(args)...);
    };
};

struct Arg_func : Allow_complex
{
    static std::string name() { return "arg"; };

    template<class T1, class ... Args>
    static auto eval(const T1& a1, Args&& ... args) -> decltype(arg(std::declval<T1>()))
    {
        return arg(a1, std::forward<Args>(args)...);
    };
};

struct Angle_func : Allow_complex
{
    static std::string name() { return "angle"; };

    template<class T1, class ... Args>
    static auto eval(const T1& a1, Args&& ... args) -> decltype(angle(std::declval<T1>()))
    {
        return angle(a1, std::forward<Args>(args)...);
    };
};

struct Inv_func : Allow_complex
{
    static std::string name() { return "inv"; };

    template<class T1, class ... Args>
    static auto eval(const T1& a1, Args&& ... args) -> decltype(inv(std::declval<T1>()))
    {
        return inv(a1, std::forward<Args>(args)...);
    };
};
struct Invs_func : Allow_complex
{
    static std::string name() { return "invs"; };

    template<class T1, class ... Args>
    static auto eval(const T1& a1, Args&& ... args) -> decltype(invs(std::declval<T1>()))
    {
        return invs(a1, std::forward<Args>(args)...);
    };
};
struct Op_not_func : Allow_complex
{
    static std::string name() { return "op_not"; };

    template<class T1, class ... Args>
    static auto eval(const T1& a1, Args&& ...) -> decltype(!(std::declval<T1>()))
    {
        return !(a1);
    };
    template<class ... Args>
    static auto eval(const dynamic::object& a1, Args&& ...) -> dynamic::object
    {
        return dynamic::object(OBool(operator!(a1)));
    };
    template<class T, class ... Args>
    static auto eval(const dynamic::object_type<T>& a1, Args&& ...) -> OBool
    {
        return OBool(operator!(a1));
    };
};
struct Op_bool_func : Allow_complex
{
    static std::string name() { return "op_bool"; };

    template<class T1, class ... Args>
    static auto eval(const T1& a1, Args&& ...) -> decltype((bool)(std::declval<T1>()))
    {
        return (bool)(a1);
    };
    template<class ... Args>
    static auto eval(const dynamic::object& a1, Args&& ...) -> dynamic::object
    {
        return dynamic::object(OBool((bool)(a1)));
    };
    template<class T, class ... Args>
    static auto eval(const dynamic::object_type<T>& a1, Args&& ...) -> OBool
    {
        return OBool((bool)(a1));
    };
};
struct Op_neg_func : Allow_complex
{
    static std::string name() { return "neg"; };

    template<class T1, class ... Args>
    static auto eval(const T1& a1, Args&& ... args) -> decltype(neg(std::declval<T1>()))
    {
        return neg(a1, std::forward<Args>(args)...);
    };
};
struct Op_true_func : Allow_complex
{
    static std::string name() { return "is_true"; };

    template<class T1, class ... Args>
    static auto eval(const T1& a1, Args&& ... args) -> decltype(is_true(std::declval<T1>()))
    {
        return is_true(a1, std::forward<Args>(args)...);
    };
};
struct Is_nan_func : Allow_complex
{
    static std::string name() { return "is_nan"; };

    template<class T1, class ... Args>
    static auto eval(const T1& a1, Args&& ... args) -> decltype(is_nan(std::declval<T1>()))
    {
        return is_nan(a1, std::forward<Args>(args)...);
    };
};
struct Is_inf_func : Allow_complex
{
    static std::string name() { return "is_inf"; };

    template<class T1, class ... Args>
    static auto eval(const T1& a1, Args&& ... args) -> decltype(is_inf(std::declval<T1>()))
    {
        return is_inf(a1, std::forward<Args>(args)...);
    };
};
struct Is_regular_func : Allow_complex
{
    static std::string name() { return "is_regular"; };

    template<class T1, class ... Args>
    static auto eval(const T1& a1, Args&& ... args) -> decltype(is_regular(std::declval<T1>()))
    {
        return is_regular(a1, std::forward<Args>(args)...);
    };
};
struct Is_normal_func : Allow_complex
{
    static std::string name() { return "is_normal"; };

    template<class T1, class ... Args>
    static auto eval(const T1& a1, Args&& ... args) -> decltype(is_normal(std::declval<T1>()))
    {
        return is_normal(a1, std::forward<Args>(args)...);
    };
};
struct Is_finite_func : Allow_complex
{
    static std::string name() { return "is_finite"; };

    template<class T1, class ... Args>
    static auto eval(const T1& a1, Args&& ... args) -> decltype(is_finite(std::declval<T1>()))
    {
        return is_finite(a1, std::forward<Args>(args)...);
    };
};
struct Is_int_func : Allow_complex
{
    static std::string name() { return "is_int"; };

    template<class T1, class ... Args>
    static auto eval(const T1& a1, Args&& ... args) -> decltype(is_int(std::declval<T1>()))
    {
        return is_int(a1, std::forward<Args>(args)...);
    };
};
struct Is_real_func : Allow_complex
{
    static std::string name() { return "is_real"; };

    template<class T1, class ... Args>
    static auto eval(const T1& a1, Args&& ... args) -> decltype(is_real(std::declval<T1>()))
    {
        return is_real(a1, std::forward<Args>(args)...);
    };
};
struct Is_zero_func : Allow_complex
{
    static std::string name() { return "is_zero"; };

    template<class T1, class ... Args>
    static auto eval(const T1& a1, Args&& ... args) -> decltype(is_zero(std::declval<T1>()))
    {
        return is_zero(a1, std::forward<Args>(args)...);
    };
};
struct Is_one_func : Allow_complex
{
    static std::string name() { return "is_one"; };

    template<class T1, class ... Args>
    static auto eval(const T1& a1, Args&& ... args) -> decltype(is_one(std::declval<T1>()))
    {
        return is_one(a1, std::forward<Args>(args)...);
    };
};
struct Eps_func : Allow_complex
{
    static std::string name() { return "eps"; };

    template<class T1, class ... Args>
    static auto eval(const T1& a1, Args&& ... args) -> decltype(eps(std::declval<T1>()))
    {
        return eps(a1, std::forward<Args>(args)...);
    };
};
struct Asin_func : Allow_complex
{
    static std::string name() { return "asin"; };

    template<class T1, class ... Args>
    static auto eval(const T1& a1, Args&& ... args) -> decltype(asin(std::declval<T1>()))
    {
        return asin(a1, std::forward<Args>(args)...);
    };
};
struct Asin_c_func : Allow_complex
{
    static std::string name() { return "asin_c"; };

    template<class T1, class ... Args>
    static auto eval(const T1& a1, Args&& ... args) -> decltype(asin_c(std::declval<T1>()))
    {
        return asin_c(a1, std::forward<Args>(args)...);
    };
};

struct Acos_func : Allow_complex
{
    static std::string name() { return "acos"; };

    template<class T1, class ... Args>
    static auto eval(const T1& a1, Args&& ... args) -> decltype(acos(std::declval<T1>()))
    {
        return acos(a1, std::forward<Args>(args)...);
    };
};
struct Acos_c_func : Allow_complex
{
    static std::string name() { return "acos_c"; };

    template<class T1, class ... Args>
    static auto eval(const T1& a1, Args&& ... args) -> decltype(acos_c(std::declval<T1>()))
    {
        return acos_c(a1, std::forward<Args>(args)...);
    };
};
struct Asinh_func : Allow_complex
{
    static std::string name() { return "asinh"; };

    template<class T1, class ... Args>
    static auto eval(const T1& a1, Args&& ... args) -> decltype(asinh(std::declval<T1>()))
    {
        return asinh(a1, std::forward<Args>(args)...);
    };
};

struct Acosh_func : Allow_complex
{
    static std::string name() { return "acosh"; };

    template<class T1, class ... Args>
    static auto eval(const T1& a1, Args&& ... args) -> decltype(acosh(std::declval<T1>()))
    {
        return acosh(a1, std::forward<Args>(args)...);
    };
};
struct Acosh_c_func : Allow_complex
{
    static std::string name() { return "acosh_c"; };

    template<class T1, class ... Args>
    static auto eval(const T1& a1, Args&& ... args) -> decltype(acosh_c(std::declval<T1>()))
    {
        return acosh_c(a1, std::forward<Args>(args)...);
    };
};

struct Atanh_func : Allow_complex
{
    static std::string name() { return "atanh"; };

    template<class T1, class ... Args>
    static auto eval(const T1& a1, Args&& ... args) -> decltype(atanh(std::declval<T1>()))
    {
        return atanh(a1, std::forward<Args>(args)...);
    };
};
struct Atanh_c_func : Allow_complex
{
    static std::string name() { return "atanh_c"; };

    template<class T1, class ... Args>
    static auto eval(const T1& a1, Args&& ... args) -> decltype(atanh_c(std::declval<T1>()))
    {
        return atanh_c(a1, std::forward<Args>(args)...);
    };
};

struct Atan_func : Allow_complex
{
    static std::string name() { return "atan"; };

    template<class T1, class ... Args>
    static auto eval(const T1& a1, Args&& ... args) -> decltype(atan(std::declval<T1>()))
    {
        return atan(a1, std::forward<Args>(args)...);
    };
};
struct Acot_func : Allow_complex
{
    static std::string name() { return "acot"; };

    template<class T1, class ... Args>
    static auto eval(const T1& a1, Args&& ... args) -> decltype(acot(std::declval<T1>()))
    {
        return acot(a1, std::forward<Args>(args)...);
    };
};
struct Asec_func : Allow_complex
{
    static std::string name() { return "asec"; };

    template<class T1, class ... Args>
    static auto eval(const T1& a1, Args&& ... args) -> decltype(asec(std::declval<T1>()))
    {
        return asec(a1, std::forward<Args>(args)...);
    };
};
struct Asec_c_func : Allow_complex
{
    static std::string name() { return "asec_c"; };

    template<class T1, class ... Args>
    static auto eval(const T1& a1, Args&& ... args) -> decltype(asec_c(std::declval<T1>()))
    {
        return asec_c(a1, std::forward<Args>(args)...);
    };
};
struct Acsc_func : Allow_complex
{
    static std::string name() { return "acsc"; };

    template<class T1, class ... Args>
    static auto eval(const T1& a1, Args&& ... args) -> decltype(acsc(std::declval<T1>()))
    {
        return acsc(a1, std::forward<Args>(args)...);
    };
};
struct Acsc_c_func : Allow_complex
{
    static std::string name() { return "acsc_c"; };

    template<class T1, class ... Args>
    static auto eval(const T1& a1, Args&& ... args) -> decltype(acsc_c(std::declval<T1>()))
    {
        return acsc_c(a1, std::forward<Args>(args)...);
    };
};

struct Acoth_func : Allow_complex
{
    static std::string name() { return "acoth"; };

    template<class T1, class ... Args>
    static auto eval(const T1& a1, Args&& ... args) -> decltype(acoth(std::declval<T1>()))
    {
        return acoth(a1, std::forward<Args>(args)...);
    };
};
struct Acoth_c_func : Allow_complex
{
    static std::string name() { return "acoth_c"; };

    template<class T1, class ... Args>
    static auto eval(const T1& a1, Args&& ... args) -> decltype(acoth_c(std::declval<T1>()))
    {
        return acoth_c(a1, std::forward<Args>(args)...);
    };
};
struct Asech_func : Allow_complex
{
    static std::string name() { return "asech"; };

    template<class T1, class ... Args>
    static auto eval(const T1& a1, Args&& ... args) -> decltype(asech(std::declval<T1>()))
    {
        return asech(a1, std::forward<Args>(args)...);
    };
};
struct Asech_c_func : Allow_complex
{
    static std::string name() { return "asech_c"; };

    template<class T1, class ... Args>
    static auto eval(const T1& a1, Args&& ... args) -> decltype(asech_c(std::declval<T1>()))
    {
        return asech_c(a1, std::forward<Args>(args)...);
    };
};
struct Acsch_func : Allow_complex
{
    static std::string name() { return "acsch"; };

    template<class T1, class ... Args>
    static auto eval(const T1& a1, Args&& ... args) -> decltype(acsch(std::declval<T1>()))
    {
        return acsch(a1, std::forward<Args>(args)...);
    };
};

struct Sqrt_func : Allow_complex
{
    static std::string name() { return "sqrt"; };

    template<class T1, class ... Args>
    static auto eval(const T1& a1, Args&& ... args) -> decltype(sqrt(std::declval<T1>()))
    {
        return sqrt(a1, std::forward<Args>(args)...);
    };
};

struct Sqrt1pm1_func : Allow_complex
{
    static std::string name() { return "sqrt1pm1"; };

    template<class T1, class ... Args>
    static auto eval(const T1& a1, Args&& ... args) -> decltype(sqrt1pm1(std::declval<T1>()))
    {
        return sqrt1pm1(a1, std::forward<Args>(args)...);
    };
};

struct Cbrt_func : Notcomplex_func
{
    using Notcomplex_func::eval;

    static std::string name() { return "cbrt"; };

    template<class T1, class ... Args>
    static auto eval(const T1& a1, Args&& ... args) -> decltype(cbrt(std::declval<T1>()))
    {
        return cbrt(a1, std::forward<Args>(args)...);
    };
};

struct Sqrt_c_func : Allow_complex
{
    static std::string name() { return "sqrt_c"; };

    template<class T1, class ... Args>
    static auto eval(const T1& a1, Args&& ... args) -> decltype(sqrt_c(std::declval<T1>()))
    {
        return sqrt_c(a1, std::forward<Args>(args)...);
    };
};

struct Plus_func : Allow_complex
{
    static std::string name() { return "plus"; };

    template<class T1, class T2, class ... Args>
    static auto eval(const T1& a1, const T2& a2, Args&& ... args) 
        -> decltype(plus(a1, a2))
    {
        return plus(a1, a2, std::forward<Args>(args)...);
    };
};

struct Minus_func : Allow_complex
{
    static std::string name() { return "minus"; };

    template<class T1, class T2, class ... Args>
    static auto eval(const T1& a1, const T2& a2, Args&& ... args) 
        -> decltype(minus(a1, a2))
    {
        return minus(a1, a2, std::forward<Args>(args)...);
    };
};

struct Mult_func : Allow_complex
{
    static std::string name() { return "mul"; };

    template<class T1, class T2, class ... Args>
    static auto eval(const T1& a1, const T2& a2, Args&& ... args) 
        -> decltype(mul(a1, a2))
    {
        return mul(a1, a2, std::forward<Args>(args)...);
    };
};

struct MMult_func : Allow_complex
{
    static std::string name() { return "mmul"; };

    template<class T1, class T2, class ... Args>
    static auto eval(const T1& a1, const T2& a2, Args&& ... args) 
        -> decltype(mmul(a1, a2))
    {
        return mmul(a1, a2, std::forward<Args>(args)...);
    };
};

struct Div_func : Allow_complex
{
    static std::string name() { return "div"; };

    template<class T1, class T2, class ... Args>
    static auto eval(const T1& a1, const T2& a2, Args&& ... args) 
        -> decltype(div(a1, a2))
    {
        return div(a1, a2, std::forward<Args>(args)...);
    };
};

struct Idiv_func : Allow_complex
{
    static std::string name() { return "idiv"; };

    template<class T1, class T2, class ... Args>
    static auto eval(const T1& a1, const T2& a2, Args&& ... args) 
        -> decltype(idiv(a1, a2))
    {
        return idiv(a1, a2, std::forward<Args>(args)...);
    };
};

struct Div_0_func : Allow_complex
{
    static std::string name() { return "div_0"; };

    template<class T1, class T2, class ... Args>
    static auto eval(const T1& a1, const T2& a2, Args&& ... args) 
        -> decltype(div_0(a1, a2))
    {
        return div_0(a1, a2, std::forward<Args>(args)...);
    };
};

struct Div_1_func : Allow_complex
{
    static std::string name() { return "div_1"; };

    template<class T1, class T2, class ... Args>
    static auto eval(const T1& a1, const T2& a2, Args&& ... args) 
        -> decltype(div_1(a1, a2))
    {
        return div_1(a1, a2, std::forward<Args>(args)...);
    };
};


struct Hypot_func : Allow_complex
{
    static std::string name() { return "hypot"; };

    template<class T1, class T2, class ... Args>
    static auto eval(const T1& a1, const T2& a2, Args&& ... args) 
        -> decltype(hypot(a1, a2))
    {
        return hypot(a1, a2, std::forward<Args>(args)...);
    };
};

struct Pow_func : Allow_complex
{
    static std::string name() { return "pow"; };

    template<class T1, class T2, class ... Args>
    static auto eval(const T1& a1, const T2& a2, Args&& ... args) 
        -> decltype(pow(a1, a2))
    {
        return pow(a1, a2, std::forward<Args>(args)...);
    };
};

struct Pow_c_func : Allow_complex
{
    static std::string name() { return "pow_c"; };

    template<class T1, class T2, class ... Args>
    static auto eval(const T1& a1, const T2& a2, Args&& ... args) 
        -> decltype(pow_c(a1, a2))
    {
        return pow_c(a1, a2, std::forward<Args>(args)...);
    };
};

struct Pow_rc_func : Allow_complex
{
    static std::string name() { return "pow rc"; };

    template<class T1, class T2, class ... Args>
    static auto eval(const T1& a1, const T2& a2, Args&& ... args) 
        -> decltype(pow(real(a1), a2))
    {
        return pow(real(a1), a2, std::forward<Args>(args)...);
    };
};

struct Pow_cr_func : Allow_complex
{
    static std::string name() { return "pow cr"; };

    template<class T1, class T2, class ... Args>
    static auto eval(const T1& a1, const T2& a2, Args&& ... args) 
        -> decltype(pow(a1, real(a2)))
    {
        return pow(a1, real(a2), std::forward<Args>(args)...);
    };
};

struct Atan2_func : Notcomplex_func
{
    static std::string name() { return "atan2"; };

    template<class T1, class T2, class ... Args>
    static auto eval(const T1& a1, const T2& a2, Args&& ... args) 
        -> decltype(atan2(a1, a2))
    {
        return atan2(a1, a2, std::forward<Args>(args)...);
    };
};

struct Mod_func : Notcomplex_func
{
    static std::string name() { return "mod"; };

    template<class T1, class T2, class ... Args>
    static auto eval(const T1& a1, const T2& a2, Args&& ... args) 
        -> decltype(mod(a1, a2))
    {
        return mod(a1, a2, std::forward<Args>(args)...);
    };
};

struct Rem_func : Notcomplex_func
{
    static std::string name() { return "rem"; };

    template<class T1, class T2, class ... Args>
    static auto eval(const T1& a1, const T2& a2, Args&& ... args) 
        -> decltype(rem(a1, a2))
    {
        return rem(a1, a2, std::forward<Args>(args)...);
    };
};
struct Ldexp3p_func_t : Allow_complex
{
    static std::string name() { return "ldexp +3"; };

    template<class T1, class ... Args>
    static auto eval(const T1& a1, Args&& ... args) 
        -> decltype(ldexp(a1, Integer()))
    {
        return ldexp(a1, 3);
    };
};

struct Ldexp3m_func_t : Allow_complex
{
    static std::string name() { return "ldexp -3"; };

    template<class T1, class ... Args>
    static auto eval(const T1& a1, Args&& ... args) 
        -> decltype(ldexp(a1, Integer()))
    {
        return ldexp(a1, -3);
    };
};

struct Scalbn3p_func_t : Allow_complex
{
    static std::string name() { return "scalbn +3"; };

    template<class T1, class ... Args>
    static auto eval(const T1& a1, Args&& ... args) 
        -> decltype(scalbn(a1, Integer()))
    {
        return scalbn(a1, 3);
    };
};

struct Scalbn3m_func_t : Allow_complex
{
    static std::string name() { return "scalbn -3"; };

    template<class T1, class ... Args>
    static auto eval(const T1& a1, Args&& ... args) 
        -> decltype(scalbn(a1, Integer()))
    {
        return scalbn(a1, -3);
    };
};

struct Frexp_func_t : Notcomplex_func
{
    using Notcomplex_func::eval;

    static std::string name() { return "frexp"; };

    template<class T1, class ... Args>
    static auto eval(const T1& a1, Args&& ... args) 
        -> decltype(frexp(a1, std::declval<Integer&>()))
    {
        Integer exp;
        return frexp(a1, exp);
    };
};

template<class T>
struct modf_return_type
{
    using type0 = typename matcl::details::real_type<T>::type;
    using type  = typename matcl::details::unify_types<type0, Float>::type;
};
template<class T>
struct modf_return_type<dynamic::object_type<T>>
{
    using type0 = typename matcl::details::real_type<T>::type;
    using type1 = typename matcl::details::unify_types<type0, Float>::type;
    using type  = dynamic::object_type<type1>;
};
template<>
struct modf_return_type<dynamic::object>
{
    using type  = dynamic::object;
};
struct Modf_frac_func_t : Notcomplex_func
{
    using Notcomplex_func::eval;

    static std::string name() { return "modf frac"; };

    template<class T1, class ... Args>
    static auto eval(const T1& a1, Args&& ... args) 
        -> typename modf_return_type<T1>::type
    {
        using Ret   = typename modf_return_type<T1>::type;
        Ret int_part;
        return modf(a1, int_part);
    };
};

struct Modf_int_func_t : Notcomplex_func
{
    using Notcomplex_func::eval;

    static std::string name() { return "modf int"; };

    template<class T1, class ... Args>
    static auto eval(const T1& a1, Args&& ... args) 
        -> typename modf_return_type<T1>::type
    {
        using Ret   = typename modf_return_type<T1>::type;
        Ret int_part;
        modf(a1, int_part);

        return int_part;
    };
};

struct Logb_func : Allow_complex
{
    static std::string name() { return "logb"; };

    template<class T1, class ... Args>
    static auto eval(const T1& a1, Args&& ... args) 
        -> decltype(logb(a1))
    {
        return logb(a1);
    };
};

struct Ilogb_func : Allow_complex
{
    static std::string name() { return "ilogb"; };

    template<class T1, class ... Args>
    static auto eval(const T1& a1, Args&& ... args) 
        -> decltype(ilogb(a1))
    {
        return ilogb(a1);
    };
};

struct Eeq_func : Allow_complex
{
    static std::string name() { return "eeq"; };

    template<class T1, class T2, class ... Args>
    static auto eval(const T1& a1, const T2& a2, Args&& ... args) 
        -> decltype(eeq(a1, a2))
    {
        return eeq(a1, a2, std::forward<Args>(args)...);
    };
};

struct Neq_func : Allow_complex
{
    static std::string name() { return "neq"; };

    template<class T1, class T2, class ... Args>
    static auto eval(const T1& a1, const T2& a2, Args&& ... args) 
        -> decltype(neq(a1, a2))
    {
        return neq(a1, a2, std::forward<Args>(args)...);
    };
};

struct Leq_func : Allow_complex
{
    static std::string name() { return "leq"; };

    template<class T1, class T2, class ... Args>
    static auto eval(const T1& a1, const T2& a2, Args&& ... args) 
        -> decltype(leq(a1, a2))
    {
        return leq(a1, a2, std::forward<Args>(args)...);
    };
};

struct Geq_func : Allow_complex
{
    static std::string name() { return "geq"; };

    template<class T1, class T2, class ... Args>
    static auto eval(const T1& a1, const T2& a2, Args&& ... args) 
        -> decltype(geq(a1, a2))
    {
        return geq(a1, a2, std::forward<Args>(args)...);
    };
};

struct Gt_func : Allow_complex
{
    static std::string name() { return "gt"; };

    template<class T1, class T2, class ... Args>
    static auto eval(const T1& a1, const T2& a2, Args&& ... args) 
        -> decltype(gt(a1, a2))
    {
        return gt(a1, a2, std::forward<Args>(args)...);
    };
};

struct Lt_func : Allow_complex
{
    static std::string name() { return "lt"; };

    template<class T1, class T2, class ... Args>
    static auto eval(const T1& a1, const T2& a2, Args&& ... args) 
        -> decltype(lt(a1, a2))
    {
        return lt(a1, a2, std::forward<Args>(args)...);
    };
};

struct Min_func : Allow_complex
{
    static std::string name() { return "min"; };

    template<class T1, class T2, class ... Args>
    static auto eval(const T1& a1, const T2& a2, Args&& ... args) 
        -> decltype(min(a1, a2))
    {
        return min(a1, a2, std::forward<Args>(args)...);
    };
};

struct Max_func : Allow_complex
{
    static std::string name() { return "max"; };

    template<class T1, class T2, class ... Args>
    static auto eval(const T1& a1, const T2& a2, Args&& ... args) 
        -> decltype(max(a1, a2))
    {
        return max(a1, a2, std::forward<Args>(args)...);
    };
};

struct Elem_and_func : Allow_complex
{
    static std::string name() { return "elem_and"; };

    template<class T1, class T2, class ... Args>
    static auto eval(const T1& a1, const T2& a2, Args&& ... args) 
        -> decltype(elem_and(a1, a2))
    {
        return elem_and(a1, a2, std::forward<Args>(args)...);
    };
};

struct Elem_or_func : Allow_complex
{
    static std::string name() { return "elem_or"; };

    template<class T1, class T2, class ... Args>
    static auto eval(const T1& a1, const T2& a2, Args&& ... args) 
        -> decltype(elem_or(a1, a2))
    {
        return elem_or(a1, a2, std::forward<Args>(args)...);
    };
};

struct Elem_xor_func : Allow_complex
{
    static std::string name() { return "elem_xor"; };

    template<class T1, class T2, class ... Args>
    static auto eval(const T1& a1, const T2& a2, Args&& ... args) 
        -> decltype(elem_xor(a1, a2))
    {
        return elem_xor(a1, a2, std::forward<Args>(args)...);
    };
};

struct Op_or_func : Allow_complex
{
    static std::string name() { return "op_or"; };

    template<class T1, class T2, class ... Args>
    static auto eval(const T1& a1, const T2& a2, Args&& ... args) 
        -> decltype(op_or(a1, a2))
    {
        return op_or(a1, a2, std::forward<Args>(args)...);
    };

    template<class ... Args>
    static auto eval(const dynamic::object& a1, const dynamic::object& a2, Args&& ...) 
        -> dynamic::object
    {
        return dynamic::object(OBool(op_or(a1,a2)));
    };
    template<class T1, class T2, class ... Args>
    static auto eval(const dynamic::object_type<T1>& a1, const dynamic::object_type<T2>& a2, Args&& ...) 
        -> OBool
    {
        return OBool(op_or(a1,a2));
    };
};

struct Op_xor_func : Allow_complex
{
    static std::string name() { return "op_xor"; };

    template<class T1, class T2, class ... Args>
    static auto eval(const T1& a1, const T2& a2, Args&& ... args) 
        -> decltype(op_xor(a1, a2))
    {
        return op_xor(a1, a2, std::forward<Args>(args)...);
    };

    template<class ... Args>
    static auto eval(const dynamic::object& a1, const dynamic::object& a2, Args&& ...) 
        -> dynamic::object
    {
        return dynamic::object(OBool(op_xor(a1,a2)));
    };
    template<class T1, class T2, class ... Args>
    static auto eval(const dynamic::object_type<T1>& a1, const dynamic::object_type<T2>& a2, Args&& ...) 
        -> OBool
    {
        return OBool(op_xor(a1,a2));
    };

};

struct Op_and_func : Allow_complex
{
    static std::string name() { return "op_and"; };

    template<class T1, class T2, class ... Args>
    static auto eval(const T1& a1, const T2& a2, Args&& ... args) 
        -> decltype(op_and(a1, a2))
    {
        return op_and(a1, a2, std::forward<Args>(args)...);
    };

    template<class ... Args>
    static auto eval(const dynamic::object& a1, const dynamic::object& a2, Args&& ...) 
        -> dynamic::object
    {
        return dynamic::object(OBool(op_and(a1,a2)));
    };
    template<class T1, class T2, class ... Args>
    static auto eval(const dynamic::object_type<T1>& a1, const dynamic::object_type<T2>& a2, Args&& ...) 
        -> OBool
    {
        return OBool(op_and(a1,a2));
    };

};

struct Eeq_nan_func : Allow_complex
{
    static std::string name() { return "eeq_nan"; };

    template<class T1, class T2, class ... Args>
    static auto eval(const T1& a1, const T2& a2, Args&& ... args) 
        -> decltype(eeq_nan(a1, a2))
    {
        return eeq_nan(a1, a2, std::forward<Args>(args)...);
    };
};

struct Neq_nan_func : Allow_complex
{
    static std::string name() { return "neq_nan"; };

    template<class T1, class T2, class ... Args>
    static auto eval(const T1& a1, const T2& a2, Args&& ... args) 
        -> decltype(neq_nan(a1, a2))
    {
        return neq_nan(a1, a2, std::forward<Args>(args)...);
    };
};

}};
