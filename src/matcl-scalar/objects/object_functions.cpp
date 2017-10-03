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

#include "matcl-scalar/objects/object_functions.h"
#include "matcl-dynamic/matcl_dynamic.h"
#include "matcl-dynamic/matcl_function_names.h"
#include "matcl-scalar/details/scalfunc_helpers.h"

namespace matcl { namespace dynamic { namespace functions
{

struct validators
{
    static bool check_ret_unit(const function& f)
    {
        return f.return_type() == predefined::type_unit();
    };

    static void print_ret_unit(std::ostream& os)
    {
        os << "function must return " << predefined::type_unit().to_string();
    };

    static bool check_ret_string(const function& f)
    {
        return f.return_type() == predefined::type_string();
    };

    static void print_ret_string(std::ostream& os)
    {
        os << "function must return " << predefined::type_string().to_string();
    };
};

static function_validator validator_disp()
{
    return function_validator(&validators::check_ret_unit,
                              &validators::print_ret_unit);
};

static function_validator validator_to_string()
{
    return function_validator(&validators::check_ret_string,
                              &validators::print_ret_string);
};

function_name disp::eval()
{
    static function_name f("disp", validator_disp());
    return f;
}

function_name to_string::eval()
{
    static function_name f("to_string", validator_to_string());
    return f;
}

}}};

namespace matcl { namespace dynamic
{

object dynamic::op_neg(const object& a)
{
    using func  = functions::op_neg;
    return eval_function::eval(func::eval(), a);
}

object dynamic::op_true(const object& a)
{
    using func  = functions::op_true;
    return eval_function::eval(func::eval(), a);
}

object dynamic::inv(const object& a)
{
    using func  = functions::inv;
    return eval_function::eval(func::eval(), a);
}

object dynamic::invs(const object& a)
{
    using func  = functions::invs;
    return eval_function::eval(func::eval(), a);
}

object dynamic::pow(const object& a, const object& b)
{
    using func  = functions::pow;
    return eval_function::eval(func::eval(), a, b);
};
object dynamic::pow_c(const object& a, const object& b)
{
    using func  = functions::pow_c;
    return eval_function::eval(func::eval(), a, b);
};

object dynamic::div(const object& a, const object& b)
{
    return a / b;
};

object dynamic::div_0(const object& a, const object& b)
{
    using func  = functions::div_0;
    return eval_function::eval(func::eval(), a, b);
};

object dynamic::div_1(const object& a, const object& b)
{    
    using func  = functions::div_1;
    return eval_function::eval(func::eval(), a, b);
};

object dynamic::arg(const object& a)
{
    using func  = functions::arg;
    return eval_function::eval(func::eval(), a);
};

object dynamic::conj(const object& a)
{
    using func  = functions::conj;
    return eval_function::eval(func::eval(), a);
};

object dynamic::abs(const object& a)
{
    using func  = functions::abs;
    return eval_function::eval(func::eval(), a);
};

object dynamic::abs2(const object& a)
{
    using func  = functions::abs2;
    return eval_function::eval(func::eval(), a);
};

object dynamic::is_nan(const object& a)
{
    using func  = functions::is_nan;
    return eval_function::eval(func::eval(), a);
};

object dynamic::is_finite(const object& a)
{
    using func  = functions::is_finite;
    return eval_function::eval(func::eval(), a);
};

object dynamic::is_inf(const object& a)
{
    using func  = functions::is_inf;
    return eval_function::eval(func::eval(), a);
};

object dynamic::is_regular(const object& a)
{
    using func  = functions::is_regular;
    return eval_function::eval(func::eval(), a);
};

object dynamic::is_normal(const object& a)
{
    using func  = functions::is_normal;
    return eval_function::eval(func::eval(), a);
};

object dynamic::is_int(const object& a)
{
    using func  = functions::is_int;
    return eval_function::eval(func::eval(), a);
};

object dynamic::is_real(const object& a)
{
    using func  = functions::is_real;
    return eval_function::eval(func::eval(), a);
};

fp_type dynamic::fpclassify(const object& a)
{
    using func  = functions::fpclassify;
    object ret  = eval_function::eval(func::eval(), a);

    //it is ensured, that this function returns OInteger
    Integer v   = static_cast<const details::object_data<Integer>*>(ret.get_data())->get();
    return matcl::raw::details::int_to_fptype(v);
};

object dynamic::copysign(const object& x, const object& y)
{
    using func  = functions::copysign;
    object ret  = eval_function::eval(func::eval(), x, y);
    return ret;
};

object dynamic::fma(const object& a, const object& b, const object& c)
{
    using func  = functions::fma;
    return eval_function::eval(func::eval(), a, b, c);
}

object dynamic::fms(const object& a, const object& b, const object& c)
{
    using func  = functions::fms;
    return eval_function::eval(func::eval(), a, b, c);
}

object dynamic::dot2_ac(const object& a, const object& b, const object& c, const object& d)
{
    using func  = functions::dot2_ac;
    return eval_function::eval(func::eval(), a, b, c, d);
}

object dynamic::eeq_nan(const object& x, const object& y)
{
    using func  = functions::eeq_nan;
    object ret  = eval_function::eval(func::eval(), x, y);
    return ret;
};

object dynamic::neq_nan(const object& x, const object& y)
{
    using func  = functions::neq_nan;
    object ret  = eval_function::eval(func::eval(), x, y);
    return ret;
};

object dynamic::elem_mul(const object& x, const object& y)
{
    using func  = functions::elem_mul;
    object ret  = eval_function::eval(func::eval(), x, y);
    return ret;
};
object dynamic::hypot(const object& x, const object& y)
{
    using func  = functions::hypot;
    object ret  = eval_function::eval(func::eval(), x, y);
    return ret;
};
object dynamic::atan2(const object& x, const object& y)
{
    using func  = functions::atan2;
    object ret  = eval_function::eval(func::eval(), x, y);
    return ret;
};
object dynamic::mod(const object& x, const object& y)
{
    using func  = functions::mod;
    object ret  = eval_function::eval(func::eval(), x, y);
    return ret;
};
object dynamic::rem(const object& x, const object& y)
{
    using func  = functions::rem;
    object ret  = eval_function::eval(func::eval(), x, y);
    return ret;
};

object dynamic::nextafter(const object& x, const object& y)
{
    using func  = functions::nextafter;
    object ret  = eval_function::eval(func::eval(), x, y);
    return ret;
};

object dynamic::nextabove(const object& x)
{
    using func  = functions::nextabove;
    object ret  = eval_function::eval(func::eval(), x);
    return ret;
};

object dynamic::nextbelow(const object& x)
{
    using func  = functions::nextbelow;
    object ret  = eval_function::eval(func::eval(), x);
    return ret;
};

object dynamic::signbit(const object& x)
{
    using func  = functions::signbit;
    return eval_function::eval(func::eval(), x);
};

object dynamic::sign(const object& x)
{
    using func  = functions::sign;
    return eval_function::eval(func::eval(), x);
};

object dynamic::isign(const object& x)
{
    using func  = functions::isign;
    return eval_function::eval(func::eval(), x);
};

object dynamic::ldexp(const object& x, const object& exp)
{
    using func  = functions::ldexp;
    return eval_function::eval(func::eval(), x, exp);
};

object dynamic::scalbn(const object& x, const object& exp)
{
    using func  = functions::ldexp;
    return eval_function::eval(func::eval(), x, exp);
};

object dynamic::frexp(const object& x, object& exp)
{
    using func  = functions::frexp;
    return eval_function::eval(func::eval(), x, exp);
};

object dynamic::modf(const object& x, object& exp)
{
    using func  = functions::modf;
    return eval_function::eval(func::eval(), x, exp);
};

object dynamic::fdim(const object& x, const object& exp)
{
    using func  = functions::fdim;
    return eval_function::eval(func::eval(), x, exp);
};

object dynamic::logb(const object& x)
{
    using func  = functions::logb;
    return eval_function::eval(func::eval(), x);
};

object dynamic::ilogb(const object& x)
{
    using func  = functions::ilogb;
    return eval_function::eval(func::eval(), x);
};

object dynamic::eps(const object& x)
{
    using func  = functions::eps;
    object ret  = eval_function::eval(func::eval(), x);
    return ret;
};

object dynamic::sqrt(const object& x)
{
    using func  = functions::sqrt;
    object ret  = eval_function::eval(func::eval(), x);
    return ret;
};

object dynamic::cbrt(const object& x)
{
    using func  = functions::cbrt;
    object ret  = eval_function::eval(func::eval(), x);
    return ret;
};

object dynamic::sqrt_c(const object& x)
{
    using func  = functions::sqrt_c;
    object ret  = eval_function::eval(func::eval(), x);
    return ret;
};

object dynamic::exp(const object& x)
{
    using func  = functions::exp;
    object ret  = eval_function::eval(func::eval(), x);
    return ret;
};

object dynamic::expm1(const object& x)
{
    using func  = functions::expm1;
    object ret  = eval_function::eval(func::eval(), x);
    return ret;
};

object dynamic::expi(const object& x)
{
    using func  = functions::expi;
    object ret  = eval_function::eval(func::eval(), x);
    return ret;
};

object dynamic::exp2(const object& x)
{
    using func  = functions::exp2;
    object ret  = eval_function::eval(func::eval(), x);
    return ret;
};

object dynamic::exp10(const object& x)
{
    using func  = functions::exp10;
    object ret  = eval_function::eval(func::eval(), x);
    return ret;
};

object dynamic::log(const object& x)
{
    using func  = functions::log;
    object ret  = eval_function::eval(func::eval(), x);
    return ret;
};

object dynamic::log_c(const object& x)
{
    using func  = functions::log_c;
    object ret  = eval_function::eval(func::eval(), x);
    return ret;
};

object dynamic::log1p(const object& x)
{
    using func  = functions::log1p;
    object ret  = eval_function::eval(func::eval(), x);
    return ret;
};

object dynamic::log1p_c(const object& x)
{
    using func  = functions::log1p_c;
    object ret  = eval_function::eval(func::eval(), x);
    return ret;
};

object dynamic::log2(const object& x)
{
    using func  = functions::log2;
    object ret  = eval_function::eval(func::eval(), x);
    return ret;
};

object dynamic::log2_c(const object& x)
{
    using func  = functions::log2_c;
    object ret  = eval_function::eval(func::eval(), x);
    return ret;
};

object dynamic::log10(const object& x)
{
    using func  = functions::log10;
    object ret  = eval_function::eval(func::eval(), x);
    return ret;
};

object dynamic::log10_c(const object& x)
{
    using func  = functions::log10_c;
    object ret  = eval_function::eval(func::eval(), x);
    return ret;
};

object dynamic::sin(const object& x)
{
    using func  = functions::sin;
    object ret  = eval_function::eval(func::eval(), x);
    return ret;
};

object dynamic::cos(const object& x)
{
    using func  = functions::cos;
    object ret  = eval_function::eval(func::eval(), x);
    return ret;
};

object dynamic::tan(const object& x)
{
    using func  = functions::tan;
    object ret  = eval_function::eval(func::eval(), x);
    return ret;
};

object dynamic::cot(const object& x)
{
    using func  = functions::cot;
    object ret  = eval_function::eval(func::eval(), x);
    return ret;
};

object dynamic::sec(const object& x)
{
    using func  = functions::sec;
    object ret  = eval_function::eval(func::eval(), x);
    return ret;
};

object dynamic::csc(const object& x)
{
    using func  = functions::csc;
    object ret  = eval_function::eval(func::eval(), x);
    return ret;
};

object dynamic::sinh(const object& x)
{
    using func  = functions::sinh;
    object ret  = eval_function::eval(func::eval(), x);
    return ret;
};

object dynamic::cosh(const object& x)
{
    using func  = functions::cosh;
    object ret  = eval_function::eval(func::eval(), x);
    return ret;
};

object dynamic::tanh(const object& x)
{
    using func  = functions::tanh;
    object ret  = eval_function::eval(func::eval(), x);
    return ret;
};

object dynamic::coth(const object& x)
{
    using func  = functions::coth;
    object ret  = eval_function::eval(func::eval(), x);
    return ret;
};

object dynamic::sech(const object& x)
{
    using func  = functions::sech;
    object ret  = eval_function::eval(func::eval(), x);
    return ret;
};

object dynamic::csch(const object& x)
{
    using func  = functions::csch;
    object ret  = eval_function::eval(func::eval(), x);
    return ret;
};

object dynamic::asin(const object& x)
{
    using func  = functions::asin;
    object ret  = eval_function::eval(func::eval(), x);
    return ret;
};

object dynamic::asin_c(const object& x)
{
    using func  = functions::asin_c;
    object ret  = eval_function::eval(func::eval(), x);
    return ret;
};

object dynamic::acos(const object& x)
{
    using func  = functions::acos;
    object ret  = eval_function::eval(func::eval(), x);
    return ret;
};

object dynamic::acos_c(const object& x)
{
    using func  = functions::acos_c;
    object ret  = eval_function::eval(func::eval(), x);
    return ret;
};

object dynamic::atan(const object& x)
{
    using func  = functions::atan;
    object ret  = eval_function::eval(func::eval(), x);
    return ret;
};

object dynamic::acot(const object& x)
{
    using func  = functions::acot;
    object ret  = eval_function::eval(func::eval(), x);
    return ret;
};

object dynamic::asec(const object& x)
{
    using func  = functions::asec;
    object ret  = eval_function::eval(func::eval(), x);
    return ret;
};

object dynamic::asec_c(const object& x)
{
    using func  = functions::asec_c;
    object ret  = eval_function::eval(func::eval(), x);
    return ret;
};

object dynamic::acsc(const object& x)
{
    using func  = functions::acsc;
    object ret  = eval_function::eval(func::eval(), x);
    return ret;
};

object dynamic::acsc_c(const object& x)
{
    using func  = functions::acsc_c;
    object ret  = eval_function::eval(func::eval(), x);
    return ret;
};

object dynamic::asinh(const object& x)
{
    using func  = functions::asinh;
    object ret  = eval_function::eval(func::eval(), x);
    return ret;
};

object dynamic::acosh(const object& x)
{
    using func  = functions::acosh;
    object ret  = eval_function::eval(func::eval(), x);
    return ret;
};

object dynamic::acosh_c(const object& x)
{
    using func  = functions::acosh_c;
    object ret  = eval_function::eval(func::eval(), x);
    return ret;
};

object dynamic::atanh(const object& x)
{
    using func  = functions::atanh;
    object ret  = eval_function::eval(func::eval(), x);
    return ret;
};

object dynamic::atanh_c(const object& x)
{
    using func  = functions::atanh_c;
    object ret  = eval_function::eval(func::eval(), x);
    return ret;
};

object dynamic::acoth(const object& x)
{
    using func  = functions::acoth;
    object ret  = eval_function::eval(func::eval(), x);
    return ret;
};

object dynamic::acoth_c(const object& x)
{
    using func  = functions::acoth_c;
    object ret  = eval_function::eval(func::eval(), x);
    return ret;
};

object dynamic::asech(const object& x)
{
    using func  = functions::asech;
    object ret  = eval_function::eval(func::eval(), x);
    return ret;
};

object dynamic::asech_c(const object& x)
{
    using func  = functions::asech_c;
    object ret  = eval_function::eval(func::eval(), x);
    return ret;
};

object dynamic::acsch(const object& x)
{
    using func  = functions::acsch;
    object ret  = eval_function::eval(func::eval(), x);
    return ret;
};

object dynamic::floor(const object& x)
{
    using func  = functions::floor;
    object ret  = eval_function::eval(func::eval(), x);
    return ret;
};

object dynamic::ceil(const object& x)
{
    using func  = functions::ceil;
    object ret  = eval_function::eval(func::eval(), x);
    return ret;
};

object dynamic::round(const object& x)
{
    using func  = functions::round;
    object ret  = eval_function::eval(func::eval(), x);
    return ret;
};

object dynamic::trunc(const object& x)
{
    using func  = functions::trunc;
    object ret  = eval_function::eval(func::eval(), x);
    return ret;
};

object dynamic::ifloor(const object& x)
{
    using func  = functions::ifloor;
    object ret  = eval_function::eval(func::eval(), x);
    return ret;
};

object dynamic::iceil(const object& x)
{
    using func  = functions::iceil;
    object ret  = eval_function::eval(func::eval(), x);
    return ret;
};

object dynamic::iround(const object& x)
{
    using func  = functions::iround;
    object ret  = eval_function::eval(func::eval(), x);
    return ret;
};

object dynamic::itrunc(const object& x)
{
    using func  = functions::itrunc;
    object ret  = eval_function::eval(func::eval(), x);
    return ret;
};

object dynamic::sqrt1pm1_c(const object& x)
{
    dynamic::function_name fn = functions::sqrt1pm1_c::eval();
    return eval_function::eval(fn, x);
};

object dynamic::sqrt1pm1(const object& obj)
{
    dynamic::function_name fn = functions::sqrt1pm1::eval();
    return eval_function::eval(fn, obj);
};

object dynamic::powm1(const object& x, const object& y)
{
    dynamic::function_name fn = functions::powm1::eval();
    return eval_function::eval(fn, x, y);
};

object dynamic::elem_and(const object& x, const object& y)
{
    dynamic::function_name fn = functions::elem_and::eval();
    return eval_function::eval(fn, x, y);
};

object dynamic::elem_or(const object& x, const object& y)
{
    dynamic::function_name fn = functions::elem_or::eval();
    return eval_function::eval(fn, x, y);
};

object dynamic::elem_xor(const object& x, const object& y)
{
    dynamic::function_name fn = functions::elem_xor::eval();
    return eval_function::eval(fn, x, y);
};

bool dynamic::op_and(const object& x, const object& y)
{
    dynamic::function_name fn = functions::op_and::eval();
    object ret = eval_function::eval(fn, x, y);

    //it is ensured, that this function returns OBool
    return static_cast<const details::object_data<bool>*>(ret.get_data())->get();
};

bool dynamic::op_or(const object& x, const object& y)
{
    dynamic::function_name fn = functions::op_or::eval();
    object ret = eval_function::eval(fn, x, y);

    //it is ensured, that this function returns OBool
    return static_cast<const details::object_data<bool>*>(ret.get_data())->get();
};

bool dynamic::op_xor(const object& x, const object& y)
{
    dynamic::function_name fn = functions::op_xor::eval();
    object ret = eval_function::eval(fn, x, y);

    //it is ensured, that this function returns OBool
    return static_cast<const details::object_data<bool>*>(ret.get_data())->get();
};

object dynamic::max(const object& x, const object& y)
{
    dynamic::function_name fn = functions::max::eval();
    return eval_function::eval(fn, x, y);
};

object dynamic::min(const object& x, const object& y)
{
    dynamic::function_name fn = functions::min::eval();
    return eval_function::eval(fn, x, y);
};

object dynamic::rising_factorial(const object& x, Integer i)
{
    dynamic::function_name fn = functions::rising_factorial::eval();
    OInteger oi(i);
    return eval_function::eval(fn, x, object(oi));
};

object dynamic::falling_factorial(const object& x, Integer i)
{
    dynamic::function_name fn = functions::falling_factorial::eval();
    OInteger oi(i);
    return eval_function::eval(fn, x, object(oi));
};

object dynamic::bernoulli_b2n(Type t, Integer n)
{
    dynamic::function_name fn = functions::bernoulli_b2n::eval();
    OInteger oi(n);
    return eval_function_template{t}.eval(fn, object(oi));
}

Integer dynamic::max_bernoulli_b2n(Type t)
{
    dynamic::function_name fn = functions::max_bernoulli_b2n::eval();
    dynamic::object ret = eval_function_template{t}.eval(fn);

    //it is ensured, that this function returns OInteger
    return static_cast<const details::object_data<Integer>*>(ret.get_data())->get();
}

object dynamic::factorial(Type t, Integer i)
{
    dynamic::function_name fn = functions::factorial::eval();
    OInteger oi(i);
    return eval_function_template{t}.eval(fn, object(oi));
};

object dynamic::double_factorial(Type t, Integer i)
{
    dynamic::function_name fn = functions::double_factorial::eval();
    OInteger oi(i);
    return eval_function_template{t}.eval(fn, object(oi));
};

object dynamic::binomial_coefficient(Type t,  Integer n, Integer k)
{
    dynamic::function_name fn = functions::binomial_coefficient::eval();
    OInteger on(n);
    OInteger ok(k);
    return eval_function_template{t}.eval(fn, object(on), object(ok));
};

}};