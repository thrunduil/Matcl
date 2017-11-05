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

#include "matcl-dynamic/matcl_function_names.h"
#include "matcl-dynamic/function.h"

namespace matcl { namespace dynamic { namespace functions
{

struct validators
{
    static bool is_ret_bool(function f)
    {
        return f.return_type() == predefined::type_bool();
    };

    static bool is_ret_int(function f)
    {
        return f.return_type() == predefined::type_int();
    };

    static bool check_ret_bool(function f)
    {
        return is_ret_bool(f);
    };

    static void print_ret_bool(std::ostream& os)
    {
        os  << "function returns "
            << predefined::type_bool().to_string();;
    };

    static bool check_ret_int(function f)
    {
        return is_ret_int(f);
    };

    static void print_ret_int(std::ostream& os)
    {
        os  << "function returns "
            << predefined::type_int().to_string();;
    };

    static bool check_sec_ref(function f)
    {
        if (f.number_arguments() < 2)
            return false;

        Type t = f.argument_type(1);
        
        if (Type::is_reference(t) == false)
            return false;

        return true;
    };

    static void print_sec_ref(std::ostream& os)
    {
        os << "function takes at least two arguments and the second argument must be a reference";
    };

    static bool check_sec_int(function f)
    {
        if (f.number_arguments() < 2)
            return false;

        Type t = f.argument_type(1);
        
        if (t != predefined::type_int())
            return false;

        return true;
    };

    static void print_sec_int(std::ostream& os)
    {
        os << "function takes at least two arguments and the second argument must be integer";
    };

};

static function_validator validator_frexp()
{
    return function_validator(&validators::check_sec_ref, 
                              &validators::print_sec_ref);
};

static function_validator validator_modf()
{
    return function_validator(&validators::check_sec_ref, 
                              &validators::print_sec_ref);
};

static function_validator validator_ret_int()
{
    return function_validator(&validators::check_ret_int,
                              &validators::print_ret_int);
};

static function_validator validator_ret_bool()
{
    return function_validator(&validators::check_ret_bool,
                              &validators::print_ret_bool);
};

const function_name& functions::inv::eval()
{
    static function_name f("inv");
    return f;
};

const function_name& functions::invs::eval()
{
    static function_name f("invs");
    return f;
};

const function_name& functions::pow::eval()
{
    static function_name f("pow");
    return f;
};
const function_name& functions::pow_c::eval()
{
    static function_name f("pow_c");
    return f;
};

const function_name& functions::arg::eval()
{
    static function_name f("arg");
    return f;
};

const function_name& functions::conj::eval()
{
    static function_name f("conj");
    return f;
};

const function_name& functions::abs::eval()
{
    static function_name f("abs");
    return f;
};

const function_name& functions::abs2::eval()
{
    static function_name f("abs2");
    return f;
};

const function_name& functions::op_neg::eval()
{
    static function_name f("op_neg");
    return f;
};

const function_name& functions::op_true::eval()
{
    static function_name f("op_true");
    return f;
};

const function_name& functions::is_nan::eval()
{
    static function_name f("is_nan");
    return f;
};

const function_name& functions::is_finite::eval()
{
    static function_name f("is_finite");
    return f;
};

const function_name& functions::is_inf::eval()
{
    static function_name f("is_inf");
    return f;
};

const function_name& functions::is_regular::eval()
{
    static function_name f("is_regular");
    return f;
};

const function_name& functions::is_normal::eval()
{
    static function_name f("is_normal");
    return f;
};

const function_name& functions::is_int::eval()
{
    static function_name f("is_int");
    return f;
};

const function_name& functions::is_real::eval()
{
    static function_name f("is_real");
    return f;
};

const function_name& functions::fpclassify::eval()
{
    static function_name f("fpclassify", validator_ret_int());
    return f;
};

const function_name& functions::fma::eval()
{
    static function_name f("fma");
    return f;
};

const function_name& functions::fms::eval()
{
    static function_name f("fms");
    return f;
};

const function_name& functions::dot2_ac::eval()
{
    static function_name f("dot2_ac");
    return f;
};

const function_name& functions::elem_mul::eval()
{
    static function_name f("elem_mul");
    return f;
};

const function_name& functions::eeq_nan::eval()
{
    static function_name f("eeq_nan");
    return f;
};

const function_name& functions::neq_nan::eval()
{
    static function_name f("neq_nan");
    return f;
};

const function_name& functions::div_0::eval()
{
    static function_name f("div_0");
    return f;
};

const function_name& functions::div_1::eval()
{
    static function_name f("div_1");
    return f;
};

const function_name& functions::atan2::eval()
{
    static function_name f("atan2");
    return f;
};

const function_name& functions::mod::eval()
{
    static function_name f("mod");
    return f;
};

const function_name& functions::rem::eval()
{
    static function_name f("rem");
    return f;
};

const function_name& functions::hypot::eval()
{
    static function_name f("hypot");
    return f;
};

const function_name& functions::copysign::eval()
{
    static function_name f("copysign");
    return f;
};

const function_name& functions::nextafter::eval()
{
    static function_name f("nextafter");
    return f;
};

const function_name& functions::float_distance::eval()
{
    static function_name f("float_distance");
    return f;
};

const function_name& functions::nextbelow::eval()
{
    static function_name f("nextbelow");
    return f;
};

const function_name& functions::nextabove::eval()
{
    static function_name f("nextabove");
    return f;
};

const function_name& functions::signbit::eval()
{
    static function_name f("signbit");
    return f;
};

const function_name& functions::sign::eval()
{
    static function_name f("sign");
    return f;
};

const function_name& functions::isign::eval()
{
    static function_name f("isign");
    return f;
};

const function_name& functions::ldexp::eval()
{
    static function_name f("ldexp");
    return f;
};

const function_name& functions::frexp::eval()
{
    static function_name f("frexp", validator_frexp());
    return f;
};

const function_name& functions::modf::eval()
{
    static function_name f("modf", validator_modf());
    return f;
};

const function_name& functions::fdim::eval()
{
    static function_name f("fdim");
    return f;
};

const function_name& functions::logb::eval()
{
    static function_name f("logb");
    return f;
};

const function_name& functions::ilogb::eval()
{
    static function_name f("ilogb");
    return f;
};

const function_name& functions::eps::eval()
{
    static function_name f("eps");
    return f;
};

const function_name& functions::sqrt::eval()
{
    static function_name f("sqrt");
    return f;
};

const function_name& functions::cbrt::eval()
{
    static function_name f("cbrt");
    return f;
};

const function_name& functions::sqrt_c::eval()
{
    static function_name f("sqrt_c");
    return f;
};

const function_name& functions::exp::eval()
{
    static function_name f("exp");
    return f;
};

const function_name& functions::expm1::eval()
{
    static function_name f("expm1");
    return f;
};

const function_name& functions::expi::eval()
{
    static function_name f("expi");
    return f;
};


const function_name& functions::exp2::eval()
{
    static function_name f("exp2");
    return f;
};

const function_name& functions::exp10::eval()
{
    static function_name f("exp10");
    return f;
};

const function_name& functions::log::eval()
{
    static function_name f("log");
    return f;
};

const function_name& functions::log_c::eval()
{
    static function_name f("log_c");
    return f;
};

const function_name& functions::log1p::eval()
{
    static function_name f("log1p");
    return f;
};

const function_name& functions::log1p_c::eval()
{
    static function_name f("log1p_c");
    return f;
};

const function_name& functions::log2::eval()
{
    static function_name f("log2");
    return f;
};

const function_name& functions::log2_c::eval()
{
    static function_name f("log2_c");
    return f;
};

const function_name& functions::log10::eval()
{
    static function_name f("log10");
    return f;
};

const function_name& functions::log10_c::eval()
{
    static function_name f("log10_c");
    return f;
};

const function_name& functions::sin::eval()
{
    static function_name f("sin");
    return f;
};

const function_name& functions::cos::eval()
{
    static function_name f("cos");
    return f;
};

const function_name& functions::tan::eval()
{
    static function_name f("tan");
    return f;
};

const function_name& functions::cot::eval()
{
    static function_name f("cot");
    return f;
};

const function_name& functions::sec::eval()
{
    static function_name f("sec");
    return f;
};

const function_name& functions::csc::eval()
{
    static function_name f("csc");
    return f;
};

const function_name& functions::sinh::eval()
{
    static function_name f("sinh");
    return f;
};

const function_name& functions::cosh::eval()
{
    static function_name f("cosh");
    return f;
};

const function_name& functions::tanh::eval()
{
    static function_name f("tanh");
    return f;
};

const function_name& functions::coth::eval()
{
    static function_name f("coth");
    return f;
};

const function_name& functions::sech::eval()
{
    static function_name f("sech");
    return f;
};

const function_name& functions::csch::eval()
{
    static function_name f("csch");
    return f;
};

const function_name& functions::asin::eval()
{
    static function_name f("asin");
    return f;
};

const function_name& functions::asin_c::eval()
{
    static function_name f("asin_c");
    return f;
};

const function_name& functions::acos::eval()
{
    static function_name f("acos");
    return f;
};

const function_name& functions::acos_c::eval()
{
    static function_name f("acos_c");
    return f;
};

const function_name& functions::atan::eval()
{
    static function_name f("atan");
    return f;
};

const function_name& functions::acot::eval()
{
    static function_name f("acot");
    return f;
};

const function_name& functions::asec::eval()
{
    static function_name f("asec");
    return f;
};

const function_name& functions::asec_c::eval()
{
    static function_name f("asec_c");
    return f;
};

const function_name& functions::acsc::eval()
{
    static function_name f("acsc");
    return f;
};

const function_name& functions::acsc_c::eval()
{
    static function_name f("acsc_c");
    return f;
};

const function_name& functions::asinh::eval()
{
    static function_name f("asinh");
    return f;
};

const function_name& functions::acosh::eval()
{
    static function_name f("acosh");
    return f;
};

const function_name& functions::acosh_c::eval()
{
    static function_name f("acosh_c");
    return f;
};

const function_name& functions::atanh::eval()
{
    static function_name f("atanh");
    return f;
};

const function_name& functions::atanh_c::eval()
{
    static function_name f("atanh_c");
    return f;
};

const function_name& functions::acoth::eval()
{
    static function_name f("acoth");
    return f;
};

const function_name& functions::acoth_c::eval()
{
    static function_name f("acoth_c");
    return f;
};

const function_name& functions::asech::eval()
{
    static function_name f("asech");
    return f;
};

const function_name& functions::asech_c::eval()
{
    static function_name f("asech_c");
    return f;
};

const function_name& functions::acsch::eval()
{
    static function_name f("acsch");
    return f;
};

const function_name& functions::floor::eval()
{
    static function_name f("floor");
    return f;
};

const function_name& functions::ceil::eval()
{
    static function_name f("ceil");
    return f;
};

const function_name& functions::round::eval()
{
    static function_name f("round");
    return f;
};

const function_name& functions::trunc::eval()
{
    static function_name f("trunc");
    return f;
};

const function_name& functions::ifloor::eval()
{
    static function_name f("ifloor");
    return f;
};

const function_name& functions::iceil::eval()
{
    static function_name f("iceil");
    return f;
};

const function_name& functions::iround::eval()
{
    static function_name f("iround");
    return f;
};

const function_name& functions::itrunc::eval()
{
    static function_name f("itrunc");
    return f;
};

const function_name& functions::op_and::eval()
{
    static function_name f("op_and", validator_ret_bool());
    return f;
};

const function_name& functions::op_or::eval()
{
    static function_name f("op_or", validator_ret_bool());
    return f;
};

const function_name& functions::op_xor::eval()
{
    static function_name f("op_xor", validator_ret_bool());
    return f;
};

const function_name& functions::elem_and::eval()
{
    static function_name f("elem_and");
    return f;
};

const function_name& functions::elem_or::eval()
{
    static function_name f("elem_or");
    return f;
};

const function_name& functions::elem_xor::eval()
{
    static function_name f("elem_xor");
    return f;
};

const function_name& functions::max::eval()
{
    static function_name f("max");
    return f;
};

const function_name& functions::min::eval()
{
    static function_name f("min");
    return f;
};

const function_name& functions::sqrt1pm1_c::eval()
{
    static function_name f("sqrt1pm1_c");
    return f;
};
const function_name& functions::sqrt1pm1::eval()
{
    static function_name f("sqrt1pm1");
    return f;
};
const function_name& functions::powm1::eval()
{
    static function_name f("powm1");
    return f;
};
const function_name& functions::rising_factorial::eval()
{
    static function_name f("rising_factorial");
    return f;
};
const function_name& functions::falling_factorial::eval()
{
    static function_name f("falling_factorial");
    return f;
};
const function_name& functions::bernoulli_b2n::eval()
{
    static function_name f("bernoulli_b2n");
    return f;
};
const function_name& functions::max_bernoulli_b2n::eval()
{
    static function_name f("max_bernoulli_b2n", validator_ret_int());
    return f;
};

const function_name& functions::factorial::eval()
{
    static function_name f("factorial");
    return f;
};

const function_name& functions::double_factorial::eval()
{
    static function_name f("double_factorial");
    return f;
};

const function_name& functions::binomial_coefficient::eval()
{
    static function_name f("binomial_coefficient");
    return f;
};

};};};
