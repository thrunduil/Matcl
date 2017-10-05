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

#include "matcl-dynamic/predefined_functions.h"
#include "matcl-dynamic/details/object.inl"
#include "matcl-dynamic/details/object_data.inl"
#include "matcl-core/details/stack_array.h"
#include "matcl-dynamic/initialization.h"
#include "type_table.h"

namespace matcl { namespace dynamic { namespace functions
{

//---------------------------------------------------------------
//                  validators
//---------------------------------------------------------------
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

    static bool check_2args_ret_bool(function f)
    {
        return check_2args(f) && is_ret_bool(f);
    };

    static void print_2args_ret_bool(std::ostream& os)
    {
        os << "function takes two arguments and returns " 
           << predefined::type_bool().to_string();
    };

    static bool check_2args(function f)
    {
        return f.number_arguments() == 2;
    };
    
    static void print_2args(std::ostream& os)
    {
        os << "function takes two arguments";
    };

    static bool check_1args_ret_bool(function f)
    {
        return check_1args(f) && is_ret_bool(f);
    };

    static void print_1args_ret_bool(std::ostream& os)
    {
        os  << "function takes one argument and returns "
            << predefined::type_bool().to_string();;
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

    static bool check_1args(function f)
    {
        return f.number_arguments() == 1;
    };
    
    static void print_1args(std::ostream& os)
    {
        os << "function takes one argument";
    };

    static bool check_2args_sec_ref(function f)
    {
        if (f.number_arguments() != 2)
            return false;

        Type t = f.argument_type(1);
        
        if (Type::is_reference(t) == false)
            return false;

        return true;
    };

    static void print_2args_sec_ref(std::ostream& os)
    {
        os << "function takes two arguments and the second argument must be a reference";
    };
};

static function_validator validator_1arg_ret_bool()
{
    return function_validator(&validators::check_1args_ret_bool,
                              &validators::print_1args_ret_bool);
};

//---------------------------------------------------------------
//                  functions
//---------------------------------------------------------------
enum class func_code : size_t
{
    op_eeq,     op_neq,     op_lt,      op_leq,     op_geq,
    op_gt,      op_uminus,  op_plus,    op_minus,   op_mul,
    op_div,     idiv,       op_bool,    op_not,     real,
    imag,       is_zero,    is_one,

    size
};

using function_name_pod = md::pod_type<function_name>;

// function names must be initialized very early; dynamic initialization
// phase is too late; we could do initialization using static variables
// in functions, but overhead in this case is very high.
// function names are initialized by matcl_dynamic initializer

static
function_name_pod func_name[(size_t)func_code::size];

static inline
const function_name* cast(const function_name_pod* arr)
{
    return reinterpret_cast<const function_name*>(arr);
};

function_name functions::op_eeq::eval()
{
    return cast(func_name)[(size_t)func_code::op_eeq];
};

function_name functions::op_neq::eval()
{
    return cast(func_name)[(size_t)func_code::op_neq];
};

function_name functions::op_lt::eval()
{
    return cast(func_name)[(size_t)func_code::op_lt];
};

function_name functions::op_leq::eval()
{
    return cast(func_name)[(size_t)func_code::op_leq];
};

function_name functions::op_geq::eval()
{
    return cast(func_name)[(size_t)func_code::op_geq];
};

function_name functions::op_gt::eval()
{
    return cast(func_name)[(size_t)func_code::op_gt];
};

function_name functions::op_uminus::eval()
{
    return cast(func_name)[(size_t)func_code::op_uminus];
};

function_name functions::op_plus::eval()
{
    return cast(func_name)[(size_t)func_code::op_plus];
};

function_name functions::op_minus::eval()
{
    return cast(func_name)[(size_t)func_code::op_minus];
};

function_name functions::op_mul::eval()
{
    return cast(func_name)[(size_t)func_code::op_mul];
};

function_name functions::op_div::eval()
{
    return cast(func_name)[(size_t)func_code::op_div];
};

function_name functions::idiv::eval()
{
    return cast(func_name)[(size_t)func_code::idiv];
};

function_name functions::op_bool::eval()
{
    return cast(func_name)[(size_t)func_code::op_bool];
};

function_name functions::op_not::eval()
{
    return cast(func_name)[(size_t)func_code::op_not];
};

function_name functions::real::eval()
{
    return cast(func_name)[(size_t)func_code::real];
};

function_name functions::imag::eval()
{
    return cast(func_name)[(size_t)func_code::imag];
};

function_name functions::is_zero::eval()
{
    return cast(func_name)[(size_t)func_code::is_zero];
};

function_name functions::is_one::eval()
{
    return cast(func_name)[(size_t)func_code::is_one];
};

};};};

namespace matcl { namespace dynamic
{

#define init_func(name)                                         \
new(functions::func_name +(size_t)functions::func_code::name)   \
    function_name(#name);

void details::initialize_funcions()
{
    init_func(op_eeq)
    init_func(op_neq)
    init_func(op_lt)
    init_func(op_leq)
    init_func(op_geq)
    init_func(op_gt)
    init_func(op_uminus)
    init_func(op_plus)
    init_func(op_minus)
    init_func(op_mul)
    init_func(op_div)
    init_func(idiv)
    init_func(real)
    init_func(imag)

    init_func(op_bool)
    init_func(op_not)
    init_func(is_zero)
    init_func(is_one)

    new (functions::func_name + (size_t)functions::func_code::op_bool) 
        function_name("op_bool", functions::validator_1arg_ret_bool());

    new (functions::func_name + (size_t)functions::func_code::op_not) 
        function_name("op_not", functions::validator_1arg_ret_bool());

    new (functions::func_name + (size_t)functions::func_code::is_zero) 
        function_name("is_zero", functions::validator_1arg_ret_bool());

    new (functions::func_name + (size_t)functions::func_code::is_one)
        function_name("is_one", functions::validator_1arg_ret_bool());
};

}}

namespace matcl { namespace dynamic
{

object dynamic::operator==(const object& a, const object& b)
{
    using func  = functions::op_eeq;
    return eval_function::eval(func::eval(), a,b);
};

object dynamic::operator!=(const object& a, const object& b)
{
    using func  = functions::op_neq;
    return eval_function::eval(func::eval(), a,b);
};

object dynamic::operator<(const object& a, const object& b)
{
    using func  = functions::op_lt;
    return eval_function::eval(func::eval(), a,b);
}

object dynamic::operator<=(const object& a, const object& b)
{
    using func  = functions::op_leq;
    return eval_function::eval(func::eval(), a,b);
};

object dynamic::operator>=(const object& a, const object& b)
{
    using func  = functions::op_geq;
    return eval_function::eval(func::eval(), a,b);
}

object dynamic::operator>(const object& a, const object& b)
{
    using func  = functions::op_gt;
    return eval_function::eval(func::eval(), a,b);
}

bool dynamic::operator!(const object& a)
{
    if (a.is_null() == true)
        return true;

    using func  = functions::op_not;
    object ret = eval_function::eval(func::eval(), a);

    //it is ensured, that this function returns OBool
    return static_cast<const details::object_data<bool>*>(ret.get_data())->get();
}

bool dynamic::cast_bool(const object& a)
{
    if (a.is_null() == true)
        return false;

    using func  = functions::op_bool;
    object ret = eval_function::eval(func::eval(), a);

    //it is ensured, that this function returns OBool
    return static_cast<const details::object_data<bool>*>(ret.get_data())->get();
}

object dynamic::operator-(const object& a)
{
    using func  = functions::op_uminus;
    return eval_function::eval(func::eval(), a);
}

object dynamic::operator+(const object& a, const object& b)
{
    using func  = functions::op_plus;
    return eval_function::eval(func::eval(), a,b);
};

object dynamic::operator-(const object& a, const object& b)
{
    using func  = functions::op_minus;
    return eval_function::eval(func::eval(), a,b);
};

object dynamic::operator*(const object& a, const object& b)
{
    using func  = functions::op_mul;
    return eval_function::eval(func::eval(), a,b);
};

object dynamic::operator/(const object& a, const object& b)
{
    using func  = functions::op_div;
    return eval_function::eval(func::eval(), a,b);
};

object dynamic::real(const object& a)
{
    using func  = functions::real;
    return eval_function::eval(func::eval(), a);
};

object dynamic::idiv(const object& a, const object& b)
{
    using func  = functions::idiv;
    return eval_function::eval(func::eval(), a, b);
};

object dynamic::imag(const object& a)
{
    using func  = functions::imag;
    return eval_function::eval(func::eval(), a);
};

bool dynamic::is_zero(const object& a)
{
    return a.is_zero();
};

bool dynamic::is_one(const object& a)
{
    return a.is_one();
};

}};