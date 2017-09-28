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
#include "matcl-dynamic/object.h"
#include "matcl-dynamic/details/object_data.h"
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
function_name functions::op_eeq::eval()
{
    static function_name f("op_eeq");
    return f;
};

function_name functions::op_neq::eval()
{
    static function_name f("op_neq");
    return f;
};

function_name functions::op_lt::eval()
{
    static function_name f("op_lt");
    return f;
};

function_name functions::op_leq::eval()
{
    static function_name f("op_leq");
    return f;
};

function_name functions::op_geq::eval()
{
    static function_name f("op_geq");
    return f;
};

function_name functions::op_gt::eval()
{
    static function_name f("op_gt");
    return f;
};

function_name functions::op_uminus::eval()
{
    static function_name f("op_uminus");
    return f;
};

function_name functions::op_plus::eval()
{
    static function_name f("op_plus");
    return f;
};

function_name functions::op_minus::eval()
{
    static function_name f("op_minus");
    return f;
};

function_name functions::op_mul::eval()
{
    static function_name f("op_mul");
    return f;
};

function_name functions::op_div::eval()
{
    static function_name f("op_div");
    return f;
};

function_name functions::idiv::eval()
{
    static function_name f("idiv");
    return f;
};

function_name functions::op_bool::eval()
{
    static function_name f("op_bool", validator_1arg_ret_bool());
    return f;
};

function_name functions::op_not::eval()
{
    static function_name f("op_not", validator_1arg_ret_bool());
    return f;
};

function_name functions::real::eval()
{
    static function_name f("real");
    return f;
};

function_name functions::imag::eval()
{
    static function_name f("imag");
    return f;
};

function_name functions::is_zero::eval()
{
    static function_name f("is_zero", validator_1arg_ret_bool());
    return f;
};

function_name functions::is_one::eval()
{
    static function_name f("is_one", validator_1arg_ret_bool());
    return f;
};

};};};

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