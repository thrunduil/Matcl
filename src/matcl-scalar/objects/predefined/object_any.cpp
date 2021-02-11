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

#include "matcl-dynamic/register_function.h"
#include "matcl-dynamic/object_type.h"
#include "matcl-dynamic/details/register_function_macro.h"
#include "matcl-scalar/matcl_scalar.h"
#include "matcl-scalar/IO/scalar_io.h"

namespace matcl { namespace dynamic
{

struct f_disp : register_function<f_disp, functions::disp>
{
	static Unit eval(const Any& v)
	{
		matcl::disp(v.get().get_stored());
        return Unit();
	};
};

struct f_to_string : register_function<f_to_string, functions::to_string>
{
	static String eval(const Any& v)
	{
		String ret(matcl::to_string(v.get().get_stored()));
        return ret;
	};
};

struct f_neg : register_function<f_neg, functions::op_neg>
{
	static OBool eval(const Any& v)
	{
		return OBool(!(v.get().get_stored()));
	};
};

struct f_true : register_function<f_true, functions::op_true>
{
	static OBool eval(const Any& v)
	{
		return OBool(cast_bool(v.get().get_stored()));
	};
};

struct f_mul : dynamic::register_function_template_return<f_mul, functions::elem_mul>
{
    static dynamic::Template eval(const dynamic::OType& t_ret, const dynamic::Template& x, 
                             const dynamic::Template& y) 
    {
        const object& o1    = x.get();
        const object& o2    = y.get();
        object ret          = o1 * o2;
        Type type_ret       = t_ret.get();
        return object(type_ret, ret);
    }
    static dynamic::Type eval_return(int n_template, const dynamic::Type templates[], int n_arg, 
                                     const dynamic::Type args[])
    {
        (void)n_template;
        (void)templates;

        Type t  = operations::return_type(functions::op_mul::eval(), n_arg, args);
        return t;
    };
};

struct f_div_0 : dynamic::register_function_template_return<f_div_0, functions::div_0>
{
    static dynamic::Template eval(const dynamic::OType& t_ret, const dynamic::Template& x, 
                             const dynamic::Template& y) 
    {
        const object& o1    = x.get();
        const object& o2    = y.get();

        Type type_ret       = t_ret.get();

        if (o1.is_zero() == true && o2.is_zero() == true)
            return object(type_ret);
    
        object ret          = o1 / o2;        
        return object(type_ret, ret);
    }
    static dynamic::Type eval_return(int n_template, const dynamic::Type templates[], int n_arg, 
                                     const dynamic::Type args[])
    {
        (void)n_template;
        (void)templates;

        Type t  = operations::return_type(functions::op_div::eval(), n_arg, args);
        return t;
    };
};

struct f_div_1 : dynamic::register_function_template_return<f_div_1, functions::div_1>
{
    static dynamic::Template eval(const dynamic::OType& t_ret, const dynamic::Template& x, 
                             const dynamic::Template& y) 
    {
        const object& o1    = x.get();
        const object& o2    = y.get();

        Type type_ret       = t_ret.get();

        if (o1.is_zero() == true && o2.is_zero() == true)
            return object::make_one(type_ret);

        object ret          = o1 / o2;        
        return object(type_ret, ret);
    }
    static dynamic::Type eval_return(int n_template, const dynamic::Type templates[], int n_arg, 
                                     const dynamic::Type args[])
    {
        (void)n_template;
        (void)templates;

        Type t      = operations::return_type(functions::op_div::eval(), n_arg, args);
        bool has1   = operations::has_one(t);
        return has1 ? t : Type();
    };
};

struct f_op_and : dynamic::register_function_template_return<f_op_and, functions::op_and>
{
    static dynamic::OBool eval(const dynamic::Template& x, const dynamic::Template& y) 
    {
        const object& o1    = x.get();
        const object& o2    = y.get();

        bool ret            = cast_bool(o1) && cast_bool(o2);        
        return OBool(ret);
    }
    static dynamic::Type eval_return(int, const dynamic::Type[], int, const dynamic::Type args[])
    {
        if (args[0] == Any::get_static_type() || args[1] == Any::get_static_type())
            return Type();
        else
            return OBool::get_static_type();
    };
};

struct f_op_or : dynamic::register_function_template_return<f_op_or, functions::op_or>
{
    static dynamic::OBool eval(const dynamic::Template& x, const dynamic::Template& y) 
    {
        const object& o1    = x.get();
        const object& o2    = y.get();

        bool ret            = cast_bool(o1) || cast_bool(o2);        
        return OBool(ret);
    }
    static dynamic::Type eval_return(int, const dynamic::Type[], int, const dynamic::Type args[])
    {
        if (args[0] == Any::get_static_type() || args[1] == Any::get_static_type())
            return Type();
        else
            return OBool::get_static_type();
    };
};

struct f_op_xor : dynamic::register_function_template_return<f_op_xor, functions::op_xor>
{
    static dynamic::OBool eval(const dynamic::Template& x, const dynamic::Template& y) 
    {
        const object& o1    = x.get();
        const object& o2    = y.get();

        bool ret            = cast_bool(o1) ^ cast_bool(o2);        
        return OBool(ret);
    }
    static dynamic::Type eval_return(int, const dynamic::Type[], int, const dynamic::Type args[])
    {
        if (args[0] == Any::get_static_type() || args[1] == Any::get_static_type())
            return Type();
        else
            return OBool::get_static_type();
    };
};

struct f_elem_and : dynamic::register_function_template_return<f_elem_and, functions::elem_and>
{
    static dynamic::OBool eval(const dynamic::Template& x, const dynamic::Template& y) 
    {
        const object& o1    = x.get();
        const object& o2    = y.get();

        bool ret            = cast_bool(o1) && cast_bool(o2);        
        return OBool(ret);
    }
    static dynamic::Type eval_return(int, const dynamic::Type[], int, const dynamic::Type args[])
    {
        if (args[0] == Any::get_static_type() || args[1] == Any::get_static_type())
            return Type();
        else
            return OBool::get_static_type();
    };
};

struct f_elem_or : dynamic::register_function_template_return<f_elem_or, functions::elem_or>
{
    static dynamic::OBool eval(const dynamic::Template& x, const dynamic::Template& y) 
    {
        const object& o1    = x.get();
        const object& o2    = y.get();

        bool ret            = cast_bool(o1) || cast_bool(o2);        
        return OBool(ret);
    }
    static dynamic::Type eval_return(int, const dynamic::Type[], int, const dynamic::Type args[])
    {
        if (args[0] == Any::get_static_type() || args[1] == Any::get_static_type())
            return Type();
        else
            return OBool::get_static_type();
    };
};

struct f_elem_xor : dynamic::register_function_template_return<f_elem_xor, functions::elem_xor>
{
    static dynamic::OBool eval(const dynamic::Template& x, const dynamic::Template& y) 
    {
        const object& o1    = x.get();
        const object& o2    = y.get();

        bool ret            = cast_bool(o1) ^ cast_bool(o2);        
        return OBool(ret);
    }
    static dynamic::Type eval_return(int, const dynamic::Type[], int, const dynamic::Type args[])
    {
        if (args[0] == Any::get_static_type() || args[1] == Any::get_static_type())
            return Type();
        else
            return OBool::get_static_type();
    };
};

struct f_max : dynamic::register_function_template_return<f_max, functions::max>
{
    static dynamic::Template eval(const dynamic::OType& t_ret, const dynamic::Template& x, 
                             const dynamic::Template& y) 
    {
        const object& o1    = x.get();
        const object& o2    = y.get();

        return (bool)(o1 > o2) ? object(t_ret.get(), o1) : object(t_ret.get(), o2);
    }
    static dynamic::Type eval_return(int, const dynamic::Type[], int, const dynamic::Type args[])
    {
        Type t1 = args[0];
        Type t2 = args[1];

        if (t1 == Any::get_static_type() || t2 == Any::get_static_type())
            return Type();

        //we require gt function to exist and return bool
        Type gt     = operations::return_type(functions::op_gt::eval(), 2, args);

        if (gt != OBool::get_static_type())
            return Type();

        Type ret    = operations::unify_types(t1, t2);
        return ret;
    };
};

struct f_min : dynamic::register_function_template_return<f_min, functions::min>
{
    static dynamic::Template eval(const dynamic::OType& t_ret, const dynamic::Template& x, 
                             const dynamic::Template& y) 
    {
        const object& o1    = x.get();
        const object& o2    = y.get();

        return (bool)(o1 < o2) ? object(t_ret.get(), o1) : object(t_ret.get(), o2);
    }

    static dynamic::Type eval_return(int, const dynamic::Type[], int, const dynamic::Type args[])
    {
        Type t1 = args[0];
        Type t2 = args[1];

        if (t1 == Any::get_static_type() || t2 == Any::get_static_type())
            return Type();

        //we require gt function to exist and return bool
        Type gt     = operations::return_type(functions::op_gt::eval(), 2, args);

        if (gt != OBool::get_static_type())
            return Type();

        Type ret    = operations::unify_types(t1, t2);
        return ret;
    };
};

struct f_eeq_nan : dynamic::register_function_template_return<f_eeq_nan, functions::eeq_nan>
{
    static dynamic::Template eval(const dynamic::OType& t_ret, const dynamic::Template& x, 
                             const dynamic::Template& y) 
    {
        const object& o1    = x.get();
        const object& o2    = y.get();

        return object(t_ret.get(), o1 == o2);
    }
    static dynamic::Type eval_return(int, const dynamic::Type[], int n_arg, const dynamic::Type args[])
    {
        Type t1 = args[0];
        Type t2 = args[1];

        if (t1 == Any::get_static_type() || t2 == Any::get_static_type())
            return Type();

        Type t      = operations::return_type(functions::op_eeq::eval(), n_arg, args);
        return t;
    };
};

struct f_neq_nan : dynamic::register_function_template_return<f_neq_nan, functions::neq_nan>
{
    static dynamic::Template eval(const dynamic::OType& t_ret, const dynamic::Template& x, 
                             const dynamic::Template& y) 
    {
        const object& o1    = x.get();
        const object& o2    = y.get();

        return object(t_ret.get(), o1 != o2);
    }
    static dynamic::Type eval_return(int, const dynamic::Type[], int n_arg, const dynamic::Type args[])
    {
        Type t1 = args[0];
        Type t2 = args[1];

        if (t1 == Any::get_static_type() || t2 == Any::get_static_type())
            return Type();

        Type t      = operations::return_type(functions::op_neq::eval(), n_arg, args);
        return t;
    };
};

};};
