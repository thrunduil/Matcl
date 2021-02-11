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

#pragma warning(push)
#pragma warning(disable: 4244)
#pragma warning(disable: 4800) //forcing value to bool 'true' or 'false'

#include "test_dynamic.h"

#include "eval_cons.h"

#include "matcl-dynamic/matcl_dynamic.h"
#include "matcl-scalar/objects/typed_object_functions.h"
#include "matcl-scalar/lib_functions/manip.h"
#include "matcl-scalar/IO/scalar_io.h"

#include "rand_scalars.h"

#pragma warning(pop)

#include "matcl-core/IO/archive.h"
#include <iostream>

namespace matcl { namespace test
{

void test_dynamic::make()
{
    test_basic();
    test_object();
    test_one();
    test_object_type();
    test_special_types();
    test_cons();
    test_cons_val();
    test_object_func();
    test_object_func_val();
    test_func_compare();
    test_func_uminus();
    test_func_reim();
    test_predefined_func(dynamic::functions::op_plus::eval(), false);
    test_predefined_func(dynamic::functions::op_minus::eval(), false);
    test_predefined_func(dynamic::functions::op_mul::eval(), false);
    test_predefined_func(dynamic::functions::op_div::eval(), true);
    test_predefined_func(dynamic::functions::idiv::eval(), false);
    test_predefined_unify();    

    test_compite_func();
    //value tests for compare, uminus, reim are performed in test_gmp_object
};

void test_dynamic::test_basic()
{
    mdy::object a = mdy::object(mdy::OComplex(1.0));
    mdy::object b = mdy::object(mdy::OFloat_complex(2.0f));

    b = a;
    a = b;
};

bool test_dynamic::test_predefined_func_impl(const dynamic::function_name& func, bool promote_int)
{
    mdy::Type t_int      = (promote_int == false) ? mdy::predefined::type_int() : mdy::predefined::type_real();
    mdy::Type t_float    = mdy::predefined::type_float();
    mdy::Type t_real     = mdy::predefined::type_real();
    mdy::Type t_compl    = mdy::predefined::type_complex();
    mdy::Type t_fcompl   = mdy::predefined::type_float_complex();

    mdy::object obj_int      = mdy::object(mdy::OInteger(1));
    mdy::object obj_float    = mdy::object(t_float);
    mdy::object obj_real     = mdy::object(t_real);
    mdy::object obj_compl    = mdy::object(t_compl);
    mdy::object obj_fcompl   = mdy::object(t_fcompl);

    mdy::object o_ii;
    mdy::object o_if;
    mdy::object o_ir;
    mdy::object o_ic;
    mdy::object o_ifc;

    mdy::eval_function::eval(func, o_ii, obj_int, obj_int);
    mdy::eval_function::eval(func, o_if, obj_int, obj_float);
    mdy::eval_function::eval(func, o_ir, obj_int, obj_real);
    mdy::eval_function::eval(func, o_ic, obj_int, obj_compl);
    mdy::eval_function::eval(func, o_ifc, obj_int, obj_fcompl);

    mdy::Type t_ii       = mdy::operations::unify_types(t_int,t_int);
    mdy::Type t_if       = mdy::operations::unify_types(t_int,t_float);
    mdy::Type t_ir       = mdy::operations::unify_types(t_int,t_real);
    mdy::Type t_ic       = mdy::operations::unify_types(t_int,t_compl);
    mdy::Type t_ifc      = mdy::operations::unify_types(t_int,t_fcompl);

    if (t_ii != o_ii.get_type() || t_if != o_if.get_type() || t_ir != o_ir.get_type() 
        || t_ic != o_ic.get_type() ||t_ifc!= o_ifc.get_type())
    {
        return false;
    };

    mdy::object o_fi;
    mdy::object o_ff;
    mdy::object o_fr;
    mdy::object o_fc;
    mdy::object o_ffc;

    mdy::eval_function::eval(func, o_fi, obj_float, obj_int);
    mdy::eval_function::eval(func, o_ff, obj_float, obj_float);
    mdy::eval_function::eval(func, o_fr, obj_float, obj_real);
    mdy::eval_function::eval(func, o_fc, obj_float, obj_compl);
    mdy::eval_function::eval(func, o_ffc,obj_float, obj_fcompl);

    mdy::Type t_fi       = mdy::operations::unify_types(t_float,t_int);
    mdy::Type t_ff       = mdy::operations::unify_types(t_float,t_float);
    mdy::Type t_fr       = mdy::operations::unify_types(t_float,t_real);
    mdy::Type t_fc       = mdy::operations::unify_types(t_float,t_compl);
    mdy::Type t_ffc      = mdy::operations::unify_types(t_float,t_fcompl);

    if (t_fi != o_fi.get_type() || t_ff != o_ff.get_type() || t_fr != o_fr.get_type() 
        || t_fc != o_fc.get_type() || t_ffc!= o_ffc.get_type())
    {
        return false;
    };

    mdy::object o_ri;
    mdy::object o_rf;
    mdy::object o_rr;
    mdy::object o_rc;
    mdy::object o_rfc;

    mdy::eval_function::eval(func, o_ri,  obj_real, obj_int);
    mdy::eval_function::eval(func, o_rf,  obj_real, obj_float);
    mdy::eval_function::eval(func, o_rr,  obj_real, obj_real);
    mdy::eval_function::eval(func, o_rc,  obj_real, obj_compl);
    mdy::eval_function::eval(func, o_rfc, obj_real, obj_fcompl);

    mdy::Type t_ri       = mdy::operations::unify_types(t_real,t_int);
    mdy::Type t_rf       = mdy::operations::unify_types(t_real,t_float);
    mdy::Type t_rr       = mdy::operations::unify_types(t_real,t_real);
    mdy::Type t_rc       = mdy::operations::unify_types(t_real,t_compl);
    mdy::Type t_rfc      = mdy::operations::unify_types(t_real,t_fcompl);

    if (t_ri != o_ri.get_type() || t_rf != o_rf.get_type() || t_rr != o_rr.get_type() 
        || t_rc != o_rc.get_type() ||t_rfc!= o_rfc.get_type())
    {
        return false;
    };

    mdy::object o_ci;
    mdy::object o_cf;
    mdy::object o_cr;
    mdy::object o_cc;
    mdy::object o_cfc;

    mdy::eval_function::eval(func, o_ci,  obj_compl, obj_int);
    mdy::eval_function::eval(func, o_cf,  obj_compl, obj_float);
    mdy::eval_function::eval(func, o_cr,  obj_compl, obj_real);
    mdy::eval_function::eval(func, o_cc,  obj_compl, obj_compl);
    mdy::eval_function::eval(func, o_cfc, obj_compl, obj_fcompl);

    mdy::Type t_ci       = mdy::operations::unify_types(t_compl,t_int);
    mdy::Type t_cf       = mdy::operations::unify_types(t_compl,t_float);
    mdy::Type t_cr       = mdy::operations::unify_types(t_compl,t_real);
    mdy::Type t_cc       = mdy::operations::unify_types(t_compl,t_compl);
    mdy::Type t_cfc      = mdy::operations::unify_types(t_compl,t_fcompl);

    if (t_ci != o_ci.get_type() || t_cf != o_cf.get_type() || t_cr != o_cr.get_type() 
        || t_cc != o_cc.get_type() ||t_cfc!= o_cfc.get_type())
    {
        return false;
    };

    mdy::object o_fci;
    mdy::object o_fcf;
    mdy::object o_fcr;
    mdy::object o_fcc;
    mdy::object o_fcfc;

    mdy::eval_function::eval(func, o_fci,  obj_fcompl, obj_int);
    mdy::eval_function::eval(func, o_fcf,  obj_fcompl, obj_float);
    mdy::eval_function::eval(func, o_fcr,  obj_fcompl, obj_real);
    mdy::eval_function::eval(func, o_fcc,  obj_fcompl, obj_compl);
    mdy::eval_function::eval(func, o_fcfc, obj_fcompl, obj_fcompl);

    mdy::Type t_fci       = mdy::operations::unify_types(t_fcompl,t_int);
    mdy::Type t_fcf       = mdy::operations::unify_types(t_fcompl,t_float);
    mdy::Type t_fcr       = mdy::operations::unify_types(t_fcompl,t_real);
    mdy::Type t_fcc       = mdy::operations::unify_types(t_fcompl,t_compl);
    mdy::Type t_fcfc      = mdy::operations::unify_types(t_fcompl,t_fcompl);

    if (t_fci != o_fci.get_type() || t_fcf != o_fcf.get_type() || t_fcr != o_fcr.get_type() 
        || t_fcc != o_fcc.get_type() || t_fcfc!= o_fcfc.get_type())
    {
        return false;
    };

    return true;
};

bool test_dynamic::test_predefined_unify_impl()
{
    mdy::Type t_int      = mdy::predefined::type_int();
    mdy::Type t_float    = mdy::predefined::type_float();
    mdy::Type t_real     = mdy::predefined::type_real();
    mdy::Type t_compl    = mdy::predefined::type_complex();
    mdy::Type t_fcompl   = mdy::predefined::type_float_complex();
    mdy::Type t_unit     = mdy::predefined::type_unit();
    mdy::Type t_any      = mdy::predefined::type_any();
    mdy::Type t_string   = mdy::predefined::type_string();

    mdy::Type t_ii       = mdy::operations::unify_types(t_int,t_int);
    mdy::Type t_if       = mdy::operations::unify_types(t_int,t_float);
    mdy::Type t_ir       = mdy::operations::unify_types(t_int,t_real);
    mdy::Type t_ic       = mdy::operations::unify_types(t_int,t_compl);
    mdy::Type t_ifc      = mdy::operations::unify_types(t_int,t_fcompl);
    mdy::Type t_iu       = mdy::operations::unify_types(t_int,t_unit);
    mdy::Type t_ia       = mdy::operations::unify_types(t_int,t_any);
    mdy::Type t_is       = mdy::operations::unify_types(t_int,t_string);

    if (t_ii != t_int || t_if != t_real || t_ir != t_real || t_ic != t_compl)
    {
        return false;
    };
    if (t_ifc!= t_compl || t_iu != t_any || t_ia != t_any || t_is != t_any)
    {
        return false;
    };

    mdy::Type t_fi       = mdy::operations::unify_types(t_float,t_int);
    mdy::Type t_ff       = mdy::operations::unify_types(t_float,t_float);
    mdy::Type t_fr       = mdy::operations::unify_types(t_float,t_real);
    mdy::Type t_fc       = mdy::operations::unify_types(t_float,t_compl);
    mdy::Type t_ffc      = mdy::operations::unify_types(t_float,t_fcompl);
    mdy::Type t_fu       = mdy::operations::unify_types(t_float,t_unit);
    mdy::Type t_fa       = mdy::operations::unify_types(t_float,t_any);
    mdy::Type t_fs       = mdy::operations::unify_types(t_float,t_string);

    if ( t_fi != t_real || t_ff != t_float || t_fr != t_real || t_fc != t_compl)
        return false;
    if ( t_ffc!= t_fcompl || t_fu != t_any || t_fa != t_any || t_fs != t_any)
        return false;

    mdy::Type t_ri       = mdy::operations::unify_types(t_real,t_int);
    mdy::Type t_rf       = mdy::operations::unify_types(t_real,t_float);
    mdy::Type t_rr       = mdy::operations::unify_types(t_real,t_real);
    mdy::Type t_rc       = mdy::operations::unify_types(t_real,t_compl);
    mdy::Type t_rfc      = mdy::operations::unify_types(t_real,t_fcompl);
    mdy::Type t_ru       = mdy::operations::unify_types(t_real,t_unit);
    mdy::Type t_ra       = mdy::operations::unify_types(t_real,t_any);
    mdy::Type t_rs       = mdy::operations::unify_types(t_real,t_string);

    if ( t_ri != t_real || t_rf != t_real || t_rr != t_real || t_rc != t_compl)
        return false;
    if ( t_rfc!= t_compl || t_ru != t_any || t_ra != t_any || t_rs != t_any)
        return false;

    mdy::Type t_ci       = mdy::operations::unify_types(t_compl,t_int);
    mdy::Type t_cf       = mdy::operations::unify_types(t_compl,t_float);
    mdy::Type t_cr       = mdy::operations::unify_types(t_compl,t_real);
    mdy::Type t_cc       = mdy::operations::unify_types(t_compl,t_compl);
    mdy::Type t_cfc      = mdy::operations::unify_types(t_compl,t_fcompl);
    mdy::Type t_cu       = mdy::operations::unify_types(t_compl,t_unit);
    mdy::Type t_ca       = mdy::operations::unify_types(t_compl,t_any);
    mdy::Type t_cs       = mdy::operations::unify_types(t_compl,t_string);

    if ( t_ci != t_compl || t_cf != t_compl || t_cr != t_compl || t_cc != t_compl)
        return false;
    if ( t_cfc!= t_compl || t_cu != t_any || t_ca != t_any || t_cs != t_any)
        return false;

    mdy::Type t_fci      = mdy::operations::unify_types(t_fcompl,t_int);
    mdy::Type t_fcf      = mdy::operations::unify_types(t_fcompl,t_float);
    mdy::Type t_fcr      = mdy::operations::unify_types(t_fcompl,t_real);
    mdy::Type t_fcc      = mdy::operations::unify_types(t_fcompl,t_compl);
    mdy::Type t_fcfc     = mdy::operations::unify_types(t_fcompl,t_fcompl);
    mdy::Type t_fcu      = mdy::operations::unify_types(t_fcompl,t_unit);
    mdy::Type t_fca      = mdy::operations::unify_types(t_fcompl,t_any);
    mdy::Type t_fcs      = mdy::operations::unify_types(t_fcompl,t_string);

    if ( t_fci != t_compl || t_fcf != t_fcompl || t_fcr != t_compl || t_fcc != t_compl)
        return false;
    if ( t_fcfc!= t_fcompl || t_fcu != t_any || t_fca != t_any || t_fcs != t_any)
        return false;

    mdy::Type t_ui       = mdy::operations::unify_types(t_unit,t_int);
    mdy::Type t_uf       = mdy::operations::unify_types(t_unit,t_float);
    mdy::Type t_ur       = mdy::operations::unify_types(t_unit,t_real);
    mdy::Type t_uc       = mdy::operations::unify_types(t_unit,t_compl);
    mdy::Type t_ufc      = mdy::operations::unify_types(t_unit,t_fcompl);
    mdy::Type t_uu       = mdy::operations::unify_types(t_unit,t_unit);
    mdy::Type t_ua       = mdy::operations::unify_types(t_unit,t_any);
    mdy::Type t_us       = mdy::operations::unify_types(t_unit,t_string);

    if ( t_ui != t_any || t_uf != t_any || t_ur != t_any || t_uc != t_any)
        return false;
    if ( t_ufc!= t_any || t_uu != t_unit || t_ua != t_any || t_us != t_any)
        return false;

    mdy::Type t_ai       = mdy::operations::unify_types(t_any,t_int);
    mdy::Type t_af       = mdy::operations::unify_types(t_any,t_float);
    mdy::Type t_ar       = mdy::operations::unify_types(t_any,t_real);
    mdy::Type t_ac       = mdy::operations::unify_types(t_any,t_compl);
    mdy::Type t_afc      = mdy::operations::unify_types(t_any,t_fcompl);
    mdy::Type t_au       = mdy::operations::unify_types(t_any,t_unit);
    mdy::Type t_aa       = mdy::operations::unify_types(t_any,t_any);
    mdy::Type t_as       = mdy::operations::unify_types(t_any,t_string);

    if ( t_ai != t_any || t_af != t_any || t_ar != t_any || t_ac != t_any)
        return false;
    if ( t_afc!= t_any || t_au != t_any || t_aa != t_any || t_as != t_any)
        return false;

    mdy::Type t_si       = mdy::operations::unify_types(t_string,t_int);
    mdy::Type t_sf       = mdy::operations::unify_types(t_string,t_float);
    mdy::Type t_sr       = mdy::operations::unify_types(t_string,t_real);
    mdy::Type t_sc       = mdy::operations::unify_types(t_string,t_compl);
    mdy::Type t_sfc      = mdy::operations::unify_types(t_string,t_fcompl);
    mdy::Type t_su       = mdy::operations::unify_types(t_string,t_unit);
    mdy::Type t_sa       = mdy::operations::unify_types(t_string,t_any);
    mdy::Type t_ss       = mdy::operations::unify_types(t_string,t_string);

    if ( t_si != t_any || t_sf != t_any || t_sr != t_any || t_sc != t_any)
        return false;
    if ( t_sfc!= t_any || t_su != t_any || t_sa != t_any || t_ss != t_string)
        return false;

    return true;
};

void test_dynamic::test_predefined_unify()
{
    bool res = false;

    try
    {
        res = test_predefined_unify_impl();
    }
    catch(std::exception& ex)
    {
        out_stream << "predefined unify EXCEPTION: " << ex.what() << "\n";
    };

    if (res == true)
        out_stream << "predefined unify: " << "ok" << "\n";
    else
        out_stream << "predefined unify: " << "FAILED" << "\n";
};
void test_dynamic::test_object()
{
    bool res = false;

    try
    {
        res = test_object_impl();
    }
    catch(std::exception& ex)
    {
        out_stream << "test_object EXCEPTION: " << ex.what() << "\n";
    };

    if (res == true)
        out_stream << "object: " << "ok" << "\n";
    else
        out_stream << "object: " << "FAILED" << "\n";
};

void test_dynamic::test_one()
{
    bool res = false;

    try
    {
        res = test_one_impl();
    }
    catch(std::exception& ex)
    {
        out_stream << "test_one EXCEPTION: " << ex.what() << "\n";
    };

    if (res == true)
        out_stream << "one: " << "ok" << "\n";
    else
        out_stream << "one: " << "FAILED" << "\n";
};

void test_dynamic::test_object_type()
{
    bool res = false;

    try
    {
        res = test_object_type_impl();
    }
    catch(std::exception& ex)
    {
        out_stream << "test_object_type EXCEPTION: " << ex.what() << "\n";
    };

    if (res == true)
        out_stream << "typed object: " << "ok" << "\n";
    else
        out_stream << "typed object: " << "FAILED" << "\n";
};
void test_dynamic::test_special_types()
{
    bool res = false;

    try
    {
        res = test_special_types_impl();
    }
    catch(std::exception& ex)
    {
        out_stream << "test_special_types EXCEPTION: " << ex.what() << "\n";
    };

    if (res == true)
        out_stream << "special types: " << "ok" << "\n";
    else
        out_stream << "special types: " << "FAILED" << "\n";
};

void test_dynamic::test_cons()
{
    bool res = false;

    try
    {
        res = test_cons_impl();
    }
    catch(std::exception& ex)
    {
        out_stream << "test_cons EXCEPTION: " << ex.what() << "\n";
    };

    if (res == true)
        out_stream << "constructors: " << "ok" << "\n";
    else
        out_stream << "constructors: " << "FAILED" << "\n";
};

void test_dynamic::test_cons_val()
{
    bool res = false;

    try
    {
        res = test_cons_val_impl();
    }
    catch(std::exception& ex)
    {
        out_stream << "test_cons_val EXCEPTION: " << ex.what() << "\n";
    };

    if (res == true)
        out_stream << "constructors val: " << "ok" << "\n";
    else
        out_stream << "constructors val: " << "FAILED" << "\n";
};

void test_dynamic::test_object_func()
{
    bool res = false;

    try
    {
        res = test_object_func_impl();
    }
    catch(std::exception& ex)
    {
        out_stream << "test_object_func EXCEPTION: " << ex.what() << "\n";
    };

    if (res == true)
        out_stream << "object_func: " << "ok" << "\n";
    else
        out_stream << "object_func: " << "FAILED" << "\n";
};

void test_dynamic::test_object_func_val()
{
    bool res = false;

    try
    {
        res = test_object_func_val_impl();
    }
    catch(std::exception& ex)
    {
        out_stream << "test_object_func_val EXCEPTION: " << ex.what() << "\n";
    };

    if (res == true)
        out_stream << "object_func val: " << "ok" << "\n";
    else
        out_stream << "object_func val: " << "FAILED" << "\n";
};

void test_dynamic::test_func_uminus()
{
    bool res = false;

    try
    {
        res = test_func_uminus_impl();
    }
    catch(std::exception& ex)
    {
        out_stream << "test_func_uminus EXCEPTION: " << ex.what() << "\n";
    };

    if (res == true)
        out_stream << "func_uminus: " << "ok" << "\n";
    else
        out_stream << "func_uminus: " << "FAILED" << "\n";
};

void test_dynamic::test_func_reim()
{
    bool res = false;

    try
    {
        res = test_func_reim_impl();
    }
    catch(std::exception& ex)
    {
        out_stream << "test_func_reim EXCEPTION: " << ex.what() << "\n";
    };

    if (res == true)
        out_stream << "reim: " << "ok" << "\n";
    else
        out_stream << "reim: " << "FAILED" << "\n";
};

void test_dynamic::test_func_compare()
{
    bool res = false;

    try
    {
        res = test_func_compare_impl();
    }
    catch(std::exception& ex)
    {
        out_stream << "test_func_compare EXCEPTION: " << ex.what() << "\n";
    };

    if (res == true)
        out_stream << "compare_func: " << "ok" << "\n";
    else
        out_stream << "compare_func: " << "FAILED" << "\n";
};

void test_dynamic::test_predefined_func(const dynamic::function_name& func, bool promote_int)
{
    bool res = false;

    try
    {
        res = test_predefined_func_impl(func, promote_int);
    }
    catch(std::exception& ex)
    {
        out_stream << "predefined function "
                << func.to_string() << " EXCEPTION: " << ex.what() << "\n";
    };

    if (res == true)
        out_stream << "predefined function: " 
                  << func.to_string() << " ok" << "\n";
    else
        out_stream << "predefined function: " 
                  << func.to_string() << " FAILED" << "\n";
};

bool test_dynamic::test_object_impl()
{
    double res  = 0;

    // constructors
    {
        mdy::object a;
        if (a.is_null() == false)
            res += 1;
    }
    {
        mdy::Type t;

        mdy::object a0;
        mdy::object a(t,a0);

        if (a.is_null() == false)
            res += 1;
    };

    try
    {
        mdy::Type t = mdy::predefined::type_float();

        mdy::object a0;
        mdy::object a(t,a0);
        res += 1;
    }
    catch(...){};

    try
    {
        mdy::Type t;

        mdy::object a0(1.1);

        mdy::object a(t,a0);
        res += 1;
    }
    catch(...){};

    {
        mdy::Type t = mdy::predefined::type_float();

        mdy::object a0(1.1);

        mdy::object a(t,a0);
        if (a.get_type() != t)
            res += 1;
    }
    try
    {
        mdy::Type t;

        mdy::object a0(1.1);

        mdy::object a(t,std::move(a0));
        res += 1;
    }
    catch(...){};

    {
        mdy::Type t = mdy::predefined::type_float();

        mdy::object a0(1.1);

        mdy::object a(t,std::move(a0));
        if (a.get_type() != t)
            res += 1;
    }
    {
        mdy::Type t = mdy::predefined::type_float();

        mdy::object a0(1.1f);

        mdy::object a(t,a0);
        if (a.get_type() != t)
            res += 1;
    }
    {
        mdy::Type t = mdy::predefined::type_float();

        mdy::object a0(1.1f);

        mdy::object a(t,std::move(a0));
        if (a.get_type() != t)
            res += 1;
    }
    {
        mdy::Type t;

        mdy::object a(t);

        if (a.get_type() != t)
            res += 1;
    }
    {
        mdy::Type t = mdy::predefined::type_float();

        mdy::object a(t);

        if (a.get_type() != t)
            res += 1;
    }
    {
        mdy::object a0(1.1);
        mdy::object a(a0);

        if (a.is_unique() == true)
            res += 1;
    }
    {
        mdy::object a0(1.1);
        mdy::object a(std::move(a0));

        if (a.is_unique() == false)
            res += 1;
    }

    {
        mdy::object a(true);

        if (a.get_type() != mdy::predefined::type_bool())
            res += 1;
    }
    {
        mdy::object a(1);

        if (a.get_type() != mdy::predefined::type_int())
            res += 1;
    }
    {
        mdy::object a(1.1f);

        if (a.get_type() != mdy::predefined::type_float())
            res += 1;
    }
    {
        mdy::object a(1.1);

        if (a.get_type() != mdy::predefined::type_real())
            res += 1;
    }
    {
        mdy::object a(Float_complex(1.1f, 1.f));

        if (a.get_type() != mdy::predefined::type_float_complex())
            res += 1;
    }
    {
        mdy::object a(Complex(1.1f, 1.f));

        if (a.get_type() != mdy::predefined::type_complex())
            res += 1;
    }
    {
        std::string s("abc");
        mdy::object a(s);

        if (a.get_type() != mdy::predefined::type_string())
            res += 1;
    }
    {
        std::string s("abc");
        mdy::object a(std::move(s));

        if (a.get_type() != mdy::predefined::type_string())
            res += 1;
    }

    //assignment
    {
        mdy::object a;
        a = a;

        if (a.get_type() != mdy::Type())
            res += 1;
    }
    {
        mdy::object a(Complex(1.1));
        a = a;

        if (a.get_type() != mdy::predefined::type_complex())
            res += 1;
    }
    try
    {
        mdy::object a(Complex(1.1));
        mdy::object b;
        a = b;
        res += 1;
    }
    catch(...){};

    try
    {
        mdy::object a;
        mdy::object b(1.1f);
        a = b;
        res += 1;
    }
    catch(...){};

    {
        mdy::object a;
        mdy::object b;
        a = b;

        if (a.get_type() != mdy::Type())
            res += 1;
    }
    {
        mdy::object a(Complex(1.1));
        mdy::object b(1.1f);
        a = b;

        if (a.get_type() != mdy::predefined::type_complex())
            res += 1;
    }
    {
        mdy::object a(Complex(1.1));
        mdy::object b(1.1f);
        a = std::move(b);

        if (a.get_type() != mdy::predefined::type_complex())
            res += 1;
    }
    {
        mdy::object a(Complex(1.1));
        mdy::object b(Complex(1.1,2));
        a = b;

        if (a.get_type() != mdy::predefined::type_complex())
            res += 1;
    }
    {
        mdy::object a(Complex(1.1));
        mdy::object b(Complex(1.1,2));
        a = std::move(b);

        if (a.get_type() != mdy::predefined::type_complex())
            res += 1;
    }
    try
    {
        mdy::object a;
        mdy::object b(Complex(1.1,2));
        a = b;
        res += 1;
    }
    catch(...){};

    try
    {
        mdy::object a(Complex(1.1,2));
        mdy::object b("abc");
        a = b;
        res += 1;
    }
    catch(...){};

    {
        mdy::object a;
        a.reset(a);

        if (a.get_type() != mdy::Type())
            res += 1;
    }
    {
        mdy::object a(Complex(1.1));
        a.reset(a);

        if (a.get_type() != mdy::predefined::type_complex())
            res += 1;
    }
    {
        mdy::object a(Complex(1.1));
        mdy::object b(1.1f);
        a.reset(b);

        if (a.get_type() != mdy::predefined::type_float())
            res += 1;
    }
    {
        mdy::object a(Complex(1.1));
        mdy::object b(1.1f);
        a.reset(std::move(b));

        if (a.get_type() != mdy::predefined::type_float())
            res += 1;
    }
    {
        mdy::object a(Complex(1.1));
        mdy::object b(Complex(1.1,2));
        a.reset(b);

        if (a.get_type() != mdy::predefined::type_complex())
            res += 1;
    }
    {
        mdy::object a(Complex(1.1));
        mdy::object b(Complex(1.1,2));
        a.reset(std::move(b));

        if (a.get_type() != mdy::predefined::type_complex())
            res += 1;
    }
    {
        mdy::object a;
        mdy::object b(Complex(1.1,2));
        a.reset(b);
        if (a.get_type() != mdy::predefined::type_complex())
            res += 1;
    }
    {
        mdy::object a;
        mdy::object b(Complex(1.1,2));
        a.reset(std::move(b));
        if (a.get_type() != mdy::predefined::type_complex())
            res += 1;
    }
    {
        mdy::object a(Complex(1.1,2));
        mdy::object b;
        a.reset(b);
        if (a.get_type() != mdy::Type())
            res += 1;
    }
    {
        mdy::object a(Complex(1.1,2));
        mdy::object b;
        a.reset(std::move(b));
        if (a.get_type() != mdy::Type())
            res += 1;
    }

    {
        mdy::object a(Complex(1.1,2));
        mdy::object b("abc");

        a.reset(b);

        if (a.get_type() != mdy::predefined::type_string())
            res += 1;
    }

    // functions
    {
        mdy::object a;
        mdy::object b(a);
        mdy::object c = b.clone();

        if (c.is_unique() == false)
            res += 1;
    }
    {
        mdy::object a(1);
        mdy::object b(a);
        mdy::object c = b.clone();

        if (c.is_unique() == false)
            res += 1;
    }
    {
        mdy::object a;
        mdy::object b(a);
        b.make_unique();

        if (b.is_unique() == false)
            res += 1;
    }
    {
        mdy::object a(1);
        mdy::object b(a);
        b.make_unique();

        if (b.is_unique() == false)
            res += 1;
    }
    {
        mdy::object a;
        mdy::object b(1.1f);
        swap(a,b);

        if (a.get_type() != mdy::predefined::type_float())
            res += 1;
    }
    {
        mdy::object a(1);
        mdy::object b(1.1f);
        swap(a,b);

        if (a.get_type() != mdy::predefined::type_float())
            res += 1;
        if (b.get_type() != mdy::predefined::type_int())
            res += 1;
    }
    {
        mdy::object a;

        if (a.get_type() != mdy::Type())
            res += 1;
        if (a.get_static_type() != mdy::predefined::type_any())
            res += 1;
    }
    {
        mdy::object a(1);

        if (a.get_type() != mdy::predefined::type_int())
            res += 1;
        if (a.get_static_type() != mdy::predefined::type_any())
            res += 1;
    }
    {
        mdy::object a(1);
        mdy::object b;

        if (a.is_null() == true)
            res += 1;
        if (b.is_null() == false)
            res += 1;
    }

    {
        mdy::object a;

        if (a.is_zero() == false)
            res += 1;
    }
    {
        mdy::object a(1);
        mdy::object b(0);

        if (a.is_zero() == true)
            res += 1;
        if (b.is_zero() == false)
            res += 1;
    }
    {
        mdy::object a;

        if (a.is_one() == true)
            res += 1;
    }
    {
        mdy::object a(1.0);
        mdy::object b(0.0);

        if (a.is_one() == false)
            res += 1;
        if (b.is_one() == true)
            res += 1;
    }
    {
        mdy::object a;

        if ((bool)a == true)
            res += 1;
    }
    {
        mdy::object a(1.0);
        mdy::object b(0.0);

        if ((bool)a == false)
            res += 1;
        if ((bool)b == true)
            res += 1;
    }
    try
    {
        mdy::object a;
        mdy::Type t = mdy::predefined::type_int();
        mdy::object b = convert(t,a);
        res += 1;
    }
    catch(...){};

    try
    {
        mdy::object a;
        mdy::Type t = mdy::predefined::type_int();
        mdy::object b = cast(t,a);
        res += 1;
    }
    catch(...){};

    try
    {
        mdy::object a(Complex(1));
        mdy::Type t = mdy::predefined::type_int();
        mdy::object b = convert(t,a);
        res += 1;
    }
    catch(...){};

    {
        mdy::object a(1.0);
        mdy::Type t = mdy::predefined::type_float_complex();
        mdy::object b = convert(t,a);

        if (b.get_type() != t)
            res += 1;
    }
    {
        mdy::object a(1.0);
        mdy::Type t = mdy::predefined::type_int();
        mdy::object b = cast(t,a);

        if (b.get_type() != t)
            res += 1;
    }
    {
        std::ostringstream ss;
        oarchive ia(ss);

        mdy::object a;
        a.serialize(ia.get(), 0);

        std::istringstream ss2(ss.str());
        iarchive ia2(ss2);

        mdy::object b;
        b.serialize(ia2.get(), 0);

        if (a != b)
            res += 1;
    };
    {
        std::ostringstream ss;
        oarchive ia(ss);

        mdy::object a(1);
        a.serialize(ia.get(), 0);

        std::istringstream ss2(ss.str());
        iarchive ia2(ss2);

        mdy::object b;
        b.serialize(ia2.get(), 0);

        if (a != b)
            res += 1;
    };
    {
        std::ostringstream ss;

        mdy::object a;
        save_data(ss, a);

        std::istringstream ss2(ss.str());

        mdy::object b;
        load_data(ss2, b);

        if (a != b)
            res += 1;
    };
    {
        std::ostringstream ss;

        mdy::object a(1);
        save_data(ss, a);

        std::istringstream ss2(ss.str());

        mdy::object b(mdy::predefined::type_int());
        load_data(ss2, b);

        if (a != b)
            res += 1;
    };

    {
        std::ostringstream ss;

        mdy::object a;
        ss << a;

        std::istringstream ss2(ss.str());

        mdy::object b;
        ss2 >> b;

        if (a != b)
            res += 1;
    };
    {
        std::ostringstream ss;

        mdy::object a(1);
        ss << a;

        std::istringstream ss2(ss.str());

        mdy::object b;
        ss2 >> b;

        if (a != b)
            res += 1;
    };

    return (res == 0) ? true : false;
}

bool test_dynamic::test_object_type_impl()
{
    double res = 0;

    {
        OFloat f;
        if (f.is_zero() == false)
            res += 1;
    }

    {
        float val   = 1.1f;
        mdy::object a(val);

        OFloat f(a,dynamic::from_object());
        if (f.get() != val)
            res += 1;
        if (f.get_type() != mdy::predefined::type_float())
            res += 1;
    }
    {
        float val   = 1.1f;
        mdy::object a(val);

        OFloat f(std::move(a),dynamic::from_object());
        if (f.get() != val)
            res += 1;
        if (f.get_type() != mdy::predefined::type_float())
            res += 1;
    }
    {
        double val   = 1.1;
        mdy::object a(val);

        OFloat f(a,dynamic::from_object());

        if (f.get() != (float)val)
            res += 1;
        if (f.get_type() != mdy::predefined::type_float())
            res += 1;
    }
    {
        double val   = 1.1;
        mdy::object a(val);

        OFloat f(std::move(a),dynamic::from_object());
        if (f.get() != (float)val)
            res += 1;
        if (f.get_type() != mdy::predefined::type_float())
            res += 1;
    }
    {
        double val   = 1.1;
        mdy::object a = mdy::object(Complex(val));

        OComplex f(a,dynamic::from_object());
        if (f.get() != (Complex)val)
            res += 1;
        if (f.get_type() != mdy::predefined::type_complex())
            res += 1;
    }
    {
        double val   = 1.1;
        mdy::object a = mdy::object(Complex(val));

        OComplex f(std::move(a),dynamic::from_object());
        if (f.get() != (Complex)val)
            res += 1;
        if (f.get_type() != mdy::predefined::type_complex())
            res += 1;
    }

    {
        OFloat a(1.1f);
        OReal f(a);
        if (f.get() != a.get())
            res += 1;
        if (f.get_type() != mdy::predefined::type_real())
            res += 1;
    }
    {
        OComplex a = OComplex(3,2);
        OFloat_complex f(a);
        if (real(f.get()) != real(a.get()))
            res += 1;
        if (imag(f.get()) != imag(a.get()))
            res += 1;
        if (f.get_type() != mdy::predefined::type_float_complex())
            res += 1;
    }
    {
        Integer a(3);
        OFloat f(a);
        if (f.get() != a)
            res += 1;
        if (f.get_type() != mdy::predefined::type_float())
            res += 1;
    }
    {
        Float a(3.f);
        OReal f(a);
        if (f.get() != a)
            res += 1;
        if (f.get_type() != mdy::predefined::type_real())
            res += 1;
    }
    {
        Complex a = Complex(3,2);
        OFloat_complex f(a);
        if (real(f.get()) != real(a))
            res += 1;
        if (imag(f.get()) != imag(a))
            res += 1;
        if (f.get_type() != mdy::predefined::type_float_complex())
            res += 1;
    }

    {
        Complex a(3.1, -2.1);
        OComplex f(real(a), imag(a));
        if (f.get() != a)
            res += 1;
        if (f.get_type() != mdy::predefined::type_complex())
            res += 1;
    }

    {
        Float v(1.1f);
        OFloat f1(v);
        OFloat f2(f1);

        if (f2.get() != v)
            res += 1;
        if (f2.get_type() != mdy::predefined::type_float())
            res += 1;
    }
    {
        Float v(1.1f);
        OFloat f1(v);
        OFloat f2(std::move(f1));

        if (f2.get() != v)
            res += 1;
        if (f2.get_type() != mdy::predefined::type_float())
            res += 1;
    }
    {
        Real r(1.1);
        Float v(2.1f);

        OFloat f1(v);
        OReal r1(r);
        f1 = r1;

        if (f1.get() != (float)r)
            res += 1;
        if (f1.get_type() != mdy::predefined::type_float())
            res += 1;
    };
    {
        Real r(1.1);
        Float v(2.1f);

        OFloat f1(v);
        OReal r1(r);

        r1 = f1;

        if (r1.get() != v)
            res += 1;
        if (r1.get_type() != mdy::predefined::type_real())
            res += 1;
    };
    {
        Real r(1.1);
        Float v(2.1f);
        OReal r1(r);

        r1 = OFloat(v);

        if (r1.get() != v)
            res += 1;
        if (r1.get_type() != mdy::predefined::type_real())
            res += 1;
    };
    {
        Real r(1.1);
        Float v(2.1f);
        OFloat f1(v);

        f1 = OReal(r);

        if (f1.get() != (float)r)
            res += 1;
        if (f1.get_type() != mdy::predefined::type_float())
            res += 1;
    };

    {
        Float r(1.1f);
        Float v(2.1f);
        OFloat f1(v);
        OFloat f2(r);

        f1 = f2;

        if (f1.get() != r)
            res += 1;
        if (f1.get_type() != mdy::predefined::type_float())
            res += 1;
    };

    {
        Float v(2.1f);
        OFloat f1(v);

        if (f1.get() != v)
            res += 1;
        if (f1.get_type() != mdy::predefined::type_float())
            res += 1;
    };
    {
        Float v1(2.1f);
        Float v2(-2.1f);

        OFloat f1(v1);
        f1.get_unique() = v2;

        if (f1.get() != v2)
            res += 1;
        if (f1.get_type() != mdy::predefined::type_float())
            res += 1;
    };

    {
        Float v(2.1f);
        OFloat f1(v);

        if (f1.get_static_type() != mdy::predefined::type_float())
            res += 1;
        if (f1.get_type() != mdy::predefined::type_float())
            res += 1;
        if (f1.get_type() != mdy::predefined::type_float())
            res += 1;
    };
    {
        Float v(2.1f);
        OFloat f1(v);
        OFloat f2(v);

        swap(f1, f2);

        if (f1.get_static_type() != mdy::predefined::type_float())
            res += 1;
        if (f1.get_type() != mdy::predefined::type_float())
            res += 1;
    };

    {
        Float v1(2.1f);
        Float v2(3.1f);
        OFloat f1(v1);
        OFloat f2(v2);

        f1.reset(f2);

        if (f1.get() != v2)
            res += 1;
        if (f1.get_type() != mdy::predefined::type_float())
            res += 1;
    };
    {
        Float v1(2.1f);
        Float v2(3.1f);
        OFloat f1(v1);
        OFloat f2(v2);

        f1.reset(std::move(f2));

        if (f1.get() != v2)
            res += 1;
        if (f1.get_type() != mdy::predefined::type_float())
            res += 1;
    };
    {
        Float v1(2.1f);
        OFloat f1(v1);

        OFloat f2 = f1.clone();

        if (f2.get() != v1)
            res += 1;
        if (f2.get_type() != mdy::predefined::type_float())
            res += 1;
    };

    {
        Float v1(2.1f);
        OFloat f1(v1);

        if (f1.is_unique() == false)
            res += 1;
    };
    {
        Float v1(2.1f);
        OFloat f1(v1);
        OFloat f2(f1);
        f2.make_unique();

        if (f2.is_unique() == false)
            res += 1;
    };

    {
        Float v1(0.f);
        OFloat f1(v1);

        if (f1.is_null() == true)
            res += 1;
    };
    {
        Float v1(0.f);
        OFloat f1(v1);

        if (f1.is_zero() == false)
            res += 1;
    };
    {
        Float v1(1.f);
        OFloat f1(v1);

        if (f1.is_one() == false)
            res += 1;
    };
    {
        Float v1(0.f);
        OFloat f1(v1);

        if ((bool)f1 == true)
            res += 1;
    };
    {
        Float v1(1.f);
        OFloat f1(v1);
        OReal f2 = f1.convert<Real>();

        if (f2.get() != Real(v1))
            res += 1;
    };

    {
        Complex v1(1.f);
        OComplex f1(v1);
        OFloat_complex f2 = f1.convert<Float_complex>();
    }

    {
        Complex v1(1.f);
        OComplex f1(v1);
        OReal f2 = f1.cast<Real>();

        if (f2.get() != real(v1))
            res += 1;
    };

    {
        Float v1(1.f);
        OFloat f1(v1);
        OReal f2 = f1.cast<Real>();

        if (f2.get() != Real(v1))
            res += 1;
    };
    {
        Float v1(1.f);
        OFloat f1(v1);
        mdy::object f2 = mdy::object(f1);

        if (OFloat(f2,dynamic::from_object()).get() != v1)
            res += 1;
    };
    {
        Float v1(1.f);
        Float v2(2.f);
        OFloat f1(v1);
        OFloat f2(v2);
        swap(f1, f2);

        if (f2.get() != v1)
            res += 1;
        if (f1.get() != v2)
            res += 1;
    };

    {
        Float v1(1.f);
        OFloat f1(v1);

        std::ostringstream ss;

        ss << f1;

        std::istringstream ss2(ss.str());

        OFloat f2;
        ss2 >> f2;

        if (f1.get() != f2.get())
            res += 1;
    };

    return (res == 0.0) ? true : false;
}

bool test_dynamic::test_special_types_impl()
{
    double res = 0;

    //unit type
    {
        mdy::unit_type u;
        mdy::unit_type u2(1);
        mdy::unit_type u3(u);
        u2 = u3;
        u3 = 1;
    }

    {
        mdy::Unit u;
        if (u.is_null() == true)
            res += 1;
        if (u.is_zero() == false)
            res += 1;
        if (u.is_one() == true)
            res += 1;
    }
    {
        mdy::Unit u1;
        mdy::Unit u2 = u1;
        if (u1.get_type() != mdy::predefined::type_unit())
            res += 1;
        if (u1.get_static_type() != mdy::predefined::type_unit())
            res += 1;
    }

    {
        mdy::Unit u1("ans");
        mdy::Unit u2 = u1;
        u2 = OInteger(123);

        mdy::object o;
        mdy::Unit u3(o);

        if (u3.get_type() != mdy::predefined::type_unit())
            res += 1;
    }
    {
        std::ostringstream ss;

        mdy::Unit a;
        ss << a;

        std::istringstream ss2(ss.str());

        mdy::Unit b;
        ss2 >> b;

        if (b.get_type() != mdy::predefined::type_unit())
            res += 1;
    };
    {
        std::ostringstream ss;
        oarchive ia(ss);

        mdy::Unit a;
        a.serialize(ia.get(), 0);

        std::istringstream ss2(ss.str());
        iarchive ia2(ss2);

        mdy::object b;
        b.serialize(ia2.get(), 0);

        if (b.get_type() != mdy::predefined::type_unit())
            res += 1;
    };

    //unit type
    {
        mdy::null_type u;
        mdy::null_type u2(u);
        mdy::null_type u3;
        u2 = u3;
    }
    {
        mdy::ONull u;
        if (u.is_null() == false)
            res += 1;
        if (u.is_zero() == false)
            res += 1;
        if (u.is_one() == true)
            res += 1;
    }
    {
        mdy::ONull u1;
        mdy::ONull u2 = u1;
        if (u1.get_type() != mdy::predefined::type_null())
            res += 1;
        if (u1.get_static_type() != mdy::predefined::type_null())
            res += 1;
    }
    {
        std::ostringstream ss;

        mdy::ONull a;
        ss << a;

        std::istringstream ss2(ss.str());

        mdy::ONull b;
        ss2 >> b;

        if (a != b)
            res += 1;
    };

    {
        std::ostringstream ss;
        oarchive ia(ss);

        mdy::ONull a;
        a.serialize(ia.get(), 0);

        std::istringstream ss2(ss.str());
        iarchive ia2(ss2);

        mdy::object b;
        b.serialize(ia2.get(), 0);

        if (mdy::object(a) != b)
            res += 1;
    };

    //any type
    {
        mdy::any_type a;
        if (a.get_stored().is_null() != true)
            res += 1;
        if ((bool)a == true)
            res += 1;
        if (!a == false)
            res += 1;
    }
    {
        mdy::any_type a = 1;
        if (a.get_stored().is_null() != false)
            res += 1;
        if ((bool)a == false)
            res += 1;
        if (!a == true)
            res += 1;
    }
    {
        mdy::object i(1);
        mdy::any_type a(i);

        if (a.get_stored().get_type() != mdy::predefined::type_int())
            res += 1;
    }
    {
        mdy::object i(1);
        mdy::any_type a(std::move(i));

        if (a.get_stored().get_type() != mdy::predefined::type_int())
            res += 1;
    }
    {
        mdy::object i(1);
        mdy::Any a(i);
        mdy::any_type a2 = mdy::any_type(mdy::object(a));

        if (a2.get_stored().get_type() != mdy::predefined::type_int())
            res += 1;
    }

    {
        mdy::object i(1);
        mdy::Any a(i);
        mdy::any_type a2(a);

        if (a2.get_stored().get_type() != mdy::predefined::type_int())
            res += 1;
    }
    {
        mdy::object i(1);
        mdy::Any a(i);
        mdy::any_type a2(std::move(a));

        if (a2.get_stored().get_type() != mdy::predefined::type_int())
            res += 1;
    }
    {
        mdy::object i(1);
        mdy::Unit a(i);
        mdy::any_type a2(a);

        if (a2.get_stored().get_type() != mdy::predefined::type_unit())
            res += 1;
    }
    {
        mdy::OInteger a(1);
        mdy::any_type a2(a);

        if (a2.get_stored().get_type() != mdy::predefined::type_int())
            res += 1;
    }
    {
        mdy::OInteger a(1);
        mdy::any_type a2(std::move(a));

        if (a2.get_stored().get_type() != mdy::predefined::type_int())
            res += 1;
    }
    {
        Integer i = 1;
        mdy::any_type a2(i);

        if (a2.get_stored().get_type() != mdy::predefined::type_int())
            res += 1;
    }
    {
        std::string s = "abc";
        mdy::any_type a2(s);

        if (a2.get_stored().get_type() != mdy::predefined::type_string())
            res += 1;
    }
    {
        std::string s = "abc";
        mdy::any_type a2(std::move(s));

        if (a2.get_stored().get_type() != mdy::predefined::type_string())
            res += 1;
    }
    {
        Integer i = 1;
        mdy::any_type a2(i);
        mdy::any_type a3(a2);

        if (a3.get_stored().get_type() != mdy::predefined::type_int())
            res += 1;
    }
    {
        Integer i = 1;
        mdy::any_type a2(i);
        mdy::any_type a3(std::move(a2));

        if (a3.get_stored().get_type() != mdy::predefined::type_int())
            res += 1;
    }
    {
        mdy::object i(1);
        mdy::any_type a;
        a = i;

        if (a.get_stored().get_type() != mdy::predefined::type_int())
            res += 1;
    }
    {
        mdy::object i(1);
        mdy::any_type a;
        a = std::move(i);

        if (a.get_stored().get_type() != mdy::predefined::type_int())
            res += 1;
    }
    {
        mdy::object i(1);
        mdy::Any a(i);
        mdy::any_type a2;
        a2 = mdy::any_type(mdy::object(a));

        if (a2.get_stored().get_type() != mdy::predefined::type_int())
            res += 1;
    }
    {
        mdy::object i(1);
        mdy::Any a(i);
        mdy::any_type a2;
        a2 = a;

        if (a2.get_stored().get_type() != mdy::predefined::type_int())
            res += 1;
    }
    {
        mdy::object i(1);
        mdy::Any a(i);
        mdy::any_type a2;
        a2 = std::move(a);

        if (a2.get_stored().get_type() != mdy::predefined::type_int())
            res += 1;
    }
    {
        mdy::object i(1);
        mdy::Unit a(i);
        mdy::any_type a2;
        a2 = a;

        if (a2.get_stored().get_type() != mdy::predefined::type_unit())
            res += 1;
    }
    {
        mdy::OInteger a(1);
        mdy::any_type a2;
        a2 = a;

        if (a2.get_stored().get_type() != mdy::predefined::type_int())
            res += 1;
    }
    {
        mdy::OInteger a(1);
        mdy::any_type a2;
        a2 = std::move(a);

        if (a2.get_stored().get_type() != mdy::predefined::type_int())
            res += 1;
    }
    {
        Integer i = 1;
        mdy::any_type a2;
        a2 = i;

        if (a2.get_stored().get_type() != mdy::predefined::type_int())
            res += 1;
    }
    {
        std::string s = "abc";
        mdy::any_type a2;
        a2 = s;

        if (a2.get_stored().get_type() != mdy::predefined::type_string())
            res += 1;
    }
    {
        std::string s = "abc";
        mdy::any_type a2;
        a2 = std::move(s);

        if (a2.get_stored().get_type() != mdy::predefined::type_string())
            res += 1;
    }
    {
        Integer i = 1;
        mdy::any_type a2(i);
        mdy::any_type a3;
        a3 = a2;

        if (a3.get_stored().get_type() != mdy::predefined::type_int())
            res += 1;
    }
    {
        Integer i = 1;
        mdy::any_type a2(i);
        mdy::any_type a3;
        a3 = std::move(a2);

        if (a3.get_stored().get_type() != mdy::predefined::type_int())
            res += 1;
    }

    {
        Integer i = 1;
        mdy::any_type a2(i);

        if (a2.is_zero() == true)
            res += 1;
    }
    {
        mdy::any_type a2 = mdy::any_type(mdy::object());

        if (a2.is_zero() == false)
            res += 1;
    }
    {
        mdy::any_type a2;

        if (a2.is_zero() == false)
            res += 1;
    }
    {
        mdy::any_type a2;
        mdy::any_type a3 = a2.clone();

        if (a2.is_zero() == false)
            res += 1;
    }

    {
        Integer i = 1;
        mdy::any_type a2(i);
        mdy::any_type a3(1.1);

        if (a2 == a3)
            res += 1;
        if ((a2 != a3) == false)
            res += 1;
    }
    {
        Integer i = 1;
        mdy::any_type a2(i);
        mdy::any_type a3 = mdy::any_type(OInteger(i));

        if ((a2 == a3) == true)
            res += 1;
        if ((a2 != a3) == false)
            res += 1;
    }
    {
        OInteger i(1);
        mdy::any_type a2(i);
        mdy::any_type a3(i);

        if ((a2 == a3) == false)
            res += 1;
        if ((a2 != a3) == true)
            res += 1;
    }

    {
        std::ostringstream ss;

        mdy::any_type a(1);
        ss << a;

        std::istringstream ss2(ss.str());

        mdy::any_type b;
        ss2 >> b;

        if (a.get_stored() != b.get_stored())
            res += 1;
    };
    {
        std::ostringstream ss;
        oarchive ia(ss);

        mdy::any_type a(std::string("abc"));
        a.serialize(ia.get(), 0);

        std::istringstream ss2(ss.str());
        iarchive ia2(ss2);

        mdy::object b;
        b.serialize(ia2.get(), 0);

        if (a.get_stored() != b)
            res += 1;
    };
    {
        std::ostringstream ss;
        oarchive ia(ss);

        mdy::any_type a(mdy::Any(std::string("abc")));
        a.serialize(ia.get(), 0);

        std::istringstream ss2(ss.str());
        iarchive ia2(ss2);

        mdy::object b;
        b.serialize(ia2.get(), 0);

        if (a.get_stored() != b)
            res += 1;
    };

    {
        mdy::Any a;
        if (a.is_null() != false)
            res += 1;
        //if ((bool)a != false)
        //    res += 1;
        if (a.get_type() != mdy::predefined::type_any())
            res += 1;
        if (a.get().is_zero() != true)
            res += 1;
    }
    {
        mdy::object i(1);
        mdy::Any a(i);

        //if ((bool)a != true)
        //    res += 1;
        if (a.get_type() != mdy::predefined::type_any())
            res += 1;
        if (a.get().get_stored().get_type() != mdy::predefined::type_int())
            res += 1;
    }
    {
        mdy::object i(1);
        mdy::Any a(std::move(i));

        if (a.get_type() != mdy::predefined::type_any())
            res += 1;
        if (a.get().get_stored().get_type() != mdy::predefined::type_int())
            res += 1;
    }
    {
        mdy::object i(1);
        mdy::Any a(i);
        mdy::Any a2 = mdy::Any(mdy::any_type(mdy::object(a)));

        if (a2.get_type() != mdy::predefined::type_any())
            res += 1;
        if (a2.get().get_stored().get_type() != mdy::predefined::type_int())
            res += 1;
    }
    {
        mdy::object i(1);
        mdy::Any a(i);
        mdy::Any a2(a);

        if (a2.get().get_stored().get_type() != mdy::predefined::type_int())
            res += 1;
    }
    {
        mdy::object i(1);
        mdy::Any a(i);
        mdy::Any a2(std::move(a));

        if (a2.get().get_stored().get_type() != mdy::predefined::type_int())
            res += 1;
    }
    {
        mdy::object i(1);
        mdy::Unit a(i);
        mdy::Any a2(a);

        if (a2.get().get_stored().get_type() != mdy::predefined::type_unit())
            res += 1;
    }
    {
        mdy::OInteger a(1);
        mdy::Any a2(a);

        if (a2.get().get_stored().get_type() != mdy::predefined::type_int())
            res += 1;
    }
    {
        mdy::OInteger a(1);
        mdy::Any a2(std::move(a));

        if (a2.get().get_stored().get_type() != mdy::predefined::type_int())
            res += 1;
    }
    {
        Integer i = 1;
        mdy::Any a2(i);

        if (a2.get().get_stored().get_type() != mdy::predefined::type_int())
            res += 1;
    }
    {
        std::string s = "abc";
        mdy::Any a2(s);

        if (a2.get().get_stored().get_type() != mdy::predefined::type_string())
            res += 1;
    }
    {
        std::string s = "abc";
        mdy::Any a2(std::move(s));

        if (a2.get().get_stored().get_type() != mdy::predefined::type_string())
            res += 1;
    }
    {
        Integer i = 1;
        mdy::Any a2(i);
        mdy::Any a3(a2);

        if (a3.get().get_stored().get_type() != mdy::predefined::type_int())
            res += 1;
    }
    {
        Integer i = 1;
        mdy::Any a2(i);
        mdy::Any a3(std::move(a2));

        if (a3.get().get_stored().get_type() != mdy::predefined::type_int())
            res += 1;
    }
    {
        mdy::object i(1);
        mdy::Any a;
        a = i;

        if (a.get().get_stored().get_type() != mdy::predefined::type_int())
            res += 1;
    }
    {
        mdy::object i(1);
        mdy::Any a;
        a = std::move(i);

        if (a.get().get_stored().get_type() != mdy::predefined::type_int())
            res += 1;
    }
    {
        mdy::object i(1);
        mdy::Any a(i);
        mdy::Any a2;
        a2 = mdy::Any(mdy::any_type(mdy::object(a)));

        if (a2.get().get_stored().get_type() != mdy::predefined::type_int())
            res += 1;
    }
    {
        mdy::object i(1);
        mdy::Any a(i);
        mdy::Any a2;
        a2 = a;

        if (a2.get().get_stored().get_type() != mdy::predefined::type_int())
            res += 1;
    }
    {
        mdy::object i(1);
        mdy::Any a(i);
        mdy::Any a2;
        a2 = std::move(a);

        if (a2.get().get_stored().get_type() != mdy::predefined::type_int())
            res += 1;
    }
    {
        mdy::object i(1);
        mdy::Unit a(i);
        mdy::Any a2;
        a2 = a;

        if (a2.get().get_stored().get_type() != mdy::predefined::type_unit())
            res += 1;
    }
    {
        mdy::OInteger a(1);
        mdy::Any a2;
        a2 = a;

        if (a2.get().get_stored().get_type() != mdy::predefined::type_int())
            res += 1;
    }
    {
        mdy::OInteger a(1);
        mdy::Any a2;
        a2 = std::move(a);

        if (a2.get().get_stored().get_type() != mdy::predefined::type_int())
            res += 1;
    }
    {
        Integer i = 1;
        mdy::Any a2;
        a2 = OInteger(i);

        if (a2.get().get_stored().get_type() != mdy::predefined::type_int())
            res += 1;
    }
    {
        std::string s = "abc";
        mdy::Any a2;
        a2 = String(s);

        if (a2.get().get_stored().get_type() != mdy::predefined::type_string())
            res += 1;
    }
    {
        std::string s = "abc";
        mdy::Any a2;
        a2 = String(std::move(s));

        if (a2.get().get_stored().get_type() != mdy::predefined::type_string())
            res += 1;
    }
    {
        Integer i = 1;
        mdy::Any a2(i);
        mdy::Any a3;
        a3 = a2;

        if (a3.get().get_stored().get_type() != mdy::predefined::type_int())
            res += 1;
    }
    {
        Integer i = 1;
        mdy::Any a2(i);
        mdy::Any a3;
        a3 = std::move(a2);

        if (a3.get().get_stored().get_type() != mdy::predefined::type_int())
            res += 1;
    }

    {
        Integer i = 1;
        mdy::Any a2(i);

        if (a2.is_zero() == true)
            res += 1;
        if (a2.is_one() == true)
            res += 1;
    }
    {
        mdy::Any a2;

        if (a2.is_zero() == false)
            res += 1;
    }
    {
        mdy::Any a2 = mdy::Any(mdy::any_type(mdy::object()));

        if (a2.is_zero() == false)
            res += 1;
    }
    {
        mdy::Any a2;
        mdy::Any a3 = a2.clone();

        if (a2.is_zero() == false)
            res += 1;
    }
    {
        Integer i = 1;
        mdy::Any a2(i);
        mdy::Any a3(1.1);

        if (a2.get() == a3.get())
            res += 1;
        if ((a2.get() != a3.get()) == false)
            res += 1;
    }
    {
        Integer i = 1;
        mdy::Any a2(i);
        mdy::Any a3 = mdy::Any(mdy::any_type(OInteger(i)));

        if ((a2.get() == a3.get()) == true)
            res += 1;
        if ((a2.get() != a3.get()) == false)
            res += 1;
    }
    {
        OInteger i(1);
        mdy::Any a2(i);
        mdy::Any a3(i);

        if ((a2.get() == a3.get()) == false)
            res += 1;
        if ((a2.get() != a3.get()) == true)
            res += 1;
    }
    {
        OInteger i(1);
        mdy::object ii(i);
        mdy::Any a2(ii);
        mdy::Any a3(ii);

        if ((a2.get() == a3.get()) == false)
            res += 1;
        if ((a2.get() != a3.get()) == true)
            res += 1;
    }

    {
        OInteger i(1);
        mdy::Any a2;
        mdy::Any a3;

        a2 = i;
        a3 = i;

        if ((a2.get() == a3.get()) == false)
            res += 1;
        if ((a2.get() != a3.get()) == true)
            res += 1;
    }
    {
        OInteger i(1);
        mdy::Any a2;
        mdy::Any a3;

        a2 = mdy::object(i);
        a3 = mdy::object(i);

        if ((a2.get() == a3.get()) == false)
            res += 1;
        if ((a2.get() != a3.get()) == true)
            res += 1;
    }
    {
        std::ostringstream ss;

        mdy::Any a(1);
        ss << a;

        std::istringstream ss2(ss.str());

        mdy::Any b;
        ss2 >> b;

        if (a.get().get_stored() != b.get().get_stored())
            res += 1;
    };
    {
        std::ostringstream ss;
        oarchive ia(ss);

        mdy::Any a(std::string("abc"));
        a.serialize(ia.get(), 0);

        std::istringstream ss2(ss.str());
        iarchive ia2(ss2);

        mdy::Any b;
        b.serialize(ia2.get(), 0);

        if (a.get().get_stored() != b.get().get_stored())
            res += 1;
    };
    {
        std::ostringstream ss;
        oarchive ia(ss);

        mdy::Any a(mdy::Any(std::string("abc")));
        a.serialize(ia.get(), 0);

        std::istringstream ss2(ss.str());
        iarchive ia2(ss2);

        mdy::object b;
        b.serialize(ia2.get(), 0);

        if (a.get().get_stored() != mdy::Any(b).get().get_stored())
            res += 1;
    };
    {
        std::ostringstream ss;
        oarchive ia(ss);

        mdy::Any a(mdy::Any(std::string("abc")));
        a.serialize(ia.get(), 0);

        std::istringstream ss2(ss.str());
        iarchive ia2(ss2);

        mdy::Any b;
        b.serialize(ia2.get(), 0);

        if (a.get().get_stored() != b.get().get_stored())
            res += 1;
    };

    return (res == 0.0) ? true : false;
}

std::set<test_dynamic::type2> test_dynamic::get_missing_cons()
{
    std::set<type2> missing;

    mdy::Type tu = mdy::predefined::type_unit();
    mdy::Type tb = mdy::predefined::type_bool();
    mdy::Type ta = mdy::predefined::type_any();
    mdy::Type ti = mdy::predefined::type_int();
    mdy::Type tr = mdy::predefined::type_real();
    mdy::Type tf = mdy::predefined::type_float();
    mdy::Type tc = mdy::predefined::type_complex();
    mdy::Type tcf = mdy::predefined::type_float_complex();
    mdy::Type ts = mdy::predefined::type_string();

    missing.insert(type2(tb,tu));   //Unit -> bool
    missing.insert(type2(tb,ta));   //Any -> bool
    missing.insert(type2(tb,ti));   //Int -> bool
    missing.insert(type2(tb,tr));   //Real -> bool
    missing.insert(type2(tb,tf));   //Float -> bool
    missing.insert(type2(tb,tc));   //Compl -> bool
    missing.insert(type2(tb,tcf));  //FCompl -> bool
    missing.insert(type2(tb,ts));   //String -> bool

    missing.insert(type2(ti,tu));   //Unit -> Int
    missing.insert(type2(ti,ta));   //Any -> Int
    //missing.insert(type2(ti,tb));   //Bool -> Int
    missing.insert(type2(ti,ts));   //String -> Int

    missing.insert(type2(tr,tu));   //Unit -> Real
    missing.insert(type2(tr,ta));   //Any -> Real
    //missing.insert(type2(tr,tb));   //Bool -> Real
    missing.insert(type2(tr,ts));   //String -> Real

    missing.insert(type2(tf,tu));   //Unit -> Float
    missing.insert(type2(tf,ta));   //Any -> Float
    //missing.insert(type2(tf,tb));   //Bool -> Float
    missing.insert(type2(tf,ts));   //String -> Float

    missing.insert(type2(tc,tu));   //Unit -> Compl
    missing.insert(type2(tc,ta));   //Any -> Compl
    //missing.insert(type2(tc,tb));   //Bool -> Compl
    missing.insert(type2(tc,ts));   //String -> Compl

    missing.insert(type2(tcf,tu));   //Unit -> FCompl
    missing.insert(type2(tcf,ta));   //Any -> FCompl
    //missing.insert(type2(tcf,tb));   //Bool -> FCompl
    missing.insert(type2(tcf,ts));   //String -> FCompl

    missing.insert(type2(ts,tu));   //Unit -> String
    missing.insert(type2(ts,ta));   //Any -> String
    missing.insert(type2(ts,tb));   //Bool -> String
    missing.insert(type2(ts,ti));   //Int -> String
    missing.insert(type2(ts,tr));   //Real -> String
    missing.insert(type2(ts,tf));   //Float -> String
    missing.insert(type2(ts,tc));   //Complex -> String
    missing.insert(type2(ts,tcf));  //Fcompl -> String

    missing.insert(type2(ti, tr));  //Real->Int
    missing.insert(type2(ti, tf));  //Float->Int
    missing.insert(type2(ti, tc));  //Compl->Int
    missing.insert(type2(ti, tcf)); //FCompl->Int
    missing.insert(type2(tr, tc));  //Complex->Real
    missing.insert(type2(tr, tcf)); //FComplex->Real
    missing.insert(type2(tf, tc));  //Complex->Float
    missing.insert(type2(tf, tcf)); //FComplex->Float

    return missing;
}
std::set<test_dynamic::type2> test_dynamic::get_missing_assign()
{
    std::set<test_dynamic::type2> missing_as = get_missing_cons();

    mdy::Type tu = mdy::predefined::type_unit();
    mdy::Type tb = mdy::predefined::type_bool();
    mdy::Type ta = mdy::predefined::type_any();
    mdy::Type ti = mdy::predefined::type_int();
    mdy::Type tr = mdy::predefined::type_real();
    mdy::Type tf = mdy::predefined::type_float();
    mdy::Type tc = mdy::predefined::type_complex();
    mdy::Type tcf = mdy::predefined::type_float_complex();
    mdy::Type ts = mdy::predefined::type_string();

    missing_as.insert(type2(ti,tr));
    missing_as.insert(type2(ti,tf));
    missing_as.insert(type2(ti,tc));
    missing_as.insert(type2(ti,tcf));
    missing_as.insert(type2(tr,tc));
    missing_as.insert(type2(tr,tcf));
    missing_as.insert(type2(tf,tc));
    missing_as.insert(type2(tf,tcf));

    missing_as.insert(type2(ti,tb));
    missing_as.insert(type2(tr,tb));
    missing_as.insert(type2(tf,tb));
    missing_as.insert(type2(tc,tb));
    missing_as.insert(type2(tcf,tb));

    return missing_as;
};

bool test_dynamic::test_cons_impl()
{
    mdy::Type tu = mdy::predefined::type_unit();
    mdy::Type tb = mdy::predefined::type_bool();
    mdy::Type ta = mdy::predefined::type_any();
    mdy::Type ti = mdy::predefined::type_int();
    mdy::Type tr = mdy::predefined::type_real();
    mdy::Type tf = mdy::predefined::type_float();
    mdy::Type tc = mdy::predefined::type_complex();
    mdy::Type tcf = mdy::predefined::type_float_complex();
    mdy::Type ts = mdy::predefined::type_string();

    std::vector<mdy::Type> types;
    types.push_back(tu);
    types.push_back(ta);
    types.push_back(tb);
    types.push_back(ti);
    types.push_back(tr);
    types.push_back(tf);
    types.push_back(tc);
    types.push_back(tcf);
    types.push_back(ts);

    using type2 = std::pair<mdy::Type, mdy::Type>;

    std::set<type2> missing     = get_missing_cons();
    std::set<type2> missing_as  = get_missing_assign();

    double res = 0;
    for (const auto& pos1 : types)
    for (const auto& pos2 : types)
    {
        res += test_cons_1(pos1, pos2, missing, missing_as);
    };

    return (res == 0) ? true : false;
}

double test_dynamic::test_cons_1(mdy::Type t1, mdy::Type t2, const type2_set& missing, 
                                 const type2_set& missing_as)
{
    mdy::object rhs(t2);
    mdy::object lhs(t1);

    double res = 0;

    try
    {        
        mdy::object t(t1, rhs);

        if (missing.find(type2(t1, t2)) != missing.end())
            res += 1;
    }
    catch(std::exception& ex)
    {
        if (missing.find(type2(t1, t2)) == missing.end())
        {
            disp(ex.what());
            res += 1;
        };
    }

    try
    {
        lhs = rhs;

        if (missing_as.find(type2(t1, t2)) != missing_as.end())
            res += 1;
    }
    catch(std::exception& ex)
    {
        if (missing_as.find(type2(t1, t2)) == missing_as.end())
        {
            disp(ex.what());
            res += 1;
        };
    }

    return res;
};

bool test_dynamic::test_cons_val_impl()
{
    std::vector<Scalar> scalars1;
    std::vector<Scalar> scalars2;

    Integer N   = 1000;
    rand_scalars::make(scalars1, N, true);
    rand_scalars::make(scalars2, N, true);

    double res = 0.0;

    for (Integer i = 0; i < 1000; ++i)
        res     += test_cons_val_1(scalars1[i],scalars2[i], i);

    return (res == 0) ? true : false;
}

double test_dynamic::test_cons_val_1(const Scalar& s1, const Scalar& s2, Integer code)
{
    eval_cons_val test1(code);
    double res1 = test1.make(s1, s2);

    eval_assign_val test2(code);
    double res2 = test2.make(s1, s2);
    return res1 + res2;
};

bool test_dynamic::test_func_compare_impl()
{
    mdy::Type tu = mdy::predefined::type_unit();
    mdy::Type tb = mdy::predefined::type_bool();
    mdy::Type ta = mdy::predefined::type_any();
    mdy::Type ti = mdy::predefined::type_int();
    mdy::Type tr = mdy::predefined::type_real();
    mdy::Type tf = mdy::predefined::type_float();
    mdy::Type tc = mdy::predefined::type_complex();
    mdy::Type tcf = mdy::predefined::type_float_complex();
    mdy::Type ts = mdy::predefined::type_string();

    std::vector<mdy::Type> types;
    //types.push_back(tu);
    //types.push_back(ta);
    types.push_back(tb);
    types.push_back(ti);
    types.push_back(tr);
    types.push_back(tf);
    types.push_back(tc);
    types.push_back(tcf);
    types.push_back(ts);

    using type2 = std::pair<mdy::Type, mdy::Type>;

    std::set<type2> missing;

    missing.insert(type2(tb,ti));
    missing.insert(type2(tb,tr));
    missing.insert(type2(tb,tf));
    missing.insert(type2(tb,tc));
    missing.insert(type2(tb,tcf));
    missing.insert(type2(tb,ts));

    missing.insert(type2(ti,tb));
    missing.insert(type2(tr,tb));
    missing.insert(type2(tf,tb));
    missing.insert(type2(tc,tb));
    missing.insert(type2(tcf,tb));
    missing.insert(type2(ts,tb));

    missing.insert(type2(ts,ti));
    missing.insert(type2(ts,tr));
    missing.insert(type2(ts,tf));
    missing.insert(type2(ts,tc));
    missing.insert(type2(ts,tcf));

    missing.insert(type2(ti,ts));
    missing.insert(type2(tr,ts));
    missing.insert(type2(tf,ts));
    missing.insert(type2(tc,ts));
    missing.insert(type2(tcf,ts));

    double res = 0;
    for (const auto& pos1 : types)
    for (const auto& pos2 : types)
    {
        res += test_compare_1(pos1, pos2, missing);
    };

    return (res == 0) ? true : false;
}

double test_dynamic::test_compare_1(mdy::Type t1, mdy::Type t2, const type2_set& missing)
{
    mdy::object rhs(t2);
    mdy::object lhs(t1);

    double res = 0;

    try
    {        
        bool res1 = (bool)(lhs == rhs);
        bool res2 = (bool)(lhs != rhs);
        bool res3 = (bool)(lhs < rhs);
        bool res4 = (bool)(lhs > rhs);
        bool res5 = (bool)(lhs <= rhs);
        bool res6 = (bool)(lhs >= rhs);

        (void)res1;
        (void)res2;
        (void)res3;
        (void)res4;
        (void)res5;
        (void)res6;

        if (missing.find(type2(t1, t2)) != missing.end())
            res += 1;
    }
    catch(std::exception& ex)
    {
        if (missing.find(type2(t1, t2)) == missing.end())
        {
            disp(ex.what());
            res += 1;
        };
    }

    return res;
};

bool test_dynamic::test_object_func_impl()
{
    mdy::Type tu = mdy::predefined::type_unit();
    mdy::Type tb = mdy::predefined::type_bool();
    mdy::Type ta = mdy::predefined::type_any();
    mdy::Type ti = mdy::predefined::type_int();
    mdy::Type tr = mdy::predefined::type_real();
    mdy::Type tf = mdy::predefined::type_float();
    mdy::Type tc = mdy::predefined::type_complex();
    mdy::Type tcf = mdy::predefined::type_float_complex();
    mdy::Type ts = mdy::predefined::type_string();

    std::vector<mdy::Type> types;
    //types.push_back(tu);
    //types.push_back(ta);
    types.push_back(mdy::Type());
    types.push_back(tb);
    types.push_back(ti);
    types.push_back(tr);
    types.push_back(tf);
    types.push_back(tc);
    types.push_back(tcf);
    types.push_back(ts);

    double res = 0;
    for (const auto& pos1 : types)
        res += test_func_1(pos1);

    return (res == 0) ? true : false;
}

double test_dynamic::test_func_1(mdy::Type t1)
{
    mdy::object lhs(t1);

    double res = 0;

    try
    {        
        (bool)lhs;
        !lhs;
        lhs.clone();
    }
    catch(std::exception& ex)
    {
        disp(ex.what());
        res += 1;
    }

    if (t1 != mdy::Type())
    {
        if (lhs.is_zero() == false)
            res += 1;
        if (lhs.is_null() == true)
            res += 1;
    };

    return res;
};

bool test_dynamic::test_object_func_val_impl()
{
    std::vector<Scalar> scalars1;

    Integer N   = 1000;
    rand_scalars::make(scalars1, N, true);

    double res = 0.0;

    for (Integer i = 0; i < 1000; ++i)
        res     += test_func_1_val(scalars1[i], i);

    return (res == 0) ? true : false;
}

double test_dynamic::test_func_1_val(const Scalar& s1, Integer code)
{
    eval_object_func_val test1(code);
    double res1 = test1.make(s1);

    return res1;
};

bool test_dynamic::test_func_uminus_impl()
{
    mdy::Type tu = mdy::predefined::type_unit();
    mdy::Type tb = mdy::predefined::type_bool();
    mdy::Type ta = mdy::predefined::type_any();
    mdy::Type ti = mdy::predefined::type_int();
    mdy::Type tr = mdy::predefined::type_real();
    mdy::Type tf = mdy::predefined::type_float();
    mdy::Type tc = mdy::predefined::type_complex();
    mdy::Type tcf = mdy::predefined::type_float_complex();
    mdy::Type ts = mdy::predefined::type_string();

    std::vector<mdy::Type> types;
    //types.push_back(tu);
    //types.push_back(ta);
    //types.push_back(mdy::Type());
    //types.push_back(tb);
    types.push_back(ti);
    types.push_back(tr);
    types.push_back(tf);
    types.push_back(tc);
    types.push_back(tcf);
    //types.push_back(ts);

    double res = 0;
    for (const auto& pos1 : types)
    {
        res += test_func_uminus_1(pos1);
    };

    return (res == 0) ? true : false;
}

double test_dynamic::test_func_uminus_1(mdy::Type t1)
{
    mdy::object lhs(t1);

    double res = 0;

    try
    {        
        mdy::object r = -lhs;

        if (r.get_type() != lhs.get_type())
            res += 1;
    }
    catch(std::exception& ex)
    {
        disp(ex.what());
        res += 1;
    }

    return res;
};

bool test_dynamic::test_func_reim_impl()
{
    mdy::Type tu = mdy::predefined::type_unit();
    mdy::Type tb = mdy::predefined::type_bool();
    mdy::Type ta = mdy::predefined::type_any();
    mdy::Type ti = mdy::predefined::type_int();
    mdy::Type tr = mdy::predefined::type_real();
    mdy::Type tf = mdy::predefined::type_float();
    mdy::Type tc = mdy::predefined::type_complex();
    mdy::Type tcf = mdy::predefined::type_float_complex();
    mdy::Type ts = mdy::predefined::type_string();

    std::vector<mdy::Type> types;
    //types.push_back(tu);
    //types.push_back(ta);
    //types.push_back(mdy::Type());
    //types.push_back(tb);
    types.push_back(ti);
    types.push_back(tr);
    types.push_back(tf);
    types.push_back(tc);
    types.push_back(tcf);
    //types.push_back(ts);

    double res = 0;
    for (const auto& pos1 : types)
    {
        res += test_func_reim_1(pos1);
    };

    return (res == 0) ? true : false;
}

static bool is_complex_type(mdy::Type t)
{
    return t == mdy::predefined::type_complex() || t == mdy::predefined::type_float_complex();
}

double test_dynamic::test_func_reim_1(mdy::Type t1)
{
    mdy::object lhs(t1);

    double res = 0;

    try
    {        
        mdy::object res1 = real(lhs);
        mdy::object res2 = imag(lhs);

        if (res1.get_type() != res2.get_type())
            res += 1;

        if (res1.get_type() != lhs.get_type() && is_complex_type(t1) == false)
            res += 1;
    }
    catch(std::exception& ex)
    {
        disp(ex.what());
        res += 1;
    }

    return res;
};

template<class T>
double test_dynamic::test_one_type(bool has_one)
{
    using OT = mdy::object_type<T>;
    double res = 0;

    if (OT::has_one != has_one)
        res += 1;

    if (OT::is_zero(T()) == false)
        res += 1;

    OT one = OT::make_one();
    if (OT::is_one(one.get()) == false && has_one == true)
        res += 1;
    if (one.is_one() == false && has_one == true)
        res += 1;

    return res;
};

bool test_dynamic::test_one_impl()
{
    mdy::Type t_unit     = mdy::predefined::type_unit();
    mdy::Type t_any      = mdy::predefined::type_any();
    mdy::Type t_string   = mdy::predefined::type_string();
    mdy::Type t_null     = mdy::Type();

    double res          = 0.0;    

    std::vector<mdy::Type> types_0{t_unit, t_any, t_string, t_null};
    
    for (const auto& pos : types_0)
    {
        if (mdy::operations::has_one(pos) == true)
            res             += 1;

        try
        {
            mdy::object::make_one(pos);
            res += 1;
        }
        catch(...){};
    };    

    mdy::Type t_int      = mdy::predefined::type_int();
    mdy::Type t_float    = mdy::predefined::type_float();
    mdy::Type t_real     = mdy::predefined::type_real();
    mdy::Type t_compl    = mdy::predefined::type_complex();
    mdy::Type t_fcompl   = mdy::predefined::type_float_complex();
    mdy::Type t_mi       = MP_int::get_static_type();
    mdy::Type t_mf       = MP_float::get_static_type();
    mdy::Type t_mc       = MP_complex::get_static_type();
    mdy::Type t_mr       = MP_rational::get_static_type();

    std::vector<mdy::Type> types_1{t_int,t_float,t_real,t_compl,t_fcompl,t_mi,t_mf,t_mc,t_mr};

    mdy::object o1(1);

    for (const auto& pos : types_1)
    {
        if (mdy::operations::has_one(pos) == false)
            res             += 1;

        try
        {
            mdy::object obj = mdy::object::make_one(pos);
            if (obj != o1)
                res         += 1;
            if (obj.is_one() == false)
                res         += 1;
        }
        catch(...)
        {
            res += 1;
        };
    };  

    res += test_one_type<Integer>(true);
    res += test_one_type<Float>(true);
    res += test_one_type<Real>(true);
    res += test_one_type<Complex>(true);
    res += test_one_type<Float_complex>(true);
    res += test_one_type<mp_int>(true);
    res += test_one_type<mp_float>(true);
    res += test_one_type<mp_complex>(true);
    res += test_one_type<mp_rational>(true);
    res += test_one_type<std::string>(false);
    res += test_one_type<mdy::any_type>(false);

    return res == 0.0;
};

#pragma warning(push)
#pragma warning (disable : 4800) // forcing value to bool 'true' or 'false' (performance warning)

template<class T>
void test_dynamic::test_compile_func_unary()
{
    T val;
    Integer iexp= 0;

    bool x  = cast_bool(val);
    x       = !val;
    T y     = -val;
    real(val);
    imag(val);
    is_zero(val);
    is_one(val);
    fpclassify(val);
    is_nan(val);
    is_finite(val);
    is_inf(val);
    is_regular(val);
    is_normal(val);
    is_int(val);
    is_real(val);
    neg(val);    
    is_false(val);
    is_true(val);
    uminus(val);
    y = +val;
    inv(val);
    invs(val);
    abs(val);
    abs2(val);
    arg(val);
    angle(val);
    conj(val);
    sqrt(val);
    sqrt_c(val);    
    exp(val);
    exp2(val);
    pow2(val);
    exp10(val);
    pow10(val);
    expm1(val);
    expi(val);
    ldexp(val, iexp);
    scalbn(val, iexp);
    log(val);
    log_c(val);
    log1p(val);
    log1p_c(val);
    log2(val);
    log2_c(val);
    log10(val);
    log10_c(val);
    logb(val);
    ilogb(val);
    sin(val);
    cos(val);
    tan(val);
    cot(val);
    sec(val);
    csc(val);
    sinh(val);
    cosh(val);
    tanh(val);
    coth(val);
    sech(val);
    csch(val);
    asin(val);
    acos(val);
    atan(val);
    asinh(val);
    acosh(val);
    atanh(val);
    asin_c(val);
    acos_c(val);
    acosh_c(val);
    atanh_c(val);
    eps(val);
    floor(val);
    ceil(val);
    trunc(val);
    fix(val);
    round(val);
    sign(val);        
};

template<class T>
void test_dynamic::test_compile_func_unary_notcompl()
{
    using TF    = typename md::unify_types<T,Float>::type;
    T val;
    TF sec;
    Integer isec;

    frexp(val,isec);
    modf(val,sec);
    cbrt(val);
    ifloor(val);
    iceil(val);
    itrunc(val);
    ifix(val);
    iround(val);
    isign(val);
    signbit(val);
};
template<class T>
void test_dynamic::test_compile_func_unary_notmp()
{
    T val;

    sqrt1pm1(val);
    sqrt1pm1_c(val);

    acot(val);
    asec(val);
    acsc(val);
    acoth(val);
    asech(val);
    acsch(val);
    asec_c(val);
    acsc_c(val);
    acoth_c(val);
    asech_c(val);
};

template<class T>
void test_dynamic::test_compile_func_unary_str()
{
    T val;

    bool x  = cast_bool(val);
    x       = !val;
    is_zero(val);
    is_one(val);
    neg(val);    
    is_false(val);
    is_true(val);

    T y = +val;
};

template<class T1, class T2>
void test_dynamic::test_compile_func_binary()
{
    T1 x = T1();
    T2 y = T2(1);

    bool res;

    res = (bool)(x == y);
    res = (bool)(x != y);
    res = (bool)(x <= y);
    res = (bool)(x >= y);
    res = (bool)(x > y);
    res = (bool)(x < y);
    auto v1 = x + y;
    auto v2 = x - y;
    auto v3 = x * y;
    auto v4 = x / y;

    (void)v1;
    (void)v2;
    (void)v3;
    (void)v4;

    idiv(x,y);
    plus(x,y);
    minus(x,y);
    mul(x,y);
    div(x,y);
    div_0(x,y);
    div_1(x,y);
    hypot(x,y);
    pow(x,y);
    pow_c(x,y);
    x&&y;
    x||y;
    op_and(x,y);
    op_or(x,y);
    op_xor(x,y);
    elem_and(x,y);
    elem_or(x,y);
    elem_xor(x,y);
    eeq(x,y);
    neq(x,y);
    eeq_nan(x,y);
    neq_nan(x,y);
    leq(x,y);
    geq(x,y);
    lt(x,y);
    gt(x,y);
    max(x,y);
    min(x,y);
};

template<class T1, class T2>
void test_dynamic::test_compile_func_binary_notcompl()
{
    T1 x;
    T2 y;

    mod(x,y);
    rem(x,y);
    atan2(x,y);
    copysign(x,y);
    nextafter(x,y);
    float_distance(x,y);
    nextabove(x);
    nextbelow(x);
    fdim(x,y);
};

template<class T1, class T2>
void test_dynamic::test_compile_func_binary_notmp()
{
    T1 x;
    T2 y;

    powm1(x,y);
};

template<class T1, class T2>
void test_dynamic::test_compile_func_binary_str()
{
    T1 x;
    T2 y;

    bool res;

    res = (bool)(x == y);
    res = (bool)(x != y);
    res = (bool)(x <= y);
    res = (bool)(x >= y);
    res = (bool)(x > y);
    res = (bool)(x < y);

    auto v1 = x + y;
    plus(x,y);
    x&&y;
    x||y;
    x&y;
    x|y;
    op_and(x,y);
    op_or(x,y);
    op_xor(x,y);
    elem_and(x,y);
    elem_or(x,y);
    elem_xor(x,y);
    eeq(x,y);
    neq(x,y);
    eeq_nan(x,y);
    neq_nan(x,y);
    leq(x,y);
    geq(x,y);
    lt(x,y);
    gt(x,y);

    auto r1 = max(x,y);
    auto r2 = min(x,y);
};

#pragma warning(pop)

void test_dynamic::test_compite_func()
{
    test_compile_func_unary_str<std::string>();
    test_compile_func_unary<Integer>();
    test_compile_func_unary<Real>();
    test_compile_func_unary<Float>();
    test_compile_func_unary<Complex>();
    test_compile_func_unary<Float_complex>();
    test_compile_func_unary<mp_int>();
    test_compile_func_unary<mp_float>();
    test_compile_func_unary<mp_complex>();
    test_compile_func_unary<mp_rational>();

    test_compile_func_unary_notcompl<Integer>();
    test_compile_func_unary_notcompl<Real>();
    test_compile_func_unary_notcompl<Float>();
    test_compile_func_unary_notcompl<mp_int>();
    test_compile_func_unary_notcompl<mp_float>();
    test_compile_func_unary_notcompl<mp_rational>();
    test_compile_func_unary_notmp<Integer>();
    test_compile_func_unary_notmp<Real>();
    test_compile_func_unary_notmp<Float>();
    test_compile_func_unary_notmp<Float_complex>();
    test_compile_func_unary_notmp<Complex>();

    //test_compile_func_unary<String>();
    test_compile_func_unary<OInteger>();
    test_compile_func_unary<OReal>();
    test_compile_func_unary<OFloat>();
    test_compile_func_unary<OComplex>();
    test_compile_func_unary<OFloat_complex>();
    test_compile_func_unary<MP_int>();
    test_compile_func_unary<MP_float>();
    test_compile_func_unary<MP_complex>();
    test_compile_func_unary<MP_rational>();
    test_compile_func_unary_notcompl<OInteger>();
    test_compile_func_unary_notcompl<OReal>();
    test_compile_func_unary_notcompl<OFloat>();
    test_compile_func_unary_notcompl<MP_int>();
    test_compile_func_unary_notcompl<MP_float>();
    test_compile_func_unary_notcompl<MP_rational>();
    test_compile_func_unary_notmp<OInteger>();
    test_compile_func_unary_notmp<OReal>();
    test_compile_func_unary_notmp<OFloat>();
    test_compile_func_unary_notmp<OFloat_complex>();
    test_compile_func_unary_notmp<OComplex>();

    test_compile_func_binary_str<std::string,std::string>();
    test_compile_func_binary<Integer,Integer>();
    test_compile_func_binary<Real,Real>();
    test_compile_func_binary<Float,Float>();
    test_compile_func_binary<Complex,Complex>();
    test_compile_func_binary<Float_complex,Float_complex>();
    test_compile_func_binary<mp_int,mp_int>();
    test_compile_func_binary<mp_float,mp_float>();
    test_compile_func_binary<mp_rational,mp_rational>();
    test_compile_func_binary<mp_complex,mp_complex>();
    test_compile_func_binary_notcompl<Integer,Integer>();
    test_compile_func_binary_notcompl<Real,Real>();
    test_compile_func_binary_notcompl<Float,Float>();
    test_compile_func_binary_notcompl<mp_int,mp_int>();
    test_compile_func_binary_notcompl<mp_float,mp_float>();
    test_compile_func_binary_notcompl<mp_rational,mp_rational>();
    test_compile_func_binary_notmp<Integer,Integer>();
    test_compile_func_binary_notmp<Real,Real>();
    test_compile_func_binary_notmp<Float,Float>();

    test_compile_func_binary_str<String,std::string>();
    test_compile_func_binary<OInteger,Integer>();
    test_compile_func_binary<OReal,Real>();
    test_compile_func_binary<OFloat,Float>();
    test_compile_func_binary<OComplex,Complex>();
    test_compile_func_binary<OFloat_complex,Float_complex>();
    test_compile_func_binary<MP_int,mp_int>();
    test_compile_func_binary<MP_float,mp_float>();
    test_compile_func_binary<MP_rational,mp_rational>();
    test_compile_func_binary<MP_complex,mp_complex>();
    test_compile_func_binary_notcompl<OInteger,Integer>();
    test_compile_func_binary_notcompl<OReal,Real>();
    test_compile_func_binary_notcompl<OFloat,Float>();
    test_compile_func_binary_notcompl<MP_int,mp_int>();
    test_compile_func_binary_notcompl<MP_float,mp_float>();
    test_compile_func_binary_notcompl<MP_rational,mp_rational>();
    test_compile_func_binary_notmp<OInteger,Integer>();
    test_compile_func_binary_notmp<OReal,Real>();
    test_compile_func_binary_notmp<OFloat,Float>();

    test_compile_func_binary_str<std::string,String>();
    test_compile_func_binary<Integer,OInteger>();
    test_compile_func_binary<Real,OReal>();
    test_compile_func_binary<Float,OFloat>();
    test_compile_func_binary<Complex,OComplex>();
    test_compile_func_binary<Float_complex,OFloat_complex>();
    test_compile_func_binary<mp_int,MP_int>();
    test_compile_func_binary<mp_float,MP_float>();
    test_compile_func_binary<mp_rational,MP_rational>();
    test_compile_func_binary<mp_complex,MP_complex>();
    test_compile_func_binary_notcompl<Integer,OInteger>();
    test_compile_func_binary_notcompl<Real,OReal>();
    test_compile_func_binary_notcompl<Float,OFloat>();
    test_compile_func_binary_notcompl<mp_int,MP_int>();
    test_compile_func_binary_notcompl<mp_float,MP_float>();
    test_compile_func_binary_notcompl<mp_rational,MP_rational>();
    test_compile_func_binary_notmp<Integer,OInteger>();
    test_compile_func_binary_notmp<Real,OReal>();
    test_compile_func_binary_notmp<Float,OFloat>();

    test_compile_func_binary_str<String,String>();
    test_compile_func_binary<OInteger,OInteger>();
    test_compile_func_binary<OReal,OReal>();
    test_compile_func_binary<OFloat,OFloat>();
    test_compile_func_binary<OComplex,OComplex>();
    test_compile_func_binary<OFloat_complex,OFloat_complex>();
    test_compile_func_binary<MP_int,MP_int>();
    test_compile_func_binary<MP_float,MP_float>();
    test_compile_func_binary<MP_rational,MP_rational>();
    test_compile_func_binary<MP_complex,MP_complex>();
    test_compile_func_binary_notcompl<OInteger,OInteger>();
    test_compile_func_binary_notcompl<OReal,OReal>();
    test_compile_func_binary_notcompl<OFloat,OFloat>();
    test_compile_func_binary_notcompl<MP_int,MP_int>();
    test_compile_func_binary_notcompl<MP_float,MP_float>();
    test_compile_func_binary_notcompl<MP_rational,MP_rational>();
    test_compile_func_binary_notmp<OInteger,OInteger>();
    test_compile_func_binary_notmp<OReal,OReal>();
    test_compile_func_binary_notmp<OFloat,OFloat>();
};

}};
