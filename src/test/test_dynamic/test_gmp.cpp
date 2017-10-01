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

#include "test_gmp.h"
#include "matcl-mp/matcl_mp.h"
#include "matcl-core/IO/archive.h"
#include "matcl-core/lib_functions/constants.h"
#include "matcl-scalar/objects/object_functions.h"
#include "matcl-scalar/objects/typed_object_functions.h"
#include "rand_scalars.h"
#include "utils.h"
#include "matcl-scalar/lib_functions/func_binary.h"
#include "matcl-scalar/lib_functions/manip.h"
#include "matcl-scalar/lib_functions/func_forwarding.h"

#include "matcl-scalar/matcl_scalar.h"

#include <iostream>

namespace matcl
{

static bool is_nan(bool)
{
    return false;
};

static std::ostream& operator<<(std::ostream& os, fp_type t)
{
    os << (Integer)t;
    return os;
};

};

namespace matcl { namespace test
{

void test_gmp()
{
    gmp_tester test;
    test.make();
};

void gmp_tester::make()
{   
    test_inv();  

    test_ldexp_p3();
    test_ldexp_m3();
    test_scalbn_p3();
    test_scalbn_m3();
    test_frexp();
    test_modf_frac();
    test_modf_int();
    test_logb();
    test_ilogb();

    test_sqrt(); 
    test_cbrt(); 
    test_sqrt_c(); 
    test_exp(); 
    test_expm1(); 
    test_expi(); 
    test_exp2();
    test_exp10();
    test_log(); 
    test_log1p();    
    test_log2(); 
    test_log10(); 
    test_log_c();
    test_log1p_c();
    test_log2_c(); 
    test_log10_c(); 

    test_sin();
    test_cos();
    test_tan();    
    test_cot();
    test_sec();
    test_csc();
    test_sinh();
    test_cosh();
    test_tanh();
    test_coth();
    test_sech();
    test_csch();
    test_asin();
    test_asin_c();
    test_acos();
    test_acos_c();
    test_atan();
    test_asinh();
    test_acosh();
    test_acosh_c();
    test_atanh();
    test_atanh_c();
    test_inv();  
    test_invs();  
    test_floor();
    test_ceil();
    test_round();
    test_trunc();
    test_ifloor();
    test_iceil();
    test_iround();
    test_itrunc();
    test_sign();
    test_isign();
    test_fpclassify();

    test_func(test_mp_int(), "mp_int");
    test_func(test_mp_float(), "mp_float");
    test_func(test_mp_complex(), "mp_complex");
    test_func(test_mp_rational(), "mp_rational");
    
    test_uminus();
    test_reim();
    test_is();
    test_next();
    test_signbit();    
    test_eps();

    test_combinatorics();
    test_combinatorics2();
};

void gmp_tester::test_reim()
{
    bool res;

    try
    {
        res = test_reim_impl();
    }
    catch(std::exception& ex)
    {
        std::cout << "mp reim EXCEPTION: " << ex.what() << "\n";
    };

    if (res == true)
        std::cout << "mp reim: " << "ok" << "\n";
    else
        std::cout << "mp reim: " << "FAILED" << "\n";
};

void gmp_tester::test_uminus()
{
    bool res;

    try
    {
        res = test_uminus_impl();
    }
    catch(std::exception& ex)
    {
        std::cout << "mp uminus EXCEPTION: " << ex.what() << "\n";
    };

    if (res == true)
        std::cout << "mp uminus: " << "ok" << "\n";
    else
        std::cout << "mp uminus: " << "FAILED" << "\n";
};

void gmp_tester::test_is()
{
    bool res;

    try
    {
        res = test_is_impl();
    }
    catch(std::exception& ex)
    {
        std::cout << "mp is EXCEPTION: " << ex.what() << "\n";
    };

    if (res == true)
        std::cout << "mp is: " << "ok" << "\n";
    else
        std::cout << "mp is: " << "FAILED" << "\n";
};
void gmp_tester::test_next()
{
    bool res;

    try
    {
        res = test_next_impl();
    }
    catch(std::exception& ex)
    {
        std::cout << "mp next EXCEPTION: " << ex.what() << "\n";
    };

    if (res == true)
        std::cout << "mp next: " << "ok" << "\n";
    else
        std::cout << "mp next: " << "FAILED" << "\n";
};
void gmp_tester::test_signbit()
{
    bool res;

    try
    {
        res = test_sign_impl();
    }
    catch(std::exception& ex)
    {
        std::cout << "mp sign EXCEPTION: " << ex.what() << "\n";
    };

    if (res == true)
        std::cout << "mp sign: " << "ok" << "\n";
    else
        std::cout << "mp sign: " << "FAILED" << "\n";
};
void gmp_tester::test_eps()
{
    bool res;

    try
    {
        res = test_eps_impl();
    }
    catch(std::exception& ex)
    {
        std::cout << "mp eps EXCEPTION: " << ex.what() << "\n";
    };

    if (res == true)
        std::cout << "mp eps: " << "ok" << "\n";
    else
        std::cout << "mp eps: " << "FAILED" << "\n";
};

void gmp_tester::test_func(double dif, const std::string& test_name)
{
    if (dif == 0.0)
        std::cout << test_name << ": ok" << "\n";
    else
        std::cout << test_name << ": FAILED" << "\n";
};

double gmp_tester::test_mp_int()
{
    double res = 0;

    try
    {
        {
            mp_int z;
            if (z != 0)
                res += 1.0;
        };
        {
            int val = 1;

            mp_int z(val);
            if (z != val)
                res += 1.0;
        };

        {
            int val = 1;

            mp_int z = val;
            if (z != val)
                res += 1.0;
        };

        {
            mp_int val = -1;

            mp_int z(val);
            if (z != val)
                res += 1.0;
        };
        {
            mp_int val = -1;
            mp_int val2 = val;

            mp_int z(std::move(val2));
            if (z != val)
                res += 1.0;
        };
        {
            mp_int val = -1;

            mp_int z = val;
            if (z != val)
                res += 1.0;
        };
        {
            mp_int val = -1;

            mp_int z = mp_int(val);
            if (z != val)
                res += 1.0;
        };

        {
            mp_int z("123");
            if (z != 123)
                res += 1.0;
        };
        {
            mp_int z("-123");
            if (z != -123)
                res += 1.0;
        };
        {
            mp_int z(0);
            if (z)
                res += 1.0;
        };
        {
            mp_int z(-1);
            if (!z)
                res += 1.0;
        };

        {
            int val = 123;
            mp_int z(val);
            mp_float f = z.cast_mp_float();
            mp_rational q = z.cast_mp_rational();
            mp_complex zz = z.cast_mp_complex();

            if (f != val)
                res += 1;
            if (q != val)
                res += 1;
            if (zz != val)
                res += 1;
        };
        {
            std::ostringstream ss;
            oarchive ia(ss);

            mp_int z(-1234567);
            z.serialize(ia.get(), 0);

            std::istringstream ss2(ss.str());
            iarchive ia2(ss2);

            mp_int z2;
            z2.serialize(ia2.get(), 0);

            if (z != z2)
                res += 1;
        };

        {
            mp_int z(-123);
            if (z.to_string() != "-123")
                res += 1.0;
        };

        return res;
    }
    catch(...)
    {
        return 1.0;
    };
};

double gmp_tester::test_mp_float()
{
    double res = 0;

    precision old_prec = mp_float::get_default_precision();
    precision default_prec = precision(118);
    mp_float::set_default_precision(default_prec);

    try
    {
        //constructors and conversions
        {
            mp_float z;              
            if (z != 0.)
                res += 1.0;
            if (z.get_precision() != default_prec)
                res += 1;
        };
        {
            Integer val = 1;

            mp_float z(val);
            if (z != val)
                res += 1.0;
            if (z.get_precision() != default_prec)
                res += 1;
        };
        {
            double val = 1;

            mp_float z(val);
            if (z != val)
                res += 1.0;
            if (z.get_precision() != default_prec)
                res += 1;
        };
        {
            float val = 1;

            mp_float z(val);
            if (z != val)
                res += 1.0;
            if (z.get_precision() != default_prec)
                res += 1;
        };

        {
            Integer val = 1;
            precision prec = precision(500);
            mp_float z(val, prec);

            if (z != val)
                res += 1.0;
            if (z.get_precision() != prec)
                res += 1;
        };
        {
            Integer val = 1;
            precision prec = precision(0);
            mp_float z(val, prec);

            if (z != val)
                res += 1.0;
            if (z.get_precision() != default_prec)
                res += 1;
        };
        {
            double val = 1;
            precision prec = precision(500);

            mp_float z(val, prec);
            if (z != val)
                res += 1.0;
            if (z.get_precision() != prec)
                res += 1;
        };
        {
            double val = 1;
            precision prec = precision(0);

            mp_float z(val, prec);
            if (z != val)
                res += 1.0;
            if (z.get_precision() != default_prec)
                res += 1;
        };
        {
            float val   = 1;
            precision prec = precision(500);

            mp_float z(val, prec);
            if (z != val)
                res += 1.0;
            if (z.get_precision() != prec)
                res += 1;
        };
        {
            float val   = 1;
            precision prec = precision(0);

            mp_float z(val, prec);
            if (z != val)
                res += 1.0;
            if (z.get_precision() != default_prec)
                res += 1;
        };
        {
            double val = matcl::constants::nan();

            mp_float z(val);
            if (is_nan(z) == false)
                res += 1.0;
            if (z.get_precision() != default_prec)
                res += 1;
        };
        {
            float val = matcl::constants::f_nan();

            mp_float z(val);
            if (is_nan(z) == false)
                res += 1.0;
            if (z.get_precision() != default_prec)
                res += 1;
        };
        {
            double val = matcl::constants::inf();

            mp_float z(val);
            if (z != val)
                res += 1.0;
            if (z.get_precision() != default_prec)
                res += 1;
        };
        {
            float val = matcl::constants::f_inf();

            mp_float z(val);
            if (z != val)
                res += 1.0;
            if (z.get_precision() != default_prec)
                res += 1;
        };
        {
            double val = -matcl::constants::inf();

            mp_float z(val);
            if (z != val)
                res += 1.0;
            if (z.get_precision() != default_prec)
                res += 1;
        };
        {
            float val = -matcl::constants::f_inf();

            mp_float z(val);
            if (z != val)
                res += 1.0;
            if (z.get_precision() != default_prec)
                res += 1;
        };

        //conversions
        {
            Integer tmp   = 15;
            mp_int val(tmp);
            precision prec = precision(500);

            mp_float z(val, prec);
            if (z != tmp)
                res += 1.0;
            if (z.get_precision() != prec)
                res += 1;
        };
        {
            Integer tmp   = 15;
            mp_int val(tmp);
            precision prec = precision(0);

            mp_float z(val, prec);
            if (z != tmp)
                res += 1.0;
            if (z.get_precision() != default_prec)
                res += 1;
        };

        {
            double tmp   = 1./2.;
            mp_rational val(1,2);
            precision prec = precision(500);

            mp_float z(val, prec);
            if (z != tmp)
                res += 1.0;
            if (z.get_precision() != prec)
                res += 1;
        };
        {
            double tmp   = 1./2.;
            mp_rational val(1,2);
            precision prec = precision(0);

            mp_float z(val, prec);
            if (z != tmp)
                res += 1.0;
            if (z.get_precision() != default_prec)
                res += 1;
        };

        //copy constructors
        {
            mp_float val = -1.0f;

            mp_float z(val);

            if (z != val)
                res += 1.0;
            if (z.get_precision() != default_prec)
                res += 1;
        };
        
        {
            mp_float val = -1.0;
            mp_float val2 = val;

            mp_float z(std::move(val2));
            if (z != val)
                res += 1.0;
            if (z.get_precision() != default_prec)
                res += 1;
        };

        {
            mp_float val = -1.0f;
            precision prec = precision(10);

            mp_float z(val, prec);

            if (z != val)
                res += 1.0;
            if (z.get_precision() != prec)
                res += 1;
        };
        
        {
            mp_float val = -1.0;
            mp_float val2 = val;
            precision prec = precision(10);

            mp_float z(std::move(val2), prec);
            if (z != val)
                res += 1.0;
            if (z.get_precision() != prec)
                res += 1;
        };

        //string constructors
        {
            mp_float z("123.5");
            if (z != 123.5)
                res += 1.0;
            if (z.get_precision() != default_prec)
                res += 1;
        };
        {
            mp_float z("-123.5");   
            if (z != -123.5)
                res += 1.0;
            if (z.get_precision() != default_prec)
                res += 1;
        };
        {
            mp_float z("-1.235e2");   
            if (z != -123.5)
                res += 1.0;
            if (z.get_precision() != default_prec)
                res += 1;
        };

        {
            precision prec = precision(300);
            mp_float z("-123.5", prec);   
            if (z != -123.5)
                res += 1.0;
            if (z.get_precision() != prec)
                res += 1.0;
        };
        {
            precision prec = precision(300);
            mp_float z("-1.235e2", prec);   
            if (z != -123.5)
                res += 1.0;
            if (z.get_precision() != prec)
                res += 1.0;
        };
        {
            mp_float z("-inf");            
            if (z != -constants::inf())
                res += 1.0;
            if (z.get_precision() != default_prec)
                res += 1;
        };
        {
            mp_float z("inf");            
            if (z != constants::inf())
                res += 1.0;
            if (z.get_precision() != default_prec)
                res += 1;
        };
        {
            mp_float z("nan");            
            if (is_nan(z) == false)
                res += 1.0;
            if (z.get_precision() != default_prec)
                res += 1;
        };

        //assignments
        {
            Integer val = 1;

            mp_float z;
            z = val;
            if (z != val)
                res += 1.0;
            if (z.get_precision() != default_prec)
                res += 1;
        };
        {
            float val = 1;

            mp_float z;
            z = val;
            if (z != val)
                res += 1.0;
            if (z.get_precision() != default_prec)
                res += 1;
        };
        {
            double val = 1;

            mp_float z;
            z = val;
            if (z != val)
                res += 1.0;
            if (z.get_precision() != default_prec)
                res += 1;
        };

        {
            precision prec = precision(30);
            mp_float val(-1.0f, prec);

            mp_float z;
            z = val;
            if (z != val)
                res += 1.0;
            if (z.get_precision() != prec)
                res += 1;
        };
        {
            precision prec = precision(30);
            mp_float val(-1.0f, prec);
            mp_float val2 = val;            

            mp_float z;
            z = std::move(val2);
            if (z != val)
                res += 1.0;
            if (z.get_precision() != prec)
                res += 1;
        };

        //casts
        {
            int val = 123;
            mp_float z(val);
            mp_int i = z.cast_mp_int();
            mp_rational q = z.cast_mp_rational();
            mp_complex zz = z.cast_mp_complex();

            if (i != val)
                res += 1;
            if (q != val)
                res += 1;
            if (zz != val)
                res += 1;
        };

        //member func
        {
            mp_float z(0.0);
            if (z)
                res += 1.0;
        };
        {
            mp_float z(-1.01);
            if (!z)
                res += 1.0;
        };

        {
            mp_float z(-123.234);
            precision prec = z.get_precision();
            z.set_precision(precision(2 * prec));

            precision prec2 = z.get_precision();

            if (prec2 != 2 * prec)
                res += 1;
        }
        {
            mp_float z(-123.5);

            if (z.to_string() != "-123.5")
                res += 1.0;
        };

        {
            std::ostringstream ss;
            oarchive ia(ss);

            mp_float z(-1.23);

            for(int i = 0; i < 0; ++i)
                z   = z*z;

            z.serialize(ia.get(), 0);

            std::istringstream ss2(ss.str());
            iarchive ia2(ss2);

            mp_float z2;
            z2.serialize(ia2.get(), 0);

            if (z != z2)
                res += 1;
            if (z.get_precision() != z2.get_precision())
                res += 1;
            if (z.get_precision() != default_prec)
                res += 1;
        };

        {
            std::ostringstream ss;
            oarchive ia(ss);

            precision prec = precision(12);
            mp_float z(-1.23, prec);

            z.serialize(ia.get(), 0);

            std::istringstream ss2(ss.str());
            iarchive ia2(ss2);

            mp_float z2;
            z2.serialize(ia2.get(), 0);

            if (z != z2)
                res += 1;
            if (z.get_precision() != z2.get_precision())
                res += 1;
            if (z.get_precision() != prec)
                res += 1;
        };        
    }
    catch(...)
    {
        res += 1;
    };

    mp_float::set_default_precision(old_prec);

    return res;
};

static double test_prec(const mp_complex& z, precision prec)
{
    if (z.get_precision() != prec)
        return 1.0;
    if (z.real().get_precision() != prec)
        return 1.0;
    if (z.imag().get_precision() != prec)
        return 1.0;

    return 0.0;
};
static double test_prec(const mp_float& z, precision prec)
{
    if (z.get_precision() != prec)
        return 1.0;

    return 0.0;
};

double gmp_tester::test_mp_complex()
{
    double res = 0;

    precision old_prec      = mp_float::get_default_precision();
    precision default_prec  = precision(118);
    mp_float::set_default_precision(default_prec);

    try
    {        
        //constructors and conversions
        {
            mp_complex z;              
            if (z != 0.)
                res += 1.0;
            res += test_prec(z, default_prec);
        };
        {
            Integer val = 1;

            mp_complex z(val);
            if (z != val)
                res += 1.0;
            res += test_prec(z, default_prec);
        };
        {
            double val = 1;

            mp_complex z(val);
            if (z != val)
                res += 1.0;
            res += test_prec(z, default_prec);
        };
        {
            float val = 1;

            mp_complex z(val);
            if (z != val)
                res += 1.0;
            res += test_prec(z, default_prec);
        };
        {
            Integer val = 1;
            precision prec = precision(500);
            mp_complex z(val, prec);

            if (z != val)
                res += 1.0;
            res += test_prec(z, prec);
        };
        {
            Integer val = 1;
            precision prec = precision(0);
            mp_complex z(val, prec);

            if (z != val)
                res += 1.0;
            res += test_prec(z, default_prec);
        };

        {
            double val = 1;
            precision prec = precision(500);
            mp_complex z(val, prec);

            if (z != val)
                res += 1.0;
            res += test_prec(z, prec);
        };
        {
            double val = 1;
            precision prec = precision(0);
            mp_complex z(val, prec);

            if (z != val)
                res += 1.0;
            res += test_prec(z, default_prec);
        };
        {
            float val = 1;
            precision prec = precision(500);
            mp_complex z(val, prec);

            if (z != val)
                res += 1.0;
            res += test_prec(z, prec);
        };
        {
            float val = 1;
            precision prec = precision(0);
            mp_complex z(val, prec);

            if (z != val)
                res += 1.0;
            res += test_prec(z, default_prec);
        };
        {
            double val = matcl::constants::nan();

            mp_complex z(val);
            if (is_nan(z.real()) == false)
                res += 1.0;
            if (is_zero(z.imag()) == false)
                res += 1.0;
            res += test_prec(z, default_prec);
        };
        {
            double val = matcl::constants::inf();

            mp_complex z(val);
            if (is_inf(z.real()) == false)
                res += 1.0;
            if (is_zero(z.imag()) == false)
                res += 1.0;
            res += test_prec(z, default_prec);
        };

        {
            double val  = matcl::constants::nan();
            precision prec = precision(40);
            mp_complex z(val,prec);
            if (is_nan(z.real()) == false)
                res += 1.0;
            if (is_zero(z.imag()) == false)
                res += 1.0;
            res += test_prec(z, prec);
        };
        {
            double val = matcl::constants::inf();
            precision prec = precision(40);
            mp_complex z(val,prec);
            if (is_inf(z.real()) == false)
                res += 1.0;
            if (is_zero(z.imag()) == false)
                res += 1.0;
            res += test_prec(z, prec);
        };

        {
            Complex val = 1;
            precision prec = precision(500);
            mp_complex z(val, prec);

            if (z != val)
                res += 1.0;
            res += test_prec(z, prec);
        };
        {
            Complex val = 1;
            precision prec = precision(0);
            mp_complex z(val, prec);

            if (z != val)
                res += 1.0;
            res += test_prec(z, default_prec);
        };
        {
            Complex val(1,2);
            mp_complex z(val);

            if (z != val)
                res += 1.0;
            res += test_prec(z, default_prec);
        };

        {
            Float_complex val(1, constants::f_inf());
            precision prec = precision(0);
            mp_complex z(val, prec);

            if (z != val)
                res += 1.0;
            res += test_prec(z, default_prec);
        };
        {
            Float_complex val(-1,2);
            mp_complex z(val);

            if (z != val)
                res += 1.0;
            res += test_prec(z, default_prec);
        };

        {
            mp_float val(1,precision(20));
            precision prec = precision(500);
            mp_complex z(val, prec);

            if (z != val)
                res += 1.0;
            res += test_prec(z, prec);
        };
        {
            mp_float val(1,precision(20));
            precision prec = precision(0);
            mp_complex z(val, prec);

            if (z != val)
                res += 1.0;
            res += test_prec(z, default_prec);
        };
        {
            mp_float val(1,precision(20));
            mp_complex z(val);

            if (z != val)
                res += 1.0;
            res += test_prec(z, precision(20));
        };

        {
            mp_float val(1,precision(20));
            mp_float val2= val;
            precision prec = precision(500);
            mp_complex z(std::move(val2), prec);

            if (z != val)
                res += 1.0;
            res += test_prec(z, prec);
        };
        {
            mp_float val(1,precision(20));
            mp_float val2= val;
            precision prec = precision(0);
            mp_complex z(std::move(val2), prec);

            if (z != val)
                res += 1.0;
            res += test_prec(z, default_prec);
        };
        {
            mp_float val(1,precision(20));
            mp_float val2= val;
            mp_complex z(std::move(val2));

            if (z != val)
                res += 1.0;
            res += test_prec(z, precision(20));
        };

        {
            double re = 1;
            double im = 2;
            precision prec = precision(500);
            mp_complex z(re, im, prec);

            if (z.real() != re)
                res += 1.0;
            if (z.imag() != im)
                res += 1.0;
            res += test_prec(z, prec);
        };
        {
            double re = 1;
            double im = 2;
            precision prec = precision(0);
            mp_complex z(re, im, prec);

            if (z.real() != re)
                res += 1.0;
            if (z.imag() != im)
                res += 1.0;
            res += test_prec(z, default_prec);
        };        
        {
            double re = 1;
            double im = 2;
            precision prec = precision(0);
            mp_complex z(re, im);

            if (z.real() != re)
                res += 1.0;
            if (z.imag() != im)
                res += 1.0;
            res += test_prec(z, default_prec);
        };        

        {
            float re = 1;
            float im = 2;
            precision prec = precision(500);
            mp_complex z(re, im, prec);

            if (z.real() != re)
                res += 1.0;
            if (z.imag() != im)
                res += 1.0;
            res += test_prec(z, prec);
        };
        {
            float re = 1;
            float im = 2;
            precision prec = precision(0);
            mp_complex z(re, im, prec);

            if (z.real() != re)
                res += 1.0;
            if (z.imag() != im)
                res += 1.0;
            res += test_prec(z, default_prec);
        };        
        {
            float re = 1;
            float im = 2;
            precision prec = precision(0);
            mp_complex z(re, im);

            if (z.real() != re)
                res += 1.0;
            if (z.imag() != im)
                res += 1.0;
            res += test_prec(z, default_prec);
        };        

        {
            mp_float re(1, precision(20));
            mp_float im(2, precision(40));

            mp_float re2 = re;
            mp_float im2 = im;

            precision prec = precision(500);
            mp_complex z(re, im, prec);

            if (z.real() != re)
                res += 1.0;
            if (z.imag() != im)
                res += 1.0;
            res += test_prec(z, prec);
        };
        {
            mp_float re(1, precision(20));
            mp_float im(2, precision(40));

            mp_float re2 = re;
            mp_float im2 = im;

            precision prec = precision(0);
            mp_complex z(re, im, prec);

            if (z.real() != re)
                res += 1.0;
            if (z.imag() != im)
                res += 1.0;
            res += test_prec(z, default_prec);
        };        
        {
            mp_float re(1, precision(20));
            mp_float im(2, precision(40));

            mp_float re2 = re;
            mp_float im2 = im;

            precision prec = precision(0);
            mp_complex z(re, im);

            if (z.real() != re)
                res += 1.0;
            if (z.imag() != im)
                res += 1.0;
            res += test_prec(z, precision(40));
        };        

        {
            mp_float re(1, precision(20));
            mp_float im(2, precision(40));

            mp_float re2 = re;
            mp_float im2 = im;

            precision prec = precision(500);
            mp_complex z(std::move(re2), im, prec);

            if (z.real() != re)
                res += 1.0;
            if (z.imag() != im)
                res += 1.0;
            res += test_prec(z, prec);
        };
        {
            mp_float re(1, precision(20));
            mp_float im(2, precision(40));

            mp_float re2 = re;
            mp_float im2 = im;

            precision prec = precision(0);
            mp_complex z(std::move(re2), im, prec);

            if (z.real() != re)
                res += 1.0;
            if (z.imag() != im)
                res += 1.0;
            res += test_prec(z, default_prec);
        };        
        {
            mp_float re(1, precision(20));
            mp_float im(2, precision(40));

            mp_float re2 = re;
            mp_float im2 = im;

            precision prec = precision(0);
            mp_complex z(std::move(re2), im);

            if (z.real() != re)
                res += 1.0;
            if (z.imag() != im)
                res += 1.0;
            res += test_prec(z, precision(40));
        };  

        {
            mp_float re(1, precision(20));
            mp_float im(2, precision(40));

            mp_float re2 = re;
            mp_float im2 = im;

            precision prec = precision(500);
            mp_complex z(re, std::move(im2), prec);

            if (z.real() != re)
                res += 1.0;
            if (z.imag() != im)
                res += 1.0;
            res += test_prec(z, prec);
        };
        {
            mp_float re(1, precision(20));
            mp_float im(2, precision(40));

            mp_float re2 = re;
            mp_float im2 = im;

            precision prec = precision(0);
            mp_complex z(re, std::move(im2), prec);

            if (z.real() != re)
                res += 1.0;
            if (z.imag() != im)
                res += 1.0;
            res += test_prec(z, default_prec);
        };        
        {
            mp_float re(1, precision(20));
            mp_float im(2, precision(40));

            mp_float re2 = re;
            mp_float im2 = im;

            precision prec = precision(0);
            mp_complex z(re, std::move(im2));

            if (z.real() != re)
                res += 1.0;
            if (z.imag() != im)
                res += 1.0;
            res += test_prec(z, precision(40));
        };  

        {
            mp_float re(1, precision(20));
            mp_float im(2, precision(40));

            mp_float re2 = re;
            mp_float im2 = im;

            precision prec = precision(500);
            mp_complex z(std::move(re2), std::move(im2), prec);

            if (z.real() != re)
                res += 1.0;
            if (z.imag() != im)
                res += 1.0;
            res += test_prec(z, prec);
        };
        {
            mp_float re(1, precision(20));
            mp_float im(2, precision(40));

            mp_float re2 = re;
            mp_float im2 = im;

            precision prec = precision(0);
            mp_complex z(std::move(re2), std::move(im2), prec);

            if (z.real() != re)
                res += 1.0;
            if (z.imag() != im)
                res += 1.0;
            res += test_prec(z, default_prec);
        };        
        {
            mp_float re(1, precision(20));
            mp_float im(2, precision(40));

            mp_float re2 = re;
            mp_float im2 = im;

            precision prec = precision(0);
            mp_complex z(std::move(re2), std::move(im2));

            if (z.real() != re)
                res += 1.0;
            if (z.imag() != im)
                res += 1.0;
            res += test_prec(z, precision(40));
        };  

        //conversions
        {
            mp_int re(112);
            precision prec = precision(500);
            mp_complex z(re, prec);

            if (z.real() != mp_float(re))
                res += 1.0;
            res += test_prec(z, prec);
        };
        {
            mp_int re(112);
            precision prec = precision(0);
            mp_complex z(re, prec);

            if (z.real() != mp_float(re))
                res += 1.0;
            res += test_prec(z, default_prec);
        };        
        {
            mp_int re(112);
            precision prec = precision(0);
            mp_complex z(re);

            if (z.real() != mp_float(re))
                res += 1.0;
            res += test_prec(z, default_prec);
        };       

        {
            mp_rational re(112,113);
            precision prec = precision(500);
            mp_complex z(re, prec);

            if (z.real() != mp_float(re, prec))
                res += 1.0;
            res += test_prec(z, prec);
        };
        {
            mp_rational re(112,113);
            precision prec = precision(0);
            mp_complex z(re, prec);

            if (z.real() != mp_float(re,prec))
                res += 1.0;
            res += test_prec(z, default_prec);
        };        
        {
            mp_rational re(112,113);
            precision prec = precision(0);
            mp_complex z(re);

            if (z.real() != mp_float(re))
                res += 1.0;
            res += test_prec(z, default_prec);
        };

        //copy constructors
        {
            mp_complex val(-1.0f, precision(20));

            mp_complex z(val);

            if (z != val)
                res += 1.0;
            res += test_prec(z, precision(20));
        };
        {
            mp_complex val(-1.0f, precision(20));
            precision prec = precision(60);
            mp_complex z(val, prec);

            if (z != val)
                res += 1.0;
            res += test_prec(z, prec);
        };
        {
            mp_complex val(-1.0f, precision(20));
            precision prec = precision(0);
            mp_complex z(val, prec);

            if (z != val)
                res += 1.0;
            res += test_prec(z, default_prec);
        };

        {
            mp_complex val(-1.0f, precision(20));
            mp_complex val2(val);

            mp_complex z(std::move(val2));

            if (z != val)
                res += 1.0;
            res += test_prec(z, precision(20));
        };
        {
            mp_complex val(-1.0f, precision(20));
            mp_complex val2(val);
            precision prec = precision(60);
            mp_complex z(std::move(val2), prec);

            if (z != val)
                res += 1.0;
            res += test_prec(z, prec);
        };
        {
            mp_complex val(-1.0f, precision(20));
            mp_complex val2(val);
            precision prec = precision(0);
            mp_complex z(std::move(val2), prec);

            if (z != val)
                res += 1.0;
            res += test_prec(z, default_prec);
        };

        //assignments
        {
            Integer val = 1;

            mp_complex z(0, precision(20));
            z = val;
            if (z != val)
                res += 1.0;
            res += test_prec(z, default_prec);
        };
        {
            float val = 1;

            mp_complex z(0, precision(20));
            z = val;
            if (z != val)
                res += 1.0;
            res += test_prec(z, default_prec);
        };
        {
            double val = 1;

            mp_complex z(0, precision(20));
            z = val;
            if (z != val)
                res += 1.0;
            res += test_prec(z, default_prec);
        };
        {
            Complex val = 1;

            mp_complex z(0, precision(20));
            z = val;
            if (z != val)
                res += 1.0;
            res += test_prec(z, default_prec);
        };
        {
            Float_complex val = 1;

            mp_complex z(0, precision(20));
            z = val;
            if (z != val)
                res += 1.0;
            res += test_prec(z, default_prec);
        };
        {
            precision prec = precision(30);
            mp_float val(-1.0f, prec);

            mp_complex z(0, precision(20));
            z = val;
            if (z != val)
                res += 1.0;
            res += test_prec(z, prec);
        };
        {
            precision prec = precision(30);
            mp_float val(-1.0f, prec);
            mp_float val2 = val;            

            mp_complex z(0, precision(20));
            z = std::move(val2);
            if (z != val)
                res += 1.0;
            res += test_prec(z, prec);
        };

        {
            precision prec = precision(30);
            mp_complex val(-1.0f, prec);

            mp_complex z(0, precision(20));
            z = val;
            if (z != val)
                res += 1.0;
            res += test_prec(z, prec);
        };
        {
            precision prec = precision(30);
            mp_complex val(-1.0f, prec);
            mp_complex val2 = val;            

            mp_complex z(0, precision(20));
            z = std::move(val2);
            if (z != val)
                res += 1.0;
            res += test_prec(z, prec);
        };

        //cast
        {
            int val = 123;
            mp_complex z(val);
            mp_float f = z.cast_mp_float();
            mp_rational q = z.cast_mp_rational();
            mp_complex i = z.cast_mp_int();

            if (f != val)
                res += 1;
            if (q != val)
                res += 1;
            if (i != val)
                res += 1;
        };

        //member func
        {
            precision prec = precision(30);
            float re = -1.0f;
            float im = -2.9f;
            mp_complex val(re, im, prec);     

            mp_float vre = val.real();
            mp_float vim = val.imag();

            if (vre != re)
                res += 1.0;
            if (vim != im)
                res += 1.0;
            res += test_prec(vre, prec);
            res += test_prec(vim, prec);
        };
        {
            mp_complex z(1.0, 2.0, precision(40));
            precision prec(10);
            mp_float re(4.0, prec);

            z.set_real(re);

            if (z.real() != re)
                res += 1.0;
            res += test_prec(z, prec);
        }
        {
            mp_complex z(1.0, 2.0, precision(40));
            precision prec(100);
            mp_float re(4.0, prec);
            mp_float re2(re);

            z.set_real(std::move(re2));

            if (z.real() != re)
                res += 1.0;
            res += test_prec(z, prec);
        }
        {
            mp_complex z(1.0, 2.0, precision(40));
            precision prec(10);
            mp_float im(4.0, prec);

            z.set_imag(im);

            if (z.imag() != im)
                res += 1.0;
            res += test_prec(z, prec);
        }
        {
            mp_complex z(1.0, 2.0, precision(40));
            precision prec(100);
            mp_float im(4.0, prec);
            mp_float im2(im);

            z.set_imag(std::move(im2));

            if (z.imag() != im)
                res += 1.0;
            res += test_prec(z, prec);
        }

        {
            mp_complex z(0.0);
            if (z)
                res += 1.0;
        };
        {
            mp_complex z(-1.01);
            if (!z)
                res += 1.0;
        };

        {
            mp_complex z(-123.234, 2.991);
            precision prec = z.get_precision();
            z.set_precision(precision(2 * prec));

            precision prec2 = z.get_precision();

            if (prec2 != 2 * prec)
                res += 1;
        }
        {
            mp_complex z(-123.5, 123.5);

            if (z.to_string() != "-123.5+123.5i")
                res += 1.0;
        };

        {
            std::ostringstream ss;
            oarchive ia(ss);

            mp_complex z(-1.23, 1.23);

            z.serialize(ia.get(), 0);

            std::istringstream ss2(ss.str());
            iarchive ia2(ss2);

            mp_complex z2;
            z2.serialize(ia2.get(), 0);

            if (z != z2)
                res += 1;
            if (z.get_precision() != z2.get_precision())
                res += 1;

            res += test_prec(z2, default_prec);
        };

        {
            std::ostringstream ss;
            oarchive ia(ss);

            mp_complex z(-1.23, constants::inf());

            z.serialize(ia.get(), 0);

            std::istringstream ss2(ss.str());
            iarchive ia2(ss2);

            mp_complex z2;
            z2.serialize(ia2.get(), 0);

            if (z != z2)
                res += 1;
            if (z.get_precision() != z2.get_precision())
                res += 1;

            res += test_prec(z2, default_prec);
        };
        {
            std::ostringstream ss;
            oarchive ia(ss);

            mp_complex z(-1.23, constants::nan());

            z.serialize(ia.get(), 0);

            std::istringstream ss2(ss.str());
            iarchive ia2(ss2);

            mp_complex z2;
            z2.serialize(ia2.get(), 0);

            if (z.real() != z2.real())
                res += 1;
            if (is_nan(z2.imag()) == false)
                res += 1;
            if (z.get_precision() != z2.get_precision())
                res += 1;

            res += test_prec(z2, default_prec);
        };
        {
            std::ostringstream ss;
            oarchive ia(ss);

            precision prec = precision(12);
            mp_complex z(-1.23, 1.23, prec);

            z.serialize(ia.get(), 0);

            std::istringstream ss2(ss.str());
            iarchive ia2(ss2);

            mp_complex z2;
            z2.serialize(ia2.get(), 0);

            if (z != z2)
                res += 1;
            if (z.get_precision() != z2.get_precision())
                res += 1;
            res += test_prec(z2, prec);
        };
        {
            std::ostringstream ss;
            oarchive ia(ss);

            precision prec = precision(12);
            mp_complex z(constants::inf(), 1.23, prec);

            z.serialize(ia.get(), 0);

            std::istringstream ss2(ss.str());
            iarchive ia2(ss2);

            mp_complex z2;
            z2.serialize(ia2.get(), 0);

            if (z != z2)
                res += 1;
            if (z.get_precision() != z2.get_precision())
                res += 1;
            res += test_prec(z2, prec);
        };
        {
            std::ostringstream ss;
            oarchive ia(ss);

            precision prec = precision(12);
            mp_complex z(constants::nan(), 1.23, prec);

            z.serialize(ia.get(), 0);

            std::istringstream ss2(ss.str());
            iarchive ia2(ss2);

            mp_complex z2;
            z2.serialize(ia2.get(), 0);

            if (z.imag() != z2.imag())
                res += 1;
            if (is_nan(z2.real()) == false)
                res += 1;
            if (z.get_precision() != z2.get_precision())
                res += 1;
            res += test_prec(z2, prec);
        };
    }
    catch(...)
    {
        res += 1;
    };

    mp_float::set_default_precision(old_prec);
    return res;
};

double gmp_tester::test_mp_rational()
{
    double res = 0;

    try
    {
        {
            mp_rational z;
            if (z != 0)
                res += 1.0;
        };
        {
            int val = 1;

            mp_rational z(val);
            if (z != val)
                res += 1.0;
        };
        {
            mp_int val = 1;

            mp_rational z(val);
            if (z != val)
                res += 1.0;
        };

        //conversions
        {
            Float val = 0.25;
            mp_rational z(val);

            if (z*4 != 1)
                res += 1;
        };
        {
            Real val = 0.25;
            mp_rational z(val);

            if (z*4 != 1)
                res += 1;
        };
        {
            mp_float val(0.25);
            mp_rational z(val);

            if (z*4 != 1)
                res += 1;
        };

        {
            Integer num = 2;
            Integer den = 5;

            mp_rational z(num, den);
            if (z * den != num)
                res += 1.0;
        };
        {
            Integer num = 2;
            Integer den = -5;

            mp_rational z(num, den);
            if (z * den != num)
                res += 1.0;
        };
        {
            Integer num = -2;
            Integer den = 5;

            mp_rational z(num, den);
            if (z * den != num)
                res += 1.0;
        };
        {
            Integer num = -2;
            Integer den = -5;

            mp_rational z(num, den);
            if (z * den != num)
                res += 1.0;
        };
        {
            Integer num = 2;
            Integer den = 4;

            mp_rational z(num, den);
            if (z * den != num)
                res += 1.0;
        };
        {
            Integer num = 4;
            Integer den = 2;

            mp_rational z(num, den);
            if (z * den != num)
                res += 1.0;
        };

        {
            mp_int num = 2;
            mp_int den = 5;

            mp_rational z(num, den);
            if (z * den != num)
                res += 1.0;
        };
        {
            mp_int num = 2;
            mp_int den = 4;

            mp_rational z(num, den);
            if (z * den != num)
                res += 1.0;
        };
        {
            mp_int num = -2;
            mp_int den = 4;

            mp_rational z(num, den);
            if (z * den != num)
                res += 1.0;
        };
        {
            mp_int num = 2;
            mp_int den = -4;

            mp_rational z(num, den);
            if (z * den != num)
                res += 1.0;
        };
        {
            mp_int num = -2;
            mp_int den = -4;

            mp_rational z(num, den);
            if (z * den != num)
                res += 1.0;
        };

        {
            mp_int num = 4;
            mp_int den = 2;

            mp_rational z(num, den);
            if (z * den != num)
                res += 1.0;
        };

        {
            mp_rational val(16,331);
            mp_rational val2(val);

            if (val != val2)
                res += 1.0;
        };
        {
            mp_rational val(16,331);
            mp_rational val2(val);
            mp_rational val3(std::move(val2));

            if (val != val3)
                res += 1.0;
        };

        {
            mp_rational z("123/17");
            if (z != mp_rational(123,17))
                res += 1.0;
        };
        {
            mp_rational z("-123/17");
            if (z != mp_rational(123,-17))
                res += 1.0;
        };
        {
            mp_rational z("2/4");
            if (z != mp_rational(1,2))
                res += 1.0;
        };
        {
            mp_rational z("-2/4");
            if (z != mp_rational(-1,2))
                res += 1.0;
        };
        {
            mp_rational z("-2/-4");
            if (z != mp_rational(1,2))
                res += 1.0;
        };
        {
            mp_rational z("2/-4");
            if (z != mp_rational(-1,2))
                res += 1.0;
        };

        {
            Integer val = 1;

            mp_rational z;
            z = val;
            if (z != val)
                res += 1.0;
        };
        {
            mp_int val = 1;

            mp_rational z;
            z = val;
            if (z != val)
                res += 1.0;
        };
        {
            mp_int val = 1;
            mp_int val2 = val;

            mp_rational z;
            z = std::move(val2);
            if (z != val)
                res += 1.0;
        };
        {
            mp_rational val(1, -1);

            mp_rational z;
            z = val;
            if (z != val)
                res += 1.0;
        };
        {
            mp_rational val(-1, -1);
            mp_rational val2(val);

            mp_rational z;
            z = std::move(val2);

            if (z != val)
                res += 1.0;
        };

        {
            mp_rational z(0);
            if (z)
                res += 1.0;
        };
        {
            mp_rational z(-1, 1);
            if (!z)
                res += 1.0;
        };

        try
        {
            mp_rational z(-1, 0);            
        }
        catch(...)
        {
            res += 1.0;
        };

        try
        {
            mp_rational z = mp_rational(1, 0);            
        }
        catch(...)
        {
            res += 1.0;
        };

        try
        {
            mp_rational z = mp_rational(0, 0);            
        }
        catch(...)
        {
            res += 1.0;
        };

        try
        {
            mp_rational z = mp_rational(1, -0);            
        }
        catch(...)
        {
            res += 1.0;
        };

        try
        {
            mp_rational z = mp_rational(1, +0);            
        }
        catch(...)
        {
            res += 1.0;
        };

        try
        {
            mp_rational z("-2/-.4");
            res += 1;
        }
        catch(...)
        {};

        {
            mp_int num(5);
            mp_int den(7);
            mp_rational val(num, den);

            if (val.numerator() != num)
                res += 1.0;
            if (val.denominator() != den)
                res += 1.0;
        };
        {
            mp_int num(5);
            mp_int den(7);

            mp_rational val(num * 2, den * 2);

            if (val.numerator() != num)
                res += 1.0;
            if (val.denominator() != den)
                res += 1.0;
        };

        //casts
        {
            int val = 123;
            mp_rational z(val);
            mp_float f = z.cast_mp_float();
            mp_rational i = z.cast_mp_int();
            mp_complex q = z.cast_mp_complex();

            if (f != val)
                res += 1;
            if (i != val)
                res += 1;
            if (z != val)
                res += 1;
        };
        {
            std::ostringstream ss;
            oarchive ia(ss);

            mp_rational z(123,321);
            z.serialize(ia.get(), 0);

            std::istringstream ss2(ss.str());
            iarchive ia2(ss2);

            mp_rational z2;
            z2.serialize(ia2.get(), 0);

            if (z != z2)
                res += 1;
        };

        {
            std::ostringstream ss;
            oarchive ia(ss);

            mp_rational z(123,-321);
            z.serialize(ia.get(), 0);

            std::istringstream ss2(ss.str());
            iarchive ia2(ss2);

            mp_rational z2;
            z2.serialize(ia2.get(), 0);

            if (z != z2)
                res += 1;
        };

        {
            mp_rational z(-2,4);
            if (z.to_string() != "-1/2")
                res += 1.0;
        };

        return res;
    }
    catch(std::exception& ex)
    {
        std::string err = ex.what();
        std::cout << err << "\n";
        return 1.0;
    }
    catch(...)
    {
        return 1.0;
    };
};

bool gmp_tester::test_reim_impl()
{
    std::vector<Scalar> scalars1;

    Integer N   = 1000;
    rand_scalars::make(scalars1, N, true);

    double res = 0.0;

    for (Integer i = 0; i < 1000; ++i)
        res     += test_reim(scalars1[i], i);

    return res == 0.0 ? true : false;
};

bool gmp_tester::test_uminus_impl()
{
    std::vector<Scalar> scalars1;

    Integer N   = 1000;
    rand_scalars::make(scalars1, N, true);

    double res = 0.0;

    for (Integer i = 0; i < 1000; ++i)
        res     += test_uminus(scalars1[i], i);

    return res == 0.0 ? true : false;
};

bool gmp_tester::test_is_impl()
{
    std::vector<Scalar> scalars1;

    Integer N   = 1000;
    rand_scalars::make(scalars1, N, true);

    double res = 0.0;

    for (Integer i = 0; i < 1000; ++i)
        res     += test_is(scalars1[i], i);

    return res == 0.0 ? true : false;
};
bool gmp_tester::test_next_impl()
{
    std::vector<Scalar> scalars1;

    Integer N   = 1000;
    rand_scalars::make(scalars1, N, true);

    double res = 0.0;

    for (Integer i = 0; i < 1000; ++i)
        res     += test_next(scalars1[i], i);

    return res == 0.0 ? true : false;
};
bool gmp_tester::test_sign_impl()
{
    std::vector<Scalar> scalars1;

    Integer N   = 1000;
    rand_scalars::make(scalars1, N, true);

    double res = 0.0;

    for (Integer i = 0; i < 1000; ++i)
        res     += test_sign(scalars1[i], i);

    return res == 0.0 ? true : false;
};
bool gmp_tester::test_eps_impl()
{
    std::vector<Scalar> scalars1;

    Integer N   = 1000;
    rand_scalars::make(scalars1, N, true);

    double res = 0.0;

    for (Integer i = 0; i < 1000; ++i)
        res     += test_eps(scalars1[i], i);

    return res == 0.0 ? true : false;
};

template<class Func>
void gmp_tester::test_scalar_func()
{
    std::vector<Scalar> scalars1;

    Integer N   = 1000;
    rand_scalars::make(scalars1, N, true);

    double res = 0.0;

    for (Integer i = 0; i < 1000; ++i)
        res     += test_scalar<Func>(scalars1[i], i);

    if (res != 0.0)
        std::cout << Func::name() << ": FAILED" << "\n";
    else
        std::cout << Func::name() << ": ok" << "\n";        
};

template<class Derived>
struct eval_function : eval_scalars_1<eval_function<Derived>>
{
    Integer code;
    bool    test_fin;
    bool    zeroinf_allowed;

    eval_function(Integer c, bool test_fin_, bool zeroinf_allowed_) 
        : code(c), test_fin(test_fin_), zeroinf_allowed(zeroinf_allowed_)
    {};

    template<class T1>
    auto eval_op(const T1& a1) -> decltype(Derived::eval(std::declval<T1>()))
    {
        return Derived::eval(a1);
    };

    template<class T1, class T2>
    bool different(const T1& a, const T2& b)
    {
        if (a == b)
            return false;

        if (is_nan(a) == true && is_nan(b) == true)
            return false;

        mp_complex v1(a);
        mp_complex v2(b, get_prec<T2>());

        mp_float dif    = abs(v1 - v2);
        auto tol        = 2*eps(v2);

        if (dif < tol)
            return false;

        return true;
    }

    template<class T1>
    double eval_scal_func(const T1& s1)
    {
        if (code == -1)
            std::cout << "break\n";

        auto res        = eval_op(s1);

        Integer     v0_1    = convert_scalar<Integer>(s1);
        Real        v0_2    = convert_scalar<Real>(s1);
        Complex     v0_3    = convert_scalar<Complex>(s1);

        bool zeroinf_1      = is_zero(v0_1) || is_finite(v0_1) == false;
        bool zeroinf_2      = is_zero(v0_2) || is_finite(v0_2) == false;
        bool zeroinf_3      = is_zero(real(v0_3)) || is_finite(real(v0_3)) == false;

        bool is_nan         = (test_fin == true)  && matcl::is_nan(s1);
        bool is_fin         = (test_fin == false) || matcl::is_finite(s1);

        mp_int      v1_1(v0_1);
        mp_float    v1_2(v0_2);
        mp_rational v1_3(v0_1);
        mp_complex  v1_4(v0_3);

        double out = 0;

        auto res1_1     = eval_op(v1_1);
        auto res1_2     = eval_op(v0_1);
        auto res2_1     = eval_op(v1_2) ;
        auto res2_2     = eval_op(v0_2);
        auto res3_1     = eval_op(v1_3);
        auto res3_2     = eval_op(v0_1);
        auto res4_1     = eval_op(v1_4);
        auto res4_2     = eval_op(v0_3);

        bool test_1     = zeroinf_allowed == true || zeroinf_1 == false;
        bool test_2     = zeroinf_allowed == true || zeroinf_2 == false;
        bool test_3     = zeroinf_allowed == true || zeroinf_1 == false;
        bool test_4     = zeroinf_allowed == true || zeroinf_3 == false;

        if (different(res1_1, res1_2) && is_fin == true && test_1)
            out     += 1;
        if (different(res2_1, res2_2) && is_nan == false  && test_2)
            out     += 1;
        if (different(res3_1, res3_2) && is_fin == true  && test_3)
            out     += 1;
        if (different(res4_1, res4_2) && is_nan == false  && test_4)
            out     += 1;

        if (out != 0)
        {
            std::cout << code << " " << s1 << "\n";
            std::cout << res1_1 << " " << res1_2 << "\n";
            std::cout << res2_1 << " " << res2_2 << "\n";
            std::cout << res3_1 << " " << res3_2 << "\n";
            std::cout << res4_1 << " " << res4_2 << "\n";
        }

        return out;
    };
};

template<class Derived>
struct eval_function_num : eval_scalars_1<eval_function_num<Derived>>
{
    Integer code;
    Real    mult_compl;
    Real    mult_real;
    bool    naninf_allowed;

    eval_function_num(Integer c, Real mult_compl_, Real mult_real_, bool naninf_allowed_) 
        : code(c), mult_compl(mult_compl_), mult_real(mult_real_), naninf_allowed(naninf_allowed_)
    {};

    template<class T1>
    auto eval_op(const T1& a1) -> decltype(Derived::eval(std::declval<T1>()))
    {
        return Derived::eval(a1);
    };

    template<class T1, class T2>
    bool different_re(const T1& a, const T2& b)
    {
        if (std::is_same<T1,mp_rational>::value == true)
            return different_rational(a,b);

        if (convert_scalar<Real>(a) == convert_scalar<Real>(b))
            return false;

        if (is_nan(a) == true && is_nan(b) == true)
            return false;

        Real v1         = convert_scalar<Real>(a);
        Real v2         = convert_scalar<Real>(b);

        Real dif        = (Real)abs(v1 - v2);
        Real mult       = 4 * mult_real;
        Real tol        = mult * eps(b);

        if (dif < tol)
            return false;
        return true;
    }
    template<class T1, class T2>
    bool different_rational(const T1& a, const T2& b)
    {
        if (convert_scalar<Real>(a) == convert_scalar<Real>(b))
            return false;

        if (is_nan(a) == true && is_nan(b) == true)
            return false;

        if (is_finite(b) == false)
            return false;

        Real v1         = convert_scalar<Real>(a);
        Real v2         = convert_scalar<Real>(b);

        Real dif        = (Real)abs(v1 - v2);
        Real mult       = 4*mult_real;
        Real tol        = mult * eps(b);

        if (dif < tol)
            return false;
        return true;
    }
    template<class T1, class T2>
    bool different_compl(const T1& a, const T2& b)
    {
        if (convert_scalar<Complex>(a) == convert_scalar<Complex>(b))
            return false;

        if (is_nan(b) == true && is_nan(a) == true)
            return false;

        if (is_finite(a) == false && is_finite(b) == false)
        {
            if (real(a) == real(b))
                return different_re(imag(a), imag(b));

            if (imag(a) == imag(a))
                return different_re(real(a), real(b));
        };

        Complex v1  = convert_scalar<Complex>(a);
        Complex v2  = convert_scalar<Complex>(b);

        Real dif        = (Real)abs(v1 - v2);
        Real mult       = 10*mult_compl;
        Real tol        = mult*eps(b);

        if (dif < tol)
            return false;

        if (code == -1)
        {
            std::cout << a << " " << b << "\n";
            std::cout << dif << " " << tol << "\n";
        }
        return true;
    }
    template<class T1>
    double eval_scal_func(const T1& s1)
    {
        if (code == -1)
            std::cout << "break\n";

        if (naninf_allowed == false && (is_inf(s1) || is_nan(s1)))
            return 0.0;

        auto res        = eval_op(s1);

        Integer     si  = convert_scalar<Integer>(s1);
        Real        sr  = convert_scalar<Real>(s1);
        Complex     sc  = convert_scalar<Complex>(s1);

        bool test_fin   = false;
        bool is_fin     = matcl::is_finite(s1);

        mp_int      vi(si);
        mp_float    vf(sr);
        mp_rational vr(si);
        mp_complex  vc(sc);

        double out = 0;

        auto res_i_1    = eval_op(vi);
        auto res_i_2    = eval_op(si);
        auto res_f_1    = eval_op(vf) ;
        auto res_f_2    = eval_op(sr);
        auto res_r_1    = eval_op(vr);
        auto res_r_2    = eval_op(si);
        auto res_c_1    = eval_op(vc);
        auto res_c_2    = eval_op(sc);

        bool test_1     = true;
        bool test_2     = true;
        bool test_3     = true;
        bool test_4     = true;

        if (different_re(res_i_1, res_i_2) && is_fin)
            out     += 1;
        if (different_re(res_f_1, res_f_2))
            out     += 1;
        if (different_re(res_r_1, res_r_2) && is_fin)
            out     += 1;
        if (different_compl(res_c_1, res_c_2))
            out     += 1;

        if (out != 0)
        {
            std::cout << code << " " << s1 << "\n";
            std::cout << res_i_1 << " " << res_i_2 << " " << res_i_1 - res_i_2 << "\n";
            std::cout << res_f_1 << " " << res_f_2 << " " << res_f_1 - res_f_2 << "\n";
            std::cout << res_r_1 << " " << res_r_2 << " " << res_r_1 - res_r_2 << "\n";
            std::cout << res_c_1 << " " << res_c_2 << " " << res_c_1 - res_c_2 << "\n";
        }

        return out;
    };
};

struct eval_uminus : eval_function<eval_uminus>
{
    eval_uminus(Integer c) : eval_function(c,true,true){};

    template<class T1>
    static T1 eval(const T1& a1)
    {
        return -a1;
    };
};

struct eval_re : eval_function<eval_re>
{
    eval_re(Integer c) : eval_function(c,true,true){};

    template<class T1>
    static auto eval(const T1& a1) -> decltype(real(std::declval<T1>()))
    {
        return real(a1);
    };
};

struct eval_im : eval_function<eval_im>
{
    eval_im(Integer c) : eval_function(c,true,true){};

    template<class T1>
    static auto eval(const T1& a1) -> decltype(imag(std::declval<T1>()))
    {
        return imag(a1);
    };
};

struct eval_conj : eval_function<eval_conj>
{
    eval_conj(Integer c) : eval_function(c,true,true){};

    template<class T1>
    static auto eval(const T1& a1) -> decltype(conj(std::declval<T1>()))
    {
        return conj(a1);
    };
};

struct eval_abs : eval_function<eval_abs>
{
    eval_abs(Integer c) : eval_function(c,true,true){};

    template<class T1>
    static auto eval(const T1& a1) -> decltype(abs(std::declval<T1>()))
    {
        return abs(a1);
    };
};

struct eval_abs2 : eval_function<eval_abs>
{
    eval_abs2(Integer c) : eval_function(c,true,true){};

    template<class T1>
    static auto eval(const T1& a1) -> decltype(abs2(std::declval<T1>()))
    {
        return abs2(a1);
    };
};

struct eval_arg : eval_function<eval_arg>
{
    eval_arg(Integer c) : eval_function(c,true,true){};

    template<class T1>
    static auto eval(const T1& a1) -> decltype(arg(std::declval<T1>()))
    {
        return arg(a1);
    };
};

struct eval_angle : eval_function<eval_angle>
{
    eval_angle(Integer c) : eval_function(c,true,true){};

    template<class T1>
    static auto eval(const T1& a1) -> decltype(angle(std::declval<T1>()))
    {
        return angle(a1);
    };
};

double gmp_tester::test_uminus(const Scalar& s1, Integer code)
{
    eval_uminus test(code);
    return test.make(s1);
};

double gmp_tester::test_reim(const Scalar& s1, Integer code)
{
    eval_re test1(code);
    double res1 = test1.make(s1);

    eval_im test2(code);
    double res2 = test2.make(s1);

    eval_conj test3(code);
    double res3 = test3.make(s1);

    eval_abs test4(code);
    double res4 = test4.make(s1);

    eval_abs2 test5(code);
    double res5 = test5.make(s1);

    eval_arg test6(code);
    double res6 = test6.make(s1);

    eval_angle test7(code);
    double res7 = test7.make(s1);

    return res1 + res2 + res3 + res4 + res5 + res6 + res7;
};

struct eval_isnan : eval_function<eval_isnan>
{
    eval_isnan(Integer c) : eval_function(c,false,true){};

    template<class T1>
    static auto eval(const T1& a1) -> bool
    {
        return is_nan(a1);
    };
};

struct eval_isinf : eval_function<eval_isinf>
{
    eval_isinf(Integer c) : eval_function(c,false,true){};

    template<class T1>
    static auto eval(const T1& a1) -> bool
    {
        return is_inf(a1);
    };
};

struct eval_isfin : eval_function<eval_isfin>
{
    eval_isfin(Integer c) : eval_function(c,false,true){};

    template<class T1>
    static auto eval(const T1& a1) -> bool
    {
        return is_finite(a1);
    };
};

struct eval_isreg : eval_function<eval_isreg>
{
    eval_isreg(Integer c) : eval_function(c,false,true){};

    template<class T1>
    static auto eval(const T1& a1) -> bool
    {
        return is_regular(a1);
    };
};

struct eval_isnorm : eval_function<eval_isnorm>
{
    eval_isnorm(Integer c) : eval_function(c,false,true){};

    template<class T1>
    static auto eval(const T1& a1) -> bool
    {
        return is_normal(a1);
    };
};

struct eval_isint : eval_function<eval_isint>
{
    eval_isint(Integer c) : eval_function(c,false,true){};

    template<class T1>
    static auto eval(const T1& a1) -> bool
    {
        return is_int(a1);
    };
};

struct eval_isreal : eval_function<eval_isreal>
{
    eval_isreal(Integer c) : eval_function(c,false,true){};

    template<class T1>
    static auto eval(const T1& a1) -> bool
    {
        return is_real(a1);
    };
};

struct eval_iszero : eval_function<eval_iszero>
{
    eval_iszero(Integer c) : eval_function(c,false,true){};

    template<class T1>
    static auto eval(const T1& a1) -> bool
    {
        return is_zero(a1);
    };
};

struct eval_op_neg : eval_function<eval_op_neg>
{
    eval_op_neg(Integer c) : eval_function(c,false,true){};

    template<class T1>
    static auto eval(const T1& a1) -> bool
    {
        return neg(a1);
    };
};

struct eval_op_true : eval_function<eval_op_true>
{
    eval_op_true(Integer c) : eval_function(c,false,true){};

    template<class T1>
    static auto eval(const T1& a1) -> bool
    {
        return is_true(a1);
    };
};

struct eval_isone : eval_function<eval_isone>
{
    eval_isone(Integer c) : eval_function(c,false,true){};

    template<class T1>
    static auto eval(const T1& a1) -> bool
    {
        return is_one(a1);
    };
};

double gmp_tester::test_is(const Scalar& s1, Integer code)
{
    eval_isnan test1(code);
    double res1 = test1.make(s1);

    eval_isinf test2(code);
    double res2 = test2.make(s1);

    eval_isfin test3(code);
    double res3 = test3.make(s1);

    eval_isreg test4(code);
    double res4 = test4.make(s1);

    eval_isint test5(code);
    double res5 = test5.make(s1);

    eval_isreal test6(code);
    double res6 = test6.make(s1);

    eval_iszero test7(code);
    double res7 = test7.make(s1);

    eval_isone test8(code);
    double res8 = test8.make(s1);

    eval_isnorm test9(code);
    double res9 = test9.make(s1);

    eval_op_neg test10(code);
    double res10 = test10.make(s1);

    eval_op_true test11(code);
    double res11 = test11.make(s1);

    return res1 + res2 + res3 + res4 + res5 + res6 + res7 + res8 + res9
            + res10 + res11;
};

struct eval_nextabove   : eval_function<eval_nextabove>
                        , Notcomplex_func
{
    using Notcomplex_func::eval;

    eval_nextabove(Integer c) : eval_function(c,true,false){};

    template<class T1>
    static auto eval(const T1& a1) -> decltype(nextabove(a1))
    {
        return nextabove(a1);
    };
};
struct eval_nextbelow   : eval_function<eval_nextbelow>
                        , Notcomplex_func
{
    using Notcomplex_func::eval;

    eval_nextbelow(Integer c) : eval_function(c,true,false){};

    template<class T1>
    static auto eval(const T1& a1) -> decltype(nextbelow(a1))
    {
        return nextbelow(a1);
    };
};
struct eval_signbit2 : eval_function<eval_signbit2>
{
    eval_signbit2(Integer c) : eval_function(c,false,true){};

    template<class T1>
    static auto eval(const T1& a1) -> decltype(signbit(a1))
    {
        return signbit(a1);
    };

    static Real eval(const Complex&)        { return 0.0; };
    static Real eval(const Float_complex&)  { return 0.0; };
    static Real eval(const mp_complex&)     { return 0.0; };
};
struct eval_eps : eval_function<eval_eps>
{
    eval_eps(Integer c) : eval_function(c,true,false){};

    template<class T1>
    static auto eval(const T1& a1) -> decltype(eps(a1))
    {
        return eps(a1);
    };
};

double gmp_tester::test_next(const Scalar& s1, Integer code)
{
    eval_nextabove test1(code);
    double res1 = test1.make(s1);

    eval_nextbelow test2(code);
    double res2 = test2.make(s1);

    return res1 + res2;
};

double gmp_tester::test_sign(const Scalar& s1, Integer code)
{
    eval_signbit2 test1(code);
    double res1 = test1.make(s1);

    return res1;
};

double gmp_tester::test_eps(const Scalar& s1, Integer code)
{
    eval_eps test1(code);
    double res1 = test1.make(s1);

    return res1;
};

template<class Func>
struct eval_scalar_func_templ : eval_function_num<eval_scalar_func_templ<Func>>
{
    eval_scalar_func_templ(Integer c, Real mult_c, Real mult_r, bool naninf_allowed)
        : eval_function_num(c, mult_c, mult_r, naninf_allowed){};

    template<class T1>
    static auto eval(const T1& a1) -> decltype(Func::eval(a1))
    {
        return Func::eval(a1);
    };
};

template<class Func>
double gmp_tester::test_scalar(const Scalar& s, Integer code)
{
    bool naninf_allowed
                = (std::is_same<Func,IFloor_func>::value || std::is_same<Func,ICeil_func>::value
                || std::is_same<Func,ITrunc_func>::value || std::is_same<Func,IRound_func>::value) == false;

    bool is_inacc_compl
                = std::is_same<Func,Log_func>::value || std::is_same<Func,Log_c_func>::value
                ||std::is_same<Func,Log2_func>::value || std::is_same<Func,Log2_c_func>::value
                || std::is_same<Func,Log10_func>::value || std::is_same<Func,Log10_c_func>::value
                || std::is_same<Func,Log1p_func>::value || std::is_same<Func,Log1p_c_func>::value
                || std::is_same<Func,Exp2_func>::value || std::is_same<Func,Exp10_func>::value
                || std::is_same<Func,Expm1_func>::value|| std::is_same<Func,Acot_func>::value
                || std::is_same<Func,Asec_func>::value  || std::is_same<Func,Acoth_func>::value
                || std::is_same<Func,Acoth_c_func>::value|| std::is_same<Func,Asech_func>::value
                || std::is_same<Func,Asech_c_func>::value;

    bool is_inacc_compl2
                = std::is_same<Func,Asec_func>::value || std::is_same<Func,Asech_func>::value
                || std::is_same<Func,Asech_c_func>::value
                || std::is_same<Func,Acoth_c_func>::value   || std::is_same<Func,Acoth_func>::value;

    bool is_inacc_re
                =  std::is_same<Func,Exp2_func>::value || std::is_same<Func,Exp10_func>::value
                || std::is_same<Func,Acoth_c_func>::value   || std::is_same<Func,Acoth_func>::value
                || std::is_same<Func,Asech_c_func>::value;

    bool is_exp10   = std::is_same<Func,Exp10_func>::value;

    Real mult_re    = 1;
    Real mult_compl = 1;
    
    if (is_inacc_compl)
        mult_compl  = 10;
    if (is_inacc_compl2)
        mult_compl  = mult_compl * 10;

    if (is_inacc_re)
        mult_re     = 10;

    if (is_exp10)
    {
        //highly inaccurate if imaginary part is large
        mult_compl  = 5000;
        mult_re     = 100;
    };

    eval_scalar_func_templ<Func> test1(code, mult_compl, mult_re, naninf_allowed);
    double res1 = test1.make(s);
    return res1;
};

void gmp_tester::test_ldexp_p3()
{
    return test_scalar_func<Ldexp3p_func_t>();
}
void gmp_tester::test_ldexp_m3()
{
    return test_scalar_func<Ldexp3m_func_t>();
}
void gmp_tester::test_scalbn_p3()
{
    return test_scalar_func<Scalbn3p_func_t>();
}
void gmp_tester::test_scalbn_m3()
{
    return test_scalar_func<Scalbn3m_func_t>();
}
void gmp_tester::test_frexp()
{
    return test_scalar_func<Frexp_func_t>();
}
void gmp_tester::test_modf_frac()
{
    return test_scalar_func<Modf_frac_func_t>();
}
void gmp_tester::test_modf_int()
{
    return test_scalar_func<Modf_int_func_t>();
}
void gmp_tester::test_logb()
{
    return test_scalar_func<Logb_func>();
}
void gmp_tester::test_ilogb()
{
    return test_scalar_func<Ilogb_func>();
}

void gmp_tester::test_inv()
{
    return test_scalar_func<Inv_func>();
}
void gmp_tester::test_invs()
{
    return test_scalar_func<Invs_func>();
}

void gmp_tester::test_exp()
{
    return test_scalar_func<Exp_func>();
}
void gmp_tester::test_expm1()
{
    return test_scalar_func<Expm1_func>();
}
void gmp_tester::test_expi()
{
    return test_scalar_func<Expi_func>();
}

void gmp_tester::test_sqrt()
{
    return test_scalar_func<Sqrt_func>();
}
void gmp_tester::test_cbrt()
{
    return test_scalar_func<Cbrt_func>();
}

void gmp_tester::test_sqrt_c()
{
    return test_scalar_func<Sqrt_c_func>();
}

void gmp_tester::test_exp2()
{
    return test_scalar_func<Exp2_func>();
}
void gmp_tester::test_exp10()
{
    return test_scalar_func<Exp10_func>();
}
void gmp_tester::test_log()
{
    return test_scalar_func<Log_func>();
}
void gmp_tester::test_log1p()
{
    return test_scalar_func<Log1p_func>();
}
void gmp_tester::test_log2()
{
    return test_scalar_func<Log2_func>();
}
void gmp_tester::test_log10()
{
    return test_scalar_func<Log10_func>();
}
void gmp_tester::test_log_c()
{
    return test_scalar_func<Log_c_func>();
}
void gmp_tester::test_log1p_c()
{
    return test_scalar_func<Log1p_c_func>();
}
void gmp_tester::test_log2_c()
{
    return test_scalar_func<Log2_c_func>();
}
void gmp_tester::test_log10_c()
{
    return test_scalar_func<Log10_c_func>();
}

void gmp_tester::test_sin()
{
    return test_scalar_func<Sin_func>();
}
void gmp_tester::test_cos()
{
    return test_scalar_func<Cos_func>();
}
void gmp_tester::test_tan()
{
    return test_scalar_func<Tan_func>();
}
void gmp_tester::test_cot()
{
    return test_scalar_func<Cot_func>();
}
void gmp_tester::test_sec()
{
    return test_scalar_func<Sec_func>();
}
void gmp_tester::test_csc()
{
    return test_scalar_func<Csc_func>();
}
void gmp_tester::test_sinh()
{
    return test_scalar_func<Sinh_func>();
}
void gmp_tester::test_cosh()
{
    return test_scalar_func<Cosh_func>();
}
void gmp_tester::test_tanh()
{
    return test_scalar_func<Tanh_func>();
}
void gmp_tester::test_coth()
{
    return test_scalar_func<Coth_func>();
}
void gmp_tester::test_sech()
{
    return test_scalar_func<Sech_func>();
}
void gmp_tester::test_csch()
{
    return test_scalar_func<Csch_func>();
}
void gmp_tester::test_asin()
{
    return test_scalar_func<Asin_func>();
}
void gmp_tester::test_asin_c()
{
    return test_scalar_func<Asin_c_func>();
}
void gmp_tester::test_acos()
{
    return test_scalar_func<Acos_func>();
}
void gmp_tester::test_acos_c()
{
    return test_scalar_func<Acos_c_func>();
}

void gmp_tester::test_atan()
{
    return test_scalar_func<Atan_func>();
}

void gmp_tester::test_asinh()
{
    return test_scalar_func<Asinh_func>();
}
void gmp_tester::test_acosh()
{
    return test_scalar_func<Acosh_func>();
}
void gmp_tester::test_acosh_c()
{
    return test_scalar_func<Acosh_c_func>();
}
void gmp_tester::test_atanh()
{
    return test_scalar_func<Atanh_func>();
}
void gmp_tester::test_atanh_c()
{
    return test_scalar_func<Atanh_c_func>();
}

void gmp_tester::test_floor()
{
    return test_scalar_func<Floor_func>();
}
void gmp_tester::test_ceil()
{
    return test_scalar_func<Ceil_func>();
}
void gmp_tester::test_round()
{
    return test_scalar_func<Round_func>();
}
void gmp_tester::test_combinatorics()
{
    test_scalar_func<Factorial_func>();
    test_scalar_func<Double_factorial_func>();
    test_scalar_func<Binomial_coefficient_func>();
}
void gmp_tester::test_combinatorics2()
{
    double res = 0;

    {
        Real res1   = rising_factorial(5,3);
        Real res2   = rising_factorial(Real(5),3);
        Float res3  = rising_factorial(Float(5),3);

        if (res1 != res2 || abs(res2 - res3) > eps(res3))
            res     += 1;

        Real res1_2 = dynamic::rising_factorial(OInteger(5), 3).get();
        Real res2_2 = dynamic::rising_factorial(OReal(Real(5)), 3).get();
        Real res3_2 = dynamic::rising_factorial(OFloat(Float(5)), 3).get();

        if (res1_2 != res1 || res2_2 != res2 || res3_2 != res3)
            res     += 1;

        Object ores1_3 = dynamic::rising_factorial(Object(OInteger(5)), 3);
        Object ores2_3 = dynamic::rising_factorial(Object(OReal(Real(5))), 3);
        Object ores3_3 = dynamic::rising_factorial(Object(OFloat(Float(5))), 3);

        Real res1_3     = OReal(ores1_3, dynamic::from_object()).get();
        Real res2_3     = OReal(ores2_3, dynamic::from_object()).get();
        Real res3_3     = OReal(ores3_3, dynamic::from_object()).get();

        if (res1_3 != res1 || res2_3 != res2 || res3_3 != res3)
            res     += 1;
    }
    {
        Real res1   = falling_factorial(5,3);
        Real res2   = falling_factorial(Real(5),3);
        Real res3   = falling_factorial(Float(5),3);

        if (res1 != res2 || (Float)res1 != (Float)res3)
            res     += 1;

        Real res1_2 = dynamic::falling_factorial(OInteger(5), 3).get();
        Real res2_2 = dynamic::falling_factorial(OReal(Real(5)), 3).get();
        Real res3_2 = dynamic::falling_factorial(OFloat(Float(5)), 3).get();

        if (res1_2 != res1 || res2_2 != res2 || res3_2 != res3)
            res     += 1;

        Object ores1_3 = dynamic::falling_factorial(Object(OInteger(5)), 3);
        Object ores2_3 = dynamic::falling_factorial(Object(OReal(Real(5))), 3);
        Object ores3_3 = dynamic::falling_factorial(Object(OFloat(Float(5))), 3);

        Real res1_3     = OReal(ores1_3, dynamic::from_object()).get();
        Real res2_3     = OReal(ores2_3, dynamic::from_object()).get();
        Real res3_3     = OReal(ores3_3, dynamic::from_object()).get();

        if (res1_3 != res1 || res2_3 != res2 || res3_3 != res3)
            res     += 1;
    }
    {
        Real res1   = bernoulli_b2n<Integer>(6);
        Real res2   = bernoulli_b2n<Real>(6);
        Real res3   = bernoulli_b2n<Float>(6);

        if (res1 != res2 || (Float)res1 != (Float)res3)
            res     += 1;

        Real res1_2 = dynamic::bernoulli_b2n<OInteger>(6).get();
        Real res2_2 = dynamic::bernoulli_b2n<OReal>(6).get();
        Real res3_2 = dynamic::bernoulli_b2n<OFloat>(6).get();

        if (res1_2 != res1 || res2_2 != res2 || res3_2 != res3)
            res     += 1;

        Object ores1_3 = dynamic::bernoulli_b2n(OInteger::get_static_type(), 6);
        Object ores2_3 = dynamic::bernoulli_b2n(OReal::get_static_type(), 6);
        Object ores3_3 = dynamic::bernoulli_b2n(OFloat::get_static_type(), 6);

        Real res1_3     = OReal(ores1_3, dynamic::from_object()).get();
        Real res2_3     = OReal(ores2_3, dynamic::from_object()).get();
        Real res3_3     = OReal(ores3_3, dynamic::from_object()).get();

        if (res1_3 != res1 || res2_3 != res2 || res3_3 != res3)
            res     += 1;
    }
    {
        Integer res1   = max_bernoulli_b2n<Integer>();
        Integer res2   = max_bernoulli_b2n<Real>();
        Integer res3   = max_bernoulli_b2n<Float>();

        if (res1 != res2 || res1 < res3)
            res     += 1;

        Integer res1_2 = dynamic::max_bernoulli_b2n<OInteger>();
        Integer res2_2 = dynamic::max_bernoulli_b2n<OReal>();
        Integer res3_2 = dynamic::max_bernoulli_b2n<OFloat>();

        if (res1_2 != res1 || res2_2 != res2 || res3_2 != res3)
            res     += 1;

        Integer res1_3 = dynamic::max_bernoulli_b2n(OInteger::get_static_type());
        Integer res2_3 = dynamic::max_bernoulli_b2n(OReal::get_static_type());
        Integer res3_3 = dynamic::max_bernoulli_b2n(OFloat::get_static_type());

        if (res1_3 != res1 || res2_3 != res2 || res3_3 != res3)
            res     += 1;
    }

    if (res == 0.0)
        std::cout << "combinatorics2: ok" << "\n";
    else
        std::cout << "combinatorics2" << ": FAILED" << "\n";
};

void gmp_tester::test_trunc()
{
    return test_scalar_func<Trunc_func>();
}
void gmp_tester::test_ifloor()
{
    return test_scalar_func<IFloor_func>();
}
void gmp_tester::test_iceil()
{
    return test_scalar_func<ICeil_func>();
}
void gmp_tester::test_iround()
{
    return test_scalar_func<IRound_func>();
}
void gmp_tester::test_itrunc()
{
    return test_scalar_func<ITrunc_func>();
}
void gmp_tester::test_sign()
{
    return test_scalar_func<Sign_func>();
}
void gmp_tester::test_isign()
{
    return test_scalar_func<ISign_func>();
}
void gmp_tester::test_fpclassify()
{
    return test_scalar_func<Fpclassify_func>();
}

}};
