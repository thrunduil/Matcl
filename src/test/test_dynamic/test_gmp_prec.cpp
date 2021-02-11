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

#include "test_gmp_prec.h"
#include "matcl-scalar/lib_functions/func_forwarding.h"

#include <iostream>
#include <iomanip>
#include "boost/io/ios_state.hpp"
#include "matcl-scalar/lib_functions/scalar_gen.h"
#include "matcl-scalar/lib_functions/func_unary.h"
#include "matcl-scalar/lib_functions/func_binary.h"
#include "matcl-scalar/IO/scalar_io.h"

#pragma warning(push)
#pragma warning(disable:4127) // conditional expression is constant

namespace matcl { namespace test
{

void test_gmp_prec(std::ostream& os)
{
    gmp_tester_prec test;
    test.make(os);
};

struct Default_bin_rand
{
    using prec_vec  = std::vector<precision>;
    using max_vec   = std::vector<double>;

    static void rand(const prec_vec& prec, const max_vec& max, mp_complex& s1,
                     mp_complex& s2, mp_float& s3)
    {
        s1 = gmp_tester_prec::rand_scalar_c(prec, max);
        s2 = gmp_tester_prec::rand_scalar_c(prec, max);
        s3 = gmp_tester_prec::rand_scalar(prec, max);
    };
};

struct Rand_pow_cc
{
    using prec_vec  = std::vector<precision>;
    using max_vec   = std::vector<double>;

    static void rand(const prec_vec& prec, const max_vec& max, mp_complex& s1,
                     mp_complex& s2, mp_float& s3)
    {
        return gmp_tester_prec::rand_scalar_c_pow_cc(prec, max, s1, s2, s3);
    };
};

struct Rand_pow_cr
{
    using prec_vec  = std::vector<precision>;
    using max_vec   = std::vector<double>;

    static void rand(const prec_vec& prec, const max_vec& max, mp_complex& s1,
                     mp_complex& s2, mp_float& s3)
    {
        return gmp_tester_prec::rand_scalar_c_pow_cr(prec, max, s1, s2, s3);
    };
};

struct Rand_pow_rc
{
    using prec_vec  = std::vector<precision>;
    using max_vec   = std::vector<double>;

    static void rand(const prec_vec& prec, const max_vec& max, mp_complex& s1,
                     mp_complex& s2, mp_float& s3)
    {
        return gmp_tester_prec::rand_scalar_c_pow_rc(prec, max, s1, s2, s3);
    };
};

struct Rand_log
{
    static mp_complex rand_scalar_c(precision prec, double max_re, double max_im)
    {
        (void)max_im;
        return gmp_tester_prec::rand_scalar_c_log(prec, max_re);
    };
};
struct Rand_log1p
{
    static mp_complex rand_scalar_c(precision prec, double max_re, double max_im)
    {
        (void)max_im;
        return gmp_tester_prec::rand_scalar_c_log1p(prec, max_re);
    };
};
struct Rand_expm1
{
    static mp_complex rand_scalar_c(precision prec, double max_re, double max_im)
    {
        (void)max_im;
        return gmp_tester_prec::rand_scalar_c_expm1(prec, max_re);
    };
};

struct Default_crnd
{
    static mp_complex rand_scalar_c(precision prec, double max_re, double max_im)
    {
        return gmp_tester_prec::rand_scalar_c(prec, max_re, max_im);
    };
};

struct Exp2_crnd
{
    static mp_complex rand_scalar_c(precision prec, double max_re, double max_im)
    {
        return gmp_tester_prec::rand_scalar_c_exp2(prec, max_re, max_im);
    };
};

struct Exp10_crnd
{
    static mp_complex rand_scalar_c(precision prec, double max_re, double max_im)
    {
        return gmp_tester_prec::rand_scalar_c_exp10(prec, max_re, max_im);
    };
};

struct Pi_const
{
    static std::string  name()              { return "pi"; };
    static mp_float     eval(precision p)   { return constants::mp_pi(p); };
};
struct Pi_2_const
{
    static std::string  name()              { return "pi_2"; };
    static mp_float     eval(precision p)   { return constants::mp_pi_2(p); };
};

struct Pi_4_const
{
    static std::string  name()              { return "pi_4"; };
    static mp_float     eval(precision p)   { return constants::mp_pi_4(p); };
};

struct E_const
{
    static std::string  name()              { return "e"; };
    static mp_float     eval(precision p)   { return constants::mp_e(p); };
};

struct Ln2_const
{
    static std::string  name()              { return "ln2"; };
    static mp_float     eval(precision p)   { return constants::mp_ln2(p); };
};

struct Ln10_const
{
    static std::string  name()              { return "ln10"; };
    static mp_float     eval(precision p)   { return constants::mp_ln10(p); };
};

struct Log2e_const
{
    static std::string  name()              { return "log2e"; };
    static mp_float     eval(precision p)   { return constants::mp_log2e(p); };
};

struct Log10e_const
{
    static std::string  name()              { return "log10e"; };
    static mp_float     eval(precision p)   { return constants::mp_log10e(p); };
};

void gmp_tester_prec::make(std::ostream& os)
{
    Integer N       = 500;    
    Integer Nb      = 500;

  #ifdef _DEBUG
    Integer Nc      = 10;
  #else
    Integer Nc      = 100;
  #endif

    prec_vec prec   = {precision(10), precision(23), precision(53), 
                        precision(100), precision(500)};
    max_vec max     = {1.0, 1.0e-1, 1.0e1, 1.0e-2, 1.0e2, 1.0e-4, 1.0e4, 1.0e-8, 1.0e8,
                        1.0e-16, 1.0e16, 1.0e-32, 1.0e32};

    max_vec max_b   = {1.0, 1.0, 1.0e-1, 1.0e1, 1.0e-1, 1.0e1, 1.0e-2, 1.0e2, 1.0e-2, 1.0e2, 
                        1.0e-4, 1.0e4, 1.0e-8, 1.0e8, 1.0e-16, 1.0e16, 1.0e-32, 1.0e32};

    test_complex_bin_special(os, N, prec);    
    test_complex_special(os, N, prec);    
    test_bin_real(os, Nb, prec, max_b);
    test_bin_complex(os, Nb, prec, max_b);
    test_constants(os);    
    test_complex(os, Nc, prec, max);            
    test_real(os, N, prec, max);      
};

void gmp_tester_prec::test_constants(std::ostream& os)
{
    for (size_t prec = 10; prec < 500; prec += 5)
        test_constants(os, precision(prec));

    os << "test constants finished" << "\n";
};

void gmp_tester_prec::test_real(std::ostream& os, Integer N, const prec_vec& prec_v, 
                              const max_vec& max_v)
{
    for (const auto& prec : prec_v)
    for (const auto& max : max_v)
        test_real(os, N, prec, max);

    os << "test real finished" << "\n";
};

void gmp_tester_prec::test_bin_real(std::ostream& os, Integer N, const prec_vec& prec_v, 
                              const max_vec& max_v)
{
    for (const auto& prec : prec_v)
        test_bin_real(os, N, prec, max_v);

    os << "test bin real finished" << "\n";
};

void gmp_tester_prec::test_bin_complex(std::ostream& os, Integer N, const prec_vec& prec_v, 
                              const max_vec& max_v)
{
    for (const auto& prec : prec_v)
        test_bin_complex(os, N, prec, prec_v, max_v);

    os << "test bin complex finished" << "\n";
};

void gmp_tester_prec::test_bin_complex(std::ostream& os_res, Integer N, precision prec,
                                       const prec_vec& prec_v, const max_vec& max)
{
    std::ostringstream os;

    os << "\n";
    os << std::string(50,'-') << "\n";
    os << "COMPLEX BIN  PRECISION: " << prec << "\n";
    os << std::string(50,'-') << "\n";
    os << "\n";

    bool error = false;
    error |= test_bin_func_c<Plus_func, Default_bin_rand>(os, N, prec, prec_v, max, 0.5);
    error |= test_bin_func_c<Minus_func, Default_bin_rand>(os, N, prec, prec_v, max, 0.5);    

    error |= test_bin_func_c<Mult_func, Default_bin_rand>(os, N, prec, prec_v, max, 1.0);
    error |= test_bin_func_c<Div_func, Default_bin_rand>(os, N, prec, prec_v, max, 1.0);
    error |= test_bin_func_c<Hypot_func, Default_bin_rand>(os, N, prec, prec_v, max, 1.0);
    error |= test_bin_func_c<Pow_func, Default_bin_rand>(os, N, prec, prec_v, max, 1.0);
    error |= test_bin_func_c<Pow_rc_func, Default_bin_rand>(os, N, prec, prec_v, max, 1.0);
    error |= test_bin_func_c<Pow_cr_func, Default_bin_rand>(os, N, prec, prec_v, max, 1.0);

    if (error)
        os_res << os.str();
};

void gmp_tester_prec::test_complex(std::ostream& os, Integer N, const prec_vec& prec_v, 
                              const max_vec& max_v)
{
    for (const auto& prec : prec_v)
    for (const auto& max_re : max_v)
    for (const auto& max_im : max_v)
        test_complex(os, N, prec, max_re, max_im);

    os << "test complex finished" << "\n";
};

void gmp_tester_prec::test_complex_special(std::ostream& os, Integer N, const prec_vec& prec_v)
{
    max_vec max     = {1.0, 1.0e-1, 1.0e-2, 1.0e-4, 1.0e-8, 1.0e-16, 1.0e-32, 1.0e-64, 1.0e-128};

    for (const auto& prec : prec_v)
    for (const auto& max_v : max)
        test_complex_special(os, N, prec, max_v);

    os << "test complex special finished" << "\n";
};

void gmp_tester_prec::test_complex_bin_special(std::ostream& os, Integer N, const prec_vec& prec_v)
{
    max_vec max     = {1.0, 1.0e-1, 1.0e-2, 1.0e-4, 1.0e-8, 1.0e-16, 1.0e-32, 1.0e-64, 1.0e-128};

    for (const auto& prec : prec_v)
        test_complex_bin_special(os, N, prec, prec_v, max);

    os << "test complex bin special finished" << "\n";
};

bool gmp_tester_prec::test_real(std::ostream& os_res, Integer N, precision prec, double max)
{
    std::ostringstream os;

    os << "\n";
    os << std::string(50,'-') << "\n";
    os << "REAL     PRECISION: " << prec << ", MAX: " << max << "\n";
    os << std::string(50,'-') << "\n";
    os << "\n";

    bool error = false;
    error |= test_scalar_func<Ldexp3p_func_t>(os, N, prec, max, 0.0);
    error |= test_scalar_func<Ldexp3m_func_t>(os, N, prec, max, 0.0);
    error |= test_scalar_func<Frexp_func_t>(os, N, prec, max, 0.0);
    error |= test_scalar_func<Modf_frac_func_t>(os, N, prec, max, 0.0);
    error |= test_scalar_func<Modf_int_func_t>(os, N, prec, max, 0.0);
    error |= test_scalar_func<Sign_func>(os, N, prec, max, 0.0);

    error |= test_scalar_func<Sqrt_func>(os, N, prec, max, 0.5);
    error |= test_scalar_func<Cbrt_func>(os, N, prec, max, 0.5);
    error |= test_scalar_func<Exp_func>(os, N, prec, max, 0.5);
    error |= test_scalar_func<Exp2_func>(os, N, prec, max, 0.5);
    error |= test_scalar_func<Exp10_func>(os, N, prec, max, 0.5);
    error |= test_scalar_func<Expm1_func>(os, N, prec, max, 0.5);
    error |= test_scalar_func<Log_func>(os, N, prec, max, 0.5);
    error |= test_scalar_func<Log2_func>(os, N, prec, max, 0.5);
    error |= test_scalar_func<Log10_func>(os, N, prec, max, 0.5);
    error |= test_scalar_func<Log1p_func>(os, N, prec, max, 0.5);
    error |= test_scalar_func<Sin_func>(os, N, prec, max, 0.5);
    error |= test_scalar_func<Cos_func>(os, N, prec, max, 0.5);
    error |= test_scalar_func<Tan_func>(os, N, prec, max, 0.5);
    error |= test_scalar_func<Cot_func>(os, N, prec, max, 0.5);
    error |= test_scalar_func<Sec_func>(os, N, prec, max, 0.5);
    error |= test_scalar_func<Csc_func>(os, N, prec, max, 0.5);
    error |= test_scalar_func<Sinh_func>(os, N, prec, max, 0.5);
    error |= test_scalar_func<Cosh_func>(os, N, prec, max, 0.5);
    error |= test_scalar_func<Tanh_func>(os, N, prec, max, 0.5);
    error |= test_scalar_func<Coth_func>(os, N, prec, max, 0.5);
    error |= test_scalar_func<Sech_func>(os, N, prec, max, 0.5);
    error |= test_scalar_func<Csch_func>(os, N, prec, max, 0.5);
    error |= test_scalar_func<Asin_func>(os, N, prec, max, 0.5);
    error |= test_scalar_func<Acos_func>(os, N, prec, max, 0.5);
    error |= test_scalar_func<Atan_func>(os, N, prec, max, 0.5);
    error |= test_scalar_func<Asinh_func>(os, N, prec, max, 0.5);
    error |= test_scalar_func<Acosh_func>(os, N, prec, max, 0.5);    
    error |= test_scalar_func<Atanh_func>(os, N, prec, max, 0.5);    
    error |= test_scalar_func<Uminus_func>(os, N, prec, max, 0.5);
    error |= test_scalar_func<Abs_func>(os, N, prec, max, 0.5);
    error |= test_scalar_func<Abs2_func>(os, N, prec, max, 0.5);
    error |= test_scalar_func<Arg_func>(os, N, prec, max, 0.5);
    error |= test_scalar_func<Inv_func>(os, N, prec, max, 0.5);

    error |= test_scalar_func<Real_func>(os, N, prec, max, 0.0);
    error |= test_scalar_func<Imag_func>(os, N, prec, max, 0.0);
    error |= test_scalar_func<Conj_func>(os, N, prec, max, 0.0);
    error |= test_scalar_func<Ldexp3p_func_t>(os, N, prec, max, 0.0);
    error |= test_scalar_func<Ldexp3m_func_t>(os, N, prec, max, 0.0);
    error |= test_scalar_func<Frexp_func_t>(os, N, prec, max, 0.0);

    if (error)
        os_res << os.str();

    return error;
};

bool gmp_tester_prec::test_constants(std::ostream& os_res, precision prec)
{
    std::ostringstream os;

    os << "\n";
    os << std::string(50,'-') << "\n";
    os << "REAL     PRECISION: " << prec << "\n";
    os << std::string(50,'-') << "\n";
    os << "\n";

    bool error = false;
    error |= test_constant<Pi_const>(os, prec, 0.5);
    error |= test_constant<Pi_2_const>(os, prec, 0.5);
    error |= test_constant<Pi_4_const>(os, prec, 0.5);
    error |= test_constant<E_const>(os, prec, 0.5);
    error |= test_constant<Ln2_const>(os, prec, 0.5);

    error |= test_constant<Ln10_const>(os, prec, 1.0);
    error |= test_constant<Log2e_const>(os, prec, 1.0);
    error |= test_constant<Log10e_const>(os, prec, 1.0);

    if (error)
        os_res << os.str();

    return error;
};

bool gmp_tester_prec::test_bin_real(std::ostream& os_res, Integer N, precision prec, const max_vec& max)
{
    std::ostringstream os;

    os << "\n";
    os << std::string(50,'-') << "\n";
    os << "REAL BIN, PRECISION: " << prec << "\n";
    os << std::string(50,'-') << "\n";
    os << "\n";

    bool error = false;
    error |= test_bin_func<Plus_func>(os, N, prec, max, 0.5);
    error |= test_bin_func<Minus_func>(os, N, prec, max, 0.5);
    error |= test_bin_func<Mult_func>(os, N, prec, max, 0.5);
    error |= test_bin_func<Div_func>(os, N, prec, max, 0.5);
    error |= test_bin_func<Hypot_func>(os, N, prec, max, 0.5);
    error |= test_bin_func<Atan2_func>(os, N, prec, max, 0.5);    
    error |= test_bin_func<Rem_func>(os, N, prec, max, 0.5);
    error |= test_bin_func<Mod_func>(os, N, prec, max, 1.0);
    error |= test_bin_func<Pow_func>(os, N, prec, max, 1.0);

    if (error)
        os_res << os.str();

    return error;
};

bool gmp_tester_prec::test_complex(std::ostream& os_res, Integer N, precision prec, 
                                   double max_re, double max_im)
{
    std::ostringstream os;

    os << "\n";
    os << std::string(50,'-') << "\n";
    os << "COMPLEX      PRECISION: " << prec << ", MAX RE: " << max_re<< ", MAX IM: " << max_im << "\n";
    os << std::string(50,'-') << "\n";
    os << "\n";

    bool error = false;

    error |= test_scalar_func_c<Real_func,Default_crnd>(os, N, prec, max_re, max_im, 0.0);
    error |= test_scalar_func_c<Imag_func,Default_crnd>(os, N, prec, max_re, max_im, 0.0);
    error |= test_scalar_func_c<Conj_func,Default_crnd>(os, N, prec, max_re, max_im, 0.0);
    error |= test_scalar_func_c<Ldexp3p_func_t,Default_crnd>(os, N, prec, max_re, max_im, 0.0);
    error |= test_scalar_func_c<Ldexp3m_func_t,Default_crnd>(os, N, prec, max_re, max_im, 0.0);
    error |= test_scalar_func_c<Frexp_func_t,Default_crnd>(os, N, prec, max_re, max_im, 0.0);
    error |= test_scalar_func_c<Modf_frac_func_t,Default_crnd>(os, N, prec, max_re, max_im, 0.0);
    error |= test_scalar_func_c<Modf_int_func_t,Default_crnd>(os, N, prec, max_re, max_im, 0.0);

    error |= test_scalar_func_c<Uminus_func,Default_crnd>(os, N, prec, max_re, max_im, 0.5);
    error |= test_scalar_func_c<Abs_func,Default_crnd>(os, N, prec, max_re, max_im, 0.5);
    error |= test_scalar_func_c<Asin_func,Default_crnd>(os, N, prec, max_re, max_im, 1.0);
    error |= test_scalar_func_c<Asinh_func,Default_crnd>(os, N, prec, max_re, max_im, 1.0);
    error |= test_scalar_func_c<Acos_func,Default_crnd>(os, N, prec, max_re, max_im, 1.0);
    error |= test_scalar_func_c<Acosh_func,Default_crnd>(os, N, prec, max_re, max_im, 1.0);
    error |= test_scalar_func_c<Atanh_func,Default_crnd>(os, N, prec, max_re, max_im, 1.0);
    error |= test_scalar_func_c<Atan_func,Default_crnd>(os, N, prec, max_re, max_im, 1.0);
    error |= test_scalar_func_c<Sign_func,Default_crnd>(os, N, prec, max_re, max_im, 1.0);
    error |= test_scalar_func_c<Sqrt_func,Default_crnd>(os, N, prec, max_re, max_im, 1.0);
    error |= test_scalar_func_c<Exp_func,Default_crnd>(os, N, prec, max_re, max_im, 1.0);
    error |= test_scalar_func_c<Abs2_func,Default_crnd>(os, N, prec, max_re, max_im, 1.0);
    error |= test_scalar_func_c<Arg_func,Default_crnd>(os, N, prec, max_re, max_im, 1.0);
    error |= test_scalar_func_c<Inv_func,Default_crnd>(os, N, prec, max_re, max_im, 1.0);    
    error |= test_scalar_func_c<Sin_func,Default_crnd>(os, N, prec, max_re, max_im, 1.0);
    error |= test_scalar_func_c<Cos_func,Default_crnd>(os, N, prec, max_re, max_im, 1.0);
    error |= test_scalar_func_c<Sinh_func,Default_crnd>(os, N, prec, max_re, max_im, 1.0);
    error |= test_scalar_func_c<Cosh_func,Default_crnd>(os, N, prec, max_re, max_im, 1.0);    
    error |= test_scalar_func_c<Sec_func,Default_crnd>(os, N, prec, max_re, max_im, 1.0);
    error |= test_scalar_func_c<Csc_func,Default_crnd>(os, N, prec, max_re, max_im, 1.0);
    error |= test_scalar_func_c<Sech_func,Default_crnd>(os, N, prec, max_re, max_im, 1.0);
    error |= test_scalar_func_c<Csch_func,Default_crnd>(os, N, prec, max_re, max_im, 1.0);
    error |= test_scalar_func_c<Tanh_func,Default_crnd>(os, N, prec, max_re, max_im, 1.0);
    error |= test_scalar_func_c<Tan_func,Default_crnd>(os, N, prec, max_re, max_im, 1.0);
    error |= test_scalar_func_c<Cot_func,Default_crnd>(os, N, prec, max_re, max_im, 1.0);    
    error |= test_scalar_func_c<Coth_func,Default_crnd>(os, N, prec, max_re, max_im, 1.0);
    error |= test_scalar_func_c<Log_func,Default_crnd>(os, N, prec, max_re, max_im, 1.0);
    error |= test_scalar_func_c<Log1p_func,Default_crnd>(os, N, prec, max_re, max_im, 1.0);
    error |= test_scalar_func_c<Log2_func,Default_crnd>(os, N, prec, max_re, max_im, 1.0);
    error |= test_scalar_func_c<Log10_func,Default_crnd>(os, N, prec, max_re, max_im, 1.0);

    error |= test_scalar_func_c<Expm1_func,Default_crnd>(os, N, prec, max_re, max_im, 1.0);
    error |= test_scalar_func_c<Expi_func,Default_crnd>(os, N, prec, max_re, max_im, 1.0);
    error |= test_scalar_func_c<Exp2_func,Default_crnd>(os, N, prec, max_re, max_im, 1.0);
    error |= test_scalar_func_c<Exp2_func,Exp2_crnd>(os, N, prec, max_re, max_im, 1.0);
    error |= test_scalar_func_c<Exp10_func,Default_crnd>(os, N, prec, max_re, max_im, 1.0);
    error |= test_scalar_func_c<Exp10_func,Exp10_crnd>(os, N, prec, max_re, max_im, 1.0);

    if (error)
        os_res << os.str();

    return error;
};

bool gmp_tester_prec::test_complex_special(std::ostream& os_res, Integer N, precision prec, double max_v)
{
    std::ostringstream os;

    os << "\n";
    os << std::string(50,'-') << "\n";
    os << "COMPLEX SPEC  PRECISION: " << prec << ", MAX: " << max_v << "\n";
    os << std::string(50,'-') << "\n";
    os << "\n";

    bool error = false;
    error   |= test_scalar_func_c<Log_func, Rand_log>(os, N, prec, max_v, max_v, 1.0);
    error   |= test_scalar_func_c<Log1p_func, Rand_log1p>(os, N, prec, max_v, max_v, 1.0);
    error   |= test_scalar_func_c<Atanh_func, Rand_log>(os, N, prec, max_v, max_v, 1.0);
    error   |= test_scalar_func_c<Atan_func, Rand_log>(os, N, prec, max_v, max_v, 1.0);
    error   |= test_scalar_func_c<Expm1_func, Rand_expm1>(os, N, prec, max_v, max_v, 1.0);

    if (error)
        os_res << os.str();

    return error;
};

bool gmp_tester_prec::test_complex_bin_special(std::ostream& os_res, Integer N, precision prec,
                                               const prec_vec& prec_v, const max_vec& max_v)
{
    std::ostringstream os;

    os << "\n";
    os << std::string(50,'-') << "\n";
    os << "COMPLEX BIN SPEC  PRECISION: " << prec << "\n";
    os << std::string(50,'-') << "\n";
    os << "\n";

    bool error = false;
    error   |= test_bin_func_c<Pow_func, Rand_pow_cc>(os, N, prec, prec_v, max_v, 1.0);
    error   |= test_bin_func_c<Pow_rc_func, Rand_pow_rc>(os, N, prec, prec_v, max_v, 1.0);
    error   |= test_bin_func_c<Pow_cr_func, Rand_pow_cr>(os, N, prec, prec_v, max_v, 1.0);

    if (error)
        os_res << os.str();

    return error;
};

mp_float gmp_tester_prec::rand_scalar(precision prec, const max_vec& max)
{
    Integer n   = (Integer)max.size();
    Integer ver = matcl::abs(matcl::irand()) % n;
    return mp_float(matcl::randn() * max[ver], precision(prec + 4));
};

mp_float gmp_tester_prec::rand_scalar(const prec_vec& prec, const max_vec& max)
{
    Integer nm  = (Integer)max.size();
    Integer np  = (Integer)prec.size();
    Integer vm  = matcl::abs(matcl::irand()) % nm;
    Integer vp  = matcl::abs(matcl::irand()) % np;
    return mp_float(matcl::randn() * max[vm], prec[vp]);
};

void gmp_tester_prec::rand_scalar_c_bin(const prec_vec& prec, const max_vec& max, mp_complex& s1,
                    mp_complex& s2, mp_float& s3)
{
    s1 = gmp_tester_prec::rand_scalar_c(prec, max);
    s2 = gmp_tester_prec::rand_scalar_c(prec, max);
    s3 = gmp_tester_prec::rand_scalar(prec, max);
};
void gmp_tester_prec::rand_scalar_c_pow_cc(const prec_vec& prec, const max_vec& max, mp_complex& s1,
                    mp_complex& s2, mp_float& s3)
{
    Integer c       = abs(matcl::irand()) % 20;

    if (c == 0)
    {
        //a = -1 + 0i for x < 0;
        s2          = gmp_tester_prec::rand_scalar_c(prec, max);
        s1          = mp_complex(-1.0, 0.0, s2.get_precision());        
    }
    else if (c == 1)
    {
        //a == i
        s2          = gmp_tester_prec::rand_scalar_c(prec, max);
        s1          = mp_complex(0.0, 1.0, s2.get_precision());
    }
    else if (c == 2)
    {
        //a == -i
        s2          = gmp_tester_prec::rand_scalar_c(prec, max);
        s1          = mp_complex(0.0, -1.0, s2.get_precision());
    }
    else
    {
        //rand scalars such that log(z) * phi ~ pi/2 * k
        mp_complex z    = gmp_tester_prec::rand_scalar_c(prec, max);
        precision p1    = z.get_precision() + 10;
        Integer k       = matcl::irand() % 1000;
        mp_complex exp  = div(constants::mp_pi(p1) / 2 * k, log(z,p1), p1) * mp_complex(0.0,1.0,p1);

        s1 = z;
        s2 = exp;
    };

    s3 = gmp_tester_prec::rand_scalar(prec, max);
};

void gmp_tester_prec::rand_scalar_c_pow_rc(const prec_vec& prec, const max_vec& max, mp_complex& s1,
                    mp_complex& s2, mp_float& s3)
{
    //rand scalars such that log(z) * phi ~ pi/2 * k
    mp_float z      = abs(gmp_tester_prec::rand_scalar(prec, max));
    precision p1    = z.get_precision() + 10;
    Integer k       = matcl::irand() % 1000;    
    mp_complex phi  = div(constants::mp_pi(p1) * k, 2*log(z,p1), p1) * mp_complex(0.0, 1.0, p1);
    phi             = phi + gmp_tester_prec::rand_scalar_c(prec, max) * 1e-6;

    s1 = z;
    s2 = phi;
    s3 = gmp_tester_prec::rand_scalar(prec, max);
};
void gmp_tester_prec::rand_scalar_c_pow_cr(const prec_vec& prec, const max_vec& max, mp_complex& s1,
                    mp_complex& s2, mp_float& s3)
{
    Integer c       = abs(matcl::irand()) % 20;

    //s2 is int, imag(s1) = 0; real(s1) < 0
    //s2 is int, real(s1) = 0
    //s2 is int, abs(real(s1)) = abs(imag(s1))
    //s2*2 is int imag(s1) = 0

    if (c == 0)
    {
        s2          = Real(matcl::irand() % 20) * 0.5;
        s1          = gmp_tester_prec::rand_scalar(prec, max);
        s1          = mp_complex(0,real(s1),s1.get_precision());
    }
    else if (c == 1)
    {
        s2          = Real(matcl::irand() % 20) * 0.5;
        s1          = gmp_tester_prec::rand_scalar(prec, max);
    }
    else if (c == 2)
    {
        s2          = Real(matcl::irand() % 20) * 0.5;
        Real sign1  = matcl::irand() < 0 ? -1.0 : 1.0;
        Real sign2  = matcl::irand() < 0 ? -1.0 : 1.0;

        (void)sign1;
        (void)sign2;

        s1          = gmp_tester_prec::rand_scalar(prec, max);
        s1          = mp_complex(sign1 * real(s1), sign1 * real(s1),s1.get_precision());
    }
    else
    {
        //rand scalars such that log(z) * phi ~ pi/2 * k
        mp_complex z    = gmp_tester_prec::rand_scalar_c(prec, max);
        precision p1    = z.get_precision() + 10;
        Integer k       = matcl::irand() % 1000;
        mp_complex exp  = div(constants::mp_pi(p1) / 2 * k, imag(log(z,p1)), p1);

        s1 = z;
        s2 = exp;
    };
    s3          = gmp_tester_prec::rand_scalar(prec, max);
};

mp_complex gmp_tester_prec::rand_scalar_c(const prec_vec& prec, const max_vec& max)
{
    Integer np  = (Integer)prec.size();
    Integer vp  = matcl::abs(matcl::irand()) % np;
    return mp_complex(rand_scalar(prec[vp], max), rand_scalar(prec[vp], max));
};
mp_float gmp_tester_prec::rand_scalar(precision prec, double max)
{
    return mp_float(matcl::randn() * max, precision(prec + 4));
};
mp_complex  gmp_tester_prec::rand_scalar_c(precision prec, double max_re, double max_im)
{
    return mp_complex(rand_scalar(prec, max_re), rand_scalar(prec, max_im));
};
mp_complex  gmp_tester_prec::rand_scalar_c_log(precision prec, double max)
{
    //mp_complex z = mp_complex(1.0, 1e-1 * matcl::randn(), prec);
    //return z;

    //rand complex value z such that |z| = 1
    mp_float x = rand_scalar(prec, max);
    mp_float y = sqrt(1-x*x);

    int case_t  = matcl::abs(matcl::irand()) % 2;
    switch(case_t)
    {
        case 0:     return mp_complex(x,y);
        default:    return mp_complex(y,x);
    }
};
mp_complex  gmp_tester_prec::rand_scalar_c_exp2(precision prec, double max_re, double max_im)
{
    //rand complex value z such that im(ln(2) * z) = pi/2*k
    mp_float x  = rand_scalar(prec, max_re);
    mp_float kr = rand_scalar(precision(10), max(max_im,1.0));
    mp_float k;
    modf(kr, k);

    if (is_zero(k))
        k           = 1;

    precision p2    = precision(prec + 5);
    mp_float y      = mul(constants::mp_pi(p2)/constants::mp_ln2(p2)/2, k, prec);

    return mp_complex(x,y);
};
mp_complex  gmp_tester_prec::rand_scalar_c_exp10(precision prec, double max_re, double max_im)
{
    //rand complex value z such that im(ln(10) * z) = pi/2*k
    mp_float x  = rand_scalar(prec, max_re);
    mp_float kr = rand_scalar(precision(10), max(max_im,1.0));
    mp_float k;
    modf(kr, k);

    if (is_zero(k))
        k           = 1;

    precision p2    = precision(prec + 5);
    mp_float y      = mul(constants::mp_pi(p2)/constants::mp_ln10(p2)/2, k, prec);

    return mp_complex(x,y);
};

mp_complex  gmp_tester_prec::rand_scalar_c_log1p(precision prec, double max)
{
    mp_float x = rand_scalar(prec, max);
    mp_float y = rand_scalar(prec, max);
    return mp_complex(x,y);
};

mp_complex  gmp_tester_prec::rand_scalar_c_expm1(precision prec, double max)
{
    //numbers very close to zero
    mp_float x = rand_scalar(prec, 1e-100*max);
    mp_float y = rand_scalar(prec, 1e-100*max);
    return mp_complex(x,y);
};

Real gmp_tester_prec::calc_ulp(const mp_float& res, const mp_float& res_ext)
{
    mp_float dif    = abs(res - res_ext);
    mp_float tol    = eps(res);
    Real ulp        = (dif / tol).cast_float();

    return ulp;
};

template<class Derived>
struct eval_scalar_func_prec
{
    Integer code;

    eval_scalar_func_prec(Integer c) 
        : code(c)
    {};

    template<class T1>
    auto eval_fun(const T1& a1, precision p) -> decltype(Derived::eval(std::declval<T1>()))
    {
        return Derived::eval(a1, p);
    };

    template<class T1>
    auto eval_fun(const T1& a1) -> decltype(Derived::eval(std::declval<T1>()))
    {
        return Derived::eval(a1);
    };

    double make(const mp_float& s, precision prec)
    {
        if (code == -1)
            disp("break");

        precision prec_ext  = precision(prec + 50);

        mp_float res        = eval_fun(s, prec);
        mp_float res_ext    = eval_fun(mp_float(s,prec_ext), prec_ext);
        Real ulp            = gmp_tester_prec::calc_ulp(res, res_ext);
        return ulp;
    };
};
template<class Derived>
struct eval_bin_func_prec
{
    Integer code;

    eval_bin_func_prec(Integer c) 
        : code(c)
    {};

    template<class T1>
    auto eval_fun(const T1& a1, const T1& a2, precision p) 
        -> decltype(Derived::eval(a1, a2))
    {
        return Derived::eval(a1, a2, p);
    };

    template<class T1>
    auto eval_fun(const T1& a1, const T1& a2) -> decltype(Derived::eval(a1, a2))
    {
        return Derived::eval(a1, a2);
    };

    Real calc_ulp(const mp_float& res, const mp_float& res_ext)
    {
        mp_float dif    = abs(res - res_ext);
        mp_float tol    = eps(res);
        Real ulp        = (dif / tol).cast_float();

        return ulp;
    };

    double make(const mp_float& s1, const mp_float& s2, precision prec)
    {
        if (code == -1)
            disp("break");

        precision prec_ext  = precision(prec + 50);

        mp_float res        = eval_fun(s1, s2, prec);
        mp_float res_ext    = eval_fun(s1, mp_float(s2, prec_ext), prec_ext);
        Real ulp            = calc_ulp(res, res_ext);
        return ulp;
    };
};

template<class Derived>
struct has_real_return
{
    using ret_type = decltype(Derived::eval(std::declval<mp_complex>()));
    static const bool value = std::is_same<ret_type, mp_float>::value;
};
template<class Derived, class T1, class T2>
struct has_real_return2
{
    using ret_type = decltype(Derived::eval(T1(), T2()));
    static const bool value = std::is_same<ret_type, mp_float>::value;
};

template<class Derived>
struct eval_scalar_func_prec_c
{
    using dtup_2    = std::tuple<double,double>;
    static const bool real_return   = has_real_return<Derived>::value;

    Integer code;

    eval_scalar_func_prec_c(Integer c) 
        : code(c)
    {};

    template<class T1>
    auto eval_fun(const T1& a1) -> decltype(Derived::eval(std::declval<T1>()))
    {
        return Derived::eval(a1);
    };
    template<class T1>
    auto eval_fun(const T1& a1, precision p) -> decltype(Derived::eval(std::declval<T1>()))
    {
        return Derived::eval(a1, p);
    };

    dtup_2 calc_ulp(const mp_complex& res, const mp_complex& res_ext)
    {
        mp_complex dif      = res - res_ext;
        mp_float dif_re     = real(dif);
        mp_float dif_im     = imag(dif);
        mp_float dif_tot    = abs(dif);

        bool is_z_re        = is_zero(real(res));
        bool is_z_im        = is_zero(imag(res));

        mp_float tol_re     = is_z_re ? matcl::eps(0) : eps(real(res));
        mp_float tol_im     = is_z_im ? matcl::eps(0) : eps(imag(res));

        Real ulp_re         = (abs(dif_re) / tol_re).cast_float();
        Real ulp_im         = (abs(dif_im) / tol_im).cast_float();

        if (real_return == true)
            return dtup_2(ulp_re, 0.0);
        else
            return dtup_2(ulp_re, ulp_im);
    };

    dtup_2 make(const mp_complex& s, precision prec)
    {
        if (code == -1)
            disp("break");

        precision prec_ext  = precision(prec + 50);

        mp_complex res_ext  = eval_fun(mp_complex(s, prec_ext), prec_ext);
        mp_complex res      = eval_fun(s, prec);        
        dtup_2 ulp          = calc_ulp(res, res_ext);
        return ulp;
    };
};

template<class Derived, class T1, class T2>
struct eval_bin_func_prec_c
{
    using dtup_2    = std::tuple<double,double>;
    static const bool real_return   = has_real_return2<Derived,T1,T2>::value;

    Integer code;

    eval_bin_func_prec_c(Integer c) 
        : code(c)
    {};

    template<class T1, class T2>
    auto eval_fun(const T1& a1, const T2& a2) -> decltype(Derived::eval(a1, a2))
    {
        return Derived::eval(a1,a2);
    };
    template<class T1, class T2>
    auto eval_fun(const T1& a1, const T2& a2, precision p) -> decltype(Derived::eval(a1, a2))
    {
        return Derived::eval(a1, a2, p);
    };

    dtup_2 calc_ulp(const mp_complex& res, const mp_complex& res_ext)
    {
        mp_complex dif      = res - res_ext;
        mp_float dif_re     = real(dif);
        mp_float dif_im     = imag(dif);
        mp_float dif_tot    = abs(dif);

        bool is_z_re        = is_zero(real(res));
        bool is_z_im        = is_zero(imag(res));

        mp_float tol_re     = is_z_re ? matcl::eps(0) : eps(real(res));
        mp_float tol_im     = is_z_im ? matcl::eps(0) : eps(imag(res));

        Real ulp_re         = (abs(dif_re) / tol_re).cast_float();
        Real ulp_im         = (abs(dif_im) / tol_im).cast_float();

        if (real_return == true)
            return dtup_2(ulp_re, 0.0);
        else
            return dtup_2(ulp_re, ulp_im);
    };

    dtup_2 make(const T1& s10, const T2& s20, precision prec)
    {
        if (code == -1)
            disp("break");

        precision prec_ext  = precision(prec + 50);

        T1 s1(s10, prec_ext);
        T2 s2(s20, prec_ext);

        mp_complex res_ext  = eval_fun(s1, s2, prec_ext);
        mp_complex res      = eval_fun(s1, s2, prec);        
        dtup_2 ulp          = calc_ulp(res, res_ext);
        return ulp;
    };
};

template<class Func>
double gmp_tester_prec::test_scalar(const mp_float& s, precision p, Integer code)
{
    eval_scalar_func_prec<Func> test1(code);
    double res1 = test1.make(s, p);
    return res1;
};

template<class Func>
double gmp_tester_prec::test_constant(precision p)
{
    precision prec_ext  = precision(p + 20);

    mp_float res        = Func::eval(p);
    mp_float res_ext    = Func::eval(prec_ext);
    Real ulp            = calc_ulp(res, res_ext);
    return ulp;
};

template<class Func>
double gmp_tester_prec::test_bin(const mp_float& s1, const mp_float& s2, precision p, Integer code)
{
    eval_bin_func_prec<Func> test1(code);
    double res1 = test1.make(s1, s2, p);
    return res1;
};

template<class Func>
auto gmp_tester_prec::test_scalar_c(const mp_complex& s, precision p, Integer code) -> dtup_2
{
    eval_scalar_func_prec_c<Func> test1(code);
    dtup_2 res1 = test1.make(s, p);
    return res1;
};

template<class Func, class T1, class T2>
auto gmp_tester_prec::test_bin_c(const T1& s1, const T2& s2, precision p, Integer code) -> dtup_2
{
    eval_bin_func_prec_c<Func, T1, T2> test1(code);
    dtup_2 res1 = test1.make(s1, s2, p);
    return res1;
};

template<class Func>
bool gmp_tester_prec::test_scalar_func(std::ostream& os, Integer N, precision prec, double max,
                                       double accuracy)
{
    std::vector<mp_float> scalars;
    
    for (Integer i = 0; i < N; ++i)
        scalars.push_back(rand_scalar(prec, max));

    double res_mean = 0.0;
    double res_max  = 0.0;
    Integer N_res   = 0;

    for (Integer i = 0; i < N; ++i)
    {
        double res  = test_scalar<Func>(scalars[i], prec, i);

        if (matcl::is_finite(res) == false)
            continue;

        if (res > accuracy)
            res  = test_scalar<Func>(scalars[i], prec, i);

        res_max     = std::max(res_max, res);
        res_mean    += res;
        N_res       += 1;
    };

    N_res           = std::max(N_res, 1);
    res_mean        = res_mean / N_res;

    bool error      =  res_max > accuracy;

    if (error == false)
        return error;

    boost::io::ios_flags_saver old_flags(os);
    boost::io::ios_precision_saver old_prec(os);

    os  << std::setprecision(3);
    os  << std::setw(7)     << Func::name() 
        << ", max ulp: "    << std::fixed << res_max
        << ", mean ulp: "   << std::fixed << res_mean
        << "\n";

    return true;
};

template<class Func>
bool gmp_tester_prec::test_constant(std::ostream& os, precision prec, double accuracy)
{
    double res  = test_constant<Func>(prec);

    if (res > accuracy)
    {
        res     = test_constant<Func>(prec);
    }

    bool error  =  res > accuracy;

    if (error == false)
        return error;

    boost::io::ios_flags_saver old_flags(os);
    boost::io::ios_precision_saver old_prec(os);

    os  << std::setprecision(3);
    os  << std::setw(7)     << Func::name() 
        << ", ulp: "    << std::fixed << res
        << "\n";

    return true;
};

template<class Func>
bool gmp_tester_prec::test_bin_func(std::ostream& os, Integer N, precision prec, const max_vec& max,
                                    double accuracy)
{
    std::vector<mp_float> scalars1;
    std::vector<mp_float> scalars2;
    
    for (Integer i = 0; i < N; ++i)
    {
        scalars1.push_back(rand_scalar(prec, max));
        scalars2.push_back(rand_scalar(prec, max));
    };

    double res_mean = 0.0;
    double res_max  = 0.0;
    Integer N_res   = 0;

    for (Integer i = 0; i < N; ++i)
    {
        double res  = test_bin<Func>(scalars1[i], scalars2[i], prec, i);

        if (matcl::is_finite(res) == false)
            continue;

        res_max     = std::max(res_max, res);
        res_mean    += res;
        N_res       += 1;
    };

    N_res           = std::max(N_res, 1);
    res_mean        = res_mean / N_res;

    bool error      =  res_max > accuracy;

    if (error == false)
        return error;

    boost::io::ios_flags_saver old_flags(os);
    boost::io::ios_precision_saver old_prec(os);

    os  << std::setprecision(3);
    os  << std::setw(7)     << Func::name() 
        << ", max ulp: "    << std::fixed << res_max
        << ", mean ulp: "   << std::fixed << res_mean
        << "\n";

    return true;
};

template<class Func, class Rand>
bool gmp_tester_prec::test_bin_func_c(std::ostream& os, Integer N, precision p,
                            const prec_vec& prec, const max_vec& max, Real accuracy)
{
 
    double res_mean_re  = 0.0;
    double res_mean_im  = 0.0;
    double res_max_re   = 0.0;
    double res_max_im   = 0.0;

    Integer N_res       = 0;

    for (Integer i = 0; i < N; ++i)
    {
        mp_complex s1, s2;
        mp_float s3;
        Rand::rand(prec, max, s1, s2, s3);
 
        double res_re, res_im;
        std::tie(res_re, res_im) = test_bin_c<Func>(s1, s2, p, i);

        if (is_finite(res_re) == false || is_finite(res_im) == false)
            continue;

        if (res_re > accuracy || res_im > accuracy)
        {
            std::tie(res_re, res_im) = test_bin_c<Func>(s1, s2, p, i);
        };

        res_max_re      = std::max(res_max_re, res_re);
        res_max_im      = std::max(res_max_im, res_im);

        res_mean_re     += res_re;
        res_mean_im     += res_im;

        N_res           += 1;
    };
    for (Integer i = 0; i < N; ++i)
    {
        mp_complex s1, s2;
        mp_float s3;
        Rand::rand(prec, max, s1, s2, s3);

        double res_re, res_im;
        std::tie(res_re, res_im) = test_bin_c<Func>(s1, s3, p, i);

        if (is_finite(res_re) == false || is_finite(res_im) == false)
            continue;

        if (res_re > accuracy || res_im > accuracy)
        {
            std::tie(res_re, res_im) = test_bin_c<Func>(s1, s3, p, i);
        };

        res_max_re      = std::max(res_max_re, res_re);
        res_max_im      = std::max(res_max_im, res_im);

        res_mean_re     += res_re;
        res_mean_im     += res_im;

        N_res           += 1;
    };
    for (Integer i = 0; i < N; ++i)
    {
        mp_complex s1, s2;
        mp_float s3;
        Rand::rand(prec, max, s1, s2, s3);

        double res_re, res_im;
        std::tie(res_re, res_im) = test_bin_c<Func>(s3, s2, p, i);

        if (is_finite(res_re) == false || is_finite(res_im) == false)
            continue;

        if (res_re > accuracy || res_im > accuracy)
        {
            double res_re1, res_im1;
            std::tie(res_re1, res_im1) = test_bin_c<Func>(s3, s2, p, i);
        };

        res_max_re      = std::max(res_max_re, res_re);
        res_max_im      = std::max(res_max_im, res_im);

        res_mean_re     += res_re;
        res_mean_im     += res_im;

        N_res           += 1;
    };

    N_res               = std::max(N_res, 1);

    res_mean_re         = res_mean_re / N_res;
    res_mean_im         = res_mean_im / N_res;

    bool is_add         = std::is_same<Func,Plus_func>::value || std::is_same<Func,Minus_func>::value;
    bool error          =  res_max_re > accuracy || res_max_im > accuracy;

    (void)is_add;

    if (error == false)
        return error;

    boost::io::ios_flags_saver old_flags(os);
    boost::io::ios_precision_saver old_prec(os);

    os  << std::setprecision(3);
    os  << std::setw(7) << Func::name() 
        << ", max ulp (re,im): "
            << std::fixed << res_max_re << " , "
            << std::fixed << res_max_im << " , "
        << ", mean ulp (re,im): "
            << std::fixed << res_mean_re << " , "
            << std::fixed << res_mean_im << " , "
        << "\n";

    return true;
};

template<class Func, class Rand>
bool gmp_tester_prec::test_scalar_func_c(std::ostream& os, Integer N, precision prec, 
                                         double max_re, double max_im, double accuracy)
{
    std::vector<mp_complex> scalars;
    
    for (Integer i = 0; i < N; ++i)
        scalars.push_back(Rand::rand_scalar_c(prec, max_re, max_im));

    double res_mean_re  = 0.0;
    double res_mean_im  = 0.0;
    double res_max_re   = 0.0;
    double res_max_im   = 0.0;

    Integer N_res       = 0;

    for (Integer i = 0; i < N; ++i)
    {
        double res_re, res_im;
        std::tie(res_re, res_im) = test_scalar_c<Func>(scalars[i], prec, i);

        if (is_finite(res_re) == false || is_finite(res_im) == false)
            continue;

        if (res_re > accuracy || res_im > accuracy)
        {
            std::tie(res_re, res_im) = test_scalar_c<Func>(scalars[i], prec, i);
        };
        res_max_re      = std::max(res_max_re, res_re);
        res_max_im      = std::max(res_max_im, res_im);

        res_mean_re     += res_re;
        res_mean_im     += res_im;

        N_res           += 1;
    };

    N_res               = std::max(N_res, 1);

    res_mean_re         = res_mean_re / N_res;
    res_mean_im         = res_mean_im / N_res;

    bool error          =  res_max_re > accuracy 
                        || res_max_im > accuracy;

    if (error == false)
        return error;

    boost::io::ios_flags_saver old_flags(os);
    boost::io::ios_precision_saver old_prec(os);

    os  << std::setprecision(3);
    os  << std::setw(7) << Func::name() 
        << ", max ulp (re,im): "
            << std::fixed << res_max_re << " , "
            << std::fixed << res_max_im << " , "
        << ", mean ulp (re,im): "
            << std::fixed << res_mean_re << " , "
            << std::fixed << res_mean_im << " , "
        << "\n";

    return true;
};

}};

#pragma warning(pop)