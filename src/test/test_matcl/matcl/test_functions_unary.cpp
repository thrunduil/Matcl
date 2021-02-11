/*
 *  This file is a part of Matrix Computation Library (MATCL)
 *
 *  Copyright (c) Pawe³ Kowal 2018 - 2021
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

#include "test_set.h"
#include "test_functions_unary.h"

#include "test/test_matcl/framework/matrix_set/matrix_set_1.h"
#include "matcl-core/IO/logger.h"
#include "test/test_matcl/framework/data_struct/test_selector.h"

#include <boost/thread.hpp>

namespace gd = matcl::details;

namespace matcl { namespace test
{

class test_unary
{
    unary_functions_list&		tf;
    const test::options&		opts;
    Integer                     thread_id;

    public:
        test_unary(unary_functions_list& tf, const test::options& opts, Integer id)
            :tf(tf),opts(opts), thread_id(id)
        {};

        test_unary(const test_unary& tu)
            :tf(tu.tf),opts(tu.opts), thread_id(tu.thread_id)
        {};

        void make()
        {   
            /*
            Integer code = 130;
            Matrix mat = tf.get_matrix(code);
            matcl::disp(mat);
            */

            tf.make(opts);
        };

        void operator()()
        {
            make();
        };

    private:		
        test_unary& operator=(const test_unary&) = delete;
};

void test_unary_st(const rand_matrix_ptr& rand)
{
    test::options opts;

    try
    {
        opts.show_partial_res   = 0;
        opts.first_matrix_code  = 0;
        opts.show_memleaks      = true;

        test::mat_set_1 ms1(rand);
        unary_functions_list tf(ms1,rand);
        
        test_unary tu(tf,opts,0);
        tu.make();
    }
    catch(const std::exception& ex)
    {
        matcl::out_stream <<ex.what();
    };
};

void test_unary_mt(const rand_matrix_ptr& rand)
{
    test::options opts;
    try
    {
        opts.show_partial_res   = 0;
        opts.first_matrix_code  = 0;
        opts.show_memleaks      = false;

        test::mat_set_2 ms1(rand);
        unary_functions_list tf(ms1,rand);

        boost::thread_group tg;

        for (int i = 0; i < 10; i++)
            tg.create_thread(test_unary(tf,opts,i));

        tg.join_all();
    }
    catch(const std::exception& ex)
    {
        matcl::out_stream <<ex.what();
    };
};

static const Integer max_int = 1000;

unary_functions_list::unary_functions_list(const matrix_set& ms,rand_matrix_ptr rand)
:m_tests(ms), m_rand(rand)
{};

void unary_functions_list::make(options opts)
{
    m_options = opts;

    SELECT_TEST (3, test_ldexp());
    SELECT_TEST (3, test_scalbn());
    SELECT_TEST (3, test_sqrt1pm1());
    SELECT_TEST (3, test_cbrt());
    SELECT_TEST (3, test_log1p());
    SELECT_TEST (3, test_sqrt1pm1_c());
    SELECT_TEST (3, test_log1p_c());
    SELECT_TEST (3, test_op_plus());
    SELECT_TEST (3, test_expm1());
    SELECT_TEST (3, test_signbit());
    SELECT_TEST (3, test_logb());
    SELECT_TEST (3, test_ilogb());
    SELECT_TEST (3, test_exp2());
    SELECT_TEST (3, test_exp10());
    SELECT_TEST (3, test_expi());
    SELECT_TEST (3, test_frexp());
    SELECT_TEST (3, test_modf());
    SELECT_TEST (3, test_fpclassify());
    SELECT_TEST (3, test_nextbelow());
    SELECT_TEST (3, test_nextabove());
    SELECT_TEST (3, test_eps());
    SELECT_TEST (3, test_bernoulli_b2n());
    SELECT_TEST (3, test_max_bernoulli_b2n());
    SELECT_TEST (3, test_prime());
    SELECT_TEST (3, test_prime_max_count());
    SELECT_TEST (3, test_factorial());
    SELECT_TEST (3, test_double_factorial());
    SELECT_TEST (3, test_rising_factorial());
    SELECT_TEST (3, test_binomial_coefficient());
    SELECT_TEST (3, test_falling_factorial());        
    SELECT_TEST (3, test_eval_scalar_func());

    SELECT_TEST (3, test_real());
    SELECT_TEST (3, test_imag());
    SELECT_TEST (3, test_abs());
    SELECT_TEST (3, test_abs2());
    SELECT_TEST (3, test_arg());
    SELECT_TEST (3, test_angle());
    SELECT_TEST (3, test_conj());

    SELECT_TEST (3, test_sqrt());
    SELECT_TEST (3, test_pow2());
    SELECT_TEST (3, test_exp());
    SELECT_TEST (3, test_log());
    SELECT_TEST (3, test_log2());
    SELECT_TEST (3, test_log10());
    SELECT_TEST (3, test_floor());
    SELECT_TEST (3, test_ceil());
    SELECT_TEST (3, test_round());
    SELECT_TEST (3, test_fix());
    SELECT_TEST (3, test_trunc());
    SELECT_TEST (3, test_sign());
    SELECT_TEST (3, test_isign());
    SELECT_TEST (3, test_sin());
    SELECT_TEST (3, test_cos());
    SELECT_TEST (3, test_tan());
    SELECT_TEST (3, test_cot());
    SELECT_TEST (3, test_sec());
    SELECT_TEST (3, test_csc());
    SELECT_TEST (3, test_asin());
    SELECT_TEST (3, test_acos());
    SELECT_TEST (3, test_atan());
    SELECT_TEST (3, test_acot());
    SELECT_TEST (3, test_asec());
    SELECT_TEST (3, test_acsc());
    SELECT_TEST (3, test_sinh());
    SELECT_TEST (3, test_cosh());
    SELECT_TEST (3, test_tanh());
    SELECT_TEST (3, test_coth());
    SELECT_TEST (3, test_sech());
    SELECT_TEST (3, test_csch());
    SELECT_TEST (3, test_asinh());
    SELECT_TEST (3, test_acosh());
    SELECT_TEST (3, test_atanh());
    SELECT_TEST (3, test_acoth());
    SELECT_TEST (3, test_asech());
    SELECT_TEST (3, test_acsch());

    SELECT_TEST (3, test_sqrt_c());
    SELECT_TEST (3, test_log_c());
    SELECT_TEST (3, test_log2_c());
    SELECT_TEST (3, test_log10_c());
    SELECT_TEST (3, test_asin_c());
    SELECT_TEST (3, test_acos_c());
    SELECT_TEST (3, test_asec_c());
    SELECT_TEST (3, test_acsc_c());
    SELECT_TEST (3, test_acosh_c());
    SELECT_TEST (3, test_atanh_c());
    SELECT_TEST (3, test_acoth_c());
    SELECT_TEST (3, test_asech_c());

    SELECT_TEST (3, test_is_true());
    SELECT_TEST (3, test_is_false());
    SELECT_TEST (3, test_ifloor());
    SELECT_TEST (3, test_iceil());
    SELECT_TEST (3, test_iround());
    SELECT_TEST (3, test_ifix());
    SELECT_TEST (3, test_itrunc());

    SELECT_TEST (3, test_is_inf());
    SELECT_TEST (3, test_is_nan());
    SELECT_TEST (3, test_is_finite());
    SELECT_TEST (3, test_is_regular());
    SELECT_TEST (3, test_is_normal());
    SELECT_TEST (3, test_is_int());
    SELECT_TEST (3, test_is_real());

    SELECT_TEST (3, test_is_scalar_true());
    SELECT_TEST (3, test_is_scalar_false());	

    SELECT_TEST (3, test_op_un_minus());
    SELECT_TEST (3, test_op_neg());
    SELECT_TEST (3, test_inv());
    SELECT_TEST (3, test_invs());
};

Matrix unary_functions_list::get_matrix(int code) const
{
    return m_tests.get_matrix(code);
};

//
void unary_functions_list::test_is_inf()
{
    return test_function<function_isinf>();
};
void unary_functions_list::test_is_nan()
{
    return test_function<function_isnan>();
};
void unary_functions_list::test_is_finite()
{
    return test_function<function_isfin>();
};
void unary_functions_list::test_is_regular()
{
    return test_function<function_isregul>();
};
void unary_functions_list::test_is_normal()
{
    return test_function<function_isnorm>();
};
void unary_functions_list::test_is_int()
{
    return test_function<function_isint>();
};
void unary_functions_list::test_is_real()
{
    return test_function<function_isreal>();
};

void unary_functions_list::test_is_scalar_true()
{
    return test_function<function_isstrue>();
};
void unary_functions_list::test_is_scalar_false()
{
    return test_function<function_issfalse>();
};	

void unary_functions_list::test_sqrt1pm1()
{
    return test_function<function_sqrt1pm1>();
};;
void unary_functions_list::test_cbrt()
{
    return test_function<function_cbrt>();
};;
void unary_functions_list::test_log1p()
{
    return test_function<function_log1p>();
};
void unary_functions_list::test_sqrt1pm1_c()
{
    return test_function<function_sqrt1pm1_c>();
};
void unary_functions_list::test_log1p_c()
{
    return test_function<function_log1p_c>();
};
void unary_functions_list::test_op_plus()
{
    return test_function<function_op_plus>();
};
void unary_functions_list::test_expm1()
{
    return test_function<function_expm1>();
};

void unary_functions_list::test_logb()
{
    return test_function<function_logb>();
};
void unary_functions_list::test_exp2()
{
    return test_function<function_exp2>();
};
void unary_functions_list::test_exp10()
{
    return test_function<function_exp10>();
};
void unary_functions_list::test_expi()
{
    return test_function<function_expi>();
};

void unary_functions_list::test_real()
{
    return test_function<function_real>();
};;
void unary_functions_list::test_imag()
{
    return test_function<function_imag>();
};;			
void unary_functions_list::test_abs()
{
    return test_function<function_abs>();
};

void unary_functions_list::test_arg()
{
    return test_function<function_arg>();
};;
void unary_functions_list::test_angle()
{
    return test_function<function_angle>();
};;
void unary_functions_list::test_conj()
{
    return test_function<function_conj>();
};;			
void unary_functions_list::test_sqrt()
{
    return test_function<function_sqrt>();
};;
void unary_functions_list::test_pow2()
{
    return test_function<function_pow2>();
};;			
void unary_functions_list::test_exp()
{
    return test_function<function_exp>();
};;
void unary_functions_list::test_log()
{
    return test_function<function_log>();
};;
void unary_functions_list::test_log2()
{
    return test_function<function_log2>();
};;
void unary_functions_list::test_log10()
{
    return test_function<function_log10>();
};;			
void unary_functions_list::test_floor()
{
    return test_function<function_floor>();
};;
void unary_functions_list::test_ceil()
{
    return test_function<function_ceil>();
};;			
void unary_functions_list::test_round()
{
    return test_function<function_round>();
};;
void unary_functions_list::test_fix()
{
    return test_function<function_fix>();
};;
void unary_functions_list::test_trunc()
{
    return test_function<function_trunc>();
};;
void unary_functions_list::test_ifloor()
{
    return test_function<function_ifloor>();
};;
void unary_functions_list::test_iceil()
{
    return test_function<function_iceil>();
};;
void unary_functions_list::test_iround()
{
    return test_function<function_iround>();
};;
void unary_functions_list::test_ifix()
{
    return test_function<function_ifix>();
};;
void unary_functions_list::test_itrunc()
{
    return test_function<function_itrunc>();
};;

void unary_functions_list::test_sign()
{
    return test_function<function_sign>();
};;
void unary_functions_list::test_isign()
{
    return test_function<function_isign>();
};;			

void unary_functions_list::test_sin()
{
    return test_function<function_sin>();
};;
void unary_functions_list::test_cos()
{
    return test_function<function_cos>();
};;
void unary_functions_list::test_tan()
{
    return test_function<function_tan>();
};;
void unary_functions_list::test_cot()
{
    return test_function<function_cot>();
};;
void unary_functions_list::test_sec()
{
    return test_function<function_sec>();
};;
void unary_functions_list::test_csc()
{
    return test_function<function_csc>();
};;
void unary_functions_list::test_asin()
{
    return test_function<function_asin>();
};;
void unary_functions_list::test_acos()
{
    return test_function<function_acos>();
};;			
void unary_functions_list::test_atan()
{
    return test_function<function_atan>();
};;
void unary_functions_list::test_acot()
{
    return test_function<function_acot>();
};;			
void unary_functions_list::test_asec()
{
    return test_function<function_asec>();
};;
void unary_functions_list::test_acsc()
{
    return test_function<function_acsc>();
};;			
void unary_functions_list::test_sinh()
{
    return test_function<function_sinh>();
};;
void unary_functions_list::test_cosh()
{
    return test_function<function_cosh>();
};;			
void unary_functions_list::test_tanh()
{
    return test_function<function_tanh>();
};;
void unary_functions_list::test_coth()
{
    return test_function<function_coth>();
};;			
void unary_functions_list::test_sech()
{
    return test_function<function_sech>();
};;
void unary_functions_list::test_csch()
{
    return test_function<function_csch>();
};;			
void unary_functions_list::test_asinh()
{
    return test_function<function_asinh>();
};;
void unary_functions_list::test_acosh()
{
    return test_function<function_acosh>();
};			
void unary_functions_list::test_atanh()
{
    return test_function<function_atanh>();
};
void unary_functions_list::test_acoth()
{
    return test_function<function_acoth>();
};			
void unary_functions_list::test_asech()
{
    return test_function<function_asech>();
};
void unary_functions_list::test_acsch()
{
    return test_function<function_acsch>();
};		
void unary_functions_list::test_sqrt_c()
{
    return test_function<function_sqrt_c>();
};		
void unary_functions_list::test_log_c()
{
    return test_function<function_log_c>();
};		
void unary_functions_list::test_log2_c()
{
    return test_function<function_log2_c>();
};		
void unary_functions_list::test_log10_c()
{
    return test_function<function_log10_c>();
};		
void unary_functions_list::test_asin_c()
{
    return test_function<function_asin_c>();
};		
void unary_functions_list::test_acos_c()
{
    return test_function<function_acos_c>();
};		
void unary_functions_list::test_asec_c()
{
    return test_function<function_asec_c>();
};		
void unary_functions_list::test_acsc_c()
{
    return test_function<function_acsc_c>();
};		
void unary_functions_list::test_acosh_c()
{
    return test_function<function_acosh_c>();
};		
void unary_functions_list::test_atanh_c()
{
    return test_function<function_atanh_c>();
};		
void unary_functions_list::test_acoth_c()
{
    return test_function<function_acoth_c>();
};		
void unary_functions_list::test_asech_c()
{
    return test_function<function_asech_c>();
};		

void unary_functions_list::test_abs2()
{
    Real out = 0.;
    test_function_abs2 tf;
    out += m_tests.make(&tf,m_options);

    if (out == 0.)
        matcl::out_stream << std::string() +   "abs2: OK" + "\n";
    else
        matcl::out_stream << std::string() +   "abs2: FAILED"  + "\n";
};

void unary_functions_list::test_is_true()
{
    Real out = 0.;
    test_function_is_true tf;
    out += m_tests.make(&tf,m_options);

    if (out == 0.)
        matcl::out_stream << std::string() +   "is_true: OK" + "\n";
    else
        matcl::out_stream << std::string() +   "is_true: FAILED"  + "\n";
};

void unary_functions_list::test_signbit()
{
    Real out = 0.;
    test_function_signbit tf;
    out += m_tests.make(&tf,m_options);

    if (out == 0.)
        matcl::out_stream << std::string() +   "signbit: OK" + "\n";
    else
        matcl::out_stream << std::string() +   "signbit: FAILED"  + "\n";
};

void unary_functions_list::test_ilogb()
{
    Real out = 0.;
    test_function_ilogb tf;
    out += m_tests.make(&tf,m_options);

    if (out == 0.)
        matcl::out_stream << std::string() +   "ilogb: OK" + "\n";
    else
        matcl::out_stream << std::string() +   "ilogb: FAILED"  + "\n";
};

void unary_functions_list::test_fpclassify()
{
    Real out = 0.;
    test_function_fpclassify tf;
    out += m_tests.make(&tf,m_options);

    if (out == 0.)
        matcl::out_stream << std::string() +   "fpclassify: OK" + "\n";
    else
        matcl::out_stream << std::string() +   "fpclassify: FAILED"  + "\n";
};

void unary_functions_list::test_nextabove()
{
    Real out = 0.;
    test_function_nextabove tf;
    out += m_tests.make(&tf,m_options);

    if (out == 0.)
        matcl::out_stream << std::string() +   "nextabove: OK" + "\n";
    else
        matcl::out_stream << std::string() +   "nextabove: FAILED"  + "\n";
};

void unary_functions_list::test_nextbelow()
{
    Real out = 0.;
    test_function_nextbelow tf;
    out += m_tests.make(&tf,m_options);

    if (out == 0.)
        matcl::out_stream << std::string() +   "nextbelow: OK" + "\n";
    else
        matcl::out_stream << std::string() +   "nextbelow: FAILED"  + "\n";
};

void unary_functions_list::test_frexp()
{
    Real out = 0.;
    test_function_frexp tf;
    out += m_tests.make(&tf,m_options);

    if (out == 0.)
        matcl::out_stream << std::string() +   "frexp: OK" + "\n";
    else
        matcl::out_stream << std::string() +   "frexp: FAILED"  + "\n";
};

void unary_functions_list::test_modf()
{
    Real out = 0.;
    test_function_modf tf;
    out += m_tests.make(&tf,m_options);

    if (out == 0.)
        matcl::out_stream << std::string() +   "modf: OK" + "\n";
    else
        matcl::out_stream << std::string() +   "modf: FAILED"  + "\n";
};

void unary_functions_list::test_eps()
{
    Real out = 0.;
    test_function_eps tf;
    out += m_tests.make(&tf,m_options);

    if (out == 0.)
        matcl::out_stream << std::string() +   "eps: OK" + "\n";
    else
        matcl::out_stream << std::string() +   "eps: FAILED"  + "\n";
};

void unary_functions_list::test_ldexp()
{
    Real out = 0.;
    test_function_ldexp tf;
    out += m_tests.make(&tf,m_options);

    if (out == 0.)
        matcl::out_stream << std::string() +   "ldexp: OK" + "\n";
    else
        matcl::out_stream << std::string() +   "ldexp: FAILED"  + "\n";
};
void unary_functions_list::test_scalbn()
{
    Real out = 0.;
    test_function_scalbn tf;
    out += m_tests.make(&tf,m_options);

    if (out == 0.)
        matcl::out_stream << std::string() +   "scalbn: OK" + "\n";
    else
        matcl::out_stream << std::string() +   "scalbn: FAILED"  + "\n";
};

void unary_functions_list::test_bernoulli_b2n()
{
    Real out = 0.;
    test_function_bernoulli_b2n tf;
    out += m_tests.make(&tf,m_options);

    if (out == 0.)
        matcl::out_stream << std::string() +   "bernoulli_b2n: OK" + "\n";
    else
        matcl::out_stream << std::string() +   "bernoulli_b2n: FAILED"  + "\n";
};
void unary_functions_list::test_max_bernoulli_b2n()
{
    Real out = 0.;
    test_function_max_bernoulli_b2n tf;
    out += m_tests.make(&tf,m_options);

    if (out == 0.)
        matcl::out_stream << std::string() +   "max_bernoulli_b2n: OK" + "\n";
    else
        matcl::out_stream << std::string() +   "max_bernoulli_b2n: FAILED"  + "\n";
};
void unary_functions_list::test_prime()
{
    Real out = 0.;
    test_function_prime tf;
    out += m_tests.make(&tf,m_options);

    if (out == 0.)
        matcl::out_stream << std::string() +   "prime: OK" + "\n";
    else
        matcl::out_stream << std::string() +   "prime: FAILED"  + "\n";
};
void unary_functions_list::test_prime_max_count()
{
    Real out = 0.;
    test_function_prime_max_count tf;
    out += m_tests.make(&tf,m_options);

    if (out == 0.)
        matcl::out_stream << std::string() +   "prime_max_count: OK" + "\n";
    else
        matcl::out_stream << std::string() +   "prime_max_count: FAILED"  + "\n";
};
void unary_functions_list::test_factorial()
{
    Real out = 0.;
    test_function_factorial tf;
    out += m_tests.make(&tf,m_options);

    if (out == 0.)
        matcl::out_stream << std::string() +   "factorial: OK" + "\n";
    else
        matcl::out_stream << std::string() +   "factorial: FAILED"  + "\n";
};
void unary_functions_list::test_double_factorial()
{
    Real out = 0.;
    test_function_double_factorial tf;
    out += m_tests.make(&tf,m_options);

    if (out == 0.)
        matcl::out_stream << std::string() +   "double_factorial: OK" + "\n";
    else
        matcl::out_stream << std::string() +   "double_factorial: FAILED"  + "\n";
};
void unary_functions_list::test_rising_factorial()
{
    Real out = 0.;
    test_function_rising_factorial tf;
    out += m_tests.make(&tf,m_options);

    if (out == 0.)
        matcl::out_stream << std::string() +   "rising_factorial: OK" + "\n";
    else
        matcl::out_stream << std::string() +   "rising_factorial: FAILED"  + "\n";
};
void unary_functions_list::test_binomial_coefficient()
{
    Real out = 0.;
    test_function_binomial_coefficient tf;
    out += m_tests.make(&tf,m_options);

    if (out == 0.)
        matcl::out_stream << std::string() +   "binomial_coefficient: OK" + "\n";
    else
        matcl::out_stream << std::string() +   "binomial_coefficient: FAILED"  + "\n";
};
void unary_functions_list::test_falling_factorial()
{
    Real out = 0.;
    test_function_falling_factorial tf;
    out += m_tests.make(&tf,m_options);

    if (out == 0.)
        matcl::out_stream << std::string() +   "falling_factorial: OK" + "\n";
    else
        matcl::out_stream << std::string() +   "falling_factorial: FAILED"  + "\n";
};
void unary_functions_list::test_eval_scalar_func()
{
    Real out = 0.;
    test_function_eval_scalar_func tf;
    out += m_tests.make(&tf,m_options);

    if (out == 0.)
        matcl::out_stream << std::string() +   "eval_scalar_func: OK" + "\n";
    else
        matcl::out_stream << std::string() +   "eval_scalar_func: FAILED"  + "\n";
};

void unary_functions_list::test_is_false()
{
    Real out = 0.;
    test_function_is_false tf;
    out += m_tests.make(&tf,m_options);

    if (out == 0.)
        matcl::out_stream << std::string() +   "is_false: OK" + "\n";
    else
        matcl::out_stream << std::string() +   "is_false: FAILED"  + "\n";
};

template<class Func>
void unary_functions_list::test_function()
{
    Real out = 0.;
    Func func;
    out = m_tests.make(func.function(),m_options);

    if (out == 0.)
        matcl::out_stream << std::string() +   func.name() + ": OK" + "\n";
    else
        matcl::out_stream << std::string() +   func.name() + ": FAILED"  + "\n";
};

void unary_functions_list::test_op_un_minus()
{
    Real out = 0.;
    test_function_op_un_minus tf;
    out += m_tests.make(&tf,m_options);

    if (out == 0.)
        matcl::out_stream << std::string() +   "operator-: OK" + "\n";
    else
        matcl::out_stream << std::string() +   "operator-: FAILED"  + "\n";
};		

void unary_functions_list::test_invs()
{
    Real out = 0.;
    test_function_invs tf;
    out += m_tests.make(&tf,m_options);

    if (out == 0.)
        matcl::out_stream << std::string() +   "invs: OK" + "\n";
    else
        matcl::out_stream << std::string() +   "invs: FAILED"  + "\n";
};		

void unary_functions_list::test_inv()
{
    Real out = 0.;
    test_function_inv tf;
    out += m_tests.make(&tf,m_options);

    if (out == 0.)
        matcl::out_stream << std::string() +   "inv: OK" + "\n";
    else
        matcl::out_stream << std::string() +   "inv: FAILED"  + "\n";
};		

void unary_functions_list::test_op_neg()
{
    Real out = 0.;
    test_function_op_neg tf;
    out += m_tests.make(&tf,m_options);

    if (out == 0.)
        matcl::out_stream << std::string() +   "operator~: OK" + "\n";
    else
        matcl::out_stream << std::string() +   "operator~: FAILED"  + "\n";
};			

//
Real test_function_op_un_minus::eval_mat(const Matrix& mat,bool,int code )
{
    (void)code;
    Matrix out_full     = -full(mat);	
    try
    {		
        Matrix out		= -mat;
        check_struct(out);

        Matrix out2     = uminus(mat);

        Real dif		= norm_1(out - out_full);
        dif				= dif + norm_1(-out - mat);
        dif				+= norm_1(out - out2);

        return dif;
    }
    catch(const std::exception& ex)
    {
        m_is_error  = true;
        m_error     = ex.what();
        return 1.;
    }
    catch(...)
    {
        m_is_error  = true;
        m_error     = "unknown error";
        return 1.;
    };
};
Real test_function_op_un_minus::eval_scalar(const Scalar& s,bool,int code  )
{
    (void)code;
    return eval_scalar_impl<test_function_op_un_minus>(s);
};

Real test_function_inv::eval_mat(const Matrix& mat,bool,int code )
{
    (void)mat;
    (void)code;
    return 0.0;
};
Real test_function_inv::eval_scalar(const Scalar& s,bool,int code  )
{
    (void)code;
    return eval_scalar_impl<test_function_inv>(s);
};

Real test_function_invs::eval_mat(const Matrix& mat,bool,int code )
{
    (void)code;
    Matrix out_full     = invs(full(mat));	
    try
    {		
        Matrix out		= invs(mat);
        check_struct(out);

        Matrix out2     = div(1,mat);

        Real dif		= norm_1(out - out_full);
        dif				+= norm_1(out - out2);

        if (abs(dif) < norm_1(out) * constants::eps(out.get_value_code()) * 10.0)
            dif = 0.0;

        return dif;
    }
    catch(const std::exception& ex)
    {
        m_is_error  = true;
        m_error     = ex.what();
        return 1.;
    }
    catch(...)
    {
        m_is_error  = true;
        m_error     = "unknown error";
        return 1.;
    };
};
Real test_function_invs::eval_scalar(const Scalar& s,bool,int code  )
{
    (void)code;
    return eval_scalar_impl<test_function_invs>(s);
};

Real test_function_op_neg::eval_mat(const Matrix& mat,bool,int code )
{
    (void)code;
    Matrix out_full = ~full(mat);	
    try
    {		
        Matrix out		= ~mat;
        Matrix out2		= neg(mat);
        check_struct(out);

        Real dif		= norm_1(out - out_full);
        dif             +=norm_1(out - out2);
        dif				= dif + norm_1(~out - is_true(mat));

        return dif;
    }
    catch(const std::exception& ex)
    {
        m_is_error  = true;
        m_error     = ex.what();
        return 1.;
    }
    catch(...)
    {
        m_is_error  = true;
        m_error     = "unknown error";
        return 1.;
    };
};
Real test_function_op_neg::eval_scalar(const Scalar& s,bool,int code  )
{
    (void)code;
    return eval_scalar_impl<test_function_op_neg>(s);
};

Real test_function_signbit::eval_mat(const Matrix& mat,bool,int code )
{
    (void)code;

    try
    {		
        bool is_complex = matrix_traits::is_float_complex(mat.get_value_code());

        (void)is_complex;

        Matrix out_full = signbit(full(mat));
        Matrix out		= signbit(mat);

        check_struct(out);

        Real dif        = norm_1(out - out_full)/(norm_1(out) + constants::eps());

        return dif;
    }
    catch(const error::function_not_defined_for_complex&)
    {
        if (any_vec(imag(mat)) == true)
            return 0.0;
        else
            return 1.0;
    }
    catch(const std::exception& ex)
    {
        m_is_error  = true;
        m_error     = ex.what();
        return 1.;
    }
    catch(...)
    {
        m_is_error  = true;
        m_error     = "unknown error";
        return 1.;
    };
};
Real test_function_signbit::eval_scalar(const Scalar& s,bool,int code  )
{
    (void)code;
    return eval_scalar_impl<test_function_signbit>(s);
};

Real test_function_ilogb::eval_mat(const Matrix& mat,bool,int code )
{
    (void)code;

    try
    {		
        Matrix out_full = ilogb(full(mat));
        Matrix out		= ilogb(mat);

        check_struct(out);

        Real dif        = norm_1(out - out_full)/(norm_1(out) + constants::eps());

        return dif;
    }
    catch(const error::function_not_defined_for_complex&)
    {
        if (any_vec(imag(mat)) == true)
            return 0.0;
        else
            return 1.0;
    }
    catch(const std::exception& ex)
    {
        m_is_error  = true;
        m_error     = ex.what();
        return 1.;
    }
    catch(...)
    {
        m_is_error  = true;
        m_error     = "unknown error";
        return 1.;
    };
};
Real test_function_ilogb::eval_scalar(const Scalar& s,bool,int code  )
{
    (void)code;
    return eval_scalar_impl<test_function_ilogb>(s);
};

Real test_function_fpclassify::eval_mat(const Matrix& mat,bool,int code )
{
    (void)mat;
    (void)code;
    return 0.0;
};

Real test_function_fpclassify::eval_scalar(const Scalar& s,bool,int code  )
{
    (void)code;
    return eval_scalar_impl<test_function_fpclassify>(s);
};

Real test_function_ldexp::eval_mat(const Matrix& mat,bool,int code )
{
    (void)mat;
    (void)code;
    return 0.0;
};
Real test_function_ldexp::eval_scalar(const Scalar& s,bool,int code  )
{
    (void)code;
    return eval_scalar_impl<test_function_ldexp>(s);
};

Real test_function_nextabove::eval_mat(const Matrix& mat,bool,int code )
{
    (void)mat;
    (void)code;
    return 0.0;
};
Real test_function_nextabove::eval_scalar(const Scalar& s,bool,int code  )
{
    (void)code;
    return eval_scalar_impl<test_function_nextabove>(s);
};

Real test_function_eps::eval_mat(const Matrix& mat,bool,int code )
{
    (void)mat;
    (void)code;
    return 0.0;
};
Real test_function_eps::eval_scalar(const Scalar& s,bool,int code  )
{
    (void)code;
    return eval_scalar_impl<test_function_eps>(s);
};

Real test_function_modf::eval_mat(const Matrix& mat,bool,int code )
{
    (void)mat;
    (void)code;
    return 0.0;
};
Real test_function_modf::eval_scalar(const Scalar& s,bool,int code  )
{
    (void)code;
    return eval_scalar_impl<test_function_modf>(s);
};

Real test_function_frexp::eval_mat(const Matrix& mat,bool,int code )
{
    (void)code;
    (void)mat;
    return 0.0;
};
Real test_function_frexp::eval_scalar(const Scalar& s,bool,int code  )
{
    (void)code;
    return eval_scalar_impl<test_function_frexp>(s);
};

Real test_function_nextbelow::eval_mat(const Matrix& mat,bool,int code )
{
    (void)code;
    (void)mat;
    return 0.0;
};
Real test_function_nextbelow::eval_scalar(const Scalar& s,bool,int code  )
{
    (void)code;
    return eval_scalar_impl<test_function_nextbelow>(s);
};

Real test_function_scalbn::eval_mat(const Matrix& mat,bool,int code )
{
    (void)mat;
    (void)code;
    return 0.0;
};

Real test_function_scalbn::eval_scalar(const Scalar& s,bool,int code  )
{
    (void)code;
    return eval_scalar_impl<test_function_scalbn>(s);
};

Real test_function_bernoulli_b2n::eval_mat(const Matrix& mat,bool,int code )
{
    (void)mat;
    (void)code;
    return 0.0;
};
Real test_function_bernoulli_b2n::eval_scalar(const Scalar& s,bool,int code  )
{
    (void)s;
    (void)code;
    Real v1     = bernoulli_b2n(0);
    Real v2     = bernoulli_b2n(1);

    Float f1    = fbernoulli_b2n(0);
    Float f2    = fbernoulli_b2n(1);

    Real w1     = bernoulli_b2n<Real>(0);
    Real w2     = bernoulli_b2n<Real>(1);
    Real w3     = bernoulli_b2n<Float>(0);
    Real w4     = bernoulli_b2n<Float>(1);

    bernoulli_b2n<Integer>(0);
    bernoulli_b2n<Complex>(0);
    bernoulli_b2n<Float_complex>(0);

    Real out    = norm_1(v1 - f1);
    out         += norm_1(v2 - f2);
    out         += norm_1(v1 - w1);
    out         += norm_1(v2 - w2);
    out         += norm_1(f1 - w3);
    out         += norm_1(f2 - w4);

    if (out < constants::f_eps() * 10.0)
        out     = 0.0;

    return out;
};

Real test_function_max_bernoulli_b2n::eval_mat(const Matrix& mat,bool,int code )
{
    (void)mat;
    (void)code;
    return 0.0;
};
Real test_function_max_bernoulli_b2n::eval_scalar(const Scalar& s,bool,int code  )
{
    (void)s;
    (void)code;

    //just test, that the fuction compiles and not throws
    Integer v1  = max_bernoulli_b2n();
    Integer v2  = fmax_bernoulli_b2n();
    Integer w1  = max_bernoulli_b2n<Real>();
    Integer w2  = max_bernoulli_b2n<Float>();

    Real out    = norm_1(v1 - w1);
    out         +=norm_1(v2 - w2);

    max_bernoulli_b2n<Integer>();
    max_bernoulli_b2n<Complex>();
    max_bernoulli_b2n<Float_complex>();

    return out;
};

Real test_function_prime::eval_mat(const Matrix& mat,bool,int code )
{
    (void)mat;
    (void)code;
    return 0.0;
};
Real test_function_prime::eval_scalar(const Scalar& s,bool,int code  )
{
    (void)s;
    (void)code;
    //just test, that the fuction compiles and not throws
    prime(1);
    prime(2);
    prime(prime_max_count());

    return 0.0;
};

Real test_function_prime_max_count::eval_mat(const Matrix& mat,bool,int code )
{
    (void)mat;
    (void)code;
    return 0.0;
};
Real test_function_prime_max_count::eval_scalar(const Scalar& s,bool,int code  )
{
    (void)s;
    (void)code;
    //just test, that the fuction compiles and not throws
    prime_max_count();

    return 0.0;
};

Real test_function_factorial::eval_mat(const Matrix& mat,bool,int code )
{
    (void)mat;
    (void)code;
    return 0.0;
};
Real test_function_factorial::eval_scalar(const Scalar& s,bool,int code  )
{
    (void)s;
    (void)code;
    Real v1     = factorial(0);
    Real v2     = factorial(1);

    Float f1    = ffactorial(0);
    Float f2    = ffactorial(0);

    Real w1     = factorial<Real>(0);
    Real w2     = factorial<Real>(1);
    Real w3     = factorial<Float>(0);
    Real w4     = factorial<Float>(1);

    factorial<Integer>(0);
    factorial<Complex>(0);
    factorial<Float_complex>(0);

    Real out    = norm_1(v1 - f1);
    out         += norm_1(v2 - f2);
    out         += norm_1(v1 - w1);
    out         += norm_1(v2 - w2);
    out         += norm_1(f1 - w3);
    out         += norm_1(f2 - w4);

    return out;
};

Real test_function_double_factorial::eval_mat(const Matrix& mat,bool,int code )
{
    (void)mat;
    (void)code;
    return 0.0;
};
Real test_function_double_factorial::eval_scalar(const Scalar& s,bool,int code  )
{
    (void)s;
    (void)code;
    Real v1     = double_factorial(0);
    Real v2     = double_factorial(1);

    Float f1    = fdouble_factorial(0);
    Float f2    = fdouble_factorial(0);

    Real w1     = double_factorial<Real>(0);
    Real w2     = double_factorial<Real>(1);
    Real w3     = double_factorial<Float>(0);
    Real w4     = double_factorial<Float>(1);

    double_factorial<Integer>(0);
    double_factorial<Complex>(0);
    double_factorial<Float_complex>(0);

    Real out    = norm_1(v1 - f1);
    out         += norm_1(v2 - f2);
    out         += norm_1(v1 - w1);
    out         += norm_1(v2 - w2);
    out         += norm_1(f1 - w3);
    out         += norm_1(f2 - w4);

    return out;
};

Real test_function_rising_factorial::eval_mat(const Matrix& mat,bool,int code )
{
    (void)mat;
    (void)code;
    return 0.0;
};
Real test_function_rising_factorial::eval_scalar(const Scalar& s,bool,int code  )
{
    (void)code;
    return eval_scalar_impl<test_function_rising_factorial>(s);
};

Real test_function_falling_factorial::eval_mat(const Matrix& mat,bool,int code )
{
    (void)mat;
    (void)code;
    return 0.0;
};
Real test_function_falling_factorial::eval_scalar(const Scalar& s,bool,int code  )
{
    (void)code;
    return eval_scalar_impl<test_function_falling_factorial>(s);
};

Real test_function_binomial_coefficient::eval_mat(const Matrix& mat,bool,int code )
{
    (void)mat;
    (void)code;
    return 0.0;
};
Real test_function_binomial_coefficient::eval_scalar(const Scalar& s,bool,int code  )
{
    (void)s;
    (void)code;
    Real v1     = binomial_coefficient(0, 0);
    Real v2     = binomial_coefficient(1, 0);

    Float f1    = fbinomial_coefficient(0, 0);
    Float f2    = fbinomial_coefficient(0, 0);

    Real w1     = binomial_coefficient<Real>(0, 0);
    Real w2     = binomial_coefficient<Real>(1, 0);
    Real w3     = binomial_coefficient<Float>(0, 0);
    Real w4     = binomial_coefficient<Float>(1, 0);

    binomial_coefficient<Integer>(0, 0);
    binomial_coefficient<Complex>(0, 0);
    binomial_coefficient<Float_complex>(0, 0);

    Real out    = norm_1(v1 - f1);
    out         += norm_1(v2 - f2);
    out         += norm_1(v1 - w1);
    out         += norm_1(v2 - w2);
    out         += norm_1(f1 - w3);
    out         += norm_1(f2 - w4);

    return out;
};

Real test_function_eval_scalar_func::eval_mat(const Matrix& mat,bool,int code )
{
    (void)code;
    try
    {		
        Real dif = 0;
        {
            func_sin fs;
            Matrix out		= eval_scalar_func(mat,fs);
            check_struct(out);

            dif				+= norm_1(out - sin(mat));
        };
        {
            func_sin fs;
            Matrix out		= eval_scalar_func(mat,fs,fs);
            check_struct(out);

            dif				+= norm_1(out - sin(mat));
        };
        {
            func_cos fc;
            Matrix out		= eval_scalar_func(mat,fc);
            check_struct(out);

            dif				+= norm_1(out - cos(mat));
        };
        {
            func_cos fc;
            Matrix out		= eval_scalar_func(mat,fc,fc);
            check_struct(out);

            dif				+= norm_1(out - cos(mat));
        };
        return dif;
    }
    catch(const std::exception& ex)
    {
        m_is_error  = true;
        m_error     = ex.what();
        return 1.;
    }
    catch(...)
    {
        m_is_error  = true;
        m_error     = "unknown error";
        return 1.;
    };
};
Real test_function_eval_scalar_func::eval_scalar(const Scalar& s,bool,int code )
{
    (void)code;
    return eval_scalar_impl<test_function_eval_scalar_func>(s);
};


Real test_function_is_true::eval_mat(const Matrix& mat,bool,int code )
{
    (void)code;
    Matrix out_full = is_true(full(mat));	
    try
    {		
        Matrix out		= is_true(mat);
        check_struct(out);

        Matrix out2		= (mat != 0);
        check_struct(out2);
        
        Real dif		= norm_1(out - out_full);
        dif				= dif + norm_1(out - out2);
        return dif;
    }
    catch(const std::exception& ex)
    {
        m_is_error  = true;
        m_error     = ex.what();
        return 1.;
    }
    catch(...)
    {
        m_is_error  = true;
        m_error     = "unknown error";
        return 1.;
    };
};
Real test_function_is_true::eval_scalar(const Scalar& s,bool,int code  )
{
    (void)code;
    return eval_scalar_impl<test_function_is_true>(s);
};

Real test_function_is_false::eval_mat(const Matrix& mat,bool,int code )
{
    (void)code;
    Matrix out_full = is_false(full(mat));	
    try
    {		
        Matrix out		= is_false(mat);
        check_struct(out);

        Matrix out2		= (mat == 0);
        check_struct(out2);

        Real dif		= norm_1(out - out_full);
        dif				= dif + norm_1(out - out2);
        dif				= dif + norm_1(is_false(out) - is_true(mat));
        return dif;
    }
    catch(const std::exception& ex)
    {
        m_is_error  = true;
        m_error     = ex.what();
        return 1.;
    }
    catch(...)
    {
        m_is_error  = true;
        m_error     = "unknown error";
        return 1.;
    };
};
Real test_function_is_false::eval_scalar(const Scalar& s,bool,int code  )
{
    (void)code;
    return eval_scalar_impl<test_function_is_false>(s);
};

Real test_function_abs2::eval_mat(const Matrix& mat,bool,int code )
{
    (void)code;
    try
    {		
        Matrix res1 = abs(mat);
        Matrix res2 = abs2(mat);

        res1        = matcl::mul(res1,res1);

        Real dif	= norm_1(res1 - res2);

        if (dif < norm_1(res2) * constants::eps(mat.get_value_code()) * 10.)
            dif     = 0.0;

        check_struct(res2);
        return dif;
    }
    catch(const std::exception& ex)
    {
        m_is_error  = true;
        m_error     = ex.what();
        return 1.;
    }
    catch(...)
    {
        m_is_error  = true;
        m_error     = "unknown error";
        return 1.;
    };
};
Real test_function_abs2::eval_scalar(const Scalar& s,bool,int code  )
{
    (void)code;
    return eval_scalar_impl<test_function_abs2>(s);
};

};};
