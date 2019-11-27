/*
 *  This file is a part of Matrix Computation Library (MATCL)
 *
 *  Copyright (c) Pawe³ Kowal 2018
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
#include "test_functions_bin.h"

#include "test/test_matcl/framework/matrix_set/matrix_set_bin_1.h"
#include "test/test_matcl/framework/matrix_set/test_options.h"
#include "test/test_matcl/framework/data_struct/test_selector.h"
#include "matcl-core/IO/logger.h"
#include "matcl-core/matrix/enums.h"
#include "matcl-specfunc/matcl_specfunc.h"
#include <boost/thread.hpp>

namespace matcl { namespace test
{

class test_binary
{
    bin_functions_list&				tf;
    const test::options&			opts;
    Integer                         thread_id;

    public:
        test_binary(bin_functions_list& tf, const test::options& opts, Integer id)
            :tf(tf),opts(opts), thread_id(id)
        {};

        test_binary(const test_binary& tb)
            :tf(tb.tf),opts(tb.opts), thread_id(tb.thread_id)
        {};

        void make()
        {
            /*
            Integer code    = 1233;
            Matrix mat1     = tf.get_matrix(code).first;
            Matrix mat2     = tf.get_matrix(code).second;

            matcl::disp(mat1);
            matcl::disp(mat2);
            */

            tf.make(opts,thread_id);
        };

        void operator()()
        {
            make();
        };

    private:		
        test_binary& operator=(const test_binary&) = delete;
};

class test_mult
{
    bin_functions_list&				tf;
    const test::options&			opts;
    Integer                         thread_id;

    public:
        test_mult(bin_functions_list& tf, const test::options& opts, Integer id)
            :tf(tf),opts(opts), thread_id(id)
        {};

        test_mult(const test_mult& tb)
            :tf(tb.tf),opts(tb.opts), thread_id(tb.thread_id)
        {};

        void make()
        {
            /*
            using matrix_pair = test::matrix_set_bin::matrix_pair;            

            Integer code    = 75875;
            matrix_pair mp  = tf.get_matrix(code);

            dynamic_mat_set&    m_ds = tf.m_ds;

            Matrix mat10    = mp.first;
            Matrix mat20    = mp.second;

            trans_type m_ta = trans_type::conj_trans;
            trans_type m_tb = trans_type::no_trans;

            Matrix mat1     = trans(mat10, m_ta);
            Matrix mat2     = trans(mat20, m_tb);    

            disp(mat1);
            disp(mat2);

            Matrix alpha    = test_function_gemm::rand_scalar(code);
            Matrix beta     = test_function_gemm::rand_scalar(code + 1);
            value_code vc_C = test_function_gemm::get_value_code(alpha, beta, mat1, mat2);

            Integer r       = (m_ta == trans_type::no_trans) ? mat1.rows() : mat1.cols();
            Integer c       = (m_tb == trans_type::no_trans) ? mat2.cols() : mat2.rows();

            if (mat1.is_scalar() && mat2.is_scalar())
            {
                r           += code % 3;
                c           += code % 3;
            };

            Integer r2      = r + 2 + (code % 5);
            Integer c2      = c + 2 + (code % 5);

            Integer rs      = 1 + (code % 2);
            Integer cs      = 1 + (code % 3);

            Integer re      = rs + r - 1;
            Integer ce      = cs + c - 1;

            Matrix C        = test_function_gemm::rand_mat_C(r2, c2, vc_C, code + 2, m_ds);
            vc_C            = C.get_value_code();

            colon cr(rs, re);
            colon cc(cs, ce);

            Matrix alpha_c  = convert_value(alpha, vc_C);
            Matrix out_full = (alpha_c * trans(mat1, m_ta)) * trans(mat2, m_tb) + beta * C(cr, cc);

            disp(C);
            disp(out_full);

            Matrix out      = C;
            Integer ver     = test_function_gemm_sub::rand_colon_ver(rs, re, cs, ce, code + 3);

            switch (ver)
            {
                case 0:
                    gemm(alpha, mat1, mat2, m_ta, m_tb, beta, out(rs, cs));
                    break;
                case 1:
                    gemm(alpha, mat1, mat2, m_ta, m_tb, beta, out(rs, cc));
                    break;
                case 2:
                    gemm(alpha, mat1, mat2, m_ta, m_tb, beta, out(cr, cs));
                    break;
                case 3:        
                    gemm(alpha, mat1, mat2, m_ta, m_tb, beta, out(cr, cc));
                    break;
            };

            Matrix out2     = out(cr, cc);
            disp(out2);
            */

            tf.make_mult(opts);
        };
        void operator()()
        {
            make();
        };

    private:		
        test_mult& operator=(const test_mult&) = delete;
};

class test_kron
{
    bin_functions_list&				tf;
    const test::options&			opts;
    Integer                         thread_id;

    public:
        test_kron(bin_functions_list& tf, const test::options& opts, Integer id)
            :tf(tf),opts(opts), thread_id(id)
        {};

        test_kron(const test_kron& tb)
            :tf(tb.tf),opts(tb.opts), thread_id(tb.thread_id)
        {};

        void make()
        {
            /*
            using matrix_pair   = test::matrix_set_bin::matrix_pair;
            matrix_pair mp      = tf.get_matrix(3);

            Matrix mat1         = mp.first;
            Matrix mat2         = mp.second;

            disp(mat1);
            disp(mat2);
            */

            tf.make_kron(opts);
        };
        void operator()()
        {
            make();
        };

    private:		
        test_kron& operator=(const test_kron&) = delete;
};

void test_bin_st(const rand_matrix_ptr& rand)
{
    test::options opts;
    try
    {
        opts.show_partial_res   = false;
        opts.first_matrix_code  = 0;
        opts.show_memleaks      = true;

        test::mat_set_bin_1 ms(rand);
        dynamic_mat_set ds(rand);
        bin_functions_list tf(ms, ds);

        test_binary tbin (tf,opts,0);

        tbin.make();
    }
    catch(const std::exception& ex)
    {
        matcl::out_stream <<ex.what();
    };
};

void test_bin_mt(const rand_matrix_ptr& rand)
{
    test::options opts;
    try
    {
        opts.show_partial_res   = 0;
        opts.first_matrix_code  = 0;
        opts.show_memleaks      = false;

        test::mat_set_bin_2 ms1(rand);
        dynamic_mat_set ds(rand);
        bin_functions_list tf(ms1, ds);

        boost::thread_group tg;

        for (int i = 0; i < 20; i++)
            tg.create_thread(test_binary(tf,opts,i));

        tg.join_all();
    }
    catch(const std::exception& ex)
    {
        matcl::out_stream <<ex.what();
    };
};

void test_mult_st(const rand_matrix_ptr& rand)
{
    test::options opts;
    try
    {
        opts.show_partial_res   = false;
        opts.first_matrix_code  = 0;
        opts.show_memleaks      = true;

        test::mat_set_bin_mult_1 ms(rand);
        dynamic_mat_set ds(rand);

        bin_functions_list tf(ms, ds);

        test_mult tbin (tf,opts,0);

        tbin.make();
    }
    catch(const std::exception& ex)
    {
        matcl::out_stream <<ex.what();
    };
};

void test_mult_mt(const rand_matrix_ptr& rand)
{
    test::options opts;
    try
    {
        opts.show_partial_res   = 0;
        opts.first_matrix_code  = 0;
        opts.show_memleaks      = false;

        test::mat_set_bin_mult_2 ms1(rand);
        dynamic_mat_set ds(rand);

        bin_functions_list tf(ms1, ds);

        boost::thread_group tg;

        for (int i = 0; i < 20; i++)
            tg.create_thread(test_mult(tf,opts,i));

        tg.join_all();
    }
    catch(const std::exception& ex)
    {
        matcl::out_stream <<ex.what();
    };
};

void test_kron_st(const rand_matrix_ptr& rand)
{
    test::options opts;
    try
    {
        opts.show_partial_res   = false;
        opts.first_matrix_code  = 0;
        opts.show_memleaks      = true;

        test::mat_set_bin_kron_1 ms(rand);
        dynamic_mat_set ds(rand);

        bin_functions_list tf(ms, ds);

        test_kron tbin (tf,opts,0);

        tbin.make();
    }
    catch(const std::exception& ex)
    {
        matcl::out_stream <<ex.what();
    };
};

void test_kron_mt(const rand_matrix_ptr& rand)
{
    test::options opts;
    try
    {
        opts.show_partial_res   = 0;
        opts.first_matrix_code  = 0;
        opts.show_memleaks      = false;

        test::mat_set_bin_kron_2 ms1(rand);
        dynamic_mat_set ds(rand);

        bin_functions_list tf(ms1, ds);

        boost::thread_group tg;

        for (int i = 0; i < 20; i++)
            tg.create_thread(test_kron(tf,opts,i));

        tg.join_all();
    }
    catch(const std::exception& ex)
    {
        matcl::out_stream <<ex.what();
    };
};

static const Integer max_int = 1000;

void bin_functions_list::make(options opts, Integer thread_id)
{
    (void)thread_id;

    m_options = opts;       

    SELECT_TEST (3, test_op_plus());
    SELECT_TEST (3, test_op_minus());
    SELECT_TEST (3, test_op_mult());    
    SELECT_TEST (3, test_op_div());

    SELECT_TEST (3, test_plus());
    SELECT_TEST (3, test_minus());
    SELECT_TEST (3, test_mul());
    SELECT_TEST (3, test_div());
    SELECT_TEST (3, test_div_0());
    SELECT_TEST (3, test_div_1());
    SELECT_TEST (3, test_idiv());
    SELECT_TEST (3, test_pow());
    SELECT_TEST (3, test_pow_c());

    SELECT_TEST (3, test_max_bin());
    SELECT_TEST (3, test_min_bin());
    SELECT_TEST (3, test_xor());    

    SELECT_TEST (3, test_rem());
    SELECT_TEST (3, test_mod());	
    SELECT_TEST (3, test_op_or());
    SELECT_TEST (3, test_op_and());
    SELECT_TEST (3, test_op_eeq());
    SELECT_TEST (3, test_op_neq());
    SELECT_TEST (3, test_eeq_nan());
    SELECT_TEST (3, test_neq_nan());
    SELECT_TEST (3, test_op_lt());
    SELECT_TEST (3, test_op_leq());
    SELECT_TEST (3, test_op_gt());
    SELECT_TEST (3, test_op_geq());	

    SELECT_TEST (3, test_atan2());
    SELECT_TEST (3, test_hypot());

    SELECT_TEST (3, test_copysign());
    SELECT_TEST (3, test_fdim());
    SELECT_TEST (3, test_nextafter());
    SELECT_TEST (3, test_powm1());

    SELECT_TEST (3, test_kron());
    SELECT_TEST (3, test_eval_bin());
    
    SELECT_TEST (3, test_chain_mult());    

    SELECT_TEST (3, test_fma());
    SELECT_TEST (3, test_fms());
    SELECT_TEST (3, test_dot2_ac());

    //TODO: remove
    SELECT_TEST (3, test_beta());        
};

void bin_functions_list::make_mult(options opts)
{
    m_options = opts;
        
    SELECT_TEST(3, test_mmul_abs());
    SELECT_TEST(3, test_gemm());
    SELECT_TEST(3, test_gemm_sub());
    SELECT_TEST(3, test_mmul());
    SELECT_TEST(3, test_mmul_ext());

    SELECT_TEST /**/(3, test_op_mult());	  //disabled scons parsing
    SELECT_TEST /**/(3, test_chain_mult());  //disabled scons parsing
    SELECT_TEST /**/(3, test_kron());  //disabled scons parsing
};

void bin_functions_list::make_kron(options opts)
{
    m_options = opts;

    SELECT_TEST /**/ (3, test_kron());  //disabled scons parsing
};

bin_functions_list::matrix_pair bin_functions_list::get_matrix(int code) const
{
    return m_tests.get_matrix(code);
};
bin_functions_list::scalar_pair bin_functions_list::get_scalar(int code) const
{
    return m_tests.get_scalar(code);
};

void bin_functions_list::test_op_plus()
{
    Real out = 0.;
    test_function_op_plus tf;
    out += m_tests.make(&tf,m_options);

    if (out == 0.)
        matcl::out_stream << std::string() +   "operator+: OK" + "\n";
    else
        matcl::out_stream << std::string() +   "operator+: FAILED"  + "\n";
};

void bin_functions_list::test_op_minus()
{
    Real out = 0.;
    test_function_op_minus tf;
    out += m_tests.make(&tf,m_options);

    if (out == 0.)
        matcl::out_stream << std::string() +   "operator-: OK" + "\n";
    else
        matcl::out_stream << std::string() +   "operator-: FAILED"  + "\n";
};

void bin_functions_list::test_op_mult()
{
    Real out = 0.;
    test_function_op_mult tf;
    out += m_tests.make(&tf,m_options);

    if (out == 0.)
        matcl::out_stream << std::string() +   "operator*: OK" + "\n";
    else
        matcl::out_stream << std::string() +   "operator*: FAILED"  + "\n";
};

void bin_functions_list::test_mmul()
{
    Real out = 0.;

    {
        test_function_mmul tf(trans_type::no_trans, trans_type::no_trans);
        out += m_tests.make(&tf,m_options);
    }
    {
        test_function_mmul tf(trans_type::no_trans, trans_type::trans);
        out += m_tests.make(&tf,m_options);
    }    
    {
        test_function_mmul tf(trans_type::no_trans, trans_type::conj_trans);
        out += m_tests.make(&tf,m_options);
    }    

    {
        test_function_mmul tf(trans_type::trans, trans_type::no_trans);
        out += m_tests.make(&tf,m_options);
    }
    {
        test_function_mmul tf(trans_type::trans, trans_type::trans);
        out += m_tests.make(&tf,m_options);
    }
    {
        test_function_mmul tf(trans_type::trans, trans_type::conj_trans);
        out += m_tests.make(&tf,m_options);
    }    
    {
        test_function_mmul tf(trans_type::conj_trans, trans_type::no_trans);
        out += m_tests.make(&tf,m_options);
    }    
    {
        test_function_mmul tf(trans_type::conj_trans, trans_type::trans);
        out += m_tests.make(&tf,m_options);
    }
    {
        test_function_mmul tf(trans_type::conj_trans, trans_type::conj_trans);
        out += m_tests.make(&tf,m_options);
    }

    if (out == 0.)
        matcl::out_stream << std::string() +   "mmul: OK" + "\n";
    else
        matcl::out_stream << std::string() +   "mmul: FAILED"  + "\n";
};

void bin_functions_list::test_gemm()
{
    Real out = 0.;

    {
        test_function_gemm tf(trans_type::no_trans, trans_type::no_trans, m_ds);
        out += m_tests.make(&tf,m_options);
    } 
    {
        test_function_gemm tf(trans_type::no_trans, trans_type::trans, m_ds);
        out += m_tests.make(&tf,m_options);
    }    
    {
        test_function_gemm tf(trans_type::no_trans, trans_type::conj_trans, m_ds);
        out += m_tests.make(&tf,m_options);
    }    

    {
        test_function_gemm tf(trans_type::trans, trans_type::no_trans, m_ds);
        out += m_tests.make(&tf,m_options);
    }
    {
        test_function_gemm tf(trans_type::trans, trans_type::trans, m_ds);
        out += m_tests.make(&tf,m_options);
    }
    {
        test_function_gemm tf(trans_type::trans, trans_type::conj_trans, m_ds);
        out += m_tests.make(&tf,m_options);
    }    
    {
        test_function_gemm tf(trans_type::conj_trans, trans_type::no_trans, m_ds);
        out += m_tests.make(&tf,m_options);
    }    
    {
        test_function_gemm tf(trans_type::conj_trans, trans_type::trans, m_ds);
        out += m_tests.make(&tf,m_options);
    }
    {
        test_function_gemm tf(trans_type::conj_trans, trans_type::conj_trans, m_ds);
        out += m_tests.make(&tf,m_options);
    }

    if (out == 0.)
        matcl::out_stream << std::string() +   "gemm: OK" + "\n";
    else
        matcl::out_stream << std::string() +   "gemm: FAILED"  + "\n";
};

void bin_functions_list::test_gemm_sub()
{
    Real out = 0.;

    {
        test_function_gemm_sub tf(trans_type::no_trans, trans_type::no_trans, m_ds);
        out += m_tests.make(&tf,m_options);
    } 
    {
        test_function_gemm_sub tf(trans_type::no_trans, trans_type::trans, m_ds);
        out += m_tests.make(&tf,m_options);
    }    
    {
        test_function_gemm_sub tf(trans_type::no_trans, trans_type::conj_trans, m_ds);
        out += m_tests.make(&tf,m_options);
    }    

    {
        test_function_gemm_sub tf(trans_type::trans, trans_type::no_trans, m_ds);
        out += m_tests.make(&tf,m_options);
    }
    {
        test_function_gemm_sub tf(trans_type::trans, trans_type::trans, m_ds);
        out += m_tests.make(&tf,m_options);
    }
    {
        test_function_gemm_sub tf(trans_type::trans, trans_type::conj_trans, m_ds);
        out += m_tests.make(&tf,m_options);
    }    
    {
        test_function_gemm_sub tf(trans_type::conj_trans, trans_type::no_trans, m_ds);
        out += m_tests.make(&tf,m_options);
    }    
    {
        test_function_gemm_sub tf(trans_type::conj_trans, trans_type::trans, m_ds);
        out += m_tests.make(&tf,m_options);
    }
    {
        test_function_gemm_sub tf(trans_type::conj_trans, trans_type::conj_trans, m_ds);
        out += m_tests.make(&tf,m_options);
    }

    if (out == 0.)
        matcl::out_stream << std::string() +   "gemm_sub: OK" + "\n";
    else
        matcl::out_stream << std::string() +   "gemm_sub: FAILED"  + "\n";
};

void bin_functions_list::test_mmul_abs()
{
    Real out = 0.;

    {
        test_function_mmul_abs tf(trans_type::no_trans, m_ds);
        out += m_tests.make(&tf,m_options);
    } 
    {
        test_function_mmul_abs tf(trans_type::trans, m_ds);
        out += m_tests.make(&tf,m_options);
    } 
    {
        test_function_mmul_abs tf(trans_type::conj_trans, m_ds);
        out += m_tests.make(&tf,m_options);
    }    

    if (out == 0.)
        matcl::out_stream << std::string() +   "mmul_abs: OK" + "\n";
    else
        matcl::out_stream << std::string() +   "mmul_abs: FAILED"  + "\n";
};


void bin_functions_list::test_mmul_ext()
{
    Real out = 0.;
    {
        test_function_mmul_ext tf(trans_type_ext::no_trans, trans_type_ext::conj);
        out += m_tests.make(&tf,m_options);
    }
    {
        test_function_mmul_ext tf(trans_type_ext::trans, trans_type_ext::conj);
        out += m_tests.make(&tf,m_options);
    }
    {
        test_function_mmul_ext tf(trans_type_ext::conj_trans, trans_type_ext::conj);
        out += m_tests.make(&tf,m_options);
    }

    {
        test_function_mmul_ext tf(trans_type_ext::conj, trans_type_ext::no_trans);
        out += m_tests.make(&tf,m_options);
    }
    {
        test_function_mmul_ext tf(trans_type_ext::conj, trans_type_ext::trans);
        out += m_tests.make(&tf,m_options);
    }
    {
        test_function_mmul_ext tf(trans_type_ext::conj, trans_type_ext::conj_trans);
        out += m_tests.make(&tf,m_options);
    }
    {
        test_function_mmul_ext tf(trans_type_ext::conj, trans_type_ext::conj);
        out += m_tests.make(&tf,m_options);
    }

    if (out == 0.)
        matcl::out_stream << std::string() +   "mmul_ext: OK" + "\n";
    else
        matcl::out_stream << std::string() +   "mmul_ext: FAILED"  + "\n";
};
void bin_functions_list::test_chain_mult()
{
    Real out = 0.;
    test_function_chain_mult tf;
    out += m_tests.make(&tf,m_options);

    if (out == 0.)
        matcl::out_stream << std::string() +   "chain mult: OK" + "\n";
    else
        matcl::out_stream << std::string() +   "chain mult: FAILED"  + "\n";
};

void bin_functions_list::test_op_div()
{
    Real out = 0.;
    test_function_op_div tf;
    out += m_tests.make(&tf,m_options);

    if (out == 0.)
        matcl::out_stream << std::string() +   "operator/: OK" + "\n";
    else
        matcl::out_stream << std::string() +   "operator/: FAILED"  + "\n";
};

void bin_functions_list::test_plus()
{
    Real out = 0.;
    test_function_plus tf;
    out += m_tests.make(&tf,m_options);

    if (out == 0.)
        matcl::out_stream << std::string() +   "plus: OK" + "\n";
    else
        matcl::out_stream << std::string() +   "plus: FAILED"  + "\n";
};

void bin_functions_list::test_minus()
{
    Real out = 0.;
    test_function_minus tf;
    out += m_tests.make(&tf,m_options);

    if (out == 0.)
        matcl::out_stream << std::string() +   "minus: OK" + "\n";
    else
        matcl::out_stream << std::string() +   "minus: FAILED"  + "\n";
};

void bin_functions_list::test_mul()
{
    Real out = 0.;
    test_function_mul tf;
    out += m_tests.make(&tf,m_options);

    if (out == 0.)
        matcl::out_stream << std::string() +   "mul: OK" + "\n";
    else
        matcl::out_stream << std::string() +   "mul: FAILED"  + "\n";
};

void bin_functions_list::test_div()
{
    Real out = 0.;
    test_function_div tf;
    out += m_tests.make(&tf,m_options);

    if (out == 0.)
        matcl::out_stream << std::string() +   "div: OK" + "\n";
    else
        matcl::out_stream << std::string() +   "div: FAILED"  + "\n";
};

void bin_functions_list::test_div_0()
{
    Real out = 0.;
    test_function_div_0 tf;
    out += m_tests.make(&tf,m_options);

    if (out == 0.)
        matcl::out_stream << std::string() +   "div_0: OK" + "\n";
    else
        matcl::out_stream << std::string() +   "div_0: FAILED"  + "\n";
};

void bin_functions_list::test_div_1()
{
    Real out = 0.;
    test_function_div_1 tf;
    out += m_tests.make(&tf,m_options);

    if (out == 0.)
        matcl::out_stream << std::string() +   "div_1: OK" + "\n";
    else
        matcl::out_stream << std::string() +   "div_1: FAILED"  + "\n";
};

void bin_functions_list::test_idiv()
{
    Real out = 0.;
    test_function_idiv tf;
    out += m_tests.make(&tf,m_options);

    if (out == 0.)
        matcl::out_stream << std::string() +   "idiv: OK" + "\n";
    else
        matcl::out_stream << std::string() +   "idiv: FAILED"  + "\n";
};

void bin_functions_list::test_beta()
{
    Real out = 0.;
    test_function_beta tf;
    out += m_tests.make(&tf,m_options);

    if (out == 0.)
        matcl::out_stream << std::string() +   "beta: OK" + "\n";
    else
        matcl::out_stream << std::string() +   "beta: FAILED"  + "\n";
};

void bin_functions_list::test_eval_bin()
{
    Real out = 0.;
    test_function_eval_bin tf;
    out += m_tests.make(&tf,m_options);

    if (out == 0.)
        matcl::out_stream << std::string() +   "eval_bin: OK" + "\n";
    else
        matcl::out_stream << std::string() +   "eval_bin: FAILED"  + "\n";
};

void bin_functions_list::test_atan2()
{
    Real out = 0.;
    test_function_atan2 tf;
    out += m_tests.make(&tf,m_options);

    if (out == 0.)
        matcl::out_stream << std::string() +   "atan2: OK" + "\n";
    else
        matcl::out_stream << std::string() +   "atan2: FAILED"  + "\n";
};

void bin_functions_list::test_hypot()
{
    Real out = 0.;
    test_function_hypot tf;
    out += m_tests.make(&tf,m_options);

    if (out == 0.)
        matcl::out_stream << std::string() +   "hypot: OK" + "\n";
    else
        matcl::out_stream << std::string() +   "hypot: FAILED"  + "\n";
};

void bin_functions_list::test_copysign()
{
    Real out = 0.;
    test_function_copysign tf;
    out += m_tests.make(&tf,m_options);

    if (out == 0.)
        matcl::out_stream << std::string() +   "copysign: OK" + "\n";
    else
        matcl::out_stream << std::string() +   "copysign: FAILED"  + "\n";
};
void bin_functions_list::test_fdim()
{
    Real out = 0.;
    test_function_fdim tf;
    out += m_tests.make(&tf,m_options);

    if (out == 0.)
        matcl::out_stream << std::string() +   "fdim: OK" + "\n";
    else
        matcl::out_stream << std::string() +   "fdim: FAILED"  + "\n";
};
void bin_functions_list::test_nextafter()
{
    Real out = 0.;
    test_function_nextafter tf;
    out += m_tests.make(&tf,m_options);

    if (out == 0.)
        matcl::out_stream << std::string() +   "nextafter: OK" + "\n";
    else
        matcl::out_stream << std::string() +   "nextafter: FAILED"  + "\n";
};
void bin_functions_list::test_fma()
{
    Real out = 0.;
    test_function_fma tf;
    out += m_tests.make(&tf,m_options);

    if (out == 0.)
        matcl::out_stream << std::string() +   "fma: OK" + "\n";
    else
        matcl::out_stream << std::string() +   "fma: FAILED"  + "\n";
};
void bin_functions_list::test_fms()
{
    Real out = 0.;
    test_function_fms tf;
    out += m_tests.make(&tf,m_options);

    if (out == 0.)
        matcl::out_stream << std::string() +   "fms: OK" + "\n";
    else
        matcl::out_stream << std::string() +   "fms: FAILED"  + "\n";
};
void bin_functions_list::test_dot2_ac()
{
    Real out = 0.;
    test_function_dot2_ac tf;
    out += m_tests.make(&tf,m_options);

    if (out == 0.)
        matcl::out_stream << std::string() +   "dot2_ac: OK" + "\n";
    else
        matcl::out_stream << std::string() +   "dot2_ac: FAILED"  + "\n";
};

void bin_functions_list::test_powm1()
{
    Real out = 0.;
    test_function_powm1 tf;
    out += m_tests.make(&tf,m_options);

    if (out == 0.)
        matcl::out_stream << std::string() +   "powm1: OK" + "\n";
    else
        matcl::out_stream << std::string() +   "powm1: FAILED"  + "\n";
};

void bin_functions_list::test_pow()
{
    Real out = 0.;
    test_function_pow tf;
    out += m_tests.make(&tf,m_options);

    if (out == 0.)
        matcl::out_stream << std::string() +   "pow: OK" + "\n";
    else
        matcl::out_stream << std::string() +   "pow: FAILED"  + "\n";
};

void bin_functions_list::test_pow_c()
{
    Real out = 0.;
    test_function_pow_c tf;
    out += m_tests.make(&tf,m_options);

    if (out == 0.)
        matcl::out_stream << std::string() +   "pow_c: OK" + "\n";
    else
        matcl::out_stream << std::string() +   "pow_c: FAILED"  + "\n";
};

void bin_functions_list::test_max_bin()
{
    Real out = 0.;
    test_function_max_bin tf;
    out += m_tests.make(&tf,m_options);

    if (out == 0.)
        matcl::out_stream << std::string() +   "bin max: OK" + "\n";
    else
        matcl::out_stream << std::string() +   "bin max: FAILED"  + "\n";
};

void bin_functions_list::test_min_bin()
{
    Real out = 0.;
    test_function_min_bin tf;
    out += m_tests.make(&tf,m_options);

    if (out == 0.)
        matcl::out_stream << std::string() +   "bin min: OK" + "\n";
    else
        matcl::out_stream << std::string() +   "bin min: FAILED"  + "\n";
};

void bin_functions_list::test_xor()
{
    Real out = 0.;
    test_function_xor tf;
    out += m_tests.make(&tf,m_options);

    if (out == 0.)
        matcl::out_stream << std::string() +   "xor: OK" + "\n";
    else
        matcl::out_stream << std::string() +   "xor: FAILED"  + "\n";
};

void bin_functions_list::test_rem()
{
    Real out = 0.;
    test_function_rem tf;
    out += m_tests.make(&tf,m_options);

    if (out == 0.)
        matcl::out_stream << std::string() +   "rem: OK" + "\n";
    else
        matcl::out_stream << std::string() +   "rem: FAILED"  + "\n";
};

void bin_functions_list::test_mod()
{
    Real out = 0.;
    test_function_mod tf;
    out += m_tests.make(&tf,m_options);

    if (out == 0.)
        matcl::out_stream << std::string() +   "mod: OK" + "\n";
    else
        matcl::out_stream << std::string() +   "mod: FAILED"  + "\n";
};

void bin_functions_list::test_kron()
{
    Real out = 0.;
    test_function_kron tf;
    out += m_tests.make(&tf,m_options);

    if (out == 0.)
        matcl::out_stream << std::string() +   "kron: OK" + "\n";
    else
        matcl::out_stream << std::string() +   "kron: FAILED"  + "\n";
};

void bin_functions_list::test_op_or()
{
    Real out = 0.;
    test_function_op_or tf;
    out += m_tests.make(&tf,m_options);

    if (out == 0.)
        matcl::out_stream << std::string() +   "operator|: OK" + "\n";
    else
        matcl::out_stream << std::string() +   "operator|: FAILED"  + "\n";
};

void bin_functions_list::test_op_and()
{
    Real out = 0.;
    test_function_op_and tf;
    out += m_tests.make(&tf,m_options);

    if (out == 0.)
        matcl::out_stream << std::string() +   "operator&: OK" + "\n";
    else
        matcl::out_stream << std::string() +   "operator&: FAILED"  + "\n";
};

void bin_functions_list::test_op_eeq()
{
    Real out = 0.;
    test_function_op_eeq tf;
    out += m_tests.make(&tf,m_options);

    if (out == 0.)
        matcl::out_stream << std::string() +   "operator==: OK" + "\n";
    else
        matcl::out_stream << std::string() +   "operator==: FAILED"  + "\n";
};

void bin_functions_list::test_op_neq()
{
    Real out = 0.;
    test_function_op_neq tf;
    out += m_tests.make(&tf,m_options);

    if (out == 0.)
        matcl::out_stream << std::string() +   "operator!=: OK" + "\n";
    else
        matcl::out_stream << std::string() +   "operator!=: FAILED"  + "\n";
};

void bin_functions_list::test_eeq_nan()
{
    Real out = 0.;
    test_function_eeq_nan tf;
    out += m_tests.make(&tf,m_options);

    if (out == 0.)
        matcl::out_stream << std::string() +   "eeq_nan: OK" + "\n";
    else
        matcl::out_stream << std::string() +   "eeq_nan: FAILED"  + "\n";
};

void bin_functions_list::test_neq_nan()
{
    Real out = 0.;
    test_function_neq_nan tf;
    out += m_tests.make(&tf,m_options);

    if (out == 0.)
        matcl::out_stream << std::string() +   "neq_nan: OK" + "\n";
    else
        matcl::out_stream << std::string() +   "neq_nan: FAILED"  + "\n";
};

void bin_functions_list::test_op_lt()
{
    Real out = 0.;
    test_function_op_lt tf;
    out += m_tests.make(&tf,m_options);

    if (out == 0.)
        matcl::out_stream << std::string() +   "operator<: OK" + "\n";
    else
        matcl::out_stream << std::string() +   "operator<: FAILED"  + "\n";
};

void bin_functions_list::test_op_leq()
{
    Real out = 0.;
    test_function_op_leq tf;
    out += m_tests.make(&tf,m_options);

    if (out == 0.)
        matcl::out_stream << std::string() +   "operator<=: OK" + "\n";
    else
        matcl::out_stream << std::string() +   "operator<=: FAILED"  + "\n";
};

void bin_functions_list::test_op_gt()
{
    Real out = 0.;
    test_function_op_gt tf;
    out += m_tests.make(&tf,m_options);

    if (out == 0.)
        matcl::out_stream << std::string() +   "operator>: OK" + "\n";
    else
        matcl::out_stream << std::string() +   "operator>: FAILED"  + "\n";
};

void bin_functions_list::test_op_geq()
{
    Real out = 0.;
    test_function_op_geq tf;
    out += m_tests.make(&tf,m_options);

    if (out == 0.)
        matcl::out_stream << std::string() +   "operator>=: OK" + "\n";
    else
        matcl::out_stream << std::string() +   "operator>=: FAILED"  + "\n";
};

Real test_function_op_plus::eval_mat(const Matrix& mat1,const Matrix& mat2, int )
{
    if (mat1.is_scalar() == false || mat2.is_scalar() == false)
    {
        if (mat1.rows() != mat2.rows() || mat1.cols() != mat2.cols())
            return 0.0;    
    };

    try
    {		
        Matrix out_full = full(mat1) + full(mat2);	

        Matrix out		= mat1+mat2;
        check_struct(out);

        Real dif		= norm_1(out - out_full);
        dif				= dif + norm_1(mat2+mat1-out);        

        if (dif < constants::eps(out.get_value_code()) * norm_1(out) * 10.)
            dif = 0;

        dif				+= check_value_type(mat1,mat2,out,false);

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
Real test_function_op_plus::eval_scalar(const Scalar& s1, const Scalar& s2, int )
{
    Real out = eval_scalar_impl<test_function_op_plus>(s1, s2);
    return out;
};

Real test_function_op_minus::eval_mat(const Matrix& mat1,const Matrix& mat2, int )
{
    if (mat1.is_scalar() == false || mat2.is_scalar() == false)
    {
        if (mat1.rows() != mat2.rows() || mat1.cols() != mat2.cols())
            return 0.0;    
    };

    try
    {	
        Matrix out_full = full(mat1) - full(mat2);	
        Matrix out		= mat1 - mat2;
        check_struct(out);
        Real dif		= norm_1(out - out_full);

        Matrix out2		= mat2 - mat1;
        check_struct(out2);
        dif				+= norm_1(out - (-out2));

        if (dif < constants::eps(out.get_value_code()) * norm_1(out) * 10.)
            dif = 0;

        dif				+= check_value_type(mat1,mat2,out,false);
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
Real test_function_op_minus::eval_scalar(const Scalar& s1, const Scalar& s2, int )
{
    Real out = eval_scalar_impl<test_function_op_minus>(s1, s2);
    return out;
};

Real test_function_op_mult::eval_mat(const Matrix& mat1,const Matrix& mat2, int )
{
    if (mat1.is_scalar() == false || mat2.is_scalar() == false)
    {
        if (mat1.cols() != mat2.rows())
            return 0.0;    
    };

    Matrix out_full = full(mat1) * full(mat2);	

    try
    {		
        Matrix out		= mat1 * mat2;
        check_struct(out);

        Real dif		= norm_1(out - out_full);        
        Real nrm        = norm_1(out) + norm_1(mat1) * norm_1(mat2);
        Integer K       = mat1.cols();

        if (dif < constants::eps(out.get_value_code()) * nrm * K * 10.)
            dif = 0;

        dif				+= check_value_type(mat1,mat2,out,false);

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
Real test_function_op_mult::eval_scalar(const Scalar& s1, const Scalar& s2, int )
{
    Real out = eval_scalar_impl<test_function_op_mult>(s1, s2);
    return out;
};

Real test_function_mmul::eval_mat(const Matrix& mat10,const Matrix& mat20, int )
{
    if (mat10.is_scalar() == false || mat20.is_scalar() == false)
    {
        if (mat10.cols() != mat20.rows())
            return 0.0;    
    };

    Matrix mat1     = trans(mat10, m_ta);
    Matrix mat2     = trans(mat20, m_tb);

    Matrix out_full = trans(mat1, m_ta) * trans(mat2, m_tb);

    try
    {		
        Matrix out		= mmul(mat1, mat2, m_ta, m_tb);
        check_struct(out);

        Real dif		= norm_1(out - out_full);        
        Real nrm        = norm_1(out) + norm_1(mat1) * norm_1(mat2);
        Integer K       = mat1.cols();

        if (dif < constants::eps(out.get_value_code()) * nrm * K * 10.)
            dif = 0;

        dif				+= check_value_type(mat1,mat2,out,false);

        out             = Matrix();
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
Real test_function_mmul::eval_scalar(const Scalar& s1, const Scalar& s2, int )
{
    Real out = eval_scalar_impl<test_function_mmul>(s1, s2);
    return out;
};

Real test_function_mmul_ext::eval_mat(const Matrix& mat10,const Matrix& mat20, int )
{
    if (mat10.is_scalar() == false || mat20.is_scalar() == false)
    {
        if (mat10.cols() != mat20.rows())
            return 0.0;    
    };

    bool is_tr_1    = m_ta == trans_type_ext::trans || m_ta == trans_type_ext::conj_trans;
    bool is_tr_2    = m_tb == trans_type_ext::trans || m_tb == trans_type_ext::conj_trans;

    Matrix mat1     = is_tr_1 ? trans(mat10) : mat10;
    Matrix mat2     = is_tr_2 ? trans(mat20) : mat20;

    Matrix out_full = trans(mat1, m_ta) * trans(mat2, m_tb);

    try
    {		
        Matrix out		= mmul(mat1, mat2, m_ta, m_tb);
        check_struct(out);

        Real dif		= norm_1(out - out_full);        
        Real nrm        = norm_1(out) + norm_1(mat1) * norm_1(mat2);
        Integer K       = mat1.cols();

        if (dif < constants::eps(out.get_value_code()) * nrm * K * 10.)
            dif = 0;

        dif				+= check_value_type(mat1,mat2,out,false);

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
Real test_function_mmul_ext::eval_scalar(const Scalar& s1, const Scalar& s2, int )
{
    Real out = eval_scalar_impl<test_function_mmul_ext>(s1, s2);
    return out;
};

size_t test_function_gemm::rand_int(Integer seed)
{
    size_t val  = test_options::get_seed();
    boost::hash_combine(val, seed);
    return val;
};
Matrix test_function_gemm::rand_scalar(Integer seed)
{
    Integer type    = rand_int(seed) % 22;

    switch(type)
    {
        case 0:
            return Integer(0);
        case 1:
            return Integer(1);
        case 2:
            return Integer(-1);
        case 3:
            return Integer(3);

        case 4:
            return Float(0);
        case 5:
            return Float(1);
        case 6:
            return Float(-1);
        case 7:
            return Float(3);

        case 8:
            return Real(0);
        case 9:
            return Real(1);
        case 10:
            return Real(-1);
        case 11:
            return Real(3);

        case 12:
            return Float_complex(0);
        case 13:
            return Float_complex(1);
        case 14:
            return Float_complex(-1);
        case 15:
            return Float_complex(3);
        case 16:
            return Float_complex(3, 1);

        case 17:
            return Complex(0);
        case 18:
            return Complex(1);
        case 19:
            return Complex(-1);
        case 20:
            return Complex(3);
        case 21:
            return Complex(3, 1);
    }

    return 1.0;
};
value_code test_function_gemm::get_value_code(const Matrix& alpha, const Matrix& beta, const Matrix& mat1, 
                                const Matrix& mat2)
{
    value_code vc1  = matrix_traits::unify_value_types(alpha.get_value_code(), beta.get_value_code());
    value_code vc2  = matrix_traits::unify_value_types(mat1.get_value_code(), mat2.get_value_code());
    value_code vc   = matrix_traits::unify_value_types(vc1, vc2);

    return vc;
};

Matrix test_function_gemm::rand_mat_C(Integer r, Integer c, value_code vc_C, Integer code,
                                      dynamic_mat_set& ds)
{
    Integer seed= (Integer)rand_int(code);
    Matrix C    = ds.rand_dense(r, c, seed);

    value_code vc   = matrix_traits::unify_value_types(vc_C, C.get_value_code());
    C           = convert_value(C, vc);
    return C;
};

Real test_function_gemm::eval_mat(const Matrix& mat10,const Matrix& mat20, int code)
{
    if (mat10.is_scalar() == false || mat20.is_scalar() == false)
    {
        if (mat10.cols() != mat20.rows())
            return 0.0;    
    };

    Matrix mat1     = trans(mat10, m_ta);
    Matrix mat2     = trans(mat20, m_tb);    

    try
    {		
        Matrix alpha    = rand_scalar(code);
        Matrix beta     = rand_scalar(code + 1);
        value_code vc_C = get_value_code(alpha, beta, mat1, mat2);

        Integer r       = (m_ta == trans_type::no_trans) ? mat1.rows() : mat1.cols();
        Integer c       = (m_tb == trans_type::no_trans) ? mat2.cols() : mat2.rows();

        if (mat1.is_scalar() && mat2.is_scalar())
        {
            r           += code % 3;
            c           += code % 3;
        };

        Matrix C        = rand_mat_C(r, c, vc_C, code + 2, m_ds);
        vc_C            = C.get_value_code();

        Matrix alpha_c  = convert_value(alpha, vc_C);
        Matrix out_full = (alpha_c * trans(mat1, m_ta)) * trans(mat2, m_tb) + beta * C;
        
        Matrix out      = C;
        gemm(alpha, mat1, mat2, m_ta, m_tb, beta, out);
        check_struct(out);

        Real dif		= norm_1(out - out_full);        
        Real nrm        = norm_1(out) + norm_1(mat1) * norm_1(mat2) * norm_1(alpha);
        Integer K       = C.cols();

        if (dif < constants::eps(out.get_value_code()) * nrm * K * 10.)
            dif = 0;

        dif				+= check_value_type(out_full,out);

        out             = Matrix();
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
Real test_function_gemm::eval_scalar(const Scalar&, const Scalar&, int )
{
    return 0.0;
};

Integer test_function_gemm_sub::rand_colon_ver(Integer rs, Integer re, Integer cs, Integer ce,
                                               Integer seed)
{
    Integer val = test_function_gemm::rand_int(seed) % 4;

    if (rs == re)
    {
        if (cs == ce)
        {
            //possible version: 0, 1, 2, 3
            return val;
        }
        else
        {            
            //possible version: 1, 3
            return (val % 2) == 0 ? 1 : 3;
        }        
    }
    else
    {
        if (cs == ce)
        {
            //possible version: 2, 3
            return (val % 2) == 0 ? 2 : 3;
        }
        else
        {
            //possible version: 3
            return 3;
        }
    };
};

Real test_function_gemm_sub::eval_mat(const Matrix& mat10,const Matrix& mat20, int code)
{
    if (mat10.is_scalar() == false || mat20.is_scalar() == false)
    {
        if (mat10.cols() != mat20.rows())
            return 0.0;    
    };

    Matrix mat1     = trans(mat10, m_ta);
    Matrix mat2     = trans(mat20, m_tb);    

    try
    {		
        Matrix alpha    = test_function_gemm::rand_scalar(code);
        Matrix beta     = test_function_gemm::rand_scalar(code + 1);
        value_code vc_C = test_function_gemm::get_value_code(alpha, beta, mat1, mat2);

        Integer r       = (m_ta == trans_type::no_trans) ? mat1.rows() : mat1.cols();
        Integer c       = (m_tb == trans_type::no_trans) ? mat2.cols() : mat2.rows();

        if (mat1.is_scalar() && mat2.is_scalar())
        {
            r           += code % 3;
            c           += code % 3;
        };

        Integer r2      = r + 2 + (code % 5);
        Integer c2      = c + 2 + (code % 5);

        Integer rs      = 1 + (code % 2);
        Integer cs      = 1 + (code % 3);

        Integer re      = rs + r - 1;
        Integer ce      = cs + c - 1;

        Matrix C        = test_function_gemm::rand_mat_C(r2, c2, vc_C, code + 2, m_ds);
        vc_C            = C.get_value_code();

        colon cr(rs, re);
        colon cc(cs, ce);

        Matrix alpha_c  = convert_value(alpha, vc_C);
        Matrix out_full = (alpha_c * trans(mat1, m_ta)) * trans(mat2, m_tb) + beta * C(cr, cc);
        
        Matrix out      = C;
        Integer ver     = rand_colon_ver(rs, re, cs, ce, code + 3);

        switch (ver)
        {
            case 0:
                gemm(alpha, mat1, mat2, m_ta, m_tb, beta, out(rs, cs));
                break;
            case 1:
                gemm(alpha, mat1, mat2, m_ta, m_tb, beta, out(rs, cc));
                break;
            case 2:
                gemm(alpha, mat1, mat2, m_ta, m_tb, beta, out(cr, cs));
                break;
            case 3:        
                gemm(alpha, mat1, mat2, m_ta, m_tb, beta, out(cr, cc));
                break;
        };

        check_struct(out);

        Matrix out2     = out(cr, cc);
        Real dif		= norm_1(out2 - out_full);        
        Real nrm        = norm_1(out2) + norm_1(mat1) * norm_1(mat2) * norm_1(alpha);
        Integer K       = C.cols();

        if (dif < constants::eps(out2.get_value_code()) * nrm * K * 10.)
            dif = 0;

        dif				+= check_value_type(out_full,out2);

        out             = Matrix();
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
Real test_function_gemm_sub::eval_scalar(const Scalar&, const Scalar&, int )
{
    return 0.0;
};

Matrix test_function_mmul_abs::rand_mat_C(Integer r, Integer c, Integer code, dynamic_mat_set& ds)
{
    Integer seed= (Integer)test_function_gemm::rand_int(code);
    Matrix C    = ds.rand_dense(r, c, seed);
    return C;
};

Real test_function_mmul_abs::eval_mat(const Matrix& mat10,const Matrix& mat20, int code)
{
    if (mat10.is_scalar() == false || mat20.is_scalar() == false)
    {
        if (mat10.cols() != mat20.rows())
            return 0.0;    
    };

    Matrix mat1     = trans(mat10, m_ta);
    Matrix mat2     = mat20;

    try
    {		
        Integer r       = (m_ta == trans_type::no_trans) ? mat1.rows() : mat1.cols();
        Integer c       = mat2.cols();

        if (mat1.is_scalar() && mat2.is_scalar())
        {
            r           += code % 3;
            c           += code % 3;
        };

        Matrix C        = rand_mat_C(r, c, code + 2, m_ds);
        value_code vc_C = test_function_gemm::get_value_code(mat1, mat1, mat2, C);

        Matrix C2       = convert_value(C, vc_C);
        Matrix mat22    = convert_value(mat2, vc_C);
        Matrix mat12    = convert_value(mat1, vc_C);
        Matrix out_full = abs(trans(mat12, m_ta)) * abs(mat22) + abs(C2);        

        Matrix out      = mmul_abs(mat1, mat2, m_ta, C);        
        check_struct(out);

        Real dif		= norm_1(out - out_full);        
        Real nrm        = norm_1(out) + norm_1(mat1) * norm_1(mat2);
        Integer K       = C.cols();

        if (dif < constants::eps(out.get_value_code()) * nrm * K * 10.)
            dif = 0;

        dif				+= check_value_type(out_full,out);

        out             = Matrix();
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
Real test_function_mmul_abs::eval_scalar(const Scalar&, const Scalar&, int )
{
    return 0.0;
};

Real test_function_chain_mult::eval_mat(const Matrix& mat0,const Matrix& mat2, int )
{
    if (mat0.is_scalar() == false || mat2.is_scalar() == false)
    {
        if (mat0.cols() != mat2.rows())
            return 0.0;    
    };

    try
    {
        Matrix mat1     = mat0 + 0.f;
        int dim_approx  = mat1.rows() * mat1.cols() * mat2.rows() * mat2.cols();
        (void)dim_approx;

        Matrix out_full_1 = full(mat1) * full(mat2);
        Matrix out_full_2 = out_full_1 * full(trans(mat2)) * full(trans(mat1));
        Matrix out_full_3 = out_full_2 * full(mat1) * full(mat2);

        Matrix out_1	= chain_mult(mat1,mat2);
        Matrix out_2	= chain_mult(mat1,mat2,trans(mat2),trans(mat1));
        Matrix out_3	= chain_mult(mat1,mat2,trans(mat2),trans(mat1),mat1,mat2);

        check_struct(out_1);
        check_struct(out_2);
        check_struct(out_3);

        bool is_sing    = matrix_traits::is_single_precision(mat1.get_value_code())
                        ||matrix_traits::is_single_precision(mat2.get_value_code());

        Real eps        = is_sing ? constants::f_eps() : constants::eps();

        Real dif1		= norm_1(out_1 - out_full_1);
        Real nr_1       = norm_1(mat1);
        Real nr_2       = norm_1(mat2);
        Real nr         = nr_1 + nr_2 + nr_1 * nr_2;
        Real nrm1       = norm_1(out_1) + nr; 

        if (dif1 < eps * nrm1 * 100.)
            dif1        = 0;

        Real dif2		= norm_1(out_2 - out_full_2);
        Real nrm2       = norm_1(out_2) + nr * nr; 
        if (dif2 < eps * nrm2 * 100.)
            dif2        = 0;

        Real dif3       = norm_1(out_3 - out_full_3);
        Real nrm3       = norm_1(out_3) + nr * nr * nr; 
        if (dif3 < eps * nrm3 * 100.)
            dif3        = 0;

        Real dif        = dif1 + dif2 + dif3;
        dif				+= check_value_type(out_1,out_full_1);
        dif		    	+= check_value_type(out_2,out_full_2);
        dif		    	+= check_value_type(out_3,out_full_3);

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

Real test_function_chain_mult::eval_scalar(const Scalar& , const Scalar& , int )
{
    return 0;
};

Real test_function_op_div::eval_mat(const Matrix& mat1,const Matrix& mat2, int )
{
    if (mat1.is_scalar() == false || mat2.is_scalar() == false)
    {
        if (mat1.rows() != mat2.rows() || mat1.cols() != mat2.cols())
            return 0.0;    
    };

    Matrix out_full = full(mat1) / full(mat2);	

    try
    {		
        Matrix out		= mat1 / mat2;
        check_struct(out);
        
        Matrix out2     = div(mat1, mat2);

        Real dif		= norm_1(out - out_full);        
        dif             += norm_1(out - out2);

        if (dif < constants::eps(out.get_value_code()) * norm_1(out) * 10.)
            dif         = 0;

        dif				+= check_value_type(mat1,mat2,out,true);
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
Real test_function_op_div::eval_scalar(const Scalar& s1, const Scalar& s2, int )
{	
    Real out = eval_scalar_impl<test_function_op_div>(s1, s2);
    return out;
};

Real test_function_plus::eval_mat(const Matrix& mat1,const Matrix& mat2, int )
{
    if (mat1.is_scalar() == false || mat2.is_scalar() == false)
    {
        if (mat1.rows() != mat2.rows() || mat1.cols() != mat2.cols())
            return 0.0;    
    }; 

    try
    {		
        Matrix out_full = plus(full(mat1),full(mat2));
        Matrix out		= plus(mat1,mat2);
        check_struct(out);

        Real dif		= norm_1(out - out_full);
        dif				+=norm_1(out - (mat1+mat2));

        if (dif < constants::eps(out.get_value_code()) * norm_1(out) * 10.)
            dif         = 0;

        dif				+= check_value_type(mat1,mat2,out,false);
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
Real test_function_plus::eval_scalar(const Scalar& s1, const Scalar& s2, int )
{
    Real out = eval_scalar_impl<test_function_plus>(s1, s2);
    return out;
};

Real test_function_minus::eval_mat(const Matrix& mat1,const Matrix& mat2, int )
{
    if (mat1.is_scalar() == false || mat2.is_scalar() == false)
    {
        if (mat1.rows() != mat2.rows() || mat1.cols() != mat2.cols())
            return 0.0;    
    };
    
    try
    {		
        Matrix out_full = minus(full(mat1),full(mat2));
        Matrix out		= minus(mat1,mat2);
        check_struct(out);

        Real dif		= norm_1(out - out_full);
        dif				+=norm_1(out - (mat1-mat2));

        if (dif < constants::eps(out.get_value_code()) * norm_1(out) * 10.)
            dif         = 0;

        dif				+= check_value_type(mat1,mat2,out,false);
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
Real test_function_minus::eval_scalar(const Scalar& s1, const Scalar& s2, int )
{
    Real out = eval_scalar_impl<test_function_minus>(s1, s2);
    return out;
};

Real test_function_mul::eval_mat(const Matrix& mat1,const Matrix& mat2, int )
{
    if (mat1.is_scalar() == false || mat2.is_scalar() == false)
    {
        if (mat1.rows() != mat2.rows() || mat1.cols() != mat2.cols())
            return 0.0;    
    };   

    try
    {
        Matrix out_full = mul(full(mat1),full(mat2));
        Matrix out		= mul(mat1,mat2);
        check_struct(out);

        Real dif 		= norm_1(out - out_full);
        dif             += norm_1(out - mul(mat2,mat1));

        if (dif < constants::eps(out.get_value_code()) * norm_1(out) * 10.)
            dif         = 0;

        dif             += check_value_type(mat1,mat2,out,false);
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
Real test_function_mul::eval_scalar(const Scalar& s1, const Scalar& s2, int )
{
    Real out = eval_scalar_impl<test_function_mul>(s1, s2);
    return out;
};

Real test_function_div::eval_mat(const Matrix& mat1,const Matrix& mat2, int )
{
    if (mat1.is_scalar() == false || mat2.is_scalar() == false)
    {
        if (mat1.rows() != mat2.rows() || mat1.cols() != mat2.cols())
            return 0.0;    
    };

    try
    {		
        Matrix out_full = div(full(mat1),full(mat2));
        Matrix out		= div(mat1,mat2);
        check_struct(out);

        Real dif		= norm_1(out - out_full);
        if (dif < constants::eps(out.get_value_code()) * norm_1(out) * 10.)
            dif         = 0;

        dif				+= check_value_type(mat1,mat2,out,true);
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
Real test_function_div::eval_scalar(const Scalar& s1, const Scalar& s2, int )
{
    Real out = eval_scalar_impl<test_function_div>(s1, s2);
    return out;
};

Real test_function_div_0::eval_mat(const Matrix& mat1,const Matrix& mat2, int )
{
    if (mat1.is_scalar() == false || mat2.is_scalar() == false)
    {
        if (mat1.rows() != mat2.rows() || mat1.cols() != mat2.cols())
            return 0.0;    
    };   

    try
    {		
        Matrix out_full = div_0(full(mat1),full(mat2));
        Matrix out		= div_0(mat1,mat2);
        check_struct(out);

        Real dif		= norm_1(out - out_full);
        if (dif < constants::eps(out.get_value_code()) * norm_1(out) * 10.)
            dif         = 0;

        dif				+= check_value_type(mat1,mat2,out,true);
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
Real test_function_div_0::eval_scalar(const Scalar& s1, const Scalar& s2, int )
{
    Real out = eval_scalar_impl<test_function_div_0>(s1, s2);
    return out;
};

Real test_function_div_1::eval_mat(const Matrix& mat1,const Matrix& mat2, int )
{
    if (mat1.is_scalar() == false || mat2.is_scalar() == false)
    {
        if (mat1.rows() != mat2.rows() || mat1.cols() != mat2.cols())
            return 0.0;    
    };       

    try
    {		
        Matrix out_full = div_1(full(mat1),full(mat2));
        Matrix out		= div_1(mat1,mat2);
        check_struct(out);

        Real dif		= norm_1(out - out_full);
        if (dif < constants::eps(out.get_value_code()) * norm_1(out) * 10.)
            dif         = 0;

        dif				+= check_value_type(mat1,mat2,out,true);
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
Real test_function_div_1::eval_scalar(const Scalar& s1, const Scalar& s2, int )
{
    Real out = eval_scalar_impl<test_function_div_1>(s1, s2);
    return out;
};

Real test_function_idiv::eval_mat(const Matrix& mat1,const Matrix& mat2, int )
{
    if (mat1.is_scalar() == false || mat2.is_scalar() == false)
    {
        if (mat1.rows() != mat2.rows() || mat1.cols() != mat2.cols())
            return 0.0;    
    };          

    try
    {		
        Matrix out_full = idiv(full(mat1),full(mat2));
        Matrix out		= idiv(mat1,mat2);

        check_struct(out);

        Real dif		= norm_1(out - out_full);
        if (dif < constants::eps(out.get_value_code()) * norm_1(out) * 10.)
            dif         = 0;

        dif				+= check_value_type(mat1,mat2,out,false);
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
Real test_function_idiv::eval_scalar(const Scalar& s1, const Scalar& s2, int )
{
    Real out = eval_scalar_impl<test_function_idiv>(s1, s2);
    return out;
};

Real test_function_pow::eval_mat(const Matrix& mat1,const Matrix& mat2, int )
{
    if (mat1.is_scalar() == false || mat2.is_scalar() == false)
    {
        if (mat1.rows() != mat2.rows() || mat1.cols() != mat2.cols())
            return 0.0;    
    };       

    try
    {		
        Matrix out_full = pow(full(mat1),full(mat2));
        Matrix out		= pow(mat1,mat2);
        check_struct(out);

        Real dif		= norm_1(out - out_full);

        if (dif < constants::eps(out.get_value_code()) * norm_1(out) * 10.)
            dif         = 0;

        dif				+= test_function_pow_c::check_value_type_pow(mat1,mat2,out);
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
Real test_function_pow::eval_scalar(const Scalar& s1, const Scalar& s2, int )
{
    Real out = eval_scalar_impl<test_function_pow>(s1, s2);
    return out;
};

Real test_function_pow_c::eval_mat(const Matrix& mat1,const Matrix& mat2, int )
{
    if (mat1.is_scalar() == false || mat2.is_scalar() == false)
    {
        if (mat1.rows() != mat2.rows() || mat1.cols() != mat2.cols())
            return 0.0;    
    };          
    
    try
    {		
        Matrix out_full = pow_c(full(mat1),full(mat2));
        Matrix out		= pow_c(mat1,mat2);
        check_struct(out);

        Real dif		= norm_1(out - out_full);

        if (dif < constants::eps(out.get_value_code()) * norm_1(out) * 10.)
            dif         = 0;

        Matrix mat1_c   = mat1 + Float_complex(0.0f, 0.0f);
        Matrix out2     = pow(mat1_c, mat2);
        dif				+= check_value_type_pow(mat1,mat2,out);
        dif				+= compare_value_type(out,out2);

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
Real test_function_pow_c::eval_scalar(const Scalar& s1, const Scalar& s2, int )
{
    Real out = eval_scalar_impl<test_function_pow_c>(s1, s2);
    return out;
};

Real test_function_pow_c::check_value_type_pow(const Matrix& mat1, const Matrix& mat2, const Matrix& out)
{
    value_code v1   = mat1.get_value_code();
    value_code v2   = mat2.get_value_code();
    value_code v3   = out.get_value_code();

    value_code v;

    if (v2  == value_code::v_integer)
    {
        v           = matrix_traits::unify_value_types(v1, value_code::v_float);
    }
    else
    {
        v           = matrix_traits::unify_value_types(v1, v2);
        v           = matrix_traits::unify_value_types(v, value_code::v_float);
    };

    v               = matrix_traits::real_value_type(v);
    v3              = matrix_traits::real_value_type(v3);

    if (v == v3)
        return 0.0;
    else 
        return 1.0;

};
Real test_function_pow_c::compare_value_type(const Matrix& mat1, const Matrix& mat2)
{
    value_code v1   = mat1.get_value_code();
    value_code v2   = mat2.get_value_code();

    v1              = matrix_traits::real_value_type(v1);
    v2              = matrix_traits::real_value_type(v2);

    if (v1 == v2)
        return 0.0;
    else 
        return 1.0;
};

Real test_function_max_bin::eval_mat(const Matrix& mat1,const Matrix& mat2, int )
{
    if (mat1.is_scalar() == false || mat2.is_scalar() == false)
    {
        if (mat1.rows() != mat2.rows() || mat1.cols() != mat2.cols())
            return 0.0;    
    };          
    
    try
    {		
        Matrix out_full =  max(full(mat1),full(mat2));	
        Matrix out		= max(mat1,mat2);
        check_struct(out);

        Real dif		= norm_1(out - out_full);
        dif				+= norm_1(out - max(mat2,mat1));
        dif				+= check_value_type(mat1,mat2,out,false);
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
Real test_function_max_bin::eval_scalar(const Scalar& s1, const Scalar& s2, int )
{
    Real out = eval_scalar_impl<test_function_max_bin>(s1, s2);
    return out;
};

Real test_function_min_bin::eval_mat(const Matrix& mat1,const Matrix& mat2, int )
{
    if (mat1.is_scalar() == false || mat2.is_scalar() == false)
    {
        if (mat1.rows() != mat2.rows() || mat1.cols() != mat2.cols())
            return 0.0;    
    };               

    try
    {		
        Matrix out_full = min(full(mat1),full(mat2));
        Matrix out		= min(mat1,mat2);
        check_struct(out);

        Real dif		= norm_1(out - out_full);
        dif				+= norm_1(out - min(mat2,mat1));
        dif				+= check_value_type(mat1,mat2,out,false);
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
Real test_function_min_bin::eval_scalar(const Scalar& s1, const Scalar& s2, int )
{
    Real out = eval_scalar_impl<test_function_min_bin>(s1, s2);
    return out;
};

Real test_function_xor::eval_mat(const Matrix& mat1,const Matrix& mat2, int )
{
    if (mat1.is_scalar() == false || mat2.is_scalar() == false)
    {
        if (mat1.rows() != mat2.rows() || mat1.cols() != mat2.cols())
            return 0.0;    
    };

    try
    {		
        Matrix out_full = elem_xor(full(mat1),full(mat2));	
        Matrix out		= elem_xor(mat1,mat2);
        check_struct(out);

        Matrix out2     = mat1 ^ mat2;

        Real dif		= norm_1(out - out_full);
        dif				+= norm_1(out - elem_xor(mat2,mat1));
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
Real test_function_xor::eval_scalar(const Scalar& s1, const Scalar& s2, int )
{
    Real out = eval_scalar_impl<test_function_xor>(s1, s2);
    return out;
};

Real test_function_rem::eval_mat(const Matrix& mat1,const Matrix& mat2, int )
{
    if (mat1.is_scalar() == false || mat2.is_scalar() == false)
    {
        if (mat1.rows() != mat2.rows() || mat1.cols() != mat2.cols())
            return 0.0;    
    };   

    value_code vc1  = mat1.get_value_code();
    value_code vc2  = mat2.get_value_code();

    if (matrix_traits::is_float_complex(vc1) == true || matrix_traits::is_float_complex(vc2) == true)
    {
        return 0.0;
    }

    try
    {		
        Matrix out_full = rem(full(mat1),full(mat2));	
        Matrix out		= rem(mat1,mat2);
        check_struct(out);

        Real dif		= norm_1(out - out_full);
        dif				+= check_value_type(mat1,mat2,out,false);
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
Real test_function_rem::eval_scalar(const Scalar& s1, const Scalar& s2, int )
{
    Real out = eval_scalar_impl<test_function_rem>(s1, s2);
    return out;
};

Real test_function_mod::eval_mat(const Matrix& mat1,const Matrix& mat2, int )
{
    if (mat1.is_scalar() == false || mat2.is_scalar() == false)
    {
        if (mat1.rows() != mat2.rows() || mat1.cols() != mat2.cols())
            return 0.0;    
    };       

    value_code vc1  = mat1.get_value_code();
    value_code vc2  = mat2.get_value_code();

    if (matrix_traits::is_float_complex(vc1) == true || matrix_traits::is_float_complex(vc2) == true)
    {
        return 0.0;
    }

    try
    {		
        Matrix out_full = mod(full(mat1),full(mat2));
        Matrix out		= mod(mat1,mat2);
        check_struct(out);

        Real dif		= norm_1(out - out_full);
        dif				+= check_value_type(mat1,mat2,out,false);
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
Real test_function_mod::eval_scalar(const Scalar& s1, const Scalar& s2, int )
{
    Real out = eval_scalar_impl<test_function_mod>(s1, s2);
    return out;
};

Real test_function_kron::eval_mat(const Matrix& mat1,const Matrix& mat2, int )
{
    Matrix out_full = kron(full(mat1),full(mat2));	

    try
    {		
        Matrix out		= kron(mat1,mat2);
        check_struct(out);

        Real dif		= norm_1(out - out_full);
        dif				+= check_value_type(mat1,mat2,out,false);

        if (abs(dif) < norm_1(out) * constants::eps(out.get_value_code()) * 10.0)
            dif         = 0.0;

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
Real test_function_kron::eval_scalar(const Scalar& s1, const Scalar& s2, int )
{
    Real out = eval_scalar_impl<test_function_kron>(s1, s2);
    return out;
};

Real test_function_op_or::eval_mat(const Matrix& mat1,const Matrix& mat2, int )
{
    if (mat1.is_scalar() == false || mat2.is_scalar() == false)
    {
        if (mat1.rows() != mat2.rows() || mat1.cols() != mat2.cols())
            return 0.0;    
    };           

    try
    {		
        Matrix out_full =(full(mat1)|full(mat2));	
        Matrix out		= (mat1|mat2);
        check_struct(out);

        Matrix out2     = elem_or(mat1, mat2);

        Real dif		= norm_1(out - out_full);
        dif				+= norm_1(out - (mat2|mat1));
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
Real test_function_op_or::eval_scalar(const Scalar& s1, const Scalar& s2, int )
{
    Real out = eval_scalar_impl<test_function_op_or>(s1, s2);
    return out;
};

Real test_function_op_and::eval_mat(const Matrix& mat1,const Matrix& mat2, int )
{
    if (mat1.is_scalar() == false || mat2.is_scalar() == false)
    {
        if (mat1.rows() != mat2.rows() || mat1.cols() != mat2.cols())
            return 0.0;    
    };            

    try
    {		
        Matrix out_full = (full(mat1)&full(mat2));	
        Matrix out		= (mat1&mat2);
        check_struct(out);

        Matrix out2     = elem_and(mat1, mat2);

        Real dif		= norm_1(out - out_full);
        dif				+= norm_1(out -(mat2&mat1));
        dif             += norm_1(out - out2);

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
Real test_function_op_and::eval_scalar(const Scalar& s1, const Scalar& s2, int )
{
    Real out = eval_scalar_impl<test_function_op_and>(s1, s2);
    return out;
};

Real test_function_op_eeq::eval_mat(const Matrix& mat1,const Matrix& mat2, int )
{
    if (mat1.is_scalar() == false || mat2.is_scalar() == false)
    {
        if (mat1.rows() != mat2.rows() || mat1.cols() != mat2.cols())
            return 0.0;    
    };
    
    try
    {		
        Matrix out_full = (full(mat1)==full(mat2));	
        Matrix out		= (mat1 == mat2);
        check_struct(out);

        Matrix out2     = eeq(mat1, mat2);

        Real dif		= norm_1(out - out_full);
        dif				+= norm_1(out - (mat2==mat1));
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
Real test_function_op_eeq::eval_scalar(const Scalar& s1, const Scalar& s2, int )
{
    Real out = eval_scalar_impl<test_function_op_eeq>(s1, s2);
    return out;
};

Real test_function_op_neq::eval_mat(const Matrix& mat1,const Matrix& mat2, int )
{
    if (mat1.is_scalar() == false || mat2.is_scalar() == false)
    {
        if (mat1.rows() != mat2.rows() || mat1.cols() != mat2.cols())
            return 0.0;    
    };

    try
    {		
        Matrix out_full = (full(mat1)!=full(mat2));	
        Matrix out		= (mat1 != mat2);
        check_struct(out);

        Matrix out2     = neq(mat1, mat2);

        Real dif		= norm_1(out - out_full);
        dif				+= norm_1(out - (mat2!=mat1));
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
Real test_function_op_neq::eval_scalar(const Scalar& s1, const Scalar& s2, int )
{
    Real out = eval_scalar_impl<test_function_op_neq>(s1, s2);
    return out;
};

Real test_function_eeq_nan::eval_mat(const Matrix& mat1,const Matrix& mat2, int )
{
    if (mat1.is_scalar() == false || mat2.is_scalar() == false)
    {
        if (mat1.rows() != mat2.rows() || mat1.cols() != mat2.cols())
            return 0.0;    
    };

    try
    {		
        Matrix out_full = (eeq_nan(full(mat1), full(mat2)));	
        Matrix out		= eeq_nan(mat1, mat2);
        check_struct(out);

        Real dif		= norm_1(out - out_full);

        Matrix I1       = is_nan(mat1);
        Matrix I2       = is_nan(mat2);
        Matrix I12      = mul(I1, I2);

        Matrix I1_a     = is_nan(real(mat1)) & is_nan(imag(mat1));
        Matrix I2_a     = is_nan(real(mat2)) & is_nan(imag(mat2));
        Matrix I12_a    = mul(I1_a, I2_a);

        Matrix sel_fin  = find(I12 == 0);
        Matrix sel_nan2 = find(I12_a == 1);
        Matrix sel_nan1 = find((I1 == 1 | I2 == 1) & (I12 == 0));

        Matrix out2     = (mat2 == mat1);
        dif				+= norm_1(out(sel_fin) - out2(sel_fin));
        dif				+= norm_1(out(sel_nan2) - 1);
        dif				+= norm_1(out(sel_nan1) - 0);

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
Real test_function_eeq_nan::eval_scalar(const Scalar& s1, const Scalar& s2, int )
{
    Real out = eval_scalar_impl<test_function_eeq_nan>(s1, s2);
    return out;
};

Real test_function_neq_nan::eval_mat(const Matrix& mat1,const Matrix& mat2, int )
{
    if (mat1.is_scalar() == false || mat2.is_scalar() == false)
    {
        if (mat1.rows() != mat2.rows() || mat1.cols() != mat2.cols())
            return 0.0;    
    };

    try
    {		
        Matrix out_full = (neq_nan(full(mat1), full(mat2)));	
        Matrix out		= neq_nan(mat1, mat2);
        check_struct(out);

        Real dif		= norm_1(out - out_full);

        Matrix I1       = is_nan(mat1);
        Matrix I2       = is_nan(mat2);
        Matrix I12      = mul(I1, I2);

        Matrix I1_a     = is_nan(real(mat1)) & is_nan(imag(mat1));
        Matrix I2_a     = is_nan(real(mat2)) & is_nan(imag(mat2));
        Matrix I12_a    = mul(I1_a, I2_a);

        Matrix sel_fin  = find(I12 == 0);
        Matrix sel_nan2 = find(I12_a == 1);
        Matrix sel_nan1 = find((I1 == 1 | I2 == 1) & (I12 == 0));

        Matrix out2     = (mat2 != mat1);
        dif				+= norm_1(out(sel_fin) - out2(sel_fin));

        dif				+= norm_1(out(sel_nan2) - 0);
        dif				+= norm_1(out(sel_nan1) - 1);

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
Real test_function_neq_nan::eval_scalar(const Scalar& s1, const Scalar& s2, int )
{
    Real out = eval_scalar_impl<test_function_neq_nan>(s1, s2);
    return out;
};

Real test_function_op_lt::eval_mat(const Matrix& mat1,const Matrix& mat2, int )
{
    if (mat1.is_scalar() == false || mat2.is_scalar() == false)
    {
        if (mat1.rows() != mat2.rows() || mat1.cols() != mat2.cols())
            return 0.0;    
    };
    
    try
    {		
        Matrix out_full = (full(mat1) < full(mat2));	
        Matrix out		= (mat1 < mat2);
        check_struct(out);

        Matrix out2     = lt(mat1, mat2);

        Real dif1		= norm_1(out - out_full);
        Matrix non_nan_elements = 1 - is_nan(mat1 + mat2);
        Real dif2       = norm_1(mul(out - ~(mat1>=mat2), non_nan_elements));
        Real dif3       = norm_1(mul(out - ~(mat2<=mat1), non_nan_elements));
        Real dif4       = norm_1(out - out2);

        Real dif        = dif1 + dif2 + dif3 + dif4;

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
Real test_function_op_lt::eval_scalar(const Scalar& s1, const Scalar& s2, int )
{
    Real out = eval_scalar_impl<test_function_op_lt>(s1, s2);
    return out;
};

Real test_function_op_leq::eval_mat(const Matrix& mat1,const Matrix& mat2, int )
{
    if (mat1.is_scalar() == false || mat2.is_scalar() == false)
    {
        if (mat1.rows() != mat2.rows() || mat1.cols() != mat2.cols())
            return 0.0;    
    };   

    try
    {		
        Matrix out_full =  (full(mat1) <= full(mat2));	
        Matrix out		= (mat1 <= mat2);
        check_struct(out);

        Matrix out2     = leq(mat1, mat2);

        Real dif1		= norm_1(out - out_full);
        Matrix non_nan_elements = 1 - is_nan(mat1 + mat2);
        Real dif2       = norm_1(mul(out - ~(mat1>mat2), non_nan_elements));
        Real dif3       = norm_1(mul(out - ~(mat2<mat1), non_nan_elements));
        Real dif4       = norm_1(out - out2);

        Real dif        = dif1 + dif2 + dif3 + dif4;
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
Real test_function_op_leq::eval_scalar(const Scalar& s1, const Scalar& s2, int )
{
    Real out = eval_scalar_impl<test_function_op_leq>(s1, s2);
    return out;
};

Real test_function_op_gt::eval_mat(const Matrix& mat1,const Matrix& mat2, int )
{
    if (mat1.is_scalar() == false || mat2.is_scalar() == false)
    {
        if (mat1.rows() != mat2.rows() || mat1.cols() != mat2.cols())
            return 0.0;    
    };     

    try
    {		
        Matrix out_full =  (full(mat1) > full(mat2));
        Matrix out		= (mat1 > mat2);
        check_struct(out);

        Matrix out2     = gt(mat1, mat2);

        Real dif1		= norm_1(out - out_full);
        Real dif2       = norm_1(out - (mat2<mat1));
        Matrix non_nan_elements = 1 - is_nan(mat1 + mat2);
        Real dif3       = norm_1(mul(out - ~(mat1<=mat2), non_nan_elements));
        Real dif4       = norm_1(mul(out - ~(mat2>=mat1), non_nan_elements));
        Real dif5       = norm_1(out - out2);

        Real dif        = dif1 + dif2 + dif3 + dif4 + dif5;

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
Real test_function_op_gt::eval_scalar(const Scalar& s1, const Scalar& s2, int )
{
    Real out = eval_scalar_impl<test_function_op_gt>(s1, s2);
    return out;
};

Real test_function_op_geq::eval_mat(const Matrix& mat1,const Matrix& mat2, int )
{
    if (mat1.is_scalar() == false || mat2.is_scalar() == false)
    {
        if (mat1.rows() != mat2.rows() || mat1.cols() != mat2.cols())
            return 0.0;    
    };           

    try
    {		
        Matrix out_full = (full(mat1) >= full(mat2));	
        Matrix out		= (mat1 >= mat2);
        check_struct(out);

        Matrix out2     = geq(mat1, mat2);

        Real dif1		= norm_1(out - out_full);
        Real dif2       = norm_1(out - (mat2<=mat1));
        Matrix non_nan_elements = 1 - is_nan(mat1 + mat2);
        Real dif3       = norm_1(mul(out - ~(mat1<mat2), non_nan_elements));
        Real dif4       = norm_1(mul(out - ~(mat2>mat1), non_nan_elements));
        Real dif5       = norm_1(out - out2);

        Real dif        = dif1 + dif2 + dif3 + dif4 + dif5;

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
Real test_function_op_geq::eval_scalar(const Scalar& s1, const Scalar& s2, int )
{
    Real out = eval_scalar_impl<test_function_op_geq>(s1, s2);
    return out;
};

Real test_function_beta::eval_mat(const Matrix& mat1,const Matrix& mat2, int )
{
    if (mat1.is_scalar() == false || mat2.is_scalar() == false)
    {
        if (mat1.rows() != mat2.rows() || mat1.cols() != mat2.cols())
            return 0.0;    
    };

    try
    {
        Matrix out_full;
        out_full    = beta(full(mat1),full(mat2)); 
        Matrix out	= beta(mat2,mat1); // Note reversed order of args. Beta should be symmetric.
        check_struct(out);

        return norm_1(out - out_full);
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

Real test_function_beta::eval_scalar(const Scalar& s1, const Scalar& s2, int )
{
    Real out = eval_scalar_impl<test_function_beta>(s1, s2);
    return out;
}

Real test_function_eval_bin::eval_mat(const Matrix& mat1,const Matrix& mat2, int )
{
    if (mat1.is_scalar() == false || mat2.is_scalar() == false)
    {
        if (mat1.rows() != mat2.rows() || mat1.cols() != mat2.cols())
            return 0.0;    
    };

    try
    {
        Matrix out1 = eval_binary_func(mat1, mat2, binfunc_plus());
        Matrix out3 = eval_binary_func(mat1, mat2, binfunc_plus(), binfunc_plus());
        Matrix out2 = plus(mat1, mat2); 
        check_struct(out1);

        Real dif1   = norm_1(out1 - out2);
        dif1        += norm_1(out1 - out3);
        Real tol0   = norm_1(mat1) + norm_1(mat2);
        Real tol1   = tol0 + norm_1(out1);

        if (dif1 < tol1 * constants::eps(out1.get_value_code()) * 10.0)
            dif1    = 0;

        Real dif    = dif1;
        dif         += check_value_type(out1, out2);

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

Real test_function_eval_bin::eval_scalar(const Scalar&, const Scalar&, int )
{
    return 0.0;
}


Real test_function_atan2::eval_mat(const Matrix& mat1,const Matrix& mat2, int )
{
    if (mat1.is_scalar() == false || mat2.is_scalar() == false)
    {
        if (mat1.rows() != mat2.rows() || mat1.cols() != mat2.cols())
            return 0.0;    
    };

    value_code vc1  = mat1.get_value_code();
    value_code vc2  = mat2.get_value_code();

    if (matrix_traits::is_float_complex(vc1) == true || matrix_traits::is_float_complex(vc2) == true)
    {
        return 0.0;
    }

    try
    {		
        Matrix out_full;
        out_full        = atan2(full(mat1),full(mat2));
        Matrix out		= atan2(mat1,mat2);
        check_struct(out);

        Real dif		= norm_1(out - out_full);
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
Real test_function_atan2::eval_scalar(const Scalar& s1, const Scalar& s2, int )
{
    Real out = eval_scalar_impl<test_function_atan2>(s1, s2);
    return out;
};

Real test_function_hypot::eval_mat(const Matrix& mat1,const Matrix& mat2, int )
{
    if (mat1.is_scalar() == false || mat2.is_scalar() == false)
    {
        if (mat1.rows() != mat2.rows() || mat1.cols() != mat2.cols())
            return 0.0;    
    };

    try
    {		
        Matrix out_full;
        out_full        = hypot(full(mat1),full(mat2));
        Matrix out		= hypot(mat1,mat2);
        check_struct(out);

        Real dif		= norm_1(out - out_full);
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
Real test_function_hypot::eval_scalar(const Scalar& s1, const Scalar& s2, int )
{
    Real out = eval_scalar_impl<test_function_hypot>(s1, s2);
    return out;
};

Real test_function_copysign::eval_mat(const Matrix&,const Matrix&, int )
{
    return 0.0;
};
Real test_function_copysign::eval_scalar(const Scalar& s1, const Scalar& s2, int )
{
    Real out = eval_scalar_impl<test_function_copysign>(s1, s2);
    return out;
};

Real test_function_fdim::eval_mat(const Matrix&,const Matrix&, int )
{
    return 0.0;
};
Real test_function_fdim::eval_scalar(const Scalar& s1, const Scalar& s2, int )
{
    Real out = eval_scalar_impl<test_function_fdim>(s1, s2);
    return out;
};

Real test_function_fma::eval_mat(const Matrix&,const Matrix&, int )
{
    return 0.0;
};
Real test_function_fma::eval_scalar(const Scalar& s1, const Scalar& s2, int )
{
    Real out = eval_scalar_impl<test_function_fma>(s1, s2);
    return out;
};

Real test_function_fms::eval_mat(const Matrix&,const Matrix&, int )
{
    return 0.0;
};
Real test_function_fms::eval_scalar(const Scalar& s1, const Scalar& s2, int )
{
    Real out = eval_scalar_impl<test_function_fms>(s1, s2);
    return out;
};

Real test_function_dot2_ac::eval_mat(const Matrix&,const Matrix&, int )
{
    return 0.0;
};
Real test_function_dot2_ac::eval_scalar(const Scalar& s1, const Scalar& s2, int )
{
    Real out = eval_scalar_impl<test_function_dot2_ac>(s1, s2);
    return out;
};

Real test_function_nextafter::eval_mat(const Matrix&,const Matrix&, int )
{
    return 0.0;
};
Real test_function_nextafter::eval_scalar(const Scalar& s1, const Scalar& s2, int )
{
    Real out = eval_scalar_impl<test_function_nextafter>(s1, s2);
    return out;
};

Real test_function_powm1::eval_mat(const Matrix&,const Matrix&, int )
{
    return 0.0;
};
Real test_function_powm1::eval_scalar(const Scalar& s1, const Scalar& s2, int )
{
    Real out = eval_scalar_impl<test_function_powm1>(s1, s2);
    return out;
};

};};
