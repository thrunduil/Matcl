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
#include "test_functions_vecfunc.h"

#include "test/test_matcl/framework/matrix_set/matrix_set_1.h"
#include "matcl-core/IO/logger.h"
#include "test/test_matcl/framework/data_struct/test_selector.h"

#include <boost/thread.hpp>


namespace gd = matcl::details;

namespace matcl { namespace test
{

class test_vecfunc
{
    vecfunc_functions_list&		tf;
    const test::options&		opts;
    Integer                     thread_id;

    public:
        test_vecfunc(vecfunc_functions_list& tf, const test::options& opts, Integer id)
            :tf(tf),opts(opts), thread_id(id)
        {};

        test_vecfunc(const test_vecfunc& tu)
            :tf(tu.tf),opts(tu.opts), thread_id(tu.thread_id)
        {};

        void make()
        {          
            /*
            Integer code = 499;
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
        test_vecfunc& operator=(const test_vecfunc&) = delete;
};

void test_vecfunc_st(const rand_matrix_ptr& rand)
{
    test::options opts;

    try
    {
        opts.show_partial_res   = 0;
        opts.first_matrix_code  = 0;
        opts.show_memleaks      = true;

        test::mat_set_1 ms1(rand);
        vecfunc_functions_list tf(ms1,rand);
        
        test_vecfunc tu(tf,opts,0);
        tu.make();
    }
    catch(const std::exception& ex)
    {
        matcl::out_stream <<ex.what();
    };
};

void test_vecfunc_mt(const rand_matrix_ptr& rand)
{
    test::options opts;
    try
    {
        opts.show_partial_res   = 0;
        opts.first_matrix_code  = 0;
        opts.show_memleaks      = false;

        test::mat_set_2 ms1(rand);
        vecfunc_functions_list tf(ms1,rand);

        boost::thread_group tg;

        for (int i = 0; i < 10; i++)
            tg.create_thread(test_vecfunc(tf,opts,i));

        tg.join_all();
    }
    catch(const std::exception& ex)
    {
        matcl::out_stream <<ex.what();
    };
};

vecfunc_functions_list::vecfunc_functions_list(const matrix_set& ms,rand_matrix_ptr rand)
:m_tests(ms), m_rand(rand)
{};

void vecfunc_functions_list::make(options opts)
{
    m_options = opts;

    SELECT_TEST (3, test_nnz());
    SELECT_TEST (3, test_sum());
    SELECT_TEST (3, test_cumsum());
    SELECT_TEST (3, test_prod());
    SELECT_TEST (3, test_cumprod());
    //SELECT_TEST (3, test_sumsq());    
    SELECT_TEST (3, test_min());
    SELECT_TEST (3, test_max());
    SELECT_TEST (3, test_min_abs());
    SELECT_TEST (3, test_max_abs());
    SELECT_TEST (3, test_mean());
    SELECT_TEST (3, test_std());

    SELECT_TEST (3, test_min2());
    SELECT_TEST (3, test_max2());
    SELECT_TEST (3, test_min_abs2());
    SELECT_TEST (3, test_max_abs2());
    SELECT_TEST (3, test_any());
    SELECT_TEST (3, test_all());

    SELECT_TEST (3, test_nnz_vec());
    SELECT_TEST (3, test_all_vec());
    SELECT_TEST (3, test_any_vec());
    SELECT_TEST (3, test_sum_vec());
    SELECT_TEST (3, test_prod_vec());
    SELECT_TEST (3, test_mean_vec());
    SELECT_TEST (3, test_std_vec());
    SELECT_TEST (3, test_min_vec());
    SELECT_TEST (3, test_max_vec());
    SELECT_TEST (3, test_min_abs_vec());
    SELECT_TEST (3, test_max_abs_vec());
    SELECT_TEST (3, test_min2_vec());
    SELECT_TEST (3, test_max2_vec());
    SELECT_TEST (3, test_min_abs2_vec());
    SELECT_TEST (3, test_max_abs2_vec());
};

Matrix vecfunc_functions_list::get_matrix(int code) const
{
    return m_tests.get_matrix(code);
};

//
void vecfunc_functions_list::test_nnz()
{
    Real out = 0.;
    for (int i = 1; i <= 2; i++)
    {
        test_function_nnz tf(i);
        out += m_tests.make(&tf,m_options);
    };

    if (out == 0.)
        matcl::out_stream << std::string() + "nnz: OK" + "\n";
    else
        matcl::out_stream << std::string() + "nnz: FAILED"  + "\n";   
};
void vecfunc_functions_list::test_nnz_vec()
{
    Real out = 0.;
    test_function_nnz_vec tf;
    out += m_tests.make(&tf,m_options);

    if (out == 0.)
        matcl::out_stream << std::string() + "nnz_vec: OK" + "\n";
    else
        matcl::out_stream << std::string() + "nnz_vec: FAILED"  + "\n";   
};

void vecfunc_functions_list::test_sum()
{
    Real out = 0.;
    for (int i = 1; i <= 2; i++)
    {
        test_function_sum tf(i);
        out += m_tests.make(&tf,m_options);
    };

    if (out == 0.)
        matcl::out_stream << std::string() + "sum: OK" + "\n";
    else
        matcl::out_stream << std::string() + "sum: FAILED"  + "\n";
};
void vecfunc_functions_list::test_sum_vec()
{
    Real out = 0.;

    test_function_sum_vec tf;
    out += m_tests.make(&tf,m_options);

    if (out == 0.)
        matcl::out_stream << std::string() + "sum_vec: OK" + "\n";
    else
        matcl::out_stream << std::string() + "sum_vec: FAILED"  + "\n";
};

void vecfunc_functions_list::test_cumsum()
{
    Real out = 0.;
    for (int i = 1; i <= 2; i++)
    {
        test_function_cumsum tf(i);
        out += m_tests.make(&tf,m_options);
    };

    if (out == 0.)
        matcl::out_stream << std::string() + "cumsum: OK" + "\n";
    else
        matcl::out_stream << std::string() + "cumsum: FAILED"  + "\n";
};

void vecfunc_functions_list::test_prod()
{
    Real out = 0.;
    for (int i = 1; i <= 2; i++)
    {
        test_function_prod tf(i);
        out += m_tests.make(&tf,m_options);
    };

    if (out == 0.)
        matcl::out_stream << std::string() + "prod: OK" + "\n";
    else
        matcl::out_stream << std::string() + "prod: FAILED"  + "\n";
};
void vecfunc_functions_list::test_prod_vec()
{
    Real out = 0.;

    test_function_prod_vec tf;
    out += m_tests.make(&tf,m_options);

    if (out == 0.)
        matcl::out_stream << std::string() + "prod_vec: OK" + "\n";
    else
        matcl::out_stream << std::string() + "prod_vec: FAILED"  + "\n";
};

void vecfunc_functions_list::test_cumprod()
{
    Real out = 0.;
    for (int i = 1; i <= 2; i++)
    {
        test_function_cumprod tf(i);
        out += m_tests.make(&tf,m_options);
    };

    if (out == 0.)
        matcl::out_stream << std::string() + "cumprod: OK" + "\n";
    else
        matcl::out_stream << std::string() + "cumprod: FAILED"  + "\n";
};

void vecfunc_functions_list::test_min()
{
    Real out = 0.;
    for (int i = 1; i <= 2; i++)
    {
        test_function_min tf(i);
        out += m_tests.make(&tf,m_options);
    };

    if (out == 0.)
        matcl::out_stream << std::string() + "min: OK" + "\n";
    else
        matcl::out_stream << std::string() + "min: FAILED"  + "\n";
};
void vecfunc_functions_list::test_min_vec()
{
    Real out = 0.;
    for (int i = 1; i <= 2; i++)
    {
        test_function_min_vec tf(i);
        out += m_tests.make(&tf,m_options);
    };

    if (out == 0.)
        matcl::out_stream << std::string() + "min_vec: OK" + "\n";
    else
        matcl::out_stream << std::string() + "min_vec: FAILED"  + "\n";
};

void vecfunc_functions_list::test_min_abs()
{
    Real out = 0.;
    for (int i = 1; i <= 2; i++)
    {
        test_function_min_abs tf(i);
        out += m_tests.make(&tf,m_options);
    };

    if (out == 0.)
        matcl::out_stream << std::string() + "min_abs: OK" + "\n";
    else
        matcl::out_stream << std::string() + "min_abs: FAILED"  + "\n";
};
void vecfunc_functions_list::test_min_abs_vec()
{
    Real out = 0.;

    test_function_min_abs_vec tf;
    out += m_tests.make(&tf,m_options);

    if (out == 0.)
        matcl::out_stream << std::string() + "min_abs_vec: OK" + "\n";
    else
        matcl::out_stream << std::string() + "min_abs_vec: FAILED"  + "\n";
};

void vecfunc_functions_list::test_min2()
{
    Real out = 0.;
    for (int i = 1; i <= 2; i++)
    {
        test_function_min2 tf(i);
        out += m_tests.make(&tf,m_options);
    };

    if (out == 0.)
        matcl::out_stream << std::string() + "min2: OK" + "\n";
    else
        matcl::out_stream << std::string() + "min2: FAILED"  + "\n";
};
void vecfunc_functions_list::test_min2_vec()
{
    Real out = 0.;

    test_function_min2_vec tf;
    out += m_tests.make(&tf,m_options);

    if (out == 0.)
        matcl::out_stream << std::string() + "min2_vec: OK" + "\n";
    else
        matcl::out_stream << std::string() + "min2_vec: FAILED"  + "\n";
};

void vecfunc_functions_list::test_min_abs2()
{
    Real out = 0.;
    for (int i = 1; i <= 2; i++)
    {
        test_function_min_abs2 tf(i);
        out += m_tests.make(&tf,m_options);
    };

    if (out == 0.)
        matcl::out_stream << std::string() + "min_abs2: OK" + "\n";
    else
        matcl::out_stream << std::string() + "min_abs2: FAILED"  + "\n";
};
void vecfunc_functions_list::test_min_abs2_vec()
{
    Real out = 0.;

    test_function_min_abs2_vec tf;
    out += m_tests.make(&tf,m_options);

    if (out == 0.)
        matcl::out_stream << std::string() + "min_abs2_vec: OK" + "\n";
    else
        matcl::out_stream << std::string() + "min_abs2_vec: FAILED"  + "\n";
};

void vecfunc_functions_list::test_max()
{
    Real out = 0.;
    for (int i = 1; i <= 2; i++)
    {
        test_function_max tf(i);
        out += m_tests.make(&tf,m_options);
    };

    if (out == 0.)
        matcl::out_stream << std::string() + "max: OK" + "\n";
    else
        matcl::out_stream << std::string() + "max: FAILED"  + "\n";
};
void vecfunc_functions_list::test_max_vec()
{
    Real out = 0.;

    test_function_max_vec tf;
    out += m_tests.make(&tf,m_options);

    if (out == 0.)
        matcl::out_stream << std::string() + "max_vec: OK" + "\n";
    else
        matcl::out_stream << std::string() + "max_vec: FAILED"  + "\n";
};

void vecfunc_functions_list::test_max_abs()
{
    Real out = 0.;
    for (int i = 1; i <= 2; i++)
    {
        test_function_max_abs tf(i);
        out += m_tests.make(&tf,m_options);
    };

    if (out == 0.)
        matcl::out_stream << std::string() + "max_abs: OK" + "\n";
    else
        matcl::out_stream << std::string() + "max_abs: FAILED"  + "\n";
};
void vecfunc_functions_list::test_max_abs_vec()
{
    Real out = 0.;

    test_function_max_abs_vec tf;
    out += m_tests.make(&tf,m_options);

    if (out == 0.)
        matcl::out_stream << std::string() + "max_abs_vec: OK" + "\n";
    else
        matcl::out_stream << std::string() + "max_abs_vec: FAILED"  + "\n";
};

void vecfunc_functions_list::test_max2()
{
    Real out = 0.;
    for (int i = 1; i <= 2; i++)
    {
        test_function_max2 tf(i);
        out += m_tests.make(&tf,m_options);
    };

    if (out == 0.)
        matcl::out_stream << std::string() + "max2: OK" + "\n";
    else
        matcl::out_stream << std::string() + "max2: FAILED"  + "\n";
};	
void vecfunc_functions_list::test_max2_vec()
{
    Real out = 0.;
    test_function_max2_vec tf;
    out += m_tests.make(&tf,m_options);

    if (out == 0.)
        matcl::out_stream << std::string() + "max2_vec: OK" + "\n";
    else
        matcl::out_stream << std::string() + "max2_vec: FAILED"  + "\n";
};	

void vecfunc_functions_list::test_max_abs2()
{
    Real out = 0.;
    for (int i = 1; i <= 2; i++)
    {
        test_function_max_abs2 tf(i);
        out += m_tests.make(&tf,m_options);
    };

    if (out == 0.)
        matcl::out_stream << std::string() + "max_abs2: OK" + "\n";
    else
        matcl::out_stream << std::string() + "max_abs2: FAILED"  + "\n";
};	
void vecfunc_functions_list::test_max_abs2_vec()
{
    Real out = 0.;

    test_function_max_abs2_vec tf;
    out += m_tests.make(&tf,m_options);

    if (out == 0.)
        matcl::out_stream << std::string() + "max_abs2_vec: OK" + "\n";
    else
        matcl::out_stream << std::string() + "max_abs2_vec: FAILED"  + "\n";
};	

void vecfunc_functions_list::test_mean()
{
    Real out = 0.;
    for (int i = 1; i <= 2; i++)
    {
        test_function_mean tf(i);
        out += m_tests.make(&tf,m_options);
    };

    if (out == 0.)
        matcl::out_stream << std::string() +   "mean: OK" + "\n";
    else
        matcl::out_stream << std::string() +   "mean: FAILED"  + "\n";
};
void vecfunc_functions_list::test_mean_vec()
{
    Real out = 0.;

    test_function_mean_vec tf;
    out += m_tests.make(&tf,m_options);

    if (out == 0.)
        matcl::out_stream << std::string() +   "mean_vec: OK" + "\n";
    else
        matcl::out_stream << std::string() +   "mean_vec: FAILED"  + "\n";
};

void vecfunc_functions_list::test_std()
{
    Real out = 0.;
    for (int i = 1; i <= 2; i++)
    {
        test_function_std tf(i);
        out += m_tests.make(&tf,m_options);
    };

    if (out == 0.)
        matcl::out_stream << std::string() +   "std: OK" + "\n";
    else
        matcl::out_stream << std::string() +   "std: FAILED"  + "\n";
};
void vecfunc_functions_list::test_std_vec()
{
    Real out = 0.;
    for (int i = 1; i <= 2; i++)
    {
        test_function_std_vec tf(i);
        out += m_tests.make(&tf,m_options);
    };

    if (out == 0.)
        matcl::out_stream << std::string() +   "std_vec: OK" + "\n";
    else
        matcl::out_stream << std::string() +   "std_vec: FAILED"  + "\n";
};

void vecfunc_functions_list::test_all()
{
    Real out = 0.;
    for (int i = 1; i <= 2; i++)
    {
        test_function_all tf(i);
        out += m_tests.make(&tf,m_options);
    };

    if (out == 0.)
        matcl::out_stream << std::string() +   "all: OK" + "\n";
    else
        matcl::out_stream << std::string() +   "all: FAILED"  + "\n";
};
void vecfunc_functions_list::test_all_vec()
{
    Real out = 0.;

    test_function_all_vec tf;
    out += m_tests.make(&tf,m_options);

    if (out == 0.)
        matcl::out_stream << std::string() +   "all_vec: OK" + "\n";
    else
        matcl::out_stream << std::string() +   "all_vec: FAILED"  + "\n";
};

void vecfunc_functions_list::test_any()
{
    Real out = 0.;
    for (int i = 1; i <= 2; i++)
    {
        test_function_any tf(i);
        out += m_tests.make(&tf,m_options);
    };

    if (out == 0.)
        matcl::out_stream << std::string() +   "any: OK" + "\n";
    else
        matcl::out_stream << std::string() +   "any: FAILED"  + "\n";
};
void vecfunc_functions_list::test_any_vec()
{
    Real out = 0.;
    for (int i = 1; i <= 2; i++)
    {
        test_function_any_vec tf(i);
        out += m_tests.make(&tf,m_options);
    };

    if (out == 0.)
        matcl::out_stream << std::string() +   "any_vec: OK" + "\n";
    else
        matcl::out_stream << std::string() +   "any_vec: FAILED"  + "\n";
};

//-----------------------------------------------------------------
Real test_function_sum::eval_mat(const Matrix& mat,bool,int code )
{
    (void)code;
    Matrix out_full = sum(full(mat),d);	
    try
    {		
        Matrix out		= sum(mat,d);
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
Real test_function_sum::eval_scalar(const Scalar& s,bool,int code  )
{
    (void)code;
    return eval_scalar_impl<test_function_sum>(s);
};

Real test_function_sum_vec::eval_mat(const Matrix& mat,bool,int code )
{
    (void)code;
    Matrix out_full = sum_vec(full(mat));	
    try
    {		
        Matrix out		= sum_vec(mat);
        Matrix out2     = sum(sum(mat,1),2);

        check_struct(out);
        Real dif        = norm_1(out - out_full);
        dif             += norm_1(out - out2);

        if (dif < sum_vec(abs(mat)).get_scalar<Real>() * constants::eps(mat.get_value_code()) * 10.0)
            dif         = 0;

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
Real test_function_sum_vec::eval_scalar(const Scalar& s,bool,int code  )
{
    (void)code;
    return eval_scalar_impl<test_function_sum_vec>(s);
};

Real test_function_cumsum::eval_mat(const Matrix& mat,bool,int code )
{
    (void)code;
    Matrix out_full = cumsum(full(mat),d);

    // test by matrix multiplication
    // for d==1 (by cols) cumsum=cumsum_factor*mat, 
    // where cumsum_factor has 0 above diag, 1 elsewhere
    Integer csf_dim         = d==1 ? mat.rows() : mat.cols();
    Matrix cumsum_factor    = tril(iones(csf_dim,csf_dim));    
    Matrix cumsum_by_mult   = d==1 ? cumsum_factor*mat : mat*trans(cumsum_factor);

    try
    {
        Matrix out		= cumsum(mat,d);

        Real dif        =  0;
        dif             += norm_1(out - cumsum_by_mult);
        //take roundoff error of testing by matrix multiplication into account
        if (abs(dif) < ::powl(::sqrtl(csf_dim),3) * norm_1(mat) * constants::eps(mat.get_value_code()))
            dif = 0;

        dif             += norm_1(out - out_full);
        check_struct(out);
        return dif ;
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
Real test_function_cumsum::eval_scalar(const Scalar& s,bool,int code  )
{
    (void)code;
    return eval_scalar_impl<test_function_cumsum>(s);
};

Real test_function_prod::eval_mat(const Matrix& mat,bool,int code )
{
    (void)code;
    Matrix out_full = prod(full(mat),d);	
    try
    {		
        Matrix out		= prod(mat,d);
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
Real test_function_prod::eval_scalar(const Scalar& s,bool,int code  )
{
    (void)code;
    return eval_scalar_impl<test_function_prod>(s);
};

Real test_function_prod_vec::eval_mat(const Matrix& mat,bool,int code )
{
    (void)code;
    Matrix out_full = prod_vec(full(mat));	
    try
    {		
        Matrix out		= prod_vec(mat);
        Matrix out2     = prod(prod(mat,1),2);
        check_struct(out);

        Real dif        = norm_1(out - out_full);
        dif             += norm_1(out - out2);

        Real scal       = prod_vec(abs(mat)).get_scalar<Real>();
        if (dif < scal * constants::eps(mat.get_value_code()) * mat.numel() * 10.0)
            dif         = 0;

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
Real test_function_prod_vec::eval_scalar(const Scalar& s,bool,int code  )
{
    (void)code;
    return eval_scalar_impl<test_function_prod_vec>(s);
};

Real test_function_cumprod::eval_mat(const Matrix& mat,bool,int code )
{
    (void)code;
    Matrix out_full = cumprod(full(mat),d);	
    try
    {		
        Matrix out		= cumprod(mat,d);
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

Real test_function_cumprod::eval_scalar(const Scalar& s,bool,int code  )
{
    (void)code;
    return eval_scalar_impl<test_function_cumprod>(s);
};

Real test_function_nnz::eval_mat(const Matrix& mat,bool,int code )
{
    (void)code;
    Matrix out_full = nnz(full(mat),d);	
    try
    {		
        Matrix out		= nnz(mat,d);
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

Real test_function_nnz::eval_scalar(const Scalar& s,bool,int code  )
{
    (void)code;
    return eval_scalar_impl<test_function_nnz>(s);
};

Real test_function_nnz_vec::eval_mat(const Matrix& mat,bool,int code )
{
    (void)code;
    Matrix out_full = nnz_vec(full(mat));	
    try
    {		
        Matrix out		= nnz_vec(mat);
        check_struct(out);
        Matrix out2     = sum_vec(nnz(mat,1));

        Real dif        = norm_1(out - out_full);
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

Real test_function_nnz_vec::eval_scalar(const Scalar& s,bool,int code  )
{
    (void)code;
    return eval_scalar_impl<test_function_nnz_vec>(s);
};

Real test_function_min::eval_mat(const Matrix& mat,bool,int code )
{
    (void)code;
    Matrix out_full = min_d(full(mat),d);	
    try
    {		
        Matrix out		= min_d(mat,d);
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

Real test_function_min::eval_scalar(const Scalar& s,bool,int code  )
{
    (void)code;
    return eval_scalar_impl<test_function_min>(s);
};

Real test_function_min_vec::eval_mat(const Matrix& mat,bool,int code )
{
    (void)code;
    Matrix out_full = min_vec(full(mat));	
    try
    {		
        Matrix out		= min_vec(mat);
        Matrix out2     = min_d(min_d(mat,1),2);

        check_struct(out);
        Real dif        = norm_1(out - out_full);
        dif             +=norm_1(out - out2);
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

Real test_function_min_vec::eval_scalar(const Scalar& s,bool,int code  )
{
    (void)code;
    return eval_scalar_impl<test_function_min_vec>(s);
};

Real test_function_min_abs::eval_mat(const Matrix& mat,bool,int code )
{
    (void)code;
    Matrix out_full = min_abs_d(full(mat),d);	
    try
    {		
        Matrix out1	= min_abs_d(mat,d);
        check_struct(out1);
        Matrix out2	= min_d(abs(mat),d);
        check_struct(out2);
        Real nrm    = norm_1(out1 - out_full);
        nrm         = nrm + norm_1(out2 - out_full);
        return nrm;
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

Real test_function_min_abs::eval_scalar(const Scalar& s,bool,int code  )
{
    (void)code;
    return eval_scalar_impl<test_function_min_abs>(s);
};

Real test_function_min_abs_vec::eval_mat(const Matrix& mat,bool,int code )
{
    (void)code;
    Matrix out_full = min_abs_vec(full(mat));	
    try
    {		
        Matrix out1	= min_abs_vec(mat);
        check_struct(out1);
        Matrix out2	= min_d(vec(abs(mat)), 1);
        check_struct(out2);

        Real nrm    = norm_1(out1 - out_full);
        nrm         = nrm + norm_1(out2 - out_full);
        return nrm;
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

Real test_function_min_abs_vec::eval_scalar(const Scalar& s,bool,int code  )
{
    (void)code;
    return eval_scalar_impl<test_function_min_abs_vec>(s);
};

Real test_function_max::eval_mat(const Matrix& mat,bool,int code )
{
    (void)code;
    Matrix out_full = max_d(full(mat),d);	
    try
    {		
        Matrix out		= max_d(mat,d);
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

Real test_function_max::eval_scalar(const Scalar& s,bool,int code  )
{
    (void)code;
    return eval_scalar_impl<test_function_max>(s);
};

Real test_function_max_vec::eval_mat(const Matrix& mat,bool,int code )
{
    (void)code;
    Matrix out_full = max_vec(full(mat));	
    try
    {		
        Matrix out		= max_vec(mat);
        Matrix out2     = max_d(max_d(mat,1),2);

        check_struct(out);
        Real dif        = norm_1(out - out_full);
        dif             +=norm_1(out - out2);
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

Real test_function_max_vec::eval_scalar(const Scalar& s,bool,int code  )
{
    (void)code;
    return eval_scalar_impl<test_function_max_vec>(s);
};

Real test_function_max_abs::eval_mat(const Matrix& mat,bool,int code )
{
    (void)code;
    Matrix out_full = max_abs_d(full(mat),d);	
    try
    {		
        Matrix out1	= max_abs_d(mat,d);
        Matrix out2	= max_d(abs(mat),d);
        check_struct(out1);
        check_struct(out2);
        Real nrm    = norm_1(out1 - out_full);
        nrm         = nrm + norm_1(out2 - out_full);
        return nrm;
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

Real test_function_max_abs::eval_scalar(const Scalar& s,bool,int code  )
{
    (void)code;
    return eval_scalar_impl<test_function_max_abs>(s);
};

Real test_function_max_abs_vec::eval_mat(const Matrix& mat,bool,int code )
{
    (void)code;
    Matrix out_full = max_abs_vec(full(mat));	
    try
    {		
        Matrix out1	= max_abs_vec(mat);
        Matrix out2	= max_d(max_d(abs(mat),1),2);
        check_struct(out1);

        Real nrm    = norm_1(out1 - out_full);
        nrm         = nrm + norm_1(out2 - out_full);
        return nrm;
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

Real test_function_max_abs_vec::eval_scalar(const Scalar& s,bool,int code  )
{
    (void)code;
    return eval_scalar_impl<test_function_max_abs_vec>(s);
};

Real test_function_min2::eval_mat(const Matrix& mat,bool,int code )
{
    (void)code;
    Real dif = 0;

    Matrix out_full, iout_full;
    tie(out_full,iout_full) = min2(full(mat),d);	
    
    try
    {		
        Matrix out, iout;
        tie(out,iout) = min2(mat,d);
        check_struct(out);
        check_struct(iout);

        dif = norm_1(out - out_full);
        dif = dif + norm_1(iout - iout_full);

        Matrix pos;
        if (d == 1)
        {
            pos = iout;
            if (mat.rows()>0)
                pos = pos + irange(0,mat.cols()-1)*mat.rows();
        }
        else
        {
            pos = (iout-1)*mat.rows();
            if (mat.cols()>0)
                pos = pos + trans(irange(1,mat.rows()));
        };

        dif = dif + norm_1(out - mat(pos));
        dif = dif + norm_1(out - min_d(mat,d));
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

    return dif;
};

Real test_function_min2::eval_scalar(const Scalar& s,bool,int code  )
{
    (void)code;
    return eval_scalar_impl<test_function_min2>(s);
};

Real test_function_min2_vec::eval_mat(const Matrix& mat,bool,int code )
{
    (void)code;
    Real dif = 0;

    Matrix out_full;
    Integer rout_full, cout_full;
    tie(out_full,rout_full, cout_full) = min2_vec(full(mat));	
    
    try
    {		
        Matrix out;
        Integer rout, cout;
        tie(out,rout,cout) = min2_vec(mat);
        check_struct(out);

        dif += norm_1(out - out_full);

        if (mat.rows() > 0 && mat.cols() > 0)
            dif += norm_1(out - mat(rout,cout));

        Matrix out2, iout;
        tie(out2,iout) = min2(min2(mat,1).get<1>(),2);
        
        dif = dif + norm_1(out - out2);
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

    return dif;
};

Real test_function_min2_vec::eval_scalar(const Scalar& s,bool,int code  )
{
    (void)code;
    return eval_scalar_impl<test_function_min2_vec>(s);
};

Real test_function_min_abs2_vec::eval_mat(const Matrix& mat,bool,int code )
{
    (void)code;
    Real dif = 0;

    Matrix out_full;
    Integer rout_full, cout_full;
    tie(out_full,rout_full, cout_full) = min_abs2_vec(full(mat));	
    
    try
    {		
        Matrix out;
        Integer rout, cout;
        tie(out,rout,cout) = min_abs2_vec(mat);
        check_struct(out);

        dif += norm_1(out - out_full);

        if (mat.rows() > 0 && mat.cols() > 0)
            dif += norm_1(out - abs(mat(rout,cout)));

        Matrix out2, iout;
        tie(out2,iout) = min_abs2(min_abs2(mat,1).get<1>(),2);
        
        dif = dif + norm_1(out - out2);
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

    return dif;
};

Real test_function_min_abs2_vec::eval_scalar(const Scalar& s,bool,int code  )
{
    (void)code;
    return eval_scalar_impl<test_function_min_abs2_vec>(s);
};

Real test_function_max2_vec::eval_mat(const Matrix& mat,bool,int code )
{
    (void)code;
    Real dif = 0;

    Matrix out_full;
    Integer rout_full, cout_full;
    tie(out_full,rout_full, cout_full) = max2_vec(full(mat));	
    
    try
    {		
        Matrix out;
        Integer rout, cout;
        tie(out,rout,cout) = max2_vec(mat);
        check_struct(out);

        dif += norm_1(out - out_full);

        if (mat.rows() > 0 && mat.cols() > 0)
            dif += norm_1(out - mat(rout,cout));

        Matrix out2, iout;
        tie(out2,iout) = max2(max2(mat,1).get<1>(),2);
        
        dif = dif + norm_1(out - out2);
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

    return dif;
};

Real test_function_max2_vec::eval_scalar(const Scalar& s,bool,int code  )
{
    (void)code;
    return eval_scalar_impl<test_function_max2_vec>(s);
};

Real test_function_max_abs2_vec::eval_mat(const Matrix& mat,bool,int code )
{
    (void)code;
    Real dif = 0;

    Matrix out_full;
    Integer rout_full, cout_full;
    tie(out_full,rout_full, cout_full) = max_abs2_vec(full(mat));	
    
    try
    {		
        Matrix out;
        Integer rout, cout;
        tie(out,rout,cout) = max_abs2_vec(mat);
        check_struct(out);

        dif += norm_1(out - out_full);

        if (mat.rows() > 0 && mat.cols() > 0)
            dif += norm_1(out - abs(mat(rout,cout)));

        Matrix out2, iout;
        tie(out2,iout) = max_abs2(max_abs2(mat,1).get<1>(),2);
        
        dif = dif + norm_1(out - out2);
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

    return dif;
};

Real test_function_max_abs2_vec::eval_scalar(const Scalar& s,bool,int code  )
{
    (void)code;
    return eval_scalar_impl<test_function_max_abs2_vec>(s);
};

Real test_function_min_abs2::eval_mat(const Matrix& mat,bool,int code )
{
    (void)code;
    Real dif = 0;

    Matrix out_full, iout_full;
    tie(out_full,iout_full) = min_abs2(full(mat),d);	
    
    try
    {		
        Matrix out1, iout1;
        tie(out1,iout1) = min_abs2(mat,d);
        check_struct(out1);
        check_struct(iout1);

        Matrix out2, iout2;
        tie(out2,iout2) = min2(abs(mat),d);
        check_struct(out2);
        check_struct(iout2);

        dif += norm_1(out1 - out_full);
        dif += norm_1(iout1 - iout_full);
        dif += norm_1(out2 - out_full);
        dif += norm_1(iout2 - iout_full);

        Matrix pos;
        if (d == 1)
        {
            pos = iout1;
            if (mat.rows()>0)
                pos = pos + irange(0,mat.cols()-1)*mat.rows();
        }
        else
        {
            pos = (iout1-1)*mat.rows();
            if (mat.cols()>0)
                pos = pos + trans(irange(1,mat.rows()));
        };

        dif = dif + norm_1(out1 - abs(mat(pos)));
        dif = dif + norm_1(out1 - min_abs_d(mat,d));
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

    return dif;
};

Real test_function_min_abs2::eval_scalar(const Scalar& s,bool,int code  )
{
    (void)code;
    return eval_scalar_impl<test_function_min_abs2>(s);
};

Real test_function_max2::eval_mat(const Matrix& mat,bool,int code )
{
    (void)code;
    Matrix out_full, iout_full;
    tie(out_full,iout_full) = max2(full(mat),d);	
    try
    {		
        Matrix out, iout;
        tie(out,iout) = max2(mat,d);
        check_struct(out);
        check_struct(iout);

        Real dif = norm_1(out - out_full);
        dif = dif + norm_1(iout - iout_full);

        Matrix pos;
        if (d == 1)
        {
            pos = iout;
            if (mat.rows()>0)
                pos = pos + irange(0,mat.cols()-1)*mat.rows();
        }
        else
        {
            pos = (iout-1)*mat.rows();
            if (mat.cols()>0)
                pos = pos + trans(irange(1,mat.rows()));
        };

        dif = dif + norm_1(out - mat(pos));
        dif = dif + norm_1(out - max_d(mat,d));

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

Real test_function_max2::eval_scalar(const Scalar& s,bool,int code  )
{
    (void)code;
    return eval_scalar_impl<test_function_max2>(s);
};

Real test_function_max_abs2::eval_mat(const Matrix& mat,bool,int code )
{
    (void)code;
    Real dif = 0;

    Matrix out_full, iout_full;
    tie(out_full,iout_full) = max_abs2(full(mat),d);	
    
    try
    {		
        Matrix out1, iout1;
        tie(out1,iout1) = max_abs2(mat,d);
        check_struct(out1);
        check_struct(iout1);

        Matrix out2, iout2;
        tie(out2,iout2) = max2(abs(mat),d);
        check_struct(out2);
        check_struct(iout2);

        dif += norm_1(out1 - out_full);
        dif += norm_1(iout1 - iout_full);
        dif += norm_1(out2 - out_full);
        dif += norm_1(iout2 - iout_full);

        Matrix pos;
        if (d == 1)
        {
            pos = iout1;
            if (mat.rows()>0)
                pos = pos + irange(0,mat.cols()-1)*mat.rows();
        }
        else
        {
            pos = (iout1-1)*mat.rows();
            if (mat.cols()>0)
                pos = pos + trans(irange(1,mat.rows()));
        };

        dif = dif + norm_1(out1 - abs(mat(pos)));
        dif = dif + norm_1(out1 - max_abs_d(mat,d));
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

    return dif;
};

Real test_function_max_abs2::eval_scalar(const Scalar& s,bool,int code  )
{
    (void)code;
    return eval_scalar_impl<test_function_max_abs2>(s);
};

Real test_function_mean::eval_mat(const Matrix& mat,bool,int code )
{
    (void)code;

    Matrix out_full = mean(full(mat),d);	
    try
    {		
        Matrix out		= mean(mat,d);
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

Real test_function_mean::eval_scalar(const Scalar& s,bool,int code  )
{
    (void)code;
    return eval_scalar_impl<test_function_mean>(s);
};

Real test_function_mean_vec::eval_mat(const Matrix& mat,bool,int code )
{
    (void)code;

    Matrix out_full = mean_vec(full(mat));	
    try
    {		
        Matrix out		= mean_vec(mat);
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

Real test_function_mean_vec::eval_scalar(const Scalar& s,bool,int code  )
{
    (void)code;
    return eval_scalar_impl<test_function_mean_vec>(s);
};

Real test_function_std::eval_mat(const Matrix& mat,bool,int code )
{
    (void)code;

    Matrix out_full_1   = std(full(mat),d, true);	
    Matrix out_full_2   = std(full(mat),d, false);	
    
    try
    {		
        Matrix out_1	= std(mat,d, true);
        Matrix out_2    = std(mat,d, false);

        check_struct(out_1);
        check_struct(out_2);

        Real dif        = norm_1(out_1 - out_full_1);
        dif             += norm_1(out_2 - out_full_2);
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

Real test_function_std::eval_scalar(const Scalar& s,bool,int code  )
{
    (void)code;
    return eval_scalar_impl<test_function_std>(s);
};

Real test_function_std_vec::eval_mat(const Matrix& mat,bool,int code )
{
    (void)code;

    Matrix out_full_1   = std_vec(full(mat), true);	
    Matrix out_full_2   = std_vec(full(mat), false);	
    
    try
    {		
        Matrix out_1	= std_vec(mat, true);
        Matrix out_2    = std_vec(mat, false);

        check_struct(out_1);
        check_struct(out_2);

        Real dif        = norm_1(out_1 - out_full_1);
        dif             += norm_1(out_2 - out_full_2);
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

Real test_function_std_vec::eval_scalar(const Scalar& s,bool,int code  )
{
    (void)code;
    return eval_scalar_impl<test_function_std_vec>(s);
};

class Test_all : public matcl::test_function
{
    public:
        virtual bool eval(Integer v) const 
        {
            return v != 0;
        };
        virtual bool eval(Real v) const
        {
            return v != 0.;
        };
        virtual bool eval(Float v) const
        {
            return v != 0.;
        };
        virtual bool eval(const Complex& v) const
        {
            return v != 0.;
        };
        virtual bool eval(const Float_complex& v) const
        {
            return v != 0.;
        };
        virtual bool eval(const Object& v) const
        {
            return (bool)(v != 0.);
        };
};

Real test_function_all::eval_mat(const Matrix& mat,bool,int code )
{
    (void)code;
    Matrix out_full = all(full(mat),d);	

    try
    {		
        Matrix out		= all(mat,d);
        check_struct(out);

        Test_all ta;
        Matrix out2		= all(mat,ta,d);
        check_struct(out2);

        Real dif		= norm_1(out - out_full);
        dif				= dif + norm_1(out - out2);
        dif				= dif + norm_1(out - ~any(~mat,d));
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

Real test_function_all::eval_scalar(const Scalar& s,bool,int code  )
{
    (void)code;
    return eval_scalar_impl<test_function_all>(s);
};

Real test_function_all_vec::eval_mat(const Matrix& mat,bool,int code )
{
    (void)code;
    bool out_full = all_vec(full(mat));	

    try
    {		
        bool out		= all_vec(mat);

        Test_all ta;
        bool out2		= all_vec(mat,ta);
        Matrix out3		= all(all(mat,1),2);
        Matrix out4		= all(all(mat,ta,1),ta,2);

        Real dif		= norm_1(out - out_full);
        dif				= dif + norm_1(out - out2);
        dif				= dif + norm_1(Matrix(out) - out3);
        dif				= dif + norm_1(Matrix(out2) - out4);
        dif				= dif + norm_1(out - !any_vec(~mat));
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

Real test_function_all_vec::eval_scalar(const Scalar& s,bool,int code  )
{
    (void)code;
    return eval_scalar_impl<test_function_all_vec>(s);
};

Real test_function_any::eval_mat(const Matrix& mat,bool,int code )
{
    (void)code;
    Matrix out_full = any(full(mat),d);	

    try
    {		
        Matrix out		= any(mat,d);
        check_struct(out);

        Test_all ta;
        Matrix out2		= any(mat,ta,d);
        check_struct(out2);

        Real dif		= norm_1(out - out_full);
        dif				= dif + norm_1(out - out2);
        dif				= dif + norm_1(out - ~all(~mat,d));
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

Real test_function_any::eval_scalar(const Scalar& s,bool,int code  )
{
    (void)code;
    return eval_scalar_impl<test_function_any>(s);
};

Real test_function_any_vec::eval_mat(const Matrix& mat,bool,int code )
{
    (void)code;
    bool out_full = any_vec(full(mat));	

    try
    {		
        bool out		= any_vec(mat);

        Test_all ta;
        bool out2		= any_vec(mat,ta);

        Matrix out3     = any(any(mat,1),2);
        Matrix out4     = any(any(mat,ta,1),ta,2);

        Real dif		= norm_1(out - out_full);
        dif				= dif + norm_1(out - out2);
        dif				= dif + norm_1(Matrix(out) - out3);
        dif				= dif + norm_1(Matrix(out2) - out4);
        dif				= dif + norm_1(out - !all_vec(~mat));
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

Real test_function_any_vec::eval_scalar(const Scalar& s,bool,int code  )
{
    (void)code;
    return eval_scalar_impl<test_function_any_vec>(s);
};

};};
