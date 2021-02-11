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
#include "test_functions_matfunc.h"

#include "test/test_matcl/framework/matrix_set/matrix_set_1.h"
#include "matcl-core/IO/logger.h"
#include "test/test_matcl/framework/data_struct/test_selector.h"

#include <boost/thread.hpp>


namespace gd = matcl::details;

namespace matcl { namespace test
{

class test_matfunc
{
    matfunc_functions_list&		tf;
    const test::options&		opts;
    Integer                     thread_id;

    public:
        test_matfunc(matfunc_functions_list& tf, const test::options& opts, Integer id)
            :tf(tf),opts(opts), thread_id(id)
        {};

        test_matfunc(const test_matfunc& tu)
            :tf(tu.tf),opts(tu.opts), thread_id(tu.thread_id)
        {};

        void make()
        {   
            /*
            Integer code = 169;
            Matrix mat = tf.get_matrix(code);
            matcl::disp(mat);
	
            Matrix res1     = mat * ctrans(mat);
            Matrix res2     = herprod(mat, false);

            Matrix res3     = ctrans(mat) * mat;
            Matrix res4     = herprod(mat, true);

            Real dif1		= norm_1(res1 - res2);
            Real dif2		= norm_1(res3 - res4);

            disp(res1);
            disp(res2);
            disp(res4);

            check_struct(res2);
            check_struct(res4);

            Real dif        = dif1 + dif2;
            Real tol        = matcl::norm_1(mat);

            if (dif < 10.0 * mat.numel() * tol * tol * constants::eps(mat.get_value_code()))
                dif         = 0.;
            */

            tf.make(opts);
        };

        void operator()()
        {
            make();
        };

    private:		
        test_matfunc& operator=(const test_matfunc&) = delete;
};

void test_matfunc_st(const rand_matrix_ptr& rand)
{
    test::options opts;

    try
    {
        opts.show_partial_res   = 0;
        opts.first_matrix_code  = 0;
        opts.show_memleaks      = true;

        test::mat_set_1 ms1(rand);
        dynamic_mat_set ms(rand);
        matfunc_functions_list tf(ms1,rand,ms);
        
        test_matfunc tu(tf,opts,0);
        tu.make();
    }
    catch(const std::exception& ex)
    {
        matcl::out_stream <<ex.what();
    };
};

void test_matfunc_mt(const rand_matrix_ptr& rand)
{
    test::options opts;
    try
    {
        opts.show_partial_res   = 0;
        opts.first_matrix_code  = 0;
        opts.show_memleaks      = false;

        test::mat_set_2 ms1(rand);
        dynamic_mat_set ms(rand);
        matfunc_functions_list tf(ms1,rand,ms);

        boost::thread_group tg;

        for (int i = 0; i < 10; i++)
            tg.create_thread(test_matfunc(tf,opts,i));

        tg.join_all();
    }
    catch(const std::exception& ex)
    {
        matcl::out_stream <<ex.what();
    };
};

static const Integer max_int = 1000;

matfunc_functions_list::matfunc_functions_list(const matrix_set& ms,rand_matrix_ptr rand, dynamic_mat_set& msd)
:m_tests(ms), m_rand(rand), m_ms(msd)
{};

void matfunc_functions_list::make(options opts)
{
    m_options = opts;

    SELECT_TEST (3, test_scale_rowscols());
    SELECT_TEST (3, test_scale_cols());
    SELECT_TEST (3, test_scale_rows());        

    SELECT_TEST (3, test_symsum());
    SELECT_TEST (3, test_hersum());

    SELECT_TEST (3, test_symprod());
    SELECT_TEST (3, test_herprod());
};

Matrix matfunc_functions_list::get_matrix(int code) const
{
    return m_tests.get_matrix(code);
};

void matfunc_functions_list::test_symprod()
{
    Real out = 0.;
    test_function_symprod tf;
    out += m_tests.make(&tf,m_options);

    if (out == 0.)
        matcl::out_stream << std::string() +   "symprod: OK" + "\n";
    else
        matcl::out_stream << std::string() +   "symprod: FAILED"  + "\n";
};

void matfunc_functions_list::test_herprod()
{
    Real out = 0.;
    test_function_herprod tf;
    out += m_tests.make(&tf,m_options);

    if (out == 0.)
        matcl::out_stream << std::string() +   "herprod: OK" + "\n";
    else
        matcl::out_stream << std::string() +   "herprod: FAILED"  + "\n";
};

void matfunc_functions_list::test_symsum()
{
    Real out = 0.;
    test_function_symsum tf;
    out += m_tests.make(&tf,m_options);

    if (out == 0.)
        matcl::out_stream << std::string() +   "symsum: OK" + "\n";
    else
        matcl::out_stream << std::string() +   "symsum: FAILED"  + "\n";
};

void matfunc_functions_list::test_hersum()
{
    Real out = 0.;
    test_function_hersum tf;
    out += m_tests.make(&tf,m_options);

    if (out == 0.)
        matcl::out_stream << std::string() +   "hersum: OK" + "\n";
    else
        matcl::out_stream << std::string() +   "hersum: FAILED"  + "\n";
};

void matfunc_functions_list::test_scale_rows()
{
    Real out = 0.;
    test_function_scale_rows tf(m_ms);
    out += m_tests.make(&tf,m_options);

    if (out == 0.)
        matcl::out_stream << std::string() +   "scale_rows: OK" + "\n";
    else
        matcl::out_stream << std::string() +   "scale_rows: FAILED"  + "\n";
};

void matfunc_functions_list::test_scale_cols()
{
    Real out = 0.;
    test_function_scale_cols tf(m_ms);
    out += m_tests.make(&tf,m_options);

    if (out == 0.)
        matcl::out_stream << std::string() +   "scale_cols: OK" + "\n";
    else
        matcl::out_stream << std::string() +   "scale_cols: FAILED"  + "\n";
};

void matfunc_functions_list::test_scale_rowscols()
{
    Real out = 0.;
    test_function_scale_rowscols tf(m_ms);
    out += m_tests.make(&tf,m_options);

    if (out == 0.)
        matcl::out_stream << std::string() +   "scale_rowscols: OK" + "\n";
    else
        matcl::out_stream << std::string() +   "scale_rowscols: FAILED"  + "\n";
};

//-------------------------------------------------------------
Real test_function_symprod::eval_mat(const Matrix& mat,bool,int code )
{
    (void)code;

    try
    {		
        Matrix res1 = mat * trans(mat);
        Matrix res2 = symprod(mat, false);

        Matrix res3 = trans(mat) * mat;
        Matrix res4 = symprod(mat, true);

        Real dif1	= norm_1(res1 - res2);
        Real dif2	= norm_1(res3 - res4);

        check_struct(res2);
        check_struct(res4);

        Real tol1 = matcl::norm_1(res1) * constants::eps(mat.get_value_code());
        Real tol2 = matcl::norm_1(res3) * constants::eps(mat.get_value_code());

        if (dif1 < 10.0 * tol1 * mat.numel())
            dif1 = 0;

        if (dif2 < 10.0 * tol2 * mat.numel())
            dif2 = 0;

        Real dif = dif1 + dif2;        

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
Real test_function_symprod::eval_scalar(const Scalar& s,bool,int code  )
{
    (void)code;
    return eval_scalar_impl<test_function_symprod>(s);
};

Real test_function_herprod::eval_mat(const Matrix& mat,bool,int code )
{
    (void)code;
    try
    {		
        Matrix res1     = mat * ctrans(mat);
        Matrix res2     = herprod(mat, false);

        Matrix res3     = ctrans(mat) * mat;
        Matrix res4     = herprod(mat, true);

        Real dif1		= norm_1(res1 - res2);
        Real dif2		= norm_1(res3 - res4);

        check_struct(res2);
        check_struct(res4);

        Real dif        = dif1 + dif2;
        Real tol        = matcl::norm_1(mat);

        if (dif < 10.0 * mat.numel() * tol * tol * constants::eps(mat.get_value_code()))
            dif         = 0.;

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

Real test_function_herprod::eval_scalar(const Scalar& s,bool,int code  )
{
    (void)code;
    return eval_scalar_impl<test_function_herprod>(s);
};

Real test_function_symsum::eval_mat(const Matrix& mat,bool,int code )
{
    (void)code;

    try
    {		
        if (mat.is_square() == false)
            return 0.0;

        Matrix res1 = mat + trans(mat);
        Matrix res2 = symsum(mat);

        Real dif1	= norm_1(res1 - res2);

        check_struct(res2);

        Real dif = dif1;
        Real tol = matcl::norm_1(mat);

        if (dif < 10.0 * mat.numel() * tol * constants::eps(mat.get_value_code()))
            dif = 0;

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
Real test_function_symsum::eval_scalar(const Scalar& s,bool,int code  )
{
    (void)code;
    return eval_scalar_impl<test_function_symsum>(s);
};

Real test_function_hersum::eval_mat(const Matrix& mat,bool,int code )
{
    (void)code;

    try
    {		
        if (mat.is_square() == false)
            return 0.0;

        Matrix res1     = mat + ctrans(mat);
        Matrix res2     = hersum(mat);

        Real dif1		= norm_1(res1 - res2);

        check_struct(res2);

        Real dif        = dif1;
        Real tol        = matcl::norm_1(mat);

        if (dif < 10.0 * mat.numel() * tol * constants::eps(mat.get_value_code()))
            dif         = 0.;

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
Real test_function_hersum::eval_scalar(const Scalar& s,bool,int code  )
{
    (void)code;
    return eval_scalar_impl<test_function_hersum>(s);
};

Real test_function_scale_rows::eval_mat(const Matrix& mat,bool,int code )
{
    (void)code;

    try
    {		
        if (mat.is_square() == false)
            return 0.0;

        Integer r       = mat.rows();

        using container = dynamic_mat_set::container;
        long n1 = matcl::details::no_existing_objects();
        const container& mc = m_ms.get(r,1);
        n1 = matcl::details::no_existing_objects() - n1;
        m_new_objects = n1;

        size_t size = mc.size();

        Real dif            = 0.0;

        for (size_t i = 0; i < size ; i ++ )
        {
            Matrix D        = mc[i];
            Matrix res1     = scale_rows(mat, D);
            Matrix res2     = bdiag(D) * mat;
                
            check_struct(res1);
            
            if (check_value_code(res1.get_value_code(),res2.get_value_code()) == false)
                dif         += 1;

            Real loc_dif    = norm_1(res1 - res2);

            if (loc_dif < norm_1(res1) * constants::eps(res1.get_value_code()) * res1.numel() * 10.)
                loc_dif     = 0.0;

            dif             += loc_dif;
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
Real test_function_scale_rows::eval_scalar(const Scalar& s,bool,int code  )
{
    (void)code;
    (void)s;

    return 0.0;
};

bool test_function_scale_rows::check_value_code(value_code v1, value_code v2) const
{
    value_code v1_r = matrix_traits::real_value_type(v1);
    value_code v2_r = matrix_traits::real_value_type(v2);

    return v1_r == v2_r;
}

Real test_function_scale_cols::eval_mat(const Matrix& mat,bool,int code )
{
    (void)code;

    try
    {		
        if (mat.is_square() == false)
            return 0.0;

        Integer c       = mat.cols();

        using container = dynamic_mat_set::container;
        long n1 = matcl::details::no_existing_objects();
        const container& mc = m_ms.get(c,1);
        n1 = matcl::details::no_existing_objects() - n1;
        m_new_objects = n1;

        size_t size = mc.size();

        Real dif            = 0.0;

        for (size_t i = 0; i < size ; i ++ )
        {
            Matrix D        = mc[i];
            Matrix res1     = scale_cols(mat, D);
            Matrix res2     = mat * bdiag(D);
                
            check_struct(res1);
            
            if (check_value_code(res1.get_value_code(),res2.get_value_code()) == false)
                dif         += 1;

            Real loc_dif    = norm_1(res1 - res2);

            if (loc_dif < norm_1(res1) * constants::eps(res1.get_value_code()) * res1.numel() * 10.)
                loc_dif     = 0.0;

            dif             += loc_dif;
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
Real test_function_scale_cols::eval_scalar(const Scalar& s,bool,int code  )
{
    (void)code;
    (void)s;

    return 0.0;
};
bool test_function_scale_cols::check_value_code(value_code v1, value_code v2) const
{
    value_code v1_r = matrix_traits::real_value_type(v1);
    value_code v2_r = matrix_traits::real_value_type(v2);

    return v1_r == v2_r;
}

Real test_function_scale_rowscols::eval_mat(const Matrix& mat,bool,int code )
{
    (void)code;

    try
    {		
        if (mat.is_square() == false)
            return 0.0;

        Integer r       = mat.rows();
        Integer c       = mat.cols();

        using container = dynamic_mat_set::container;
        long n1 = matcl::details::no_existing_objects();
        
        const container& mc_c = m_ms.get(c,1);
        const container& mc_r = m_ms.get(r,1);

        n1 = matcl::details::no_existing_objects() - n1;
        m_new_objects = n1;

        size_t size_r   = mc_r.size();
        size_t size_c   = mc_c.size();
        size_t size     = std::min(size_r, size_c);

        Real dif            = 0.0;

        for (size_t i = 0; i < size ; i ++ )
        {
            Matrix D_r      = mc_r[i];
            Matrix D_c      = mc_c[i];
            Matrix res1     = scale_rowscols(mat, D_r, D_c);
            Matrix res2     = bdiag(D_r) * convert_value(mat,res1.get_value_code()) * bdiag(D_c);
                
            check_struct(res1);
            
            if (check_value_code(res1.get_value_code(),res2.get_value_code()) == false)
                dif         += 1;

            Real loc_dif    = norm_1(res1 - res2);

            if (loc_dif < norm_1(res1) * constants::eps(res1.get_value_code()) * res1.numel() * 10.)
                loc_dif     = 0.0;

            dif             += loc_dif;
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
Real test_function_scale_rowscols::eval_scalar(const Scalar& s,bool,int code  )
{
    (void)code;
    (void)s;
    return 0.0;
};

bool test_function_scale_rowscols::check_value_code(value_code v1, value_code v2) const
{
    value_code v1_r = matrix_traits::real_value_type(v1);
    value_code v2_r = matrix_traits::real_value_type(v2);

    return v1_r == v2_r;
}

};};
