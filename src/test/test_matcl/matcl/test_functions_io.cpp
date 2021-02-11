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
#include "test_functions_io.h"

#include "test/test_matcl/framework/matrix_set/matrix_set_bin_1.h"
#include "test/test_matcl/framework/data_struct/test_selector.h"
#include "matcl-core/IO/logger.h"

#include "matcl-matrep/matcl_matrep.h"
#include "matcl-file/matcl_file.h"

#include <boost/thread.hpp>

namespace matcl { namespace test
{

class test_io
{
    io_functions_list&			    tf;
    const test::options&			opts;
    Integer                         thread_id;

    public:
        test_io(io_functions_list& tf, const test::options& opts, Integer id)
            :tf(tf),opts(opts),thread_id(id)
        {};

        test_io(const test_io& ta)
            :tf(ta.tf),opts(ta.opts), thread_id(ta.thread_id)
        {};

        void make()
        {
            using matrix_pair   = test::matrix_set_bin::matrix_pair;
            using scalar_pair   = test::matrix_set_bin::scalar_pair;

            /*
            matrix_pair mp = tf.get_matrix(335);  
            
            Matrix mat1     = mp.first;
            Matrix mat2     = mp.second;

            disp(mat1);
            disp(mat2);

            Real tol1   = norm_1(mat1) * constants::eps(mat1.get_value_code());
            Real tol2   = norm_1(mat2) * constants::eps(mat2.get_value_code());
            Real dif    = 0;

            {
                std::stringstream ss;
                mm_save(ss, mat1);
                mm_save(ss, mat2);

                ss.seekg(std::ios_base::beg);

                disp(ss.str());

                Matrix tmp1,tmp2;
                mm_load(ss, tmp1);
                mm_load(ss, tmp2);

                disp(tmp1);
                disp(tmp2);

                check_struct(tmp1);
                check_struct(tmp2);

                Real dif1 = norm_1(mat1 - tmp1);
                Real dif2 = norm_1(mat2 - tmp2);

                if (dif1 < tol1)
                    dif1    = 0;
                if (dif2 < tol2)
                    dif2    = 0;

                dif         += dif1;
                dif         += dif2;

                tmp1        = Matrix();
                tmp2        = Matrix();
            };

            {
                std::stringstream ss;
                std::string comm = "test";
            
                mm_save(ss, mat1, comm);
                mm_save(ss, mat2, comm);
                mm_save(ss, mat2);

                ss.seekg(std::ios_base::beg);

                Matrix tmp1, tmp2, tmp3;
                std::string com1, com2;
            
                mm_load(ss, tmp1, com1);
                mm_load(ss, tmp2, com2);
                mm_load(ss, tmp3);

                check_struct(tmp1);
                check_struct(tmp2);
                check_struct(tmp3);

                Real dif1   = norm_1(mat1 - tmp1);
                Real dif2   = norm_1(mat2 - tmp2);
                Real dif3   = norm_1(mat2 - tmp3);

                if (dif1 < tol1)
                    dif1    = 0;
                if (dif2 < tol2)
                    dif2    = 0;
                if (dif3 < tol2)
                    dif3    = 0;

                dif         += dif1;
                dif         += dif2;
                dif         += dif3;

                if (com1 != comm)
                    dif += 1;
                if (com2 != comm)
                    dif += 1;

                {
                    ss.seekg(std::ios_base::beg);
                    std::stringstream os_matcl;

                    disp(ss.str());
                    convert_mm_to_matcl(ss, os_matcl);
                    convert_mm_to_matcl(ss, os_matcl);
                    convert_mm_to_matcl(ss, os_matcl);

                    Matrix tmp1, tmp2, tmp3;
                    std::string com1, com2;

                    os_matcl.seekg(std::ios_base::beg);
                    disp(os_matcl.str());

                    load(os_matcl, tmp1, com1);
                    load(os_matcl, tmp2, com2);
                    load(os_matcl, tmp3);

                    check_struct(tmp1);
                    check_struct(tmp2);
                    check_struct(tmp3);

                    Real dif1   = norm_1(mat1 - tmp1);
                    Real dif2   = norm_1(mat2 - tmp2);
                    Real dif3   = norm_1(mat2 - tmp3);

                    if (dif1 < tol1)
                        dif1    = 0;
                    if (dif2 < tol2)
                        dif2    = 0;
                    if (dif3 < tol2)
                        dif3    = 0;

                    dif         += dif1;
                    dif         += dif2;
                    dif         += dif3;

                    if (com1 != comm)
                        dif += 1;
                    if (com2 != comm)
                        dif += 1;
                };
            }

            {
                std::stringstream ss;
                std::string comm = "test";
            
                save(ss, mat1, comm);
                save(ss, mat2, comm);
                save(ss, mat2);

                ss.seekg(std::ios_base::beg);

                Matrix tmp1, tmp2, tmp3;
                std::string com1, com2;

                std::stringstream ss3;
                convert_matcl_to_mm(ss, ss3);
                convert_matcl_to_mm(ss, ss3);
                convert_matcl_to_mm(ss, ss3);

                ss3.seekg(std::ios_base::beg);

                mm_load(ss3, tmp1, com1);
                mm_load(ss3, tmp2, com2);
                mm_load(ss3, tmp3);

                check_struct(tmp1);
                check_struct(tmp2);
                check_struct(tmp3);

                Real dif1   = norm_1(mat1 - tmp1);
                Real dif2   = norm_1(mat2 - tmp2);
                Real dif3   = norm_1(mat2 - tmp3);

                if (dif1 < tol1)
                    dif1    = 0;
                if (dif2 < tol2)
                    dif2    = 0;
                if (dif3 < tol2)
                    dif3    = 0;

                dif         += dif1;
                dif         += dif2;
                dif         += dif3;

                if (com1 != comm)
                    dif += 1;
                if (com2 != comm)
                    dif += 1;
            };
            */

            tf.make(opts,thread_id);
        };
        void operator()()
        {
            make();
        };

    private:		
        test_io& operator=(const test_io&) = delete;
};

void test_io_st(const rand_matrix_ptr& rand)
{
    test::options opts;

    try
    {
        opts.show_partial_res   = 0;
        opts.first_matrix_code  = 0;
        opts.show_memleaks      = true;

        mat_set_bin_1 ms1(rand);
        io_functions_list tf(ms1);

        test_io ta(tf,opts,0);
        ta.make();
    }
    catch(const std::exception& ex)
    {
        matcl::out_stream <<ex.what();
    };
};

void test_io_mt(const rand_matrix_ptr& rand)
{
    test::options opts;
    try
    {
        opts.show_partial_res   = 0;
        opts.first_matrix_code  = 0;
        opts.show_memleaks      = false;

        mat_set_bin_2 ms1(rand);
        io_functions_list tf(ms1);

        boost::thread_group tg;

        for (int i = 0; i < 10; i++)
            tg.create_thread(test_io(tf,opts, i));

        tg.join_all();
    }
    catch(const std::exception& ex)
    {
        matcl::out_stream <<ex.what();
    };
};

void io_functions_list::make(options opts, Integer thread_id)
{
    (void)thread_id;

    m_options = opts;
    
    SELECT_TEST (3, test_io_mm());
    SELECT_TEST (3, test_io());
    SELECT_TEST (3, test_io2());
    SELECT_TEST (3, test_serialize());    
};

io_functions_list::matrix_pair io_functions_list::get_matrix(int code) const
{
    return m_tests.get_matrix(code);
};
io_functions_list::scalar_pair io_functions_list::get_scalar(int code) const
{
    return m_tests.get_scalar(code);
};

void io_functions_list::test_io()
{
    Real out = 0.;
    test_function_io_formatted tf;

    out += m_tests.make(&tf,m_options);

    if (out == 0.)
        matcl::out_stream << std::string() +   "io_format: OK" + "\n";
    else
        matcl::out_stream << std::string() +   "io_format: FAILED"  + "\n";
};

void io_functions_list::test_io2()
{
    Real out = 0.;
    test_function_io_formatted2 tf;

    out += m_tests.make(&tf,m_options);

    if (out == 0.)
        matcl::out_stream << std::string() +   "io_format2: OK" + "\n";
    else
        matcl::out_stream << std::string() +   "io_format2: FAILED"  + "\n";
};

void io_functions_list::test_io_mm()
{
    Real out = 0.;
    test_function_io_mm tf;

    out += m_tests.make(&tf,m_options);

    if (out == 0.)
        matcl::out_stream << std::string() +   "io_mm: OK" + "\n";
    else
        matcl::out_stream << std::string() +   "io_mm: FAILED"  + "\n";
};

void io_functions_list::test_serialize()
{
    Real out = 0.;
    test_function_serialize tf;
    out += m_tests.make(&tf,m_options);

    if (out == 0.)
        matcl::out_stream << std::string() +   "serialize: OK" + "\n";
    else
        matcl::out_stream << std::string() +   "serialize: FAILED"  + "\n";
};

Real test_function_serialize::eval_mat(const Matrix& mat1,const Matrix& mat2, int )
{
    try
    {	
        Real dif = 0;

        {
            std::ostringstream ss;
            oarchive ia(ss);

            save(ia,mat1);
            save(ia,mat2);

            Matrix tmp1, tmp2;
            std::istringstream ss2(ss.str());
            iarchive ia2(ss2);

            load(ia2,tmp1);
            load(ia2,tmp2);

            check_struct(tmp1);
            check_struct(tmp2);
            
            dif += norm_1(mat1 - tmp1);
            dif += norm_1(mat2 - tmp2);
        };
        {
            std::ostringstream ss;
            oarchive ia(ss);

            save(ia,mat1);
            save(ia,mat2);
            save(ia,mat1);
            save(ia,mat2);

            Matrix tmp1, tmp2, tmp3, tmp4;
            std::istringstream ss2(ss.str());
            iarchive ia2(ss2);

            load(ia2,tmp1);
            load(ia2,tmp2);
            load(ia2,tmp3);
            load(ia2,tmp4);

            check_struct(tmp1);
            check_struct(tmp2);
            check_struct(tmp3);
            check_struct(tmp4);
            
            dif += norm_1(mat1 - tmp1);
            dif += norm_1(mat2 - tmp2);
            dif += norm_1(mat1 - tmp3);
            dif += norm_1(mat2 - tmp4);
        };

        //if (abs(dif)<1e-10)
        //    dif = 0;

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
Real test_function_serialize::eval_scalar(const Scalar& s1, const Scalar& s2, int )
{
    Real out = eval_scalar_impl<test_function_serialize>(s1, s2);
    return out;
};

Real test_function_io_formatted::eval_mat(const Matrix& mat1,const Matrix& mat2, int )
{
    try
    {	
        Real dif = 0;

        {
            std::ostringstream ss;
            ss << mat1;
            ss << mat2;

            Matrix tmp1,tmp2;
            std::istringstream ss2(ss.str());
            ss2 >> tmp1;
            ss2 >> tmp2;

            check_struct(tmp1);
            check_struct(tmp2);

            dif += norm_1(mat1 - tmp1);
            dif += norm_1(mat2 - tmp2);
        };

        //if (abs(dif)<1e-10)
        //    dif = 0;

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

Real test_function_io_formatted::eval_scalar(const Scalar& s1, const Scalar& s2, int )
{
    Real out = eval_scalar_impl<test_function_io_formatted>(s1, s2);
    return out;
};

Real test_function_io_formatted2::eval_mat(const Matrix& mat1,const Matrix& mat2, int )
{
    try
    {	
        Real dif = 0;

        {
            std::ostringstream ss;
            save(ss, mat1);
            save(ss, mat2);

            Matrix tmp1,tmp2;
            std::istringstream ss2(ss.str());
            load(ss2, tmp1);
            load(ss2, tmp2);

            check_struct(tmp1);
            check_struct(tmp2);

            dif += norm_1(mat1 - tmp1);
            dif += norm_1(mat2 - tmp2);
        };

        {
            std::ostringstream ss;
            std::string comm = "test";
            save(ss, mat1, comm);
            save(ss, mat2, comm);
            save(ss, mat2);

            Matrix tmp1, tmp2, tmp3;
            std::string com1, com2;
            std::istringstream ss2(ss.str());
            load(ss2, tmp1, com1);
            load(ss2, tmp2, com2);
            load(ss2, tmp3);

            check_struct(tmp1);
            check_struct(tmp2);
            check_struct(tmp3);

            dif += norm_1(mat1 - tmp1);
            dif += norm_1(mat2 - tmp2);
            dif += norm_1(mat2 - tmp3);

            if (com1 != comm)
                dif += 1;
            if (com2 != comm)
                dif += 1;
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

Real test_function_io_formatted2::eval_scalar(const Scalar&, const Scalar&, int )
{
    return 0.0;
};

Real test_function_io_mm::eval_mat(const Matrix& mat1,const Matrix& mat2, int )
{
    try
    {	
        Real dif = 0;

        Real tol1   = norm_1(mat1) * constants::eps(mat1.get_value_code());
        Real tol2   = norm_1(mat2) * constants::eps(mat2.get_value_code());

        {
            std::stringstream ss;
            mm_save(ss, mat1);
            mm_save(ss, mat2);

            ss.seekg(std::ios_base::beg);

            Matrix tmp1,tmp2;
            mm_load(ss, tmp1);
            mm_load(ss, tmp2);

            check_struct(tmp1);
            check_struct(tmp2);

            Real dif1 = norm_1(mat1 - tmp1);
            Real dif2 = norm_1(mat2 - tmp2);

            if (dif1 < tol1)
                dif1    = 0;
            if (dif2 < tol2)
                dif2    = 0;

            dif         += dif1;
            dif         += dif2;
        };

        {
            std::stringstream ss;
            std::string comm = "test";
            
            mm_save(ss, mat1, comm);
            mm_save(ss, mat2, comm);
            mm_save(ss, mat2);

            ss.seekg(std::ios_base::beg);

            Matrix tmp01, tmp02, tmp03;
            std::string com01, com02;
            
            mm_load(ss, tmp01, com01);
            mm_load(ss, tmp02, com02);
            mm_load(ss, tmp03);

            check_struct(tmp01);
            check_struct(tmp02);
            check_struct(tmp03);

            Real dif01   = norm_1(mat1 - tmp01);
            Real dif02   = norm_1(mat2 - tmp02);
            Real dif03   = norm_1(mat2 - tmp03);

            if (dif01 < tol1)
                dif01    = 0;
            if (dif02 < tol2)
                dif02    = 0;
            if (dif03 < tol2)
                dif03    = 0;

            dif         += dif01;
            dif         += dif02;
            dif         += dif03;

            if (com01 != comm)
                dif += 1;
            if (com02 != comm)
                dif += 1;

            {
                ss.seekg(std::ios_base::beg);
                std::stringstream os_matcl;

                convert_mm_to_matcl(ss, os_matcl);
                convert_mm_to_matcl(ss, os_matcl);
                convert_mm_to_matcl(ss, os_matcl);

                Matrix tmp1, tmp2, tmp3;
                std::string com1, com2;

                os_matcl.seekg(std::ios_base::beg);
                load(os_matcl, tmp1, com1);
                load(os_matcl, tmp2, com2);
                load(os_matcl, tmp3);

                check_struct(tmp1);
                check_struct(tmp2);
                check_struct(tmp3);

                Real dif1   = norm_1(mat1 - tmp1);
                Real dif2   = norm_1(mat2 - tmp2);
                Real dif3   = norm_1(mat2 - tmp3);

                if (dif1 < tol1)
                    dif1    = 0;
                if (dif2 < tol2)
                    dif2    = 0;
                if (dif3 < tol2)
                    dif3    = 0;

                dif         += dif1;
                dif         += dif2;
                dif         += dif3;

                if (com1 != comm)
                    dif += 1;
                if (com2 != comm)
                    dif += 1;
            };
        }

        {
            std::stringstream ss;
            std::string comm = "test";
            
            save(ss, mat1, comm);
            save(ss, mat2, comm);
            save(ss, mat2);

            ss.seekg(std::ios_base::beg);

            Matrix tmp1, tmp2, tmp3;
            std::string com1, com2;

            std::stringstream ss3;
            convert_matcl_to_mm(ss, ss3);
            convert_matcl_to_mm(ss, ss3);
            convert_matcl_to_mm(ss, ss3);

            ss3.seekg(std::ios_base::beg);

            mm_load(ss3, tmp1, com1);
            mm_load(ss3, tmp2, com2);
            mm_load(ss3, tmp3);

            check_struct(tmp1);
            check_struct(tmp2);
            check_struct(tmp3);

            Real dif1   = norm_1(mat1 - tmp1);
            Real dif2   = norm_1(mat2 - tmp2);
            Real dif3   = norm_1(mat2 - tmp3);

            if (dif1 < tol1)
                dif1    = 0;
            if (dif2 < tol2)
                dif2    = 0;
            if (dif3 < tol2)
                dif3    = 0;

            dif         += dif1;
            dif         += dif2;
            dif         += dif3;

            if (com1 != comm)
                dif += 1;
            if (com2 != comm)
                dif += 1;
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

Real test_function_io_mm::eval_scalar(const Scalar&, const Scalar&, int )
{
    return 0.0;
};

};};
