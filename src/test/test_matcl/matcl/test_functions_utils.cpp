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
#include "test_functions_utils.h"

#include "test/test_matcl/framework/matrix_set/matrix_set_1.h"
#include "matcl-core/IO/logger.h"
#include "test/test_matcl/framework/data_struct/test_selector.h"

#include <boost/thread.hpp>


namespace matcl { namespace test
{

void test_utils_st(const rand_matrix_ptr& rand)
{
    (void)rand;

    test::options opts;

    try
    {
        opts.show_partial_res   = 0;
        opts.first_matrix_code  = 0;
        opts.show_memleaks      = true;

        utils_functions_list tf;
        
        tf.make();
    }
    catch(const std::exception& ex)
    {
        matcl::out_stream <<ex.what();
    };
};

void test_utils_mt(const rand_matrix_ptr& rand)
{
    (void)rand;

    test::options opts;
    try
    {
        opts.show_partial_res   = 0;
        opts.first_matrix_code  = 0;
        opts.show_memleaks      = false;

        boost::thread_group tg;

        for (int i = 0; i < 10; i++)
            tg.create_thread(utils_functions_list());

        tg.join_all();
    }
    catch(const std::exception& ex)
    {
        matcl::out_stream <<ex.what();
    };
};

utils_functions_list::utils_functions_list()
{};

void utils_functions_list::make()
{
    SELECT_TEST /**/ (3, test_thread());
    SELECT_TEST /**/ (3, test_tuple());
};

void utils_functions_list::test_thread()
{
    Real out = 0.;
    test_function_thread tf;
    out += tf.make();

    if (out == 0.)
        matcl::out_stream << std::string() +   "threads: OK" + "\n";
    else
        matcl::out_stream << std::string() +   "threads: FAILED"  + "\n";
};

void utils_functions_list::test_tuple()
{
    Real out = 0.;
    test_function_tuple tf;
    out += tf.make();

    if (out == 0.)
        matcl::out_stream << std::string() +   "tuple: OK" + "\n";
    else
        matcl::out_stream << std::string() +   "tuple: FAILED"  + "\n";
};

//
class test_thread_cl
{
    private:
        const std::vector<Matrix>&  mat_vec;
        const std::vector<Matrix>&  mat_val;
        Integer                     thread_id;
        Real*                       result;

    public:
        test_thread_cl(const std::vector<Matrix>& mat, const std::vector<Matrix>& val, 
                       Integer id, Real& res)	
            : mat_vec(mat), mat_val(val), thread_id(id), result(&res)
        {} ;

        test_thread_cl(test_thread_cl&& other)	
            : mat_vec(other.mat_vec), mat_val(other.mat_val), thread_id(other.thread_id) 
            , result(other.result)
        {} ;
        test_thread_cl(const test_thread_cl& other)	
            : mat_vec(other.mat_vec), mat_val(other.mat_val), thread_id(other.thread_id) 
            , result(other.result)
        {} ;

        void make()
        {
            *result = 0;
            for (int i = 0; i < 10000; i++)
            {
                Integer pos = abs(irand())%mat_vec.size();			
                Matrix tmp = Matrix(mat_vec[pos]);
                make(tmp,100);
                Matrix out = sum(sum(tmp,1),2);
                *result += norm_1(out - mat_val[pos]);
            };

            matcl::mutex                m1;
            matcl::spinlock_mutex       m2;
            matcl::nonblocking_mutex    m3;

            {
                std::lock_guard<matcl::mutex> lock(m1);
            }
            {
                std::lock_guard<spinlock_mutex> lock(m2);
            }
            {
                std::lock_guard<matcl::nonblocking_mutex> lock(m3);
            };
        };

        void operator()()
        {
            make();
        };

    private:
        void make(Matrix tmp, int n)
        {
            if (n > 0)
            {
                Matrix tmp2 = Matrix(tmp);
                make(tmp2,n-1);
            };
        };
    
        test_thread_cl& operator=(const test_thread_cl&) = delete;
};

Real test_function_thread::make()
{
    try
    {
        boost::thread_group tg;

        std::vector<Matrix> mat_vec(3);
        std::vector<Matrix> mat_val(3);
      
        std::vector<Real> results(20);

        for (int i = 0; i < 3; i++)
        {
            mat_vec[i] = randn(4,4);
            mat_val[i] = sum(sum(mat_vec[i],1),2);
        };

        for (int i = 0; i < 20; i++)
        {
            test_thread_cl thr(mat_vec,mat_val,i, results[i]);
            tg.create_thread(std::move(thr));
        };

        tg.join_all();

        Real dif    = 0;
        for (int i = 0; i < 20; ++i)
            dif     += results[i];

        return dif;
    }
    catch(const std::exception& ex)
    {
        (void)ex;
        return 1.0;
    }
    catch(...)
    {
        return 2.;
    };
};

Real test_function_tuple::make()
{
    try
    {
        using tuple     = matcl::tuple<int,int>;
        using std_tuple = std::tuple<int,int>;

        tuple t(1,2);
        std_tuple st1   = t.to_standard_tuple();
        std_tuple&& st2 = tuple(1,2).to_standard_tuple();

        tuple t1        = matcl::from_standard_tuple(st1);
        tuple t2        = matcl::from_standard_tuple(std::move(st2));

        return 0.0;
    }
    catch(std::exception& ex)
    {
        (void)ex;
        return 1.0;
    }
    catch(...)
    {
        return 2.0;
    }
};

};};
