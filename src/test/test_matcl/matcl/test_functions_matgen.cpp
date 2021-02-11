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
#include "test_functions_matgen.h"

#include "test/test_matcl/framework/matrix_set/matrix_set_1.h"
#include "matcl-core/IO/logger.h"
#include "test/test_matcl/framework/data_struct/test_selector.h"

#include <boost/thread.hpp>


namespace matcl { namespace test
{

void test_matgen_st(const rand_matrix_ptr& rand)
{
    (void)rand;

    test::options opts;

    try
    {
        opts.show_partial_res   = 0;
        opts.first_matrix_code  = 0;
        opts.show_memleaks      = true;

        matgen_functions_list tf;
        
        tf.make();
    }
    catch(const std::exception& ex)
    {
        matcl::out_stream <<ex.what();
    };
};

void test_matgen_mt(const rand_matrix_ptr& rand)
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
            tg.create_thread(matgen_functions_list());

        tg.join_all();
    }
    catch(const std::exception& ex)
    {
        matcl::out_stream <<ex.what();
    };
};

matgen_functions_list::matgen_functions_list()
{};

void matgen_functions_list::make()
{
    SELECT_TEST (3, test_seq());
    SELECT_TEST (3, test_zeros());
    SELECT_TEST (3, test_ones());
    SELECT_TEST (3, test_eye());
    SELECT_TEST (3, test_rand());
    SELECT_TEST (3, test_constructors());
};

void matgen_functions_list::test_seq()
{
    Real out = 0.;
    test_function_seq tf;
    out += tf.make();

    if (out == 0.)
        matcl::out_stream << std::string() +   "sequences: OK" + "\n";
    else
        matcl::out_stream << std::string() +   "sequences: FAILED"  + "\n";
};

void matgen_functions_list::test_zeros()
{
    Real out = 0.;
    test_function_zeros tf;
    out += tf.make();

    if (out == 0.)
        matcl::out_stream << std::string() +   "zeros: OK" + "\n";
    else
        matcl::out_stream << std::string() +   "zeros: FAILED"  + "\n";
};

void matgen_functions_list::test_ones()
{
    Real out = 0.;
    test_function_ones tf;
    out += tf.make();

    if (out == 0.)
        matcl::out_stream << std::string() +   "ones: OK" + "\n";
    else
        matcl::out_stream << std::string() +   "ones: FAILED"  + "\n";
};

void matgen_functions_list::test_eye()
{
    Real out = 0.;
    test_function_eye tf;
    out += tf.make();

    if (out == 0.)
        matcl::out_stream << std::string() +   "eye: OK" + "\n";
    else
        matcl::out_stream << std::string() +   "eye: FAILED"  + "\n";
};

void matgen_functions_list::test_rand()
{
    Real out = 0.;
    test_function_rand tf;
    out += tf.make();

    if (out == 0.)
        matcl::out_stream << std::string() +   "rand: OK" + "\n";
    else
        matcl::out_stream << std::string() +   "rand: FAILED"  + "\n";
};

void matgen_functions_list::test_constructors()
{
    Real out = 0.;
    test_function_construct tf;
    out += tf.make();

    if (out == 0.)
        matcl::out_stream << std::string() +   "constructors: OK" + "\n";
    else
        matcl::out_stream << std::string() +   "constructors: FAILED"  + "\n";
};

Real test_function_construct::make()
{
    Real out = 0.;
    try
    {
        Complex I	        = constants::i();
        Float_complex If	= constants::f_i();

        //
        out			+= norm_1(make_integer_dense(2,2)-zeros(2,2));
        out			+= norm_1(make_integer_dense(0,2)-zeros(0,2));
        out			+= norm_1(make_integer_dense(2,0)-zeros(2,0));
        out			+= norm_1(make_integer_dense(0,0)-zeros(0,0));
        out			+= norm_1(make_integer_dense(1,1)-zeros(1,1));

        out			+= norm_1(make_real_dense(2,2)-zeros(2,2));
        out			+= norm_1(make_real_dense(0,2)-zeros(0,2));
        out			+= norm_1(make_real_dense(2,0)-zeros(2,0));
        out			+= norm_1(make_real_dense(0,0)-zeros(0,0));
        out			+= norm_1(make_real_dense(1,1)-zeros(1,1));

        out			+= norm_1(make_float_dense(2,2)-zeros(2,2));
        out			+= norm_1(make_float_dense(0,2)-zeros(0,2));
        out			+= norm_1(make_float_dense(2,0)-zeros(2,0));
        out			+= norm_1(make_float_dense(0,0)-zeros(0,0));
        out			+= norm_1(make_float_dense(1,1)-zeros(1,1));

        out			+= norm_1(make_complex_dense(2,2)-zeros(2,2));
        out			+= norm_1(make_complex_dense(0,2)-zeros(0,2));
        out			+= norm_1(make_complex_dense(2,0)-zeros(2,0));
        out			+= norm_1(make_complex_dense(0,0)-zeros(0,0));
        out			+= norm_1(make_complex_dense(1,1)-zeros(1,1));

        out			+= norm_1(make_float_complex_dense(2,2)-zeros(2,2));
        out			+= norm_1(make_float_complex_dense(0,2)-zeros(0,2));
        out			+= norm_1(make_float_complex_dense(2,0)-zeros(2,0));
        out			+= norm_1(make_float_complex_dense(0,0)-zeros(0,0));
        out			+= norm_1(make_float_complex_dense(1,1)-zeros(1,1));

        //
        out			+= norm_1(make_integer_dense(1,2,2)-ones(2,2));
        out			+= norm_1(make_integer_dense(1,0,2)-ones(0,2));
        out			+= norm_1(make_integer_dense(1,2,0)-ones(2,0));
        out			+= norm_1(make_integer_dense(1,0,0)-ones(0,0));
        out			+= norm_1(make_integer_dense(1,1,1)-ones(1,1));

        out			+= norm_1(make_real_dense(1.,2,2)-ones(2,2));
        out			+= norm_1(make_real_dense(1.,0,2)-ones(0,2));
        out			+= norm_1(make_real_dense(1.,2,0)-ones(2,0));
        out			+= norm_1(make_real_dense(1.,0,0)-ones(0,0));
        out			+= norm_1(make_real_dense(1.,1,1)-ones(1,1));

        out			+= norm_1(make_float_dense(1.,2,2)-ones(2,2));
        out			+= norm_1(make_float_dense(1.,0,2)-ones(0,2));
        out			+= norm_1(make_float_dense(1.,2,0)-ones(2,0));
        out			+= norm_1(make_float_dense(1.,0,0)-ones(0,0));
        out			+= norm_1(make_float_dense(1.,1,1)-ones(1,1));

        out			+= norm_1(make_complex_dense(1+0*I,2,2)-ones(2,2));
        out			+= norm_1(make_complex_dense(1+0*I,0,2)-ones(0,2));
        out			+= norm_1(make_complex_dense(1+0*I,2,0)-ones(2,0));
        out			+= norm_1(make_complex_dense(1+0*I,0,0)-ones(0,0));
        out			+= norm_1(make_complex_dense(1+0*I,1,1)-ones(1,1));

        out			+= norm_1(make_float_complex_dense(1.f+0.f*If,2,2)-ones(2,2));
        out			+= norm_1(make_float_complex_dense(1.f+0.f*If,0,2)-ones(0,2));
        out			+= norm_1(make_float_complex_dense(1.f+0.f*If,2,0)-ones(2,0));
        out			+= norm_1(make_float_complex_dense(1.f+0.f*If,0,0)-ones(0,0));
        out			+= norm_1(make_float_complex_dense(1.f+0.f*If,1,1)-ones(1,1));

        for (Integer nz = 0; nz < 3; ++nz)
        {
            out			+= norm_1(make_integer_sparse(2,2,nz)-zeros(2,2));
            out			+= norm_1(make_integer_sparse(0,2,nz)-zeros(0,2));
            out			+= norm_1(make_integer_sparse(2,0,nz)-zeros(2,0));
            out			+= norm_1(make_integer_sparse(0,0,nz)-zeros(0,0));
            out			+= norm_1(make_integer_sparse(1,1,nz)-zeros(1,1));

            out			+= norm_1(make_real_sparse(2,2,nz)-zeros(2,2));
            out			+= norm_1(make_real_sparse(0,2,nz)-zeros(0,2));
            out			+= norm_1(make_real_sparse(2,0,nz)-zeros(2,0));
            out			+= norm_1(make_real_sparse(0,0,nz)-zeros(0,0));
            out			+= norm_1(make_real_sparse(1,1,nz)-zeros(1,1));

            out			+= norm_1(make_float_sparse(2,2,nz)-zeros(2,2));
            out			+= norm_1(make_float_sparse(0,2,nz)-zeros(0,2));
            out			+= norm_1(make_float_sparse(2,0,nz)-zeros(2,0));
            out			+= norm_1(make_float_sparse(0,0,nz)-zeros(0,0));
            out			+= norm_1(make_float_sparse(1,1,nz)-zeros(1,1));

            out			+= norm_1(make_complex_sparse(2,2,nz)-zeros(2,2));
            out			+= norm_1(make_complex_sparse(0,2,nz)-zeros(0,2));
            out			+= norm_1(make_complex_sparse(2,0,nz)-zeros(2,0));
            out			+= norm_1(make_complex_sparse(0,0,nz)-zeros(0,0));
            out			+= norm_1(make_complex_sparse(1,1,nz)-zeros(1,1));

            out			+= norm_1(make_float_complex_sparse(2,2,nz)-zeros(2,2));
            out			+= norm_1(make_float_complex_sparse(0,2,nz)-zeros(0,2));
            out			+= norm_1(make_float_complex_sparse(2,0,nz)-zeros(2,0));
            out			+= norm_1(make_float_complex_sparse(0,0,nz)-zeros(0,0));
            out			+= norm_1(make_float_complex_sparse(1,1,nz)-zeros(1,1));
        };

        //
        Matrix B,F,d1,d2;
        B			= make_integer_band(0,2,2,0,0);
        out			+= norm_1(B-zeros(2,2));
        B			= make_integer_band(0,2,2,-1,0);
        out			+= norm_1(B-zeros(2,2));
        B			= make_integer_band(0,2,2,0,1);
        out			+= norm_1(B-zeros(2,2));
        B			= make_integer_band(0,2,2,-1,1);
        out			+= norm_1(B-zeros(2,2));
        B			= make_integer_band(0,2,2,-2,2);
        out			+= norm_1(B-zeros(2,2));

        B			= make_integer_band(1,2,2,0,0);	tie(d1,d2,F) = find3(B);
        out			+= norm_1(F-1);
        B			= make_integer_band(1,2,2,-1,0);	tie(d1,d2,F) = find3(B);
        out			+= norm_1(F-1);
        B			= make_integer_band(1,2,2,0,1);	tie(d1,d2,F) = find3(B);
        out			+= norm_1(F-1);
        B			= make_integer_band(1,2,2,-1,1);	tie(d1,d2,F) = find3(B);
        out			+= norm_1(F-1);
        out			+= norm_1(B-1);
        B			= make_integer_band(1,2,2,-3,3);	tie(d1,d2,F) = find3(B);
        out			+= norm_1(F-1);

        //
        B			= make_real_band(0,2,2,0,0);
        out			+= norm_1(B-zeros(2,2));
        B			= make_real_band(0,2,2,-1,0);
        out			+= norm_1(B-zeros(2,2));
        B			= make_real_band(0,2,2,0,1);
        out			+= norm_1(B-zeros(2,2));
        B			= make_real_band(0,2,2,-1,1);
        out			+= norm_1(B-zeros(2,2));
        B			= make_real_band(0,2,2,-3,3);
        out			+= norm_1(B-zeros(2,2));

        B			= make_real_band(1,2,2,0,0);	tie(d1,d2,F) = find3(B);
        out			+= norm_1(F-1);
        B			= make_real_band(1,2,2,-1,0);	tie(d1,d2,F) = find3(B);
        out			+= norm_1(F-1);
        B			= make_real_band(1,2,2,0,1);	tie(d1,d2,F) = find3(B);
        out			+= norm_1(F-1);
        B			= make_real_band(1,2,2,-1,1);	tie(d1,d2,F) = find3(B);
        out			+= norm_1(F-1);
        out			+= norm_1(B-1);
        B			= make_real_band(1,2,2,-3,3);	tie(d1,d2,F) = find3(B);
        out			+= norm_1(F-1);

        //
        B			= make_float_band(0,2,2,0,0);
        out			+= norm_1(B-zeros(2,2));
        B			= make_float_band(0,2,2,-1,0);
        out			+= norm_1(B-zeros(2,2));
        B			= make_float_band(0,2,2,0,1);
        out			+= norm_1(B-zeros(2,2));
        B			= make_float_band(0,2,2,-1,1);
        out			+= norm_1(B-zeros(2,2));
        B			= make_float_band(0,2,2,-3,3);
        out			+= norm_1(B-zeros(2,2));

        B			= make_float_band(1,2,2,0,0);	tie(d1,d2,F) = find3(B);
        out			+= norm_1(F-1);
        B			= make_float_band(1,2,2,-1,0);	tie(d1,d2,F) = find3(B);
        out			+= norm_1(F-1);
        B			= make_float_band(1,2,2,0,1);	tie(d1,d2,F) = find3(B);
        out			+= norm_1(F-1);
        B			= make_float_band(1,2,2,-1,1);	tie(d1,d2,F) = find3(B);
        out			+= norm_1(F-1);
        out			+= norm_1(B-1);
        B			= make_float_band(1,2,2,-3,3);	tie(d1,d2,F) = find3(B);
        out			+= norm_1(F-1);

        //
        B			= make_complex_band(0,2,2,0,0);
        out			+= norm_1(B-zeros(2,2));
        B			= make_complex_band(0,2,2,-1,0);
        out			+= norm_1(B-zeros(2,2));
        B			= make_complex_band(0,2,2,0,1);
        out			+= norm_1(B-zeros(2,2));
        B			= make_complex_band(0,2,2,-1,1);
        out			+= norm_1(B-zeros(2,2));
        B			= make_complex_band(0,2,2,-3,3);
        out			+= norm_1(B-zeros(2,2));

        B			= make_complex_band(1,2,2,0,0);	tie(d1,d2,F) = find3(B);
        out			+= norm_1(F-1);
        B			= make_complex_band(1,2,2,-1,0);	tie(d1,d2,F) = find3(B);
        out			+= norm_1(F-1);
        B			= make_complex_band(1,2,2,0,1);	tie(d1,d2,F) = find3(B);
        out			+= norm_1(F-1);
        B			= make_complex_band(1,2,2,-1,1);	tie(d1,d2,F) = find3(B);
        out			+= norm_1(F-1);
        out			+= norm_1(B-1);
        B			= make_complex_band(1,2,2,-3,3);	tie(d1,d2,F) = find3(B);
        out			+= norm_1(F-1);

        //
        B			= make_float_complex_band(0,2,2,0,0);
        out			+= norm_1(B-zeros(2,2));
        B			= make_float_complex_band(0,2,2,-1,0);
        out			+= norm_1(B-zeros(2,2));
        B			= make_float_complex_band(0,2,2,0,1);
        out			+= norm_1(B-zeros(2,2));
        B			= make_float_complex_band(0,2,2,-1,1);
        out			+= norm_1(B-zeros(2,2));
        B			= make_float_complex_band(0,2,2,-3,3);
        out			+= norm_1(B-zeros(2,2));

        B			= make_float_complex_band(1,2,2,0,0);	tie(d1,d2,F) = find3(B);
        out			+= norm_1(F-1);
        B			= make_float_complex_band(1,2,2,-1,0);	tie(d1,d2,F) = find3(B);
        out			+= norm_1(F-1);
        B			= make_float_complex_band(1,2,2,0,1);	tie(d1,d2,F) = find3(B);
        out			+= norm_1(F-1);
        B			= make_complex_band(1,2,2,-1,1);	tie(d1,d2,F) = find3(B);
        out			+= norm_1(F-1);
        out			+= norm_1(B-1);
        B			= make_float_complex_band(1,2,2,-3,3);	tie(d1,d2,F) = find3(B);
        out			+= norm_1(F-1);

        //
        B			= make_integer_band(2,2,0,0);
        out			+= norm_1(B-zeros(2,2));
        B			= make_integer_band(2,2,-1,0);
        out			+= norm_1(B-zeros(2,2));
        B			= make_integer_band(2,2,0,1);
        out			+= norm_1(B-zeros(2,2));
        B			= make_integer_band(2,2,-1,1);
        out			+= norm_1(B-zeros(2,2));
        B			= make_integer_band(2,2,-3,3);
        out			+= norm_1(B-zeros(2,2));

        //
        B			= make_real_band(2,2,0,0);
        out			+= norm_1(B-zeros(2,2));
        B			= make_real_band(2,2,-1,0);
        out			+= norm_1(B-zeros(2,2));
        B			= make_real_band(2,2,0,1);
        out			+= norm_1(B-zeros(2,2));
        B			= make_real_band(2,2,-1,1);
        out			+= norm_1(B-zeros(2,2));
        B			= make_real_band(2,2,-3,3);
        out			+= norm_1(B-zeros(2,2));

        //
        B			= make_float_band(2,2,0,0);
        out			+= norm_1(B-zeros(2,2));
        B			= make_float_band(2,2,-1,0);
        out			+= norm_1(B-zeros(2,2));
        B			= make_float_band(2,2,0,1);
        out			+= norm_1(B-zeros(2,2));
        B			= make_float_band(2,2,-1,1);
        out			+= norm_1(B-zeros(2,2));
        B			= make_float_band(2,2,-3,3);
        out			+= norm_1(B-zeros(2,2));

        //
        B			= make_complex_band(2,2,0,0);
        out			+= norm_1(B-zeros(2,2));
        B			= make_complex_band(2,2,-1,0);
        out			+= norm_1(B-zeros(2,2));
        B			= make_complex_band(2,2,0,1);
        out			+= norm_1(B-zeros(2,2));
        B			= make_complex_band(2,2,-1,1);
        out			+= norm_1(B-zeros(2,2));
        B			= make_complex_band(2,2,-3,3);
        out			+= norm_1(B-zeros(2,2));

        //
        B			= make_float_complex_band(2,2,0,0);
        out			+= norm_1(B-zeros(2,2));
        B			= make_float_complex_band(2,2,-1,0);
        out			+= norm_1(B-zeros(2,2));
        B			= make_float_complex_band(2,2,0,1);
        out			+= norm_1(B-zeros(2,2));
        B			= make_float_complex_band(2,2,-1,1);
        out			+= norm_1(B-zeros(2,2));
        B			= make_float_complex_band(2,2,-3,3);
        out			+= norm_1(B-zeros(2,2));

        //
        Matrix A;
        A			= irand(2,2);
        out			+= norm_1(make_integer_dense(2,2,A.get_array<Integer>())-A);
        A			= irand(0,2);
        out			+= norm_1(make_integer_dense(0,2,A.get_array<Integer>())-A);
        A			= irand(2,0);
        out			+= norm_1(make_integer_dense(2,0,A.get_array<Integer>())-A);
        A			= irand(0,0);
        out			+= norm_1(make_integer_dense(0,0,A.get_array<Integer>())-A);
        A			= irand(1,1);
        out			+= norm_1(make_integer_dense(1,1,A.get_array<Integer>())-A);

        //
        A			= randn(2,2);
        out			+= norm_1(make_real_dense(2,2,A.get_array<Real>())-A);
        A			= randn(0,2);
        out			+= norm_1(make_real_dense(0,2,A.get_array<Real>())-A);
        A			= randn(2,0);
        out			+= norm_1(make_real_dense(2,0,A.get_array<Real>())-A);
        A			= randn(0,0);
        out			+= norm_1(make_real_dense(0,0,A.get_array<Real>())-A);
        A			= randn(1,1);
        out			+= norm_1(make_real_dense(1,1,A.get_array<Real>())-A);

        //
        A			= frandn(2,2);
        out			+= norm_1(make_float_dense(2,2,A.get_array<Float>())-A);
        A			= frandn(0,2);
        out			+= norm_1(make_float_dense(0,2,A.get_array<Float>())-A);
        A			= frandn(2,0);
        out			+= norm_1(make_float_dense(2,0,A.get_array<Float>())-A);
        A			= frandn(0,0);
        out			+= norm_1(make_float_dense(0,0,A.get_array<Float>())-A);
        A			= frandn(1,1);
        out			+= norm_1(make_float_dense(1,1,A.get_array<Float>())-A);

        //
        A			= crandn(2,2);
        out			+= norm_1(make_complex_dense(2,2,A.get_array<Complex>())-A);
        A			= crandn(0,2);
        out			+= norm_1(make_complex_dense(0,2,A.get_array<Complex>())-A);
        A			= crandn(2,0);
        out			+= norm_1(make_complex_dense(2,0,A.get_array<Complex>())-A);
        A			= crandn(0,0);
        out			+= norm_1(make_complex_dense(0,0,A.get_array<Complex>())-A);
        A			= crandn(1,1);
        out			+= norm_1(make_complex_dense(1,1,A.get_array<Complex>())-A);

        A			= randn(2,2);
        out			+= norm_1(make_complex_dense(2,2,A.get_array<Real>(),A.get_array<Real>())-A- I*A);
        A			= randn(0,2);
        out			+= norm_1(make_complex_dense(0,2,A.get_array<Real>(),A.get_array<Real>())-A- I*A);
        A			= randn(2,0);
        out			+= norm_1(make_complex_dense(2,0,A.get_array<Real>(),A.get_array<Real>())-A- I*A);
        A			= randn(0,0);
        out			+= norm_1(make_complex_dense(0,0,A.get_array<Real>(),A.get_array<Real>())-A- I*A);
        A			= randn(1,1);
        out			+= norm_1(make_complex_dense(1,1,A.get_array<Real>(),A.get_array<Real>())-A- I*A);

        //
        A			= fcrandn(2,2);
        out			+= norm_1(make_float_complex_dense(2,2,A.get_array<Float_complex>())-A);
        A			= fcrandn(0,2);
        out			+= norm_1(make_float_complex_dense(0,2,A.get_array<Float_complex>())-A);
        A			= fcrandn(2,0);
        out			+= norm_1(make_float_complex_dense(2,0,A.get_array<Float_complex>())-A);
        A			= fcrandn(0,0);
        out			+= norm_1(make_float_complex_dense(0,0,A.get_array<Float_complex>())-A);
        A			= fcrandn(1,1);
        out			+= norm_1(make_float_complex_dense(1,1,A.get_array<Float_complex>())-A);

        A			= frandn(2,2);
        out			+= norm_1(make_float_complex_dense(2,2,A.get_array<Float>(),A.get_array<Float>())-A- If*A);
        A			= frandn(0,2);
        out			+= norm_1(make_float_complex_dense(0,2,A.get_array<Float>(),A.get_array<Float>())-A- If*A);
        A			= frandn(2,0);
        out			+= norm_1(make_float_complex_dense(2,0,A.get_array<Float>(),A.get_array<Float>())-A- If*A);
        A			= frandn(0,0);
        out			+= norm_1(make_float_complex_dense(0,0,A.get_array<Float>(),A.get_array<Float>())-A- If*A);
        A			= frandn(1,1);
        out			+= norm_1(make_float_complex_dense(1,1,A.get_array<Float>(),A.get_array<Float>())-A- If*A);

        //
        Matrix Ir,Ic,Ix,C,D;
        A			= isprand(5,5,.1);
        tie(Ir,Ic,Ix) = find3(A);
        Ir			= full(Ir);
        Ic			= full(Ic);
        Ix			= full(Ix);
        B			= make_sparse_matrix(Ir,Ic,Ix,5,5);
        C			= make_integer_sparse(Ir.get_array<Integer>(),Ic.get_array<Integer>(),
                            Ix.get_array<Integer>(),5,5,Ir.rows());
        out			+= norm_1(A-B);
        out			+= norm_1(A-C);

        //
        A			= sprandn(5,5,.1);
        tie(Ir,Ic,Ix) = find3(A);
        Ir			= full(Ir);
        Ic			= full(Ic);
        Ix			= full(Ix);
        B			= make_sparse_matrix(Ir,Ic,Ix,5,5);
        C			= make_real_sparse(Ir.get_array<Integer>(),Ic.get_array<Integer>(),
                            Ix.get_array<Real>(),5,5,Ir.rows());
        out			+= norm_1(A-B);
        out			+= norm_1(A-C);

        //
        A			= fsprandn(5,5,.1);
        tie(Ir,Ic,Ix) = find3(A);
        Ir			= full(Ir);
        Ic			= full(Ic);
        Ix			= full(Ix);
        B			= make_sparse_matrix(Ir,Ic,Ix,5,5);
        C			= make_float_sparse(Ir.get_array<Integer>(),Ic.get_array<Integer>(),
                            Ix.get_array<Float>(),5,5,Ir.rows());
        out			+= norm_1(A-B);
        out			+= norm_1(A-C);

        //
        A			= csprandn(5,5,.1);
        tie(Ir,Ic,Ix) = find3(A);
        Ir			= full(Ir);
        Ic			= full(Ic);
        Ix			= full(Ix);
        Matrix Ixr  = full(real(Ix));
        Matrix Ixi  = full(imag(Ix));
        B			= make_sparse_matrix(Ir,Ic,Ix,5,5);
        C			= make_complex_sparse(Ir.get_array<Integer>(),Ic.get_array<Integer>(),
                            Ix.get_array<Complex>(),5,5,Ir.rows());
        D			= make_complex_sparse(Ir.get_array<Integer>(),Ic.get_array<Integer>(),
                            Ixr.get_array<Real>(),Ixi.get_array<Real>(),
                            5,5,Ir.rows());
        out			+= norm_1(A-B);
        out			+= norm_1(A-C);
        out			+= norm_1(A-D);

        //
        A			= fcsprandn(5,5,.1);
        tie(Ir,Ic,Ix) = find3(A);
        Ir			= full(Ir);
        Ic			= full(Ic);
        Ix			= full(Ix);
        Ixr         = full(real(Ix));
        Ixi         = full(imag(Ix));
        B			= make_sparse_matrix(Ir,Ic,Ix,5,5);
        C			= make_float_complex_sparse(Ir.get_array<Integer>(),Ic.get_array<Integer>(),
                            Ix.get_array<Float_complex>(),5,5,Ir.rows());
        D			= make_float_complex_sparse(Ir.get_array<Integer>(),Ic.get_array<Integer>(),
                            Ixr.get_array<Float>(),Ixi.get_array<Float>(),
                            5,5,Ir.rows());

        out			+= norm_1(A-B);
        out			+= norm_1(A-C);
        out			+= norm_1(A-D);

        //
        A			= isprand(5,5,.1);
        tie(Ir,Ic,Ix) = find3(A);
        Ir			= full(Ir);
        Ic			= full(Ic);
        Ix			= full(Ix);
        B			= make_sparse_matrix(Ir,Ic,Ix,5,5,5);
        C			= make_integer_sparse(Ir.get_array<Integer>(),Ic.get_array<Integer>(),
                            Ix.get_array<Integer>(),5,5,Ir.rows(),5);
        out			+= norm_1(A-B);
        out			+= norm_1(A-C);

        //
        A			= sprandn(5,5,.1);
        tie(Ir,Ic,Ix) = find3(A);
        Ir			= full(Ir);
        Ic			= full(Ic);
        Ix			= full(Ix);
        B			= make_sparse_matrix(Ir,Ic,Ix,5,5,5);
        C			= make_real_sparse(Ir.get_array<Integer>(),Ic.get_array<Integer>(),
                            Ix.get_array<Real>(),5,5,Ir.rows(),5);
        out			+= norm_1(A-B);
        out			+= norm_1(A-C);

        //
        A			= fsprandn(5,5,.1);
        tie(Ir,Ic,Ix) = find3(A);
        Ir			= full(Ir);
        Ic			= full(Ic);
        Ix			= full(Ix);
        B			= make_sparse_matrix(Ir,Ic,Ix,5,5,5);
        C			= make_float_sparse(Ir.get_array<Integer>(),Ic.get_array<Integer>(),
                            Ix.get_array<Float>(),5,5,Ir.rows(),5);
        out			+= norm_1(A-B);
        out			+= norm_1(A-C);

        //
        A			= csprandn(5,5,.1);
        tie(Ir,Ic,Ix) = find3(A);
        Ir			= full(Ir);
        Ic			= full(Ic);
        Ix			= full(Ix);
        Ixr         = full(real(Ix));
        Ixi         = full(imag(Ix));

        B			= make_sparse_matrix(Ir,Ic,Ix,5,5,5);
        C			= make_complex_sparse(Ir.get_array<Integer>(),Ic.get_array<Integer>(),
                            Ix.get_array<Complex>(),5,5,Ir.rows(),5);
        D			= make_complex_sparse(Ir.get_array<Integer>(),Ic.get_array<Integer>(),
                            Ixr.get_array<Real>(),Ixi.get_array<Real>(),
                            5,5,Ir.rows(),5);
        out			+= norm_1(A-B);
        out			+= norm_1(A-C);
        out			+= norm_1(A-D);

        //
        A			= fcsprandn(5,5,.1);
        tie(Ir,Ic,Ix) = find3(A);
        Ir			= full(Ir);
        Ic			= full(Ic);
        Ix			= full(Ix);
        Ixr         = full(real(Ix));
        Ixi         = full(imag(Ix));

        B			= make_sparse_matrix(Ir,Ic,Ix,5,5,5);
        C			= make_float_complex_sparse(Ir.get_array<Integer>(),Ic.get_array<Integer>(),
                            Ix.get_array<Float_complex>(),5,5,Ir.rows(),5);
        D			= make_float_complex_sparse(Ir.get_array<Integer>(),Ic.get_array<Integer>(),
                            Ixr.get_array<Float>(),Ixi.get_array<Float>(),
                            5,5,Ir.rows(),5);
        out			+= norm_1(A-B);
        out			+= norm_1(A-C);
        out			+= norm_1(A-D);

        for (Integer vc = (Integer)value_code::v_integer; vc <= (Integer)value_code::v_complex; ++vc)
        {
            out         += test_dense_ld(0,0, vc);
            out         += test_dense_ld(0,1, vc);
            out         += test_dense_ld(1,0, vc);
            out         += test_dense_ld(1,1, vc);
            out         += test_dense_ld(3,2, vc);
            out         += test_dense_ld(2,3, vc);
        };

        for (Integer vc = (Integer)value_code::v_integer; vc <= (Integer)value_code::v_complex; ++vc)
        {
            out         += test_other_dense_cons(0,0, vc);
            out         += test_other_dense_cons(0,1, vc);
            out         += test_other_dense_cons(1,0, vc);
            out         += test_other_dense_cons(1,1, vc);
            out         += test_other_dense_cons(3,2, vc);
            out         += test_other_dense_cons(2,3, vc);
            out         += test_other_dense_cons(3,3, vc);

            out         += test_other_sparse_cons(0,0, 0, vc);
            out         += test_other_sparse_cons(0,1, 0, vc);
            out         += test_other_sparse_cons(1,0, 0, vc);
            out         += test_other_sparse_cons(1,1, 0, vc);
            out         += test_other_sparse_cons(3,2, 0, vc);
            out         += test_other_sparse_cons(2,3, 0, vc);
            out         += test_other_sparse_cons(3,3, 0, vc);

            out         += test_other_sparse_cons(0,0, 3, vc);
            out         += test_other_sparse_cons(0,1, 3, vc);
            out         += test_other_sparse_cons(1,0, 3, vc);
            out         += test_other_sparse_cons(1,1, 3, vc);
            out         += test_other_sparse_cons(3,2, 3, vc);
            out         += test_other_sparse_cons(2,3, 3, vc);
            out         += test_other_sparse_cons(3,3, 3, vc);

            out         += test_other_band_cons(0, 0, 0, 0, vc);
            out         += test_other_band_cons(0, 1, -1, 1, vc);            
            out         += test_other_band_cons(2, 0, -1, 0, vc);
            out         += test_other_band_cons(3, 3, 0, 1, vc);
            out         += test_other_band_cons(4, 6, -2, 2, vc);
            out         += test_other_band_cons(7, 9, 2, 3, vc);

            out         += test_other_band_cons(7, 9, 0, 0, vc);
            out         += test_other_band_cons(7, 9, 0, 2, vc);
            out         += test_other_band_cons(7, 9, -1, 0, vc);
            out         += test_other_band_cons(7, 9, -1, -1, vc);
            out         += test_other_band_cons(7, 9, -2, 1, vc);
            out         += test_other_band_cons(7, 9, -2, 2, vc);
            out         += test_other_band_cons(7, 9, 1, 2, vc);
            out         += test_other_band_cons(7, 9, 2, 1, vc);
            out         += test_other_band_cons(7, 9, -2, -3, vc);
        };
    }
    catch(std::exception& ex)
    {
        matcl::out_stream  <<  std::string() +ex.what() + "\n";
        out = 1;
    };

    return out;
};

Real test_function_construct::test_dense_ld(Integer r, Integer c, int vc)
{
    Matrix A    = rand(r + 1, c, (value_code)vc);
    A           = A(colon(1,r), colon());
    Integer ld  = r + 1;

    (void)ld;

    Matrix B;
    Matrix C;

    Matrix As;

    switch(vc)
    {
        case value_code::v_integer:
        {
            using Mat       = raw::Matrix<Integer, struct_dense>;
            const Mat& m    = convert(A, Mat::matrix_code).get_impl<Mat>();
            As              = Matrix(m,false);
            B               = make_integer_dense(r, c, m.ptr(), m.ld());
            C               = make_dense_foreign(r, c, (Integer*)m.ptr(), m.ld());
            break;
        }
        case value_code::v_real:
        {
            using Mat       = raw::Matrix<Real, struct_dense>;
            const Mat& m    = convert(A, Mat::matrix_code).get_impl<Mat>();
            As              = Matrix(m,false);

            B               = make_real_dense(r, c, m.ptr(), m.ld());
            C               = make_dense_foreign(r, c, (Real*)m.ptr(), m.ld());
            break;
        }
        case value_code::v_float:
        {
            using Mat       = raw::Matrix<Float, struct_dense>;
            const Mat& m    = convert(A, Mat::matrix_code).get_impl<Mat>();
            As              = Matrix(m,false);

            B               = make_float_dense(r, c, m.ptr(), m.ld());
            C               = make_dense_foreign(r, c, (Float*)m.ptr(), m.ld());
            break;
        }
        case value_code::v_complex:
        {
            using Mat       = raw::Matrix<Complex, struct_dense>;
            const Mat& m    = convert(A, Mat::matrix_code).get_impl<Mat>();
            As              = Matrix(m,false);

            B               = make_complex_dense(r, c, m.ptr(), m.ld());
            C               = make_dense_foreign(r, c, (Complex*)m.ptr(), m.ld());
            break;
        }
        case value_code::v_float_complex:
        {
            using Mat       = raw::Matrix<Float_complex, struct_dense>;
            const Mat& m    = convert(A, Mat::matrix_code).get_impl<Mat>();
            As              = Matrix(m,false);

            B               = make_float_complex_dense(r, c, m.ptr(), m.ld());
            C               = make_dense_foreign(r, c, (Float_complex*)m.ptr(), m.ld());
            break;
        }
    }

    Real out1   = norm_1(A - B);
    Real out2   = norm_1(A - C);
    Real out    = out1 + out2;

    return out;
};

Real test_function_construct::test_other_dense_cons(Integer r, Integer c, int vc)
{
    Real out    = 0;
    Matrix A    = make_dense_matrix(r, c, (value_code)vc);

    if (A.get_value_code() != (value_code)vc)
        out     += 1;
    if (A.get_struct_code() != struct_code::struct_dense)
        out     += 1;

    out         += norm_1(A - zeros(r,c));

    A           = make_dense_noinit(r, c, (value_code)vc);
    if (A.get_value_code() != (value_code)vc)
        out     += 1;
    if (A.get_struct_code() != struct_code::struct_dense)
        out     += 1;

    switch(vc)
    {
        case value_code::v_integer:
        {
            Integer* ptr;
            A   = make_integer_dense_noinit(r,c, ptr);
            break;
        }
        case value_code::v_real:
        {
            Real* ptr;
            A   = make_real_dense_noinit(r,c, ptr);
            break;
        }
        case value_code::v_float:
        {
            Float* ptr;
            A   = make_float_dense_noinit(r,c, ptr);
            break;
        }
        case value_code::v_complex:
        {
            Complex* ptr;
            A   = make_complex_dense_noinit(r,c, ptr);
            break;
        }
        case value_code::v_float_complex:
        {
            Float_complex* ptr;
            A   = make_float_complex_dense_noinit(r,c, ptr);
            break;
        }
    }

    if (A.rows() != r || A.cols() != c)
        out += 1.0;
    if (A.get_value_code() != (value_code)vc)
        out += 1.0;
    if (A.get_struct_code() != struct_code::struct_dense)
        out += 1.0;

    return out;
};

Real test_function_construct::test_other_sparse_cons(Integer r, Integer c, Integer nz, int vc)
{
    Real out    = 0;
    Matrix A    = make_sparse_matrix(r, c, nz, (value_code)vc);

    if (A.get_value_code() != (value_code)vc)
        out     += 1;
    if (A.get_struct_code() != struct_code::struct_sparse)
        out     += 1;

    out         += norm_1(A - zeros(r,c));

    A           = make_sparse_noinit(r, c, nz, (value_code)vc);
    if (A.get_value_code() != (value_code)vc)
        out     += 1;
    if (A.get_struct_code() != struct_code::struct_sparse)
        out     += 1;

    switch(vc)
    {
        case value_code::v_integer:
        {
            A   = make_integer_sparse_noinit(r,c, nz);
            break;
        }
        case value_code::v_real:
        {
            A   = make_real_sparse_noinit(r,c, nz);
            break;
        }
        case value_code::v_float:
        {
            A   = make_float_sparse_noinit(r,c, nz);
            break;
        }
        case value_code::v_complex:
        {
            A   = make_complex_sparse_noinit(r,c, nz);
            break;
        }
        case value_code::v_float_complex:
        {
            A   = make_float_complex_sparse_noinit(r,c, nz);
            break;
        }
    }

    if (A.rows() != r || A.cols() != c)
        out += 1.0;
    if (A.get_value_code() != (value_code)vc)
        out += 1.0;
    if (A.get_struct_code() != struct_code::struct_sparse)
        out += 1.0;

    return out;
};

Real test_function_construct::test_other_band_cons(Integer r, Integer c, Integer fd, Integer ld, int vc)
{
    Integer ldiag   = std::max(0, -fd);
    Integer udiag   = std::max(0, ld);
    bool empty      = (r == 0 || c == 0);

    Real out    = 0;
    Matrix A    = make_band_matrix(r, c, fd, ld, (value_code)vc);

    if (A.get_value_code() != (value_code)vc)
        out     += 1;
    if (A.get_struct_code() != struct_code::struct_banded)
        out     += 1;
    if (A.structural_ldiags(false) != ldiag && empty == false)
        out += 1.0;
    if (A.structural_udiags(false) != udiag && empty == false)
        out += 1.0;

    out         += norm_1(A - zeros(r,c));

    A           = make_band_noinit(r, c, fd, ld, (value_code)vc);
    if (A.get_value_code() != (value_code)vc)
        out     += 1;
    if (A.get_struct_code() != struct_code::struct_banded)
        out     += 1;
    if (A.structural_ldiags(false) != ldiag && empty == false)
        out += 1.0;
    if (A.structural_udiags(false) != udiag && empty == false)
        out += 1.0;

    switch(vc)
    {
        case value_code::v_integer:
        {
            A   = make_integer_band_noinit(r,c, fd, ld);
            break;
        }
        case value_code::v_real:
        {
            A   = make_real_band_noinit(r,c, fd, ld);
            break;
        }
        case value_code::v_float:
        {
            A   = make_float_band_noinit(r,c, fd, ld);
            break;
        }
        case value_code::v_complex:
        {
            A   = make_complex_band_noinit(r,c, fd, ld);
            break;
        }
        case value_code::v_float_complex:
        {
            A   = make_float_complex_band_noinit(r,c, fd, ld);
            break;
        }
    }

    if (A.rows() != r || A.cols() != c)
        out += 1.0;
    if (A.get_value_code() != (value_code)vc)
        out += 1.0;
    if (A.get_struct_code() != struct_code::struct_banded)
        out += 1.0;
    if (A.structural_ldiags(false) != ldiag && empty == false)
        out += 1.0;
    if (A.structural_udiags(false) != udiag && empty == false)
        out += 1.0;

    return out;
};

Real test_function_seq::make()
{
    Real out = 0.;
    try
    {
        //range
        {
            Matrix A    = range(1.0, -1.0);
            if (A.rows() != 1 || A.cols() != 0)
                out     += 1;
            if (A.get_value_code() != value_code::v_real)
                out     += 1;

            A           = range(1.0, 5.0);
            if (A.rows() != 1 || A.cols() != 5)
                out     += 1;
            if (A.get_value_code() != value_code::v_real)
                out     += 1;

            Matrix B    = mat_row().add(1.0).add(2.0).add(3.0).add(4.0).add(5.0);
            out         += norm_1(A - B);

            Matrix C    = range(1.0, 1.0, 5.0);
            out         += norm_1(A - C);

            A           = range(1.0, 2.0, 5.0);
            if (A.rows() != 1 || A.cols() != 3)
                out     += 1;
            if (A.get_value_code() != value_code::v_real)
                out     += 1;

            B           = mat_row().add(1.0).add(3.0).add(5.0);
            out         += norm_1(A - B);

            A           = range(1.0, 2.0, 6.0);
            if (A.rows() != 1 || A.cols() != 3)
                out     += 1;
            if (A.get_value_code() != value_code::v_real)
                out     += 1;

            B           = mat_row().add(1.0).add(3.0).add(5.0);
            out         += norm_1(A - B);

            A           = range(1.0, -2.0, 2.0);
            if (A.rows() != 1 || A.cols() != 0)
                out     += 1;
            if (A.get_value_code() != value_code::v_real)
                out     += 1;

            A           = range(1.0, -2.0, -3.0);
            if (A.rows() != 1 || A.cols() != 3)
                out     += 1;
            if (A.get_value_code() != value_code::v_real)
                out     += 1;

            B           = mat_row().add(1.0).add(-1.0).add(-3.0);
            out         += norm_1(A - B);

            A           = range(1.0, -2.0, -4.0);
            if (A.rows() != 1 || A.cols() != 3)
                out     += 1;
            if (A.get_value_code() != value_code::v_real)
                out     += 1;

            out         += norm_1(A - B);
        };

        {
            Matrix A    = frange(1.0f, -1.0f);
            if (A.rows() != 1 || A.cols() != 0)
                out     += 1;
            if (A.get_value_code() != value_code::v_float)
                out     += 1;

            A           = frange(1.0f, 5.0f);
            if (A.rows() != 1 || A.cols() != 5)
                out     += 1;
            if (A.get_value_code() != value_code::v_float)
                out     += 1;

            Matrix B    = mat_row().add(1.0f).add(2.0f).add(3.0f).add(4.0f).add(5.0f);
            out         += norm_1(A - B);

            Matrix C    = frange(1.0f, 1.0f, 5.0f);
            out         += norm_1(A - C);

            A           = frange(1.0f, 2.0f, 5.0f);
            if (A.rows() != 1 || A.cols() != 3)
                out     += 1;
            if (A.get_value_code() != value_code::v_float)
                out     += 1;

            B           = mat_row().add(1.0f).add(3.0f).add(5.0f);
            out         += norm_1(A - B);

            A           = frange(1.0f, 2.0f, 6.0f);
            if (A.rows() != 1 || A.cols() != 3)
                out     += 1;
            if (A.get_value_code() != value_code::v_float)
                out     += 1;

            B           = mat_row().add(1.0f).add(3.0f).add(5.0f);
            out         += norm_1(A - B);

            A           = frange(1.0f, -2.0f, 2.0f);
            if (A.rows() != 1 || A.cols() != 0)
                out     += 1;
            if (A.get_value_code() != value_code::v_float)
                out     += 1;

            A           = frange(1.0f, -2.0f, -3.0f);
            if (A.rows() != 1 || A.cols() != 3)
                out     += 1;
            if (A.get_value_code() != value_code::v_float)
                out     += 1;

            B           = mat_row().add(1.0f).add(-1.0f).add(-3.0f);
            out         += norm_1(A - B);

            A           = frange(1.0f, -2.0f, -4.0f);
            if (A.rows() != 1 || A.cols() != 3)
                out     += 1;
            if (A.get_value_code() != value_code::v_float)
                out     += 1;

            out         += norm_1(A - B);
        };

        {
            Matrix A    = irange(1, -1);
            if (A.rows() != 1 || A.cols() != 0)
                out     += 1;
            if (A.get_value_code() != value_code::v_integer)
                out     += 1;

            A           = irange(1, 5);
            if (A.rows() != 1 || A.cols() != 5)
                out     += 1;
            if (A.get_value_code() != value_code::v_integer)
                out     += 1;

            Matrix B    = mat_row().add(1).add(2).add(3).add(4).add(5);
            out         += norm_1(A - B);

            Matrix C    = irange(1, 1, 5);
            out         += norm_1(A - C);

            A           = irange(1, 2, 5);
            if (A.rows() != 1 || A.cols() != 3)
                out     += 1;
            if (A.get_value_code() != value_code::v_integer)
                out     += 1;

            B           = mat_row().add(1).add(3).add(5);
            out         += norm_1(A - B);

            A           = irange(1, 2, 6);
            if (A.rows() != 1 || A.cols() != 3)
                out     += 1;
            if (A.get_value_code() != value_code::v_integer)
                out     += 1;

            B           = mat_row().add(1).add(3).add(5);
            out         += norm_1(A - B);

            A           = irange(1, -2, 2);
            if (A.rows() != 1 || A.cols() != 0)
                out     += 1;
            if (A.get_value_code() != value_code::v_integer)
                out     += 1;

            A           = irange(1, -2, -3);
            if (A.rows() != 1 || A.cols() != 3)
                out     += 1;
            if (A.get_value_code() != value_code::v_integer)
                out     += 1;

            B           = mat_row().add(1).add(-1).add(-3);
            out         += norm_1(A - B);

            A           = irange(1, -2, -4);
            if (A.rows() != 1 || A.cols() != 3)
                out     += 1;
            if (A.get_value_code() != value_code::v_integer)
                out     += 1;

            out         += norm_1(A - B);
        };

        //linspace
        {
            Matrix A;

            for (int i = 0; i <= 1; ++i)
            {
                value_code vc   = (i == 0)? value_code::v_real : value_code::v_float;

                out         += test_linspace(0.0, 0.0, -1, vc, A);
                out         += test_linspace(0.0, 0.0, 0, vc, A);
                out         += test_linspace(0.0, 0.0, 1, vc, A);
                out         += test_linspace(0.0, 0.0, 3, vc, A);

                out         += test_linspace(5.0, 0.0, -1, vc, A);
                out         += test_linspace(5.0, 0.0, 0, vc, A);
                out         += test_linspace(5.0, 0.0, 1, vc, A);
                out         += test_linspace(5.0, 0.0, 3, vc, A);

                out         += test_linspace(0.0, 5.0, -1, vc, A);
                out         += test_linspace(0.0, 5.0, 0, vc, A);
                out         += test_linspace(0.0, 5.0, 1, vc, A);
                out         += test_linspace(0.0, 5.0, 3, vc, A);

                out         += test_linspace(0.0, constants::inf(), -1, vc, A);
                out         += test_linspace(0.0, constants::inf(), 0, vc, A);
                out         += test_linspace(0.0, constants::inf(), 1, vc, A);
                out         += test_linspace(0.0, constants::inf(), 5, vc, A);

                out         += test_linspace(constants::inf(), 0.0, -1, vc, A);
                out         += test_linspace(constants::inf(), 0.0, 0, vc, A);
                out         += test_linspace(constants::inf(), 0.0, 1, vc, A);
                out         += test_linspace(constants::inf(), 0.0, 5, vc, A);

                out         += test_linspace(0.0, constants::nan(), -1, vc, A);
                out         += test_linspace(0.0, constants::nan(), 0, vc, A);
                out         += test_linspace(0.0, constants::nan(), 1, vc, A);
                out         += test_linspace(0.0, constants::nan(), 5, vc, A);

                out         += test_linspace(constants::nan(), 0.0, -1, vc, A);
                out         += test_linspace(constants::nan(), 0.0, 0, vc, A);
                out         += test_linspace(constants::nan(), 0.0, 1, vc, A);
                out         += test_linspace(constants::nan(), 0.0, 5, vc, A);
            };
        };

        //logspace
        {
            Matrix A;

            for (int i = 0; i <= 1; ++i)
            {
                value_code vc   = (i == 0)? value_code::v_real : value_code::v_float;

                out         += test_logspace(0.0, 0.0, -1, vc, A);
                out         += test_logspace(0.0, 0.0, 0, vc, A);
                out         += test_logspace(0.0, 0.0, 1, vc, A);
                out         += test_logspace(0.0, 0.0, 3, vc, A);

                out         += test_logspace(5.0, 0.0, -1, vc, A);
                out         += test_logspace(5.0, 0.0, 0, vc, A);
                out         += test_logspace(5.0, 0.0, 1, vc, A);
                out         += test_logspace(5.0, 0.0, 3, vc, A);

                out         += test_logspace(0.0, 5.0, -1, vc, A);
                out         += test_logspace(0.0, 5.0, 0, vc, A);
                out         += test_logspace(0.0, 5.0, 1, vc, A);
                out         += test_logspace(0.0, 5.0, 3, vc, A);

                out         += test_logspace(0.0, constants::inf(), -1, vc, A);
                out         += test_logspace(0.0, constants::inf(), 0, vc, A);
                out         += test_logspace(0.0, constants::inf(), 1, vc, A);
                out         += test_logspace(0.0, constants::inf(), 5, vc, A);

                out         += test_logspace(constants::inf(), 0.0, -1, vc, A);
                out         += test_logspace(constants::inf(), 0.0, 0, vc, A);
                out         += test_logspace(constants::inf(), 0.0, 1, vc, A);
                out         += test_logspace(constants::inf(), 0.0, 5, vc, A);

                out         += test_logspace(0.0, constants::nan(), -1, vc, A);
                out         += test_logspace(0.0, constants::nan(), 0, vc, A);
                out         += test_logspace(0.0, constants::nan(), 1, vc, A);
                out         += test_logspace(0.0, constants::nan(), 5, vc, A);

                out         += test_logspace(constants::nan(), 0.0, -1, vc, A);
                out         += test_logspace(constants::nan(), 0.0, 0, vc, A);
                out         += test_logspace(constants::nan(), 0.0, 1, vc, A);
                out         += test_logspace(constants::nan(), 0.0, 5, vc, A);
            };
        };
    }
    catch(std::exception& ex)
    {
        matcl::out_stream  <<  std::string() +ex.what() + "\n";
        out = 1;
    };

    return out;
};

Real test_function_seq::test_linspace(Real s, Real e, Integer n, value_code vc, Matrix& A)
{
    Real out    = 0.0;

    if (vc == value_code::v_real)
        A       = linspace(s, e, n);
    else
        A       = flinspace(Float(s), Float(e), n);

    Integer n2  = std::max(n, 0);

    if (A.get_value_code() != vc)
        out     += 1;
    if (A.rows() != 1 || A.cols() != n2)
        out     += 1;

    if (n < 1)
        return out;

    Real sA     = A(1).get_scalar<Real>();
    Real eA     = A(end).get_scalar<Real>();

    if (n == 1)
    {
        out         += norm_1(e - eA);
    }
    else
    {
        out         += norm_1(s - sA);
        out         += norm_1(e - eA);
    };

    return out;
};

Real test_function_seq::test_logspace(Real s, Real e, Integer n, value_code vc, Matrix& A)
{
    Real out    = 0.0;

    if (vc == value_code::v_real)
        A       = logspace(s, e, n);
    else
        A       = flogspace(Float(s), Float(e), n);

    Integer n2  = std::max(n, 0);

    if (A.get_value_code() != vc)
        out     += 1;
    if (A.rows() != 1 || A.cols() != n2)
        out     += 1;

    if (n < 1)
        return out;

    Real sA     = A(1).get_scalar<Real>();
    Real eA     = A(end).get_scalar<Real>();
    Real s10    = std::pow(10., s);
    Real e10    = std::pow(10., e);

    if (n == 1)
    {
        out         += norm_1(e10 - eA);
    }
    else
    {
        out         += norm_1(s10 - sA);
        out         += norm_1(e10 - eA);
    };

    return out;
};

Real test_function_zeros::make()
{
    Real out = 0.;
    try
    {
        Matrix A;

        for (int i = 0; i <= 4; ++i)
        {
            value_code vc = value_code::v_integer;
            switch(i)
            {
                case 0:
                    vc  = value_code::v_integer;
                    break;
                case 1:
                    vc  = value_code::v_real;
                    break;
                case 2:
                    vc  = value_code::v_float;
                    break;
                case 3:
                    vc  = value_code::v_complex;
                    break;
                case 4:
                    vc  = value_code::v_float_complex;
                    break;
            };

            out         += test_zeros(0, 0, vc, A);
            out         += test_zeros(0, 1, vc, A);
            out         += test_zeros(2, 0, vc, A);
            out         += test_zeros(3, 3, vc, A);
            out         += test_zeros(4, 6, vc, A);
            out         += test_zeros(7, 9, vc, A);

            out         += test_spzeros(0, 0, 0, vc, A);
            out         += test_spzeros(0, 1, 0, vc, A);
            out         += test_spzeros(2, 0, 0, vc, A);
            out         += test_spzeros(3, 3, 0, vc, A);
            out         += test_spzeros(4, 6, 0, vc, A);
            out         += test_spzeros(7, 9, 0, vc, A);

            out         += test_spzeros(0, 0, 10, vc, A);
            out         += test_spzeros(0, 1, 10, vc, A);
            out         += test_spzeros(2, 0, 10, vc, A);
            out         += test_spzeros(3, 3, 10, vc, A);
            out         += test_spzeros(4, 6, 10, vc, A);
            out         += test_spzeros(7, 9, 10, vc, A);

            out         += test_bzeros(0, 0, 0, 0, vc, A);
            out         += test_bzeros(0, 1, -1, 1, vc, A);            
            out         += test_bzeros(2, 0, -1, 0, vc, A);
            out         += test_bzeros(3, 3, 0, 1, vc, A);
            out         += test_bzeros(4, 6, -2, 2, vc, A);
            out         += test_bzeros(7, 9, 2, 3, vc, A);

            out         += test_bzeros(7, 9, 0, 0, vc, A);
            out         += test_bzeros(7, 9, 0, 2, vc, A);
            out         += test_bzeros(7, 9, -1, 0, vc, A);
            out         += test_bzeros(7, 9, -1, -1, vc, A);
            out         += test_bzeros(7, 9, -2, 1, vc, A);
            out         += test_bzeros(7, 9, -2, 2, vc, A);
            out         += test_bzeros(7, 9, 1, 2, vc, A);
            out         += test_bzeros(7, 9, 2, 1, vc, A);
            out         += test_bzeros(7, 9, -2, -3, vc, A);
        };
    }
    catch(std::exception& ex)
    {
        matcl::out_stream  <<  std::string() +ex.what() + "\n";
        out = 1;
    };

    return out;
};

Real test_function_zeros::test_zeros(Integer r, Integer c, value_code vc, Matrix& A)
{
    switch(vc)
    {
        case value_code::v_integer:
            A   = izeros(r,c);
            break;
        case value_code::v_real:
            A   = zeros(r,c);
            break;
        case value_code::v_float:
            A   = fzeros(r,c);
            break;
        case value_code::v_complex:
            A   = czeros(r,c);
            break;
        case value_code::v_float_complex:
            A   = fczeros(r,c);
            break;
    }

    Real out = 0.0;

    if (A.rows() != r || A.cols() != c)
        out += 1.0;
    if (A.get_value_code() != vc)
        out += 1.0;
    if (A.get_struct_code() != struct_code::struct_dense)
        out += 1.0;

    check_struct(A);

    out += norm_1(A);

    Matrix B    = zeros(r,c,vc);

    if (B.rows() != r || B.cols() != c)
        out += 1.0;
    if (B.get_value_code() != vc)
        out += 1.0;
    if (B.get_struct_code() != struct_code::struct_dense)
        out += 1.0;

    check_struct(B);

    out += norm_1(B);
    return out;
};

Real test_function_zeros::test_spzeros(Integer r, Integer c, Integer nz, value_code vc, Matrix& A)
{
    switch(vc)
    {
        case value_code::v_integer:
            A   = ispzeros(r,c,nz);
            break;
        case value_code::v_real:
            A   = spzeros(r,c,nz);
            break;
        case value_code::v_float:
            A   = fspzeros(r,c,nz);
            break;
        case value_code::v_complex:
            A   = cspzeros(r,c,nz);
            break;
        case value_code::v_float_complex:
            A   = fcspzeros(r,c,nz);
            break;
    }

    Real out = 0.0;

    if (A.rows() != r || A.cols() != c)
        out += 1.0;
    if (A.get_value_code() != vc)
        out += 1.0;
    if (A.get_struct_code() != struct_code::struct_sparse)
        out += 1.0;

    check_struct(A);

    out += norm_1(A);

    Matrix B    = spzeros(r,c,nz,vc);

    if (B.rows() != r || B.cols() != c)
        out += 1.0;
    if (B.get_value_code() != vc)
        out += 1.0;
    if (B.get_struct_code() != struct_code::struct_sparse)
        out += 1.0;

    check_struct(B);

    out += norm_1(B);
    return out;
};

Real test_function_zeros::test_bzeros(Integer r, Integer c, Integer fd, Integer ld, value_code vc, Matrix& A)
{
    switch(vc)
    {
        case value_code::v_integer:
            A   = ibzeros(r,c,fd,ld);
            break;
        case value_code::v_real:
            A   = bzeros(r,c,fd,ld);
            break;
        case value_code::v_float:
            A   = fbzeros(r,c,fd,ld);
            break;
        case value_code::v_complex:
            A   = cbzeros(r,c,fd,ld);
            break;
        case value_code::v_float_complex:
            A   = fcbzeros(r,c,fd,ld);
            break;
    }

    Real out = 0.0;

    if (A.rows() != r || A.cols() != c)
        out += 1.0;
    if (A.get_value_code() != vc)
        out += 1.0;
    if (A.get_struct_code() != struct_code::struct_banded)
        out += 1.0;

    check_struct(A);

    out += norm_1(A);

    Matrix B    = bzeros(r,c,fd,ld,vc);

    if (B.rows() != r || B.cols() != c)
        out += 1.0;
    if (B.get_value_code() != vc)
        out += 1.0;
    if (B.get_struct_code() != struct_code::struct_banded)
        out += 1.0;

    check_struct(B);

    out += norm_1(B);

    if (r == 0 || c == 0)
        return out;

    Integer ldiag   = std::max(0, -fd);
    Integer udiag   = std::max(0, ld);

    if (A.structural_ldiags(false) != ldiag)
        out += 1.0;
    if (A.structural_udiags(false) != udiag)
        out += 1.0;
    if (B.structural_ldiags(false) != ldiag)
        out += 1.0;
    if (B.structural_udiags(false) != udiag)
        out += 1.0;

    return out;
};

Real test_function_ones::make()
{
    Real out = 0.;
    try
    {
        Matrix A;

        for (int i = 0; i <= 4; ++i)
        {
            value_code vc = value_code::v_integer;

            switch(i)
            {
                case 0:
                    vc  = value_code::v_integer;
                    break;
                case 1:
                    vc  = value_code::v_real;
                    break;
                case 2:
                    vc  = value_code::v_float;
                    break;
                case 3:
                    vc  = value_code::v_complex;
                    break;
                case 4:
                    vc  = value_code::v_float_complex;
                    break;
            };

            out         += test_ones(0, 0, vc, A);
            out         += test_ones(0, 1, vc, A);
            out         += test_ones(2, 0, vc, A);
            out         += test_ones(3, 3, vc, A);
            out         += test_ones(4, 6, vc, A);
            out         += test_ones(7, 9, vc, A);

            out         += test_spones(0, 0, vc, A);
            out         += test_spones(0, 1, vc, A);
            out         += test_spones(2, 0, vc, A);
            out         += test_spones(3, 3, vc, A);
            out         += test_spones(4, 6, vc, A);
            out         += test_spones(7, 9, vc, A);

            out         += test_bones(0, 0, vc, A);
            out         += test_bones(0, 1, vc, A);
            out         += test_bones(2, 0, vc, A);
            out         += test_bones(3, 3, vc, A);
            out         += test_bones(4, 6, vc, A);
            out         += test_bones(7, 9, vc, A);
        };
    }
    catch(std::exception& ex)
    {
        matcl::out_stream  <<  std::string() +ex.what() + "\n";
        out = 1;
    };

    return out;
};

Real test_function_ones::test_ones(Integer r, Integer c, value_code vc, Matrix& A)
{
    switch(vc)
    {
        case value_code::v_integer:
            A   = iones(r,c);
            break;
        case value_code::v_real:
            A   = ones(r,c);
            break;
        case value_code::v_float:
            A   = fones(r,c);
            break;
        case value_code::v_complex:
            A   = cones(r,c);
            break;
        case value_code::v_float_complex:
            A   = fcones(r,c);
            break;
    }

    Real out = 0.0;

    if (A.rows() != r || A.cols() != c)
        out += 1.0;
    if (A.get_value_code() != vc)
        out += 1.0;
    if (A.get_struct_code() != struct_code::struct_dense)
        out += 1.0;

    check_struct(A);

    out += norm_1(A - 1);

    Matrix B    = ones(r,c,vc);

    if (B.rows() != r || B.cols() != c)
        out += 1.0;
    if (B.get_value_code() != vc)
        out += 1.0;
    if (B.get_struct_code() != struct_code::struct_dense)
        out += 1.0;

    check_struct(B);

    out += norm_1(B - 1);
    return out;
};

Real test_function_ones::test_spones(Integer r, Integer c, value_code vc, Matrix& A)
{
    switch(vc)
    {
        case value_code::v_integer:
            A   = ispones(r,c);
            break;
        case value_code::v_real:
            A   = spones(r,c);
            break;
        case value_code::v_float:
            A   = fspones(r,c);
            break;
        case value_code::v_complex:
            A   = cspones(r,c);
            break;
        case value_code::v_float_complex:
            A   = fcspones(r,c);
            break;
    }

    Real out = 0.0;

    if (A.rows() != r || A.cols() != c)
        out += 1.0;
    if (A.get_value_code() != vc)
        out += 1.0;
    if (A.get_struct_code() != struct_code::struct_sparse)
        out += 1.0;

    check_struct(A);

    out += norm_1(A - 1);

    Matrix B    = spones(r,c,vc);

    if (B.rows() != r || B.cols() != c)
        out += 1.0;
    if (B.get_value_code() != vc)
        out += 1.0;
    if (B.get_struct_code() != struct_code::struct_sparse)
        out += 1.0;

    check_struct(B);

    out += norm_1(B - 1);
    return out;
};

Real test_function_ones::test_bones(Integer r, Integer c, value_code vc, Matrix& A)
{
    switch(vc)
    {
        case value_code::v_integer:
            A   = ibones(r,c);
            break;
        case value_code::v_real:
            A   = bones(r,c);
            break;
        case value_code::v_float:
            A   = fbones(r,c);
            break;
        case value_code::v_complex:
            A   = cbones(r,c);
            break;
        case value_code::v_float_complex:
            A   = fcbones(r,c);
            break;
    }

    Real out = 0.0;

    if (A.rows() != r || A.cols() != c)
        out += 1.0;
    if (A.get_value_code() != vc)
        out += 1.0;
    if (A.get_struct_code() != struct_code::struct_banded)
        out += 1.0;

    check_struct(A);

    out += norm_1(A - 1);

    Matrix B    = bones(r,c,vc);

    if (B.rows() != r || B.cols() != c)
        out += 1.0;
    if (B.get_value_code() != vc)
        out += 1.0;
    if (B.get_struct_code() != struct_code::struct_banded)
        out += 1.0;

    check_struct(B);

    out += norm_1(B - 1);

    return out;
};

Real test_function_eye::make()
{
    Real out = 0.;
    try
    {
        Matrix A;

        for (int i = 0; i <= 4; ++i)
        {
            value_code vc = value_code::v_integer;
            switch(i)
            {
                case 0:
                    vc  = value_code::v_integer;
                    break;
                case 1:
                    vc  = value_code::v_real;
                    break;
                case 2:
                    vc  = value_code::v_float;
                    break;
                case 3:
                    vc  = value_code::v_complex;
                    break;
                case 4:
                    vc  = value_code::v_float_complex;
                    break;
            };

            out         += test_eye(0, 0, vc, A);
            out         += test_eye(0, 1, vc, A);
            out         += test_eye(2, 0, vc, A);
            out         += test_eye(3, 3, vc, A);
            out         += test_eye(4, 6, vc, A);
            out         += test_eye(7, 9, vc, A);

            out         += test_speye(0, 0, vc, A);
            out         += test_speye(0, 1, vc, A);
            out         += test_speye(2, 0, vc, A);
            out         += test_speye(3, 3, vc, A);
            out         += test_speye(4, 6, vc, A);
            out         += test_speye(7, 9, vc, A);

            out         += test_beye(0, 0, 0, 0, vc, A);
            out         += test_beye(0, 1, -1, 1, vc, A);            
            out         += test_beye(2, 0, -1, 0, vc, A);
            out         += test_beye(3, 3, 0, 1, vc, A);
            out         += test_beye(4, 6, -2, 2, vc, A);
            out         += test_beye(7, 9, 2, 3, vc, A);

            out         += test_beye(9, 9, 0, 0, vc, A);
            out         += test_beye(9, 9, 0, 2, vc, A);
            out         += test_beye(9, 9, -1, 0, vc, A);
            out         += test_beye(9, 9, -1, -1, vc, A);
            out         += test_beye(9, 9, -2, 1, vc, A);
            out         += test_beye(9, 9, -2, 2, vc, A);
            out         += test_beye(9, 9, 1, 2, vc, A);
            out         += test_beye(9, 9, 2, 1, vc, A);
            out         += test_beye(9, 9, -2, -3, vc, A);
        };
    }
    catch(std::exception& ex)
    {
        matcl::out_stream  <<  std::string() +ex.what() + "\n";
        out = 1;
    };

    return out;
};

Real test_function_eye::test_eye(Integer r, Integer c, value_code vc, Matrix& A)
{
    switch(vc)
    {
        case value_code::v_integer:
            A   = ieye(r,c);
            break;
        case value_code::v_real:
            A   = eye(r,c);
            break;
        case value_code::v_float:
            A   = feye(r,c);
            break;
        case value_code::v_complex:
            A   = ceye(r,c);
            break;
        case value_code::v_float_complex:
            A   = fceye(r,c);
            break;
    }

    Real out = 0.0;

    if (A.rows() != r || A.cols() != c)
        out += 1.0;
    if (A.get_value_code() != vc)
        out += 1.0;
    if (A.get_struct_code() != struct_code::struct_dense)
        out += 1.0;

    check_struct(A);

    Matrix D    = A.diag(0);
    out         += norm_1(D - 1);

    Matrix AA   = A;
    AA.diag(0)   = 0;
    out         += norm_1(AA);

    Matrix B    = eye(r,c,vc);

    if (B.rows() != r || B.cols() != c)
        out += 1.0;
    if (B.get_value_code() != vc)
        out += 1.0;
    if (B.get_struct_code() != struct_code::struct_dense)
        out += 1.0;

    check_struct(B);

    D           = B.diag(0);
    out         += norm_1(D - 1);

    B.diag(0)   = 0;
    out         += norm_1(B);

    return out;
};

Real test_function_eye::test_speye(Integer r, Integer c, value_code vc, Matrix& A)
{
    switch(vc)
    {
        case value_code::v_integer:
            A   = ispeye(r,c);
            break;
        case value_code::v_real:
            A   = speye(r,c);
            break;
        case value_code::v_float:
            A   = fspeye(r,c);
            break;
        case value_code::v_complex:
            A   = cspeye(r,c);
            break;
        case value_code::v_float_complex:
            A   = fcspeye(r,c);
            break;
    }

    Real out = 0.0;

    if (A.rows() != r || A.cols() != c)
        out += 1.0;
    if (A.get_value_code() != vc)
        out += 1.0;
    if (A.get_struct_code() != struct_code::struct_sparse)
        out += 1.0;

    check_struct(A);

    out += norm_1(A - eye(r,c));

    Matrix B    = speye(r,c,vc);

    if (B.rows() != r || B.cols() != c)
        out += 1.0;
    if (B.get_value_code() != vc)
        out += 1.0;
    if (B.get_struct_code() != struct_code::struct_sparse)
        out += 1.0;

    check_struct(B);

    out += norm_1(B - eye(r,c));
    return out;
};

Real test_function_eye::test_beye(Integer r, Integer c, Integer fd, Integer ld, value_code vc, Matrix& A)
{
    switch(vc)
    {
        case value_code::v_integer:
            A   = ibeye(r,c,fd,ld);
            break;
        case value_code::v_real:
            A   = beye(r,c,fd,ld);
            break;
        case value_code::v_float:
            A   = fbeye(r,c,fd,ld);
            break;
        case value_code::v_complex:
            A   = cbeye(r,c,fd,ld);
            break;
        case value_code::v_float_complex:
            A   = fcbeye(r,c,fd,ld);
            break;
    }

    Real out = 0.0;

    if (A.rows() != r || A.cols() != c)
        out += 1.0;
    if (A.get_value_code() != vc)
        out += 1.0;
    if (A.get_struct_code() != struct_code::struct_banded)
        out += 1.0;

    check_struct(A);

    out += norm_1(A - eye(r,c));

    Matrix B    = beye(r,c,fd,ld,vc);

    if (B.rows() != r || B.cols() != c)
        out += 1.0;
    if (B.get_value_code() != vc)
        out += 1.0;
    if (B.get_struct_code() != struct_code::struct_banded)
        out += 1.0;

    check_struct(B);

    out += norm_1(B - eye(r,c));

    if (r == 0 || c == 0)
        return out;

    Integer ldiag   = std::max(0, -fd);
    Integer udiag   = std::max(0, ld);

    if (A.structural_ldiags(false) != ldiag)
        out += 1.0;
    if (A.structural_udiags(false) != udiag)
        out += 1.0;
    if (B.structural_ldiags(false) != ldiag)
        out += 1.0;
    if (B.structural_udiags(false) != udiag)
        out += 1.0;

    return out;
};

Real test_function_rand::make()
{
    Real out = 0.;
    try
    {
        Matrix A;

        for (int i = 0; i <= 4; ++i)
        {
            value_code vc = value_code::v_integer;

            switch(i)
            {
                case 0:
                    vc  = value_code::v_integer;
                    break;
                case 1:
                    vc  = value_code::v_real;
                    break;
                case 2:
                    vc  = value_code::v_float;
                    break;
                case 3:
                    vc  = value_code::v_complex;
                    break;
                case 4:
                    vc  = value_code::v_float_complex;
                    break;
            };

            out         += test_dense(1, 1, vc, 0, A);
            out         += test_dense(1, 1, vc, 1, A);

            out         += test_dense(0, 0, vc, 2, A);
            out         += test_dense(0, 1, vc, 2, A);
            out         += test_dense(1, 1, vc, 2, A);
            out         += test_dense(2, 0, vc, 2, A);
            out         += test_dense(3, 3, vc, 2, A);
            out         += test_dense(4, 6, vc, 2, A);
            out         += test_dense(7, 9, vc, 2, A);

            out         += test_dense(0, 0, vc, 3, A);
            out         += test_dense(0, 1, vc, 3, A);
            out         += test_dense(1, 1, vc, 3, A);
            out         += test_dense(2, 0, vc, 3, A);
            out         += test_dense(3, 3, vc, 3, A);
            out         += test_dense(4, 6, vc, 3, A);
            out         += test_dense(7, 9, vc, 3, A);

            out         += test_sparse(0, 0, vc, 2, A);
            out         += test_sparse(0, 1, vc, 2, A);
            out         += test_sparse(1, 1, vc, 2, A);
            out         += test_sparse(2, 0, vc, 2, A);
            out         += test_sparse(3, 3, vc, 2, A);
            out         += test_sparse(4, 6, vc, 2, A);
            out         += test_sparse(7, 9, vc, 2, A);

            out         += test_sparse(0, 0, vc, 3, A);
            out         += test_sparse(0, 1, vc, 3, A);
            out         += test_sparse(1, 1, vc, 3, A);
            out         += test_sparse(2, 0, vc, 3, A);
            out         += test_sparse(3, 3, vc, 3, A);
            out         += test_sparse(4, 6, vc, 3, A);
            out         += test_sparse(7, 9, vc, 3, A);

            for (int type = 2; type <= 3; ++type)
            {
                out         += test_band(0, 0, 0, 0, vc, type, A);
                out         += test_band(0, 1, -1, 1, vc, type, A);            
                out         += test_band(2, 0, -1, 0, vc, type, A);
                out         += test_band(3, 3, 0, 1, vc, type, A);
                out         += test_band(4, 6, -2, 2, vc, type, A);
                out         += test_band(7, 9, 2, 3, vc, type, A);

                out         += test_band(9, 9, 0, 0, vc, type, A);
                out         += test_band(9, 9, 0, 2, vc, type, A);
                out         += test_band(9, 9, -1, 0, vc, type, A);
                out         += test_band(9, 9, -1, -1, vc, type, A);
                out         += test_band(9, 9, -2, 1, vc, type, A);
                out         += test_band(9, 9, -2, 2, vc, type, A);
                out         += test_band(9, 9, 1, 2, vc, type, A);
                out         += test_band(9, 9, 2, 1, vc, type, A);
                out         += test_band(9, 9, -2, -3, vc, type, A);
            };
        };

        out += test_randperm(0, A);
        out += test_randperm(1, A);
        out += test_randperm(10, A);
    }
    catch(std::exception& ex)
    {
        matcl::out_stream  <<  std::string() +ex.what() + "\n";
        out = 1;
    };

    return out;
};

Real test_function_rand::test_dense(Integer r, Integer c, value_code vc, int type, Matrix& A)
{
    Matrix B;

    if (type == 0)
    {
        //scalar rand
        switch(vc)
        {
            case value_code::v_integer:
                A   = irand();
                break;
            case value_code::v_real:
                A   = rand();
                break;
            case value_code::v_float:
                A   = frand();
                break;
            case value_code::v_complex:
                A   = crand();
                break;
            case value_code::v_float_complex:
                A   = fcrand();
                break;
        }

        B = A;
    }
    else if (type == 1)
    {
        //scalar randn
        switch(vc)
        {
            case value_code::v_integer:
                A   = irand();
                break;
            case value_code::v_real:
                A   = randn();
                break;
            case value_code::v_float:
                A   = frandn();
                break;
            case value_code::v_complex:
                A   = crandn();
                break;
            case value_code::v_float_complex:
                A   = fcrandn();
                break;
        }

        B = A;
    }
    else if (type == 2)
    {
        //matrix rand
        switch(vc)
        {
            case value_code::v_integer:
                A   = irand(r,c);
                break;
            case value_code::v_real:
                A   = rand(r,c);
                break;
            case value_code::v_float:
                A   = frand(r,c);
                break;
            case value_code::v_complex:
                A   = crand(r,c);
                break;
            case value_code::v_float_complex:
                A   = fcrand(r,c);
                break;
        }

        B = rand(r,c,vc);
    }
    else if (type == 3)
    {
        //matrix randn
        switch(vc)
        {
            case value_code::v_integer:
                A   = irand(r,c);
                break;
            case value_code::v_real:
                A   = randn(r,c);
                break;
            case value_code::v_float:
                A   = frandn(r,c);
                break;
            case value_code::v_complex:
                A   = crandn(r,c);
                break;
            case value_code::v_float_complex:
                A   = fcrandn(r,c);
                break;
        }

        if (vc != value_code::v_integer)
            B = randn(r,c,vc);
        else
            B = A;
    };

    Real out = 0.0;

    if (A.rows() != r || A.cols() != c)
        out += 1.0;
    if (A.get_value_code() != vc)
        out += 1.0;
    if (A.get_struct_code() != struct_code::struct_dense && type > 1)
        out += 1.0;

    check_struct(A);

    //values should be small
    Real n_out      = norm_1(A);
    if (n_out > 1e3 && vc != value_code::v_integer)
        out         += 1.0;

    if (B.rows() != r || B.cols() != c)
        out += 1.0;
    if (B.get_value_code() != vc)
        out += 1.0;
    if (B.get_struct_code() != struct_code::struct_dense && type > 1)
        out += 1.0;

    check_struct(B);

    n_out           = norm_1(B);
    if (n_out > 1e3 && vc != value_code::v_integer)
        out         += 1.0;

    return out;
};

Real test_function_rand::test_sparse(Integer r, Integer c, value_code vc, int type, Matrix& A)
{
    Matrix B;

    Real d  = 0.3;

    if (type == 2)
    {
        //matrix rand
        switch(vc)
        {
            case value_code::v_integer:
                A   = isprand(r,c,d);
                break;
            case value_code::v_real:
                A   = sprand(r,c,d);
                break;
            case value_code::v_float:
                A   = fsprand(r,c,d);
                break;
            case value_code::v_complex:
                A   = csprand(r,c,d);
                break;
            case value_code::v_float_complex:
                A   = fcsprand(r,c,d);
                break;
        }

        B = sprand(r,c,d,vc);
    }
    else if (type == 3)
    {
        //matrix randn
        switch(vc)
        {
            case value_code::v_integer:
                A   = isprand(r,c,d);
                break;
            case value_code::v_real:
                A   = sprandn(r,c,d);
                break;
            case value_code::v_float:
                A   = fsprandn(r,c,d);
                break;
            case value_code::v_complex:
                A   = csprandn(r,c,d);
                break;
            case value_code::v_float_complex:
                A   = fcsprandn(r,c,d);
                break;
        }

        if (vc != value_code::v_integer)
            B = sprandn(r,c,d,vc);
        else
            B = A;
    };

    Real out = 0.0;

    if (A.rows() != r || A.cols() != c)
        out += 1.0;
    if (A.get_value_code() != vc)
        out += 1.0;
    if (A.get_struct_code() != struct_code::struct_sparse)
        out += 1.0;

    check_struct(A);

    //values should be small
    Real n_out      = norm_1(A);
    if (n_out > 1e3 && vc != value_code::v_integer)
        out         += 1.0;

    if (B.rows() != r || B.cols() != c)
        out += 1.0;
    if (B.get_value_code() != vc)
        out += 1.0;
    if (B.get_struct_code() != struct_code::struct_sparse)
        out += 1.0;

    check_struct(B);

    n_out           = norm_1(B);
    if (n_out > 1e3 && vc != value_code::v_integer)
        out         += 1.0;

    return out;
};

Real test_function_rand::test_band(Integer r, Integer c, Integer fd, Integer ld, value_code vc, int type, Matrix& A)
{
    Matrix B;

    if (type == 2)
    {
        //matrix rand
        switch(vc)
        {
            case value_code::v_integer:
                A   = irand_band(r,c,fd,ld);
                break;
            case value_code::v_real:
                A   = rand_band(r,c,fd,ld);
                break;
            case value_code::v_float:
                A   = frand_band(r,c,fd,ld);
                break;
            case value_code::v_complex:
                A   = crand_band(r,c,fd,ld);
                break;
            case value_code::v_float_complex:
                A   = fcrand_band(r,c,fd,ld);
                break;
        }

        B = rand_band(r,c,fd,ld,vc);
    }
    else if (type == 3)
    {
        //matrix randn
        switch(vc)
        {
            case value_code::v_integer:
                A   = irand_band(r,c,fd,ld);
                break;
            case value_code::v_real:
                A   = randn_band(r,c,fd,ld);
                break;
            case value_code::v_float:
                A   = frandn_band(r,c,fd,ld);
                break;
            case value_code::v_complex:
                A   = crandn_band(r,c,fd,ld);
                break;
            case value_code::v_float_complex:
                A   = fcrandn_band(r,c,fd,ld);
                break;
        }

        if (vc != value_code::v_integer)
            B = randn_band(r,c,fd,ld,vc);
        else
            B = A;
    };

    Real out = 0.0;

    if (A.rows() != r || A.cols() != c)
        out += 1.0;
    if (A.get_value_code() != vc)
        out += 1.0;
    if (A.get_struct_code() != struct_code::struct_banded)
        out += 1.0;

    check_struct(A);

    //values should be small
    Real n_out      = norm_1(A);
    if (n_out > 1e3 && vc != value_code::v_integer)
        out         += 1.0;

    if (B.rows() != r || B.cols() != c)
        out += 1.0;
    if (B.get_value_code() != vc)
        out += 1.0;
    if (B.get_struct_code() != struct_code::struct_banded)
        out += 1.0;

    check_struct(B);

    n_out           = norm_1(B);
    if (n_out > 1e3 && vc != value_code::v_integer)
        out         += 1.0;

    if (r == 0 || c == 0)
        return out;

    Integer ldiag   = std::max(0, -fd);
    Integer udiag   = std::max(0, ld);

    if (A.structural_ldiags(false) != ldiag)
        out += 1.0;
    if (A.structural_udiags(false) != udiag)
        out += 1.0;
    if (B.structural_ldiags(false) != ldiag)
        out += 1.0;
    if (B.structural_udiags(false) != udiag)
        out += 1.0;

    return out;
};

Real test_function_rand::test_randperm(Integer n, Matrix& A)
{
    permvec p   = randperm(n);
    A           = p.to_matrix();

    if (p.length() != n)
        return 1.0;

    if (any_vec(A < 0) || any_vec(A > n))
        return 1.0;

    if (A.rows() != n || A.cols() != 1)
        return 1.0;

    return 0;
};

};};
