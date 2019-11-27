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

#include "test_main_body.h"

#include <iostream>
#include "matcl-matrep/matcl_matrep.h"

#pragma warning( push )
#pragma warning(disable:4702)	//  unreachable code
#include <boost/lexical_cast.hpp>
#pragma warning( pop )

#include "matcl-core/IO/archive.h"

#include "data_struct/test_selector.h"
#include "matcl-core/IO/logger.h"
#include "crash_behaviour.h"

#include "matcl-linalg/matcl_linalg.h"
#include <fstream>

#include <array> 

using namespace std;
namespace matcl { namespace test
{

void test_setups(int argc, const char* argv[])
{
    disable_crash_dialog();
    report_crash_to_stderr();

    //plot(eye(9)); // no plotting in tester by default

    bool open_log = true;
    for (int i = 1; i < argc; i++)
        open_log = std::string(argv[i]) != "--no-log" && open_log;

    if (open_log)
    {
        std::shared_ptr<std::ofstream> log(new std::ofstream("matcl-test.log"));
        matcl::set_logger(log);
    }

    for (int i = 1; i < argc; i++)
    {
        string arg = std::string(argv[i]);
        size_t pos = arg.find("=");
        if (pos != string::npos)
        {
            stringstream ss(arg.substr(0, 1));
            int layer;
            ss >> layer ? layer : 0;
            test::test_selector().add_selection(layer, arg.substr(pos + 1));
        }
    }
}

int test_main_body(void (*tester)(const rand_matrix_ptr& ))
{
    //TODO: new tests
    //      - NaN/Inf * Matrix  -> dense NaN matrix
    //      - NAN/Matrix, Matrix/NaN -> dense NaN matrix
    //      - new functions
    //      - test inplace functions
    //      - str_matrices with NaNs
    try
    {
        if (test::test_selector().is_selected(1, "full_matrices1"))   
        {
            matcl::out_stream << std::string(50, '-') + "\n";
            matcl::out_stream << std::string() + "\t\t\tfull_matrices1" + "\n";
            matcl::out_stream << std::string(50, '-') + "\n";

            test::rand_matrix_ptr rand(new test::rand_matrix_1());
            tester(rand);
        }
        //{
        //       test::rand_matrix_ptr rand(new test::rand_matrix_obj());
        //       tester(rand);
        //   }
        if (test::test_selector().is_selected(1, "full_matrices2"))   
        {
            matcl::out_stream << std::string(50, '-') + "\n";
            matcl::out_stream << std::string() + "\t\t\tfull_matrices2" + "\n";
            matcl::out_stream << std::string(50, '-') + "\n";

            test::rand_matrix_ptr rand(new test::rand_matrix_1());
            tester(rand);
        }
        //{
        //       test::rand_matrix_ptr rand(new test::rand_matrix_obj());
        //       tester(rand);
        //   }
     
        if (test::test_selector().is_selected(1, "sparse_matrices_true"))   
        {
            matcl::out_stream << std::string(50, '-') + "\n";
            matcl::out_stream << std::string() + "\t\t\tsparse_matrices_true" + "\n";
            matcl::out_stream << std::string(50, '-') + "\n";

            test::rand_matrix_ptr rand(new test::rand_matrix_1_sp(true));
            tester(rand);
        }
        if (test::test_selector().is_selected(1, "sparse_matrices_false"))   
        {
            matcl::out_stream << std::string(50, '-') + "\n";
            matcl::out_stream << std::string() + "\t\t\tsparse_matrices_false" + "\n";
            matcl::out_stream << std::string(50, '-') + "\n";

            test::rand_matrix_ptr rand(new test::rand_matrix_1_sp(false));
            tester(rand);
        }
        if (test::test_selector().is_selected(1, "str_matrices"))   
        {
            matcl::out_stream << std::string(50, '-') + "\n";
            matcl::out_stream << std::string() + "\t\t\tstr_matrices" + "\n";
            matcl::out_stream << std::string(50, '-') + "\n";

            test::rand_matrix_ptr rand(new test::rand_matrix_1_str());
            tester(rand);
        }
        if (test::test_selector().is_selected(1, "dense_matrices_rand"))   
        {
            matcl::out_stream << std::string(50, '-') + "\n";
            matcl::out_stream << std::string() + "\t\t\tdense_matrices_rand" + "\n";
            matcl::out_stream << std::string(50, '-') + "\n";

            test::rand_matrix_ptr base(new test::rand_matrix_1());
            test::rand_matrix_ptr rand(new test::rand_matrix_dense(base));
            tester(rand);
        }
        if (test::test_selector().is_selected(1, "sparse_matrices_rand"))   
        {
            matcl::out_stream << std::string(50, '-') + "\n";
            matcl::out_stream << std::string() + "\t\t\tsparse_matrices_rand" + "\n";
            matcl::out_stream << std::string(50, '-') + "\n";

            test::rand_matrix_ptr base(new test::rand_matrix_1());
            test::rand_matrix_ptr rand(new test::rand_matrix_sparse(base));
            tester(rand);
        }
    }
    catch(const std::exception& ex)
    {
        matcl::out_stream << ex.what();
        matcl::out_stream << "\nTesting terminated because of errors\n";
        return 1;
    }
    catch(...)
    {
        matcl::out_stream << "\nTesting terminated because of UNKNOWN ERROR\n";
        return 1;
    }

    matcl::out_stream << std::string() + "\nTesting finished \n";
    return 0;

}

}}
