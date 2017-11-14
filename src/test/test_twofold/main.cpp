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

#include "test_twofold.h"
#include "test_twofold_simd.h"
#include "test_poly.h"
#include "matcl-core/IO/logger.h"
#include "matcl-core/float/twofold.h"

#include <iostream>
#include <fstream>

using namespace matcl;

int main(int argc, const char* argv[])
{
    (void)argc;
    (void)argv;

    using log_ptr   = std::shared_ptr<std::ofstream>;

    try
    {         
        {
            std::string log_file_name   = std::string("log_test_twofold.txt");
            log_ptr log = log_ptr(new std::ofstream(log_file_name));
            set_logger(log);
        };

        matcl::test::test_poly_dyn();
        matcl::test::test_poly();

        matcl::test::test_functions_simd();
        matcl::test::test_functions();                
        matcl::test::test_fma();        

        matcl::test::test_error();
        matcl::test::test_error_simd();
        
        matcl::test::test_io();
        matcl::test::test_io_simd();

        std::cout << "\n";
        std::cout << "finished" << "\n";
    }
    catch(std::exception& ex)
    {
        std::cout << ex.what() << "\n";
        return 1;
    };

    return 0;
}
