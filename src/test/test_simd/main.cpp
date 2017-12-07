/*
 *  This file is a part of Matrix Computation Library (MATCL)
 *
 *  Copyright (c) Pawe� Kowal 2017
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

#include "test_simd_config.h"
#include "test_simd.h"
#include "test_simd_compl.h"
#include "matcl-core/IO/logger.h"

#include <iostream>
#include <fstream>

using namespace matcl;

int main(int argc, const char* argv[])
{
    (void)argc;
    (void)argv;

    using log_ptr   = std::shared_ptr<std::ofstream>;

    //TODO: is_nan, pow2k, shift_left, shift_right, _arithmetic, if_nan_else, if_zero_else, if_then_else
    // casts, gather, extract_low, extract_high, bitwise_not
    try
    {         
        {
            std::string log_file_name   = std::string("log_test_simd_") + MATCL_TEST_SIMD_TAG + ".txt";
            log_ptr log = log_ptr(new std::ofstream(log_file_name));
            set_logger(log);
        };

        matcl::test::test_performance_real();        
        matcl::test::test_performance_complex();        

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
