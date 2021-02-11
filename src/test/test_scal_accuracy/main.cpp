/*
 *  This file is a part of Matrix Computation Library (MATCL)
 *
 *  Copyright (c) Pawe³ Kowal 2017 - 2021
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

#include "test_scal_accuracy.h"
#include "matcl-mp/matcl_mp.h"
#include "matcl-core/IO/logger.h"

#include <iostream>

int main(int argc, const char* argv[])
{
    (void)argc;
    (void)argv;

    using log_ptr   = std::shared_ptr<std::ofstream>;

    try
    {         
        {
            std::string log_file_name   = std::string("log_test_scal_accuracy.txt");
            log_ptr log = log_ptr(new std::ofstream(log_file_name));
            matcl::set_logger(log);
        };

        matcl::test::test_scal_accuracy(matcl::out_stream);
        matcl::out_stream << "finished" << "\n";
    }
    catch(std::exception& ex)
    {
        matcl::out_stream << ex.what() << "\n";
        return 1;
    };

    return 0;
}
