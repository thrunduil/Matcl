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

#include "test_dynamic.h"
#include "test_gmp.h"
#include "test_gmp_prec.h"
#include "test_gmp_object.h"
#include "matcl-mp/matcl_mp.h"
#include "matcl-core/options/matcl_options.h"

#include "matcl-scalar/IO/scalar_io.h"
#include "matcl-scalar/IO/formatted_disp.h"
#include "matcl-core/IO/logger.h"

#include <iostream>
#include <iomanip>

namespace mdy = matcl::dynamic;
using namespace matcl;

void example_object();

/*
static void break_func()
{
    disp("break");
};
*/

#include <iostream>
#include <sstream>

int main(int argc, const char* argv[])
{
    (void)argc;
    (void)argv;

    using log_ptr   = std::shared_ptr<std::ofstream>;

    try
    {         
        {
            std::string log_file_name   = std::string("log_test_dynamic.txt");
            log_ptr log = log_ptr(new std::ofstream(log_file_name));
            set_logger(log);
        };

        //details::leak_detector::break_at_codes({37204, 51566, 52566, 25195, 36023, 40397},
        //                                       &break_func);

        example_object();

        matcl::test::test_gmp_object();
        matcl::test::test_gmp();
        matcl::test::test_gmp_bin();        

        matcl::test::test_dynamic test;
        test.make();  

        matcl::test::test_gmp_prec(out_stream);

        disp("finished");
    }
    catch(std::exception& ex)
    {
        disp(ex.what());
        return 1;
    };

    return 0;
}
