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
#include "test_performance.h"
#include "matcl-mp/matcl_mp.h"
#include "matcl-core/options/matcl_options.h"

#include "matcl-scalar/IO/scalar_io.h"
#include "matcl-scalar/IO/formatted_disp.h"

#include <iostream>
#include <iomanip>

namespace mdy = matcl::dynamic;
using namespace matcl;

int main(int argc, const char* argv[])
{
    try
    { 
        matcl::test::test_performance();

        /*
        matcl::test::test_gmp();        
        matcl::test::test_gmp_bin();        
        matcl::test::test_gmp_object();                     

        matcl::test::test_dynamic test;
        test.make();  

        matcl::test::test_gmp_prec(std::cout);
        matcl::test::test_performance();
        */

        matcl::free_caches();

      #if MATCL_DEBUG_MEMORY
        matcl::details::leak_detector::report_leaks(std::cout);
      #endif

        std::cout << "finished" << "\n";
    }
    catch(std::exception& ex)
    {
        std::cout << ex.what() << "\n";
        return 1;
    };

    return 0;
}
