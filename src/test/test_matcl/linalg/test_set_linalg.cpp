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

#include "test/test_matcl/framework/data_struct/test_selector.h"
#include "test_set_linalg.h"
#include "matcl-matrep/details/refcount.h"

#include <iostream>

namespace matcl { namespace test
{

void test_all_linalg(const rand_matrix_ptr& rand)
{
    long N = 0;
    
    if (test::test_selector().is_selected(2,"single_thread"))   
    {
        matcl::out_stream << std::string(50, '-') + "\n";
        matcl::out_stream << std::string() + "\t\t\tsingle_thread" + "\n";
        matcl::out_stream << std::string(50, '-') + "\n";

        if (test::test_selector().is_selected_no_scons_parse(4,"unary"))   
            test::test_linalg_st(rand);     N = matcl::details::no_existing_objects();

        if (test::test_selector().is_selected_no_scons_parse(4,"binary"))   
            test::test_linalg_bin_st(rand); N = matcl::details::no_existing_objects();               
    }
    
    if (test::test_selector().is_selected(2,"multi_thread"))   
    {
        matcl::out_stream << std::string(50, '-') + "\n";
        matcl::out_stream << std::string() + "\t\t\tmulti_thread" + "\n";
        matcl::out_stream << std::string(50, '-') + "\n";

        if (test::test_selector().is_selected_no_scons_parse(4,"unary"))   
            test::test_linalg_st(rand);     N = matcl::details::no_existing_objects();

        if (test::test_selector().is_selected_no_scons_parse(4,"binary"))   
            test::test_linalg_bin_st(rand); N = matcl::details::no_existing_objects();
    }
};

};};
