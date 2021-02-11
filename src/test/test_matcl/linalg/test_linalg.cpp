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


#include "test/test_matcl/framework/test_main_body.h"
#include "test_set_linalg.h"
#include "test/test_matcl/framework/data_struct/test_selector.h"

using namespace matcl;

int main(int argc, const char* argv[])
{
    init_genrand(2);
    
    //test::test_selector().add_selection(1,"full_matrices1");
    //test::test_selector().add_selection(1,"full_matrices2");
    //test::test_selector().add_selection(1,"sparse_matrices_true");    
    //test::test_selector().add_selection(1,"sparse_matrices_false");  
    //test::test_selector().add_selection(1,"str_matrices");      
    //test::test_selector().add_selection(1,"dense_matrices_rand");          
    
    //test::test_selector().add_selection(2,"multi_thread");          

    //test::test_selector().add_selection(3,"test_eigs()");

    //test::test_selector().add_selection(4,"unary");
    //test::test_selector().add_selection(4,"binary");

    test::test_setups(argc, argv);
    
    int ret = test::test_main_body(&test::test_all_linalg);
    
    if (ret != 0) 
        return ret;
}
