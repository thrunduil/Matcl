/*
 *  This file is a part of Matrix Computation Library (MATCL)
 *
 *  Copyright (c) Pawe� Kowal 2018 - 2021
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
#include "test_set.h"
#include "test/test_matcl/framework/data_struct/test_selector.h"
#include "test/test_matcl/framework/matrix_set/matrix_set.h"
#include "test/test_matcl/framework/matrix_set/test_options.h"

#include "matcl-linalg/matcl_linalg.h"
#include "matcl-blas-lapack/blas/blas.h"
#include "matcl-blas-lapack/blas_loader/blas_loader.h"

#include <iostream>

using namespace matcl;

int main(int argc, const char* argv[])
{
    init_genrand(4);

    matcl::test::test_options::set_seed(3);
    matcl::test::test_options::set_binmat_prob(0.1);
    matcl::test::test_options::set_num_colons(5);
    matcl::test::test_options::set_num_scalar_groups(1);
    matcl::test::test_options::set_num_matrix_groups(1);

    //matcl::test::test_options::set_num_colons(20);
    
    test::test_setups(argc, argv);

    #ifdef WIN64
        std::string plugin_name  = "matcl-openblas-plugin-x64-Release.dll";
        //std::string plugin_name  = "matcl-clapack-plugin-x64-Release.dll";
        //std::string plugin_name  = "matcl-mkl-plugin-x64-Release.dll";
    #else
        //std::string plugin_name  = "matcl-clapack-plugin-Win32-Release.dll";
        std::string plugin_name  = "matcl-openblas-plugin-Win32-Release.dll";
        //std::string plugin_name  = "matcl-mkl-plugin-Win32-Release.dll";
    #endif

    //matcl::lapack::load_blas_plugin(plugin_name);
    //matcl::lapack::initialize_plugin();
    //raw_blas_lapack::set_num_threads_blas_kernel(1);

    //test::test_selector().add_selection(3,"test_herprod()");
    
    //test::test_selector().add_selection(4,"utils");
    //test::test_selector().add_selection(4,"matgen");
    //test::test_selector().add_selection(4,"matrix");
    //test::test_selector().add_selection(4,"manip"); 
    //test::test_selector().add_selection(4,"assign");     
    //test::test_selector().add_selection(4,"matfunc");     
    //test::test_selector().add_selection(4,"unary");     
        
    //test::test_selector().add_selection(2,"multi_thread");
    test::test_selector().add_selection(2,"single_thread");
    
    test::test_selector().add_selection(1,"full_matrices1");
    //test::test_selector().add_selection(1,"full_matrices2");
    //test::test_selector().add_selection(1,"sparse_matrices_true");
    //test::test_selector().add_selection(1,"sparse_matrices_false");
    //test::test_selector().add_selection(1,"str_matrices");
    //test::test_selector().add_selection(1,"dense_matrices_rand");
    //test::test_selector().add_selection(1,"sparse_matrices_rand");   

    int ret = test::test_main_body(&test::test_all);
    
    if (ret != 0) 
        return ret;
}
