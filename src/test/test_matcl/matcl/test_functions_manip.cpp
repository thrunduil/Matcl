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

#pragma once

#include "test_set.h"
#include "test_functions_manip.h"

#include "test/test_matcl/framework/matrix_set/matrix_set_1.h"
#include "test/test_matcl/framework/matrix_set/test_options.h"
#include "matcl-core/IO/logger.h"
#include "test/test_matcl/framework/data_struct/test_selector.h"

#include <boost/thread.hpp>


namespace gd = matcl::details;

namespace matcl { namespace test
{

class test_manip
{
    manip_functions_list&		tf;
    const test::options&		opts;
    Integer                     thread_id;

    public:
        test_manip(manip_functions_list& tf, const test::options& opts, Integer id)
            :tf(tf),opts(opts), thread_id(id)
        {};

        test_manip(const test_manip& tu)
            :tf(tu.tf),opts(tu.opts), thread_id(tu.thread_id)
        {};

        void make()
        {   
            /*
            Integer code = 139;
            Matrix mat = tf.get_matrix(code);
            matcl::disp(mat);

            double tol      = 0.0;

            if (mat.rows() == mat.cols())
            {
                bool is     = is_her(mat, tol);
                bool is2    = is_her(mat + ctrans(mat), tol);

                disp(mat - ctrans(mat));

                (void)is2;
                (void)is;

                //if (is == true && norm_1(mat - ctrans(mat)) > 0.0)
                //    out     += 1;
            }
            */

            tf.make(opts);
        };

        void operator()()
        {
            make();
        };

    private:		
        test_manip& operator=(const test_manip&) = delete;
};

void test_manip_st(const rand_matrix_ptr& rand)
{
    test::options opts;

    try
    {
        opts.show_partial_res   = 0;
        opts.first_matrix_code  = 0;
        opts.show_memleaks      = true;

        test::mat_set_1 ms1(rand);
        manip_functions_list tf(ms1,rand);
        
        test_manip tu(tf,opts,0);
        tu.make();
    }
    catch(const std::exception& ex)
    {
        matcl::out_stream <<ex.what();
    };
};

void test_manip_mt(const rand_matrix_ptr& rand)
{
    test::options opts;
    try
    {
        opts.show_partial_res   = 0;
        opts.first_matrix_code  = 0;
        opts.show_memleaks      = false;

        test::mat_set_2 ms1(rand);
        manip_functions_list tf(ms1,rand);

        boost::thread_group tg;

        for (int i = 0; i < 10; i++)
            tg.create_thread(test_manip(tf,opts,i));

        tg.join_all();
    }
    catch(const std::exception& ex)
    {
        matcl::out_stream <<ex.what();
    };
};

manip_functions_list::manip_functions_list(const matrix_set& ms,rand_matrix_ptr rand)
:m_tests(ms), m_rand(rand)
{};

void manip_functions_list::make(options opts)
{
    m_options = opts;

    SELECT_TEST (3, test_diag());
    SELECT_TEST (3, test_bdiag());
    SELECT_TEST (3, test_spdiag());
    SELECT_TEST (3, test_diags());
    SELECT_TEST (3, test_bdiags());
    SELECT_TEST (3, test_spdiags());

    SELECT_TEST (3, test_find());
    SELECT_TEST (3, test_find2());
    SELECT_TEST (3, test_find3());

    SELECT_TEST (3, test_sort());
    SELECT_TEST (3, test_sort2());
    SELECT_TEST (3, test_sortrows());
    SELECT_TEST (3, test_sortrows2());
    SELECT_TEST (3, test_sortcols());
    SELECT_TEST (3, test_sortcols2());

    SELECT_TEST (3, test_sortrows_dim());
    SELECT_TEST (3, test_sortrows2_dim());
    SELECT_TEST (3, test_sortcols_dim());
    SELECT_TEST (3, test_sortcols2_dim());
    SELECT_TEST (3, test_is_sorted());
    SELECT_TEST (3, test_is_sorted_cols());
    SELECT_TEST (3, test_is_sorted_rows());

    SELECT_TEST (3, test_cat());
    SELECT_TEST (3, test_get_diag());
    SELECT_TEST (3, test_tril());
    SELECT_TEST (3, test_triu());  
    SELECT_TEST (3, test_reshape());
    SELECT_TEST (3, test_rot90());
    SELECT_TEST (3, test_vec());
    SELECT_TEST (3, test_trans());
    SELECT_TEST (3, test_ctrans());
    SELECT_TEST (3, test_trans_t());
    SELECT_TEST (3, test_flipud());
    SELECT_TEST (3, test_fliplr());
    SELECT_TEST (3, test_full());
    SELECT_TEST (3, test_sparse());
    SELECT_TEST (3, test_band());
    SELECT_TEST (3, test_clone());

    SELECT_TEST (3, test_get_lu());
    SELECT_TEST (3, test_repmat());
    SELECT_TEST (3, test_convert());

    SELECT_TEST (3, test_matrix_fwd());    
    SELECT_TEST (3, test_convert_val());        
    SELECT_TEST (3, test_horzcat());        
    SELECT_TEST (3, test_vertcat());        
    SELECT_TEST (3, test_blkdiag());        
    SELECT_TEST (3, test_select_band());
    SELECT_TEST (3, test_checks());    
    SELECT_TEST (3, test_nnz_m());    
    SELECT_TEST (3, test_drop_sparse());
};

Matrix manip_functions_list::get_matrix(int code) const
{
    return m_tests.get_matrix(code);
};

//
void manip_functions_list::test_diag()
{
    Real out = 0.;
    {
        test_function_diag tf(0);
        out += m_tests.make(&tf,m_options);
    };
    {
        test_function_diag tf(1);
        out += m_tests.make(&tf,m_options);
    };
    {
        test_function_diag tf(-1);
        out += m_tests.make(&tf,m_options);
    };

    if (out == 0.)
        matcl::out_stream << std::string() +   "diag: OK" + "\n";
    else
        matcl::out_stream << std::string() +   "diag: FAILED"  + "\n";
};

void manip_functions_list::test_bdiag()
{
    Real out = 0.;
    {
        test_function_bdiag tf(0);
        out += m_tests.make(&tf,m_options);
    };
    {
        test_function_bdiag tf(1);
        out += m_tests.make(&tf,m_options);
    };
    {
        test_function_bdiag tf(-1);
        out += m_tests.make(&tf,m_options);
    };

    if (out == 0.)
        matcl::out_stream << std::string() +   "bdiag: OK" + "\n";
    else
        matcl::out_stream << std::string() +   "bdiag: FAILED"  + "\n";
};

void manip_functions_list::test_spdiag()
{
    Real out = 0.;
    {
        test_function_spdiag tf(0);
        out += m_tests.make(&tf,m_options);
    };
    {
        test_function_spdiag tf(1);
        out += m_tests.make(&tf,m_options);
    };
    {
        test_function_spdiag tf(-1);
        out += m_tests.make(&tf,m_options);
    };

    if (out == 0.)
        matcl::out_stream << std::string() +   "spdiag: OK" + "\n";
    else
        matcl::out_stream << std::string() +   "spdiag: FAILED"  + "\n";
};

void manip_functions_list::test_diags()
{
    Real out = 0.;
    {
        test_function_diags tf;
        out += m_tests.make(&tf,m_options);
    };

    if (out == 0.)
        matcl::out_stream << std::string() +   "diags: OK" + "\n";
    else
        matcl::out_stream << std::string() +   "diags: FAILED"  + "\n";
};

void manip_functions_list::test_bdiags()
{
    Real out = 0.;
    {
        test_function_bdiags tf;
        out += m_tests.make(&tf,m_options);
    };

    if (out == 0.)
        matcl::out_stream << std::string() +   "bdiags: OK" + "\n";
    else
        matcl::out_stream << std::string() +   "bdiags: FAILED"  + "\n";
};

void manip_functions_list::test_spdiags()
{
    Real out = 0.;
    {
        test_function_spdiags tf;
        out += m_tests.make(&tf,m_options);
    };

    if (out == 0.)
        matcl::out_stream << std::string() +   "spdiags: OK" + "\n";
    else
        matcl::out_stream << std::string() +   "spdiags: FAILED"  + "\n";
};

void manip_functions_list::test_find()
{
    Real out = 0.;
    test_function_find tf;
    out += m_tests.make(&tf,m_options);

    if (out == 0.)
        matcl::out_stream << std::string() +   "find: OK" + "\n";
    else
        matcl::out_stream << std::string() +   "find: FAILED"  + "\n";
};

void manip_functions_list::test_find2()
{
    Real out = 0.;
    test_function_find2 tf;
    out += m_tests.make(&tf,m_options);

    if (out == 0.)
        matcl::out_stream << std::string() +   "find2: OK" + "\n";
    else
        matcl::out_stream << std::string() +   "find2: FAILED"  + "\n";
};

void manip_functions_list::test_find3()
{
    Real out = 0.;
    test_function_find3 tf;
    out += m_tests.make(&tf,m_options);

    if (out == 0.)
        matcl::out_stream << std::string() +   "find3: OK" + "\n";
    else
        matcl::out_stream << std::string() +   "find3: FAILED"  + "\n";
};

void manip_functions_list::test_sort()
{
    Real out = 0.;
    for (int j = 1; j <= 2; j++)
    {
        for (int i = 1; i <= 2; i++)
        {
            bool is_asc = (j == 1);
            test_function_sort tf(i,is_asc);
            out += m_tests.make(&tf,m_options);
        };
    };

    if (out == 0.)
        matcl::out_stream << std::string() +   "sort: OK" + "\n";
    else
        matcl::out_stream << std::string() +   "sort: FAILED"  + "\n";
};

void manip_functions_list::test_sort2()
{
    Real out = 0.;
    for (int j = 1; j <= 2; j++)
    {
        for (int i = 1; i <= 2; i++)
        {
            bool is_asc = (j == 1);
            test_function_sort2 tf(i,is_asc);
            out += m_tests.make(&tf,m_options);
        };
    };

    if (out == 0.)
        matcl::out_stream << std::string() +   "sort2: OK" + "\n";
    else
        matcl::out_stream << std::string() +   "sort2: FAILED"  + "\n";
};

void manip_functions_list::test_sortrows()
{
    Real out = 0.;
    test_function_sortrows tf;
    out += m_tests.make(&tf,m_options);

    if (out == 0.)
        matcl::out_stream << std::string() +   "sortrows: OK" + "\n";
    else
        matcl::out_stream << std::string() +   "sortrows: FAILED"  + "\n";    
};

void manip_functions_list::test_sortrows2()
{
    Real out = 0.;
    test_function_sortrows2 tf;
    out += m_tests.make(&tf,m_options);

    if (out == 0.)
        matcl::out_stream << std::string() +   "sortrows2: OK" + "\n";
    else
        matcl::out_stream << std::string() +   "sortrows2: FAILED"  + "\n";
};

void manip_functions_list::test_sortrows_dim()
{
    Real out = 0.;

    {
        Matrix dim = mat_row();
        test_function_sortrows_dim tf(dim);
        out += m_tests.make(&tf,m_options);
    };
    {
        Matrix dim = (mat_row(), 1);
        test_function_sortrows_dim tf(dim);
        out += m_tests.make(&tf,m_options);
    };
    {
        Matrix dim = (mat_row(), -1);
        test_function_sortrows_dim tf(dim);
        out += m_tests.make(&tf,m_options);
    };
    {
        Matrix dim = (mat_row(), 3);
        test_function_sortrows_dim tf(dim);
        out += m_tests.make(&tf,m_options);
    };
    {
        Matrix dim = (mat_row(), -3);
        test_function_sortrows_dim tf(dim);
        out += m_tests.make(&tf,m_options);
    };
    {
        Matrix dim = (mat_row(), 1, 2, 3);
        test_function_sortrows_dim tf(dim);
        out += m_tests.make(&tf,m_options);
    };
    {
        Matrix dim = (mat_row(), 1, -2, 3);
        test_function_sortrows_dim tf(dim);
        out += m_tests.make(&tf,m_options);
    };
    {
        Matrix dim = (mat_row(), -1, -2, -3);
        test_function_sortrows_dim tf(dim);
        out += m_tests.make(&tf,m_options);
    };
    {
        Matrix dim = (mat_row(), 1, 3, 2);
        test_function_sortrows_dim tf(dim);
        out += m_tests.make(&tf,m_options);
    };
    {
        Matrix dim = (mat_row(), 1, 3, -2);
        test_function_sortrows_dim tf(dim);
        out += m_tests.make(&tf,m_options);
    };
    {
        Matrix dim = (mat_row(), -1, -3, -2);
        test_function_sortrows_dim tf(dim);
        out += m_tests.make(&tf,m_options);
    };
    {
        Matrix dim = (mat_row(), 3, 2, 1);
        test_function_sortrows_dim tf(dim);
        out += m_tests.make(&tf,m_options);
    };
    {
        Matrix dim = (mat_row(), -3, 2, 1);
        test_function_sortrows_dim tf(dim);
        out += m_tests.make(&tf,m_options);
    };
    {
        Matrix dim = (mat_row(), -3, -2, -1);
        test_function_sortrows_dim tf(dim);
        out += m_tests.make(&tf,m_options);
    };

    if (out == 0.)
        matcl::out_stream << std::string() +   "sortrows_dim: OK" + "\n";
    else
        matcl::out_stream << std::string() +   "sortrows_dim: FAILED"  + "\n";
};

void manip_functions_list::test_sortrows2_dim()
{
    Real out = 0.;
    {
        Matrix dim = mat_row();
        test_function_sortrows_dim2 tf(dim);
        out += m_tests.make(&tf,m_options);
    };
    {
        Matrix dim = (mat_row(), 1);
        test_function_sortrows_dim2 tf(dim);
        out += m_tests.make(&tf,m_options);
    };
    {
        Matrix dim = (mat_row(), -1);
        test_function_sortrows_dim2 tf(dim);
        out += m_tests.make(&tf,m_options);
    };
    {
        Matrix dim = (mat_row(), 3);
        test_function_sortrows_dim2 tf(dim);
        out += m_tests.make(&tf,m_options);
    };
    {
        Matrix dim = (mat_row(), -3);
        test_function_sortrows_dim2 tf(dim);
        out += m_tests.make(&tf,m_options);
    };
    {
        Matrix dim = (mat_row(), 1, 2, 3);
        test_function_sortrows_dim2 tf(dim);
        out += m_tests.make(&tf,m_options);
    };
    {
        Matrix dim = (mat_row(), 1, -2, 3);
        test_function_sortrows_dim2 tf(dim);
        out += m_tests.make(&tf,m_options);
    };
    {
        Matrix dim = (mat_row(), -1, -2, -3);
        test_function_sortrows_dim2 tf(dim);
        out += m_tests.make(&tf,m_options);
    };
    {
        Matrix dim = (mat_row(), 1, 3, 2);
        test_function_sortrows_dim2 tf(dim);
        out += m_tests.make(&tf,m_options);
    };
    {
        Matrix dim = (mat_row(), 1, 3, -2);
        test_function_sortrows_dim2 tf(dim);
        out += m_tests.make(&tf,m_options);
    };
    {
        Matrix dim = (mat_row(), -1, -3, -2);
        test_function_sortrows_dim2 tf(dim);
        out += m_tests.make(&tf,m_options);
    };
    {
        Matrix dim = (mat_row(), 3, 2, 1);
        test_function_sortrows_dim2 tf(dim);
        out += m_tests.make(&tf,m_options);
    };
    {
        Matrix dim = (mat_row(), -3, 2, 1);
        test_function_sortrows_dim2 tf(dim);
        out += m_tests.make(&tf,m_options);
    };
    {
        Matrix dim = (mat_row(), -3, -2, -1);
        test_function_sortrows_dim2 tf(dim);
        out += m_tests.make(&tf,m_options);
    };

    if (out == 0.)
        matcl::out_stream << std::string() +   "sortrows_dim2: OK" + "\n";
    else
        matcl::out_stream << std::string() +   "sortrows_dim2: FAILED"  + "\n";
};

void manip_functions_list::test_sortcols()
{
    Real out = 0.;

    test_function_sortcols tf;
    out += m_tests.make(&tf,m_options);

    if (out == 0.)
        matcl::out_stream << std::string() +   "sortcols: OK" + "\n";
    else
        matcl::out_stream << std::string() +   "sortcols: FAILED"  + "\n";
};

void manip_functions_list::test_sortcols2()
{
    Real out = 0.;

    test_function_sortcols2 tf;
    out += m_tests.make(&tf,m_options);

    if (out == 0.)
        matcl::out_stream << std::string() +   "sortcols2: OK" + "\n";
    else
        matcl::out_stream << std::string() +   "sortcols2: FAILED"  + "\n";
};

void manip_functions_list::test_sortcols_dim()
{
    Real out = 0.;
    {
        Matrix dim = mat_row();
        test_function_sortcols_dim tf(dim);
        out += m_tests.make(&tf,m_options);
    };
    {
        Matrix dim = (mat_row(), 1);
        test_function_sortcols_dim tf(dim);
        out += m_tests.make(&tf,m_options);
    };
    {
        Matrix dim = (mat_row(), -1);
        test_function_sortcols_dim tf(dim);
        out += m_tests.make(&tf,m_options);
    };
    {
        Matrix dim = (mat_row(), 3);
        test_function_sortcols_dim tf(dim);
        out += m_tests.make(&tf,m_options);
    };
    {
        Matrix dim = (mat_row(), -3);
        test_function_sortcols_dim tf(dim);
        out += m_tests.make(&tf,m_options);
    };
    {
        Matrix dim = (mat_row(), 1, 2, 3);
        test_function_sortcols_dim tf(dim);
        out += m_tests.make(&tf,m_options);
    };
    {
        Matrix dim = (mat_row(), 1, -2, 3);
        test_function_sortcols_dim tf(dim);
        out += m_tests.make(&tf,m_options);
    };
    {
        Matrix dim = (mat_row(), -1, -2, -3);
        test_function_sortcols_dim tf(dim);
        out += m_tests.make(&tf,m_options);
    };
    {
        Matrix dim = (mat_row(), 1, 3, 2);
        test_function_sortcols_dim tf(dim);
        out += m_tests.make(&tf,m_options);
    };
    {
        Matrix dim = (mat_row(), 1, 3, -2);
        test_function_sortcols_dim tf(dim);
        out += m_tests.make(&tf,m_options);
    };
    {
        Matrix dim = (mat_row(), -1, -3, -2);
        test_function_sortcols_dim tf(dim);
        out += m_tests.make(&tf,m_options);
    };
    {
        Matrix dim = (mat_row(), 3, 2, 1);
        test_function_sortcols_dim tf(dim);
        out += m_tests.make(&tf,m_options);
    };
    {
        Matrix dim = (mat_row(), -3, 2, 1);
        test_function_sortcols_dim tf(dim);
        out += m_tests.make(&tf,m_options);
    };
    {
        Matrix dim = (mat_row(), -3, -2, -1);
        test_function_sortcols_dim tf(dim);
        out += m_tests.make(&tf,m_options);
    };

    if (out == 0.)
        matcl::out_stream << std::string() +   "sortcols_dim: OK" + "\n";
    else
        matcl::out_stream << std::string() +   "sortcols_dim: FAILED"  + "\n";
};

void manip_functions_list::test_sortcols2_dim()
{
    Real out = 0.;
    {
        Matrix dim = mat_row();
        test_function_sortcols_dim2 tf(dim);
        out += m_tests.make(&tf,m_options);
    };
    {
        Matrix dim = (mat_row(), 1);
        test_function_sortcols_dim2 tf(dim);
        out += m_tests.make(&tf,m_options);
    };
    {
        Matrix dim = (mat_row(), -1);
        test_function_sortcols_dim2 tf(dim);
        out += m_tests.make(&tf,m_options);
    };
    {
        Matrix dim = (mat_row(), 3);
        test_function_sortcols_dim2 tf(dim);
        out += m_tests.make(&tf,m_options);
    };
    {
        Matrix dim = (mat_row(), -3);
        test_function_sortcols_dim2 tf(dim);
        out += m_tests.make(&tf,m_options);
    };
    {
        Matrix dim = (mat_row(), 1, 2, 3);
        test_function_sortcols_dim2 tf(dim);
        out += m_tests.make(&tf,m_options);
    };
    {
        Matrix dim = (mat_row(), 1, -2, 3);
        test_function_sortcols_dim2 tf(dim);
        out += m_tests.make(&tf,m_options);
    };
    {
        Matrix dim = (mat_row(), -1, -2, -3);
        test_function_sortcols_dim2 tf(dim);
        out += m_tests.make(&tf,m_options);
    };
    {
        Matrix dim = (mat_row(), 1, 3, 2);
        test_function_sortcols_dim2 tf(dim);
        out += m_tests.make(&tf,m_options);
    };
    {
        Matrix dim = (mat_row(), 1, 3, -2);
        test_function_sortcols_dim2 tf(dim);
        out += m_tests.make(&tf,m_options);
    };
    {
        Matrix dim = (mat_row(), -1, -3, -2);
        test_function_sortcols_dim2 tf(dim);
        out += m_tests.make(&tf,m_options);
    };
    {
        Matrix dim = (mat_row(), 3, 2, 1);
        test_function_sortcols_dim2 tf(dim);
        out += m_tests.make(&tf,m_options);
    };
    {
        Matrix dim = (mat_row(), -3, 2, 1);
        test_function_sortcols_dim2 tf(dim);
        out += m_tests.make(&tf,m_options);
    };
    {
        Matrix dim = (mat_row(), -3, -2, -1);
        test_function_sortcols_dim2 tf(dim);
        out += m_tests.make(&tf,m_options);
    };

    if (out == 0.)
        matcl::out_stream << std::string() +   "sortcols_dim2: OK" + "\n";
    else
        matcl::out_stream << std::string() +   "sortcols_dim2: FAILED"  + "\n";
};

void manip_functions_list::test_is_sorted()
{
    Real out = 0.;
    for (int j = 1; j <= 2; j++)
    {
        for (int i = 1; i <= 2; i++)
        {
            bool is_asc = (j == 1);
            test_function_is_sorted tf(i,is_asc);
            out += m_tests.make(&tf,m_options);
        };
    };

    if (out == 0.)
        matcl::out_stream << std::string() +   "is_sorted: OK" + "\n";
    else
        matcl::out_stream << std::string() +   "is_sorted: FAILED"  + "\n";
};

void manip_functions_list::test_is_sorted_rows()
{
    Real out = 0.;
    test_function_is_sorted_rows tf;
    out += m_tests.make(&tf,m_options);

    if (out == 0.)
        matcl::out_stream << std::string() +   "is_sorted_rows: OK" + "\n";
    else
        matcl::out_stream << std::string() +   "is_sorted_rows: FAILED"  + "\n";
};

void manip_functions_list::test_is_sorted_cols()
{
    Real out = 0.;
    test_function_is_sorted_cols tf;
    out += m_tests.make(&tf,m_options);

    if (out == 0.)
        matcl::out_stream << std::string() +   "is_sorted_cols: OK" + "\n";
    else
        matcl::out_stream << std::string() +   "is_sorted_cols: FAILED"  + "\n";
};

void manip_functions_list::test_get_diag()
{
    Real out = 0.;
    {
        test_function_get_diag tf(0);
        out += m_tests.make(&tf,m_options);
    };
    {
        test_function_get_diag tf(1);
        out += m_tests.make(&tf,m_options);
    };
    {
        test_function_get_diag tf(-1);
        out += m_tests.make(&tf,m_options);
    };
    {
        test_function_get_diag tf(3);
        out += m_tests.make(&tf,m_options);
    };
    {
        test_function_get_diag tf(-3);
        out += m_tests.make(&tf,m_options);
    };

    if (out == 0.)
        matcl::out_stream << std::string() + "get_diag: OK" + "\n";
    else
        matcl::out_stream << std::string() + "get_diag: FAILED"  + "\n";
};

void manip_functions_list::test_tril()
{
    Real out = 0.;
    {
        test_function_tril tf(0);
        out += m_tests.make(&tf,m_options);
    };
    {
        test_function_tril tf(1);
        out += m_tests.make(&tf,m_options);
    };
    {
        test_function_tril tf(-1);
        out += m_tests.make(&tf,m_options);
    };
    {
        test_function_tril tf(3);
        out += m_tests.make(&tf,m_options);
    };
    {
        test_function_tril tf(-3);
        out += m_tests.make(&tf,m_options);
    };

    if (out == 0.)
        matcl::out_stream << std::string() + "tril: OK" + "\n";
    else
        matcl::out_stream << std::string() + "tril: FAILED"  + "\n";
};

void manip_functions_list::test_triu()
{
    Real out = 0.;

    {
        test_function_triu tf(0);
        out += m_tests.make(&tf,m_options);
    };
    {
        test_function_triu tf(1);
        out += m_tests.make(&tf,m_options);
    };
    {
        test_function_triu tf(-1);
        out += m_tests.make(&tf,m_options);
    };
    {
        test_function_triu tf(3);
        out += m_tests.make(&tf,m_options);
    };
    {
        test_function_triu tf(-3);
        out += m_tests.make(&tf,m_options);
    };

    if (out == 0.)
        matcl::out_stream << std::string() + "triu: OK" + "\n";
    else
        matcl::out_stream << std::string() + "triu: FAILED"  + "\n";
};

void manip_functions_list::test_rot90()
{
    Real out = 0.;
    {
        test_function_rot90 tf(0);
        out += m_tests.make(&tf,m_options);
    };
    {
        test_function_rot90 tf(1);
        out += m_tests.make(&tf,m_options);
    };
    {
        test_function_rot90 tf(2);
        out += m_tests.make(&tf,m_options);
    };
    {
        test_function_rot90 tf(3);
        out += m_tests.make(&tf,m_options);
    };
    {
        test_function_rot90 tf(4);
        out += m_tests.make(&tf,m_options);
    };

    if (out == 0.)
        matcl::out_stream << std::string() + "rot90: OK" + "\n";
    else
        matcl::out_stream << std::string() + "rot90: FAILED"  + "\n";
};

void manip_functions_list::test_reshape()
{
    Real out = 0.;

    test_function_reshape tf;
    out += m_tests.make(&tf,m_options);

    if (out == 0.)
        matcl::out_stream << std::string() + "reshape: OK" + "\n";
    else
        matcl::out_stream << std::string() + "reshape: FAILED"  + "\n";
};

void manip_functions_list::test_cat()
{
    Real out = 0.;
    Integer seed = test_options::get_seed();

    for (int i = 0; i <= 1; i++)
    {
        {
            test_function_cat tf(0,0,0,i,m_rand,seed + i);
            out += m_tests.make(&tf,m_options);
        };
        {
            test_function_cat tf(1,1,1,i,m_rand,seed + i);
            out += m_tests.make(&tf,m_options);
        };
        {
            test_function_cat tf(1,2,2,i,m_rand,seed + i);
            out += m_tests.make(&tf,m_options);
        };
        {
            test_function_cat tf(1,3,3,i,m_rand,seed + i);
            out += m_tests.make(&tf,m_options);
        };
        {
            test_function_cat tf(2,1,1,i,m_rand,seed + i);
            out += m_tests.make(&tf,m_options);
        };
        {
            test_function_cat tf(2,2,2,i,m_rand,seed + i);
            out += m_tests.make(&tf,m_options);
        };
        {
            test_function_cat tf(2,3,3,i,m_rand,seed + i);
            out += m_tests.make(&tf,m_options);
        };
        {
            test_function_cat tf(2,3,6,i,m_rand,seed + i);
            out += m_tests.make(&tf,m_options);
        };
    };

    if (out == 0.)
        matcl::out_stream << std::string() + "cat: OK" + "\n";
    else
        matcl::out_stream << std::string() + "cat: FAILED"  + "\n";
};

void manip_functions_list::test_trans_t()
{
    Real out = 0.;
    test_function_trans_t tf;
    out += m_tests.make(&tf,m_options);

    if (out == 0.)
        matcl::out_stream << std::string() +   "trans_t: OK" + "\n";
    else
        matcl::out_stream << std::string() +   "trans_t: FAILED"  + "\n";
};


void manip_functions_list::test_vec()
{
    return test_function<function_vec>();
};
void manip_functions_list::test_trans()
{
    return test_function<function_trans>();
};
void manip_functions_list::test_ctrans()
{
    return test_function<function_ctrans>();
};;			
void manip_functions_list::test_flipud()
{
    return test_function<function_flipud>();
};;
void manip_functions_list::test_fliplr()
{
    return test_function<function_fliplr>();
};;			

void manip_functions_list::test_sparse()
{
    Real out = 0.;
    test_function_sparse tf;
    out += m_tests.make(&tf,m_options);

    if (out == 0.)
        matcl::out_stream << std::string() +   "sparse: OK" + "\n";
    else
        matcl::out_stream << std::string() +   "sparse: FAILED"  + "\n";
};		

void manip_functions_list::test_full()
{
    Real out = 0.;
    test_function_full tf;
    out += m_tests.make(&tf,m_options);

    if (out == 0.)
        matcl::out_stream << std::string() +   "full: OK" + "\n";
    else
        matcl::out_stream << std::string() +   "full: FAILED"  + "\n";
};		

void manip_functions_list::test_band()
{
    Real out = 0.;
    test_function_band tf;
    out += m_tests.make(&tf,m_options);

    if (out == 0.)
        matcl::out_stream << std::string() +   "band: OK" + "\n";
    else
        matcl::out_stream << std::string() +   "band: FAILED"  + "\n";
};	

void manip_functions_list::test_clone()
{
    Real out = 0.;
    test_function_clone tf;
    out += m_tests.make(&tf,m_options);

    if (out == 0.)
        matcl::out_stream << std::string() +   "clone: OK" + "\n";
    else
        matcl::out_stream << std::string() +   "clone: FAILED"  + "\n";
}

template<class Func>
void manip_functions_list::test_function()
{
    Real out = 0.;
    Func func;
    out = m_tests.make(func.function(),m_options);

    if (out == 0.)
        matcl::out_stream << std::string() +   func.name() + ": OK" + "\n";
    else
        matcl::out_stream << std::string() +   func.name() + ": FAILED"  + "\n";
};

void manip_functions_list::test_repmat()
{
    Real out = 0.;
    for (int j = 1; j <= 2; j++)
    {
        for (int i = 1; i <= 2; i++)
        {
            test_function_repmat tf(i,j);
            out += m_tests.make(&tf,m_options);
        };
    };

    if (out == 0.)
        matcl::out_stream << std::string() +   "repmat: OK" + "\n";
    else
        matcl::out_stream << std::string() +   "repmat: FAILED"  + "\n";
};

void manip_functions_list::test_horzcat()
{
    Real out = 0.;

    test_function_horzcat tf;
    out += m_tests.make(&tf,m_options);

    if (out == 0.)
        matcl::out_stream << std::string() +   "horzcat: OK" + "\n";
    else
        matcl::out_stream << std::string() +   "horzcat: FAILED"  + "\n";
};

void manip_functions_list::test_vertcat()
{
    Real out = 0.;

    test_function_vertcat tf;
    out += m_tests.make(&tf,m_options);

    if (out == 0.)
        matcl::out_stream << std::string() +   "vertcat: OK" + "\n";
    else
        matcl::out_stream << std::string() +   "vertcat: FAILED"  + "\n";
};

void manip_functions_list::test_blkdiag()
{
    Real out = 0.;

    test_function_blkdiag tf;
    out += m_tests.make(&tf,m_options);

    if (out == 0.)
        matcl::out_stream << std::string() +   "blkdiag: OK" + "\n";
    else
        matcl::out_stream << std::string() +   "blkdiag: FAILED"  + "\n";
};

void manip_functions_list::test_checks()
{
    Real out = 0.;

    test_function_checks tf;
    out += m_tests.make(&tf,m_options);

    if (out == 0.)
        matcl::out_stream << std::string() +   "test_checks: OK" + "\n";
    else
        matcl::out_stream << std::string() +   "test_checks: FAILED"  + "\n";
};

void manip_functions_list::test_select_band()
{
    Real out = 0.;

    {
        test_function_select_band tf(-2,-1);
        out += m_tests.make(&tf,m_options);
    }
    {
        test_function_select_band tf(-1,-1);
        out += m_tests.make(&tf,m_options);
    }
    {
        test_function_select_band tf(-1,0);
        out += m_tests.make(&tf,m_options);
    }
    {
        test_function_select_band tf(-1,1);
        out += m_tests.make(&tf,m_options);
    }
    {
        test_function_select_band tf(0,0);
        out += m_tests.make(&tf,m_options);
    }
    {
        test_function_select_band tf(0,1);
        out += m_tests.make(&tf,m_options);
    }
    {
        test_function_select_band tf(1,1);
        out += m_tests.make(&tf,m_options);
    }
    {
        test_function_select_band tf(1,2);
        out += m_tests.make(&tf,m_options);
    }

    if (out == 0.)
        matcl::out_stream << std::string() +   "select_band: OK" + "\n";
    else
        matcl::out_stream << std::string() +   "select_band: FAILED"  + "\n";
};

void manip_functions_list::test_get_lu()
{
    Real out = 0.;
    test_function_get_lu tf;
    out += m_tests.make(&tf,m_options);

    if (out == 0.)
        matcl::out_stream << std::string() +   "get_lu: OK" + "\n";
    else
        matcl::out_stream << std::string() +   "get_lu: FAILED"  + "\n";
};

void manip_functions_list::test_nnz_m()
{
    Real out = 0.;
    test_function_nnz_m tf;
    out += m_tests.make(&tf,m_options);

    if (out == 0.)
        matcl::out_stream << std::string() +   "nnz_m: OK" + "\n";
    else
        matcl::out_stream << std::string() +   "nnz_m: FAILED"  + "\n";
};

void manip_functions_list::test_drop_sparse()
{
    Real out = 0.;
    test_function_drop_sparse tf;
    out += m_tests.make(&tf,m_options);

    if (out == 0.)
        matcl::out_stream << std::string() +   "drop_sparse: OK" + "\n";
    else
        matcl::out_stream << std::string() +   "drop_sparse: FAILED"  + "\n";
};

void manip_functions_list::test_convert()
{
    Real out = 0.;
    for (Integer i = (Integer)mat_code::integer_scalar; i <= (Integer)mat_code::object_band;i++)
    {
        //ignore object value type
        if (i % 6 == 5)
            continue;

        test_function_convert tf(i);
        out += m_tests.make(&tf,m_options);
    };

    if (out == 0.)
        matcl::out_stream << std::string() +   "convert: OK" + "\n";
    else
        matcl::out_stream << std::string() +   "convert: FAILED"  + "\n";
};

void manip_functions_list::test_convert_val()
{
    Real out = 0.;
    for (Integer i = (Integer)value_code::v_integer; i <= (Integer)value_code::v_complex;i++)
    {
        test_function_convert_val tf(i);
        out += m_tests.make(&tf,m_options);
    };

    if (out == 0.)
        matcl::out_stream << std::string() +   "convert_value: OK" + "\n";
    else
        matcl::out_stream << std::string() +   "convert_value: FAILED"  + "\n";
};

void manip_functions_list::test_matrix_fwd()
{
    Real out = 0.;
    test_function_matrix_fwd tf;
    out += m_tests.make(&tf,m_options);

    if (out == 0.)
        matcl::out_stream << std::string() +   "matrix_fwd: OK" + "\n";
    else
        matcl::out_stream << std::string() +   "matrix_fwd: FAILED"  + "\n";
};

//
Real test_function_diag::eval_mat(const Matrix& mat,bool,int code )
{
    (void) code;
    try
    {
        Real dif		= 0;
        Matrix tmp		= vec(mat);
        // make dense diagonal only if mat not too big
        if (Real(tmp.numel()) * Real(tmp.numel()) < 10000)
        {
            Matrix out	= diag(tmp,d);
                
            if (d > 0 && d > out.cols())
                return dif;

            if (d < 0 && -d > out.rows())
                return dif;

            Matrix out2	= get_diag(out,d);

            check_struct(out);
            check_struct(out2);

            dif			+= norm_1(out2 - tmp);
        }
        return dif;
    }
    catch(const std::exception& ex)
    {
        m_is_error  = true;
        m_error     = ex.what();
        return 1.;
    }
    catch(...)
    {
        m_is_error  = true;
        m_error     = "unknown error";
        return 1.;
    };
};
Real test_function_diag::eval_scalar(const Scalar& ,bool,int code )
{
    (void) code;
    return 0;
};

Real test_function_bdiag::eval_mat(const Matrix& mat,bool,int code )
{
    (void) code;
    try
    {
        Real dif		= 0;
        Matrix tmp		= vec(mat);
        {
            Matrix out	= bdiag(tmp,d);
            // make dense diagonal only if mat not too big
            if (Real(tmp.numel()) * Real(tmp.numel()) < 10000)
            {
                Matrix outd	= diag(tmp,d);
                dif			+= norm_1(out - outd);
            }

            if (d > 0 && d > out.cols())
                return dif;

            if (d < 0 && -d > out.rows())
                return dif;

            Matrix out2	= get_diag(out,d);

            check_struct(out);
            check_struct(out2);

            dif			+= norm_1(out2 - tmp);
        };

        return dif;
    }
    catch(const std::exception& ex)
    {
        m_is_error  = true;
        m_error     = ex.what();
        return 1.;
    }
    catch(...)
    {
        m_is_error  = true;
        m_error     = "unknown error";
        return 1.;
    };
};
Real test_function_bdiag::eval_scalar(const Scalar& ,bool,int code )
{
    (void) code;
    return 0;
};

Real test_function_spdiag::eval_mat(const Matrix& mat,bool,int code )
{
    (void) code;
    try
    {
        Real dif		= 0;
        Matrix tmp		= vec(mat);
        {
            Matrix out	= spdiag(tmp,d);
            // make dense diagonal only if mat not too big
            if (Real(tmp.numel()) * Real(tmp.numel()) < 10000)
            {
                Matrix outd	= diag(tmp,d);
                dif			+= norm_1(out - outd);
            }

            if (d > 0 && d > out.cols())
                return dif;

            if (d < 0 && -d > out.rows())
                return dif;

            Matrix out2	= get_diag(out,d);

            check_struct(out);
            check_struct(out2);

            dif			+= norm_1(out2 - tmp);
        };
        return dif;
    }
    catch(const std::exception& ex)
    {
        m_is_error  = true;
        m_error     = ex.what();
        return 1.;
    }
    catch(...)
    {
        m_is_error  = true;
        m_error     = "unknown error";
        return 1.;
    };
};
Real test_function_spdiag::eval_scalar(const Scalar& ,bool,int code )
{
    (void) code;
    return 0;
};

Real test_function_diags::eval_mat(const Matrix& mat,bool,int code )
{
    if (mat.rows() * mat.cols() > 200)
    {
        //Matrix too large
        return 0.;
    };

    (void) code;
    try
    {
        Real dif		= 0;
        Matrix tmp		= vec(mat);
        Integer m       = mat.rows();
        Integer n       = mat.cols();

        // make dense diagonal only if mat not too big
        if (Real(tmp.numel()) * Real(tmp.numel()) < 10000)
        {
            Matrix out	= diags(tmp,(mat_row(), 0),m*n,m*n);
            check_struct(out);

            dif			+= norm_1(out - diag(tmp,0));

            if (1 < m*n)
            {
                Matrix out2	= diags(tmp,(mat_row(), 1),m*n,m*n);
                check_struct(out2);

                dif			+= norm_1(out2 - diag(tmp,1)(colon(1,m*n),colon(1,m*n)) );
            };
            
            if (1 < m*n)
            {
                Matrix out2	= diags(tmp,(mat_row(), -1),m*n,m*n);
                check_struct(out2);

                dif			+= norm_1(out2 - diag(tmp,-1)(colon(1,m*n),colon(1,m*n)) );
            };
        }

        if (n > 1 && m > 1)
        {
            {
                Matrix out	= diags(mat(colon(),colon(1,2)),(mat_row(), -1, 0),m,m);
                check_struct(out);

                dif			+= norm_1(get_diag(out,-1) - mat(colon(1,m-1),1));
                dif			+= norm_1(get_diag(out,0) - mat(colon(),2));
            }
            {
                Matrix out	= diags(mat(colon(),colon(1,2)),(mat_row(), 0, 1),m,m);
                check_struct(out);

                dif			+= norm_1(get_diag(out,0) - mat(colon(),1));
                dif			+= norm_1(get_diag(out,1) - mat(colon(1,m-1),2));
            }
            {
                Matrix out	= diags(mat(colon(),colon(1,2)),(mat_row(), -1, 1),m,m);
                check_struct(out);

                dif			+= norm_1(get_diag(out,-1) - mat(colon(1,m-1),1));
                dif			+= norm_1(get_diag(out,1) - mat(colon(1,m-1),2));
            };
        };

        return dif;
    }
    catch(const std::exception& ex)
    {
        m_is_error  = true;
        m_error     = ex.what();
        return 1.;
    }
    catch(...)
    {
        m_is_error  = true;
        m_error     = "unknown error";
        return 1.;
    };
};
Real test_function_diags::eval_scalar(const Scalar& ,bool,int code )
{
    (void) code;
    return 0;
};

Real test_function_bdiags::eval_mat(const Matrix& mat,bool,int code )
{
    (void) code;
    try
    {
        Real dif		= 0;
        Matrix tmp		= vec(mat);
        Integer m = mat.rows(), n = mat.cols(); 

        // make dense diagonal only if mat not too big
        if (Real(tmp.numel()) * Real(tmp.numel()) < 10000)
        {
            Matrix out0	    = bdiags(tmp,(mat_row(), 0),m*n,m*n);
            Matrix out02	= spdiags(tmp,(mat_row(), 0),m*n,m*n);
            check_struct(out0);

            dif			+= norm_1(out0 - out02);

            if (1 < m*n)
            {
                Matrix out	= bdiags(tmp,(mat_row(), 1),m*n,m*n);
                Matrix out2	= spdiags(tmp,(mat_row(), 1),m*n,m*n);
                check_struct(out);

                dif			+= norm_1(out - out2);
            };
        
            if (1 < m*n)
            {
                Matrix out	= bdiags(tmp,(mat_row(), -1),m*n,m*n);
                Matrix out2	= spdiags(tmp,(mat_row(), -1),m*n,m*n);
                check_struct(out);

                dif			+= norm_1(out - out2);
            };		
        }
        if (n > 1 && m > 1)
        {
            {
                Matrix out	= bdiags(mat(colon(),colon(1,2)),(mat_row(), -1, 0),m,m);
                Matrix out2	= spdiags(mat(colon(),colon(1,2)),(mat_row(), -1, 0),m,m);
                check_struct(out);

                dif			+= norm_1(out - out2);
            }
            {
                Matrix out	= bdiags(mat(colon(),colon(1,2)),(mat_row(), 0, 1),m,m);
                Matrix out2	= spdiags(mat(colon(),colon(1,2)),(mat_row(), 0, 1),m,m);
                check_struct(out);

                dif			+= norm_1(out - out2);
            }
            {
                Matrix out	= bdiags(mat(colon(),colon(1,2)),(mat_row(), -1, 1),m,m);
                Matrix out2	= spdiags(mat(colon(),colon(1,2)),(mat_row(), -1, 1),m,m);
                check_struct(out);

                dif			+= norm_1(out - out2);
            };
        };

        return dif;
    }
    catch(const std::exception& ex)
    {
        m_is_error  = true;
        m_error     = ex.what();
        return 1.;
    }
    catch(...)
    {
        m_is_error  = true;
        m_error     = "unknown error";
        return 1.;
    };
};
Real test_function_bdiags::eval_scalar(const Scalar& ,bool,int code )
{
    (void) code;
    return 0;
};

Real test_function_spdiags::eval_mat(const Matrix& mat,bool,int code )
{
    (void) code;
    try
    {
        Real dif		= 0;
        Matrix tmp		= vec(mat);
        Integer m       = mat.rows();
        Integer n       = mat.cols();

        // make dense diagonal only if mat not too big
        if (Real(tmp.numel()) * Real(tmp.numel()) < 10000)
        {
            Matrix out0	    = spdiags(tmp,(mat_row(), 0),m*n,m*n);
            Matrix out02	= bdiags(tmp,(mat_row(), 0),m*n,m*n);
            check_struct(out0);

            dif			+= norm_1(out0 - out02);

            if (1 < m*n)
            {
                Matrix out	= spdiags(tmp,(mat_row(), 1),m*n,m*n);
                Matrix out2	= bdiags(tmp,(mat_row(), 1),m*n,m*n);
                check_struct(out);

                dif			+= norm_1(out - out2);
            };
            
            if (1 < m*n)
            {
                Matrix out	= spdiags(tmp,(mat_row(), -1),m*n,m*n);
                Matrix out2	= bdiags(tmp,(mat_row(), -1),m*n,m*n);
                check_struct(out);

                dif			+= norm_1(out - out2);
            };		
        }
        if (n > 1 && m > 1)
        {
            {
                Matrix out	= spdiags(mat(colon(),colon(1,2)),(mat_row(), -1, 0),m,m);
                Matrix out2	= bdiags(mat(colon(),colon(1,2)),(mat_row(), -1, 0),m,m);
                check_struct(out);

                dif			+= norm_1(out - out2);
            }
            {
                Matrix out	= spdiags(mat(colon(),colon(1,2)),(mat_row(), 0, 1),m,m);
                Matrix out2	= bdiags(mat(colon(),colon(1,2)),(mat_row(), 0, 1),m,m);
                check_struct(out);

                dif			+= norm_1(out - out2);
            }
            {
                Matrix out	= spdiags(mat(colon(),colon(1,2)),(mat_row(), -1, 1),m,m);
                Matrix out2	= bdiags(mat(colon(),colon(1,2)),(mat_row(), -1, 1),m,m);
                check_struct(out);

                dif			+= norm_1(out - out2);
            };
        };

        return dif;
    }
    catch(const std::exception& ex)
    {
        m_is_error  = true;
        m_error     = ex.what();
        return 1.;
    }
    catch(...)
    {
        m_is_error  = true;
        m_error     = "unknown error";
        return 1.;
    };
};
Real test_function_spdiags::eval_scalar(const Scalar& ,bool,int code )
{
    (void) code;
    return 0;
};

class Test_g0 : public matcl::test_function
{
    public:
        virtual bool eval(Integer v) const 
        {
            return v > 0;
        };
        virtual bool eval(Real v) const
        {
            return v > 0.;
        };
        virtual bool eval(Float v) const
        {
            return v > 0.;
        };
        virtual bool eval(const Complex& v) const
        {
            return v > 0.;
        };
        virtual bool eval(const Float_complex& v) const
        {
            return v > 0.;
        };
        virtual bool eval(const Object& v) const
        {
            return (bool)(v > 0.);
        };
};

class Test_e0 : public matcl::test_function
{
    public:
        virtual bool eval(Integer v) const 
        {
            return v == 0;
        };
        virtual bool eval(Real v) const
        {
            return v == 0.;
        };
        virtual bool eval(Float v) const
        {
            return v == 0.;
        };
        virtual bool eval(const Complex& v) const
        {
            return v == 0.;
        };
        virtual bool eval(const Float_complex& v) const
        {
            return v == 0.;
        };
        virtual bool eval(const Object& v) const
        {
            return (bool)(v == 0.);
        };
};

Real test_function_find::eval_mat(const Matrix& mat,bool,int code )
{
    (void)code;
    Matrix out_full     = find(full(mat));	
    Matrix out_full2    = find(full(mat)>0);	

    Matrix out_fulle    = find(full(mat)==0);
    try
    {		
        Matrix out		= find(mat);
        Matrix out2		= find(mat>0);
        Matrix oute		= find(mat==0);
        check_struct(out);
        check_struct(out2);
        check_struct(oute);

        Real dif		= norm_1(out - out_full);
        dif				+= norm_1(out2 - out_full2);
        dif				+= norm_1(oute - out_fulle);
        dif				+= (bool(all(mat(out))) == false);
        dif				+= (bool(all(mat(out2)>0)) == false);
        dif				+= (bool(all(mat(oute)==0)) == false);

        Matrix m1		= mat;
        Matrix m2		= mat;
        Matrix me		= mat;
        m1(out)			= zeros(0,0);
        m2(out2)		= zeros(0,0);
        me(oute)		= zeros(0,0);

        check_struct(m1);
        check_struct(m2);
        check_struct(me);

        dif				+= (!(any_vec(m1)) == false);
        dif				+= (!(any_vec(m2>0)) == false);
        dif				+= (!(any_vec(me==0)) == false);

        Test_g0 ta;
        Matrix out3		= find(mat,ta);
        check_struct(out3);
        Test_e0 tae;
        Matrix out3e		= find(mat,tae);
        check_struct(out3e);

        dif				+= norm_1(out2 - out3);
        dif				+= norm_1(oute - out3e);
        return dif;
    }
    catch(const std::exception& ex)
    {
        m_is_error  = true;
        m_error     = ex.what();
        return 1.;
    }
    catch(...)
    {
        m_is_error  = true;
        m_error     = "unknown error";
        return 1.;
    };
};
Real test_function_find::eval_scalar(const Scalar& ,bool,int code  )
{
    (void)code;
    return 0.;
};

Real test_function_find2::eval_mat(const Matrix& mat,bool,int code )
{
    (void)code;
    mat_tup_2 out_full  = find2(full(mat));	
    mat_tup_2 out_full2 = find2(full(mat)>0);	
    
    mat_tup_2 out_fulle = find2(full(mat)==0);
    try
    {		
        mat_tup_2 out	= find2(mat);
        mat_tup_2 out2	= find2(mat>0);
        mat_tup_2 oute	= find2(mat==0);

        check_struct(out.get<1>());
        check_struct(out.get<2>());
        
        check_struct(out2.get<1>());
        check_struct(out2.get<2>());
        
        check_struct(oute.get<1>());
        check_struct(oute.get<2>());

        Real dif		= norm_1(out.get<1>() - out_full.get<1>());
        dif				+= norm_1(out.get<2>() - out_full.get<2>());
        dif				+= norm_1(out2.get<1>() - out_full2.get<1>());
        dif				+= norm_1(out2.get<2>() - out_full2.get<2>());
        dif				+= norm_1(oute.get<1>() - out_fulle.get<1>());
        dif				+= norm_1(oute.get<2>() - out_fulle.get<2>());

        Test_g0 ta;
        mat_tup_2 out5	= find2(mat,ta);
        check_struct(out5.get<1>());
        check_struct(out5.get<2>());
        
        Test_e0 tae;
        mat_tup_2 out5e	= find2(mat,tae);
        check_struct(out5e.get<1>());
        check_struct(out5e.get<2>());

        dif				+= norm_1(out2.get<1>() - out5.get<1>());
        dif				+= norm_1(out2.get<2>() - out5.get<2>());
        dif				+= norm_1(oute.get<1>() - out5e.get<1>());
        dif				+= norm_1(oute.get<2>() - out5e.get<2>());
        return dif;
    }
    catch(const std::exception& ex)
    {
        m_is_error  = true;
        m_error     = ex.what();
        return 1.;
    }
    catch(...)
    {
        m_is_error  = true;
        m_error     = "unknown error";
        return 1.;
    };
};

Real test_function_find2::eval_scalar(const Scalar& ,bool,int code  )
{
    (void)code;
    return 0.;
};

Real test_function_find3::eval_mat(const Matrix& mat,bool,int code )
{
    (void)code;
    mat_tup_3 out_full  = find3(full(mat));	
    mat_tup_3 out_full2 = find3(full(mat)>0);	
    
    mat_tup_3 out_fulle = find3(full(mat)==0);

    try
    {		
        mat_tup_3 out	= find3(mat);
        mat_tup_3 out2	= find3(mat>0);
        mat_tup_3 oute	= find3(mat==0);

        Matrix out3		= mat(find(mat));
        Matrix out4		= find(mat>0);
        Matrix out4e	= find(mat==0);

        Real dif		= norm_1(out.get<1>() - out_full.get<1>());
        dif				+= norm_1(out.get<2>() - out_full.get<2>());
        dif				+= norm_1(out.get<3>() - out_full.get<3>());
        dif				+= norm_1(out2.get<1>() - out_full2.get<1>());
        dif				+= norm_1(out2.get<2>() - out_full2.get<2>());
        dif				+= norm_1(out2.get<3>() - out_full2.get<3>());
        dif				+= norm_1(out.get<3>() - out3);
        dif				+= norm_1(oute.get<1>() - out_fulle.get<1>());
        dif				+= norm_1(oute.get<2>() - out_fulle.get<2>());
        dif				+= norm_1(oute.get<3>() - out_fulle.get<3>());

        Test_g0 ta;
        mat_tup_3 out5	= find3(mat,ta);
        dif				+= norm_1(out2.get<1>() - out5.get<1>());
        dif				+= norm_1(out2.get<2>() - out5.get<2>());
        dif				+= norm_1(mat(out4) - out5.get<3>());
        Test_e0 tae;
        mat_tup_3 out5e	= find3(mat,tae);
        dif				+= norm_1(oute.get<1>() - out5e.get<1>());
        dif				+= norm_1(oute.get<2>() - out5e.get<2>());
        dif				+= norm_1(mat(out4e) - out5e.get<3>());

        check_struct(out.get<1>());
        check_struct(out.get<2>());
        check_struct(out2.get<1>());
        check_struct(out2.get<2>());
        check_struct(oute.get<1>());
        check_struct(oute.get<2>());
        check_struct(out3);
        check_struct(out4);
        check_struct(out4e);
        check_struct(out5.get<1>());
        check_struct(out5.get<2>());
        check_struct(out5e.get<1>());
        check_struct(out5e.get<2>());

        return dif;
    }
    catch(const std::exception& ex)
    {
        m_is_error  = true;
        m_error     = ex.what();
        return 1.;
    }
    catch(...)
    {
        m_is_error  = true;
        m_error     = "unknown error";
        return 1.;
    };
};
Real test_function_find3::eval_scalar(const Scalar& ,bool,int code  )
{
    (void)code;
    return 0.;
};

Real test_function_sort::eval_mat(const Matrix& mat,bool,int code )
{
    (void)code;
    Matrix out_full = sort(full(mat),d,is_asc);	

    try
    {		
        Matrix out		= sort(mat,d,is_asc);
        check_struct(out);

        Real dif		= norm_1(out - out_full);
        return dif;
    }
    catch(const std::exception& ex)
    {
        m_is_error  = true;
        m_error     = ex.what();
        return 1.;
    }
    catch(...)
    {
        m_is_error  = true;
        m_error     = "unknown error";
        return 1.;
    };
};
Real test_function_sort::eval_scalar(const Scalar& ,bool,int code  )
{
    (void)code;
    return 0.;
};

Real test_function_sort2::eval_mat(const Matrix& mat,bool,int code )
{
    (void)code;
    mat_tup_2 out_full = sort2(full(mat),d,is_asc);	

    try
    {		
        mat_tup_2 out	= sort2(mat,d,is_asc);
        Real dif		= norm_1(out.get<1>() - out_full.get<1>());
        dif				+=norm_1(out.get<2>() - out_full.get<2>());

        check_struct(out.get<1>());
        check_struct(out.get<2>());

        Matrix pos;
        Matrix iout     = out.get<2>();

        if (d == 1)
        {
            pos = iout;
            if (mat.cols() > 0)
                pos = pos + iones(mat.rows(),1)*(irange(0,mat.cols()-1)*mat.rows());
        }
        else
        {
            pos = (iout-1)*mat.rows();
            if (mat.rows() > 0)
                pos = pos + trans(irange(1,mat.rows()))*iones(1,mat.cols());
        };

        dif = dif + norm_1(out.get<1>() - mat(pos));
        dif = dif + norm_1(out.get<1>() - sort(mat,d,is_asc));

        return dif;
    }
    catch(const std::exception& ex)
    {
        m_is_error  = true;
        m_error     = ex.what();
        return 1.;
    }
    catch(...)
    {
        m_is_error  = true;
        m_error     = "unknown error";
        return 1.;
    };
};
Real test_function_sort2::eval_scalar(const Scalar& ,bool,int code)
{
    (void)code;
    return 0.;
};

Real test_function_sortrows::eval_mat(const Matrix& mat,bool,int code )
{
    (void)code;
    Matrix out_full     = sortrows(full(mat));	

    try
    {		
        Matrix out		= sortrows(mat);
        check_struct(out);
        Real dif		= norm_1(out - out_full);
        return dif;
    }
    catch(const std::exception& ex)
    {
        m_is_error  = true;
        m_error     = ex.what();
        return 1.;
    }
    catch(...)
    {
        m_is_error  = true;
        m_error     = "unknown error";
        return 1.;
    };
};
Real test_function_sortrows::eval_scalar(const Scalar& ,bool,int code )
{
    (void)code;
    return 0.;
};

Real test_function_sortrows_dim::eval_mat(const Matrix& mat,bool,int code )
{
    (void)code;
    Integer c		= mat.cols();
    Matrix out1		= sortrows(mat,irange(1,c));
    Matrix out2		= sortrows(mat);

    check_struct(out1);
    check_struct(out2);

    Real dif		= norm_1(out1 - out2);

    Matrix out_full;
    try
    {
        out_full = sortrows(full(mat),m);	
    }
    catch(...)
    {
        return dif;
    };

    try
    {				
        Matrix out		= sortrows(mat,m);		
        check_struct(out);
        dif				+= norm_1(out - out_full); 
        return dif;
    }
    catch(const std::exception& ex)
    {
        m_is_error  = true;
        m_error     = ex.what();
        return 1.;
    }
    catch(...)
    {
        m_is_error  = true;
        m_error     = "unknown error";
        return 1.;
    };
};
Real test_function_sortrows_dim::eval_scalar(const Scalar& ,bool,int code )
{
    (void)code;
    return 0.;
};

Real test_function_sortrows2::eval_mat(const Matrix& mat,bool,int code )
{
    (void)code;
    mat_tup_2 out_full = sortrows2(full(mat));	

    try
    {		
        mat_tup_2 out	= sortrows2(mat);
        check_struct(out.get<1>());
        check_struct(out.get<2>());

        Real dif		= norm_1(out.get<1>() - out_full.get<1>());
        dif				+=norm_1(out.get<2>() - out_full.get<2>());

        Matrix iout		= out.get<2>();
        dif				+= norm_1(out.get<1>() - mat(iout,colon()));
        dif				+= norm_1(out.get<1>() - sortrows(mat));

        return dif;
    }
    catch(const std::exception& ex)
    {
        m_is_error  = true;
        m_error     = ex.what();
        return 1.;
    }
    catch(...)
    {
        m_is_error  = true;
        m_error     = "unknown error";
        return 1.;
    };
};
Real test_function_sortrows2::eval_scalar(const Scalar& ,bool,int code )
{
    (void)code;
    return 0.;
};

Real test_function_sortrows_dim2::eval_mat(const Matrix& mat,bool,int code )
{
    (void)code;
    Integer c		= mat.cols();
    mat_tup_2 out1	= sortrows2(mat,irange(1,c));
    mat_tup_2 out2	= sortrows2(mat);
    Real dif		= norm_1(out1.get<1>() - out2.get<1>());
    dif				+= norm_1(out1.get<2>() - out2.get<2>());

    check_struct(out1.get<1>());
    check_struct(out1.get<2>());
    check_struct(out2.get<1>());
    check_struct(out2.get<2>());

    mat_tup_2 out_full;
    try
    {
        out_full	= sortrows2(full(mat),m);	
    }
    catch(...)
    {
        return dif;
    };

    try
    {		
        mat_tup_2 out	= sortrows2(mat,m);

        check_struct(out.get<1>());
        check_struct(out.get<2>());

        dif				+= norm_1(out.get<1>() - out_full.get<1>());
        dif				+=norm_1(out.get<2>() - out_full.get<2>());

        Matrix iout		= out.get<2>();
        dif				+= norm_1(out.get<1>() - mat(iout,colon()));
        dif				+= norm_1(out.get<1>() - sortrows(mat,m));

        return dif;
    }
    catch(const std::exception& ex)
    {
        m_is_error  = true;
        m_error     = ex.what();
        return 1.;
    }
    catch(...)
    {
        m_is_error  = true;
        m_error     = "unknown error";
        return 1.;
    };
};
Real test_function_sortrows_dim2::eval_scalar(const Scalar& ,bool,int code )
{
    (void)code;
    return 0.;
};

Real test_function_sortcols::eval_mat(const Matrix& mat,bool,int code )
{
    (void)code;
    Matrix out_full     = sortcols(full(mat));	

    try
    {		
        Matrix out		= sortcols(mat);
        check_struct(out);
        Real dif		= norm_1(out - out_full);
        return dif;
    }
    catch(const std::exception& ex)
    {
        m_is_error  = true;
        m_error     = ex.what();
        return 1.;
    }
    catch(...)
    {
        m_is_error  = true;
        m_error     = "unknown error";
        return 1.;
    };

};
Real test_function_sortcols::eval_scalar(const Scalar& ,bool,int code )
{
    (void)code;
    return 0.;
};

Real test_function_sortcols_dim::eval_mat(const Matrix& mat,bool,int code )
{
    (void)code;
    Integer r		= mat.rows();
    Matrix out1		= sortcols(mat,irange(1,r));
    Matrix out2		= sortcols(mat);

    check_struct(out1);
    check_struct(out2);

    Real dif		= norm_1(out1 - out2);

    Matrix out_full;
    try
    {
        out_full    = sortcols(full(mat),m);	
    }
    catch(...)
    {
        return dif;
    };

    try
    {				
        Matrix out		= sortcols(mat,m);		
        check_struct(out);
        dif				+= norm_1(out - out_full); 
        return dif;
    }
    catch(const std::exception& ex)
    {
        m_is_error  = true;
        m_error     = ex.what();
        return 1.;
    }
    catch(...)
    {
        m_is_error  = true;
        m_error     = "unknown error";
        return 1.;
    };
};
Real test_function_sortcols_dim::eval_scalar(const Scalar& ,bool,int code )
{
    (void)code;
    return 0.;
};

Real test_function_sortcols2::eval_mat(const Matrix& mat,bool,int code )
{
    (void)code;
    mat_tup_2 out_full = sortcols2(full(mat));	

    try
    {		
        mat_tup_2 out	= sortcols2(mat);
        check_struct(out.get<1>());
        check_struct(out.get<2>());

        Real dif		= norm_1(out.get<1>() - out_full.get<1>());
        dif				+=norm_1(out.get<2>() - out_full.get<2>());

        Matrix iout		= out.get<2>();
        dif				+= norm_1(out.get<1>() - mat(colon(),iout));
        dif				+= norm_1(out.get<1>() - sortcols(mat));

        return dif;
    }
    catch(const std::exception& ex)
    {
        m_is_error  = true;
        m_error     = ex.what();
        return 1.;
    }
    catch(...)
    {
        m_is_error  = true;
        m_error     = "unknown error";
        return 1.;
    };
};
Real test_function_sortcols2::eval_scalar(const Scalar& ,bool,int code )
{
    (void)code;
    return 0.;
};

Real test_function_sortcols_dim2::eval_mat(const Matrix& mat,bool,int code )
{
    (void)code;
    Integer r		= mat.rows();
    mat_tup_2 out1	= sortcols2(mat,irange(1,r));
    mat_tup_2 out2	= sortcols2(mat);

    check_struct(out1.get<1>());
    check_struct(out1.get<2>());
    check_struct(out2.get<1>());
    check_struct(out2.get<2>());

    Real dif		= norm_1(out1.get<1>() - out2.get<1>());
    dif				+= norm_1(out1.get<2>() - out2.get<2>());

    mat_tup_2 out_full;
    try
    {
        out_full	= sortcols2(full(mat),m);	
    }
    catch(...)
    {
        return dif;
    };

    try
    {		
        mat_tup_2 out	= sortcols2(mat,m);

        check_struct(out.get<1>());
        check_struct(out.get<2>());

        dif				+= norm_1(out.get<1>() - out_full.get<1>());
        dif				+=norm_1(out.get<2>() - out_full.get<2>());

        Matrix iout		= out.get<2>();
        dif			    += norm_1(out.get<1>() - mat(colon(),iout));
        dif				+= norm_1(out.get<1>() - sortcols(mat,m));

        return dif;
    }
    catch(const std::exception& ex)
    {
        m_is_error  = true;
        m_error     = ex.what();
        return 1.;
    }
    catch(...)
    {
        m_is_error  = true;
        m_error     = "unknown error";
        return 1.;
    };
};
Real test_function_sortcols_dim2::eval_scalar(const Scalar& ,bool,int code )
{
    (void)code;
    return 0.;
};

Real test_function_is_sorted::eval_mat(const Matrix& mat,bool,int code )
{
    (void)code;
    Matrix out_full = issorted(full(mat),d,is_asc);	

    try
    {		
        Matrix out	= issorted(mat,d,is_asc);
        check_struct(out);

        Real dif		= norm_1(out - out_full);

        Matrix ms		= sort(mat,d,is_asc);
        Matrix outs		= issorted(ms,d,is_asc);
        check_struct(outs);
        check_struct(ms);

        bool iss		= true;
        if (mat.rows() > 0 && mat.cols() > 0)
            iss			= bool(all(all(outs,1),2));

        dif				+= (iss == false);

        return dif;
    }
    catch(const std::exception& ex)
    {
        m_is_error  = true;
        m_error     = ex.what();
        return 1.;
    }
    catch(...)
    {
        m_is_error  = true;
        m_error     = "unknown error";
        return 1.;
    };
};
Real test_function_is_sorted::eval_scalar(const Scalar& ,bool,int code )
{
    (void)code;
    return 0.;
};

Real test_function_is_sorted_cols::eval_mat(const Matrix& mat,bool,int code )
{
    (void)code;
    bool out_full	= issorted_cols(full(mat));	

    try
    {		
        bool out		= issorted_cols(mat);
        Real dif		= norm_1(out - out_full);

        Matrix ms		= sortcols(mat);
        check_struct(ms);
        bool outs		= issorted_cols(ms);

        dif				+= (outs == false);

        return dif;
    }
    catch(const std::exception& ex)
    {
        m_is_error  = true;
        m_error     = ex.what();
        return 1.;
    }
    catch(...)
    {
        m_is_error  = true;
        m_error     = "unknown error";
        return 1.;
    };
};
Real test_function_is_sorted_cols::eval_scalar(const Scalar& ,bool,int code )
{
    (void)code;
    return 0.;
};

Real test_function_is_sorted_rows::eval_mat(const Matrix& mat,bool,int code )
{
    (void)code;
    bool out_full	    = issorted_rows(full(mat));	

    try
    {		
        bool out		= issorted_rows(mat);
        Real dif		= norm_1(out - out_full);

        Matrix ms		= sortrows(mat);
        check_struct(ms);
        bool outs		= issorted_rows(ms);

        dif				+= (outs == false);

        return dif;
    }
    catch(const std::exception& ex)
    {
        m_is_error  = true;
        m_error     = ex.what();
        return 1.;
    }
    catch(...)
    {
        m_is_error  = true;
        m_error     = "unknown error";
        return 1.;
    };
};
Real test_function_is_sorted_rows::eval_scalar(const Scalar& ,bool,int code )
{
    (void)code;
    return 0.;
};

Real test_function_cat::eval_mat(const Matrix& mat,bool,int code)
{
    dynamic_mat_set ds(m_rand);
    try
    {		
        Integer p = 1;
        mat_col mc;
        mat_col mc_full;
        for (Integer j = 1; j <= m; j++)
        {
            mat_row mr;
            mat_row mr_full;
            for (Integer i = 1; i <= n; i++)
            {
                if (p == pos)
                {
                    (mr, mat);
                    (mr_full, full(mat));
                }
                else
                {
                    Matrix tmp;
                    tmp = ds.rand(mat.rows(), mat.cols(),seed+code);
                    (mr, tmp);
                    (mr_full, full(tmp));
                };
                p++;
            };

            (mc, mr);
            (mc_full, mr_full);
        };
        
        Matrix out		= mc;
        Matrix out_full = mc_full;
        check_struct(out);

        return norm_1(out - out_full);
    }
    catch(const std::exception& ex)
    {
        m_is_error  = true;
        m_error     = ex.what();
        return 1.;
    }
    catch(...)
    {
        m_is_error  = true;
        m_error     = "unknown error";
        return 1.;
    };
};

Real test_function_cat::eval_scalar(const Scalar& ,bool,int code  )
{
    (void)code;
    return 0;
};

Real test_function_get_diag::eval_mat(const Matrix& mat,bool,int code )
{
    (void)code;
    try
    {		
        Matrix out_full = get_diag(full(mat),d);	
        Matrix out		= get_diag(mat,d);
        check_struct(out);
        return norm_1(out - out_full);
    }
    catch(const error::invalid_diag&)
    {
        if (mat.rows()*mat.cols() == 0 || d + 1 > mat.cols() || 1 - d > mat.rows())
            return 0.;
        else
            return 1.;
    }
    catch(const std::exception& ex)
    {
        m_is_error  = true;
        m_error     = ex.what();
        return 1.;
    }
    catch(...)
    {
        m_is_error  = true;
        m_error     = "unknown error";
        return 1.;
    };
};

Real test_function_get_diag::eval_scalar(const Scalar& s,bool,int code  )
{
    (void)code;
    switch(s.get_value_code())
    {
        case value_code::v_integer:
        {
            Integer value = s.get_int();
            return norm_1(diag(value,d) - diag(full(value),d));
        }
        case value_code::v_float:
        {
            Float value		= s.get_float();
            Float_complex val_c	= value;
            Real dif		= norm_1(diag(value,d) - diag(full(value),d));
            dif				+=norm_1(diag(value,d) - diag(val_c,d));
            return dif;
        }
        case value_code::v_real:
        {
            Real value		= s.get_real();
            Complex val_c	= value;
            Real dif		= norm_1(diag(value,d) - diag(full(value),d));
            dif				+=norm_1(diag(value,d) - diag(val_c,d));
            return dif;
        }
        case value_code::v_float_complex:
        {
            Float_complex value = s.get_fcomplex();
            return norm_1(diag(value,d) - diag(full(value),d));
        }
        case value_code::v_complex:
        {
            Complex value = s.get_complex();
            return norm_1(diag(value,d) - diag(full(value),d));
        }
        default:
        {
            return 0;
        }
    };
};

Real test_function_tril::eval_mat(const Matrix& mat,bool,int code )
{
    (void)code;
    Matrix out_full = tril(full(mat),d);	
    try
    {		
        Matrix out		= tril(mat,d);
        check_struct(out);
        return norm_1(out - out_full);
    }
    catch(const std::exception& ex)
    {
        m_is_error  = true;
        m_error     = ex.what();
        return 1.;
    }
    catch(...)
    {
        m_is_error  = true;
        m_error     = "unknown error";
        return 1.;
    };
};
Real test_function_tril::eval_scalar(const Scalar& s,bool,int code  )
{
    (void)code;
    try
    {
        tril(full(1),d);
    }
    catch(...)
    {
        return 0;
    };

    switch(s.get_value_code())
    {
        case value_code::v_integer:
        {
            Integer value = s.get_int();
            return norm_1(tril(value,d) - tril(full(value),d));
        }
        case value_code::v_float:
        {
            Float value		= s.get_float();
            Float_complex val_c	= value;
            Real dif		= norm_1(tril(value,d) - tril(full(value),d));
            dif				+=norm_1(tril(value,d) - tril(val_c,d));
            return dif;
        }
        case value_code::v_real:
        {
            Real value		= s.get_real();
            Complex val_c	= value;
            Real dif		= norm_1(tril(value,d) - tril(full(value),d));
            dif				+=norm_1(tril(value,d) - tril(val_c,d));
            return dif;
        }
        case value_code::v_float_complex:
        {
            Float_complex value = s.get_fcomplex();
            return norm_1(tril(value,d) - tril(full(value),d));
        }
        case value_code::v_complex:
        {
            Complex value = s.get_complex();
            return norm_1(tril(value,d) - tril(full(value),d));
        }
        default:
        {
            return 0;
        }
    };
};

Real test_function_triu::eval_mat(const Matrix& mat,bool,int code )
{
    (void)code;
    Matrix out_full = triu(full(mat),d);
    try
    {		
        Matrix out		= triu(mat,d);        
        check_struct(out);
        return norm_1(out - out_full);
    }
    catch(const std::exception& ex)
    {
        m_is_error  = true;
        m_error     = ex.what();
        return 1.;
    }
    catch(...)
    {
        m_is_error  = true;
        m_error     = "unknown error";
        return 1.;
    };
};
Real test_function_triu::eval_scalar(const Scalar& s,bool,int code  )
{
    (void)code;
    try
    {
        triu(full(1),d);
    }
    catch(...)
    {
        return 0;
    };

    switch(s.get_value_code())
    {
        case value_code::v_integer:
        {
            Integer value = s.get_int();
            return norm_1(triu(value,d) - triu(full(value),d));
        }
        case value_code::v_float:
        {
            Float value		= s.get_float();
            Float_complex val_c	= value;
            Real dif		= norm_1(triu(value,d) - triu(full(value),d));
            dif				+=norm_1(triu(value,d) - triu(val_c,d));
            return dif;
        }
        case value_code::v_real:
        {
            Real value		= s.get_real();
            Complex val_c	= value;
            Real dif		= norm_1(triu(value,d) - triu(full(value),d));
            dif				+=norm_1(triu(value,d) - triu(val_c,d));
            return dif;
        }
        case value_code::v_float_complex:
        {
            Float_complex value = s.get_fcomplex();
            return norm_1(triu(value,d) - triu(full(value),d));
        }
        case value_code::v_complex:
        {
            Complex value = s.get_complex();
            return norm_1(triu(value,d) - triu(full(value),d));
        }
        default:
        {
            return 0;
        }
    };
};

Real test_function_rot90::eval_mat(const Matrix& mat,bool,int code )
{
    (void)code;
    Matrix out_full = rot90(full(mat),d);	
    try
    {		
        Matrix out		= rot90(mat,d);
        check_struct(out);
        return norm_1(out - out_full);
    }
    catch(const std::exception& ex)
    {
        m_is_error  = true;
        m_error     = ex.what();
        return 1.;
    }
    catch(...)
    {
        m_is_error  = true;
        m_error     = "unknown error";
        return 1.;
    };
};

Real test_function_rot90::eval_scalar(const Scalar& s,bool,int code  )
{
    (void)code;
    switch(s.get_value_code())
    {
        case value_code::v_integer:
        {
            Integer value = s.get_int();
            return norm_1(rot90(value,d) - rot90(full(value),d));
        }
        case value_code::v_float:
        {
            Float value		= s.get_float();
            Float_complex val_c	= value;
            Real dif		= norm_1(rot90(value,d) - rot90(full(value),d));
            dif				+=norm_1(rot90(value,d) - rot90(val_c,d));
            return dif;
        }
        case value_code::v_float_complex:
        {
            Float_complex value = s.get_fcomplex();
            return norm_1(rot90(value,d) - rot90(full(value),d));
        }
        case value_code::v_real:
        {
            Real value		= s.get_real();
            Complex val_c	= value;
            Real dif		= norm_1(rot90(value,d) - rot90(full(value),d));
            dif				+=norm_1(rot90(value,d) - rot90(val_c,d));
            return dif;
        }
        case value_code::v_complex:
        {
            Complex value = s.get_complex();
            return norm_1(rot90(value,d) - rot90(full(value),d));
        }
        default:
        {
            return 0;
        }
    };
};

Real test_function_reshape::eval_mat(const Matrix& mat,bool,int code )
{
    (void)code;
    Integer r = mat.rows();
    Integer c = mat.cols();
    Matrix mf = full(mat);
    Matrix out_full_1 = reshape(mf,r,c);	
    Matrix out_full_2 = reshape(mf,c,r);	
    Matrix out_full_3 = reshape(mf,r*c,1);	
    Matrix out_full_4 = reshape(mf,1,r*c);	
    Matrix out_full_5;
    Matrix out_full_6;

    if (r % 2 == 0)
    {
        out_full_5 = reshape(mf,r/2,c*2);	
    }
    if (c%2 == 0)
    {
        out_full_6 = reshape(mf,r*2,c/2);	
    }

    try
    {		
        Matrix out_1 = reshape(mat,r,c);	
        Matrix out_2 = reshape(mat,c,r);	
        Matrix out_3 = reshape(mat,r*c,1);	
        Matrix out_4 = reshape(mat,1,r*c);	

        Matrix out_12 = reshape(out_1,r,c);	
        Matrix out_22 = reshape(out_2,r,c);	
        Matrix out_32 = reshape(out_3,r,c);	
        Matrix out_42 = reshape(out_4,r,c);	
        Matrix out_5, out_52;
        Matrix out_6, out_62;

        if (r % 2 == 0)
        {
            out_5   = reshape(mat,r/2,c*2);	
            out_52  = reshape(out_5,r,c);
        }
        if (c % 2 == 0)
        {
            out_6 = reshape(mat,r*2,c/2);	
            out_62  = reshape(out_6,r,c);
        }

        check_struct(out_1);
        check_struct(out_2);
        check_struct(out_3);
        check_struct(out_4);
        check_struct(out_5);
        check_struct(out_6);

        Real dif		= 0;
        dif				+=norm_1(out_1 - out_full_1);
        dif				+=norm_1(out_2 - out_full_2);
        dif				+=norm_1(out_3 - out_full_3);
        dif				+=norm_1(out_4 - out_full_4);
        dif			    +=norm_1(out_5 - out_full_5);
        dif			    +=norm_1(out_6 - out_full_6);

        dif				+=norm_1(mat - out_12);
        dif				+=norm_1(mat - out_22);
        dif				+=norm_1(mat - out_32);
        dif				+=norm_1(mat - out_42);

        if (r % 2 == 0)
            dif			+=norm_1(mat - out_52);

        if (c % 2 == 0)
            dif		    +=norm_1(mat - out_62);

        return dif;
    }
    catch(const std::exception& ex)
    {
        m_is_error  = true;
        m_error     = ex.what();
        return 1.;
    }
    catch(...)
    {
        m_is_error  = true;
        m_error     = "unknown error";
        return 1.;
    };
};
Real test_function_reshape::eval_scalar(const Scalar& ,bool,int code  )
{
    (void)code;
    return 0.;
};

Real test_function_sparse::eval_mat(const Matrix& mat,bool,int code )
{
    (void)code;
    try
    {		
        Real dif		= norm_1(mat - sparse(mat));
        return dif;
    }
    catch(const std::exception& ex)
    {
        m_is_error  = true;
        m_error     = ex.what();
        return 1.;
    }
    catch(...)
    {
        m_is_error  = true;
        m_error     = "unknown error";
        return 1.;
    };
};
Real test_function_sparse::eval_scalar(const Scalar& ,bool,int code )
{
    (void)code;
    return 0.;
};

Real test_function_full::eval_mat(const Matrix& mat,bool,int code )
{
    (void)code;
    try
    {		
        Matrix mat_f    = full(mat);
        Real dif		= norm_1(mat - mat_f);
        if (mat_f.get_struct_code() != struct_code::struct_dense)
            dif         += 1.0;

        return dif;
    }
    catch(const std::exception& ex)
    {
        m_is_error  = true;
        m_error     = ex.what();
        return 1.;
    }
    catch(...)
    {
        m_is_error  = true;
        m_error     = "unknown error";
        return 1.;
    };
};
Real test_function_full::eval_scalar(const Scalar& ,bool,int code )
{
    (void)code;
    return 0.;
};
Real test_function_band::eval_mat(const Matrix& mat,bool,int code )
{
    (void)code;
    try
    {		
        Real dif		= norm_1(mat - band(mat));
        return dif;
    }
    catch(const std::exception& ex)
    {
        m_is_error  = true;
        m_error     = ex.what();
        return 1.;
    }
    catch(...)
    {
        m_is_error  = true;
        m_error     = "unknown error";
        return 1.;
    };
};
Real test_function_band::eval_scalar(const Scalar& ,bool,int code )
{
    (void)code;
    return 0.;
};

Real test_function_clone::eval_mat(const Matrix& mat,bool,int code )
{
    (void)code;
    try
    {	
        Real dif		= 0;
        {
            Matrix A		= clone(mat);
            check_struct(A);

            dif				= norm_1(mat - A);
            dif				+= (A.is_unique()==false);
        };
        {
            Matrix B		= mat;
            Matrix A		= clone(mat);
            check_struct(A);
            dif				= norm_1(mat - A);
            dif				+= (A.is_unique()==false);
        };
        {
            Matrix B		= mat;
            Matrix A		= clone(B);
            check_struct(A);
            dif				= norm_1(mat - A);
            dif				+= (A.is_unique()==false);
        };
        return dif;
    }
    catch(const std::exception& ex)
    {
        m_is_error  = true;
        m_error     = ex.what();
        return 1.;
    }
    catch(...)
    {
        m_is_error  = true;
        m_error     = "unknown error";
        return 1.;
    };
};
Real test_function_clone::eval_scalar(const Scalar& ,bool,int code )
{
    (void) code;
    return 0.;
};

Real test_function_repmat::eval_mat(const Matrix& mat,bool,int code )
{
    (void)code;
    Matrix out_full = repmat(full(mat),m,n);	
    try
    {		
        Matrix out		= repmat(mat,m,n);
        check_struct(out);
        return norm_1(out - out_full);
    }
    catch(const std::exception& ex)
    {
        m_is_error  = true;
        m_error     = ex.what();
        return 1.;
    }
    catch(...)
    {
        m_is_error  = true;
        m_error     = "unknown error";
        return 1.;
    };
};
Real test_function_repmat::eval_scalar(const Scalar& ,bool,int code  )
{
    (void)code;
    return 0.;
};

Real test_function_convert::eval_mat(const Matrix& mat,bool,int )
{
    auto ew = matcl::error::enable_warnings(false);

    Matrix out_full;
    try
    {
        Matrix mf   = full(mat);
        out_full    = convert(mf,(matcl::mat_code)code);	
    }
    catch(...)
    {
        return 0.;
    };

    try
    {		
        Matrix out		= convert(mat,(matcl::mat_code)code);
        check_struct(out);

        Real dif		= norm_1(out - out_full);

        value_code v1   = matcl::matrix_traits::get_value_type((matcl::mat_code)code);
        value_code v2   = matcl::matrix_traits::get_value_type(mat.get_matrix_code());
        value_code v12  = matcl::matrix_traits::unify_value_types(v1, v2);

        if (v1 == v12)
            dif			+=norm_1(out - mat);

        return dif;
    }
    catch(const std::exception& ex)
    {
        m_is_error  = true;
        m_error     = ex.what();
        return 1.;
    }
    catch(...)
    {
        m_is_error  = true;
        m_error     = "unknown error";
        return 1.;
    };
};
Real test_function_convert::eval_scalar(const Scalar& ,bool,int code_arg )
{
    (void)code_arg;
    return 0.;
};

Real test_function_convert_val::eval_mat(const Matrix& mat,bool,int )
{
    auto ew = matcl::error::enable_warnings(false);

    try
    {		
        Matrix out		= convert_value(mat, (matcl::value_code)code);
        check_struct(out);

        mat_code mc     = matrix_traits::get_matrix_type((matcl::value_code)code, mat.get_struct_code());
        Matrix out2     = convert(mat, mc);

        Real dif		= norm_1(out - out2);
        if (out.get_value_code() != (matcl::value_code)code)
            dif         += 1.0;
        if (out.get_struct_code() != mat.get_struct_code())
            dif         += 1.0;

        return dif;
    }
    catch(const std::exception& ex)
    {
        m_is_error  = true;
        m_error     = ex.what();
        return 1.;
    }
    catch(...)
    {
        m_is_error  = true;
        m_error     = "unknown error";
        return 1.;
    };
};
Real test_function_convert_val::eval_scalar(const Scalar& ,bool,int code_arg )
{
    (void)code_arg;
    return 0.;
};

Real test_function_get_lu::eval_mat(const Matrix& mat,bool,int code )
{
    (void)code;
    try
    {
        //Integer ld = mat.ldiags();
        //Integer ud = mat.udiags();

        Matrix mat_sp   = sparse(mat);
        Matrix mat_d    = full(mat);
        Matrix mat_bd   = band(mat);

        Integer ld_sp   = get_ld(mat_sp, -1);
        Integer ud_sp   = get_ud(mat_sp, -1);    

        Integer ld_d    = get_ld(mat_d, -1);
        Integer ud_d    = get_ud(mat_d, -1);

        Integer ld_bd   = get_ld(mat_bd, -1);
        Integer ud_bd   = get_ud(mat_bd, -1);

        if (ld_sp != ld_d || ld_bd != ld_d)
            return 1.;

        if (ud_sp != ud_d || ud_bd != ud_d)
            return 1.;

        return 0.;
    }
    catch(const std::exception& ex)
    {
        m_is_error  = true;
        m_error     = ex.what();
        return 1.;
    }
    catch(...)
    {
        m_is_error  = true;
        m_error     = "unknown error";
        return 1.;
    };    
};
Real test_function_get_lu::eval_scalar(const Scalar& ,bool,int code )
{
    (void)code;
    return 0.;
};

Real test_function_trans_t::eval_mat(const Matrix& mat,bool,int code )
{
    (void)code;
    try
    {
        Matrix m1       = trans(mat, trans_type::no_trans);
        Matrix m2       = trans(mat, trans_type::trans);
        Matrix m3       = trans(mat, trans_type::conj_trans);

        Real out        = norm_1(m1 - mat);
        out             += norm_1(m2 - trans(mat));
        out             += norm_1(m3 - ctrans(mat));

        return out;
    }
    catch(const std::exception& ex)
    {
        m_is_error  = true;
        m_error     = ex.what();
        return 1.;
    }
    catch(...)
    {
        m_is_error  = true;
        m_error     = "unknown error";
        return 1.;
    };    
};
Real test_function_trans_t::eval_scalar(const Scalar& ,bool,int code )
{
    (void)code;
    return 0.;
};

Real test_function_matrix_fwd::eval_mat(const Matrix& mat,bool,int code )
{
    (void)code;
    try
    {
        Real out        = 0.0;
        out             += norm_1(mat.rows() - rows(mat));
        out             += norm_1(mat.cols() - cols(mat));
        out             += norm_1(mat.structural_ldiags(true) - structural_ldiags(mat, true));
        out             += norm_1(mat.structural_ldiags(false) - structural_ldiags(mat, false));
        out             += norm_1(mat.structural_udiags(true) - structural_udiags(mat, true));
        out             += norm_1(mat.structural_udiags(false) - structural_udiags(mat, false));
        out             += norm_1(mat.structural_nnz() - structural_nnz(mat));
        out             += norm_1(mat.length() - length(mat));
        out             += norm_1(mat.numel() - numel(mat));
        out             += norm_1(mat.all_finite() - all_finite(mat));        
        
        return out;
    }
    catch(const std::exception& ex)
    {
        m_is_error  = true;
        m_error     = ex.what();
        return 1.;
    }
    catch(...)
    {
        m_is_error  = true;
        m_error     = "unknown error";
        return 1.;
    };    
};
Real test_function_matrix_fwd::eval_scalar(const Scalar& ,bool,int code )
{
    (void)code;
    return 0.;
};

Real test_function_horzcat::eval_mat(const Matrix& mat,bool,int code )
{
    (void)code;
    try
    {
        Matrix out1     = horzcat(mat,mat);
        Matrix out2     = horzcat({mat,mat,mat});
        Matrix out3     = horzcat({mat,mat});
        Matrix out4     = horzcat({mat});

        Matrix t2       = mat_row().add(mat).add(mat);
        Matrix t3       = mat_row().add(mat).add(mat).add(mat);

        Real out        = 0.0;
        out             += norm_1(out1 - t2);
        out             += norm_1(out2 - t3);
        out             += norm_1(out3 - t2);
        out             += norm_1(out4 - mat);

        check_struct(out1);
        check_struct(out2);
        check_struct(out3);
        check_struct(out4);

        return out;
    }
    catch(const std::exception& ex)
    {
        m_is_error  = true;
        m_error     = ex.what();
        return 1.;
    }
    catch(...)
    {
        m_is_error  = true;
        m_error     = "unknown error";
        return 1.;
    };    
};
Real test_function_horzcat::eval_scalar(const Scalar& ,bool,int code )
{
    (void)code;
    return 0.;
};

Real test_function_vertcat::eval_mat(const Matrix& mat,bool,int code )
{
    (void)code;
    try
    {
        Matrix out1     = vertcat(mat,mat);
        Matrix out2     = vertcat({mat,mat,mat});
        Matrix out3     = vertcat({mat,mat});
        Matrix out4     = vertcat({mat});

        Matrix t2       = mat_col().add(mat).add(mat);
        Matrix t3       = mat_col().add(mat).add(mat).add(mat);

        Real out        = 0.0;
        out             += norm_1(out1 - t2);
        out             += norm_1(out2 - t3);
        out             += norm_1(out3 - t2);
        out             += norm_1(out4 - mat);

        check_struct(out1);
        check_struct(out2);
        check_struct(out3);
        check_struct(out4);

        return out;
    }
    catch(const std::exception& ex)
    {
        m_is_error  = true;
        m_error     = ex.what();
        return 1.;
    }
    catch(...)
    {
        m_is_error  = true;
        m_error     = "unknown error";
        return 1.;
    };    
};
Real test_function_vertcat::eval_scalar(const Scalar& ,bool,int code )
{
    (void)code;
    return 0.;
};

Real test_function_blkdiag::eval_mat(const Matrix& mat,bool,int code )
{
    (void)code;
    try
    {
        Matrix out1     = blkdiag(mat,mat);
        Matrix out2     = blkdiag({mat,mat,mat});
        Matrix out3     = blkdiag({mat,mat});
        Matrix out4     = blkdiag({mat});

        Integer r       = mat.rows();
        Integer c       = mat.cols();
        Matrix t2       = mat_col().add(mat_row().add(mat).add(zeros(r,c)))
                                   .add(mat_row().add(zeros(r,c)).add(mat));
        Matrix t3       = blkdiag(t2, full(mat));

        Real out        = 0.0;
        out             += norm_1(out1 - t2);
        out             += norm_1(out2 - t3);
        out             += norm_1(out3 - t2);
        out             += norm_1(out4 - mat);

        check_struct(out1);
        check_struct(out2);
        check_struct(out3);
        check_struct(out4);

        return out;
    }
    catch(const std::exception& ex)
    {
        m_is_error  = true;
        m_error     = ex.what();
        return 1.;
    }
    catch(...)
    {
        m_is_error  = true;
        m_error     = "unknown error";
        return 1.;
    };    
};
Real test_function_blkdiag::eval_scalar(const Scalar& ,bool,int code )
{
    (void)code;
    return 0.;
};

Real test_function_select_band::eval_mat(const Matrix& mat,bool,int code )
{
    (void)code;
    try
    {
        Matrix out1     = select_band(mat,m_fd, m_ld);
        Matrix out2     = tril(triu(mat, m_fd), m_ld);

        Real out        = 0.0;
        out             += norm_1(out1 - out2);

        check_struct(out1);
        check_struct(out2);

        return out;
    }
    catch(const std::exception& ex)
    {
        m_is_error  = true;
        m_error     = ex.what();
        return 1.;
    }
    catch(...)
    {
        m_is_error  = true;
        m_error     = "unknown error";
        return 1.;
    };    
};
Real test_function_select_band::eval_scalar(const Scalar& ,bool,int code )
{
    (void)code;
    return 0.;
};

Real test_function_checks::eval_mat(const Matrix& mat,bool,int code )
{
    double tol  = 0.0;

    (void)code;
    try
    {
        Real out        = 0.0;

        {
            bool is     = is_tril(mat);
            bool is2    = is_tril(tril(mat,0));

            if (is2 == false)
                out     += 1;

            if (is == true && norm_1(mat - tril(mat,0)) > 0.0)
                out     += 1;
        }

        {
            bool is     = is_triu(mat);
            bool is2    = is_triu(triu(mat,0));

            if (is2 == false)
                out     += 1;

            if (is == true && norm_1(mat - triu(mat,0)) > 0.0)
                out     += 1;
        }
        {
            bool is     = is_diag(mat);
            bool is2    = is_diag(select_band(mat,0,0));

            if (is2 == false)
                out     += 1;

            if (is == true && norm_1(mat - select_band(mat,0,0)) > 0.0)
                out     += 1;
        }
        if (mat.rows() == mat.cols())
        {
            bool is     = is_sym(mat, tol);
            bool is2    = is_sym(mat + trans(mat), tol);

            if (is2 == false && mat.all_finite() == true)
                out     += 1;

            if (is == true && norm_1(mat - trans(mat)) > 0.0)
                out     += 1;
        }
        else
        {
            bool is     = is_sym(mat, tol);

            if (is)
                out     += 1;
        };
        
        if (mat.rows() == mat.cols())
        {
            bool is     = is_her(mat, tol);
            bool is2    = is_her(mat + ctrans(mat), tol);

            if (is2 == false && mat.all_finite() == true)
                out     += 1;

            if (is == true && norm_1(mat - ctrans(mat)) > 0.0)
                out     += 1;
        }
        else
        {
            bool is     = is_her(mat, tol);

            if (is)
                out     += 1;
        };
        {
            bool is     = is_real_matrix(mat);

            if (is == true)
            {
                if (matrix_traits::is_real(mat.get_value_code()) == false)
                    out += 1;
            }
            else
            {
                if (matrix_traits::is_real(mat.get_value_code()) == true)
                    out += 1;
            }
        }
        {
            bool is     = is_real_float_matrix(mat);

            if (is == true)
            {
                if (matrix_traits::is_float_real(mat.get_value_code()) == false)
                    out += 1;
            }
            else
            {
                if (matrix_traits::is_float_real(mat.get_value_code()) == true)
                    out += 1;
            }
        }
        {
            bool is     = is_complex_matrix(mat);

            if (is == true)
            {
                if (matrix_traits::is_float_complex(mat.get_value_code()) == false)
                    out += 1;
            }
            else
            {
                if (matrix_traits::is_float_complex(mat.get_value_code()) == true)
                    out += 1;
            }
        }
        {
            bool is     = is_integer_matrix(mat);

            if (is == true)
            {
                if (mat.get_value_code() != value_code::v_integer)
                    out += 1;
            }
            else
            {
                if (mat.get_value_code() == value_code::v_integer)
                    out += 1;
            }
        }
        {
            bool is     = is_object_matrix(mat);

            if (is == true)
            {
                if (mat.get_value_code() != value_code::v_object)
                    out += 1;
            }
            else
            {
                if (mat.get_value_code() == value_code::v_object)
                    out += 1;
            }
        }
        {
            bool is     = is_dense_matrix(mat);

            if (is == true)
            {
                if (mat.get_struct_code() != struct_code::struct_dense)
                    out += 1;
            }
            else
            {
                if (mat.get_struct_code() == struct_code::struct_dense)
                    out += 1;
            }
        }

        {
            bool is     = is_sparse_matrix(mat);

            if (is == true)
            {
                if (mat.get_struct_code() != struct_code::struct_sparse)
                    out += 1;
            }
            else
            {
                if (mat.get_struct_code() == struct_code::struct_sparse)
                    out += 1;
            }
        }
        {
            bool is     = is_band_matrix(mat);

            if (is == true)
            {
                if (mat.get_struct_code() != struct_code::struct_banded)
                    out += 1;
            }
            else
            {
                if (mat.get_struct_code() == struct_code::struct_banded)
                    out += 1;
            }
        }
        {
            bool is     = is_scalar_matrix(mat);

            if (is == true)
            {
                if (mat.get_struct_code() != struct_code::struct_scalar)
                    out += 1;
            }
            else
            {
                if (mat.get_struct_code() == struct_code::struct_scalar)
                    out += 1;
            }
        }

        if (mat.is_scalar_type() == false)
        {
            bool is     = is_same_matrix(mat,mat);

            if (is == false)
                out += 1;
        }

        return out;
    }
    catch(const std::exception& ex)
    {
        m_is_error  = true;
        m_error     = ex.what();
        return 1.;
    }
    catch(...)
    {
        m_is_error  = true;
        m_error     = "unknown error";
        return 1.;
    };    
};
Real test_function_checks::eval_scalar(const Scalar& ,bool,int code )
{
    (void)code;
    return 0.;
};

Real test_function_nnz_m::eval_mat(const Matrix& mat,bool,int code )
{
    (void)code;
    try
    {
        Integer val     = nnz(mat);
        
        Matrix Ir, Ic, Ix;
        tie(Ir, Ic, Ix) = find3(mat);
        Integer val2    = Ir.length();

        if (val != val2)
            return 1.0;
        else
            return 0.0;
    }
    catch(const std::exception& ex)
    {
        m_is_error  = true;
        m_error     = ex.what();
        return 1.;
    }
    catch(...)
    {
        m_is_error  = true;
        m_error     = "unknown error";
        return 1.;
    };    
};
Real test_function_nnz_m::eval_scalar(const Scalar& ,bool,int code )
{
    (void)code;
    return 0.;
};

Real test_function_drop_sparse::eval_mat(const Matrix& mat,bool,int code )
{
    (void)code;
    try
    {
        Real tol        = 1e-1;
        bool is_sparse  = is_sparse_matrix(mat);
        Matrix out      = drop_sparse(mat, tol);

        Real dif        = 0.0;

        if (is_sparse == false)
        {
            if (norm_1(mat - out) != 0.0)
                dif     += 1;

            return dif;
        };
        
        Matrix Ir, Ic, Ix;
        tie(Ir, Ic, Ix) = find3(abs(mat) > tol);
        Matrix out2     = make_sparse_matrix(Ir, Ic, Ix, mat.rows(), mat.cols());
        Matrix out3     = mul(mat, out2);

        dif             += norm_1(out3 - out);
        return dif;
    }
    catch(const std::exception& ex)
    {
        m_is_error  = true;
        m_error     = ex.what();
        return 1.;
    }
    catch(...)
    {
        m_is_error  = true;
        m_error     = "unknown error";
        return 1.;
    };    
};
Real test_function_drop_sparse::eval_scalar(const Scalar& ,bool,int code )
{
    (void)code;
    return 0.;
};

};};
