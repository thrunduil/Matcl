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

#pragma once

#include "test_set.h"
#include "test_functions_assign.h"

#include "test/test_matcl/framework/matrix_set/matrix_set_1.h"
#include "matcl-core/IO/logger.h"
#include "test/test_matcl/framework/data_struct/test_selector.h"
#include "test/test_matcl/framework/matrix_set/colon_set.h"
#include "test/test_matcl/framework/matrix_set/test_options.h"

#include "matcl-matrep/matcl_matrep.h"

#include <boost/thread.hpp>

#pragma warning(push)
#pragma warning(disable:4702)
#include <boost/lexical_cast.hpp>
#pragma warning(pop)

namespace matcl { namespace test
{

static Integer num_colons()
{
    return test_options::num_colons();
};

static bool is_valid_diag(Integer d, const Matrix& mat)
{
    Integer r   = mat.rows();
    Integer c   = mat.cols();

    if (d == 0)
        return true;

    if (r == 0 || c == 0)
        return false;

    if (d + 1 <= c && 1 - d <= r)
        return true;

    return false;
};

static bool is_colon_valid(const colon& c1, Integer size)
{
    Integer rows    = size;
    Integer cols    = 1;

    Integer st      = 0;
    Integer en      = 0;

    switch (c1.m_flag)
    {
        case colon::colon_type::t_all:            
            return true;
        case colon::colon_type::t_last:            
            return size == 0 ? false : true;
        case colon::colon_type::t_end_end:
            st = c1.m_s + size;
            en = c1.m_e + size;
            break;
        case colon::colon_type::t_end:
            st = c1.m_s;
            en = c1.m_e + size;
            break;
        case colon::colon_type::t_rev_end:
            st = c1.m_s + size;
            en = c1.m_e;
            break;
        case colon::colon_type::t_one:
        case colon::colon_type::t_range_simple:
        case colon::colon_type::t_range_mat:
            st = c1.m_s;
            en = c1.m_e;
            break;
        case colon::colon_type::t_matrix1:            
        {
            if (c1.m_mat_size.m_imat_12->is_empty())
                return true;

            return (max_d(max_d(*c1.m_mat_size.m_imat_12, 1),2).get_scalar<Integer>() <= size &&
                    min_d(min_d(*c1.m_mat_size.m_imat_12, 1),2).get_scalar<Integer>() > 0) ? true : false;
        }
        case colon::colon_type::t_matrix2:
        {
            bool valid_1    = false;
            bool valid_2    = false;

            if (c1.m_mat_size.m_imat_12[0].is_empty())
            {
                valid_1 = true;
            }
            else
            {
                valid_1 = max_d(max_d(c1.m_mat_size.m_imat_12[0], 1),2).get_scalar<Integer>() <= rows 
                        && min_d(min_d(c1.m_mat_size.m_imat_12[0], 1),2).get_scalar<Integer>() > 0
                        ? true : false;
            };

            if (c1.m_mat_size.m_imat_12[1].is_empty())
            {
                valid_2 = true;
            }
            else
            {
                valid_1 = (max_d(max_d(c1.m_mat_size.m_imat_12[1], 1),2).get_scalar<Integer>() <= cols
                    && min_d(min_d(c1.m_mat_size.m_imat_12[1], 1),2).get_scalar<Integer>() > 0) 
                    ? true : false;
            };
            
            return valid_1 && valid_2;
        }
        default:
            assert(0 && "Impossible colon_type");
    }

    if ((c1.m_i > 0 && st <= en) || (c1.m_i < 0 && st >= en))
        return (std::max(st, en) <= size && std::min(st, en) > 0) ? true : false;
    else 
        return true; 
};

static bool is_colon_all(const colon& c1)
{
    switch (c1.m_flag)
    {
        case colon::colon_type::t_all:            
            return true;
        default:
            return false;
    }
};

class test_assign
{
    assign_functions_list&			tf;
    const test::options&			opts;
    Integer                         thread_id;

    public:
        test_assign(assign_functions_list& tf, const test::options& opts, Integer id)
            :tf(tf),opts(opts), thread_id(id)
        {};

        test_assign(const test_assign& ta)
            :tf(ta.tf),opts(ta.opts), thread_id(ta.thread_id)
        {};

        void make()
        {            
            /*
            Matrix mat = tf.get_matrix(281);
            disp(mat);

            Integer d   = -1;
            mat.diag(d).add_sparse();
            disp(mat);

            check_struct(mat);
            */
            tf.make(opts);
        };
        void operator()()
        {
            make();
        };

    private:		
        test_assign& operator=(const test_assign&) = delete;
};

void test_assign_st(const rand_matrix_ptr& rand)
{
    test::options opts;

    try
    {
        opts.show_partial_res   = 0;
        opts.first_matrix_code  = 0;
        opts.show_memleaks      = true;

        mat_set_1 ms1(rand);
        dynamic_mat_set ms(rand);
        assign_functions_list tf(ms1,ms);

        test_assign ta(tf,opts,0);
        ta.make();
    }
    catch(const std::exception& ex)
    {
        matcl::out_stream <<ex.what();
    };
};

void test_assign_mt(const rand_matrix_ptr& rand)
{
    test::options opts;
    try
    {
        opts.show_partial_res   = 0;
        opts.first_matrix_code  = 0;
        opts.show_memleaks      = false;

        mat_set_2 ms1(rand);
        dynamic_mat_set ms(rand);
        assign_functions_list tf(ms1,ms);

        boost::thread_group tg;

        for (int i = 0; i < 10; i++)
            tg.create_thread(test_assign(tf,opts,i));

        tg.join_all();
    }
    catch(const std::exception& ex)
    {
        matcl::out_stream <<ex.what();
    };
};


static const Integer max_int = 1000;

assign_functions_list::assign_functions_list(const matrix_set& ms1, dynamic_mat_set& ms)
:ms(ms),m_tests(ms1)
{};

void assign_functions_list::make(options opts)
{
    m_options = opts;    

    SELECT_TEST (3, test_get_int());
    SELECT_TEST (3, test_get_int_int());
    SELECT_TEST (3, test_get_colon());
    SELECT_TEST (3, test_get_colon_colon());    
    SELECT_TEST (3, test_get_diag());

    SELECT_TEST (3, test_delcols());
    SELECT_TEST (3, test_delrows());
    SELECT_TEST (3, test_delrowscols());

    SELECT_TEST (3, test_del_int());
    SELECT_TEST (3, test_del_int_int());
    SELECT_TEST (3, test_del_colon());
    SELECT_TEST (3, test_del_colon_colon());    
    SELECT_TEST (3, test_del_diag());
    
    SELECT_TEST (3, test_drop_int());
    SELECT_TEST (3, test_drop_int_int());
    SELECT_TEST (3, test_drop_colon());
    SELECT_TEST (3, test_drop_colon_colon());
    SELECT_TEST (3, test_drop_diag());        

    SELECT_TEST (3, test_set_int_scal());
    SELECT_TEST (3, test_set_int2_scal());
    SELECT_TEST (3, test_set_col_scal());
    SELECT_TEST (3, test_set_col2_scal());
    SELECT_TEST (3, test_set_diag_scal());

    SELECT_TEST (3, test_set_int_mat());
    SELECT_TEST (3, test_set_int2_mat());    
    SELECT_TEST (3, test_set_col_mat());
    SELECT_TEST (3, test_set_col2_mat());    
    SELECT_TEST (3, test_set_diag_mat());    
    
    SELECT_TEST (3, test_drop_sparse_int());
    SELECT_TEST (3, test_drop_sparse_int_int());
    SELECT_TEST (3, test_drop_sparse_col());
    SELECT_TEST (3, test_drop_sparse_col2());
    SELECT_TEST (3, test_drop_sparse_diag());

    SELECT_TEST (3, test_add_sparse_int());
    SELECT_TEST (3, test_add_sparse_int_int());
    SELECT_TEST (3, test_add_sparse_col());
    SELECT_TEST (3, test_add_sparse_col2());
    SELECT_TEST (3, test_add_sparse_diag());
};

Matrix assign_functions_list::get_matrix(int code) const
{
    return m_tests.get_matrix(code);
};

void assign_functions_list::test_get_int()
{
    Real out = 0.;
    {
        test_function_get_int tf(1);
        out += m_tests.make(&tf,m_options);
    };
    {
        test_function_get_int tf(5);
        out += m_tests.make(&tf,m_options);
    };

    if (out == 0.)
        matcl::out_stream << std::string() +   "get_int: OK" + "\n";
    else
        matcl::out_stream << std::string() +   "get_int: FAILED"  + "\n";
};

void assign_functions_list::test_get_int_int()
{
    Real out = 0.;
    {
        test_function_get_int_int tf(1,1);
        out += m_tests.make(&tf,m_options);
    };
    {
        test_function_get_int_int tf(1,3);
        out += m_tests.make(&tf,m_options);
    };
    {
        test_function_get_int_int tf(3,1);
        out += m_tests.make(&tf,m_options);
    };
    {
        test_function_get_int_int tf(3,3);
        out += m_tests.make(&tf,m_options);
    };

    if (out == 0.)
        matcl::out_stream << std::string() +   "get_int_int: OK" + "\n";
    else
        matcl::out_stream << std::string() +   "get_int_int: FAILED"  + "\n";
};

void assign_functions_list::test_get_colon()
{
    Real out = 0.;

    test_function_get_colon tf;
    out += m_tests.make(&tf,m_options);

    if (out == 0.)
        matcl::out_stream << std::string() +   "get_colon: OK" + "\n";
    else
        matcl::out_stream << std::string() +   "get_colon: FAILED"  + "\n";
};

void assign_functions_list::test_get_colon_colon()
{
    Real out = 0.;

    test_function_get_colon_colon tf;
    out += m_tests.make(&tf,m_options);

    if (out == 0.)
        matcl::out_stream << std::string() +   "get_colon_colon: OK" + "\n";
    else
        matcl::out_stream << std::string() +   "get_colon_colon: FAILED"  + "\n";
};

void assign_functions_list::test_get_diag()
{
    Real out = 0.;

    {
        test_function_get_diag_sub tf(-1);
        out += m_tests.make(&tf,m_options);
    }

    {
        test_function_get_diag_sub tf(0);
        out += m_tests.make(&tf,m_options);
    }

    {
        test_function_get_diag_sub tf(1);
        out += m_tests.make(&tf,m_options);
    }

    if (out == 0.)
        matcl::out_stream << std::string() +   "diag_method: OK" + "\n";
    else
        matcl::out_stream << std::string() +   "diag_method: FAILED"  + "\n";
};

void assign_functions_list::test_delcols()
{
    Real out = 0.;

    test_function_delcols tf;
    out += m_tests.make(&tf,m_options);

    if (out == 0.)
        matcl::out_stream << std::string() +   "delcols: OK" + "\n";
    else
        matcl::out_stream << std::string() +   "delcols: FAILED"  + "\n";
};

void assign_functions_list::test_delrows()
{
    Real out = 0.;

    test_function_delrows tf;
    out += m_tests.make(&tf,m_options);

    if (out == 0.)
        matcl::out_stream << std::string() +   "delrows: OK" + "\n";
    else
        matcl::out_stream << std::string() +   "delrows: FAILED"  + "\n";
};

void assign_functions_list::test_delrowscols()
{
    Real out = 0.;

    test_function_delrowscols tf;
    out += m_tests.make(&tf,m_options);

    if (out == 0.)
        matcl::out_stream << std::string() +   "delrowscols: OK" + "\n";
    else
        matcl::out_stream << std::string() +   "delrowscols: FAILED"  + "\n";
};

void assign_functions_list::test_drop_int()
{
    Real out = 0.;
    {
        test_function_drop_int tf(1);
        out += m_tests.make(&tf,m_options);
    };
    {
        test_function_drop_int tf(5);
        out += m_tests.make(&tf,m_options);
    };

    if (out == 0.)
        matcl::out_stream << std::string() +   "drop_int: OK" + "\n";
    else
        matcl::out_stream << std::string() +   "drop_int: FAILED"  + "\n";
};

void assign_functions_list::test_drop_diag()
{
    Real out = 0.;
    {
        test_function_drop_diag tf(-1);
        out += m_tests.make(&tf,m_options);
    };
    {
        test_function_drop_diag tf(0);
        out += m_tests.make(&tf,m_options);
    };
    {
        test_function_drop_diag tf(1);
        out += m_tests.make(&tf,m_options);
    };

    if (out == 0.)
        matcl::out_stream << std::string() +   "drop_diag: OK" + "\n";
    else
        matcl::out_stream << std::string() +   "drop_diag: FAILED"  + "\n";
};

void assign_functions_list::test_drop_int_int()
{
    Real out = 0.;
    {
        test_function_drop_int_int tf(1,1);
        out += m_tests.make(&tf,m_options);
    };
    {
        test_function_drop_int_int tf(1,3);
        out += m_tests.make(&tf,m_options);
    };
    {
        test_function_drop_int_int tf(3,1);
        out += m_tests.make(&tf,m_options);
    };
    {
        test_function_drop_int_int tf(3,3);
        out += m_tests.make(&tf,m_options);
    };
    if (out == 0.)
        matcl::out_stream << std::string() +   "drop_int_int: OK" + "\n";
    else
        matcl::out_stream << std::string() +   "drop_int_int: FAILED"  + "\n";
};

void assign_functions_list::test_drop_colon()
{
    Real out = 0.;

    test_function_drop_colon tf;
    out += m_tests.make(&tf,m_options);

    if (out == 0.)
        matcl::out_stream << std::string() +   "drop_colon: OK" + "\n";
    else
        matcl::out_stream << std::string() +   "drop_colon: FAILED"  + "\n";
};

void assign_functions_list::test_drop_colon_colon()
{
    Real out = 0.;

    test_function_drop_colon_colon tf;
    out += m_tests.make(&tf,m_options);

    if (out == 0.)
        matcl::out_stream << std::string() +   "drop_colon_colon: OK" + "\n";
    else
        matcl::out_stream << std::string() +   "drop_colon_colon: FAILED"  + "\n";
};

//
void assign_functions_list::test_drop_sparse_int()
{
    Real out = 0.;

    for (int i = 0; i <= 2; ++i)
    {
        Real tol = Real(i) - 1.0;
        {
            test_function_drop_sparse_int tf(1, tol);
            out += m_tests.make(&tf,m_options);
        };
        {
            test_function_drop_sparse_int tf(5, tol);
            out += m_tests.make(&tf,m_options);
        };
    };

    if (out == 0.)
        matcl::out_stream << std::string() +   "drop_sparse_int: OK" + "\n";
    else
        matcl::out_stream << std::string() +   "drop_sparse_int: FAILED"  + "\n";
};

void assign_functions_list::test_drop_sparse_diag()
{
    Real out = 0.;

    for (int i = 0; i <= 2; ++i)
    {
        Real tol = Real(i) - 1.0;
        {
            test_function_drop_sparse_diag tf(-1, tol);
            out += m_tests.make(&tf,m_options);
        };
        {
            test_function_drop_sparse_diag tf(0, tol);
            out += m_tests.make(&tf,m_options);
        };
        {
            test_function_drop_sparse_diag tf(1, tol);
            out += m_tests.make(&tf,m_options);
        };
    };

    if (out == 0.)
        matcl::out_stream << std::string() +   "drop_sparse_diag: OK" + "\n";
    else
        matcl::out_stream << std::string() +   "drop_sparse_diag: FAILED"  + "\n";
};

void assign_functions_list::test_drop_sparse_int_int()
{
    Real out = 0.;

    for (int i = 0; i <= 2; ++i)
    {
        Real tol    = Real(i) - 1.0;
        {
            test_function_drop_sparse_int_int tf(1,1, tol);
            out += m_tests.make(&tf,m_options);
        };
        {
            test_function_drop_sparse_int_int tf(1,3, tol);
            out += m_tests.make(&tf,m_options);
        };
        {
            test_function_drop_sparse_int_int tf(3,1, tol);
            out += m_tests.make(&tf,m_options);
        };
        {
            test_function_drop_sparse_int_int tf(3,3, tol);
            out += m_tests.make(&tf,m_options);
        };
    };

    if (out == 0.)
        matcl::out_stream << std::string() +   "drop_sparse_int_int: OK" + "\n";
    else
        matcl::out_stream << std::string() +   "drop_sparse_int_int: FAILED"  + "\n";
};

void assign_functions_list::test_drop_sparse_col()
{
    Real out = 0.;

    {
        test_function_drop_sparse_colon tf(0.0);
        out += m_tests.make(&tf,m_options);
    };
    {
        test_function_drop_sparse_colon tf(1.0);
        out += m_tests.make(&tf,m_options);
    };
    {
        test_function_drop_sparse_colon tf(-1.0);
        out += m_tests.make(&tf,m_options);
    };

    if (out == 0.)
        matcl::out_stream << std::string() +   "drop_sparse_colon: OK" + "\n";
    else
        matcl::out_stream << std::string() +   "drop_sparse_colon: FAILED"  + "\n";
};

void assign_functions_list::test_drop_sparse_col2()
{
    Real out = 0.;

    {
        test_function_drop_sparse_colon_colon tf(0.0);
        out += m_tests.make(&tf,m_options);
    }
    {
        test_function_drop_sparse_colon_colon tf(1.0);
        out += m_tests.make(&tf,m_options);
    }
    {
        test_function_drop_sparse_colon_colon tf(-1.0);
        out += m_tests.make(&tf,m_options);
    }

    if (out == 0.)
        matcl::out_stream << std::string() +   "drop_sparse_colon_colon: OK" + "\n";
    else
        matcl::out_stream << std::string() +   "drop_sparse_colon_colon: FAILED"  + "\n";
};

//
void assign_functions_list::test_add_sparse_int()
{
    Real out = 0.;
    {
        test_function_add_sparse_int tf(1);
        out += m_tests.make(&tf,m_options);
    };
    {
        test_function_add_sparse_int tf(5);
        out += m_tests.make(&tf,m_options);
    };

    if (out == 0.)
        matcl::out_stream << std::string() +   "add_sparse_int: OK" + "\n";
    else
        matcl::out_stream << std::string() +   "drop_sparse_int: FAILED"  + "\n";
};

void assign_functions_list::test_add_sparse_diag()
{
    Real out = 0.;
    {
        test_function_add_sparse_diag tf(-1);
        out += m_tests.make(&tf,m_options);
    };
    {
        test_function_add_sparse_diag tf(0);
        out += m_tests.make(&tf,m_options);
    };
    {
        test_function_add_sparse_diag tf(1);
        out += m_tests.make(&tf,m_options);
    };

    if (out == 0.)
        matcl::out_stream << std::string() +   "add_sparse_diag: OK" + "\n";
    else
        matcl::out_stream << std::string() +   "add_sparse_diag: FAILED"  + "\n";
};

void assign_functions_list::test_add_sparse_int_int()
{
    Real out = 0.;
    {
        test_function_add_sparse_int_int tf(1,1);
        out += m_tests.make(&tf,m_options);
    };
    {
        test_function_add_sparse_int_int tf(1,3);
        out += m_tests.make(&tf,m_options);
    };
    {
        test_function_add_sparse_int_int tf(3,1);
        out += m_tests.make(&tf,m_options);
    };
    {
        test_function_add_sparse_int_int tf(3,3);
        out += m_tests.make(&tf,m_options);
    };
    if (out == 0.)
        matcl::out_stream << std::string() +   "add_sparse_int_int: OK" + "\n";
    else
        matcl::out_stream << std::string() +   "add_sparse_int_int: FAILED"  + "\n";
};

void assign_functions_list::test_add_sparse_col()
{
    Real out = 0.;

    test_function_add_sparse_colon tf;
    out += m_tests.make(&tf,m_options);

    if (out == 0.)
        matcl::out_stream << std::string() +   "add_sparse_colon: OK" + "\n";
    else
        matcl::out_stream << std::string() +   "add_sparse_colon: FAILED"  + "\n";
};

void assign_functions_list::test_add_sparse_col2()
{
    Real out = 0.;

    test_function_add_sparse_colon_colon tf;
    out += m_tests.make(&tf,m_options);

    if (out == 0.)
        matcl::out_stream << std::string() +   "add_sparse_colon_colon: OK" + "\n";
    else
        matcl::out_stream << std::string() +   "add_sparse_colon_colon: FAILED"  + "\n";
};

void assign_functions_list::test_del_int()
{
    Real out = 0.;
    {
        test_function_del_int tf(1);
        out += m_tests.make(&tf,m_options);
    };
    {
        test_function_del_int tf(5);
        out += m_tests.make(&tf,m_options);
    };

    if (out == 0.)
        matcl::out_stream << std::string() +   "del_int: OK" + "\n";
    else
        matcl::out_stream << std::string() +   "del_int: FAILED"  + "\n";
};

void assign_functions_list::test_del_diag()
{
    Real out = 0.;
    {
        test_function_del_diag tf(-1);
        out += m_tests.make(&tf,m_options);
    };
    {
        test_function_del_diag tf(0);
        out += m_tests.make(&tf,m_options);
    };
    {
        test_function_del_diag tf(1);
        out += m_tests.make(&tf,m_options);
    };

    if (out == 0.)
        matcl::out_stream << std::string() +   "del_diag: OK" + "\n";
    else
        matcl::out_stream << std::string() +   "del_diag: FAILED"  + "\n";
};

void assign_functions_list::test_del_int_int()
{
    Real out = 0.;
    {
        test_function_del_int_int tf(1,1);
        out += m_tests.make(&tf,m_options);
    };
    {
        test_function_del_int_int tf(1,3);
        out += m_tests.make(&tf,m_options);
    };
    {
        test_function_del_int_int tf(3,1);
        out += m_tests.make(&tf,m_options);
    };
    {
        test_function_del_int_int tf(3,3);
        out += m_tests.make(&tf,m_options);
    };

    if (out == 0.)
        matcl::out_stream << std::string() +   "del_int_int: OK" + "\n";
    else
        matcl::out_stream << std::string() +   "del_int_int: FAILED"  + "\n";
};

void assign_functions_list::test_del_colon()
{
    Real out = 0.;

    test_function_del_colon tf;
    out += m_tests.make(&tf,m_options);

    if (out == 0.)
        matcl::out_stream << std::string() +   "del_colon: OK" + "\n";
    else
        matcl::out_stream << std::string() +   "del_colon: FAILED"  + "\n";
};

void assign_functions_list::test_del_colon_colon()
{
    Real out = 0.;

    test_function_del_colon_colon tf;
    out += m_tests.make(&tf,m_options);

    if (out == 0.)
        matcl::out_stream << std::string() +   "del_colon_colon: OK" + "\n";   
    else
        matcl::out_stream << std::string() +   "del_colon_colon: FAILED"  + "\n";
};

void assign_functions_list::test_set_int_scal()
{
    Real out = 0.;
    {
        test_function_set_int_scal tf(1);
        out += m_tests.make(&tf,m_options);
    };
    {
        test_function_set_int_scal tf(5);
        out += m_tests.make(&tf,m_options);
    };

    if (out == 0.)
        matcl::out_stream << std::string() +   "set_int_scal: OK" + "\n";
    else
        matcl::out_stream << std::string() +   "set_int_scal: FAILED"  + "\n";
};

void assign_functions_list::test_set_int2_scal()
{
    Real out = 0.;
    {
        test_function_set_int2_scal tf(1,1);
        out += m_tests.make(&tf,m_options);
    };
    {
        test_function_set_int2_scal tf(1,3);
        out += m_tests.make(&tf,m_options);
    };
    {
        test_function_set_int2_scal tf(3,1);
        out += m_tests.make(&tf,m_options);
    };
    {
        test_function_set_int2_scal tf(3,3);
        out += m_tests.make(&tf,m_options);
    };

    if (out == 0.)
        matcl::out_stream << std::string() +   "set_int2_scal: OK" + "\n";
    else
        matcl::out_stream << std::string() +   "set_int2_scal: FAILED"  + "\n";
};

void assign_functions_list::test_set_col_scal()
{
    Real out = 0.;

    test_function_set_col_scal tf;
    out += m_tests.make(&tf,m_options);

    if (out == 0.)
        matcl::out_stream << std::string() +   "set_col_scal: OK" + "\n";
    else
        matcl::out_stream << std::string() +   "set_col_scal: FAILED"  + "\n";
};

void assign_functions_list::test_set_col2_scal()
{
    Real out = 0.;

    test_function_set_col2_scal tf;
    out += m_tests.make(&tf,m_options);

    if (out == 0.)
        matcl::out_stream << std::string() +   "set_col2_scal: OK" + "\n";
    else
        matcl::out_stream << std::string() +   "set_col2_scal: FAILED"  + "\n";
};

void assign_functions_list::test_set_int_mat()
{
    Real out = 0.;
    {
        test_function_set_int_mat tf(1);
        out += m_tests.make(&tf,m_options);
    };
    {
        test_function_set_int_mat tf(5);
        out += m_tests.make(&tf,m_options);
    };

    if (out == 0.)
        matcl::out_stream << std::string() +   "set_int_mat: OK" + "\n";
    else
        matcl::out_stream << std::string() +   "set_int_mat: FAILED"  + "\n";
};

void assign_functions_list::test_set_int2_mat()
{
    Real out = 0.;
    {
        test_function_set_int2_mat tf(1,1);
        out += m_tests.make(&tf,m_options);
    };
    {
        test_function_set_int2_mat tf(1,3);
        out += m_tests.make(&tf,m_options);
    };
    {
        test_function_set_int2_mat tf(3,1);
        out += m_tests.make(&tf,m_options);
    };
    {
        test_function_set_int2_mat tf(3,3);
        out += m_tests.make(&tf,m_options);
    };

    if (out == 0.)
        matcl::out_stream << std::string() +   "set_int2_mat: OK" + "\n";
    else
        matcl::out_stream << std::string() +   "set_int2_mat: FAILED"  + "\n";
};

void assign_functions_list::test_set_col_mat()
{
    Real out = 0.;

    test_function_set_col_mat tf(ms);
    out += m_tests.make(&tf,m_options);

    if (out == 0.)
        matcl::out_stream << std::string() +   "set_col_mat: OK" + "\n";
    else
        matcl::out_stream << std::string() +   "set_col_mat: FAILED"  + "\n";
};

void assign_functions_list::test_set_diag_scal()
{
    Real out = 0.;

    {
        test_function_set_diag_scal tf(-1);
        out += m_tests.make(&tf,m_options);
    }

    {
        test_function_set_diag_scal tf(0);
        out += m_tests.make(&tf,m_options);
    }

    {
        test_function_set_diag_scal tf(1);
        out += m_tests.make(&tf,m_options);
    }

    if (out == 0.)
        matcl::out_stream << std::string() +   "set_diag_scal: OK" + "\n";
    else
        matcl::out_stream << std::string() +   "set_diag_scal: FAILED"  + "\n";
};

void assign_functions_list::test_set_col2_mat()
{
    Real out = 0.;

    test_function_set_col2_mat tf(ms);
    out += m_tests.make(&tf,m_options);

    if (out == 0.)
        matcl::out_stream << std::string() +   "set_col2_mat: OK" + "\n";
    else
        matcl::out_stream << std::string() +   "set_col2_mat: FAILED"  + "\n";
};

void assign_functions_list::test_set_diag_mat()
{
    Real out = 0.;

    {
        test_function_set_diag_mat tf(-1,ms);
        out += m_tests.make(&tf,m_options);
    }

    {
        test_function_set_diag_mat tf(0,ms);
        out += m_tests.make(&tf,m_options);
    }

    {
        test_function_set_diag_mat tf(1,ms);
        out += m_tests.make(&tf,m_options);
    }

    if (out == 0.)
        matcl::out_stream << std::string() +   "set_diag_mat: OK" + "\n";
    else
        matcl::out_stream << std::string() +   "set_diag_mat: FAILED"  + "\n";
};

Real test_function_delcols::eval_mat(const Matrix& mat,bool,int code )
{
    (void)code;

    colon_set ts(num_colons(), mat.cols(), code);
    Integer size    = ts.size();
    Real dif        = 0;

    for (Integer i = 0; i < size; ++i)
    {
        colon c     = ts.get(i);
        dif         += eval_mat_col(mat, c);
    }

    return dif;
};

Real test_function_delcols::eval_mat_col(const Matrix& mat, const colon& c1)
{
    try
    {		
        Matrix out_full;
        {
            Matrix mf = full(mat);
            out_full = mf.delcols(c1);	
        }

        Matrix out		= mat.delcols(c1);
        check_struct(out);

        Matrix c2		= irange(1,mat.cols());
        c2				= c2.delcols(c1);
        Real dif		= norm_1(out - out_full);

        Matrix out2		= mat(colon(),c2);
        check_struct(out2);

        dif				+=norm_1(out - out2);

        return dif;
    }
    catch(const std::exception& ex)
    {
        bool valid_1    = is_colon_valid(c1, mat.cols());

        if (valid_1 == false)
            return 0.0;

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

Real test_function_delcols::eval_scalar(const Scalar& ,bool,int code )
{
    (void)code;
    return 0.;
};

Real test_function_delrows::eval_mat(const Matrix& mat,bool,int code )
{
    (void)code;

    colon_set ts(num_colons(), mat.rows(), code);
    Integer size    = ts.size();
    Real dif        = 0;

    for (Integer i = 0; i < size; ++i)
    {
        colon c     = ts.get(i);
        dif         += eval_mat_col(mat, c);
    }

    return dif;
};

Real test_function_delrows::eval_mat_col(const Matrix& mat, const colon& c1)
{
    try
    {		
        Matrix out_full;
        {
            Matrix mf   = full(mat);
            out_full    = mf.delrows(c1);	
        }

        Matrix out		= mat.delrows(c1);
        check_struct(out);

        Real dif		= norm_1(out - out_full);
        Matrix c2		= irange(1,mat.rows());
        c2				= c2.delcols(c1);

        Matrix out2		= mat(c2,colon());
        check_struct(out2);

        dif				+=norm_1(out - out2);
        return dif;
    }
    catch(const std::exception& ex)
    {
        bool valid_1    = is_colon_valid(c1, mat.rows());

        if (valid_1 == false)
            return 0.0;

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

Real test_function_delrows::eval_scalar(const Scalar& ,bool,int code )
{
    (void)code;
    return 0.;
};

Real test_function_delrowscols::eval_mat(const Matrix& mat,bool,int code )
{
    (void)code;

    colon_set ts1(num_colons(), mat.rows(), code);
    colon_set ts2(num_colons(), mat.cols(), code+1);
    Integer size1   = ts1.size();
    Integer size2   = ts2.size();
    Real dif        = 0;

    for (Integer i = 0; i < size1; ++i)
    for (Integer j = 0; j < size2; ++j)
    {
        colon c1    = ts1.get(i);
        colon c2    = ts2.get(i);
        dif         += eval_mat_col(mat, c1, c2);
    }

    return dif;
};

Real test_function_delrowscols::eval_mat_col(const Matrix& mat, const colon& c1, const colon& c2)
{
    try
    {		
        Matrix out_full;
        {
            Matrix mf   = full(mat);
            out_full    = mf.delrowscols(c1,c2);	
        }

        Matrix out		= mat.delrowscols(c1,c2);
        check_struct(out);

        Real dif		= norm_1(out - out_full);

        Matrix out2     = mat.delrows(c1).delcols(c2);
        dif				+=norm_1(out - out2);

        out_full        = Matrix();
        out             = Matrix();
        out2            = Matrix();
        return dif;
    }
    catch(const std::exception& ex)
    {
        bool valid_1    = is_colon_valid(c1, mat.rows());
        bool valid_2    = is_colon_valid(c2, mat.cols());

        if (valid_1 == false || valid_2 == false)
            return 0.0;

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

Real test_function_delrowscols::eval_scalar(const Scalar& ,bool,int code )
{
    (void)code;
    return 0.;
};

Real test_function_get_colon::eval_mat(const Matrix& mat,bool,int code )
{
    (void)code;

    colon_set ts(num_colons(), mat.rows() * mat.cols(), code);
    Integer size    = ts.size();
    Real dif        = 0;

    for (Integer i = 0; i < size; ++i)
    {
        colon c     = ts.get(i);
        dif         += eval_mat_col(mat, c);
    }

    return dif;
};

Real test_function_get_colon::eval_mat_col(const Matrix& mat, const colon& c1)
{
    Matrix out_full;
    try
    {
        out_full = full(mat)(c1);	
    }
    catch(error::invalid_single_index&)
    {
        bool is_valid   = is_colon_valid(c1, mat.rows() * mat.cols());
        return is_valid ? 1.0 : 0.0;
    }

    try
    {		
        Matrix out		= mat(c1);
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
Real test_function_get_colon::eval_scalar(const Scalar& ,bool,int code )
{
    (void)code;
    return 0.;
};

Real test_function_get_colon_colon::eval_mat(const Matrix& mat,bool,int code )
{
    (void)code;

    colon_set ts1(num_colons(), mat.rows(), code);
    colon_set ts2(num_colons(), mat.cols(), code + 1);

    Integer size1   = ts1.size();
    Integer size2   = ts2.size();

    Real dif        = 0;

    for (Integer i = 0; i < size1; ++i)
    for (Integer j = 0; j < size2; ++j)
    {
        colon c1    = ts1.get(i);
        colon c2    = ts2.get(j);
        dif         += eval_mat_col(mat, c1, c2);
    }

    return dif;
};

Real test_function_get_colon_colon::eval_mat_col(const Matrix& mat, const colon& c1, const colon& c2)
{
    try
    {		
        Matrix out_full;
        {
            out_full = full(mat)(c1,c2);	
        }

        Matrix out		= mat(c1,c2);
        check_struct(out);

        Real dif		= norm_1(out - out_full);

        return dif;
    }
    catch(const std::exception& ex)
    {
        bool valid_1    = is_colon_valid(c1, mat.rows());
        bool valid_2    = is_colon_valid(c2, mat.cols());

        if (valid_1 == false || valid_2 == false)
            return 0.0;

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
Real test_function_get_colon_colon::eval_scalar(const Scalar& ,bool,int code )
{
    (void)code;
    return 0.;
};

Real test_function_get_int::eval_mat(const Matrix& mat,bool,int code )
{
    (void)code;
    Matrix out_full;
    try
    {
        out_full = full(mat)(i);	
    }
    catch(...)
    {
        return 0.;
    };
    try
    {		
        Matrix out		= mat(i);
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
Real test_function_get_int::eval_scalar(const Scalar& ,bool,int code )
{
    (void)code;
    return 0.;
};

Real test_function_get_int_int::eval_mat(const Matrix& mat,bool,int code )
{
    (void)code;
    Matrix out_full;
    try
    {
        out_full = full(mat)(i,j);	
    }
    catch(...)
    {
        return 0.;
    };
    try
    {		
        Matrix out		= mat(i,j);
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
Real test_function_get_int_int::eval_scalar(const Scalar& ,bool,int code )
{
    (void)code;
    return 0.;
};

Real test_function_del_colon::eval_mat(const Matrix& mat,bool,int code )
{
    (void)code;

    colon_set ts(num_colons(), mat.rows() * mat.cols(), code);
    Integer size    = ts.size();
    Real dif        = 0;

    for (Integer i = 0; i < size; ++i)
    {
        colon c     = ts.get(i);
        dif         += eval_mat_col(mat, c);
    }

    return dif;
};
Real test_function_del_colon::eval_mat_col(const Matrix& mat, const colon& c1)
{
    try
    {		
        Matrix out_full;
        {
            out_full    = full(mat);	
            out_full(c1)= zeros(0,0);
        }

        Matrix out		= mat;
        out(c1)			= zeros(0,0);
        check_struct(out);

        Real dif		= norm_1(out - out_full);

        return dif;
    }
    catch(const std::exception& ex)
    {
        bool valid_1    = is_colon_valid(c1, mat.rows() * mat.cols());

        if (valid_1 == false)
            return 0.0;

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
Real test_function_del_colon::eval_scalar(const Scalar& ,bool,int code )
{
    (void)code;
    return 0.;
};

Real test_function_del_colon_colon::eval_mat(const Matrix& mat,bool,int code )
{
    (void)code;

    colon_set ts1(num_colons(), mat.rows(), code);
    colon_set ts2(num_colons(), mat.cols(), code + 1);

    Integer size1   = ts1.size();
    Integer size2   = ts2.size();

    Real dif        = 0;

    for (Integer i = 0; i < size1; ++i)
    for (Integer j = 0; j < size2; ++j)
    {
        colon c1    = ts1.get(i);
        colon c2    = ts2.get(j);
        dif         += eval_mat_col(mat, c1, c2);
    }

    return dif;
};
Real test_function_del_colon_colon::eval_mat_col(const Matrix& mat, const colon& c1, const colon& c2)
{
    try
    {		
        Matrix out_full;
        {
            out_full        = full(mat);	
            out_full(c1,c2) = izeros(0,0);
        }

        Matrix out		= mat;
        out(c1,c2)		= izeros(0,0);
        check_struct(out);

        Real dif		= norm_1(out - out_full);

        return dif;
    }
    catch(const std::exception& ex)
    {
        bool valid_1    = is_colon_valid(c1, mat.rows());
        bool valid_2    = is_colon_valid(c2, mat.cols());

        if (valid_1 == false || valid_2 == false)
            return 0.0;

        bool is_all_1   = is_colon_all(c1);
        bool is_all_2   = is_colon_all(c2);
        
        if (is_all_1 == false && is_all_2 == false)
            return 0.0;

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

Real test_function_del_colon_colon::eval_scalar(const Scalar& ,bool,int code )
{
    (void)code;
    return 0.;
};

Real test_function_del_int::eval_mat(const Matrix& mat,bool,int code )
{
    (void)code;
    Matrix out_full;
    try
    {
        out_full    = full(mat);	
        out_full(i) = zeros(0,0);
    }
    catch(...)
    {
        return 0.;
    };

    try
    {		
        Matrix out		= mat;
        out(i)			= zeros(0,0);
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
Real test_function_del_int::eval_scalar(const Scalar& ,bool,int code )
{
    (void)code;
    return 0.;
};

Real test_function_del_diag::eval_mat(const Matrix& mat,bool,int code )
{
    (void)code;
    (void)mat;

    //mat.diag(d) = [] is interpreted as assignment not as deletion of a diagonal
    return 0.0;
};
Real test_function_del_diag::eval_scalar(const Scalar& ,bool,int code )
{
    (void)code;
    return 0.;
};

Real test_function_del_int_int::eval_mat(const Matrix& mat,bool,int code )
{
    (void)code;
    Matrix out_full;
    try
    {
        out_full        = full(mat);	
        out_full(i,j)   = zeros(0,0);
    }
    catch(...)
    {
        return 0.;
    };

    try
    {		
        Matrix out		= mat;
        out(i,j)		= zeros(0,0);
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
Real test_function_del_int_int::eval_scalar(const Scalar& ,bool,int code )
{
    (void)code;
    return 0.;
};

Real test_function_drop_colon::eval_mat(const Matrix& mat,bool,int code )
{
    (void)code;

    colon_set ts(num_colons(), mat.rows() * mat.cols(), code);
    Integer size    = ts.size();
    Real dif        = 0;

    for (Integer i = 0; i < size; ++i)
    {
        colon c     = ts.get(i);
        dif         += eval_mat_col(mat, c);
    }

    return dif;
};
Real test_function_drop_colon::eval_mat_col(const Matrix& mat, const colon& c1)
{
    try
    {		
        Matrix mat2 = mat;
        {
            Matrix mf   = full(mat);
            mf(c1)      = 0;		
            mat2(c1)    = 1;
        }

        Matrix mf1_0 = full(mat);
        Matrix mf2_0 = full(mat);
        Matrix mf3_0 = full(mat);
        Matrix mf4_0 = full(mat);
        Matrix mf5_0 = full(mat);

        mf1_0(c1) = Integer();
        mf2_0(c1) = Real();
        mf3_0(c1) = Complex();
        mf4_0(c1) = Float();
        mf5_0(c1) = Float_complex();

        Matrix out1_0	= mat;
        Matrix out2_0	= mat;
        Matrix out3_0	= mat;
        Matrix out4_0	= mat;
        Matrix out5_0	= mat;

        Matrix out1_1	= mat;
        Matrix out2_1	= mat;
        Matrix out3_1	= mat;
        Matrix out4_1	= mat;
        Matrix out5_1	= mat;

        out1_0(c1) = Integer();
        out2_0(c1) = Real();
        out3_0(c1) = Complex();
        out4_0(c1) = Float();
        out5_0(c1) = Float_complex();

        out1_1(c1) = full(Integer());
        out2_1(c1) = full(Real());
        out3_1(c1) = full(Complex());        
        out4_1(c1) = full(Float());        
        out5_1(c1) = full(Float_complex());        

        check_struct(out1_0);
        check_struct(out2_0);
        check_struct(out3_0);
        check_struct(out4_0);
        check_struct(out5_0);

        check_struct(out1_1);
        check_struct(out2_1);
        check_struct(out3_1);
        check_struct(out4_1);
        check_struct(out5_1);

        Real dif		= norm_1(out1_0 - mf1_0);
        dif				+= norm_1(out2_0 - mf2_0);
        dif				+= norm_1(out3_0 - mf3_0);
        dif				+= norm_1(out4_0 - mf4_0);
        dif				+= norm_1(out5_0 - mf5_0);

        dif				+= norm_1(out1_0(c1));
        dif				+= norm_1(out2_0(c1));
        dif				+= norm_1(out3_0(c1));
        dif				+= norm_1(out4_0(c1));
        dif				+= norm_1(out5_0(c1));

        dif				+= norm_1(out1_1(c1));
        dif				+= norm_1(out2_1(c1));
        dif				+= norm_1(out3_1(c1));
        dif				+= norm_1(out4_1(c1));
        dif				+= norm_1(out5_1(c1));

        dif				+= norm_1(mat2(c1)-1);
        return dif;
    }
    catch(const std::exception& ex)
    {
        bool valid_1    = is_colon_valid(c1, mat.rows() * mat.cols());

        if (valid_1 == false)
            return 0.0;

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
Real test_function_drop_colon::eval_scalar(const Scalar& ,bool,int code )
{
    (void)code;
    return 0.;
};

Real test_function_drop_colon_colon::eval_mat(const Matrix& mat,bool,int code )
{
    (void)code;

    colon_set ts1(num_colons(), mat.rows(), code);
    colon_set ts2(num_colons(), mat.cols(), code + 1);

    Integer size1   = ts1.size();
    Integer size2   = ts2.size();

    Real dif        = 0;

    for (Integer i = 0; i < size1; ++i)
    for (Integer j = 0; j < size2; ++j)
    {
        colon c1    = ts1.get(i);
        colon c2    = ts2.get(j);
        dif         += eval_mat_col(mat, c1, c2);
    }

    return dif;
};
Real test_function_drop_colon_colon::eval_mat_col(const Matrix& mat, const colon& c1, const colon& c2)
{
    try
    {		
        Matrix mat2 = mat;
        Matrix mat3 = mat;
        {
            Matrix mf   = full(mat);
            mf(c1,c2)   = 0;		
            mat2(c1,c2) = 1;
            mat3(c1,c2) = 1.f;
        }

        Matrix mf1_0 = full(mat);
        Matrix mf2_0 = full(mat);
        Matrix mf3_0 = full(mat);
        Matrix mf4_0 = full(mat);
        Matrix mf5_0 = full(mat);

        mf1_0(c1,c2) = Integer();
        mf2_0(c1,c2) = Real();
        mf3_0(c1,c2) = Complex();
        mf4_0(c1,c2) = Float();
        mf5_0(c1,c2) = Float_complex();

        Matrix out1_0	= mat;
        Matrix out2_0	= mat;
        Matrix out3_0	= mat;
        Matrix out4_0	= mat;
        Matrix out5_0	= mat;

        Matrix out1_1	= mat;
        Matrix out2_1	= mat;
        Matrix out3_1	= mat;
        Matrix out4_1	= mat;
        Matrix out5_1	= mat;

        out1_0(c1,c2) = Integer();
        out2_0(c1,c2) = Real();
        out3_0(c1,c2) = Complex();
        out4_0(c1,c2) = Float();
        out5_0(c1,c2) = Float_complex();

        out1_1(c1,c2) = full(Integer());
        out2_1(c1,c2) = full(Real());
        out3_1(c1,c2) = full(Complex());
        out4_1(c1,c2) = full(Float());
        out5_1(c1,c2) = full(Float_complex());

        check_struct(out1_0);
        check_struct(out2_0);
        check_struct(out3_0);
        check_struct(out4_0);
        check_struct(out5_0);

        check_struct(out1_1);
        check_struct(out2_1);
        check_struct(out3_1);
        check_struct(out4_1);
        check_struct(out5_1);

        check_struct(mat2);
        check_struct(mat3);

        Real dif		= norm_1(out1_0 - mf1_0);
        dif				+= norm_1(out2_0 - mf2_0);
        dif				+= norm_1(out3_0 - mf3_0);
        dif				+= norm_1(out4_0 - mf4_0);
        dif				+= norm_1(out5_0 - mf5_0);

        dif				+= norm_1(out1_0(c1,c2));
        dif				+= norm_1(out2_0(c1,c2));
        dif				+= norm_1(out3_0(c1,c2));
        dif				+= norm_1(out4_0(c1,c2));
        dif				+= norm_1(out5_0(c1,c2));

        dif				+= norm_1(out1_1(c1,c2));
        dif				+= norm_1(out2_1(c1,c2));
        dif				+= norm_1(out3_1(c1,c2));
        dif				+= norm_1(out4_1(c1,c2));
        dif				+= norm_1(out5_1(c1,c2));

        dif				+= norm_1(mat2(c1,c2)-1);
        dif				+= norm_1(mat3(c1,c2)-1.f);

        return dif;
    }
    catch(const std::exception& ex)
    {
        bool valid_1    = is_colon_valid(c1, mat.rows());
        bool valid_2    = is_colon_valid(c2, mat.cols());

        if (valid_1 == false || valid_2 == false)
            return 0.0;

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
Real test_function_drop_colon_colon::eval_scalar(const Scalar& ,bool,int code )
{
    (void)code;
    return 0.;
};

Real test_function_drop_int::eval_mat(const Matrix& mat,bool,int code )
{
    (void)code;
    Matrix mat2     = mat;
    Matrix mat3     = mat;
    try
    {
        Matrix mf   = full(mat);
        mf(i)       = 0;		
        mat2(i)     = 1;
        mat3(i)     = 1.f;
    }
    catch(...)
    {
        return 0.;
    };

    try
    {		
        Matrix mf1_0    = full(mat);
        Matrix mf2_0    = full(mat);
        Matrix mf3_0    = full(mat);
        Matrix mf4_0    = full(mat);
        Matrix mf5_0    = full(mat);

        mf1_0(i)        = Integer();
        mf2_0(i)        = Real();
        mf3_0(i)        = Complex();
        mf4_0(i)        = Float();
        mf5_0(i)        = Float_complex();

        Matrix out1_0	= mat;
        Matrix out2_0	= mat;
        Matrix out3_0	= mat;
        Matrix out4_0	= mat;
        Matrix out5_0	= mat;

        Matrix out1_1	= mat;
        Matrix out2_1	= mat;
        Matrix out3_1	= mat;
        Matrix out4_1	= mat;
        Matrix out5_1	= mat;

        out1_0(i)       = Integer();
        out2_0(i)       = Real();
        out3_0(i)       = Complex();
        out4_0(i)       = Float();
        out5_0(i)       = Float_complex();

        out1_1(i)       = full(Integer());
        out2_1(i)       = full(Real());
        out3_1(i)       = full(Complex());
        out4_1(i)       = full(Float());
        out5_1(i)       = full(Float_complex());

        check_struct(out1_0);
        check_struct(out2_0);
        check_struct(out3_0);
        check_struct(out4_0);
        check_struct(out5_0);

        check_struct(out1_1);
        check_struct(out2_1);
        check_struct(out3_1);
        check_struct(out4_1);
        check_struct(out5_1);

        check_struct(mat2);
        check_struct(mat3);

        Real dif		= norm_1(out1_0 - mf1_0);
        dif				+= norm_1(out2_0 - mf2_0);
        dif				+= norm_1(out3_0 - mf3_0);
        dif				+= norm_1(out4_0 - mf4_0);
        dif				+= norm_1(out5_0 - mf5_0);

        dif				+= norm_1(out1_0(i));
        dif				+= norm_1(out2_0(i));
        dif				+= norm_1(out3_0(i));
        dif				+= norm_1(out4_0(i));
        dif				+= norm_1(out5_0(i));

        dif				+= norm_1(out1_1(i));
        dif				+= norm_1(out2_1(i));
        dif				+= norm_1(out3_1(i));
        dif				+= norm_1(out4_1(i));
        dif				+= norm_1(out5_1(i));

        dif				+= norm_1(mat2(i)-1);
        dif				+= norm_1(mat3(i)-1.f);

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
Real test_function_drop_int::eval_scalar(const Scalar& ,bool,int code )
{
    (void)code;
    return 0.;
};

Real test_function_drop_diag::eval_mat(const Matrix& mat,bool,int code )
{
    (void)code;

    try
    {		
        Matrix mat2     = mat;
        Matrix mat3     = mat;
        Matrix mf       = full(mat);
        mf.diag(d)      = 0;		
        mat2.diag(d)    = 1;
        mat3.diag(d)    = 1.f;

        Matrix mf1_0    = full(mat);
        Matrix mf2_0    = full(mat);
        Matrix mf3_0    = full(mat);
        Matrix mf4_0    = full(mat);
        Matrix mf5_0    = full(mat);

        mf1_0.diag(d)   = Integer();
        mf2_0.diag(d)   = Real();
        mf3_0.diag(d)   = Complex();
        mf4_0.diag(d)   = Float();
        mf5_0.diag(d)   = Float_complex();

        Matrix out1_0	= mat;
        Matrix out2_0	= mat;
        Matrix out3_0	= mat;
        Matrix out4_0	= mat;
        Matrix out5_0	= mat;

        Matrix out1_1	= mat;
        Matrix out2_1	= mat;
        Matrix out3_1	= mat;
        Matrix out4_1	= mat;
        Matrix out5_1	= mat;

        out1_0.diag(d)  = Integer();
        out2_0.diag(d)  = Real();
        out3_0.diag(d)  = Complex();
        out4_0.diag(d)  = Float();
        out5_0.diag(d)  = Float_complex();

        out1_1.diag(d)  = full(Integer());
        out2_1.diag(d)  = full(Real());
        out3_1.diag(d)  = full(Complex());
        out4_1.diag(d)  = full(Float());
        out5_1.diag(d)  = full(Float_complex());

        check_struct(out1_0);
        check_struct(out2_0);
        check_struct(out3_0);
        check_struct(out4_0);
        check_struct(out5_0);

        check_struct(out1_1);
        check_struct(out2_1);
        check_struct(out3_1);
        check_struct(out4_1);
        check_struct(out5_1);

        check_struct(mat2);
        check_struct(mat3);

        Real dif		= norm_1(out1_0 - mf1_0);
        dif				+= norm_1(out2_0 - mf2_0);
        dif				+= norm_1(out3_0 - mf3_0);
        dif				+= norm_1(out4_0 - mf4_0);
        dif				+= norm_1(out5_0 - mf5_0);

        dif				+= norm_1(out1_0.diag(d));
        dif				+= norm_1(out2_0.diag(d));
        dif				+= norm_1(out3_0.diag(d));
        dif				+= norm_1(out4_0.diag(d));
        dif				+= norm_1(out5_0.diag(d));

        dif				+= norm_1(out1_1.diag(d));
        dif				+= norm_1(out2_1.diag(d));
        dif				+= norm_1(out3_1.diag(d));
        dif				+= norm_1(out4_1.diag(d));
        dif				+= norm_1(out5_1.diag(d));

        dif				+= norm_1(mat2.diag(d)-1);
        dif				+= norm_1(mat3.diag(d)-1.f);

        return dif;
    }
    catch(const std::exception& ex)
    {
        bool valid  = is_valid_diag(d, mat);

        if (valid == false)
            return 0.0;

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
Real test_function_drop_diag::eval_scalar(const Scalar& ,bool,int code )
{
    (void)code;
    return 0.;
};


Real test_function_drop_int_int::eval_mat(const Matrix& mat,bool,int code )
{
    (void)code;
    Matrix mat2 = mat;
    Matrix mat3 = mat;
    try
    {
        Matrix mf   = full(mat);
        mf(i,j)     = 0;		
        mat2(i,j)   = 1;
        mat3(i,j)   = 1.f;
    }
    catch(...)
    {
        return 0.;
    };

    try
    {		
        Matrix mf1_0 = full(mat);
        Matrix mf2_0 = full(mat);
        Matrix mf3_0 = full(mat);
        Matrix mf4_0 = full(mat);
        Matrix mf5_0 = full(mat);

        mf1_0(i,j)  = Integer();
        mf2_0(i,j)  = Real();
        mf3_0(i,j)  = Complex();
        mf4_0(i,j)  = Float();
        mf5_0(i,j)  = Float_complex();

        Matrix out1_0	= mat;
        Matrix out2_0	= mat;
        Matrix out3_0	= mat;
        Matrix out4_0	= mat;
        Matrix out5_0	= mat;

        Matrix out1_1	= mat;
        Matrix out2_1	= mat;
        Matrix out3_1	= mat;
        Matrix out4_1	= mat;
        Matrix out5_1	= mat;

        out1_0(i,j)     = Integer();
        out2_0(i,j)     = Real();
        out3_0(i,j)     = Complex();
        out4_0(i,j)     = Float();
        out5_0(i,j)     = Float_complex();

        out1_1(i,j)     = full(Integer());
        out2_1(i,j)     = full(Real());
        out3_1(i,j)     = full(Complex());
        out4_1(i,j)     = full(Float());
        out5_1(i,j)     = full(Float_complex());

        check_struct(out1_0);
        check_struct(out2_0);
        check_struct(out3_0);
        check_struct(out4_0);
        check_struct(out5_0);

        check_struct(out1_1);
        check_struct(out2_1);
        check_struct(out3_1);
        check_struct(out4_1);
        check_struct(out5_1);

        check_struct(mat2);
        check_struct(mat3);

        Real dif		= norm_1(out1_0 - mf1_0);
        dif				+= norm_1(out2_0 - mf2_0);
        dif				+= norm_1(out3_0 - mf3_0);
        dif				+= norm_1(out4_0 - mf4_0);
        dif				+= norm_1(out5_0 - mf5_0);

        dif				+= norm_1(out1_0(i,j));
        dif				+= norm_1(out2_0(i,j));
        dif				+= norm_1(out3_0(i,j));
        dif				+= norm_1(out4_0(i,j));
        dif				+= norm_1(out5_0(i,j));

        dif				+= norm_1(out1_1(i,j));
        dif				+= norm_1(out2_1(i,j));
        dif				+= norm_1(out3_1(i,j));
        dif				+= norm_1(out4_1(i,j));
        dif				+= norm_1(out5_1(i,j));

        dif				+= norm_1(mat2(i,j)-1);
        dif				+= norm_1(mat3(i,j)-1.f);

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
Real test_function_drop_int_int::eval_scalar(const Scalar& ,bool,int code )
{
    (void)code;
    return 0.;
};

Real test_function_set_int_scal::eval_mat(const Matrix& mat,bool,int code )
{
    (void)code;
    try
    {
        Matrix mf   = full(mat);
        mf(i)       = 1;
    }
    catch(...)
    {
        return 0.;
    };

    try
    {		
        Matrix mf1_0 = full(mat);
        Matrix mf2_0 = full(mat);
        Matrix mf3_0 = full(mat);
        Matrix mf4_0 = full(mat);
        Matrix mf5_0 = full(mat);

        Matrix mf1  = full(mat);
        Matrix mf2  = full(mat);
        Matrix mf3  = full(mat);
        Matrix mf4  = full(mat);
        Matrix mf5  = full(mat);

        Matrix mf_self = full(mat);

        mf1_0(i)    = Integer();
        mf2_0(i)    = Real();
        mf3_0(i)    = Complex();
        mf4_0(i)    = Float();
        mf5_0(i)    = Float_complex();

        mf1(i)      = Integer(1);
        mf2(i)      = Real(1);
        mf3(i)      = Complex(1);
        mf4(i)      = Float(1);
        mf5(i)      = Float_complex(1);

        mf_self(i)  = mf_self(i);

        Matrix out1_0	= mat;
        Matrix out2_0	= mat;
        Matrix out3_0	= mat;
        Matrix out4_0	= mat;
        Matrix out5_0	= mat;

        Matrix out1		= mat;
        Matrix out2		= mat;
        Matrix out3		= mat;
        Matrix out4		= mat;
        Matrix out5		= mat;

        Matrix out_self = mat.clone();

        out1_0(i)		= Integer();
        out2_0(i)		= Real();
        out3_0(i)		= Complex();
        out4_0(i)		= Float();
        out5_0(i)		= Float_complex();

        out1(i)			= Integer(1);
        out2(i)			= Real(1);
        out3(i)			= Complex(1);
        out4(i)			= Float(1);
        out5(i)			= Float_complex(1);

        out_self(i)     = out_self(i);

        check_struct(out1_0);
        check_struct(out2_0);
        check_struct(out3_0);
        check_struct(out4_0);
        check_struct(out5_0);

        check_struct(out1);
        check_struct(out2);
        check_struct(out3);
        check_struct(out4);
        check_struct(out5);

        check_struct(out_self);

        Real dif		= norm_1(out1_0 - mf1_0);
        dif				+= norm_1(out2_0 - mf2_0);
        dif				+= norm_1(out3_0 - mf3_0);
        dif				+= norm_1(out4_0 - mf4_0);
        dif				+= norm_1(out5_0 - mf5_0);

        dif				+= norm_1(out1 - mf1);
        dif				+= norm_1(out2 - mf2);
        dif				+= norm_1(out3 - mf3);
        dif				+= norm_1(out4 - mf4);
        dif				+= norm_1(out5 - mf5);

        dif				+= norm_1(out1_0(i));
        dif				+= norm_1(out2_0(i));
        dif				+= norm_1(out3_0(i));
        dif				+= norm_1(out4_0(i));
        dif				+= norm_1(out5_0(i));

        dif				+= norm_1(out1(i)-1);
        dif				+= norm_1(out2(i)-1);
        dif				+= norm_1(out3(i)-1);
        dif				+= norm_1(out4(i)-1);
        dif				+= norm_1(out5(i)-1);

        dif				+= norm_1(out_self - mf_self);
        dif				+= norm_1(out_self - mat);

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
Real test_function_set_int_scal::eval_scalar(const Scalar& ,bool,int code )
{
    (void)code;
    return 0.;
};

Real test_function_set_int2_scal::eval_mat(const Matrix& mat,bool,int code )
{
    (void)code;
    try
    {
        Matrix mf   = full(mat);
        mf(i,j)     = 1;
    }
    catch(...)
    {
        return 0.;
    };

    try
    {		
        Matrix mf1_0 = full(mat);
        Matrix mf2_0 = full(mat);
        Matrix mf3_0 = full(mat);
        Matrix mf4_0 = full(mat);
        Matrix mf5_0 = full(mat);

        Matrix mf1  = full(mat);
        Matrix mf2  = full(mat);
        Matrix mf3  = full(mat);
        Matrix mf4  = full(mat);
        Matrix mf5  = full(mat);

        Matrix mf_self = full(mat);

        mf1_0(i,j)  = Integer();
        mf2_0(i,j)  = Real();
        mf3_0(i,j)  = Complex();
        mf4_0(i,j)  = Float();
        mf5_0(i,j)  = Float_complex();

        mf1(i,j)    = Integer(1);
        mf2(i,j)    = Real(1);
        mf3(i,j)    = Complex(1);
        mf4(i,j)    = Float(1);
        mf5(i,j)    = Float_complex(1);

        mf_self(i,j) = mf_self(i,j);

        Matrix out1_0	= mat;
        Matrix out2_0	= mat;
        Matrix out3_0	= mat;
        Matrix out4_0	= mat;
        Matrix out5_0	= mat;

        Matrix out1		= mat;
        Matrix out2		= mat;
        Matrix out3		= mat;
        Matrix out4		= mat;
        Matrix out5		= mat;

        Matrix out_self = mat.clone();

        out1_0(i,j)		= Integer();
        out2_0(i,j)		= Real();
        out3_0(i,j)		= Complex();
        out4_0(i,j)		= Float();
        out5_0(i,j)		= Float_complex();

        out1(i,j)		= Integer(1);
        out2(i,j)		= Real(1);
        out3(i,j)		= Complex(1);
        out4(i,j)		= Float(1);
        out5(i,j)		= Float_complex(1);

        out_self(i,j)   = out_self(i,j);

        check_struct(out1_0);
        check_struct(out2_0);
        check_struct(out3_0);
        check_struct(out4_0);
        check_struct(out5_0);

        check_struct(out1);
        check_struct(out2);
        check_struct(out3);
        check_struct(out4);
        check_struct(out5);

        check_struct(out_self);

        Real dif		= norm_1(out1_0 - mf1_0);
        dif				+= norm_1(out2_0 - mf2_0);
        dif				+= norm_1(out3_0 - mf3_0);
        dif				+= norm_1(out4_0 - mf4_0);
        dif				+= norm_1(out5_0 - mf5_0);

        dif				+= norm_1(out1 - mf1);
        dif				+= norm_1(out2 - mf2);
        dif				+= norm_1(out3 - mf3);
        dif				+= norm_1(out4 - mf4);
        dif				+= norm_1(out5 - mf5);

        dif				+= norm_1(out1_0(i,j));
        dif				+= norm_1(out2_0(i,j));
        dif				+= norm_1(out3_0(i,j));
        dif				+= norm_1(out4_0(i,j));
        dif				+= norm_1(out5_0(i,j));

        dif				+= norm_1(out1(i,j)-1);
        dif				+= norm_1(out2(i,j)-1);
        dif				+= norm_1(out3(i,j)-1);
        dif				+= norm_1(out4(i,j)-1);
        dif				+= norm_1(out5(i,j)-1);
        
        dif				+= norm_1(out_self - mf_self);
        dif				+= norm_1(out_self - mat);

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
Real test_function_set_int2_scal::eval_scalar(const Scalar& ,bool,int code )
{
    (void)code;
    return 0.;
};

Real test_function_set_col_scal::eval_mat(const Matrix& mat,bool,int code )
{
    (void)code;

    colon_set ts(num_colons(), mat.rows() * mat.cols(), code);
    Integer size    = ts.size();
    Real dif        = 0;

    for (Integer i = 0; i < size; ++i)
    {
        colon c     = ts.get(i);
        dif         += eval_mat_col(mat, c);
    }

    return dif;
};
Real test_function_set_col_scal::eval_mat_col(const Matrix& mat, const colon& c1)
{
    try
    {		
        Matrix mf1		= full(mat);
        Matrix mf2		= full(mat);
        Matrix mf3		= full(mat);
        Matrix mf4		= full(mat);
        Matrix mf5		= full(mat);

        Matrix mf_self  = full(mat);

        mf1(c1)			= Integer(1);
        mf2(c1)			= Real(1);
        mf3(c1)			= Complex(1);
        mf4(c1)			= Float(1);
        mf5(c1)			= Float_complex(1);

        mf_self(c1)     = mf_self(c1);

        Matrix out1		= mat;
        Matrix out2		= mat;
        Matrix out3		= mat;
        Matrix out4		= mat;
        Matrix out5		= mat;

        Matrix out_self = mat.clone();

        out1(c1)		= Integer(1);
        out2(c1)		= Real(1);
        out3(c1)		= Complex(1);
        out4(c1)		= Float(1);
        out5(c1)		= Float_complex(1);

        out_self(c1)    = out_self(c1);

        check_struct(out1);
        check_struct(out2);
        check_struct(out3);
        check_struct(out4);
        check_struct(out5);

        check_struct(out_self);

        Real dif		= 0;
        dif				+= norm_1(out1 - mf1);
        dif				+= norm_1(out2 - mf2);
        dif				+= norm_1(out3 - mf3);
        dif				+= norm_1(out4 - mf4);
        dif				+= norm_1(out5 - mf5);

        bool empty      = out1(c1).rows() == 0 || out1(c1).cols() == 0;

        if (empty == false)
        {
            dif		    += norm_1(out1(c1)-1);
            dif		    += norm_1(out2(c1)-1);
            dif		    += norm_1(out3(c1)-1);
            dif		    += norm_1(out4(c1)-1);
            dif		    += norm_1(out5(c1)-1);
        };        

        dif				+= norm_1(out_self - mf_self);
        dif				+= norm_1(out_self - mat);

        return dif;
    }
    catch(const std::exception& ex)
    {
        bool valid_1    = is_colon_valid(c1, mat.rows() * mat.cols());

        if (valid_1 == false)
            return 0.0;

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
Real test_function_set_col_scal::eval_scalar(const Scalar& ,bool,int code )
{
    (void)code;
    return 0;
};

Real test_function_set_col2_scal::eval_mat(const Matrix& mat,bool,int code )
{
    (void)code;

    colon_set ts1(num_colons(), mat.rows(), code);
    colon_set ts2(num_colons(), mat.cols(), code + 1);

    Integer size1   = ts1.size();
    Integer size2   = ts2.size();

    Real dif        = 0;

    for (Integer i = 0; i < size1; ++i)
    for (Integer j = 0; j < size2; ++j)
    {
        colon c1    = ts1.get(i);
        colon c2    = ts2.get(j);
        dif         += eval_mat_col(mat, c1, c2);
    }

    return dif;
};
Real test_function_set_col2_scal::eval_mat_col(const Matrix& mat, const colon& c1, const colon& c2)
{
    try
    {		
        Matrix mf1		= full(mat);
        Matrix mf2		= full(mat);
        Matrix mf3		= full(mat);
        Matrix mf4		= full(mat);
        Matrix mf5		= full(mat);

        Matrix mf_self  = full(mat);

        mf1(c1,c2)		= Integer(1);
        mf2(c1,c2)		= Real(1);
        mf3(c1,c2)		= Complex(1);
        mf4(c1,c2)		= Float(1);
        mf5(c1,c2)		= Float_complex(1);

        mf_self(c1,c2)  = mf_self(c1,c2);

        Matrix out1		= mat;
        Matrix out2		= mat;
        Matrix out3		= mat;
        Matrix out4		= mat;
        Matrix out5		= mat;

        Matrix out_self = mat.clone();

        out1(c1,c2)		= Integer(1);
        out2(c1,c2)		= Real(1);
        out3(c1,c2)		= Complex(1);
        out4(c1,c2)		= Float(1);
        out5(c1,c2)		= Float_complex(1);

        out_self(c1,c2) = out_self(c1,c2);

        check_struct(out1);
        check_struct(out2);
        check_struct(out3);
        check_struct(out4);
        check_struct(out5);

        check_struct(out_self);

        Real dif		= 0;
        dif				+= norm_1(out1 - mf1);
        dif				+= norm_1(out2 - mf2);
        dif				+= norm_1(out3 - mf3);
        dif				+= norm_1(out4 - mf4);
        dif				+= norm_1(out5 - mf5);

        dif				+= norm_1(out1(c1,c2)-1);
        dif				+= norm_1(out2(c1,c2)-1);
        dif				+= norm_1(out3(c1,c2)-1);
        dif				+= norm_1(out4(c1,c2)-1);
        dif				+= norm_1(out5(c1,c2)-1);
        
        dif				+= norm_1(out_self - mf_self);
        dif				+= norm_1(out_self - mat);

        return dif;
    }
    catch(const std::exception& ex)
    {
        bool valid_1    = is_colon_valid(c1, mat.rows());
        bool valid_2    = is_colon_valid(c2, mat.cols());

        if (valid_1 == false || valid_2 == false)
            return 0.0;

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
Real test_function_set_col2_scal::eval_scalar(const Scalar& ,bool,int code )
{
    (void)code;
    return 0.;
};

Real test_function_set_diag_scal::eval_mat(const Matrix& mat,bool,int code )
{
    (void)code;

    try
    {		
        Matrix mf   = full(mat);
        mf.diag(d)  = 1;

        Matrix mf1		= full(mat);
        Matrix mf2		= full(mat);
        Matrix mf3		= full(mat);
        Matrix mf4		= full(mat);
        Matrix mf5		= full(mat);

        Matrix mf_self  = full(mat);

        mf1.diag(d)		= Integer(1);
        mf2.diag(d)		= Real(1);
        mf3.diag(d)		= Complex(1);
        mf4.diag(d)		= Float(1);
        mf5.diag(d)		= Float_complex(1);
        
        mf_self.diag(d) = mf_self.diag(d);

        Matrix out1		= mat;
        Matrix out2		= mat;
        Matrix out3		= mat;
        Matrix out4		= mat;
        Matrix out5		= mat;

        Matrix out_self = mat.clone();

        out1.diag(d)	= Integer(1);
        out2.diag(d)	= Real(1);
        out3.diag(d)	= Complex(1);
        out4.diag(d)	= Float(1);
        out5.diag(d)	= Float_complex(1);

        out_self.diag(d)= out_self.diag(d);

        check_struct(out1);
        check_struct(out2);
        check_struct(out3);
        check_struct(out4);
        check_struct(out5);

        check_struct(out_self);

        Real dif		= 0;
        dif				+= norm_1(out1 - mf1);
        dif				+= norm_1(out2 - mf2);
        dif				+= norm_1(out3 - mf3);
        dif				+= norm_1(out4 - mf4);
        dif				+= norm_1(out5 - mf5);

        dif				+= norm_1(out1.diag(d)-1);
        dif				+= norm_1(out2.diag(d)-1);
        dif				+= norm_1(out3.diag(d)-1);
        dif				+= norm_1(out4.diag(d)-1);
        dif				+= norm_1(out5.diag(d)-1);
        
        dif				+= norm_1(out_self - mf_self);
        dif				+= norm_1(out_self - mat);

        return dif;
    }
    catch(const std::exception& ex)
    {
        bool valid  = is_valid_diag(d, mat);

        if (valid == false)
            return 0.0;

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
Real test_function_set_diag_scal::eval_scalar(const Scalar& ,bool,int code )
{
    (void)code;
    return 0;
};

Real test_function_set_col_mat::eval_mat(const Matrix& mat,bool show_res,int code)
{
    (void)code;
    (void)show_res;

    colon_set ts(num_colons(), mat.rows() * mat.cols(), code);
    Integer size    = ts.size();
    Real dif        = 0;

    for (Integer i = 0; i < size; ++i)
    {
        colon c     = ts.get(i);
        dif         += eval_mat_col(mat, c);
    }

    return dif;
};
Real test_function_set_col_mat::eval_mat_col(const Matrix& mat, const colon& c1)
{
    m_new_objects = 0;  

    Real dif		= 0;

    try
    {				
        Integer r, c;
        {
            Matrix tmp = mat(c1);
            r = tmp.rows();
            c = tmp.cols();
        }

        using container = dynamic_mat_set::container;

        long n1 = matcl::details::no_existing_objects();
        // make sure mat is not too big for test
        if (Real(r*c) * Real(r*c) * 64 < Real(constants::max_int()))
        {
            const container& mc = ms.get(r*c);
            n1 = matcl::details::no_existing_objects() - n1;
            m_new_objects = n1;

            size_t size = mc.size();

            for (size_t i = 0; i < size ; i ++ )
            {
                Matrix mf1		= full(mat);
                Matrix out1		= mat;
                Matrix out2		= mat;

                mf1(c1)		    = mc[i];

                out1(c1)		= mc[i];
                out2(c1)		= full(mc[i]);
                Matrix mat2		= mc[i];
                
                check_struct(out1);
                check_struct(out2);

                dif				+= norm_1(out1 - mf1);
                dif				+= norm_1(out1 - out2);
            };
        }
    }
    catch(const std::exception& ex)
    {
        bool valid_1    = is_colon_valid(c1, mat.rows() * mat.cols());

        if (valid_1 == false)
            return 0.0;

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

    return dif;
};
Real test_function_set_col_mat::eval_scalar(const Scalar& ,bool,int code )
{
    (void)code;
    return 0.;
};

Real test_function_set_col2_mat::eval_mat(const Matrix& mat,bool show_res,int code)
{
    (void)code;
    (void)show_res;

    colon_set ts1(num_colons(), mat.rows(), code);
    colon_set ts2(num_colons(), mat.cols(), code + 1);

    Integer size1   = ts1.size();
    Integer size2   = ts2.size();

    Real dif        = 0;

    for (Integer i = 0; i < size1; ++i)
    for (Integer j = 0; j < size2; ++j)
    {
        colon c1    = ts1.get(i);
        colon c2    = ts2.get(j);
        dif         += eval_mat_col(mat, c1, c2);
    }

    return dif;
};
Real test_function_set_col2_mat::eval_mat_col(const Matrix& mat, const colon& c1, const colon& c2)
{
    m_new_objects = 0;

    try
    {		
        Integer r, c;
        {
            Matrix tmp = mat(c1,c2);
            r = tmp.rows();
            c = tmp.cols();
        }
        Real dif		= 0;

        using container = dynamic_mat_set::container;
        
        long n1 = matcl::details::no_existing_objects();
        const container& mc = ms.get(r,c);
        n1 = matcl::details::no_existing_objects() - n1;
        m_new_objects = n1;

        size_t size = mc.size();

        for (size_t i = 0; i < size ; i ++ )
        {
            Matrix mf1		= full(mat);
            Matrix out1		= mat;
            Matrix out2		= mat;

            mf1(c1,c2)	    = mc[i];

            out1(c1,c2)		= mc[i];
            out2(c1,c2)		= full(mc[i]);
            Matrix mat2		= mc[i];
            
            check_struct(out1);
            check_struct(out2);

            dif				+= norm_1(out1 - mf1);
            dif				+= norm_1(out1 - out2);
        };

        return dif;
    }
    catch(const std::exception& ex)
    {
        bool valid_1    = is_colon_valid(c1, mat.rows());
        bool valid_2    = is_colon_valid(c2, mat.cols());

        if (valid_1 == false || valid_2 == false)
            return 0.0;

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
Real test_function_set_col2_mat::eval_scalar(const Scalar& ,bool,int code )
{
    (void)code;
    return 0.;
};

Real test_function_set_diag_mat::eval_mat(const Matrix& mat,bool show_res,int code)
{
    (void)code;
    m_new_objects = 0;  

    Real dif		= 0;

    try
    {				
        Integer r, c;
        Matrix tmp = mat.diag(d);
        r = tmp.rows();
        c = tmp.cols();

        using container = dynamic_mat_set::container;
        long n1 = matcl::details::no_existing_objects();
        const container& mc = ms.get(r,c);
        n1 = matcl::details::no_existing_objects() - n1;
        m_new_objects = n1;

        size_t size = mc.size();

        for (size_t i = 0; i < size ; i ++ )
        {
            Matrix mf1		= full(mat);
            Matrix out1		= mat;
            Matrix out2		= mat;

            try
            {
                mf1.diag(d)	= mc[i];
            }
            catch(...)
            {
                continue;
            };
            out1.diag(d)	= mc[i];
            out2.diag(d)	= full(mc[i]);
            Matrix mat2		= mc[i];
                
            check_struct(out1);
            check_struct(out2);

            dif				+= norm_1(out1 - mf1);
            dif				+= norm_1(out1 - out2);
            dif				+= norm_1(out1.diag(d) - mat2);
            if (show_res)
            {
                matcl::out_stream  << boost::lexical_cast<std::string>(i) + " " 
                            + boost::lexical_cast<std::string>(dif) + "\n";
            };
        };		
    }
    catch(const std::exception& ex)
    {
        bool valid  = is_valid_diag(d, mat);

        if (valid == false)
            return 0.0;

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

    return dif;
};
Real test_function_set_diag_mat::eval_scalar(const Scalar& ,bool,int code )
{
    (void)code;
    return 0.;
};

Real test_function_set_int_mat::eval_mat(const Matrix& mat,bool,int code )
{
    (void)code;
    try
    {
        Matrix mf   = full(mat);
        mf(i)       = full(1);
    }
    catch(...)
    {
        return 0.;
    };

    try
    {		
        Matrix mf1_0 = full(mat);
        Matrix mf2_0 = full(mat);
        Matrix mf3_0 = full(mat);
        Matrix mf4_0 = full(mat);
        Matrix mf5_0 = full(mat);

        Matrix mf1  = full(mat);
        Matrix mf2  = full(mat);
        Matrix mf3  = full(mat);
        Matrix mf4  = full(mat);
        Matrix mf5  = full(mat);

        mf1_0(i)    = full(Integer());
        mf2_0(i)    = full(Real());
        mf3_0(i)    = full(Complex());
        mf4_0(i)    = full(Float());
        mf5_0(i)    = full(Float_complex());

        mf1(i)  = full(Integer(1));
        mf2(i)  = full(Real(1));
        mf3(i)  = full(Complex(1));
        mf4(i)  = full(Float(1));
        mf5(i)  = full(Float_complex(1));

        Matrix out1_0	= mat;
        Matrix out2_0	= mat;
        Matrix out3_0	= mat;
        Matrix out4_0	= mat;
        Matrix out5_0	= mat;

        Matrix out1		= mat;
        Matrix out2		= mat;
        Matrix out3		= mat;
        Matrix out4		= mat;
        Matrix out5		= mat;

        out1_0(i)		= full(Integer());
        out2_0(i)		= full(Real());
        out3_0(i)		= full(Complex());
        out4_0(i)		= full(Float());
        out5_0(i)		= full(Float_complex());

        out1(i)			= full(Integer(1));
        out2(i)			= full(Real(1));
        out3(i)			= full(Complex(1));
        out4(i)			= full(Float(1));
        out5(i)			= full(Float_complex(1));

        check_struct(out1_0);
        check_struct(out2_0);
        check_struct(out3_0);
        check_struct(out4_0);
        check_struct(out5_0);

        check_struct(out1);
        check_struct(out2);
        check_struct(out3);
        check_struct(out4);
        check_struct(out5);

        Real dif		= norm_1(out1_0 - mf1_0);
        dif				+= norm_1(out2_0 - mf2_0);
        dif				+= norm_1(out3_0 - mf3_0);
        dif				+= norm_1(out4_0 - mf4_0);
        dif				+= norm_1(out5_0 - mf5_0);

        dif				+= norm_1(out1 - mf1);
        dif				+= norm_1(out2 - mf2);
        dif				+= norm_1(out3 - mf3);
        dif				+= norm_1(out4 - mf4);
        dif				+= norm_1(out5 - mf5);

        dif				+= norm_1(out1_0(i));
        dif				+= norm_1(out2_0(i));
        dif				+= norm_1(out3_0(i));
        dif				+= norm_1(out4_0(i));
        dif				+= norm_1(out5_0(i));

        dif				+= norm_1(out1(i)-1);
        dif				+= norm_1(out2(i)-1);
        dif				+= norm_1(out3(i)-1);
        dif				+= norm_1(out4(i)-1);
        dif				+= norm_1(out5(i)-1);

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
Real test_function_set_int_mat::eval_scalar(const Scalar& ,bool,int code )
{
    (void)code;
    return 0.;
};

Real test_function_set_int2_mat::eval_mat(const Matrix& mat,bool,int code )
{
    (void)code;
    try
    {
        Matrix mf   = full(mat);
        mf(i,j)     = full(1);
    }
    catch(...)
    {
        return 0.;
    };

    try
    {		
        Matrix mf1_0 = full(mat);
        Matrix mf2_0 = full(mat);
        Matrix mf3_0 = full(mat);
        Matrix mf4_0 = full(mat);
        Matrix mf5_0 = full(mat);

        Matrix mf1  = full(mat);
        Matrix mf2  = full(mat);
        Matrix mf3  = full(mat);
        Matrix mf4  = full(mat);
        Matrix mf5  = full(mat);

        mf1_0(i,j)  = full(Integer());
        mf2_0(i,j)  = full(Real());
        mf3_0(i,j)  = full(Complex());
        mf4_0(i,j)  = full(Float());
        mf5_0(i,j)  = full(Float_complex());

        mf1(i,j)    = full(Integer(1));
        mf2(i,j)    = full(Real(1));
        mf3(i,j)    = full(Complex(1));
        mf4(i,j)    = full(Float(1));
        mf5(i,j)    = full(Float_complex(1));

        Matrix out1_0	= mat;
        Matrix out2_0	= mat;
        Matrix out3_0	= mat;
        Matrix out4_0	= mat;
        Matrix out5_0	= mat;

        Matrix out1		= mat;
        Matrix out2		= mat;
        Matrix out3		= mat;
        Matrix out4		= mat;
        Matrix out5		= mat;

        out1_0(i,j)		= full(Integer());
        out2_0(i,j)		= full(Real());
        out3_0(i,j)		= full(Complex());
        out4_0(i,j)		= full(Float());
        out5_0(i,j)		= full(Float_complex());

        out1(i,j)		= full(Integer(1));
        out2(i,j)		= full(Real(1));
        out3(i,j)		= full(Complex(1));
        out4(i,j)		= full(Float(1));
        out5(i,j)		= full(Float_complex(1));

        check_struct(out1_0);
        check_struct(out2_0);
        check_struct(out3_0);
        check_struct(out4_0);
        check_struct(out5_0);

        check_struct(out1);
        check_struct(out2);
        check_struct(out3);
        check_struct(out4);
        check_struct(out5);

        Real dif		= norm_1(out1_0 - mf1_0);
        dif				+= norm_1(out2_0 - mf2_0);
        dif				+= norm_1(out3_0 - mf3_0);
        dif				+= norm_1(out4_0 - mf4_0);
        dif				+= norm_1(out5_0 - mf5_0);

        dif				+= norm_1(out1 - mf1);
        dif				+= norm_1(out2 - mf2);
        dif				+= norm_1(out3 - mf3);
        dif				+= norm_1(out4 - mf4);
        dif				+= norm_1(out5 - mf5);

        dif				+= norm_1(out1_0(i,j));
        dif				+= norm_1(out2_0(i,j));
        dif				+= norm_1(out3_0(i,j));
        dif				+= norm_1(out4_0(i,j));
        dif				+= norm_1(out5_0(i,j));

        dif				+= norm_1(out1(i,j)-1);
        dif				+= norm_1(out2(i,j)-1);
        dif				+= norm_1(out3(i,j)-1);
        dif				+= norm_1(out4(i,j)-1);
        dif				+= norm_1(out5(i,j)-1);

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
Real test_function_set_int2_mat::eval_scalar(const Scalar& ,bool,int code )
{
    (void)code;
    return 0.;
};

Real test_function_get_diag_sub::eval_mat(const Matrix& mat,bool,int code )
{
    (void)code;	
    try
    {		
        Matrix out_full     = full(mat).diag(d);
        Matrix out		    = mat.diag(d);

        Integer st, diagsize;
        Integer r = mat.rows();
        Integer c = mat.cols();

        if (d >= 0)
        {
            st          = d * mat.rows() + 1;
            diagsize    = (r + d >= c) ? c - d : r;
        }
        else
        {
            st          = - d + 1;
            diagsize    = (r + d >= c) ? c : r + d;
        }

        Matrix out_test = mat(colon(st,mat.rows() + 1, end));

        if (diagsize > 0 ) 
            out_test = out_test(colon(1,diagsize));

        check_struct(out);

        Real diff       = norm_1(out - out_full);
        diff            += norm_1(out - get_diag(mat,d));
        diff            += norm_1(out - out_test);
        return diff;
    }
    catch(const std::exception& ex)
    {
        bool valid  = is_valid_diag(d, mat);

        if (valid == false)
            return 0.0;

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
Real test_function_get_diag_sub::eval_scalar(const Scalar& ,bool,int code  )
{
    (void)code;
    return 0;
};

Real test_function_drop_sparse_int::eval_mat(const Matrix& mat,bool,int code )
{
    (void)code;
    Matrix out_full;
    try
    {
        out_full    = full(mat);	
        out_full(i).drop_sparse(tol);
    }
    catch(...)
    {
        return 0.;
    };

    try
    {		
        Matrix out		= mat;
        out(i).drop_sparse(tol);
        check_struct(out);

        bool is_sp      = is_sparse_matrix(mat);
        Real dif        = 0.0;

        //do nothing for nonsparse matrices
        if (is_sp == false)
        {
            dif         += norm_1(out_full - mat);
            return dif;
        };

        Matrix mat_old  = out;
        mat_old(i)      = mat(i);
        dif             += norm_1(mat - mat_old);

        Matrix mat_new  = mat(i);
        Matrix I        = find(abs(mat_new) <= tol);
        mat_new(I)      = 0;
        
        dif		        += norm_1(out(i) - mat_new);

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
Real test_function_drop_sparse_int::eval_scalar(const Scalar& ,bool,int code )
{
    (void)code;
    return 0.;
};

Real test_function_drop_sparse_int_int::eval_mat(const Matrix& mat,bool,int code )
{
    (void)code;
    Matrix out_full;
    try
    {
        out_full    = full(mat);	
        out_full(i,j).drop_sparse(tol);
    }
    catch(...)
    {
        return 0.;
    };

    try
    {		
        Matrix out		= mat;
        out(i,j).drop_sparse(tol);
        check_struct(out);

        bool is_sp      = is_sparse_matrix(mat);
        Real dif        = 0.0;

        //do nothing for nonsparse matrices
        if (is_sp == false)
        {
            dif         += norm_1(out_full - mat);
            return dif;
        };

        Matrix mat_old  = out;
        mat_old(i,j)    = mat(i,j);
        dif             += norm_1(mat - mat_old);

        Matrix mat_new  = mat(i,j);
        Matrix I        = find(abs(mat_new) <= tol);
        mat_new(I)      = 0;
        
        dif		        += norm_1(out(i,j) - mat_new);

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
Real test_function_drop_sparse_int_int::eval_scalar(const Scalar& ,bool,int code )
{
    (void)code;
    return 0.;
};

Real test_function_drop_sparse_diag::eval_mat(const Matrix& mat,bool,int code )
{
    (void)code;

    try
    {		
        Matrix out_full;
        out_full    = full(mat);	
        out_full.diag(d).drop_sparse(tol);

        Matrix out		= mat;
        out.diag(d).drop_sparse(tol);
        check_struct(out);

        bool is_sp      = is_sparse_matrix(mat);
        Real dif        = 0.0;

        //do nothing for nonsparse matrices
        if (is_sp == false)
        {
            dif         += norm_1(out_full - mat);
            return dif;
        };

        Matrix mat_old  = out;
        mat_old.diag(d) = mat.diag(d);
        dif             += norm_1(mat - mat_old);

        Matrix mat_new  = mat.diag(d);
        Matrix I        = find(abs(mat_new) <= tol);
        mat_new(I)      = 0;
        
        dif		        += norm_1(out.diag(d) - mat_new);

        return dif;
    }
    catch(const std::exception& ex)
    {
        bool valid  = is_valid_diag(d, mat);

        if (valid == false)
            return 0.0;

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
Real test_function_drop_sparse_diag::eval_scalar(const Scalar& ,bool,int code )
{
    (void)code;
    return 0.;
};

Real test_function_drop_sparse_colon::eval_mat(const Matrix& mat,bool,int code )
{
    (void)code;

    colon_set ts(num_colons(), mat.rows() * mat.cols(), code);
    Integer size    = ts.size();
    Real dif        = 0;

    for (Integer i = 0; i < size; ++i)
    {
        colon c     = ts.get(i);
        dif         += eval_mat_col(mat, c);
    }

    return dif;
};
Real test_function_drop_sparse_colon::eval_mat_col(const Matrix& mat, const colon& c)
{
    Matrix out_full;
    try
    {
        out_full    = full(mat);	
        out_full(c).drop_sparse(tol);
    }
    catch(...)
    {
        return 0.;
    };

    try
    {		
        Matrix out		= mat;
        out(c).drop_sparse(tol);
        check_struct(out);

        bool is_sp      = is_sparse_matrix(mat);
        Real dif        = 0.0;

        //do nothing for nonsparse matrices
        if (is_sp == false)
        {
            dif         += norm_1(out_full - mat);
            return dif;
        };

        Matrix mat_old  = out;
        mat_old(c)      = mat(c);
        dif             += norm_1(mat - mat_old);

        Matrix mat_new  = mat(c);
        Matrix I        = find(abs(mat_new) <= tol);
        mat_new(I)      = 0;
        
        dif		        += norm_1(out(c) - mat_new);

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

Real test_function_drop_sparse_colon::eval_scalar(const Scalar& ,bool,int code )
{
    (void)code;
    return 0.;
};

Real test_function_drop_sparse_colon_colon::eval_mat(const Matrix& mat,bool,int code )
{
    (void)code;

    colon_set ts1(num_colons(), mat.rows(), code);
    colon_set ts2(num_colons(), mat.cols(), code + 1);

    Integer size1   = ts1.size();
    Integer size2   = ts2.size();

    Real dif        = 0;

    for (Integer i = 0; i < size1; ++i)
    for (Integer j = 0; j < size2; ++j)
    {
        colon c1    = ts1.get(i);
        colon c2    = ts2.get(j);
        dif         += eval_mat_col(mat, c1, c2);
    }

    return dif;
};
Real test_function_drop_sparse_colon_colon::eval_mat_col(const Matrix& mat, const colon& c1, const colon& c2)
{
    Matrix out_full;
    try
    {
        out_full    = full(mat);	
        out_full(c1,c2).drop_sparse(tol);
    }
    catch(...)
    {
        return 0.;
    };

    try
    {		
        Matrix out		= mat;
        out(c1,c2).drop_sparse(tol);
        check_struct(out);

        bool is_sp      = is_sparse_matrix(mat);
        Real dif        = 0.0;

        //do nothing for nonsparse matrices
        if (is_sp == false)
        {
            dif         += norm_1(out_full - mat);
            return dif;
        };

        Matrix mat_old  = out;
        mat_old(c1,c2)  = mat(c1,c2);
        dif             += norm_1(mat - mat_old);

        Matrix mat_new  = mat(c1,c2);
        Matrix I        = find(abs(mat_new) <= tol);
        mat_new(I)      = 0;
        
        dif		        += norm_1(out(c1,c2) - mat_new);

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

Real test_function_drop_sparse_colon_colon::eval_scalar(const Scalar& ,bool,int code )
{
    (void)code;
    return 0.;
};

Real test_function_add_sparse_int::eval_mat(const Matrix& mat,bool,int code )
{
    (void)code;
    Matrix out_full;
    struct_flag sf_full;

    try
    {
        out_full    = full(mat);	
        sf_full     = out_full.get_struct();
        out_full(i).add_sparse();
    }
    catch(...)
    {
        return 0.;
    };

    try
    {		
        Matrix out		= mat;
        out(i).add_sparse();
        check_struct(out);

        bool is_sp      = is_sparse_matrix(mat);
        Real dif        = 0.0;

        //do nothing for nonsparse matrices
        if (is_sp == false)
        {
            dif         += norm_1(out_full - mat);
            if (out_full.get_struct() != sf_full)
                dif     += 1;

            return dif;
        };

        dif             += norm_1(mat - out);

        if (mat.get_struct() != out.get_struct())
            dif     += 1;

        //TODO: check pattern

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
Real test_function_add_sparse_int::eval_scalar(const Scalar& ,bool,int code )
{
    (void)code;
    return 0.;
};

Real test_function_add_sparse_int_int::eval_mat(const Matrix& mat,bool,int code )
{
    (void)code;
    Matrix out_full;
    struct_flag sf_out;

    try
    {
        out_full    = full(mat);	
        sf_out      = out_full.get_struct();
        out_full(i,j).add_sparse();
    }
    catch(...)
    {
        return 0.;
    };

    try
    {		
        Matrix out		= mat;
        out(i,j).add_sparse();
        check_struct(out);

        bool is_sp      = is_sparse_matrix(mat);
        Real dif        = 0.0;

        //do nothing for nonsparse matrices
        if (is_sp == false)
        {
            dif         += norm_1(out_full - mat);
            if (out_full.get_struct() != sf_out)
                dif     += 1;

            return dif;
        };

        dif             += norm_1(mat - out);

        if (mat.get_struct() != out.get_struct())
            dif     += 1;

        //TODO: check pattern

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
Real test_function_add_sparse_int_int::eval_scalar(const Scalar& ,bool,int code )
{
    (void)code;
    return 0.;
};

Real test_function_add_sparse_diag::eval_mat(const Matrix& mat,bool,int code )
{
    (void)code;
    
    try
    {		
        Matrix out_full     = full(mat);
        struct_flag sf_out  = out_full.get_struct();

        out_full.diag(d).add_sparse();
        check_struct(out_full);
        
        Matrix out		= mat;
        out.diag(d).add_sparse();
        check_struct(out);

        bool is_sp      = is_sparse_matrix(mat);
        Real dif        = 0.0;

        //do nothing for nonsparse matrices
        if (is_sp == false)
        {
            dif         += norm_1(out_full - mat);
            if (out_full.get_struct() != sf_out)
                dif     += 1;

            return dif;
        };

        dif             += norm_1(mat - out);

        if (mat.get_struct() != out.get_struct())
            dif     += 1;

        Matrix D        = out.diag(d);

        if (D.length() != D.structural_nnz())
            dif         += 1;

        return dif;
    }
    catch(const std::exception& ex)
    {
        bool valid  = is_valid_diag(d, mat);

        if (valid == false)
            return 0.0;

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
Real test_function_add_sparse_diag::eval_scalar(const Scalar& ,bool,int code )
{
    (void)code;
    return 0.;
};
Real test_function_add_sparse_colon::eval_mat(const Matrix& mat,bool,int code )
{
    (void)code;

    colon_set ts(num_colons(), mat.rows() * mat.cols(), code);
    Integer size    = ts.size();
    Real dif        = 0;

    for (Integer i = 0; i < size; ++i)
    {
        colon c     = ts.get(i);
        dif         += eval_mat_col(mat, c);
    }

    return dif;
};
Real test_function_add_sparse_colon::eval_mat_col(const Matrix& mat,const colon& c)
{
    Matrix out_full;
    struct_flag sf_out;
    try
    {
        out_full    = full(mat);	
        sf_out      = out_full.get_struct();
        out_full(c).add_sparse();
    }
    catch(...)
    {
        return 0.;
    };

    try
    {		
        Matrix out		= mat;
        out(c).add_sparse();
        check_struct(out);

        bool is_sp      = is_sparse_matrix(mat);
        Real dif        = 0.0;

        //do nothing for nonsparse matrices
        if (is_sp == false)
        {
            dif         += norm_1(out_full - mat);
            if (out_full.get_struct() != sf_out)
                dif     += 1;

            return dif;
        };

        dif             += norm_1(mat - out);

        if (mat.get_struct() != out.get_struct())
            dif     += 1;

        Matrix S        = out(c);

        if (S.length() != S.structural_nnz() && is_sparse_matrix(S) == true)
            dif         += 1;

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

Real test_function_add_sparse_colon::eval_scalar(const Scalar& ,bool,int code )
{
    (void)code;
    return 0.;
};

Real test_function_add_sparse_colon_colon::eval_mat(const Matrix& mat,bool,int code )
{
    (void)code;

    colon_set ts1(num_colons(), mat.rows(), code);
    colon_set ts2(num_colons(), mat.cols(), code + 1);

    Integer size1   = ts1.size();
    Integer size2   = ts2.size();

    Real dif        = 0;

    for (Integer i = 0; i < size1; ++i)
    for (Integer j = 0; j < size2; ++j)
    {
        colon c1    = ts1.get(i);
        colon c2    = ts2.get(j);
        dif         += eval_mat_col(mat, c1, c2);
    }

    return dif;
};
Real test_function_add_sparse_colon_colon::eval_mat_col(const Matrix& mat, const colon& c1, const colon& c2)
{
    Matrix out_full;
    struct_flag sf_out;

    try
    {
        out_full    = full(mat);	
        sf_out      = out_full.get_struct();
        out_full(c1,c2).add_sparse();
    }
    catch(...)
    {
        return 0.;
    };

    try
    {		
        Matrix out		= mat;
        out(c1,c2).add_sparse();
        check_struct(out);

        bool is_sp      = is_sparse_matrix(mat);
        Real dif        = 0.0;

        //do nothing for nonsparse matrices
        if (is_sp == false)
        {
            dif         += norm_1(out_full - mat);
            if (out_full.get_struct() != sf_out)
                dif     += 1;

            return dif;
        };

        dif             += norm_1(mat - out);

        if (mat.get_struct() != out.get_struct())
            dif     += 1;

        Matrix S        = out(c1, c2);

        if (S.rows() * S.cols() != S.structural_nnz() && is_sparse_matrix(S) == true)
            dif         += 1;

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

Real test_function_add_sparse_colon_colon::eval_scalar(const Scalar& ,bool,int code )
{
    (void)code;
    return 0.;
};

};};
