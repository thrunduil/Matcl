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

#include "test_set_linalg.h"
#include <vector>
#include "matcl-matrep/matcl_matrep.h"

#include "matcl-core/IO/logger.h"
#include "test/test_matcl/framework/matrix_set/test_options.h"
#include "test_functions_linalg.h"
#include "test/test_matcl/framework/matrix_set/matrix_set_1.h"
#include "test/test_matcl/framework/matrix_utils.h"
#include "matcl-linalg/matcl_linalg.h"

#include "test_functions_linalg_bin.h"
#include "test/test_matcl/framework/data_struct/test_selector.h"
#include "test/test_matcl/framework/matrix_set/matrix_set_bin_1.h"

#include <boost/thread.hpp>

namespace matcl { namespace test
{

class test_linalg_bin
{
    linalg_bin_functions_list&		tf;
    const test::options&			opts;
    Integer                         thread_id;

    public:
        test_linalg_bin(linalg_bin_functions_list& tf, const test::options& opts, Integer id)
            :tf(tf),opts(opts), thread_id(id)
        {};

        test_linalg_bin(const test_linalg_bin& tb)
            :tf(tb.tf),opts(tb.opts), thread_id(tb.thread_id)
        {};

        void make()
        {
            /*
            {
                Integer code    = 153748;
                Matrix mat10    = tf.get_matrix(code).first;
                Matrix mat20    = tf.get_matrix(code).second;

                matcl::disp(mat10);
                matcl::disp(mat20);

                mat10           = randn(100, 100);
                mat20           = randn(100, 10);

                Matrix matA     = mat10;
                Matrix matB     = herprod(mat20);
                //Matrix matB   = mat20;

                //make mat2 posdef
                Integer size    = matB.rows();
                //value_code vc1  = matA.get_value_code();
                //value_code vc2  = matB.get_value_code();

                //Real diag       = sqrt(constants::eps(vc2)) * norm_1(matB);
                //Real diagB      = 1.e-10 * norm_1(matB);
                //matB            = matB + Float(diagB) * speye(size, size, vc2);

                //Real diagA      = 1.e-1 * norm_1(matA);
                //matA            = matA + Float(diagA) * speye(size, size, vc1);

                // max_trials to accommodate non-determinism of ARPACK
                Integer max_trials = 5;
                Integer k       = std::max(0, std::min(5, size - 3));
    
                Matrix U;
                Matrix T;
                bool conv = false;

                //Matrix OP       = linsolve(matB, matA);

                for (int trial = 0; trial < max_trials; trial++)
                {
                    matcl::options opts2{opt::speigs::return_nonconvergent(true)};
                    tie(U,T,conv) = pbschur(matA, matB, false, k, cluster_type::LM, opts2);

                    if (conv == true)
                        break;
                }
        
                if (conv == false)
                {
                    matcl::options opts2{opt::speigs::return_nonconvergent(false)};
                    tie(U,T,conv) = pbschur(matA, matB, false, k, cluster_type::LM, opts2);
                }

                disp(U);
                disp(T);

                for (int i = 0; i < 1; ++i)
                {
                    check_struct(U);
                    check_struct(T);

                    // check finitness
                    if (matA.all_finite() == true && matB.all_finite() == true)
                    {
                        if (U.all_finite() == false)
                            break;

                        if (T.all_finite() == false)
                            break;
                    }
                    else
                    {
                        if (U.all_finite() == true)
                            break;

                        if (T.all_finite() == true)
                            break;
                    };

                    if (has_struct_qtriu(T) == false)
                        break;

                    Real dif    = 0.;

                    // check correctness
                    Matrix mt_1     = matA * U;
                    Matrix mt_2     = U * T;

                    dif             += norm_1(mt_1 - mt_2);
                    Real tol1       = error_mult(500.0, matA, U) + error_mult(500.0, U, T);

                    if (dif < tol1)
                        dif = 0.0;

                    Matrix mt3      = ctrans(U) * matB * U;
                    dif             += norm_1(mt3 - eye(U.cols()));
                    Real tol3       = error_mult(1000.0, ctrans(U), matB, U);

                    if (dif < tol3)
                        dif         = 0.0;
                }
            };
            */

            tf.make(opts,thread_id);
        };

        void operator()()
        {
            make();
        };

    private:		
        test_linalg_bin& operator=(const test_linalg_bin&) = delete;
};

void test_linalg_bin_st(const rand_matrix_ptr& rand)
{
    test::options opts;
    try
    {
        opts.show_partial_res = false;
        opts.first_matrix_code = 0;
        opts.show_memleaks = true;

        test::mat_set_bin_1 ms(rand);
        dynamic_mat_set ms_dyn(rand);
        linalg_bin_functions_list tf(ms,ms_dyn,rand->is_nan_allowed());

        test_linalg_bin tbin (tf,opts,0);

        tbin.make();
    }
    catch(const std::exception& ex)
    {
        matcl::out_stream <<ex.what();
        matcl::out_stream.flush();
    };
};

void test_linalg_bin_mt(const rand_matrix_ptr& rand)
{
    test::options opts;
    try
    {
        opts.show_partial_res = 0;
        opts.first_matrix_code = 0;
        opts.show_memleaks = false;

        test::mat_set_bin_2 ms1(rand);
        dynamic_mat_set ms_dyn(rand);
        linalg_bin_functions_list tf(ms1,ms_dyn,rand->is_nan_allowed());

        boost::thread_group tg;

        for (int i = 0; i < 20; i++)
        {
            tg.create_thread(test_linalg_bin(tf,opts,i));
        };

        tg.join_all();
    }
    catch(const std::exception& ex)
    {
        matcl::out_stream <<ex.what();
    };
};

linalg_bin_functions_list::matrix_pair 
linalg_bin_functions_list::get_matrix(int code) const
{
    return m_tests.get_matrix(code);
};

void linalg_bin_functions_list::make(options opts, Integer thread_id)
{
    (void)thread_id;

    m_options = opts;       
    
    SELECT_TEST (3, test_gschur());    
    SELECT_TEST (3, test_gen_sym_eigen());
    SELECT_TEST (3, test_geigs());
};

//---------------------------------------------------------------
void linalg_bin_functions_list::test_gschur()
{
    test_function_gschur tf;
    Real out = m_tests.make(&tf,m_options);

    if (out == 0.)
        matcl::out_stream << std::string() + "gschur: OK" + "\n";
    else
        matcl::out_stream << std::string() + "gschur: FAILED"  + "\n";
};

void linalg_bin_functions_list::test_geigs()
{
    test_function_geigs tf;
    Real out = m_tests.make(&tf,m_options);

    if (out == 0.)
        matcl::out_stream << std::string() + "geigs: OK" + "\n";
    else
        matcl::out_stream << std::string() + "geigs: FAILED"  + "\n";
};

void linalg_bin_functions_list::test_gen_sym_eigen()
{
    test_function_gen_sym_eigen tf;
    Real out = m_tests.make(&tf,m_options);

    if (out == 0.)
        matcl::out_stream << std::string() + "gen_sym_eigen: OK" + "\n";
    else
        matcl::out_stream << std::string() + "gen_sym_eigen: FAILED"  + "\n";
};

//----------------------------------------------------------------------------
//                  test_function_gschur
//----------------------------------------------------------------------------
Real test_function_gschur::eval_mat(const Matrix& mat1,const Matrix& mat2,int code)
{
    (void)code;    

    if (mat1.is_square() == false)
        return 0.;

    if (mat2.is_square() == false)
        return 0.;

    if (mat1.rows() != mat2.rows())
        return 0.;

    Matrix aa1, bb1, q1, z1;
    Matrix aa1_compl, bb1_compl, q1_compl, z1_compl;
    Matrix aa1_sel, bb1_sel, q1_sel, z1_sel;

    Matrix temp_select;
    
    int size    = mat1.rows();

    try
    {
        tie(q1,z1,aa1,bb1)                          = gschur(mat1, mat2);
        tie(q1_compl,z1_compl,aa1_compl,bb1_compl)  = gschur_compl(mat1, mat2);
    
        //prepare valid selection vector to maintain paired eigenvalues
        {
            gschur_decomposition obj(mat1, mat2);

            temp_select = mod(floor(real(obj.eig()) * 123), 2);

            if (size > 0 && temp_select.all_finite() == false)
                temp_select = zeros(aa1.rows(), 1);
        }

        tie(q1_sel, z1_sel, aa1_sel, bb1_sel)   = ordgschur(q1, z1, aa1, bb1, temp_select);
    }
    catch(const error::error_size_geig&)
    {
        m_is_error  = true;
        m_error     = "error_size_geigr should not appear here - matrices were conformant & square";
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

    try
    {
        // check struct
        check_struct(aa1);
        check_struct(bb1);
        check_struct(q1);
        check_struct(z1);
        check_struct(aa1_compl);
        check_struct(bb1_compl);
        check_struct(q1_compl);
        check_struct(z1_compl);
        check_struct(aa1_sel);
        check_struct(bb1_sel);
        check_struct(q1_sel);
        check_struct(z1_sel);

        // check finitness
        if (mat1.all_finite() == true && mat2.all_finite() == true)
        {
            if (aa1.all_finite() == false)
                return 1.;

            if (bb1.all_finite() == false)
                return 1.;

            if (q1.all_finite() == false)
                return 1.;

            if (z1.all_finite() == false)
                return 1.;

            if (aa1_compl.all_finite() == false)
                return 1.;

            if (bb1_compl.all_finite() == false)
                return 1.;

            if (q1_compl.all_finite() == false)
                return 1.;

            if (z1_compl.all_finite() == false)
                return 1.;

            if (aa1_sel.all_finite() == false)
                return 1.;

            if (bb1_sel.all_finite() == false)
                return 1.;

            if (q1_sel.all_finite() == false)
                return 1.;

            if (z1_sel.all_finite() == false)
                return 1.;
        }
        else
        {
            if (aa1.all_finite() == true)
                return 1.;

            if (bb1.all_finite() == true)
                return 1.;

            if (q1.all_finite() == true)
                return 1.;

            if (z1.all_finite() == true)
                return 1.;

            if (aa1_compl.all_finite() == true)
                return 1.;

            if (bb1_compl.all_finite() == true)
                return 1.;

            if (q1_compl.all_finite() == true)
                return 1.;

            if (z1_compl.all_finite() == true)
                return 1.;

            if (aa1_sel.all_finite() == true)
                return 1.;

            if (bb1_sel.all_finite() == true)
                return 1.;

            if (q1_sel.all_finite() == true)
                return 1.;

            if (z1_sel.all_finite() == true)
                return 1.;
        }

        // check dimmensions
        if (aa1.is_square() == false || bb1.is_square() == false 
            || q1.is_square() == false || z1.is_square() == false 
            || aa1_compl.is_square() == false || bb1_compl.is_square() == false
            || q1_compl.is_square() == false || z1_compl.is_square() == false
            || aa1_sel.is_square() == false || bb1_sel.is_square() == false
            || q1_sel.is_square() == false || z1_sel.is_square() == false)
        {
            return 1.0;
        };

        // check if flags are set
        if (has_struct_unitary(q1) == false)
            return 1.;

        if (has_struct_unitary(z1) == false)
            return 1.;

        if (has_struct_unitary(q1_compl) == false)
            return 1.;

        if (has_struct_unitary(z1_compl) == false)
            return 1.;

        if (has_struct_unitary(q1_sel) == false)
            return 1.;

        if (has_struct_unitary(z1_sel) == false)
            return 1.;

        if (has_struct_qtriu(aa1) == false)
            return 1.0;

        if (has_struct_triu(bb1) == false)
            return 1.0;

        if (has_struct_qtriu(aa1_sel) == false)
            return 1.0;

        if (has_struct_triu(bb1_sel) == false)
            return 1.0;

        if (has_struct_triu(aa1_compl) == false)
            return 1.0;

        if (has_struct_triu(bb1_compl) == false)
            return 1.0;

        Real dif = 0.;

        //check unitary matrices
        dif     +=norm_1(q1 * ctrans(q1) - eye(size));
        dif     +=norm_1(z1 * ctrans(z1) - eye(size));
        dif     +=norm_1(q1_compl * ctrans(q1_compl) - eye(size));
        dif     +=norm_1(z1_compl * ctrans(z1_compl) - eye(size));
        dif     +=norm_1(q1_sel * ctrans(q1_sel) - eye(size));
        dif     +=norm_1(z1_sel * ctrans(z1_sel) - eye(size));

        value_code vc1  = mat1.get_value_code();
        value_code vc2  = mat2.get_value_code();
        value_code vc   = matrix_traits::unify_value_types(vc1, vc2);

        Real tol_U      = 100.0 * constants::eps(vc) * mat1.rows();

        if (dif < tol_U)
            dif         = 0.0;
        else
            return dif;

        //check if factorization is correct
        Matrix mat1_check1          = q1 * aa1 * ctrans(z1);            
        Matrix mat1_check1_compl    = q1_compl * aa1_compl * ctrans(z1_compl);            
        Matrix mat1_check1_sel      = q1_sel * aa1_sel * ctrans(z1_sel);    
        
        dif     += norm_1(mat1 - mat1_check1);
        dif     += norm_1(mat1 - mat1_check1_compl);
        dif     += norm_1(mat1 - mat1_check1_sel);

        if (dif < error_tolerance2(500.0, mat1, mat2) )
            dif = 0.0;
        else
            return dif;

        Matrix mat2_check1          = q1 * bb1 * ctrans(z1);   
        Matrix mat2_check1_compl    = q1_compl * bb1_compl * ctrans(z1_compl);   
        Matrix mat2_check1_sel      = q1_sel * bb1_sel * ctrans(z1_sel);   

        dif     += norm_1(mat2 - mat2_check1);
        dif     += norm_1(mat2 - mat2_check1_compl);
        dif     += norm_1(mat2 - mat2_check1_sel);

        if (dif < error_tolerance2(500.0, mat1, mat2) )
            dif = 0.0;
        else
            return dif;

        return dif;
    }
    catch(const std::exception& ex)
    {
        m_is_error = true;
        m_error = ex.what();
        return 1.;
    }
    catch(...)
    {
        m_is_error = true;
        m_error = "unknown error";
        return 1.;
    };
};

Real test_function_gschur::eval_scalar(const Scalar&,const Scalar&, int)
{
    return 0;
};

//----------------------------------------------------------------------------
//                  test_function_gen_sym_eigen
//----------------------------------------------------------------------------
Real test_function_gen_sym_eigen::eval_mat(const Matrix& mat1_in, const Matrix& mat2_in, 
                                           int code)
{
    (void)code;    

    if (mat1_in.is_square() == false)
        return 0.;

    if (mat2_in.is_square() == false)
        return 0.;

    if (mat1_in.rows() != mat2_in.rows())
        return 0.;

    Matrix v1, d1, v2, d2, v3, d3, mat1, mat2;

    try
    {
        mat1        = hersum(mat1_in);
        mat2        = herprod(mat2_in);

        //make mat2 posdef
        value_code vc2  = mat2.get_value_code();
        Real diag   = sqrt(constants::eps(vc2)) * norm_1(mat2);
        mat2        = mat2 + Float(diag) * speye(mat2.rows(), mat2.rows(), mat2.get_value_code());

        tie(v1, d1) = gen_sym_eigen(mat1, mat2, gschur_sym_type::A_B);
        tie(v2, d2) = gen_sym_eigen(mat1, mat2, gschur_sym_type::AB);
        tie(v3, d3) = gen_sym_eigen(mat1, mat2, gschur_sym_type::BA);
    }
    catch(const error::error_nonposdef&)
    {
        // check if matrix has in fact been non positive definite
        try
        {
            auto ret    = chol(full(mat2), false);

            // it was not
            return 1.;
        }
        catch(const error::error_nonposdef&)
        {
            return 0.;
        }

        throw;
    }
    catch(const error::error_size_geig&)
    {
        m_is_error  = true;
        m_error     = "error_size_geig should not appear here - matrices were conformant & square";
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

    int size    = mat1.rows();

    try
    {
        // check struct
        check_struct(d1);
        check_struct(d2);
        check_struct(d3);
        check_struct(v1);
        check_struct(v2);
        check_struct(v3);

        // check finitness
        if (mat1.all_finite() == true && mat2.all_finite() == true)
        {
            if (d1.all_finite() == false)
                return 1.;

            if (d2.all_finite() == false)
                return 1.;

            if (d3.all_finite() == false)
                return 1.;

            if (v1.all_finite() == false)
                return 1.;

            if (v2.all_finite() == false)
                return 1.;

            if (v3.all_finite() == false)
                return 1.;
        }
        else
        {
            if (d1.all_finite() == true)
                return 1.;

            if (d2.all_finite() == true)
                return 1.;

            if (d3.all_finite() == true)
                return 1.;

            if (v1.all_finite() == true)
                return 1.;

            if (v2.all_finite() == true)
                return 1.;

            if (v3.all_finite() == true)
                return 1.;
        }

        // check dimmensions
        if (d1.is_square() == false || v1.is_square() == false 
            || d2.is_square() == false || v2.is_square() == false 
            || d3.is_square() == false || v3.is_square() == false)
        {
            return 1.0;
        };

        // check if flags are set
        if (has_struct_diag(d1) == false)
            return 1.;

        if (has_struct_diag(d2) == false)
            return 1.;

        if (has_struct_diag(d3) == false)
            return 1.;

        Real dif    = 0.;
        Matrix one  = ones(1, 1, v3.get_value_code());

        // check correct factorizations of all variants

        // variant 1: Ax = lBx
        Matrix mat1_check1  = (one * mat2) * v1 * d1;   

        dif                 += norm_1(mat1_check1 - (one * mat1) * v1);
        Real tol1           = error_mult(500.0, mat2, v1, d1) 
                            + error_mult(500.0, mat1, v1);

        if (dif < tol1)
            dif             = 0.;
        else
            return dif;

        Matrix mat2_check1  = ctrans(v1) * mat2 * v1;   
                        
        dif                 += norm_1(mat2_check1 - eye(size));
        Real tol2           = error_mult(500.0, v1, mat2, v1);

        if (dif < tol2)
            dif             = 0.;
        else
            return dif;

        // variant 2: ABx = lx
        Matrix mat1_check2  = (one * mat1) * mat2 * v2;

        dif                 += norm_1(mat1_check2 - v2 * d2);
        Real tol3           = error_mult(500.0, mat1, mat2, v2) + error_mult(500.0, v2, d2);

        if (dif < tol3)
            dif             = 0.;
        else
            return dif;

        Matrix mat2_check2  = ctrans(v2) * mat2 * v2;        

        dif                 += norm_1(mat2_check2 - eye(size));
        Real tol4           = error_mult(500.0, v2, mat2, v2);

        if (dif < tol4)
            dif             = 0.;
        else
            return dif;

        // variant 3: BAx = lx
        Matrix mat1_check3  = (one * mat2) * mat1 * v3;

        dif                 += norm_1(mat1_check3 - v3 * d3);
        Real tol5           = error_mult(500.0, mat2, mat1, v3) + error_mult(500.0, v3, d3);

        if (dif < tol5)
            dif             = 0.;
        else
            return dif;

        Matrix chol_mat2;
        permvec perm;
                    
        Matrix mat2_check3i = linsolve(one * mat2, v3);
        Matrix mat2_check3  = ctrans(v3) * mat2_check3i;

        dif                 += norm_1(mat2_check3 - eye(size));
        Real tol6           = error_mult(10000.0, v3, mat2_check3i);

        if (dif < tol6)
            dif             = 0.;

        return dif;
    }
    catch(const std::exception& ex)
    {
        m_is_error = true;
        m_error = ex.what();
        return 1.;
    }
    catch(...)
    {
        m_is_error = true;
        m_error = "unknown error";
        return 1.;
    };
};
Real test_function_gen_sym_eigen::eval_scalar(const Scalar&,const Scalar&, int)
{
    return 0;
};

//----------------------------------------------------------------------------
//                  test_function_geigs
//----------------------------------------------------------------------------
Real test_function_geigs::eval_mat(const Matrix& mat1_in,const Matrix& mat2_in,int code)
{
    (void)code;

    if (mat1_in.is_square() == false)
        return 0.;

    if (mat2_in.is_square() == false)
        return 0.;

    if (mat1_in.rows() != mat2_in.rows())
        return 0.;

    Matrix matA     = mat1_in;
    Matrix matB     = herprod(mat2_in);

    //make mat2 posdef
    Integer size    = matB.rows();
    
    value_code vc1  = matA.get_value_code();
    value_code vc2  = matB.get_value_code();

    (void)vc1;
    (void)vc2;

    Real diagB      = sqrt(constants::eps(vc2)) * norm_1(matB);
    //Real diagB    = 1.e-1 * norm_1(matB);
    matB            = matB + Float(diagB) * speye(size, size, vc2);

    //Real diagA    = 1.e-1 * norm_1(matA);
    //matA          = matA + Float(diagA) * speye(size, size, vc1);

    // max_trials to accommodate non-determinism of ARPACK
    Integer max_trials = 5;
    Integer k       = std::max(0, std::min(5, size - 3));
    
    Matrix U;
    Matrix T;
    bool conv = false;

    try
    {
        for (int trial = 0; trial < max_trials; trial++)
        {
            matcl::options opts{opt::speigs::return_nonconvergent(true)};
            tie(U,T,conv) = pbschur(matA, matB, false, k, cluster_type::LM, opts);

            if (conv == true)
                break;
        }
        
        if (conv == false)
        {
            matcl::options opts{opt::speigs::return_nonconvergent(false)};
            tie(U,T,conv) = pbschur(matA, matB, false, k, cluster_type::LM, opts);
        }
    }
    catch(const error::invalid_speigs_k&)
    {
        //problem is too small
        return 0.0;
    }
    catch(const error::error_nonposdef&)
    {
        return 0.;
    }
    catch(const error::error_singular&)
    {
        return 0.;
    }
    catch(const error::error_size_geig& ex)
    {
        m_is_error  = true;
        m_error     = static_cast<const error::matcl_exception_linalg&>(ex).what();
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
        
    try
    {
        check_struct(U);
        check_struct(T);

        // check finitness
        if (matA.all_finite() == true && matB.all_finite() == true)
        {
            if (U.all_finite() == false)
                return 1.0;

            if (T.all_finite() == false)
                return 1.0;
        }
        else
        {
            if (U.all_finite() == true)
                return 1.0;

            if (T.all_finite() == true)
                return 1.0;
        };

        if (has_struct_qtriu(T) == false)
            return 1.0;

        Real dif    = 0.;

        // check correctness
        Matrix mt_1     = matA * U;
        Matrix mt_2     = U * T;

        dif             += norm_1(mt_1 - mt_2);
        Real tol1       = error_mult(5000.0, matA, U) + error_mult(5000.0, U, T);
                 
        // Schur vectors are highly inaccurate, we need to add quite substantial tolerance
        tol1            += error_tolerance(1000.0, matA);

        // 
        if (dif < tol1)
            dif = 0.0;
        else
            return dif;

        Matrix mt3      = ctrans(U) * matB * U;
        dif             += norm_1(mt3 - eye(U.cols()));
        Real tol3       = error_mult(5000.0, ctrans(U), matB, U);

        if (dif < tol3)
            dif         = 0.0;
        else
            return dif;

        return dif;
    }
    catch(const std::exception& ex)
    {
        m_is_error = true;
        m_error = ex.what();
        return 1.;
    }
    catch(...)
    {
        m_is_error = true;
        m_error = "unknown error";
        return 1.;
    };

    return 0;
};

Real test_function_geigs::eval_scalar(const Scalar&,const Scalar&, int)
{
    return 0;
};

};};
