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
            Integer code = 8713;
            Matrix mat1     = tf.get_matrix(code).first;
            Matrix mat2     = tf.get_matrix(code).second;

            matcl::disp(mat1);
            matcl::disp(mat2);

            Matrix out_full = kron(full(mat1) , full(mat2));	
            Matrix out      = kron(mat1 , mat2);	

            disp(out_full);
            disp(out);
            */

            tf.make(opts,thread_id);
        };

        void operator()()
        {
            make();
        };

    private:		
        test_linalg_bin& operator=(const test_linalg_bin&);
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

//
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
    {
        matcl::out_stream << std::string() + "gschur: OK" + "\n";
    }
    else
    {
        matcl::out_stream << std::string() + "gschur: FAILED"  + "\n";
    };
};

void linalg_bin_functions_list::test_geigs()
{
    test_function_geigs tf;
    Real out = m_tests.make(&tf,m_options);

    if (out == 0.)
    {
        matcl::out_stream << std::string() + "geigs: OK" + "\n";
    }
    else
    {
        matcl::out_stream << std::string() + "geigs: FAILED"  + "\n";
    };
};

void linalg_bin_functions_list::test_gen_sym_eigen()
{
    test_function_gen_sym_eigen tf;
    Real out = m_tests.make(&tf,m_options);

    if (out == 0.)
    {
        matcl::out_stream << std::string() + "gen_sym_eigen: OK" + "\n";
    }
    else
    {
        matcl::out_stream << std::string() + "gen_sym_eigen: FAILED"  + "\n";
    };
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
        Real dif = 0;  
        //check if any selection was done
        if ((sum(temp_select,1) > 0) &&  (temp_select(1,1) == 0) && 
            all(all(is_finite(q1_sel + z1_sel), 1), 2))
        {
            dif += (norm_1(q1 - q1_sel) == 0);
            dif += (norm_1(z1 - z1_sel) == 0);
            // not tested = problems with id and id-like matrices
            //dif += (norm_1(aa1 - aa1_sel) == 0);
            //dif += (norm_1(bb1 - bb1_sel) == 0);
        }
        // check finitenes of all results
        if (any(any(neg(is_finite(aa1       + bb1       + q1       + z1 +
                                  aa1_compl + bb1_compl + q1_compl + z1_compl +
                                  aa1_sel   + bb1_sel   + q1_sel   + z1_sel 
                                  )),1),2))
        {
            if (all(all(is_finite(mat1 + mat2),1),2))
            {
                dif += 1000;
            }
        }
        //check unitariness of qs, z
        dif +=norm_1(q1 * ctrans(q1) - eye(size));
        dif +=norm_1(z1 * ctrans(z1) - eye(size));
        dif +=norm_1(q1_compl * ctrans(q1_compl) - eye(size));
        dif +=norm_1(z1_compl * ctrans(z1_compl) - eye(size));
        dif +=norm_1(q1_sel * ctrans(q1_sel) - eye(size));
        dif +=norm_1(z1_sel * ctrans(z1_sel) - eye(size));
        // check correct factorizations
        Matrix mat1_check1 = q1 * aa1 * ctrans(z1);    
        dif += norm_1(mat1 - mat1_check1);
        Matrix mat1_check1_compl = q1_compl * aa1_compl * ctrans(z1_compl);    
        dif += norm_1(mat1 - mat1_check1_compl);
        Matrix mat1_check1_sel = q1_sel * aa1_sel * ctrans(z1_sel);    
        dif += norm_1(mat1 - mat1_check1_sel);
        //take roundoff error of testing by multiplication into account
        if (dif < 18 *
            ::powl(size, 3./2.) * constants::eps()
            * max(max(norm_1(mat1), norm_1(mat2)),1))
        {
            dif = 0;
        }
        
        Matrix mat2_check1 = q1 * bb1 * ctrans(z1);   
        dif += norm_1(mat2 - mat2_check1);
        Matrix mat2_check1_compl = q1_compl * bb1_compl * ctrans(z1_compl);   
        dif += norm_1(mat2 - mat2_check1_compl);
        Matrix mat2_check1_sel = q1_sel * bb1_sel * ctrans(z1_sel);   
        dif += norm_1(mat2 - mat2_check1_sel);
        //take roundoff error of testing by multiplication into account
        if (dif < 18 *
            ::powl(size, 3./2.) * constants::eps() 
            * max(max(norm_1(mat1), norm_1(mat2)),1))
        {
            dif = 0;
        }

        // check qtriu aa and triu bb
        dif += norm_1(tril(aa1,-2));
        dif += norm_1(tril(bb1,-1));
        dif += norm_1(tril(aa1_sel,-2));
        dif += norm_1(tril(bb1_sel,-1));
        // check triu for complex
        dif += norm_1(tril(aa1_compl,-1));
        dif += norm_1(tril(bb1_compl,-1));
        
        // check dimmensions (only two...)
        if (aa1.is_square() == false || bb1_sel.is_square() == false 
            || q1_compl.is_square() == false)
        {
            dif += 1;
        };
        
        if(size > 1 && mat1.structural_nnz() > 0 && mat2.structural_nnz() > 0 
            && all(all(is_finite(aa1 + bb1), 1), 2)
            && (!aa1.get_struct().is_qtriu()        || !bb1.get_struct().is_triu()
             || !is_unitary(q1.get_struct())        || !is_unitary(z1.get_struct())
             || !aa1_compl.get_struct().is_triu()   || !bb1_compl.get_struct().is_triu()
             || !is_unitary(q1_compl.get_struct())  || !is_unitary(z1_compl.get_struct())
             || !aa1_sel.get_struct().is_qtriu()    || !bb1_sel.get_struct().is_triu()
             || !is_unitary(q1_sel.get_struct())    || !is_unitary(z1_sel.get_struct())))
        {
            dif += 43247895;
        }
        
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

Real test_function_geigs::eval_mat(const Matrix& matA,const Matrix& matB_in,int code)
{
    (void)code;
    (void)matA;
    (void)matB_in;

    //TODO:
    /*
    (void)code;
    Matrix matB = matB_in;
    if ( ! matB.get_struct().is_id()) matB = ctrans(matB) * matB;

    Matrix U;
    Matrix T;
    bool conv;
    // max_trials to accommodate non-determinism of ARPACK
    Integer max_trials = 5;
    
    Integer size = matA.rows();
    Integer k = std::max(0, std::min(5, size - 3));
    int trial;
    try
    {
        for (trial = 0; trial < max_trials; trial++)
        {
            matcl::options opts{opt::speigs::return_nonconvergent(true)};
            tie(U,T,conv) = itschur(matA, matB, k, cluster_type::LM, opts);
            if (conv == true)
                break;
        }

        if (conv == false) 
        {
            matcl::options opts{opt::speigs::return_nonconvergent(false)};
            tie(U,T,conv) = itschur(matA, matB, k, cluster_type::LM, opts);
        }
    }
    catch(const error::error_nonposdef_chol&)
    {
        return 0.;
    }
    catch(const error::error_singular&)
    {
        return 0.;
    }
    catch(const error::error_size_geig& )
    {
        if (matA.is_square() && matB.is_square() && 
            matA.rows() == matB.cols() &&
            is_sym(matB)) throw;
        else return 0.;
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
    try
    {
        Real dif = 0;
        dif += norm_1(matA * V - matB * V * (size > 0 ? diag(D) : matA));
        Matrix s = svd1(matB,svd_algorithm::dc);
        Real condB = size > 1 ? (max2(s).get<1>() / min2(s).get<1>()).get_scalar<Real>() : 1.0;
        Real thresh = 1000 * constants::eps() 
                    * max(norm_1(matA), 1.0)
                    * condB
                    * (( condB > 1e7 ) ? condB : 1);
        //roundoff error of testing into account
        if (dif < thresh)
        {
            dif = 0;
        }
        // check finiteness of results
        if (any(any(neg(is_finite(V)),1),2) ||
            any(any(neg(is_finite(D)),1),2))
        {
            if ( all(all(is_finite(matA),1),2) && all(all(is_finite(matB),1),2))
            {
                return 1000;
            }
        }
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
    */
    return 0;
};
Real test_function_geigs::eval_scalar(const Scalar&,const Scalar&, int)
{
    return 0;
};

Real test_function_gen_sym_eigen::eval_mat(const Matrix& randmat1,const Matrix& randmat2,int code)
{
    (void)code;    
    Matrix  v1, d1, v2, d2, v3, d3, mat1, mat2;
    int     size    = randmat1.rows();
    try
    {
        if (randmat1.get_struct().is_id())
        {
            mat1 = randmat1;
        }
        else
        {
            mat1 = randmat1 + ctrans(randmat1);
            if(randmat1.get_value_code() == value_code::v_complex)
            {
                mat1.add_struct(predefined_struct_type::her);
            }
            else
            {
                mat1.add_struct(predefined_struct_type::sym);
            }
        }
        mat2 = randmat2 * ctrans(randmat2); // semdefinite matrix also accepted for testing!
    }
    catch(const error::invalid_eeop&)
    {
        return 0.;
    }
    try
    {
        tie(v1, d1)  =   gen_sym_eigen(mat1, mat2, gschur_sym_type::A_B);
        tie(v2, d2)  =   gen_sym_eigen(mat1, mat2, gschur_sym_type::AB);
        tie(v3, d3)  =   gen_sym_eigen(mat1, mat2, gschur_sym_type::BA);
    }
    catch(const error::error_nonposdef&)
    {
        // correct exception!
        if(mat2.structural_nnz()  == 0)
        {
            return 0.;
        }
        try
        {
            chol(full(mat2),false); // check if matrix has in fact been non positive definite
            return 1.; // it was not
        }
        catch(const error::error_nonposdef&)
        {
            return 0.;
        }
        throw;
    }
    catch(const error::error_size_geig&)
    {
        if (!randmat1.is_square() || !randmat2.is_square() || mat1.cols() != mat2.cols())
        {
            return 0.;
        }
        else
        {
            m_is_error = true;
            m_error = "error_size_geig should not appear here - matrices were conformant & square";
            return 1.;
        }
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

    try
    {
        Real dif = 0;  
        // check finitenes of all results
        if (any(any(neg(is_finite(v1 + v2 + v3 + d1 + d2 + d3)),1),2))
        {
            if (all(all(is_finite(mat1 + mat2),1),2))
            {
                dif += 1000;
            }
        }
        // check correct factorizations of all variants
        // variant 1: Ax = lBx
        Matrix mat1_check1 = mat2 * v1 * d1;    
        Matrix mat2_check1 = ctrans(v1) * full(mat2) * v1;
        dif += norm_1(mat1_check1 - mat1 * v1);
        if(all(all(is_finite(v1),1),2))
        {
            dif += norm_1(mat2_check1 - eye(size));
        }
        Real cond_mat2, inv_norm_mat_2;
        if (size > 0)
        {
            Matrix svd_mat2 = svd1(mat2,svd_algorithm::dc);
            cond_mat2 = (svd_mat2(1) /svd_mat2(size)).get_scalar<Real>();
            inv_norm_mat_2 = (1. / svd_mat2(size)).get_scalar<Real>();
        }
        else
        {
            cond_mat2 = 0;
            inv_norm_mat_2 = 0;
        }
        if (dif < 400 * 
            (max(inv_norm_mat_2, max(cond_mat2, max(norm_1(mat1), norm_1(mat2)))) + constants::eps()) *
            ::powl(size, 3./2.) * constants::eps()
            || cond_mat2 > 1. / ::sqrt(constants::eps()) )
        {
            dif = 0;
        }
        // variant 2: ABx = lx
        Matrix mat1_check2 = mat1 * 1. * mat2 * v2;
        Matrix mat2_check2 = ctrans(v2) * mat2 * v2;
        dif += norm_1(mat1_check2 - v2 * d2);
        if(all(all(is_finite(v2),1),2))
        {
            dif += norm_1(mat2_check2 - eye(size));
        }
        if (dif < 6 * 
            (cond_mat2 * max(norm_1(mat1), norm_1(mat2)) + constants::eps()) *
            ::powl(size, 3./2.) * constants::eps())
        {
            dif = 0;
        }
        // variant 3: BAx = lx
        Matrix mat1_check3 = mat2 * 1. * mat1 * v3;    
        Matrix mat2_check3;
        dif += norm_1(mat1_check3 - v3 * d3);
        if (size != 0 && all(all(is_finite(v3),1),2))
        {
            try
            {
                if (cond_mat2 < 10e13)
                {
                    Matrix chol_mat2;
                    permvec perm;
                    
                    std::tie(chol_mat2, perm) = chol(mat2,true);
                    Matrix igrek = linsolve(ctrans(chol_mat2), v3);
                    mat2_check3 = ctrans(v3) * linsolve(chol_mat2, igrek);
                    dif += norm_1(mat2_check3 - eye(size));
                }
            }
            catch (std::exception&)
            {
                disp(std::string( "Unexpected error during testing gen_sym_eigen, chol(mat2) failed") );
                throw;
            }
        }
        //take roundoff error of testing by multiplication into account
        if (dif < 400 * 
            (cond_mat2 * max(norm_1(mat1) * norm_1(mat2), max(norm_1(mat1), norm_1(mat2))) + constants::eps()) *
            ::powl(size, 3./2.) * constants::eps())
        {
            dif = 0;
        }
        // check diag d's
        dif += norm_1(tril(d1,-1) + triu(d1,1));
        dif += norm_1(tril(d2,-1) + triu(d2,1));
        dif += norm_1(tril(d3,-1) + triu(d3,1));
        // check dimmensions (only one...)
        if (d1.is_square() == false)
        {
            dif += 1;
        };
        if(size > 1 && mat1.structural_nnz() > 0 && mat2.structural_nnz() > 0 &&
            all(all(is_finite(d1), 1), 2)
            && (!d1.get_struct().is_diag() || !d2.get_struct().is_diag() || !d3.get_struct().is_diag()))
        {
            dif += 43247895;
        }
        check_struct(d1);
        check_struct(d2);
        check_struct(d3);
        check_struct(v1);
        check_struct(v2);
        check_struct(v3);
        if (dif > 10e-14)
        {
            disp(std::string( "Problematyczne macierze, A:"));
            disp(mat1);
            disp(std::string( "B:"));
            disp(mat2);
        }
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

};};
