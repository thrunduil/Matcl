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


#include "test_set_linalg.h"
#include <vector>
#include "matcl-matrep/matcl_matrep.h"

#include "matcl-core/IO/logger.h"
#include "test_functions_linalg.h"
#include "test/test_matcl/framework/matrix_set/matrix_set_1.h"
#include "test/test_matcl/framework/matrix_utils.h"
#include "matcl-linalg/matcl_linalg.h"

#include "test_functions_linalg.h"
#include "test/test_matcl/framework/data_struct/test_selector.h"

#include <boost/thread.hpp>

namespace matcl { namespace test
{

class test_linalg
{
    linalg_functions_list&		tf;
    const test::options&		opts;
    Integer                     thread_id;

    public:
        test_linalg(linalg_functions_list& tf, const test::options& opts, Integer id)
            :tf(tf),opts(opts), thread_id(id)
        {};

        test_linalg(const test_linalg& tu)
            :tf(tu.tf),opts(tu.opts), thread_id(tu.thread_id)
        {};

        void make()
        {   
            /*
            {
                int code    = 278;
                Matrix mat  = tf.get_matrix(code);

                matcl::disp(mat);
            };
            */

             tf.make(opts);
        };

        void operator()()
        {
            make();
        };

    private:		
        test_linalg& operator=(const test_linalg&) = delete;
};

void test_linalg_st(const rand_matrix_ptr& rand)
{
    test::options opts;

    try
    {
        opts.show_partial_res = 0;
        opts.first_matrix_code = 0;
        opts.show_memleaks = true;

        test::mat_set_1 ms1(rand);
        dynamic_mat_set ms(rand);
        linalg_functions_list tf(ms1,ms);
        
        test_linalg tu(tf,opts,0);
        tu.make();
    }
    catch(const std::exception& ex)
    {
        matcl::out_stream <<ex.what();
        matcl::out_stream.flush();
    };
};

void test_linalg_mt(const rand_matrix_ptr& rand)
{
    test::options opts;
    try
    {
        opts.show_partial_res = 0;
        opts.first_matrix_code = 0;
        opts.show_memleaks = false;

        test::mat_set_2 ms1(rand);
        dynamic_mat_set ms(rand);
        linalg_functions_list tf(ms1,ms);

        boost::thread_group tg;

        for (int i = 0; i < 10; i++)
        {
            tg.create_thread(test_linalg(tf,opts,i));
        };

        tg.join_all();
    }
    catch(const std::exception& ex)
    {
        matcl::out_stream <<ex.what();
    };
};

linalg_functions_list::linalg_functions_list(const matrix_set& t,dynamic_mat_set& ms)
:m_tests(t),ms(ms)
{};

void linalg_functions_list::make(options opts)
{
    m_options = opts;
    
    SELECT_TEST (3, test_svd());
    SELECT_TEST (3, test_norm());
    SELECT_TEST (3, test_cond());
    SELECT_TEST (3, test_lu_rook());    
    SELECT_TEST (3, test_lu_partial());   
    SELECT_TEST (3, test_lu_complete());    
    SELECT_TEST (3, test_chol());
    SELECT_TEST (3, test_chol_rr());    
    SELECT_TEST (3, test_cholmod());
    SELECT_TEST (3, test_qr());
    SELECT_TEST (3, test_ldl());
    
    SELECT_TEST (3, test_linsolve_0_NT());
    SELECT_TEST (3, test_linsolve_0_T());
    SELECT_TEST (3, test_linsolve_0_CT());
    SELECT_TEST (3, test_linsolve_1_NT());
    SELECT_TEST (3, test_linsolve_1_T());
    SELECT_TEST (3, test_linsolve_1_CT());
    SELECT_TEST (3, test_linsolve_3_NT());
    SELECT_TEST (3, test_linsolve_3_T());
    SELECT_TEST (3, test_linsolve_3_CT());
    SELECT_TEST (3, test_linsolve_rev());    
    SELECT_TEST (3, test_linsolve_rev2_0());       
    SELECT_TEST (3, test_linsolve_rev2_1());       
    SELECT_TEST (3, test_linsolve_rev2_3());       
    
    SELECT_TEST (3, test_hess());        
    SELECT_TEST (3, test_schur());
    SELECT_TEST (3, test_eigs());    
};

Matrix linalg_functions_list::get_matrix(int code) const
{
    return m_tests.get_matrix(code);
};

void linalg_functions_list::test_lu_partial()
{
    Real out = 0.;

    test_function_lu_partial tf;
    out += m_tests.make(&tf,m_options);

    if (out == 0.)
        matcl::out_stream << std::string() + "lu partial: OK" + "\n";
    else
        matcl::out_stream << std::string() + "lu partial: FAILED"  + "\n";
};

void linalg_functions_list::test_lu_rook()
{
    Real out = 0.;

    test_function_lu_rook tf;
    out += m_tests.make(&tf,m_options);

    if (out == 0.)
        matcl::out_stream << std::string() + "lu rook: OK" + "\n";
    else
        matcl::out_stream << std::string() + "lu rook: FAILED"  + "\n";
};

void linalg_functions_list::test_norm()
{
    Real out = 0.;

    {
        test_function_norm tf(1);
        out += m_tests.make(&tf,m_options);
    }
    {
        test_function_norm tf(2);
        out += m_tests.make(&tf,m_options);
    }
    {
        test_function_norm tf(constants::inf());
        out += m_tests.make(&tf,m_options);
    }
    {
        test_function_norm tf(-1);
        out += m_tests.make(&tf,m_options);
    }
    {
        test_function_norm tf(-2);
        out += m_tests.make(&tf,m_options);
    }

    if (out == 0.)
        matcl::out_stream << std::string() + "norm: OK" + "\n";
    else
        matcl::out_stream << std::string() + "norm: FAILED"  + "\n";
};

void linalg_functions_list::test_svd()
{
    test_function_svd tf;
    Real out = m_tests.make(&tf,m_options);

    if (out == 0.)
        matcl::out_stream << std::string() + "svd: OK" + "\n";
    else
        matcl::out_stream << std::string() + "svd: FAILED"  + "\n";
};

void linalg_functions_list::test_chol()
{
    test_function_chol tf;
    Real out = m_tests.make(&tf,m_options);

    if (out == 0.)
        matcl::out_stream << std::string() + "chol: OK" + "\n";
    else
        matcl::out_stream << std::string() + "chol: FAILED"  + "\n";
};

void linalg_functions_list::test_chol_rr()
{
    test_function_chol_rr tf;
    Real out = m_tests.make(&tf,m_options);

    if (out == 0.)
        matcl::out_stream << std::string() + "chol_rr: OK" + "\n";
    else
        matcl::out_stream << std::string() + "chol_rr: FAILED"  + "\n";
};

void linalg_functions_list::test_cholmod()
{
    test_function_cholmod tf;
    Real out = m_tests.make(&tf,m_options);

    if (out == 0.)
        matcl::out_stream << std::string() + "cholmod: OK" + "\n";
    else
        matcl::out_stream << std::string() + "cholmod: FAILED"  + "\n";
};

void linalg_functions_list::test_lu_complete()
{
    Real out = 0.;

    test_function_lu_complete tf;
    out += m_tests.make(&tf,m_options);

    if (out == 0.)
        matcl::out_stream << std::string() + "lu complete: OK" + "\n";
    else
        matcl::out_stream << std::string() + "lu complete: FAILED"  + "\n";
};

void linalg_functions_list::test_linsolve_0_NT()
{
    Real out = 0.;

    {
        test_function_linsolve tf(0, trans_type::no_trans, ms);
        out += m_tests.make(&tf,m_options);
    }

    if (out == 0.)
        matcl::out_stream << std::string() + "linsolve 0 NT: OK" + "\n";
    else
        matcl::out_stream << std::string() + "linsolve 0 NT: FAILED"  + "\n";
};

void linalg_functions_list::test_linsolve_0_T()
{
    Real out = 0.;

    {
        test_function_linsolve tf(0, trans_type::trans, ms);
        out += m_tests.make(&tf,m_options);
    }

    if (out == 0.)
        matcl::out_stream << std::string() + "linsolve 0 T: OK" + "\n";
    else
        matcl::out_stream << std::string() + "linsolve 0 T: FAILED"  + "\n";
};

void linalg_functions_list::test_linsolve_0_CT()
{
    Real out = 0.;

    {
        test_function_linsolve tf(0, trans_type::conj_trans, ms);
        out += m_tests.make(&tf,m_options);
    }

    if (out == 0.)
        matcl::out_stream << std::string() + "linsolve 0 CT: OK" + "\n";
    else
        matcl::out_stream << std::string() + "linsolve 0 CT: FAILED"  + "\n";
};

void linalg_functions_list::test_linsolve_1_NT()
{
    Real out = 0.;

    {
        test_function_linsolve tf(1, trans_type::no_trans, ms);
        out += m_tests.make(&tf,m_options);
    }

    if (out == 0.)
        matcl::out_stream << std::string() + "linsolve 1 NT: OK" + "\n";
    else
        matcl::out_stream << std::string() + "linsolve 1 NT: FAILED"  + "\n";
};

void linalg_functions_list::test_linsolve_1_T()
{
    Real out = 0.;

    {
        test_function_linsolve tf(1, trans_type::trans, ms);
        out += m_tests.make(&tf,m_options);
    }

    if (out == 0.)
        matcl::out_stream << std::string() + "linsolve 1 T: OK" + "\n";
    else
        matcl::out_stream << std::string() + "linsolve 1 T: FAILED"  + "\n";
};

void linalg_functions_list::test_linsolve_1_CT()
{
    Real out = 0.;

    {
        test_function_linsolve tf(1, trans_type::conj_trans, ms);
        out += m_tests.make(&tf,m_options);
    }

    if (out == 0.)
        matcl::out_stream << std::string() + "linsolve 1 CT: OK" + "\n";
    else
        matcl::out_stream << std::string() + "linsolve 1 CT: FAILED"  + "\n";
};

void linalg_functions_list::test_linsolve_3_NT()
{
    Real out = 0.;

    {
        test_function_linsolve tf(3, trans_type::no_trans, ms);
        out += m_tests.make(&tf,m_options);
    }

    if (out == 0.)
        matcl::out_stream << std::string() + "linsolve 3 NT: OK" + "\n";
    else
        matcl::out_stream << std::string() + "linsolve 3 NT: FAILED"  + "\n";
};

void linalg_functions_list::test_linsolve_3_T()
{
    Real out = 0.;

    {
        test_function_linsolve tf(3, trans_type::trans, ms);
        out += m_tests.make(&tf,m_options);
    }

    if (out == 0.)
        matcl::out_stream << std::string() + "linsolve 3 T: OK" + "\n";
    else
        matcl::out_stream << std::string() + "linsolve 3 T: FAILED"  + "\n";
};

void linalg_functions_list::test_linsolve_3_CT()
{
    Real out = 0.;

    {
        test_function_linsolve tf(3, trans_type::conj_trans, ms);
        out += m_tests.make(&tf,m_options);
    }

    if (out == 0.)
        matcl::out_stream << std::string() + "linsolve 3 CT: OK" + "\n";
    else
        matcl::out_stream << std::string() + "linsolve 3 CT: FAILED"  + "\n";
};

void linalg_functions_list::test_linsolve_rev()
{
    Real out = 0.;

    {
        test_function_linsolve_rev tf(0,ms);
        out += m_tests.make(&tf,m_options);
    }
    {
        test_function_linsolve_rev tf(1,ms);
        out += m_tests.make(&tf,m_options);
    }
    {
        test_function_linsolve_rev tf(3,ms);
        out += m_tests.make(&tf,m_options);
    }

    if (out == 0.)
        matcl::out_stream << std::string() + "linsolve_rev: OK" + "\n";
    else
        matcl::out_stream << std::string() + "linsolve_rev: FAILED"  + "\n";
};

void linalg_functions_list::test_linsolve_rev2_0()
{
    Real out = 0.;

    {
        test_function_linsolve_rev2 tf(0,ms);
        out += m_tests.make(&tf,m_options);
    }

    if (out == 0.)
        matcl::out_stream << std::string() + "linsolve_rev2 0: OK" + "\n";
    else
        matcl::out_stream << std::string() + "linsolve_rev2 0: FAILED"  + "\n";
};

void linalg_functions_list::test_linsolve_rev2_1()
{
    Real out = 0.;

    {
        test_function_linsolve_rev2 tf(1,ms);
        out += m_tests.make(&tf,m_options);
    }

    if (out == 0.)
        matcl::out_stream << std::string() + "linsolve_rev2 1: OK" + "\n";
    else
        matcl::out_stream << std::string() + "linsolve_rev2 1: FAILED"  + "\n";
};

void linalg_functions_list::test_linsolve_rev2_3()
{
    Real out = 0.;

    {
        test_function_linsolve_rev2 tf(3,ms);
        out += m_tests.make(&tf,m_options);
    }

    if (out == 0.)
        matcl::out_stream << std::string() + "linsolve_rev2 3: OK" + "\n";
    else
        matcl::out_stream << std::string() + "linsolve_rev2 3: FAILED"  + "\n";
};

void linalg_functions_list::test_qr()
{
    test_function_qr tf;
    Real out = m_tests.make(&tf,m_options);

    if (out == 0.)
        matcl::out_stream << std::string() + "qr: OK" + "\n";
    else
        matcl::out_stream << std::string() + "qr: FAILED"  + "\n";
};

void linalg_functions_list::test_hess()
{
    test_function_hess tf;
    Real out = m_tests.make(&tf,m_options);

    if (out == 0.)
        matcl::out_stream << std::string() + "hess: OK" + "\n";
    else
        matcl::out_stream << std::string() + "hess: FAILED"  + "\n";
};

void linalg_functions_list::test_schur()
{
    test_function_schur tf;
    Real out = m_tests.make(&tf,m_options);

    if (out == 0.)
        matcl::out_stream << std::string() + "schur: OK" + "\n";
    else
        matcl::out_stream << std::string() + "schur: FAILED"  + "\n";
};

void linalg_functions_list::test_eigs()
{
    test_function_eigs tf;
    Real out = m_tests.make(&tf,m_options);

    if (out == 0.)
        matcl::out_stream << std::string() + "eigs: OK" + "\n";
    else
        matcl::out_stream << std::string() + "eigs: FAILED"  + "\n";
};

void linalg_functions_list::test_ldl()
{
    test_function_ldl tf;
    Real out = m_tests.make(&tf,m_options);

    if (out == 0.)
        matcl::out_stream << std::string() + "ldl: OK" + "\n";
    else
        matcl::out_stream << std::string() + "ldl: FAILED"  + "\n";
};

void linalg_functions_list::test_cond()
{
    test_function_cond tf;
    Real out = m_tests.make(&tf,m_options);

    if (out == 0.)
        matcl::out_stream << std::string() + "cond: OK" + "\n";
    else
        matcl::out_stream << std::string() + "cond: FAILED"  + "\n";
};

//----------------------------------------------------------------------------
//                  test_function_norm
//----------------------------------------------------------------------------
Real test_function_norm::eval_mat(const Matrix& mat,bool ,int code)
{
    (void)code;

    Real val_f;

    try
    {		
        val_f = norm(full(mat),m_type);
    }
    catch(const std::exception& )
    {
        return 0;
    }
    
    try
    {		
        Real val    = norm(mat,m_type);
        Real dif    = 0;

        if (is_finite(val + val_f))
        {
            dif     += abs(val_f - val);

            if (dif < error_tolerance(100.0, mat))
                dif = 0.;
            else
                return dif;
        }

        dif     += (is_nan(val) == is_nan(val_f)) ? 0 : 1;
        dif     += (is_inf(val) == is_inf(val_f)) ? 0 : 1;

        if (is_inf(val) || is_inf(val_f))
            dif += (sign(val) == sign(val_f)) ? 0 : 1;

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

Real test_function_norm::eval_scalar(const Scalar&, bool,int)
{
    return 0;
};

//----------------------------------------------------------------------------
//                  test_function_lu_partial
//----------------------------------------------------------------------------
Real test_function_lu_partial::eval_mat(const Matrix& mat,bool show, int code)
{
    (void) code;
    (void)show;

    try
    {		
        Matrix L, U;
        permvec p,q;

        matcl::options opts{opt::linsolve::pivot(opt::linsolve::pivot_type::partial)};
        std::tie(L,U,p,q) = lu(mat,opts);

        Matrix P = p.to_matrix();
        Matrix Q = q.to_matrix();

        check_struct(L);
        check_struct(U);
        check_struct(P);
        check_struct(Q);

        Matrix dif = L*U - mat(colon(P),colon(Q));
        check_struct(dif);

        Real d = 0;
        
        bool fin    = mat.all_finite();
    
        if (fin == true)
            d += norm_1(dif);

        //take roundoff error into account
        if (d < error_tolerance(100.0, mat))
            d = 0;

        if (has_struct_tril(L) == false)
            d += 1;

        if (has_struct_triu(U) == false)
            d += 1;

        return d;
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

Real test_function_lu_partial::eval_scalar(const Scalar& ,bool , int code)
{
    (void) code;
    return 0;
};

//----------------------------------------------------------------------------
//                  test_function_lu_rook
//----------------------------------------------------------------------------
Real test_function_lu_rook::eval_mat(const Matrix& mat,bool , int code)
{
    (void) code;
    try
    {		
        Matrix L, U;
        permvec p,q;

        matcl::options opts{opt::linsolve::pivot(opt::linsolve::pivot_type::rook)};
        std::tie(L,U,p,q) = lu(mat,opts);

        Matrix P = p.to_matrix();
        Matrix Q = q.to_matrix();

        check_struct(L);
        check_struct(U);
        check_struct(P);
        check_struct(Q);

        Integer M = mat.rows();
        Integer N = mat.cols();
        Integer K = min(M,N);

        if (L.cols() != K)
            return 1;

        Matrix dif = L*U - mat(colon(P),colon(Q));
        check_struct(dif);

        Real d = 0;

        bool fin    = mat.all_finite();
    
        if (fin == true)
            d += norm_1(dif);
        
        //take roundoff error into account
        if (d < error_tolerance(100.0, mat))
            d = 0;

        if (has_struct_tril(L) == false)
            d += 1;

        if (has_struct_triu(U) == false)
            d += 1;

        return d;
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

Real test_function_lu_rook::eval_scalar(const Scalar& ,bool , int code)
{
    (void) code;
    return 0;
};

//----------------------------------------------------------------------------
//                  test_function_lu_complete
//----------------------------------------------------------------------------
Real test_function_lu_complete::eval_mat(const Matrix& mat,bool , int code)
{
    (void) code;
    try
    {		
        Matrix L, U;
        permvec p,q;

        matcl::options opts{opt::linsolve::pivot(opt::linsolve::pivot_type::complete)};
        std::tie(L,U,p,q) = lu(mat,opts);

        Matrix P = p.to_matrix();
        Matrix Q = q.to_matrix();

        check_struct(L);
        check_struct(U);
        check_struct(P);
        check_struct(Q);

        Integer M = mat.rows();
        Integer N = mat.cols();
        Integer K = min(M,N);

        if (L.cols() != K)
            return 1;

        Matrix dif = L*U - mat(colon(P),colon(Q));
        check_struct(dif);

        Real d = 0;

        bool fin    = mat.all_finite();
    
        if (fin == true)
            d += norm_1(dif);

        //take roundoff error into account
        if (d < error_tolerance(100.0, mat))
            d = 0;

        if (has_struct_tril(L) == false)
            d += 1;

        if (has_struct_triu(U) == false)
            d += 1;

        return d;
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

Real test_function_lu_complete::eval_scalar(const Scalar& ,bool, int code )
{
    (void) code;
    return 0;
};

//----------------------------------------------------------------------------
//                  test_function_svd
//----------------------------------------------------------------------------

Real test_function_svd::eval_mat(const Matrix& mat,bool ,int code)
{
    Real dif = 0.0;

    dif += eval_1(mat,code,true, svd_algorithm::qr);
    dif += eval_1(mat,code,false, svd_algorithm::qr);
    dif += eval_1(mat,code,true, svd_algorithm::dc);
    dif += eval_1(mat,code,false, svd_algorithm::dc);
    
    return dif;
};

Real test_function_svd::eval_1(const Matrix& mat,int code, bool economy, svd_algorithm alg)
{
    (void)code;

    Matrix U_f, S_f, V_f, S_f1;

    try
    {		
        tie(U_f,S_f,V_f)    = svd(mat,economy,alg);
        S_f1                = bdiags(svd1(mat,alg), 0, S_f.rows(), S_f.cols());
    }
    catch(...)
    {
        m_is_error = true;
        m_error = "unknown error";
        return 1.;
    };

    try
    {
        Matrix M = U_f * S_f * ctrans(V_f);

        check_struct(U_f);
        check_struct(S_f);
        check_struct(S_f1);
        check_struct(V_f);
                
        Real dif = norm_1(M - mat);
        
        // test singular values-only version
        dif += norm_1(S_f - S_f1);
        
        dif += norm_1(ctrans(U_f)*U_f - speye(U_f.cols()));
        dif += norm_1(ctrans(V_f)*V_f - speye(V_f.cols()));

        if (economy == false)
        {
            dif += norm_1(U_f*ctrans(U_f) - speye(U_f.rows()));
            dif += norm_1(V_f*ctrans(V_f) - speye(V_f.rows()));
        };

        // take roundoff differences between different algorithms used into account
        if (dif < error_tolerance(100.0, mat))
            dif = 0;
        else
            return dif;

        // check finitenes of all results
        if (mat.all_finite() == true)
        {
            if (U_f.all_finite() == false)
                dif += 1;

            if (V_f.all_finite() == false)
                dif += 1;

            if (S_f.all_finite() == false)
                dif += 1;

            if (S_f1.all_finite() == false)
                dif += 1;
        };
        
        if (U_f.is_square() == false && economy == false)
            dif += 1;

        if (V_f.is_square() == false && economy == false)
            dif += 1;

        if (economy == false && U_f.all_finite() == true &&
            U_f.rows() > 1 && is_unitary(U_f.get_struct()) == false && U_f.get_struct().is_id() == false)
        {
            dif += 1;
        };

        if (economy == false && V_f.all_finite() == true &&
            V_f.rows() > 1 && is_unitary(V_f.get_struct()) == false && V_f.get_struct().is_id() == false)
        {
            dif += 1;
        };        

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

Real test_function_svd::eval_scalar(const Scalar&, bool,int)
{
    return 0;
};

//----------------------------------------------------------------------------
//                  test_function_cond
//----------------------------------------------------------------------------
Real test_function_cond::eval_mat(const Matrix& mat,bool ,int code)
{
    (void)code;
    Real res;

    try
    {
        res = cond(mat);
    }
    catch(const error::square_matrix_required& )
    {
        if (mat.is_square() ) 
        {
            m_is_error = true;
            m_error = "error::square_matrix_required should not appear here - matrix was square";
            return 1.;
        }
        else
        {
            return 0.;
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
        Real dif        = 0;
        Matrix s        = svd1(mat,svd_algorithm::dc);
        Real test_cond  = (mat.rows() > 1 ) ? (s(1) / s(matcl::end)).get_scalar<Real>() : 1;
        dif             += norm_1(test_cond - res);
        
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

Real test_function_cond::eval_scalar(const Scalar&, bool,int)
{
    return 0;
};

//----------------------------------------------------------------------------
//                  test_function_chol
//----------------------------------------------------------------------------
Real test_function_chol::eval_mat(const Matrix& mat,bool ,int code)
{
    (void)code;

    Real out    = 0.;
    out         += eval(mat,true);
    out         += eval(mat,false);

    return out;
};

Real test_function_chol::eval(const Matrix& mat, bool upper)
{
    Matrix S;
    permvec p;

    Matrix M0 = mat;

    if (mat.get_value_code() == matcl::value_code::v_integer)
    {
        M0 = matcl::convert(M0,matcl::matrix_traits::get_matrix_type(matcl::value_code::v_real, 
                                                                     mat.get_struct_code()));
    };

    Matrix M        = herprod(M0, true);
    value_code vc   = M.get_value_code();
    Matrix r        = make_scalar(100.0f * epsilon_mat(M), vc) * speye(M.rows(), vc);

    M = M + r;

    try
    {		
        std::tie(S,p) = chol(M, upper);
    }
    catch(const std::exception& )
    {
        return 1;
    }
    catch(...)
    {
        m_is_error = true;
        m_error = "unknown error";
        return 1.;
    };

    try
    {
        Matrix M1   = upper? ctrans(S) * S : S * ctrans(S);
        check_struct(S);

        Real dif = norm_1(M(p, p) - M1);

        if (dif < error_tolerance(100.0, M1))
            dif = 0.0;
        else
            return dif;

        if (M1.is_square() == false)
            dif += 1;
    
        if (S.is_square() == false)
            dif += 1;

        if (upper == true && has_struct_triu(S) == false)
            dif += 1.0;

        if (upper == false && has_struct_tril(S) == false)
            dif += 1.0;

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

Real test_function_chol::eval_scalar(const Scalar&, bool,int)
{
    return 0;
};

//----------------------------------------------------------------------------
//                  test_function_chol_rr
//----------------------------------------------------------------------------
Real test_function_chol_rr::eval_mat(const Matrix& mat,bool ,int code)
{
    (void)code;

    Real out = 0.;

    out += eval(mat,true);
    out += eval(mat,false);

    return out;
};

Real test_function_chol_rr::eval(const Matrix& mat,bool upper)
{
    Matrix S;
    permvec p;
    Integer rank;
    Matrix M0 = mat;

    if (mat.get_value_code() == matcl::value_code::v_integer)
    {
        M0 = matcl::convert(M0,matcl::matrix_traits::get_matrix_type(matcl::value_code::v_real, mat.get_struct_code()));
    };
    
    Matrix M = herprod(M0, true);

    try
    {		
        std::tie(S,p,rank) = chol_rr(M,upper);
    }
    catch(const std::exception& )
    {
        return 1;
    }
    catch(...)
    {
        m_is_error = true;
        m_error = "unknown error";
        return 1.;
    };

    try
    {
        Matrix M1 = upper ? ctrans(S) * S : S * ctrans(S);
        check_struct(S);

        Real dif = norm_1(M(p, p) - M1);

        if (dif < error_tolerance(300.0, M1))
            dif = 0.0;
        else
            return dif;

        if (M1.is_square() == false)
            dif += 1;

        if (upper == true && has_struct_triu(S) == false)
            dif += 1.0;

        if (upper == false && has_struct_tril(S) == false)
            dif += 1.0;

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

Real test_function_chol_rr::eval_scalar(const Scalar&, bool,int)
{
    return 0;
};

//----------------------------------------------------------------------------
//                  test_function_cholmod
//----------------------------------------------------------------------------
Real test_function_cholmod::eval_mat(const Matrix& mat,bool show,int code)
{
    (void)code;

    Real out = 0.;
    out += eval(mat,true,show);
    out += eval(mat,false,show);

    return out;
};

Real test_function_cholmod::eval(const Matrix& mat,bool upper, bool show)
{
    Matrix S, est_rank, norm_E;
    permvec p;
    Matrix M;

    try
    {		 
        M = hersum(mat);
        std::tie(S,p,est_rank,norm_E) = cholmod(M,upper);
    }
    catch(const std::exception& )
    {
        return 0;
    }
    catch(...)
    {
        m_is_error = true;
        m_error = "unknown error";
        return 1.;
    };

    try
    {
        check_struct(S);
        Real dif = 0;

        // check factorization M(p,p) + E = S'S
        Matrix M1   = upper ? ctrans(S) * S : S * ctrans(S);
        int N       = M.rows();

        Matrix E = M(p, p) - M1;
        dif = norm_1(norm(E, 2.) - norm_E.get_scalar<Real>());
        
        if (dif < error_tolerance(10.0, M1))
            dif = 0.;
        else
            return dif;

        // check estimated rank is plausible
        if (est_rank.get_scalar<Integer>() < 0 || est_rank.get_scalar<Integer>() > N)
            dif += 1000.0;

        // check bound on ||E||
        Matrix eig_M = schur_decomposition(M).eig();

        Real tol_E  = error_tolerance(10.0, eig_M);

        if (N > 0)
            tol_E   += 100.0 * abs(min_d(eig_M)).get_scalar<Real>() * N; 

        if (N > 0 && norm_E > tol_E)
        {
            dif += 10000;

            if(show)
                disp(std::string("Bound on ||E|| not satisfied for diagonal matrix"));
        }

        // check finitenes of all results
        if (mat.all_finite() == true)
        {
            if (S.all_finite() == false)
                dif     += 1.0;

            if (est_rank.all_finite() == false)
                dif     += 1.0;

            if (norm_E.all_finite() == false)
                dif     += 1.0;
        }

        if (M1.is_square() == false)
            dif += 1;

        if (upper == true && has_struct_triu(S) == false)
            dif += 1.0;

        if (upper == false && has_struct_tril(S) == false)
            dif += 1.0;

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

Real test_function_cholmod::eval_scalar(const Scalar&, bool,int)
{
    return 0;
};

//----------------------------------------------------------------------------
//                  test_function_qr
//----------------------------------------------------------------------------
Real test_function_qr::eval_mat(const Matrix& mat,bool ,int code)
{
    (void)code;

    // number x means result of tuple-x version
    unitary_matrix q2;
    unitary_matrix q3;
    unitary_matrix q2eco;
    unitary_matrix q3eco;
    Matrix r1_internal;
    Matrix r1;
    Matrix r2;
    Matrix r3;
    Matrix r2eco;
    Matrix r3eco;
    permvec e3;
    permvec e3eco;

    try
    {		
        r1_internal             = qr_internal(mat);
        r1                      = qr(mat);
        tie(q2, r2)             = qr2(mat);
        tie(q3, r3, e3)         = qr3(mat);
        tie(q2eco, r2eco)       = qr2(mat, true);
        tie(q3eco, r3eco, e3eco)= qr3 (mat, true);
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
        Matrix Q2   = q2.to_matrix();
        Matrix Q2_e = q2eco.to_matrix();
        Matrix Q3   = q3.to_matrix();
        Matrix Q3_e = q3eco.to_matrix();

        Real dif = 0;

        // check finitenes of all results

        if (mat.all_finite() == true)
        {
            if (r1_internal.all_finite() == false)
                dif += 1.0;

            if (r1.all_finite() == false)
                dif += 1.0;

            if (q2.all_finite() == false)
                dif += 1.0;

            if (r2.all_finite() == false)
                dif += 1.0;

            if (q3.all_finite() == false)
                dif += 1.0;

            if (r3.all_finite() == false)
                dif += 1.0;

            if (q2eco.all_finite() == false)
                dif += 1.0;

            if (r2eco.all_finite() == false)
                dif += 1.0;

            if (q3eco.all_finite() == false)
                dif += 1.0;

            if (r3eco.all_finite() == false)
                dif += 1.0;

            if (Q2.all_finite() == false)
                dif += 1.0;

            if (Q2_e.all_finite() == false)
                dif += 1.0;

            if (Q3.all_finite() == false)
                dif += 1.0;

            if (Q3_e.all_finite() == false)
                dif += 1.0;

            if (dif != 0.0)
                return dif;
        };

        if (mat.all_finite() == false)
        {
            if (r1_internal.all_finite() == true)
                dif += 1.0;

            if (r1.all_finite() == true)
                dif += 1.0;

            if (q2.all_finite() == true)
                dif += 1.0;

            if (r2.all_finite() == true)
                dif += 1.0;

            if (q3.all_finite() == true)
                dif += 1.0;

            if (r3.all_finite() == true)
                dif += 1.0;

            if (q2eco.all_finite() == true)
                dif += 1.0;

            if (r2eco.all_finite() == true)
                dif += 1.0;

            if (q3eco.all_finite() == true)
                dif += 1.0;

            if (r3eco.all_finite() == true)
                dif += 1.0;

            if (Q2.all_finite() == true)
                dif += 1.0;

            if (Q2_e.all_finite() == true)
                dif += 1.0;

            if (Q3.all_finite() == true)
                dif += 1.0;

            if (Q3_e.all_finite() == true)
                dif += 1.0;

            if (dif != 0.0)
                return dif;
        }

        // check dimmensions
        if (q2.is_square() == false)
            dif += 1;

        if (q3.is_square() == false)
            dif += 1;

        if (q2eco.cols() != min(mat.cols(), mat.rows()))
            dif +=1;

        if (q3eco.cols() != min(mat.cols(), mat.rows()))
            dif +=1;

        if (r2eco.rows() != min(mat.cols(), mat.rows()))
            dif +=1;

        if (r3eco.rows() != min(mat.cols(), mat.rows()))
            dif +=1;

        // check uppertriangular r1 r2 r3 (r1_internal is NOT uppertriangular)
        if (has_struct_triu(r1) == false)
            dif     += 1.0;

        if (has_struct_triu(r2) == false)
            dif     += 1.0;

        if (has_struct_triu(r3) == false)
            dif     += 1.0;

        if (has_struct_triu(r2eco) == false)
            dif     += 1.0;

        if (has_struct_triu(r3eco) == false)
            dif     += 1.0;

        if (has_struct_unitary(Q2) == false)
            dif     += 1.0;

        if (has_struct_unitary(Q3) == false)
            dif     += 1.0;

        if (Q2_e.is_square() && has_struct_unitary(Q2_e) == false)
            dif     += 1.0;

        if (Q3_e.is_square() && has_struct_unitary(Q3_e) == false)
            dif     += 1.0;

        if (dif != 0.0)
            return dif;

        // check struct
        check_struct(r1_internal);

        check_struct(r1);
        check_struct(r2);
        check_struct(r3);
        check_struct(Q2);
        check_struct(Q3);
        check_struct(r2eco);
        check_struct(r3eco);
        check_struct(Q2_e);
        check_struct(Q3_e);

        if (mat.all_finite() == false)
            return dif;

        //check unitariness of q
        dif         +=norm_1(ctrans(Q2) * Q2 - eye(q2.cols()));
        dif         +=norm_1(ctrans(Q3) * Q3 - eye(q3.cols()));
        dif         +=norm_1(ctrans(Q2_e) * Q2_e - eye(q2eco.cols()));
        dif         +=norm_1(ctrans(Q3_e) * Q3_e - eye(q3eco.cols()));

        value_code vc   = mat.get_value_code();
        Real tol_U      = 10.0 * constants::eps(vc) * Q2.rows();

        if (dif < tol_U)
            dif     = 0.0;
        else
            return dif;

        // check correct factorizations
        Matrix mat_check1_internal = q2 * triu(r1_internal);
        Matrix mat_check1       = q2 * r1;
        Matrix mat_check2       = q2 * r2;
        Matrix mat_check3p      = q3 * r3;
        Matrix mat_check2eco    = q2eco * r2eco;
        Matrix mat_check3ecop   = q3eco * r3eco;

        dif         += norm_1(mat - mat_check1_internal);
        dif         += norm_1(mat - mat_check1);
        dif         += norm_1(mat - mat_check2);
        dif         += norm_1(mat(colon(), e3) - mat_check3p);
        dif         += norm_1(mat - mat_check2eco);
        dif         += norm_1(mat(colon(),e3eco) - mat_check3ecop);

        if (dif < error_tolerance(600.0, mat))
            dif     = 0.0;

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
Real test_function_qr::eval_scalar(const Scalar&, bool,int)
{
    return 0;
};

//----------------------------------------------------------------------------
//                  test_function_ldl
//----------------------------------------------------------------------------
Real test_function_ldl::eval_mat(const Matrix& mat_temp,bool ,int code)
{
    Real dif = 0;
    dif     += eval_1(mat_temp, code, true);
    dif     += eval_1(mat_temp, code, false);
    return dif;
};

Real test_function_ldl::eval_1(const Matrix& mat_temp,int code, bool upper)
{
    (void)code;
    Matrix mat, matherm;
    try
    {
        if (mat_temp.get_struct().is_id())
        {
            mat = mat_temp; 
            matherm = mat_temp;
        }
        else
        {
            if (mat_temp.rows() == mat_temp.cols())
            {
                mat     = symsum(mat_temp);
                matherm = hersum(mat_temp);
            }
            else
            {
                mat     = symprod(mat_temp);
                matherm = herprod(mat_temp);
            };
        }
    }
    catch(const error::invalid_eeop&)
    {
        return 0.;
    }

    Matrix l,d;
    Matrix lh,dh;
    permvec p, ph;

    try
    {
        tie(l,d,p)            =   ldl(mat, upper);
        tie(lh,dh,ph)         =   ldl_herm(matherm, upper);
    }
    catch(const error::error_size_ldl&)
    {
        if (!mat.is_square())
        {
            return 0.;
        }
        else
        {
            m_is_error = true;
            m_error = "error_size_ldl should not appear here - matrix was square";
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

        //check unit lower-triangularness of l
        if (upper == false && has_struct_tril(l) == false)
            dif     += 1.0;

        if (upper == false && has_struct_tril(lh) == false)
            dif     += 1.0;

        if (upper == true && has_struct_triu(l) == false)
            dif     += 1.0;

        if (upper == true && has_struct_triu(lh) == false)
            dif     += 1.0;

        if (has_struct_qtril(d) == false)
            dif     += 1.0;

        if (has_struct_qtriu(d) == false)
            dif     += 1.0;

        if (has_struct_qtril(dh) == false)
            dif     += 1.0;

        if (has_struct_qtriu(dh) == false)
            dif     += 1.0;

        dif += norm_1(get_diag(l) - 1.0);
        dif += norm_1(get_diag(lh) - 1.0);

        if (dif != 0.0)
            return dif;

        check_struct(l);
        check_struct(d);
        check_struct(lh);
        check_struct(dh);

        // check finitenes of all results
        if (mat.all_finite() == true)
        {
            if (l.all_finite() == false)
                dif     += 1.0;

            if (d.all_finite() == false)
                dif     += 1.0;

            if (lh.all_finite() == false)
                dif     += 1.0;

            if (dh.all_finite() == false)
                dif     += 1.0;
        };

        if (mat.all_finite() == false)
        {
            if (l.all_finite() == true)
                dif     += 1.0;

            if (d.all_finite() == true)
                dif     += 1.0;

            if (lh.all_finite() == true)
                dif     += 1.0;

            if (dh.all_finite() == true)
                dif     += 1.0;
        };

        // check dimmensions
        if (d.is_square() == false  || l.is_square() == false
            || dh.is_square() == false || lh.is_square() == false)
        {
            dif += 1;
        };

        if (dif != 0.0)
            return dif;

        if (mat.all_finite() == false)
            return dif;

        // check correctness of factorization   
        Matrix mat_check_tr;
        Matrix mat_checkh_tr;

        mat_check_tr    = l * d * trans(l);  
        mat_checkh_tr   = lh * dh * ctrans(lh);

        dif             += norm_1(mat(p,p) - mat_check_tr);

        if (dif < error_tolerance(100.0, mat))
            dif         = 0.0;         
        else
            return dif;

        dif             += norm_1(matherm(ph,ph) - mat_checkh_tr);               

        if (dif < error_tolerance(100.0, matherm))
            dif         = 0.0;         

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

Real test_function_ldl::eval_scalar(const Scalar&, bool,int)
{
    return 0;
};

//----------------------------------------------------------------------------
//                  test_function_linsolve
//----------------------------------------------------------------------------
Real test_function_linsolve::eval_mat(const Matrix& mat,bool , int code)
{
    (void)code;
    m_new_objects = 0;

    if (mat.is_square() == false)
        return 0.;

    Integer r = mat.rows();

    try
    {	
        Real dif		= 0;
        using container = dynamic_mat_set::container;
        
        long n1         = matcl::details::no_existing_objects();
        const container& mc = ms.get(r,K);

        n1              = matcl::details::no_existing_objects() - n1;
        m_new_objects   = n1;

        size_t size = mc.size();

        for (size_t i = 0; i < size ; i ++ )
            dif += eval_mat_impl(mat,mc[i]);

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

Real test_function_linsolve::eval_mat_impl(const Matrix& A,const Matrix& B)
{
    Real out    = 0.;

    permvec p    = randperm(A.rows());
    permvec q    = randperm(A.cols());
    permvec e    = permvec::identity(A.rows());

    if (m_trans == trans_type::no_trans)
    {
        out         = out + test_notrans(A,B,e,e);
        out         = out + test_notrans(A,B,p,e);
        out         = out + test_notrans(A,B,e,q);
        out         = out + test_notrans(A,B,p,q);
    };

    if (m_trans == trans_type::trans)
    {
        out         = out + test_trans(A,B,e,e);
        out         = out + test_trans(A,B,p,e);
        out         = out + test_trans(A,B,e,q);
        out         = out + test_trans(A,B,p,q);
    };

    if (m_trans == trans_type::conj_trans)
    {
        out         = out + test_ctrans(A,B,e,e);
        out         = out + test_ctrans(A,B,p,e);
        out         = out + test_ctrans(A,B,e,q);
        out         = out + test_ctrans(A,B,p,q);
    };

    return out;
};

Real test_function_linsolve::test_notrans(const Matrix& A,const Matrix& B, permvec p, permvec q)
{
    try
    {	
        Matrix X       = linsolve(A, p, q, B);

        check_struct(X);

        Real dif		= 0;              

        if (A.all_finite() == true && B.all_finite() == true)
        {
            // overflow is possible
            if (X.all_finite() == false)
            {
                Matrix smin     = svd1(A, svd_algorithm::dc)(end);
                Real eps        = constants::eps(X.get_value_code());

                if (smin < eps)
                    return 0.0;
                else
                    return 1.;
            };
        }
        else
        {
            if (X.all_finite() == true)
                return 1.;
        }

        Matrix test     = A(invperm(p).to_matrix(), invperm(q).to_matrix())*X - B;

        dif             += norm_1(test);
        Real tol        = norm(A, 1) * norm(X,1) + norm(B, 1);
        Real eps        = constants::eps(X.get_value_code());

        if (dif < 100. * tol * eps)
            dif         = 0.;

        return dif;                
    }
    catch(error::error_lsolve& ex)
    {
        if (A.is_square() && A.rows() == B.rows())
        {
            m_is_error = true;

            m_error = static_cast<error::matcl_exception_linalg&>(ex).what();
            return 1.;
        };

        return 0.;
    }
    catch(error::error_singular&)
    {
        Integer size    = min(A.rows(),A.cols());

        if (size == 0)
            return 0.;

        Matrix s        = svd1(A, svd_algorithm::dc);

        if (s(size) == 0)
            return 0.;

        Matrix cond     = s(1) / s(size);
        Real tol        = 1.0 / constants::eps(A.get_value_code());

        if (cond > 1.e-2 * tol)
            return 0;

        m_is_error = true;
        m_error = "singular matrix";
        return 1.;
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

Real test_function_linsolve::test_trans(const Matrix& A,const Matrix& B, permvec p, permvec q)
{
    try
    {	
        Matrix Xf       = linsolve(A, p, q, B,trans_type::trans);
        check_struct(Xf);

        Real dif		= 0;
        
        if (A.all_finite() == true && B.all_finite() == true)
        {
            // overflow is possible
            if (Xf.all_finite() == false)
            {
                Matrix smin     = svd1(A, svd_algorithm::dc)(end);
                Real eps        = constants::eps(Xf.get_value_code());

                if (smin < eps)
                    return 0.0;
                else
                    return 1.;
            };
        }
        else
        {
            if (Xf.all_finite() == true)
                return 1.;
        }
       
        Matrix test     = matcl::trans(A(invperm(p), invperm(q))) * Xf - B;

        dif             += norm_1(test);
        Real tol        = norm(A, 1) * norm(Xf,1) + norm(B, 1);
        Real eps        = constants::eps(Xf.get_value_code());

        if (dif < 100. * tol * eps)
            dif         = 0.;

        return dif;
    }
    catch(error::error_lsolve& ex)
    {
        if (A.is_square() && A.rows() == B.rows())
        {
            m_is_error = true;
            m_error = static_cast<error::matcl_exception_linalg&>(ex).what();
            return 1.;
        };

        return 0.;
    }
    catch(error::error_singular&)
    {
        Integer size    = min(A.rows(),A.cols());

        if (size == 0)
            return 0.;

        Matrix s        = svd1(A, svd_algorithm::dc);

        if (s(size) == 0)
            return 0.;

        Matrix cond     = s(1) / s(size);
        Real tol        = 1.0 / constants::eps(A.get_value_code());

        if (cond > 1.e-2 * tol)
            return 0;

        m_is_error = true;
        m_error = "singular matrix";
        return 1.;
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

Real test_function_linsolve::test_ctrans(const Matrix& A,const Matrix& B, permvec p, permvec q)
{
    try
    {	
        Matrix Xf       = linsolve(A, p, q, B,trans_type::conj_trans);
        check_struct(Xf);

        Real dif		= 0;
        
        if (A.all_finite() == true && B.all_finite() == true)
        {
            // overflow is possible
            if (Xf.all_finite() == false)
            {
                Matrix smin     = svd1(A, svd_algorithm::dc)(end);
                Real eps        = constants::eps(Xf.get_value_code());

                if (smin < eps)
                    return 0.0;
                else
                    return 1.;
            };
        }
        else
        {
            if (Xf.all_finite() == true)
                return 1.;
        }
       
        Matrix test     = matcl::ctrans(A(invperm(p), invperm(q)))*Xf - B;

        dif             += norm_1(test);
        Real tol        = norm(A, 1) * norm(Xf,1) + norm(B, 1);
        Real eps        = constants::eps(Xf.get_value_code());

        if (dif < 100. * tol * eps)
            dif         = 0.;
        else
            return dif;
                
        return dif;
    }
    catch(error::error_lsolve& ex)
    {
        if (A.is_square() && A.rows() == B.rows())
        {
            m_is_error = true;
            m_error = static_cast<error::matcl_exception_linalg&>(ex).what();
            return 1.;
        };

        return 0.;
    }
    catch(error::error_singular&)
    {
        Integer size    = min(A.rows(),A.cols());

        if (size == 0)
            return 0.;

        Matrix s        = svd1(A, svd_algorithm::dc);

        if (s(size) == 0)
            return 0.;

        Matrix cond     = s(1) / s(size);
        Real tol        = 1.0 / constants::eps(A.get_value_code());

        if (cond > 1.e-2 * tol)
            return 0;

        m_is_error = true;
        m_error = "singular matrix";
        return 1.;
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

Real test_function_linsolve::eval_scalar(const Scalar& ,bool, int code )
{
    (void) code;
    return 0;
};

//----------------------------------------------------------------------------
//                  test_function_linsolve_rev
//----------------------------------------------------------------------------
Real test_function_linsolve_rev::eval_mat(const Matrix& mat,bool , int code)
{
    (void)code;
    m_new_objects = 0;

    if (mat.is_square() == false)
        return 0.;

    Integer r = mat.rows();

    try
    {	
        Real dif		= 0;
        using container = dynamic_mat_set::container;
        
        long n1             = matcl::details::no_existing_objects();
        const container& mc = ms.get(r,K);
        n1                  = matcl::details::no_existing_objects() - n1;
        m_new_objects       = n1;

        size_t size = mc.size();

        for (size_t i = 0; i < size ; i ++ )
            dif += eval_mat_impl(mat, mc[i]);

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

Real test_function_linsolve_rev::eval_mat_impl(const Matrix& A,const Matrix& B)
{
    Real out    = 0.;

    permvec p    = randperm(A.rows());
    permvec q    = randperm(A.cols());
    permvec e    = permvec::identity(A.rows());

    out         = out + test_notrans(A,B,e,e);
    out         = out + test_notrans(A,B,p,e);
    out         = out + test_notrans(A,B,e,q);
    out         = out + test_notrans(A,B,p,q);

    return out;
};

Real test_function_linsolve_rev::test_notrans(const Matrix& A,const Matrix& B, permvec p, permvec q)
{
    try
    {	
        Matrix Xf       = linsolve_rev(A, p, q, trans(B));

        check_struct(Xf);

        Real dif		= 0;              

        if (A.all_finite() == true && B.all_finite() == true)
        {
            // overflow is possible
            if (Xf.all_finite() == false)
            {
                Matrix smin     = svd1(A, svd_algorithm::dc)(end);
                Real eps        = constants::eps(Xf.get_value_code());

                if (smin < eps)
                    return 0.0;
                else
                    return 1.;
            };
        }
        else
        {
            if (Xf.all_finite() == true)
                return 1.;
        }

        Matrix test     = Xf * A (invperm(p), invperm(q))- trans(B);

        dif             += norm_1(test);
        Real tol        = norm(A, 1) * norm(Xf,1) + norm(B, 1);
        Real eps        = constants::eps(Xf.get_value_code());

        if (dif < 100. * tol * eps)
            dif         = 0.;

        return dif;                
    }
    catch(error::error_lsolve& ex)
    {
        if (A.is_square() && A.rows() == B.rows())
        {
            m_is_error = true;

            m_error = static_cast<error::matcl_exception_linalg&>(ex).what();
            return 1.;
        };

        return 0.;
    }
    catch(error::error_singular&)
    {
        Integer size    = min(A.rows(),A.cols());

        if (size == 0)
            return 0.;

        Matrix s        = svd1(A, svd_algorithm::dc);

        if (s(size) == 0)
            return 0.;

        Matrix cond     = s(1) / s(size);
        Real tol        = 1.0 / constants::eps(A.get_value_code());

        if (cond > 1.e-2 * tol)
            return 0;

        m_is_error = true;
        m_error = "singular matrix";
        return 1.;
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

Real test_function_linsolve_rev::eval_scalar(const Scalar& ,bool, int code )
{
    (void) code;
    return 0;
};

//----------------------------------------------------------------------------
//                  test_function_linsolve_rev2
//----------------------------------------------------------------------------
Real test_function_linsolve_rev2::eval_mat(const Matrix& mat,bool , int code)
{
    (void)code;
    m_new_objects = 0;

    if (mat.is_square() == false)
        return 0.;

    Integer r = mat.rows();

    try
    {	
        Real dif		= 0;
        using container = dynamic_mat_set::container;
        
        long n1             = matcl::details::no_existing_objects();
        const container& mc = ms.get(r,K);
        n1                  = matcl::details::no_existing_objects() - n1;
        m_new_objects       = n1;

        size_t size = mc.size();

        for (size_t i = 0; i < size ; i ++ )
            dif += eval_mat_impl(mat,mc[i]);

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

Real test_function_linsolve_rev2::eval_mat_impl(const Matrix& A,const Matrix& B)
{
    Real out    = 0.;

    permvec p    = randperm(A.rows());
    permvec q    = randperm(A.cols());
    permvec e    = permvec::identity(A.rows());

    out         = out + test_notrans(A,B,e,e);
    out         = out + test_notrans(A,B,p,e);
    out         = out + test_notrans(A,B,e,q);
    out         = out + test_notrans(A,B,p,q);

    out         = out + test_trans(A,B,e,e);
    out         = out + test_trans(A,B,p,e);
    out         = out + test_trans(A,B,e,q);
    out         = out + test_trans(A,B,p,q);

    out         = out + test_ctrans(A,B,e,e);
    out         = out + test_ctrans(A,B,p,e);
    out         = out + test_ctrans(A,B,e,q);
    out         = out + test_ctrans(A,B,p,q);

    return out;
};

Real test_function_linsolve_rev2::test_notrans(const Matrix& A,const Matrix& B, permvec p, permvec q)
{
    try
    {	
        Matrix Xf       = linsolve_rev2(A, p, q, trans(B), trans_type::no_trans);

        check_struct(Xf);

        Real dif		= 0;              

        if (A.all_finite() == true && B.all_finite() == true)
        {
            // overflow is possible
            if (Xf.all_finite() == false)
            {
                Matrix smin     = svd1(A, svd_algorithm::dc)(end);
                Real eps        = constants::eps(Xf.get_value_code());

                if (smin < eps)
                    return 0.0;
                else
                    return 1.;
            };
        }
        else
        {
            if (Xf.all_finite() == true)
                return 1.;
        }

        Matrix test     = Xf * A(invperm(p), invperm(q)) - trans(B);

        dif             += norm_1(test);
        Real tol        = norm(A, 1) * norm(Xf,1) + norm(B, 1);
        Real eps        = constants::eps(Xf.get_value_code());

        if (dif < 100. * tol * eps)
            dif         = 0.;
                
        return dif;
    }
    catch(error::error_lsolve& ex)
    {
        if (A.is_square() && A.rows() == B.rows())
        {
            m_is_error = true;

            m_error = static_cast<error::matcl_exception_linalg&>(ex).what();
            return 1.;
        };

        return 0.;
    }
    catch(error::error_singular&)
    {
        Integer size    = min(A.rows(),A.cols());

        if (size == 0)
            return 0.;

        Matrix s        = svd1(A, svd_algorithm::dc);

        if (s(size) == 0)
            return 0.;

        Matrix cond     = s(1) / s(size);
        Real tol        = 1.0 / constants::eps(A.get_value_code());

        if (cond > 1.e-2 * tol)
            return 0;

        m_is_error = true;
        m_error = "singular matrix";
        return 1.;
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

Real test_function_linsolve_rev2::test_trans(const Matrix& A,const Matrix& B, permvec p, permvec q)
{
    try
    {	
        Matrix Xf       = linsolve_rev2(A,p,q, trans(B), trans_type::trans);
        check_struct(Xf);

        Real dif		= 0;              

        if (A.all_finite() == true && B.all_finite() == true)
        {
            // overflow is possible
            if (Xf.all_finite() == false)
            {
                Matrix smin     = svd1(A, svd_algorithm::dc)(end);
                Real eps        = constants::eps(Xf.get_value_code());

                if (smin < eps)
                    return 0.0;
                else
                    return 1.;
            };
        }
        else
        {
            if (Xf.all_finite() == true)
                return 1.;
        }

        Matrix test     = Xf * matcl::trans(A(invperm(p), invperm(q))) - trans(B);

        dif             += norm_1(test);
        Real tol        = norm(A, 1) * norm(Xf,1) + norm(B, 1);
        Real eps        = constants::eps(Xf.get_value_code());

        if (dif < 100. * tol * eps)
            dif         = 0.;
                
        return dif;
    }
    catch(error::error_lsolve& ex)
    {
        if (A.is_square() && A.rows() == B.rows())
        {
            m_is_error = true;

            m_error = static_cast<error::matcl_exception_linalg&>(ex).what();
            return 1.;
        };

        return 0.;
    }
    catch(error::error_singular&)
    {
        Integer size    = min(A.rows(),A.cols());

        if (size == 0)
            return 0.;

        Matrix s        = svd1(A, svd_algorithm::dc);

        if (s(size) == 0)
            return 0.;

        Matrix cond     = s(1) / s(size);
        Real tol        = 1.0 / constants::eps(A.get_value_code());

        if (cond > 1.e-2 * tol)
            return 0;

        m_is_error = true;
        m_error = "singular matrix";
        return 1.;
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

Real test_function_linsolve_rev2::test_ctrans(const Matrix& A,const Matrix& B, permvec p, permvec q)
{
    try
    {	
        Matrix Xf       = linsolve_rev2(A, p, q, trans(B), trans_type::conj_trans);
        check_struct(Xf);

        Real dif		= 0;              

        if (A.all_finite() == true && B.all_finite() == true)
        {
            // overflow is possible
            if (Xf.all_finite() == false)
            {
                Matrix smin     = svd1(A, svd_algorithm::dc)(end);
                Real eps        = constants::eps(Xf.get_value_code());

                if (smin < eps)
                    return 0.0;
                else
                    return 1.;
            };
        }
        else
        {
            if (Xf.all_finite() == true)
                return 1.;
        }

        Matrix test     = Xf * matcl::ctrans(A(invperm(p), invperm(q))) - trans(B);

        dif             += norm_1(test);
        Real tol        = norm(A, 1) * norm(Xf,1) + norm(B, 1);
        Real eps        = constants::eps(Xf.get_value_code());

        if (dif < 100. * tol * eps)
            dif         = 0.;
                
        return dif;
    }
    catch(error::error_lsolve& ex)
    {
        if (A.is_square() && A.rows() == B.rows())
        {
            m_is_error = true;

            m_error = static_cast<error::matcl_exception_linalg&>(ex).what();
            return 1.;
        };

        return 0.;
    }
    catch(error::error_singular&)
    {
        Integer size    = min(A.rows(),A.cols());

        if (size == 0)
            return 0.;

        Matrix s        = svd1(A, svd_algorithm::dc);

        if (s(size) == 0)
            return 0.;

        Matrix cond     = s(1) / s(size);
        Real tol        = 1.0 / constants::eps(A.get_value_code());

        if (cond > 1.e-2 * tol)
            return 0;

        m_is_error = true;
        m_error = "singular matrix";
        return 1.;
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

Real test_function_linsolve_rev2::eval_scalar(const Scalar& ,bool, int code )
{
    (void) code;
    return 0;
};

//----------------------------------------------------------------------------
//                  test_function_linsolve_rev2
//----------------------------------------------------------------------------

Real test_function_hess::eval_mat(const Matrix& mat,bool ,int code)
{
    (void)code;

    if (mat.is_square() == false)
        return 0.;

    // number x means result of tuple-x version
    unitary_matrix u2;
    Matrix h1;
    Matrix h2;

    try
    {
        h1          = hess(mat);
        tie(u2, h2) = hess2(mat);
    }
    catch(const error::error_hess_nonsq&)
    {
        m_is_error  = true;
        m_error     = "error_hess_nonsq should not appear here - matrix was square";
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
        check_struct(h1);
        check_struct(h2);

        if (mat.all_finite() == true)
        {
            // overflow is possible
            if (h1.all_finite() == false)
                return 1.;

            if (h2.all_finite() == false)
                return 1.;

            if (u2.all_finite() == false)
                return 1.;
        }
        else
        {
            if (h1.all_finite() == true)
                return 1.;

            if (h2.all_finite() == true)
                return 1.;

            if (u2.all_finite() == true)
                return 1.;
        }

        Real dif = 0;

        //check unitariness of u
        dif     += norm_1((ctrans(u2) * u2).to_matrix() - eye(u2.cols()));
        dif     += norm_1((u2 * ctrans(u2)).to_matrix() - eye(u2.cols()));

        value_code vc   = mat.get_value_code();
        Real tol_U      = 10.0 * constants::eps(vc) * u2.rows();

        if (dif < tol_U)
            dif         = 0.0;
        else
            return dif;

        // check correct factorizations
        Matrix mat_check1 = u2 * h1 * ctrans(u2);
        Matrix mat_check2 = u2 * h2 * ctrans(u2);

        dif     += norm_1(mat - mat_check1);
        dif     += norm_1(mat - mat_check2);

        if (dif < error_tolerance(100.0, mat) )
            dif = 0.0;
        else
            return dif;

        // check dimmensions
        if (h2.is_square() == false)
            dif += 1;

        if (h1.is_square() == false)
            dif += 1;

        if (has_struct_hessu(h1) == false)
            dif += 1;

        if (has_struct_hessu(h2) == false)
            dif += 1;

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

Real test_function_hess::eval_scalar(const Scalar&, bool,int)
{
    return 0;
};

//----------------------------------------------------------------------------
//                  test_function_schur
//----------------------------------------------------------------------------
Real test_function_schur::eval_mat(const Matrix& mat,bool ,int code)
{
    (void)code;

    if (mat.is_square() == false)
        return 0.;

    Matrix q, t;       
    Matrix q_compl, t_compl, temp_select;
    Matrix q_sel, t_sel;
    Matrix mat_sym;
    Matrix v_sym_rr, d_sym_rr;
    Matrix v_sym_qr, d_sym_qr;
    Matrix v_sym_dc, d_sym_dc;
    
    int     size    = mat.rows();

    try
    {
        tie(q,t)                = schur(mat);
        tie(q_compl,t_compl)    = schur_compl(mat);

        //prepare valid selection vector to maintain paired eigenvalues
        {
            temp_select         = mod(floor(real(get_diag(t)) * 123.), 2);

            if (temp_select.all_finite() == false)
                temp_select     = zeros(t.rows(), 1);
        }

        tie(q_sel, t_sel)       = ordschur(q, t, temp_select);
        
        // build matrix and call symmetric case
        mat_sym                 = mat + ctrans(mat);
        value_code vc           = mat_sym.get_value_code();
                    
        if (matrix_traits::is_real(vc) == true)
            mat_sym.set_struct(predefined_struct_type::sym);                        
        else
            mat_sym.set_struct(predefined_struct_type::her);

        tie(v_sym_qr, d_sym_qr)   = schur(mat_sym, schur_sym_alg::qr);
        tie(v_sym_dc, d_sym_dc)   = schur(mat_sym, schur_sym_alg::dc);
        tie(v_sym_rr, d_sym_rr)   = schur(mat_sym, schur_sym_alg::rrr);
    }
    catch(const error::error_size_eig&)
    {
        m_is_error  = true;
        m_error     = "error_size_eig should not appear here - matrix was square";
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
        check_struct(q);
        check_struct(t);
        check_struct(q_compl);
        check_struct(t_compl);
        check_struct(q_sel);
        check_struct(t_sel);
        check_struct(mat_sym);
        check_struct(v_sym_rr);
        check_struct(d_sym_rr);
        check_struct(v_sym_qr);
        check_struct(d_sym_qr);
        check_struct(v_sym_dc);
        check_struct(d_sym_dc);

        // check finitness
        if (mat.all_finite() == true)
        {
            if (q.all_finite() == false)
                return 1.;

            if (t.all_finite() == false)
                return 1.;

            if (q_compl.all_finite() == false)
                return 1.;

            if (t_compl.all_finite() == false)
                return 1.;

            if (q_sel.all_finite() == false)
                return 1.;

            if (t_sel.all_finite() == false)
                return 1.;

            if (mat_sym.all_finite() == false)
                return 1.;

            if (v_sym_rr.all_finite() == false)
                return 1.;

            if (d_sym_rr.all_finite() == false)
                return 1.;

            if (v_sym_qr.all_finite() == false)
                return 1.;

            if (d_sym_qr.all_finite() == false)
                return 1.;

            if (v_sym_dc.all_finite() == false)
                return 1.;

            if (d_sym_dc.all_finite() == false)
                return 1.;
        };

        if (mat.all_finite() == false)
        {
            if (q.all_finite() == true)
                return 1.;

            if (t.all_finite() == true)
                return 1.;

            if (q_compl.all_finite() == true)
                return 1.;

            if (t_compl.all_finite() == true)
                return 1.;

            if (q_sel.all_finite() == true)
                return 1.;

            if (t_sel.all_finite() == true)
                return 1.;

            if (mat_sym.all_finite() == true)
                return 1.;

            if (v_sym_rr.all_finite() == true)
                return 1.;

            if (d_sym_rr.all_finite() == true)
                return 1.;

            if (v_sym_qr.all_finite() == true)
                return 1.;

            if (d_sym_qr.all_finite() == true)
                return 1.;

            if (v_sym_dc.all_finite() == true)
                return 1.;

            if (d_sym_dc.all_finite() == true)
                return 1.;
        };

        // check dimmensions
        if (q.is_square() == false || t.is_square() == false 
            || q_compl.is_square() == false || t_compl.is_square() == false 
            || q_sel.is_square() == false   || t_sel.is_square() == false)
        {
            return 1.0;
        };

        // check if flags are set

        if (has_struct_unitary(v_sym_dc) == false)
            return 1.0;

        if (has_struct_unitary(v_sym_qr) == false)
            return 1.0;

        // rr algorithm does not set unitary flag
        //if (has_struct_unitary(v_sym_rr) == false)
        //    return 1.0;

        if (has_struct_unitary(q) == false)
            return 1.0;

        if (has_struct_unitary(q_compl) == false)
            return 1.0;

        if (has_struct_unitary(q_sel) == false)
            return 1.0;

        if (has_struct_qtriu(t) == false)
            return 1.0;

        if (has_struct_qtriu(t_sel) == false)
            return 1.0;

        if (has_struct_triu(t_compl) == false)
            return 1.0;

        if (has_struct_diag(d_sym_rr) == false)
            return 1.0;

        if (has_struct_diag(d_sym_qr) == false)
            return 1.0;

        if (has_struct_diag(d_sym_dc) == false)
            return 1.0;

        // check if symmetric alg was used in symmetric case
        if ( mat.get_struct().is_hermitian(mat.is_square(), is_real_matrix(mat)) == true)
        {
            if (has_struct_diag(t) == false)
                return 1.0;

            if (has_struct_diag(t_sel) == false)
                return 1.0;

            if (has_struct_diag(t_compl) == false)
                return 1.0;
        }

        Real dif = 0;

        //check unitary matrices
        
        dif     += norm_1(ctrans(v_sym_dc) * v_sym_dc - eye(size));
        dif     += norm_1(ctrans(v_sym_qr) * v_sym_qr - eye(size));
        dif     += norm_1(ctrans(v_sym_rr) * v_sym_rr - eye(size));
        dif     += norm_1(q * ctrans(q) - eye(size));
        dif     += norm_1(q_compl * ctrans(q_compl) - eye(size));
        dif     += norm_1(q_sel * ctrans(q_sel) - eye(size));

        value_code vc   = mat.get_value_code();
        Real tol_U      = 100.0 * constants::eps(vc) * mat.rows();

        if (dif < tol_U)
            dif         = 0.0;
        else
            return dif;

        //check if factorization is correct
        Matrix mat_check_sym_dc = v_sym_dc * d_sym_dc * ctrans(v_sym_dc);
        Matrix mat_check_sym_qr = v_sym_qr * d_sym_qr * ctrans(v_sym_qr);
        Matrix mat_check_sym_rr = v_sym_rr * d_sym_rr * ctrans(v_sym_rr);

        Matrix mat_check        = q * t * ctrans(q);    
        Matrix mat_check_compl  = q_compl * t_compl * ctrans(q_compl);    
        Matrix mat_check_sel    = q_sel * t_sel * ctrans(q_sel);   

        dif     += norm_1(mat_check_sym_dc - mat_sym);
        dif     += norm_1(mat_check_sym_qr - mat_sym);
        dif     += norm_1(mat_check_sym_rr - mat_sym);

        if (dif < error_tolerance(500.0, mat_sym) )
            dif = 0.0;
        else
            return dif;

        dif     += norm_1(mat - mat_check); 
        dif     += norm_1(mat - mat_check_compl);
        dif     += norm_1(mat - mat_check_sel);

        if (dif < error_tolerance(1000.0, mat) )
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

Real test_function_schur::eval_scalar(const Scalar&, bool,int)
{
    return 0;
};

//----------------------------------------------------------------------------
//                  test_function_eigs
//----------------------------------------------------------------------------
Real test_function_eigs::eval_mat(const Matrix& mat,bool ,int code)
{
    (void)code; 

    if (mat.is_square() == false)
        return 0.0;

    Matrix U;
    Matrix T;
    bool conv = false;
    
    // max_trials to accommodate non-determinism of ARPACK
    Integer max_trials = 5;

    Integer size    = mat.rows();
    Integer k       = std::max(0, std::min(5, size - 3));

    try
    {
        for (int trial = 0; trial < max_trials; trial++)
        {
            tie(U,T,conv) = pschur(mat, k, cluster_type::LM, 
                                   matcl::options{opt::speigs::return_nonconvergent(true)});
    
            if (conv == true)
                break;
        }
        
        if (conv == false)
        {
            tie(U,T,conv) = pschur(mat, k, cluster_type::LM, 
                                   matcl::options{opt::speigs::return_nonconvergent(false)});
        }
    }
    catch(const error::invalid_speigs_k&)
    {
        //problem is too small
        return 0.0;
    }
    catch(const error::error_size_eig& ex)
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
        if (mat.all_finite() == true)
        {
            if (U.all_finite() == false)
                return 1.0;

            if (T.all_finite() == false)
                return 1.0;
        };

        if (mat.all_finite() == false)
        {
            if (U.all_finite() == true)
                return 1.0;

            if (T.all_finite() == true)
                return 1.0;
        };

        if (has_struct_qtriu(T) == false)
            return 1.0;

        Real dif    = 0;

        // check orthogonality
        dif     += norm_1(ctrans(U) * U - eye(U.cols()));

        value_code vc   = mat.get_value_code();
        Real tol_U      = 1000.0 * constants::eps(vc) * U.cols();

        if (dif < tol_U)
            dif         = 0.0;
        else
            return dif;

        // check correctness
        dif         += norm_1(mat * U - U * T);

        if (dif < error_tolerance(1000.0, mat) )
            dif = 0.0;
        else
            return dif;
        
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

Real test_function_eigs::eval_scalar(const Scalar&, bool,int)
{
    return 0;
};

};};
