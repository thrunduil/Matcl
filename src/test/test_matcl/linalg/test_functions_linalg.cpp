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
#include "test_functions_linalg.h"
#include "test/test_matcl/framework/matrix_set/matrix_set_1.h"
#include "test/test_matcl/framework/matrix_utils.h"
#include "matcl-linalg/matcl_linalg.h"

#include "test_functions_linalg.h"
#include "test/test_matcl/framework/data_struct/test_selector.h"

#include <boost/thread.hpp>

namespace matcl { namespace test
{

//TODO

// return eps * norm_1(mat)
static Real epsilon_mat(const Matrix& mat)
{
    value_code vc  = mat.get_value_code();

    // add some small value when result is close to zero
    Real mr         = std::sqrt(constants::min_real(vc));

    Real ret        = constants::eps(vc) * norm_1(mat) + mr;
    return ret;
};

static Real error_tolerance(Real mult, const Matrix& mat)
{
    Real tol = mult * (min(mat.rows(), mat.cols())) * epsilon_mat(mat);

    return tol;
};

template<class V>
static Matrix make_scalar(const V& v, value_code vc)
{
    switch (vc)
    {
        case value_code::v_integer:
            return Matrix(convert_scalar<Integer>(v));
        case value_code::v_float:
            return Matrix(convert_scalar<Float>(v));
        case value_code::v_real:
            return Matrix(convert_scalar<Real>(v));
        case value_code::v_float_complex:
            return Matrix(convert_scalar<Float_complex>(v));
        case value_code::v_complex:
            return Matrix(convert_scalar<Complex>(v));
        case value_code::v_object:
            return Matrix(Object(v));
        default:
            throw std::runtime_error("invalid value code");
    };
}

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
            {
                int code = 279;

                Matrix mat = tf.get_matrix(code);

                matcl::disp(full(mat));

                Matrix S, est_rank, norm_E;
                permvec p;
                Matrix M;

                M           = hersum(mat);
                bool upper  = true;

                disp(M);
                std::tie(S,p,est_rank,norm_E) = cholmod(M,upper);

                check_struct(S);
                Real dif = 0;

                disp(S);

                // check factorization M(p,p) + E = S'S
                Matrix M1   = upper ? ctrans(S) * S : S * ctrans(S);
                int N       = M.rows();

                Matrix E = M(p, p) - M1;
                dif = norm_1(norm(E, 2.) - norm_E.get_scalar<Real>());
        
                if (dif < error_tolerance(10.0, M1))
                    dif = 0.;

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
            };

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
    
    //TODO

    //SELECT_TEST (3, test_svd());
    //SELECT_TEST (3, test_norm());
    //SELECT_TEST (3, test_cond());
    //SELECT_TEST (3, test_lu_rook());    
    //SELECT_TEST (3, test_lu_partial());   
    //SELECT_TEST (3, test_lu_complete());    
    //SELECT_TEST (3, test_chol());
    //SELECT_TEST (3, test_chol_rr());
    
    SELECT_TEST (3, test_cholmod());
    
    SELECT_TEST (3, test_linsolve());
    SELECT_TEST (3, test_linsolve_rev());
    SELECT_TEST (3, test_linsolve_rev2());        
    
    SELECT_TEST (3, test_qr());
    SELECT_TEST (3, test_hess());
    SELECT_TEST (3, test_schur());
    SELECT_TEST (3, test_eigs());
    SELECT_TEST (3, test_ldl());    
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

void linalg_functions_list::test_linsolve()
{
    Real out = 0.;

    {
        test_function_linsolve tf(0,ms);
        out += m_tests.make(&tf,m_options);
    }
    {
        test_function_linsolve tf(1,ms);
        out += m_tests.make(&tf,m_options);
    }
    {
        test_function_linsolve tf(3,ms);
        out += m_tests.make(&tf,m_options);
    }

    if (out == 0.)
        matcl::out_stream << std::string() + "linsolve: OK" + "\n";
    else
        matcl::out_stream << std::string() + "linsolve: FAILED"  + "\n";
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

void linalg_functions_list::test_linsolve_rev2()
{
    Real out = 0.;

    {
        test_function_linsolve_rev2 tf(0,ms);
        out += m_tests.make(&tf,m_options);
    }
    {
        test_function_linsolve_rev2 tf(1,ms);
        out += m_tests.make(&tf,m_options);
    }
    {
        test_function_linsolve_rev2 tf(3,ms);
        out += m_tests.make(&tf,m_options);
    }

    if (out == 0.)
        matcl::out_stream << std::string() + "linsolve_rev2: OK" + "\n";
    else
        matcl::out_stream << std::string() + "linsolve_rev2: FAILED"  + "\n";
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
//                  test_function_linsolve
//----------------------------------------------------------------------------
Real test_function_linsolve::eval_mat(const Matrix& mat,bool , int code)
{
    (void)code;
    m_new_objects = 0;

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

//TODO: rewrite from this point to the end


Real test_function_linsolve::test_notrans(const Matrix& A,const Matrix& B, permvec p, permvec q)
{
    Matrix X, X_f;
    try
    {		
        X       = linsolve(A,p,q, B);
    }
    catch(error::error_lsolve&)
    {
        if (A.is_square() && A.cols() == B.cols())
            throw;

        return 0.;
    }
    catch(error::error_singular&)
    {
        Integer size    = min(A.rows(),A.cols());
        Matrix s        = svd1(A,svd_algorithm::dc);

        if (size == 0 || !(s(1) / s(size) < 1e+14))
        {
            return 0;
        }
        if (norm(B,1) == 0)
        {
            return 0;
        };
        throw;
    }
    catch(std::exception&)
    {
        throw;
    };

    try
    {		        
        Real dif		= 0;
        check_struct(X);      

        dif             += norm_1(A(invperm(p).to_matrix(), invperm(q).to_matrix())*X - B);
        Real n          = 1000*(norm(A)*norm(X) + norm(B) + constants::eps());

        if (is_nan(n) == false && is_inf(n) == false)
        {
            dif             = dif/n;
        }
        else
        {
            dif         = 0;
        };
        if (dif < constants::eps())
        {
            dif = 0;
        };
        return dif;
    }
    catch(error::error_singular&)
    {
        if (norm(B,1) == 0)
        {
            return 0;
        };
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
        Matrix X    = linsolve(full(A),B, trans_type::trans);
    }
    catch(error::error_lsolve&)
    {
        return 0.;
    }
    catch(error::error_singular&)
    {
        return 0.;
    }
    catch(std::exception&)
    {
        throw;
    };

    try
    {	
        Matrix Xf       = linsolve(A,p,q, B,trans_type::trans);

        Real dif		= 0;
        check_struct(Xf);
       
        dif             += norm_1(matcl::trans(A(invperm(p).to_matrix(), invperm(q).to_matrix()))*Xf - B);
        Real n          = 1000*(norm(A)*norm(Xf) + norm(B) + constants::eps());

        if (is_nan(n) == false && is_inf(n) == false)
        {
            dif             = dif/n;
        }
        else
        {
            dif         = 0;
        };
        if (dif < constants::eps())
        {
            dif = 0;
        };
        return dif;
    }
    catch(error::error_singular&)
    {
        return 0;
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
        Matrix X    = linsolve(full(A),B, trans_type::conj_trans);
    }
    catch(error::error_lsolve&)
    {
        return 0.;
    }
    catch(error::error_singular&)
    {
        return 0.;
    }
    catch(std::exception&)
    {
        throw;
    };

    try
    {	
        Matrix Xf       = linsolve(A,p,q, B,trans_type::conj_trans);

        Real dif		= 0;
        check_struct(Xf);
       
        dif             += norm_1(matcl::ctrans(A(invperm(p).to_matrix(),invperm(q).to_matrix()))*Xf - B);
        Real n          = 1000*(norm(A)*norm(Xf) + norm(B) + constants::eps());

        if (is_nan(n) == false && is_inf(n) == false)
        {
            dif             = dif/n;
        }
        else
        {
            dif         = 0;
        };
        if (dif < constants::eps())
        {
            dif = 0;
        };
        return dif;
    }
    catch(error::error_singular&)
    {
        return 0.;
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

Real test_function_linsolve_rev::eval_mat(const Matrix& mat,bool , int code)
{
    (void)code;
    m_new_objects = 0;

    Integer r = mat.rows();

    try
    {	
        Real dif		= 0;
        using container = dynamic_mat_set::container;
        
        long n1 = matcl::details::no_existing_objects();
        const container& mc = ms.get(r,K);
        n1 = matcl::details::no_existing_objects() - n1;
        m_new_objects = n1;

        size_t size = mc.size();

        for (size_t i = 0; i < size ; i ++ )
        {
            dif += eval_mat_impl(mat,mc[i]);
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
        Matrix X    = linsolve(A, B, trans_type::trans);
    }
    catch(error::error_lsolve&)
    {
        return 0.;
    }
    catch(error::error_singular&)
    {
        return 0.;
    }
    catch(std::exception&)
    {
        throw;
    };

    try
    {	
        Matrix Xf       = linsolve_rev(A,p,q, trans(B));

        Real dif		= 0;
        check_struct(Xf);       

        dif             += norm_1(Xf*A (invperm(p).to_matrix(),invperm(q).to_matrix())- trans(B));
        Real n          = 1000*(norm(A)*norm(Xf) + norm(B) + constants::eps());

        if (is_nan(n) == false && is_inf(n) == false)
        {
            dif             = dif/n;
        }
        else
        {
            dif         = 0;
        };
        if (dif < constants::eps())
        {
            dif = 0;
        };
        return dif;
    }
    catch(error::error_singular&)
    {
        return 0;
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

    Integer r = mat.rows();

    try
    {	
        Real dif		= 0;
        using container = dynamic_mat_set::container;
        
        long n1 = matcl::details::no_existing_objects();
        const container& mc = ms.get(r,K);
        n1 = matcl::details::no_existing_objects() - n1;
        m_new_objects = n1;

        size_t size = mc.size();

        for (size_t i = 0; i < size ; i ++ )
        {
            dif += eval_mat_impl(mat,mc[i]);
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
        Matrix X    = linsolve(A, B, trans_type::trans);
    }
    catch(error::error_lsolve&)
    {
        return 0.;
    }
    catch(error::error_singular&)
    {
        return 0.;
    }
    catch(std::exception&)
    {
        throw;
    };

    try
    {	
        Matrix Xf       = linsolve_rev2(A,p,q, trans(B), trans_type::no_trans);

        Real dif		= 0;
        check_struct(Xf);      

        dif             += norm_1(Xf*A(invperm(p).to_matrix(),invperm(q).to_matrix()) - trans(B));
        Real n          = 1000*(norm(A)*norm(Xf) + norm(B) + constants::eps());

        if (is_nan(n) == false && is_inf(n) == false)
        {
            dif             = dif/n;
        }
        else
        {
            dif         = 0;
        };
        if (dif < constants::eps())
        {
            dif = 0;
        };
        return dif;
    }
    catch(error::error_singular&)
    {
        return 0;
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
        Matrix X    = linsolve(A, B, trans_type::no_trans);
    }
    catch(error::error_lsolve&)
    {
        return 0.;
    }
    catch(error::error_singular&)
    {
        return 0.;
    }
    catch(std::exception&)
    {
        throw;
    };

    try
    {	
        Matrix Xf       = linsolve_rev2(A,p,q, trans(B), trans_type::trans);

        Real dif		= 0;
        check_struct(Xf);

        dif             += norm_1(Xf*matcl::trans(A(invperm(p).to_matrix(), invperm(q).to_matrix())) - trans(B));
        Real n          = 1000*(norm(A)*norm(Xf) + norm(B) + constants::eps());

        if (is_nan(n) == false && is_inf(n) == false)
        {
            dif             = dif/n;
        }
        else
        {
            dif         = 0;
        };
        if (dif < constants::eps())
        {
            dif = 0;
        };
        return dif;
    }
    catch(error::error_singular&)
    {
        return 0;
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
        Matrix X    = linsolve(A, B, trans_type::no_trans);
    }
    catch(error::error_lsolve&)
    {
        return 0.;
    }
    catch(error::error_singular&)
    {
        return 0.;
    }
    catch(std::exception&)
    {
        throw;
    };

    try
    {	
        Matrix Xf       = linsolve_rev2(A,p,q, trans(B), trans_type::conj_trans);

        Real dif		= 0;
        check_struct(Xf);
       
        dif             += norm_1(Xf*matcl::ctrans(A(invperm(p).to_matrix(),invperm(q).to_matrix())) - trans(B));
        Real n          = 1000*(norm(A)*norm(Xf) + norm(B) + constants::eps());

        if (is_nan(n) == false && is_inf(n) == false)
        {
            dif             = dif/n;
        }
        else
        {
            dif         = 0;
        };
        if (dif < constants::eps())
        {
            dif = 0;
        };
        return dif;
    }
    catch(error::error_singular&)
    {
        return 0;
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
        r1_internal = qr_internal(mat);
        r1 = qr(mat);
        tie (q2, r2) = qr2(mat);
        tie (q3, r3, e3) = qr3(mat);
        tie (q2eco, r2eco) = qr2(mat, true);
        tie (q3eco, r3eco, e3eco) = qr3 (mat, true);
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
        if (any(any(neg(is_finite(
            r1_internal + r1 + q2 * r2 + q3 * r3 + q2eco * r2eco + q3eco * r3eco)),1),2))
        {
            if (all(all(is_finite(mat),1),2))
            {
                dif += 1000;
            }
        }

        Matrix Q2   = q2.to_matrix();
        Matrix Q2_e = q2eco.to_matrix();
        Matrix Q3   = q3.to_matrix();
        Matrix Q3_e = q3eco.to_matrix();

        //check unitariness of q
        dif +=norm_1(ctrans(Q2) * Q2 - eye(q2.cols()));
        dif +=norm_1(ctrans(Q3) * Q3 - eye(q3.cols()));
        dif +=norm_1(ctrans(Q2_e) * Q2_e - eye(q2eco.cols()));
        dif +=norm_1(ctrans(Q3_e) * Q3_e - eye(q3eco.cols()));
        // check correct factorizations
        Matrix mat_check1_internal = q2 * triu(r1_internal);
        Matrix mat_check1 = q2 * r1;
        Matrix mat_check2 = q2 * r2;
        Matrix mat_check3p = q3 * r3;
        Matrix mat_check2eco = q2eco * r2eco;
        Matrix mat_check3ecop = q3eco * r3eco;
        if(any(any(is_nan(mat) + is_inf(mat),1),2))
        {
            if (! (any(any(is_nan(r1_internal),1),2)) && 
                  (any(any(is_nan(r1),1),2)) && 
                  (any(any(is_nan(Q2),1),2)) && 
                  (any(any(is_nan(r2),1),2)) && 
                  (any(any(is_nan(Q3),1),2)) && 
                  (any(any(is_nan(r3),1),2)) &&  
                  (any(any(is_nan(Q2_e),1),2)) &&  
                  (any(any(is_nan(r2eco),1),2)) &&  
                  (any(any(is_nan(Q3_e),1),2)) &&  
                  (any(any(is_nan(r3eco),1),2))
               )
            {
                dif += 4321654783;
            }
        }
        else
        {
            dif += norm_1(mat - mat_check1_internal) / 
                (norm_1(mat_check1_internal) + constants::eps());
            dif += norm_1(mat - mat_check1) / (norm_1(mat_check1) + constants::eps());
            dif += norm_1(mat - mat_check2) / (norm_1(mat_check2) + constants::eps());
            dif += norm_1(mat(colon(),e3) - mat_check3p) / (norm_1(mat_check3p) + constants::eps());
            dif += norm_1(mat - mat_check2eco) / (norm_1(mat_check2eco) + constants::eps());
            dif += norm_1(mat(colon(),e3eco) - mat_check3ecop) / (norm_1(mat_check3ecop) + constants::eps());
        }
        // check uppertriangular r1 r2 r3 (r1_internal is NOT uppertriangular)
        dif += norm_1(tril(r1, -1));
        dif += norm_1(tril(r2, -1));
        dif += norm_1(tril(r3, -1));
        dif += norm_1(tril(r2eco, -1));
        dif += norm_1(tril(r3eco, -1));
        // check dimmensions
        if (q2.is_square() == false)
        {
            dif += 1;
        };
        if (q3.is_square() == false)
        {
            dif += 1;
        };
        if (q2eco.cols() != min(mat.cols(), mat.rows()))
        {
            dif+=1;
        };
        if (q3eco.cols() != min(mat.cols(), mat.rows()))
        {
            dif+=1;
        };
        if (r2eco.rows() != min(mat.cols(), mat.rows()))
        {
            dif+=1;
        };
        if (r3eco.rows() != min(mat.cols(), mat.rows()))
        {
            dif+=1;
        };
        if(mat.cols() > 1 && mat.structural_nnz() > 0 &&
            all(all(is_finite(mat), 1), 2)
            && (  !r1.get_struct().is_triu() 
             ||   !r2.get_struct().is_triu()           || !is_unitary(Q2.get_struct())
             ||   !r3.get_struct().is_triu()           || !is_unitary(Q3.get_struct())
             ||   !r2eco.get_struct().is_triu()        
             ||   !r3eco.get_struct().is_triu()        
              ) )
        {
            dif += 43247895;
        }

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

        if (abs(dif) < 1e-13) // TODO: was 1e-14 before - make sure it is OK
        {
            dif = 0.;
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
Real test_function_qr::eval_scalar(const Scalar&, bool,int)
{
    return 0;
};

Real test_function_hess::eval_mat(const Matrix& mat,bool ,int code)
{
    (void)code;
    
    // number x means result of tuple-x version
    unitary_matrix u2;
    Matrix h1;
    Matrix h2;

    try
    {
        h1 = hess(mat);
        tie (u2, h2) = hess2(mat);
    }
    catch(const error::error_hess_nonsq&)
    {
        m_is_error = true;
        if (!mat.is_square())
        {
            m_is_error = false;
            return 0.;
        }
        m_error = "error_hess_nonsq should not appear here - matrix was square";
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

    try
    {
        Real dif = 0;
        // check finitenes of all results
        if (any(any(neg(is_finite(h1 + h2)),1),2))
        {
            if (all(all(is_finite(mat),1),2))
            {
                dif += 1000;
            }
        }
        //check unitariness of u
        dif +=norm_1((ctrans(u2) * u2).to_matrix() - eye(u2.cols()));
        dif +=norm_1((u2 * ctrans(u2)).to_matrix() - eye(u2.cols()));
        // check correct factorizations
        Matrix mat_check1 = u2*h1*ctrans(u2);
        Matrix mat_check2 = u2*h2*ctrans(u2);
        dif += norm_1(mat - mat_check1) / (norm_1(mat_check1) + constants::eps());
        dif += norm_1(mat - mat_check2) / (norm_1(mat_check2) + constants::eps());
        // check hessenberg h2 h1
        dif += norm_1(tril(h1,-2));
        dif += norm_1(tril(h2,-2));
        // check dimmensions
        if (h2.is_square() == false)
        {
            dif += 1;
        };
        if (h1.is_square() == false)
        {
            dif += 1;
        };
        if(mat.numel() > 1 && mat.structural_nnz() > 0 &&
            all(all(is_finite(mat), 1), 2))
        {
            dif += 43247895;
        }
        check_struct(h1);
        check_struct(h2);
        
        if (abs(dif) < 1e-14)
        {
            dif = 0.;
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
Real test_function_hess::eval_scalar(const Scalar&, bool,int)
{
    return 0;
};

Real test_function_schur::eval_mat(const Matrix& mat,bool ,int code)
{
    (void)code;
    Matrix  q;        Matrix  t;       
    Matrix  q_compl;  Matrix  t_compl;       Matrix temp_select;
    Matrix  q_sel;    Matrix  t_sel;
    Matrix  v_sym;    Matrix  d_sym;         Matrix  mat_sym;
    Matrix  v_sym_rr; Matrix  d_sym_rr;
    Matrix  v_sym_qr; Matrix  d_sym_qr;
    Matrix  v_sym_dc; Matrix  d_sym_dc;
    int     size    = mat.rows();
    try
    {
        tie(q,t)            =   schur(mat);
        tie(q_compl,t_compl)=   schur_compl(mat);
        //prepare valid selection vector to maintain paired eigenvalues
        {
            temp_select         =   mod(floor(real(get_diag(t)) * 123.), 2);
            if ((bool)all(is_finite(temp_select),1) == false)
            {
                temp_select     =   zeros(t.rows(), 1);
            }
        }
        tie(q_sel, t_sel)   =   ordschur(q, t, temp_select);
        // build matrix and call symmetric case
        mat_sym = mat + ctrans(mat);
        if (mat_sym.get_value_code() == value_code::v_complex)
        {
            mat_sym.set_struct(predefined_struct_type::her);
        }
        else
        {
            mat_sym.set_struct(predefined_struct_type::sym);
        }
        tie(v_sym_qr, d_sym_qr)     =   schur(mat_sym, schur_sym_alg::qr);
        tie(v_sym_dc, d_sym_dc)     =   schur(mat_sym, schur_sym_alg::dc);
        tie(v_sym_rr, d_sym_rr)     =   schur(mat_sym, schur_sym_alg::rrr);
    }
    catch(const error::error_size_eig&)
    {
        if (!mat.is_square())
        {
            return 0.;
        }
        else
        {
            m_is_error = true;
            m_error = "error_size_eig should not appear here - matrix was square";
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
        // check symmetric case
        //check unitary v
        dif += norm_1(ctrans(v_sym) * v_sym - eye(size));
        dif += norm_1(ctrans(v_sym_dc) * v_sym_dc - eye(size));
        dif += norm_1(ctrans(v_sym_qr) * v_sym_qr - eye(size));
        dif += norm_1(ctrans(v_sym_rr) * v_sym_rr - eye(size));

        //check correct factorization
        Matrix mat_check_sym    = v_sym * d_sym * ctrans(v_sym);
        Matrix mat_check_sym_dc = v_sym_dc * d_sym_dc * ctrans(v_sym_dc);
        Matrix mat_check_sym_qr = v_sym_qr * d_sym_qr * ctrans(v_sym_qr);
        Matrix mat_check_sym_rr = v_sym_rr * d_sym_rr * ctrans(v_sym_rr);

        dif += norm_1(mat_check_sym - mat_sym);
        dif += norm_1(mat_check_sym_dc - mat_sym);
        dif += norm_1(mat_check_sym_qr - mat_sym);
        dif += norm_1(mat_check_sym_rr - mat_sym);

        //roundoff error of testing by multiplication into account
        if (dif < 10 * ::powl(size, 3./2.) * constants::eps() * (norm_1(mat_sym) + constants::eps()))
        {
            dif = 0;
        }
        //check diagonal
        if (all(all(is_finite(d_sym), 1), 2) 
            && size > 1 && (!d_sym.get_struct().is_diag() || !d_sym_dc.get_struct().is_diag()
                            || !d_sym_qr.get_struct().is_diag() || !d_sym_rr.get_struct().is_diag()))
        {
            dif += 1000000;
        }
        check_struct(d_sym);
        check_struct(d_sym_qr);
        check_struct(d_sym_dc);
        check_struct(d_sym_rr);

        // check general case
        //check if any selection was done
        if ((sum(temp_select,1) > 0)
            && (temp_select(1,1) == 0)
            && all(all(is_finite(q_sel + t_sel), 1), 2))
        {
            dif += (norm_1(q - q_sel) == 0);
            dif += (norm_1(t - t_sel) == 0);
        }
        //check unitariness of q's
        dif +=norm_1(q * ctrans(q) - eye(size));
        dif +=norm_1(q_compl * ctrans(q_compl) - eye(size));
        dif +=norm_1(q_sel * ctrans(q_sel) - eye(size));
        // check correct factorization
        Matrix mat_check = q * t * ctrans(q);    
        Matrix mat_check_compl = q_compl * t_compl * ctrans(q_compl);    
        Matrix mat_check_sel = q_sel * t_sel * ctrans(q_sel);   
        dif += norm_1(mat - mat_check); 
        dif += norm_1(mat - mat_check_compl);
        dif += norm_1(mat - mat_check_sel);
        //take roundoff error of testing by multiplication into account
        if (dif < 9 *
            ::powl(size, 3./2.) * constants::eps()
            * (norm_1(mat) + constants::eps()))
        {
            dif = 0;
        }
        // check finitenes of all results
        if (any(any(neg(is_finite(q + t + q_compl + t_compl + 
                                  q_sel + t_sel + v_sym + d_sym + d_sym_rr + d_sym_dc + d_sym_qr)),1),2))
        {
            if (all(all(is_finite(mat),1),2))
            {
                dif += 1000;
            }
        }
        // check qtriu t's
        dif += norm_1(tril(t,-2));
        dif += norm_1(tril(t_sel,-2));
        // check triu of explicit complex Schur factorization
        dif += norm_1(tril(t_compl,-1));
        // check dimmensions
        if (q.is_square() == false       || t.is_square() == false || 
            q_compl.is_square() == false || t_compl.is_square() == false || 
            q_sel.is_square() == false   || t_sel.is_square() == false)
        {
            dif += 1;
        };
        // check if symmetric alg was used in symmetric case
        if (   (mat.get_struct().is_hermitian(mat.is_square(), is_real_matrix(mat))) )
        {
            if (size > 1 && !t.get_struct().is_diag() && !t_sel.get_struct().is_diag())
            {
                dif += 10000;
            }
        }
        check_struct(q);
        check_struct(t);
        check_struct(q_compl);
        check_struct(t_compl);
        check_struct(q_sel);
        check_struct(t_sel);
        check_struct(v_sym);
        check_struct(v_sym_dc);
        check_struct(v_sym_qr);
        check_struct(v_sym_rr);
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
Real test_function_eigs::eval_mat(const Matrix& mat,bool ,int code)
{
    (void)code; 
    Matrix U;
    Matrix T;
    bool conv = false;
    // max_trials to accommodate non-determinism of ARPACK
    Integer max_trials = 5;

    Integer size = mat.rows();
    Integer k = std::max(0, std::min(5, size - 3));
    try
    {
        for (int trial = 0; trial < max_trials; trial++)
        {
            tie(U,T,conv) = pschur(mat, k, cluster_type::LM, matcl::options{opt::speigs::return_nonconvergent(true)});
            if (conv == true)
                break;
        }
        
        if (conv == false)
        {
            tie(U,T,conv) = pschur(mat, k, cluster_type::LM, matcl::options{opt::speigs::return_nonconvergent(false)});
        }
    }
    catch(const error::error_size_eig& )
    {
        if (mat.is_square()) throw;
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
        dif += norm_1(mat * U - U * T);
        //roundoff error of testing by multiplication into account
        if (dif < 3 * ::powl(size, 3./2.) * constants::eps() * (norm_1(mat) + constants::eps()))
        {
            dif = 0;
        }
        // check finiteness of results
        if (any(any(neg(is_finite(U)),1),2) ||
            any(any(neg(is_finite(T)),1),2))
        {
            if (   all(all(is_finite(mat),1),2))
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
};
Real test_function_eigs::eval_scalar(const Scalar&, bool,int)
{
    return 0;
};
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
            if (mat.rows() == mat.cols())
            {
                mat     = mat_temp + trans(mat_temp);
                matherm = mat_temp + ctrans(mat_temp);
            }
            else
            {
                mat     = mat_temp * trans(mat_temp);
                matherm = mat_temp * ctrans(mat_temp);
            };
            mat.set_struct(predefined_struct_type::sym);
            matherm.set_struct(predefined_struct_type::her);
        }
    }
    catch(const error::invalid_eeop&)
    {
        return 0.;
    }

    int size = mat.rows();
    Matrix l,d;
    Matrix lh,dh;
    permvec p, ph;

    try
    {
        tie(l,d,p)            =   ldl(mat,upper);
        tie(lh,dh,ph)         =   ldl_herm(matherm,upper);
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
        if (upper == false)
        {
            dif += norm_1(triu(l,1));
            dif += norm_1(triu(lh,1));
        }
        else
        {
            dif += norm_1(tril(l,1));
            dif += norm_1(tril(lh,1));
        };
        dif += norm_1(diag(get_diag(l)) - speye(size,size));
        dif += norm_1(diag(get_diag(lh)) - speye(size,size));
        
        // check correct factorization   
        Matrix mat_check_tr;
        Matrix mat_checkh_tr;

        mat_check_tr    = l * d * trans(l);  
        mat_checkh_tr   = lh * dh * ctrans(lh);

        dif += norm_1(mat(p,p) - mat_check_tr);
        dif += norm_1(matherm(p,p) - mat_checkh_tr);
        
        //take roundoff error of testing by multiplication into account
        if (dif < 8 *
            (::powl(size, 3./2.) + size) * constants::eps()
            * (::pow(norm_1(l), 2.) * norm_1(d) + constants::eps()))
        {
            dif = 0;
        }
        
        // check finitenes of all results
        if (any(any(neg(is_finite(l + d + lh + dh)),1),2))
        {
            if (all(all(is_finite(mat),1),2))
            {
                dif += 1000;
            }
        }
        
        // check dimmensions
        if (d.is_square() == false  || l.is_square() == false
         || dh.is_square() == false || lh.is_square() == false)
        {
            dif += 1;
        };
        
        //check block diagonal d (coarse check, TODO: improve)
        //will work only with present inexact qtriu/l check struct
        if(all(all(is_finite(mat),1),2))
        {
            Matrix block_diag_help = d.clone();
            d.set_struct(predefined_struct_type::qtril);
            check_struct(d);
            d.set_struct(predefined_struct_type::qtriu);
            check_struct(d);
        }
        check_struct(l);
        check_struct(d);
        check_struct(lh);
        check_struct(dh);

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

};};
