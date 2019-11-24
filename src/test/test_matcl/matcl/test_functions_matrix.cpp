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
#include "test_functions_matrix.h"

#include "test/test_matcl/framework/matrix_set/matrix_set_1.h"
#include "matcl-core/IO/logger.h"
#include "test/test_matcl/framework/data_struct/test_selector.h"

#include <boost/thread.hpp>

namespace matcl { namespace test
{

namespace md = matcl::details;

class test_matrix
{
    matrix_functions_list&		tf;
    const test::options&		opts;
    Integer                     thread_id;

    public:
        test_matrix(matrix_functions_list& tf, const test::options& opts, Integer id)
            :tf(tf),opts(opts), thread_id(id)
        {};

        test_matrix(const test_matrix& tu)
            :tf(tu.tf),opts(tu.opts), thread_id(tu.thread_id)
        {};

        void make()
        {   
            /*
            Integer code = 151;
            Matrix mat = tf.get_matrix(code);
            matcl::disp(mat);
            */

            tf.make(opts);
        };

        void operator()()
        {
            make();
        };

    private:		
        test_matrix& operator=(const test_matrix&) = delete;
};

void test_matrix_st(const rand_matrix_ptr& rand)
{
    test::options opts;

    try
    {
        opts.show_partial_res   = 0;
        opts.first_matrix_code  = 0;
        opts.show_memleaks      = true;

        test::mat_set_1 ms1(rand);
        matrix_functions_list tf(ms1,rand);
        
        test_matrix tu(tf,opts,0);
        tu.make();
    }
    catch(const std::exception& ex)
    {
        matcl::out_stream <<ex.what();
    };
};

void test_matrix_mt(const rand_matrix_ptr& rand)
{
    test::options opts;
    try
    {
        opts.show_partial_res   = 0;
        opts.first_matrix_code  = 0;
        opts.show_memleaks      = false;

        test::mat_set_2 ms1(rand);
        matrix_functions_list tf(ms1,rand);

        boost::thread_group tg;

        for (int i = 0; i < 10; i++)
            tg.create_thread(test_matrix(tf,opts,i));

        tg.join_all();
    }
    catch(const std::exception& ex)
    {
        matcl::out_stream <<ex.what();
    };
};

matrix_functions_list::matrix_functions_list(const matrix_set& ms,rand_matrix_ptr rand)
:m_tests(ms), m_rand(rand)
{};

void matrix_functions_list::make(options opts)
{
    m_options = opts;
    
    SELECT_TEST (3, test_mat_compile());
    SELECT_TEST (3, test_resize_reserve());
    SELECT_TEST (3, test_resize_reserve_b());    

    SELECT_TEST (3, test_get_scalar());
    SELECT_TEST (3, test_get_scalar_unique());
    
    SELECT_TEST (3, test_get_array());
    SELECT_TEST (3, test_get_const_array());
    SELECT_TEST (3, test_get_impl_unique());    
};

Matrix matrix_functions_list::get_matrix(int code) const
{
    return m_tests.get_matrix(code);
};


void matrix_functions_list::test_mat_compile()
{
    Real out = 0.;
    test_function_mat_compile tf;
    out += m_tests.make(&tf,m_options);

    if (out == 0.)
        matcl::out_stream << std::string() +   "mat_compile: OK" + "\n";
    else
        matcl::out_stream << std::string() +   "mat_compile: FAILED"  + "\n";
};

void matrix_functions_list::test_resize_reserve()
{
    Real out = 0.;
    test_function_resize_reserve tf;
    out += m_tests.make(&tf,m_options);

    if (out == 0.)
        matcl::out_stream << std::string() +   "resize_reserve: OK" + "\n";
    else
        matcl::out_stream << std::string() +   "resize_reserve: FAILED"  + "\n";
};

void matrix_functions_list::test_resize_reserve_b()
{
    Real out = 0.;
    test_function_resize_reserve_b tf;
    out += m_tests.make(&tf,m_options);

    if (out == 0.)
        matcl::out_stream << std::string() +   "resize_reserve_band: OK" + "\n";
    else
        matcl::out_stream << std::string() +   "resize_reserve_band: FAILED"  + "\n";
};

void matrix_functions_list::test_get_scalar()
{
    Real out = 0.;
    test_function_get_scal tf;
    out += m_tests.make(&tf,m_options);

    if (out == 0.)
        matcl::out_stream << std::string() +   "get_scal: OK" + "\n";
    else
        matcl::out_stream << std::string() +   "get_scal: FAILED"  + "\n";
};

void matrix_functions_list::test_get_scalar_unique()
{
    Real out = 0.;
    test_function_get_scal_unique tf;
    out += m_tests.make(&tf,m_options);

    if (out == 0.)
        matcl::out_stream << std::string() +   "get_scal_unique: OK" + "\n";
    else
        matcl::out_stream << std::string() +   "get_scal_unique: FAILED"  + "\n";
};

void matrix_functions_list::test_get_array()
{
    Real out = 0.;
    test_function_get_array tf;
    out += m_tests.make(&tf,m_options);

    if (out == 0.)
        matcl::out_stream << std::string() +   "get_array: OK" + "\n";
    else
        matcl::out_stream << std::string() +   "get_array: FAILED"  + "\n";
};

void matrix_functions_list::test_get_const_array()
{
    Real out = 0.;
    test_function_get_const_array tf;
    out += m_tests.make(&tf,m_options);

    if (out == 0.)
        matcl::out_stream << std::string() +   "get_const_array: OK" + "\n";
    else
        matcl::out_stream << std::string() +   "get_const_array: FAILED"  + "\n";
};

void matrix_functions_list::test_get_impl_unique()
{
    Real out = 0.;
    test_function_get_impl_unique tf;
    out += m_tests.make(&tf,m_options);

    if (out == 0.)
        matcl::out_stream << std::string() +   "get_impl_unique: OK" + "\n";
    else
        matcl::out_stream << std::string() +   "get_impl_unique: FAILED"  + "\n";
};

//--------------------------------------------------------------------

Real test_function_get_array::eval_mat(const Matrix& mat,bool,int code )
{
    (void)code;
    try
    {
        Matrix A = full(mat);
        matcl::mat_code	mt = A.get_matrix_code();

        if (A.rows() == 0 || A.cols() == 0)
            return 0;

        Matrix A2           = A;

        {
            auto warn = error::enable_warnings(false);
            {
                Matrix A3 = A2;
                A3.get_array_unique<Integer>();
            }
            {
                Matrix A3 = A2;
                A3.get_array_unique<Real>();
            }
            {
                Matrix A3 = A2;
                A3.get_array_unique<Float>();
            }
            {
                Matrix A3 = A2;
                A3.get_array_unique<Complex>();
            }
            {
                Matrix A3 = A2;
                A3.get_array_unique<Float_complex>();
            }
        };

        switch (mt)
        {
            case mat_code::integer_dense:
            {
                using value_type    = Integer;
                value_type* ptr     = A.get_array_unique<value_type>();

                if (!ptr)
                    return 0;

                Real dif    = norm_1(A(1)-ptr[0]);
                ptr[0]      = 1;
                dif         += norm_1(A(1)-1);

                return dif; 
            }
            case mat_code::real_dense:
            {
                using value_type    = Real;
                value_type* ptr     = A.get_array_unique<value_type>();

                if (!ptr)
                    return 0;

                Real dif    = norm_1(A(1)-ptr[0]);
                ptr[0]      = 1;
                dif         += norm_1(A(1)-1);

                return dif;
            }
            case mat_code::float_dense:
            {
                using value_type    = Float;
                value_type* ptr     = A.get_array_unique<value_type>();

                if (!ptr)
                    return 0;

                Real dif    = norm_1(A(1)-ptr[0]);
                ptr[0]      = 1;
                dif         += norm_1(A(1)-1);

                return dif;
            }
            case mat_code::complex_dense:
            {
                using value_type    = Complex;
                value_type* ptr     = A.get_array_unique<value_type>();

                if (!ptr)
                    return 0;

                Real dif    = norm_1(A(1)-ptr[0]);
                ptr[0]      = 1;
                dif         += norm_1(A(1)-1);

                return dif;
            };
            case mat_code::float_complex_dense:
            {
                using value_type    = Float_complex;
                value_type* ptr     = A.get_array_unique<value_type>();

                if (!ptr)
                    return 0;

                Real dif    = norm_1(A(1)-ptr[0]);
                ptr[0]      = 1;
                dif         += norm_1(A(1)-1);

                return dif;
            };
        };

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
Real test_function_get_array::eval_scalar(const Scalar& ,bool,int code )
{
    (void)code;
    return 0.;
};

Real test_function_get_const_array::eval_mat(const Matrix& mat,bool,int code )
{
    (void)code;
    try
    {
        Matrix A            = full(mat);
        matcl::mat_code	mt  = A.get_matrix_code();

        if (A.rows() == 0 || A.cols() == 0)
            return 0;

        Matrix A2   = A;

        {
            auto warn = error::enable_warnings(false);
            {
                Matrix A3 = A2;
                A3.get_array<Integer>();
            }
            {
                Matrix A3 = A2;
                A3.get_array<Real>();
            }
            {
                Matrix A3 = A2;
                A3.get_array<Float>();
            }
            {
                Matrix A3 = A2;
                A3.get_array<Complex>();
            }
            {
                Matrix A3 = A2;
                A3.get_array<Float_complex>();
            }
        };

        switch (mt)
        {
            case mat_code::integer_dense:
            {
                using value_type        = Integer;
                const value_type* ptr   = A.get_array<value_type>();

                if (!ptr)
                    return 0;

                Real dif    = norm_1(A(1)-ptr[0]);
                return dif; 
            }
            case mat_code::real_dense:
            {
                using value_type        = Real;
                const value_type* ptr   = A.get_array<value_type>();

                if (!ptr)
                    return 0;

                Real dif    = norm_1(A(1)-ptr[0]);
                return dif;
            }
            case mat_code::float_dense:
            {
                using value_type        = Float;
                const value_type* ptr   = A.get_array<value_type>();

                if (!ptr)
                    return 0;

                Real dif    = norm_1(A(1)-ptr[0]);
                return dif;
            }
            case mat_code::complex_dense:
            {
                using value_type        = Complex;
                const value_type* ptr   = A.get_array<value_type>();

                if (!ptr)
                    return 0;

                Real dif    = norm_1(A(1)-ptr[0]);
                return dif;
            }
            case mat_code::float_complex_dense:
            {
                using value_type        = Float_complex;
                const value_type* ptr   = A.get_array<value_type>();

                if (!ptr)
                    return 0;

                Real dif    = norm_1(A(1)-ptr[0]);
                return dif;
            }
        };

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
Real test_function_get_const_array::eval_scalar(const Scalar& ,bool,int code )
{
    (void)code;
    return 0.;
};
Real test_function_get_impl_unique::eval_mat(const Matrix& mat,bool,int code )
{
    // test whether get_impl_unique clears struct_flag correcty
    (void)code;
    try
    {
        Matrix A            = mat;
        struct_flag mat_sf  = mat.get_struct();
        mat_code	mt      = A.get_matrix_code();
        
        if (mt < mat_code::integer_dense) 
            return 0.; // scalars

        switch (mt)
        {
            case mat_code::integer_dense:
                A.get_impl_unique<raw::Matrix<Integer, struct_dense>>();
                break;
            case mat_code::real_dense:
                A.get_impl_unique<raw::Matrix<Real, struct_dense>>();
                break;
            case mat_code::float_dense:
                A.get_impl_unique<raw::Matrix<Float, struct_dense>>();
                break;
            case mat_code::complex_dense:
                A.get_impl_unique<raw::Matrix<Complex, struct_dense>>();
                break;
            case mat_code::float_complex_dense:
                A.get_impl_unique<raw::Matrix<Float_complex, struct_dense>>();
                break;
            case mat_code::integer_sparse:
                A.get_impl_unique<raw::Matrix<Integer, struct_sparse>>();
                break;
            case mat_code::real_sparse:
                A.get_impl_unique<raw::Matrix<Real, struct_sparse>>();
                break;
            case mat_code::float_sparse:
                A.get_impl_unique<raw::Matrix<Float, struct_sparse>>();
                break;
            case mat_code::complex_sparse:
                A.get_impl_unique<raw::Matrix<Complex, struct_sparse>>();
                break;
            case mat_code::float_complex_sparse:
                A.get_impl_unique<raw::Matrix<Float_complex, struct_sparse>>();
                break;
            case mat_code::integer_band:
                A.get_impl_unique<raw::Matrix<Integer, struct_banded>>();
                break;
            case mat_code::real_band:
                A.get_impl_unique<raw::Matrix<Real, struct_banded>>();
                break;
            case mat_code::float_band:
                A.get_impl_unique<raw::Matrix<Float, struct_banded>>();
                break;
            case mat_code::complex_band:
                A.get_impl_unique<raw::Matrix<Complex, struct_banded>>();
                break;
            case mat_code::float_complex_band:
                A.get_impl_unique<raw::Matrix<Float_complex, struct_banded>>();
                break;
            default:
                return 1.;
        }

        if (mat_sf != mat.get_struct()) 
            return 1.; // original struct flag changed
        
        return 0;
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
Real test_function_get_impl_unique::eval_scalar(const Scalar& ,bool,int code )
{
    (void)code;
    return 0.;
};

Real test_function_get_scal::eval_mat(const Matrix& A,bool,int code )
{
    (void)code;
    try
    {
        if (A.rows() != 1 || A.cols() != 1)
            return 0;

        matcl::value_code vt    = A.get_value_code();

        Matrix A2   = A;

        {
            auto warn = error::enable_warnings(false);
            {
                Matrix A3 = A2;
                A3.get_scalar<Integer>();
            }
            {
                Matrix A3 = A2;
                A3.get_scalar<Real>();
            }
            {
                Matrix A3 = A2;
                A3.get_scalar<Float>();
            }
            {
                Matrix A3 = A2;
                A3.get_scalar<Complex>();
            }
            {
                Matrix A3 = A2;
                A3.get_scalar<Float_complex>();
            }
        };

        switch (vt)
        {
            case value_code::v_integer:
            {
                Real dif = norm_1(A(1)-A.get_scalar<Integer>());
                return dif; 
            }
            case value_code::v_float:
            {
                Real dif = norm_1(A(1)-A.get_scalar<Float>());
                return dif; 
            }
            case value_code::v_float_complex:
            {
                Real dif = norm_1(A(1)-A.get_scalar<Float_complex>());
                return dif; 
            }
            case value_code::v_real:
            {
                Real dif = norm_1(A(1)-A.get_scalar<Real>());
                return dif; 
            }
            case value_code::v_complex:
            {
                Real dif = norm_1(A(1)-A.get_scalar<Complex>());
                return dif; 
            };
        };
        return 0;
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
Real test_function_get_scal::eval_scalar(const Scalar& ,bool,int code )
{
    (void)code;
    return 0.;
};

Real test_function_get_scal_unique::eval_mat(const Matrix& mat,bool,int code )
{
    (void)code;
    try
    {
        if (mat.rows() != 1 || mat.cols() != 1)
            return 0;

        Matrix A                = mat;
        matcl::value_code vt    = A.get_value_code();

        Matrix A2   = A;

        {
            auto warn = error::enable_warnings(false);
            {
                Matrix A3 = A2;
                A3.get_scalar_unique<Integer>();
            }
            {
                Matrix A3 = A2;
                A3.get_scalar_unique<Real>();
            }
            {
                Matrix A3 = A2;
                A3.get_scalar_unique<Float>();
            }
            {
                Matrix A3 = A2;
                A3.get_scalar_unique<Complex>();
            }
            {
                Matrix A3 = A2;
                A3.get_scalar_unique<Float_complex>();
            }
        };

        switch (vt)
        {
            case value_code::v_integer:
            {
                using value_type    = Integer;
                value_type& val     = A.get_scalar_unique<value_type>();

                Real dif    = norm_1(A(1)-val);
                val         = 1;
                dif         += norm_1(A(1)-1);

                return dif; 
            }
            case value_code::v_float:
            {
                using value_type    = Float;
                value_type& val     = A.get_scalar_unique<value_type>();

                Real dif    = norm_1(A(1)-val);
                val         = 1;
                dif         += norm_1(A(1)-1);

                return dif; 
            }
            case value_code::v_float_complex:
            {
                using value_type    = Float_complex;
                value_type& val     = A.get_scalar_unique<value_type>();

                Real dif    = norm_1(A(1)-val);
                val         = 1;
                dif         += norm_1(A(1)-1);

                return dif; 
            }
            case value_code::v_real:
            {
                using value_type    = Real;
                value_type& val     = A.get_scalar_unique<value_type>();

                Real dif    = norm_1(A(1)-val);
                val         = 1;
                dif         += norm_1(A(1)-1);

                return dif; 
            }
            case value_code::v_complex:
            {
                using value_type    = Complex;
                value_type& val     = A.get_scalar_unique<value_type>();

                Real dif    = norm_1(A(1)-val);
                val         = 1;
                dif         += norm_1(A(1)-1);

                return dif; 
            }
        };
        return 0;
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
Real test_function_get_scal_unique::eval_scalar(const Scalar& ,bool,int code )
{
    (void)code;
    return 0.;
};

Real test_function_mat_compile::eval_mat(const Matrix& mat,bool,int code )
{
    (void)code;
    //just test if given functions are defined

    try
    {
        //constructurs
        {
            Matrix m;
            m       = Matrix(true);
            m       = Matrix(m);
            m       = Matrix(std::move(m));
            m       = Matrix(1);

            m       = Matrix(unsigned char(1));
            m       = Matrix(signed char(1));
            m       = Matrix(unsigned short(1));
            m       = Matrix(signed short(1));
            m       = Matrix(unsigned int(1));
            m       = Matrix(signed int(1));
            m       = Matrix(unsigned long(1));
            m       = Matrix(signed long(1));
            m       = Matrix(unsigned long long(1));
            m       = Matrix(signed long long(1));
            m       = Matrix(float(1.));
            m       = Matrix(double(1.));
            m       = Matrix(long double(1.));
            m       = Matrix(Complex(1.));
            m       = Matrix(Float_complex(1.));
            m       = Matrix(Float(1.));
            m       = Matrix(Real(1.));
            m       = Matrix(Integer(1.));
        };

        Real dif    = 0;

        Integer r   = mat.rows();
        Integer c   = mat.cols();
        Integer l   = mat.length();
        Integer nz  = mat.structural_nnz();
        Integer ld1 = mat.structural_ldiags(false);
        Integer ld2 = mat.structural_ldiags(true);
        Integer ud1 = mat.structural_udiags(false);
        Integer ud2 = mat.structural_udiags(true);
        Real num    = mat.numel();
        
        bool af     = mat.all_finite();
        bool is_e   = mat.is_empty();
        bool is_s   = mat.is_scalar();
        bool is_sq  = mat.is_square();
        bool is_v   = mat.is_vector();
        bool is_mt  = mat.is_matrix_type();
        bool is_st  = mat.is_scalar_type();

        if (Real(r) * Real(c) != num)
            dif     += 1;

        if (is_e == true && l != 0)
            dif     += 1;
        if (is_e == false && l != std::max(r,c))
            dif     += 1;
        
        if (nz < 0 || Real(nz) > num)
            dif     += 1;

        if (ld1 < 0 || ld1 > r)
            dif     += 1;
        if (ld2 < 0 || ld2 > r)
            dif     += 1;

        if (ud1 < 0 || ud1 > c)
            dif     += 1;
        if (ud2 < 0 || ud2 > c)
            dif     += 1;
        if (ld1 < ld2)
            dif     += 1;
        if (ud1 < ud2)
            dif     += 1;

        if (is_e == true && (r != 0 && c != 0))
            dif     += 1;
        if (is_e == false && (r == 0 || c == 0))
            dif     += 1;

        if (is_s == true && (r != 1 || c != 1))
            dif     += 1;
        if (is_s == false && (r == 1 && c == 1))
            dif     += 1;

        if (is_sq == true && (r != c))
            dif     += 1;
        if (is_sq == false && (r == c))
            dif     += 1;

        if (is_v == true && (r != 1 && c != 1 && is_e == false))
            dif     += 1;
        if (is_v == false && (r == 1 || c == 1 || is_e == true))
            dif     += 1;

        if (is_st == true && is_s == false)
            dif     += 1;

        Matrix sum  = sum_vec(mat);

        if (af == true && (bool)is_finite(sum) == false)
            dif     += 1;
        if (af == false && (bool)is_finite(sum) == true)
            dif     += 1;

        value_code vc   = mat.get_value_code();
        struct_code sc  = mat.get_struct_code();
        mat_code mc     = mat.get_matrix_code();
        mat_code mc2    = matrix_traits::get_matrix_type(vc,sc);

        if (mc != mc2)
            dif         += 1;

        {
            Matrix m2   = mat;
            bool is_un  = m2.is_unique();

            if (is_un == true && is_mt == true)
                dif     += 1;

            m2.make_unique();

            is_un       = m2.is_unique();

            if (is_un == false && is_mt == true)
                dif     += 1;
        };

        ti::ti_object ti = mat.get_type();
        (void)ti;

        struct_flag sf  = mat.get_struct();
        mat.set_struct(sf);
        mat.add_struct(sf);

        struct_flag sf2 = mat.get_struct();

        if (sf2 != sf)
            dif         += 1;

        Matrix m2       = mat.clone();
        if (m2.is_unique() == false)
            dif         += 1;
        if (m2.get_matrix_code() != mc)
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
Real test_function_mat_compile::eval_scalar(const Scalar& ,bool,int code )
{
    (void)code;
    return 0.;
};

Real test_function_resize_reserve::check_struct_eq(struct_flag s1, struct_flag s2)
{
    if (s1 == s2)
        return 0.0;
    else
        return 1.0;
}

Real test_function_resize_reserve::eval_mat(const Matrix& mat1,bool,int code )
{
    (void)code;
    // crude tests against reserve-related heap coruptions
    {
        Matrix test_big_reserve1 = mat1;
        test_big_reserve1.reserve(test_big_reserve1.rows(), test_big_reserve1.cols() + 10000);
    }
    {
        Matrix test_big_reserve2 = mat1;
        test_big_reserve2.reserve(test_big_reserve2.rows() + 10000, test_big_reserve2.cols());
    }
    {
        Matrix test_big_resize1 = mat1;
        test_big_resize1.resize(test_big_resize1.rows(), test_big_resize1.cols() + 10000);
    }
    {
        Matrix test_big_resize2 = mat1;
        test_big_resize2.resize(test_big_resize2.rows() + 10000, test_big_resize2.cols());
    }

    try
    {
        Matrix mat  = mat1;
        Matrix B    = mat;

        Real dif    = 0;

        struct_flag old_struct  = mat.get_struct();

        mat.reserve(mat.rows() + 4, mat.cols() + 3);

        if (has_struct(mat, old_struct) == false)
            dif         += 1.0;

        Matrix mat_f    = full(mat1);
        Matrix B_f      = mat_f;
        mat_f.reserve(mat_f.rows() + 4, mat_f.cols() + 3);
        check_struct(mat_f);

        check_struct(mat);
        dif             += norm_1(B - mat);
        dif             += norm_1(mat - mat_f);
        
        B.resize(mat.rows() + 2, mat.cols() + 2);
        check_struct(B);

        mat.resize(mat.rows() + 2, mat.cols() + 2);
        check_struct(mat);

        B_f.resize(mat_f.rows() + 2, mat_f.cols() + 2);
        check_struct(B_f);

        mat_f.resize(mat_f.rows() + 2, mat_f.cols() + 2);
        check_struct(mat_f);

        dif     += norm_1(B - mat);
        dif     += norm_1(mat - mat_f);
        dif     += norm_1(B - B_f);

        Matrix C = mat;
        mat.resize(mat.rows() - 1, mat.cols() - 1);
        check_struct(mat);

        Matrix tmp = C(colon(1,C.rows()-1),colon(1,C.cols()-1));
        check_struct(tmp);

        Matrix C_f = mat_f;
        mat_f.resize(mat_f.rows() - 1, mat_f.cols() - 1);
        check_struct(mat_f);

        Matrix tmp_f = C_f(colon(1,C_f.rows()-1),colon(1,C_f.cols()-1));

        dif     += norm_1(mat - tmp);
        dif     += norm_1(mat - mat_f);
        dif     += norm_1(tmp - tmp_f);

        mat = C;
        mat.resize(mat.rows(), mat.cols() - 1);
        check_struct(mat);

        tmp = C(colon(),colon(1,C.cols()-1));
        check_struct(tmp);

        mat_f = C_f;
        mat_f.resize(mat_f.rows(), mat_f.cols() - 1);
        check_struct(mat_f);

        tmp_f = C_f(colon(),colon(1,C_f.cols()-1));

        dif     += norm_1(mat - tmp);
        dif     += norm_1(mat - mat_f);
        dif     += norm_1(tmp - tmp_f);

        return dif ;
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
Real test_function_resize_reserve::eval_scalar(const Scalar& ,bool,int code )
{
    (void)code;
    return 0.;
};

Real test_function_resize_reserve_b::eval_mat(const Matrix& mat1,bool,int code )
{
    (void)code;
    try
    {
        if (mat1.get_struct_code() != struct_code::struct_banded)
            return 0;

        struct_flag old_struct  = mat1.get_struct();

        Real diff   = 0.0;

        Matrix mat  = mat1;
        Matrix B    = mat;
        mat.reserve_band(mat.rows() + 4, mat.cols() + 3, -(mat.structural_ldiags(false) + 2), 
                         mat.structural_udiags(false) + 2);
        check_struct(mat);

        Matrix mat_f = full(mat1);
        Matrix B_f  = mat_f;
        mat_f.reserve_band(mat_f.rows() + 4, mat_f.cols() + 3, -(mat_f.structural_ldiags(false) + 2), 
                           mat_f.structural_udiags(false) + 2);
        check_struct(mat_f);

        B.resize_band(mat.rows() + 3, mat.cols() + 3, -(mat.structural_ldiags(false) + 1), 
                      mat.structural_udiags(false) + 1);
        check_struct(B);  

        mat.resize_band(mat.rows() + 3, mat.cols() + 3, -(mat.structural_ldiags(false) + 1), 
                        mat.structural_udiags(false) + 1);
        check_struct(mat);

        B_f.resize_band(mat_f.rows() + 3, mat_f.cols() + 3, -(mat_f.structural_ldiags(false) + 1), 
                        mat_f.structural_udiags(false) + 1);
        check_struct(B_f);

        mat_f.resize_band(mat_f.rows() + 3, mat_f.cols() + 3, -(mat_f.structural_ldiags(false) + 1), 
                          mat_f.structural_udiags(false) + 1);
        check_struct(mat_f);

        diff        += norm_1(B - mat);
        diff        += norm_1(B - B_f);
        diff        += norm_1(mat - mat_f);

        Matrix C    = mat;
        mat.resize_band(mat.rows() - 1, mat.cols() - 1, -(mat.structural_ldiags(false) - 1), 
                        mat.structural_udiags(false) - 1);
        check_struct(mat);

        Matrix tmp  = C(colon(1,C.rows()-1),colon(1,C.cols()-1));
        {
            Integer l   = std::max(C.structural_ldiags(false) - 1,0);
            Integer u   = std::max(C.structural_udiags(false) - 1,0);
            tmp         = tril(triu(tmp,-l),u);
        };
        check_struct(tmp);

        Matrix C_f = mat_f;
        mat_f.resize_band(mat_f.rows() - 1, mat_f.cols() - 1, -(mat_f.structural_ldiags(false) - 1), 
                          mat_f.structural_udiags(false) - 1);
        check_struct(mat_f);

        Matrix tmp_f = C_f(colon(1,C_f.rows()-1),colon(1,C_f.cols()-1));
        {
            Integer l   = mat_f.structural_ldiags(false);
            Integer u   = mat_f.structural_udiags(false);
            tmp_f       = tril(triu(tmp_f,-l),u);
        }
        
        diff    += norm_1(mat - tmp);
        diff    += norm_1(mat - mat_f);
        diff    += norm_1(tmp - tmp_f);

        mat = C;
        mat.resize_band(C.rows(), C.cols() - 1, -(C.structural_ldiags(false) - 1), C.structural_udiags(false) - 1);
        check_struct(mat);

        tmp = C(colon(),colon(1,C.cols()-1));
        {
            Integer l   = std::max(C.structural_ldiags(false) - 1,0);
            Integer u   = std::max(C.structural_udiags(false) - 1,0);
            tmp         = tril(triu(tmp,-l),u);
        };
        check_struct(tmp);

        mat_f = C_f;
        mat_f.resize_band(C_f.rows(), C_f.cols() - 1, -(C_f.structural_ldiags(false) - 1), 
                          C_f.structural_udiags(false) - 1);
        check_struct(mat_f);

        tmp_f = C_f(colon(),colon(1,C_f.cols()-1));
        {
            Integer l   = mat_f.structural_ldiags(false);
            Integer u   = mat_f.structural_udiags(false);
            tmp_f       = tril(triu(tmp_f,-l),u);
        };

        diff    += norm_1(mat - tmp);
        diff    += norm_1(mat - mat_f);
        diff    += norm_1(tmp - tmp_f);
        return diff;
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
Real test_function_resize_reserve_b::eval_scalar(const Scalar& ,bool,int code )
{
    (void)code;
    return 0.;
};

};};
