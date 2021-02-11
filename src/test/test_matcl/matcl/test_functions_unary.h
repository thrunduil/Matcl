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

#include "matcl-matrep/matcl_matrep.h"
#include "test/test_matcl/framework/matrix_set/matrix_set.h"
#include "test/test_matcl/framework/matrix_utils.h"
#include "matcl-linalg/matcl_linalg.h"

namespace matcl { namespace test
{

class unary_functions_list
{
    private:
        const matrix_set&	m_tests;        
        options				m_options;

    public:
        rand_matrix_ptr     m_rand;

    public:
        void		make(options opts);

        Matrix		get_matrix(int code) const;

    private:
        void    test_real();            void    test_conj();		
        void    test_imag();			void    test_abs();
        void    test_arg();				void    test_angle();    
        void    test_abs2();

        void    test_sqrt();            void    test_acsch();
        void    test_pow2();			void    test_exp();
        void    test_log();				void    test_log2();
        void    test_log10();			void    test_floor();
        void    test_ceil();			void    test_round();
        void    test_fix();				void    test_trunc();
        void    test_sign();			void    test_sin();
        void    test_cos();				void    test_tan();
        void    test_cot();				void    test_sec();
        void    test_csc();				void    test_asin();
        void    test_acos();			void    test_atan();
        void    test_acot();			void    test_asec();
        void    test_acsc();			void    test_sinh();
        void    test_cosh();			void    test_tanh();
        void    test_coth();			void    test_sech();
        void    test_csch();			void    test_asinh();
        void    test_acosh();			void    test_atanh();
        void    test_acoth();			void    test_asech();        

        void    test_sqrt_c();          void    test_log_c();
        void    test_log2_c();          void    test_log10_c();
        void    test_asin_c();          void    test_acos_c();
        void    test_asec_c();          void    test_acsc_c();
        void    test_acosh_c();         void    test_atanh_c();
        void    test_acoth_c();         void    test_asech_c();

        void    test_sqrt1pm1();        void    test_cbrt();
        void    test_log1p();           void    test_sqrt1pm1_c();
        void    test_log1p_c();         void    test_op_plus();
        void    test_expm1();           void    test_exp10();
        void    test_signbit();         void    test_logb();
        void    test_ilogb();           void    test_exp2();
        void    test_fpclassify();      void    test_nextabove();
        void    test_nextbelow();       void    test_eps();
        void    test_expi();            void    test_frexp();
        void    test_modf();

        void    test_bernoulli_b2n();   void    test_max_bernoulli_b2n();
        void    test_prime();           void    test_prime_max_count();
        void    test_factorial();       void    test_double_factorial();
        void    test_rising_factorial();void    test_binomial_coefficient();
        void    test_falling_factorial();
        
        void    test_eval_scalar_func();

        void    test_is_true();			void    test_is_false();		
        void    test_ifloor();			void    test_iceil();
        void    test_iround();			void    test_ifix();
        void    test_itrunc();			void    test_is_inf();
        void    test_is_nan();			void    test_is_finite();
        void    test_is_regular();      void    test_is_int();
        void    test_is_real();         void    test_is_normal();
        void    test_is_scalar_true();	void    test_is_scalar_false();
        void    test_isign();   

        void    test_op_un_minus();		void    test_op_neg();
        void    test_ldexp();           void    test_scalbn();
        void    test_inv();             void    test_invs();

        template<class Func>
        void    test_function();

    public:
        unary_functions_list(const matrix_set&,rand_matrix_ptr rand);

    private:
        unary_functions_list(const unary_functions_list&) = delete;
        unary_functions_list& operator=(const unary_functions_list&) = delete;
};

class test_function_abs2 : public unary_function
{
    public:
        test_function_abs2(){};

        virtual Real eval_mat(const Matrix& mat,bool show_res,int code) override;
        virtual Real eval_scalar(const Scalar& s, bool show_res,int code) override;

        template<class Val>
        Real eval_scal_func(const Val& value) const
        {
            using Val_c = typename md::complex_type<Val>::type;
            Val_c val_c = Val_c(value);

            auto res1   = abs(value);
            auto res2   = abs2(value);
            auto res3   = abs2(val_c);

            Real dif    = norm_1(res2 - res1 * res1);
            dif         += norm_1(res2 - res3);

            value_code vc   = matrix_traits::value_code<Val>::value;            

            if (dif < res2 * constants::eps(vc) * 10.0)
                dif     = 0.0;

            return dif;
        }
};

class test_function_is_true : public unary_function
{
    public:
        test_function_is_true()	{};

        virtual Real eval_mat(const Matrix& mat,bool show_res,int code) override;
        virtual Real eval_scalar(const Scalar& s, bool show_res,int code) override;

        template<class Val>
        Real eval_scal_func(const Val& value) const
        {
            using Val_c = typename md::complex_type<Val>::type;
            Val_c val_c = Val_c(value);

            Real dif		= norm_1(Matrix(is_true(value)) - is_true(full(value)));
            dif				+=norm_1(Matrix(is_true(value)) - Matrix(is_true(val_c)));

            return dif;
        }
};

class test_function_is_false : public unary_function
{
    public:
        test_function_is_false()	{};

        virtual Real eval_mat(const Matrix& mat,bool show_res,int code) override;
        virtual Real eval_scalar(const Scalar& s, bool show_res,int code) override;

        template<class Val>
        Real eval_scal_func(const Val& value) const
        {
            using Val_c = typename md::complex_type<Val>::type;
            Val_c val_c = Val_c(value);

            Real dif		= norm_1(Matrix(is_false(value)) - is_false(full(value)));
            dif				+=norm_1(Matrix(is_false(value)) - Matrix(is_false(val_c)));

            return dif;
        }
};

class test_function_op_un_minus : public unary_function
{
    private:
    public:
        test_function_op_un_minus()	{};

        virtual Real eval_mat(const Matrix& mat,bool show_res,int code) override;
        virtual Real eval_scalar(const Scalar& s, bool show_res,int code) override;

        template<class Val>
        Real eval_scal_func(const Val& value) const
        {
            using Val_c = typename md::complex_type<Val>::type;
            Val_c val_c = Val_c(value);

            Real dif		= norm_1(-value - (-full(value)));
            dif				+=norm_1(-value - (-val_c));
            dif				+=norm_1(-value - uminus(value));

            return dif;
        }
};

class test_function_inv : public unary_function
{
    private:
    public:
        test_function_inv()	{};

        virtual Real eval_mat(const Matrix& mat,bool show_res,int code) override;
        virtual Real eval_scalar(const Scalar& s, bool show_res,int code) override;

        template<class Val>
        Real eval_scal_func(const Val& value) const
        {
            bool is_zero    = matcl::is_zero(value);

            if (is_zero == true)
                return 0.0;

            using Val_c = typename md::complex_type<Val>::type;
            Val_c val_c = Val_c(value);

            Real dif		= norm_1(inv(value) - inv(full(value)));
            dif				+=norm_1(inv(value) - inv(val_c));
            dif				+=norm_1(inv(value) - div(Val(1), value));

            return dif;
        }
};

class test_function_invs : public unary_function
{
    private:
    public:
        test_function_invs()	{};

        virtual Real eval_mat(const Matrix& mat,bool show_res,int code) override;
        virtual Real eval_scalar(const Scalar& s, bool show_res,int code) override;

        template<class Val>
        Real eval_scal_func(const Val& value) const
        {
            using Val_c = typename md::complex_type<Val>::type;
            Val_c val_c = Val_c(value);

            Real dif		= norm_1(invs(value) - invs(full(value)));
            dif				+=norm_1(invs(value) - invs(val_c));
            dif				+=norm_1(invs(value) - div(Val(1), value));

            return dif;
        }
};

class test_function_op_neg : public unary_function
{
    public:
        test_function_op_neg()	{};

        virtual Real eval_mat(const Matrix& mat,bool show_res,int code) override;
        virtual Real eval_scalar(const Scalar& s, bool show_res,int code) override;

        template<class Val>
        Real eval_scal_func(const Val& value) const
        {
            using Val_c = typename md::complex_type<Val>::type;
            Val_c val_c = Val_c(value);

            Real dif		= norm_1(Matrix(!value) - (~full(value)));
            dif				+=norm_1(Matrix(!value) - Matrix(!val_c));
            dif				+=norm_1(Matrix(!value) - Matrix(neg(val_c)));

            return dif;
        }
};

class test_function_fpclassify : public unary_function
{
    public:
        test_function_fpclassify()	{};

        virtual Real eval_mat(const Matrix& mat,bool show_res,int code) override;
        virtual Real eval_scalar(const Scalar& s, bool show_res,int code) override;

        template<class Val>
        Real eval_scal_func(const Val& value) const
        {
            fpclassify(real(value));
            return 0.0;
        }
};

class test_function_nextabove : public unary_function
{
    public:
        test_function_nextabove()	{};

        virtual Real eval_mat(const Matrix& mat,bool show_res,int code) override;
        virtual Real eval_scalar(const Scalar& s, bool show_res,int code) override;

        template<class Val>
        Real eval_scal_func(const Val& value) const
        {
            using Val_c = typename md::complex_type<Val>::type;
            using Val_r = typename md::real_type<Val_c>::type;
            Val_c val_c = Val_c(value);

            Val_r v1    = nextabove(real(value));
            Val_r v3    = nextafter(real(value), constants::inf());

            Real dif    = norm_1(Matrix(v1) - Matrix(v3));

            return dif;
        }
};

class test_function_eps : public unary_function
{
    public:
        test_function_eps()	{};

        virtual Real eval_mat(const Matrix& mat,bool show_res,int code) override;
        virtual Real eval_scalar(const Scalar& s, bool show_res,int code) override;

        template<class Val>
        Real eval_scal_func(const Val& value) const
        {
            using Val_c = typename md::complex_type<Val>::type;
            using Val_r = typename md::real_type<Val_c>::type;
            Val_c val_c = Val_c(value);

            Val_r v1    = eps(value);
            Val_r v2    = eps(val_c);

            Real dif    = norm_1(Matrix(v1) - eps(Matrix(value)));
            dif         +=norm_1(Matrix(v1) - Matrix(v2));

            return dif;
        }
};

class test_function_nextbelow : public unary_function
{
    public:
        test_function_nextbelow()	{};

        virtual Real eval_mat(const Matrix& mat,bool show_res,int code) override;
        virtual Real eval_scalar(const Scalar& s, bool show_res,int code) override;

        template<class Val>
        Real eval_scal_func(const Val& value) const
        {
            using Val_c = typename md::complex_type<Val>::type;
            using Val_r = typename md::real_type<Val_c>::type;
            Val_c val_c = Val_c(value);

            Val_r v1    = nextbelow(real(value));
            Val_r v3    = nextafter(real(value), -constants::inf());

            Real dif    = norm_1(Matrix(v1) - Matrix(v3));

            return dif;
        }
};

class test_function_ldexp : public unary_function
{
    public:
        test_function_ldexp()	{};

        virtual Real eval_mat(const Matrix& mat,bool show_res,int code) override;
        virtual Real eval_scalar(const Scalar& s, bool show_res,int code) override;

        template<class Val>
        Real eval_scal_func(const Val& value) const
        {
            using Val_c = typename md::complex_type<Val>::type;
            Val_c val_c = Val_c(value);

            Real dif    = 0;

            for (Integer exp = -2; exp <= 2; ++exp)
            {
                auto res    = ldexp(value, exp);
                Real difl   = 0.0;
                difl        += norm_1(res - value * pow2(exp)); 
                difl        += norm_1(res - ldexp(val_c, exp));

                using T         = decltype(res);
                value_code vc   = matrix_traits::value_code<T>::value;
                Real tol        = norm_1(res) * constants::eps(vc) * 10.0;

                if (difl < tol)
                    difl    = 0;

                dif         += difl;
            };
            
            return dif;
        }
};

class test_function_scalbn : public unary_function
{
    public:
        test_function_scalbn()	{};

        virtual Real eval_mat(const Matrix& mat,bool show_res,int code) override;
        virtual Real eval_scalar(const Scalar& s, bool show_res,int code) override;

        template<class Val>
        Real eval_scal_func(const Val& value) const
        {
            using Val_c = typename md::complex_type<Val>::type;
            Val_c val_c = Val_c(value);

            Real dif    = 0;

            for (Integer exp = -2; exp <= 2; ++exp)
            {
                auto res    = scalbn(value, exp);
                Real difl   = 0.0;
                difl        += norm_1(res - value * pow(2, exp)); 
                difl        += norm_1(res - scalbn(val_c, exp));

                using T         = decltype(res);
                value_code vc   = matrix_traits::value_code<T>::value;
                Real tol        = norm_1(res) * constants::eps(vc) * 10.0;

                if (difl < tol)
                    difl    = 0;

                dif         += difl;
            };
            
            return dif;
        }
};

class test_function_modf : public unary_function
{
    public:
        test_function_modf()	{};

        virtual Real eval_mat(const Matrix& mat,bool show_res,int code) override;
        virtual Real eval_scalar(const Scalar& s, bool show_res,int code) override;

        template<class Val>
        Real eval_scal_func(const Val& value) const
        {
            using Val_r = typename md::unify_types<Float, typename md::real_type<Val>::type>::type;

            Val_r frac;
            Val_r res   =  matcl::modf(real(value), frac);
            Val_r val2  = res + frac;

            Real dif    = norm_1(val2 - real(value));
            Real tol    = matcl::eps(val2) * 10.0;

            if (dif < tol)
                dif    = 0;

            return dif;
        }
};

class test_function_frexp : public unary_function
{
    public:
        test_function_frexp()	{};

        virtual Real eval_mat(const Matrix& mat,bool show_res,int code) override;
        virtual Real eval_scalar(const Scalar& s, bool show_res,int code) override;

        template<class Val>
        Real eval_scal_func(const Val& value) const
        {
            using Val_r = typename md::unify_types<Float, typename md::real_type<Val>::type>::type;

            Integer exp;
            Val_r res   = matcl::frexp(real(value), exp);
            Val_r val2  = ldexp(res, exp);

            Real dif    = norm_1(val2 - real(value));
            Real tol    = matcl::eps(val2) * 10.0;

            if (dif < tol)
                dif    = 0;

            return dif;
        }
};

class test_function_bernoulli_b2n : public unary_function
{
    public:
        test_function_bernoulli_b2n()	{};

        virtual Real eval_mat(const Matrix& mat,bool show_res,int code) override;
        virtual Real eval_scalar(const Scalar& s, bool show_res,int code) override;
};

class test_function_max_bernoulli_b2n : public unary_function
{
    public:
        test_function_max_bernoulli_b2n()	{};

        virtual Real eval_mat(const Matrix& mat,bool show_res,int code) override;
        virtual Real eval_scalar(const Scalar& s, bool show_res,int code) override;
};

class test_function_prime : public unary_function
{
    public:
        test_function_prime()	{};

        virtual Real eval_mat(const Matrix& mat,bool show_res,int code) override;
        virtual Real eval_scalar(const Scalar& s, bool show_res,int code) override;
};
class test_function_prime_max_count : public unary_function
{
    public:
        test_function_prime_max_count()	{};

        virtual Real eval_mat(const Matrix& mat,bool show_res,int code) override;
        virtual Real eval_scalar(const Scalar& s, bool show_res,int code) override;
};
class test_function_factorial : public unary_function
{
    public:
        test_function_factorial()	{};

        virtual Real eval_mat(const Matrix& mat,bool show_res,int code) override;
        virtual Real eval_scalar(const Scalar& s, bool show_res,int code) override;
};
class test_function_double_factorial : public unary_function
{
    public:
        test_function_double_factorial()	{};

        virtual Real eval_mat(const Matrix& mat,bool show_res,int code) override;
        virtual Real eval_scalar(const Scalar& s, bool show_res,int code) override;
};
class test_function_rising_factorial : public unary_function
{
    public:
        test_function_rising_factorial()	{};

        virtual Real eval_mat(const Matrix& mat,bool show_res,int code) override;
        virtual Real eval_scalar(const Scalar& s, bool show_res,int code) override;

        template<class Val>
        Real eval_scal_func(const Val& value_0) const
        {
            using Val_rc    = typename md::unify_types<Val, Float>::type;
            using Val_r     = typename md::real_type<Val_rc>::type;

            Val_r value     = Val_r(real(value_0));

            Val_r res_11    = rising_factorial(real(value),1);
            Val_r res_12    = value;
            Val_r res_21    = rising_factorial(real(value),2);
            Val_r res_22    = value * (value + Val_r(1));
            Val_r res_31    = rising_factorial(real(value),3);
            Val_r res_32    = value * (value + Val_r(1)) * (value + Val_r(2));

            Real dif    = norm_1(res_11 - res_12);
            dif         += norm_1(res_21 - res_22);
            
            dif         += norm_1(res_31 - res_32);

            if (abs(dif) < 10.0 * (eps(res_31) + eps(value)))
                dif     = 0.0;

            return dif;
        }
};

class test_function_falling_factorial : public unary_function
{
    public:
        test_function_falling_factorial()	{};

        virtual Real eval_mat(const Matrix& mat,bool show_res,int code) override;
        virtual Real eval_scalar(const Scalar& s, bool show_res,int code) override;

        template<class Val>
        Real eval_scal_func(const Val& value_0) const
        {
            using Val_rc    = typename md::unify_types<Val, Float>::type;
            using Val_r     = typename md::real_type<Val_rc>::type;

            Val_r value     = Val_r(real(value_0));

            Val_r res_11    = falling_factorial(real(value),1);
            Val_r res_12    = value;
            Val_r res_21    = falling_factorial(real(value),2);
            Val_r res_22    = value * (value - Val_r(1));
            Val_r res_31    = falling_factorial(real(value),3);
            Val_r res_32    = value * (value - Val_r(1)) * (value - Val_r(2));

            Real dif    = norm_1(res_11 - res_12);
            dif         += norm_1(res_21 - res_22);
            
            dif         += norm_1(res_31 - res_32);

            if (abs(dif) < 10.0 * (eps(res_31) + eps(value)))
                dif     = 0.0;

            return dif;
        }
};

class test_function_binomial_coefficient : public unary_function
{
    public:
        test_function_binomial_coefficient()	{};

        virtual Real eval_mat(const Matrix& mat,bool show_res,int code) override;
        virtual Real eval_scalar(const Scalar& s, bool show_res,int code) override;
};

class func_sin
{
    public:
        template<class T>
        auto eval(const T& v) const -> decltype(sin(v))
        {
            return sin(v);
        };
};
class func_cos
{
    public:
        template<class T>
        auto eval(const T& v) const -> decltype(cos(v))
        {
            return cos(v);
        };
};

class test_function_eval_scalar_func : public unary_function
{
    public:
        test_function_eval_scalar_func()	{};

        virtual Real eval_mat(const Matrix& mat,bool show_res,int code) override;
        virtual Real eval_scalar(const Scalar& s, bool show_res,int code) override;

        template<class Val>
        Real eval_scal_func(const Val& value) const
        {
            using Val_c = typename md::complex_type<Val>::type;
            Val_c val_c = Val_c(value);

            func_cos fc;

            Real dif		= norm_1(eval_scalar_func(value, fc) - eval_scalar_func(value, fc, fc));
            dif		        += norm_1(eval_scalar_func(value, fc) - eval_scalar_func(Matrix(value), fc));
            dif				+=norm_1(eval_scalar_func(value, fc) - eval_scalar_func(val_c, fc));

            return dif;
        }
};

class test_function_signbit : public unary_function
{
    public:
        test_function_signbit()	{};

        virtual Real eval_mat(const Matrix& mat,bool show_res,int code) override;
        virtual Real eval_scalar(const Scalar& s, bool show_res,int code) override;

        template<class Val>
        Real eval_scal_func(const Val& value) const
        {
            using Val_c = typename md::complex_type<Val>::type;
            using Val_r = typename md::real_type<Val_c>::type;

            Val_r v1    = signbit(real(value));

            Real dif    =  norm_1(Matrix(v1) - signbit(full(real(value))));
            return dif;
        }
};

class test_function_ilogb : public unary_function
{
    public:
        test_function_ilogb()	{};

        virtual Real eval_mat(const Matrix& mat,bool show_res,int code) override;
        virtual Real eval_scalar(const Scalar& s, bool show_res,int code) override;

        template<class Val>
        Real eval_scal_func(const Val& value) const
        {
            try
            {
                Real dif		= norm_1(Matrix(ilogb(value)) - ilogb(full(value)));
                return dif;
            }
            catch(const error::function_not_defined_for_complex&)
            {
                if (any_vec(imag(value)) == true)
                    return 0.0;
                else
                    return 1.0;
            }
        }
};

template<class func_helper>
class scal_func : public unary_function
{
    public:
        virtual Real eval_mat(const Matrix& mat,bool,int code )
        {
            (void)code;
            try
            {		
                Matrix out_full = func_helper::eval(full(mat));
                value_code vc   = matrix_traits::complex_value_type(mat.get_value_code());
                matcl::mat_code nt = matrix_traits::get_matrix_type(vc,mat.get_struct_code());

                Matrix mat_c    = convert(mat,nt);
                Matrix out		= func_helper::eval(mat);
                Matrix out_c	= func_helper::eval(mat_c);
                check_struct(out);

                Real dif        = norm_1(out - out_full)/(norm_1(out) + constants::eps());
                dif             += norm_1(out - out_c)/(norm_1(out) + 1.);

                if (abs(dif) < 1000 * out.length() * constants::eps(mat.get_value_code()))
                    dif         = 0.;

                return dif;
            }
            catch(error::scalar_required)
            {
                return mat.cols() == 1 && mat.rows() == 1 ? 1. : 0.;
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

        virtual Real eval_scalar(const Scalar& mat,bool,int code )
        {
            (void)code;

            Real dif = 0;

            switch (mat.get_value_code())
            {
                case value_code::v_integer:
                {
                    Integer val = mat.get_int();
                    Matrix res1 = func_helper::eval(val);
                    dif         = norm_1(func_helper::eval_scal(val) - res1);

                    if(abs(dif) < 1000.*constants::eps()) 
                        dif = 0.;

                    dif = dif / (norm_1(res1) + constants::eps()*(1. + norm_1(val)));

                    if(abs(dif) < 10.*constants::eps()) 
                        dif = 0.;

                    break;
                }
                case value_code::v_float:
                {
                    Float val       = mat.get_float();
                    Matrix res1     = func_helper::eval(val);
                    Float_complex val_c	= val;
                    dif		        = norm_1(func_helper::eval_scal(val) - res1);
                    dif				+=norm_1(func_helper::eval_scal(val) - func_helper::eval_scal(val_c));

                    if(abs(dif) < 1000.*constants::f_eps()) 
                        dif = 0.; // Absolute tolerance

                    dif = dif / (norm_1(res1) + constants::f_eps()*(1. + norm_1(val)));

                    if(abs(dif) < 10.*constants::f_eps()) 
                        dif = 0.;

                    break;
                }
                case value_code::v_real:
                {
                    Real val		= mat.get_real();
                    Matrix res1     = func_helper::eval(val);
                    Complex val_c	= val;
                    dif		        = norm_1(func_helper::eval_scal(val) - res1);
                    dif				+=norm_1(func_helper::eval_scal(val) - func_helper::eval_scal(val_c));

                    if(abs(dif) < 1000.*constants::eps()) 
                        dif = 0.;

                    dif = dif / (norm_1(res1) + constants::eps()*(1. + norm_1(val)));

                    if(abs(dif) < 10.*constants::eps()) 
                        dif = 0.;
                    break;
                }
                case value_code::v_float_complex:
                {
                    Float_complex val   = mat.get_fcomplex();
                    Matrix res1         = func_helper::eval(val);
                    dif                 = norm_1(func_helper::eval_scal(val) - res1);

                    if(abs(dif) < 1000.*constants::f_eps()) 
                        dif = 0.;

                    dif = dif / (norm_1(res1) + constants::f_eps()*(1. + norm_1(val)));

                    if(abs(dif) < 10.*constants::f_eps()) 
                        dif = 0.;

                    break;
                }
                case value_code::v_complex:
                {
                    Complex val     = mat.get_complex();
                    Matrix res1     = func_helper::eval(val);
                    dif             = norm_1(func_helper::eval_scal(val) - res1);

                    if(abs(dif) < 1000.*constants::eps())
                        dif = 0.;

                    dif = dif / (norm_1(res1) + constants::eps()*(1. + norm_1(val)));

                    if(abs(dif) < 10.*constants::eps()) 
                        dif = 0.;
                    break;
                }
                default:
                {
                    return 0;
                }
            };

            return dif;
        };
};

struct h_ifloor
{ 
    static Matrix eval(const Matrix& m)			{ return ifloor(real(m)); }; 
    template<class T>
    static Matrix eval_scal(T m)				{ return ifloor(real(m)); }; 
};
struct h_iceil
{ 
    static Matrix eval(const Matrix& m)			{ return iceil(real(m)); }; 
    template<class T>
    static Matrix eval_scal(T m)				{ return iceil(real(m)); }; 
};
struct h_iround
{ 
    static Matrix eval(const Matrix& m)			{ return iround(real(m)); }; 
    template<class T>
    static Matrix eval_scal(T m)				{ return iround(real(m)); }; 
};
struct h_ifix
{ 
    static Matrix eval(const Matrix& m)			{ return ifix(real(m)); }; 
    template<class T>
    static Matrix eval_scal(T m)				{ return ifix(real(m)); }; 
};
struct h_itrunc
{ 
    static Matrix eval(const Matrix& m)			{ return itrunc(real(m)); }; 
    template<class T>
    static Matrix eval_scal(T m)				{ return itrunc(real(m)); }; 
};
struct h_sign
{ 
    static Matrix eval(const Matrix& m)			{ return sign(m); }; 
    template<class T>
    static Matrix eval_scal(T m)				{ return sign(m); }; 
};
struct h_isign
{ 
    static Matrix eval(const Matrix& m)			{ return isign(real(m)); }; 
    template<class T>
    static Matrix eval_scal(T m)				{ return isign(real(m)); }; 
};
struct h_isinf
{ 
    static Matrix eval(const Matrix& m)			{ return is_inf(m); }; 
    template<class T>
    static Matrix eval_scal(T m)				{ return Matrix(is_inf(m)); }; 
};
struct h_isnan
{ 
    static Matrix eval(const Matrix& m)			{ return is_nan(m); }; 
    template<class T>
    static Matrix eval_scal(T m)				{ return Matrix(is_nan(m)); }; 
};
struct h_isfin
{ 
    static Matrix eval(const Matrix& m)			{ return is_finite(m); }; 
    template<class T>
    static Matrix eval_scal(T m)				{ return Matrix(is_finite(m)); }; 
};
struct h_isreg
{ 
    static Matrix eval(const Matrix& m)			{ return is_regular(m); }; 
    template<class T>
    static Matrix eval_scal(T m)				{ return Matrix(is_regular(m)); }; 
};
struct h_isnorm
{ 
    static Matrix eval(const Matrix& m)			{ return is_normal(m); }; 
    template<class T>
    static Matrix eval_scal(T m)				{ return Matrix(is_normal(m)); }; 
};

struct h_isint
{ 
    static Matrix eval(const Matrix& m)			{ return is_int(m); }; 
    template<class T>
    static Matrix eval_scal(T m)				{ return Matrix(is_int(m)); }; 
};
struct h_isreal
{ 
    static Matrix eval(const Matrix& m)			{ return is_real(m); }; 
    template<class T>
    static Matrix eval_scal(T m)				{ return Matrix(is_real(m)); }; 
};

#pragma warning(push)
#pragma warning(disable: 4800)//forcing value to bool 'true' or 'false' (performance warning)
struct h_isstrue
{ 
    static Matrix eval(const Matrix& m)			{ return Matrix((bool)m); }; 
    template<class T>
    static Matrix eval_scal(T m)				{ return Matrix(bool(m)); }; 
};
struct h_issfalse
{ 
    static Matrix eval(const Matrix& m)			{ return Matrix(!(m)); }; 
    template<class T>
    static Matrix eval_scal(T m)				{ return Matrix(!(m)); }; 
};
#pragma warning(pop)
struct h_sin
{ 
    static Matrix eval(const Matrix& m)			{ return sin(m); }; 
    template<class T>
    static Matrix eval_scal(T m)				{ return sin(m); }; 
};
struct h_cos
{ 
    static Matrix eval(const Matrix& m)			{ return cos(m); }; 
    template<class T>
    static Matrix eval_scal(T m)				{ return cos(m); }; 
};
struct h_tan
{ 
    static Matrix eval(const Matrix& m)			{ return tan(m); }; 
    template<class T>
    static Matrix eval_scal(T m)				{ return tan(m); }; 
};
struct h_cot
{ 
    static Matrix eval(const Matrix& m)			{ return cot(m); }; 
    template<class T>
    static Matrix eval_scal(T m)				{ return cot(m); }; 
};
struct h_sec
{ 
    static Matrix eval(const Matrix& m)			{ return sec(m); }; 
    template<class T>
    static Matrix eval_scal(T m)				{ return sec(m); }; 
};
struct h_csc
{ 
    static Matrix eval(const Matrix& m)			{ return csc(m); }; 
    template<class T>
    static Matrix eval_scal(T m)				{ return csc(m); }; 
};
struct h_asin
{ 
    static Matrix eval(const Matrix& m)			{ return asin(m); }; 
    template<class T>
    static Matrix eval_scal(T m)				{ return asin(m); }; 
};
struct h_asin_c
{ 
    static Matrix eval(const Matrix& m)			{ return asin_c(m); }; 
    template<class T>
    static Matrix eval_scal(T m)				{ return asin_c(m); }; 
};
struct h_acos
{ 
    static Matrix eval(const Matrix& m)			{ return acos(m); }; 
    template<class T>
    static Matrix eval_scal(T m)				{ return acos(m); }; 
};
struct h_acos_c
{ 
    static Matrix eval(const Matrix& m)			{ return acos_c(m); }; 
    template<class T>
    static Matrix eval_scal(T m)				{ return acos_c(m); }; 
};
struct h_atan
{ 
    static Matrix eval(const Matrix& m)			{ return atan(m); }; 
    template<class T>
    static Matrix eval_scal(T m)				{ return atan(m); }; 
};
struct h_acot
{ 
    static Matrix eval(const Matrix& m)			{ return acot(m); }; 
    template<class T>
    static Matrix eval_scal(T m)				{ return acot(m); }; 
};
struct h_asec
{ 
    static Matrix eval(const Matrix& m)			{ return asec(m); }; 
    template<class T>
    static Matrix eval_scal(T m)				{ return asec(m); }; 
};
struct h_asec_c
{ 
    static Matrix eval(const Matrix& m)			{ return asec_c(m); }; 
    template<class T>
    static Matrix eval_scal(T m)				{ return asec_c(m); }; 
};
struct h_acsc
{ 
    static Matrix eval(const Matrix& m)			{ return acsc(m); }; 
    template<class T>
    static Matrix eval_scal(T m)				{ return acsc(m); }; 
};
struct h_acsc_c
{ 
    static Matrix eval(const Matrix& m)			{ return acsc_c(m); }; 
    template<class T>
    static Matrix eval_scal(T m)				{ return acsc_c(m); }; 
};
struct h_sinh
{ 
    static Matrix eval(const Matrix& m)			{ return sinh(m); }; 
    template<class T>
    static Matrix eval_scal(T m)				{ return sinh(m); }; 
};
struct h_cosh
{ 
    static Matrix eval(const Matrix& m)			{ return cosh(m); }; 
    template<class T>
    static Matrix eval_scal(T m)				{ return cosh(m); }; 
};
struct h_tanh
{ 
    static Matrix eval(const Matrix& m)			{ return tanh(m); }; 
    template<class T>
    static Matrix eval_scal(T m)				{ return tanh(m); }; 
};
struct h_coth
{ 
    static Matrix eval(const Matrix& m)			{ return coth(m); }; 
    template<class T>
    static Matrix eval_scal(T m)				{ return coth(m); }; 
};
struct h_sech
{	
    static Matrix eval(const Matrix& m)			{ return sech(m); }; 
    template<class T>
    static Matrix eval_scal(T m)				{ return sech(m); }; 
};
struct h_csch
{ 
    static Matrix eval(const Matrix& m)			{ return csch(m); }; 
    template<class T>
    static Matrix eval_scal(T m)				{ return csch(m); }; 
};
struct h_asinh
{ 
    static Matrix eval(const Matrix& m)			{ return asinh(m); }; 
    template<class T>
    static Matrix eval_scal(T m)				{ return asinh(m); }; 
};
struct h_acosh
{ 
    static Matrix eval(const Matrix& m)			{ return acosh(m); }; 
    template<class T>
    static Matrix eval_scal(T m)				{ return acosh(m); }; 
};
struct h_acosh_c
{ 
    static Matrix eval(const Matrix& m)			{ return acosh_c(m); }; 
    template<class T>
    static Matrix eval_scal(T m)				{ return acosh_c(m); }; 
};
struct h_atanh
{ 
    static Matrix eval(const Matrix& m)			{ return atanh(m); }; 
    template<class T>
    static Matrix eval_scal(T m)				{ return atanh(m); }; 
};
struct h_atanh_c
{ 
    static Matrix eval(const Matrix& m)			{ return atanh_c(m); }; 
    template<class T>
    static Matrix eval_scal(T m)				{ return atanh_c(m); }; 
};
struct h_acoth
{ 
    static Matrix eval(const Matrix& m)			{ return acoth(m); }; 
    template<class T>
    static Matrix eval_scal(T m)				{ return acoth(m); }; 
};
struct h_acoth_c
{ 
    static Matrix eval(const Matrix& m)			{ return acoth_c(m); }; 
    template<class T>
    static Matrix eval_scal(T m)				{ return acoth_c(m); }; 
};
struct h_asech
{ 
    static Matrix eval(const Matrix& m)			{ return asech(m); }; 
    template<class T>
    static Matrix eval_scal(T m)				{ return asech(m); }; 
};
struct h_asech_c
{ 
    static Matrix eval(const Matrix& m)			{ return asech_c(m); }; 
    template<class T>
    static Matrix eval_scal(T m)				{ return asech_c(m); }; 
};
struct h_acsch
{ 
    static Matrix eval(const Matrix& m)			{ return acsch(m); }; 
    template<class T>
    static Matrix eval_scal(T m)				{ return acsch(m); }; 
};
struct h_real
{ 
    static Matrix eval(const Matrix& m)			{ return real(m); }; 
    template<class T>
    static Matrix eval_scal(T m)				{ return real(m); }; 
};
struct h_imag
{ 
    static Matrix eval(const Matrix& m)			{ return imag(m); }; 
    template<class T>
    static Matrix eval_scal(T m)				{ return imag(m); }; 
};
struct h_abs
{ 
    static Matrix eval(const Matrix& m)			{ return abs(m); }; 
    template<class T>
    static Matrix eval_scal(T m)				{ return abs(m); }; 
};
struct h_arg
{ 
    static Matrix eval(const Matrix& m)			{ return arg(m); }; 
    template<class T>
    static Matrix eval_scal(T m)				{ return arg(m); }; 
};
struct h_angle
{ 
    static Matrix eval(const Matrix& m)			{ return angle(m); }; 
    template<class T>
    static Matrix eval_scal(T m)				{ return angle(m); }; 
};
struct h_conj
{ 
    static Matrix eval(const Matrix& m)			{ return conj(m); }; 
    template<class T>
    static Matrix eval_scal(T m)				{ return conj(m); }; 
};
struct h_sqrt
{ 
    static Matrix eval(const Matrix& m)			{ return sqrt(m); }; 
    template<class T>
    static Matrix eval_scal(T m)				{ return sqrt(m); }; 
};
struct h_sqrt_c
{ 
    static Matrix eval(const Matrix& m)			{ return sqrt_c(m); }; 
    template<class T>
    static Matrix eval_scal(T m)				{ return sqrt_c(m); }; 
};
struct h_pow2
{ 
    static Matrix eval(const Matrix& m)			{ return pow2(m); }; 
    template<class T>
    static Matrix eval_scal(T m)				{ return pow2(m); }; 
};
struct h_exp
{ 
    static Matrix eval(const Matrix& m)			{ return exp(m); }; 
    template<class T>
    static Matrix eval_scal(T m)				{ return exp(m); }; 
};
struct h_log
{ 
    static Matrix eval(const Matrix& m)			{ return log(m); }; 
    template<class T>
    static Matrix eval_scal(T m)				{ return log(m); }; 
};
struct h_log_c
{ 
    static Matrix eval(const Matrix& m)			{ return log_c(m); }; 
    template<class T>
    static Matrix eval_scal(T m)				{ return log_c(m); }; 
};
struct h_log2
{ 
    static Matrix eval(const Matrix& m)			{ return log2(m); }; 
    template<class T>
    static Matrix eval_scal(T m)				{ return log2(m); }; 
};
struct h_log2_c
{ 
    static Matrix eval(const Matrix& m)			{ return log2_c(m); }; 
    template<class T>
    static Matrix eval_scal(T m)				{ return log2_c(m); }; 
};
struct h_log10
{ 
    static Matrix eval(const Matrix& m)			{ return log10(m); }; 
    template<class T>
    static Matrix eval_scal(T m)				{ return log10(m); }; 
};
struct h_log10_c
{ 
    static Matrix eval(const Matrix& m)			{ return log10_c(m); }; 
    template<class T>
    static Matrix eval_scal(T m)				{ return log10_c(m); }; 
};
struct h_floor
{ 
    static Matrix eval(const Matrix& m)			{ return floor(m); }; 
    template<class T>
    static Matrix eval_scal(T m)				{ return floor(m); }; 
};
struct h_ceil
{ 
    static Matrix eval(const Matrix& m)			{ return ceil(m); }; 
    template<class T>
    static Matrix eval_scal(T m)				{ return ceil(m); }; 
};
struct h_round
{ 
    static Matrix eval(const Matrix& m)			{ return round(m); }; 
    template<class T>
    static Matrix eval_scal(T m)				{ return round(m); }; 
};
struct h_fix
{ 
    static Matrix eval(const Matrix& m)			{ return fix(m); }; 
    template<class T>
    static Matrix eval_scal(T m)				{ return fix(m); }; 
};
struct h_trunc
{ 
    static Matrix eval(const Matrix& m)			{ return trunc(m); }; 
    template<class T>
    static Matrix eval_scal(T m)				{ return trunc(m); }; 
};
struct h_sqrt1pm1
{ 
    static Matrix eval(const Matrix& m)			{ return sqrt1pm1(m); }; 
    template<class T>
    static Matrix eval_scal(T m)				{ return sqrt1pm1(m); }; 
};
struct h_cbrt
{ 
    static Matrix eval(const Matrix& m)			{ return cbrt(real(m)); }; 
    template<class T>
    static Matrix eval_scal(T m)				{ return cbrt(real(m)); }; 
};
struct h_log1p
{ 
    static Matrix eval(const Matrix& m)			{ return log1p(m); }; 
    template<class T>
    static Matrix eval_scal(T m)				{ return log1p(m); }; 
};
struct h_sqrt1pm1_c
{ 
    static Matrix eval(const Matrix& m)			{ return sqrt1pm1_c(m); }; 
    template<class T>
    static Matrix eval_scal(T m)				{ return sqrt1pm1_c(m); }; 
};
struct h_log1p_c
{ 
    static Matrix eval(const Matrix& m)			{ return log1p_c(m); }; 
    template<class T>
    static Matrix eval_scal(T m)				{ return log1p_c(m); }; 
};
struct h_op_plus
{ 
    static Matrix eval(const Matrix& m)			{ return +(m); }; 
    template<class T>
    static Matrix eval_scal(T m)				{ return +(m); }; 
};
struct h_is_normal
{ 
    static Matrix eval(const Matrix& m)			{ return is_normal(m); }; 
    template<class T>
    static Matrix eval_scal(T m)				{ return Matrix(is_normal(m)); }; 
};
struct h_expm1
{ 
    static Matrix eval(const Matrix& m)			{ return expm1(m); }; 
    template<class T>
    static Matrix eval_scal(T m)				{ return expm1(m); }; 
};
struct h_signbit
{ 
    static Matrix eval(const Matrix& m)			{ return signbit(real(m)); }; 
    template<class T>
    static Matrix eval_scal(T m)				{ return Matrix(signbit(real(m))); }; 
};
struct h_logb
{ 
    static Matrix eval(const Matrix& m)			{ return logb(m); }; 
    template<class T>
    static Matrix eval_scal(T m)				{ return logb(m); }; 
};
struct h_ilogb
{ 
    static Matrix eval(const Matrix& m)			{ return ilogb(m); }; 
    template<class T>
    static Matrix eval_scal(T m)				{ return ilogb(m); }; 
};
struct h_exp2
{ 
    static Matrix eval(const Matrix& m)			{ return exp2(m); }; 
    template<class T>
    static Matrix eval_scal(T m)				{ return exp2(m); }; 
};
struct h_expi
{ 
    static Matrix eval(const Matrix& m)			{ return expi(m); }; 
    template<class T>
    static Matrix eval_scal(T m)				{ return expi(m); }; 
};

struct h_exp10
{ 
    static Matrix eval(const Matrix& m)			{ return exp10(m); }; 
    template<class T>
    static Matrix eval_scal(T m)				{ return exp10(m); }; 
};

using tf    = unary_function;
using ccp   = const char*;

static const char* str_isstrue = "is_scalar_true";
static const char* str_issfals = "is_scalar_false";

struct function_sqrt1pm1{ ccp name() { return "sqrt1pm1";};	scal_func<h_sqrt1pm1> f;tf* function() { return &f;}; };
struct function_cbrt    { ccp name() { return "cbrt";};	    scal_func<h_cbrt> f;    tf* function() { return &f;}; };
struct function_log1p   { ccp name() { return "log1p";};	scal_func<h_log1p> f;    tf* function() { return &f;}; };
struct function_sqrt1pm1_c{ ccp name() { return "sqrt1pm1_c";};scal_func<h_sqrt1pm1_c> f;tf* function() { return &f;}; };
struct function_log1p_c { ccp name() { return "log1p_c";}; scal_func<h_log1p_c> f;tf* function() { return &f;}; };
struct function_op_plus { ccp name() { return "op_plus";};  scal_func<h_op_plus> f; tf* function() { return &f;}; };
struct function_is_normal{ ccp name(){ return "is_normal";};scal_func<h_is_normal> f;tf* function() { return &f;}; };
struct function_expm1   { ccp name() { return "expm1";};    scal_func<h_expm1> f;   tf* function() { return &f;}; };
struct function_signbit { ccp name() { return "signbit";};  scal_func<h_signbit> f; tf* function() { return &f;}; };
struct function_logb    { ccp name() { return "logb";};     scal_func<h_logb> f;    tf* function() { return &f;}; };
struct function_ilogb   { ccp name() { return "ilogb";};    scal_func<h_ilogb> f;   tf* function() { return &f;}; };
struct function_exp2    { ccp name() { return "exp2";};     scal_func<h_exp2> f;    tf* function() { return &f;}; };
struct function_exp10   { ccp name() { return "exp10";};    scal_func<h_exp10> f;   tf* function() { return &f;}; };
struct function_expi    { ccp name() { return "expi";};     scal_func<h_expi> f;    tf* function() { return &f;}; };

struct function_real	{ ccp name() { return "real";};		scal_func<h_real> f;	tf* function() { return &f;}; };
struct function_imag	{ ccp name() { return "imag";};		scal_func<h_imag> f;	tf* function() { return &f;}; };			
struct function_abs		{ ccp name() { return "abs";};		scal_func<h_abs> f;		tf* function() { return &f;}; };
struct function_arg		{ ccp name() { return "arg";};		scal_func<h_arg> f;		tf* function() { return &f;}; };
struct function_angle	{ ccp name() { return "angle";};	scal_func<h_angle> f;	tf* function() { return &f;}; };
struct function_conj	{ ccp name() { return "conj";};		scal_func<h_conj> f;	tf* function() { return &f;}; };			
struct function_sqrt	{ ccp name() { return "sqrt";};		scal_func<h_sqrt> f;	tf* function() { return &f;}; };
struct function_pow2	{ ccp name() { return "pow2";};		scal_func<h_pow2> f;	tf* function() { return &f;}; };			
struct function_exp		{ ccp name() { return "exp";};		scal_func<h_exp> f;		tf* function() { return &f;}; };
struct function_log		{ ccp name() { return "log";};		scal_func<h_log> f;		tf* function() { return &f;}; };
struct function_log2	{ ccp name() { return "log2";};		scal_func<h_log2> f;	tf* function() { return &f;}; };
struct function_log10	{ ccp name() { return "log10";};	scal_func<h_log10> f;	tf* function() { return &f;}; };			
struct function_floor	{ ccp name() { return "floor";};	scal_func<h_floor> f;	tf* function() { return &f;}; };
struct function_ceil	{ ccp name() { return "ceil";};		scal_func<h_ceil> f;	tf* function() { return &f;}; };			
struct function_round	{ ccp name() { return "round";};	scal_func<h_round> f;	tf* function() { return &f;}; };
struct function_fix		{ ccp name() { return "fix";};		scal_func<h_fix> f;		tf* function() { return &f;}; };
struct function_trunc	{ ccp name() { return "trunc";};	scal_func<h_trunc> f;	tf* function() { return &f;}; };

struct function_sign	{ ccp name() { return "sign";};		scal_func<h_sign> f;	tf* function() { return &f;}; };			
struct function_isign	{ ccp name() { return "isign";};	scal_func<h_isign> f;	tf* function() { return &f;}; };
struct function_sin		{ ccp name() { return "sin";};		scal_func<h_sin> f;		tf* function() { return &f;}; };
struct function_cos		{ ccp name() { return "cos";};		scal_func<h_cos> f;		tf* function() { return &f;}; };
struct function_tan		{ ccp name() { return "tan";};		scal_func<h_tan> f;		tf* function() { return &f;}; };
struct function_cot		{ ccp name() { return "cot";};		scal_func<h_cot> f;		tf* function() { return &f;}; };
struct function_sec		{ ccp name() { return "sec";};		scal_func<h_sec> f;		tf* function() { return &f;}; };
struct function_csc		{ ccp name() { return "csc";};		scal_func<h_csc> f;		tf* function() { return &f;}; };
struct function_asin	{ ccp name() { return "asin";};		scal_func<h_asin> f;	tf* function() { return &f;}; };
struct function_acos	{ ccp name() { return "acos";};		scal_func<h_acos> f;	tf* function() { return &f;}; };			
struct function_atan	{ ccp name() { return "atan";};		scal_func<h_atan> f;	tf* function() { return &f;}; };
struct function_acot	{ ccp name() { return "acot";};		scal_func<h_acot> f;	tf* function() { return &f;}; };			
struct function_asec	{ ccp name() { return "asec";};		scal_func<h_asec> f;	tf* function() { return &f;}; };
struct function_acsc	{ ccp name() { return "acsc";};		scal_func<h_acsc> f;	tf* function() { return &f;}; };			
struct function_sinh	{ ccp name() { return "sinh";};		scal_func<h_sinh> f;	tf* function() { return &f;}; };
struct function_cosh	{ ccp name() { return "cosh";};		scal_func<h_cosh> f;	tf* function() { return &f;}; };			
struct function_tanh	{ ccp name() { return "tanh";};		scal_func<h_tanh> f;	tf* function() { return &f;}; };
struct function_coth	{ ccp name() { return "coth";};		scal_func<h_coth> f;	tf* function() { return &f;}; };			
struct function_sech	{ ccp name() { return "sech";};		scal_func<h_sech> f;	tf* function() { return &f;}; };
struct function_csch	{ ccp name() { return "csch";};		scal_func<h_csch> f;	tf* function() { return &f;}; };			
struct function_asinh	{ ccp name() { return "asinh";};	scal_func<h_asinh> f;	tf* function() { return &f;}; };
struct function_acosh	{ ccp name() { return "acosh";};	scal_func<h_acosh> f;	tf* function() { return &f;}; };			
struct function_atanh	{ ccp name() { return "atanh";};	scal_func<h_atanh> f;	tf* function() { return &f;}; };
struct function_acoth	{ ccp name() { return "acoth";};	scal_func<h_acoth> f;	tf* function() { return &f;}; };			
struct function_asech	{ ccp name() { return "asech";};	scal_func<h_asech> f;	tf* function() { return &f;}; };
struct function_acsch	{ ccp name() { return "acsch";};	scal_func<h_acsch> f;	tf* function() { return &f;}; };

struct function_ifloor	{ ccp name() { return "ifloor";};	scal_func<h_ifloor> f;	tf* function() { return &f;}; };
struct function_iceil	{ ccp name() { return "iceil";};	scal_func<h_iceil> f;	tf* function() { return &f;}; };			
struct function_iround	{ ccp name() { return "iround";};	scal_func<h_iround> f;	tf* function() { return &f;}; };
struct function_ifix	{ ccp name() { return "ifix";};		scal_func<h_ifix> f;	tf* function() { return &f;}; };
struct function_itrunc	{ ccp name() { return "itrunc";};	scal_func<h_itrunc> f;	tf* function() { return &f;}; };
struct function_isinf	{ ccp name() { return "is_inf";};	scal_func<h_isinf> f;	tf* function() { return &f;}; };
struct function_isnan	{ ccp name() { return "is_nan";};	scal_func<h_isnan> f;	tf* function() { return &f;}; };
struct function_isfin	{ ccp name() { return "is_finite";};scal_func<h_isfin> f;	tf* function() { return &f;}; };
struct function_isregul { ccp name() { return "is_regular";};scal_func<h_isreg> f;	tf* function() { return &f;}; };
struct function_isnorm  { ccp name() { return "is_normal";};scal_func<h_isnorm> f;	tf* function() { return &f;}; };
struct function_isint   { ccp name() { return "is_int";};   scal_func<h_isint> f;	tf* function() { return &f;}; };
struct function_isreal  { ccp name() { return "is_real";};  scal_func<h_isreal> f;	tf* function() { return &f;}; };
struct function_isstrue	{ ccp name() { return str_isstrue;};scal_func<h_isstrue> f;	tf* function() { return &f;}; };
struct function_issfalse{ ccp name() { return str_issfals;};scal_func<h_issfalse> f;tf* function() { return &f;}; };
struct function_sqrt_c	{ ccp name() { return "sqrt_c";};	scal_func<h_sqrt_c> f;	tf* function() { return &f;}; };
struct function_log_c	{ ccp name() { return "log_c";};	scal_func<h_log_c> f;	tf* function() { return &f;}; };
struct function_log2_c	{ ccp name() { return "log2_c";};	scal_func<h_log2_c> f;	tf* function() { return &f;}; };
struct function_log10_c { ccp name() { return "log10_c";};	scal_func<h_log10_c> f;tf* function() { return &f;}; };			
struct function_asin_c	{ ccp name() { return "asin_c";};	scal_func<h_asin_c> f;	tf* function() { return &f;}; };
struct function_acos_c	{ ccp name() { return "acos_c";};	scal_func<h_acos_c> f;	tf* function() { return &f;}; };			
struct function_asec_c	{ ccp name() { return "asec_c";};	scal_func<h_asec_c> f;	tf* function() { return &f;}; };
struct function_acsc_c	{ ccp name() { return "acsc_c";};	scal_func<h_acsc_c> f;	tf* function() { return &f;}; };			
struct function_acosh_c { ccp name() { return "acosh_c";};	scal_func<h_acosh_c> f;tf* function() { return &f;}; };			
struct function_atanh_c { ccp name() { return "atanh_c";};	scal_func<h_atanh_c> f;tf* function() { return &f;}; };
struct function_acoth_c { ccp name() { return "acoth_c";};	scal_func<h_acoth_c> f;tf* function() { return &f;}; };			
struct function_asech_c { ccp name() { return "asech_c";};	scal_func<h_asech_c> f;tf* function() { return &f;}; };

};};
