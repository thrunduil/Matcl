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

namespace 
{
    Real to_real(Integer val)
    {
        return val;
    };

    template<class T>
    T to_real(T val)
    {
        return val;
    };
}

class bin_functions_list
{
    public:
        using matrix_pair   = matrix_set_bin::matrix_pair;
        using scalar_pair   = matrix_set_bin::scalar_pair;

    private:
        const matrix_set_bin&	m_tests;
        options					m_options;

    public:
        dynamic_mat_set&        m_ds;

    public:
        bin_functions_list(const matrix_set_bin& t, dynamic_mat_set& ds) 
                : m_tests(t), m_ds(ds)
        {};

        void			make(options opts, Integer thread_id);
        void			make_mult(options opts);
        void			make_kron(options opts);

        matrix_pair		get_matrix(int code) const;
        scalar_pair		get_scalar(int code) const;

    private:

        void    test_op_plus();
        void    test_op_minus();		void    test_op_mult();
        void    test_op_div();			void    test_plus();
        void    test_minus();			void    test_mul();
        void    test_div();				void    test_pow();
        void    test_max_bin();			void    test_min_bin();
        void    test_xor();				void    test_rem();
        void    test_mod();				void    test_kron();
        void    test_pow_c();           void    test_div_0();
        void    test_div_1();

        void    test_fma();             void    test_fms();
        void    test_dot2_ac();
        
        void    test_op_or();
        void    test_op_and();			void    test_op_eeq();
        void    test_op_neq();			void    test_op_lt();
        void    test_op_leq();			void    test_op_gt();
        void    test_op_geq();			void    test_idiv();
        void    test_beta();            void    test_eval_bin();
        void    test_atan2();           void    test_hypot();
        void    test_copysign();        void    test_fdim();
        void    test_nextafter();       void    test_powm1();
        void    test_eeq_nan();         void    test_neq_nan();

        void    test_chain_mult();
        void    test_mmul();            void    test_mmul_ext();
        void    test_gemm();            void    test_mmul_abs();
        void    test_gemm_sub();
        
    private:
        bin_functions_list(const bin_functions_list&) = delete;
        bin_functions_list& operator=(const bin_functions_list&) = delete;
};

class test_function_op_plus : public bin_function
{
    public:
        template<class T1,class T2>
        static Real eval_scal_func(T1 a, T2 b)
        {
            using Compl_1   = typename md::complex_type<T1>::type;
            using Compl_2   = typename md::complex_type<T2>::type;

            Real dif	= norm_1((a+b) - (full(a)+full(b)));
            dif			+=norm_1((a+b) - (Compl_1(a)+(b)));
            dif			+=norm_1((a+b) - ((a)+Compl_2(b)));
            dif			+=norm_1((a+b) - (Compl_1(a)+Compl_2(b)));

            Matrix res  = abs(a) + abs(b);

            if (dif < constants::eps(res.get_value_code()) * norm_1(res) * 10.)
                dif     = 0.0;

            return dif;
        };

    public:
        virtual Real	eval_mat(const Matrix& mat1,const Matrix& mat2,int code);
        virtual Real	eval_scalar(const Scalar& s1, const Scalar& s2,int code);
};

class test_function_op_minus : public bin_function
{
    public:
        template<class T1,class T2>
        static Real eval_scal_func(T1 a, T2 b)
        {
            using Compl_1   = typename md::complex_type<T1>::type;
            using Compl_2   = typename md::complex_type<T2>::type;

            Real dif	= norm_1((a-b) - (full(a)-full(b)));
            dif			+=norm_1((a-b) - (Compl_1(a)-(b)));
            dif			+=norm_1((a-b) - ((a)-Compl_2(b)));
            dif			+=norm_1((a-b) - (Compl_1(a)-Compl_2(b)));

            Matrix res  = abs(a) + abs(b);

            if (dif < constants::eps(res.get_value_code()) * norm_1(res) * 10.)
                dif     = 0.0;

            return dif;
        };

    public:
        virtual Real	eval_mat(const Matrix& mat1,const Matrix& mat2,int code);
        virtual Real	eval_scalar(const Scalar& s1, const Scalar& s2,int code);
};

class test_function_op_mult : public bin_function
{
    public:
        template<class T1,class T2>
        static Real eval_scal_func(T1 a, T2 b)
        {
            using Compl_1   = typename md::complex_type<T1>::type;
            using Compl_2   = typename md::complex_type<T2>::type;

            Real dif	= norm_1((a*b) - (full(a)*full(b)));
            dif			+=norm_1((a*b) - (Compl_1(a)*(b)));
            dif			+=norm_1((a*b) - ((a)*Compl_2(b)));
            dif			+=norm_1((a*b) - (Compl_1(a)*Compl_2(b)));

            Matrix res  = abs(a) + abs(b) + abs(a*b);

            if (dif < constants::eps(res.get_value_code()) * norm_1(res) * 10.)
                dif     = 0.0;

            return dif;
        };

    public:
        virtual Real	eval_mat(const Matrix& mat1,const Matrix& mat2,int code);
        virtual Real	eval_scalar(const Scalar& s1, const Scalar& s2,int code);
};

class test_function_mmul : public bin_function
{
    private:
        trans_type      m_ta;
        trans_type      m_tb;

    public:
        test_function_mmul(trans_type ta, trans_type tb)
            :m_ta(ta), m_tb(tb)
        {};

        template<class T1,class T2>
        Real eval_scal_func(T1 a, T2 b)
        {
            using Compl_1   = typename md::complex_type<T1>::type;
            using Compl_2   = typename md::complex_type<T2>::type;

            Matrix out_full = trans(a, m_ta) * trans(b, m_tb);
            Matrix out		= mmul(a, b, m_ta, m_tb);

            Real dif	= norm_1(out_full - out);
            Matrix res  = abs(a) + abs(b) + abs(a*b);

            if (dif < constants::eps(res.get_value_code()) * norm_1(res) * 10.)
                dif     = 0.0;

            return dif;
        };
    public:
        virtual Real	eval_mat(const Matrix& mat1,const Matrix& mat2,int code);
        virtual Real	eval_scalar(const Scalar& s1, const Scalar& s2,int code);
};

class test_function_mmul_abs: public bin_function
{
    private:
        trans_type          m_ta;
        dynamic_mat_set&    m_ds;

    public:
        test_function_mmul_abs(trans_type ta, dynamic_mat_set& ds)
            :m_ta(ta), m_ds(ds)
        {};

    public:
        virtual Real	eval_mat(const Matrix& mat1,const Matrix& mat2,int code);
        virtual Real	eval_scalar(const Scalar& s1, const Scalar& s2,int code);

        static Matrix   rand_mat_C(Integer r, Integer c, Integer seed, dynamic_mat_set& ds);
};

class test_function_gemm : public bin_function
{
    private:
        trans_type          m_ta;
        trans_type          m_tb;
        dynamic_mat_set&    m_ds;

    public:
        test_function_gemm(trans_type ta, trans_type tb, dynamic_mat_set& ds)
            :m_ta(ta), m_tb(tb), m_ds(ds)
        {};

    public:
        virtual Real	eval_mat(const Matrix& mat1,const Matrix& mat2,int code);
        virtual Real	eval_scalar(const Scalar& s1, const Scalar& s2,int code);

        static Matrix       rand_scalar(Integer seed);
        static value_code   get_value_code(const Matrix& alpha, const Matrix& beta, const Matrix& mat1, 
                                       const Matrix& mat2);
        static Matrix       rand_mat_C(Integer r, Integer c, value_code vc_C, Integer seed,
                                       dynamic_mat_set& ds);
        static size_t       rand_int(Integer seed);
};

class test_function_gemm_sub : public bin_function
{
    private:
        trans_type          m_ta;
        trans_type          m_tb;
        dynamic_mat_set&    m_ds;

    public:
        test_function_gemm_sub(trans_type ta, trans_type tb, dynamic_mat_set& ds)
            :m_ta(ta), m_tb(tb), m_ds(ds)
        {};

    public:
        virtual Real	eval_mat(const Matrix& mat1,const Matrix& mat2,int code);
        virtual Real	eval_scalar(const Scalar& s1, const Scalar& s2,int code);

        static Integer  rand_colon_ver(Integer rs, Integer re, Integer cs, Integer ce,
                                       Integer seed);
};

class test_function_mmul_ext : public bin_function
{
    private:
        trans_type_ext  m_ta;
        trans_type_ext  m_tb;

    public:
        test_function_mmul_ext(trans_type_ext ta, trans_type_ext tb)
            :m_ta(ta), m_tb(tb)
        {};

        template<class T1,class T2>
        Real eval_scal_func(T1 a, T2 b)
        {
            using Compl_1   = typename md::complex_type<T1>::type;
            using Compl_2   = typename md::complex_type<T2>::type;

            Matrix out_full = trans(a, m_ta) * trans(b, m_tb);
            Matrix out		= mmul(a, b, m_ta, m_tb);

            Real dif	= norm_1(out_full - out);
            Matrix res  = abs(a) + abs(b) + abs(a*b);

            if (dif < constants::eps(res.get_value_code()) * norm_1(res) * 10.)
                dif     = 0.0;

            return dif;
        };

    public:
        virtual Real	eval_mat(const Matrix& mat1,const Matrix& mat2,int code);
        virtual Real	eval_scalar(const Scalar& s1, const Scalar& s2,int code);
};
class test_function_chain_mult : public bin_function
{
    public:
        virtual Real	eval_mat(const Matrix& mat1,const Matrix& mat2,int code);
        virtual Real	eval_scalar(const Scalar& s1, const Scalar& s2,int code);
};

class test_function_op_div : public bin_function
{
    public:
        template<class T1,class T2>
        static Real eval_scal_func(T1 a, T2 b)
        {
            using Compl_1   = typename md::complex_type<T1>::type;
            using Compl_2   = typename md::complex_type<T2>::type;

            Real dif	= norm_1(div(a,b) - div(full(a),full(b)));;
            dif			+=norm_1(div(a,b) - div(Compl_1(a),(b)));;
            dif			+=norm_1(div(a,b) - div((a),Compl_2(b)));;
            dif			+=norm_1(div(a,b) - div(Compl_1(a),Compl_2(b)));;

            Matrix res  = abs(a) + abs(b) + abs(div(a,b));

            if (dif < constants::eps(res.get_value_code()) * norm_1(res) * 10.)
                dif     = 0.0;

            return dif;
        };

    public:
        virtual Real	eval_mat(const Matrix& mat1,const Matrix& mat2,int code);
        virtual Real	eval_scalar(const Scalar& s1, const Scalar& s2,int code);
};

class test_function_plus : public bin_function
{
    public:
        template<class T1,class T2>
        static Real eval_scal_func(T1 a, T2 b)
        {
            using Compl_1   = typename md::complex_type<T1>::type;
            using Compl_2   = typename md::complex_type<T2>::type;

            Real dif	= norm_1(plus(a,b) - plus(full(a),full(b)));;
            dif			+=norm_1(plus(a,b) - plus(Compl_1(a),(b)));;
            dif			+=norm_1(plus(a,b) - plus((a),Compl_2(b)));;
            dif			+=norm_1(plus(a,b) - plus(Compl_1(a),Compl_2(b)));;

            Matrix res  = abs(a) + abs(b);

            if (dif < constants::eps(res.get_value_code()) * norm_1(res) * 10.)
                dif     = 0.0;

            return dif;
        };

    public:
        virtual Real	eval_mat(const Matrix& mat1,const Matrix& mat2,int code);
        virtual Real	eval_scalar(const Scalar& s1, const Scalar& s2,int code);
};

class test_function_minus : public bin_function
{
    public:
        template<class T1,class T2>
        static Real eval_scal_func(T1 a, T2 b)
        {
            using Compl_1   = typename md::complex_type<T1>::type;
            using Compl_2   = typename md::complex_type<T2>::type;

            Real dif	= norm_1(minus(a,b) - minus(full(a),full(b)));;
            dif			+=norm_1(minus(a,b) - minus(Compl_1(a),(b)));;
            dif			+=norm_1(minus(a,b) - minus((a),Compl_2(b)));;
            dif			+=norm_1(minus(a,b) - minus(Compl_1(a),Compl_2(b)));;

            Matrix res  = abs(a) + abs(b);

            if (dif < constants::eps(res.get_value_code()) * norm_1(res) * 10.)
                dif     = 0.0;

            return dif;
        };

    public:
        virtual Real	eval_mat(const Matrix& mat1,const Matrix& mat2,int code);
        virtual Real	eval_scalar(const Scalar& s1, const Scalar& s2,int code);
};

class test_function_mul : public bin_function
{
    public:
        template<class T1,class T2>
        static Real eval_scal_func(T1 a, T2 b)
        {
            using Compl_1   = typename md::complex_type<T1>::type;
            using Compl_2   = typename md::complex_type<T2>::type;

            Real dif	= norm_1(mul(a,b) - mul(full(a),full(b)));;
            dif			+=norm_1(mul(a,b) - mul(Compl_1(a),(b)));;
            dif			+=norm_1(mul(a,b) - mul((a),Compl_2(b)));;
            dif			+=norm_1(mul(a,b) - mul(Compl_1(a),Compl_2(b)));;

            Matrix res  = abs(a) + abs(b) + abs(a*b);

            if (dif < constants::eps(res.get_value_code()) * norm_1(res) * 10.)
                dif     = 0.0;

            return dif;
        };

    public:
        virtual Real	eval_mat(const Matrix& mat1,const Matrix& mat2,int code);
        virtual Real	eval_scalar(const Scalar& s1, const Scalar& s2,int code);
};

class test_function_div : public bin_function
{
    public:
        template<class T1,class T2>
        static Real eval_scal_func(T1 a, T2 b)
        {
            using Compl_1   = typename md::complex_type<T1>::type;
            using Compl_2   = typename md::complex_type<T2>::type;

            Real dif	= norm_1(div(a,b) - div(full(a),full(b)));;
            dif			+=norm_1(div(a,b) - div(Compl_1(a),(b)));;
            dif			+=norm_1(div(a,b) - div((a),Compl_2(b)));;
            dif			+=norm_1(div(a,b) - div(Compl_1(a),Compl_2(b)));;

            Matrix res  = abs(a) + abs(b) + abs(div(a,b));

            if (dif < constants::eps(res.get_value_code()) * norm_1(res) * 10.)
                dif     = 0.0;

            return dif;
        };

    public:
        virtual Real	eval_mat(const Matrix& mat1,const Matrix& mat2,int code);
        virtual Real	eval_scalar(const Scalar& s1, const Scalar& s2,int code);
};

class test_function_div_0 : public bin_function
{
    public:
        template<class T1,class T2>
        static Real eval_scal_func(T1 a, T2 b)
        {
            using Compl_1   = typename md::complex_type<T1>::type;
            using Compl_2   = typename md::complex_type<T2>::type;

            Real dif	= norm_1(div_0(a,b) - div_0(full(a),full(b)));;
            dif			+=norm_1(div_0(a,b) - div_0(Compl_1(a),(b)));;
            dif			+=norm_1(div_0(a,b) - div_0((a),Compl_2(b)));;
            dif			+=norm_1(div_0(a,b) - div_0(Compl_1(a),Compl_2(b)));;

            Matrix res  = abs(a) + abs(b) + abs(div_0(a,b));

            if (dif < constants::eps(res.get_value_code()) * norm_1(res) * 10.)
                dif     = 0.0;

            return dif;
        };

    public:
        virtual Real	eval_mat(const Matrix& mat1,const Matrix& mat2,int code);
        virtual Real	eval_scalar(const Scalar& s1, const Scalar& s2,int code);
};

class test_function_div_1 : public bin_function
{
    public:
        template<class T1,class T2>
        static Real eval_scal_func(T1 a, T2 b)
        {
            using Compl_1   = typename md::complex_type<T1>::type;
            using Compl_2   = typename md::complex_type<T2>::type;

            Real dif	= norm_1(div_1(a,b) - div_1(full(a),full(b)));;
            dif			+=norm_1(div_1(a,b) - div_1(Compl_1(a),(b)));;
            dif			+=norm_1(div_1(a,b) - div_1((a),Compl_2(b)));;
            dif			+=norm_1(div_1(a,b) - div_1(Compl_1(a),Compl_2(b)));;

            Matrix res  = abs(a) + abs(b) + abs(div_1(a,b));

            if (dif < constants::eps(res.get_value_code()) * norm_1(res) * 10.)
                dif     = 0.0;

            return dif;
        };

    public:
        virtual Real	eval_mat(const Matrix& mat1,const Matrix& mat2,int code);
        virtual Real	eval_scalar(const Scalar& s1, const Scalar& s2,int code);
};

class test_function_idiv : public bin_function
{
    public:
        template<class T1,class T2>
        static Real eval_scal_func(T1 a, T2 b)
        {
            using Compl_1   = typename md::complex_type<T1>::type;
            using Compl_2   = typename md::complex_type<T2>::type;

            Real dif	= norm_1(idiv(a,b) - idiv(full(a),full(b)));;
            dif			+=norm_1(idiv(to_real(a),b) - idiv(Compl_1(a),(b)));;
            dif			+=norm_1(idiv(a,to_real(b)) - idiv((a),Compl_2(b)));;
            dif			+=norm_1(idiv(to_real(a),to_real(b)) - idiv(Compl_1(a),Compl_2(b)));;

            Matrix res  = abs(a) + abs(b) + abs(idiv(a,b));

            if (dif < constants::eps(res.get_value_code()) * norm_1(res) * 10.)
                dif     = 0.0;

            return dif;
        };

    public:
        virtual Real	eval_mat(const Matrix& mat1,const Matrix& mat2,int code);
        virtual Real	eval_scalar(const Scalar& s1, const Scalar& s2,int code);
};

class test_function_beta : public bin_function
{
    public:
        template<class T1,class T2>
        static Real eval_scal_func(T1 a, T2 b)
        {
            using Compl_1   = typename md::complex_type<T1>::type;
            using Compl_2   = typename md::complex_type<T2>::type;

            Real dif = 0.;
            // Note reversed order of args (test for symmetricity).
            dif	+= norm_1(matcl::beta(b,a) - matcl::beta(full(a),full(b)));
            dif	+= norm_1(matcl::beta(b,a) - matcl::beta(Compl_1(a),(b)));
            dif	+= norm_1(matcl::beta(b,a) - matcl::beta((a),Compl_2(b)));
            dif	+= norm_1(matcl::beta(b,a) - matcl::beta(Compl_1(a),Compl_2(b)));

            Matrix res  = abs(a) + abs(b) + abs(beta(a,b));

            if (dif < constants::eps(res.get_value_code()) * res * 1.0e4)
                dif     = 0.0;

            return dif;
        };

    public:
        virtual Real	eval_mat(const Matrix& mat1,const Matrix& mat2,int code);
        virtual Real	eval_scalar(const Scalar& s1, const Scalar& s2,int code);
};

struct binfunc_plus
{
    template<class T1, class T2>
    static auto eval(const T1& A, const T2& B)
                    -> typename md::unify_types<T1,T2>::type
    {
        using ret = typename md::unify_types<T1,T2>::type;
        return ret(A) + ret(B);
    };
};

class test_function_eval_bin : public bin_function
{
    public:
        virtual Real	eval_mat(const Matrix& mat1,const Matrix& mat2,int code);
        virtual Real	eval_scalar(const Scalar& s1, const Scalar& s2,int code);
};

class test_function_atan2 : public bin_function
{
    public:
        template<class T1,class T2>
        static Real eval_scal_func(T1 a, T2 b)
        {
            using Compl_1   = typename md::complex_type<T1>::type;
            using Compl_2   = typename md::complex_type<T2>::type;

            Real dif = 0.;
            dif			+=norm_1(matcl::atan2(real(a),real(b)) 
                                 - matcl::atan2(full(real(a)),full(real(b))));;

            Matrix res  = abs(a) + abs(b) + abs(atan2(real(a),real(b)));

            if (dif < constants::eps(res.get_value_code()) * norm_1(res) * 10.)
                dif     = 0.0;

            return dif;
        };

        template<class T1>
        static bool is_complex(const T1&)  
        { 
            return false; 
        };

        static bool is_complex(const Complex& val)
        { 
            return imag(val) != 0; 
        };
        static bool is_complex(const Float_complex& val)
        { 
            return imag(val) != 0; 
        };

    public:
        virtual Real	eval_mat(const Matrix& mat1,const Matrix& mat2,int code);
        virtual Real	eval_scalar(const Scalar& s1, const Scalar& s2,int code);
};

class test_function_hypot : public bin_function
{
    public:
        template<class T1,class T2>
        static Real eval_scal_func(T1 a, T2 b)
        {
            using Compl_1   = typename md::complex_type<T1>::type;
            using Compl_2   = typename md::complex_type<T2>::type;

            Real dif = 0.;
            dif			+=norm_1(matcl::hypot(a,b) - matcl::hypot(full(a),full(b)));;
            dif			+=norm_1(matcl::hypot(a,b) - matcl::hypot(Compl_1(a),(b)));;
            dif			+=norm_1(matcl::hypot(a,b) - matcl::hypot((a),Compl_2(b)));;
            dif			+=norm_1(matcl::hypot(a,b) - matcl::hypot(Compl_1(a),Compl_2(b)));;

            Matrix res  = abs(a) + abs(b) + abs(hypot(a,b));

            if (dif < constants::eps(res.get_value_code()) * norm_1(res) * 10.)
                dif     = 0.0;

            return dif;
        };

    public:
        virtual Real	eval_mat(const Matrix& mat1,const Matrix& mat2,int code);
        virtual Real	eval_scalar(const Scalar& s1, const Scalar& s2,int code);
};

class test_function_pow : public bin_function
{
    public:
        template<class T1,class T2>
        static Real eval_scal_func(T1 a, T2 b)
        {
            using Compl_1   = typename md::complex_type<T1>::type;
            using Compl_2   = typename md::complex_type<T2>::type;

            Matrix base = matcl::pow(a,b);
            Matrix res1 = matcl::pow(full(a),full(b));
            Matrix res2 = matcl::pow(Compl_1(a),(b));

            Real dif	= norm_1(base - res1);;
            dif			+=norm_1(base - res2);

            Real nr     = norm_1(res1)  + abs(a) + abs(b);

            if (dif < constants::eps(base.get_value_code()) * nr * 10000.0)
                dif     = 0.0;

            Matrix res3 = matcl::pow(Compl_1(a),Compl_2(b));
            Matrix res4 = matcl::pow((a),Compl_2(b));

            Real dif2   = norm_1(base - res3);
            dif2		+=norm_1(base - res4);
            Real nr2    = norm_1(res3) + abs(a) + abs(b);

            if (dif2 < constants::eps(base.get_value_code()) * nr2 * 10000.)
                dif2     = 0.0;

            dif         += dif2;
            dif         += test_function_pow_c::compare_value_type(base, res1);

            return dif;
        };

    public:
        virtual Real	eval_mat(const Matrix& mat1,const Matrix& mat2,int code);
        virtual Real	eval_scalar(const Scalar& s1, const Scalar& s2,int code);
};

class test_function_pow_c : public bin_function
{
    public:
        template<class T1,class T2>
        static Real eval_scal_func(T1 a, T2 b)
        {
            using Compl_1   = typename md::complex_type<T1>::type;
            using Compl_2   = typename md::complex_type<T2>::type;

            Matrix base = matcl::pow_c(a,b);
            Matrix res1 = matcl::pow_c(full(a),full(b));
            Matrix res2 = matcl::pow_c(Compl_1(a),(b));

            Real dif	= norm_1(base - res1);;
            dif			+=norm_1(base - res2);

            Real nr     = norm_1(res1)  + abs(a) + abs(b);

            if (dif < constants::eps(base.get_value_code()) * nr * 10000.0)
                dif     = 0.0;

            Matrix res3 = matcl::pow_c(Compl_1(a),Compl_2(b));
            Matrix res4 = matcl::pow_c((a),Compl_2(b));

            Real dif2   = norm_1(base - res3);
            dif2		+=norm_1(base - res4);
            Real nr2    = norm_1(res3) + abs(a) + abs(b);

            if (dif2 < constants::eps(base.get_value_code()) * nr2 * 10000.)
                dif2     = 0.0;

            dif         += dif2;
            dif         += compare_value_type(base, res1);

            return dif;
        };

    public:
        virtual Real	eval_mat(const Matrix& mat1,const Matrix& mat2,int code);
        virtual Real	eval_scalar(const Scalar& s1, const Scalar& s2,int code);

    public:
        static Real     check_value_type_pow(const Matrix& mat1, const Matrix& mat2, const Matrix& out);
        static Real     compare_value_type(const Matrix& mat1, const Matrix& mat2);
};

class test_function_max_bin : public bin_function
{
    public:
        template<class T1,class T2>
        static Real eval_scal_func(T1 a, T2 b)
        {
            using Compl_1   = typename md::complex_type<T1>::type;
            using Compl_2   = typename md::complex_type<T2>::type;

            Real dif	= norm_1(matcl::max(a,b) - matcl::max(full(a),full(b)));;
            dif			+=norm_1(matcl::max(a,b) - matcl::max(Compl_1(a),(b)));;
            dif			+=norm_1(matcl::max(a,b) - matcl::max((a),Compl_2(b)));;
            dif			+=norm_1(matcl::max(a,b) - matcl::max(Compl_1(a),Compl_2(b)));;

            Matrix res  = abs(a) + abs(b);

            if (dif < constants::eps(res.get_value_code()) * norm_1(res) * 10.)
                dif     = 0.0;

            return dif;
        };

    public:
        virtual Real	eval_mat(const Matrix& mat1,const Matrix& mat2,int code);
        virtual Real	eval_scalar(const Scalar& s1, const Scalar& s2,int code);
};

class test_function_min_bin : public bin_function
{
    public:
        template<class T1,class T2>
        static Real eval_scal_func(T1 a, T2 b)
        {
            using Compl_1   = typename md::complex_type<T1>::type;
            using Compl_2   = typename md::complex_type<T2>::type;

            Real dif	= norm_1(matcl::min(a,b) - matcl::min(full(a),full(b)));;
            dif			+=norm_1(matcl::min(a,b) - matcl::min(Compl_1(a),(b)));;
            dif			+=norm_1(matcl::min(a,b) - matcl::min((a),Compl_2(b)));;
            dif			+=norm_1(matcl::min(a,b) - matcl::min(Compl_1(a),Compl_2(b)));;

            Matrix res  = abs(a) + abs(b);

            if (dif < constants::eps(res.get_value_code()) * norm_1(res) * 10.)
                dif     = 0.0;

            return dif;
        };

    public:
        virtual Real	eval_mat(const Matrix& mat1,const Matrix& mat2,int code);
        virtual Real	eval_scalar(const Scalar& s1, const Scalar& s2,int code);
};

class test_function_xor : public bin_function
{
    public:
        template<class T1,class T2>
        static Real eval_scal_func(T1 a, T2 b)
        {
            using Compl_1   = typename md::complex_type<T1>::type;
            using Compl_2   = typename md::complex_type<T2>::type;

            Real dif	= norm_1(Matrix(matcl::elem_xor(a,b)) - matcl::elem_xor(full(a),full(b)));;
            dif			+=norm_1(Matrix(matcl::elem_xor(a,b)) - Matrix(matcl::elem_xor(Compl_1(a),(b))));;
            dif			+=norm_1(Matrix(matcl::elem_xor(a,b)) - Matrix(matcl::elem_xor((a),Compl_2(b))));;
            dif			+=norm_1(Matrix(matcl::elem_xor(a,b)) - Matrix(matcl::elem_xor(Compl_1(a),Compl_2(b))));;

            return dif;
        };

    public:
        virtual Real	eval_mat(const Matrix& mat1,const Matrix& mat2,int code);
        virtual Real	eval_scalar(const Scalar& s1, const Scalar& s2,int code);
};

class test_function_rem : public bin_function
{
    public:
        template<class T1,class T2>
        static Real eval_scal_func(T1 a, T2 b)
        {
            using Compl_1   = typename md::complex_type<T1>::type;
            using Compl_2   = typename md::complex_type<T2>::type;

            Real dif	= norm_1(matcl::rem(real(a),real(b)) 
                        - matcl::rem(full(real(a)),full(real(b))));;

            Matrix res  = abs(a) + abs(b);

            if (dif < constants::eps(res.get_value_code()) * norm_1(res) * 10.)
                dif     = 0.0;

            return dif;
        };

    public:
        virtual Real	eval_mat(const Matrix& mat1,const Matrix& mat2,int code);
        virtual Real	eval_scalar(const Scalar& s1, const Scalar& s2,int code);
};

class test_function_mod : public bin_function
{
    public:
        template<class T1,class T2>
        static Real eval_scal_func(T1 a, T2 b)
        {
            using Compl_1   = typename md::complex_type<T1>::type;
            using Compl_2   = typename md::complex_type<T2>::type;

            Real dif	= norm_1(matcl::mod(real(a),real(b)) 
                                 - matcl::mod(full(real(a)),full(real(b))));;

            Matrix res  = abs(a) + abs(b);

            if (dif < constants::eps(res.get_value_code()) * norm_1(res) * 10.)
                dif     = 0.0;

            return dif;
        };

    public:
        virtual Real	eval_mat(const Matrix& mat1,const Matrix& mat2,int code);
        virtual Real	eval_scalar(const Scalar& s1, const Scalar& s2,int code);
};

class test_function_kron : public bin_function
{
    public:
        template<class T1,class T2>
        static Real eval_scal_func(T1 a, T2 b)
        {
            using Compl_1   = typename md::complex_type<T1>::type;
            using Compl_2   = typename md::complex_type<T2>::type;

            Real dif	= norm_1(matcl::kron(a,b) - matcl::kron(full(a),full(b)));;
            dif			+=norm_1(matcl::kron(a,b) - matcl::kron(Compl_1(a),(b)));;
            dif			+=norm_1(matcl::kron(a,b) - matcl::kron((a),Compl_2(b)));;
            dif			+=norm_1(matcl::kron(a,b) - matcl::kron(Compl_1(a),Compl_2(b)));;

            Matrix res  = abs(a) + abs(b) + abs(kron(a,b));

            if (dif < constants::eps(res.get_value_code()) * norm_1(res) * 10.)
                dif     = 0.0;


            return dif;
        };

    public:
        virtual Real	eval_mat(const Matrix& mat1,const Matrix& mat2,int code);
        virtual Real	eval_scalar(const Scalar& s1, const Scalar& s2,int code);
};

class test_function_op_or : public bin_function
{
    public:
        template<class T1,class T2>
        static Real eval_scal_func(T1 a, T2 b)
        {
            using Compl_1   = typename md::complex_type<T1>::type;
            using Compl_2   = typename md::complex_type<T2>::type;

            Real dif = 0;
            dif			+=norm_1(Matrix(matcl::elem_or(a,b)) - matcl::elem_or(full(a),full(b)));;
            dif			+=norm_1(Matrix(a||b) - (full(a)|full(b)));;
            dif			+=norm_1(Matrix(matcl::elem_or(a,b)) - Matrix(matcl::elem_or(Compl_1(a),(b))));;
            dif			+=norm_1(Matrix(matcl::elem_or(a,b)) - Matrix(matcl::elem_or((a),Compl_2(b))));;
            dif			+=norm_1(Matrix(matcl::elem_or(a,b)) - Matrix(matcl::elem_or(Compl_1(a),Compl_2(b))));;

            return dif;
        };

    public:
        virtual Real	eval_mat(const Matrix& mat1,const Matrix& mat2,int code);
        virtual Real	eval_scalar(const Scalar& s1, const Scalar& s2,int code);
};

class test_function_op_and : public bin_function
{
    public:
        template<class T1,class T2>
        static Real eval_scal_func(T1 a, T2 b)
        {
            using Compl_1   = typename md::complex_type<T1>::type;
            using Compl_2   = typename md::complex_type<T2>::type;

            Real dif = 0;
            dif			+=norm_1(Matrix(matcl::elem_and(a,b)) - matcl::elem_and(full(a),full(b)));;
            dif			+=norm_1(Matrix(a&&b) - (full(a)&full(b)));;
            dif			+=norm_1(Matrix(matcl::elem_and(a,b)) - Matrix(matcl::elem_and(Compl_1(a),(b))));;
            dif			+=norm_1(Matrix(matcl::elem_and(a,b)) - Matrix(matcl::elem_and((a),Compl_2(b))));;
            dif			+=norm_1(Matrix(matcl::elem_and(a,b)) - Matrix(matcl::elem_and(Compl_1(a),Compl_2(b))));;

            return dif;
        };

    public:
        virtual Real	eval_mat(const Matrix& mat1,const Matrix& mat2,int code);
        virtual Real	eval_scalar(const Scalar& s1, const Scalar& s2,int code);
};

class test_function_op_eeq : public bin_function
{
    public:
        template<class T1,class T2>
        static Real eval_scal_func(T1 a, T2 b)
        {
            using Compl_1   = typename md::complex_type<T1>::type;
            using Compl_2   = typename md::complex_type<T2>::type;

            Real dif = 0;
            dif			+=norm_1(Matrix(matcl::eeq(a,b)) - matcl::eeq(full(a),full(b)));;
            dif			+=norm_1(Matrix(a==b) - (full(a)==full(b)));;
            dif			+=norm_1(Matrix(matcl::eeq(a,b)) - Matrix(matcl::eeq(Compl_1(a),(b))));;
            dif			+=norm_1(Matrix(matcl::eeq(a,b)) - Matrix(matcl::eeq((a),Compl_2(b))));;
            dif			+=norm_1(Matrix(matcl::eeq(a,b)) - Matrix(matcl::eeq(Compl_1(a),Compl_2(b))));;

            return dif;
        };

    public:
        virtual Real	eval_mat(const Matrix& mat1,const Matrix& mat2,int code);
        virtual Real	eval_scalar(const Scalar& s1, const Scalar& s2,int code);
};

class test_function_op_neq : public bin_function
{
    public:
        template<class T1,class T2>
        static Real eval_scal_func(T1 a, T2 b)
        {
            using Compl_1   = typename md::complex_type<T1>::type;
            using Compl_2   = typename md::complex_type<T2>::type;

            Real dif = 0;
            dif			+=norm_1(Matrix(matcl::neq(a,b)) - matcl::neq(full(a),full(b)));;
            dif			+=norm_1(Matrix(a!=b) - (full(a)!=full(b)));;
            dif			+=norm_1(Matrix(matcl::neq(a,b)) - Matrix(matcl::neq(Compl_1(a),(b))));;
            dif			+=norm_1(Matrix(matcl::neq(a,b)) - Matrix(matcl::neq((a),Compl_2(b))));;
            dif			+=norm_1(Matrix(matcl::neq(a,b)) - Matrix(matcl::neq(Compl_1(a),Compl_2(b))));;

            return dif;
        };

    public:
        virtual Real	eval_mat(const Matrix& mat1,const Matrix& mat2,int code);
        virtual Real	eval_scalar(const Scalar& s1, const Scalar& s2,int code);
};

class test_function_eeq_nan : public bin_function
{
    public:
        template<class T1,class T2>
        static Real eval_scal_func(T1 a, T2 b)
        {
            using Compl_1   = typename md::complex_type<T1>::type;
            using Compl_2   = typename md::complex_type<T2>::type;

            Real dif = 0;
            dif			+=norm_1(Matrix(matcl::eeq_nan(a,b)) - matcl::eeq_nan(full(a),full(b)));;
            dif			+=norm_1(Matrix(eeq_nan(a,b)) - (eeq_nan(full(a),full(b))));;
            dif			+=norm_1(Matrix(matcl::eeq_nan(a,b)) - Matrix(matcl::eeq_nan(Compl_1(a),(b))));;
            dif			+=norm_1(Matrix(matcl::eeq_nan(a,b)) - Matrix(matcl::eeq_nan((a),Compl_2(b))));;
            dif			+=norm_1(Matrix(matcl::eeq_nan(a,b)) - Matrix(matcl::eeq_nan(Compl_1(a),Compl_2(b))));;

            return dif;
        };

    public:
        virtual Real	eval_mat(const Matrix& mat1,const Matrix& mat2,int code);
        virtual Real	eval_scalar(const Scalar& s1, const Scalar& s2,int code);
};

class test_function_neq_nan : public bin_function
{
    public:
        template<class T1,class T2>
        static Real eval_scal_func(T1 a, T2 b)
        {
            using Compl_1   = typename md::complex_type<T1>::type;
            using Compl_2   = typename md::complex_type<T2>::type;

            Real dif = 0;
            dif			+=norm_1(Matrix(matcl::neq_nan(a,b)) - matcl::neq_nan(full(a),full(b)));;
            dif			+=norm_1(Matrix(neq_nan(a,b)) - (neq_nan(full(a),full(b))));;
            dif			+=norm_1(Matrix(matcl::neq_nan(a,b)) - Matrix(matcl::neq_nan(Compl_1(a),(b))));;
            dif			+=norm_1(Matrix(matcl::neq_nan(a,b)) - Matrix(matcl::neq_nan((a),Compl_2(b))));;
            dif			+=norm_1(Matrix(matcl::neq_nan(a,b)) - Matrix(matcl::neq_nan(Compl_1(a),Compl_2(b))));;

            return dif;
        };

    public:
        virtual Real	eval_mat(const Matrix& mat1,const Matrix& mat2,int code);
        virtual Real	eval_scalar(const Scalar& s1, const Scalar& s2,int code);
};

class test_function_op_lt : public bin_function
{
    public:
        template<class T1,class T2>
        static Real eval_scal_func(T1 a, T2 b)
        {
            using Compl_1   = typename md::complex_type<T1>::type;
            using Compl_2   = typename md::complex_type<T2>::type;

            Real dif = 0;
            dif			+=norm_1(Matrix(matcl::lt(a,b)) - matcl::lt(full(a),full(b)));;
            dif			+=norm_1(Matrix(a<b) - (full(a)<full(b)));;
            dif			+=norm_1(Matrix(matcl::lt(a,b)) - Matrix(matcl::lt(Compl_1(a),(b))));;
            dif			+=norm_1(Matrix(matcl::lt(a,b)) - Matrix(matcl::lt((a),Compl_2(b))));;
            dif			+=norm_1(Matrix(matcl::lt(a,b)) - Matrix(matcl::lt(Compl_1(a),Compl_2(b))));;

            return dif;
        };

    public:
        virtual Real	eval_mat(const Matrix& mat1,const Matrix& mat2,int code);
        virtual Real	eval_scalar(const Scalar& s1, const Scalar& s2,int code);
};

class test_function_op_leq : public bin_function
{
    public:
        template<class T1,class T2>
        static Real eval_scal_func(T1 a, T2 b)
        {
            using Compl_1   = typename md::complex_type<T1>::type;
            using Compl_2   = typename md::complex_type<T2>::type;

            Real dif = 0;
            dif			+=norm_1(Matrix(matcl::leq(a,b)) - matcl::leq(full(a),full(b)));;
            dif			+=norm_1(Matrix(a<=b) - (full(a)<=full(b)));;
            dif			+=norm_1(Matrix(matcl::leq(a,b)) - Matrix(matcl::leq(Compl_1(a),(b))));;
            dif			+=norm_1(Matrix(matcl::leq(a,b)) - Matrix(matcl::leq((a),Compl_2(b))));;
            dif			+=norm_1(Matrix(matcl::leq(a,b)) - Matrix(matcl::leq(Compl_1(a),Compl_2(b))));;

            return dif;
        };

    public:
        virtual Real	eval_mat(const Matrix& mat1,const Matrix& mat2,int code);
        virtual Real	eval_scalar(const Scalar& s1, const Scalar& s2,int code);
};

class test_function_op_gt : public bin_function
{
    public:
        template<class T1,class T2>
        static Real eval_scal_func(T1 a, T2 b)
        {
            using Compl_1   = typename md::complex_type<T1>::type;
            using Compl_2   = typename md::complex_type<T2>::type;

            Real dif = 0;
            dif			+=norm_1(Matrix(matcl::gt(a,b)) - matcl::gt(full(a),full(b)));;
            dif			+=norm_1(Matrix(a>b) - (full(a)>full(b)));;
            dif			+=norm_1(Matrix(matcl::gt(a,b)) - Matrix(matcl::gt(Compl_1(a),(b))));;
            dif			+=norm_1(Matrix(matcl::gt(a,b)) - Matrix(matcl::gt((a),Compl_2(b))));;
            dif			+=norm_1(Matrix(matcl::gt(a,b)) - Matrix(matcl::gt(Compl_1(a),Compl_2(b))));;
            return dif;
        };

    public:
        virtual Real	eval_mat(const Matrix& mat1,const Matrix& mat2,int code);
        virtual Real	eval_scalar(const Scalar& s1, const Scalar& s2,int code);
};

class test_function_op_geq : public bin_function
{
    public:
        template<class T1,class T2>
        static Real eval_scal_func(T1 a, T2 b)
        {
            using Compl_1   = typename md::complex_type<T1>::type;
            using Compl_2   = typename md::complex_type<T2>::type;

            Real dif = 0;
            dif			+=norm_1(Matrix(matcl::geq(a,b)) - matcl::geq(full(a),full(b)));;
            dif			+=norm_1(Matrix(a>=b) - (full(a)>=full(b)));;
            dif			+=norm_1(Matrix(matcl::geq(a,b)) - Matrix(matcl::geq(Compl_1(a),(b))));;
            dif			+=norm_1(Matrix(matcl::geq(a,b)) - Matrix(matcl::geq((a),Compl_2(b))));;
            dif			+=norm_1(Matrix(matcl::geq(a,b)) - Matrix(matcl::geq(Compl_1(a),Compl_2(b))));;
            return dif;
        };

    public:
        virtual Real	eval_mat(const Matrix& mat1,const Matrix& mat2,int code);
        virtual Real	eval_scalar(const Scalar& s1, const Scalar& s2,int code);
};

class test_function_copysign : public bin_function
{
    public:
        template<class T1,class T2>
        static Real eval_scal_func(T1 a, T2 b)
        {
            using Compl_1   = typename md::complex_type<T1>::type;
            using Compl_2   = typename md::complex_type<T2>::type;

            Real dif = 0.;
            dif			+=norm_1(real(a) - matcl::copysign(abs(real(a)), real(a)));;
            dif			+=norm_1(real(b) - matcl::copysign(abs(real(b)), real(b)));;

            Real sign   = signbit(real(b)) ? -1.0 : 1.0;
            dif			+=norm_1(matcl::copysign(real(a),real(b)) - abs(real(a)) * sign);;

            Matrix res  = 0.0;

            if (dif < constants::eps(res.get_value_code()) * norm_1(res) * 10.)
                dif     = 0.0;

            return dif;
        };

    public:
        virtual Real	eval_mat(const Matrix& mat1,const Matrix& mat2,int code);
        virtual Real	eval_scalar(const Scalar& s1, const Scalar& s2,int code);
};

class test_function_fdim : public bin_function
{
    public:
        template<class T1,class T2>
        static Real eval_scal_func(T1 a, T2 b)
        {
            using Compl_1   = typename md::complex_type<T1>::type;
            using Compl_2   = typename md::complex_type<T2>::type;

            auto res1   = matcl::fdim(real(a),real(b));
            auto res2   = matcl::max(Matrix(real(a)) - real(b), 0.0f);
            Real dif    = norm_1(res1 - res2);

            if (dif < matcl::eps(res1) * 10.)
                dif     = 0.0;

            return dif;
        };

    public:
        virtual Real	eval_mat(const Matrix& mat1,const Matrix& mat2,int code);
        virtual Real	eval_scalar(const Scalar& s1, const Scalar& s2,int code);
};

class test_function_nextafter : public bin_function
{
    public:
        template<class T1,class T2>
        static Real eval_scal_func(T1 a, T2 b)
        {
            using Compl_1   = typename md::complex_type<T1>::type;
            using Compl_2   = typename md::complex_type<T2>::type;

            using T1f       = typename md::unify_types<T1, Float>::type;
            using T2f       = typename md::unify_types<T2, Float>::type;

            using T1r       = typename md::real_type<T1f>::type;

            T1f af          = T1f(a);
            T2f bf          = T2f(b);

            (void)bf;

            T1r da          = abs(real(af)) + 1;

            Real dif = 0.;
            dif			+=norm_1(matcl::nextafter(real(af),real(af) + da) - matcl::nextabove(real(af)));;
            dif			+=norm_1(matcl::nextafter(real(af),real(af) - da) - matcl::nextbelow(real(af)));;

            return dif;
        };

    public:
        virtual Real	eval_mat(const Matrix& mat1,const Matrix& mat2,int code);
        virtual Real	eval_scalar(const Scalar& s1, const Scalar& s2,int code);
};

class test_function_fma : public bin_function
{
    public:
        template<class T1,class T2>
        static Real eval_scal_func(T1 a, T2 b)
        {
            using Real_1    = typename md::real_type<T1>::type;
            using Real_2    = typename md::real_type<T2>::type;
            using Real_t    = typename md::unify_types2<Real_1, Real_2, Float>::type;

            Real_t af       = Real_t(real(a));
            Real_t bf       = Real_t(real(b));
            Real_t cf       = bf;            

            //TODO fma_a
            Real_t val1     = matcl::fma_f(af, bf, cf);
            Real_t val2     = (af * bf) + cf;

            Real_t val_a    = (abs(af) * abs(bf)) + abs(cf);
            Real dif        = norm_1(val1 - val2);;

            if (dif < eps(val_a) * 100.)
                dif     = 0.0;

            return dif;
        };

    public:
        virtual Real	eval_mat(const Matrix& mat1,const Matrix& mat2,int code);
        virtual Real	eval_scalar(const Scalar& s1, const Scalar& s2,int code);
};

class test_function_fms : public bin_function
{
    public:
        template<class T1,class T2>
        static Real eval_scal_func(T1 a, T2 b)
        {
            using Real_1    = typename md::real_type<T1>::type;
            using Real_2    = typename md::real_type<T2>::type;
            using Real_t    = typename md::unify_types2<Real_1, Real_2, Float>::type;

            Real_t af       = Real_t(real(a));
            Real_t bf       = Real_t(real(b));
            Real_t cf       = bf;            

            //TODO: fms_a
            Real_t val1     = matcl::fms_f(af, bf, cf);
            Real_t val2     = (af * bf) - cf;
            Real_t val_a    = (abs(af) * abs(bf)) + abs(cf);

            Real dif        = norm_1(val1 - val2);;

            if (dif < eps(val_a) * 100.)
                dif     = 0.0;

            return dif;
        };

    public:
        virtual Real	eval_mat(const Matrix& mat1,const Matrix& mat2,int code);
        virtual Real	eval_scalar(const Scalar& s1, const Scalar& s2,int code);
};

class test_function_dot2_ac : public bin_function
{
    public:
        template<class T1,class T2>
        static Real eval_scal_func(T1 a, T2 b)
        {
            using Real_1    = typename md::real_type<T1>::type;
            using Real_2    = typename md::real_type<T2>::type;
            using Real_t    = typename md::unify_types2<Real_1, Real_2, Float>::type;

            Real_t af       = Real_t(real(a));
            Real_t bf       = Real_t(real(b));
            Real_t cf       = bf;            
            Real_t df       = af;            

            Real_t val1     = matcl::dot2_a(af, bf, cf, df);
            Real_t val2     = (af * bf) + (cf * df);

            Real dif        = norm_1(val1 - val2);;

            if (dif < eps(val1) * 100.)
                dif     = 0.0;

            return dif;
        };

    public:
        virtual Real	eval_mat(const Matrix& mat1,const Matrix& mat2,int code);
        virtual Real	eval_scalar(const Scalar& s1, const Scalar& s2,int code);
};

class test_function_powm1 : public bin_function
{
    public:
        template<class T1,class T2>
        static Real eval_scal_func(T1 a, T2 b)
        {
            using Compl_1   = typename md::complex_type<T1>::type;
            using Compl_2   = typename md::complex_type<T2>::type;

            Real dif = 0.;
            auto base   = matcl::powm1(real(a),real(b));
            auto base1  = matcl::pow(real(a),real(b)) - 1.0f;
            dif			+=norm_1(base - base1);;

            Real nr     = norm_1(base)  + abs(real(a)) + abs(real(b));

            if (dif < constants::eps<decltype(base1)>() * nr * 10000.0)
                dif     = 0.0;

            return dif;
        };

    public:
        virtual Real	eval_mat(const Matrix& mat1,const Matrix& mat2,int code);
        virtual Real	eval_scalar(const Scalar& s1, const Scalar& s2,int code);
};

};};
