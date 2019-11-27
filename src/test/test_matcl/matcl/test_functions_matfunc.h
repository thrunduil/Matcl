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

#include "matcl-matrep/matcl_matrep.h"
#include "test/test_matcl/framework/matrix_set/matrix_set.h"
#include "test/test_matcl/framework/matrix_utils.h"
#include "matcl-linalg/matcl_linalg.h"

namespace matcl { namespace test
{

class matfunc_functions_list
{
    private:
        const matrix_set&	m_tests;        
        options				m_options;

    public:
        rand_matrix_ptr     m_rand;
        dynamic_mat_set&    m_ms;

    public:
        void		make(options opts);

        Matrix		get_matrix(int code) const;

    private:
        void    test_symprod();
        void    test_herprod();
        void    test_symsum();
        void    test_hersum();
        void    test_scale_rows();
        void    test_scale_cols();
        void    test_scale_rowscols();

    public:
        matfunc_functions_list(const matrix_set&, rand_matrix_ptr rand, dynamic_mat_set& msd);

    private:
        matfunc_functions_list(const matfunc_functions_list&) = delete;
        matfunc_functions_list& operator=(const matfunc_functions_list&) = delete;
};

class test_function_symprod : public unary_function
{
    public:
        test_function_symprod(){};

        virtual Real eval_mat(const Matrix& mat,bool show_res,int code) override;
        virtual Real eval_scalar(const Scalar& s, bool show_res,int code) override;

        template<class Val>
        Real eval_scal_func(const Val& value) const
        {
            using Val_c = typename md::complex_type<Val>::type;
            Val_c val_c = Val_c(value);

            auto res1   = value * trans(value);
            auto res2   = symprod(value, false);

            auto res3   = trans(value) * value;
            auto res4   = symprod(value, true);

            auto res5   = symprod(val_c, false);
            auto res6   = symprod(val_c, true);

            Real dif	= norm_1(res1 - res2);
            dif	        += norm_1(res3 - res4);
            dif	        += norm_1(res2 - res5);
            dif	        += norm_1(res4 - res6);

            return dif;
        }
};

class test_function_herprod : public unary_function
{
    public:
        test_function_herprod(){};

        virtual Real eval_mat(const Matrix& mat,bool show_res,int code) override;
        virtual Real eval_scalar(const Scalar& s, bool show_res,int code) override;

        template<class Val>
        Real eval_scal_func(const Val& value) const
        {
            using Val_c = typename md::complex_type<Val>::type;
            Val_c val_c = Val_c(value);

            auto res1   = value * ctrans(value);
            auto res2   = herprod(value, false);

            auto res3   = ctrans(value) * value;
            auto res4   = herprod(value, true);

            auto res5   = herprod(val_c, false);
            auto res6   = herprod(val_c, true);

            Real dif	= norm_1(res1 - res2);
            dif	        += norm_1(res3 - res4);
            dif	        += norm_1(res2 - res5);
            dif	        += norm_1(res4 - res6);

            if (abs(dif) < eps(res2) * 10.0)
                dif     = 0.0;

            return dif;
        }
};

class test_function_symsum : public unary_function
{
    public:
        test_function_symsum(){};

        virtual Real eval_mat(const Matrix& mat,bool show_res,int code) override;
        virtual Real eval_scalar(const Scalar& s, bool show_res,int code) override;

        template<class Val>
        Real eval_scal_func(const Val& value) const
        {
            using Val_c = typename md::complex_type<Val>::type;
            Val_c val_c = Val_c(value);

            auto res1   = value + trans(value);
            auto res2   = symsum(value);

            auto res5   = symsum(val_c);

            Real dif	= norm_1(res1 - res2);
            dif	        += norm_1(res2 - res5);

            return dif;
        }
};

class test_function_hersum : public unary_function
{
    public:
        test_function_hersum(){};

        virtual Real eval_mat(const Matrix& mat,bool show_res,int code) override;
        virtual Real eval_scalar(const Scalar& s, bool show_res,int code) override;

        template<class Val>
        Real eval_scal_func(const Val& value) const
        {
            using Val_c = typename md::complex_type<Val>::type;
            Val_c val_c = Val_c(value);

            auto res1   = value + ctrans(value);
            auto res2   = hersum(value);

            auto res5   = hersum(val_c);

            Real dif	= norm_1(res1 - res2);
            dif	        += norm_1(res2 - res5);

            return dif;
        }
};

class test_function_scale_rows : public unary_function
{
    private:
        dynamic_mat_set&    m_ms;
        long                m_new_objects;

    public:
        test_function_scale_rows(dynamic_mat_set& ms)  :m_ms(ms), m_new_objects(0){};

        virtual Real eval_mat(const Matrix& mat,bool show_res,int code) override;
        virtual Real eval_scalar(const Scalar& s, bool show_res,int code) override;

        virtual long n_new_objects() override   { return m_new_objects; };

    private:
        bool    check_value_code(value_code v1, value_code v2) const;
};

class test_function_scale_cols : public unary_function
{
    private:
        dynamic_mat_set&    m_ms;
        long                m_new_objects;

    public:
        test_function_scale_cols(dynamic_mat_set& ms)  :m_ms(ms), m_new_objects(0){};

        virtual Real eval_mat(const Matrix& mat,bool show_res,int code) override;
        virtual Real eval_scalar(const Scalar& s, bool show_res,int code) override;

        virtual long n_new_objects() override   { return m_new_objects; };

    private:
        bool    check_value_code(value_code v1, value_code v2) const;
};

class test_function_scale_rowscols : public unary_function
{
    private:
        dynamic_mat_set&    m_ms;
        long                m_new_objects;

    public:
        test_function_scale_rowscols(dynamic_mat_set& ms)  :m_ms(ms), m_new_objects(0){};

        virtual Real eval_mat(const Matrix& mat,bool show_res,int code) override;
        virtual Real eval_scalar(const Scalar& s, bool show_res,int code) override;

        virtual long n_new_objects() override   { return m_new_objects; };

    private:
        bool    check_value_code(value_code v1, value_code v2) const;
};
};};
