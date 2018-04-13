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
#include "scal_func_helper.h"

namespace matcl { namespace test
{

class vecfunc_functions_list
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
        void        test_nnz();         void    test_cumsum();
        void        test_prod();        void    test_cumprod();			
        void        test_min();         void    test_max();
        void        test_mean();        void    test_min_abs();
        void        test_max_abs();     void    test_std();        
        void        test_min2();	    void    test_max2();
        void        test_min_abs2();    void    test_max_abs2();
        void        test_all();		    void    test_any();
        void        test_sum();

        void        test_nnz_vec();     void    test_all_vec();
        void        test_any_vec();     void    test_sum_vec();
        void        test_prod_vec();    void    test_mean_vec();
        void        test_std_vec();     void    test_min_vec();
        void        test_max_vec();     void    test_min_abs_vec();
        void        test_max_abs_vec(); void    test_min2_vec();
        void        test_max2_vec();    void    test_min_abs2_vec();
        void        test_max_abs2_vec();

    public:
        vecfunc_functions_list(const matrix_set&,rand_matrix_ptr rand);

    private:
        vecfunc_functions_list(const vecfunc_functions_list&) = delete;
        vecfunc_functions_list& operator=(const vecfunc_functions_list&) = delete;
};

class test_function_nnz : public unary_function
{
    private:
        Integer d;

    public:
        test_function_nnz(Integer d) : d(d){};

        virtual Real eval_mat(const Matrix& mat,bool show_res,int code) override;
        virtual Real eval_scalar(const Scalar& s, bool show_res,int code) override;

        template<class Val>
        Real eval_scal_func(const Val& value) const
        {
            using Val_c = typename md::complex_type<Val>::type;
            Val_c val_c = Val_c(value);

            Real dif		= norm_1(nnz(value,d) - nnz(full(value),d));
            dif				+=norm_1(nnz(value,d) - nnz(val_c,d));
            return dif;
        }
};
class test_function_nnz_vec : public unary_function
{
    public:
        test_function_nnz_vec() {};

        virtual Real eval_mat(const Matrix& mat,bool show_res,int code) override;
        virtual Real eval_scalar(const Scalar& s, bool show_res,int code) override;

        template<class Val>
        Real eval_scal_func(const Val& value) const
        {
            using Val_c = typename md::complex_type<Val>::type;
            Val_c val_c = Val_c(value);

            Real dif		= norm_1(nnz_vec(value) - nnz_vec(full(value)));
            dif				+=norm_1(nnz_vec(value) - nnz_vec(val_c));
            return dif;
        }
};

class test_function_sum : public unary_function
{
    private:
        Integer d;

    public:
        test_function_sum(Integer d) : d(d){};

        template<class Val>
        Real eval_scal_func(const Val& value) const
        {
            using Val_c = typename md::complex_type<Val>::type;
            Val_c val_c = Val_c(value);

            Real dif		= norm_1(sum(value,d) - sum(full(value),d));
            dif				+=norm_1(sum(value,d) - sum(val_c,d));
            return dif;
        }

        virtual Real eval_mat(const Matrix& mat,bool show_res,int code) override;
        virtual Real eval_scalar(const Scalar& s, bool show_res,int code) override;
};

class test_function_sum_vec : public unary_function
{
    public:
        test_function_sum_vec() {};

        template<class Val>
        Real eval_scal_func(const Val& value) const
        {
            using Val_c = typename md::complex_type<Val>::type;
            Val_c val_c = Val_c(value);

            Real dif		= norm_1(sum_vec(value) - sum_vec(full(value)));
            dif				+=norm_1(sum_vec(value) - sum_vec(val_c));
            return dif;
        }

        virtual Real eval_mat(const Matrix& mat,bool show_res,int code) override;
        virtual Real eval_scalar(const Scalar& s, bool show_res,int code) override;
};

class test_function_cumsum : public unary_function
{
    private:
        Integer d;

    public:
        test_function_cumsum(Integer d) : d(d){};

        virtual Real eval_mat(const Matrix& mat,bool show_res,int code) override;
        virtual Real eval_scalar(const Scalar& s, bool show_res,int code) override;

        template<class Val>
        Real eval_scal_func(const Val& value) const
        {
            using Val_c = typename md::complex_type<Val>::type;
            Val_c val_c = Val_c(value);

            Real dif		= norm_1(cumsum(value,d) - cumsum(full(value),d));
            dif				+=norm_1(cumsum(value,d) - cumsum(val_c,d));
            return dif;
        }
};

class test_function_prod : public unary_function
{
    private:
        Integer d;

    public:
        test_function_prod(Integer d) : d(d){};

        virtual Real eval_mat(const Matrix& mat,bool show_res,int code) override;
        virtual Real eval_scalar(const Scalar& s, bool show_res,int code) override;

        template<class Val>
        Real eval_scal_func(const Val& value) const
        {
            using Val_c = typename md::complex_type<Val>::type;
            Val_c val_c = Val_c(value);

            Real dif		= norm_1(prod(value,d) - prod(full(value),d));
            dif				+=norm_1(prod(value,d) - prod(val_c,d));
            return dif;
        }
};

class test_function_prod_vec : public unary_function
{
    public:
        test_function_prod_vec() {};

        virtual Real eval_mat(const Matrix& mat,bool show_res,int code) override;
        virtual Real eval_scalar(const Scalar& s, bool show_res,int code) override;

        template<class Val>
        Real eval_scal_func(const Val& value) const
        {
            using Val_c = typename md::complex_type<Val>::type;
            Val_c val_c = Val_c(value);

            Real dif		= norm_1(prod_vec(value) - prod_vec(full(value)));
            dif				+=norm_1(prod_vec(value) - prod_vec(val_c));
            return dif;
        }
};

class test_function_cumprod : public unary_function
{
    private:
        Integer d;

    public:
        test_function_cumprod(Integer d) : d(d){};

        virtual Real eval_mat(const Matrix& mat,bool show_res,int code) override;
        virtual Real eval_scalar(const Scalar& s, bool show_res,int code) override;

        template<class Val>
        Real eval_scal_func(const Val& value) const
        {
            using Val_c = typename md::complex_type<Val>::type;
            Val_c val_c = Val_c(value);

            Real dif		= norm_1(cumprod(value,d) - cumprod(full(value),d));
            dif				+=norm_1(cumprod(value,d) - cumprod(val_c,d));
            return dif;
        }
};

class test_function_min : public unary_function
{
    private:
        Integer d;

    public:
        test_function_min(Integer d) : d(d){};

        virtual Real eval_mat(const Matrix& mat,bool show_res,int code) override;
        virtual Real eval_scalar(const Scalar& s, bool show_res,int code) override;

        template<class Val>
        Real eval_scal_func(const Val& value) const
        {
            using Val_c = typename md::complex_type<Val>::type;
            Val_c val_c = Val_c(value);

            Real dif		= norm_1(min_d(value,d) - min_d(full(value),d));
            dif				+=norm_1(min_d(value,d) - min_d(val_c,d));
            return dif;
        }
};

class test_function_min_vec : public unary_function
{
    private:
        Integer d;

    public:
        test_function_min_vec(Integer d) : d(d){};

        virtual Real eval_mat(const Matrix& mat,bool show_res,int code) override;
        virtual Real eval_scalar(const Scalar& s, bool show_res,int code) override;

        template<class Val>
        Real eval_scal_func(const Val& value) const
        {
            using Val_c = typename md::complex_type<Val>::type;
            Val_c val_c = Val_c(value);

            Real dif		= norm_1(min_vec(value) - min_vec(full(value)));
            dif				+=norm_1(min_vec(value) - min_vec(val_c));
            return dif;
        }
};

class test_function_min_abs : public unary_function
{
    private:
        Integer d;

    public:
        test_function_min_abs(Integer d) : d(d){};

        virtual Real eval_mat(const Matrix& mat,bool show_res,int code) override;
        virtual Real eval_scalar(const Scalar& s, bool show_res,int code) override;

        template<class Val>
        Real eval_scal_func(const Val& value) const
        {
            using Val_c = typename md::complex_type<Val>::type;
            Val_c val_c = Val_c(value);

            Real dif		= norm_1(min_abs_d(value,d) - min_abs_d(full(value),d));
            dif				+=norm_1(min_abs_d(value,d) - min_abs_d(val_c,d));
            return dif;
        }
};

class test_function_min_abs_vec : public unary_function
{
    public:
        test_function_min_abs_vec() {};

        virtual Real eval_mat(const Matrix& mat,bool show_res,int code) override;
        virtual Real eval_scalar(const Scalar& s, bool show_res,int code) override;

        template<class Val>
        Real eval_scal_func(const Val& value) const
        {
            using Val_c = typename md::complex_type<Val>::type;
            Val_c val_c = Val_c(value);

            Real dif		= norm_1(min_abs_vec(value) - min_abs_vec(full(value)));
            dif				+=norm_1(min_abs_vec(value) - min_abs_vec(val_c));
            return dif;
        }
};

class test_function_min2 : public unary_function
{
    private:
        Integer d;

    public:
        test_function_min2(Integer d) : d(d){};

        virtual Real eval_mat(const Matrix& mat,bool show_res,int code) override;
        virtual Real eval_scalar(const Scalar& s, bool show_res,int code) override;

        template<class Val>
        Real eval_scal_func(const Val& value) const
        {
            using Val_c = typename md::complex_type<Val>::type;
            Val_c val_c = Val_c(value);

            Real dif		= norm_1(min2(value,d).get<1>() - min2(full(value),d).get<1>());
            dif				+=norm_1(min2(value,d).get<1>() - min2(val_c,d).get<1>());
            return dif;
        }
};

class test_function_min2_vec : public unary_function
{
    public:
        test_function_min2_vec() {};

        virtual Real eval_mat(const Matrix& mat,bool show_res,int code) override;
        virtual Real eval_scalar(const Scalar& s, bool show_res,int code) override;

        template<class Val>
        Real eval_scal_func(const Val& value) const
        {
            using Val_c = typename md::complex_type<Val>::type;
            Val_c val_c = Val_c(value);

            Real dif		= norm_1(min2_vec(value).get<1>() - min2_vec(full(value)).get<1>());
            dif				+=norm_1(min2_vec(value).get<1>() - min2_vec(val_c).get<1>());
            return dif;
        }
};

class test_function_min_abs2 : public unary_function
{
    private:
        Integer d;

    public:
        test_function_min_abs2(Integer d) : d(d){};

        virtual Real eval_mat(const Matrix& mat,bool show_res,int code) override;
        virtual Real eval_scalar(const Scalar& s, bool show_res,int code) override;

        template<class Val>
        Real eval_scal_func(const Val& value) const
        {
            using Val_c = typename md::complex_type<Val>::type;
            Val_c val_c = Val_c(value);

            Real dif		= norm_1(min_abs2(value,d).get<1>() - min_abs2(full(value),d).get<1>());
            dif				+=norm_1(min_abs2(value,d).get<1>() - min_abs2(val_c,d).get<1>());
            return dif;
        }
};

class test_function_min_abs2_vec : public unary_function
{
    public:
        test_function_min_abs2_vec() {};

        virtual Real eval_mat(const Matrix& mat,bool show_res,int code) override;
        virtual Real eval_scalar(const Scalar& s, bool show_res,int code) override;

        template<class Val>
        Real eval_scal_func(const Val& value) const
        {
            using Val_c = typename md::complex_type<Val>::type;
            Val_c val_c = Val_c(value);

            Real dif		= norm_1(min_abs2_vec(value).get<1>() - min_abs2_vec(full(value)).get<1>());
            dif				+=norm_1(min_abs2_vec(value).get<1>() - min_abs2_vec(val_c).get<1>());
            return dif;
        }
};

class test_function_max : public unary_function
{
    private:
        Integer d;

    public:
        test_function_max(Integer d) : d(d){};

        virtual Real eval_mat(const Matrix& mat,bool show_res,int code) override;
        virtual Real eval_scalar(const Scalar& s, bool show_res,int code) override;

        template<class Val>
        Real eval_scal_func(const Val& value) const
        {
            using Val_c = typename md::complex_type<Val>::type;
            Val_c val_c = Val_c(value);

            Real dif		= norm_1(max_d(value,d) - max_d(full(value),d));
            dif				+=norm_1(max_d(value,d) - max_d(val_c,d));
            return dif;
        }
};
class test_function_max_vec : public unary_function
{
    public:
        test_function_max_vec() {};

        virtual Real eval_mat(const Matrix& mat,bool show_res,int code) override;
        virtual Real eval_scalar(const Scalar& s, bool show_res,int code) override;

        template<class Val>
        Real eval_scal_func(const Val& value) const
        {
            using Val_c = typename md::complex_type<Val>::type;
            Val_c val_c = Val_c(value);

            Real dif		= norm_1(max_vec(value) - max_vec(full(value)));
            dif				+=norm_1(max_vec(value) - max_vec(val_c));
            return dif;
        }
};

class test_function_max_abs : public unary_function
{
    private:
        Integer d;

    public:
        test_function_max_abs(Integer d) : d(d){};

        virtual Real eval_mat(const Matrix& mat,bool show_res,int code) override;
        virtual Real eval_scalar(const Scalar& s, bool show_res,int code) override;

        template<class Val>
        Real eval_scal_func(const Val& value) const
        {
            using Val_c = typename md::complex_type<Val>::type;
            Val_c val_c = Val_c(value);

            Real dif		= norm_1(max_abs_d(value,d) - max_abs_d(full(value),d));
            dif				+=norm_1(max_abs_d(value,d) - max_abs_d(val_c,d));
            return dif;
        }
};

class test_function_max_abs_vec : public unary_function
{
    public:
        test_function_max_abs_vec() {};

        virtual Real eval_mat(const Matrix& mat,bool show_res,int code) override;
        virtual Real eval_scalar(const Scalar& s, bool show_res,int code) override;

        template<class Val>
        Real eval_scal_func(const Val& value) const
        {
            using Val_c = typename md::complex_type<Val>::type;
            Val_c val_c = Val_c(value);

            Real dif		= norm_1(max_abs_vec(value) - max_abs_vec(full(value)));
            dif				+=norm_1(max_abs_vec(value) - max_abs_vec(val_c));
            return dif;
        }
};

class test_function_max2 : public unary_function
{
    private:
        Integer d;

    public:
        test_function_max2(Integer d) : d(d){};

        virtual Real eval_mat(const Matrix& mat,bool show_res,int code) override;
        virtual Real eval_scalar(const Scalar& s, bool show_res,int code) override;

        template<class Val>
        Real eval_scal_func(const Val& value) const
        {
            using Val_c = typename md::complex_type<Val>::type;
            Val_c val_c = Val_c(value);

            Real dif		= norm_1(max2(value,d).get<1>() - max2(full(value),d).get<1>());
            dif				+=norm_1(max2(value,d).get<1>() - max2(val_c,d).get<1>());
            return dif;
        }
};

class test_function_max2_vec : public unary_function
{
    public:
        test_function_max2_vec() {};

        virtual Real eval_mat(const Matrix& mat,bool show_res,int code) override;
        virtual Real eval_scalar(const Scalar& s, bool show_res,int code) override;

        template<class Val>
        Real eval_scal_func(const Val& value) const
        {
            using Val_c = typename md::complex_type<Val>::type;
            Val_c val_c = Val_c(value);

            Real dif		= norm_1(max2_vec(value).get<1>() - max2_vec(full(value)).get<1>());
            dif				+=norm_1(max2_vec(value).get<1>() - max2_vec(val_c).get<1>());
            return dif;
        }
};

class test_function_max_abs2 : public unary_function
{
    private:
        Integer d;

    public:
        test_function_max_abs2(Integer d) : d(d){};

        virtual Real eval_mat(const Matrix& mat,bool show_res,int code) override;
        virtual Real eval_scalar(const Scalar& s, bool show_res,int code) override;

        template<class Val>
        Real eval_scal_func(const Val& value) const
        {
            using Val_c = typename md::complex_type<Val>::type;
            Val_c val_c = Val_c(value);

            Real dif		= norm_1(max_abs2(value,d).get<1>() - max_abs2(full(value),d).get<1>());
            dif				+=norm_1(max_abs2(value,d).get<1>() - max_abs2(val_c,d).get<1>());
            return dif;
        }
};

class test_function_max_abs2_vec : public unary_function
{
    public:
        test_function_max_abs2_vec() {};

        virtual Real eval_mat(const Matrix& mat,bool show_res,int code) override;
        virtual Real eval_scalar(const Scalar& s, bool show_res,int code) override;

        template<class Val>
        Real eval_scal_func(const Val& value) const
        {
            using Val_c = typename md::complex_type<Val>::type;
            Val_c val_c = Val_c(value);

            Real dif		= norm_1(max_abs2_vec(value).get<1>() - max_abs2_vec(full(value)).get<1>());
            dif				+=norm_1(max_abs2_vec(value).get<1>() - max_abs2_vec(val_c).get<1>());
            return dif;
        }
};

class test_function_mean : public unary_function
{
    private:
        Integer d;

    public:
        test_function_mean(Integer d) : d(d){};

        virtual Real eval_mat(const Matrix& mat,bool show_res,int code) override;
        virtual Real eval_scalar(const Scalar& s, bool show_res,int code) override;

        template<class Val>
        Real eval_scal_func(const Val& value) const
        {
            using Val_c = typename md::complex_type<Val>::type;
            Val_c val_c = Val_c(value);

            Real dif		= norm_1(mean(value,d) - mean(full(value),d));
            dif				+=norm_1(mean(value,d) - mean(val_c,d));
            return dif;
        }
};

class test_function_mean_vec : public unary_function
{
    public:
        test_function_mean_vec() {};

        virtual Real eval_mat(const Matrix& mat,bool show_res,int code) override;
        virtual Real eval_scalar(const Scalar& s, bool show_res,int code) override;

        template<class Val>
        Real eval_scal_func(const Val& value) const
        {
            using Val_c = typename md::complex_type<Val>::type;
            Val_c val_c = Val_c(value);

            Real dif		= norm_1(mean_vec(value) - mean_vec(full(value)));
            dif				+=norm_1(mean_vec(value) - mean_vec(val_c));
            return dif;
        }
};

class test_function_std: public unary_function
{
    private:
        Integer d;

    public:
        test_function_std(Integer d) : d(d){};

        virtual Real eval_mat(const Matrix& mat,bool show_res,int code) override;
        virtual Real eval_scalar(const Scalar& s, bool show_res,int code) override;

        template<class Val>
        Real eval_scal_func(const Val& value) const
        {
            using Val_c = typename md::complex_type<Val>::type;
            Val_c val_c = Val_c(value);

            Real dif		= norm_1(std(value,d,true) - std(full(value),d,true));
            dif		        = norm_1(std(value,d,false) - std(full(value),d,false));
            dif				+=norm_1(std(value,d,true) - std(val_c,d,true));
            dif				+=norm_1(std(value,d,false) - std(val_c,d,false));
            return dif;
        }
};

class test_function_std_vec: public unary_function
{
    private:
        Integer d;

    public:
        test_function_std_vec(Integer d) : d(d){};

        virtual Real eval_mat(const Matrix& mat,bool show_res,int code) override;
        virtual Real eval_scalar(const Scalar& s, bool show_res,int code) override;

        template<class Val>
        Real eval_scal_func(const Val& value) const
        {
            using Val_c = typename md::complex_type<Val>::type;
            Val_c val_c = Val_c(value);

            Real dif		= norm_1(std_vec(value,true) - std_vec(full(value),true));
            dif		        = norm_1(std_vec(value,false) - std_vec(full(value),false));
            dif				+=norm_1(std_vec(value,true) - std_vec(val_c,true));
            dif				+=norm_1(std_vec(value,false) - std_vec(val_c,false));
            return dif;
        }
};

class test_function_all : public unary_function
{
    private:
        Integer d;

    public:
        test_function_all(Integer d) : d(d){};

        virtual Real eval_mat(const Matrix& mat,bool show_res,int code) override;
        virtual Real eval_scalar(const Scalar& s, bool show_res,int code) override;

        template<class Val>
        Real eval_scal_func(const Val& value) const
        {
            using Val_c = typename md::complex_type<Val>::type;
            Val_c val_c = Val_c(value);

            Real dif		= norm_1(Matrix(all(value,d)) - all(full(value),d));
            dif				+=norm_1(all(value,d) - all(val_c,d));
            return dif;
        }
};

class test_function_all_vec : public unary_function
{
    public:
        test_function_all_vec() {};

        virtual Real eval_mat(const Matrix& mat,bool show_res,int code) override;
        virtual Real eval_scalar(const Scalar& s, bool show_res,int code) override;

        template<class Val>
        Real eval_scal_func(const Val& value) const
        {
            using Val_c = typename md::complex_type<Val>::type;
            Val_c val_c = Val_c(value);

            Real dif		= norm_1(Matrix(all_vec(value)) - Matrix(all_vec(full(value))));
            dif				+=norm_1(all_vec(value) - all_vec(val_c));
            return dif;
        }
};

class test_function_any : public unary_function
{
    private:
        Integer d;

    public:
        test_function_any(Integer d) : d(d){};

        virtual Real eval_mat(const Matrix& mat,bool show_res,int code) override;
        virtual Real eval_scalar(const Scalar& s, bool show_res,int code) override;

        template<class Val>
        Real eval_scal_func(const Val& value) const
        {
            using Val_c = typename md::complex_type<Val>::type;
            Val_c val_c = Val_c(value);

            Real dif		= norm_1(Matrix(any(value,d)) - any(full(value),d));
            dif				+=norm_1(Matrix(any(value,d)) - Matrix(any(val_c,d)));
            return dif;
        }
};

class test_function_any_vec : public unary_function
{
    private:
        Integer d;

    public:
        test_function_any_vec(Integer d) : d(d){};

        virtual Real eval_mat(const Matrix& mat,bool show_res,int code) override;
        virtual Real eval_scalar(const Scalar& s, bool show_res,int code) override;

        template<class Val>
        Real eval_scal_func(const Val& value) const
        {
            using Val_c = typename md::complex_type<Val>::type;
            Val_c val_c = Val_c(value);

            Real dif		= norm_1(Matrix(any_vec(value)) - Matrix(any_vec(full(value))));
            dif				+=norm_1(Matrix(any_vec(value)) - Matrix(any_vec(val_c)));
            return dif;
        }
};

};};
