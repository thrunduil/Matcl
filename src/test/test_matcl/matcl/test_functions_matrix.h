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

class matrix_functions_list
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
        void    test_get_scalar();      void    test_get_scalar_unique();
        void    test_get_array();		void    test_get_const_array();
        void    test_get_impl_unique();
        void    test_mat_compile();

        void    test_resize_reserve();  void    test_resize_reserve_b();
        
    public:
        matrix_functions_list(const matrix_set&,rand_matrix_ptr rand);

    private:
        matrix_functions_list(const matrix_functions_list&) = delete;
        matrix_functions_list& operator=(const matrix_functions_list&) = delete;
};

class test_function_mat_compile : public unary_function
{
    public:
        test_function_mat_compile(){};

        virtual Real eval_mat(const Matrix& mat,bool show_res,int code) override;
        virtual Real eval_scalar(const Scalar& s, bool show_res,int code) override;
};

class test_function_resize_reserve : public unary_function
{
    public:
        test_function_resize_reserve(){};

        virtual Real eval_mat(const Matrix& mat,bool show_res,int code) override;
        virtual Real eval_scalar(const Scalar& s, bool show_res,int code) override;

        static Real    check_struct_eq(struct_flag s1, struct_flag s2);
};

class test_function_resize_reserve_b : public unary_function
{
    public:
        test_function_resize_reserve_b(){};

        virtual Real eval_mat(const Matrix& mat,bool show_res,int code) override;
        virtual Real eval_scalar(const Scalar& s, bool show_res,int code) override;
};

class test_function_get_scal : public unary_function
{
    public:
        test_function_get_scal(){};

        virtual Real eval_mat(const Matrix& mat,bool show_res,int code) override;
        virtual Real eval_scalar(const Scalar& s, bool show_res,int code) override;
};

class test_function_get_scal_unique : public unary_function
{
    public:
        test_function_get_scal_unique(){};

        virtual Real eval_mat(const Matrix& mat,bool show_res,int code) override;
        virtual Real eval_scalar(const Scalar& s, bool show_res,int code) override;
};

class test_function_get_array : public unary_function
{
    public:
        test_function_get_array(){};

        virtual Real eval_mat(const Matrix& mat,bool show_res,int code) override;
        virtual Real eval_scalar(const Scalar& s, bool show_res,int code) override;
};

class test_function_get_const_array : public unary_function
{
    public:
        test_function_get_const_array(){};

        virtual Real eval_mat(const Matrix& mat,bool show_res,int code) override;
        virtual Real eval_scalar(const Scalar& s, bool show_res,int code) override;
};

class test_function_get_impl_unique : public unary_function
{
    public:
        test_function_get_impl_unique(){};

        virtual Real eval_mat(const Matrix& mat,bool show_res,int code) override;
        virtual Real eval_scalar(const Scalar& s, bool show_res,int code) override;
};

};};
