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

class assign_functions_list
{
    private:
        const matrix_set&		m_tests;
        options					m_options;

    public:
        dynamic_mat_set&		ms;

    public:
        void		make(options opts);
        Matrix		get_matrix(int code) const;

    private:
        void    test_get_int();			    void    test_get_int_int();
        void    test_get_colon();		    void    test_get_colon_colon();
        void    test_get_diag();

        void    test_delcols();			    void    test_delrows();
        void    test_delrowscols();

        void    test_drop_int();		    void    test_drop_int_int();
        void    test_drop_colon();		    void    test_drop_colon_colon();
        void    test_drop_diag();

        void    test_del_int();			    void    test_del_int_int();
        void    test_del_colon();		    void    test_del_colon_colon();
        void    test_del_diag();

        void    test_set_int_scal();	    void    test_set_int2_scal();
        void    test_set_col_scal();	    void    test_set_col2_scal();
        void    test_set_diag_scal();

        void    test_set_int_mat();		    void    test_set_int2_mat();
        void    test_set_col_mat();		    void    test_set_col2_mat();        
        void    test_set_diag_mat();

        void    test_drop_sparse_int();     void    test_drop_sparse_int_int();
        void    test_drop_sparse_col();	    void    test_drop_sparse_col2();        
        void    test_drop_sparse_diag();

        void    test_add_sparse_int();      void    test_add_sparse_int_int();
        void    test_add_sparse_col();	    void    test_add_sparse_col2();        
        void    test_add_sparse_diag();

    public:
        assign_functions_list(const matrix_set& ms_bin, dynamic_mat_set& ms);

    private:
        assign_functions_list(const assign_functions_list&) = delete;
        assign_functions_list& operator=(const assign_functions_list&) = delete;
};

class test_function_get_int : public unary_function
{
    private:
        Integer i;

    public:
        test_function_get_int(Integer i)
            :i(i){};

        virtual Real eval_mat(const Matrix& mat,bool show_res,int code);
        virtual Real eval_scalar(const Scalar& s, bool show_res,int code);
};

class test_function_get_int_int : public unary_function
{
    private:
        Integer i;
        Integer j;

    public:
        test_function_get_int_int(Integer i,Integer j)
            :i(i),j(j){};

        virtual Real eval_mat(const Matrix& mat,bool show_res,int code);
        virtual Real eval_scalar(const Scalar& s, bool show_res,int code);
};

class test_function_get_colon : public unary_function
{
    public:
        test_function_get_colon()	{};

        virtual Real eval_mat(const Matrix& mat,bool show_res,int code);
        virtual Real eval_scalar(const Scalar& s, bool show_res,int code);

        Real eval_mat_col(const Matrix& mat, const colon& c1);
};
class test_function_get_colon_colon : public unary_function
{
    public:
        test_function_get_colon_colon(){};

        virtual Real eval_mat(const Matrix& mat,bool show_res,int code);
        virtual Real eval_scalar(const Scalar& s, bool show_res,int code);

        Real eval_mat_col(const Matrix& mat, const colon& c1, const colon& c2);
};

class test_function_get_diag_sub : public unary_function
{
    private:
        Integer d;

    public:
        test_function_get_diag_sub(Integer d) 
            : d(d) {};

        virtual Real eval_mat(const Matrix& mat,bool show_res,int code);
        virtual Real eval_scalar(const Scalar& s, bool show_res,int code);
};

class test_function_delcols : public unary_function
{
    public:
        test_function_delcols(){};

        virtual Real eval_mat(const Matrix& mat,bool show_res,int code);
        virtual Real eval_scalar(const Scalar& s, bool show_res,int code);

        Real eval_mat_col(const Matrix& mat, const colon& c);
};

class test_function_delrows : public unary_function
{
    public:
        test_function_delrows(){};

        virtual Real eval_mat(const Matrix& mat,bool show_res,int code);
        virtual Real eval_scalar(const Scalar& s, bool show_res,int code);

        Real eval_mat_col(const Matrix& mat, const colon& c);
};

class test_function_delrowscols : public unary_function
{
    public:
        test_function_delrowscols(){};

        virtual Real eval_mat(const Matrix& mat,bool show_res,int code);
        virtual Real eval_scalar(const Scalar& s, bool show_res,int code);

        Real eval_mat_col(const Matrix& mat, const colon& c1, const colon& c2);
};

class test_function_drop_int : public unary_function
{
    private:
        Integer i;

    public:
        test_function_drop_int(Integer i)
            :i(i){};

        virtual Real eval_mat(const Matrix& mat,bool show_res,int code);
        virtual Real eval_scalar(const Scalar& s, bool show_res,int code);
};

class test_function_drop_diag : public unary_function
{
    private:
        Integer d;

    public:
        test_function_drop_diag(Integer d)
            :d(d){};

        virtual Real eval_mat(const Matrix& mat,bool show_res,int code);
        virtual Real eval_scalar(const Scalar& s, bool show_res,int code);
};

class test_function_drop_int_int : public unary_function
{
    private:
        Integer i;
        Integer j;

    public:
        test_function_drop_int_int(Integer i,Integer j)
            :i(i),j(j){};

        virtual Real eval_mat(const Matrix& mat,bool show_res,int code);
        virtual Real eval_scalar(const Scalar& s, bool show_res,int code);
};

class test_function_drop_colon : public unary_function
{
    public:
        test_function_drop_colon(){};

        virtual Real eval_mat(const Matrix& mat,bool show_res,int code);
        virtual Real eval_scalar(const Scalar& s, bool show_res,int code);

        Real eval_mat_col(const Matrix& mat, const colon& c1);
};

class test_function_drop_colon_colon : public unary_function
{
    public:
        test_function_drop_colon_colon(){};

        virtual Real eval_mat(const Matrix& mat,bool show_res,int code);
        virtual Real eval_scalar(const Scalar& s, bool show_res,int code);

        Real eval_mat_col(const Matrix& mat, const colon& c1, const colon& c2);
};

class test_function_del_int : public unary_function
{
    private:
        Integer i;

    public:
        test_function_del_int(Integer i)
            :i(i){};

        virtual Real eval_mat(const Matrix& mat,bool show_res,int code);
        virtual Real eval_scalar(const Scalar& s, bool show_res,int code);
};

class test_function_del_int_int : public unary_function
{
    private:
        Integer i;
        Integer j;

    public:
        test_function_del_int_int(Integer i,Integer j)
            :i(i),j(j){};

        virtual Real eval_mat(const Matrix& mat,bool show_res,int code);
        virtual Real eval_scalar(const Scalar& s, bool show_res,int code);
};

class test_function_del_colon: public unary_function
{
    public:
        test_function_del_colon(){};

        virtual Real eval_mat(const Matrix& mat,bool show_res,int code);
        virtual Real eval_scalar(const Scalar& s, bool show_res,int code);

        Real eval_mat_col(const Matrix& mat, const colon& c1);
};

class test_function_del_colon_colon : public unary_function
{
    public:
        test_function_del_colon_colon(){};

        virtual Real eval_mat(const Matrix& mat,bool show_res,int code);
        virtual Real eval_scalar(const Scalar& s, bool show_res,int code);

        Real eval_mat_col(const Matrix& mat, const colon& c1, const colon& c2);
};

class test_function_del_diag : public unary_function
{
    private:
        Integer d;

    public:
        test_function_del_diag(Integer d)
            :d(d){};

        virtual Real eval_mat(const Matrix& mat,bool show_res,int code);
        virtual Real eval_scalar(const Scalar& s, bool show_res,int code);
};

class test_function_set_int_scal : public unary_function
{
    private:
        Integer i;

    public:
        test_function_set_int_scal(Integer i)
            :i(i){};

        virtual Real eval_mat(const Matrix& mat,bool show_res,int code);
        virtual Real eval_scalar(const Scalar& s, bool show_res,int code);
};

class test_function_set_int2_scal : public unary_function
{
    private:
        Integer i;
        Integer j;

    public:
        test_function_set_int2_scal(Integer i,Integer j)
            :i(i),j(j){};

        virtual Real eval_mat(const Matrix& mat,bool show_res,int code);
        virtual Real eval_scalar(const Scalar& s, bool show_res,int code);
};

class test_function_set_col_scal : public unary_function
{
    public:
        test_function_set_col_scal(){};

        virtual Real eval_mat(const Matrix& mat,bool show_res,int code);
        virtual Real eval_scalar(const Scalar& s, bool show_res,int code);

        Real eval_mat_col(const Matrix& mat, const colon& c1);
};

class test_function_set_col2_scal : public unary_function
{
    public:
        test_function_set_col2_scal(){};

        virtual Real eval_mat(const Matrix& mat,bool show_res,int code);
        virtual Real eval_scalar(const Scalar& s, bool show_res,int code);

        Real eval_mat_col(const Matrix& mat, const colon& c1, const colon& c2);
};

class test_function_set_int_mat : public unary_function
{
    private:
        Integer i;

    public:
        test_function_set_int_mat(Integer i)
            :i(i){};

        virtual Real eval_mat(const Matrix& mat,bool show_res,int code);
        virtual Real eval_scalar(const Scalar& s, bool show_res,int code);
};

class test_function_set_int2_mat : public unary_function
{
    private:
        Integer i;
        Integer j;

    public:
        test_function_set_int2_mat(Integer i,Integer j)
            :i(i),j(j){};

        virtual Real eval_mat(const Matrix& mat,bool show_res,int code);
        virtual Real eval_scalar(const Scalar& s, bool show_res,int code);
};

class test_function_set_col_mat : public unary_function
{
    private:
        dynamic_mat_set&	ms;
        long                m_new_objects;

    public:
        test_function_set_col_mat(dynamic_mat_set& ms)	
            :ms(ms),m_new_objects(0)
        {};

        virtual long n_new_objects()             { return m_new_objects; };
        virtual Real eval_mat(const Matrix& mat,bool show_res,int code);
        virtual Real eval_scalar(const Scalar& s, bool show_res,int code);

        Real eval_mat_col(const Matrix& mat, const colon& c1);

    private:
        test_function_set_col_mat(const test_function_set_col_mat&) = delete;
        test_function_set_col_mat& operator=(const test_function_set_col_mat&) = delete;
};

class test_function_set_diag_scal : public unary_function
{
    private:
        Integer d;

    public:
        test_function_set_diag_scal(Integer d)
            :d(d){};

        virtual Real eval_mat(const Matrix& mat,bool show_res,int code);
        virtual Real eval_scalar(const Scalar& s, bool show_res,int code);
};

class test_function_set_col2_mat : public unary_function
{
    private:
        dynamic_mat_set&	ms;
        long                m_new_objects;

    public:
        test_function_set_col2_mat(dynamic_mat_set& ms)
            :ms(ms)
        {};

        virtual long n_new_objects()             { return m_new_objects; };

        virtual Real eval_mat(const Matrix& mat,bool show_res,int code);
        virtual Real eval_scalar(const Scalar& s, bool show_res,int code);

        Real eval_mat_col(const Matrix& mat, const colon& c1, const colon& c2);

    private:
        test_function_set_col2_mat(const test_function_set_col2_mat&) = delete;
        test_function_set_col2_mat& operator=(const test_function_set_col2_mat&) = delete;
};

class test_function_set_diag_mat : public unary_function
{
    private:
        Integer		        d;
        dynamic_mat_set&	ms;
        long                m_new_objects;

    public:
        test_function_set_diag_mat(Integer d,dynamic_mat_set& ms)	
            :d(d),ms(ms),m_new_objects(0)
        {};

        virtual long n_new_objects()             { return m_new_objects; };

        virtual Real eval_mat(const Matrix& mat,bool show_res,int code);
        virtual Real eval_scalar(const Scalar& s, bool show_res,int code);

    private:
        test_function_set_diag_mat(const test_function_set_diag_mat&) = delete;
        test_function_set_diag_mat& operator=(const test_function_set_diag_mat&) = delete;
};

class test_function_drop_sparse_int : public unary_function
{
    private:
        Integer i;
        Real    tol;

    public:
        test_function_drop_sparse_int(Integer i, Real tol)
            :i(i), tol(tol){};

        virtual Real eval_mat(const Matrix& mat,bool show_res,int code);
        virtual Real eval_scalar(const Scalar& s, bool show_res,int code);
};

class test_function_drop_sparse_int_int : public unary_function
{
    private:
        Integer i;
        Integer j;
        Real    tol;

    public:
        test_function_drop_sparse_int_int(Integer i,Integer j, Real tol)
            :i(i),j(j),tol(tol){};

        virtual Real eval_mat(const Matrix& mat,bool show_res,int code);
        virtual Real eval_scalar(const Scalar& s, bool show_res,int code);
};

class test_function_drop_sparse_colon: public unary_function
{
    private:
        Real    tol;

    public:
        test_function_drop_sparse_colon(Real tol)
            :tol(tol){};

        virtual Real eval_mat(const Matrix& mat,bool show_res,int code);
        virtual Real eval_scalar(const Scalar& s, bool show_res,int code);

        Real eval_mat_col(const Matrix& mat, const colon& c1);
};

class test_function_drop_sparse_colon_colon : public unary_function
{
    private:
        Real    tol;

    public:
        test_function_drop_sparse_colon_colon(Real tol)
            :tol(tol){};

        virtual Real eval_mat(const Matrix& mat,bool show_res,int code);
        virtual Real eval_scalar(const Scalar& s, bool show_res,int code);

        Real eval_mat_col(const Matrix& mat, const colon& c1, const colon& c2);
};

class test_function_drop_sparse_diag : public unary_function
{
    private:
        Integer d;
        Real    tol;

    public:
        test_function_drop_sparse_diag(Integer d, Real tol)
            :d(d), tol(tol){};

        virtual Real eval_mat(const Matrix& mat,bool show_res,int code);
        virtual Real eval_scalar(const Scalar& s, bool show_res,int code);
};

//
class test_function_add_sparse_int : public unary_function
{
    private:
        Integer i;

    public:
        test_function_add_sparse_int(Integer i)
            :i(i){};

        virtual Real eval_mat(const Matrix& mat,bool show_res,int code);
        virtual Real eval_scalar(const Scalar& s, bool show_res,int code);
};

class test_function_add_sparse_int_int : public unary_function
{
    private:
        Integer i;
        Integer j;

    public:
        test_function_add_sparse_int_int(Integer i,Integer j)
            :i(i),j(j){};

        virtual Real eval_mat(const Matrix& mat,bool show_res,int code);
        virtual Real eval_scalar(const Scalar& s, bool show_res,int code);
};

class test_function_add_sparse_colon: public unary_function
{
    public:
        test_function_add_sparse_colon(){};

        virtual Real eval_mat(const Matrix& mat,bool show_res,int code);
        virtual Real eval_scalar(const Scalar& s, bool show_res,int code);

        Real eval_mat_col(const Matrix& mat, const colon& c1);
};

class test_function_add_sparse_colon_colon : public unary_function
{
    public:
        test_function_add_sparse_colon_colon(){};

        virtual Real eval_mat(const Matrix& mat,bool show_res,int code);
        virtual Real eval_scalar(const Scalar& s, bool show_res,int code);

        Real eval_mat_col(const Matrix& mat, const colon& c1, const colon& c2);
};

class test_function_add_sparse_diag : public unary_function
{
    private:
        Integer d;

    public:
        test_function_add_sparse_diag(Integer d)
            :d(d){};

        virtual Real eval_mat(const Matrix& mat,bool show_res,int code);
        virtual Real eval_scalar(const Scalar& s, bool show_res,int code);
};

};};
