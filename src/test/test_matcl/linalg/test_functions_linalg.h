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

#include <vector>
#include <string>
#include "matcl-matrep/matcl_matrep.h"
#include "test/test_matcl/framework/matrix_set/matrix_set.h"
#include "matcl-linalg/matcl_linalg.h"

namespace matcl { namespace test
{

class linalg_functions_list
{
    private:
        const matrix_set&	m_tests;
        options				m_options;
        dynamic_mat_set&	ms;

    public:
        void		        make(options opts);
        Matrix		        get_matrix(int code) const;

        dynamic_mat_set&	get_ms()    { return ms; };

    private:
        void		test_lu_partial();
        void		test_lu_rook();
        void		test_lu_complete();

        void        test_linsolve_0_NT();
        void        test_linsolve_0_T();
        void        test_linsolve_0_CT();
        void        test_linsolve_1_NT();
        void        test_linsolve_1_T();
        void        test_linsolve_1_CT();
        void        test_linsolve_3_NT();
        void        test_linsolve_3_T();
        void        test_linsolve_3_CT();

        void        test_linsolve_rev();
        void        test_linsolve_rev2_0();
        void        test_linsolve_rev2_1();
        void        test_linsolve_rev2_3();

        void		test_norm();
        void		test_svd();
        void		test_chol();
        void		test_chol_rr();
        void		test_cholmod();
        void        test_qr();
        void        test_hess();
        void        test_schur();
        void        test_eigs();
        void        test_ldl();
        void        test_cond();

    public:
        linalg_functions_list(const matrix_set&,dynamic_mat_set& ms);

    private:
        linalg_functions_list(const linalg_functions_list&) = delete;
        linalg_functions_list& operator=(const linalg_functions_list&) = delete;
};

//----------------------------------------------------------------------
class test_function_lu_partial : public unary_function
{
    public:
        test_function_lu_partial() {};

        virtual Real eval_mat(const Matrix& mat,bool show_res,int code);
        virtual Real eval_scalar(const Scalar& s, bool show_res,int code);
};
class test_function_lu_rook : public unary_function
{
    public:
        test_function_lu_rook() {};

        virtual Real eval_mat(const Matrix& mat,bool show_res,int code);
        virtual Real eval_scalar(const Scalar& s, bool show_res,int code);
};
class test_function_lu_complete : public unary_function
{
    public:
        test_function_lu_complete() {};

        virtual Real eval_mat(const Matrix& mat,bool show_res,int code);
        virtual Real eval_scalar(const Scalar& s, bool show_res,int code);
};

class test_function_linsolve : public unary_function
{
    private:
        dynamic_mat_set&	ms;
        long                m_new_objects;
        long                K;
        trans_type          m_trans;

    public:
        test_function_linsolve(long k, trans_type t, dynamic_mat_set& ms) 
            :ms(ms), m_trans(t), K(k)
        {};

        virtual long n_new_objects()    { return m_new_objects; };

        virtual Real eval_mat(const Matrix& mat,bool show_res,int code);
        virtual Real eval_scalar(const Scalar& s, bool show_res,int code);

    private:
        test_function_linsolve(const test_function_linsolve&);
        test_function_linsolve& operator=(const test_function_linsolve&);

        Real        eval_mat_impl(const Matrix& A,const Matrix& B);
        Real        test_notrans(const Matrix& A,const Matrix& B, permvec p, permvec q);
        Real        test_trans(const Matrix& A,const Matrix& B, permvec p, permvec q);
        Real        test_ctrans(const Matrix& A,const Matrix& B, permvec p, permvec q);
};

class test_function_linsolve_rev : public unary_function
{
    private:
        dynamic_mat_set&	ms;
        long                m_new_objects;
        long                K;

    public:
        test_function_linsolve_rev(long k, dynamic_mat_set& ms) :ms(ms),K(k){};

        virtual long n_new_objects()                        { return m_new_objects; };

        virtual Real eval_mat(const Matrix& mat,bool show_res,int code);
        virtual Real eval_scalar(const Scalar& s, bool show_res,int code);

    private:
        test_function_linsolve_rev(const test_function_linsolve_rev&) = delete;
        test_function_linsolve_rev& operator=(const test_function_linsolve_rev&) = delete;

        Real        eval_mat_impl(const Matrix& A,const Matrix& B);
        Real        test_notrans(const Matrix& A,const Matrix& B, permvec p, permvec q);
};

class test_function_linsolve_rev2 : public unary_function
{
    private:
        dynamic_mat_set&	ms;
        long                m_new_objects;
        long                K;

    public:
        test_function_linsolve_rev2(long k, dynamic_mat_set& ms) :ms(ms),K(k){};

        virtual long n_new_objects()    { return m_new_objects; };

        virtual Real eval_mat(const Matrix& mat,bool show_res,int code);
        virtual Real eval_scalar(const Scalar& s, bool show_res,int code);

    private:
        test_function_linsolve_rev2(const test_function_linsolve_rev2&) = delete;
        test_function_linsolve_rev2& operator=(const test_function_linsolve_rev2&) = delete;

        Real        eval_mat_impl(const Matrix& A,const Matrix& B);
        Real        test_notrans(const Matrix& A,const Matrix& B, permvec p, permvec q);
        Real        test_trans(const Matrix& A,const Matrix& B, permvec p, permvec q);
        Real        test_ctrans(const Matrix& A,const Matrix& B, permvec p, permvec q);
};

class test_function_norm : public unary_function
{
    private:
        Real    m_type;

    public:
        test_function_norm(Real t)  :m_type(t){};

        virtual Real eval_mat(const Matrix& mat,bool show_res,int code);
        virtual Real eval_scalar(const Scalar& s, bool show_res,int code);
};

class test_function_svd : public unary_function
{
    public:
        test_function_svd()         {};

        virtual Real eval_mat(const Matrix& mat,bool show_res,int code);
        virtual Real eval_scalar(const Scalar& s, bool show_res,int code);

        Real eval_1(const Matrix& mat,int code, bool economy, svd_algorithm alg);
};

class test_function_chol : public unary_function
{
    public:
        test_function_chol()         {};

        virtual Real eval_mat(const Matrix& mat,bool show_res,int code);
        virtual Real eval_scalar(const Scalar& s, bool show_res,int code);

    private:
        Real         eval(const Matrix& mat, bool upper);
};
class test_function_chol_rr : public unary_function
{
    public:
        test_function_chol_rr()         {};

        virtual Real eval_mat(const Matrix& mat,bool show_res,int code);
        virtual Real eval_scalar(const Scalar& s, bool show_res,int code);

    private:
        Real         eval(const Matrix& mat, bool upper);
};

class test_function_cholmod : public unary_function
{
    public:
        test_function_cholmod()         {};

        virtual Real eval_mat(const Matrix& mat,bool show_res,int code);
        virtual Real eval_scalar(const Scalar& s, bool show_res,int code);

        Real eval(const Matrix& mat,bool upper, bool show);
};

class test_function_qr : public unary_function
{
    public:
        test_function_qr()         {};

        virtual Real eval_mat(const Matrix& mat,bool show_res,int code);
        virtual Real eval_scalar(const Scalar& s, bool show_res,int code);
};

class test_function_hess : public unary_function
{
    public:
        test_function_hess()         {};

        virtual Real eval_mat(const Matrix& mat,bool show_res,int code);
        virtual Real eval_scalar(const Scalar& s, bool show_res,int code);
};

class test_function_schur : public unary_function
{
    public:
        test_function_schur()         {};

        virtual Real eval_mat(const Matrix& mat,bool show_res,int code);
        virtual Real eval_scalar(const Scalar& s, bool show_res,int code);
};

class test_function_eigs : public unary_function
{
    public:
        test_function_eigs()         {};

        virtual Real eval_mat(const Matrix& mat,bool show_res,int code);
        virtual Real eval_scalar(const Scalar& s, bool show_res,int code);
};

class test_function_ldl : public unary_function
{
    public:
        test_function_ldl()         {};

        virtual Real eval_mat(const Matrix& mat,bool show_res,int code);
        virtual Real eval_scalar(const Scalar& s, bool show_res,int code);

        Real eval_1(const Matrix& mat, int code, bool upper);
};

class test_function_cond : public unary_function
{
    public:
        test_function_cond()         {};

        virtual Real eval_mat(const Matrix& mat,bool show_res,int code);
        virtual Real eval_scalar(const Scalar& s, bool show_res,int code);
};

};};
