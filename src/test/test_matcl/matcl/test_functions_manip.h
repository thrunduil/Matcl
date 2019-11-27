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

class manip_functions_list
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
        void    test_diag();			void    test_bdiag();
        void    test_spdiag();			void    test_bdiags();
        void    test_spdiags();			void    test_diags();

        void    test_find();
        void    test_find2();			void    test_find3();

        void    test_sort();			void    test_sort2();
        void    test_sortrows();		void    test_sortrows2();
        void    test_sortrows_dim();	void    test_sortrows2_dim();
        void    test_sortcols();		void    test_sortcols2();
        void    test_sortcols_dim();	void    test_sortcols2_dim();
        void    test_is_sorted();		void    test_is_sorted_rows();
        void    test_is_sorted_cols();  

        void    test_tril();
        void    test_triu();			void    test_rot90();
        void    test_vec();				void    test_trans();
        void    test_ctrans();          void    test_flipud();
        void    test_fliplr();          void    test_reshape();
        void    test_repmat();			void    test_convert();			
        void    test_trans_t();

        void    test_sparse();			void    test_clone();
        void    test_band();            void    test_cat();
        void    test_get_diag();		void    test_get_lu();
        void    test_full();            void    test_matrix_fwd();
        void    test_convert_val();     void    test_horzcat();
        void    test_vertcat();         void    test_blkdiag();
        void    test_select_band();     void    test_checks();
        void    test_nnz_m();           void    test_drop_sparse();

        template<class Func>
        void    test_function();

    public:
        manip_functions_list(const matrix_set&,rand_matrix_ptr rand);

    private:
        manip_functions_list(const manip_functions_list&) = delete;
        manip_functions_list& operator=(const manip_functions_list&) = delete;
};

class test_function_diag: public unary_function
{
    private:
        int d;

    public:
        test_function_diag(int d)	:d(d){};

        virtual Real eval_mat(const Matrix& mat,bool show_res,int code) override;
        virtual Real eval_scalar(const Scalar& s, bool show_res,int code) override;
};

class test_function_bdiag: public unary_function
{
    private:
        int d;

    public:
        test_function_bdiag(int d)	:d(d){};

        virtual Real eval_mat(const Matrix& mat,bool show_res,int code) override;
        virtual Real eval_scalar(const Scalar& s, bool show_res,int code) override;
};

class test_function_spdiag: public unary_function
{
    private:
        int d;

    public:
        test_function_spdiag(int d)	:d(d){};

        virtual Real eval_mat(const Matrix& mat,bool show_res,int code) override;
        virtual Real eval_scalar(const Scalar& s, bool show_res,int code) override;
};

class test_function_diags: public unary_function
{
    public:
        test_function_diags()	{};

        virtual Real eval_mat(const Matrix& mat,bool show_res,int code) override;
        virtual Real eval_scalar(const Scalar& s, bool show_res,int code) override;
};

class test_function_bdiags: public unary_function
{
    public:
        test_function_bdiags()	{};

        virtual Real eval_mat(const Matrix& mat,bool show_res,int code) override;
        virtual Real eval_scalar(const Scalar& s, bool show_res,int code) override;
};

class test_function_spdiags: public unary_function
{
    public:
        test_function_spdiags()	{};

        virtual Real eval_mat(const Matrix& mat,bool show_res,int code) override;
        virtual Real eval_scalar(const Scalar& s, bool show_res,int code) override;
};

class test_function_find : public unary_function
{
    public:
        test_function_find(){};

        virtual Real eval_mat(const Matrix& mat,bool show_res,int code) override;
        virtual Real eval_scalar(const Scalar& s, bool show_res,int code) override;
};

class test_function_find2 : public unary_function
{
    public:
        test_function_find2(){};

        virtual Real eval_mat(const Matrix& mat,bool show_res,int code) override;
        virtual Real eval_scalar(const Scalar& s, bool show_res,int code) override;
};

class test_function_find3 : public unary_function
{
    public:
        test_function_find3(){};

        virtual Real eval_mat(const Matrix& mat,bool show_res,int code) override;
        virtual Real eval_scalar(const Scalar& s, bool show_res,int code) override;
};

class test_function_sort : public unary_function
{
    private:
        Integer d;
        bool	is_asc;

    public:
        test_function_sort(Integer d, bool is_asc)
            :d(d),is_asc(is_asc){};

        virtual Real eval_mat(const Matrix& mat,bool show_res,int code) override;
        virtual Real eval_scalar(const Scalar& s, bool show_res,int code) override;
};

class test_function_sort2 : public unary_function
{
    private:
        Integer d;
        bool	is_asc;

    public:
        test_function_sort2(Integer d, bool is_asc)	
            :d(d),is_asc(is_asc){};

        virtual Real eval_mat(const Matrix& mat,bool show_res,int code) override;
        virtual Real eval_scalar(const Scalar& s, bool show_res,int code) override;
};

class test_function_sortrows : public unary_function
{
    public:
        test_function_sortrows(){};

        virtual Real eval_mat(const Matrix& mat,bool show_res,int code) override;
        virtual Real eval_scalar(const Scalar& s, bool show_res,int code) override;
};

class test_function_sortrows2 : public unary_function
{
    public:
        test_function_sortrows2(){};

        virtual Real eval_mat(const Matrix& mat,bool show_res,int code) override;
        virtual Real eval_scalar(const Scalar& s, bool show_res,int code) override;
};

class test_function_sortrows_dim : public unary_function
{
    private:
        Matrix m;

    public:
        test_function_sortrows_dim(const Matrix& m)	:m(m){};

        virtual Real eval_mat(const Matrix& mat,bool show_res,int code) override;
        virtual Real eval_scalar(const Scalar& s, bool show_res,int code) override;
};

class test_function_sortrows_dim2 : public unary_function
{
    private:
        Matrix	m;

    public:
        test_function_sortrows_dim2(const Matrix& m)	:m(m){};

        virtual Real eval_mat(const Matrix& mat,bool show_res,int code) override;
        virtual Real eval_scalar(const Scalar& s, bool show_res,int code) override;
};

class test_function_sortcols : public unary_function
{
    public:
        test_function_sortcols(){};

        virtual Real eval_mat(const Matrix& mat,bool show_res,int code) override;
        virtual Real eval_scalar(const Scalar& s, bool show_res,int code) override;
};

class test_function_sortcols2 : public unary_function
{
    public:
        test_function_sortcols2(){};

        virtual Real eval_mat(const Matrix& mat,bool show_res,int code) override;
        virtual Real eval_scalar(const Scalar& s, bool show_res,int code) override;
};

class test_function_sortcols_dim : public unary_function
{
    private:
        Matrix	m;

    public:
        test_function_sortcols_dim(const Matrix& m)	:m(m){};

        virtual Real eval_mat(const Matrix& mat,bool show_res,int code) override;
        virtual Real eval_scalar(const Scalar& s, bool show_res,int code) override;
};

class test_function_sortcols_dim2 : public unary_function
{
    private:
        Matrix	m;

    public:
        test_function_sortcols_dim2(const Matrix& m)	:m(m){};

        virtual Real eval_mat(const Matrix& mat,bool show_res,int code) override;
        virtual Real eval_scalar(const Scalar& s, bool show_res,int code) override;
};

class test_function_is_sorted : public unary_function
{
    private:
        Integer d;
        bool is_asc;

    public:
        test_function_is_sorted(Integer d,bool is_asc)
            :d(d),is_asc(is_asc){};

        virtual Real eval_mat(const Matrix& mat,bool show_res,int code) override;
        virtual Real eval_scalar(const Scalar& s, bool show_res,int code) override;
};

class test_function_is_sorted_rows : public unary_function
{
    public:
        test_function_is_sorted_rows(){};

        virtual Real eval_mat(const Matrix& mat,bool show_res,int code) override;
        virtual Real eval_scalar(const Scalar& s, bool show_res,int code) override;
};

class test_function_is_sorted_cols : public unary_function
{
    public:
        test_function_is_sorted_cols(){};

        virtual Real eval_mat(const Matrix& mat,bool show_res,int code) override;
        virtual Real eval_scalar(const Scalar& s, bool show_res,int code) override;
};

class test_function_get_diag : public unary_function
{
    private:
        Integer d;

    public:
        test_function_get_diag(Integer d)
            : d(d) {};

        virtual Real eval_mat(const Matrix& mat,bool show_res,int code) override;
        virtual Real eval_scalar(const Scalar& s, bool show_res,int code) override;
};

class test_function_tril : public unary_function
{
    private:
        Integer d;

    public:
        test_function_tril(Integer d) : d(d) {};

        virtual Real eval_mat(const Matrix& mat,bool show_res,int code) override;
        virtual Real eval_scalar(const Scalar& s, bool show_res,int code) override;
};

class test_function_triu : public unary_function
{
    private:
        Integer d;

    public:
        test_function_triu(Integer d) : d(d) {};

        virtual Real eval_mat(const Matrix& mat,bool show_res,int code) override;
        virtual Real eval_scalar(const Scalar& s, bool show_res,int code) override;
};

class test_function_rot90 : public unary_function
{
    private:
        Integer d;

    public:
        test_function_rot90(Integer d) : d(d) {};

        virtual Real eval_mat(const Matrix& mat,bool show_res,int code) override;
        virtual Real eval_scalar(const Scalar& s, bool show_res,int code) override;
};

class test_function_reshape : public unary_function
{
    public:
        test_function_reshape(){};

        virtual Real eval_mat(const Matrix& mat,bool show_res,int code) override;
        virtual Real eval_scalar(const Scalar& s, bool show_res,int code) override;
};

class test_function_cat : public unary_function
{
    private:
        Integer         m;
        Integer         n;
        Integer         pos;
        Integer         type;
        Integer         seed;
        rand_matrix_ptr m_rand;

    public:
        test_function_cat(Integer m, Integer n, Integer pos, Integer type,
                        rand_matrix_ptr rand, Integer s) 
            : m(m), n(n), pos(pos), type(type), m_rand(rand), seed(s)
        {};

        virtual Real eval_mat(const Matrix& mat,bool show_res,int code) override;
        virtual Real eval_scalar(const Scalar& s, bool show_res,int code) override;
};

class test_function_sparse : public unary_function
{
    public:
        test_function_sparse()	{};

        virtual Real eval_mat(const Matrix& mat,bool show_res,int code) override;
        virtual Real eval_scalar(const Scalar& s, bool show_res,int code) override;
};

class test_function_full : public unary_function
{
    public:
        test_function_full()	{};

        virtual Real eval_mat(const Matrix& mat,bool show_res,int code) override;
        virtual Real eval_scalar(const Scalar& s, bool show_res,int code) override;
};

class test_function_band : public unary_function
{
    public:
        test_function_band()	{};

        virtual Real eval_mat(const Matrix& mat,bool show_res,int code) override;
        virtual Real eval_scalar(const Scalar& s, bool show_res,int code) override;
};

class test_function_clone : public unary_function
{
    public:
        test_function_clone()	{};

        virtual Real eval_mat(const Matrix& mat,bool show_res,int code) override;
        virtual Real eval_scalar(const Scalar& s, bool show_res,int code) override;
};

class test_function_repmat : public unary_function
{
    private:
        Integer m;
        Integer n;

    public:
        test_function_repmat(Integer m, Integer n) : m(m), n(n){};

        virtual Real eval_mat(const Matrix& mat,bool show_res,int code) override;
        virtual Real eval_scalar(const Scalar& s, bool show_res,int code) override;
};

class test_function_get_lu : public unary_function
{
    public:
        test_function_get_lu(){};

        virtual Real eval_mat(const Matrix& mat,bool show_res,int code) override;
        virtual Real eval_scalar(const Scalar& s, bool show_res,int code) override;
};

class test_function_trans_t : public unary_function
{
    public:
        test_function_trans_t(){};

        virtual Real eval_mat(const Matrix& mat,bool show_res,int code) override;
        virtual Real eval_scalar(const Scalar& s, bool show_res,int code) override;
};

class test_function_convert : public unary_function
{
    private:
        Integer code;

    public:
        test_function_convert(Integer code)	:code(code){};

        virtual Real eval_mat(const Matrix& mat,bool show_res,int code) override;
        virtual Real eval_scalar(const Scalar& s, bool show_res,int code) override;
};

class test_function_convert_val : public unary_function
{
    private:
        Integer code;

    public:
        test_function_convert_val(Integer code)	:code(code){};

        virtual Real eval_mat(const Matrix& mat,bool show_res,int code) override;
        virtual Real eval_scalar(const Scalar& s, bool show_res,int code) override;
};

class test_function_matrix_fwd : public unary_function
{
    public:
        test_function_matrix_fwd(){};

        virtual Real eval_mat(const Matrix& mat,bool show_res,int code) override;
        virtual Real eval_scalar(const Scalar& s, bool show_res,int code) override;
};

class test_function_horzcat : public unary_function
{
    public:
        test_function_horzcat(){};

        virtual Real eval_mat(const Matrix& mat,bool show_res,int code) override;
        virtual Real eval_scalar(const Scalar& s, bool show_res,int code) override;
};

class test_function_vertcat : public unary_function
{
    public:
        test_function_vertcat(){};

        virtual Real eval_mat(const Matrix& mat,bool show_res,int code) override;
        virtual Real eval_scalar(const Scalar& s, bool show_res,int code) override;
};

class test_function_blkdiag : public unary_function
{
    public:
        test_function_blkdiag(){};

        virtual Real eval_mat(const Matrix& mat,bool show_res,int code) override;
        virtual Real eval_scalar(const Scalar& s, bool show_res,int code) override;
};

class test_function_checks : public unary_function
{
    public:
        test_function_checks(){};

        virtual Real eval_mat(const Matrix& mat,bool show_res,int code) override;
        virtual Real eval_scalar(const Scalar& s, bool show_res,int code) override;
};

class test_function_select_band : public unary_function
{
    private:
        Integer m_fd;
        Integer m_ld;

    public:
        test_function_select_band(Integer fd, Integer ld)
            :m_fd(fd), m_ld(ld)
        {};

        virtual Real eval_mat(const Matrix& mat,bool show_res,int code) override;
        virtual Real eval_scalar(const Scalar& s, bool show_res,int code) override;
};

class test_function_nnz_m : public unary_function
{
    public:
        test_function_nnz_m(){};

        virtual Real eval_mat(const Matrix& mat,bool show_res,int code) override;
        virtual Real eval_scalar(const Scalar& s, bool show_res,int code) override;
};

class test_function_drop_sparse : public unary_function
{
    public:
        test_function_drop_sparse(){};

        virtual Real eval_mat(const Matrix& mat,bool show_res,int code) override;
        virtual Real eval_scalar(const Scalar& s, bool show_res,int code) override;
};

struct h_vec		
{ 
    static Matrix eval(const Matrix& m)			{ return vec(m); }; 
    template<class T>
    static Matrix eval_scal(T m)				{ return vec(m); }; 
};

struct h_trans
{ 
    static Matrix eval(const Matrix& m)			{ return trans(m); }; 
    template<class T>
    static Matrix eval_scal(T m)				{ return trans(m); }; 
};
struct h_ctrans
{ 
    static Matrix eval(const Matrix& m)			{ return ctrans(m); }; 
    template<class T>
    static Matrix eval_scal(T m)				{ return ctrans(m); }; 
};
struct h_flipud
{ 
    static Matrix eval(const Matrix& m)			{ return flipud(m); }; 
    template<class T>
    static Matrix eval_scal(T m)				{ return flipud(m); }; 
};
struct h_fliplr
{ 
    static Matrix eval(const Matrix& m)			{ return fliplr(m); }; 
    template<class T>
    static Matrix eval_scal(T m)				{ return fliplr(m); }; 
};

using tf    = unary_function;
using ccp   = const char*;

struct function_vec		{ ccp name() { return "vec";};		scal_func<h_vec> f;		tf* function() { return &f;}; };
struct function_trans	{ ccp name() { return "trans";};	scal_func<h_trans> f;	tf* function() { return &f;}; };
struct function_ctrans	{ ccp name() { return "ctrans";};	scal_func<h_ctrans> f;	tf* function() { return &f;}; };			
struct function_flipud	{ ccp name() { return "flipud";};	scal_func<h_flipud> f;	tf* function() { return &f;}; };
struct function_fliplr	{ ccp name() { return "fliplr";};	scal_func<h_fliplr> f;	tf* function() { return &f;}; };			

};};
