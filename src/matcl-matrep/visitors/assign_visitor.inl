/*
 *  This file is a part of Matrix Computation Library (MATCL)
 *
 *  Copyright (c) Pawe³ Kowal 2017 - 2021
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

#include "matcl-matrep/visitors/assign_visitor_utils.h"
#include "matcl-matrep/visitors/assign_visitor_mat.inl"
#include "matcl-matrep/algs/dense_algs.h"
#include "matcl-matrep/algs/sparse_algs.h"
#include "matcl-matrep/algs/banded_algs.h"
#include "matcl-internals/func/converter.h"
#include "matcl-internals/base/optim_params.h"
#include "matcl-scalar/details/matfunc_helpers.h"

#pragma warning( push )
#pragma warning(disable:4127)	// conditional expression is constant

namespace matcl { namespace details
{
    namespace mr = matcl::raw;
    namespace mrd = matcl::raw::details;

template<class M1>
struct assign_functor_drop
{
	static Matrix& eval(Matrix& A, M1& mat, Integer i,Integer j, Real tol);
	static Matrix& eval(Matrix& A,M1& mat, Integer i, Real tol);
	static Matrix& eval(Matrix& A, M1& mat, const colon& c1, Real tol);
	static Matrix& eval(Matrix& A, M1& mat, const colon& c1,const colon& c2, Real tol);
	static Matrix& eval_diag(Matrix& A,M1& mat, Integer d, Real tol);
};

template<class M1>
struct assign_functor_add
{
	static Matrix& eval(Matrix& A, M1& mat, Integer i,Integer j);
	static Matrix& eval(Matrix& A,M1& mat, Integer i);
	static Matrix& eval(Matrix& A, M1& mat, const colon& c1);
	static Matrix& eval(Matrix& A, M1& mat, const colon& c1,const colon& c2);
	static Matrix& eval_diag(Matrix& A,M1& mat, Integer d);
};

template<class M1,class val_type, bool req_promo, class struct_type>
struct assign_impl
{};

template<class V, class S>
struct assign_drop_impl
{
    using M1    = raw::Matrix<V,S>;

    static Matrix& eval(Matrix& A, M1& mat, Integer i, Integer j, Real)
    {
        error::check_index(i, j, mat.rows(), mat.cols());
        return A; 
    };

    static Matrix& eval(Matrix& A, M1& mat, Integer i, Real)
    { 
        error::check_single_index(i, mat.rows(), mat.cols());
        return A; 
    };

    static Matrix& eval(Matrix& A, M1& mat, const colon& c1, Real)
    { 
        colon_info ci;
        make_index(mat.rows(),mat.cols(),c1,ci);		

        return A; 
    };

    static Matrix& eval(Matrix& A, M1& mat, const colon& c1, const colon& c2, Real)
    { 
        colon_info ci;
        make_index(mat.rows(),mat.cols(),c1,c2,ci);

        return A; 
    };

    static Matrix& eval_diag(Matrix& A, M1& mat, Integer d, Real)
    { 
        error::check_diag(d, mat.rows(), mat.cols());
        return A; 
    };
};

template<class V, class S>
struct assign_add_impl
{
    using M1    = raw::Matrix<V,S>;

    static Matrix& eval(Matrix& A, M1& mat, Integer i,Integer j)
    { 
        error::check_index(i, j, mat.rows(), mat.cols());
        return A; 
    };

    static Matrix& eval(Matrix& A, M1& mat, Integer i)
    { 
        error::check_single_index(i, mat.rows(), mat.cols());
        return A; 
    };

    static Matrix& eval(Matrix& A, M1& mat, const colon& c1)
    {
        colon_info ci;
        make_index(mat.rows(),mat.cols(),c1,ci);	

        return A; 
    };

    static Matrix& eval(Matrix& A, M1& mat, const colon& c1,const colon& c2)
    { 
        colon_info ci;
        make_index(mat.rows(),mat.cols(),c1,c2,ci);

        return A; 
    };

    static Matrix& eval_diag(Matrix& A, M1& mat, Integer d)
    { 
        error::check_diag(d, mat.rows(), mat.cols());
        return A; 
    };
};

template<class V>
struct assign_drop_impl<V, struct_sparse>
{
    using M1    = raw::Matrix<V,struct_sparse>;

	static Matrix& eval(Matrix& A, M1& mat, Integer i,Integer j, Real tol)
    {
        error::check_index(i, j, mat.rows(), mat.cols());

        mrd::sparse_ccs<V>& rep = mat.rep();

        --i;
        --j;

        Integer k;
        bool has_element = rep.has_element(i,j,k);

        if (has_element == false)
            return A;

        const V& el     = rep.ptr_x()[k];

        auto aelem      = mrd::abs_helper<V>::eval(el);
        using AW        = decltype(aelem);
        bool need_drop  = (bool)mrd::leq_helper<AW,Real>::eval(aelem, tol); 

        if (need_drop == false)
            return A;

        rep.remove_element(i,j,k);
        return A;
    };

	static Matrix& eval(Matrix& A,M1& mat, Integer i, Real tol)
    {
        Integer r = mat.rows();
        if (r == 0)
            error::check_index(i, 0);

        Integer j;
        pos2ind(i,r,i,j);

        return eval(A,mat,i+1,j+1,tol);
    };

	static Matrix& eval(Matrix& A, M1& mat, const colon& c1, Real tol)
    {
        colon_info ci;
        make_index(mat.rows(),mat.cols(),c1,ci);		

        if (ci.rows() == 0)
            return A;

        if (tol < 0)
            return A;

        //struct optimization should not be allowed
        A = Matrix(algorithm::drop_entries_2(mat,ci, tol),false);
        return A;
    };

	static Matrix& eval(Matrix& A, M1& mat, const colon& c1,const colon& c2, Real tol)
    {
        colon_info ci;
        make_index(mat.rows(),mat.cols(),c1,c2,ci);

        Integer nr = ci.rows(), nc = ci.cols();

        if (nr == 0 || nc == 0 || tol < 0)
            return A;

        //struct optimization should not be allowed
        A = Matrix(algorithm::drop_entries(mat,ci, tol),false);
        return A;
    };

	static Matrix& eval_diag(Matrix& A,M1& mat, Integer d, Real tol)
    {
        //struct optimization should not be allowed
        A = Matrix(algorithm::sparse_drop_diag(mat,d,tol),false);
        return A;
    }
};

template<class V>
struct assign_add_impl<V, struct_sparse>
{
    using M1    = raw::Matrix<V,struct_sparse>;

	static Matrix& eval(Matrix& A, M1& mat, Integer i,Integer j)
    {
        error::check_index(i, j, mat.rows(), mat.cols());

        mrd::sparse_ccs<V>& rep = mat.rep();

        --i;
        --j;

        Integer k;
        bool has_element = rep.has_element(i,j,k);

        if (has_element == true)
            return A;

        rep.add_element(i,j,k);
        return A;
    };

	static Matrix& eval(Matrix& A,M1& mat, Integer i)
    {
        Integer r = mat.rows();
        if (r == 0)
            error::check_index(i, 0);

        Integer j;
        pos2ind(i,r,i,j);

        return eval(A,mat,i+1,j+1);
    }

	static Matrix& eval(Matrix& A, M1& mat, const colon& c1)
    {
        colon_info ci;
        make_index(mat.rows(),mat.cols(),c1,ci);		

        if (ci.rows() == 0)
            return A;

        //struct optimization should not be allowed
        A = Matrix(algorithm::add_entries_2(mat,ci),false);
        return A;
    };

	static Matrix& eval(Matrix& A, M1& mat, const colon& c1,const colon& c2)
    {
        colon_info ci;
        make_index(mat.rows(),mat.cols(),c1,c2,ci);

        Integer nr = ci.rows(), nc = ci.cols();

        if (nr == 0 || nc == 0)
            return A;
        
        //struct optimization should not be allowed
        A = Matrix(algorithm::add_entries(mat,ci),false);
        return A;
    }

	static Matrix& eval_diag(Matrix& A,M1& mat, Integer d)
    {
        //struct optimization should not be allowed
        A = Matrix(algorithm::sparse_add_diag(mat,d),false);
        return A;
    };
};

template<class M1,class val_type, class struct_type>
struct assign_impl<M1,val_type,true,struct_type>
{
    using val_type_1    = typename M1::value_type;
    using val_type_ret  = typename unify_types<val_type,val_type_1>::type;
    using new_mat_type  = raw::Matrix<val_type_ret,struct_type>;

    static Matrix& eval(Matrix& A,M1& mat,const val_type& val,Integer i,Integer j)
    {
        ti::ti_type<val_type_1> ti_1        = get_matrix_ti<val_type_1>::eval(A);
        ti::ti_type<val_type_ret> ret_ti    = link_ret_ti<val_type_ret,val_type_1,val_type>
                                                ::eval(ti_1,val);

        A                   = Matrix(raw::converter<new_mat_type,M1>::eval(mat,ret_ti),false);
        val_type_ret val2   = raw::converter<val_type_ret,val_type>::eval(val);

        return assign_impl<new_mat_type,val_type_ret,false,struct_type>
                    ::eval(A,A.get_impl_unique<new_mat_type>(),val2,i,j);
    };

    static Matrix& eval(Matrix& A,M1& mat,const val_type& val,Integer i)
    {
        ti::ti_type<val_type_1> ti_1        = get_matrix_ti<val_type_1>::eval(A);
        ti::ti_type<val_type_ret> ret_ti    = link_ret_ti<val_type_ret,val_type_1,val_type>
                                                ::eval(ti_1,val);
        A                   = Matrix(raw::converter<new_mat_type,M1>::eval(mat,ret_ti),false);
        val_type_ret val2   = raw::converter<val_type_ret,val_type>::eval(val);

        return assign_impl<new_mat_type,val_type_ret,false,struct_type>
                    ::eval(A,A.get_impl_unique<new_mat_type>(),val2,i);
    };

    static Matrix& eval(Matrix& A,M1& mat,const val_type& val,const colon& c1)
    {
        ti::ti_type<val_type_1> ti_1        = get_matrix_ti<val_type_1>::eval(A);
        ti::ti_type<val_type_ret> ret_ti    = link_ret_ti<val_type_ret,val_type_1,val_type>
                                                ::eval(ti_1,val);

        A                   = Matrix(raw::converter<new_mat_type,M1>::eval(mat,ret_ti),false);
        val_type_ret val2   = raw::converter<val_type_ret,val_type>::eval(val);

        return assign_impl<new_mat_type,val_type_ret,false,struct_type>
                    ::eval(A,A.get_impl_unique<new_mat_type>(),val2,c1);
    };

    static Matrix& eval(Matrix& A,M1& mat,const val_type& val,const colon& c1,const colon& c2)
    {
        ti::ti_type<val_type_1> ti_1        = get_matrix_ti<val_type_1>::eval(A);
        ti::ti_type<val_type_ret> ret_ti    = link_ret_ti<val_type_ret,val_type_1,val_type>
                                                ::eval(ti_1,val);
        A                   = Matrix(raw::converter<new_mat_type,M1>::eval(mat,ret_ti),false);
        val_type_ret val2   = raw::converter<val_type_ret,val_type>::eval(val);

        return assign_impl<new_mat_type,val_type_ret,false,struct_type>
                    ::eval(A,A.get_impl_unique<new_mat_type>(),val2,c1,c2);
    };

    static Matrix& eval_diag(Matrix& A,M1& mat,const val_type& val,Integer d)
    {
        ti::ti_type<val_type_1> ti_1        = get_matrix_ti<val_type_1>::eval(A);
        ti::ti_type<val_type_ret> ret_ti    = link_ret_ti<val_type_ret,val_type_1,val_type>
                                                ::eval(ti_1,val);
        A                   = Matrix(raw::converter<new_mat_type,M1>::eval(mat,ret_ti),false);
        val_type_ret val2   = raw::converter<val_type_ret,val_type>::eval(val);

        return assign_impl<new_mat_type,val_type_ret,false,struct_type>
                    ::eval_diag(A,A.get_impl_unique<new_mat_type>(),val2,d);
    };
};

template<class M1,class val_type>
struct assign_impl<M1,val_type,false,struct_dense>
{
    using val_type_ret  = typename M1::value_type;

    static Matrix& eval(Matrix& A,M1& mat, const val_type& val,Integer i,Integer j)
    {
        Integer r = mat.rows();
        Integer c = mat.cols();

        if (i < 1 || j < 1 || i > r || j > c)
            throw error::invalid_index(i, j, r, c);

        mrd::assign_helper(mat.ptr()[i-1+(j-1)*mat.ld()],val);

        mat.get_struct().reset();
        return A;
    };

    static Matrix& eval(Matrix& A, M1& mat, const val_type& val, Integer i)
    {
        Integer r = mat.rows();
        Integer c = mat.cols();

        if (r == 0)
            error::check_index(i,0);

        Integer j;
        pos2ind(i,r,i,j);

        if (i < 0 || j < 0 || i >= r || j >= c)
            throw error::invalid_index(i+1, j+1, r, c);

        mrd::assign_helper(mat.ptr()[i+j*mat.ld()],val);

        mat.get_struct().reset();
        return A;
    };

    static Matrix& eval(Matrix& A,M1& mat, const val_type& val,const colon& c1)
    {
        colon_info c_info;
        make_index(mat.rows(),mat.cols(), c1, c_info);        

        Integer s = c_info.rows();
        Integer r = mat.rows();

        if (s == 0)
            return A;

        typename M1::value_type* ptr_mat    = mat.ptr();
        Integer mat_ld                      = mat.ld();        

        if (c_info.r_flag == 0)
        {
            bool single             = c_info.is_double_mat_colon() == false;

            if (single)
            {
                mr::integer_dense ci_ri = c_info.get_rim_1();
                const Integer* ptr_ri   = ci_ri.ptr();

                for (Integer i = 0; i < s; ++i)
                {
                    Integer pos = ptr_ri[i];

                    Integer im, jm;
                    pos2ind(pos,r,im,jm);
                    mrd::assign_helper(ptr_mat[im+jm*mat_ld],val);
                };
            }
            else
            {
                mr::integer_dense ri    = c_info.get_rim_r();
                mr::integer_dense ci    = c_info.get_rim_c();
                const Integer* ptr_ri   = ri.ptr();
                const Integer* ptr_ci   = ci.ptr();

                for (Integer i = 0; i < s; ++i)
                {
                    Integer im  = ptr_ri[i] - 1;
                    Integer jm  = ptr_ci[i] - 1;
                    mrd::assign_helper(ptr_mat[im+jm*mat_ld],val);
                };
            };
        }
        else
        {
            if (mat.ld() == mat.rows())
            {
                for (Integer i = 0, pos = c_info.r_start-1; i < c_info.r_size; ++i, pos += c_info.r_step)
                    mrd::assign_helper(ptr_mat[pos],val);
            }
            else
            {
                for (Integer i = 0, pos = c_info.r_start-1; i < c_info.r_size; ++i, pos += c_info.r_step)
                {
                    while (pos >= r)
                    {
                        pos -= r;
                        ptr_mat += mat_ld;
                    };
                    mrd::assign_helper(ptr_mat[pos],val);
                };
            };
        };

        mat.get_struct().reset();
        return A;		
    };

    static Matrix& eval(Matrix& A,M1& mat, const val_type& val,const colon& c1,const colon& c2)
    {
        colon_info ci;
        make_index(mat.rows(),mat.cols(),c1,c2,ci);

        Integer nr = ci.rows(), nc = ci.cols();

        if (nr == 0 || nc == 0)
            return A;

        algorithm::dense_change_entries(mat,ci,val);
        return A;
    };

    static Matrix& eval_diag(Matrix& A,M1& mat, const val_type& val,Integer d)
    {
        algorithm::dense_change_diag<val_type>(mat,d,val);
        return A;
    };
};

template<class M1,class val_type>
struct assign_impl<M1,val_type,false,struct_sparse>
{
    using val_type_ret  = typename M1::value_type;

    static Matrix& eval(Matrix& A,M1& mat,const val_type& val,Integer i,Integer j)
    {
        error::check_index(i, j, mat.rows(), mat.cols());

        using value_type = typename M1::value_type;
        raw::details::sparse_ccs<value_type>& rep = mat.rep();

        --i;
        --j;

        Integer k;
        bool has_element = rep.has_element(i,j,k);

        if (mrd::is_zero(val))
        {
            if (has_element == false)
            {
                return A;
            }
            else
            {
                rep.remove_element(i,j,k);
                mat.get_struct().reset();
                return A;
            }
        }

        if (has_element == false)
        {
            rep.add_element(i,j,k);
            mrd::assign_helper(rep.ptr_x()[k],val);
            mat.get_struct().reset();
            return A;
        }
        else
        {
            mrd::assign_helper(rep.ptr_x()[k],val);
            mat.get_struct().reset();
            return A;
        };
    };

    static Matrix& eval(Matrix& A,M1& mat,const val_type& val,Integer i)
    {
        Integer r = mat.rows();
        if (r == 0)
            error::check_index(i, 0);

        Integer j;
        pos2ind(i,r,i,j);
        return eval(A,mat,val,i+1,j+1);
    };

    static Matrix& eval(Matrix& A,M1& mat,const val_type& val,const colon& c1)
    {
        colon_info ci;
        make_index(mat.rows(),mat.cols(),c1,ci);		

        if (ci.rows() == 0)
            return A;

        if (mrd::is_zero(val) == true)
        {
            A = Matrix(algorithm::zero_entries_2(mat,ci),true);
            return A;
        }

        Real density = (Real(mat.nnz()) + Real(ci.rows())) / (mat.rows()+1.) / (mat.cols()+1.);
        if (density > optim_params::max_sparse_density_max)
        {
            using value_type    = typename M1::value_type;
            using matrix_type   = raw::Matrix<value_type,struct_dense>;

            A               = matcl::full(Matrix(mat,false));
            matrix_type& B  = A.get_impl_unique<matrix_type>();

            algorithm::dense_change_entries_2(B,ci,val);
            return A;
        };

        A = Matrix(algorithm::change_entries_2(mat,ci,val),true);
        return A;
    };

    static Matrix& eval(Matrix& A, M1& mat,const val_type& val,const colon& c1,const colon& c2)
    {
        colon_info ci;
        make_index(mat.rows(),mat.cols(),c1,c2,ci);

        Integer nr = ci.rows(), nc = ci.cols();

        if (nr == 0 || nc == 0)
            return A;

        if (mrd::is_zero(val) == true)
        {
            A = Matrix(algorithm::zero_entries(mat,ci),true);
            return A;
        }
        
        Real density = (Real(mat.nnz()) + Real(ci.rows()) *Real(ci.cols()))
                    /(mat.rows()+1.)/(mat.cols()+1.);
        if (density > optim_params::max_sparse_density_max)
        {
            using value_type    = typename M1::value_type;
            using matrix_type   = raw::Matrix<value_type,struct_dense>;

            A               = matcl::full(Matrix(mat,false));
            matrix_type& B  = A.get_impl_unique<matrix_type>();

            algorithm::dense_change_entries(B,ci,val);
            return A;
        };

        A = Matrix(algorithm::change_entries(mat,ci,val),true);
        return A;
    };

    static Matrix& eval_diag(Matrix& A,M1& mat,const val_type& val,Integer d)
    {
        if (mrd::is_zero(val) == true)
        {
            A = Matrix(algorithm::zero_entries_diag(mat,d),true);
            return A;
        }

        A = Matrix(algorithm::sparse_change_diag(mat,d,val),true);
        return A;
    };
};

template<class M1,class val_type>
struct assign_impl<M1,val_type,false,struct_banded>
{
    using val_type_ret  = typename M1::value_type;

    static Matrix& eval(Matrix& A, M1& mat, const val_type& val, Integer i)
    {
        Integer r = mat.rows();
        if (r == 0)
            error::check_index(i, 0);

        Integer j;
        pos2ind(i,r,i,j);
        return eval(A,mat,val,i+1,j+1);
    };

    static Matrix& eval(Matrix& A,M1& mat,const val_type& val,Integer i,Integer j)
    {
        Integer r = mat.rows();
        Integer c = mat.cols();

        --i;
        --j;

        if (i < 0 || j < 0 || i >= r || j >= c)
            throw error::invalid_index(i+1, j+1, r, c);

        Integer first_row = mat.first_row(j);
        Integer last_row = mat.last_row(j);
        
        if ( i < first_row || i > last_row )
        {
            using value_type    = typename M1::value_type;
            using matrix_type   = raw::Matrix<value_type,struct_dense>;

            A = Matrix(raw::converter<matrix_type,M1>::eval(mat,get_matrix_ti<value_type>::eval(A)),false);
            matrix_type& mat2 = A.get_impl_unique<matrix_type>();

            mrd::assign_helper(mat2.ptr()[i+j*mat2.ld()],val);

            mat2.get_struct().reset();
            return A;
        }
      
        mrd::assign_helper(mat.rep_ptr()[mat.element_pos(i,j)], val);

        mat.get_struct().reset();
        return A;
    };

    static Matrix& eval(Matrix& A,M1& mat,const val_type& val,const colon& c1)
    {
        colon_info ci;

        make_index(mat.rows(),mat.cols(),c1,ci);
        Integer s   = ci.rows();
        Integer r   = mat.rows();

        if (s == 0)
            return A;

        Integer i   = 0;
        bool b_cont = true;
        bool single = ci.is_double_mat_colon() == false;

        if (ci.r_flag == 0)
        {            
            if (single == true)
            {
                mr::integer_dense ci_ri = ci.get_rim_1();
                const Integer* ptr_ri = ci_ri.ptr();

                while (b_cont && i<s)
                {				
                    b_cont = set(mat,ptr_ri[i],val);
                    ++i;
                };
            }
            else
            {
                mr::integer_dense rim   = ci.get_rim_r();
                mr::integer_dense cim   = ci.get_rim_c();
                const Integer* ptr_ri   = rim.ptr();
                const Integer* ptr_ci   = cim.ptr();

                while (b_cont && i < s)
                {				
                    b_cont = set(mat,ptr_ri[i],ptr_ci[i],val);
                    ++i;
                };
            }
        }
        else
        {
            Integer pos = ci.r_start;
            while (b_cont && i<s)
            {				
                b_cont = set(mat,pos,val);
                pos += ci.r_step;
                ++i;
            };
        };

        if (b_cont == true)
        {
            mat.get_struct().reset();
            return A;
        };

        using value_type    = typename M1::value_type;
        using matrix_type   = raw::Matrix<value_type,struct_dense>;

        A = Matrix(raw::converter<matrix_type,M1>::eval(A.get_impl<M1>()),false);
        matrix_type& B      = A.get_impl_unique<matrix_type>();
        value_type* ptr_B   = B.ptr();
        Integer B_ld        = B.ld();

        --i;

        if (ci.r_flag == 0)
        {
            if (single == true)
            {
                mr::integer_dense ri    = ci.get_rim_1();
                const Integer* ptr_ri   = ri.ptr();

                for (; i < s; ++i)
                {
                    Integer pos = ptr_ri[i];

                    Integer im, jm;
                    pos2ind(pos,r,im,jm);
                    mrd::assign_helper(ptr_B[im+jm*B_ld],val);
                };		
            }
            else
            {
                mr::integer_dense rim   = ci.get_rim_r();
                mr::integer_dense cim   = ci.get_rim_c();
                const Integer* ptr_ri   = rim.ptr();
                const Integer* ptr_ci   = cim.ptr();

                for (; i < s; ++i)
                {
                    Integer im          = ptr_ri[i] - 1;
                    Integer jm          = ptr_ci[i] - 1;
                    mrd::assign_helper(ptr_B[im+jm*B_ld],val);
                };		
            }
        }
        else
        {            
            Integer pos = ci.r_start - 1 + imult(i,ci.r_step);
            for (; i<s; ++i)
            {
                mrd::assign_helper(ptr_B[pos],val);
                pos += ci.r_step;
            };	
        };

        B.get_struct().reset();
        return A;
    };

    static Matrix& eval(Matrix& A,M1& mat,const val_type& val,const colon& c1,const colon& c2)
    {
        using t_val_type    = typename M1::value_type;

        colon_info ci;
        make_index(mat.rows(),mat.cols(),c1,c2,ci);

        Integer nr = ci.rows(), nc = ci.cols();

        if (nr == 0 || nc == 0)
            return A;

        Integer i = 0, j = 0;
        Integer mat_ld  = mat.ld();

        if (ci.r_flag == 0 && ci.c_flag == 0)
        {			
            mr::integer_dense ci_ci = ci.get_cim_2();
            mr::integer_dense ci_ri = ci.get_rim_2();

            const Integer* ptr_ri = ci_ri.ptr();
            const Integer* ptr_ci = ci_ci.ptr();

            for (j = 0; j < nc; ++j)
            {
                Integer c = ptr_ci[j]-1;
                Integer first_row = mat.first_row(c);
                Integer last_row = mat.last_row(c);

                Integer pos_c = imult(c,mat_ld);
                t_val_type* mat_ptr = mat.rep_ptr() + mat.first_elem_pos(c) - first_row + pos_c;

                for (i = 0; i <nr; ++i)
                {
                    Integer r = ptr_ri[i]-1;

                    if (r < first_row || r > last_row)
                        goto exit_label;

                    mrd::assign_helper(mat_ptr[r],val);
                };
            };
        }
        else if (ci.r_flag == 0 && ci.c_flag == 1)
        {
            Integer pos_c = imult(ci.c_start-1,mat.ld());
            Integer dpos_c = imult(ci.c_step,mat.ld());
            Integer c;

            mr::integer_dense ci_ri = ci.get_rim_2();
            const Integer* ptr_ri = ci_ri.ptr();

            for (j = 0, c = ci.c_start-1; j < nc; ++j, c += ci.c_step)
            {
                Integer first_row   = mat.first_row(c);
                Integer last_row    = mat.last_row(c);
                t_val_type* mat_ptr = mat.rep_ptr()+ pos_c + mat.first_elem_pos(c) - first_row;

                for (i = 0; i < nr; ++i)
                {
                    Integer r = ptr_ri[i] - 1;
                    if (r < first_row || r > last_row)
                        goto exit_label;

                    mrd::assign_helper(mat_ptr[r],val);
                };

                pos_c += dpos_c;
            };
        }
        else if (ci.r_flag == 1 && ci.c_flag == 0)
        {
            mr::integer_dense ci_ci = ci.get_cim_2();
            const Integer* ptr_ci = ci_ci.ptr();

            for (j = 0; j < nc; ++j)
            {
                Integer c = ptr_ci[j] - 1;
                Integer pos_c = imult(c,mat_ld);

                Integer first_row   = mat.first_row(c);
                Integer last_row    = mat.last_row(c);

                t_val_type* mat_ptr = mat.rep_ptr()+ pos_c + mat.first_elem_pos(c) - first_row;

                Integer r;

                for (i = 0, r = ci.r_start - 1; i < nr; ++i, r += ci.r_step)
                {
                    if (r < first_row || r > last_row)
                        goto exit_label;

                    mrd::assign_helper(mat_ptr[r],val);
                };
            };
        }
        else if (ci.r_flag == 1 && ci.c_flag == 1)
        {
            Integer c = ci.c_start - 1;
            Integer pos_c = imult(c,mat.ld());
            Integer dpos_c = imult(ci.c_step,mat.ld());

            for (j = 0; j < nc; ++j, c += ci.c_step)
            {
                Integer first_row   = mat.first_row(c);
                Integer last_row    = mat.last_row(c);
                t_val_type* mat_ptr = mat.rep_ptr()+ pos_c + mat.first_elem_pos(c) - first_row;

                Integer r;

                for (i = 0, r = ci.r_start - 1; i < nr; ++i, r += ci.r_step)
                {
                    if (r < first_row || r > last_row)
                        goto exit_label;

                    mrd::assign_helper(mat_ptr[r],val);
                };

                pos_c += dpos_c;
            };
        };

    exit_label:

        if ( i >= nr && j >= nc)
        {
            mat.get_struct().reset();
            return A;
        };

        using value_type    = typename M1::value_type;
        using matrix_type   = raw::Matrix<value_type,struct_dense>;
        A = Matrix(raw::converter<matrix_type,M1>::eval(A.get_impl<M1>()),false);
        matrix_type& B      = A.get_impl_unique<matrix_type>();
        value_type* ptr_B   = B.ptr();
        Integer B_ld        = B.ld();

        if (ci.r_flag == 0 && ci.c_flag == 0)
        {	
            mr::integer_dense ci_ci = ci.get_cim_2();
            mr::integer_dense ci_ri = ci.get_rim_2();

            const Integer* ptr_ri = ci_ri.ptr();
            const Integer* ptr_ci = ci_ci.ptr();

            for (; j < nc; ++j)
            {
                Integer pos_c = imult(ptr_ci[j]-1,B_ld);
                ptr_B = B.ptr() + pos_c;

                for (; i < nr; ++i)
                    mrd::assign_helper(ptr_B[ptr_ri[i] -1],val);

                i = 0;
            };
        }
        else if (ci.r_flag == 0 && ci.c_flag == 1)
        {
            Integer c = ci.c_start + imult(j,ci.c_step);
            Integer pos_c = imult(c-1,B.ld());
            Integer dpos_c = imult(ci.c_step,B.ld());

            mr::integer_dense ci_ri = ci.get_rim_2();
            const Integer* ptr_ri = ci_ri.ptr();

            ptr_B += pos_c;

            for (; j < nc; ++j)
            {		
                for (; i < nr; ++i)
                    mrd::assign_helper(ptr_B[ptr_ri[i] -1],val);

                i = 0;

                ptr_B += dpos_c;
            };
        }
        else if (ci.r_flag == 1 && ci.c_flag == 0)
        {	
            Integer pos_0 = ci.r_start - 1 + imult(i,ci.r_step);

            mr::integer_dense ci_ci = ci.get_cim_2();
            const Integer* ptr_ci = ci_ci.ptr();

            for (; j < nc; ++j)
            {
                Integer pos_c = imult(ptr_ci[j]-1,B_ld);
                Integer pos = pos_0;

                ptr_B = B.ptr() + pos_c;

                for (; i < nr; ++i, pos += ci.r_step)
                    mrd::assign_helper(ptr_B[pos],val);

                i = 0;
                pos_0 = ci.r_start - 1;
            };
        }
        else if (ci.r_flag == 1 && ci.c_flag == 1)
        {	
            Integer c = ci.c_start + imult(j,ci.c_step);
            Integer pos_0 = ci.r_start - 1 + imult(i,ci.r_step);

            Integer pos_s = imult(c-1,B.ld());	
            Integer dpos = imult(ci.c_step,B.ld());

            ptr_B = B.ptr() + pos_s;

            for (; j < nc; ++j)
            {								
                Integer pos = pos_0;
                for (; i < nr; ++i, pos += ci.r_step)
                    mrd::assign_helper(ptr_B[pos],val);

                i = 0;
                pos_0 = ci.r_start - 1;
                ptr_B += dpos;
            };
        };

        B.get_struct().reset();
        return A;
    };

    static bool set(M1& mat,Integer pos, const val_type& val)
    {
        Integer i,j;
        pos2ind(pos,mat.rows(),i,j);

        Integer first_row = mat.first_row(j);
        Integer last_row = mat.last_row(j);

        if ( i < first_row || i > last_row )
            return false;

        mrd::assign_helper(mat.rep_ptr()[mat.element_pos(i,j)], val);
        return true;
    };

    static bool set(M1& mat, Integer r1, Integer c1, const val_type& val)
    {
        Integer i           = r1 - 1;
        Integer j           = c1 - 1;

        Integer first_row = mat.first_row(j);
        Integer last_row = mat.last_row(j);

        if ( i < first_row || i > last_row )
            return false;

        mrd::assign_helper(mat.rep_ptr()[mat.element_pos(i,j)], val);
        return true;
    };

    static Matrix& eval_diag(Matrix& A,M1& mat,const val_type& val,Integer d)
    {
        A = Matrix(algorithm::band_change_diag(mat,d,val),true);
        return A;
    };
};

template<class M1,class val_type>
struct assign_functor_types
{
    using value_type    = typename M1::value_type;
    using struct_type   = typename M1::struct_type;

    using val_type_pr   = typename unify_types<value_type,val_type>::type;

    static const int code_1 = (int)details::type_to_code<value_type>::value;
    static const int code_2 = (int)details::type_to_code<val_type>::value;
    static const int code_ret = (int)details::type_to_code<val_type_pr>::value;
    static const bool req_promo = (code_1 != code_ret) || (code_ret == (int)mat_code::object_scalar);
};

template<class M1,class val_type>
struct assign_functor_scal_types
{
    using value_type = typename M1::value_type;
    static const bool req_promo = details::type_to_code<value_type>::value == mat_code::object_scalar;
};

template<class val_type>
struct assign_functor_scal_types<Object, val_type>
{
    static const bool req_promo = true;
};

template<class val_type>
struct assign_functor_scal_types<Real, val_type>
{
    static const bool req_promo = false;
};

template<class val_type>
struct assign_functor_scal_types<Float, val_type>
{
    static const bool req_promo = false;
};

template<class val_type>
struct assign_functor_scal_types<Integer, val_type>
{
    static const bool req_promo = false;
};

template<class val_type>
struct assign_functor_scal_types<Complex, val_type>
{
    static const bool req_promo = false;
};

template<class val_type>
struct assign_functor_scal_types<Float_complex, val_type>
{
    static const bool req_promo = false;
};

template<class M1,class val_type>
Matrix& assign_functor<M1,val_type>::eval(Matrix& A, M1& mat, const val_type& val,Integer i,Integer j)
{
    using AFT           = assign_functor_types<M1,val_type>;
    using val_type_pr   = typename AFT::val_type_pr;
    using struct_type   = typename AFT::struct_type;

    return assign_impl<M1,val_type_pr,AFT::req_promo,struct_type>
            ::eval(A,mat,raw::converter<val_type_pr,val_type>::eval(val),i,j);
};

template<class M1,class val_type>
Matrix& assign_functor<M1,val_type>::eval(Matrix& A,M1& mat, const val_type& val,Integer i)
{
    using AFT           = assign_functor_types<M1,val_type>;
    using val_type_pr   = typename AFT::val_type_pr;
    using struct_type   = typename AFT::struct_type;

    return assign_impl<M1,val_type_pr,AFT::req_promo, struct_type>
            ::eval(A,mat,raw::converter<val_type_pr,val_type>::eval(val),i);
};

template<class M1,class val_type>
Matrix& assign_functor<M1,val_type>::eval(Matrix& A, M1& mat,const val_type& val,const colon& c1)
{
    using AFT   = assign_functor_types<M1,val_type>;    

    bool tz     = md::has_trivial_assignment<M1,val_type>::eval(mat,val);
    (void)tz;

    bool is_dense = (A.get_struct_code() != struct_code::struct_dense);

    if (c1.m_flag == colon::t_all && is_dense == false)
    {
        using val_type_2    = typename md::real_type<val_type>::type;

        if (md::is_complex<val_type>::value == true 
            &&AFT::code_1 != (int)mat_code::object_scalar
            && tz && mrd::is_zero(mrd::imag_helper<val_type>::eval(val)))
        {
            using real_t    = typename md::real_type<val_type>::type;
            real_t val_r    = mrd::real_helper<val_type>::eval(val);
            raw::Matrix<real_t,struct_dense> tmp(get_matrix_ti<real_t>::eval(A),val_r,mat.rows(),mat.cols());
            A = Matrix(tmp,true);
            return A;
        }
        else
        {            
            using first_type    = typename AFT::value_type;
            using val_ret       = typename AFT::val_type_pr;
            using matrix_type   = raw::Matrix<val_ret,struct_dense>;

            ti::ti_type<first_type> ti_1    = get_matrix_ti<first_type>::eval(A);
            ti::ti_type<val_ret> ret_ti     = link_ret_ti<val_ret,first_type,val_type>
                                                ::eval(ti_1,val);

            val_ret val_init    = raw::converter<val_ret,val_type>::eval(val,ret_ti);
            matrix_type tmp(ret_ti,val_init,mat.rows(),mat.cols());
            A = Matrix(tmp,true);
            return A;
        };
    };    

    if (md::is_complex<val_type>::value == true
        && tz && mrd::is_zero(mrd::imag_helper<val_type>::eval(val)))
    {            
        using real_t    = typename md::real_type<val_type>::type;
        real_t val_r    = mrd::real_helper<val_type>::eval(val);

        return assign_functor<M1,real_t>::eval(A,A.get_impl_unique<M1>(),val_r,c1);
    }
    else
    {
        using val_ret       = typename AFT::val_type_pr;
        val_ret val_conv    = raw::converter<val_ret,val_type>::eval(val);
        return assign_impl<M1,val_ret,AFT::req_promo,typename AFT::struct_type>
                        ::eval(A,A.get_impl_unique<M1>(), val_conv,c1);
    };
}

template<class M1,class val_type>
Matrix& assign_functor<M1,val_type>::eval(Matrix& A, M1& mat,const val_type& val,const colon& c1,const colon& c2)
{
    if (c1.m_flag == colon::t_all && c2.m_flag == colon::t_all)
        return eval(A,mat,val,colon());

    using AFT   = assign_functor_types<M1,val_type>;

    bool tz     = md::has_trivial_assignment<M1,val_type>::eval(mat,val);
    (void)tz;

    if (md::is_complex<val_type>::value == true
        && tz && mrd::is_zero(mrd::imag_helper<val_type>::eval(val)))
    {            
        using real_t    = typename md::real_type<val_type>::type;
        real_t val_r    = mrd::real_helper<val_type>::eval(val);
        return assign_functor<M1,real_t>::eval(A,mat,val_r,c1,c2);
    }
    else
    {
        using val_ret       = typename AFT::val_type_pr;
        using impl_type     = assign_impl<M1,val_ret,AFT::req_promo,typename AFT::struct_type>;
        return impl_type::eval(A,mat,raw::converter<val_ret,val_type>::eval(val),c1,c2);
    };
};

template<class M1,class val_type>
Matrix& assign_functor<M1,val_type>::eval_diag(Matrix& A,M1& mat, const val_type& val,Integer d)
{
    using AFT   = assign_functor_types<M1,val_type>;

    return assign_impl<M1,typename AFT::val_type_pr,AFT::req_promo,typename AFT::struct_type>
                ::eval_diag(A,mat,raw::converter<typename AFT::val_type_pr,val_type>::eval(val),d);
};

template<class M1>
Matrix& assign_functor_drop<M1>::eval(Matrix& A, M1& mat, Integer i,Integer j, Real tol)
{
    using VT    = typename M1::value_type;
    using ST    = typename M1::struct_type;
    return assign_drop_impl<VT, ST>::eval(A,mat,i,j,tol);
}

template<class M1>
Matrix& assign_functor_drop<M1>::eval(Matrix& A,M1& mat, Integer i, Real tol)
{
    using VT    = typename M1::value_type;
    using ST    = typename M1::struct_type;
    return assign_drop_impl<VT, ST>::eval(A,mat,i,tol);
}

template<class M1>
Matrix& assign_functor_drop<M1>::eval(Matrix& A, M1& mat, const colon& c1, Real tol)
{
    using VT    = typename M1::value_type;
    using ST    = typename M1::struct_type;
    return assign_drop_impl<VT, ST>::eval(A,mat,c1,tol);
}

template<class M1>
Matrix& assign_functor_drop<M1>::eval(Matrix& A, M1& mat, const colon& c1,const colon& c2, Real tol)
{
    if (c1.m_flag == colon::t_all && c2.m_flag == colon::t_all)
        return eval(A,mat,colon(),tol);

    using VT    = typename M1::value_type;
    using ST    = typename M1::struct_type;
    return assign_drop_impl<VT, ST>::eval(A,mat,c1,c2,tol);
}

template<class M1>
Matrix& assign_functor_drop<M1>::eval_diag(Matrix& A,M1& mat, Integer d, Real tol)
{
    using VT    = typename M1::value_type;
    using ST    = typename M1::struct_type;
    return assign_drop_impl<VT, ST>::eval_diag(A,mat,d,tol);
}

template<class M1>
Matrix& assign_functor_add<M1>::eval(Matrix& A, M1& mat, Integer i,Integer j)
{
    using VT    = typename M1::value_type;
    using ST    = typename M1::struct_type;
    return assign_add_impl<VT, ST>::eval(A,mat,i,j);
}

template<class M1>
Matrix& assign_functor_add<M1>::eval(Matrix& A,M1& mat, Integer i)
{
    using VT    = typename M1::value_type;
    using ST    = typename M1::struct_type;
    return assign_add_impl<VT, ST>::eval(A,mat,i);
}

template<class M1>
Matrix& assign_functor_add<M1>::eval(Matrix& A, M1& mat, const colon& c1)
{
    using VT    = typename M1::value_type;
    using ST    = typename M1::struct_type;
    return assign_add_impl<VT, ST>::eval(A,mat,c1);
}

template<class M1>
Matrix& assign_functor_add<M1>::eval(Matrix& A, M1& mat, const colon& c1,const colon& c2)
{
    if (c1.m_flag == colon::t_all && c2.m_flag == colon::t_all)
        return eval(A,mat,colon());

    using VT    = typename M1::value_type;
    using ST    = typename M1::struct_type;
    return assign_add_impl<VT, ST>::eval(A,mat,c1,c2);
};

template<class M1>
Matrix& assign_functor_add<M1>::eval_diag(Matrix& A,M1& mat, Integer d)
{
    using VT    = typename M1::value_type;
    using ST    = typename M1::struct_type;
    return assign_add_impl<VT, ST>::eval_diag(A,mat,d);
}

template<class M1,class val_type, bool req_promo>
struct assign_scal_impl{};

template<class val_type>
struct assign_scal_impl<Object, val_type, true>
{
    using val_type_1    = Object;
    using val_type_ret  = typename unify_types<val_type,val_type_1>::type;

    static Matrix& eval(Matrix& A ,const val_type& val,Integer i,Integer j)
    {
        ti::ti_type<val_type_1> ti_1 = get_matrix_ti<val_type_1>::eval(A);
        ti::ti_type<val_type_ret> ret_ti = link_ret_ti<val_type_ret, Object, val_type>
                                            ::eval(ti_1,val);

        A = Matrix(raw::converter<val_type_ret, Object>::eval(A.get_scalar<Object>(), ret_ti), false);

        return assign_scal_impl<Object, val_type, false>::eval(A, val, i, j);
    };

    static Matrix& eval(Matrix& A, const val_type& val,Integer i)
    {	
        ti::ti_type<val_type_1> ti_1 = get_matrix_ti<val_type_1>::eval(A);
        ti::ti_type<val_type_ret> ret_ti = link_ret_ti<val_type_ret, Object, val_type>
                                            ::eval(ti_1,val);
        A = Matrix(raw::converter<val_type_ret, Object>::eval(A.get_scalar<Object>(), ret_ti), false);

        return assign_scal_impl<Object, val_type, false>::eval(A, val, i);
    }

    static Matrix& eval(Matrix& A, const val_type& val,const colon& c1)
    {
        ti::ti_type<val_type_1> ti_1 = get_matrix_ti<val_type_1>::eval(A);
        ti::ti_type<val_type_ret> ret_ti = link_ret_ti<val_type_ret, Object, val_type>
                                            ::eval(ti_1,val);

        A = Matrix(raw::converter<val_type_ret, Object>::eval(A.get_scalar<Object>(), ret_ti), false);

        return assign_scal_impl<Object, val_type, false>::eval(A, val, c1);
    }

    static Matrix& eval(Matrix& A, const val_type& val,const colon& c1,const colon& c2)
    {
        ti::ti_type<val_type_1> ti_1 = get_matrix_ti<val_type_1>::eval(A);
        ti::ti_type<val_type_ret> ret_ti = link_ret_ti<val_type_ret, Object, val_type>
                                            ::eval(ti_1,val);

        A = Matrix(raw::converter<val_type_ret, Object>::eval(A.get_scalar<Object>(), ret_ti), false);

        return assign_scal_impl<Object, val_type, false>::eval(A, val, c1, c2);
    };

};

template<class M1,class val_type>
struct assign_scal_impl<M1, val_type, false>
{
    static Matrix& eval(Matrix& A, const val_type& val,Integer i,Integer j)
    {
        error::check_index(i, j, 1, 1);
        A = val;
        return A;
    };

    static Matrix& eval(Matrix& A, const val_type& val,Integer i)
    {
        error::check_index(i, 1);
        A = val;
        return A;
    }

    static Matrix& eval(Matrix& A,const val_type& val,const colon& c1)
    {
        colon_info c_info;
        make_index(1,1,c1,c_info);

        Integer s = c_info.rows();

        if (s == 0)
        {
            return A;
        };
        A = val;
        return A;
    }

    static Matrix& eval(Matrix& A,const val_type& val,const colon& c1,const colon& c2)
    {
        colon_info ci;
        make_index(1,1,c1,c2,ci);

        if (ci.rows() == 0 || ci.cols() == 0)
        {
            return A;
        };		

        A = val;
        return A;
    };

};

template<class M1,class val_type>
Matrix& assign_functor_scal<M1,val_type>::eval(Matrix& A, const val_type& val,Integer i,Integer j)
{
    using AFST  = assign_functor_scal_types<M1, val_type>;
    return assign_scal_impl<M1, val_type, AFST::req_promo>::eval(A, val, i, j);
};

template<class M1,class val_type>
Matrix& assign_functor_scal<M1,val_type>::eval(Matrix& A,const val_type& val,Integer i)
{
    using AFST  = assign_functor_scal_types<M1, val_type>;
    return assign_scal_impl<M1, val_type, AFST::req_promo>::eval(A, val, i);
};

template<class M1,class val_type>
Matrix& assign_functor_scal<M1,val_type>::eval(Matrix& A,const val_type& val,const colon& c1)
{
    using AFST  = assign_functor_scal_types<M1, val_type>;
    return assign_scal_impl<M1, val_type, AFST::req_promo>::eval(A, val, c1);
};

template<class M1,class val_type>
Matrix& assign_functor_scal<M1,val_type>::eval(Matrix& A,const val_type& val,const colon& c1,const colon& c2)
{
    using AFST  = assign_functor_scal_types<M1, val_type>;
    return assign_scal_impl<M1, val_type, AFST::req_promo>::eval(A, val, c1, c2);
};

};};

#pragma warning( pop )