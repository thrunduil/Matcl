/*
 *  This file is a part of Matrix Computation Library (MATCL)
 *
 *  Copyright (c) Pawe³ Kowal 2017 - 2018
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

#include "matcl-matrep/visitors/insert_visitor.h"

namespace matcl { namespace details
{

template<class out_type,class out_str,class in_type,class struct_type,bool is_valid>
struct insert_visitor_impl{};

template<class out_type,class out_str,class in_type,bool is_valid>
struct insert_visitor_scal_impl{};

template<class out_type,class out_str,class in_type,class struct_type>
struct insert_visitor_impl<out_type,out_str,in_type,struct_type,false>
{
    static void eval(out_type& , Integer , Integer , const in_type& )
    {
        matcl_assert(0,"incompatible value types");
        throw;
    };
};

template<class out_type,class out_str,class in_type>
struct insert_visitor_scal_impl<out_type,out_str,in_type,false>
{
    static void eval(out_type& , Integer , Integer , const in_type& )
    {
        matcl_assert(0,"incompatible value types");
        throw;
    };
};

template<class out_type,class in_type>
struct insert_visitor_impl<out_type,struct_dense,in_type,struct_dense,true>
{
    static void eval(out_type& out, Integer row_st, Integer col_st, const in_type& mat)
    {
        Integer r = mat.rows(), c = mat.cols();
        using value_type    = typename out_type::value_type;
        using value_type_in = typename in_type::value_type;

        const value_type_in* ptr_mat = mat.ptr();
        value_type* ptr_out = out.ptr();
        Integer mat_ld      = mat.ld();
        Integer out_ld      = out.ld();

        ptr_out += col_st + row_st;
        for(Integer j = 0; j < c; ++j)
        {
            for(Integer i = 0; i < r; ++i)
            {
                value_type tmp = matcl::raw::converter_scalar<value_type,value_type_in>
                                    ::eval(ptr_mat[i],out.get_type());
                mrd::reset_helper(ptr_out[i],tmp);
            };
            ptr_mat += mat_ld;
            ptr_out += out_ld;
        };
    };;
};

template<class out_type,class in_type>
struct insert_visitor_scal_impl<out_type,struct_dense,in_type,true>
{
    static void eval(out_type& out, Integer row_st, Integer col_st, const in_type& mat)
    {
        col_st += row_st;
        using value_type = typename out_type::value_type;
        mrd::reset_helper(out.ptr()[col_st], matcl::raw::converter_scalar<value_type,in_type>::eval(mat));
    };
};

template<class out_type,class in_type>
struct insert_visitor_impl<out_type,struct_sparse,in_type,struct_dense,true>
{
    static void eval(out_type& out, Integer row_st, Integer col_st, const in_type& mat)
    {
        using value_type    = typename out_type::value_type;
        using value_type_in = typename in_type::value_type;

        if (mat.cols() == 0)
            return;

        raw::details::sparse_ccs<value_type>& d = out.rep();
        const value_type_in* ptr_mat = mat.ptr();

        Integer * d_c = d.ptr_c();
        Integer * d_r = d.ptr_r();
        value_type* d_x = d.ptr_x();

        Integer r       = mat.rows();
        Integer c       = mat.cols();		
        Integer mat_ld  = mat.ld();
        Integer k       = d_c[col_st+1];

        for(Integer j = 0; j < c; ++j)
        {			
            for(Integer i = 0, rr = row_st; i < r; ++i, ++rr)
            {
                value_type tmp = matcl::raw::converter_scalar<value_type,value_type_in>
                                    ::eval(ptr_mat[i],out.get_type());

                if (mrd::is_zero(tmp) == false)
                {					
                    d_r[k] = rr;
                    mrd::reset_helper(d_x[k],tmp);
                    ++k;
                };
            };
            ptr_mat += mat_ld;
            d_c[col_st + j+1] = k;
        };
    };;
};

template<class out_type,class in_type>
struct insert_visitor_scal_impl<out_type,struct_sparse,in_type,true>
{
    static void eval(out_type& out, Integer row_st, Integer col_st, const in_type& mat)
    {
        using value_type = typename out_type::value_type;

        raw::details::sparse_ccs<value_type>& d = out.rep();

        Integer * d_c   = d.ptr_c();
        Integer * d_r   = d.ptr_r();
        value_type* d_x = d.ptr_x();
        Integer k       = d_c[col_st+1];
        d_r[k]          = row_st;

        const value_type& conv  = matcl::raw::converter_scalar<value_type, in_type>
                                        ::eval(mat,out.get_type());

        mrd::reset_helper(d_x[k], conv);
        ++d_c[col_st + 1];
    };
};

template<class out_type,class in_type>
struct insert_visitor_impl<out_type,struct_dense,in_type,struct_banded,true>
{
    static void eval(out_type& out, Integer row_st, Integer col_st, const in_type& mat)
    {
        Integer r = mat.rows(), c = mat.cols();
        using value_type    = typename out_type::value_type;
        using value_type_in = typename in_type::value_type;

        if (mat.cols() == 0)
            return;

        value_type* ptr_out = out.ptr();
        ptr_out         += col_st+row_st;
        Integer mat_ld  = mat.ld();
        Integer out_ld  = out.ld();

        value_type* ptr_out_S = ptr_out;

        value_type Z = default_value<value_type>(out.get_type());
        
        for(Integer j = 0; j < c; ++j)
        {
            for(Integer i = 0; i < r; ++i)
                mrd::reset_helper(ptr_out[i],Z);

            ptr_out += out_ld;
        };

        ptr_out = ptr_out_S;
        
        if (mat.first_diag() == mat.last_diag())
        {
            Integer rc                      = mat.diag_length(mat.first_diag());
            const value_type_in* mat_ptr    = mat.rep_ptr() + mat.first_elem_diag(mat.first_diag());

            for (Integer j = 0; j < rc; ++j, ptr_out += out_ld+1)
            {
                value_type tmp = matcl::raw::converter_scalar<value_type,value_type_in>
                                    ::eval(mat_ptr[0],out.get_type());

                mrd::reset_helper(*ptr_out,tmp);
                mat_ptr += mat_ld;
            };
        }
        else
        {
            const value_type_in* mat_ptr = mat.rep_ptr();
            Integer matc        = mat.cols();

            for (Integer j = 0; j < matc; ++j)
            {
                Integer row_f   = mat.first_row(j);
                Integer row_l   = mat.last_row(j);
                Integer row_p   = mat.first_elem_pos(j);

                for (Integer i = row_f; i <= row_l; ++i, ++row_p)
                {
                    value_type tmp = matcl::raw::converter_scalar<value_type,value_type_in>
                                    ::eval(mat_ptr[row_p],out.get_type());
                    mrd::reset_helper(ptr_out[i],tmp);
                };		

                mat_ptr         += mat_ld;
                ptr_out         += out_ld;
            };
        };
    };
};

template<class out_type,class in_type>
struct insert_visitor_impl<out_type,struct_sparse,in_type,struct_banded,true>
{
    static void eval(out_type& out, Integer row_st, Integer col_st, const in_type& mat)
    {
        using value_type    = typename out_type::value_type;
        using value_type_in = typename in_type::value_type;

        if (mat.cols() == 0)
            return;

        raw::details::sparse_ccs<value_type>& d = out.rep();

        Integer * d_c   = d.ptr_c();
        Integer * d_r   = d.ptr_r();
        value_type* d_x = d.ptr_x();

        Integer c       = mat.cols();	
        Integer k       = d_c[col_st+1];
        Integer mat_ld  = mat.ld();

        if (mat.first_diag() == mat.last_diag())
        {
            Integer rc  = mat.diag_length(mat.first_diag());

            const value_type_in* mat_ptr = mat.rep_ptr() + mat.first_elem_diag(mat.first_diag());

            for(Integer j = 1, row = row_st; j <= rc; ++j, ++row)
            {			
                value_type tmp = matcl::raw::converter_scalar<value_type,value_type_in>
                                    ::eval(*mat_ptr,out.get_type());

                d_r[k] = row;
                mrd::reset_helper(d_x[k],tmp);
                ++k;

                d_c[col_st + j] = k;
                mat_ptr += mat_ld;
            };

            for (Integer j = rc + 1; j <= c; ++j)
                d_c[col_st + j] = k;
        }
        else
        {
            const value_type_in* mat_ptr = mat.rep_ptr();

            for(Integer j = 0; j < c; ++j)
            {			
                Integer row_f = mat.first_row(j);
                Integer row_l = mat.last_row(j);
                Integer ii    = mat.first_elem_pos(j);

                for(Integer i = row_f, rr = row_st + row_f; i <= row_l; ++i, ++ii, ++rr)
                {
                    value_type tmp = matcl::raw::converter_scalar<value_type,value_type_in>
                                        ::eval(mat_ptr[ii],out.get_type());
                    d_r[k] = rr;
                    mrd::reset_helper(d_x[k],tmp);
                    ++k;
                };

                mat_ptr  += mat_ld;
                d_c[col_st + j + 1] = k;
            };
        };
    };
};

template<class out_type,class in_type>
struct insert_visitor_impl<out_type,struct_dense,in_type,struct_sparse,true>
{
    static void eval(out_type& out, Integer row_st, Integer col_st, const in_type& mat)
    {
        Integer r = mat.rows(), c = mat.cols();

        if (c == 0)
            return;

        using val_type_out      = typename out_type::value_type;
        val_type_out Z          = default_value<val_type_out>(out.get_type());
        val_type_out* ptr_out   = out.ptr();
        Integer out_ld          = out.ld();

        ptr_out += col_st + row_st;

        val_type_out* ptr_out_s = ptr_out;

        for(Integer j = 0; j < c; ++j)
        {
            for(Integer i = 0; i < r; ++i)
                mrd::reset_helper(ptr_out[i],Z);

            ptr_out+= out_ld;
        };

        ptr_out = ptr_out_s;

        using value_type    = typename in_type::value_type;

        if (mat.nnz() > 0)
        {
            const raw::details::sparse_ccs<value_type>& Ad = mat.rep();
            const Integer* Ad_c		= Ad.ptr_c();
            const Integer* Ad_r		= Ad.ptr_r();
            const value_type* Ad_x	= Ad.ptr_x();

            for (Integer j = 0; j < c; ++j)
            {
                Integer k;
                for (k = Ad_c[j]; k < Ad_c[j + 1] ; ++k)
                {
                    val_type_out tmp = matcl::raw::converter_scalar<val_type_out,value_type>
                                        ::eval(Ad_x[k],out.get_type());
                    mrd::reset_helper(ptr_out[Ad_r[k]],tmp);
                };

                ptr_out += out_ld;
            };
        };

        return;
    };;
};

template<class out_type,class in_type>
struct insert_visitor_impl<out_type,struct_sparse,in_type,struct_sparse,true>
{
    static void eval(out_type& out, Integer row_st, Integer col_st, const in_type& mat)
    {
        using value_type    = typename out_type::value_type;
        using value_type_in = typename in_type::value_type;

        if (mat.cols() == 0)
            return;

        raw::details::sparse_ccs<value_type>& d = out.rep();

        Integer * d_c = d.ptr_c();
        Integer * d_r = d.ptr_r();
        value_type* d_x = d.ptr_x();

        const raw::details::sparse_ccs<value_type_in>& Ad = mat.rep();
        const Integer* Ad_c		= Ad.ptr_c();
        const Integer* Ad_r		= Ad.ptr_r();
        const value_type_in* Ad_x	= Ad.ptr_x();

        Integer c = mat.cols();		
        Integer p = d_c[col_st+1];

        if (mat.nnz() == 0)
        {
            for(Integer j = 0; j < c; ++j)
                d_c[col_st + 1 + j] = p;

            return;
        };

        for(Integer j = 0; j < c; ++j)
        {		
            for (Integer k = Ad_c[j]; k < Ad_c[j + 1] ; ++k)
            {
                value_type tmp = matcl::raw::converter_scalar<value_type,value_type_in>
                                    ::eval(Ad_x[k],out.get_type());

                d_r[p]  = row_st + Ad_r[k];
                mrd::reset_helper(d_x[p],tmp);
                ++p;
            };

            d_c[col_st + 1 + j] = p;
        };
    };
};

template<class Matrix_type, class Const_matrix_type = mr::const_matrix<Matrix_type>>
void insert_sparse_cols(Matrix_type& out, const std::list<Const_matrix_type>& mat_vec)
{
    matcl_assert(out.rep().offset() == 0,"");

    using value_type    = typename Matrix_type::value_type;

    raw::details::sparse_ccs<value_type>& d = out.rep();

    Integer * d_c   = d.ptr_c();
    Integer * d_r   = d.ptr_r();
    value_type* d_x = d.ptr_x();

    using iterator  = typename std::list<Const_matrix_type>::const_iterator;
    Integer c = out.cols(), nz = 0;		

    for(Integer j = 0; j < c; ++j)
    {		
        d_c[j] = nz;

        Integer row = 0;
        for (iterator I = mat_vec.begin(), E = mat_vec.end(); I!= E; ++I)
        {
            const raw::details::sparse_ccs<value_type>& Ad = I->get().rep();
            
            if (Ad.nnz() == 0)
            {
                row += I->get().rows();
                continue;
            };

            const Integer* Ad_c		= Ad.ptr_c();
            const Integer* Ad_r		= Ad.ptr_r();
            const value_type* Ad_x	= Ad.ptr_x();

            for (Integer k = Ad_c[j]; k < Ad_c[j + 1] ; ++k)
            {
                const value_type& tmp = Ad_x[k];
                d_r[nz] = row + Ad_r[k];
                mrd::reset_helper(d_x[nz],tmp);
                ++nz;
            };

            row += I->get().rows();
        };		
    };

    d_c[c] = nz;
    d.add_memory(-1);
};

};};