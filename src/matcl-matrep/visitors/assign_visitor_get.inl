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

#include "matcl-core/error/exception_classes.h"
#include "matcl-scalar/details/scalfunc_helpers.h"

namespace matcl { namespace details
{

namespace mr = matcl::raw;

template<class M1,class struct_type>
struct get_impl
{};

template<class M1>
struct get_impl<M1,struct_dense>
{
    using value_type = typename M1::value_type;

    static void eval(const M1& mat,Integer i,Integer j, Matrix& ret)
    {
        ret = mat(i,j);
    };

    static void eval(const M1& mat,Integer i, Matrix& ret)
    {
        Integer r = mat.rows();
        Integer j;
        pos2ind(i,r,i,j);
        ret = mat(i+1,j+1);
    };

    static void eval(const Matrix& A,const M1& mat,const colon& c1, Matrix& ret)
    {
        colon_info c_info;
        make_index(mat.rows(),mat.cols(), c1, c_info);

        Integer s = c_info.rows();

        if (s == 1)
        {
            eval(mat,c_info.r_start,ret);
            return;
        };

        if (c_info.r_step == 1)
        {
            M1 tmp = mat.make_view(c_info.r_start, c_info.r_end);
            ret = Matrix(tmp.reshape(c_info.rep_rows(), c_info.rep_cols()),false);
            return;
        };

        using matrix_type = raw::Matrix<value_type,struct_dense>;
        const typename M1::value_type* ptr_mat = mat.ptr();
        Integer mat_ld  = mat.ld();
        Integer r       = mat.rows();

        ti::ti_type<value_type> ti_A = get_matrix_ti<value_type>::eval(A);

        matrix_type out(ti_A);		

        if (c_info.r_flag == 0)
        {            
            out.reset_unique(c_info.rep_rows(),c_info.rep_cols());

            value_type* ptr_out = out.ptr();

            bool single             = c_info.is_double_mat_colon() == false;

            if (single == true)
            {
                const mr::integer_dense& ci_ri  = c_info.get_rim_1();
                const Integer* ptr_ri           = ci_ri.ptr();

                for (Integer i = 0; i < s; ++i)
                {
                    Integer pos = ptr_ri[i];
                    Integer r2, c2;
                    pos2ind(pos,r,r2,c2);

                    mrd::reset_helper(ptr_out[i],ptr_mat[r2+c2*mat_ld]);
                };
            }
            else
            {
                const mr::integer_dense& rim   = c_info.get_rim_r();
                const mr::integer_dense& cim   = c_info.get_rim_c();

                const Integer* ptr_ri   = rim.ptr();
                const Integer* ptr_ci   = cim.ptr();

                for (Integer i = 0; i < s; ++i)
                {
                    Integer r2          = ptr_ri[i] - 1;
                    Integer c2          = ptr_ci[i] - 1;

                    mrd::reset_helper(ptr_out[i], ptr_mat[r2 + c2 * mat_ld]);
                };
            };
        }
        else
        {
            out.reset_unique(c_info.rep_rows(),c_info.rep_cols());

            if (s == 0)
            {
                ret = Matrix(std::move(out),true);
                return;
            };

            value_type* ptr_out = out.ptr();

            if (mat.ld() == mat.rows())
            {
                for (Integer i = 0, pos = c_info.r_start-1; i < s; ++i, pos += c_info.r_step)
                {
                    mrd::reset_helper(ptr_out[i],ptr_mat[pos]);
                };
            }
            else
            {
                Integer pos = c_info.r_start-1;
                ptr_mat     = mat.ptr() + (pos/r)*mat.ld();
                pos         = pos % r;

                if (c_info.r_step > 0)
                {
                    for (Integer i = 0; i < s; ++i, pos += c_info.r_step)
                    {
                        ptr_mat     += (pos/r)*mat_ld;
                        pos         = pos%r;
                        mrd::reset_helper(ptr_out[i],ptr_mat[pos]);
                    };
                }
                else
                {
                    for (Integer i = 0; i < s; ++i, pos += c_info.r_step)
                    {
                        ptr_mat     += (pos/r)*mat_ld;
                        pos         = pos%r;
                        if (pos < 0)
                        {
                            pos     += r;
                            ptr_mat -= mat_ld;
                        };
                        mrd::reset_helper(ptr_out[i],ptr_mat[pos]);
                    };
                };
            };
        };
        
        ret = Matrix(out,true);
    };

    static void eval(const Matrix& A,const M1& mat,const colon& c1,const colon& c2,  Matrix& ret)
    {
        colon_info c_info;
        make_index(mat.rows(),mat.cols(), c1, c2, c_info);

        using matrix_type = raw::Matrix<value_type,struct_dense>;

        const Integer* ptr_ri   = c_info.get_ri_2_ptr();
        const Integer* ptr_ci   = c_info.get_ci_2_ptr();
        Integer mat_ld          = mat.ld();

        const typename M1::value_type* ptr_mat = mat.ptr();

        if (c_info.rows() == 1 && c_info.cols() == 1)
        {
            eval(mat,c_info.r_start,c_info.c_start,ret);
            return;
        };

        ti::ti_type<value_type> ti_A = get_matrix_ti<value_type>::eval(A);

        if (c_info.r_flag == 0 && c_info.c_flag == 0)
        {
            Integer r = c_info.rows(), c = c_info.cols();

            error::check_size(r,c);

            matrix_type out(ti_A,r,c);
            value_type* ptr_out = out.ptr();
            Integer out_ld      = out.ld();

            if (r == 0 || c == 0)
            {
                ret = Matrix(std::move(out),false);
                return;
            };

            for (Integer i = 0; i < c; ++i)
            {
                Integer pos_mat = imult(ptr_ci[i]-1,mat_ld);
                ptr_mat = mat.ptr() + pos_mat;

                for (Integer j = 0; j < r; ++j)
                {
                    mrd::reset_helper(ptr_out[j],ptr_mat[ptr_ri[j]-1]);
                };
                ptr_out += out_ld;
            };		
            ret = Matrix(std::move(out),true);
            return;
        }
        else if (c_info.r_flag == 0 && c_info.c_flag == 1)
        {
            Integer r = c_info.rows();

            matrix_type out(ti_A,r,c_info.c_size);
            value_type* ptr_out = out.ptr();
            Integer out_ld      = out.ld();

            if (r == 0 || c_info.c_size == 0)
            {
                ret = Matrix(std::move(out),false);
                return;
            };

            Integer pos_mat = imult(c_info.c_start-1,mat.ld());
            Integer dpos = imult(c_info.c_step,mat.ld());

            ptr_mat = mat.ptr() + pos_mat;

            for (Integer i = 0; i < c_info.c_size; ++i)
            {
                for (Integer j = 0; j < r; ++j)
                {
                    mrd::reset_helper(ptr_out[j],ptr_mat[ptr_ri[j]-1]);
                };
                ptr_out += out_ld;
                ptr_mat += dpos;
            };		
            ret = Matrix(std::move(out),true);
            return;
        }
        else if (c_info.r_flag == 1 && c_info.c_flag == 0)
        {
            Integer c = c_info.cols();

            matrix_type out(ti_A,c_info.r_size,c);
            value_type* ptr_out = out.ptr();
            Integer out_ld      = out.ld();

            if (c_info.r_size == 0 || c == 0)
            {
                ret = Matrix(std::move(out),false);
                return;
            };
            Integer r_start = c_info.r_start;
            Integer r_step = c_info.r_step;
            Integer r_size = c_info.r_size;

            for (Integer i = 0; i < c; ++i)
            {
                Integer pos_mat = imult(ptr_ci[i]-1,mat_ld);
                ptr_mat = mat.ptr() + pos_mat;
                for (Integer j = 0, jj = r_start-1; j < r_size; ++j, jj+= r_step)
                {
                    mrd::reset_helper(ptr_out[j],ptr_mat[jj]);
                };
                ptr_out += out_ld;
            };		
            ret = Matrix(std::move(out),true);
            return;
        }
        else if (c_info.r_flag == 1 && c_info.c_flag == 1)			
        {
            if (c_info.r_step == 1 && c_info.c_step == 1)
            {
                ret = Matrix(mat.make_view(c_info.r_start, c_info.r_end, c_info.c_start, 
                                c_info.c_end),false);
                return;
            };

            matrix_type out(ti_A,c_info.r_size,c_info.c_size);
            value_type* ptr_out = out.ptr();
            Integer out_ld      = out.ld();

            if (c_info.r_size == 0 || c_info.c_size == 0)
            {
                ret = Matrix(std::move(out),false);
                return;
            };
            
            Integer r_start = c_info.r_start;
            Integer r_step = c_info.r_step;
            Integer r_size = c_info.r_size;

            Integer pos_mat = imult(c_info.c_start-1,mat.ld());
            Integer dpos    = imult(c_info.c_step,mat.ld());

            ptr_mat = mat.ptr() + pos_mat;
            for (Integer i = 0; i < c_info.c_size; ++i)
            {				
                for (Integer j = 0, jj = r_start-1; j < r_size; ++j, jj+= r_step)
                {
                    mrd::reset_helper(ptr_out[j],ptr_mat[jj]);
                };
                ptr_out += out_ld;
                ptr_mat += dpos;
            };		
            ret = Matrix(std::move(out),true);
            return;
        }
        matcl_assert(0,"unknown case");
        throw;
    };
};

template<class M1>
struct get_impl<M1,struct_sparse>
{
    using value_type = typename M1::value_type;

    static void eval(const M1& mat,Integer i,Integer j, Matrix& ret)
    {
        ret = mat(i,j);
    };

    static void eval(const M1& mat,Integer i, Matrix& ret)
    {
        Integer r = mat.rows();
        Integer j;
        pos2ind(i,r,i,j);
        eval(mat,i+1,j+1,ret);
    };

    static void eval(const Matrix& ,const M1& mat,const colon& c1, Matrix& ret)
    {
        colon_info c_info;
        make_index(mat.rows(),mat.cols(), c1, c_info);

        if (c_info.rows() == 1)
        {
            eval(mat,c_info.r_start,ret);
            return;
        };

        algorithm::get_submatrix_2(ret, mat,c_info);
    };

    static void eval(const Matrix& ,const M1& mat,const colon& c1,const colon& c2, Matrix& ret)
    {
        colon_info c_info;
        make_index(mat.rows(),mat.cols(), c1, c2, c_info);

        if (c_info.rows() == 1 && c_info.cols() == 1)
        {
            eval(mat,c_info.r_start,c_info.c_start,ret);
            return;
        };

        if (c_info.r_step == 1 && c_info.c_step == 1 && c_info.r_start == 1 && c_info.r_end == mat.rows())
        {
            ret = Matrix((M1)(mat.make_view(c_info.c_start, c_info.c_end)),false);
            return;
        };

        algorithm::get_submatrix(ret, mat, c_info);
        return;
    };
};

template<class M1>
struct get_impl<M1,struct_banded>
{
    using value_type = typename M1::value_type;

    static void eval(const M1& mat,Integer i, Matrix& ret)
    {
        Integer r = mat.rows();
        Integer j;
        pos2ind(i,r,i,j);
        eval(mat,i+1,j+1,ret);
    };

    static void eval(const M1& mat,Integer i,Integer j, Matrix& ret)
    {
        Integer r = mat.rows();
        Integer c = mat.cols();

        --i;
        --j;

        if (i < 0 || j < 0 || i >= r || j >= c)
            throw error::invalid_index(i+1, j+1, r, c);

        if ( i < mat.first_row(j) || i > mat.last_row(j) )
        {
            ret = md::default_value<value_type>(mat.get_type());
            return;
        }
        ret = mat.rep_ptr()[mat.element_pos(i,j)];
    };

    static void eval(const Matrix& A, const M1& mat,const colon& c1, Matrix& ret)
    {
        colon_info c_info;
        make_index(mat.rows(),mat.cols(), c1, c_info);

        Integer s = c_info.rows();		
        using matrix_type = raw::Matrix<value_type,struct_sparse>;
        
        if (s == 1)
            eval(mat,c_info.r_start,ret);

        ti::ti_type<value_type> ti_A = get_matrix_ti<value_type>::eval(A);
        

        if (c_info.r_flag == 0)
        {
            bool single             = c_info.is_double_mat_colon() == false;

            Integer mr = c_info.rep_rows();
            Integer mc = c_info.rep_cols();            

            if (mat.rows() == 0 || mat.cols() == 0 || s == 0)
            {
                ret = Matrix(matrix_type(ti_A,mr,mc),true);
                return;
            };

            Integer nz_est = icast(Real(mat.nnz()*s)/Real(mat.rows())/Real(mat.cols()));

            matrix_type out(ti_A,mr,mc,nz_est);

            Integer nz      = 0;
            Integer nz_max  = nz_est;

            raw::details::sparse_ccs<value_type>& d = out.rep();
            Integer* d_c		= d.ptr_c();
            Integer* d_r		= d.ptr_r();
            value_type * d_x	= d.ptr_x();
            Integer matr        = mat.rows();

            if (single == true)
            {
                const mr::integer_dense& ci_ri  = c_info.get_rim_1();
                const Integer* ptr_ri           = ci_ri.ptr();

                for (Integer j = 1, p = 1; j <= mc ; ++j)
                {
                    for (Integer i = 1; i <= mr; ++i, ++p)
                    {
                        Integer pos = ptr_ri[p-1];
                        Integer r, c;
                        pos2ind(pos,matr,r,c);
                    
                        if ( r >= mat.first_row(c) && r <= mat.last_row(c) )
                        {
                            const value_type& val = mat.rep_ptr()[mat.element_pos(r,c)];

                            if (nz >= nz_max)
                            {
                                d.add_memory( d.nzmax() + 1);
                                d_r		= d.ptr_r();
                                d_x		= d.ptr_x();
                                nz_max	= d.nzmax();
                            };

                            d_r[nz]	= i-1;
                            mrd::reset_helper(d_x[nz],val);
                            ++d_c[j];
                            ++nz;
                        };	
                    };
                };
            }
            else
            {
                const mr::integer_dense& rim   = c_info.get_rim_r();
                const mr::integer_dense& cim   = c_info.get_rim_c();

                const Integer* ptr_ri   = rim.ptr();
                const Integer* ptr_ci   = cim.ptr();

                for (Integer j = 1, p = 1; j <= mc ; ++j)
                {
                    for (Integer i = 1; i <= mr; ++i, ++p)
                    {
                        Integer r       = ptr_ri[p-1] - 1;
                        Integer c       = ptr_ci[p-1] - 1;
                    
                        if ( r >= mat.first_row(c) && r <= mat.last_row(c) )
                        {
                            const value_type& val = mat.rep_ptr()[mat.element_pos(r,c)];

                            if (nz >= nz_max)
                            {
                                d.add_memory( d.nzmax() + 1);
                                d_r		= d.ptr_r();
                                d_x		= d.ptr_x();
                                nz_max	= d.nzmax();
                            };

                            d_r[nz]	= i-1;
                            mrd::reset_helper(d_x[nz],val);
                            ++d_c[j];
                            ++nz;
                        };	
                    };
                };
            };

            for (Integer i = 1; i <= mc; ++i)
                d_c[i] += d_c[i-1];
            
            ret = Matrix(std::move(out),true);
            return;
        }
        else
        {
            if (mat.rows() == 0 || mat.cols() == 0 || s == 0)
            {
                matrix_type out(ti_A,c_info.rep_rows(),c_info.rep_cols());
                out.get_struct().set(predefined_struct_type::diag);
                ret = Matrix(std::move(out),true);
                return;
            };

            Integer nz_est = icast(Real(mat.nnz()*s)/Real(mat.rows())/Real(mat.cols()));
            matrix_type out(ti_A,s,1,nz_est);
            Integer nz = 0;
            Integer nz_max = nz_est;

            raw::details::sparse_ccs<value_type>& d = out.rep();
            Integer* d_c		= d.ptr_c();
            Integer* d_r		= d.ptr_r();
            value_type * d_x	= d.ptr_x();
            Integer matr        = mat.rows();

            for (Integer i = 1, pos = c_info.r_start; i <= s; ++i, pos += c_info.r_step)
            {
                Integer r, c;
                pos2ind(pos,matr,r,c);

                if ( r >= mat.first_row(c) && r <= mat.last_row(c) )
                {
                    const value_type& val =  mat.rep_ptr()[mat.element_pos(r,c)];

                    if (nz >= nz_max)
                    {
                        d.add_memory( d.nzmax() + 1);
                        d_r		= d.ptr_r();
                        d_x		= d.ptr_x();
                        nz_max	= d.nzmax();
                    };

                    d_r[nz]	= i-1;
                    mrd::reset_helper(d_x[nz],val);
                    ++d_c[1];
                    ++nz;
                };
            };

            if (out.rows() == c_info.rep_rows())
            {
                ret = Matrix(std::move(out),true);
                return;
            }
            else
            {
                 mrd::manip_reshape_helper<matrix_type>
                     ::eval_reshape(ret, out,c_info.rep_rows(), c_info.rep_cols());
                 return;
            };
        };
    };

    static void eval(const Matrix& A, const M1& mat,const colon& c1,const colon& c2, Matrix& ret)
    {
        colon_info c_info;
        make_index(mat.rows(),mat.cols(), c1, c2, c_info);

        if (c_info.rows() == 1 && c_info.cols() == 1)
        {
            eval(mat,c_info.r_start,c_info.c_start,ret);
            return;
        };

        using matrix_type = raw::Matrix<value_type,struct_sparse>;

        ti::ti_type<value_type> ti_A = get_matrix_ti<value_type>::eval(A);
        Integer mat_ld  = mat.ld();

        if (c_info.r_flag == 0 && c_info.c_flag == 0)
        {
            const mr::integer_dense& ci_ri = c_info.get_rim_2();
            const mr::integer_dense& ci_ci = c_info.get_cim_2();
            const Integer* ptr_ri = ci_ri.ptr();
            const Integer* ptr_ci = ci_ci.ptr();

            Integer mr = ci_ri.size(), mc = ci_ci.size();
            error::check_size(mr,mc);

            if (mr == 0 || mc == 0)
            {
                ret = Matrix(matrix_type(ti_A,mr,mc),false);
                return;
            };

            Integer nz_est = Integer(mat.nnz()*Real(mr)/Real(mat.rows())*Real(mc)/Real(mat.cols()));

            matrix_type out(ti_A,mr,mc,nz_est);
            Integer nz = 0;
            Integer nz_max = nz_est;

            raw::details::sparse_ccs<value_type>& d = out.rep();
            Integer* d_c		= d.ptr_c();
            Integer* d_r		= d.ptr_r();
            value_type * d_x	= d.ptr_x();

            for (Integer i = 0; i < mc; ++i)
            {
                Integer c           = ptr_ci[i] - 1;
                Integer first_row   = mat.first_row(c);
                Integer last_row    = mat.last_row(c);

                Integer pos_B       = mat.first_elem_pos(c) - first_row + imult(c,mat_ld);

                const value_type* ptr_B = mat.rep_ptr() + pos_B;

                for (Integer j = 0; j < mr; ++j)
                {
                    Integer r = ptr_ri[j] - 1;			

                    if ( r >= first_row && r <= last_row )
                    {
                        const value_type& val =  ptr_B[r];

                        if (nz >= nz_max)
                        {
                            d.add_memory( d.nzmax() + 1);
                            d_r		= d.ptr_r();
                            d_x		= d.ptr_x();						
                            nz_max	= d.nzmax();
                        };

                        d_r[nz]	= j;
                        mrd::reset_helper(d_x[nz],val);
                        ++d_c[i+1];
                        ++nz;
                    };
                };
            };		

            for (Integer i = 0; i < mc; ++i)
                d_c[i+1] += d_c[i];

            ret = Matrix(std::move(out),true);
            return;
        }
        else if (c_info.r_flag == 0 && c_info.c_flag == 1)
        {
            const mr::integer_dense& ci_ri  = c_info.get_rim_2();
            const Integer* ptr_ri           = ci_ri.ptr();

            Integer mr = ci_ri.size(), mc = c_info.c_size;
            
            if (mr == 0 || mc == 0)
            {
                ret = Matrix(matrix_type(ti_A,mr,mc),false);
                return;
            };

            Integer nz_est = Integer(mat.nnz()*Real(mr)/Real(mat.rows())*Real(mc)/Real(mat.cols()));

            matrix_type out(ti_A,mr,mc,nz_est);
            Integer nz = 0;
            Integer nz_max = nz_est;

            raw::details::sparse_ccs<value_type>& d = out.rep();
            Integer* d_c		= d.ptr_c();
            Integer* d_r		= d.ptr_r();
            value_type * d_x	= d.ptr_x();

            Integer pos_c       = imult(c_info.c_start-1,mat.ld());
            Integer dpos_c      = imult(c_info.c_step,mat.ld());

            for (Integer i = 0, c = c_info.c_start-1; i < c_info.c_size; ++i, c+= c_info.c_step)
            {
                Integer first_row   = mat.first_row(c);
                Integer last_row    = mat.last_row(c);

                Integer pos_B       = mat.first_elem_pos(c) - first_row + pos_c;
                const value_type* ptr_B = mat.rep_ptr() + pos_B;

                for (Integer j = 0; j < mr; ++j)
                {
                    Integer r = ptr_ri[j]-1;

                    if ( r >= first_row && r <= last_row )
                    {
                        const value_type& val =  ptr_B[r];

                        if (nz >= nz_max)
                        {
                            d.add_memory( d.nzmax() + 1);
                            d_r		= d.ptr_r();
                            d_x		= d.ptr_x();						
                            nz_max	= d.nzmax();
                        };
                        d_r[nz]	= j;
                        mrd::reset_helper(d_x[nz],val);
                        ++d_c[i+1];
                        ++nz;
                    };
                };

                pos_c += dpos_c;
            };		

            for (Integer i = 0; i < mc; ++i)
                d_c[i+1] += d_c[i];

            ret = Matrix(std::move(out),true);
            return;
        }
        else if (c_info.r_flag == 1 && c_info.c_flag == 0)
        {
            const mr::integer_dense& ci_ci  = c_info.get_cim_2();
            const Integer* ptr_ci           = ci_ci.ptr();

            Integer mr = c_info.r_size, mc = ci_ci.size();
            
            if (mr == 0 || mc == 0)
            {
                ret = Matrix(matrix_type(ti_A,mr,mc),false);
                return;
            };

            Integer nz_est = Integer(mat.nnz()*Real(mr)/Real(mat.rows())*Real(mc)/Real(mat.cols()));

            matrix_type out(ti_A,mr,mc,nz_est);
            Integer nz = 0;
            Integer nz_max = nz_est;

            raw::details::sparse_ccs<value_type>& d = out.rep();
            Integer* d_c		= d.ptr_c();
            Integer* d_r		= d.ptr_r();
            value_type * d_x	= d.ptr_x();

            for (Integer i = 0; i < mc; ++i)
            {
                d_c[i]          = nz;
                Integer c       = ptr_ci[i]-1;
                Integer fr      = mat.first_row(c);
                Integer lr      = mat.last_row(c);

                if (lr < fr)
                    continue;

                Integer pos_c   = imult(c,mat_ld);
                Integer pos_B   = pos_c +  mat.first_elem_pos(c) - fr;
                const value_type* ptr_B = mat.rep_ptr() + pos_B;

                for (Integer j = 0, r = c_info.r_start - 1; j < c_info.r_size; ++j, r+= c_info.r_step)
                {
                    if ( r >= fr && r <= lr )
                    {
                        const value_type& val =  ptr_B[r];

                        if (nz >= nz_max)
                        {
                            d.add_memory( d.nzmax() + 1);
                            d_r		= d.ptr_r();
                            d_x		= d.ptr_x();				
                            nz_max	= d.nzmax();
                        };
                        d_r[nz]	    = j;
                        mrd::reset_helper(d_x[nz],val);						
                        ++nz;
                    };
                };
            };		

            d_c[mc]                 = nz;
            ret = Matrix(std::move(out),true);
            return;
        }
        else if (c_info.r_flag == 1 && c_info.c_flag == 1)			
        {
            Integer mr = c_info.r_size, mc = c_info.c_size;

            if (c_info.r_step == 1 && c_info.c_step == 1 && (c_info.c_start == c_info.r_start))
            {
                ret = Matrix(mat.make_view(c_info.c_start, c_info.r_end, c_info.c_end),false);
                return;
            };

            if (mr == 0 || mc == 0)
            {
                ret = Matrix(matrix_type(ti_A,mr,mc),false);
                return;
            };

            Integer nz_est = Integer(mat.nnz()*Real(mr)/Real(mat.rows())*Real(mc)/Real(mat.cols()));

            matrix_type out(ti_A,mr,mc,nz_est);
            Integer nz = 0;
            Integer nz_max = nz_est;

            raw::details::sparse_ccs<value_type>& d = out.rep();
            Integer* d_c		= d.ptr_c();
            Integer* d_r		= d.ptr_r();
            value_type * d_x	= d.ptr_x();

            Integer pos_c       = imult(c_info.c_start-1,mat.ld());
            Integer dpos_c      = imult(c_info.c_step,mat.ld());

            for (Integer i = 0, c = c_info.c_start-1; i < c_info.c_size; ++i, c+= c_info.c_step)
            {				
                Integer fr      = mat.first_row(c);
                Integer lr      = mat.last_row(c);
                Integer pos_B   = pos_c + mat.first_elem_pos(c) - fr;
                const value_type* ptr_B = mat.rep_ptr() + pos_B;

                for (Integer j = 0, r = c_info.r_start - 1; j < c_info.r_size; ++j, r+= c_info.r_step)
                {
                    if ( r >= fr && r <= lr )
                    {
                        const value_type& val =  ptr_B[r];

                        if (nz >= nz_max)
                        {
                            d.add_memory( d.nzmax() + 1);
                            d_r		= d.ptr_r();
                            d_x		= d.ptr_x();				
                            nz_max	= d.nzmax();
                        };

                        d_r[nz]	= j;
                        mrd::reset_helper(d_x[nz],val);
                        ++d_c[i+1];
                        ++nz;
                    };		
                };
                pos_c += dpos_c;
            };

            for (Integer i = 0; i < mc; ++i)
                d_c[i+1] += d_c[i];

            ret = Matrix(std::move(out),true);
            return;
        }

        matcl_assert(0,"unknown case");
        throw;
    };

    static void make_view_2(const M1& mat, const colon_info& c, Matrix& ret)
    {
        ret = mat.make_view(c.r_end-1, c.c_end-1);
        return;
    };
};

template<class M1>
struct get_sub_functor
{
    using value_type    = typename M1::value_type;
    using struct_type   = typename M1::struct_type;

    static void eval(const M1& A,Integer i,Integer j, Matrix& ret)
    {
        get_impl<M1,struct_type>::eval(A,i,j,ret);
        return;
    };

    static void eval(const M1& A,Integer i, Matrix& ret)
    {
        get_impl<M1,struct_type>::eval(A,i,ret);
        return;
    }

    static void eval(const Matrix& A, const M1& mat,const colon& c1, Matrix& ret)
    {
        if (c1.m_flag == colon::t_all)
        {
            ret = vec(A);
            return;
        };
        get_impl<M1,struct_type>::eval(A,mat,c1,ret);
        return;
    }

    static void eval(const Matrix& A, const M1& mat,const colon& c1,const colon& c2, Matrix& ret)
    {
        if (c1.m_flag == colon::t_all && c2.m_flag == colon::t_all)
        {
            ret = Matrix(A);
            return;
        };
        get_impl<M1,struct_type>::eval(A,mat,c1,c2,ret);
        return;
    }
};

template<class M1>
struct get_sub_functor_scal
{
    static void eval(const M1& A,Integer i,Integer j, Matrix& ret)
    {
        error::check_index(i, j, 1, 1);
        ret = A;
        return;
    };

    static void eval(const M1& A,Integer i, Matrix& ret)
    {
        error::check_index(i, 1);
        ret = A;
        return;
    }

    static void eval(const Matrix& A, const M1& mat,const colon& c1, Matrix& ret)
    {
        colon_info c_info;
        make_index(1,1,c1,c_info);

        Integer s = c_info.rows();
        if (s == 1)
        {
            ret = Matrix(A);
            return;
        };
        
        ti::ti_type<M1> ti_A = get_matrix_ti<M1>::eval(A);

        if (mrd::is_zero(mat))
        {
            using SparseMatrix = raw::Matrix<M1,struct_sparse>;
            SparseMatrix out(ti_A,c_info.rep_rows(),c_info.rep_cols());
            ret = Matrix(std::move(out),false);
            return;
        }
        else
        {
            using FullMatrix = raw::Matrix<M1,struct_dense>;
            FullMatrix out(ti_A,mat,c_info.rep_rows(),c_info.rep_cols());
            ret = Matrix(std::move(out),false);
            return;
        };        
    }

    static void eval(const Matrix& A, const M1& mat,const colon& c1,const colon& c2, Matrix& ret)
    {
        colon_info ci;
        make_index(1,1,c1,c2,ci);

        Integer r = ci.rows(), c = ci.cols();
        if (r == 1 && c == 1)
        {
            ret = Matrix(A);
            return;
        };
        
        ti::ti_type<M1> ti_A = get_matrix_ti<M1>::eval(A);

        if (mrd::is_zero(mat))
        {
            using SparseMatrix = raw::Matrix<M1,struct_sparse>;
            SparseMatrix out(ti_A,r,c);
            ret = Matrix(std::move(out),false);
            return;
        }
        else
        {
            using FullMatrix = raw::Matrix<M1,struct_dense>;
            FullMatrix out(ti_A,mat,r,c);
            ret = Matrix(std::move(out),false);
            return;
        };
    }
};

};};