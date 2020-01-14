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

#include "matcl-matrep/algs/sparse_algs.h"
#include "matcl-matrep/algs/sparse_algs_utils.h"
#include "matcl-internals/base/sort.h"
#include "matcl-internals/algs/scatter.h"
#include "matcl-internals/container/mat_s.h"
#include "matcl-scalar/details/scalfunc_helpers.h"
#include "matcl-matrep/details/struct_flag_predefined.h"
#include "matcl-matrep/lib_functions/func_binary.h"

namespace matcl { namespace algorithm { namespace details
{

namespace mr = matcl::raw;
namespace mrd = matcl::raw::details;
namespace md = matcl::details;

class column_iterator
{
    private:
        Integer					m_flag;
        Integer					m_size;
        Integer					m_dpos;

        raw::integer_dense		m_vec;
        Integer					m_pos;	

    public:
        column_iterator(const md::colon_info& c_inf)
            :m_vec(ti::ti_int())
        {
            m_flag = c_inf.c_flag;
            if (m_flag == 0)
            {
                m_vec.assign_to_fresh(c_inf.get_cim_2());
    
                sort_type s_type = is_sorted(m_vec);
                m_dpos	= 1;
                if (s_type == sorted_increasing)
                {
                }
                else if (s_type == sorted_decreasing)
                {
                    m_dpos = -1;
                } 
                else
                {
                    raw::integer_dense ci2 = m_vec.copy();
                    ci2.get_struct().reset();
                    utils::sort_q(ci2.ptr(),ci2.size());
                    m_vec.assign_to_fresh(std::move(ci2));
                };

                if (m_dpos > 0)
                    m_pos = 1;
                else
                    m_pos = m_vec.size();
            }
            else
            {
                m_dpos = c_inf.c_step;
                if (m_dpos > 0)
                {
                    m_pos = c_inf.c_start;
                }
                else
                {
                    m_dpos = -m_dpos;
                    m_pos = c_inf.c_end;
                };
            };
            m_size = c_inf.cols();
        };

        bool valid() const
        {
            return m_size > 0;
        };

        Integer get() const
        {
            if (m_flag == 0)
                return m_vec.ptr()[m_pos-1];
            else
                return m_pos;
        };

        void next()
        {
            m_pos += m_dpos;
            --m_size;
        };
};

template<class SM>
struct add_entries_impl
{
    using value_type = typename SM::value_type;

    static void eval_0(Matrix& ret, const SM& A,const md::colon_info& ci)
    {
        using value_type = typename SM::value_type;

        bool single = ci.is_double_mat_colon() == false;

        if (single == false)
            return eval_0_dc(ret, A, ci);

        Integer r = A.rows();
        Integer c = A.cols();

        const raw::details::sparse_ccs<value_type>& Ad = A.rep();
        
        mr::const_matrix<raw::integer_dense> ri = ci.get_rim_1();        

        sort_type s_type    = is_sorted(ri.get());
        bool incr           = true;

        if (s_type == sorted_increasing)
        {
        }
        else if (s_type == sorted_decreasing)
        {
            incr = false;
        } 
        else if (single == true)
        {
            raw::integer_dense ri2 = ri.get().copy();
            ri2.get_struct().reset();
            utils::sort_q(ri2.ptr(),ri2.size());
            ri.rebind(std::move(ri2));
        };

        const Integer* ptr_ri = ri.get().ptr();

        raw::details::sparse_ccs<value_type> d(A.get_type(), r, c, A.nnz() + ri.get().size());

        Integer nz				= 0;

        Integer * d_c			= d.ptr_c();
        Integer * d_r			= d.ptr_r();
        value_type * d_x		= d.ptr_x();

        const Integer * Ad_c	= Ad.ptr_c();
        const Integer * Ad_r	= Ad.ptr_r();
        const value_type * Ad_x	= Ad.ptr_x();
        value_type Z            = md::default_value<value_type>(A.get_type());

        Integer step_ri;
        Integer pos, pos_row, pos_col;
        if (incr)
        {
            pos     = 1;
            step_ri = 1;
        }
        else
        {
            pos     = ri.get().size();
            step_ri = -1;
        };

        Integer or      = ri.get().size();
        Integer old_pos = ptr_ri[pos-1];
        md::pos2ind(old_pos,r,pos_row,pos_col);

        for (Integer j = 0; j < c; ++j)
        {
            d_c[j] = nz;

            if (j < pos_col)
            {
                for (Integer k = Ad_c[j]; k < Ad_c[j+1]; ++k)
                {
                    d_r[nz] = Ad_r[k];
                    mrd::assign_helper(d_x[nz],Ad_x[k]);
                    ++nz;
                };
                continue;
            };

            Integer i;
            for (i = Ad_c[j]; i < Ad_c[j+1];)
            {
                Integer p = Ad_r[i];
                if (p < pos_row)
                {
                    d_r[nz] = p;
                    mrd::assign_helper(d_x[nz],Ad_x[i]);
                    ++nz;
                    ++i;
                    continue;
                };

                while (p > pos_row)
                {
                    d_r[nz] = pos_row;
                    mrd::assign_helper(d_x[nz],Z);
                    ++nz;

                    do
                    {
                        pos+=step_ri;
                        if (pos > or || pos < 1)
                        {
                            pos_row = r;
                            pos_col = c;
                            goto exit_row_flag;
                        };

                        md::pos2ind(ptr_ri[pos-1],r,pos_row,pos_col);
                        if (j != pos_col)
                        {
                            old_pos = ptr_ri[pos-1];
                            goto exit_row_flag;
                        };
                    }
                    while (ptr_ri[pos-1] == old_pos);
                    old_pos = ptr_ri[pos-1];
                };

                if (p == pos_row)
                {
                    d_r[nz] = pos_row;
                    mrd::assign_helper(d_x[nz],Ad_x[i]);
                    ++nz;
                    ++i;

                    do
                    {
                        pos+=step_ri;
                        if (pos > or || pos < 1)
                        {
                            pos_row = r;
                            pos_col = c;
                            goto exit_row_flag;
                        };

                        md::pos2ind(ptr_ri[pos-1],r,pos_row,pos_col);
                        if (j != pos_col)
                        {
                            old_pos = ptr_ri[pos-1];
                            goto exit_row_flag;
                        };
                    }
                    while ( p == pos_row);
                    old_pos = ptr_ri[pos-1];
                };
            };

            while (j == pos_col)
            {
                d_r[nz] = pos_row;
                mrd::assign_helper(d_x[nz],Z);
                ++nz;

                do
                {
                    pos+=step_ri;
                    if (pos > or || pos < 1)
                    {
                        pos_row = r;
                        pos_col = c;
                        goto exit_row_flag;
                    };

                    md::pos2ind(ptr_ri[pos-1],r,pos_row,pos_col);
                }
                while (ptr_ri[pos-1] == old_pos);

                old_pos = ptr_ri[pos-1];
            };

            continue;

            exit_row_flag:

            for (; i < Ad_c[j+1];++i)
            {
                d_r[nz] = Ad_r[i];
                mrd::assign_helper(d_x[nz],Ad_x[i]);
                ++nz;
            };
        };

        d_c[c] = nz;

        SM ret_l = SM(raw::sparse_matrix_base<value_type>(d));
        ret_l.set_struct(A.get_struct());

        ret         = Matrix(ret_l, false);
        return;
    };

    static void eval_0_dc(Matrix& ret, const SM& A,const md::colon_info& ci)
    {
        using value_type = typename SM::value_type;

        Integer r = A.rows();
        Integer c = A.cols();

        const raw::details::sparse_ccs<value_type>& Ad = A.rep();
        
        raw::integer_dense cr   = ci.get_rim_r().copy();
        raw::integer_dense cc   = ci.get_rim_c().copy();

        sort_rows_cols(cr.ptr(), cc.ptr(), cr.size());

        const Integer* ptr_r    = cr.ptr();
        const Integer* ptr_c    = cc.ptr();

        raw::details::sparse_ccs<value_type> d(A.get_type(),r, c, A.nnz() + cr.size());

        Integer nz				= 0;

        Integer * d_c			= d.ptr_c();
        Integer * d_r			= d.ptr_r();
        value_type * d_x		= d.ptr_x();

        const Integer * Ad_c	= Ad.ptr_c();
        const Integer * Ad_r	= Ad.ptr_r();
        const value_type * Ad_x	= Ad.ptr_x();
        value_type Z            = md::default_value<value_type>(A.get_type());

        Integer pos             = 1;
        Integer pos_row         = ptr_r[0] - 1;
        Integer pos_col         = ptr_c[0] - 1;
        Integer or              = cr.size();

        Integer old_pos_row     = pos_row;
        Integer old_pos_col     = pos_col;

        for (Integer j = 0; j < c; ++j)
        {
            d_c[j] = nz;

            if (j < pos_col)
            {
                for (Integer k = Ad_c[j]; k < Ad_c[j+1]; ++k)
                {
                    d_r[nz] = Ad_r[k];
                    mrd::assign_helper(d_x[nz],Ad_x[k]);
                    ++nz;
                };
                continue;
            };

            Integer i;
            for (i = Ad_c[j]; i < Ad_c[j+1];)
            {
                Integer p = Ad_r[i];
                if (p < pos_row)
                {
                    d_r[nz] = p;
                    mrd::assign_helper(d_x[nz],Ad_x[i]);
                    ++nz;
                    ++i;
                    continue;
                };

                while (p > pos_row)
                {
                    d_r[nz] = pos_row;
                    mrd::assign_helper(d_x[nz],Z);
                    ++nz;

                    do
                    {
                        pos += 1;

                        if (pos > or || pos < 1)
                        {
                            pos_row = r;
                            pos_col = c;
                            goto exit_row_flag;
                        };

                        pos_row = ptr_r[pos - 1] - 1;
                        pos_col = ptr_c[pos - 1] - 1;

                        if (j != pos_col)
                        {
                            old_pos_row = pos_row;
                            old_pos_col = pos_col;
                            goto exit_row_flag;
                        };
                    }
                    while (pos_row == old_pos_row && pos_col == old_pos_col);

                    old_pos_row = pos_row;
                    old_pos_col = pos_col;
                };

                if (p == pos_row)
                {
                    d_r[nz] = pos_row;
                    mrd::assign_helper(d_x[nz],Ad_x[i]);
                    ++nz;
                    ++i;

                    do
                    {
                        pos += 1;
                        if (pos > or || pos < 1)
                        {
                            pos_row = r;
                            pos_col = c;
                            goto exit_row_flag;
                        };

                        pos_row = ptr_r[pos - 1] - 1;
                        pos_col = ptr_c[pos - 1] - 1;
                        
                        if (j != pos_col)
                        {
                            old_pos_row = pos_row;
                            old_pos_col = pos_col;

                            goto exit_row_flag;
                        };
                    }
                    while ( p == pos_row);
                    
                    old_pos_row = pos_row;
                    old_pos_col = pos_col;
                };
            };

            while (j == pos_col)
            {
                d_r[nz] = pos_row;
                mrd::assign_helper(d_x[nz],Z);
                ++nz;

                do
                {
                    pos += 1;
                    if (pos > or || pos < 1)
                    {
                        pos_row = r;
                        pos_col = c;
                        goto exit_row_flag;
                    };

                    pos_row = ptr_r[pos - 1] - 1;
                    pos_col = ptr_c[pos - 1] - 1;
                }
                while (pos_row == old_pos_row && pos_col == old_pos_col);

                old_pos_row = pos_row;
                old_pos_col = pos_col;
            };

            continue;

            exit_row_flag:

            for (; i < Ad_c[j+1];++i)
            {
                d_r[nz] = Ad_r[i];
                mrd::assign_helper(d_x[nz],Ad_x[i]);
                ++nz;
            };
        };

        d_c[c] = nz;

        SM ret_l    = SM(raw::sparse_matrix_base<value_type>(d));
        ret_l.set_struct(A.get_struct());

        ret         = Matrix(ret_l, false);
        return;
    };

    static void eval_1(Matrix& ret, const SM& A,const md::colon_info& ci)
    {
        using value_type = typename SM::value_type;

        Integer r   = A.rows();
        Integer c   = A.cols();
        const raw::details::sparse_ccs<value_type>& Ad = A.rep();

        Integer rs = ci.r_step,rf,rl;
        if (rs < 0)
        {
            rs = -rs;
            rf = ci.r_end;
            rl = ci.r_start;
        }
        else
        {
            rl = ci.r_end;
            rf = ci.r_start;
        };

        raw::details::sparse_ccs<value_type> d(A.get_type(),r, c, A.nnz() + ci.rows());        

        Integer nz				= 0;

        Integer * d_c			= d.ptr_c();
        Integer * d_r			= d.ptr_r();
        value_type * d_x		= d.ptr_x();

        const Integer * Ad_c	= Ad.ptr_c();
        const Integer * Ad_r	= Ad.ptr_r();
        const value_type * Ad_x	= Ad.ptr_x();

        value_type Z            = md::default_value<value_type>(A.get_type());

        Integer pos_row, pos_col, pos = rf;
        md::pos2ind(rf,r,pos_row,pos_col);

        for (Integer j = 0; j < c; ++j)
        {
            d_c[j] = nz;

            if (j < pos_col)
            {
                for (Integer k = Ad_c[j]; k < Ad_c[j+1]; ++k)
                {
                    d_r[nz] = Ad_r[k];
                    mrd::assign_helper(d_x[nz],Ad_x[k]);
                    ++nz;
                };
                continue;
            };

            Integer i;
            for (i = Ad_c[j]; i < Ad_c[j+1];)
            {
                Integer p = Ad_r[i];
                if (p < pos_row)
                {
                    d_r[nz] = p;
                    mrd::assign_helper(d_x[nz],Ad_x[i]);
                    ++nz;
                    ++i;
                    continue;
                };

                while (p > pos_row)
                {
                    d_r[nz] = pos_row;
                    mrd::assign_helper(d_x[nz],Z);
                    ++nz;

                    pos+=rs;
                    if (pos > rl)
                    {
                        pos_row = r;
                        pos_col = c;
                        goto exit_row_flag;
                    };
                    md::pos2ind(pos,r,pos_row,pos_col);

                    if (j != pos_col)
                        goto exit_row_flag;
                };
                
                if (p == pos_row)
                {
                    d_r[nz] = pos_row;
                    mrd::assign_helper(d_x[nz],Ad_x[i]);
                    ++nz;
                    ++i;

                    do
                    {
                        pos+=rs;
                        if (pos > rl)
                        {
                            pos_row = r;
                            pos_col = c;
                            goto exit_row_flag;
                        };

                        md::pos2ind(pos,r,pos_row,pos_col);
                        if (j != pos_col)
                            goto exit_row_flag;
                    }
                    while(p == pos_row);
                };
            };

            while (j == pos_col)
            {
                d_r[nz] = pos_row;
                mrd::assign_helper(d_x[nz],Z);
                ++nz;

                pos+= rs;
                if (pos > rl)
                {
                    pos_row = r;
                    pos_col = c;
                    goto exit_row_flag;
                };
                md::pos2ind(pos,r,pos_row,pos_col);
            };

            continue;

            exit_row_flag:

            for (; i < Ad_c[j+1];++i)
            {
                d_r[nz] = Ad_r[i];
                mrd::assign_helper(d_x[nz],Ad_x[i]);
                ++nz;
            };
        };

        d_c[c] = nz;

        SM ret_l    = SM(raw::sparse_matrix_base<value_type>(d));
        ret_l.set_struct(A.get_struct());

        ret         = Matrix(ret_l, false);
        return;
    };	

    static void eval_02(Matrix& ret, const SM& A,const md::colon_info& c_in)
    {
        using value_type = typename SM::value_type;

        Integer r   = A.rows();
        Integer c   = A.cols();

        const raw::details::sparse_ccs<value_type>& Ad = A.rep();

        ti::ti_type<value_type> ti = A.get_type();

        const Integer * Ad_c		= Ad.ptr_c();
        const Integer * Ad_r		= Ad.ptr_r();
        const value_type * Ad_x		= Ad.ptr_x();	

        Integer nzA = icast_c(Real(A.nnz()) + Real(c_in.rows())*Real(c_in.cols()));
        raw::details::sparse_ccs<value_type> d(A.get_type(),r, c, nzA);

        Integer * d_c				= d.ptr_c();
        Integer * d_r				= d.ptr_r();
        value_type * d_x			= d.ptr_x();

        column_iterator c_it(c_in);

        using workspace     = md::workspace2<value_type>;
        using scatter       = matcl::algorithm::scatter;

        workspace                   work_x(ti,r);
        
        scatter sc                  = scatter::get(r, c);
        const raw::integer_dense& ri= c_in.get_rim_2();
        const Integer* ptr_ri       = ri.ptr();
        value_type Z                = md::default_value<value_type>(A.get_type());

        Integer nz = 0;

        for (Integer j = 0; j < c; ++j)
        {
            d_c[j] = nz;

            if (!c_it.valid() || c_it.get() > j+1)
            {
                for (Integer k = Ad_c[j]; k < Ad_c[j+1]; ++k)
                {
                    d_r[nz] = Ad_r[k];
                    mrd::assign_helper(d_x[nz],Ad_x[k]);
                    ++nz;
                };
                continue;
            };

            while (c_it.valid() && c_it.get() == j+1)
                c_it.next();

            Integer nz_old		= nz;
            bool b_added		= false;
            auto mark           = sc.next_mark();

            for (Integer k = Ad_c[j]; k < Ad_c[j+1]; ++k)
            {
                Integer p		= Ad_r[k];
                sc[p]		    = mark;			        // p is new entry in column j //
                d_r[nz++]		= p;			        // add p to pattern of d(:,j) //
                mrd::assign_helper(work_x[p],Ad_x[k]);	// x(p) = A(p,j) //
            };

            Integer ri_size     = ri.size();
            
            for (Integer k = 0; k < ri_size; ++k)
            {
                Integer p		= ptr_ri[k]-1;
                if (sc[p] < mark)
                {
                    d_r[nz++]	= p;		            // add p to pattern of d(:,j) //
                    mrd::assign_helper(work_x[p],Z);	// x(p) = val //
                    sc[p]       = mark;
                    b_added		= true;
                };
            };
            Integer nz_new		= nz - nz_old;

            if (b_added)
                utils::sort_q(d_r+nz_old,nz_new);

            for (Integer k = nz_old; k < nz ; ++k) 
                mrd::assign_helper(d_x[k],work_x[d_r[k]]);
        };

        d_c[c] = nz;

        SM ret_l = SM(raw::sparse_matrix_base<value_type>(d));
        ret_l.set_struct(A.get_struct());

        ret     = Matrix(ret_l, false);
        return;
    };

    static void eval_12(Matrix& ret, const SM& A,const md::colon_info& c_in)
    {
        if (c_in.r_start == 1 && c_in.r_end == A.rows() && c_in.r_step == 1)
            return add_cols(ret, A, c_in);

        using value_type = typename SM::value_type;

        Integer r = A.rows(), c = A.cols();
        const raw::details::sparse_ccs<value_type>& Ad = A.rep();

        const Integer * Ad_c		= Ad.ptr_c();
        const Integer * Ad_r		= Ad.ptr_r();
        const value_type * Ad_x		= Ad.ptr_x();	

        Integer nzA = icast_c(Real(A.nnz()) + Real(c_in.rows())*Real(c_in.cols()));
        raw::details::sparse_ccs<value_type> d(A.get_type(),r, c, nzA);

        Integer * d_c				= d.ptr_c();
        Integer * d_r				= d.ptr_r();
        value_type * d_x			= d.ptr_x();
        value_type Z                = md::default_value<value_type>(A.get_type());

        column_iterator c_it(c_in);

        Integer nz					= 0;

        Integer rs = c_in.r_step, rf, rl;
        if (rs < 0)
        {
            rs = -rs;
            rf = c_in.r_end;
            rl = c_in.r_start;
        }
        else
        {
            rf = c_in.r_start;		
            rl = c_in.r_end;
        };
        --rf;
        --rl;
        
        for (Integer j = 0; j < c; ++j)
        {
            d_c[j]					= nz;

            if (!c_it.valid() || c_it.get() > j+1)
            {
                for (Integer k = Ad_c[j]; k < Ad_c[j+1]; ++k)
                {
                    d_r[nz] = Ad_r[k];
                    mrd::assign_helper(d_x[nz],Ad_x[k]);
                    ++nz;
                };
                continue;
            };

            while (c_it.valid() && c_it.get() == j+1)
                c_it.next();

            Integer la			= Ad_c[j+1];
            Integer i			= Ad_c[j];

            Integer k = rf;
            while(i < la && k <= rl)
            {
                Integer p		= Ad_r[i];

                if (p < k)
                {
                    d_r[nz]		= Ad_r[i];
                    mrd::assign_helper(d_x[nz],Ad_x[i]);
                    ++nz;
                    ++i;
                }
                else if (p > k)
                {
                    d_r[nz]		= k;
                    mrd::assign_helper(d_x[nz],Z);
                    ++nz;
                    k			+= rs;
                }
                else
                {
                    d_r[nz]		= k;
                    mrd::assign_helper(d_x[nz],Ad_x[i]);
                    ++nz;
                    k			+= rs;
                    ++i;
                };
            };

            while(i < la)
            {
                d_r[nz]			= Ad_r[i];
                mrd::assign_helper(d_x[nz],Ad_x[i]);
                ++nz;
                ++i;
            };
            
            while(k <= rl)
            {
                d_r[nz]			= k;
                mrd::assign_helper(d_x[nz],Z);
                ++nz;
                k				+= rs;
            };
        };
        d_c[c] = nz;

        SM ret_l = SM(raw::sparse_matrix_base<value_type>(d));
        ret_l.set_struct(A.get_struct());

        ret     = Matrix(ret_l, false);
        return;
    };

    static void add_cols(Matrix& ret, const SM& A,const md::colon_info& c_in)
    {
        using value_type = typename SM::value_type;

        Integer r = A.rows(), c = A.cols();
        const raw::details::sparse_ccs<value_type>& Ad = A.rep();

        const Integer * Ad_c		= Ad.ptr_c();
        const Integer * Ad_r		= Ad.ptr_r();
        const value_type * Ad_x		= Ad.ptr_x();	

        Integer nzA = icast_c(Real(A.nnz()) + Real(c_in.rows())*Real(c_in.cols()));
        raw::details::sparse_ccs<value_type> d(A.get_type(),r, c, nzA);

        Integer * d_c				= d.ptr_c();
        Integer * d_r				= d.ptr_r();
        value_type * d_x			= d.ptr_x();

        column_iterator c_it(c_in);

        Integer nz					= 0;
        value_type Z                = md::default_value<value_type>(A.get_type());

        for (Integer j = 0; j < c; ++j)
        {
            d_c[j]					= nz;

            if (!c_it.valid() || c_it.get() > j+1)
            {
                for (Integer k = Ad_c[j]; k < Ad_c[j+1]; ++k)
                {
                    d_r[nz] = Ad_r[k];
                    mrd::assign_helper(d_x[nz],Ad_x[k]);
                    ++nz;
                };
                continue;
            };

            while (c_it.valid() && c_it.get() == j+1)
                c_it.next();

            value_type * d_x_sav    = d_x + nz;

            for (Integer i = 0; i < r; ++i)
            {
                d_r[nz] = i;
                mrd::assign_helper(d_x[nz],Z);
                ++nz;
            };            

            for (Integer k = Ad_c[j]; k < Ad_c[j+1]; ++k)
            {
                Integer pos = Ad_r[k];
                mrd::assign_helper(d_x_sav[pos],Ad_x[k]);
            };
        };
        d_c[c] = nz;

        SM ret_l = SM(raw::sparse_matrix_base<value_type>(d));
        ret_l.set_struct(A.get_struct());

        ret     = Matrix(ret_l, false);
        return;
    };

};

template<class SM>
struct change_entries_impl
{
    using value_type = typename SM::value_type;

    static void eval_0(Matrix& ret, const SM& A,const md::colon_info& ci, const value_type& val)
    {
        bool single = ci.is_double_mat_colon() == false;

        if (single == false)
            return eval_0_dc(ret, A, ci, val);

        using value_type = typename SM::value_type;

        Integer r = A.rows();
        Integer c = A.cols();
        const raw::details::sparse_ccs<value_type>& Ad = A.rep();

        mr::const_matrix<raw::integer_dense> ri = ci.get_rim_1();

        sort_type s_type = is_sorted(ri.get());

        bool incr = true;
        if (s_type == sorted_increasing)
        {
        }
        else if (s_type == sorted_decreasing)
        {
            incr = false;
        } 
        else
        {
            raw::integer_dense ri2 = ri.get().copy();
            ri2.get_struct().reset();
            utils::sort_q(ri2.ptr(),ri2.size());
            ri.rebind(std::move(ri2));
        };

        const Integer* ptr_ri = ri.get().ptr();

        raw::details::sparse_ccs<value_type> d(A.get_type(),r, c, A.nnz() + ri.get().size());

        Integer nz				= 0;

        Integer * d_c			= d.ptr_c();
        Integer * d_r			= d.ptr_r();
        value_type * d_x		= d.ptr_x();

        const Integer * Ad_c	= Ad.ptr_c();
        const Integer * Ad_r	= Ad.ptr_r();
        const value_type * Ad_x	= Ad.ptr_x();

        Integer step_ri;
        Integer pos, pos_row, pos_col;
        if (incr)
        {
            pos = 1;
            step_ri = 1;
        }
        else
        {
            pos     = ri.get().size();
            step_ri = -1;
        };

        Integer or      = ri.get().size();
        Integer old_pos = ptr_ri[pos-1];
        md::pos2ind(old_pos,r,pos_row,pos_col);

        for (Integer j = 0; j < c; ++j)
        {
            d_c[j] = nz;

            if (j < pos_col)
            {
                for (Integer k = Ad_c[j]; k < Ad_c[j+1]; ++k)
                {
                    d_r[nz] = Ad_r[k];
                    mrd::assign_helper(d_x[nz],Ad_x[k]);
                    ++nz;
                };
                continue;
            };

            Integer i;
            for (i = Ad_c[j]; i < Ad_c[j+1];)
            {
                Integer p = Ad_r[i];
                if (p < pos_row)
                {
                    d_r[nz] = p;
                    mrd::assign_helper(d_x[nz],Ad_x[i]);
                    ++nz;
                    ++i;
                    continue;
                };
                while (p > pos_row)
                {
                    d_r[nz] = pos_row;
                    mrd::assign_helper(d_x[nz],val);
                    ++nz;

                    do
                    {
                        pos+=step_ri;
                        if (pos > or || pos < 1)
                        {
                            pos_row = r;
                            pos_col = c;
                            goto exit_row_flag;
                        };
                        md::pos2ind(ptr_ri[pos-1],r,pos_row,pos_col);
                        if (j != pos_col)
                        {
                            old_pos = ptr_ri[pos-1];
                            goto exit_row_flag;
                        };
                    }
                    while (ptr_ri[pos-1] == old_pos);
                    old_pos = ptr_ri[pos-1];
                };

                if ( p == pos_row)
                {
                    d_r[nz] = pos_row;
                    mrd::assign_change_helper(d_x[nz],Ad_x[i],val);
                    ++nz;
                    ++i;

                    do
                    {
                        pos+=step_ri;
                        if (pos > or || pos < 1)
                        {
                            pos_row = r;
                            pos_col = c;
                            goto exit_row_flag;
                        };
                        md::pos2ind(ptr_ri[pos-1],r,pos_row,pos_col);
                        if (j != pos_col)
                        {
                            old_pos = ptr_ri[pos-1];
                            goto exit_row_flag;
                        };
                    }
                    while ( p == pos_row);
                    old_pos = ptr_ri[pos-1];
                };
            };

            while (j == pos_col)
            {
                d_r[nz] = pos_row;
                mrd::assign_helper(d_x[nz],val);
                ++nz;

                do
                {
                    pos+=step_ri;
                    if (pos > or || pos < 1)
                    {
                        pos_row = r;
                        pos_col = c;
                        goto exit_row_flag;
                    };
                    md::pos2ind(ptr_ri[pos-1],r,pos_row,pos_col);
                }
                while (ptr_ri[pos-1] == old_pos);
                old_pos = ptr_ri[pos-1];
            };

            continue;

            exit_row_flag:

            for (; i < Ad_c[j+1];++i)
            {
                d_r[nz] = Ad_r[i];
                mrd::assign_helper(d_x[nz],Ad_x[i]);
                ++nz;
            };
        };
        d_c[c] = nz;

        ret = Matrix(raw::sparse_matrix_base<value_type>(d), false);
        return;
    };

    static void eval_0_dc(Matrix& ret, const SM& A,const md::colon_info& ci, 
                          const value_type& val)
    {
        using value_type = typename SM::value_type;

        Integer r = A.rows();
        Integer c = A.cols();
        const raw::details::sparse_ccs<value_type>& Ad = A.rep();

        raw::integer_dense cr   = ci.get_rim_r().copy();
        raw::integer_dense cc   = ci.get_rim_c().copy();

        sort_rows_cols(cr.ptr(), cc.ptr(), cr.size());

        const Integer* ptr_r    = cr.ptr();
        const Integer* ptr_c    = cc.ptr();

        raw::details::sparse_ccs<value_type> d(A.get_type(),r, c, A.nnz() + cr.size());

        Integer nz				= 0;

        Integer * d_c			= d.ptr_c();
        Integer * d_r			= d.ptr_r();
        value_type * d_x		= d.ptr_x();

        const Integer * Ad_c	= Ad.ptr_c();
        const Integer * Ad_r	= Ad.ptr_r();
        const value_type * Ad_x	= Ad.ptr_x();

        Integer pos             = 1;
        Integer or              = cr.size();

        Integer pos_row         = ptr_r[0] - 1;
        Integer pos_col         = ptr_c[0] - 1;
        Integer old_pos_row     = pos_row;
        Integer old_pos_col     = pos_col;

        for (Integer j = 0; j < c; ++j)
        {
            d_c[j] = nz;

            if (j < pos_col)
            {
                for (Integer k = Ad_c[j]; k < Ad_c[j+1]; ++k)
                {
                    d_r[nz] = Ad_r[k];
                    mrd::assign_helper(d_x[nz],Ad_x[k]);
                    ++nz;
                };
                continue;
            };

            Integer i;
            for (i = Ad_c[j]; i < Ad_c[j+1];)
            {
                Integer p = Ad_r[i];
                if (p < pos_row)
                {
                    d_r[nz] = p;
                    mrd::assign_helper(d_x[nz],Ad_x[i]);
                    ++nz;
                    ++i;
                    continue;
                };
                while (p > pos_row)
                {
                    d_r[nz] = pos_row;
                    mrd::assign_helper(d_x[nz],val);
                    ++nz;

                    do
                    {
                        pos += 1;
                        if (pos > or || pos < 1)
                        {
                            pos_row = r;
                            pos_col = c;
                            goto exit_row_flag;
                        };

                        pos_row = ptr_r[pos-1] - 1;
                        pos_col = ptr_c[pos-1] - 1;

                        if (j != pos_col)
                        {
                            old_pos_row = pos_row;
                            old_pos_col = pos_col;
                            goto exit_row_flag;
                        };
                    }
                    while (pos_row == old_pos_row && pos_col == old_pos_col);

                    old_pos_row = pos_row;
                    old_pos_col = pos_col;
                };

                if ( p == pos_row)
                {
                    d_r[nz] = pos_row;
                    mrd::assign_change_helper(d_x[nz],Ad_x[i],val);
                    ++nz;
                    ++i;

                    do
                    {
                        pos += 1;
                        if (pos > or || pos < 1)
                        {
                            pos_row = r;
                            pos_col = c;
                            goto exit_row_flag;
                        };

                        pos_row = ptr_r[pos-1] - 1;
                        pos_col = ptr_c[pos-1] - 1;

                        if (j != pos_col)
                        {
                            old_pos_row = pos_row;
                            old_pos_col = pos_col;

                            goto exit_row_flag;
                        };
                    }
                    while (p == pos_row);

                    old_pos_row = pos_row;
                    old_pos_col = pos_col;
                };
            };

            while (j == pos_col)
            {
                d_r[nz] = pos_row;
                mrd::assign_helper(d_x[nz],val);
                ++nz;

                do
                {
                    pos += 1;
                    if (pos > or || pos < 1)
                    {
                        pos_row = r;
                        pos_col = c;
                        goto exit_row_flag;
                    };

                    pos_row = ptr_r[pos-1] - 1;
                    pos_col = ptr_c[pos-1] - 1;
                }
                while (pos_row == old_pos_row && pos_col == old_pos_col);

                old_pos_row = pos_row;
                old_pos_col = pos_col;
            };

            continue;

            exit_row_flag:

            for (; i < Ad_c[j+1];++i)
            {
                d_r[nz] = Ad_r[i];
                mrd::assign_helper(d_x[nz],Ad_x[i]);
                ++nz;
            };
        };
        d_c[c] = nz;

        ret = Matrix(raw::sparse_matrix_base<value_type>(d), false);
    };

    static void eval_1(Matrix& ret, const SM& A,const md::colon_info& ci, 
                       const value_type& val)
    {
        using value_type = typename SM::value_type;

        Integer r = A.rows(), c = A.cols();
        const raw::details::sparse_ccs<value_type>& Ad = A.rep();

        Integer rs = ci.r_step,rf,rl;
        if (rs < 0)
        {
            rs = -rs;
            rf = ci.r_end;
            rl = ci.r_start;
        }
        else
        {
            rl = ci.r_end;
            rf = ci.r_start;
        };

        raw::details::sparse_ccs<value_type> d(A.get_type(),r, c, A.nnz() + ci.rows());

        Integer nz				= 0;

        Integer * d_c			= d.ptr_c();
        Integer * d_r			= d.ptr_r();
        value_type * d_x		= d.ptr_x();

        const Integer * Ad_c	= Ad.ptr_c();
        const Integer * Ad_r	= Ad.ptr_r();
        const value_type * Ad_x	= Ad.ptr_x();

        Integer pos_row, pos_col, pos = rf;
        md::pos2ind(rf,r,pos_row,pos_col);

        for (Integer j = 0; j < c; ++j)
        {
            d_c[j] = nz;

            if (j < pos_col)
            {
                for (Integer k = Ad_c[j]; k < Ad_c[j+1]; ++k)
                {
                    d_r[nz] = Ad_r[k];
                    mrd::assign_helper(d_x[nz],Ad_x[k]);
                    ++nz;
                };
                continue;
            };

            Integer i;
            for (i = Ad_c[j]; i < Ad_c[j+1];)
            {
                Integer p = Ad_r[i];
                if (p < pos_row)
                {
                    d_r[nz] = p;
                    mrd::assign_helper(d_x[nz],Ad_x[i]);
                    ++nz;
                    ++i;
                    continue;
                };
                while (p > pos_row)
                {
                    d_r[nz] = pos_row;
                    mrd::assign_helper(d_x[nz],val);
                    ++nz;

                    pos+=rs;
                    if (pos > rl)
                    {
                        pos_row = r;
                        pos_col = c;
                        goto exit_row_flag;
                    };
                    md::pos2ind(pos,r,pos_row,pos_col);
                    if (j != pos_col)
                    {
                        goto exit_row_flag;
                    };
                };

                if ( p == pos_row)
                {
                    d_r[nz] = pos_row;
                    mrd::assign_change_helper(d_x[nz],Ad_x[i],val);
                    ++nz;
                    ++i;

                    do
                    {
                        pos+=rs;
                        if (pos > rl)
                        {
                            pos_row = r;
                            pos_col = c;
                            goto exit_row_flag;
                        };

                        md::pos2ind(pos,r,pos_row,pos_col);
                        if (j != pos_col)
                            goto exit_row_flag;
                    }
                    while(p == pos_row);
                };
            };

            while (j == pos_col)
            {
                d_r[nz] = pos_row;
                mrd::assign_helper(d_x[nz],val);
                ++nz;

                pos+= rs;
                if (pos > rl)
                {
                    pos_row = r;
                    pos_col = c;
                    goto exit_row_flag;
                };
                md::pos2ind(pos,r,pos_row,pos_col);
            };

            continue;

            exit_row_flag:

            for (; i < Ad_c[j+1];++i)
            {
                d_r[nz] = Ad_r[i];
                mrd::assign_helper(d_x[nz],Ad_x[i]);
                ++nz;
            };
        };

        d_c[c] = nz;
        ret = Matrix(raw::sparse_matrix_base<value_type>(d), false);
    };	

    static void eval_02(Matrix& ret, const SM& A,const md::colon_info& c_in, 
                        const value_type& val)
    {
        if (imult(c_in.rows(),c_in.cols()) < A.rows())
            return eval_02w(ret, A, c_in, val);
        else
            return eval_02s(ret, A, c_in, val);
    };

    static void eval_02s(Matrix& ret, const SM& A,const md::colon_info& c_in, 
                         const value_type& val)
    {
        using value_type = typename SM::value_type;

        Integer r = A.rows(), c = A.cols();
        const raw::details::sparse_ccs<value_type>& Ad = A.rep();

        ti::ti_type<value_type> ti = A.get_type();

        const Integer * Ad_c		= Ad.ptr_c();
        const Integer * Ad_r		= Ad.ptr_r();
        const value_type * Ad_x		= Ad.ptr_x();	

        Integer nzA = icast_c(Real(A.nnz()) + Real(c_in.rows())*Real(c_in.cols()));
        raw::details::sparse_ccs<value_type> d(A.get_type(),r, c, nzA);

        Integer * d_c				= d.ptr_c();
        Integer * d_r				= d.ptr_r();
        value_type * d_x			= d.ptr_x();

        column_iterator c_it(c_in);

        using workspace     = md::workspace2<value_type>;
        using scatter       = matcl::algorithm::scatter;

        workspace                   work_x(ti,r);
        
        scatter sc                      = scatter::get(r, c);
        const raw::integer_dense& ri    = c_in.get_rim_2();
        const Integer* ptr_ri           = ri.ptr();

        Integer nz = 0;

        for (Integer j = 0; j < c; ++j)
        {
            d_c[j] = nz;

            if (!c_it.valid() || c_it.get() > j+1)
            {
                for (Integer k = Ad_c[j]; k < Ad_c[j+1]; ++k)
                {
                    d_r[nz] = Ad_r[k];
                    mrd::assign_helper(d_x[nz],Ad_x[k]);
                    ++nz;
                };
                continue;
            };

            while (c_it.valid() && c_it.get() == j+1)
                c_it.next();

            Integer nz_old		= nz;
            bool b_added		= false;
            auto mark           = sc.next_mark();

            for (Integer k = Ad_c[j]; k < Ad_c[j+1]; ++k)
            {
                Integer p		= Ad_r[k];
                sc[p]		    = mark;			        // p is new entry in column j //
                d_r[nz++]		= p;			        // add p to pattern of d(:,j) //
                mrd::assign_helper(work_x[p],Ad_x[k]);	// x(p) = A(p,j) //
            };

            Integer ri_size     = ri.size();
            
            for (Integer k = 0; k < ri_size; ++k)
            {
                Integer p		= ptr_ri[k]-1;
                if (sc[p] < mark)
                {
                    d_r[nz++]	= p;		            // add p to pattern of d(:,j) //
                    mrd::assign_helper(work_x[p],val);	// x(p) = val //
                    sc[p]       = mark;
                    b_added		= true;
                }
                else
                {
                    mrd::assign_helper(work_x[p],val);  // p exists in d(:,j) //
                };
            };
            Integer nz_new		= nz - nz_old;

            if (b_added)
                utils::sort_q(d_r+nz_old,nz_new);

            for (Integer k = nz_old; k < nz ; ++k) 
                mrd::assign_helper(d_x[k],work_x[d_r[k]]);
        };

        d_c[c] = nz;

        ret = Matrix(raw::sparse_matrix_base<value_type>(d), false);
    };

    static void eval_02w(Matrix& ret, const SM& A,const md::colon_info& c_in, 
                         const value_type& val)
    {
        using value_type = typename SM::value_type;

        Integer r = A.rows(), c = A.cols();
        const raw::details::sparse_ccs<value_type>& Ad = A.rep();

        const Integer * Ad_c		= Ad.ptr_c();
        const Integer * Ad_r		= Ad.ptr_r();
        const value_type * Ad_x		= Ad.ptr_x();	

        Integer nzA = icast_c(Real(A.nnz()) + Real(c_in.rows())*Real(c_in.cols()));
        raw::details::sparse_ccs<value_type> d(A.get_type(),r, c, nzA);

        Integer * d_c				= d.ptr_c();
        Integer * d_r				= d.ptr_r();
        value_type * d_x			= d.ptr_x();

        column_iterator c_it(c_in);

        Integer nz					= 0;

        mr::const_matrix<raw::integer_dense> ri = c_in.get_rim_2();        
        sort_type rs_type			            = is_sorted(ri.get());
        
        if (rs_type != sorted_increasing)
        {
            raw::integer_dense ri2 = ri.get().copy();
            ri2.get_struct().reset();
            utils::sort_q(ri2.ptr(), ri2.size());
            ri.rebind(std::move(ri2));
        };

        const Integer* ptr_ri       = ri.get().ptr();
        Integer n					= ri.get().size();

        for (Integer j = 0; j < c; ++j)
        {
            d_c[j]					= nz;

            if (!c_it.valid() || c_it.get() > j+1)
            {
                for (Integer k = Ad_c[j]; k < Ad_c[j+1]; ++k)
                {
                    d_r[nz] = Ad_r[k];
                    mrd::assign_helper(d_x[nz],Ad_x[k]);
                    ++nz;
                };
                continue;
            };

            while (c_it.valid() && c_it.get() == j+1)
                c_it.next();

            Integer la			= Ad_c[j+1];
            Integer i			= Ad_c[j];

            Integer k = 0;
            while(i < la && k < n)
            {
                Integer p		= Ad_r[i];
                Integer pr		= ptr_ri[k]-1;

                if (p < pr)
                {
                    d_r[nz]		= Ad_r[i];
                    mrd::assign_helper(d_x[nz],Ad_x[i]);
                    ++nz;
                    ++i;
                }
                else 
                {	
                    if (p > pr)
                    {
                        d_r[nz]		= pr;
                        mrd::assign_helper(d_x[nz],val);
                        ++nz;
                        ++k;
                    }
                    else
                    {
                        d_r[nz]		= pr;
                        mrd::assign_change_helper(d_x[nz],Ad_x[i],val);
                        ++nz;
                        ++k;
                        ++i;
                    };

                    while ( k < n && ptr_ri[k] == ptr_ri[k-1])
                        ++k;
                };
            };

            while(i < la)
            {
                d_r[nz]			= Ad_r[i];
                mrd::assign_helper(d_x[nz],Ad_x[i]);
                ++nz;
                ++i;
            };
            
            while(k < n)
            {
                Integer pr		= ptr_ri[k]-1;

                d_r[nz]			= pr;
                mrd::assign_helper(d_x[nz],val);
                ++nz;
                ++k;

                while ( k < n && ptr_ri[k] == ptr_ri[k-1])
                    ++k;
            };
        };
        d_c[c] = nz;

        ret = Matrix(raw::sparse_matrix_base<value_type>(d), false);
    };

    static void eval_12(Matrix& ret, const SM& A,const md::colon_info& c_in, 
                        const value_type& val)
    {
        if (c_in.r_start == 1 && c_in.r_end == A.rows() && c_in.r_step == 1)
            return change_cols(ret, A, c_in, val);

        using value_type = typename SM::value_type;

        Integer r = A.rows(), c = A.cols();
        const raw::details::sparse_ccs<value_type>& Ad = A.rep();

        const Integer * Ad_c		= Ad.ptr_c();
        const Integer * Ad_r		= Ad.ptr_r();
        const value_type * Ad_x		= Ad.ptr_x();	

        Integer nzA = icast_c(Real(A.nnz()) + Real(c_in.rows())*Real(c_in.cols()));
        raw::details::sparse_ccs<value_type> d(A.get_type(),r, c, nzA);

        Integer * d_c				= d.ptr_c();
        Integer * d_r				= d.ptr_r();
        value_type * d_x			= d.ptr_x();

        column_iterator c_it(c_in);

        Integer nz					= 0;

        Integer rs = c_in.r_step, rf, rl;
        if (rs < 0)
        {
            rs = -rs;
            rf = c_in.r_end;
            rl = c_in.r_start;
        }
        else
        {
            rf = c_in.r_start;		
            rl = c_in.r_end;
        };
        --rf;
        --rl;
        
        for (Integer j = 0; j < c; ++j)
        {
            d_c[j]					= nz;

            if (!c_it.valid() || c_it.get() > j+1)
            {
                for (Integer k = Ad_c[j]; k < Ad_c[j+1]; ++k)
                {
                    d_r[nz] = Ad_r[k];
                    mrd::assign_helper(d_x[nz],Ad_x[k]);
                    ++nz;
                };
                continue;
            };

            while (c_it.valid() && c_it.get() == j+1)
                c_it.next();

            Integer la			= Ad_c[j+1];
            Integer i			= Ad_c[j];

            Integer k = rf;
            while(i < la && k <= rl)
            {
                Integer p		= Ad_r[i];

                if (p < k)
                {
                    d_r[nz]		= Ad_r[i];
                    mrd::assign_helper(d_x[nz],Ad_x[i]);
                    ++nz;
                    ++i;
                }
                else if (p > k)
                {
                    d_r[nz]		= k;
                    mrd::assign_helper(d_x[nz],val);
                    ++nz;
                    k			+= rs;
                }
                else
                {
                    d_r[nz]		= k;
                    mrd::assign_change_helper(d_x[nz],Ad_x[i],val);
                    ++nz;
                    k			+= rs;
                    ++i;
                };
            };

            while(i < la)
            {
                d_r[nz]			= Ad_r[i];
                mrd::assign_helper(d_x[nz],Ad_x[i]);
                ++nz;
                ++i;
            };
            
            while(k <= rl)
            {
                d_r[nz]			= k;
                mrd::assign_helper(d_x[nz],val);
                ++nz;
                k				+= rs;
            };
        };
        d_c[c] = nz;

        ret = Matrix(raw::sparse_matrix_base<value_type>(d), false);
    };

    static void change_cols(Matrix& ret, const SM& A,const md::colon_info& c_in, 
                            const value_type& val)
    {
        using value_type = typename SM::value_type;

        Integer r = A.rows(), c = A.cols();
        const raw::details::sparse_ccs<value_type>& Ad = A.rep();

        const Integer * Ad_c		= Ad.ptr_c();
        const Integer * Ad_r		= Ad.ptr_r();
        const value_type * Ad_x		= Ad.ptr_x();	

        Integer nzA = icast_c(Real(A.nnz()) + Real(c_in.rows())*Real(c_in.cols()));
        raw::details::sparse_ccs<value_type> d(A.get_type(),r, c, nzA);

        Integer * d_c				= d.ptr_c();
        Integer * d_r				= d.ptr_r();
        value_type * d_x			= d.ptr_x();

        column_iterator c_it(c_in);

        Integer nz					= 0;
        
        for (Integer j = 0; j < c; ++j)
        {
            d_c[j]					= nz;

            if (!c_it.valid() || c_it.get() > j+1)
            {
                for (Integer k = Ad_c[j]; k < Ad_c[j+1]; ++k)
                {
                    d_r[nz] = Ad_r[k];
                    mrd::assign_helper(d_x[nz],Ad_x[k]);
                    ++nz;
                };
                continue;
            };

            while (c_it.valid() && c_it.get() == j+1)
                c_it.next();

            for (Integer i = 0; i < r; ++i)
            {
                d_r[nz] = i;
                mrd::assign_helper(d_x[nz],val);
                ++nz;
            };
        };
        d_c[c] = nz;

        ret = Matrix(raw::sparse_matrix_base<value_type>(d), false);
    };
};

template<class SM>
void change_entries_functor<SM>::eval(Matrix& ret, SM& A, const md::colon_info& ci, 
                                    const value_type& val)
{
    if (ci.rows() == 0 || ci.cols() == 0)
    {
        ret = Matrix(A, false);
        return;
    };

    if (ci.r_flag == 0)
        return change_entries_impl<SM>::eval_02(ret, A, ci, val);
    else
        return change_entries_impl<SM>::eval_12(ret, A, ci, val);
};

template<class SM>
void add_entries_functor<SM>::eval(Matrix& ret, SM& A, const md::colon_info& ci)
{
    if (ci.rows() == 0 || ci.cols() == 0)
    {
        ret = Matrix(A, false);
        return;
    };

    if (ci.r_flag == 0)
        return add_entries_impl<SM>::eval_02(ret, A, ci);
    else
        return add_entries_impl<SM>::eval_12(ret, A, ci);
};

template<class SM>
void change_entries_functor_2<SM>::eval(Matrix& ret, SM& A,const md::colon_info& ci, 
                                    const value_type& val)
{
    if (ci.rows() == 0)
    {
        ret = Matrix(A, false);
        return;
    };

    if (ci.r_flag == 0)
        return change_entries_impl<SM>::eval_0(ret, A, ci, val);
    else
        return change_entries_impl<SM>::eval_1(ret, A, ci, val);
};

template<class SM>
void add_entries_functor_2<SM>::eval(Matrix& ret, SM& A,const md::colon_info& ci)
{
    if (ci.rows() == 0)
    {
        ret = Matrix(A, false);
        return;
    };

    if (ci.r_flag == 0)
        return add_entries_impl<SM>::eval_0(ret, A, ci);
    else
        return add_entries_impl<SM>::eval_1(ret, A, ci);
};

template<class V>
void sparse_change_diag_functor<V>::eval(Matrix& ret, SM& mat, Integer d, const DM& B)
{
    Integer r = mat.rows();
    Integer c = mat.cols();

    error::check_diag(d,r,c);

    Integer s, r1, c1;
    if (d >= 0)
    {
        s   = (r + d >= c) ? c - d : r;
        r1  = 0;
        c1  = d;
    }
    else
    {
        s   = (r + d >= c) ? c : r + d;
        c1  = 0;
        r1  = -d;
    }

    error::check_assign_1(s,B.size(),1);

    if (s == 0)
    {
        ret = Matrix(mat, false);
        return;
    };

    raw::details::sparse_ccs<V>& rep = mat.rep();

    const Integer* d_c      = rep.ptr_c();
    const Integer* d_r      = rep.ptr_r();    
    V * d_x		            = rep.ptr_x();
    Integer off             = rep.offset();
    DM B_tmp                = B.make_explicit();
    const V* ptr_val        = B_tmp.ptr();

    bool modif_needed       = false;

    for (Integer i = 0; i < s; ++i, ++c1, ++r1, ++ptr_val)
    {        
        const V& val = *ptr_val;

        Integer pos;
        if (rep.has_element(r1,c1,pos))
        {
            mrd::assign_helper(d_x[pos],val);
        }
        else
        {
            modif_needed = true;
            break;
        };
    };

    value_struct_class vt = value_struct_class::vc_general;

    if (modif_needed == false)
    {
        mat.set_struct(md::predefined_struct::get_set_diag(mat.get_struct(), d, vt,
                            is_real_matrix(mat), r == c));

        ret     = Matrix(mat, false);
        return;
    };

    raw::details::sparse_ccs<V> out(mat.get_type(),r, c, mat.nnz() + s);

    Integer* out_r          = out.ptr_r();
    Integer* out_c          = out.ptr_c();
    V * out_x		        = out.ptr_x();

    d_r                     = rep.ptr_r() + off;
    d_x		                = rep.ptr_x() + off;

    Integer nz = d_c[c1] - off;
    for (Integer i = 0; i < c1; ++i)
        out_c[i] = d_c[i] - off;

    for (Integer i = 0; i < nz; ++i)
    {
        out_r[i] = d_r[i];
        mrd::assign_helper(out_x[i],d_x[i]);
    };

    d_r                     = rep.ptr_r();
    d_x		                = rep.ptr_x();

    for (Integer i = c1; i < c; ++i, ++r1)
    {
        out_c[i]			= nz;

        Integer k           = d_c[i];
        Integer k2          = d_c[i+1];
        for (; k < k2; ++k)
        {
            Integer row = d_r[k];
            if (row >= r1)
            {
                if (row == r1)
                    ++k;

                break;
            };

            out_r[nz] = row;
            mrd::assign_helper(out_x[nz],d_x[k]);
            ++nz;
        };

        if (r1 < r)
        {
            const V& val = *ptr_val;
            ++ptr_val;

            out_r[nz] = r1;
            mrd::assign_helper(out_x[nz],val);
            ++nz;
        }

        for (; k < k2; ++k)
        {
            out_r[nz] = d_r[k];
            mrd::assign_helper(out_x[nz],d_x[k]);
            ++nz;
        };
    };

    out_c[c] = nz;

    out.get_struct() = md::predefined_struct::get_set_diag(mat.get_struct(),d,vt,
                                is_real_matrix(out), r == c);

    ret = Matrix(raw::sparse_matrix_base<V>(out), false);
    return;
};

template<class V>
void sparse_change_diag_functor<V>::eval(Matrix& ret, SM& mat, Integer d, const V& val)
{
    Integer r = mat.rows();
    Integer c = mat.cols();

    error::check_diag(d,r,c);

    Integer s, r1, c1;
    if (d >= 0)
    {
        s = (r + d >= c) ? c - d : r;
        r1 = 0;
        c1 = d;
    }
    else
    {
        s = (r + d >= c) ? c : r + d;
        c1 = 0;
        r1 = -d;
    }

    if (s == 0)
    {
        ret = Matrix(mat, false);
        return;
    };

    raw::details::sparse_ccs<V>& rep = mat.rep();

    const Integer* d_c      = rep.ptr_c();
    const Integer* d_r      = rep.ptr_r();
    V * d_x		            = rep.ptr_x();
    Integer off             = rep.offset();

    bool modif_needed       = false;

    for (Integer i = 0; i < s; ++i, ++c1, ++r1)
    {
        Integer pos;
        if (rep.has_element(r1,c1,pos))
        {
            mrd::assign_helper(d_x[pos],val);
        }
        else
        {
            modif_needed = true;
            break;
        };
    };

    bool tz = md::has_trivial_assignment<SM,V>::eval(mat,val);
    value_struct_class vt = md::predefined_struct::get_value_type(val,tz);

    if (modif_needed == false)
    {
        mat.set_struct(md::predefined_struct::get_set_diag(mat.get_struct(),d,vt,
                                is_real_matrix(mat), r == c));

        ret = Matrix(mat, false);
        return;
    };

    raw::details::sparse_ccs<V> out(mat.get_type(),r, c, mat.nnz() + s);

    Integer* out_r          = out.ptr_r();
    Integer* out_c          = out.ptr_c();
    V * out_x		        = out.ptr_x();

    d_r                     = rep.ptr_r() + off;
    d_x		                = rep.ptr_x() + off;

    Integer nz = d_c[c1] - off;
    for (Integer i = 0; i < c1; ++i)
        out_c[i] = d_c[i] - off;

    for (Integer i = 0; i < nz; ++i)
    {
        out_r[i] = d_r[i];
        mrd::assign_helper(out_x[i],d_x[i]);
    };

    d_r                     = rep.ptr_r();
    d_x		                = rep.ptr_x();

    for (Integer i = c1; i < c; ++i, ++r1)
    {
        out_c[i]			= nz;

        Integer k           = d_c[i];
        Integer k2          = d_c[i+1];
        for (; k < k2; ++k)
        {
            Integer row = d_r[k];
            if (row >= r1)
            {
                if (row == r1)
                    ++k;

                break;
            };

            out_r[nz] = row;
            mrd::assign_helper(out_x[nz],d_x[k]);
            ++nz;
        };

        if (r1 < r)
        {
            out_r[nz] = r1;
            mrd::assign_helper(out_x[nz],val);
            ++nz;
        };

        for (; k < k2; ++k)
        {
            out_r[nz] = d_r[k];
            mrd::assign_helper(out_x[nz],d_x[k]);
            ++nz;
        };
    };

    out_c[c] = nz;
    out.get_struct() = md::predefined_struct::get_set_diag(mat.get_struct(),d,vt,
                                    is_real_matrix(out), r == c);

    ret = Matrix(raw::sparse_matrix_base<V>(out), false);
};

template<class V>
void sparse_add_diag_functor<V>::eval(Matrix& ret, SM& mat, Integer d)
{
    Integer r = mat.rows();
    Integer c = mat.cols();

    error::check_diag(d,r,c);

    Integer s, r1, c1;
    if (d >= 0)
    {
        s = (r + d >= c) ? c - d : r;
        r1 = 0;
        c1 = d;
    }
    else
    {
        s = (r + d >= c) ? c : r + d;
        c1 = 0;
        r1 = -d;
    }

    if (s == 0)
    {
        ret = Matrix(mat, false);
        return;
    };

    raw::details::sparse_ccs<V>& rep = mat.rep();

    const Integer* d_c      = rep.ptr_c();
    const Integer* d_r      = rep.ptr_r();
    V * d_x		            = rep.ptr_x();
    Integer off             = rep.offset();

    bool add_needed         = false;

    for (Integer i = 0; i < s; ++i, ++c1, ++r1)
    {
        Integer pos;
        if (rep.has_element(r1,c1,pos) == false)
        {
            add_needed      = true;
            break;
        };
    };

    if (add_needed == false)
    {
        ret = Matrix(mat, false);
        return;
    };

    raw::details::sparse_ccs<V> out(mat.get_type(),r, c, mat.nnz() + s);

    Integer* out_r          = out.ptr_r();
    Integer* out_c          = out.ptr_c();
    V * out_x		        = out.ptr_x();

    d_r                     = rep.ptr_r() + off;
    d_x		                = rep.ptr_x() + off;
    const V& Z              = md::default_value<V>(mat.get_type());

    Integer nz = d_c[c1] - off;
    for (Integer i = 0; i < c1; ++i)
        out_c[i] = d_c[i] - off;

    for (Integer i = 0; i < nz; ++i)
    {
        out_r[i] = d_r[i];
        mrd::assign_helper(out_x[i],d_x[i]);
    };

    d_r                     = rep.ptr_r();
    d_x		                = rep.ptr_x();

    for (Integer i = c1; i < c; ++i, ++r1)
    {
        out_c[i]			= nz;

        Integer k           = d_c[i];
        Integer k2          = d_c[i+1];
        Integer row         = -1;

        for (; k < k2; ++k)
        {
            row         = d_r[k];

            if (row >= r1)
                break;

            out_r[nz]   = row;
            mrd::assign_helper(out_x[nz],d_x[k]);
            ++nz;
        };

        if (r1 < r && row != r1)
        {
            out_r[nz] = r1;
            mrd::assign_helper(out_x[nz], Z);
            ++nz;
        };
        
        for (; k < k2; ++k)
        {
            out_r[nz] = d_r[k];
            mrd::assign_helper(out_x[nz],d_x[k]);
            ++nz;
        };
    };

    out_c[c] = nz;
    out.get_struct() = mat.get_struct();

    ret = Matrix(raw::sparse_matrix_base<V>(out), false);
};

template<class V>
void sparse_drop_diag_functor<V>::eval(Matrix& ret, SM& mat, Integer d, Real tol0)
{
    Integer r   = mat.rows();
    Integer c   = mat.cols();
    Real tol    = tol0;

    error::check_diag(d,r,c);

    Integer s, r1, c1;
    if (d >= 0)
    {
        s   = (r + d >= c) ? c - d : r;
        r1  = 0;
        c1  = d;
    }
    else
    {
        s   = (r + d >= c) ? c : r + d;
        c1  = 0;
        r1  = -d;
    }

    if (s == 0 || tol < 0)
    {
        ret     = Matrix(mat, false);
        return;
    };

    raw::details::sparse_ccs<V>& rep = mat.rep();

    Integer* d_c            = rep.ptr_c();
    Integer* d_r            = rep.ptr_r();
    V * d_x		            = rep.ptr_x();
    Integer off             = rep.offset();

    bool modif_needed       = false;

    for (Integer i = 0; i < s; ++i, ++c1, ++r1)
    {
        Integer pos;

        if (rep.has_element(r1,c1,pos) == false)
            continue;

        const V& tmp        = d_x[pos];

        if (mrd::abs_helper<V>::eval(tmp) <= tol)
        {
            d_r[pos]        = -1;
            modif_needed    = true;
        };
    };

    if (modif_needed == false)
    {
        ret = Matrix(mat, false);
        return;
    };

    for (Integer i = 0, pos = off; i < c; ++i)
    {
        for (Integer j = d_c[i]; j < d_c[i+1]; ++j)
        {
            if (d_r[j] == -1)
                continue;

            if (j == pos)
            {
                ++pos;
                continue;
            };

            d_r[pos] = d_r[j];
            mrd::reset_helper(d_x[pos],d_x[j]);
            ++pos;
        };

        d_c[i] = pos;
    };

    for (Integer i = c; i > 0; --i)
        d_c[i] = d_c[i-1];

    d_c[0] = off;

    rep.add_memory(-1);

    if (tol > 0)
        mat.get_struct().reset();

    ret = Matrix(mat, false);
    return;
};

template<class V>
void sparse_drop_diag_functor<V>::eval_zero(Matrix& ret, SM& mat, Integer d)
{
    Integer r   = mat.rows();
    Integer c   = mat.cols();

    error::check_diag(d,r,c);

    Integer s, r1, c1;
    if (d >= 0)
    {
        s   = (r + d >= c) ? c - d : r;
        r1  = 0;
        c1  = d;
    }
    else
    {
        s   = (r + d >= c) ? c : r + d;
        c1  = 0;
        r1  = -d;
    }

    if (s == 0)
    {
        ret = Matrix(mat, false);
        return;
    };

    raw::details::sparse_ccs<V>& rep = mat.rep();

    Integer* d_c            = rep.ptr_c();
    Integer* d_r            = rep.ptr_r();
    V * d_x		            = rep.ptr_x();
    Integer off             = rep.offset();

    bool modif_needed       = false;

    for (Integer i = 0; i < s; ++i, ++c1, ++r1)
    {
        Integer pos;

        if (rep.has_element(r1,c1,pos) == false)
            continue;

        d_r[pos]            = -1;
        modif_needed        = true;
    };

    if (modif_needed == false)
    {
        ret = Matrix(mat, false);
        return;
    };

    for (Integer i = 0, pos = off; i < c; ++i)
    {
        for (Integer j = d_c[i]; j < d_c[i+1]; ++j)
        {
            if (d_r[j] == -1)
                continue;

            if (j == pos)
            {
                ++pos;
                continue;
            };

            d_r[pos] = d_r[j];
            mrd::reset_helper(d_x[pos],d_x[j]);
            ++pos;
        };

        d_c[i] = pos;
    };

    for (Integer i = c; i > 0; --i)
        d_c[i] = d_c[i-1];

    d_c[0] = off;

    rep.add_memory(-1);
    mat.get_struct().reset();

    ret     = Matrix(mat, false);
    return;
};

};};};

template struct matcl::algorithm::details::change_entries_functor<matcl::raw::integer_sparse>;
template struct matcl::algorithm::details::change_entries_functor<matcl::raw::real_sparse>;
template struct matcl::algorithm::details::change_entries_functor<matcl::raw::float_sparse>;
template struct matcl::algorithm::details::change_entries_functor<matcl::raw::complex_sparse>;
template struct matcl::algorithm::details::change_entries_functor<matcl::raw::float_complex_sparse>;
template struct matcl::algorithm::details::change_entries_functor<matcl::raw::object_sparse>;

template struct matcl::algorithm::details::change_entries_functor_2<matcl::raw::integer_sparse>;
template struct matcl::algorithm::details::change_entries_functor_2<matcl::raw::real_sparse>;
template struct matcl::algorithm::details::change_entries_functor_2<matcl::raw::float_sparse>;
template struct matcl::algorithm::details::change_entries_functor_2<matcl::raw::complex_sparse>;
template struct matcl::algorithm::details::change_entries_functor_2<matcl::raw::float_complex_sparse>;
template struct matcl::algorithm::details::change_entries_functor_2<matcl::raw::object_sparse>;

template struct matcl::algorithm::details::add_entries_functor<matcl::raw::integer_sparse>;
template struct matcl::algorithm::details::add_entries_functor<matcl::raw::real_sparse>;
template struct matcl::algorithm::details::add_entries_functor<matcl::raw::float_sparse>;
template struct matcl::algorithm::details::add_entries_functor<matcl::raw::complex_sparse>;
template struct matcl::algorithm::details::add_entries_functor<matcl::raw::float_complex_sparse>;
template struct matcl::algorithm::details::add_entries_functor<matcl::raw::object_sparse>;

template struct matcl::algorithm::details::add_entries_functor_2<matcl::raw::integer_sparse>;
template struct matcl::algorithm::details::add_entries_functor_2<matcl::raw::real_sparse>;
template struct matcl::algorithm::details::add_entries_functor_2<matcl::raw::float_sparse>;
template struct matcl::algorithm::details::add_entries_functor_2<matcl::raw::complex_sparse>;
template struct matcl::algorithm::details::add_entries_functor_2<matcl::raw::float_complex_sparse>;
template struct matcl::algorithm::details::add_entries_functor_2<matcl::raw::object_sparse>;

template struct matcl::algorithm::details::sparse_change_diag_functor<matcl::Integer>;
template struct matcl::algorithm::details::sparse_change_diag_functor<matcl::Real>;
template struct matcl::algorithm::details::sparse_change_diag_functor<matcl::Float>;
template struct matcl::algorithm::details::sparse_change_diag_functor<matcl::Complex>;
template struct matcl::algorithm::details::sparse_change_diag_functor<matcl::Float_complex>;
template struct matcl::algorithm::details::sparse_change_diag_functor<matcl::Object>;

template struct matcl::algorithm::details::sparse_drop_diag_functor<matcl::Integer>;
template struct matcl::algorithm::details::sparse_drop_diag_functor<matcl::Real>;
template struct matcl::algorithm::details::sparse_drop_diag_functor<matcl::Float>;
template struct matcl::algorithm::details::sparse_drop_diag_functor<matcl::Complex>;
template struct matcl::algorithm::details::sparse_drop_diag_functor<matcl::Float_complex>;
template struct matcl::algorithm::details::sparse_drop_diag_functor<matcl::Object>;

template struct matcl::algorithm::details::sparse_add_diag_functor<matcl::Integer>;
template struct matcl::algorithm::details::sparse_add_diag_functor<matcl::Real>;
template struct matcl::algorithm::details::sparse_add_diag_functor<matcl::Float>;
template struct matcl::algorithm::details::sparse_add_diag_functor<matcl::Complex>;
template struct matcl::algorithm::details::sparse_add_diag_functor<matcl::Float_complex>;
template struct matcl::algorithm::details::sparse_add_diag_functor<matcl::Object>;
