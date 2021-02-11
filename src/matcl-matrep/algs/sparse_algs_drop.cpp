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

#include "matcl-matrep/algs/sparse_algs.h"
#include "matcl-matrep/algs/sparse_algs_utils.h"
#include "matcl-internals/base/sort.h"
#include "matcl-matrep/matrix/matrix.h"
#include "matcl-internals/container/mat_s.h"
#include "matcl-core/utils/workspace.h"
#include "matcl-scalar/details/scalfunc_helpers.h"
#include "matcl-matrep/lib_functions/func_binary.h"

namespace matcl { namespace algorithm { namespace details
{

namespace mrd = matcl::raw::details;

class double_colon_it
{
    private:
        raw::integer_dense  m_r;
        raw::integer_dense  m_c;
        Integer             m_pos;
        Integer             m_length;
        Integer             m_row_end;

    public:
        double_colon_it(const raw::integer_dense& ri, const raw::integer_dense& ci)
            :m_r(ri.copy()), m_c(ci.copy()), m_pos(-1), m_length(ri.length())
            ,m_row_end(0)
        {
            sort_rows_cols(m_r.ptr(), m_c.ptr(), m_length);
        };

        Integer next_column()
        {
            ++m_pos;

            if (m_pos == m_length)
                return -1;

            m_row_end = calc_row_end();

            return m_c.ptr()[m_pos] - 1;
        };

        Integer next_row()
        {
            if (m_pos < m_row_end)
                return m_r.ptr()[m_pos++] - 1;
            else
                return -1;
        };

    private:

        Integer calc_row_end() const
        {
            Integer c = m_c.ptr()[m_pos];
            Integer i = m_pos;

            for (; i < m_length; ++i)
            {
                Integer cl = m_c.ptr()[i];

                if (c != cl)
                    break;
            };

            return i + 1;
        };
};

class column_iterator_drop
{
    private:
        const matcl::details::colon_info& ci;
        Integer pos;

    public:
        column_iterator_drop(const matcl::details::colon_info& ci) : ci(ci) 
        {
            if (ci.c_flag == 0)
                pos = 0;
            else
                pos = ci.c_start - 1;
        };

        Integer size()
        {
            return ci.cols();
        };

        Integer get()
        {
            if (ci.c_flag == 0)
            {
                return ci.get_cim_2().ptr()[pos++];
            }
            else
            {
                Integer tmp = pos;
                pos += ci.c_step;
                return tmp + 1;
            };
        };

    private:
        column_iterator_drop(const column_iterator_drop&) = delete;
        column_iterator_drop& operator=(const column_iterator_drop&) = delete;
};

template<class SM>
struct drop_entries_impl
{
    static SM eval_20(SM& in,const matcl::details::colon_info& ci, Real tol0)
    {
        using value_type = typename SM::value_type;

        Real tol    = tol0;
        Integer r   = in.rows();
        Integer c   = in.cols();
        raw::details::sparse_ccs<value_type>& Ad = in.rep();

        bool any_deleted	= false;

        Integer * Ad_c		= Ad.ptr_c();
        Integer * Ad_r		= Ad.ptr_r();
        value_type * Ad_x	= Ad.ptr_x();		

        Matrix ri(ci.get_rim_2(),false);
        Integer n = ci.get_rim_2().size();

        column_iterator_drop c_it(ci);

        if (imult(n,c) < r)
        {
            raw::integer_dense ri_raw = ri.impl_unique<raw::integer_dense>();
            const Integer* ptr_ri_raw = ri_raw.ptr();

            utils::sort_q(ri_raw.ptr(),n);

            Integer size = c_it.size();
            for (Integer i = 0; i < size; ++i)
            {
                Integer col     = c_it.get() - 1;
                Integer nz_row  = Ad_c[col+1] - Ad_c[col];

                if (nz_row == 0)
                    continue;

                Integer la  = Ad_c[col+1];
                Integer j   = Ad_c[col];

                Integer k   = 0;
                while(j < la && k < n)
                {
                    Integer p = Ad_r[j];

                    if (p < ptr_ri_raw[k]-1)
                    {
                        ++j;
                    }
                    else if (p > ptr_ri_raw[k]-1)
                    {
                        ++k;
                    }
                    else
                    {
                        const value_type& tmp   = Ad_x[j];
                        
                        if (mrd::abs_helper<value_type>::eval(tmp) <= tol)
                        {
                            Ad_r[j]         = -1;
                            any_deleted     = true;
                        };

                        ++j;
                        ++k;
                    };
                };
            };
        }
        else
        {
            raw::integer_dense ri_raw = ri.impl_unique<raw::integer_dense>();
            const Integer* ptr_ri_raw = ri_raw.ptr();

            matcl::pod_workspace<Integer> v_work_ind;
            bool work_initialized = false;

            Integer r_min = 0, r_max = 0;

            Integer cit_size    = c_it.size();
            Integer ri_size     = ri_raw.size();

            for (Integer i = 1; i <= cit_size; ++i)
            {
                Integer col     = c_it.get() - 1;
                Integer nz_row  = Ad_c[col+1] - Ad_c[col];

                if (nz_row > 0 && work_initialized == false)
                {
                    init_row_selector(ri_raw,0,v_work_ind,r_min,r_max);                    

                    for (Integer k = 0; k < ri_size; ++k)
                        v_work_ind[ptr_ri_raw[k]-1-r_min] = 1;

                    work_initialized = true;
                };

                for (Integer j = Ad_c[col]; j < Ad_c[col+1]; ++j)
                {
                    Integer row = Ad_r[j];

                    if (row < r_min || row > r_max)
                        continue;

                    if (v_work_ind[row-r_min] == 1)
                    {
                        const value_type& tmp   = Ad_x[j];

                        if (mrd::abs_helper<value_type>::eval(tmp) <= tol)
                        {
                            Ad_r[j]     = -1;
                            any_deleted = true;
                        };
                    };
                };
            };
        };

        if (any_deleted == false)
            return in;

        for (Integer i = 0, pos = Ad.offset(); i < c; ++i)
        {
            for (Integer j = Ad_c[i]; j < Ad_c[i+1]; ++j)
            {
                if (Ad_r[j] == -1)
                    continue;

                if (j == pos)
                {
                    ++pos;
                    continue;
                };

                Ad_r[pos] = Ad_r[j];
                mrd::reset_helper(Ad_x[pos],Ad_x[j]);
                ++pos;
            };

            Ad_c[i] = pos;
        };

        for (Integer i = c; i >0; --i)
            Ad_c[i] = Ad_c[i-1];

        Ad_c[0] = Ad.offset();

        Ad.add_memory(-1);
        in.get_struct().reset();
        return in;
    };

    static SM eval_20_zero(SM& in,const matcl::details::colon_info& ci)
    {
        using value_type = typename SM::value_type;

        Integer r   = in.rows();
        Integer c   = in.cols();
        raw::details::sparse_ccs<value_type>& Ad = in.rep();

        Integer * Ad_c		= Ad.ptr_c();
        Integer * Ad_r		= Ad.ptr_r();
        value_type * Ad_x	= Ad.ptr_x();		

        Matrix ri(ci.get_rim_2(),false);
        Integer n           = ci.get_rim_2().size();
        bool any_modif      = false;
        const value_type& Z = md::default_value<value_type>(in.get_type());

        column_iterator_drop c_it(ci);

        if (imult(n,c) < r)
        {
            raw::integer_dense ri_raw = ri.impl_unique<raw::integer_dense>();
            const Integer* ptr_ri_raw = ri_raw.ptr();

            utils::sort_q(ri_raw.ptr(),n);

            Integer size = c_it.size();
            for (Integer i = 0; i < size; ++i)
            {
                Integer col     = c_it.get() - 1;
                Integer nz_row  = Ad_c[col+1] - Ad_c[col];

                if (nz_row == 0)
                    continue;

                Integer la  = Ad_c[col+1];
                Integer j   = Ad_c[col];

                Integer k   = 0;

                while(j < la && k < n)
                {
                    Integer p = Ad_r[j];

                    if (p < ptr_ri_raw[k]-1)
                    {
                        ++j;
                    }
                    else if (p > ptr_ri_raw[k]-1)
                    {
                        ++k;
                    }
                    else
                    {
                        Ad_x[j]         = Z;
                        any_modif       = true;
                        ++j;
                        ++k;
                    };
                };
            };
        }
        else
        {
            raw::integer_dense ri_raw = ri.impl_unique<raw::integer_dense>();
            const Integer* ptr_ri_raw = ri_raw.ptr();

            matcl::pod_workspace<Integer> v_work_ind;
            bool work_initialized = false;

            Integer r_min = 0, r_max = 0;

            Integer cit_size    = c_it.size();
            Integer ri_size     = ri_raw.size();

            for (Integer i = 1; i <= cit_size; ++i)
            {
                Integer col     = c_it.get() - 1;
                Integer nz_row  = Ad_c[col+1] - Ad_c[col];

                if (nz_row > 0 && work_initialized == false)
                {
                    init_row_selector(ri_raw,0,v_work_ind,r_min,r_max);                    

                    for (Integer k = 0; k < ri_size; ++k)
                        v_work_ind[ptr_ri_raw[k]-1-r_min] = 1;

                    work_initialized = true;
                };

                for (Integer j = Ad_c[col]; j < Ad_c[col+1]; ++j)
                {
                    Integer row = Ad_r[j];

                    if (row < r_min || row > r_max)
                        continue;

                    if (v_work_ind[row-r_min] == 1)
                    {
                        Ad_x[j]     = Z;
                        any_modif   = true;
                    };
                };
            };
        };

        if (any_modif == false)
            return in;
        
        in.get_struct().reset();
        return in;
    };

    static SM eval_21(SM& in,const matcl::details::colon_info& ci, Real tol0)
    {
        using value_type = typename SM::value_type;

        Integer c = in.cols();
        raw::details::sparse_ccs<value_type>& Ad = in.rep();

        bool any_deleted	= false;
        Real tol            = tol0;

        Integer * Ad_c		= Ad.ptr_c();
        Integer * Ad_r		= Ad.ptr_r();
        value_type * Ad_x	= Ad.ptr_x();		

        Integer rs = ci.r_step, rf, rl;

        if (rs < 0)
        {
            rs = -rs;
            rf = ci.r_end;
            rl = ci.r_start;
        }
        else
        {
            rf = ci.r_start;		
            rl = ci.r_end;
        };
        --rf;
        --rl;

        column_iterator_drop c_it(ci);

        Integer c_it_size = c_it.size();
        for (Integer i = 1; i <= c_it_size; ++i)
        {
            Integer col = c_it.get() - 1;

            for (Integer j = Ad_c[col]; j < Ad_c[col+1]; ++j)
            {
                Integer row = Ad_r[j];
                
                if (row < rf || row > rl)
                    continue;

                if ((row-rf)%rs == 0)
                {
                    const value_type& tmp   = Ad_x[j];
                    if (mrd::abs_helper<value_type>::eval(tmp) <= tol)
                    {
                        Ad_r[j] = -1;
                        any_deleted = true;
                    };
                };
            };
        };

        if (any_deleted == false)
            return in;

        for (Integer i = 0, pos = Ad.offset(); i < c; ++i)
        {
            for (Integer j = Ad_c[i]; j < Ad_c[i+1]; ++j)
            {
                if (Ad_r[j] == -1)
                    continue;

                if (j == pos)
                {
                    ++pos;
                    continue;
                };

                Ad_r[pos] = Ad_r[j];
                mrd::reset_helper(Ad_x[pos],Ad_x[j]);
                ++pos;
            };

            Ad_c[i] = pos;
        };

        for (Integer i = c; i >0; --i)
            Ad_c[i] = Ad_c[i-1];

        Ad_c[0] = Ad.offset();

        Ad.add_memory(-1);
        in.get_struct().reset();
        return in;
    };

    static SM eval_21_zero(SM& in,const matcl::details::colon_info& ci)
    {
        using value_type = typename SM::value_type;

        raw::details::sparse_ccs<value_type>& Ad = in.rep();

        bool any_modif	    = false;

        Integer * Ad_c		= Ad.ptr_c();
        Integer * Ad_r		= Ad.ptr_r();
        value_type * Ad_x	= Ad.ptr_x();		

        Integer rs = ci.r_step, rf, rl;

        if (rs < 0)
        {
            rs = -rs;
            rf = ci.r_end;
            rl = ci.r_start;
        }
        else
        {
            rf = ci.r_start;		
            rl = ci.r_end;
        };

        --rf;
        --rl;

        column_iterator_drop c_it(ci);

        const value_type& Z = md::default_value<value_type>(in.get_type());

        Integer c_it_size = c_it.size();

        for (Integer i = 1; i <= c_it_size; ++i)
        {
            Integer col = c_it.get() - 1;

            for (Integer j = Ad_c[col]; j < Ad_c[col+1]; ++j)
            {
                Integer row = Ad_r[j];
                
                if (row < rf || row > rl)
                    continue;

                if ((row-rf)%rs == 0)
                {
                    Ad_x[j]     = Z;
                    any_modif   = true;
                };
            };
        };

        if (any_modif == false)
            return in;

        in.get_struct().reset();
        return in;
    };

    static SM eval_0(SM& out,const matcl::details::colon_info& ci, Real tol0)
    {
        if (ci.is_double_mat_colon() == true)
            return eval_0_dc(out, ci, tol0);

        using value_type = typename SM::value_type;
        Integer r   = out.rows();
        Integer c   = out.cols();
        Real tol    = tol0;

        Matrix ri(ci.get_rim_1(),false);
        raw::integer_dense ri_raw = ri.impl_unique<raw::integer_dense>();

        sort_type s_type = is_sorted(ri_raw);
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
            ri.make_unique();
            raw::integer_dense ri_tmp = ri.impl_unique<raw::integer_dense>();
            utils::sort_q(ri_tmp.ptr(),ri_tmp.size());

            ri_raw.assign_to_fresh(std::move(ri_tmp));
        };
        
        bool any_deleted	= false;

        raw::details::sparse_ccs<value_type>& Ad = out.rep();
        Integer * Ad_c		= Ad.ptr_c();
        Integer * Ad_r		= Ad.ptr_r();
        value_type * Ad_x	= Ad.ptr_x();		

        Integer step_ri;
        Integer pos, pos_row, pos_col;

        if (incr)
        {
            pos = 0;
            step_ri = 1;
        }
        else
        {
            pos = ri_raw.size()-1;
            step_ri = -1;
        };

        Integer or = ri_raw.size();
        const Integer* ptr_ri_raw = ri_raw.ptr();

        matcl::details::pos2ind(ptr_ri_raw[pos],r,pos_row,pos_col);

        for (Integer j = 0; j < c; ++j)
        {
            if (j < pos_col)
                continue;

            for (Integer i = Ad_c[j]; i < Ad_c[j+1];++i)
            {
                Integer p = Ad_r[i];

                if (p < pos_row)
                    continue;

                while (p > pos_row)
                {
                    pos+=step_ri;

                    if (pos > or-1 || pos < 0)
                        goto exit_flag;

                    matcl::details::pos2ind(ptr_ri_raw[pos],r,pos_row,pos_col);

                    if (j != pos_col)
                        goto exit_row_flag;
                };

                while (p == pos_row)
                {
                    const value_type& tmp   = Ad_x[i];

                    if (mrd::abs_helper<value_type>::eval(tmp) <= tol)
                    {
                        Ad_r[i] = -1;
                        any_deleted = true;
                    };

                    pos += step_ri;
                    if (pos > or-1 || pos < 0)
                        goto exit_flag;

                    matcl::details::pos2ind(ptr_ri_raw[pos],r,pos_row,pos_col);
                    if (j != pos_col)
                        goto exit_row_flag;
                };
            };

            while (j == pos_col)
            {
                pos+= step_ri;
                if (pos > or-1 || pos < 0)
                    goto exit_flag;

                matcl::details::pos2ind(ptr_ri_raw[pos],r,pos_row,pos_col);
            };

            exit_row_flag:
            continue;
        };

    exit_flag:

        if (any_deleted == false)
            return out;

        pos = Ad.offset();
        for (Integer i = 0; i < c; ++i)
        {
            for (Integer j = Ad_c[i]; j < Ad_c[i+1]; ++j)
            {
                if (Ad_r[j] == -1)
                    continue;

                if (j == pos)
                {
                    ++pos;
                    continue;
                };
                Ad_r[pos] = Ad_r[j];
                mrd::reset_helper(Ad_x[pos],Ad_x[j]);
                ++pos;
            };
            Ad_c[i] = pos;
        };

        for (Integer i = c; i >0; --i)
            Ad_c[i] = Ad_c[i-1];

        Ad_c[0] = Ad.offset();

        Ad.add_memory(-1);
        out.get_struct().reset();
        return out;
    };

    static SM eval_0_dc(SM& out,const matcl::details::colon_info& ci, Real tol0)
    {
        using value_type = typename SM::value_type;
        Integer c   = out.cols();
        Real tol    = tol0;
        
        double_colon_it it(ci.get_rim_r(), ci.get_rim_c());

        bool any_deleted	= false;

        raw::details::sparse_ccs<value_type>& Ad = out.rep();
        Integer * Ad_c		= Ad.ptr_c();
        Integer * Ad_r		= Ad.ptr_r();
        value_type * Ad_x	= Ad.ptr_x();		

        Integer pos_row, pos_col;

        for (;;)
        {
            pos_col         = it.next_column();

            if (pos_col < 0)
                break;

            pos_row         = it.next_row();

            for (Integer i = Ad_c[pos_col]; i < Ad_c[pos_col+1];++i)
            {
                Integer p = Ad_r[i];

                if (p < pos_row)
                    continue;

                while (p > pos_row)
                {
                    pos_row = it.next_row();

                    if (pos_row < 0)
                        goto exit_row_flag;
                };

                while (p == pos_row)
                {
                    const value_type& tmp   = Ad_x[i];

                    if (mrd::abs_helper<value_type>::eval(tmp) <= tol)
                    {
                        Ad_r[i]     = -1;
                        any_deleted = true;
                    };

                    pos_row = it.next_row();

                    if (pos_row < 0)
                        goto exit_row_flag;
                };
            };

            exit_row_flag:
            continue;
        };

        if (any_deleted == false)
            return out;

        Integer pos = Ad.offset();

        for (Integer i = 0; i < c; ++i)
        {
            for (Integer j = Ad_c[i]; j < Ad_c[i+1]; ++j)
            {
                if (Ad_r[j] == -1)
                    continue;

                if (j == pos)
                {
                    ++pos;
                    continue;
                };

                Ad_r[pos] = Ad_r[j];
                mrd::reset_helper(Ad_x[pos],Ad_x[j]);
                ++pos;
            };
            Ad_c[i] = pos;
        };

        for (Integer i = c; i >0; --i)
            Ad_c[i] = Ad_c[i-1];

        Ad_c[0] = Ad.offset();

        Ad.add_memory(-1);
        out.get_struct().reset();
        return out;
    };

    static SM eval_0_zero_dc(SM& out,const matcl::details::colon_info& ci)
    {
        using value_type = typename SM::value_type;

        double_colon_it it(ci.get_rim_r(), ci.get_rim_c());
        
        bool any_modif	    = false;

        raw::details::sparse_ccs<value_type>& Ad = out.rep();
        Integer * Ad_c		= Ad.ptr_c();
        Integer * Ad_r		= Ad.ptr_r();
        value_type * Ad_x	= Ad.ptr_x();		

        const value_type& Z = md::default_value<value_type>(out.get_type());

        Integer pos_row, pos_col;

        for (;;)
        {
            pos_col         = it.next_column();

            if (pos_col < 0)
                break;

            pos_row         = it.next_row();

            for (Integer i = Ad_c[pos_col]; i < Ad_c[pos_col+1];++i)
            {
                Integer p = Ad_r[i];

                if (p < pos_row)
                    continue;

                while (p > pos_row)
                {
                    pos_row = it.next_row();

                    if (pos_row < 0)
                        goto exit_row_flag;
                };

                while (p == pos_row)
                {
                    Ad_x[i]         = Z;
                    any_modif       = true;

                    pos_row = it.next_row();

                    if (pos_row < 0)
                        goto exit_row_flag;
                };
            };

          exit_row_flag:
            continue;
        };

        if (any_modif == false)
            return out;

        out.get_struct().reset();
        return out;
    };

    static SM eval_0_zero(SM& out,const matcl::details::colon_info& ci)
    {
        if (ci.is_double_mat_colon() == true)
            return eval_0_zero_dc(out, ci);

        using value_type = typename SM::value_type;
        Integer r   = out.rows();
        Integer c   = out.cols();

        Matrix ri(ci.get_rim_1(),false);
        raw::integer_dense ri_raw = ri.impl_unique<raw::integer_dense>();

        sort_type s_type    = is_sorted(ri_raw);
        bool incr           = true;

        if (s_type == sorted_increasing)
        {
        }
        else if (s_type == sorted_decreasing)
        {
            incr = false;
        } 
        else
        {
            ri.make_unique();
            raw::integer_dense ri_tmp = ri.impl_unique<raw::integer_dense>();
            utils::sort_q(ri_tmp.ptr(),ri_tmp.size());

            ri_raw.assign_to_fresh(std::move(ri_tmp));
        };
        
        bool any_modif	    = false;

        raw::details::sparse_ccs<value_type>& Ad = out.rep();
        Integer * Ad_c		= Ad.ptr_c();
        Integer * Ad_r		= Ad.ptr_r();
        value_type * Ad_x	= Ad.ptr_x();		

        Integer step_ri;
        Integer pos, pos_row, pos_col;

        if (incr)
        {
            pos = 0;
            step_ri = 1;
        }
        else
        {
            pos = ri_raw.size()-1;
            step_ri = -1;
        };

        Integer or                  = ri_raw.size();
        const Integer* ptr_ri_raw   = ri_raw.ptr();
        const value_type& Z         = md::default_value<value_type>(out.get_type());

        matcl::details::pos2ind(ptr_ri_raw[pos],r,pos_row,pos_col);

        for (Integer j = 0; j < c; ++j)
        {
            if (j < pos_col)
                continue;

            for (Integer i = Ad_c[j]; i < Ad_c[j+1];++i)
            {
                Integer p = Ad_r[i];

                if (p < pos_row)
                    continue;

                while (p > pos_row)
                {
                    pos+=step_ri;

                    if (pos > or-1 || pos < 0)
                        goto exit_flag;

                    matcl::details::pos2ind(ptr_ri_raw[pos],r,pos_row,pos_col);

                    if (j != pos_col)
                        goto exit_row_flag;
                };

                while (p == pos_row)
                {
                    Ad_x[i]         = Z;
                    any_modif       = true;
                    pos             += step_ri;

                    if (pos > or-1 || pos < 0)
                        goto exit_flag;

                    matcl::details::pos2ind(ptr_ri_raw[pos],r,pos_row,pos_col);
                    if (j != pos_col)
                        goto exit_row_flag;
                };
            };

            while (j == pos_col)
            {
                pos+= step_ri;
                if (pos > or-1 || pos < 0)
                    goto exit_flag;

                matcl::details::pos2ind(ptr_ri_raw[pos],r,pos_row,pos_col);
            };

            exit_row_flag:
            continue;
        };

    exit_flag:

        if (any_modif == false)
            return out;

        out.get_struct().reset();
        return out;
    };

    static SM eval_1(SM& out,const matcl::details::colon_info& ci, Real tol0)
    {
        using value_type = typename SM::value_type;
        Integer r   = out.rows();
        Integer c   = out.cols();
        Real tol    = tol0;

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

        bool any_deleted	= false;

        raw::details::sparse_ccs<value_type>& Ad = out.rep();
        Integer * Ad_c		= Ad.ptr_c();
        Integer * Ad_r		= Ad.ptr_r();
        value_type * Ad_x	= Ad.ptr_x();		

        Integer pos = rf, pos_row, pos_col;
        matcl::details::pos2ind(pos,r,pos_row,pos_col);

        for (Integer j = 0; j < c; ++j)
        {
            if (j < pos_col)
                continue;

            for (Integer i = Ad_c[j]; i < Ad_c[j+1];++i)
            {
                Integer p = Ad_r[i];
                if (p < pos_row)
                    continue;

                while (p > pos_row)
                {
                    pos+=rs;
                    if (pos > rl)
                        goto exit_flag;

                    matcl::details::pos2ind(pos,r,pos_row,pos_col);
                    if (j != pos_col)
                        goto exit_row_flag;
                };

                while (p == pos_row)
                {
                    const value_type& tmp   = Ad_x[i];

                    if (mrd::abs_helper<value_type>::eval(tmp) <= tol)
                    {
                        Ad_r[i] = -1;
                        any_deleted = true;
                    };


                    pos+=rs;
                    if (pos > rl)
                        goto exit_flag;

                    matcl::details::pos2ind(pos,r,pos_row,pos_col);
                    if (j != pos_col)
                        goto exit_row_flag;
                };
            };

            while (j == pos_col)
            {
                pos+= rs;
                if (pos > rl)
                    goto exit_flag;

                matcl::details::pos2ind(pos,r,pos_row,pos_col);
            };

            exit_row_flag:
            continue;
        };

    exit_flag:

        if (any_deleted == false)
            return out;
        
        pos = Ad.offset();

        for (Integer i = 0; i < c; ++i)
        {
            for (Integer j = Ad_c[i]; j < Ad_c[i+1]; ++j)
            {
                if (Ad_r[j] == -1)
                    continue;

                if (j == pos)
                {
                    ++pos;
                    continue;
                };

                Ad_r[pos] = Ad_r[j];
                mrd::reset_helper(Ad_x[pos],Ad_x[j]);
                ++pos;
            };
            Ad_c[i] = pos;
        };

        for (Integer i = c; i > 0; --i)
            Ad_c[i] = Ad_c[i-1];

        Ad_c[0] = Ad.offset();

        Ad.add_memory(-1);
        out.get_struct().reset();
        return out;
    };

    static SM eval_1_zero(SM& out,const matcl::details::colon_info& ci)
    {
        using value_type = typename SM::value_type;
        Integer r   = out.rows();
        Integer c   = out.cols();

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

        bool any_modif	    = false;

        raw::details::sparse_ccs<value_type>& Ad = out.rep();
        Integer * Ad_c		= Ad.ptr_c();
        Integer * Ad_r		= Ad.ptr_r();
        value_type * Ad_x	= Ad.ptr_x();		
        const value_type& Z = md::default_value<value_type>(out.get_type());

        Integer pos = rf, pos_row, pos_col;
        matcl::details::pos2ind(pos,r,pos_row,pos_col);

        for (Integer j = 0; j < c; ++j)
        {
            if (j < pos_col)
                continue;

            for (Integer i = Ad_c[j]; i < Ad_c[j+1];++i)
            {
                Integer p = Ad_r[i];
                if (p < pos_row)
                    continue;

                while (p > pos_row)
                {
                    pos+=rs;
                    if (pos > rl)
                        goto exit_flag;

                    matcl::details::pos2ind(pos,r,pos_row,pos_col);
                    if (j != pos_col)
                        goto exit_row_flag;
                };

                while (p == pos_row)
                {
                    Ad_x[i]     = Z;
                    any_modif   = true;
                    pos         += rs;

                    if (pos > rl)
                        goto exit_flag;

                    matcl::details::pos2ind(pos,r,pos_row,pos_col);
                    if (j != pos_col)
                        goto exit_row_flag;
                };
            };

            while (j == pos_col)
            {
                pos+= rs;
                if (pos > rl)
                    goto exit_flag;

                matcl::details::pos2ind(pos,r,pos_row,pos_col);
            };

            exit_row_flag:
            continue;
        };

    exit_flag:

        if (any_modif == false)
            return out;
        
        out.get_struct().reset();
        return out;
    };

};

template<class SM>
SM drop_entries_functor<SM>::eval(SM& A,const matcl::details::colon_info& ci, Real tol)
{
    if (ci.rows() == 0 || ci.cols() == 0 || A.nnz() == 0)
        return A;

    if (ci.r_flag == 0)
        return drop_entries_impl<SM>::eval_20(A,ci, tol);
    else
        return drop_entries_impl<SM>::eval_21(A,ci, tol);
};

template<class SM>
SM zero_entries_functor<SM>::eval(SM& A,const matcl::details::colon_info& ci)
{
    if (ci.rows() == 0 || ci.cols() == 0 || A.nnz() == 0)
        return A;

    if (ci.r_flag == 0)
        return drop_entries_impl<SM>::eval_20_zero(A,ci);
    else
        return drop_entries_impl<SM>::eval_21_zero(A,ci);
};

template<class SM>
SM drop_entries_functor_2<SM>::eval(SM& A,const matcl::details::colon_info& ci, Real tol)
{
    if (ci.rows() == 0 || A.nnz() == 0)
        return A;

    if (ci.r_flag == 0)
        return drop_entries_impl<SM>::eval_0(A,ci, tol);
    else
        return drop_entries_impl<SM>::eval_1(A,ci, tol);
};

template<class SM>
SM zero_entries_functor_2<SM>::eval(SM& A,const matcl::details::colon_info& ci)
{
    if (ci.rows() == 0 || A.nnz() == 0)
        return A;

    if (ci.r_flag == 0)
        return drop_entries_impl<SM>::eval_0_zero(A,ci);
    else
        return drop_entries_impl<SM>::eval_1_zero(A,ci);
};

};};};

template struct matcl::algorithm::details::drop_entries_functor<matcl::raw::integer_sparse>;
template struct matcl::algorithm::details::drop_entries_functor<matcl::raw::real_sparse>;
template struct matcl::algorithm::details::drop_entries_functor<matcl::raw::float_sparse>;
template struct matcl::algorithm::details::drop_entries_functor<matcl::raw::complex_sparse>;
template struct matcl::algorithm::details::drop_entries_functor<matcl::raw::float_complex_sparse>;
template struct matcl::algorithm::details::drop_entries_functor<matcl::raw::object_sparse>;

template struct matcl::algorithm::details::drop_entries_functor_2<matcl::raw::integer_sparse>;
template struct matcl::algorithm::details::drop_entries_functor_2<matcl::raw::real_sparse>;
template struct matcl::algorithm::details::drop_entries_functor_2<matcl::raw::float_sparse>;
template struct matcl::algorithm::details::drop_entries_functor_2<matcl::raw::complex_sparse>;
template struct matcl::algorithm::details::drop_entries_functor_2<matcl::raw::float_complex_sparse>;
template struct matcl::algorithm::details::drop_entries_functor_2<matcl::raw::object_sparse>;

template struct matcl::algorithm::details::zero_entries_functor<matcl::raw::integer_sparse>;
template struct matcl::algorithm::details::zero_entries_functor<matcl::raw::real_sparse>;
template struct matcl::algorithm::details::zero_entries_functor<matcl::raw::float_sparse>;
template struct matcl::algorithm::details::zero_entries_functor<matcl::raw::complex_sparse>;
template struct matcl::algorithm::details::zero_entries_functor<matcl::raw::float_complex_sparse>;
template struct matcl::algorithm::details::zero_entries_functor<matcl::raw::object_sparse>;

template struct matcl::algorithm::details::zero_entries_functor_2<matcl::raw::integer_sparse>;
template struct matcl::algorithm::details::zero_entries_functor_2<matcl::raw::real_sparse>;
template struct matcl::algorithm::details::zero_entries_functor_2<matcl::raw::float_sparse>;
template struct matcl::algorithm::details::zero_entries_functor_2<matcl::raw::complex_sparse>;
template struct matcl::algorithm::details::zero_entries_functor_2<matcl::raw::float_complex_sparse>;
template struct matcl::algorithm::details::zero_entries_functor_2<matcl::raw::object_sparse>;
