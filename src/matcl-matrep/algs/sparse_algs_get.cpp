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
#include "matcl-matrep/func/raw/raw_manip.h"
#include "matcl-internals/container/mat_s.h"
#include "matcl-scalar/details/scalfunc_helpers.h"

namespace matcl { namespace algorithm { namespace details
{

namespace mr = matcl::raw;
namespace mrd = matcl::raw::details;

template<class SM>
SM get_submatrix_functor<SM>::eval(const SM& A,const matcl::details::colon_info& ci)
{
    if (ci.r_flag == 0 && ci.c_flag == 0)
        return eval_00(A,ci.get_rim_2(),ci.get_cim_2());

    if (ci.r_flag == 0 && ci.c_flag == 1)
        return eval_01(A,ci);

    if (ci.r_flag == 1 && ci.c_flag == 0)
        return eval_10(A,ci);

    if (ci.r_flag == 1 && ci.c_flag == 1)
        return eval_11(A,ci);

    matcl_assert(0,"invalid colon info");
    throw;
}

template<class SM>
SM get_submatrix_functor<SM>::eval_00(const SM& A,const Vector& ri,const Vector& ci)
{
    if (ri.size() == 0 || ci.size() == 0 || A.nnz() == 0)
        return SM(A.get_type(),ri.size(),ci.size());

    using value_type = SM::value_type;

    Integer r = A.rows(), c = A.cols();
    Integer or = ri.size(), oc = ci.size();
    const raw::details::sparse_ccs<value_type>& Ad = A.rep();

    Integer onnz = icast(Real(A.nnz())/Real(r)/Real(c)*Real(or)*Real(oc)) + 1 + or;	
    raw::details::sparse_ccs<value_type> d(A.get_type(),or, oc, onnz);

    Integer r_min       = 0;
    Integer r_max       = 0;
    sort_type s_type    = not_sorted;
    Integer n_dupl      = 0;

    matcl::pod_workspace<Integer> v_work_ind;
    row_map m_row_map;

    bool workspace_initialized = false;	

    Integer nz					= 0;

    const Integer * Ad_c		= Ad.ptr_c();
    const Integer * Ad_r		= Ad.ptr_r();
    const value_type * Ad_x		= Ad.ptr_x();

    Integer * d_c				= d.ptr_c();
    Integer * d_r				= d.ptr_r();
    value_type * d_x			= d.ptr_x();

    const Integer* ptr_ri       = ri.ptr();
    const Integer* ptr_ci       = ci.ptr();

    for (Integer j = 0; j < oc; ++j)
    {
        Integer col             = ptr_ci[j]-1;	
        d_c[j]				    = nz;
        Integer nnz_row			= Ad_c[col+1] - Ad_c[col];

        if (nnz_row == 0)
            continue;

        if (workspace_initialized == false)
        {
            init_row_selector(ri,-1,v_work_ind,r_min,r_max);
            s_type = is_sorted(ri);

            for (Integer i = 0; i < or; ++i)
            {
                Integer pos = ptr_ri[i]-1-r_min;

                if (v_work_ind[pos] == -1)
                {
                    v_work_ind[pos] = i;
                }
                else if (v_work_ind[pos] == -2)
                {
                    ++n_dupl;
                    m_row_map.insert(pos,i);
                }
                else
                {
                    // given row must be insterted many times
                    ++n_dupl;
                    m_row_map.insert(pos,v_work_ind[pos]);
                    m_row_map.insert(pos,i);
                    v_work_ind[pos] = -2;
                };
            };

            workspace_initialized = true;
        };

        if (nz + nnz_row + n_dupl > d.nzmax()) 
        {
            d.add_memory( d.nzmax() + nnz_row + n_dupl);
            d_c					= d.ptr_c();
            d_r					= d.ptr_r();
            d_x					= d.ptr_x();
        };	

        if (s_type == sorted_decreasing)
        {
            for (Integer k = Ad_c[col+1]-1; k >= Ad_c[col]; --k)
            {			
                Integer p	= Ad_r[k];

                if (p<r_min || p>r_max)
                    continue;

                Integer pos = v_work_ind[p-r_min];
                if (pos == -1)
                {
                    continue;
                }
                else if (pos == -2)
                {
                    m_row_map.add_rows(p-r_min,Ad_x[k],d_x,d_r,nz);
                }
                else
                {
                    mrd::reset_helper(d_x[nz],Ad_x[k]);
                    d_r[nz]		= pos;
                    ++nz;
                };
            };
        }
        else
        {
            Integer nz_old = nz;

            for (Integer k = Ad_c[col]; k < Ad_c[col+1]; ++k)
            {			
                Integer p	= Ad_r[k];
                if (p<r_min || p>r_max)
                    continue;

                Integer pos = v_work_ind[p-r_min];
                if (pos == -1)
                {
                    continue;
                }
                else if (pos == -2)
                {
                    m_row_map.add_rows(p-r_min,Ad_x[k],d_x,d_r,nz);
                }
                else
                {
                    mrd::reset_helper(d_x[nz],Ad_x[k]);
                    d_r[nz]		= pos;
                    ++nz;
                };
            };

            if (s_type == not_sorted && (nz - nz_old > 1) )
                utils::sort_q(d_r+nz_old,d_x + nz_old, nz - nz_old);
        };
    };

    d_c[oc] = nz;
    d.add_memory(-1);

    return raw::sparse_matrix_base<value_type>(d);
};

template<class SM>
SM get_submatrix_functor<SM>::eval_01(const SM& A,const matcl::details::colon_info& colon_info)
{
    if (colon_info.get_rim_2().size() == 0 || colon_info.c_size == 0 || A.nnz() == 0)
    {
        SM out(A.get_type(),colon_info.get_rim_2().size(),colon_info.c_size);
        return out;
    };

    using value_type = typename SM::value_type;

    Integer r   = A.rows();
    Integer c   = A.cols();
    Integer or  = colon_info.get_rim_2().size();
    Integer oc  = colon_info.c_size;
    const raw::details::sparse_ccs<value_type>& Ad = A.rep();

    Integer onnz = icast(Real(A.nnz())/Real(r)/Real(c)*Real(or)*Real(oc)) + 1 + or;	
    raw::details::sparse_ccs<value_type> d(A.get_type(),or, oc, onnz);

    Integer r_min       = 0;
    Integer r_max       = 0;
    sort_type s_type    = not_sorted;
    Integer n_dupl      = 0;

    matcl::pod_workspace<Integer> v_work_ind;
    row_map m_row_map;

    raw::integer_dense ri   = colon_info.get_rim_2();
    const Integer* ptr_ri   = ri.ptr();

    bool workspace_initialized = false;	

    Integer nz					= 0;

    const Integer * Ad_c		= Ad.ptr_c();
    const Integer * Ad_r		= Ad.ptr_r();
    const value_type * Ad_x		= Ad.ptr_x();

    Integer * d_c				= d.ptr_c();
    Integer * d_r				= d.ptr_r();
    value_type * d_x			= d.ptr_x();

    for (Integer j = 1, col = colon_info.c_start-1, pos_c = 0; j <= colon_info.c_size; 
                                            ++j, col += colon_info.c_step, ++pos_c)
    {
        d_c[pos_c]				= nz;

        Integer nnz_row			= Ad_c[col+1] - Ad_c[col];

        if (nnz_row == 0)
            continue;

        if (workspace_initialized == false)
        {
            init_row_selector(ri,-1,v_work_ind,r_min,r_max);
            s_type = is_sorted(ri);

            for (Integer i = 0; i < or; ++i)
            {
                Integer pos = ptr_ri[i]-1-r_min;
                if (v_work_ind[pos] == -1)
                {
                    v_work_ind[pos] = i;
                }
                else if (v_work_ind[pos] == -2)
                {
                    ++n_dupl;
                    m_row_map.insert(pos,i);
                }
                else
                {					
                    // given row must be insterted many times
                    ++n_dupl;
                    m_row_map.insert(pos,v_work_ind[pos]);
                    m_row_map.insert(pos,i);
                    v_work_ind[pos] = -2;
                };
            };

            workspace_initialized = true;
        };

        if (nz + nnz_row  + n_dupl> d.nzmax()) 
        {
            d.add_memory( d.nzmax() + nnz_row + n_dupl);
            d_c					= d.ptr_c();
            d_r					= d.ptr_r();
            d_x					= d.ptr_x();
        };	

        if (s_type == sorted_decreasing)
        {
            for (Integer k = Ad_c[col+1]-1; k >= Ad_c[col]; --k)
            {			
                Integer p	= Ad_r[k];
                if (p<r_min || p>r_max)
                    continue;

                Integer pos = v_work_ind[p-r_min];
                if (pos == -1)
                {
                    continue;
                }
                else if (pos == -2)
                {
                    m_row_map.add_rows(p-r_min,Ad_x[k],d_x,d_r,nz);
                }
                else
                {
                    mrd::reset_helper(d_x[nz],Ad_x[k]);
                    d_r[nz]		= pos;
                    ++nz;
                };
            };
        }
        else
        {
            Integer nz_old = nz;

            for (Integer k = Ad_c[col]; k < Ad_c[col+1]; ++k)
            {			
                Integer p	= Ad_r[k];
                if (p<r_min || p>r_max)
                    continue;

                Integer pos = v_work_ind[p-r_min];
                if (pos == -1)
                {
                    continue;
                }
                else if (pos == -2)
                {
                    m_row_map.add_rows(p-r_min,Ad_x[k],d_x,d_r,nz);
                }
                else
                {
                    mrd::reset_helper(d_x[nz],Ad_x[k]);
                    d_r[nz]		= pos;
                    ++nz;
                };
            };

            if (s_type == not_sorted && (nz - nz_old > 1) )
                utils::sort_q(d_r+nz_old,d_x + nz_old, nz - nz_old);
        };
    };

    d_c[oc] = nz;
    d.add_memory(-1);

    return raw::sparse_matrix_base<value_type>(d);
};

template<class SM>
SM get_submatrix_functor<SM>::eval_10(const SM& A,const matcl::details::colon_info& colon_info)
{
    if (colon_info.r_size == 0 || colon_info.get_cim_2().size() == 0 || A.nnz() == 0)
    {
        SM out(A.get_type(),colon_info.r_size,colon_info.get_cim_2().size());
        return out;
    };

    if (colon_info.r_start == 1 && colon_info.r_step == 1 && colon_info.r_end == A.rows())
        return get_cols_0(A,colon_info);

    using value_type = typename SM::value_type;

    Integer r   = A.rows();
    Integer c   = A.cols();
    Integer or  = colon_info.r_size;
    Integer oc  = colon_info.get_cim_2().size();
    const raw::details::sparse_ccs<value_type>& Ad = A.rep();

    Integer onnz = icast(Real(A.nnz())/Real(r)/Real(c)*Real(or)*Real(oc));	
    raw::details::sparse_ccs<value_type> d(A.get_type(),or, oc, onnz);

    Integer r_start = colon_info.r_start - 1;
    Integer r_end = colon_info.r_end - 1;
    Integer r_step = colon_info.r_step;
    bool b_increasing = (r_step > 0);

    if (!b_increasing)
    {
        Integer tmp = r_start;
        r_start = r_end;
        r_end = tmp;
        r_step = - r_step;
    };

    Integer nz					= 0;

    const Integer * Ad_c		= Ad.ptr_c();
    const Integer * Ad_r		= Ad.ptr_r();
    const value_type * Ad_x		= Ad.ptr_x();

    Integer * d_c				= d.ptr_c();
    Integer * d_r				= d.ptr_r();
    value_type * d_x			= d.ptr_x();	

    mr::integer_dense ci_ci     = colon_info.get_cim_2();
    const Integer* ptr_ci       = ci_ci.ptr();

    for (Integer j = 0; j < oc; ++j)
    {
        Integer col				= ptr_ci[j]-1;
        d_c[j]				    = nz;
        Integer nnz_row			= Ad_c[col+1] - Ad_c[col];

        if (nnz_row == 0)
            continue;

        if (nz + nnz_row > d.nzmax()) 
        {
            d.add_memory( d.nzmax() + nnz_row);
            d_c				= d.ptr_c();
            d_r				= d.ptr_r();
            d_x				= d.ptr_x();	
        };

        if (b_increasing == false)
        {
            if (r_step == 1)
            {
                for (Integer k = Ad_c[col+1]-1; k >= Ad_c[col]; --k)
                {			
                    Integer p	= Ad_r[k];
                    if (p<r_start || p>r_end)
                        continue;

                    Integer pos = r_end-p;
                    mrd::reset_helper(d_x[nz],Ad_x[k]);
                    d_r[nz]		= pos;
                    ++nz;
                };
            }
            else
            {
                for (Integer k = Ad_c[col+1]-1; k >= Ad_c[col]; --k)
                {			
                    Integer p	= Ad_r[k];
                    if (p<r_start || p>r_end)
                        continue;

                    if ((r_end-p)%r_step != 0)
                        continue;

                    Integer pos = (r_end-p)/r_step;
                    mrd::reset_helper(d_x[nz],Ad_x[k]);
                    d_r[nz]		= pos;
                    ++nz;
                };
            };
        }
        else
        {
            if (r_step == 1)
            {
                for (Integer k = Ad_c[col]; k < Ad_c[col+1]; ++k)
                {			
                    Integer p	= Ad_r[k];
                    if (p<r_start || p>r_end)
                        continue;

                    Integer pos = p-r_start;
                    mrd::reset_helper(d_x[nz],Ad_x[k]);
                    d_r[nz]		= pos;
                    ++nz;
                };
            }
            else
            {
                for (Integer k = Ad_c[col]; k < Ad_c[col+1]; ++k)
                {			
                    Integer p	= Ad_r[k];
                    if (p<r_start || p>r_end)
                        continue;

                    if ((p-r_start)%r_step != 0)
                        continue;

                    Integer pos = (p-r_start)/r_step;
                    mrd::reset_helper(d_x[nz],Ad_x[k]);
                    d_r[nz]		= pos;
                    ++nz;
                };
            };
        };
    };
    d_c[oc] = nz;
    d.add_memory(-1);

    return raw::sparse_matrix_base<value_type>(d);
};

template<class SM>
SM get_submatrix_functor<SM>::eval_11(const SM& A,const matcl::details::colon_info& colon_info)
{
    if (colon_info.r_size == 0 || colon_info.c_size == 0 || A.nnz() == 0)
    {
        SM out(A.get_type(),colon_info.r_size,colon_info.c_size);
        return out;
    };

    if (colon_info.r_start == 1 && colon_info.r_step == 1 && colon_info.r_end == A.rows())
        return get_cols_1(A,colon_info);

    using value_type = typename SM::value_type;

    Integer r = A.rows(), c = A.cols();
    Integer or = colon_info.r_size, oc = colon_info.c_size;
    const raw::details::sparse_ccs<value_type>& Ad = A.rep();

    Integer onnz = icast(Real(A.nnz())/Real(r)/Real(c)*Real(or)*Real(oc));	
    raw::details::sparse_ccs<value_type> d(A.get_type(),or, oc, onnz);

    Integer r_start = colon_info.r_start - 1;
    Integer r_end = colon_info.r_end - 1;
    Integer r_step = colon_info.r_step;
    bool b_increasing = (r_step > 0);

    if (!b_increasing)
    {
        Integer tmp = r_start;
        r_start = r_end;
        r_end = tmp;
        r_step = - r_step;
    };

    Integer nz					= 0;

    const Integer * Ad_c		= Ad.ptr_c();
    const Integer * Ad_r		= Ad.ptr_r();
    const value_type * Ad_x		= Ad.ptr_x();

    Integer * d_c				= d.ptr_c();
    Integer * d_r				= d.ptr_r();
    value_type * d_x			= d.ptr_x();	

    for (Integer j = 1, col = colon_info.c_start-1, pos_c = 0; j <= colon_info.c_size; 
                                                ++j, col += colon_info.c_step, ++pos_c)
    {	
        d_c[pos_c]				= nz;

        Integer nnz_row			= Ad_c[col+1] - Ad_c[col];

        if (nnz_row == 0)
            continue;

        if (nz + nnz_row > d.nzmax()) 
        {
            d.add_memory( d.nzmax() + nnz_row);
            d_c				= d.ptr_c();
            d_r				= d.ptr_r();
            d_x				= d.ptr_x();	
        };

        if (b_increasing == false)
        {
            if (r_step == 1)
            {
                for (Integer k = Ad_c[col+1]-1; k >= Ad_c[col]; --k)
                {			
                    Integer p	= Ad_r[k];
                    if (p<r_start || p>r_end)
                        continue;

                    Integer pos = r_end-p;
                    mrd::reset_helper(d_x[nz],Ad_x[k]);
                    d_r[nz]		= pos;
                    ++nz;
                };
            }
            else
            {
                for (Integer k = Ad_c[col+1]-1; k >= Ad_c[col]; --k)
                {			
                    Integer p	= Ad_r[k];
                    if (p<r_start || p>r_end)
                        continue;

                    if ((r_end-p)%r_step != 0)
                        continue;

                    Integer pos = (r_end-p)/r_step;
                    mrd::reset_helper(d_x[nz],Ad_x[k]);
                    d_r[nz]		= pos;
                    ++nz;
                };
            };
        }
        else
        {
            if (r_step == 1)
            {
                for (Integer k = Ad_c[col]; k < Ad_c[col+1]; ++k)
                {			
                    Integer p	= Ad_r[k];
                    if (p<r_start || p>r_end)
                        continue;

                    Integer pos = p-r_start;
                    mrd::reset_helper(d_x[nz],Ad_x[k]);
                    d_r[nz]		= pos;
                    ++nz;
                };
            }
            else
            {
                for (Integer k = Ad_c[col]; k < Ad_c[col+1]; ++k)
                {			
                    Integer p	= Ad_r[k];
                    if (p<r_start || p>r_end)
                        continue;

                    if ((p-r_start)%r_step != 0)
                        continue;

                    Integer pos = (p-r_start)/r_step;
                    mrd::reset_helper(d_x[nz],Ad_x[k]);
                    d_r[nz]		= pos;
                    ++nz;
                };
            };
        };
    };

    d_c[oc] = nz;
    d.add_memory(-1);

    return raw::sparse_matrix_base<value_type>(d);
};

template<class SM>
SM get_submatrix_functor<SM>::get_cols_0(const SM& A,const matcl::details::colon_info& colon_info)
{
    using value_type = typename SM::value_type;

    Integer r   = A.rows();
    Integer c   = A.cols();
    Integer or  = r;
    Integer oc  = colon_info.get_cim_2().size();
    const raw::details::sparse_ccs<value_type>& Ad = A.rep();

    Integer onnz = icast(Real(A.nnz())/Real(r)/Real(c)*Real(oc));
    raw::details::sparse_ccs<value_type> d(A.get_type(),or, oc, onnz);

    if (or == 0 || oc == 0 || A.nnz() == 0)
    {
        SM out= raw::sparse_matrix_base<value_type>(d);
        return out;
    };

    Integer nz					= 0;

    const Integer * Ad_c		= Ad.ptr_c();
    const Integer * Ad_r		= Ad.ptr_r();
    const value_type * Ad_x		= Ad.ptr_x();

    Integer * d_c				= d.ptr_c();
    Integer * d_r				= d.ptr_r();
    value_type * d_x			= d.ptr_x();	
    mr::integer_dense ci_ci     = colon_info.get_cim_2();
    const Integer* ptr_ci       = ci_ci.ptr();

    for (Integer j = 0; j < oc; ++j)
    {
        Integer col				= ptr_ci[j]-1;	
        d_c[j]				    = nz;
        Integer nnz_row			= Ad_c[col+1] - Ad_c[col];

        if (nz + nnz_row > d.nzmax()) 
        {
            d.add_memory( d.nzmax() + nnz_row);
            d_c				= d.ptr_c();
            d_r				= d.ptr_r();
            d_x				= d.ptr_x();	
        };	

        for (Integer k = Ad_c[col]; k < Ad_c[col+1]; ++k)
        {
            mrd::reset_helper(d_x[nz],Ad_x[k]);
            d_r[nz]			= Ad_r[k];
            ++nz;
        };
    };

    d_c[oc] = nz;
    d.add_memory(-1);

    return raw::sparse_matrix_base<value_type>(d);
};

template<class SM>
SM get_submatrix_functor<SM>::get_cols_1(const SM& A,const matcl::details::colon_info& colon_info)
{
    using value_type = typename SM::value_type;

    Integer r = A.rows(), c = A.cols();
    Integer or = r, oc = colon_info.c_size;
    const raw::details::sparse_ccs<value_type>& Ad = A.rep();

    Integer onnz = icast(Real(A.nnz())/Real(r)/Real(c)*Real(oc));
    raw::details::sparse_ccs<value_type> d(A.get_type(),or, oc, onnz);

    if (or == 0 || oc == 0 || A.nnz() == 0)
    {
        SM out = raw::sparse_matrix_base<value_type>(d);
        return out;
    };

    Integer nz					= 0;

    const Integer * Ad_c		= Ad.ptr_c();
    const Integer * Ad_r		= Ad.ptr_r();
    const value_type * Ad_x		= Ad.ptr_x();

    Integer * d_c				= d.ptr_c();
    Integer * d_r				= d.ptr_r();
    value_type * d_x			= d.ptr_x();	

    for (Integer j = 1, col = colon_info.c_start-1, pos_c = 0; j <= colon_info.c_size; 
                    ++j, col+= colon_info.c_step, ++pos_c)
    {
        d_c[pos_c]				= nz;

        Integer nnz_row			= Ad_c[col+1] - Ad_c[col];

        if (nz + nnz_row > d.nzmax()) 
        {
            d.add_memory( d.nzmax() + nnz_row);
            d_c				= d.ptr_c();
            d_r				= d.ptr_r();
            d_x				= d.ptr_x();	
        };	

        for (Integer k = Ad_c[col]; k < Ad_c[col+1]; ++k)
        {
            mrd::reset_helper(d_x[nz],Ad_x[k]);
            d_r[nz]			= Ad_r[k];
            ++nz;
        };
    };

    d_c[oc] = nz;
    d.add_memory(-1);

    return raw::sparse_matrix_base<value_type>(d);
};

template<class SM>
void get_submatrix_functor_2<SM>::eval(matcl::Matrix& ret, const SM& A,const matcl::details::colon_info& ci)
{
    if (ci.r_flag == 0)
    {
        if (ci.is_double_mat_colon() == true)
        {
            SM out = eval_0_dc(A,ci.get_rim_r(),ci.get_rim_c());
            return mrd::manip_reshape_helper<SM>::eval_reshape(ret, out, ci.rep_rows(), ci.rep_cols());
        }
        else
        {
            SM out = eval_0(A,ci.get_rim_1());
            return mrd::manip_reshape_helper<SM>::eval_reshape(ret, out, ci.rep_rows(), ci.rep_cols());
        }
    };

    if (ci.r_flag == 1)
    {
        SM out = eval_1(A,ci);
        if (out.rows() == ci.rep_rows())
        {
            ret = matcl::Matrix(out,true);
            return;
        }
        else
        {
            return mrd::manip_reshape_helper<SM>::eval_reshape(ret,out,ci.rep_rows(), ci.rep_cols());
        };
    };
    
    matcl_assert(0,"invalid colon info");
    throw;
}

template<class SM>
SM get_submatrix_functor_2<SM>::eval_0(const SM& A,const Vector& ci)
{
    if (ci.size() == 0 || A.nnz() == 0)
    {
        SM out(A.get_type(),ci.rows(),ci.cols());
        return out;
    };	

    sort_type s_type = is_sorted(ci);
    
    if (s_type == sorted_increasing)
    {
        Integer n_rep = number_dupl(ci);
        return eval_0_increasing(A,ci,n_rep);
    };
    
    if (s_type == sorted_decreasing)
    {
        Integer n_rep = number_dupl(ci);
        return eval_0_decreasing(A,ci,n_rep);
    };

    Vector ci2 = ci.copy();
    ci2.get_struct().reset();

    using value_type = typename SM::value_type;

    Integer or = ci2.size();

    matcl::pod_workspace<Integer>	v_work_sort_ind(or);	

    for (Integer i = 0; i < or; ++i)
        v_work_sort_ind[i] = i;

    utils::sort_q(ci2.ptr(),&v_work_sort_ind[0],or);

    Integer n_rep   = number_dupl(ci2);
    SM out          = eval_0_increasing(A,ci2,n_rep);

    Integer off     = out.rep().offset();
    Integer* d_r    = out.rep().ptr_r() + off;
    Integer d_nnz   = out.nnz();

    for (Integer i = 0; i < d_nnz; ++i)
        d_r[i] = v_work_sort_ind[d_r[i]];

    utils::sort_q(out.rep().ptr_r() + off,out.rep().ptr_x() + off,d_nnz);

    out.get_struct().reset();
    return out;
};

template<class SM>
SM get_submatrix_functor_2<SM>::eval_0_dc(const SM& A,const Vector& ri,const Vector& ci)
{
    if (ci.size() == 0 || A.nnz() == 0)
    {
        SM out(A.get_type(),ci.rows(),ci.cols());
        return out;
    };

    Vector ri2  = ri.copy();
    Vector ci2  = ci.copy();

    ri2.get_struct().reset();
    ci2.get_struct().reset();

    using value_type = typename SM::value_type;

    Integer or = ci2.size();

    matcl::pod_workspace<Integer>	v_work_sort_ind(or);	

    for (Integer i = 0; i < or; ++i)
        v_work_sort_ind[i] = i;

    sort_rows_cols(ri2.ptr(), ci2.ptr(), &v_work_sort_ind[0], or);

    Integer n_rep   = number_dupl(ri2, ci2);
    SM out          = eval_0_dc_increasing(A, ri2, ci2, n_rep);

    Integer off     = out.rep().offset();
    Integer* d_r    = out.rep().ptr_r() + off;
    Integer d_nnz   = out.nnz();

    for (Integer i = 0; i < d_nnz; ++i)
        d_r[i] = v_work_sort_ind[d_r[i]];

    utils::sort_q(out.rep().ptr_r() + off,out.rep().ptr_x() + off,d_nnz);

    out.get_struct().reset();
    return out;
};

template<class SM>
SM get_submatrix_functor_2<SM>::eval_1(const SM& A,const matcl::details::colon_info& ci)
{
    if (ci.rows() == 0 || A.nnz() == 0)
    {
        SM out(A.get_type(),ci.rows(),1);
        return out;
    };

    if (ci.r_start < ci.r_end)
        return eval_1_increasing(A,ci);
    else
        return eval_1_decreasing(A,ci);
};

template<class SM>
SM get_submatrix_functor_2<SM>::eval_1_increasing(const SM& A,const matcl::details::colon_info& ci)
{
    using value_type = typename SM::value_type;

    Integer r = A.rows(), c = A.cols();
    Integer or = ci.rows();
    const raw::details::sparse_ccs<value_type>& Ad = A.rep();

    Integer onnz = icast(Real(A.nnz())/Real(r)/Real(c)*Real(or)) + 1;

    raw::details::sparse_ccs<value_type> d(A.get_type(),or,1, onnz);	

    if (or == 0 || A.nnz() == 0)
    {
        SM out = raw::sparse_matrix_base<value_type>(d);
        return out;
    };

    const Integer * Ad_c		= Ad.ptr_c();
    const Integer * Ad_r		= Ad.ptr_r();
    const value_type * Ad_x		= Ad.ptr_x();

    Integer pos = 1, pos_val = ci.r_start, pos_row, pos_col;
    matcl::details::pos2ind(pos_val,r,pos_row,pos_col);

    Integer * d_c				= d.ptr_c();
    Integer * d_r				= d.ptr_r();
    value_type * d_x			= d.ptr_x();	

    d_c[0]						= 0;
    Integer nz					= 0;

    for (Integer j = 0; j < c; ++j)
    {
        if (j < pos_col)
            continue;

        Integer nnz_row			= Ad_c[j+1] - Ad_c[j];

        if (nz + nnz_row > d.nzmax()) 
        {
            d.add_memory( d.nzmax() + nnz_row);
        
            d_c				= d.ptr_c();
            d_r				= d.ptr_r();
            d_x				= d.ptr_x();	
        };

        for (Integer i = Ad_c[j]; i < Ad_c[j+1];++i)
        {
            Integer p = Ad_r[i];
            if (p < pos_row)
                continue;

            while (p > pos_row)
            {
                ++pos;
                pos_val += ci.r_step;
                if (pos > or)
                    goto exit_flag;

                matcl::details::pos2ind(pos_val,r,pos_row,pos_col);
                if (j != pos_col)
                    goto exit_row_flag;
            };

            while (p == pos_row)
            {
                mrd::reset_helper(d_x[nz],Ad_x[i]);
                d_r[nz]		= pos-1;
                ++nz;

                ++pos;
                pos_val += ci.r_step;
                if (pos > or)
                    goto exit_flag;

                matcl::details::pos2ind(pos_val,r,pos_row,pos_col);
                if (j != pos_col)
                    goto exit_row_flag;
            };
        };

        exit_row_flag:

        while (j == pos_col)
        {
            ++pos;
            pos_val += ci.r_step;
            if (pos > or)
                goto exit_flag;

            matcl::details::pos2ind(pos_val,r,pos_row,pos_col);
            if (j != pos_col)
                break;
        };
    };

  exit_flag:

    d_c[1] = nz;
    d.add_memory(-1);

    return raw::sparse_matrix_base<value_type>(d);
};

template<class SM>
SM get_submatrix_functor_2<SM>::eval_1_decreasing(const SM& A,const matcl::details::colon_info& ci)
{
    using value_type = typename SM::value_type;

    Integer r = A.rows(), c = A.cols();
    Integer or = ci.rows();
    const raw::details::sparse_ccs<value_type>& Ad = A.rep();

    Integer onnz = icast(Real(A.nnz())/Real(r)/Real(c)*Real(or)) + 1;

    raw::details::sparse_ccs<value_type> d(A.get_type(),ci.rows(),1, onnz);	

    if (ci.rows() == 0 || A.nnz() == 0)
    {
        SM out = raw::sparse_matrix_base<value_type>(d);
        return out;
    };

    const Integer * Ad_c		= Ad.ptr_c();
    const Integer * Ad_r		= Ad.ptr_r();
    const value_type * Ad_x		= Ad.ptr_x();

    Integer * d_c				= d.ptr_c();
    Integer * d_r				= d.ptr_r();
    value_type * d_x			= d.ptr_x();	

    d_c[0]						= 0;
    Integer nz					= 0;

    Integer pos = 1, pos_val = ci.r_start, pos_row, pos_col;
    matcl::details::pos2ind(pos_val,r,pos_row,pos_col);

    for (Integer j = c-1; j >= 0; --j)
    {
        if (j > pos_col)
            continue;

        Integer nnz_row			= Ad_c[j+1] - Ad_c[j];

        if (nz + nnz_row > d.nzmax()) 
        {
            d.add_memory( d.nzmax() + nnz_row);
            d_c				= d.ptr_c();
            d_r				= d.ptr_r();
            d_x				= d.ptr_x();	
        };

        for (Integer i = Ad_c[j+1] - 1; i >= Ad_c[j]; --i)
        {
            Integer p = Ad_r[i];
            if (p > pos_row)
                continue;

            while (p < pos_row)
            {
                ++pos;
                pos_val += ci.r_step;
                if (pos > or)
                    goto exit_flag;

                matcl::details::pos2ind(pos_val,r,pos_row,pos_col);
                if (j != pos_col)
                    goto exit_row_flag;
            };

            while (p == pos_row)
            {
                mrd::reset_helper(d_x[nz],Ad_x[i]);
                d_r[nz]		= pos-1;
                ++nz;

                ++pos;
                pos_val += ci.r_step;
                if (pos > or)
                    goto exit_flag;

                matcl::details::pos2ind(pos_val,r,pos_row,pos_col);
                if (j != pos_col)
                    goto exit_row_flag;
            };
        };

      exit_row_flag:

        while (j == pos_col)
        {
            ++pos;
            pos_val += ci.r_step;
            if (pos > or)
                goto exit_flag;

            matcl::details::pos2ind(pos_val,r,pos_row,pos_col);
            if (j != pos_col)
                break;
        };
    };

  exit_flag:

    d_c[1] = nz;
    d.add_memory(-1);

    return raw::sparse_matrix_base<value_type>(d);
};

template<class SM>
SM get_submatrix_functor_2<SM>::eval_0_increasing(const SM& A,const Vector& ci,Integer n_rep)
{
    using value_type = typename SM::value_type;

    Integer r = A.rows(), c = A.cols();
    Integer or = ci.size();
    const raw::details::sparse_ccs<value_type>& Ad = A.rep();

    if (r == 0 || c == 0 || A.nnz() == 0)
    {
        SM out(A.get_type(),ci.rows(), ci.cols());
        return out;
    };

    Integer onnz = icast(Real(A.nnz())/Real(r)/Real(c)*Real(or)) + 1;

    raw::details::sparse_ccs<value_type> d(A.get_type(),ci.size(), 1, onnz);	

    const Integer * Ad_c		= Ad.ptr_c();
    const Integer * Ad_r		= Ad.ptr_r();
    const value_type * Ad_x		= Ad.ptr_x();    

    const Integer* ptr_ci       = ci.ptr();

    Integer pos = 0, pos_row, pos_col;
    matcl::details::pos2ind(ptr_ci[pos],r,pos_row,pos_col);

    Integer * d_c		        = d.ptr_c();
    Integer * d_r				= d.ptr_r();
    value_type * d_x			= d.ptr_x();

    d_c[0]						= 0;
    Integer nz					= 0;

    for (Integer j = 0; j < c; ++j)
    {
        if (j < pos_col)
            continue;

        Integer nnz_row			= Ad_c[j+1] - Ad_c[j];

        if (nz + nnz_row + n_rep > d.nzmax()) 
        {
            d.add_memory( d.nzmax() + nnz_row + n_rep);
            d_c					= d.ptr_c();
            d_r					= d.ptr_r();
            d_x					= d.ptr_x();	
        };

        for (Integer i = Ad_c[j]; i < Ad_c[j+1];++i)
        {
            Integer p = Ad_r[i];
            if (p < pos_row)
                continue;

            while (p > pos_row)
            {
                ++pos;
                if (pos >= or)
                    goto exit_flag;

                matcl::details::pos2ind(ptr_ci[pos],r,pos_row,pos_col);
                if (j != pos_col)
                    goto exit_row_flag;
            };

            while (p == pos_row)
            {
                mrd::reset_helper(d_x[nz],Ad_x[i]);
                d_r[nz]		= pos;
                ++nz;

                ++pos;
                if (pos >= or)
                    goto exit_flag;

                matcl::details::pos2ind(ptr_ci[pos],r,pos_row,pos_col);
                if (j != pos_col)
                    goto exit_row_flag;
            };
        };

      exit_row_flag:

        while (j == pos_col)
        {
            ++pos;
            if (pos >= or)     
                goto exit_flag;

            matcl::details::pos2ind(ptr_ci[pos],r,pos_row,pos_col);
            if (j != pos_col)
                break;
        };
    };

  exit_flag:

    d_c[1] = nz;
    d.add_memory(-1);

    SM out = raw::sparse_matrix_base<value_type>(d);	
    return out;
};

template<class SM>
SM get_submatrix_functor_2<SM>::eval_0_dc_increasing(const SM& A, const Vector& cr, const Vector& cc,Integer n_rep)
{
    using value_type = typename SM::value_type;

    Integer r = A.rows(), c = A.cols();
    Integer or = cr.size();
    const raw::details::sparse_ccs<value_type>& Ad = A.rep();

    if (r == 0 || c == 0 || A.nnz() == 0)
    {
        SM out(A.get_type(),cr.rows(), cr.cols());
        return out;
    };

    Integer onnz = icast(Real(A.nnz())/Real(r)/Real(c)*Real(or)) + 1;

    raw::details::sparse_ccs<value_type> d(A.get_type(), cr.size(), 1, onnz);	

    const Integer * Ad_c		= Ad.ptr_c();
    const Integer * Ad_r		= Ad.ptr_r();
    const value_type * Ad_x		= Ad.ptr_x();

    const Integer* ptr_r        = cr.ptr();
    const Integer* ptr_c        = cc.ptr();

    Integer pos = 0, pos_row, pos_col;

    pos_row                     = ptr_r[pos] - 1;
    pos_col                     = ptr_c[pos] - 1;

    Integer * d_c		        = d.ptr_c();
    Integer * d_r				= d.ptr_r();
    value_type * d_x			= d.ptr_x();

    d_c[0]						= 0;
    Integer nz					= 0;

    for (Integer j = 0; j < c; ++j)
    {
        if (j < pos_col)
            continue;

        Integer nnz_row			= Ad_c[j+1] - Ad_c[j];

        if (nz + nnz_row + n_rep > d.nzmax()) 
        {
            d.add_memory( d.nzmax() + nnz_row + n_rep);
            d_c					= d.ptr_c();
            d_r					= d.ptr_r();
            d_x					= d.ptr_x();	
        };

        for (Integer i = Ad_c[j]; i < Ad_c[j+1];++i)
        {
            Integer p = Ad_r[i];
            if (p < pos_row)
                continue;

            while (p > pos_row)
            {
                ++pos;
                if (pos >= or)
                    goto exit_flag;

                pos_row         = ptr_r[pos] - 1;
                pos_col         = ptr_c[pos] - 1;

                if (j != pos_col)
                    goto exit_row_flag;
            };

            while (p == pos_row)
            {
                mrd::reset_helper(d_x[nz],Ad_x[i]);
                d_r[nz]		= pos;
                ++nz;

                ++pos;
                if (pos >= or)
                    goto exit_flag;

                pos_row         = ptr_r[pos] - 1;
                pos_col         = ptr_c[pos] - 1;

                if (j != pos_col)
                    goto exit_row_flag;
            };
        };

      exit_row_flag:

        while (j == pos_col)
        {
            ++pos;
            if (pos >= or)     
                goto exit_flag;

            pos_row         = ptr_r[pos] - 1;
            pos_col         = ptr_c[pos] - 1;

            if (j != pos_col)
                break;
        };
    };

  exit_flag:

    d_c[1] = nz;
    d.add_memory(-1);

    SM out = raw::sparse_matrix_base<value_type>(d);	
    return out;
};

template<class SM>
SM get_submatrix_functor_2<SM>::eval_0_decreasing(const SM& A,const Vector& ci,Integer n_rep)
{
    using value_type = typename SM::value_type;

    Integer r = A.rows(), c = A.cols();
    Integer or = ci.size();
    const raw::details::sparse_ccs<value_type>& Ad = A.rep();

    Integer onnz = icast(Real(A.nnz())/Real(r)/Real(c)*Real(or)) + 1;

    raw::details::sparse_ccs<value_type> d(A.get_type(),ci.size(), 1, onnz);	

    if (A.nnz() == 0)
    {
        SM out = raw::sparse_matrix_base<value_type>(d);
        return out;
    };

    const Integer * Ad_c		= Ad.ptr_c();
    const Integer * Ad_r		= Ad.ptr_r();
    const value_type * Ad_x		= Ad.ptr_x();

    Integer * d_c				= d.ptr_c();
    Integer * d_r				= d.ptr_r();
    value_type * d_x			= d.ptr_x();	

    d_c[0]						= 0;
    Integer nz					= 0;

    const Integer* ptr_ci       = ci.ptr();

    Integer pos = 0, pos_row, pos_col;
    matcl::details::pos2ind(ptr_ci[pos],r,pos_row,pos_col);

    for (Integer j = c-1; j >= 0; --j)
    {
        if (j > pos_col)
            continue;

        Integer nnz_row			= Ad_c[j+1] - Ad_c[j];

        if (nz + nnz_row + n_rep > d.nzmax()) 
        {
            d.add_memory( d.nzmax() + nnz_row + n_rep);
            d_c				= d.ptr_c();
            d_r				= d.ptr_r();
            d_x				= d.ptr_x();	
        };

        for (Integer i = Ad_c[j+1] - 1; i >= Ad_c[j]; --i)
        {
            Integer p = Ad_r[i];
            if (p > pos_row)
                continue;

            while (p < pos_row)
            {
                ++pos;
                if (pos >= or)
                    goto exit_flag;

                matcl::details::pos2ind(ptr_ci[pos],r,pos_row,pos_col);

                if (j != pos_col)
                    goto exit_row_flag;
            };

            while (p == pos_row)
            {
                mrd::reset_helper(d_x[nz],Ad_x[i]);
                d_r[nz]		= pos;
                ++nz;

                ++pos;
                if (pos >= or)
                    goto exit_flag;

                matcl::details::pos2ind(ptr_ci[pos],r,pos_row,pos_col);
                if (j != pos_col)
                    goto exit_row_flag;
            };
        };

        exit_row_flag:

        while (j == pos_col)
        {
            ++pos;
            if (pos >= or)
                goto exit_flag;

            matcl::details::pos2ind(ptr_ci[pos],r,pos_row,pos_col);
            if (j != pos_col)
                break;
        };
    };

  exit_flag:

    d_c[1] = nz;
    d.add_memory(-1);

    return raw::sparse_matrix_base<value_type>(d);
};

template struct get_submatrix_functor<raw::integer_sparse>;
template struct get_submatrix_functor<raw::real_sparse>;
template struct get_submatrix_functor<raw::float_sparse>;
template struct get_submatrix_functor<raw::complex_sparse>;
template struct get_submatrix_functor<raw::float_complex_sparse>;
template struct get_submatrix_functor<raw::object_sparse>;

template struct get_submatrix_functor_2<raw::integer_sparse>;
template struct get_submatrix_functor_2<raw::real_sparse>;
template struct get_submatrix_functor_2<raw::float_sparse>;
template struct get_submatrix_functor_2<raw::complex_sparse>;
template struct get_submatrix_functor_2<raw::float_complex_sparse>;
template struct get_submatrix_functor_2<raw::object_sparse>;

};};};
