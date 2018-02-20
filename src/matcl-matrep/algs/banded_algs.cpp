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
 
#include "matcl-matrep/algs/banded_algs.h"
#include "matcl-internals/func/converter.h"
#include "matcl-internals/base/sort.h"
#include "matcl-internals/base/utils.h"
#include "matcl-internals/container/mat_b.inl"
#include "matcl-internals/container/mat_d.inl"
#include "matcl-internals/container/mat_s.inl"
#include "matcl-matrep/utils/workspace.h"
#include "matcl-scalar/details/scalfunc_helpers.h"
#include "matcl-matrep/details/struct_flag_predefined.h"

namespace matcl { namespace algorithm { namespace details
{

namespace md    = matcl::details;
namespace mrd   = matcl::raw::details;

template<class value_type>
void del_rows_banded_functor<value_type>::eval(matcl::Matrix& ret, const BM& mat,
                        const matcl::details::colon_info& ci, bool)
{
    Integer s = ci.rows();
    if (s == 0)
    {
        ret = matcl::Matrix(mat,false);
        return;
    };

    Integer r0  = mat.rows();
    Integer c0  = mat.cols();

    if (ci.r_flag == 0)
    {		
        if (ci.is_double_mat_colon() == true)
            throw error::invalid_colon_too_many_mat();

        raw::integer_dense ritmp    = ci.get_rim_1();
        const Integer* ptr_r        = ritmp.ptr();
        Integer ce                  = ptr_r[0];

        error::check_row(ce, r0, c0);

        matcl::pod_workspace<Integer> iwork(r0,0);

        for (Integer i = 0; i < s; ++i)
        {
            Integer ce2     = ptr_r[i];
            error::check_row(ce2, r0, c0);
            iwork[ce2-1]    = -1;
        };

        Integer n_del       = 0;
        for (Integer i = 0; i < r0; ++i)
        {
            if (iwork[i] == -1)
                ++n_del;
            else
                iwork[i]    = i-n_del;
        };

        if (!c0)
        {
            ret = matcl::Matrix(SM(mat.get_type(),r0-n_del,0),false);
            return;
        };

        if (n_del == r0)
        {
            ret = matcl::Matrix(SM(mat.get_type(),0,c0), false);
            return;
        };

        Integer r1 = r0 - n_del;		
        SM out(mat.get_type(),r1,c0,mat.nnz());

        raw::details::sparse_ccs<value_type>& d = out.rep();

        Integer* d_c                = d.ptr_c();
        Integer* d_r                = d.ptr_r();
        value_type* d_x             = d.ptr_x();
        const value_type* mat_ptr   = mat.rep_ptr();

        Integer nz          = 0;
        Integer mat_ld      = mat.ld();

        for (Integer j = 0, jj = 0; j < c0; ++j, jj += mat_ld)
        {
            d_c[j]      = nz;

            Integer fr  = mat.first_row(j);
            Integer lr  = mat.last_row(j);
            Integer ii  = jj + mat.first_elem_pos(j);

            for (Integer i = fr; i <= lr; ++i,++ii)
            {
                Integer np  = iwork[i];
                if (np != -1)
                {
                    const value_type& tmp = mat_ptr[ii];
                    d_r[nz] = np;
                    mrd::reset_helper(d_x[nz],tmp);
                    ++nz;
                };
            };
        };

        d_c[c0] = nz;
        d.add_memory(-1);

        ret = matcl::Matrix(out,true);
        return;
    }
    else
    {
        Integer ne = s;
        if (!c0)
        {
            ret = matcl::Matrix(SM(mat.get_type(),r0-ne,0),false);
            return;
        };

        if (ne == r0)
        {
            ret = matcl::Matrix(SM(mat.get_type(),0,c0),false);
            return;
        };

        Integer r1 = r0 - ne;		
        SM out(mat.get_type(),r1,c0,mat.nnz());

        raw::details::sparse_ccs<value_type>& d = out.rep();

        Integer* d_c    = d.ptr_c();
        Integer* d_r    = d.ptr_r();
        value_type* d_x = d.ptr_x();

        const value_type* mat_ptr = mat.rep_ptr();

        Integer rs = ci.r_step, rf, rl;

        if (rs < 0)
        {
            rs  = -rs;
            rf  = ci.r_end;
            rl  = ci.r_start;
        }
        else
        {
            rf  = ci.r_start;		
            rl  = ci.r_end;
        };

        //make it 0-based
        rf  = rf - 1;
        rl  = rl - 1;

        Integer nz      = 0;
        Integer mat_ld   = mat.ld();

        for (Integer j = 0, jj = 0; j < c0; ++j, jj += mat_ld)
        {
            d_c[j]          = nz;

            Integer fr      = mat.first_row(j);
            Integer lr      = mat.last_row(j);
            Integer ii      = jj + mat.first_elem_pos(j);
            Integer n_del   = 0;

            if (fr > rf)
            {
                if (fr < rl)
                {
                    n_del = 1+(fr-rf)/rs;

                    if ((fr-rf)%rs == 0)
                        --n_del;

                }
                else if (fr == rl)
                {
                    n_del = ci.r_size - 1;
                }
                else
                {
                    n_del = ci.r_size;
                };
            };

            for (Integer i = fr; i <= lr; ++i, ++ii)
            {
                if (i >= rf && i <= rl && (i-rf)%rs == 0)
                {
                    ++n_del;
                    continue;
                };

                const value_type& tmp = mat_ptr[ii];
                d_r[nz]     = i-n_del;

                mrd::reset_helper(d_x[nz],tmp);
                ++nz;
            };
        };

        d_c[c0] = nz;
        d.add_memory(-1);

        ret = matcl::Matrix(out,true);
        return;
    };
};

template<class value_type>
void del_cols_banded_functor<value_type>::eval(matcl::Matrix& ret, const BM& mat,
                         const matcl::details::colon_info& ci, bool)
{
    Integer s = ci.rows();
    if (s == 0)
    {
        ret = matcl::Matrix(mat,false);
        return;
    };

    Integer r0  = mat.rows();
    Integer c0  = mat.cols();

    if (ci.r_flag == 0)
    {		
        if (ci.is_double_mat_colon() == true)
            throw error::invalid_colon_too_many_mat();

        raw::integer_dense ritmp    = ci.get_rim_1().copy();
        const Integer* ptr_r        = ritmp.ptr();

        ritmp.get_struct().reset();

        matcl::utils::sort_q(ritmp.ptr(),s);

        Integer ce = ptr_r[0];
        error::check_col(ce, r0, c0);
        Integer ne = 1;

        for (Integer i = 1; i < s; ++i)
        {
            Integer ce2 = ptr_r[i];
            error::check_col(ce2, r0, c0);

            if (ce2 != ce)
                ++ne;

            ce = ce2;
        }

        if (!r0)
        {
            ret = matcl::Matrix(SM(mat.get_type(),0,c0-ne),false);
            return;
        };

        if (ne == c0)
        {
            ret = matcl::Matrix(SM(mat.get_type(),r0,0),false);
            return;
        };

        Integer c1 = c0 - ne;		
        Integer nnz = mat.nnz();
        SM out(mat.get_type(),r0,c1,nnz);

        raw::details::sparse_ccs<value_type>& d = out.rep();

        Integer* d_c        = d.ptr_c();
        Integer* d_r        = d.ptr_r();
        value_type* d_x     = d.ptr_x();
        const value_type* mat_ptr = mat.rep_ptr();
        Integer mat_ld       = mat.ld();

        Integer nz = 0;
        for (Integer j = 0, k = 0, col = 0; j < c0; ++j)
        {
            if (k < s && j == ptr_r[k]-1)
            {
                while (k < s && j == ptr_r[k]-1)
                    ++k;
            }
            else
            {
                d_c[col] = nz;
                Integer fr = mat.first_row(j);
                Integer lr = mat.last_row(j);
                Integer ii = mat.first_elem_pos(j);

                for (Integer i = fr; i <= lr; ++i, ++ii)
                {			
                    const value_type& tmp = mat_ptr[ii];
                    d_r[nz] = i;
                    mrd::reset_helper(d_x[nz],tmp);
                    ++nz;
                };
                ++col;
            };

            mat_ptr += mat_ld;
        };

        d_c[c1] = nz;
        d.add_memory(-1);

        ret = matcl::Matrix(out,true);
        return;
    }
    else
    {
        Integer ne = s;

        if (!r0)
        {
            ret = matcl::Matrix(SM(mat.get_type(),0,c0-ne),false);
            return;
        };

        if (ne == c0)
        {
            ret = matcl::Matrix(SM(mat.get_type(),r0,0),false);
            return;
        };

        Integer nnz = mat.nnz();
        Integer c1 = c0 - ne;
        SM out(mat.get_type(),r0,c1,nnz);

        raw::details::sparse_ccs<value_type>& d = out.rep();
        Integer* d_c    = d.ptr_c();
        Integer* d_r    = d.ptr_r();
        value_type* d_x = d.ptr_x();
        const value_type* mat_ptr = mat.rep_ptr();

        Integer cs = ci.r_step, cf, cl;

        if (cs < 0)
        {
            cs = -cs;
            cf = ci.r_end;
            cl = ci.r_start;
        }
        else
        {
            cf = ci.r_start;		
            cl = ci.r_end;
        };

        //make it 0-based
        --cf;
        --cl;

        Integer jj = 0, nz = 0, col = 0;
        Integer mat_ld   = mat.ld();

        for (Integer j = 0; j < cf; ++j, jj += mat_ld)
        {
            d_c[col] = nz;

            Integer fr = mat.first_row(j);
            Integer lr = mat.last_row(j);
            Integer ii = jj + mat.first_elem_pos(j);

            for (Integer i = fr; i <= lr; ++i, ++ii)
            {			
                const value_type& tmp = mat_ptr[ii];
                d_r[nz] = i;
                mrd::reset_helper(d_x[nz],tmp);
                ++nz;
            };

            ++col;
        };

        if (cs == 1)
        {
            jj += imult(ne,mat_ld);
        }
        else
        {
            for (Integer j = cf; j <= cl; ++j, jj += mat_ld)
            {
                d_c[col] = nz;

                if (((j-cf) % cs) == 0)
                    continue;

                Integer fr = mat.first_row(j);
                Integer lr = mat.last_row(j);
                Integer ii = jj + mat.first_elem_pos(j);

                for (Integer i = fr; i <= lr; ++i, ++ii)
                {			
                    const value_type& tmp = mat_ptr[ii];
                    d_r[nz] = i;
                    mrd::reset_helper(d_x[nz],tmp);
                    ++nz;
                };

                ++col;
            };
        };
        for (Integer j = cl+1; j < c0; ++j, jj += mat_ld)
        {
            d_c[col] = nz;

            Integer fr = mat.first_row(j);
            Integer lr = mat.last_row(j);
            Integer ii = jj + mat.first_elem_pos(j);

            for (Integer i = fr; i <= lr; ++i, ++ii)
            {			
                const value_type& tmp = mat_ptr[ii];
                d_r[nz] = i;
                mrd::reset_helper(d_x[nz],tmp);
                ++nz;
            };

            ++col;
        };

        d_c[c1] = nz;
        d.add_memory(-1);

        ret = matcl::Matrix(out,true);
        return;
    };
};

template<class value_type>
void del_rowscols_banded_functor<value_type>::eval(matcl::Matrix& ret, const BM& mat,
                         const matcl::details::colon_info& ci, bool rvalue)
{
    Integer sr = ci.rows();
    Integer sc = ci.cols();

    if (sr == 0 && sc == 0)
    {
        ret = matcl::Matrix(mat,false);
        return;
    }
    else if (sr == 0)
    {
        return del_cols_banded_functor<value_type>::eval(ret,mat,ci,rvalue);
    }
    else if (sc == 0)
    {
        return del_rows_banded_functor<value_type>::eval(ret,mat,ci,rvalue);
    }

    if (ci.r_flag == 0 && ci.c_flag == 0)
        return eval_00(ret,mat,ci);
    else if (ci.r_flag == 0 && ci.c_flag == 1)
        return eval_01(ret,mat,ci);
    else if (ci.r_flag == 1 && ci.c_flag == 0)
        return eval_10(ret,mat,ci);
    else
        return eval_11(ret,mat,ci);
};

template<class value_type>
void del_rowscols_banded_functor<value_type>::eval_00(matcl::Matrix& ret, const BM& mat,
                         const matcl::details::colon_info& ci)
{
    Integer sr  = ci.rows();
    Integer sc  = ci.cols();
    Integer r0  = mat.rows();
    Integer c0  = mat.cols();

    //column colon
    raw::integer_dense citmp    = ci.get_cim_2().copy();
    const Integer* ptr_c        = citmp.ptr();

    citmp.get_struct().reset();

    matcl::utils::sort_q(citmp.ptr(),sc);

    Integer cec = ptr_c[0];
    error::check_col(cec, r0, c0);
    Integer nec = 1;

    for (Integer i = 1; i < sc; ++i)
    {
        Integer ce2 = ptr_c[i];
        error::check_col(ce2, r0, c0);

        if (ce2 != cec)
            ++nec;

        cec = ce2;
    }

    //row colon
    raw::integer_dense ritmp    = ci.get_rim_2();
    const Integer* ptr_r        = ritmp.ptr();
    Integer cer                 = ptr_r[0];

    error::check_row(cer, r0, c0);

    matcl::pod_workspace<Integer> iwork(r0,0);

    for (Integer i = 0; i < sr; ++i)
    {
        Integer ce2     = ptr_r[i];
        error::check_row(ce2, r0, c0);
        iwork[ce2-1]    = -1;
    };

    Integer n_del       = 0;
    for (Integer i = 0; i < r0; ++i)
    {
        if (iwork[i] == -1)
            ++n_del;
        else
            iwork[i]    = i-n_del;
    };

    Integer c1 = c0 - nec;
    Integer r1 = r0 - n_del;		

    if (r1 == 0 || c1 == 0)
    {
        ret = matcl::Matrix(SM(mat.get_type(),r1,c1),false);
        return;
    };
    
    Integer nnz = mat.nnz();
    SM out(mat.get_type(),r1,c1,nnz);

    raw::details::sparse_ccs<value_type>& d = out.rep();
    Integer* d_c        = d.ptr_c();
    Integer* d_r        = d.ptr_r();
    value_type* d_x     = d.ptr_x();
    const value_type* mat_ptr = mat.rep_ptr();
    Integer mat_ld       = mat.ld();

    Integer nz          = 0;

    for (Integer j = 0, k = 0, col = 0; j < c0; ++j)
    {
        if (k < sc && j == ptr_c[k]-1)
        {
            while (k < sc && j == ptr_c[k]-1)
                ++k;
        }
        else
        {
            d_c[col]    = nz;
            Integer fr  = mat.first_row(j);
            Integer lr  = mat.last_row(j);
            Integer ii  = mat.first_elem_pos(j);

            for (Integer i = fr; i <= lr; ++i,++ii)
            {
                Integer np  = iwork[i];
                if (np != -1)
                {
                    const value_type& tmp = mat_ptr[ii];
                    d_r[nz] = np;
                    mrd::reset_helper(d_x[nz],tmp);
                    ++nz;
                };
            };

            ++col;
        };

        mat_ptr += mat_ld;
    };

    d_c[c1] = nz;
    d.add_memory(-1);

    ret = matcl::Matrix(out,true);
    return;
}

template<class value_type>
void del_rowscols_banded_functor<value_type>::eval_10(matcl::Matrix& ret, const BM& mat,
                         const matcl::details::colon_info& ci)
{
    Integer sr  = ci.rows();
    Integer sc  = ci.cols();
    Integer r0  = mat.rows();
    Integer c0  = mat.cols();

    raw::integer_dense citmp    = ci.get_cim_2().copy();
    const Integer* ptr_c        = citmp.ptr();

    citmp.get_struct().reset();

    matcl::utils::sort_q(citmp.ptr(),sc);

    Integer cec = ptr_c[0];
    error::check_col(cec, r0, c0);
    Integer nec = 1;

    for (Integer i = 1; i < sc; ++i)
    {
        Integer ce2 = ptr_c[i];
        error::check_col(ce2, r0, c0);

        if (ce2 != cec)
            ++nec;

        cec = ce2;
    }

    Integer ner     = sr;

    Integer c1 = c0 - nec;
    Integer r1 = r0 - ner;		

    if (c1 == 0 || r1 == 0)
    {
        ret = matcl::Matrix(SM(mat.get_type(),r1,c1),false);
        return;
    };

    //row colon
    Integer rs = ci.r_step, rf, rl;

    if (rs < 0)
    {
        rs  = -rs;
        rf  = ci.r_end;
        rl  = ci.r_start;
    }
    else
    {
        rf  = ci.r_start;		
        rl  = ci.r_end;
    };

    //make it 0-based
    rf  = rf - 1;
    rl  = rl - 1;

    Integer nnz = mat.nnz();
    SM out(mat.get_type(),r1,c1,nnz);

    raw::details::sparse_ccs<value_type>& d = out.rep();
    Integer* d_c        = d.ptr_c();
    Integer* d_r        = d.ptr_r();
    value_type* d_x     = d.ptr_x();
    const value_type* mat_ptr = mat.rep_ptr();
    Integer mat_ld      = mat.ld();
    Integer nz          = 0;

    for (Integer j = 0, k = 0, col = 0; j < c0; ++j)
    {
        if (k < sc && j == ptr_c[k]-1)
        {
            while (k < sc && j == ptr_c[k]-1)
                ++k;
        }
        else
        {
            d_c[col]        = nz;
            Integer fr      = mat.first_row(j);
            Integer lr      = mat.last_row(j);
            Integer ii      = mat.first_elem_pos(j);
            Integer n_del   = 0;

            if (fr > rf)
            {
                if (fr < rl)
                {
                    n_del   = 1+(fr-rf)/rs;

                    if ((fr-rf)%rs == 0)
                        --n_del;
                }
                else if (fr == rl)
                {
                    n_del = sr - 1;
                }
                else
                {
                    n_del = sr;
                };
            };

            for (Integer i = fr; i <= lr; ++i, ++ii)
            {
                if (i >= rf && i <= rl && (i-rf)%rs == 0)
                {
                    ++n_del;
                    continue;
                };

                const value_type& tmp = mat_ptr[ii];
                d_r[nz]     = i-n_del;
                mrd::reset_helper(d_x[nz],tmp);
                ++nz;
            };
            ++col;
        };

        mat_ptr += mat_ld;
    };

    d_c[c1] = nz;
    d.add_memory(-1);

    ret = matcl::Matrix(out,true);
    return;
}

template<class value_type>
void del_rowscols_banded_functor<value_type>::eval_01(matcl::Matrix& ret, const BM& mat,
                         const matcl::details::colon_info& ci)
{
    Integer sr  = ci.rows();
    Integer sc  = ci.cols();
    Integer r0  = mat.rows();
    Integer c0  = mat.cols();

    //row colon
    raw::integer_dense ritmp    = ci.get_rim_2();
    const Integer* ptr_r        = ritmp.ptr();
    Integer cer                 = ptr_r[0];

    error::check_row(cer, r0, c0);

    matcl::pod_workspace<Integer> iwork(r0,0);

    for (Integer i = 0; i < sr; ++i)
    {
        Integer ce2     = ptr_r[i];
        error::check_row(ce2, r0, c0);
        iwork[ce2-1]    = -1;
    };

    Integer n_del       = 0;
    for (Integer i = 0; i < r0; ++i)
    {
        if (iwork[i] == -1)
            ++n_del;
        else
            iwork[i]    = i-n_del;
    };

    Integer nec = sc;

    Integer r1 = r0 - n_del;
    Integer c1 = c0 - nec;

    if (r1 == 0 || c1 == 0)
    {
        ret = matcl::Matrix(SM(mat.get_type(),r1,c1),false);
        return;
    };

    Integer nnz = mat.nnz();
    
    SM out(mat.get_type(),r1,c1,nnz);

    raw::details::sparse_ccs<value_type>& d = out.rep();
    Integer* d_c    = d.ptr_c();
    Integer* d_r    = d.ptr_r();
    value_type* d_x = d.ptr_x();
    const value_type* mat_ptr = mat.rep_ptr();

    //column colon
    Integer cs = ci.c_step, cf, cl;
    if (cs < 0)
    {
        cs = -cs;
        cf = ci.c_end;
        cl = ci.c_start;
    }
    else
    {
        cf = ci.c_start;		
        cl = ci.c_end;
    };

    //make it 0-based
    --cf;
    --cl;

    Integer jj = 0, nz = 0, col = 0;
    Integer mat_ld   = mat.ld();

    for (Integer j = 0; j < cf; ++j, jj += mat_ld)
    {
        d_c[col] = nz;

        Integer fr = mat.first_row(j);
        Integer lr = mat.last_row(j);
        Integer ii = jj + mat.first_elem_pos(j);

        for (Integer i = fr; i <= lr; ++i,++ii)
        {
            Integer np  = iwork[i];

            if (np != -1)
            {
                const value_type& tmp = mat_ptr[ii];
                d_r[nz] = np;
                mrd::reset_helper(d_x[nz],tmp);
                ++nz;
            };
        };

        ++col;
    };

    if (cs == 1)
    {
        jj += imult(nec,mat_ld);
    }
    else
    {
        for (Integer j = cf; j <= cl; ++j, jj += mat_ld)
        {
            d_c[col] = nz;

            if (((j-cf) % cs) == 0)
                continue;

            Integer fr = mat.first_row(j);
            Integer lr = mat.last_row(j);
            Integer ii = jj + mat.first_elem_pos(j);

            for (Integer i = fr; i <= lr; ++i,++ii)
            {
                Integer np  = iwork[i];

                if (np != -1)
                {
                    const value_type& tmp = mat_ptr[ii];
                    d_r[nz] = np;
                    mrd::reset_helper(d_x[nz],tmp);
                    ++nz;
                };
            };

            ++col;
        };
    };

    for (Integer j = cl+1; j < c0; ++j, jj += mat_ld)
    {
        d_c[col] = nz;

        Integer fr = mat.first_row(j);
        Integer lr = mat.last_row(j);
        Integer ii = jj + mat.first_elem_pos(j);

        for (Integer i = fr; i <= lr; ++i,++ii)
        {
            Integer np  = iwork[i];
            if (np != -1)
            {
                const value_type& tmp = mat_ptr[ii];
                d_r[nz] = np;
                mrd::reset_helper(d_x[nz],tmp);
                ++nz;
            };
        };

        ++col;
    };

    d_c[c1] = nz;
    d.add_memory(-1);

    ret = matcl::Matrix(out,true);
    return;
};

template<class value_type>
void del_rowscols_banded_functor<value_type>::eval_11(matcl::Matrix& ret, const BM& mat,
                         const matcl::details::colon_info& ci)
{
    Integer sr  = ci.rows();
    Integer sc  = ci.cols();
    Integer r0  = mat.rows();
    Integer c0  = mat.cols();
    
    Integer ner = sr;
    Integer nec = sc;

    Integer r1 = r0 - ner;
    Integer c1 = c0 - nec;

    if (r1 == 0 || c1 == 0)
    {
        ret = matcl::Matrix(SM(mat.get_type(),r1,c1),false);
        return;
    };

    Integer nnz = mat.nnz();    
    SM out(mat.get_type(),r1,c1,nnz);

    raw::details::sparse_ccs<value_type>& d = out.rep();
    Integer* d_c    = d.ptr_c();
    Integer* d_r    = d.ptr_r();
    value_type* d_x = d.ptr_x();
    const value_type* mat_ptr = mat.rep_ptr();

    //column colon
    Integer cs = ci.c_step, cf, cl;
    if (cs < 0)
    {
        cs = -cs;
        cf = ci.c_end;
        cl = ci.c_start;
    }
    else
    {
        cf = ci.c_start;		
        cl = ci.c_end;
    };

    //make it 0-based
    --cf;
    --cl;

    //row colon
    Integer rs = ci.r_step, rf, rl;

    if (rs < 0)
    {
        rs  = -rs;
        rf  = ci.r_end;
        rl  = ci.r_start;
    }
    else
    {
        rf  = ci.r_start;		
        rl  = ci.r_end;
    };

    //make it 0-based
    rf  = rf - 1;
    rl  = rl - 1;

    Integer jj = 0, nz = 0, col = 0;
    Integer mat_ld   = mat.ld();

    for (Integer j = 0; j < cf; ++j, jj += mat_ld)
    {
        d_c[col] = nz;

        Integer fr      = mat.first_row(j);
        Integer lr      = mat.last_row(j);
        Integer ii      = jj + mat.first_elem_pos(j);
        Integer n_del   = 0;

        if (fr > rf)
        {
            if (fr < rl)
            {
                n_del = 1+(fr-rf)/rs;

                if ((fr-rf)%rs == 0)
                    --n_del;

            }
            else if (fr == rl)
            {
                n_del = sr - 1;
            }
            else
            {
                n_del = sr;
            };
        };

        for (Integer i = fr; i <= lr; ++i, ++ii)
        {
            if (i >= rf && i <= rl && (i-rf)%rs == 0)
            {
                ++n_del;
                continue;
            };

            const value_type& tmp = mat_ptr[ii];
            d_r[nz]     = i-n_del;
            mrd::reset_helper(d_x[nz],tmp);
            ++nz;
        };

        ++col;
    };
    if (cs == 1)
    {
        jj += imult(nec,mat_ld);
    }
    else
    {
        for (Integer j = cf; j <= cl; ++j, jj += mat_ld)
        {
            d_c[col] = nz;

            if (((j-cf) % cs) == 0)
                continue;

            Integer fr      = mat.first_row(j);
            Integer lr      = mat.last_row(j);
            Integer ii      = jj + mat.first_elem_pos(j);
            Integer n_del   = 0;

            if (fr > rf)
            {
                if (fr < rl)
                {
                    n_del = 1+(fr-rf)/rs;

                    if ((fr-rf)%rs == 0)
                        --n_del;

                }
                else if (fr == rl)
                {
                    n_del = sr - 1;
                }
                else
                {
                    n_del = sr;
                };
            };

            for (Integer i = fr; i <= lr; ++i, ++ii)
            {
                if (i >= rf && i <= rl && (i-rf)%rs == 0)
                {
                    ++n_del;
                    continue;
                };
                
                const value_type& tmp = mat_ptr[ii];
                d_r[nz]     = i-n_del;
                mrd::reset_helper(d_x[nz],tmp);
                ++nz;
            };

            ++col;
        };
    };

    for (Integer j = cl+1; j < c0; ++j, jj += mat_ld)
    {
        d_c[col] = nz;

        Integer fr      = mat.first_row(j);
        Integer lr      = mat.last_row(j);
        Integer ii      = jj + mat.first_elem_pos(j);
        Integer n_del   = 0;

        if (fr > rf)
        {
            if (fr < rl)
            {
                n_del = 1+(fr-rf)/rs;

                if ((fr-rf)%rs == 0)
                    --n_del;

            }
            else if (fr == rl)
            {
                n_del = sr - 1;
            }
            else
            {
                n_del = sr;
            };
        };

        for (Integer i = fr; i <= lr; ++i, ++ii)
        {
            if (i >= rf && i <= rl && (i-rf)%rs == 0)
            {
                ++n_del;
                continue;
            };

            const value_type& tmp = mat_ptr[ii];
            d_r[nz]     = i-n_del;
            mrd::reset_helper(d_x[nz],tmp);
            ++nz;
        };

        ++col;
    };

    d_c[c1] = nz;
    d.add_memory(-1);

    ret = matcl::Matrix(out,true);
    return;
};

template<class V>
raw::Matrix<V,struct_banded> band_change_diag_functor<V>::eval(const BM& mat, Integer d, const DM& val)
{    
    Integer r   = mat.rows();
    Integer c   = mat.cols();
    Integer f_d = mat.first_diag();
    Integer l_d = mat.last_diag();

    error::check_diag(d,r,c);

    Integer s   = mat.diag_length(d);

    error::check_assign_1(s,val.size(),1);

    if (s == 0)
        return mat;

    BM A(mat.get_type());

    if (d > l_d)
    {
        l_d = d;
        A.assign_to_fresh(mat.resize(mat.rows(),mat.cols(),f_d,l_d));
    }
    else if (d < f_d)
    {
        f_d = d;
        A.assign_to_fresh(mat.resize(mat.rows(),mat.cols(),f_d,l_d));
    }
    else
    {
        A.assign_to_fresh(mat);
    };

    DM val_tmp          = val.make_explicit();
    const V* val_ptr    = val_tmp.ptr();

    V * ptr_this    = A.rep_ptr() + A.first_elem_diag(d);
    Integer A_ld    = A.ld();

    for (Integer i = 0; i < s; ++i)
    {
        mrd::assign_helper(*(ptr_this), (*val_ptr));
        ptr_this += A_ld;
        val_ptr  += 1;
    };

    A.set_struct(md::predefined_struct::get_set_diag(mat.get_struct(), d,value_struct_class::vc_general,
                              is_real_matrix(A), r == c));
    return A;
};

template<class V>
raw::Matrix<V,struct_banded> band_change_diag_functor<V>::eval(const BM& mat, Integer d, const V& val)
{
    using Mat_B     = raw::Matrix<V,struct_banded>;

    Integer r       = mat.rows();
    Integer c       = mat.cols();
    Integer f_d     = mat.first_diag();
    Integer l_d     = mat.last_diag();

    error::check_diag(d,r,c);
    
    Integer s       = mat.diag_length(d);

    if (s == 0)
        return mat;

    bool tz = matcl::details::has_trivial_assignment<BM,V>::eval(mat,val);

    BM A(mat.get_type());

    if (d > l_d)
    {
        l_d = d;
        A.assign_to_fresh(mat.resize(mat.rows(),mat.cols(),f_d,l_d));
    }
    else if (d < f_d)
    {
        f_d = d;
        A.assign_to_fresh(mat.resize(mat.rows(),mat.cols(),f_d,l_d));
    }
    else
    {
        A.assign_to_fresh(mat);
    };

    Integer lda  = A.ld();
    V * ptr_this = A.rep_ptr() + A.first_elem_diag(d);

    for (Integer i = 0; i < s; ++i)
    {
        mrd::assign_helper(*(ptr_this), val);
        ptr_this += lda;
    };

    value_struct_class vt = md::predefined_struct::get_value_type(val,tz);
    A.set_struct(md::predefined_struct::get_set_diag(mat.get_struct(), d,vt, 
                        is_real_matrix(A), r == c));
    return A;
};


};};};

template struct matcl::algorithm::details::del_cols_banded_functor<matcl::Integer>;
template struct matcl::algorithm::details::del_cols_banded_functor<matcl::Real>;
template struct matcl::algorithm::details::del_cols_banded_functor<matcl::Float>;
template struct matcl::algorithm::details::del_cols_banded_functor<matcl::Complex>;
template struct matcl::algorithm::details::del_cols_banded_functor<matcl::Float_complex>;
template struct matcl::algorithm::details::del_cols_banded_functor<matcl::Object>;

template struct matcl::algorithm::details::del_rows_banded_functor<matcl::Integer>;
template struct matcl::algorithm::details::del_rows_banded_functor<matcl::Real>;
template struct matcl::algorithm::details::del_rows_banded_functor<matcl::Float>;
template struct matcl::algorithm::details::del_rows_banded_functor<matcl::Complex>;
template struct matcl::algorithm::details::del_rows_banded_functor<matcl::Float_complex>;
template struct matcl::algorithm::details::del_rows_banded_functor<matcl::Object>;

template struct matcl::algorithm::details::del_rowscols_banded_functor<matcl::Integer>;
template struct matcl::algorithm::details::del_rowscols_banded_functor<matcl::Real>;
template struct matcl::algorithm::details::del_rowscols_banded_functor<matcl::Float>;
template struct matcl::algorithm::details::del_rowscols_banded_functor<matcl::Complex>;
template struct matcl::algorithm::details::del_rowscols_banded_functor<matcl::Float_complex>;
template struct matcl::algorithm::details::del_rowscols_banded_functor<matcl::Object>;

template struct matcl::algorithm::details::band_change_diag_functor<matcl::Integer>;
template struct matcl::algorithm::details::band_change_diag_functor<matcl::Real>;
template struct matcl::algorithm::details::band_change_diag_functor<matcl::Float>;
template struct matcl::algorithm::details::band_change_diag_functor<matcl::Complex>;
template struct matcl::algorithm::details::band_change_diag_functor<matcl::Float_complex>;
template struct matcl::algorithm::details::band_change_diag_functor<matcl::Object>;
