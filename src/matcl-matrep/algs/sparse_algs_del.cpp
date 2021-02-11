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
#include "matcl-internals/base/sort.h"
#include "matcl-internals/container/mat_s.h"
#include "matcl-core/utils/workspace.h"
#include "matcl-matrep/matrix/matrix.h"
#include "matcl-scalar/details/scalfunc_helpers.h"

namespace matcl { namespace algorithm { namespace details
{

namespace mr = matcl::raw;
namespace mrd = matcl::raw::details;
namespace md = matcl::details;

template<class SM>
void del_rows_sparse_functor<SM>::eval(matcl::Matrix& ret, const SM& mat,
                                 const matcl::details::colon_info& ci, bool rvalue)
{
    using value_type = typename SM::value_type;

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

        const raw::integer_dense& ritmp = ci.get_rim_1();
        const Integer* ptr_ri           = ritmp.ptr();
        Integer ce                      = ptr_ri[0];

        error::check_row(ce, r0, c0);

        matcl::pod_workspace<Integer> iwork(r0,0);

        for (Integer i = 0; i < s; ++i)
        {
            Integer ce2     = ptr_ri[i];
            error::check_row(ce2, r0, c0);
            iwork[ce2-1]    = -1;
        };

        Integer n_del = 0;
        for (Integer i = 0; i < r0; ++i)
        {
            if (iwork[i] == -1)
                ++n_del;
            else
                iwork[i] = i-n_del;
        };

        if (!c0)
        {
            ret = matcl::Matrix(SM(mat.get_type(),r0-n_del,0),false);
            return;
        };

        if (n_del == r0)
        {
            ret = matcl::Matrix(SM(mat.get_type(),0,c0),false);
            return;
        };

        Integer r1 = r0 - n_del;		

        bool inplace    = (rvalue == true && mat.is_unique() == true);

        SM out(mat.get_type());

        if (inplace == true)
            out.assign_to_fresh(mat);
        else
            out.assign_to_fresh(SM(mat.get_type(),r1,c0,mat.nnz()));

        raw::details::sparse_ccs<value_type>& d = out.rep();
        Integer* d_c        = d.ptr_c();
        Integer* d_r        = d.ptr_r();
        value_type* d_x     = d.ptr_x();

        const raw::details::sparse_ccs<value_type>& Ad = mat.rep();
        const Integer* Ad_c     = Ad.ptr_c();
        const Integer* Ad_r     = Ad.ptr_r();
        const value_type* Ad_x  = Ad.ptr_x();

        Integer nz              = d.offset();

        for (Integer j = 0; j < c0; ++j)
        {            
            Integer i_start     = Ad_c[j];
            Integer i_end       = Ad_c[j+1];

            d_c[j]              = nz;

            for (Integer i = i_start; i < i_end; ++i)
            {
                Integer p   = Ad_r[i];
                Integer np  = iwork[p];

                if (np == -1)
                    continue;

                const value_type& tmp = Ad_x[i];

                d_r[nz]     = np;
                mrd::reset_helper(d_x[nz],tmp);
                ++nz;
            };
        };

        d_c[c0] = nz;
        d.add_memory(-1);

        if (inplace == true)
        {
            out.set_struct(struct_flag());
            out.change_number_rows(r1);
            ret = matcl::Matrix(out,true);
            return;
        };

        ret = matcl::Matrix(out,true);
        return;
    }
    else
    {
        Integer ne = s;
        
        if (!c0)
        {
            ret = matcl::Matrix(SM(mat.get_type(),r0-ne,c0),false);
            return;
        };

        if (ne == r0)
        {
            ret = matcl::Matrix(SM(mat.get_type(),0,c0), false);
            return;
        };

        Integer r1 = r0 - ne;		

        bool inplace    = (rvalue == true && mat.is_unique() == true);

        SM out(mat.get_type());

        if (inplace == true)
            out.assign_to_fresh(mat);
        else
            out.assign_to_fresh(SM(mat.get_type(),r1,c0,mat.nnz()));

        raw::details::sparse_ccs<value_type>& d = out.rep();
        Integer* d_c        = d.ptr_c();
        Integer* d_r        = d.ptr_r();
        value_type* d_x     = d.ptr_x();

        const raw::details::sparse_ccs<value_type>& Ad = mat.rep();
        const Integer* Ad_c     = Ad.ptr_c();
        const Integer* Ad_r     = Ad.ptr_r();
        const value_type* Ad_x  = Ad.ptr_x();

        Integer rs              = ci.r_step, rf, rl;
        Integer r_size          = ci.r_size;
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

        Integer nz              = d.offset();

        if (rs == 1)
        {
            for (Integer j = 0; j < c0; ++j)
            {
                Integer i_start = Ad_c[j];
                Integer i_end   = Ad_c[j+1];
                d_c[j]          = nz;

                for (Integer i = i_start; i < i_end; ++i)
                {
                    Integer p = Ad_r[i];

                    if (p+1 >= rf && p+1 <= rl)
                        continue;

                    const value_type& tmp = Ad_x[i];

                    if (p  < rf - 1)
                        d_r[nz] = p;
                    else
                        d_r[nz] = p-r_size;

                    mrd::reset_helper(d_x[nz],tmp);
                    ++nz;
                };
            };
        }
        else
        {
            for (Integer j = 0; j < c0; ++j)
            {
                Integer i_start = Ad_c[j];
                Integer i_end   = Ad_c[j+1];
                d_c[j]          = nz;

                for (Integer i = i_start; i < i_end; ++i)
                {
                    Integer p = Ad_r[i];

                    if (p+1 >= rf && p+1 <= rl && (p+1-rf)%rs == 0)
                        continue;

                    const value_type& tmp = Ad_x[i];

                    if (p  < rf - 1)
                        d_r[nz] = p;
                    else if (p >= rl - 1) 
                        d_r[nz] = p-r_size;
                    else
                        d_r[nz] = p-(p-rf+1)/rs-1;					

                    mrd::reset_helper(d_x[nz],tmp);
                    ++nz;
                };
            };
        };

        d_c[c0] = nz;
        d.add_memory(-1);

        if (inplace == true)
        {
            out.set_struct(struct_flag());
            out.change_number_rows(r1);
            ret = matcl::Matrix(out,true);
            return;
        };

        ret = matcl::Matrix(out,true);
        return;
    };
};

template<class SM>
void del_cols_sparse_functor<SM>::eval(matcl::Matrix& ret, const SM& mat,
                                 const matcl::details::colon_info& ci, bool rvalue)
{
    using value_type = typename SM::value_type;

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

        raw::integer_dense ritmp = ci.get_rim_1().copy();
        ritmp.get_struct().reset();
        Integer* ptr_ri = ritmp.ptr();

        matcl::utils::sort_q(ptr_ri,s);

        Integer ce = ptr_ri[0];
        error::check_col(ce, r0, c0);
        Integer ne = 1;

        for (Integer i = 1; i < s; ++i)
        {
            Integer ce2 = ptr_ri[i];
            error::check_col(ce2, r0, c0);

            if (ce2 != ce)
                ++ne;

            ce = ce2;
        }

        if (!r0)
        {
            SM out(mat.get_type(),0,c0-ne);
            ret = matcl::Matrix(out,false);
            return;
        };

        if (ne == c0)
        {
            SM out(mat.get_type(),r0,0);
            ret = matcl::Matrix(out,false);
            return;
        };

        Integer c1 = c0 - ne;		
        Integer nnz = mat.nnz();

        bool inplace    = (rvalue == true && mat.is_unique() == true);

        SM out(mat.get_type());

        if (inplace == true)
            out.assign_to_fresh(mat);
        else
            out.assign_to_fresh(SM(mat.get_type(),r0,c1,nnz));

        raw::details::sparse_ccs<value_type>& d = out.rep();
        Integer* d_c    = d.ptr_c();
        Integer* d_r    = d.ptr_r();
        value_type* d_x = d.ptr_x();

        const raw::details::sparse_ccs<value_type>& Ad = mat.rep();
        const Integer* Ad_c     = Ad.ptr_c();
        const Integer* Ad_r     = Ad.ptr_r();
        const value_type* Ad_x  = Ad.ptr_x();

        Integer nz              = d.offset();

        for (Integer j = 0, k = 0, col = 0; j < c0; ++j)
        {
            if (k < s && j == ptr_ri[k]-1)
            {
                while(k < s && j == ptr_ri[k]-1)
                    ++k;
            }
            else
            {
                Integer i_start = Ad_c[j];
                Integer i_end   = Ad_c[j+1];
                d_c[col]        = nz;

                for (Integer i = i_start; i < i_end; ++i)
                {			
                    Integer p               = Ad_r[i];
                    const value_type& tmp   = Ad_x[i];

                    d_r[nz] = p;
                    mrd::reset_helper(d_x[nz],tmp);
                    ++nz;
                };

                ++col;
            };
        };

        d_c[c1] = nz;
        d.add_memory(-1);

        if (inplace == true)
        {
            out.set_struct(struct_flag());
            out.change_number_cols(c1);
            ret = matcl::Matrix(out,true);
            return;
        };

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
            SM out(mat.get_type(),r0,0);
            ret = matcl::Matrix(out,false);
            return;
        };

        Integer nnz = mat.nnz();
        Integer c1 = c0 - ne;

        bool inplace    = (rvalue == true && mat.is_unique() == true);

        SM out(mat.get_type());

        if (inplace == true)
            out.assign_to_fresh(mat);
        else
            out.assign_to_fresh(SM(mat.get_type(),r0,c1,nnz));

        raw::details::sparse_ccs<value_type>& d = out.rep();
        Integer* d_c    = d.ptr_c();
        Integer* d_r    = d.ptr_r();
        value_type* d_x = d.ptr_x();

        const raw::details::sparse_ccs<value_type>& Ad = mat.rep();
        const Integer* Ad_c     = Ad.ptr_c();
        const Integer* Ad_r     = Ad.ptr_r();
        const value_type* Ad_x  = Ad.ptr_x();

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

        Integer col = 0;
        Integer nz  = d.offset();

        if (inplace == false)
        {
            for (Integer j = 1; j < cf; ++j)
            {
                Integer i_start = Ad_c[j-1];
                Integer i_end   = Ad_c[j];
                d_c[col]        = nz;

                for (Integer i = i_start; i < i_end; ++i)
                {			
                    Integer p               = Ad_r[i];
                    const value_type& tmp   = Ad_x[i];

                    d_r[nz] = p;
                    mrd::reset_helper(d_x[nz],tmp);
                    ++nz;
                };

                ++col;
            };
        }
        else
        {
            //just skip columns 
            if (cf > 1)
            {
                col += (cf - 1);
                nz  += (Ad_c[cf-1] - Ad_c[0]);
            }
        };

        if (cs != 1)
        {
            for (Integer j = cf; j <= cl; ++j)
            {
                if (((j-cf) % cs) == 0)
                    continue;

                Integer i_start = Ad_c[j-1];
                Integer i_end   = Ad_c[j];
                d_c[col]        = nz;

                for (Integer i = i_start; i < i_end; ++i)
                {			
                    Integer p               = Ad_r[i];
                    const value_type& tmp   = Ad_x[i];

                    d_r[nz] = p;
                    mrd::reset_helper(d_x[nz],tmp);
                    ++nz;
                };
                ++col;
            };
        };

        if (inplace == true && Ad_c[cl] == nz)
        {
            //no elements have been removed; only update column indices
            for (Integer j = cl; j < c0; ++j)
            {
                d_c[col]        = nz;
                nz              = Ad_c[j+1];
                ++col;
            };

            d_c[c1]             = nz;
        }
        else
        {
            for (Integer j = cl; j < c0; ++j)
            {
                Integer i_start = Ad_c[j];
                Integer i_end   = Ad_c[j+1];
                d_c[col]        = nz;

                for (Integer i = i_start; i < i_end; ++i)
                {			
                    Integer p               = Ad_r[i];
                    const value_type& tmp   = Ad_x[i];

                    d_r[nz] = p;
                    mrd::reset_helper(d_x[nz],tmp);
                    ++nz;
                };
                ++col;
            };

            d_c[c1] = nz;
            d.add_memory(-1);
        };

        if (inplace == true)
        {
            out.set_struct(struct_flag());
            out.change_number_cols(c1);
            ret = matcl::Matrix(out,true);
            return;
        };

        ret = matcl::Matrix(out,true);
        return;
    };
};

template<class SM>
void del_rowscols_sparse_functor<SM>::eval(matcl::Matrix& ret, const SM& mat,
                                 const matcl::details::colon_info& ci, bool rvalue)
{
    using value_type = typename SM::value_type;

    Integer sr = ci.rows();
    Integer sc = ci.cols();

    if (sr == 0 && sc == 0)
    {
        ret = matcl::Matrix(mat,false);
        return;
    }
    else if (sr == 0)
    {
        return del_cols_sparse_functor<SM>::eval(ret,mat,ci,rvalue);
    }
    else if (sc == 0)
    {
        return del_rows_sparse_functor<SM>::eval(ret,mat,ci,rvalue);
    }

    if (ci.r_flag == 0 && ci.c_flag == 0)
        return eval_00(ret,mat,ci,rvalue);
    else if (ci.r_flag == 0 && ci.c_flag == 1)
        return eval_01(ret,mat,ci,rvalue);
    else if (ci.r_flag == 1 && ci.c_flag == 0)
        return eval_10(ret,mat,ci,rvalue);
    else
        return eval_11(ret,mat,ci,rvalue);
};

template<class SM>
void del_rowscols_sparse_functor<SM>::eval_00(matcl::Matrix& ret, const SM& mat,
                                 const matcl::details::colon_info& ci, bool rvalue)
{
    using value_type = typename SM::value_type;

    Integer sr  = ci.rows();
    Integer sc  = ci.cols();
    Integer r0  = mat.rows();
    Integer c0  = mat.cols();

    //column colon
    raw::integer_dense citmp = ci.get_cim_2().copy();
    citmp.get_struct().reset();
    Integer* ptr_ci = citmp.ptr();

    matcl::utils::sort_q(ptr_ci,sc);

    Integer cec = ptr_ci[0];
    error::check_col(cec, r0, c0);
    Integer nec = 1;

    for (Integer i = 1; i < sc; ++i)
    {
        Integer ce2 = ptr_ci[i];
        error::check_col(ce2, r0, c0);

        if (ce2 != cec)
            ++nec;

        cec = ce2;
    }

    //row colon
    const raw::integer_dense& ritmp = ci.get_rim_2();
    const Integer* ptr_ri           = ritmp.ptr();
    Integer cer                     = ptr_ri[0];

    error::check_row(cer, r0, c0);

    matcl::pod_workspace<Integer> iwork(r0,0);
    for (Integer i = 0; i < sr; ++i)
    {
        Integer ce2     = ptr_ri[i];
        error::check_row(ce2, r0, c0);
        iwork[ce2-1]    = -1;
    };

    Integer n_del = 0;
    for (Integer i = 0; i < r0; ++i)
    {
        if (iwork[i] == -1)
            ++n_del;
        else
            iwork[i] = i-n_del;
    };

    Integer r1 = r0 - n_del;
    Integer c1 = c0 - nec;

    if (r1 == 0 || c1 == 0)
    {
        SM out(mat.get_type(),r1,c1);
        ret = matcl::Matrix(out,false);
        return;
    };
    
    Integer nnz     = mat.nnz();

    bool inplace    = (rvalue == true && mat.is_unique() == true);

    SM out(mat.get_type());

    if (inplace == true)
        out.assign_to_fresh(mat);
    else
        out.assign_to_fresh(SM(mat.get_type(),r1,c1,nnz));

    raw::details::sparse_ccs<value_type>& d = out.rep();
    Integer* d_c    = d.ptr_c();
    Integer* d_r    = d.ptr_r();
    value_type* d_x = d.ptr_x();

    const raw::details::sparse_ccs<value_type>& Ad = mat.rep();
    const Integer* Ad_c     = Ad.ptr_c();
    const Integer* Ad_r     = Ad.ptr_r();
    const value_type* Ad_x  = Ad.ptr_x();

    Integer nz              = d.offset();

    for (Integer j = 0, k = 0, col = 0; j < c0; ++j)
    {
        if (k < sc && j == ptr_ci[k]-1)
        {
            while(k < sc && j == ptr_ci[k]-1)
                ++k;
        }
        else
        {
            Integer i_start = Ad_c[j];
            Integer i_end   = Ad_c[j+1];
            d_c[col]        = nz;

            for (Integer i = i_start; i < i_end; ++i)
            {
                Integer p   = Ad_r[i];
                Integer np  = iwork[p];

                if (np == -1)
                    continue;

                const value_type& tmp = Ad_x[i];

                d_r[nz]     = np;
                mrd::reset_helper(d_x[nz],tmp);
                ++nz;
            };

            ++col;
        };
    };

    d_c[c1] = nz;
    d.add_memory(-1);

    if (inplace == true)
    {
        out.set_struct(struct_flag());
        out.change_number_rows(r1);
        out.change_number_cols(c1);
        ret = matcl::Matrix(out,true);
        return;
    };

    ret = matcl::Matrix(out,true);
    return;
};

template<class SM>
void del_rowscols_sparse_functor<SM>::eval_10(matcl::Matrix& ret, const SM& mat,
                                 const matcl::details::colon_info& ci, bool rvalue)
{
    using value_type = typename SM::value_type;

    Integer sr  = ci.rows();
    Integer sc  = ci.cols();
    Integer r0  = mat.rows();
    Integer c0  = mat.cols();

    //column colon
    raw::integer_dense citmp = ci.get_cim_2().copy();
    citmp.get_struct().reset();
    Integer* ptr_ci = citmp.ptr();

    matcl::utils::sort_q(ptr_ci,sc);

    Integer cec = ptr_ci[0];
    error::check_col(cec, r0, c0);
    Integer nec = 1;    

    for (Integer i = 1; i < sc; ++i)
    {
        Integer ce2 = ptr_ci[i];
        error::check_col(ce2, r0, c0);

        if (ce2 != cec)
            ++nec;

        cec = ce2;
    }

    Integer ner = sr;

    Integer r1  = r0 - ner;
    Integer c1  = c0 - nec;

    if (r1 == 0 || c1 == 0)
    {
        SM out(mat.get_type(),r1,c1);
        ret = matcl::Matrix(out,false);
        return;
    };

    Integer nnz     = mat.nnz();
    bool inplace    = (rvalue == true && mat.is_unique() == true);

    SM out(mat.get_type());

    if (inplace == true)
        out.assign_to_fresh(mat);
    else
        out.assign_to_fresh(SM(mat.get_type(),r1,c1,nnz));

    raw::details::sparse_ccs<value_type>& d = out.rep();
    Integer* d_c    = d.ptr_c();
    Integer* d_r    = d.ptr_r();
    value_type* d_x = d.ptr_x();

    const raw::details::sparse_ccs<value_type>& Ad = mat.rep();
    const Integer* Ad_c     = Ad.ptr_c();
    const Integer* Ad_r     = Ad.ptr_r();
    const value_type* Ad_x  = Ad.ptr_x();

    //row colon
    Integer rs              = ci.r_step, rf, rl;
    Integer r_size          = ci.r_size;
    Integer nz              = d.offset();

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

    if (rs == 1)
    {
        for (Integer j = 0, k = 0, col = 0; j < c0; ++j)
        {
            if (k < sc && j == ptr_ci[k]-1)
            {
                while(k < sc && j == ptr_ci[k]-1)
                    ++k;
            }
            else
            {
                Integer i_start = Ad_c[j];
                Integer i_end   = Ad_c[j+1];
                d_c[col]        = nz;

                for (Integer i = i_start; i < i_end; ++i)
                {
                    Integer p = Ad_r[i];

                    if (p+1 >= rf && p+1 <= rl)
                        continue;

                    const value_type& tmp = Ad_x[i];

                    if (p  < rf - 1)
                        d_r[nz] = p;
                    else
                        d_r[nz] = p-r_size;

                    mrd::reset_helper(d_x[nz],tmp);
                    ++nz;
                };

                ++col;
            };
        };
    }
    else
    {
        for (Integer j = 0, k = 0, col = 0; j < c0; ++j)
        {
            if (k < sc && j == ptr_ci[k]-1)
            {
                while(k < sc && j == ptr_ci[k]-1)
                    ++k;
            }
            else
            {
                Integer i_start = Ad_c[j];
                Integer i_end   = Ad_c[j+1];
                d_c[col]        = nz;

                for (Integer i = i_start; i < i_end; ++i)
                {
                    Integer p = Ad_r[i];

                    if (p+1 >= rf && p+1 <= rl && (p+1-rf)%rs == 0)
                        continue;

                    const value_type& tmp = Ad_x[i];

                    if (p  < rf - 1)
                        d_r[nz] = p;
                    else if (p >= rl - 1) 
                        d_r[nz] = p-r_size;
                    else
                        d_r[nz] = p-(p-rf+1)/rs-1;					

                    mrd::reset_helper(d_x[nz],tmp);
                    ++nz;
                };

                ++col;
            };
        };
    };

    d_c[c1] = nz;
    d.add_memory(-1);

    if (inplace == true)
    {
        out.set_struct(struct_flag());
        out.change_number_rows(r1);
        out.change_number_cols(c1);
        ret = matcl::Matrix(out,true);
        return;
    };

    ret = matcl::Matrix(out,true);
    return;
};

template<class SM>
void del_rowscols_sparse_functor<SM>::eval_01(matcl::Matrix& ret, const SM& mat,
                                 const matcl::details::colon_info& ci, bool rvalue)
{
    using value_type = typename SM::value_type;

    Integer sr  = ci.rows();
    Integer sc  = ci.cols();
    Integer r0  = mat.rows();
    Integer c0  = mat.cols();

    Integer nec = sc;
        
    //row colon
    const raw::integer_dense& ritmp = ci.get_rim_2();
    const Integer* ptr_ri           = ritmp.ptr();
    Integer cer                     = ptr_ri[0];

    error::check_row(cer, r0, c0);

    matcl::pod_workspace<Integer> iwork(r0,0);
    for (Integer i = 0; i < sr; ++i)
    {
        Integer ce2     = ptr_ri[i];
        error::check_row(ce2, r0, c0);
        iwork[ce2-1]    = -1;
    };

    Integer n_del = 0;
    for (Integer i = 0; i < r0; ++i)
    {
        if (iwork[i] == -1)
            ++n_del;
        else
            iwork[i] = i-n_del;
    };

    Integer r1 = r0 - n_del;
    Integer c1 = c0 - nec;

    if (r1 == 0 || c1 == 0)
    {
        ret = matcl::Matrix(SM(mat.get_type(),r1,c1),false);
        return;
    };

    Integer nnz     = mat.nnz();    
    bool inplace    = (rvalue == true && mat.is_unique() == true);

    SM out(mat.get_type());

    if (inplace == true)
        out.assign_to_fresh(mat);
    else
        out.assign_to_fresh(SM(mat.get_type(),r1,c1,nnz));

    raw::details::sparse_ccs<value_type>& d = out.rep();
    Integer* d_c    = d.ptr_c();
    Integer* d_r    = d.ptr_r();
    value_type* d_x = d.ptr_x();

    const raw::details::sparse_ccs<value_type>& Ad = mat.rep();
    const Integer* Ad_c     = Ad.ptr_c();
    const Integer* Ad_r     = Ad.ptr_r();
    const value_type* Ad_x  = Ad.ptr_x();

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

    Integer col = 0;
    Integer nz  = d.offset();

    for (Integer j = 1; j < cf; ++j)
    {
        Integer i_start = Ad_c[j-1];
        Integer i_end   = Ad_c[j];
        d_c[col]        = nz;

        for (Integer i = i_start; i < i_end; ++i)
        {
            Integer p   = Ad_r[i];
            Integer np  = iwork[p];

            if (np == -1)
                continue;

            const value_type& tmp = Ad_x[i];

            d_r[nz]     = np;
            mrd::reset_helper(d_x[nz],tmp);
            ++nz;
        };

        ++col;
    };

    if (cs != 1)
    {
        for (Integer j = cf; j <= cl; ++j)
        {
            if (((j-cf) % cs) == 0)
                continue;

            Integer i_start = Ad_c[j-1];
            Integer i_end   = Ad_c[j];
            d_c[col]        = nz;

            for (Integer i = i_start; i < i_end; ++i)
            {
                Integer p   = Ad_r[i];
                Integer np  = iwork[p];

                if (np == -1)
                    continue;

                const value_type& tmp = Ad_x[i];

                d_r[nz]     = np;
                mrd::reset_helper(d_x[nz],tmp);
                ++nz;
            };

            ++col;
        };
    };

    for (Integer j = cl; j < c0; ++j)
    {
        Integer i_start = Ad_c[j];
        Integer i_end   = Ad_c[j+1];
        d_c[col]        = nz;

        for (Integer i = i_start; i < i_end; ++i)
        {
            Integer p   = Ad_r[i];
            Integer np  = iwork[p];

            if (np == -1)
                continue;

            const value_type& tmp = Ad_x[i];

            d_r[nz]     = np;
            mrd::reset_helper(d_x[nz],tmp);
            ++nz;
        };

        ++col;
    };

    d_c[c1] = nz;
    d.add_memory(-1);

    if (inplace == true)
    {
        out.set_struct(struct_flag());
        out.change_number_rows(r1);
        out.change_number_cols(c1);
        ret = matcl::Matrix(out,true);
        return;
    };

    ret = matcl::Matrix(out,true);
    return;
};

template<class SM>
void del_rowscols_sparse_functor<SM>::eval_11(matcl::Matrix& ret, const SM& mat,
                                 const matcl::details::colon_info& ci, bool rvalue)
{
    using value_type = typename SM::value_type;
    
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

    Integer nnz     = mat.nnz();    
    bool inplace    = (rvalue == true && mat.is_unique() == true);

    SM out(mat.get_type());

    if (inplace == true)
        out.assign_to_fresh(mat);
    else
        out.assign_to_fresh(SM(mat.get_type(),r1,c1,nnz));

    raw::details::sparse_ccs<value_type>& d = out.rep();
    Integer* d_c    = d.ptr_c();
    Integer* d_r    = d.ptr_r();
    value_type* d_x = d.ptr_x();

    const raw::details::sparse_ccs<value_type>& Ad = mat.rep();
    const Integer* Ad_c     = Ad.ptr_c();
    const Integer* Ad_r     = Ad.ptr_r();
    const value_type* Ad_x  = Ad.ptr_x();

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

    //row colon
    Integer rs              = ci.r_step, rf, rl;
    Integer r_size          = ci.r_size;
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

    Integer col = 0;
    Integer nz  = d.offset();

    if (rs == 1)
    {
        for (Integer j = 1; j < cf; ++j)
        {
            Integer i_start = Ad_c[j-1];
            Integer i_end   = Ad_c[j];
            d_c[col]        = nz;

            for (Integer i = i_start; i < i_end; ++i)
            {
                Integer p = Ad_r[i];

                if (p+1 >= rf && p+1 <= rl)
                    continue;

                const value_type& tmp = Ad_x[i];

                if (p  < rf - 1)
                    d_r[nz] = p;
                else
                    d_r[nz] = p-r_size;

                mrd::reset_helper(d_x[nz],tmp);
                ++nz;
            };

            ++col;
        };
    }
    else
    {
        for (Integer j = 1; j < cf; ++j)
        {
            Integer i_start = Ad_c[j-1];
            Integer i_end   = Ad_c[j];
            d_c[col]        = nz;

            for (Integer i = i_start; i < i_end; ++i)
            {
                Integer p = Ad_r[i];

                if (p+1 >= rf && p+1 <= rl && (p+1-rf)%rs == 0)
                    continue;

                const value_type& tmp = Ad_x[i];

                if (p  < rf - 1)
                    d_r[nz] = p;
                else if (p >= rl - 1) 
                    d_r[nz] = p-r_size;
                else
                    d_r[nz] = p-(p-rf+1)/rs-1;					

                mrd::reset_helper(d_x[nz],tmp);
                ++nz;
            };

            ++col;
        };
    };

    if (cs != 1)
    {
        if (rs == 1)
        {
            for (Integer j = cf; j <= cl; ++j)
            {
                if (((j-cf) % cs) == 0)
                    continue;

                Integer i_start = Ad_c[j-1];
                Integer i_end   = Ad_c[j];
                d_c[col]        = nz;

                for (Integer i = i_start; i < i_end; ++i)
                {
                    Integer p = Ad_r[i];

                    if (p+1 >= rf && p+1 <= rl)
                        continue;

                    const value_type& tmp = Ad_x[i];

                    if (p  < rf - 1)
                        d_r[nz] = p;
                    else
                        d_r[nz] = p-r_size;

                    mrd::reset_helper(d_x[nz],tmp);
                    ++nz;
                };

                ++col;
            };
        }
        else
        {
            for (Integer j = cf; j <= cl; ++j)
            {
                if (((j-cf) % cs) == 0)
                    continue;

                Integer i_start = Ad_c[j-1];
                Integer i_end   = Ad_c[j];
                d_c[col]        = nz;

                for (Integer i = i_start; i < i_end; ++i)
                {
                    Integer p = Ad_r[i];

                    if (p+1 >= rf && p+1 <= rl && (p+1-rf)%rs == 0)
                        continue;

                    const value_type& tmp = Ad_x[i];

                    if (p  < rf - 1)
                        d_r[nz] = p;
                    else if (p >= rl - 1) 
                        d_r[nz] = p-r_size;
                    else
                        d_r[nz] = p-(p-rf+1)/rs-1;					

                    mrd::reset_helper(d_x[nz],tmp);
                    ++nz;
                };

                ++col;
            };
        };
    };

    if (rs == 1)
    {
        for (Integer j = cl; j < c0; ++j)
        {
            Integer i_start = Ad_c[j];
            Integer i_end   = Ad_c[j+1];
            d_c[col]        = nz;

            for (Integer i = i_start; i < i_end; ++i)
            {
                Integer p = Ad_r[i];

                if (p+1 >= rf && p+1 <= rl)
                    continue;

                const value_type& tmp = Ad_x[i];

                if (p  < rf - 1)
                    d_r[nz] = p;
                else
                    d_r[nz] = p-r_size;

                mrd::reset_helper(d_x[nz],tmp);
                ++nz;
            };

            ++col;
        };
    }
    else
    {
        for (Integer j = cl; j < c0; ++j)
        {
            Integer i_start = Ad_c[j];
            Integer i_end   = Ad_c[j+1];
            d_c[col]        = nz;

            for (Integer i = i_start; i < i_end; ++i)
            {
                Integer p = Ad_r[i];

                if (p+1 >= rf && p+1 <= rl && (p+1-rf)%rs == 0)
                    continue;

                const value_type& tmp = Ad_x[i];

                if (p  < rf - 1)
                    d_r[nz] = p;
                else if (p >= rl - 1) 
                    d_r[nz] = p-r_size;
                else
                    d_r[nz] = p-(p-rf+1)/rs-1;					

                mrd::reset_helper(d_x[nz],tmp);
                ++nz;
            };

            ++col;
        };
    };

    d_c[c1] = nz;
    d.add_memory(-1);

    if (inplace == true)
    {
        out.set_struct(struct_flag());
        out.change_number_rows(r1);
        out.change_number_cols(c1);
        ret = matcl::Matrix(out,true);
        return;
    };

    ret = matcl::Matrix(out,true);
    return;
};

};};};

template struct matcl::algorithm::details::del_cols_sparse_functor<matcl::raw::integer_sparse>;
template struct matcl::algorithm::details::del_cols_sparse_functor<matcl::raw::real_sparse>;
template struct matcl::algorithm::details::del_cols_sparse_functor<matcl::raw::float_sparse>;
template struct matcl::algorithm::details::del_cols_sparse_functor<matcl::raw::complex_sparse>;
template struct matcl::algorithm::details::del_cols_sparse_functor<matcl::raw::float_complex_sparse>;
template struct matcl::algorithm::details::del_cols_sparse_functor<matcl::raw::object_sparse>;

template struct matcl::algorithm::details::del_rows_sparse_functor<matcl::raw::integer_sparse>;
template struct matcl::algorithm::details::del_rows_sparse_functor<matcl::raw::real_sparse>;
template struct matcl::algorithm::details::del_rows_sparse_functor<matcl::raw::float_sparse>;
template struct matcl::algorithm::details::del_rows_sparse_functor<matcl::raw::complex_sparse>;
template struct matcl::algorithm::details::del_rows_sparse_functor<matcl::raw::float_complex_sparse>;
template struct matcl::algorithm::details::del_rows_sparse_functor<matcl::raw::object_sparse>;

template struct matcl::algorithm::details::del_rowscols_sparse_functor<matcl::raw::integer_sparse>;
template struct matcl::algorithm::details::del_rowscols_sparse_functor<matcl::raw::real_sparse>;
template struct matcl::algorithm::details::del_rowscols_sparse_functor<matcl::raw::float_sparse>;
template struct matcl::algorithm::details::del_rowscols_sparse_functor<matcl::raw::complex_sparse>;
template struct matcl::algorithm::details::del_rowscols_sparse_functor<matcl::raw::float_complex_sparse>;
template struct matcl::algorithm::details::del_rowscols_sparse_functor<matcl::raw::object_sparse>;
