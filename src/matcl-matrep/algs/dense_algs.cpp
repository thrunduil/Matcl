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

#include "matcl-matrep/algs/dense_algs.h"
#include "matcl-internals/base/utils.h"
#include "matcl-internals/base/sort.h"
#include "matcl-internals/error/error_check_basic.h"
#include "matcl-scalar/details/scalfunc_helpers.h"
#include "matcl-matrep/details/struct_flag_predefined.h"

namespace matcl { namespace algorithm { namespace details
{

namespace mr = matcl::raw;
namespace mrd = matcl::raw::details;
namespace md = matcl::details;

template<class DM>
void dense_change_entries_functor<DM>::eval(DM& mat,const md::colon_info& ci, const value_type& val)
{
    Integer nr  = ci.rows();
    Integer nc  = ci.cols();

    if (nr == 0 && nc == 0)
        return;

    const Integer* ptr_ci   = ci.get_ci_2_ptr();
    const Integer* ptr_ri   = ci.get_ri_2_ptr();
    Integer mat_ld          = mat.ld();

    using value_type        = typename DM::value_type;

    if (ci.r_flag == 0 && ci.c_flag == 0)
    {			
        for (Integer j = 0; j < nc; ++j)
        {
            Integer pos_c = imult(ptr_ci[j] - 1,mat_ld);

            value_type* ptr_mat = mat.ptr() + pos_c - 1;

            for (Integer i = 0; i < nr; ++i)
                mrd::assign_helper(ptr_mat[ptr_ri[i]], val);
        };

        mat.get_struct().reset();
        return;
    };

    if (ci.r_flag == 0 && ci.c_flag == 1)
    {			
        Integer pos_c   = imult(ci.c_start - 1,mat.ld());
        Integer dpos_c  = imult(ci.c_step,mat.ld());

        value_type* ptr_mat = mat.ptr() + pos_c - 1;

        for (Integer j = 0; j < ci.c_size; ++j)
        {			
            for (Integer i = 0; i < nr; ++i)
                mrd::assign_helper(ptr_mat[ptr_ri[i]], val);

            ptr_mat += dpos_c;
        };

        mat.get_struct().reset();
        return;
    };

    if (ci.r_flag == 1 && ci.c_flag == 0)
    {
        for (Integer j = 0; j < nc; ++j)
        {
            Integer pos_c = imult(ptr_ci[j] - 1, mat_ld);

            value_type* ptr_mat = mat.ptr() + pos_c - 1;

            for (Integer i = 0, ii = ci.r_start; i < ci.r_size; ++i, ii+= ci.r_step)
                mrd::assign_helper(ptr_mat[ii], val);
        };

        mat.get_struct().reset();
        return;
    };

    if (ci.r_flag == 1 && ci.c_flag == 1)
    {
        Integer pos_c   = imult(ci.c_start - 1,mat.ld());
        Integer dpos_c  = imult(ci.c_step,mat.ld());

        value_type* ptr_mat = mat.ptr() + pos_c - 1;

        for (Integer j = 0; j < ci.c_size; ++j)
        {			
            for (Integer i = 0, ii = ci.r_start; i < ci.r_size; ++i, ii+=ci.r_step)
                mrd::assign_helper(ptr_mat[ii], val);

            ptr_mat += dpos_c;
        };

        mat.get_struct().reset();
        return;
    };
};

template<class DM>
void dense_change_entries_functor_2<DM>::eval(DM& mat,const md::colon_info& ci, const value_type& val)
{
    Integer s               = ci.rows();    

    if (s == 0)
        return;

    value_type* ptr_mat = mat.ptr();
    Integer mat_ld      = mat.ld();

    if (ci.r_flag == 0)
    {       
        bool single     = ci.is_double_mat_colon() == false;

        if (single == true && mat.ld() == mat.rows())
        {
            const Integer* ptr_ri   = ci.get_ri_1_ptr();

            for (Integer i = 0; i < s; ++i)
                mrd::assign_helper(ptr_mat[ptr_ri[i]-1], val);
        }
        else if (single == true)
        {
            const Integer* ptr_ri   = ci.get_ri_1_ptr();
            Integer r = mat.rows();

            for (Integer i = 0; i < s; ++i)
            {
                Integer pos = ptr_ri[i];
                Integer rm, cm;
                md::pos2ind(pos, r, rm, cm);
                mrd::assign_helper(ptr_mat[rm+cm*mat_ld], val);
            };
        }
        else
        {
            const Integer* ptr_ri   = ci.get_ri_r_ptr();
            const Integer* ptr_ci   = ci.get_ri_c_ptr();

            for (Integer i = 0; i < s; ++i)
            {
                Integer rm  = ptr_ri[i] - 1;
                Integer cm  = ptr_ci[i] - 1;
                mrd::assign_helper(ptr_mat[rm+cm*mat_ld], val);
            };
        };

        mat.get_struct().reset();
        return;
    }
    else
    {        
        if (mat.ld() == mat.rows())
        {
            for (Integer i = 0, pos = ci.r_start-1; i < ci.r_size; ++i, pos += ci.r_step)
                mrd::assign_helper(ptr_mat[pos], val);
        }
        else
        {
            Integer r = mat.rows();

            for (Integer i = 0, pos = ci.r_start-1; i < ci.r_size; ++i, pos += ci.r_step)
            {
                while (pos >= r)
                {
                    pos     -= r;
                    ptr_mat += mat_ld;
                };

                mrd::assign_helper(ptr_mat[pos], val);
            };
        };

        mat.get_struct().reset();
        return;
    };
};

template<class DM>
void del_rows_dense_functor<DM>::eval(Matrix& ret, const DM& mat,const md::colon_info& ci, bool rvalue)
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

        const mr::integer_dense& ci_ri  = ci.get_rim_1();
        raw::integer_dense ritmp        = ci_ri.copy();
        const Integer* ptr_ri           = ritmp.ptr();

        ritmp.get_struct().reset();

        utils::sort_q(ritmp.ptr(),s);

        Integer ce = ptr_ri[0];
        error::check_row(ce, r0, c0);
        Integer ne = 1;

        for (Integer i = 1; i < s; ++i)
        {
            Integer ce2 = ptr_ri[i];
            error::check_row(ce2, r0, c0);
    
            if (ce2 != ce)
                ++ne;

            ce = ce2;
        }

        if (!c0)
        {
            ret = matcl::Matrix(DM(mat.get_type(),r0-ne,c0), false);
            return;
        };

        if (ne == r0)
        {
            ret =  matcl::Matrix(DM(mat.get_type(),0,c0), false);
            return;
        };

        Integer r1 = r0 - ne;

        DM out(mat.get_type());

        bool inplace    = (rvalue == true && mat.is_unique() == true);

        if (inplace == true)
            out.assign_to_fresh(mat);
        else
            out.assign_to_fresh(DM(mat.get_type(),r1,c0));

        using V = typename DM::value_type;

        V * ptr_out         = out.ptr();
        const V* ptr_mat    = mat.ptr();

        Integer mat_ld      = mat.ld();
        Integer out_ld      = out.ld();

        for (Integer j = 0; j < c0; ++j)
        {
            for (Integer i = 0, ii = 0, k = 0; i < r0; ++i)
            {
                if (k < s && i+1 == ptr_ri[k])
                {
                    while(k < s && i+1 == ptr_ri[k])
                        ++k;
                }
                else
                {
                    mrd::reset_helper(ptr_out[ii],ptr_mat[i]);
                    ++ii;
                };
            };

            ptr_mat     += mat_ld;
            ptr_out     += out_ld;
        };

        if (inplace == true)
        {
            out.set_struct(struct_flag());
            ret = matcl::Matrix(out.resize(r1,c0),true);
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
            DM out(mat.get_type(),r0-ne,0);
            ret = matcl::Matrix(out,false);
            return;
        };

        if (ne == r0)
        {
            DM out(mat.get_type(),0,c0);
            ret = matcl::Matrix(out,false);
            return;
        };

        Integer r1 = r0 - ne;

        bool inplace    = (rvalue == true && mat.is_unique() == true);

        DM out(mat.get_type());

        if (inplace == true)
            out.assign_to_fresh(mat);
        else
            out.assign_to_fresh(DM(mat.get_type(),r1,c0));        

        using V = typename DM::value_type;

        V* ptr_out          = out.ptr();
        const V* ptr_mat    = mat.ptr();

        Integer mat_ld      = mat.ld();
        Integer out_ld      = out.ld();

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

        for (Integer j = 1; j <= c0; ++j)
        {
            Integer ii = 0, iii = 0;
            
            for (Integer i = 1; i < rf; ++i, ++ii)
                mrd::reset_helper(ptr_out[iii++], ptr_mat[ii]);
        
            if (rs == 1)
            {
                ii += ne;
            }
            else
            {
                for (Integer i = rf; i <= rl; ++i, ++ii)
                {
                    if (((i-rf) % rs) != 0)
                        mrd::reset_helper(ptr_out[iii++], ptr_mat[ii]);
                };
            }

            for (Integer i = rl + 1; i <= r0; ++i, ++ii)
                mrd::reset_helper(ptr_out[iii++],ptr_mat[ii]);

            ptr_mat += mat_ld;
            ptr_out += out_ld;
        };

        if (inplace == true)
        {
            out.set_struct(struct_flag());
            ret = matcl::Matrix(out.resize(r1,c0),true);
            return;
        };

        ret = matcl::Matrix(out,true);
        return;
    };
};

template<class DM>
void del_cols_dense_functor<DM>::eval(matcl::Matrix& ret, const DM& mat,const md::colon_info& ci, 
                                      bool rvalue)
{
    Integer s = ci.rows();
    if (s == 0)
    {
        ret = matcl::Matrix(mat, false);
        return;
    };    

    Integer r0  = mat.rows();
    Integer c0  = mat.cols();

    using V = typename DM::value_type;

    const V* ptr_mat = mat.ptr();

    if (ci.r_flag == 0)
    {		
        if (ci.is_double_mat_colon() == true)
            throw error::invalid_colon_too_many_mat();

        const mr::integer_dense& ci_ri  = ci.get_rim_1();
        raw::integer_dense ritmp        = ci_ri.copy();
        const Integer* ptr_ri           = ritmp.ptr();

        ritmp.get_struct().reset();

        utils::sort_q(ritmp.ptr(),s);

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
            ret = matcl::Matrix(DM(mat.get_type(),0,c0-ne), false);
            return;
        };

        if (ne == c0)
        {
            ret = matcl::Matrix(DM(mat.get_type(),r0,0), false);
            return;
        };

        Integer c1 = c0 - ne;

        bool inplace    = (rvalue == true && mat.is_unique() == true);

        DM out(mat.get_type());

        if (inplace == true)
            out.assign_to_fresh(mat);
        else
            out.assign_to_fresh(DM(mat.get_type(),r0,c1));        

        V* ptr_out              = out.ptr();
        Integer out_ld          = out.ld();
        Integer mat_ld          = mat.ld();

        for (Integer j = 1, k = 0; j <= c0; ++j)
        {
            if (k < s && j == ptr_ri[k])
            {
                while (k < s && j == ptr_ri[k])
                    ++k;
            }
            else
            {
                for (Integer i = 0; i < r0; ++i)	
                    mrd::reset_helper(ptr_out[i], ptr_mat[i]);

                ptr_out += out_ld;
            };

            ptr_mat += mat_ld;            
        };

        if (inplace == true)
        {
            out.set_struct(struct_flag());
            ret = matcl::Matrix(out.resize(r0,c1),true);
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
            ret = matcl::Matrix(DM(mat.get_type(),0,c0-ne), false);
            return;
        };

        if (ne == c0)
        {
            DM out(mat.get_type(),r0,0);
            ret = matcl::Matrix(out,false);
            return;
        };

        Integer c1 = c0 - ne;		

        bool inplace    = (rvalue == true && mat.is_unique() == true);

        DM out(mat.get_type());

        if (inplace == true)
            out.assign_to_fresh(mat);
        else
            out.assign_to_fresh(DM(mat.get_type(),r0,c1));

        V* ptr_out              = out.ptr();
        Integer mat_ld          = mat.ld();
        Integer out_ld          = out.ld();

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

        if (inplace == false)
        {
            for (Integer j = 1; j < cf; ++j)
            {
                for (Integer i = 0; i < r0; ++i)
                    mrd::reset_helper(ptr_out[i],ptr_mat[i]);

                ptr_mat += mat_ld;
                ptr_out += out_ld;
            };
        }
        else
        {
            ptr_mat += mat_ld * (cf - 1);
            ptr_out += out_ld * (cf - 1);
        };
        
        if (cs == 1)
        {
            ptr_mat += imult(ne,mat_ld);
        }
        else
        {
            for (Integer j = cf; j <= cl; ++j)
            {
                if (((j-cf) % cs) != 0)
                {
                    for (Integer i = 0; i < r0; ++i)
                        mrd::reset_helper(ptr_out[i],ptr_mat[i]);

                    ptr_out += out_ld;
                };
                ptr_mat += mat_ld;
            };
        };
        for (Integer j = cl + 1; j <= c0; ++j)
        {
            for (Integer i = 0; i < r0; ++i)
                mrd::reset_helper(ptr_out[i],ptr_mat[i]);

            ptr_mat += mat_ld;
            ptr_out += out_ld;
        };

        if (inplace == true)
        {
            out.set_struct(struct_flag());
            ret = matcl::Matrix(out.resize(r0,c1),true);
            return;
        };

        ret = matcl::Matrix(out,true);
        return;
    };
};

template<class DM>
void del_rowscols_dense_functor<DM>::eval(matcl::Matrix& ret, const DM& mat,const md::colon_info& ci, 
                                      bool rvalue)
{
    Integer sr = ci.rows();
    Integer sc = ci.cols();

    if (sr == 0 && sc == 0)
    {
        ret = matcl::Matrix(mat, false);
        return;
    }
    else if (sr == 0)
    {
        return del_cols_dense_functor<DM>::eval(ret,mat,ci,rvalue);
    }
    else if (sc == 0)
    {
        return del_rows_dense_functor<DM>::eval(ret,mat,ci,rvalue);
    };

    if (ci.r_flag == 0 && ci.c_flag == 0)
        return eval_00(ret,mat,ci,rvalue);
    else if (ci.r_flag == 0 && ci.c_flag == 1)
        return eval_01(ret,mat,ci,rvalue);
    else if (ci.r_flag == 1 && ci.c_flag == 0)
        return eval_10(ret,mat,ci,rvalue);
    else
        return eval_11(ret,mat,ci,rvalue);
};

template<class DM>
void del_rowscols_dense_functor<DM>::eval_00(matcl::Matrix& ret, const DM& mat,const md::colon_info& ci, 
                                      bool rvalue)
{
    Integer sc  = ci.cols();
    Integer sr  = ci.rows();
    Integer r0  = mat.rows();
    Integer c0  = mat.cols();

    using V     = typename DM::value_type;

    const V* ptr_mat            = mat.ptr();

    //column colon
    const mr::integer_dense& ci_ci  = ci.get_cim_2();
    raw::integer_dense citmp        = ci_ci.copy();
    const Integer* ptr_ci           = citmp.ptr();

    citmp.get_struct().reset();

    utils::sort_q(citmp.ptr(),sc);

    Integer cec  = ptr_ci[0];
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

    if (!r0)
    {
        ret = matcl::Matrix(DM(mat.get_type(),0,c0-nec), false);
        return;
    };

    //row colon
    const mr::integer_dense& ci_ri  = ci.get_rim_2();
    raw::integer_dense ritmp        = ci_ri.copy();
    const Integer* ptr_ri           = ritmp.ptr();

    ritmp.get_struct().reset();

    utils::sort_q(ritmp.ptr(),sr);

    Integer cer = ptr_ri[0];
    error::check_row(cer, r0, c0);
    Integer ner = 1;

    for (Integer i = 1; i < sr; ++i)
    {
        Integer ce2 = ptr_ri[i];
        error::check_row(ce2, r0, c0);
    
        if (ce2 != cer)
            ++ner;

        cer = ce2;
    }

    if (!c0)
    {
        ret = matcl::Matrix(DM(mat.get_type(),r0-ner,0), false);
        return;
    };

    //fast exit
    Integer c1 = c0 - nec;
    Integer r1 = r0 - ner;

    if (nec == c0 || ner == r0)
    {
        ret = matcl::Matrix(DM(mat.get_type(),r1,c1), false);
        return;
    };

    //general case
    bool inplace    = (rvalue == true && mat.is_unique() == true);

    DM out(mat.get_type());

    if (inplace == true)
        out.assign_to_fresh(mat);
    else
        out.assign_to_fresh(DM(mat.get_type(),r1,c1));

    V* ptr_out              = out.ptr();
    Integer out_ld          = out.ld();
    Integer mat_ld          = mat.ld();

    for (Integer j = 1, k = 0; j <= c0; ++j)
    {
        if (k < sc && j == ptr_ci[k])
        {
            while (k < sc && j == ptr_ci[k])
                ++k;
        }
        else
        {
            for (Integer i = 0, ii = 0, k2 = 0; i < r0; ++i)
            {
                if (k2 < sr && i+1 == ptr_ri[k2])
                {
                    while(k2 < sr && i+1 == ptr_ri[k2])
                        ++k2;
                }
                else
                {
                    mrd::reset_helper(ptr_out[ii],ptr_mat[i]);
                    ++ii;
                };
            };

            ptr_out += out_ld;
        };

        ptr_mat += mat_ld;            
    };

    if (inplace == true)
    {
        out.set_struct(struct_flag());
        ret = matcl::Matrix(out.resize(r1,c1),true);
        return;
    };

    ret = matcl::Matrix(out,true);
    return;
};

template<class DM>
void del_rowscols_dense_functor<DM>::eval_10(matcl::Matrix& ret, const DM& mat,const md::colon_info& ci, 
                                      bool rvalue)
{
    Integer sr  = ci.rows();
    Integer sc  = ci.cols();
    Integer r0  = mat.rows();
    Integer c0  = mat.cols();

    using V     = typename DM::value_type;

    const V* ptr_mat            = mat.ptr();

    //column colon
    const mr::integer_dense& ci_ci  = ci.get_cim_2();
    raw::integer_dense citmp        = ci_ci.copy();
    const Integer* ptr_ci           = citmp.ptr();

    citmp.get_struct().reset();

    utils::sort_q(citmp.ptr(),sc);

    Integer cec  = ptr_ci[0];
    error::check_col(cec, r0, c0);
    Integer nec  = 1;

    for (Integer i = 1; i < sc; ++i)
    {
        Integer ce2 = ptr_ci[i];
        error::check_col(ce2, r0, c0);
            
        if (ce2 != cec)
            ++nec;

        cec = ce2;
    }

    Integer ner = sr;

    Integer r1 = r0 - ner;
    Integer c1 = c0 - nec;

    if (r1 == 0 || c1 == 0)
    {
        ret = matcl::Matrix(DM(mat.get_type(),r1,c1), false);
        return;
    };

    bool inplace    = (rvalue == true && mat.is_unique() == true);

    DM out(mat.get_type());

    if (inplace == true)
        out.assign_to_fresh(mat);
    else
        out.assign_to_fresh(DM(mat.get_type(),r1,c1));        

    V* ptr_out              = out.ptr();
    Integer out_ld          = out.ld();
    Integer mat_ld          = mat.ld();

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

    for (Integer j = 1, k = 0; j <= c0; ++j)
    {
        if (k < sc && j == ptr_ci[k])
        {
            while (k < sc && j == ptr_ci[k])
                ++k;
        }
        else
        {
            Integer ii = 0, iii = 0;

            for (Integer i = 1; i < rf; ++i, ++ii)
            {
                mrd::reset_helper(ptr_out[iii++], ptr_mat[ii]);
            };

            if (rs == 1)
            {
                ii += ner;
            }
            else
            {
                for (Integer i = rf; i <= rl; ++i, ++ii)
                {
                    if (((i-rf) % rs) != 0)
                        mrd::reset_helper(ptr_out[iii++], ptr_mat[ii]);
                };
            }

            for (Integer i = rl + 1; i <= r0; ++i, ++ii)
                mrd::reset_helper(ptr_out[iii++],ptr_mat[ii]);

            ptr_out += out_ld;
        };

        ptr_mat += mat_ld;            
    };

    if (inplace == true)
    {
        out.set_struct(struct_flag());
        ret = matcl::Matrix(out.resize(r1,c1),true);
        return;
    };

    ret = matcl::Matrix(out,true);
    return;
};

template<class DM>
void del_rowscols_dense_functor<DM>::eval_01(matcl::Matrix& ret, const DM& mat,const md::colon_info& ci, 
                                      bool rvalue)
{
    using V     = typename DM::value_type;

    Integer sc  = ci.cols();
    Integer sr  = ci.rows();
    Integer r0  = mat.rows();
    Integer c0  = mat.cols();

    Integer nec = sc;
    
    //row colon
    const mr::integer_dense& ci_ri  = ci.get_rim_2();
    raw::integer_dense ritmp        = ci_ri.copy();
    const Integer* ptr_ri           = ritmp.ptr();

    ritmp.get_struct().reset();

    utils::sort_q(ritmp.ptr(),sr);

    Integer cer = ptr_ri[0];
    error::check_row(cer, r0, c0);
    Integer ner = 1;

    for (Integer i = 1; i < sr; ++i)
    {
        Integer ce2 = ptr_ri[i];
        error::check_row(ce2, r0, c0);
    
        if (ce2 != cer)
            ++ner;

        cer = ce2;
    }

    Integer r1 = r0 - ner;
    Integer c1 = c0 - nec;

    if (r1 == 0 || c1 == 0)
    {
        ret = matcl::Matrix(DM(mat.get_type(),r1,c1), false);
        return;
    };

    bool inplace    = (rvalue == true && mat.is_unique() == true);

    DM out(mat.get_type());

    if (inplace == true)
        out.assign_to_fresh(mat);
    else
        out.assign_to_fresh(DM(mat.get_type(),r1,c1));

    V* ptr_out              = out.ptr();
    Integer mat_ld          = mat.ld();
    Integer out_ld          = out.ld();
    const V* ptr_mat        = mat.ptr();

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

    for (Integer j = 1; j < cf; ++j)
    {
        for (Integer i = 0, ii = 0, k = 0; i < r0; ++i)
        {
            if (k < sr && i+1 == ptr_ri[k])
            {
                while(k < sr && i+1 == ptr_ri[k])
                    ++k;
            }
            else
            {
                mrd::reset_helper(ptr_out[ii],ptr_mat[i]);
                ++ii;
            };
        };

        ptr_mat += mat_ld;
        ptr_out += out_ld;
    };

    if (cs == 1)
    {
        ptr_mat += imult(nec,mat_ld);
    }
    else
    {
        for (Integer j = cf; j <= cl; ++j)
        {
            if (((j-cf) % cs) != 0)
            {
                for (Integer i = 0, ii = 0, k = 0; i < r0; ++i)
                {
                    if (k < sr && i+1 == ptr_ri[k])
                    {
                        while(k < sr && i+1 == ptr_ri[k])
                            ++k;
                    }
                    else
                    {
                        mrd::reset_helper(ptr_out[ii],ptr_mat[i]);
                        ++ii;
                    };
                };

                ptr_out += out_ld;
            };
            ptr_mat += mat_ld;
        };
    };
    for (Integer j = cl + 1; j <= c0; ++j)
    {
        for (Integer i = 0, ii = 0, k = 0; i < r0; ++i)
        {
            if (k < sr && i+1 == ptr_ri[k])
            {
                while(k < sr && i+1 == ptr_ri[k])
                    ++k;
            }
            else
            {
                mrd::reset_helper(ptr_out[ii],ptr_mat[i]);
                ++ii;
            };
        };

        ptr_mat += mat_ld;
        ptr_out += out_ld;
    };

    if (inplace == true)
    {
        out.set_struct(struct_flag());
        ret = matcl::Matrix(out.resize(r1,c1),true);
        return;
    };

    ret = matcl::Matrix(out,true);
    return;
};

template<class DM>
void del_rowscols_dense_functor<DM>::eval_11(matcl::Matrix& ret, const DM& mat,const md::colon_info& ci, 
                                      bool rvalue)
{
    using V     = typename DM::value_type;

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
        ret = matcl::Matrix(DM(mat.get_type(),r1,c1), false);
        return;
    };

    bool inplace    = (rvalue == true && mat.is_unique() == true);

    DM out(mat.get_type());

    if (inplace == true)
        out.assign_to_fresh(mat);
    else
        out.assign_to_fresh(DM(mat.get_type(),r1,c1));

    V* ptr_out              = out.ptr();
    Integer mat_ld          = mat.ld();
    Integer out_ld          = out.ld();
    const V* ptr_mat        = mat.ptr();

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

    for (Integer j = 1; j < cf; ++j)
    {
        Integer ii = 0, iii = 0;
            
        for (Integer i = 1; i < rf; ++i, ++ii)
            mrd::reset_helper(ptr_out[iii++], ptr_mat[ii]);
        
        if (rs == 1)
        {
            ii += ner;
        }
        else
        {
            for (Integer i = rf; i <= rl; ++i, ++ii)
            {
                if (((i-rf) % rs) != 0)
                    mrd::reset_helper(ptr_out[iii++], ptr_mat[ii]);
            };
        }

        for (Integer i = rl + 1; i <= r0; ++i, ++ii)
            mrd::reset_helper(ptr_out[iii++],ptr_mat[ii]);

        ptr_mat += mat_ld;
        ptr_out += out_ld;
    };
        
    if (cs == 1)
    {
        ptr_mat += imult(nec,mat_ld);
    }
    else
    {
        for (Integer j = cf; j <= cl; ++j)
        {
            if (((j-cf) % cs) != 0)
            {
                Integer ii = 0, iii = 0;
            
                for (Integer i = 1; i < rf; ++i, ++ii)
                    mrd::reset_helper(ptr_out[iii++], ptr_mat[ii]);
        
                if (rs == 1)
                {
                    ii += ner;
                }
                else
                {
                    for (Integer i = rf; i <= rl; ++i, ++ii)
                    {
                        if (((i-rf) % rs) != 0)
                            mrd::reset_helper(ptr_out[iii++], ptr_mat[ii]);
                    };
                }

                for (Integer i = rl + 1; i <= r0; ++i, ++ii)
                    mrd::reset_helper(ptr_out[iii++],ptr_mat[ii]);

                ptr_out += out_ld;
            };

            ptr_mat += mat_ld;
        };
    };

    for (Integer j = cl + 1; j <= c0; ++j)
    {
        Integer ii = 0, iii = 0;
            
        for (Integer i = 1; i < rf; ++i, ++ii)
            mrd::reset_helper(ptr_out[iii++], ptr_mat[ii]);
        
        if (rs == 1)
        {
            ii += ner;
        }
        else
        {
            for (Integer i = rf; i <= rl; ++i, ++ii)
            {
                if (((i-rf) % rs) != 0)
                    mrd::reset_helper(ptr_out[iii++], ptr_mat[ii]);
            };
        }

        for (Integer i = rl + 1; i <= r0; ++i, ++ii)
            mrd::reset_helper(ptr_out[iii++],ptr_mat[ii]);

        ptr_mat += mat_ld;
        ptr_out += out_ld;
    };

    if (inplace == true)
    {
        out.set_struct(struct_flag());
        ret = matcl::Matrix(out.resize(r1,c1),true);
        return;
    };

    ret = matcl::Matrix(out,true);
    return;
};

template<class DM>
void dense_change_diag_functor<DM>::eval(DM& mat, Integer d, const DM& val)
{
    Integer r = mat.rows();
    Integer c = mat.cols();

    error::check_diag(d,r,c);

    Integer st, s;
    if (d >= 0)
    {
        st = imult(d,r);
        s = (r + d >= c) ? c - d : r;
    }
    else
    {
        st = - d;
        s = (r + d >= c) ? c : r + d;
    }

    error::check_assign_1(s,val.size(),1);

    if (s == 0)
        return;

    value_type * ptr    = mat.ptr() + st;
    const DM& val_tmp   = val.make_explicit();
    const value_type* val_ptr   = val_tmp.ptr();
    Integer mat_ld              = mat.ld();

    for (Integer i = 0; i < s; ++i)
    {
        mrd::assign_helper(*(ptr), *(val_ptr++));
        ptr += mat_ld + 1;
    }

    mat.set_struct(md::predefined_struct::get_set_diag(mat.get_struct(),d,value_struct_class::vc_general,
                                    is_real_matrix(mat), r == c));

    return;
};

template<class DM>
void dense_change_diag_functor<DM>::eval(DM& mat, Integer d, const value_type& val)
{
    Integer r = mat.rows();
    Integer c = mat.cols();

    error::check_diag(d,r,c);

    bool tz = md::has_trivial_assignment<DM,value_type>::eval(mat,val);

    Integer st, s;
    if (d >= 0)
    {
        st = imult(d,r);
        s = (r + d >= c) ? c - d : r;
    }
    else
    {
        st = - d;
        s = (r + d >= c) ? c : r + d;
    }    

    if (s == 0)
        return;

    value_type * ptr    = mat.ptr() + st;
    Integer mat_ld      = mat.ld();

    for (Integer i = 0; i < s; ++i)
    {
        mrd::assign_helper(*(ptr), val);
        ptr += mat_ld + 1;
    }

    value_struct_class vt = md::predefined_struct::get_value_type(val,tz);
    mat.set_struct(md::predefined_struct::get_set_diag(mat.get_struct(),d,vt, 
                            is_real_matrix(mat), r == c));
    return;
};

};};};

template struct matcl::algorithm::details::dense_change_entries_functor<matcl::raw::integer_dense>;
template struct matcl::algorithm::details::dense_change_entries_functor<matcl::raw::real_dense>;
template struct matcl::algorithm::details::dense_change_entries_functor<matcl::raw::float_dense>;
template struct matcl::algorithm::details::dense_change_entries_functor<matcl::raw::complex_dense>;
template struct matcl::algorithm::details::dense_change_entries_functor<matcl::raw::float_complex_dense>;
template struct matcl::algorithm::details::dense_change_entries_functor<matcl::raw::object_dense>;

template struct matcl::algorithm::details::dense_change_entries_functor_2<matcl::raw::integer_dense>;
template struct matcl::algorithm::details::dense_change_entries_functor_2<matcl::raw::real_dense>;
template struct matcl::algorithm::details::dense_change_entries_functor_2<matcl::raw::float_dense>;
template struct matcl::algorithm::details::dense_change_entries_functor_2<matcl::raw::complex_dense>;
template struct matcl::algorithm::details::dense_change_entries_functor_2<matcl::raw::float_complex_dense>;
template struct matcl::algorithm::details::dense_change_entries_functor_2<matcl::raw::object_dense>;

template struct matcl::algorithm::details::del_cols_dense_functor<matcl::raw::integer_dense>;
template struct matcl::algorithm::details::del_cols_dense_functor<matcl::raw::real_dense>;
template struct matcl::algorithm::details::del_cols_dense_functor<matcl::raw::float_dense>;
template struct matcl::algorithm::details::del_cols_dense_functor<matcl::raw::complex_dense>;
template struct matcl::algorithm::details::del_cols_dense_functor<matcl::raw::float_complex_dense>;
template struct matcl::algorithm::details::del_cols_dense_functor<matcl::raw::object_dense>;

template struct matcl::algorithm::details::del_rows_dense_functor<matcl::raw::integer_dense>;
template struct matcl::algorithm::details::del_rows_dense_functor<matcl::raw::real_dense>;
template struct matcl::algorithm::details::del_rows_dense_functor<matcl::raw::float_dense>;
template struct matcl::algorithm::details::del_rows_dense_functor<matcl::raw::complex_dense>;
template struct matcl::algorithm::details::del_rows_dense_functor<matcl::raw::float_complex_dense>;
template struct matcl::algorithm::details::del_rows_dense_functor<matcl::raw::object_dense>;

template struct matcl::algorithm::details::del_rowscols_dense_functor<matcl::raw::integer_dense>;
template struct matcl::algorithm::details::del_rowscols_dense_functor<matcl::raw::real_dense>;
template struct matcl::algorithm::details::del_rowscols_dense_functor<matcl::raw::float_dense>;
template struct matcl::algorithm::details::del_rowscols_dense_functor<matcl::raw::complex_dense>;
template struct matcl::algorithm::details::del_rowscols_dense_functor<matcl::raw::float_complex_dense>;
template struct matcl::algorithm::details::del_rowscols_dense_functor<matcl::raw::object_dense>;

template struct matcl::algorithm::details::dense_change_diag_functor<matcl::raw::integer_dense>;
template struct matcl::algorithm::details::dense_change_diag_functor<matcl::raw::real_dense>;
template struct matcl::algorithm::details::dense_change_diag_functor<matcl::raw::float_dense>;
template struct matcl::algorithm::details::dense_change_diag_functor<matcl::raw::complex_dense>;
template struct matcl::algorithm::details::dense_change_diag_functor<matcl::raw::float_complex_dense>;
template struct matcl::algorithm::details::dense_change_diag_functor<matcl::raw::object_dense>;
