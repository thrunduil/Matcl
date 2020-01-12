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
#include "matcl-internals/container/mat_b.h"
#include "matcl-matrep/base/sort_iterator.h"
#include "matcl-matrep/func/raw/raw_manip.h"
#include "matcl-scalar/details/scalfunc_helpers.h"
#include "matcl-internals/func/converter.h"
#include "matcl-matrep/details/matrix.inl"

namespace matcl { namespace algorithm { namespace details
{

namespace mr    = matcl::raw;
namespace mrd   = matcl::raw::details;
namespace md    = matcl::details;

template<class T>
class elem_getter1
{
    private:
        using SM        = typename matcl::raw::Matrix<T,struct_sparse>;

        const SM&		A;
        const Integer * d_c;
        const Integer * d_r;
        const T *		d_x;
        Integer			step;
        Integer			pos_f;
        Integer			pos_l;
        Integer			row;
        Integer			elem;

    public:
        elem_getter1(const SM& A, bool decr)
            :A(A)
        {
            const raw::details::sparse_ccs<T>& Ad = A.rep();
            d_c	= Ad.ptr_c();
            d_r	= Ad.ptr_r();
            d_x	= Ad.ptr_x();
            
            if (decr)
                step	= -1;
            else
                step	= 1;
        };

        void set_col(Integer col)
        {		
            if (step < 0)
            {
                pos_f	= d_c[col]-1;
                pos_l	= d_c[col-1]-1;
                elem	= A.rows()-1;
            }
            else
            {
                pos_f	= d_c[col-1];
                pos_l	= d_c[col];
                elem	= 0;
            };

            if (pos_f == pos_l)
                row	    = -1;
            else
                row		= d_r[pos_f];					
        };

        T get(bool& has_nz)
        {			
            has_nz = false;

            if (pos_f == pos_l)
                return md::default_value<T>(A.get_type());

            T value = md::default_value<T>(A.get_type());

            if (elem == row)
            {
                mrd::reset_helper(value,d_x[pos_f]);
                pos_f	+= step;
                row		= d_r[pos_f];
                has_nz  = true;
            };			
            
            elem	+= step;
            return value;
        };

    private:
        elem_getter1(const elem_getter1&) = delete;
        elem_getter1& operator=(const elem_getter1&) = delete;
};

template<class T>
class elem_getter2
{
    private:
        using SM        = typename raw::Matrix<T,struct_sparse> ;

        const SM&		A;
        const Integer * d_c;
        const Integer * d_r;
        const T *		d_x;

        Integer			m_rf;
        Integer			m_rf0;
        Integer			m_rl;
        Integer			m_rs;
        Integer			step;

        Integer			pos_f;
        Integer			pos_l;

    public:
        elem_getter2(const SM& A, Integer rf, Integer rs, Integer rl)
            :A(A)
        {
            const raw::details::sparse_ccs<T>& Ad = A.rep();
            d_c	= Ad.ptr_c();
            d_r	= Ad.ptr_r();
            d_x	= Ad.ptr_x();
            
            m_rf0	= rf;
            m_rs	= rs;

            if (rs < 0)
            {
                m_rf	= rl;
                m_rl	= rf;
                step	= -1;
            }
            else
            {	
                step	= 1;
                m_rf	= rf;
                m_rl	= rl;
            };
        };

        void set_col(Integer col, Integer& k)
        {		
            if (m_rs < 0)
            {
                pos_f	= d_c[col]-1;
                pos_l	= d_c[col-1]-1;			
            }
            else
            {
                pos_f	= d_c[col-1];
                pos_l	= d_c[col];
            };

            if (pos_f == pos_l)
            {
                k			= m_rl+1;
            }
            else
            {
                Integer r	= d_r[pos_f];
                k			= m_rf0 + imult(r,m_rs);
            };		
        };

        const T& get(Integer& k)
        {			
            const T& value  = d_x[pos_f];			
            pos_f		    += step;

            if (pos_f == pos_l)
            {				
                k           = m_rl+1;
            }
            else
            {
                Integer r	= d_r[pos_f];
                k			= m_rf0 + imult(r,m_rs);
            };

            return value;
        };

    private:
        elem_getter2(const elem_getter2&) = delete;
        elem_getter2& operator=(const elem_getter2&) = delete;
};

class column_iterator_2
{
    private:
        Integer					m_flag;
        Integer					m_size;
        Integer					m_dpos;
        Integer					m_dpos_col;
        Integer					m_column;
        raw::integer_dense	    m_indices;
        bool					m_own_sort;

        raw::integer_dense	    m_vec;		
        Integer					m_pos;	

    public:
        column_iterator_2(const md::colon_info& c_inf)
            :m_indices(ti::ti_int()), m_vec(ti::ti_int())
        {
            m_flag = c_inf.c_flag;

            if (m_flag == 0)
            {
                m_vec.assign_to_fresh(c_inf.get_cim_2());
                sort_type m_sort	= is_sorted(m_vec);
                m_dpos				= 1;
                m_own_sort			= false;

                if (m_sort == sorted_increasing)
                {
                }
                else if (m_sort == sorted_decreasing)
                {
                    m_dpos			= -1;
                } 
                else
                {
                    m_own_sort		= true;
                    raw::integer_dense ci2 = m_vec.copy();
                    ci2.get_struct().reset();
                    m_indices.assign_to_fresh(raw::integer_dense(ti::ti_int(),1,ci2.size()));
                    Integer* ptr_ind = m_indices.ptr();

                    Integer size = ci2.size();
                    for (Integer i = 0; i < size; ++i)
                        ptr_ind[i] = i+1;

                    using iterator  = md::iterator_helper_2<Integer,Integer>;
                    iterator it_begin(ci2.ptr(),m_indices.ptr(),1);
                    iterator it_end(ci2.ptr()+ci2.size(),m_indices.ptr()+ci2.size(),1);

                    static const bool asceding = true;
                    std::stable_sort(it_begin,it_end,md::value_ind_compare<Integer,Integer,asceding>());

                    m_vec.assign_to_fresh(std::move(ci2));
                };

                if (m_dpos > 0)
                {
                    m_pos			= 1;
                    m_column		= 1;
                    m_dpos_col		= 1;
                }
                else
                {
                    m_pos			= m_vec.size();
                    m_column		= c_inf.cols();
                    m_dpos_col		= -1;
                };
            }
            else
            {				
                m_dpos				= c_inf.c_step;
                if (m_dpos>0)
                {
                    m_column		= 1;
                    m_pos			= c_inf.c_start;
                    m_dpos_col		= 1;
                }
                else
                {
                    m_pos			= c_inf.c_end;
                    m_dpos			= -m_dpos;
                    m_dpos_col		= -1;
                    m_column		= c_inf.c_size;
                };
            };

            m_size					= c_inf.cols();
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

        bool is_increasing() const
        {
            return m_dpos_col > 0;
        };

        void next()
        {
            m_pos		+= m_dpos;
            m_column	+= m_dpos_col;
            --m_size;
        };

        Integer column() const
        {
            if (m_flag == 0)
            {
                if (m_own_sort)
                    return m_indices.ptr()[m_pos-1];
                else
                    return m_column;			
            }
            else
            {
                return m_column;
            };
        };
};

template<class SM,class M2>
struct change_submatrix_0_impl
{
    using value_type    = typename SM::value_type;
    using value_type_B  = typename M2::value_type;

    static void eval_dense(const SM& A, const md::colon_info& ci, const M2& B,
                           Matrix& ret)
    {
        if (ci.is_double_mat_colon() == true)
            return eval_dense_dc(A, ci, B, ret);

        using value_type    = typename SM::value_type;
        using value_type_B  = typename M2::value_type;

        Integer r = A.rows(), c = A.cols();
        const raw::details::sparse_ccs<value_type>& Ad = A.rep();
        M2 vB = B;

        raw::integer_dense ri   = ci.get_rim_1();
        sort_type s_type        = is_sorted(ri);
        bool incr               = true;

        if (s_type == sorted_increasing)
        {
        }
        else if (s_type == sorted_decreasing)
        {
            incr = false;
        } 
        else
        {
            raw::integer_dense ri2 = ri.copy();
            ri2.get_struct().reset();
            vB.assign_to_fresh(B.copy());

            using iterator = md::iterator_helper_2<Integer,value_type_B>;
            iterator it_begin(ri2.ptr(),vB.ptr(),1);
            iterator it_end(ri2.ptr()+ri2.size(),vB.ptr()+ri2.size(),1);

            static const bool asceding = true;
            std::stable_sort(it_begin,it_end,md::value_ind_compare<Integer,value_type_B,asceding>());

            ri.assign_to_fresh(std::move(ri2));
        };

        Integer step_ri;
        Integer pos, pos_row, pos_col;

        if (incr)
        {
            pos = 0;
            step_ri = 1;
        }
        else
        {
            pos = ri.size()-1;
            step_ri = -1;
        };

        raw::details::sparse_ccs<value_type> d(A.get_type(),r, c, A.nnz() + ri.size());

        Integer nz				= 0;

        Integer * d_c			= d.ptr_c();
        Integer * d_r			= d.ptr_r();
        value_type * d_x		= d.ptr_x();

        const Integer * Ad_c	= Ad.ptr_c();
        const Integer * Ad_r	= Ad.ptr_r();
        const value_type * Ad_x	= Ad.ptr_x();

        vB.assign_to_fresh(vB.make_explicit());

        const value_type_B* ptr_B = vB.ptr();

        Integer or              = ri.size();
        const Integer* ptr_ri   = ri.ptr();
        Integer old_pos         = ptr_ri[pos];

        md::pos2ind(ptr_ri[pos],r,pos_row,pos_col);

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
                    Integer elem_pos = pos;
                    while (pos < or && pos >= 0 && ptr_ri[pos] == old_pos)
                    {
                        if(incr)
                            elem_pos = pos;

                        pos+=step_ri;
                    };
                    
                    md::pos2ind(ptr_ri[elem_pos],r,pos_row,pos_col);

                    value_type tmp(ptr_B[elem_pos]);

                    d_r[nz] = pos_row;
                    mrd::assign_helper(d_x[nz],tmp);
                    ++nz;

                    if (pos > or-1 || pos < 0)
                    {
                        pos_row = r;
                        pos_col = c;
                        goto exit_row_flag;
                    };

                    old_pos = ptr_ri[pos];
                    md::pos2ind(ptr_ri[pos],r,pos_row,pos_col);

                    if (j != pos_col)
                        goto exit_row_flag;
                };

                if (p == pos_row)
                {
                    Integer elem_pos = pos;

                    while (pos < or && pos >= 0 && ptr_ri[pos] == old_pos)
                    {
                        if(incr)
                            elem_pos = pos;

                        pos+=step_ri;
                    };
                    
                    md::pos2ind(ptr_ri[elem_pos],r,pos_row,pos_col);

                    value_type tmp(ptr_B[elem_pos]);

                    d_r[nz] = pos_row;
                    mrd::assign_change_helper(d_x[nz],Ad_x[i],tmp);
                    ++nz;
                    ++i;

                    if (pos > or-1 || pos < 0)
                    {
                        pos_row = r;
                        pos_col = c;
                        goto exit_row_flag;
                    };

                    old_pos = ptr_ri[pos];
                    md::pos2ind(ptr_ri[pos],r,pos_row,pos_col);

                    if (j != pos_col)
                        goto exit_row_flag;	
                };
            };

            while (j == pos_col)
            {
                Integer elem_pos = pos;
                while (pos < or && pos >= 0 && ptr_ri[pos] == old_pos)
                {
                    if(incr)
                        elem_pos = pos;

                    pos+=step_ri;
                };
                
                md::pos2ind(ptr_ri[elem_pos],r,pos_row,pos_col);

                value_type tmp(ptr_B[elem_pos]);

                d_r[nz] = pos_row;
                mrd::assign_helper(d_x[nz],tmp);
                ++nz;

                if (pos > or-1 || pos < 0)
                {
                    pos_row = r;
                    pos_col = c;
                    goto exit_row_flag;
                };

                old_pos = ptr_ri[pos];
                md::pos2ind(ptr_ri[pos],r,pos_row,pos_col);

                if (j != pos_col)
                    goto exit_row_flag;	
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

    static void eval_sparse_impl(const SM& A, const md::colon_info& ci, 
                    const value_type_B* ptr_B0, Integer* pos_B, Matrix& ret)
    {
        if (ci.is_double_mat_colon() == true)
            return eval_sparse_dc_impl(A, ci, ptr_B0, pos_B, ret);

        Integer r = A.rows(), c = A.cols();
        const raw::details::sparse_ccs<value_type>& Ad = A.rep();

        raw::integer_dense ri   = ci.get_rim_1();
        sort_type s_type        = is_sorted(ri);
        bool incr               = true;

        if (s_type == sorted_increasing)
        {
        }
        else if (s_type == sorted_decreasing)
        {
            incr = false;
        } 
        else
        {
            raw::integer_dense ri2 = ri.copy();
            ri2.get_struct().reset();

            using iterator = md::iterator_helper_2<Integer, Integer>;
            iterator it_begin(ri2.ptr(), pos_B, 1);
            iterator it_end(ri2.ptr()+ri2.size(), pos_B + ri2.size(),1);

            static const bool asceding = true;
            std::stable_sort(it_begin,it_end,md::value_ind_compare<Integer,Integer,asceding>());

            ri.assign_to_fresh(std::move(ri2));
        };

        Integer step_ri;
        Integer pos, pos_row, pos_col;
        if (incr)
        {
            pos = 0;
            step_ri = 1;
        }
        else
        {
            pos = ri.size()-1;
            step_ri = -1;
        };

        raw::details::sparse_ccs<value_type> d(A.get_type(),r, c, A.nnz() + ri.size());

        Integer nz				= 0;

        Integer * d_c			= d.ptr_c();
        Integer * d_r			= d.ptr_r();
        value_type * d_x		= d.ptr_x();

        const Integer * Ad_c	= Ad.ptr_c();
        const Integer * Ad_r	= Ad.ptr_r();
        const value_type * Ad_x	= Ad.ptr_x();

        Integer or              = ri.size();
        const Integer* ptr_ri   = ri.ptr();
        Integer old_pos         = ptr_ri[pos];

        md::pos2ind(ptr_ri[pos],r,pos_row,pos_col);

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
                    Integer elem_pos = pos;
                    while (pos < or && pos >= 0 && ptr_ri[pos] == old_pos)
                    {
                        if(incr)
                            elem_pos = pos;

                        pos+=step_ri;
                    };
                    
                    md::pos2ind(ptr_ri[elem_pos],r,pos_row,pos_col);

                    Integer pos_B_elem  = pos_B[elem_pos];
                    if (pos_B_elem >= 0)
                    {
                        value_type tmp(ptr_B0[pos_B_elem]);

                        d_r[nz] = pos_row;
                        mrd::assign_helper(d_x[nz],tmp);
                        ++nz;
                    };

                    if (pos > or-1 || pos < 0)
                    {
                        pos_row = r;
                        pos_col = c;
                        goto exit_row_flag;
                    };

                    old_pos = ptr_ri[pos];
                    md::pos2ind(ptr_ri[pos],r,pos_row,pos_col);

                    if (j != pos_col)
                        goto exit_row_flag;
                };

                if (p == pos_row)
                {
                    Integer elem_pos = pos;
                    while (pos < or && pos >= 0 && ptr_ri[pos] == old_pos)
                    {
                        if(incr)
                            elem_pos = pos;

                        pos+=step_ri;
                    };
                    
                    md::pos2ind(ptr_ri[elem_pos],r,pos_row,pos_col);

                    Integer pos_B_elem  = pos_B[elem_pos];
                    if (pos_B_elem >= 0)
                    {
                        value_type tmp(ptr_B0[pos_B_elem]);

                        d_r[nz] = pos_row;
                        mrd::assign_change_helper(d_x[nz],Ad_x[i],tmp);
                        ++nz;
                    };

                    ++i;

                    if (pos > or-1 || pos < 0)
                    {
                        pos_row = r;
                        pos_col = c;
                        goto exit_row_flag;
                    };

                    old_pos = ptr_ri[pos];
                    md::pos2ind(ptr_ri[pos],r,pos_row,pos_col);

                    if (j != pos_col)
                        goto exit_row_flag;	
                };
            };

            while (j == pos_col)
            {
                Integer elem_pos = pos;
                while (pos < or && pos >= 0 && ptr_ri[pos] == old_pos)
                {
                    if(incr)
                        elem_pos = pos;

                    pos+=step_ri;
                };
                
                md::pos2ind(ptr_ri[elem_pos],r,pos_row,pos_col);

                Integer pos_B_elem  = pos_B[elem_pos];
                if (pos_B_elem >= 0)
                {
                    value_type tmp(ptr_B0[pos_B_elem]);

                    d_r[nz] = pos_row;
                    mrd::assign_helper(d_x[nz],tmp);
                    ++nz;
                };

                if (pos > or-1 || pos < 0)
                {
                    pos_row = r;
                    pos_col = c;
                    goto exit_row_flag;
                };

                old_pos = ptr_ri[pos];
                md::pos2ind(ptr_ri[pos],r,pos_row,pos_col);

                if (j != pos_col)
                    goto exit_row_flag;	
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

    static void eval_dense_dc(const SM& A, const md::colon_info& ci, const M2& B,
                              Matrix& ret)
    {
        using value_type    = typename SM::value_type;
        using value_type_B  = typename M2::value_type;

        Integer r = A.rows();
        Integer c = A.cols();
        const raw::details::sparse_ccs<value_type>& Ad = A.rep();
        M2 vB = B.copy();

        raw::integer_dense cr = ci.get_rim_r().copy();
        raw::integer_dense cc = ci.get_rim_c().copy();

        sort_rows_cols_stable(cr.ptr(), cc.ptr(), vB.ptr(), cr.size());              

        raw::details::sparse_ccs<value_type> d(A.get_type(),r, c, A.nnz() + cr.size());

        Integer nz				= 0;

        Integer * d_c			= d.ptr_c();
        Integer * d_r			= d.ptr_r();
        value_type * d_x		= d.ptr_x();

        const Integer * Ad_c	= Ad.ptr_c();
        const Integer * Ad_r	= Ad.ptr_r();
        const value_type * Ad_x	= Ad.ptr_x();

        const value_type_B* ptr_B = vB.ptr();
        
        Integer or              = cr.size();
        const Integer* ptr_r    = cr.ptr();
        const Integer* ptr_c    = cc.ptr();

        Integer pos             = 0;  
        Integer pos_r           = ptr_r[0] - 1;
        Integer pos_c           = ptr_c[0] - 1;

        for (Integer j = 0; j < c; ++j)
        {
            d_c[j] = nz;

            if (j < pos_c)
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
                if (p < pos_r)
                {
                    d_r[nz] = p;
                    mrd::assign_helper(d_x[nz],Ad_x[i]);
                    ++nz;
                    ++i;
                    continue;
                };	

                while (p > pos_r)
                {
                    Integer elem_pos = pos;

                    while (pos < or && ptr_r[pos] - 1 == pos_r && ptr_c[pos] - 1 == pos_c)
                    {
                        elem_pos = pos;
                        pos += 1;
                    };

                    value_type tmp(ptr_B[elem_pos]);

                    d_r[nz] = pos_r;
                    mrd::assign_helper(d_x[nz],tmp);
                    ++nz;

                    if (pos > or-1)
                    {
                        pos_r   = r;
                        pos_c   = c;
                        goto exit_row_flag;
                    };

                    pos_r       = ptr_r[pos] - 1;
                    pos_c       = ptr_c[pos] - 1;

                    if (j != pos_c)
                        goto exit_row_flag;
                };

                if (p == pos_r)
                {
                    Integer elem_pos = pos;
                    while (pos < or && ptr_r[pos] - 1 == pos_r && ptr_c[pos] - 1 == pos_c)
                    {
                        elem_pos = pos;
                        pos     += 1;
                    };

                    value_type tmp(ptr_B[elem_pos]);

                    d_r[nz] = pos_r;
                    mrd::assign_change_helper(d_x[nz],Ad_x[i],tmp);
                    ++nz;
                    ++i;

                    if (pos > or-1)
                    {
                        pos_r   = r;
                        pos_c   = c;
                        goto exit_row_flag;
                    };

                    pos_r       = ptr_r[pos] - 1;
                    pos_c       = ptr_c[pos] - 1;

                    if (j != pos_c)
                        goto exit_row_flag;	
                };
            };

            while (j == pos_c)
            {
                Integer elem_pos = pos;
                while (pos < or && ptr_r[pos] - 1 == pos_r && ptr_c[pos] - 1 == pos_c)
                {
                    elem_pos = pos;
                    pos     += 1;
                };

                value_type tmp(ptr_B[elem_pos]);

                d_r[nz] = pos_r;
                mrd::assign_helper(d_x[nz],tmp);
                ++nz;

                if (pos > or-1)
                {
                    pos_r   = r;
                    pos_c   = c;
                    goto exit_row_flag;
                };

                pos_r       = ptr_r[pos] - 1;
                pos_c       = ptr_c[pos] - 1;

                if (j != pos_c)
                    goto exit_row_flag;	
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

    static void eval_sparse_dc_impl(const SM& A, const md::colon_info& ci, 
                    const value_type_B* ptr_B0, Integer* pos_B, Matrix& ret)
    {
        Integer r = A.rows();
        Integer c = A.cols();
        const raw::details::sparse_ccs<value_type>& Ad = A.rep();

        raw::integer_dense cr = ci.get_rim_r().copy();
        raw::integer_dense cc = ci.get_rim_c().copy();

        sort_rows_cols_stable(cr.ptr(), cc.ptr(), pos_B, cr.size());              

        raw::details::sparse_ccs<value_type> d(A.get_type(),r, c, A.nnz() + cr.size());

        Integer nz				= 0;

        Integer * d_c			= d.ptr_c();
        Integer * d_r			= d.ptr_r();
        value_type * d_x		= d.ptr_x();

        const Integer * Ad_c	= Ad.ptr_c();
        const Integer * Ad_r	= Ad.ptr_r();
        const value_type * Ad_x	= Ad.ptr_x();
        
        Integer or              = cr.size();
        const Integer* ptr_r    = cr.ptr();
        const Integer* ptr_c    = cc.ptr();

        Integer pos             = 0;  
        Integer pos_r           = ptr_r[0] - 1;
        Integer pos_c           = ptr_c[0] - 1;

        for (Integer j = 0; j < c; ++j)
        {
            d_c[j] = nz;

            if (j < pos_c)
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
                if (p < pos_r)
                {
                    d_r[nz] = p;
                    mrd::assign_helper(d_x[nz],Ad_x[i]);
                    ++nz;
                    ++i;
                    continue;
                };	

                while (p > pos_r)
                {
                    Integer elem_pos = pos;

                    while (pos < or && ptr_r[pos] - 1 == pos_r && ptr_c[pos] - 1 == pos_c)
                    {
                        elem_pos = pos;
                        pos += 1;
                    };

                    Integer pos_B_elem  = pos_B[elem_pos];
                    if (pos_B_elem >= 0)
                    {
                        value_type tmp(ptr_B0[pos_B_elem]);

                        d_r[nz] = pos_r;
                        mrd::assign_helper(d_x[nz],tmp);
                        ++nz;
                    };

                    if (pos > or-1)
                    {
                        pos_r   = r;
                        pos_c   = c;
                        goto exit_row_flag;
                    };

                    pos_r       = ptr_r[pos] - 1;
                    pos_c       = ptr_c[pos] - 1;

                    if (j != pos_c)
                        goto exit_row_flag;
                };

                if (p == pos_r)
                {
                    Integer elem_pos = pos;
                    while (pos < or && ptr_r[pos] - 1 == pos_r && ptr_c[pos] - 1 == pos_c)
                    {
                        elem_pos = pos;
                        pos     += 1;
                    };

                    Integer pos_B_elem  = pos_B[elem_pos];
                    if (pos_B_elem >= 0)
                    {
                        value_type tmp(ptr_B0[pos_B_elem]);

                        d_r[nz] = pos_r;
                        mrd::assign_change_helper(d_x[nz],Ad_x[i],tmp);
                        ++nz;
                    };

                    ++i;

                    if (pos > or-1)
                    {
                        pos_r   = r;
                        pos_c   = c;
                        goto exit_row_flag;
                    };

                    pos_r       = ptr_r[pos] - 1;
                    pos_c       = ptr_c[pos] - 1;

                    if (j != pos_c)
                        goto exit_row_flag;	
                };
            };

            while (j == pos_c)
            {
                Integer elem_pos = pos;
                while (pos < or && ptr_r[pos] - 1 == pos_r && ptr_c[pos] - 1 == pos_c)
                {
                    elem_pos = pos;
                    pos     += 1;
                };

                Integer pos_B_elem  = pos_B[elem_pos];
                if (pos_B_elem >= 0)
                {
                    value_type tmp(ptr_B0[pos_B_elem]);

                    d_r[nz] = pos_r;
                    mrd::assign_helper(d_x[nz],tmp);
                    ++nz;
                };

                if (pos > or-1)
                {
                    pos_r   = r;
                    pos_c   = c;
                    goto exit_row_flag;
                };

                pos_r       = ptr_r[pos] - 1;
                pos_c       = ptr_c[pos] - 1;

                if (j != pos_c)
                    goto exit_row_flag;	
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

    static void eval_sparse(const SM& A, const md::colon_info& ci, const M2& B,
                            Matrix& ret)
    {
        using Mat_I         = raw::Matrix<Integer, struct_dense>;

        const auto& Bd      = B.rep();
        const Integer* Bc   = Bd.ptr_c();
        const Integer* Br   = Bd.ptr_r();
        Integer N           = B.cols();

        //make position indices
        Mat_I ind           = Mat_I(ti::ti_empty(), -1, B.rows(), B.cols());
        Integer* ptr_ind    = ind.ptr();
        Integer ind_ld      = ind.ld();
        
        for (Integer c = 0; c < N; ++c)
        {
            for (Integer k = Bc[c]; k < Bc[c + 1]; ++k)
            {
                Integer r   = Br[k];
                ptr_ind[r]  = k;                
            };

            ptr_ind         += ind_ld;
        };

        Integer* pos_B              = ind.ptr();
        const value_type_B* ptr_B0  = Bd.ptr_x();

        return eval_sparse_impl(A, ci, ptr_B0, pos_B, ret);
    };

    static void eval_band(const SM& A, const md::colon_info& ci, const M2& B,
                          Matrix& ret)
    {
        using Mat_I         = raw::Matrix<Integer, struct_dense>;
        
        Integer N           = B.cols();
        Integer B_ld        = B.ld();

        //make position indices
        Mat_I ind           = Mat_I(ti::ti_empty(), -1, B.rows(), B.cols());
        Integer* ptr_ind    = ind.ptr();
        Integer ind_ld      = ind.ld();
        
        for (Integer j = 0, jj = 0; j < N; ++j, jj += B_ld)
        {
            Integer fr      = B.first_row(j);
            Integer lr      = B.last_row(j);
            Integer ii      = jj + B.first_elem_pos(j);

            for (Integer i = fr; i <= lr; ++i,++ii)
                ptr_ind[i]  = ii;                

            ptr_ind         += ind_ld;
        };

        Integer* pos_B              = ind.ptr();
        const value_type_B* ptr_B0  = B.rep_ptr();

        return eval_sparse_impl(A, ci, ptr_B0, pos_B, ret);
    };
};

template<class SM,class M2>
struct change_submatrix_1_impl
{
    static void eval_dense(const SM& A, const md::colon_info& ci, const M2& B,
                           Matrix& ret)
    {
        using value_type = typename SM::value_type;

        Integer r = A.rows(), c = A.cols();
        const raw::details::sparse_ccs<value_type>& Ad = A.rep();

        Integer rs = ci.r_step,rf,rl, pos_v, dpos_v;
        if (rs < 0)
        {
            rs      = -rs;
            rf      = ci.r_end;
            rl      = ci.r_start;
            pos_v   = ci.r_size-1;
            dpos_v  = -1;
        }
        else
        {
            rl      = ci.r_end;
            rf      = ci.r_start;
            pos_v   = 0;
            dpos_v  = 1;
        };

        raw::details::sparse_ccs<value_type> d(A.get_type(),r, c, A.nnz() + ci.rows());

        Integer nz				= 0;

        Integer * d_c			= d.ptr_c();
        Integer * d_r			= d.ptr_r();
        value_type * d_x		= d.ptr_x();

        const Integer * Ad_c	= Ad.ptr_c();
        const Integer * Ad_r	= Ad.ptr_r();
        const value_type * Ad_x	= Ad.ptr_x();
        const M2 B2             = B.make_explicit();
        const typename M2::value_type* ptr_B   = B2.ptr();

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
                while (p >= pos_row)
                {
                    value_type tmp(ptr_B[pos_v]);

                    d_r[nz] = pos_row ;
                        
                    if (p == pos_row)
                        mrd::assign_change_helper(d_x[nz],Ad_x[i],tmp);
                    else
                        mrd::assign_helper(d_x[nz],tmp);
                        
                    ++nz;

                    if (p == pos_row)
                        ++i;

                    pos     += rs;
                    pos_v   += dpos_v;

                    if (pos > rl)
                    {
                        pos_row = r;
                        pos_col = c;
                        break;
                    };

                    md::pos2ind(pos,r,pos_row,pos_col);
                    if (j != pos_col)
                        goto exit_row_flag;
                };
            };

            while (j == pos_col)
            {
                value_type tmp(ptr_B[pos_v]);

                d_r[nz] = pos_row;
                mrd::assign_helper(d_x[nz],tmp);
                ++nz;

                pos     += rs;
                pos_v   += dpos_v;
                if (pos > rl)
                {
                    pos_row = r;
                    pos_col = c;
                    break;
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

    static void eval_sparse(const SM& A, const md::colon_info& ci, const M2& B,
                            Matrix& ret)
    {        
        using value_type    = typename SM::value_type;
        using value_type_B  = typename M2::value_type;

        Integer r = A.rows(), c = A.cols();
        const raw::details::sparse_ccs<value_type>& Ad = A.rep();

        Integer rs = ci.r_step,rf,rl;
        bool decr;
        if (rs < 0)
        {
            rs      = -rs;
            rf      = ci.r_end;
            rl      = ci.r_start;
            decr    = true;
        }
        else
        {
            rl      = ci.r_end;
            rf      = ci.r_start;
            decr    = false;
        };

        raw::details::sparse_ccs<value_type> d(A.get_type(),r, c, A.nnz() + B.nnz());	

        Integer nz				= 0;

        Integer * d_c			= d.ptr_c();
        Integer * d_r			= d.ptr_r();
        value_type * d_x		= d.ptr_x();

        const Integer * Ad_c	= Ad.ptr_c();
        const Integer * Ad_r	= Ad.ptr_r();
        const value_type * Ad_x	= Ad.ptr_x();

        Integer pos_row, pos_col, pos = rf;
        md::pos2ind(rf,r,pos_row,pos_col);

        matcl::Matrix tmp_vB;
        mrd::manip_reshape_helper<M2>::eval_vec(tmp_vB, B);

        const M2& vB = tmp_vB.impl<M2>();

        elem_getter1<value_type_B> B_elem(vB,decr);
        B_elem.set_col(1);

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

                while (p >= pos_row)
                {
                    bool has_nz;
                    value_type tmp(B_elem.get(has_nz));

                    if (has_nz == true)
                    {
                        d_r[nz] = pos_row;
                        if (p == pos_row)
                            mrd::assign_change_helper(d_x[nz],Ad_x[i],tmp);
                        else
                            mrd::assign_helper(d_x[nz],tmp);
                        ++nz;
                    };

                    if (p == pos_row)
                        ++i;

                    pos +=rs;
                    if (pos > rl)
                    {
                        pos_row = r;
                        pos_col = c;
                        break;
                    };

                    md::pos2ind(pos,r,pos_row,pos_col);
                    if (j != pos_col)
                        goto exit_row_flag;
                };
            };

            while (j == pos_col)
            {
                bool has_nz;
                value_type tmp(B_elem.get(has_nz));

                if (has_nz == true)
                {
                    d_r[nz] = pos_row;
                    mrd::assign_helper(d_x[nz],tmp);
                    ++nz;
                };

                ++i;
                pos+= rs;
                if (pos > rl)
                {
                    pos_row = r;
                    pos_col = c;
                    break;
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

    static void eval_banded(const SM& A, const md::colon_info& ci, const M2& B,
                            Matrix& ret)
    {
        using val_type      = typename SM::value_type;
        using SparseMatrix  = raw::Matrix<val_type,struct_sparse>;

        SparseMatrix sB = raw::converter<SM,M2>::eval(B,A.get_type());
        return change_submatrix_1_impl<SM,SparseMatrix>::eval_sparse(A, ci, sB, ret);
    };
};

template<class SM, class M2>
struct change_submatrix_20_impl
{
    static void eval_dense(const SM& A, const md::colon_info& c_in, const M2& B, 
                         Matrix& ret)
    {
        if (imult(c_in.rows(),c_in.cols()) < A.rows())
            return eval_w(A, c_in, B, ret);
        else
            return eval_s(A, c_in, B, ret);
    };

    static void eval_sparse(const SM& A, const md::colon_info& c_in, const M2& B,
                            Matrix& ret)
    {
        if (imult(c_in.rows(),c_in.cols()) < A.rows())
            return eval_sparse_w(A, c_in, B, ret);
        else
            return eval_sparse_s(A, c_in, B, ret);
    };

    static void eval_banded(const SM& A, const md::colon_info& c_in, const M2& B,
                            Matrix& ret)
    {
        using val_type      = typename SM::value_type;
        using SparseMatrix  = raw::Matrix<val_type,struct_sparse>;

        SparseMatrix sB = raw::converter<SM,M2>::eval(B,A.get_type());
        return change_submatrix_20_impl<SM,SparseMatrix>::eval_sparse(A, c_in, sB, ret);
    }

    static void eval_w(const SM& A, const md::colon_info& c_in, const M2& B, Matrix& ret)
    {
        //M2 is dense
        using value_type    = typename SM::value_type;

        Integer r = A.rows(), c = A.cols();
        const raw::details::sparse_ccs<value_type>& Ad = A.rep();

        const Integer * Ad_c		= Ad.ptr_c();
        const Integer * Ad_r		= Ad.ptr_r();
        const value_type * Ad_x		= Ad.ptr_x();	
        Integer B_ld                = B.ld();

        Integer nzA = icast_c(Real(A.nnz()) + Real(c_in.rows())*Real(c_in.cols()));
        raw::details::sparse_ccs<value_type> d(A.get_type(),r, c, nzA);

        Integer * d_c				= d.ptr_c();
        Integer * d_r				= d.ptr_r();
        value_type * d_x			= d.ptr_x();

        column_iterator_2 c_it(c_in);

        Integer nz					= 0;

        raw::integer_dense ri		= c_in.get_rim_2();

        ti::ti_int ti_int;
        raw::integer_dense ind(ti_int);

        sort_type rs_type			= is_sorted(ri);
        bool incr					= true;
        if (rs_type != sorted_increasing)
        {
            raw::integer_dense ri2 = ri.copy();
            ri2.get_struct().reset();

            ind.reset_unique(1,ri.size());
            Integer* ptr_ind    = ind.ptr();
            Integer size        = ri.size();

            for(Integer i = 0; i < size; ++i)
                ptr_ind[i] = i + 1;

            using iterator = md::iterator_helper_2<Integer,Integer>;
            iterator it_begin(ri2.ptr(),ind.ptr(),1);
            iterator it_end(ri2.ptr()+ri2.size(),ind.ptr()+ri2.size(),1);

            static const bool asceding = true;
            std::stable_sort(it_begin,it_end,md::value_ind_compare<Integer,Integer,asceding>());

            ri.assign_to_fresh(std::move(ri2));
            incr = false;
        };

        Integer n = ri.size();
        const Integer* ptr_ri = ri.ptr();
        const typename M2::value_type* ptr_B = B.ptr();
        Integer* ptr_ind = ind.ptr();

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

            Integer m_column = c_it.column();
            c_it.next();
            while (c_it.valid() && c_it.get() == j+1)
            {
                if (c_it.is_increasing())
                    m_column    = c_it.column();

                c_it.next();
            };			

            Integer la			= Ad_c[j+1];
            Integer i			= Ad_c[j];
            Integer col			= imult(m_column-1,B_ld);
            ptr_B               = B.ptr() + col;

            Integer k = 0;
            while(i < la && k < n)
            {
                Integer p		= Ad_r[i];
                Integer pr		= ptr_ri[k]-1;

                while (k < n-1 && ptr_ri[k+1] == ptr_ri[k])
                    ++k;

                if (p < pr)
                {
                    d_r[nz]		= Ad_r[i];
                    mrd::assign_helper(d_x[nz],Ad_x[i]);
                    ++nz;
                    ++i;
                }
                else
                {
                    value_type tmp = md::default_value<value_type>(A.get_type());

                    if (incr)
                        mrd::assign_helper(tmp,value_type(ptr_B[k]));
                    else
                        mrd::assign_helper(tmp,value_type(ptr_B[ptr_ind[k] - 1]));

                    d_r[nz]		= pr;
                    if (p == pr)
                        mrd::assign_change_helper(d_x[nz],Ad_x[i],tmp);
                    else
                        mrd::assign_helper(d_x[nz],tmp);
                    ++nz;

                    ++k;
                    if (p == pr)
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
            
            while(k < n)
            {
                while (k < n-1 && ptr_ri[k+1] == ptr_ri[k])
                    ++k;

                value_type tmp = md::default_value<value_type>(A.get_type());;
                if (incr)
                    mrd::assign_helper(tmp,value_type(ptr_B[k]));
                else
                    mrd::assign_helper(tmp,value_type(ptr_B[ptr_ind[k] - 1]));

                Integer pr	= ptr_ri[k]-1;
                d_r[nz]		= pr;
                mrd::assign_helper(d_x[nz],tmp);
                ++nz;
                ++k;
            };
        };
        d_c[c] = nz;

        ret = Matrix(raw::sparse_matrix_base<value_type>(d), false);
    };

    static void eval_sparse_w(const SM& A, const md::colon_info& c_in, const M2& B, 
                        Matrix& ret)
    {
        using value_type    = typename SM::value_type;
        using value_type_B  = typename M2::value_type;

        Integer r = A.rows(), c = A.cols();
        const raw::details::sparse_ccs<value_type>& Ad = A.rep();
        const raw::details::sparse_ccs<value_type_B>& Bd = B.rep();

        const Integer * Ad_c		= Ad.ptr_c();
        const Integer * Ad_r		= Ad.ptr_r();
        const value_type * Ad_x		= Ad.ptr_x();	

        Integer nzA = icast_c(Real(A.nnz()) + Real(c_in.rows())*Real(c_in.cols()));
        raw::details::sparse_ccs<value_type> d(A.get_type(),r, c, nzA);

        Integer * d_c				= d.ptr_c();
        Integer * d_r				= d.ptr_r();
        value_type * d_x			= d.ptr_x();
        const value_type_B* B_x     = Bd.ptr_x();

        column_iterator_2 c_it(c_in);

        Integer nz					= 0;

        raw::integer_dense ri		= c_in.get_rim_2();
        sort_type rs_type			= is_sorted(ri);
        bool incr					= true;

        ti::ti_int ti_int;
        raw::integer_dense ind(ti_int);

        if (rs_type != sorted_increasing)
        {
            raw::integer_dense ri2 = ri.copy();
            ri2.get_struct().reset();
            
            ind.reset_unique(1,ri.size());
            Integer* ptr_ind    = ind.ptr();
            Integer size        = ri.size();

            for(Integer i = 0; i < size;++i)
                ptr_ind[i]	= i+1;

            using iterator = md::iterator_helper_2<Integer,Integer>;
            iterator it_begin(ri2.ptr(),ind.ptr(),1);
            iterator it_end(ri2.ptr()+ri2.size(),ind.ptr()+ri2.size(),1);

            static const bool asceding = true;
            std::stable_sort(it_begin,it_end,md::value_ind_compare<Integer,Integer,asceding>());

            ri.assign_to_fresh(std::move(ri2));
            incr = false;
        };

        Integer n               = ri.size();
        const Integer* ptr_ri   = ri.ptr();
        const Integer* ptr_ind  = ind.ptr();

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

            Integer m_column = c_it.column();
            c_it.next();

            while (c_it.valid() && c_it.get() == j+1)
            {
                if (c_it.is_increasing())
                    m_column	= c_it.column();

                c_it.next();
            };

            Integer la			= Ad_c[j+1];
            Integer i			= Ad_c[j];
            Integer col			= m_column-1;

            Integer k = 0;
            while(i < la && k < n)
            {
                Integer p		= Ad_r[i];
                Integer pr		= ptr_ri[k]-1;

                while (k < n-1 && ptr_ri[k+1] == ptr_ri[k])
                    ++k;

                if (p < pr)
                {
                    d_r[nz]		= Ad_r[i];
                    mrd::assign_helper(d_x[nz],Ad_x[i]);
                    ++nz;
                    ++i;
                }
                else
                {
                    bool has_nz;
                    Integer pos_k;

                    if (incr)
                        has_nz  = Bd.has_element(k,col, pos_k);
                    else
                        has_nz  = Bd.has_element(ptr_ind[k]-1,col, pos_k);

                    if (has_nz == true)
                    {
                        d_r[nz]		= pr;
                        if (p == pr)
                            mrd::assign_change_helper(d_x[nz],Ad_x[i], value_type(B_x[pos_k]));
                        else
                            mrd::assign_helper(d_x[nz], value_type(B_x[pos_k]));
                        ++nz;
                    };

                    ++k;
                    if (p == pr)
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

            while(k < n)
            {
                while (k < n-1 && ptr_ri[k+1] == ptr_ri[k])
                    ++k;

                bool has_nz;
                Integer pos_k;

                if (incr)
                    has_nz  = Bd.has_element(k,col, pos_k);
                else
                    has_nz  = Bd.has_element(ptr_ind[k]-1,col, pos_k);

                if (has_nz)
                {
                    Integer pr	= ptr_ri[k]-1;
                    d_r[nz]		= pr;
                    mrd::assign_helper(d_x[nz], value_type(B_x[pos_k]));
                    ++nz;
                };

                ++k;
            };
        };

        d_c[c] = nz;

        ret = Matrix(raw::sparse_matrix_base<value_type>(d), false);
    }

    static void eval_s(const SM& A, const md::colon_info& c_in, const M2& B, 
                       Matrix& ret)
    {
        //M2 is dense
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
        const typename M2::value_type* ptr_B = B.ptr();
        Integer B_ld                = B.ld();

        column_iterator_2 c_it(c_in);

        using workspace = md::workspace2<value_type>;
        using scatter   = matcl::algorithm::scatter;

        raw::integer_dense ri		= c_in.get_rim_2();

        Integer size                = ri.size();        
        const Integer* ptr_ri       = ri.ptr();

        workspace		            work_x(A.get_type(),r);
        scatter sc                  = scatter::get(r, c+1+size);
        Integer offset              = sc.current_mark().value() + c + 1;

        for (Integer k = 0; k < size; ++k)
        {
            Integer p			    = ptr_ri[k]-1;
            sc[p]			        = sc.create_mark(k+offset);
        };

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

            Integer m_column = c_it.column();
            c_it.next();
            while (c_it.valid() && c_it.get() == j+1)
            {
                if (c_it.is_increasing())
                    m_column		= c_it.column();

                c_it.next();
            };

            Integer nz_old		= nz;
            bool b_added		= false;
            auto mark           = sc.next_mark();

            Integer col			= imult(m_column-1,B_ld);
            ptr_B               = B.ptr() + col;
            
            for (Integer k = Ad_c[j]; k < Ad_c[j+1]; ++k)
            {
                Integer p		= Ad_r[k];
                if (sc[p] < mark)
                {
                    sc[p]	    = mark;	
                    d_r[nz++]	= p;	
                    mrd::assign_helper(work_x[p],Ad_x[k]);
                };
            };

            for (Integer k = 0; k < size; ++k, ++col)
            {
                Integer p		= ptr_ri[k]-1;
                value_type tmp(ptr_B[k]);
                if (sc[p] == sc.create_mark(k+offset))
                {
                    d_r[nz++]	= p;
                    mrd::assign_helper(work_x[p],tmp);
                    b_added		= true;
                };
            };

            Integer nz_new		= nz - nz_old;

            if (b_added && nz_new > 1)
                utils::sort_q(d_r+nz_old,nz_new);

            for (Integer k = d_c[j]; k < nz ; ++k) 
                mrd::assign_helper(d_x[k],work_x[d_r[k]]);
        };

        d_c[c] = nz;

        ret = Matrix(raw::sparse_matrix_base<value_type>(d), false);
    };
    
    static void eval_sparse_s(const SM& A, const md::colon_info& c_in, const M2& B, 
                              Matrix& ret)
    {
        using value_type    = typename SM::value_type;
        using value_type_B  = typename M2::value_type;

        Integer r = A.rows(), c = A.cols();
        const raw::details::sparse_ccs<value_type>& Ad = A.rep();
        const raw::details::sparse_ccs<value_type_B>& Bd = B.rep();

        const Integer * Ad_c		= Ad.ptr_c();
        const Integer * Ad_r		= Ad.ptr_r();
        const value_type * Ad_x		= Ad.ptr_x();	

        const Integer * Bd_c		= Bd.ptr_c();
        const Integer * Bd_r		= Bd.ptr_r();
        const value_type_B * Bd_x	= Bd.ptr_x();	

        raw::details::sparse_ccs<value_type> d(A.get_type(),r, c, A.nnz() + B.nnz());

        Integer * d_c				= d.ptr_c();
        Integer * d_r				= d.ptr_r();
        value_type * d_x			= d.ptr_x();

        column_iterator_2 c_it(c_in);
        using workspace             = md::workspace2<value_type>;
        using scatter               = matcl::algorithm::scatter;

        raw::integer_dense ri		= c_in.get_rim_2();
        Integer nz					= 0;
        Integer size                = ri.size();
        const Integer* ptr_ri       = ri.ptr();

        workspace		            work_x(A.get_type(),r);
        scatter sc                  = scatter::get(r, c+1+size);        

        Integer offset              = c+1 + sc.current_mark().value();

        for (Integer k = 0; k < size; ++k)
        {
            Integer p			    = ptr_ri[k]-1;
            sc[p]			        = sc.create_mark(offset + k);
        };

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

            Integer m_column = c_it.column();
            c_it.next();
            while (c_it.valid() && c_it.get() == j+1)
            {
                if (c_it.is_increasing())
                    m_column    = c_it.column();

                c_it.next();
            };

            Integer nz_old		= nz;
            bool b_added		= false;
            Integer col			= m_column-1;
            auto mark           = sc.next_mark();

            for (Integer k = Ad_c[j]; k < Ad_c[j+1]; ++k)
            {
                Integer p		= Ad_r[k];
                if (sc[p] < mark)
                {
                    sc[p]	    = mark;	
                    d_r[nz++]	= p;	
                    mrd::assign_helper(work_x[p],Ad_x[k]);
                };
            };

            for (Integer k = Bd_c[col]; k < Bd_c[col+1]; ++k)
            {
                Integer row		= Bd_r[k];
                Integer p		= ptr_ri[row]-1;
                value_type tmp  = value_type(Bd_x[k]);

                if (sc[p] == sc.create_mark(row + offset))
                {
                    d_r[nz++]	= p;
                    mrd::assign_helper(work_x[p],tmp);										
                    b_added		= true;
                };
            };

            Integer nz_new		= nz - nz_old;
            if (b_added && nz_new > 1)
                utils::sort_q(d_r+nz_old,nz_new);

            for (Integer k = d_c[j]; k < nz ; ++k) 
                mrd::assign_helper(d_x[k],work_x[d_r[k]]);
        };

        d_c[c] = nz;

        ret = Matrix(raw::sparse_matrix_base<value_type>(d), false);
    };
};

template<class SM, class M2>
struct change_submatrix_21_impl
{
    static void eval_dense(const SM& A, const md::colon_info& c_in, const M2& B, 
                           Matrix& ret)
    {
        if (c_in.r_start == 1 && c_in.r_end == A.rows() && c_in.r_step == 1)
            return change_cols_dense(A, c_in, B, ret);

        using value_type        = typename SM::value_type;
        using value_type_B      = typename M2::value_type;

        Integer r = A.rows(), c = A.cols();
        const raw::details::sparse_ccs<value_type>& Ad = A.rep();

        const Integer * Ad_c		= Ad.ptr_c();
        const Integer * Ad_r		= Ad.ptr_r();
        const value_type * Ad_x		= Ad.ptr_x();	

        raw::details::sparse_ccs<value_type> d(A.get_type(),r, c, A.nnz() + B.size());

        Integer * d_c				= d.ptr_c();
        Integer * d_r				= d.ptr_r();
        value_type * d_x			= d.ptr_x();
        const value_type_B* ptr_B   = B.ptr();
        Integer B_ld                = B.ld();

        column_iterator_2 c_it(c_in);

        Integer nz					= 0;

        Integer rs = c_in.r_step, rf, rl, pos_v, dpos_v;
        if (rs < 0)
        {
            rs      = -rs;
            rf      = c_in.r_end;
            rl      = c_in.r_start;
            pos_v   = c_in.r_size-1;
            dpos_v  = -1;
        }
        else
        {
            rf      = c_in.r_start;		
            rl      = c_in.r_end;
            pos_v   = 0;
            dpos_v  = 1;
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

            Integer m_column = c_it.column();
            c_it.next();
            while (c_it.valid() && c_it.get() == j+1)
            {
                if (c_it.is_increasing())
                    m_column		= c_it.column();

                c_it.next();
            };

            Integer la			= Ad_c[j+1];
            Integer i			= Ad_c[j];
            Integer col			= imult(m_column-1,B_ld);
            ptr_B               = B.ptr() + col;

            Integer k           = rf;
            Integer pos_B       = pos_v;

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
                else 
                {
                    value_type tmp(ptr_B[pos_B]);

                    d_r[nz]		= k;
                        
                    if (p == k)
                        mrd::assign_change_helper(d_x[nz],Ad_x[i],tmp);
                    else
                        mrd::assign_helper(d_x[nz],tmp);

                    ++nz;

                    if (p == k)
                        ++i;

                    k			+= rs;
                    pos_B		+=dpos_v;
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
                value_type tmp(ptr_B[pos_B]);

                d_r[nz]		= k;
                mrd::assign_helper(d_x[nz],tmp);
                ++nz;

                k				+= rs;
                pos_B			+=dpos_v;
            };
        };
        d_c[c] = nz;

        ret = Matrix(raw::sparse_matrix_base<value_type>(d), false);
    };

    static void eval_sparse(const SM& A, const md::colon_info& c_in, const M2& B,
                            Matrix& ret)
    {
        if (c_in.r_start == 1 && c_in.r_end == A.rows() && c_in.r_step == 1)
            return change_cols_sparse(A, c_in, B, ret);

        using value_type            = typename SM::value_type;
        using value_type_B          = typename M2::value_type;

        Integer r = A.rows(), c = A.cols();
        const raw::details::sparse_ccs<value_type>& Ad = A.rep();

        const Integer * Ad_c		= Ad.ptr_c();
        const Integer * Ad_r		= Ad.ptr_r();
        const value_type * Ad_x		= Ad.ptr_x();	

        raw::details::sparse_ccs<value_type> d(A.get_type(),r, c, A.nnz() + B.nnz());

        Integer * d_c				= d.ptr_c();
        Integer * d_r				= d.ptr_r();
        value_type * d_x			= d.ptr_x();

        column_iterator_2 c_it(c_in);

        Integer nz					= 0;

        Integer rs = c_in.r_step, rf, rl;
        bool decr = false;
        if (rs < 0)
        {
            rs      = -rs;
            rf      = c_in.r_end;
            rl      = c_in.r_start;
            decr    = true;
        }
        else
        {
            rf      = c_in.r_start;		
            rl      = c_in.r_end;
        };
        --rf;
        --rl;

        elem_getter2<value_type_B> eg(B,c_in.r_start-1,c_in.r_step,c_in.r_end-1);
        
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

            Integer m_column = c_it.column();
            c_it.next();
            while (c_it.valid() && c_it.get() == j+1)
            {
                if (c_it.is_increasing())
                    m_column    = c_it.column();

                c_it.next();
            };

            Integer la			= Ad_c[j+1];
            Integer i			= Ad_c[j];
            Integer col			= m_column;
            Integer k			= rf;
            eg.set_col(col,k);

            while(i < la)
            {
                Integer p		= Ad_r[i];

                if (p >= rf)
                    break;

                d_r[nz]			= Ad_r[i];
                mrd::assign_helper(d_x[nz],Ad_x[i]);
                ++nz;
                ++i;
            };			

            while(i < la && k <= rl)
            {
                Integer p		= Ad_r[i];

                if (p < k) 
                {
                    if ((k-p)% rs != 0)
                    {
                        d_r[nz]		= Ad_r[i];
                        mrd::assign_helper(d_x[nz],Ad_x[i]);
                        ++nz;
                    };
                    ++i;
                }
                else if (p == k)
                {
                    value_type tmp(eg.get(k));

                    d_r[nz]		= p;
                    mrd::assign_change_helper(d_x[nz],Ad_x[i],tmp);
                    ++nz;
                    ++i;
                }
                else 
                {
                    Integer k_old	= k;
                    value_type tmp(eg.get(k));

                    d_r[nz]		= k_old;
                    mrd::assign_helper(d_x[nz],tmp);
                    ++nz;
                };
            };

            while(i < la)
            {
                Integer p		= Ad_r[i];

                if (p >= rf && p <= rl)
                {
                    if ((p-rf)% rs != 0)
                    {
                        d_r[nz]		= p;
                        mrd::assign_helper(d_x[nz],Ad_x[i]);
                        ++nz;
                    };
                }
                else
                {
                    d_r[nz]		= p;
                    mrd::assign_helper(d_x[nz],Ad_x[i]);
                    ++nz;
                };
                ++i;
            };
            
            while(k <= rl)
            {
                Integer k_old	= k;
                value_type tmp(eg.get(k));

                d_r[nz]		= k_old;
                mrd::assign_helper(d_x[nz],tmp);
                ++nz;
            };
        };

        d_c[c] = nz;

        ret = Matrix(raw::sparse_matrix_base<value_type>(d), false);
    };

    static void eval_banded(const SM& A, const md::colon_info& ci, const M2& B,
                            Matrix& ret)
    {
        using val_type      = typename SM::value_type;
        using SparseMatrix  = raw::Matrix<val_type,struct_sparse>;

        SparseMatrix sB = raw::converter<SM,M2>::eval(B,A.get_type());
        return change_submatrix_21_impl<SM,SparseMatrix>::eval_sparse(A, ci, sB, ret);
    };
    
    static void change_cols_dense(const SM& A, const md::colon_info& c_in, const M2& B,
                                  Matrix& ret)
    {
        using value_type    = typename SM::value_type;

        Integer r = A.rows(), c = A.cols();
        const raw::details::sparse_ccs<value_type>& Ad = A.rep();

        const Integer * Ad_c		= Ad.ptr_c();
        const Integer * Ad_r		= Ad.ptr_r();
        const value_type * Ad_x		= Ad.ptr_x();	

        raw::details::sparse_ccs<value_type> d(A.get_type(),r, c, A.nnz() + B.size());

        Integer * d_c				= d.ptr_c();
        Integer * d_r				= d.ptr_r();
        value_type * d_x			= d.ptr_x();
        const typename M2::value_type* ptr_B = B.ptr();
        Integer B_ld                = B.ld();

        column_iterator_2 c_it(c_in);		

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

            Integer m_column = c_it.column();
            c_it.next();
            while (c_it.valid() && c_it.get() == j+1)
            {
                if (c_it.is_increasing())
                    m_column		= c_it.column();

                c_it.next();
            };

            Integer col				= imult(m_column-1,B_ld);
            ptr_B                   = B.ptr() + col;

            for (Integer i = 0; i < r; ++i)
            {
                value_type tmp(ptr_B[i]);

                d_r[nz] = i;
                mrd::assign_helper(d_x[nz],tmp);
                ++nz;
            };
        };
        d_c[c] = nz;

        ret = Matrix(raw::sparse_matrix_base<value_type>(d), false);
    };
    
    static void change_cols_sparse(const SM& A, const md::colon_info& c_in, const M2& B,
                                   Matrix& ret)
    {
        using value_type        = typename SM::value_type;
        using value_type_B      = typename M2::value_type;

        Integer r = A.rows(), c = A.cols();
        const raw::details::sparse_ccs<value_type>& Ad = A.rep();
        const raw::details::sparse_ccs<value_type_B>& Bd = B.rep();

        const Integer * Ad_c		= Ad.ptr_c();
        const Integer * Ad_r		= Ad.ptr_r();
        const value_type * Ad_x		= Ad.ptr_x();	

        const Integer * Bd_c		= Bd.ptr_c();
        const Integer * Bd_r		= Bd.ptr_r();
        const value_type_B * Bd_x	= Bd.ptr_x();	

        raw::details::sparse_ccs<value_type> d(A.get_type(),r, c, A.nnz() + B.nnz());

        Integer * d_c				= d.ptr_c();
        Integer * d_r				= d.ptr_r();
        value_type * d_x			= d.ptr_x();

        column_iterator_2 c_it(c_in);		

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

            Integer col = c_it.column();
            c_it.next();
            while (c_it.valid() && c_it.get() == j+1)
            {
                if (c_it.is_increasing())
                    col			= c_it.column();

                c_it.next();
            };			

            for (Integer i = Bd_c[col-1]; i < Bd_c[col]; ++i)
            {
                d_r[nz] = Bd_r[i];
                mrd::assign_helper(d_x[nz],value_type(Bd_x[i]));
                ++nz;
            };
        };
        d_c[c] = nz;

        ret = Matrix(raw::sparse_matrix_base<value_type>(d), false);
    };
};

template<class SM1,class DM2>
void change_submatrix_dense_functor<SM1,DM2>::eval(const SM1& A, const md::colon_info& ci, 
                                                   const DM2& B, Matrix& ret)
{
    if (ci.rows() == 0 || ci.cols() == 0)
    {
        ret = Matrix(A, false);
        return;
    };

    if (ci.r_flag == 0)
        return change_submatrix_20_impl<SM1,DM2>::eval_dense(A, ci, B, ret);
    else
        return change_submatrix_21_impl<SM1,DM2>::eval_dense(A, ci, B, ret);
};

template<class SM1,class SM2>
void change_submatrix_functor<SM1,SM2>::eval(const SM1& A, const md::colon_info& ci, 
                                             const SM2& B, Matrix& ret)
{
    if (ci.rows() == 0)
    {
        ret = Matrix(A, false);
        return;
    }

    if (ci.r_flag == 0)
        return change_submatrix_20_impl<SM1,SM2>::eval_sparse(A, ci, B, ret);
    else
        return change_submatrix_21_impl<SM1,SM2>::eval_sparse(A, ci, B, ret);
};

template<class SM1,class BM2>
void change_submatrix_band_functor<SM1,BM2>::eval(const SM1& A, const md::colon_info& ci, 
                                                  const BM2& B, Matrix& ret)
{
    if (ci.rows() == 0)
    {
        ret = Matrix(A, false);
        return;
    };

    if (ci.r_flag == 0)
        return change_submatrix_20_impl<SM1,BM2>::eval_banded(A, ci, B, ret);
    else
        return change_submatrix_21_impl<SM1,BM2>::eval_banded(A, ci, B, ret);
};

template<class SM1,class DM2>
void change_submatrix_dense_functor_2<SM1,DM2>::eval(const SM1& A, const md::colon_info& ci, 
                                                     const DM2& B, Matrix& ret)
{
    if (ci.rows() == 0)
    {
        ret = Matrix(A, false);
        return;
    };

    if (ci.r_flag == 0)
        return change_submatrix_0_impl<SM1,DM2>::eval_dense(A, ci, B, ret);
    else
        return change_submatrix_1_impl<SM1,DM2>::eval_dense(A, ci, B, ret);
};

template<class SM1,class BM2>
void change_submatrix_band_functor_2<SM1,BM2>::eval(const SM1& A, const md::colon_info& ci, 
                                                    const BM2& B, Matrix& ret)
{
    if (ci.rows() == 0)
    {
        ret = Matrix(A, false);
        return;
    };

    if (ci.r_flag == 0)
        return change_submatrix_0_impl<SM1,BM2>::eval_band(A, ci, B, ret);
    else
        return change_submatrix_1_impl<SM1,BM2>::eval_banded(A, ci, B, ret);
};

template<class SM1,class SM2>
void change_submatrix_functor_2<SM1,SM2>::eval(const SM1& A, const md::colon_info& ci, 
                                               const SM2& B, Matrix& ret)
{
    if (ci.rows() == 0)
    {
        ret = Matrix(A, false);
        return;
    };

    if (ci.r_flag == 0)
        return change_submatrix_0_impl<SM1,SM2>::eval_sparse(A, ci, B, ret);
    else
        return change_submatrix_1_impl<SM1,SM2>::eval_sparse(A, ci, B, ret);
};

template struct change_submatrix_functor<matcl::raw::integer_sparse,matcl::raw::integer_sparse>;
template struct change_submatrix_functor<matcl::raw::real_sparse,matcl::raw::integer_sparse>;
template struct change_submatrix_functor<matcl::raw::real_sparse,matcl::raw::real_sparse>;
template struct change_submatrix_functor<matcl::raw::real_sparse,matcl::raw::float_sparse>;
template struct change_submatrix_functor<matcl::raw::float_sparse,matcl::raw::float_sparse>;
template struct change_submatrix_functor<matcl::raw::complex_sparse,matcl::raw::integer_sparse>;
template struct change_submatrix_functor<matcl::raw::complex_sparse,matcl::raw::real_sparse>;
template struct change_submatrix_functor<matcl::raw::complex_sparse,matcl::raw::float_sparse>;
template struct change_submatrix_functor<matcl::raw::complex_sparse,matcl::raw::float_complex_sparse>;
template struct change_submatrix_functor<matcl::raw::complex_sparse,matcl::raw::complex_sparse>;
template struct change_submatrix_functor<matcl::raw::float_complex_sparse,matcl::raw::float_sparse>;
template struct change_submatrix_functor<matcl::raw::float_complex_sparse,matcl::raw::float_complex_sparse>;
template struct change_submatrix_functor<matcl::raw::object_sparse,matcl::raw::object_sparse>;

template struct change_submatrix_functor_2<matcl::raw::integer_sparse,matcl::raw::integer_sparse>;
template struct change_submatrix_functor_2<matcl::raw::real_sparse,matcl::raw::integer_sparse>;
template struct change_submatrix_functor_2<matcl::raw::real_sparse,matcl::raw::real_sparse>;
template struct change_submatrix_functor_2<matcl::raw::real_sparse,matcl::raw::float_sparse>;
template struct change_submatrix_functor_2<matcl::raw::float_sparse,matcl::raw::float_sparse>;
template struct change_submatrix_functor_2<matcl::raw::complex_sparse,matcl::raw::integer_sparse>;
template struct change_submatrix_functor_2<matcl::raw::complex_sparse,matcl::raw::real_sparse>;
template struct change_submatrix_functor_2<matcl::raw::complex_sparse,matcl::raw::float_sparse>;
template struct change_submatrix_functor_2<matcl::raw::complex_sparse,matcl::raw::complex_sparse>;
template struct change_submatrix_functor_2<matcl::raw::complex_sparse,matcl::raw::float_complex_sparse>;
template struct change_submatrix_functor_2<matcl::raw::float_complex_sparse,matcl::raw::float_sparse>;
template struct change_submatrix_functor_2<matcl::raw::float_complex_sparse,matcl::raw::float_complex_sparse>;
template struct change_submatrix_functor_2<matcl::raw::object_sparse,matcl::raw::object_sparse>;

template struct change_submatrix_dense_functor<matcl::raw::integer_sparse,matcl::raw::integer_dense>;
template struct change_submatrix_dense_functor<matcl::raw::real_sparse,matcl::raw::integer_dense>;
template struct change_submatrix_dense_functor<matcl::raw::real_sparse,matcl::raw::real_dense>;
template struct change_submatrix_dense_functor<matcl::raw::real_sparse,matcl::raw::float_dense>;
template struct change_submatrix_dense_functor<matcl::raw::float_sparse,matcl::raw::float_dense>;
template struct change_submatrix_dense_functor<matcl::raw::complex_sparse,matcl::raw::integer_dense>;
template struct change_submatrix_dense_functor<matcl::raw::complex_sparse,matcl::raw::real_dense>;
template struct change_submatrix_dense_functor<matcl::raw::complex_sparse,matcl::raw::float_dense>;
template struct change_submatrix_dense_functor<matcl::raw::complex_sparse,matcl::raw::complex_dense>;
template struct change_submatrix_dense_functor<matcl::raw::complex_sparse,matcl::raw::float_complex_dense>;
template struct change_submatrix_dense_functor<matcl::raw::float_complex_sparse,matcl::raw::float_dense>;
template struct change_submatrix_dense_functor<matcl::raw::float_complex_sparse,matcl::raw::float_complex_dense>;
template struct change_submatrix_dense_functor<matcl::raw::object_sparse,matcl::raw::object_dense>;

template struct change_submatrix_dense_functor_2<matcl::raw::integer_sparse,matcl::raw::integer_dense>;
template struct change_submatrix_dense_functor_2<matcl::raw::real_sparse,matcl::raw::integer_dense>;
template struct change_submatrix_dense_functor_2<matcl::raw::real_sparse,matcl::raw::real_dense>;
template struct change_submatrix_dense_functor_2<matcl::raw::real_sparse,matcl::raw::float_dense>;
template struct change_submatrix_dense_functor_2<matcl::raw::float_sparse,matcl::raw::float_dense>;
template struct change_submatrix_dense_functor_2<matcl::raw::complex_sparse,matcl::raw::integer_dense>;
template struct change_submatrix_dense_functor_2<matcl::raw::complex_sparse,matcl::raw::real_dense>;
template struct change_submatrix_dense_functor_2<matcl::raw::complex_sparse,matcl::raw::float_dense>;
template struct change_submatrix_dense_functor_2<matcl::raw::complex_sparse,matcl::raw::complex_dense>;
template struct change_submatrix_dense_functor_2<matcl::raw::complex_sparse,matcl::raw::float_complex_dense>;
template struct change_submatrix_dense_functor_2<matcl::raw::float_complex_sparse,matcl::raw::float_dense>;
template struct change_submatrix_dense_functor_2<matcl::raw::float_complex_sparse,matcl::raw::float_complex_dense>;
template struct change_submatrix_dense_functor_2<matcl::raw::object_sparse,matcl::raw::object_dense>;

template struct change_submatrix_band_functor<matcl::raw::integer_sparse,matcl::raw::integer_band>;
template struct change_submatrix_band_functor<matcl::raw::real_sparse,matcl::raw::integer_band>;
template struct change_submatrix_band_functor<matcl::raw::real_sparse,matcl::raw::real_band>;
template struct change_submatrix_band_functor<matcl::raw::real_sparse,matcl::raw::float_band>;
template struct change_submatrix_band_functor<matcl::raw::float_sparse,matcl::raw::float_band>;
template struct change_submatrix_band_functor<matcl::raw::complex_sparse,matcl::raw::integer_band>;
template struct change_submatrix_band_functor<matcl::raw::complex_sparse,matcl::raw::real_band>;
template struct change_submatrix_band_functor<matcl::raw::complex_sparse,matcl::raw::float_band>;
template struct change_submatrix_band_functor<matcl::raw::complex_sparse,matcl::raw::complex_band>;
template struct change_submatrix_band_functor<matcl::raw::complex_sparse,matcl::raw::float_complex_band>;
template struct change_submatrix_band_functor<matcl::raw::float_complex_sparse,matcl::raw::float_band>;
template struct change_submatrix_band_functor<matcl::raw::float_complex_sparse,matcl::raw::float_complex_band>;
template struct change_submatrix_band_functor<matcl::raw::object_sparse,matcl::raw::object_band>;

template struct change_submatrix_band_functor_2<matcl::raw::integer_sparse,matcl::raw::integer_band>;
template struct change_submatrix_band_functor_2<matcl::raw::real_sparse,matcl::raw::integer_band>;
template struct change_submatrix_band_functor_2<matcl::raw::real_sparse,matcl::raw::real_band>;
template struct change_submatrix_band_functor_2<matcl::raw::real_sparse,matcl::raw::float_band>;
template struct change_submatrix_band_functor_2<matcl::raw::float_sparse,matcl::raw::float_band>;
template struct change_submatrix_band_functor_2<matcl::raw::complex_sparse,matcl::raw::integer_band>;
template struct change_submatrix_band_functor_2<matcl::raw::complex_sparse,matcl::raw::real_band>;
template struct change_submatrix_band_functor_2<matcl::raw::complex_sparse,matcl::raw::float_band>;
template struct change_submatrix_band_functor_2<matcl::raw::complex_sparse,matcl::raw::complex_band>;
template struct change_submatrix_band_functor_2<matcl::raw::complex_sparse,matcl::raw::float_complex_band>;
template struct change_submatrix_band_functor_2<matcl::raw::float_complex_sparse,matcl::raw::float_band>;
template struct change_submatrix_band_functor_2<matcl::raw::float_complex_sparse,matcl::raw::float_complex_band>;
template struct change_submatrix_band_functor_2<matcl::raw::object_sparse,matcl::raw::object_band>;

};};};
