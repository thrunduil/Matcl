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

#include "matcl-matrep/container/matrix2.inl"
#include "matcl-matrep/visitors/assign_visitor_utils.h"
#include "matcl-internals/error/error_check_basic.h"
#include "matcl-core/error/exception_classes.h"
#include "matcl-core/lib_functions/constants.h"
#include "matcl-matrep/lib_functions/manip.h"
#include "matcl-internals/base/utils.h"
#include "matcl-matrep/details/matrix.inl"

namespace matcl { namespace details
{

struct make_index_2_helper
{
    static const raw::integer_dense&
    make_colon_matrix(Integer size, Integer row, Integer col, 
                bool is_row, const colon& c, colon_info& c_info)
    {
        (void)size;
        (void)c_info;

        const Matrix& mat_ind1          = *c.m_mat_size.m_imat_12;
        const raw::integer_dense& ri    = mat_ind1.get_impl<raw::integer_dense>();

        Integer r               = ri.size();
        const Integer* ptr_ri   = ri.ptr();

        if (is_row)
        {
            for (Integer i = 0; i < r; ++i)
                error::check_row(ptr_ri[i], row, col);
        }
        else
        {
            for (Integer i = 0; i < r; ++i)
                error::check_col(ptr_ri[i], row, col);
        };

        return ri;
    };

    static void make_colon_range(Integer size,Integer rows, Integer cols, bool is_row, 
                                const colon& c, colon_info& ci)
    {        
        Integer st  = 0;
        Integer in  = 0;
        Integer en  = 0;
        Integer r   = 0;

        if (c.m_flag == colon::t_all)
        {
            r               = size;

            if (size)
            {
                st          = 1;
                in          = 1;
                en          = size;                    
            }
        }
        else if (c.m_flag == colon::t_one)
        {
            st          = c.m_s;
            in          = 1;
            en          = c.m_s;
            r           = 1;
        }
        else
        {
            Integer     m_s = 0;
            Integer     m_i = 0;
            Integer     m_e = 0;

            switch (c.m_flag)
            {
                case colon::t_range_simple:
                case colon::t_range_mat:
                {            
                    m_s             = c.m_s; 
                    m_i             = c.m_i; 
                    m_e             = c.m_e; 
                    break;
                }
                case colon::t_end:
                {
                    m_s             = c.m_s; 
                    m_i             = c.m_i; 
                    m_e             = c.m_e + size;
                    break;
                }
                case colon::t_rev_end:
                {
                    m_s             = c.m_s + size; 
                    m_i             = c.m_i; 
                    m_e             = c.m_e; 
                    break;
                }
                case colon::t_end_end:
                {
                    m_s             = c.m_s + size; 
                    m_i             = c.m_i; 
                    m_e             = c.m_e + size; 
                    break;
                }
                case colon::t_last:
                {
                    Integer p       = c.m_e + size;
                    m_s             = p; 
                    m_i             = 1; 
                    m_e             = p; 
                    break;
                }
            };

            if ( m_i == 0 || m_e != m_s && same_sign(m_e - m_s, m_i) == false)
            {
                r           = 0;
            }
            else
            {
                st          = m_s;
                in          = m_i;

                switch (m_i)
                {
                    case 1:
                        en      = m_e;
                        r       = m_e - m_s + 1;
                        break;
                    case -1:
                        en      = m_e;
                        r       = m_s - m_e + 1;
                        break;

                    default:
                    {
                        Integer d,rem;
                        integer_div(m_e - m_s, m_i, d, rem);
                        en      = m_e - rem;
                        r       = d + 1;

                        if (d == 0 && in < 0)
                            in  = -in;
                    }
                };
            };
        };


        if (r)
        {
            if (is_row)
            {
                error::check_row(st, rows, cols);
                error::check_row(en, rows, cols);
            }
            else
            {
                error::check_col(st, rows, cols);
                error::check_col(en, rows, cols);
            };
        }
        else
        {
            st  = 1;
            in  = 1;
            en  = 0;
        };

        if (is_row)
        {
            ci.r_start      = st;
            ci.r_step       = in;
            ci.r_end        = en;
            ci.r_size       = r;
            ci.r_flag       = 1;
            ci.r_rep_size   = r;
        }
        else
        {
            ci.c_start      = st;
            ci.c_step       = in;
            ci.c_end        = en;
            ci.c_size       = r;
            ci.c_flag       = 1;
            ci.c_rep_size   = r;
        };
        return;
    };

    static void eval(Integer rows, Integer cols, const colon& c1, const colon& c2, colon_info& c_info)
    {
        switch (c1.m_flag)
        {            
            case colon::t_matrix1:
            {
                c_info.set_ri(make_colon_matrix(rows,rows,cols,true,c1, c_info));
                c_info.r_flag       = 0;            
                c_info.r_rep_size   = c_info.get_rim_2().size();
                break;
            }
            case colon::t_matrix2:
            {
                throw error::invalid_colon_too_many_mat();
            }
            default:
            {
                make_colon_range(rows,rows,cols, true,c1,c_info);
                break;
            }
        };

        switch (c2.m_flag)
        {            
            case colon::t_matrix1:
            {
                c_info.set_ci(make_colon_matrix(cols,rows,cols, false, c2, c_info));
                c_info.c_flag       = 0;
                c_info.c_rep_size   = c_info.get_cim_2().size();
                break;
            }
            case colon::t_matrix2:
            {
                throw error::invalid_colon_too_many_mat();
            }
            default:
            {
                make_colon_range(cols, rows, cols, false, c2, c_info);
                break;
            }
        };               
    };
};

struct make_index_1_helper
{
    static void make_range(Integer rows, Integer cols, const colon& c,colon_info& c_info)
    {        
        Integer size = imult_c(rows,cols);

        Integer st = 0;
        Integer in = 0;
        Integer en = 0;
        Integer s = 0;
        bool mat = false;

        if (c.m_flag == colon::t_all)
        {
            if (size)
            {
                st      = 1;
                in      = 1;
                en      = size;
                s       = size;
            }
            else
            {
                s       = 0;
            }
        }
        else if (c.m_flag == colon::t_one)
        {
            st      = c.m_s;
            in      = 1;
            en      = c.m_s;
            s       = 1;
        }
        else
        {
            st      = 0;
            in      = 0;
            en      = 0;
            s       = 0;

            Integer m_s = 0;
            Integer m_i = 0;
            Integer m_e = 0;

            switch(c.m_flag)
            {
                case colon::t_range_simple:
                {                
                    m_s     = c.m_s;
                    m_i     = c.m_i;
                    m_e     = c.m_e;
                    break;
                }
                case colon::t_range_mat:
                {                
                    m_s     = c.m_s;
                    m_i     = c.m_i;
                    m_e     = c.m_e;
                    mat     = true;
                    break;
                }
                case colon::t_end:
                {
                    m_s     = c.m_s;
                    m_i     = c.m_i;
                    m_e     = c.m_e + size;
                    break;
                }
                case colon::t_rev_end:
                {
                    m_s     = c.m_s + size;
                    m_i     = c.m_i;
                    m_e     = c.m_e;
                    break;
                }
                case colon::t_end_end:
                {
                    m_s     = c.m_s + size;
                    m_i     = c.m_i;
                    m_e     = c.m_e + size;
                    break;
                }
                case colon::t_last:
                {
                    Integer p = c.m_e + size;
                    m_s     = p;
                    m_i     = 1;
                    m_e     = p;
                    break;
                }
            };

            if ( m_i == 0 || m_e != m_s && same_sign(m_e - m_s, m_i) == false)
            {
                s = 0;
            }
            else
            {
                st          = m_s;
                in          = m_i;

                switch (m_i)
                {
                    case 1:
                        en      = m_e;
                        s       = m_e - m_s + 1;
                        break;
                    case -1:
                        en      = m_e;
                        s       = m_s - m_e + 1;
                        break;
                    default:
                    {
                        Integer d,r;
                        integer_div(m_e - m_s, m_i, d, r);
                        en      = m_e - r;
                        s       = d + 1;

                        if (d == 0 && in < 0)
                            in  = -in;
                    }
                };
            }
        };

        if (s)
        {
            // if colon is a range, then index <= MaxInt
            // in opposite case if As is too large, then exception is already thrown.
            error::check_index(st, size);
            error::check_index(en, size);
        }

        c_info.r_start  = st;
        c_info.r_step   = in;
        c_info.r_end    = en;
        c_info.r_size   = s;
        c_info.r_flag   = 1;
        c_info.c_flag   = -1;

        if (mat == false)
        {
            c_info.r_rep_size = s;
            c_info.c_rep_size = 1;
        }
        else
        {
            c_info.r_rep_size = c.get_matrix_rows();
            c_info.c_rep_size = c.get_matrix_cols();
        };

        return;
    };

    static void eval_matrix_1(Integer rows, Integer cols, const colon& c, colon_info& c_info)
    {
        const Matrix& mat_ind1          = *c.m_mat_size.m_imat_12;
        const raw::integer_dense& ri    = mat_ind1.get_impl<raw::integer_dense>();
        const Integer* ptr_ri           = ri.ptr();

        Integer r   = ri.size();    
        Integer As  = imult_s(rows,cols);
    
        for (Integer i = 0; i < r; ++i)
            error::check_index(ptr_ri[i], As);

        c_info.set_ri(ri);

        c_info.r_flag       = 0;
        c_info.c_flag       = -1;
        c_info.r_rep_size   = ri.rows();
        c_info.c_rep_size   = ri.cols();
    };

    static void eval_matrix_2(Integer rows, Integer cols, const colon& c, colon_info& c_info)
    {
        const Matrix& mat_ind_r         = c.m_mat_size.m_imat_12[0];
        const Matrix& mat_ind_c         = c.m_mat_size.m_imat_12[1];

        const raw::integer_dense& ri_r  = mat_ind_r.get_impl<raw::integer_dense>();
        const raw::integer_dense& ri_c  = mat_ind_c.get_impl<raw::integer_dense>();
        
        const Integer* ptr_ri_r         = ri_r.ptr();
        const Integer* ptr_ri_c         = ri_c.ptr();

        Integer s   = ri_r.size();

        for (Integer i = 0; i < s; ++i)
            error::check_row(ptr_ri_r[i], rows, cols);

        for (Integer i = 0; i < s; ++i)
                error::check_col(ptr_ri_c[i], rows, cols);

        c_info.set_ri_2(ri_r, ri_c);

        c_info.r_flag       = 0;
        c_info.c_flag       = -1;
        c_info.r_rep_size   = ri_r.rows();
        c_info.c_rep_size   = ri_r.cols();
    }

    static void eval(Integer rows, Integer cols, const colon& c, colon_info& c_info)
    {
        switch (c.m_flag)
        {            
            case colon::t_matrix1:
            {
                eval_matrix_1(rows, cols, c, c_info);
                break;
            }
            case colon::t_matrix2:
            {
                eval_matrix_2(rows, cols, c, c_info);
                break;
            }
            default:
            {
                make_range(rows,cols, c,c_info);
                break;
            }
        };
    };
};

void make_index(Integer rows, Integer cols, const colon& c1, const colon& c2, colon_info& ci)
{
    make_index_2_helper::eval(rows,cols, c1, c2, ci);
};

void make_index(Integer rows, Integer cols, const colon& c1, colon_info& ci)
{
    make_index_1_helper::eval(rows,cols, c1, ci);
};

};};
