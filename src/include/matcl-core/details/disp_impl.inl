/*
 *  This file is a part of Matrix Computation Library (MATCL)
 *
 *  Copyright (c) Pawe³ Kowal 2017
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

#include "matcl-core/details/disp_impl.h"

namespace matcl { namespace raw { namespace details
{

template<class V, class S>
void disp_elem_helper<V,S>::eval(md::disp_stream_impl& os, const disp_stream* user, 
                const V& elem, const std::vector<Integer>& widths, Integer col)
{
    Integer cvs     = os.do_get_col_values_separator_width(user);

    Integer n_subv  = os.do_get_number_subvalues(user);
    align_type at   = os.do_get_align_col(user,col);

    for(Integer v = 0; v < n_subv; ++v)
    {                
        Integer w	= widths[v];

        os.get_printer().disp_elem(w, elem, at, v);

        if(v < n_subv - 1)
            os.do_disp_col_values_separator(user,cvs,v);
    }
};

template<class V, class S>
void disp_elem_helper<V,S>::eval(md::disp_stream_impl& os, const disp_stream* user, 
                const std::string& sym, const V& elem, bool is_zero,
                const std::vector<Integer>& widths, Integer col)
{
    (void)elem;

    Integer cvs = os.do_get_col_values_separator_width(user);

    bool display = true;

    if (is_zero == true && os.do_display_zero(user) == false)
        display = false;

    Integer n_subv  = os.do_get_number_subvalues(user);
    align_type at   = os.do_get_align_col(user,col);

    for(Integer v = 0; v < n_subv; ++v)
    {                
        Integer w	= widths[v];                
                
        if (display == true)
            os.get_printer().disp_elem(w,sym,at,v);
        else
            os.get_printer().disp_elem(w," ",at,v);

        if(v < n_subv - 1)
            os.do_disp_col_values_separator(user,cvs,v);
    }
};

//--------------------------------------------------------------------------
//                         impl disp matrix
//--------------------------------------------------------------------------

template<class V>
void disp_matrix<V,struct_dense>::eval_matrix_body(md::disp_stream_impl& os, 
                    const disp_stream* user, matrix_provider_base<V>& m)
{
    Integer r = m.rows();
    Integer c = m.cols();	

    if (r == 0 || c == 0)
    {
        os.do_display_empty_matrix(user,r,c);
        os.do_end_display(user);
        return;
    }

    measures ms;

    Integer max_matrix_width= os.do_get_max_matrix_width(user);
    Integer max_matrix_pos  = os.do_get_max_matrix_position(user);
    bool is_symher          = m.is_symher();
    bool is_sym             = os.do_symmetric_dispay(user, is_symher);

    Integer mc = std::min(c,os.do_max_matrix_cols(user));
    Integer mr = std::min(r,os.do_max_matrix_rows(user));

    os.do_preparse_labels(user, mr, mc, m);

    Integer rhw     = os.do_get_row_headers_width(user,mr);

    for (Integer j = 0; j < mr; ++j)
    {
        Integer tmp = os.do_get_row_header_width_dense(user,j);
        rhw         = std::max(rhw, tmp);
    };

    rhw             = std::min(rhw, max_matrix_width/2);

    ms.csw          = os.do_get_col_separator_width(user);
    ms.rhw          = rhw;
    ms.fcw          = os.do_get_first_col_separator_width(user);
    ms.chw          = 0;

    Integer pos0    = os.do_get_line_prefix_width(user) + ms.rhw + ms.fcw;
    Integer pos     = pos0;            

    Integer lc      = 1;
    Integer cc      = 1;
    Integer n_sub   = os.do_get_number_subvalues(user);
    Integer rbl     = os.do_get_row_block(user);

    m.begin();
    m.hold_column();
    m.hold();            

    for (Integer j = 0; j < mc; ++j, ++cc)
    {
        std::vector<Integer> min_value_width(n_sub, 0);

        //room for column labels
        if(os.do_show_values_labels(user))
        {
            for(Integer v = 0; v < n_sub; ++v)
            {
                min_value_width[v]	= std::max(min_value_width[v],
                                                os.do_get_subvalue_label_width(user,v+1));
            }
        }

        align_type at = os.do_get_align_col(user,j);

        for (Integer i = 0; i < mr; ++i)
        {
            if (is_sym == true && i > j)
            {
                //lower triangle, no need to measure
            }
            else
            {
                bool is_zero;
                const V& elem       = m.get_value(user, 0, at, i, j, is_zero);
                m.next_row();

                for(Integer v = 0; v < n_sub; ++v)
                {
                    Integer vw			= os.do_get_value_min_width(user, elem, v);
                    vw					= std::min(vw,os.do_get_max_value_width(user,v));
                    min_value_width[v]	= std::max(min_value_width[v],vw);
                }
            }
        };

        Integer max_allowed_width   = max_matrix_width - (ms.rhw + ms.fcw);

        os.do_set_column_min_width(user,j,min_value_width, max_allowed_width);
                
        pos                         += os.do_get_column_width(user,j);

        if (pos > max_matrix_pos)
        {
            cc = std::max(cc-1,lc);
            --j;

            m.restore_column();
            m.hold();

            measure_column_header(os, user, ms, lc, cc, false, false);
            os.do_start_display_matrix_block(user, ms.chw, lc, cc);

            make_column_header(os, user, ms, lc, cc, false, false);

            //print matrix
            for (Integer i = 0; i < mr; ++i)
            {
                if (rbl > 0 && i % rbl == 0 && i > 0)
                    make_row_separator(os, user, ms);

                make_row(os, user, ms, i, lc, cc, m, true, false, false, is_sym);

                m.restore();
                m.next_row();
                m.hold();
            }
            if (mr < r)
            {
                make_row(os, user, ms, mr, lc, cc, m, false, false, false, is_sym);
            };

            os.do_end_display_matrix_block(user, ms.chw);

            m.restore_column();
            for (int k = lc; k <= cc; ++k)
            {
                m.next_column();
            };
            m.hold_column();
            m.hold();

            lc      = cc + 1;
            pos     = pos0;
        }
        else
        {
            pos     += ms.csw;

            m.restore();
            m.next_column();
            m.hold();
        }
    };

    cc = mc;
    if (lc <= cc)
    {
        pos += ms.csw;
        pos += 3;

        bool need_cont  = false;
        bool add_cont   = true;
        if (cc < c)
        {             
            need_cont   = true;
            if (pos > max_matrix_pos)
            {
                add_cont = false;
            };
        };

        measure_column_header(os, user, ms, lc, cc, need_cont, add_cont);
        os.do_start_display_matrix_block(user, ms.chw, lc, cc);

        make_column_header(os, user, ms, lc, cc, need_cont, add_cont);

        //print matrix
        m.restore_column();
        m.hold();

        for (Integer i = 0; i < mr; ++i)
        {
            if (rbl > 0 && i % rbl == 0 && i > 0)
                make_row_separator(os, user, ms);

            make_row(os, user, ms, i, lc, cc, m, true, need_cont, add_cont, is_sym);

            m.restore();
            m.next_row();
            m.hold();
        }
        if (mr < r)
        {
            make_row(os, user, ms, mr, lc, cc, m, false, need_cont, add_cont, is_sym);

            m.restore();
            m.next_row();
            m.next_row_header();
            m.hold();
        };                

        if (need_cont == true)
        {
            os.do_final_continuation(user,mr+2, align_type::left);
        };

        os.do_end_display_matrix_block(user, ms.chw);
    };

    os.do_end_display(user);
};

template<class V>
void disp_matrix<V,struct_dense>::measure_column_header(md::disp_stream_impl& os, 
                  const disp_stream* user, measures& ms, Integer lc, Integer cc, 
                  bool need_cont, bool add_cont)
{
    ms.chw = ms.rhw + ms.fcw;

    for (Integer i = lc; i < cc; ++i)
    {
        Integer cw  = os.do_get_column_width(user,i-1);
        ms.chw      += ms.csw + cw;
    };

    if (need_cont == true && add_cont == false)
    {
        ms.chw      += 3;
    }
    else
    {
        Integer cw  = os.do_get_column_width(user,cc-1);
        ms.chw      += cw;                
    };
    if (need_cont == true && add_cont == true)
    {
        ms.chw      += 3 + ms.csw;
    }
};

template<class V>
void disp_matrix<V,struct_dense>::make_column_header(md::disp_stream_impl& os, 
            const disp_stream* user, measures& ms, Integer lc, Integer cc, bool need_cont, 
            bool add_cont)
{
    bool show_row_headers = os.do_show_row_headers(user);

    if (show_row_headers == true)
        os.do_disp_new_line(user,0);

    ms.chw = ms.rhw + ms.fcw;

    os.do_disp_col_header_row(user,ms.rhw);
    os.do_disp_col_header_1sep(user,ms.fcw);

    for (Integer i = lc; i < cc; ++i)
    {
        Integer cw = os.do_get_column_width(user,i-1);
        os.do_disp_col_header_col_dense(user,i-1,cw);
        os.do_disp_col_header_col_separator(user,i-1,ms.csw);

        ms.chw += ms.csw + cw;
    };

    if (need_cont == true && add_cont == false)
    {
        os.do_disp_continuation_column(user,3);
        ms.chw += 3;
    }
    else
    {
        Integer cw = os.do_get_column_width(user,cc-1);
        os.do_disp_col_header_col_dense(user,cc-1,cw);
        ms.chw += cw;                
    };
    if (need_cont == true && add_cont == true)
    {
        os.do_disp_col_header_col_separator(user,cc,ms.csw);
        os.do_disp_continuation_column(user,3);
        ms.chw += 3 + ms.csw;
    }

    if (show_row_headers == true)
        os.do_disp_end_line(user,0);

    os.do_disp_first_row_separator(user,ms.chw);

    if(os.do_show_values_labels(user))
    {
        os.do_disp_new_line(user,0);
        os.do_disp_labels_row_id(user,ms.rhw);
        os.do_disp_first_col_separator(user,ms.fcw);

        for (Integer k = lc; k <= cc; ++k) 
        {
            os.do_disp_col_values_labels(user,k-1);
            if(k < cc)
                os.do_disp_col_separator(user,ms.csw,k-1);
        }
        if (need_cont == true && add_cont == false)
        {
            os.do_disp_col_separator(user,ms.csw,cc-1);
            os.do_disp_continuation_value(user,3);
        }

        os.do_disp_end_line(user,0);
        os.do_disp_first_row_separator(user,ms.chw);
    }
};

template<class V>
void disp_matrix<V,struct_dense>::make_row_separator(md::disp_stream_impl& os, 
                 const disp_stream* user, measures& ms)
{
    os.do_disp_first_row_separator(user,ms.chw);
}

template<class V>
void disp_matrix<V,struct_dense>::make_row(md::disp_stream_impl& os, const disp_stream* user, 
                measures& ms, Integer i, Integer lc, Integer cc, matrix_provider_base<V>& mat,
                bool value, bool need_cont, bool add_cont, bool is_sym)
{
    os.do_disp_new_line(user,i+1);

    if (value == true)
        os.do_disp_row_header_dense(user,i,ms.rhw);
    else
        os.do_disp_row_header_cont(user,ms.rhw);

    os.do_disp_first_col_separator(user,ms.fcw);
            
    for (Integer k = lc; k < cc; ++k) 
    {
        align_type at = os.do_get_align_col(user,k-1);

        if (value == true)
        {
            //currently width is taken correctly only for columns
            //with 0 subvalues, width is used only by disp_data_provider
            Integer w       = os.do_get_column_min_width(k-1)[0];
            bool is_zero    = true;

            if (is_sym == true && i > k-1)
            {
                const V& elem   = mat.get_value(user, w, at, i, k-1, is_zero);
                disp_elem_helper<V,struct_dense>::eval(os,user,"*",elem, is_zero, 
                                    os.do_get_column_min_width(k-1),k-1);
            }
            else
            {                    
                const V& elem   = mat.get_value(user, w, at, i, k-1, is_zero);
                disp_elem_helper<V,struct_dense>::eval(os,user,elem,
                                                os.do_get_column_min_width(k-1), k-1);
            };
        }
        else
        {
            Integer w = os.do_get_column_width(user,k-1);
            os.do_disp_continuation_value(user,w);
        };

        os.do_disp_col_separator(user,ms.csw,k-1);

        mat.next_column();
    };
            
    align_type at = os.do_get_align_col(user,cc-1);

    if (need_cont == true && add_cont == false)
    {
        os.do_disp_continuation_value(user,3);
    }
    else
    {
        if (value == true)
        {
            //currently width is taken correctly only for columns
            //with 0 subvalues, width is used only by disp_data_provider
            Integer w   = os.do_get_column_min_width(cc-1)[0];
            bool is_zero;

            if (is_sym == true && i > cc - 1)
            {
                const V& elem   = mat.get_value(user, w, at, i, cc-1, is_zero);
                disp_elem_helper<V,struct_dense>::eval(os,user,"*", elem, is_zero,
                                                os.do_get_column_min_width(cc-1),cc-1);
            }
            else
            {
                const V& elem   = mat.get_value(user, w, at, i, cc-1, is_zero);
                disp_elem_helper<V,struct_dense>::eval(os,user,elem,
                                                os.do_get_column_min_width(cc-1), cc-1);
            };
            mat.next_column();
        }
        else
        {
            Integer w = os.do_get_column_width(user,cc-1);
            os.do_disp_continuation_value(user,w);
        };
    };
    if (need_cont == true && add_cont == true)
    {
        os.do_disp_col_separator(user,ms.csw,cc-1);
        os.do_disp_continuation_value(user,3);
    };

    os.do_disp_end_line(user,i+1);
};

}}}
