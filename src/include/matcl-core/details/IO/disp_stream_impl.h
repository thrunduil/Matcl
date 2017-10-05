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

#include "matcl-core/config.h"
#include "matcl-core/general/fwd_decls.h"
#include "matcl-core/memory/global_objects.h"
#include "matcl-core/details/IO/printer.h"

#include <string>
#include <iosfwd>
#include <vector>

#pragma warning(push)
#pragma warning(disable:4251)   // needs to have dll-interface to be used by clients

namespace matcl { namespace details
{

struct MATCL_CORE_EXPORT label_iterators
{
    virtual void        begin_row_headers(){};
    virtual void        next_row_header(){};
    virtual void        begin_col_headers(){};
    virtual void        next_col_header(){};
};

//all indices are 0-based except new_line and end_line which are 1-base (0 - column headers)
class MATCL_CORE_EXPORT disp_stream_impl : public matcl_new_delete
{
    private:
        using string            = std::string;

    private:
        output_stream_ptr       m_stream;
        printer                 m_printer;
        matcl::value_code       m_vt;
        bufor_info              m_buf;
        std::vector<Integer>    m_column_width;
        Integer                 m_rows;
        Integer                 m_cols;

        std::vector<std::string>m_column_labels;
        std::vector<std::string>m_row_labels;

        Integer                 m_terminal_width;
        Integer                 m_max_cols;
        Integer                 m_max_rows;
        Integer                 m_precision;
        bool                    m_ignore_lt;
        bool                    m_display_zero;
        bool                    m_disp_data;
        Integer                 m_row_block;
        disp_mode               m_disp_mode;
        Integer                 m_max_nnz;
        bool                    m_do_restrict_sparse;
        Integer                 m_new_line_width;
        Integer                 m_end_line_width;

        std::vector<std::vector<Integer>>
                                m_column_values_width;

    public:
        explicit disp_stream_impl(const output_stream_ptr& os);

    public:
        //initialization, finalization
        void                do_init_display(const disp_stream* user, matcl::value_code vt);
        void                do_start_display(const disp_stream* user);
        void                do_end_display(const disp_stream* user);

        //preambule
        void                do_displaying_string(const disp_stream* user);
        void                do_displaying_scalar(const disp_stream* user, matcl::value_code vt, 
                                    const std::string& type_name);
        void                do_displaying_data(const disp_stream* user, Integer r, Integer c,
                                    const std::string& description,align_type at);
        void                do_displaying_dense_matrix(const disp_stream* user, Integer r, Integer c,
                                    matcl::value_code vt, const std::string& struct_name);
        void                do_displaying_banded_matrix2(const disp_stream* user, Integer r, Integer c, 
                                    Integer fd, Integer ld, matcl::value_code vt, 
                                    const std::string& struct_name);
        void                do_displaying_sparse_matrix(const disp_stream* user, Integer r, Integer c,
                                    Integer nz, matcl::value_code vt, const std::string& struct_name);
        void                do_display_empty_matrix(const disp_stream* user, Integer r, Integer c);
        bool                do_short_print_empty_matrix(const disp_stream* user);
        void                do_start_display_matrix_block(const disp_stream* user, Integer matrix_width,
                                    Integer col_start, Integer col_end);
        void                do_end_display_matrix_block(const disp_stream* user, Integer matrix_width);
        void                do_start_display_matrix_block_sparse(const disp_stream* user, Integer matrix_width);
        void                do_end_display_matrix_block_sparse(const disp_stream* user, Integer matrix_width);

        //atoms
        void                do_disp_new_line(const disp_stream* user, Integer line);
        void                do_disp_end_line(const disp_stream* user, Integer line);
        void                do_disp_empty_line(const disp_stream* user);
        void                do_display_special(const disp_stream* user, const std::string& str,
                                    disp_stream::line_type type, bool matrix_block, Integer preferred_width,
                                    align_type at);
        void                do_final_continuation(const disp_stream* user, Integer line, align_type at);
        void                do_disp_row_header_dense(const disp_stream* user, Integer r, Integer w);
        void                do_disp_row_header_cont(const disp_stream* user, Integer w);
        void                do_disp_first_col_separator(const disp_stream* user, Integer w);
        void                do_disp_continuation_value(const disp_stream* user, Integer w);
        void                do_disp_col_separator(const disp_stream* user, Integer w, Integer c);
        void                do_get_row_header_width(const disp_stream* user, Integer& w_min, 
                                Integer& w_max) const;
        void                do_get_col_header_width(const disp_stream* user, Integer c, Integer& w_min, 
                                Integer& w_max) const;
        void                do_disp_col_header_row(const disp_stream* user, Integer w);
        void                do_disp_col_header_1sep(const disp_stream* user, Integer w);
        void                do_disp_col_header_col_separator(const disp_stream* user, Integer c, Integer w);                
        void                do_disp_col_header_col_dense(const disp_stream* user, Integer c, Integer w);
        void                do_disp_continuation_column(const disp_stream* user, Integer w);
        void                do_disp_first_row_separator(const disp_stream* user, Integer w);        
        void                do_disp_labels_row_id(const disp_stream* user, Integer w);
        void                do_disp_col_values_labels(const disp_stream* user, Integer col);
        void                do_disp_sparse_separator(const disp_stream* user, Integer w);
        void                do_disp_begin_next_column_sparse(const disp_stream* user, int w);
        void                do_disp_col_values_separator(const disp_stream* user, Integer w, Integer c);
        align_type          do_get_align_row_header(const disp_stream* user);
        align_type          do_get_align_col(const disp_stream* user, Integer c);

        //measurement
        Integer             do_get_col_separator_width(const disp_stream* user)  const;
        void                do_get_row_headers_width(const disp_stream* user, Integer rows,
                                Integer& w_min, Integer& w_max) const;
        Integer             do_get_row_header_width_sparse(const disp_stream* user, Integer row) const;
        Integer             do_get_row_header_width_dense(const disp_stream* user, Integer row) const;
        Integer             do_get_first_col_separator_width(const disp_stream* user) const;        
        Integer             do_get_line_prefix_width(const disp_stream* user) const;
        Integer             do_get_sparse_separator_width(const disp_stream* user) const;
        Integer             do_get_max_matrix_width(const disp_stream* user) const;
        Integer             do_get_max_matrix_position(const disp_stream* user) const;

        //general config
        Integer             do_get_value_min_width(const disp_stream* user, Integer v,
                                Integer value_pos) const;
        Integer             do_get_value_min_width(const disp_stream* user, const Complex& v,
                                Integer value_pos) const;
        Integer             do_get_value_min_width(const disp_stream* user, const Float_complex& v,
                                Integer value_pos) const;
        Integer             do_get_value_min_width(const disp_stream* user, const Real& v, 
                                Integer value_pos) const;
        Integer             do_get_value_min_width(const disp_stream* user, const Float& v,
                                Integer value_pos) const;
        Integer             do_get_value_min_width(const disp_stream* user, const dynamic::object& v, 
                                Integer value_pos) const;
        Integer             do_get_value_min_width(const disp_stream* user, const std::string& v,
                                Integer value_pos) const;

        //access to global options
        //user must be queried only once at the beginning of display
        Integer             do_get_terminal_width(const disp_stream* user) const;       
        Integer             do_max_matrix_cols(const disp_stream* user) const;
        Integer             do_max_matrix_rows(const disp_stream* user) const;
        Integer             do_get_precision(const disp_stream* user) const;
        bool                do_symmetric_dispay(const disp_stream* user, bool is_symher) const;
        bool                do_ignore_lower_triangle(const disp_stream* user) const;
        bool                do_display_zero(const disp_stream* user) const;
        disp_mode           do_get_display_mode(const disp_stream* user) const;
        bool                do_display_data(const disp_stream* user) const;
        Integer             do_max_nnz(const disp_stream* user) const;
        bool                do_restrict_sparse_matrix(const disp_stream* user) const;
        Integer             do_get_row_block(const disp_stream* user) const;

        //displaying objects
        Integer             do_get_number_subvalues(const disp_stream* user) const;
        Integer             do_get_max_value_width(const disp_stream* user, Integer value_pos) const;
        bool                do_show_values_labels(const disp_stream* user) const;
        bool                do_show_column_header_line(const disp_stream* user) const;
        bool                do_show_column_header_row(const disp_stream* user) const;

        Integer             do_get_subvalue_label_width(const disp_stream* user, Integer value_pos) const;
        Integer             do_get_values_labels_width(const disp_stream* user, 
                                    const std::vector<Integer>& min_value_width) const;        
        void                do_disp_values_labels(const disp_stream* user, 
                                    const std::vector<Integer>& widths, Integer col);                
        Integer             do_get_column_label_width_sparse(const disp_stream* user, Integer c) const;
        void                do_preparse_labels(const disp_stream* user, Integer mr, Integer mc,
                                    label_iterators& it);

        //config printers
        void                do_set_column_min_width(const disp_stream* user, Integer c, 
                                    const std::vector<Integer>& width_vec, Integer min_allowed,
                                    Integer max_allowed);        
        Integer             do_get_column_width(const disp_stream* user, Integer c) const;
        bufor_info&         get_stream()            { return m_buf; };
        const std::vector<Integer>& 
                            do_get_column_min_width(Integer col) const;

        Integer             do_get_col_values_separator_width(const disp_stream* user) const;

        //sparse matrices
        Integer             do_max_matrix_cols_sparse(const disp_stream* user, Integer matrix_cols) const;
        Integer             do_max_matrix_rows_sparse(const disp_stream* user, Integer matrix_rows) const;
        Integer             do_max_matrix_nnz_sparse(const disp_stream* user, Integer matrix_nnz) const;        
        void                do_adjust_sparse_row_header(const disp_stream* user, Integer& width_row_header, 
                                    Integer& width_col_header, Integer& width_header_addin, Integer max_width) const;
        void                do_adjust_sparse_column_min_width(const disp_stream* user, 
                                    std::vector<Integer>& width_vec, Integer max_allowed_width) const;
        //if row or col == -1, then insert continuation mark
        void                do_sparse_row(const disp_stream* user, Integer row,Integer col, 
                                          Integer width_r, Integer width_c); 
        Integer             do_get_column_label_sparse_width(const disp_stream* user) const;
        Integer             get_min_width_row_header_addins_sparse(const disp_stream* user) const;

        //utilities
        std::vector<std::string> 
                            split_string(const std::string& str, Integer width) const;
        std::string         do_disp_string(Integer width, const std::string& str, align_type at) const;

    private:        
        Integer             do_get_label_row_id_width(const disp_stream* user) const;        
        void                do_disp_continuation_row(const disp_stream* user, Integer w);                
        void                init(const disp_stream* user);
        Integer             init_terminal_width(const disp_stream* user);
        void                adjust_width(std::vector<Integer>& width_vec, Integer max_width, 
                                bool must_be_exact) const;
        void                do_flush();
        void                do_display_matrix_line(const disp_stream* user, const std::string& str,
                                align_type at);
        void                do_display_matrix_name(const disp_stream* user);        

    public:    
        printer&            get_printer();
        output_stream_ptr   get_output_stream() const  { return m_stream; };

    private:
        disp_stream_impl(const disp_stream_impl&) = delete;
        disp_stream_impl& operator=(const disp_stream_impl&) = delete;
};

}};

#pragma warning(pop)