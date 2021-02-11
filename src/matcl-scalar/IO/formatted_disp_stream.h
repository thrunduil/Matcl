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

#pragma once

#include "matcl-scalar/config.h"
#include "matcl-scalar/general/fwd_decls.h"
#include "matcl-scalar/IO/formatted_disp.h"
#include "matcl-core/IO/disp_data_provider.h"
#include "matcl-dynamic/details/object.inl"

namespace matcl { namespace details
{

class formatted_disp_stream : public disp_data_provider
{
    private:
        using col_info      = std::vector<column_info>;
        using vec_data      = std::vector<Object>;

    private:
        disp_stream_ptr     m_disp_stream;
        options             m_options;
        column_info         m_row_label;
        col_info            m_cols;

        std::string         m_row_name;
        vec_data            m_data;

        bool                m_header_printed;
        int                 m_col;
        int                 m_col_V1;
        int                 m_col_V2;

    public:
        formatted_disp_stream(const disp_stream_ptr& os = default_disp_stream(),
                        const options& opts = options())
            :m_header_printed(false), m_disp_stream(os), m_options(opts)
        {};

        formatted_disp_stream(const formatted_disp_stream&) = delete;
        formatted_disp_stream& operator=(const formatted_disp_stream&) = delete;

    public:
        void                set_row_label(const column_info& format);
        void                add_column(const column_info& format);

        void                disp_header();

        void                disp_row_object(const std::string& label, size_t size, 
                                const Object* data);

    public:
        // return true if row labels will be displayed
        bool                display_row_label() const;

        // return number of columns to be displayed
        Integer             number_columns() const;

        // return format of column 'col'
        const column_info&  get_column_format(Integer col) const;

        // return format of row labels
        const column_info&  get_row_label_format() const;

    public:
        virtual Integer     cols() const override;
        virtual Integer     rows() const override;

        virtual void        begin() override;           //point to first element of the matrix
        virtual void        hold() override;            //store current point in variable V1
        virtual void        hold_column() override;     //store current point in variable V2
        virtual void        next_row() override;        //go to next row in given column
        virtual void        next_column() override;     //go to next column in given row
        virtual void        restore() override;         //go to the point stored in variable V1
        virtual void        restore_column() override;  //go to the point stored in variable V2

        virtual bool        short_print_empty_matrix(const disp_stream* orig) const override;
        virtual bool        show_matrix_header(const disp_stream* orig) const override;
        virtual bool        show_column_header_line(const disp_stream* orig) const override;
        virtual bool        show_column_header_row(const disp_stream* orig) const override;
        virtual bool        show_column_header_columns(const disp_stream* orig) const override;
        virtual bool        can_split(const disp_stream* orig) const override;
        virtual void        start_display_matrix_block(const disp_stream* orig, line_printer& p, 
                                Integer block_width, Integer first_col, Integer last_col) const override;
        virtual void        end_display_matrix_block(const disp_stream* orig, line_printer& p, 
                                    Integer block_width) const override;
        virtual void        start_display(const disp_stream* orig, line_printer& p) const override;
        virtual void        end_display(const disp_stream* orig, line_printer& p) const override;

        virtual void        begin_row_headers() override;
        virtual void        next_row_header() override;

        virtual std::string get_value(const disp_stream* ds, Integer width, align_type at, 
                                      Integer r, Integer c) const override;

        virtual std::string get_row_name(const disp_stream* orig, Integer r) const;

        virtual std::string get_rows_label(const disp_stream* orig) const override;
        virtual void        get_column_width_row(const disp_stream* orig, Integer& w_min, 
                                Integer& w_max) const override;
        virtual void        get_column_width(const disp_stream* orig, Integer c, Integer& w_min, 
                                Integer& w_max) const override;
        virtual align_type  get_align_row_header(const disp_stream* orig) const override;
        virtual std::string get_col_name(const disp_stream* orig, Integer c) const override;                
        virtual align_type  get_align_col(const disp_stream* orig, Integer c) const override;

    private:
        void                check_column(Integer col) const;
        void                check_row_size(Integer size, Integer req_size) const; 
};

};};
