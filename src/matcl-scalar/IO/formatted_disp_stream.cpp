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

#include "matcl-scalar/IO/formatted_disp_stream.h"
#include "matcl-core/IO/base_io.h"

namespace matcl { namespace details
{

void formatted_disp_stream::set_row_label(const column_info& format)
{
    m_row_label = format;
};

void formatted_disp_stream::add_column(const column_info& format)
{
    m_cols.push_back(format);
};

bool formatted_disp_stream::display_row_label() const
{
    return m_row_label.width() > 0;
};

Integer formatted_disp_stream::number_columns() const
{
    return (Integer)m_cols.size();
};

const column_info& formatted_disp_stream::get_row_label_format() const
{
    return m_row_label;
};

// return format of column 'col'
const column_info& formatted_disp_stream::get_column_format(Integer col) const
{
    check_column(col);
    return m_cols[col];
};

std::string formatted_disp_stream::get_rows_label(const disp_stream* orig) const
{
    (void)orig;
    return m_row_label.name();
};

void formatted_disp_stream::get_column_width_row(const disp_stream* orig, Integer& w_min, 
                        Integer& w_max) const
{
    (void)orig;
    w_min = m_row_label.width();
    w_max = m_row_label.width();
}

void formatted_disp_stream::get_column_width(const disp_stream* orig, Integer c, Integer& w_min, 
                        Integer& w_max) const
{
    (void)orig;
    w_min = m_cols[c].width();
    w_max = m_cols[c].width();
}

align_type formatted_disp_stream::get_align_row_header(const disp_stream* orig) const
{
    (void)orig;
    return m_row_label.align();
};

std::string formatted_disp_stream::get_col_name(const disp_stream* orig, Integer c) const
{
    (void)orig;
    return m_cols[c].name();
};

align_type formatted_disp_stream::get_align_col(const disp_stream* orig, Integer c) const
{
    (void)orig;
    return m_cols[c].align();
};

std::string formatted_disp_stream::get_row_name(const disp_stream* orig, Integer r) const
{
    (void)orig;
    (void)r;
    return m_row_name;
}

Integer formatted_disp_stream::cols() const 
{ 
    return (Integer)m_cols.size();
};

Integer formatted_disp_stream::rows() const
{
    return m_header_printed == false ? 0 : 1;
}

void formatted_disp_stream::begin()
{
    m_col   = 0;
};

void formatted_disp_stream::hold()
{
    m_col_V1    = m_col;
};

void formatted_disp_stream::hold_column()
{
    m_col_V2    = m_col;
};

void formatted_disp_stream::next_row()
{};

void formatted_disp_stream::next_column()
{
    ++m_col;
};

void formatted_disp_stream::restore()
{
    m_col       = m_col_V1;
};

void formatted_disp_stream::restore_column()
{
    m_col       = m_col_V2;
};

void formatted_disp_stream::begin_row_headers()
{}

bool formatted_disp_stream::short_print_empty_matrix(const disp_stream* orig) const
{
    (void)orig;
    return false;
};

bool formatted_disp_stream::show_matrix_header(const disp_stream* orig) const
{
    (void)orig;
    return m_header_printed == false ? true : false;
};

bool formatted_disp_stream::show_column_header_line(const disp_stream* orig) const
{
    (void)orig;
    return m_header_printed == false ? true : false;
};

bool formatted_disp_stream::show_column_header_row(const disp_stream* orig) const
{
    (void)orig;
    return this->display_row_label();
}

bool formatted_disp_stream::show_column_header_columns(const disp_stream* orig) const
{
    (void)orig;
    return true;
}

bool formatted_disp_stream::can_split(const disp_stream* orig) const
{
    (void)orig;
    return false;
}

void formatted_disp_stream::start_display_matrix_block(const disp_stream* orig, line_printer& p, 
                Integer block_width, Integer first_col, Integer last_col) const
{
    if (m_header_printed == true)
        return;

    orig->start_display_matrix_block(p, block_width, first_col, last_col);
};

void formatted_disp_stream::end_display_matrix_block(const disp_stream* orig, line_printer& p, 
                    Integer block_width) const
{
    (void)orig;
    (void)p;
    (void)block_width;
};

void formatted_disp_stream::start_display(const disp_stream* orig, line_printer& p) const
{
    if (m_header_printed == true)
        return;

    orig->start_display(p);
}

void formatted_disp_stream::end_display(const disp_stream* orig, line_printer& p) const
{
    (void)orig;
    (void)p;
}

void formatted_disp_stream::next_row_header()
{};

std::string formatted_disp_stream::get_value(const disp_stream* ds, Integer width, align_type at, 
                                Integer r, Integer c) const
{ 
    (void)r;

    const Object& val   = m_data[c];
    return this->to_string(ds, width, val, at);
};

void formatted_disp_stream::disp_header()
{
    if (m_header_printed == true)
        return;

    matcl::disp(*this, m_disp_stream, m_options);
    m_header_printed    = true;
}

void formatted_disp_stream::disp_row_object(const std::string& label, 
                        size_t size, const Object* data)
{
    size_t n    = m_cols.size();

    check_row_size((Integer)size, (Integer)n);

    m_data.resize(n);

    m_row_name  = label;

    for (size_t i = 0; i < n; ++i)
        m_data[i].reset(data[i]);

    m_header_printed = true;
    matcl::disp(*this, m_disp_stream, m_options);
};

void formatted_disp_stream::check_column(Integer col) const
{
    if (col < 0 || col >= number_columns())
        throw error::formatted_disp_invalid_column(col, number_columns());
}

void formatted_disp_stream::check_row_size(Integer size, Integer req_size) const
{
    if (size != req_size)
        throw error::formatted_disp_invalid_row_size(size, req_size);
}

};};
