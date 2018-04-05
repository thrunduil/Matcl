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

#include "matcl-scalar/IO/formatted_disp.h"
#include "matcl-scalar/IO/formatted_disp_stream.h"

#include <algorithm>

namespace matcl
{

//--------------------------------------------------------------------
//                      column_info
//--------------------------------------------------------------------
column_info::column_info()
    :m_width(0)
{}

column_info::column_info(const std::string& name, align_type align, int width)
    :m_name(name), m_align(align), m_width(std::max(width, 0))
{};

const std::string& column_info::name() const
{
    return m_name;
}

align_type column_info::align() const
{
    return m_align;
}

int column_info::width() const
{
    return m_width;
}

//--------------------------------------------------------------------
//                      formatted_disp
//--------------------------------------------------------------------
formatted_disp::formatted_disp(const disp_stream_ptr& os, const options& opts)
    :m_impl(new details::formatted_disp_stream(os, opts))
{};

void formatted_disp::set_row_label(const std::string& name, align_type al,
                        int width)
{
    m_impl->set_row_label(column_info(name, al, width));
};

void formatted_disp::set_row_label(const column_info& format)
{
    m_impl->set_row_label(format);
}

void formatted_disp::add_column(const std::string& name, align_type al,
                        int width)
{
    return m_impl->add_column(column_info(name, al, width));
}

void formatted_disp::add_column(const column_info& format)
{
    m_impl->add_column(format);
}

bool formatted_disp::display_row_label() const
{
    return m_impl->display_row_label();
}

Integer formatted_disp::number_columns() const
{
    return m_impl->number_columns();
}

const column_info& formatted_disp::get_column_format(Integer col) const
{
    return m_impl->get_column_format(col);
};

const column_info& formatted_disp::get_row_label_format() const
{
    return m_impl->get_row_label_format();
};

void formatted_disp::disp_header()
{
    return m_impl->disp_header();
};

void formatted_disp::disp_row(const std::string& label, formatted_disp_row& row)
{
    m_impl->disp_row_object(label, row.size(), row.data());
    row.clear();
}

void formatted_disp::disp_row_object(const std::string& label, 
                        std::initializer_list<Object> vals)
{
    return m_impl->disp_row_object(label, vals.size(), vals.begin());
}

//--------------------------------------------------------------------
//                      formatted_disp_row
//--------------------------------------------------------------------
void formatted_disp_row::clear()
{
    m_vector.clear();
}

size_t formatted_disp_row::size() const
{
    return m_vector.size();
}

const Object* formatted_disp_row::data() const
{
    return m_vector.data();
}

};
