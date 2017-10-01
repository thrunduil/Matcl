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

#include "matcl-core/IO/disp_data_provider.h"
#include "matcl-core/details/printer.h"

namespace matcl
{

std::string disp_data_provider::get_data_type_description() const
{
    return "";
}
void disp_data_provider::start_display(const disp_stream* orig, line_printer& p) const
{
    return orig->start_display(p);
};
void disp_data_provider::end_display(const disp_stream* orig, line_printer& p) const
{
    return orig->end_display(p);
};
void disp_data_provider::display_empty_matrix(const disp_stream* orig,line_printer& p, Integer r,
                    Integer c) const
{
    return orig->display_empty_matrix(p,r,c);
};
void disp_data_provider::display_matrix_name(const disp_stream* orig, line_printer& p) const
{
    return orig->display_matrix_name(p);
};
void disp_data_provider::start_display_matrix_block(const disp_stream* orig, line_printer& p, 
                    Integer block_width, Integer first_col, Integer last_col) const
{
    return orig->start_display_matrix_block(p,block_width,first_col,last_col);
};
void disp_data_provider::end_display_matrix_block(const disp_stream* orig, line_printer& p, 
                    Integer block_width) const
{
    return orig->end_display_matrix_block(p,block_width);
};
bool disp_data_provider::show_column_header(const disp_stream* orig) const
{
    return orig->show_column_header();
};
bool disp_data_provider::show_row_headers(const disp_stream* orig)  const
{
    return orig->show_row_headers();
};
align_type disp_data_provider::get_align_row_header(const disp_stream* orig) const
{
    return orig->get_align_row_header();
};
align_type disp_data_provider::get_align_col(const disp_stream* orig, Integer c) const
{
    return orig->get_align_col(c);
};
std::string disp_data_provider::get_row_name(const disp_stream* orig, Integer r) const
{
    return orig->get_row_name(r);
};
std::string disp_data_provider::get_rows_label(const disp_stream* orig) const
{
    return orig->get_rows_label();
};
std::string disp_data_provider::get_col_name(const disp_stream* orig, Integer c) const
{
    return orig->get_col_name(c);
};
std::string disp_data_provider::to_string(const disp_stream* user, Integer w, Integer val, 
                                          align_type at) const
{
    matcl::details::bufor_info buf;
    details::printer p(&buf);

    p.set_disp_zero(user->display_zero());
    p.set_precision(user->get_precision());

    p.disp_elem(w, val, at, 0);

    return buf.get_stream().str();
};
std::string disp_data_provider::to_string(const disp_stream* user, Integer w, Real val, align_type at) const
{
    matcl::details::bufor_info buf;
    details::printer p(&buf);

    p.set_disp_zero(user->display_zero());
    p.set_precision(user->get_precision());

    p.disp_elem(w, val, at, 0);

    return buf.get_stream().str();
};
std::string disp_data_provider::to_string(const disp_stream* user, Integer w, const Complex& val, 
                                          align_type at) const
{
    matcl::details::bufor_info buf;
    details::printer p(&buf);

    p.set_disp_zero(user->display_zero());
    p.set_precision(user->get_precision());

    p.disp_elem(w, val,at,0);

    return buf.get_stream().str();
};
std::string disp_data_provider::to_string(const disp_stream* user, Integer w, const std::string& val,
                                align_type at) const
{
    matcl::details::bufor_info buf;
    details::printer p(&buf);

    p.set_disp_zero(user->display_zero());
    p.set_precision(user->get_precision());

    p.disp_elem(w, val,at,0);

    return buf.get_stream().str();
};

};