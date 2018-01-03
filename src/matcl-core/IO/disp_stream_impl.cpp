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

#include "matcl-core/IO/disp_stream.h"
#include "matcl-core/memory/global_objects.h"
#include "matcl-core/details/IO/disp_stream_impl.h"
#include "matcl-core/error/exception_classes.h"

#include "matcl-core/details/integer.h"

#include <iomanip>
#include <algorithm>
#include <iostream>

namespace matcl { namespace details
{

namespace mrd = matcl::raw::details;

class line_printer_impl : public line_printer
{
    private:
        disp_stream_impl*       m_impl;
        disp_stream::line_type  m_lt;
        const disp_stream*      m_user;
        bool                    m_matrix_block;
        Integer                 m_width;

    public:
        line_printer_impl(const disp_stream* user, disp_stream_impl* impl, disp_stream::line_type lt,
                          bool matrix_block, Integer width);

        virtual void    disp_empty_line() override;
        virtual void    disp_raw_string(const std::string& str, align_type at) override;
        virtual void    disp_string(const std::string& str, align_type at) override;
};

line_printer_impl::line_printer_impl(const disp_stream* user, disp_stream_impl* impl, disp_stream::line_type lt,
                                     bool matrix_block, Integer width)
    :m_impl(impl), m_lt(lt), m_user(user), m_matrix_block(matrix_block), m_width(width)
{};
void line_printer_impl::disp_empty_line()
{
    m_impl->do_disp_empty_line(m_user);
};
void line_printer_impl::disp_string(const std::string& str, align_type at)
{
    m_impl->do_display_special(m_user, str, m_lt, m_matrix_block,  m_width,at);
};
void line_printer_impl::disp_raw_string(const std::string& str, align_type at)
{
    Integer width   = m_impl->do_get_terminal_width(m_user);
    if (m_width > 0)
        width       = std::min(width, m_width);

    std::vector<std::string> str_vec = m_impl->split_string(str, width);

    Integer N   = (Integer)str_vec.size();
    for (Integer i = 0; i < N-1; ++i)
    {
        m_impl->get_stream() << m_impl->do_disp_string(width, str_vec[i],at);
        m_impl->get_stream().new_line();
    };

    m_impl->get_stream() << m_impl->do_disp_string(width, str_vec[N],at);
    if (str.back() != '\n')
        m_impl->get_stream().new_line();
};

disp_stream_impl::disp_stream_impl(const output_stream_ptr& os)	
:m_stream(os), m_rows(1), m_cols(1), m_printer(&m_buf)
{
    //these values are not important
    m_terminal_width    = 20;
    m_max_cols          = 10;
    m_max_rows          = 10;
    m_precision         = 0;
    m_ignore_lt         = true;
    m_display_zero      = false;
    m_disp_mode         = disp_mode::standard;
    m_disp_data         = true;
    m_row_block         = 5;
    m_max_nnz           = 1;
    m_do_restrict_sparse= false;
};

void disp_stream_impl::do_init_display(const disp_stream* user, matcl::value_code vt)
{
    m_vt = vt;

    this->init(user);

    m_printer.set_disp_zero(do_display_zero(user));    
    m_printer.set_precision(do_get_precision(user));    
};

void disp_stream_impl::init(const disp_stream* user)
{
    m_terminal_width    = init_terminal_width(user);
    m_max_cols          = std::max(1,user->get_max_cols());
    m_max_rows          = std::max(1,user->get_max_rows());
    m_precision         = std::max(user->get_precision(), 2);
    m_ignore_lt         = user->ignore_lower_triangle();
    m_display_zero      = user->display_zero();
    m_disp_mode         = user->get_disp_mode();
    m_row_block         = user->get_row_block();
    m_disp_data         = user->disp_header_only() == false;
    m_max_nnz           = std::max(1,user->get_max_nnz());
    m_do_restrict_sparse= user->restrict_sparse_matrix_size();

    m_new_line_width    = std::max(0,user->new_line_width(m_terminal_width));
    m_end_line_width    = std::max(0,user->end_line_width(m_terminal_width));

    std::vector<Integer> w{m_new_line_width, m_end_line_width};
    adjust_width(w, m_terminal_width/2,false);

    m_new_line_width    = w[0];
    m_end_line_width    = w[1];
};

Integer disp_stream_impl::init_terminal_width(const disp_stream* user)
{
    Integer w = user->get_terminal_width();

    if (w != -1)
        return w;

    Integer size = m_stream->get_terminal_width();
    if (size <= 0)
    {
        size = m_stream->output_stream::get_terminal_width();
    };

    if (size < 20)
        size = 20;

    return size;
};
Integer disp_stream_impl::do_get_max_matrix_width(const disp_stream* user) const
{
    (void)user;
    Integer w = m_terminal_width - m_new_line_width - m_end_line_width;
    return std::max(w,1);
};
Integer disp_stream_impl::do_get_max_matrix_position(const disp_stream* user) const
{
    (void)user;
    Integer pos = m_terminal_width - m_end_line_width;
    return pos;
};

void disp_stream_impl::do_start_display(const disp_stream* user)
{
    line_printer_impl lp(user, this, disp_stream::line_type::start_display, false, -1);

    user->start_display(lp);
};
void disp_stream_impl::do_display_special(const disp_stream* user, const std::string& str, 
                        disp_stream::line_type type, bool matrix_block, Integer preferred_width,
                        align_type at)
{
    Integer nl_w;
    Integer el_w;

    if (matrix_block == false)
    {
        nl_w    = std::max(user->new_line_special_width(m_terminal_width, type),0);
        el_w    = std::max(0,user->end_line_special_width(m_terminal_width, type));
    }
    else
    {
        nl_w    = std::max(user->new_line_width(m_terminal_width),0);
        el_w    = std::max(0,user->end_line_width(m_terminal_width));
    }

    std::vector<Integer> w{nl_w, el_w};
    adjust_width(w, m_terminal_width/2,false);

    nl_w            = w[0];
    el_w            = w[1];

    Integer txt_w   = m_terminal_width - nl_w - el_w;

    if (preferred_width > 0)
        txt_w       = std::min(txt_w, preferred_width);

    std::string nl;
    std::string el;

    if (matrix_block == false)
    {
        nl      = user->new_line_special(type);
        el      = user->end_line_special(type);
    }
    else
    {
        nl      = user->new_line(0);
        el      = user->end_line(0);
    };

    std::vector<std::string> str_vec = split_string(str, txt_w);

    for (auto pos = str_vec.begin(); pos != str_vec.end(); ++pos)
    {
        get_stream()<< do_disp_string(nl_w, nl,at) 
                    << do_disp_string(txt_w, *pos,at) 
                    << do_disp_string(el_w, el,at);
        get_stream().new_line();
    };
};
std::string disp_stream_impl::do_disp_string(Integer width, const std::string& str, align_type at) const
{
    if (width == 0)
        return "";

    return printer::disp_string(width, str, at);
    /*
    Integer size = (Integer)str.size();

    if (width <= size)
        return printer::disp_string(width, str, at);

    std::string ret = str + std::string(width - size, ' ');
    return ret;
    */
};
void disp_stream_impl::do_disp_empty_line(const disp_stream* user)
{
    (void)user;
    get_stream().new_line();
};
void disp_stream_impl::do_display_matrix_line(const disp_stream* user, const std::string& str,
                                              align_type at)
{
    Integer txt_e   = m_terminal_width - m_new_line_width - (Integer)str.size();

    std::string nl  = user->new_line(0);
    std::string el  = user->end_line(0);

    get_stream()<< do_disp_string(m_new_line_width, nl,at) 
                << do_disp_string((Integer)str.size(), str,at) 
                << do_disp_string(txt_e, el,at);
    get_stream().new_line();
};

void disp_stream_impl::do_end_display(const disp_stream* user)
{
    line_printer_impl lp(user, this, disp_stream::line_type::end_display, false, -1);

    user->end_display(lp);

    do_flush();
    return;
}
void disp_stream_impl::do_flush()
{
    if (m_stream && m_buf.get_size() > 0)
        m_stream->disp(m_buf.get_stream());

    m_buf.clear();
};
void disp_stream_impl::do_displaying_string(const disp_stream* user)
{
    bool disp_header_only = (do_display_data(user) == false);

    if (user->show_matrix_header() == false && disp_header_only == false)
        return;

    do_display_matrix_name(user);

    line_printer_impl lp(user, this, disp_stream::line_type::matrix_type, false, -1);

    user->general_info_string(lp, disp_header_only);
};
void disp_stream_impl::do_display_matrix_name(const disp_stream* user)
{
    line_printer_impl lp(user, this, disp_stream::line_type::matrix_name, false, -1);
    user->display_matrix_name(lp);
};

void disp_stream_impl::do_displaying_scalar(const disp_stream* user, matcl::value_code vt,
                                            const std::string& type_name)
{
    bool disp_header_only = (do_display_data(user) == false);

    if (user->show_matrix_header() == false && disp_header_only == false)
    {
        m_display_zero  = true;
        m_printer.set_disp_zero(true);
        return;
    };

    do_display_matrix_name(user);

    bool as_matrix = false;

    switch(this->do_get_display_mode(user))
    {
        case disp_mode::all_dense:
        case disp_mode::scalar_dense:
            as_matrix = true;
            break;
        default:
        {
            m_printer.set_disp_zero(true);
            break;
        }
    };

    if (disp_header_only == true)
        as_matrix = true;

    line_printer_impl lp(user, this, disp_stream::line_type::matrix_type, false, -1);

    user->general_info_scalar(lp, vt, as_matrix, type_name);
};

void disp_stream_impl::do_displaying_dense_matrix(const disp_stream* user, Integer r, Integer c,
                        matcl::value_code vt, const std::string& struct_name)
{
    m_rows	    = r;
    m_cols	    = c;

    m_column_width.resize(do_max_matrix_cols(user));
    m_column_values_width.resize(do_max_matrix_cols(user));

    bool disp_header_only = (do_display_data(user) == false);

    if (user->show_matrix_header() == false && disp_header_only == false)
        return;

    do_display_matrix_name(user);

    line_printer_impl lp(user, this, disp_stream::line_type::matrix_type, false, -1);
    user->general_info_dense_matrix(lp, r,c,vt, struct_name);
};
void disp_stream_impl::do_displaying_data(const disp_stream* user, Integer r, Integer c,
                            const std::string& description, align_type at)
{
    m_rows	    = r;
    m_cols	    = c;

    bool disp_header_only = (do_display_data(user) == false);

    m_column_width.resize(do_max_matrix_cols(user));
    m_column_values_width.resize(do_max_matrix_cols(user));

    if (user->show_matrix_header() == false && disp_header_only == false)
        return;

    do_display_matrix_name(user);

    if (description.size() == 0)
        return;

    line_printer_impl lp(user, this, disp_stream::line_type::matrix_type, false, -1);
    lp.disp_string(description,at);
};
void disp_stream_impl::do_displaying_banded_matrix2(const disp_stream* user, Integer r, Integer c, 
                        Integer fd, Integer ld, matcl::value_code vt, const std::string& struct_name)
{
    m_rows	= r;
    m_cols	= c;

    m_column_width.resize(do_max_matrix_cols(user));
    m_column_values_width.resize(do_max_matrix_cols(user));

    bool disp_header_only = (do_display_data(user) == false);

    if (user->show_matrix_header() == false && disp_header_only == false)
        return;
    
    do_display_matrix_name(user);

    line_printer_impl lp(user, this, disp_stream::line_type::matrix_type, false, -1);
    user->general_info_banded_matrix(lp, r,c,fd,ld,vt, struct_name);
};

void disp_stream_impl::do_displaying_sparse_matrix(const disp_stream* user, Integer r, Integer c,
                        Integer nz, matcl::value_code vt, const std::string& struct_name)
{
    m_rows	= r;
    m_cols	= c;

    bool disp_header_only = (do_display_data(user) == false);

    if (user->show_matrix_header() == false && disp_header_only == false)
        return;   

    do_display_matrix_name(user);

    line_printer_impl lp(user, this, disp_stream::line_type::matrix_type, false, -1);
    user->general_info_sparse_matrix(lp, r,c,nz,vt, struct_name);
};

void disp_stream_impl::do_display_empty_matrix(const disp_stream* user, Integer r, Integer c)
{
    line_printer_impl lp(user, this, disp_stream::line_type::matrix_type, false, -1);
    user->display_empty_matrix(lp, r,c);
};

bool disp_stream_impl::do_short_print_empty_matrix(const disp_stream* user)
{
    return user->short_print_empty_matrix();
};

void disp_stream_impl::do_start_display_matrix_block(const disp_stream* user, Integer matrix_width,
                            Integer col_start, Integer col_end)
{
    line_printer_impl lp(user, this, disp_stream::line_type::matrix_type, true, matrix_width);
    user->start_display_matrix_block(lp, matrix_width, col_start, col_end);
};
void disp_stream_impl::do_end_display_matrix_block(const disp_stream* user, Integer matrix_width)
{
    line_printer_impl lp(user, this, disp_stream::line_type::matrix_type, true, matrix_width);
    user->end_display_matrix_block(lp, matrix_width);
};
void disp_stream_impl::do_start_display_matrix_block_sparse(const disp_stream* user, Integer matrix_width)
{
    line_printer_impl lp(user, this, disp_stream::line_type::matrix_type, true, matrix_width);
    user->start_display_matrix_block_sparse(lp, matrix_width);
};
void disp_stream_impl::do_end_display_matrix_block_sparse(const disp_stream* user, Integer matrix_width)
{
    line_printer_impl lp(user, this, disp_stream::line_type::matrix_type, true, matrix_width);
    user->end_display_matrix_block_sparse(lp, matrix_width);
};

Integer	disp_stream_impl::do_get_terminal_width(const disp_stream* user) const
{
    (void)user;
    return m_terminal_width;
};

Integer disp_stream_impl::do_get_col_separator_width(const disp_stream* user) const
{
    return static_cast<Integer>(user->get_col_separator().size());
};

Integer disp_stream_impl::do_get_number_subvalues(const disp_stream* user) const
{
    return std::max(0,user->get_number_subvalues());
};

Integer disp_stream_impl::do_get_column_width(const disp_stream* user, Integer c) const
{
    (void)user;
    return m_column_width[c];
};
Integer disp_stream_impl::do_get_first_col_separator_width(const disp_stream* user) const
{
    if (user->show_column_header_row())
    {
        return static_cast<Integer>(user->get_first_col_separator().size());
    }
    else
    {
        return 0;
    };
};

void disp_stream_impl::do_get_row_headers_width(const disp_stream* user, Integer rows,
                                                Integer& w_min, Integer& w_max) const
{
    if (user->show_column_header_row() == false)
    {
        w_min   = 0;
        w_max   = 0;

        return;
    }

    Integer r       = this->do_max_matrix_rows(user);
    do_get_row_header_width(user, w_min, w_max);

    if(user->show_values_labels())
    {
        Integer w   = do_get_label_row_id_width(user);
        w_min       = std::max(w, w_min);
    }

    (void)rows;
    (void) r;
};

Integer disp_stream_impl::do_get_row_header_width_sparse(const disp_stream* user, Integer row) const
{
    return (Integer)user->get_row_name(row).size();
};

Integer disp_stream_impl::do_get_row_header_width_dense(const disp_stream* user, Integer row) const
{
    (void)user;
    return (Integer)this->m_row_labels[row].size();
};

Integer	disp_stream_impl::do_max_matrix_cols(const disp_stream* user) const
{
    (void)user;
    return std::min(m_cols, m_max_cols);
};
Integer	disp_stream_impl::do_max_matrix_rows(const disp_stream* user) const
{
    (void)user;
    return std::min(m_rows, m_max_rows);
};
Integer disp_stream_impl::do_max_matrix_nnz_sparse(const disp_stream* user, Integer matrix_nnz) const
{
    Integer max_nnz = this->do_max_nnz(user);
    return std::min(max_nnz, matrix_nnz);
};
Integer disp_stream_impl::do_max_matrix_cols_sparse(const disp_stream* user, Integer matrix_cols) const
{
    bool restrict_matrix = this->do_restrict_sparse_matrix(user);
    if (restrict_matrix == false)
        return matrix_cols;
    else
        return std::min(matrix_cols, this->do_max_matrix_cols(user));
};
Integer disp_stream_impl::do_max_matrix_rows_sparse(const disp_stream* user, Integer matrix_rows) const
{
    bool restrict_matrix = this->do_restrict_sparse_matrix(user);
    if (restrict_matrix == false)
        return matrix_rows;
    else
        return std::min(matrix_rows, this->do_max_matrix_rows(user));
};

void disp_stream_impl::do_set_column_min_width(const disp_stream* user, Integer c, 
                            const std::vector<Integer>& width_vec, Integer min_allowed,
                            Integer max_allowed)
{
    Integer values_width        = 0;
    Integer separator_width     = 0;

    Integer n_sub               = do_get_number_subvalues(user);

    for(int i = 0; i < n_sub; ++i)
    {
        values_width            += width_vec[i];

        if(i < n_sub - 1)
            values_width        += do_get_col_values_separator_width(user);
    }

    Integer total_width         = values_width + separator_width;
    Integer modified_width      = 0;

    modified_width	            = std::max(total_width,(Integer)this->m_column_labels[c].size());
    
    modified_width	            = std::max(modified_width, min_allowed);
    modified_width              = std::min(modified_width, max_allowed);
    modified_width	            = std::max(modified_width,1);

    m_column_width[c]           = modified_width;

    m_column_values_width[c]	= width_vec;

    if (modified_width < total_width)
    {
        //column header is too large, we must decrease width of values.
        adjust_width(m_column_values_width[c], modified_width - separator_width, true);
    }
    else
    {
        //column header is too large, we must increase width of values.
        m_column_values_width[c][0] += m_column_width[c] - total_width;
    };
};
void disp_stream_impl::do_adjust_sparse_column_min_width(const disp_stream* user, 
                            std::vector<Integer>& width_vec, Integer max_allowed_width) const
{
    Integer values_width        = 0;
    Integer separator_width     = 0;

    Integer n_sub               = do_get_number_subvalues(user);

    for(int i = 0; i < n_sub; ++i)
    {
        values_width            += width_vec[i];

        if(i < n_sub - 1)
        {
            values_width        += do_get_col_values_separator_width(user);
        }
    }

    Integer total_width         = values_width + separator_width;
    Integer modified_width      = 0;

    modified_width	            = std::max(total_width,this->do_get_column_label_sparse_width(user));
    modified_width	            = std::max(modified_width,1);

    Integer res_width           = std::min(modified_width, max_allowed_width);

    if (res_width < total_width)
    {
        //column header is too large, we must decrease width of values.
        adjust_width(width_vec, res_width-separator_width, true);
    }
    else
    {
        //column header is too large, we must increase width of values.
        width_vec[0]            += res_width - total_width;
    };
};
void disp_stream_impl::do_adjust_sparse_row_header(const disp_stream* user, Integer& width_row_header, 
            Integer& width_col_header, Integer& width_header_addin, Integer max_width) const
{
    Integer tot_width   = width_row_header + width_header_addin + width_col_header;

    if (tot_width <= max_width)
        return;
    
    Integer min_width_header_addin  = this->get_min_width_row_header_addins_sparse(user);
    Integer tot_width2  = width_row_header + min_width_header_addin + width_col_header;    

    width_header_addin  = min_width_header_addin;

    if (tot_width2 <= max_width)
        return;

    max_width           = max_width - width_header_addin;

    std::vector<Integer> tmp{width_row_header, width_col_header};
    adjust_width(tmp, max_width, true);

    width_row_header    = tmp[0];
    width_col_header    = tmp[1];
};
void disp_stream_impl::adjust_width(std::vector<Integer>& width, Integer res_width, bool must_be_exact) const
{
    if (res_width <= (Integer)width.size())
    {
        //nothing can be done, set all to 1
        for (size_t i = 0; i < width.size(); ++i)
        {
            width[i]    = 1;
        };
        return;
    };

    Integer values_width        = 0;

    for(size_t i = 0; i < width.size(); ++i)
    {
        values_width            += width[i];
    }

    if (values_width <= res_width && must_be_exact == false)
        return;

    if (res_width >= values_width)
    {
        //we must increase one width
        width[0]                = width[0] + res_width - values_width;
        return;
    };

    if (width.size() == 1)
    {
        width[0]                = res_width;
        return;
    };

    std::vector<Integer*> wp;
    wp.resize(width.size());

    for (size_t i = 0; i < width.size(); ++i)
    {
        wp[i]                   = &width[i];
    };

    //sort from largest to smallest
    struct comparer
    {
        bool operator()(const Integer* a, const Integer* b) const
        {
            return  (*a) > (*b);
        };
    };
    std::sort(wp.begin(), wp.end(), comparer());

    Integer missing_width       = res_width - values_width;

    for (size_t pos = 0; pos < wp.size();++pos)
    {
        Real scale              = Real(res_width) / Real(values_width) / 1.5;
        Integer width_0         = *(wp[pos]);
        Real width_scaled       = std::round(width_0 * scale);
        Integer width_1         = std::max(1,static_cast<Integer>(width_scaled));
        Integer d_width         = width_1 - width_0;

        if (d_width < missing_width)
        {
            *(wp[pos])          = *(wp[pos]) + missing_width;
            break;
        };

        *(wp[pos])              = *(wp[pos]) + d_width;
        missing_width           = missing_width-d_width;
    };

    //missing_width == 0 here
};

static void split(const std::string &s, char delim, std::vector<std::string> &elems) 
{
    std::stringstream ss(s);
    std::string item;
    while (std::getline(ss, item, delim)) 
    {
        elems.push_back(item);
    }
}

static std::vector<std::string> split(const std::string &s, char delim) 
{
    std::vector<std::string> elems;
    split(s, delim, elems);
    return elems;
}

std::vector<std::string> disp_stream_impl::split_string(const std::string& str, Integer width) const
{
    std::vector<std::string> str_vec = split(str, '\n');
    std::vector<std::string> ret;

    for (auto it = str_vec.begin(); it != str_vec.end(); ++it)    
    {
        const std::string& tmp  = *it;

        if ((Integer)tmp.size() <= width)
        {
            ret.push_back(tmp);
            continue;
        };

        Integer pos         = 0;
        Integer rem_width   = (Integer)str.size();

        while (rem_width > 0)
        {
            std::string tmp2 = str.substr(pos, std::min(rem_width, width));
            ret.push_back(std::move(tmp2));

            pos         = pos + width;
            rem_width   = rem_width - width;
        };
    };

    return ret;
};

Integer disp_stream_impl::get_min_width_row_header_addins_sparse(const disp_stream* user) const
{
    (void)user;
    //opening symbol, separator, closing symbol
    return 1+1+1;
};

Integer disp_stream_impl::do_get_column_label_sparse_width(const disp_stream* user) const
{
    (void)user;
    //currently sparse matrices do not have headers
    return 0;
};
void disp_stream_impl::do_preparse_labels(const disp_stream* user, Integer mr, Integer mc,
                                             label_iterators& it)
{
    m_column_labels = std::vector<std::string>();
    m_row_labels    = std::vector<std::string>();

    m_row_labels.reserve(mr);
    m_column_labels.reserve(mc);

    it.begin_col_headers();
    for (Integer i = 0; i < mc; ++i)
    {
        std::string name    = user->get_col_name(i);
        m_column_labels.push_back(name);
        it.next_col_header();
    };

    it.begin_row_headers();
    for (Integer i = 0; i < mr; ++i)
    {        
        std::string name    = user->get_row_name(i);
        m_row_labels.push_back(name);
        it.next_row_header();
    };
};

Integer disp_stream_impl::do_get_column_label_width_sparse(const disp_stream* user, Integer c) const
{
    std::string col_name = user->get_col_name(c);
    return (Integer)col_name.size();
};
const std::vector<Integer>& disp_stream_impl::do_get_column_min_width(Integer col) const
{
    return m_column_values_width[col];
}
Integer	disp_stream_impl::do_get_subvalue_label_width(const disp_stream* user, Integer value_pos) const
{
    return static_cast<Integer>(user->get_subvalue_label(value_pos).size());
};

Integer disp_stream_impl::do_get_sparse_separator_width(const disp_stream* user) const
{
    return std::max(0,user->get_sparse_separator_width());
};
void disp_stream_impl::do_disp_sparse_separator(const disp_stream* user, Integer w)
{
    if (w == 0)
        return;

    std::string sep     = user->get_sparse_separator();

    if (sep.size() == 0)
    {
        get_stream() << std::string(w,' ');
        return;
    };

    if (w <= (int)sep.size())
    {
        get_stream() << sep.substr(0, w);
        return;
    };

    Integer ls = (w-(Integer)sep.size())/2;
    Integer rs = w - ls - (Integer)sep.size();

    get_stream() << std::string(ls,' ') << sep << std::string(rs,' ');
};
void disp_stream_impl::do_disp_new_line(const disp_stream* user, Integer line)
{
    std::string str = user->new_line(line);
    get_stream() << do_disp_string(m_new_line_width, str, align_type::left);
};
void disp_stream_impl::do_disp_end_line(const disp_stream* user, Integer line)
{
    std::string str = user->end_line(line);
    Integer w       = m_terminal_width - m_buf.get_line_pos();
    get_stream() << do_disp_string(w, str,align_type::left);
    get_stream().new_line();;
};
Integer disp_stream_impl::do_get_line_prefix_width(const disp_stream* user) const
{
    (void)user;
    return m_new_line_width;
};

Integer disp_stream_impl::do_get_label_row_id_width(const disp_stream* user) const
{
    return static_cast<Integer>(user->get_labels_row_id().size());
};
void disp_stream_impl::do_disp_labels_row_id(const disp_stream* user,Integer w)
{
    if (user->show_column_header_row() == false)
        return;

    align_type at       = this->do_get_align_row_header(user);
    m_printer.disp_elem(w, user->get_labels_row_id(),at,0);
};
void disp_stream_impl::do_disp_row_header_dense(const disp_stream* user, Integer r, Integer w)
{
    if (user->show_column_header_row())
    {
        std::string name    = m_row_labels[r];
        align_type at       = this->do_get_align_row_header(user);
        m_printer.disp_elem(w,name,at,0);
    };
};
void disp_stream_impl::do_disp_row_header_cont(const disp_stream* user, Integer w)
{
    if (user->show_column_header_row())
    {
        do_disp_continuation_row(user,w);
    };
};
void disp_stream_impl::do_disp_first_col_separator(const disp_stream* user, Integer w)
{
    if (user->show_column_header_row())
    {
        get_stream().disp(w, user->get_first_col_separator());
    };
};

void disp_stream_impl::do_disp_col_separator(const disp_stream* user, Integer w, Integer )
{
    get_stream().disp(w, user->get_col_separator());
};
void disp_stream_impl::do_disp_col_values_separator(const disp_stream* user, Integer w, Integer )
{
    get_stream().disp(w, user->get_col_values_separator());
};
align_type disp_stream_impl::do_get_align_row_header(const disp_stream* user)
{
    return user->get_align_row_header();
};
align_type disp_stream_impl::do_get_align_col(const disp_stream* user, Integer c)
{
    return user->get_align_col(c);
};

void disp_stream_impl::do_get_col_header_width(const disp_stream* user, Integer c, Integer& w_min, 
                                Integer& w_max) const
{
    user->get_column_width(c, w_min, w_max);

    w_min   = std::max(w_min, 1);
    w_max   = std::max(w_max, 1);

    Integer lab_width   = (Integer)m_column_labels[c].size();
    w_min   = std::max(w_min, lab_width);
    w_max   = std::max(w_max, lab_width);
};

void disp_stream_impl::do_get_row_header_width(const disp_stream* user,
                                    Integer& w_min, Integer& w_max) const
{
    if (user->show_column_header_row())
    {
        std::string row_header = user->get_rows_label();
        Integer w0  = (Integer)row_header.size();      
        user->get_column_width_row(w_min, w_max);

        w_min   = std::max(w0, w_min);
        w_max   = std::max(w0, w_max);

        w_min   = std::max(w_min, 1);
        w_max   = std::max(w_max, 1);
    }
    else
    {
        w_min   = 0;
        w_max   = 0;
    };
};
void disp_stream_impl::do_disp_col_header_row(const disp_stream* user, Integer w)
{
    if (user->show_column_header_row())
    {
        std::string row_header = user->get_rows_label();

        align_type at       = this->do_get_align_row_header(user);
        m_printer.disp_elem(w, row_header, at, 0);
    };
};

bool disp_stream_impl::do_show_column_header_line(const disp_stream* user) const
{
    return user->show_column_header_line();
};
bool disp_stream_impl::do_show_column_header_row(const disp_stream* user) const
{
    return user->show_column_header_row();
};
void disp_stream_impl::do_disp_col_header_1sep(const disp_stream* user, Integer w)
{
    if (user->show_column_header_row())
    {
        get_stream().disp(w, user->get_col_separator());
    };
};
void disp_stream_impl::do_disp_col_header_col_dense(const disp_stream* user, Integer c, Integer w)
{
    if (user->show_column_header_columns())
    {
        std::string name    = m_column_labels[c];
        align_type at       = this->do_get_align_col(user,c);
        m_printer.disp_elem(w,name,at, 0);
    };
};

void disp_stream_impl::do_disp_col_header_col_separator(const disp_stream* user, Integer c, Integer w)
{
    (void)c;
    if (user->show_column_header_columns())
    {    
        get_stream().disp(w,user->get_col_separator());
    };
};
void disp_stream_impl::do_disp_first_row_separator(const disp_stream* user, Integer w)
{
    if ((user->show_column_header_columns() || user->show_values_labels()) 
            && user->show_first_row_separator())
    {
        std::string sep = std::string(w,'-');

        do_display_matrix_line(user, sep,align_type::left);
    };
};
void disp_stream_impl::do_disp_col_values_labels(const disp_stream* user, Integer col)
{
    Integer n_sub   = do_get_number_subvalues(user);

    for(int i = 0; i < n_sub; ++i)
    {
        std::string label   = user->get_subvalue_label(i+1);
        Integer w           = m_column_values_width[col][i];
        align_type at       = this->do_get_align_col(user,col);

        m_printer.disp_elem(w, label, at, i);
        
        if(i < n_sub - 1)
        {
            do_disp_col_values_separator(user,do_get_col_values_separator_width(user), i);
        }
    }
};
void disp_stream_impl::do_disp_values_labels(const disp_stream* user, const std::vector<Integer>& widths,
                                             Integer col)
{
    Integer n_sub   = do_get_number_subvalues(user);

    for(int i = 0; i < n_sub; ++i)
    {
        std::string label   = user->get_subvalue_label(i+1);
        align_type at       = this->do_get_align_col(user, col);

        m_printer.disp_elem(widths[i], label, at, i);
        if(i < n_sub - 1)
        {
            do_disp_col_values_separator(user,do_get_col_values_separator_width(user), i);
        }
    }
};
Integer disp_stream_impl::do_get_values_labels_width(const disp_stream* user, const std::vector<Integer>& min_value_width) const
{
    Integer w = 0;
    Integer n_sub   = do_get_number_subvalues(user);

    for(int i = 0; i < n_sub; ++i)
    {
        w       += min_value_width[i];

        if(i < n_sub - 1)
        {
            w   += do_get_col_values_separator_width(user);
        }
    }
    return w;
};

void disp_stream_impl::do_sparse_row(const disp_stream* user, Integer r,Integer c, 
                                     Integer wr, Integer wc)
{
    if (r == -1)
    {
        std::string cd_r(std::min(3,wr),'.');

        if (c == -1)
        {            
            std::string cd_c(std::min(3,wc),'.');
            get_stream() << "(";
            get_stream().disp(wr, cd_r) << ", ";
            get_stream().disp(wc, cd_c) << ")";
        }
        else
        {
            std::string col_name = user->get_col_name(c - 1);
            get_stream() << "(";
            get_stream().disp(wr, cd_r)  << ", ";
            
            m_printer.disp_elem(wc, col_name, align_type::left, 0);

            get_stream() << ")";
        };
    }
    else
    {
        if (c == -1)
        {
            std::string row_name = user->get_row_name(r - 1);

            std::string cd_c(std::min(3,wc),'.');
            get_stream() << "(" ;

            m_printer.disp_elem(wr, row_name, align_type::left, 0);

            get_stream() << ", ";
            get_stream().disp(wc, cd_c) << ")";
        }
        else
        {
            std::string row_name = user->get_row_name(r - 1);
            std::string col_name = user->get_col_name(c - 1);

            get_stream() << "(";
                m_printer.disp_elem(wr, row_name, align_type::left, 0);
            get_stream() << ", ";
                m_printer.disp_elem(wc, col_name, align_type::left, 0);
            get_stream() << ")";
        };	    
    };
};

void disp_stream_impl::do_disp_continuation_value(const disp_stream* user, Integer w)
{
    (void)user;
    std::string cd_v(std::min(3,w),'.');
    m_printer.disp_elem(w,cd_v,align_type::left, 0);
};

void disp_stream_impl::do_disp_continuation_column(const disp_stream* user, Integer w)
{
    if (user->show_column_header_columns())
    {  
        std::string cd_v(std::min(3,w),'.');
        m_printer.disp_elem(w,cd_v,align_type::left, 0);
    };
};
void disp_stream_impl::do_disp_continuation_row(const disp_stream* user, Integer w)
{
    (void)user;
    std::string cd_r(std::min(3,w),'.');
    m_printer.disp_elem(w,cd_r,align_type::left, 0);
};

bool disp_stream_impl::do_symmetric_dispay(const disp_stream* user, bool is_symher) const
{
    if (is_symher == true)
    {
        return this->do_ignore_lower_triangle(user);
    }

    return false;
};

bool disp_stream_impl::do_ignore_lower_triangle(const disp_stream* user) const
{
    (void)user;
    return m_ignore_lt;
};
bool disp_stream_impl::do_display_zero(const disp_stream* user) const
{
    (void)user;
    return m_display_zero;
}
disp_mode disp_stream_impl::do_get_display_mode(const disp_stream* user) const
{
    (void)user;
    return m_disp_mode;
};
bool disp_stream_impl::do_display_data(const disp_stream* user) const
{
    (void)user;
    return m_disp_data;
};

Integer disp_stream_impl::do_get_row_block(const disp_stream* user) const
{
    (void)user;
    return m_row_block;
};


Integer disp_stream_impl::do_max_nnz(const disp_stream* user) const
{
    (void)user;
    return m_max_nnz;
};
bool disp_stream_impl::do_restrict_sparse_matrix(const disp_stream* user) const
{
    (void)user;
    return m_do_restrict_sparse;
};

Integer disp_stream_impl::do_get_value_min_width(const disp_stream* user, 
                    const dynamic::object& v, Integer value_pos) const
{
    (void)user;
    bufor_info os;
    printer(&os).disp_elem(0, v, align_type::left, value_pos);
    return cast_int64(os.get_stream().tellp());
};

Integer disp_stream_impl::do_get_value_min_width(const disp_stream* user, const std::string& v, 
                                                 Integer value_pos) const
{
    (void)user;
    (void)value_pos;
    return (Integer)v.size();
};

Integer disp_stream_impl::do_get_value_min_width(const disp_stream* user, Integer v, Integer value_pos) const
{
    (void)user;
    (void)value_pos;

    if (v == 0)     return 1;
    if (v > 0)
    {
        Integer d = (Integer)std::log10((double)v);
        if (v < 1e9)    return d + 1;
        else            return d + 2;
    }
    else
    {
        v = -v;
        Integer d = (Integer)std::log10((double)v);
        if (v < 1e9)    return d + 2;
        else            return d + 3;
    };
}

Integer disp_stream_impl::do_get_value_min_width(const disp_stream* user, const Float& v, Integer value_pos) const
{
    return do_get_value_min_width(user, Real(v), value_pos);
};
Integer disp_stream_impl::do_get_value_min_width(const disp_stream* user, const Real& v, Integer value_pos) const
{
    (void)value_pos;
    Integer precision = this->do_get_precision(user);

    if (std::isnan(v))     return 3;
    if (std::isinf(v))     return (v > 0)? 3 : 4;

    return printer::get_min_width(v, precision);
}

Integer disp_stream_impl::do_get_value_min_width(const disp_stream* user, const Float_complex& v, Integer value_pos) const
{
    Complex tmp(real(v), imag(v));
    return do_get_value_min_width(user, tmp, value_pos);
};

Integer disp_stream_impl::do_get_value_min_width(const disp_stream* user, const Complex& v, Integer value_pos) const
{
    Real r = real(v);
    Real i = imag(v);

    if (r == 0. && i == 0 && this->do_display_zero(user) == false)
    {
        return 1;
    }

    Integer wr = do_get_value_min_width(user,r, value_pos);
    Integer wc = do_get_value_min_width(user,i, value_pos);

    Integer w = wr + wc;

    if (i >= 0. || std::isnan(i) == true)
        w += 1;

    //additional space for better printing (not really required);
    //make place for additional space between and after +- sign
    w   += 2;

    return w + 1;
}
Integer disp_stream_impl::do_get_precision(const disp_stream* user) const
{
    (void)user;
    return m_precision;
};

printer& disp_stream_impl::get_printer()
{
    return m_printer;
};

void disp_stream_impl::do_final_continuation(const disp_stream* user,Integer line, align_type at)
{
    (void)line;
    if (user->show_final_continuation())
    {
        std::string str = user->get_final_continuation();

        if (str.size() == 0)
            return;

        do_display_matrix_line(user, str,at);
    };
};
void disp_stream_impl::do_disp_begin_next_column_sparse(const disp_stream* user, int w)
{
    if (user->separate_columns_sparse() == true)
    {
        std::string sep = std::string(w,'-');
        do_display_matrix_line(user, sep,align_type::left);
    };
};
Integer disp_stream_impl::do_get_col_values_separator_width(const disp_stream* user) const
{
    return static_cast<Integer>(user->get_col_values_separator().size());
};
bool disp_stream_impl::do_show_values_labels(const disp_stream* user) const
{
    return user->show_values_labels();
};
Integer	disp_stream_impl::do_get_max_value_width(const disp_stream* user, Integer value_pos) const
{
    Integer val = user->get_max_value_width(value_pos+1);
    if (val <= 0)
    {
        Integer w = m_terminal_width - m_new_line_width - m_end_line_width;
        return std::max(w/2, 1);
    }
    else
        return val;
};

}}