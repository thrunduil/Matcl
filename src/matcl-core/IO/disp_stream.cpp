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
#include "matcl-core/details/disp_stream_impl.h"
#include "matcl-core/general/exception.h"
#include "matcl-core/details/integer.h"
#include "matcl-core/options/matcl_options.h"
#include "matcl-core/options/options_disp.h"

#include <iomanip>
#include <boost/lexical_cast.hpp>

namespace matcl 
{

namespace mrd = matcl::raw::details;

//------------------------------------------------------------------------
//                          disp_stream
//------------------------------------------------------------------------
disp_stream::disp_stream(const output_stream_ptr& os)
    :m_impl(new details::disp_stream_impl(os))
{};
disp_stream::disp_stream(std::ostream& os)
    :m_impl(new details::disp_stream_impl(output_stream_ptr(new output_stream_from_ostream(os))))
{
};
disp_stream::impl_type disp_stream::impl() const
{
    return m_impl;
}
disp_stream::~disp_stream()
{
};

//------------------------------------------------------------------------
//                          disp_stream_global
//------------------------------------------------------------------------
class disp_stream_global : public disp_stream
{
    public:
        explicit disp_stream_global(const output_stream_ptr& os);

        virtual ~disp_stream_global();

    public:
        //initialization, finalization
        virtual void        start_display(line_printer& p) const override;
        virtual void        end_display(line_printer& p) const override;

        //preambule, new_line and end_line not called at this stage
        virtual void        general_info_string(line_printer& p, bool header_only) const override;
        virtual void        general_info_scalar(line_printer& p, matcl::value_code v, bool as_matrix,
                                    const std::string& type_name) const override;
        virtual void        general_info_dense_matrix(line_printer& p, Integer r, Integer c, 
                                    matcl::value_code vt, const std::string& struct_name) const override;
        virtual void        general_info_banded_matrix(line_printer& p, Integer r, Integer c, Integer fd,
                                    Integer ld, matcl::value_code vt, const std::string& struct_name) const override;
        virtual void        general_info_sparse_matrix(line_printer& p, Integer r, Integer c, Integer nz,
                                    matcl::value_code vt, const std::string& struct_name) const;
        virtual void        display_empty_matrix(line_printer& p, Integer r, Integer c) const override;
        virtual bool        show_matrix_header()  const override;  
        virtual void        display_matrix_name(line_printer& p)  const override;
        virtual void        start_display_matrix_block(line_printer& p, Integer block_width, Integer first_col,
                                    Integer last_col) const  override;
        virtual void        end_display_matrix_block(line_printer& p, Integer block_width) const override;
        virtual void        start_display_matrix_block_sparse(line_printer& p, Integer block_width) const override;
        virtual void        end_display_matrix_block_sparse(line_printer& p, Integer block_width) const override;

        virtual string      get_col_separator() const override;
        virtual bool        show_row_headers()  const override;
        virtual string      get_first_col_separator() const override;        
        virtual Integer     get_precision() const override;

        virtual bool        show_values_labels() const override;
        virtual std::string get_subvalue_label(Integer value_pos) const override;        
        virtual Integer     get_sparse_separator_width() const override;        
        virtual string      get_labels_row_id() const override;        
        virtual string      get_row_name(Integer r) const override;
        virtual std::string get_rows_label() const override;
        virtual string      get_col_values_separator() const override;
        virtual bool        show_column_header() const override;
        virtual string      get_col_name(Integer c) const override;
        virtual bool        show_first_row_separator() const override;
        virtual bool        show_final_continuation() const override;
        virtual std::string get_final_continuation() const override;
        virtual Integer     get_number_subvalues() const override;
        virtual align_type  get_align_row_header() const override;
        virtual align_type  get_align_col(Integer c) const override;

        virtual std::string new_line(Integer line) const override;
        virtual std::string end_line(Integer line) const override;
        virtual Integer     new_line_width(Integer terminal_width) const override;
        virtual Integer     end_line_width(Integer terminal_width) const override;
        virtual std::string new_line_special(line_type type) const override;
        virtual std::string end_line_special(line_type type) const override;
        virtual Integer     new_line_special_width(Integer terminal_width, line_type type) const override;
        virtual Integer     end_line_special_width(Integer terminal_width, line_type type) const override;
        virtual bool        separate_columns_sparse() const override;
        virtual bool        disp_header_only() const override;
        virtual string      get_sparse_separator() const override;
        virtual bool        display_zero() const override;
        virtual bool        ignore_lower_triangle() const override;
        virtual disp_mode   get_disp_mode() const override;
        virtual Integer     get_row_block() const override;

                            //-1 if unknown
        virtual Integer     get_terminal_width() const override;
                            //-1 if not set
        virtual Integer     get_max_value_width(Integer value_pos) const override;        

        virtual Integer     get_max_cols() const override;
        virtual Integer     get_max_rows() const override;
        virtual Integer     get_max_nnz() const override;
        virtual bool        restrict_sparse_matrix_size() const override;

        virtual std::string get_value_name(matcl::value_code vt) const override;
};

namespace mrd = matcl::raw::details;

disp_stream_global::disp_stream_global(const output_stream_ptr& os)
    :disp_stream(os)
{};
disp_stream_global::~disp_stream_global()
{};

std::string disp_stream_global::get_value_name(matcl::value_code vt) const
{
    switch (vt)
    {
        case value_code::v_integer:         return "integer";
        case value_code::v_float:           return "float";
        case value_code::v_real:            return "real";
        case value_code::v_float_complex:   return "float complex";
        case value_code::v_complex:         return "complex";
        case value_code::v_object:          return "object";
        default:
        {
            matcl_assert(0,"unknown case");
            throw;
        }
    };
};

bool disp_stream_global::show_matrix_header() const
{
    return true;
};

void disp_stream_global::general_info_string(line_printer& p, bool header_only) const
{
    if (header_only == true)
    {
        std::string ret = "string";

        p.disp_empty_line();
        p.disp_string(ret,align_type::left);
        return;
    }

    return;
};

void disp_stream_global::general_info_scalar(line_printer& p, matcl::value_code vt, bool as_matrix,
                                             const std::string& type_name) const
{
    if (as_matrix == true)
    {
        std::string v   = this->get_value_name(vt);
        std::string ret = v + " scalar";

        if (vt == matcl::value_code::v_object)
            ret += std::string(" ") + type_name;

        p.disp_empty_line();
        p.disp_string(ret,align_type::left);
        return;
    }
    else
    {
        return;
    };
};
void disp_stream_global::general_info_dense_matrix(line_printer& p, Integer r, Integer c, matcl::value_code vt,
                                          const std::string& struct_name) const
{
    std::string v = this->get_value_name(vt);

    std::ostringstream buf;
    buf << "dense " << v << " matrix, size: " << r << 'x' << c << ", type: " << struct_name;

    std::string str = buf.str();

    p.disp_empty_line();
    p.disp_string(str,align_type::left);

    return;
};
void disp_stream_global::general_info_banded_matrix(line_printer& p, Integer r, Integer c, Integer fd, Integer ld, 
                                           matcl::value_code vt, const std::string& struct_name) const
{
    std::string v = get_value_name(vt);

    std::ostringstream buf;
    buf << "banded " << v << " matrix, size: " << r << 'x' << c << ", fd: " << fd << ", ld: " << ld 
        << ", type: " << struct_name;

    std::string str = buf.str();

    p.disp_empty_line();
    p.disp_string(str,align_type::left);

    return;
};
void disp_stream_global::general_info_sparse_matrix(line_printer& p, Integer r, Integer c, Integer nz, 
                            matcl::value_code vt, const std::string& struct_name) const
{
    std::string v = get_value_name(vt);

    std::ostringstream buf;
    buf         << "sparse " << v << " matrix, size: " << r << 'x' << c << ", nnz: " << nz 
                << ", type: " << struct_name;

    std::string str = buf.str();

    p.disp_empty_line();
    p.disp_string(str,align_type::left);

    return;
};
void disp_stream_global::display_empty_matrix(line_printer& p, Integer r, Integer c) const
{
    std::ostringstream buf;
    buf << "[](" << r << 'x' << c << ")";

    std::string str = buf.str();

    p.disp_empty_line();
    p.disp_string(str,align_type::left);
    return;
};

bool disp_stream_global::show_final_continuation() const		
{ 
    return false; 
};
std::string disp_stream_global::get_final_continuation() const
{ 
    return ""; 
};
Integer disp_stream_global::get_number_subvalues() const
{
    return 1;
};
align_type disp_stream_global::get_align_row_header() const
{
    return align_type::right;
};
align_type disp_stream_global::get_align_col(Integer c) const
{
    (void)c;
    return align_type::right;
};

Integer	disp_stream_global::get_terminal_width() const
{
    return options().get_option<Integer>(matcl::opt::disp::terminal_width());
};
bool disp_stream_global::display_zero() const
{
    return options().get_option<bool>(matcl::opt::disp::display_zero());
};
bool disp_stream_global::ignore_lower_triangle() const
{
    return options().get_option<bool>(matcl::opt::disp::ignore_lower_triangle());
};
disp_mode disp_stream_global::get_disp_mode() const
{
    return static_cast<disp_mode>(options().get_option<Integer>(matcl::opt::disp::disp_mode()));
};
Integer disp_stream_global::get_row_block() const
{
    return options().get_option<Integer>(matcl::opt::disp::row_block());
}

std::string disp_stream_global::get_col_separator() const
{
    return " | ";
};
std::string disp_stream_global::get_col_values_separator() const
{
    return " : ";
};

Integer	disp_stream_global::get_max_value_width(Integer value_pos) const
{
    (void)value_pos;
    return -1;
};

Integer disp_stream_global::get_sparse_separator_width() const
{
    return 3;
};
std::string disp_stream_global::get_sparse_separator() const
{
    return "|";
};
void disp_stream_global::start_display(line_printer& p) const
{
    (void)p;
};
void disp_stream_global::end_display(line_printer& p) const
{
    (void)p;
    return;
}
std::string disp_stream_global::end_line(Integer) const
{
    return "";
};
std::string disp_stream_global::new_line(Integer) const
{
    return " ";
};
Integer disp_stream_global::new_line_width(Integer terminal_width) const
{
    (void)terminal_width;
    return 1;
};
Integer disp_stream_global::end_line_width(Integer terminal_width) const
{
    (void)terminal_width;
    return 0;
};
std::string disp_stream_global::new_line_special(line_type type) const
{
    (void)type;
    return " ";
};
std::string disp_stream_global::end_line_special(line_type type) const
{
    (void)type;
    return "";
};
Integer disp_stream_global::new_line_special_width(Integer terminal_width, line_type type) const
{
    (void)terminal_width;
    (void)type;
    return 1;
};
Integer disp_stream_global::end_line_special_width(Integer terminal_width, line_type type) const
{
    (void)terminal_width;
    (void)type;
    return 0;
};

disp_stream_global::string disp_stream_global::get_labels_row_id() const
{
    return "L";
};
bool disp_stream_global::show_row_headers() const
{
    return true;
};
std::string disp_stream_global::get_first_col_separator() const
{
    return " | ";
};
bool disp_stream_global::show_column_header() const
{
    return true;
};
bool disp_stream_global::show_first_row_separator() const
{
    return true;
};
std::string disp_stream_global::get_subvalue_label(Integer) const
{
    return std::string();
};
bool disp_stream_global::show_values_labels() const
{
    return false;
};

void disp_stream_global::display_matrix_name(line_printer& p) const
{
    (void)p;
    return;
};
void disp_stream_global::start_display_matrix_block(line_printer& p, Integer block_width, Integer first_col,
                            Integer last_col) const
{
    (void)block_width;
    (void)first_col;
    (void)last_col;
    p.disp_empty_line();
}
void disp_stream_global::end_display_matrix_block(line_printer& p, Integer block_width) const
{
    (void)p;
    (void)block_width;
};
void disp_stream_global::start_display_matrix_block_sparse(line_printer& p, Integer block_width) const
{
    (void)block_width;
    p.disp_empty_line();
};
void disp_stream_global::end_display_matrix_block_sparse(line_printer& p, Integer block_width) const
{
    (void)block_width;
    (void)p;
};

Integer disp_stream_global::get_max_cols() const
{
    return options().get_option<Integer>(matcl::opt::disp::max_cols());
};
Integer disp_stream_global::get_max_rows() const
{
    return options().get_option<Integer>(matcl::opt::disp::max_rows());
};
Integer disp_stream_global::get_max_nnz() const
{
    return options().get_option<Integer>(matcl::opt::disp::max_nnz());
};
bool disp_stream_global::restrict_sparse_matrix_size() const
{
    return options().get_option<bool>(matcl::opt::disp::restrict_sparse_matrix_size());
};

Integer disp_stream_global::get_precision() const
{
    return options().get_option<Integer>(matcl::opt::disp::precision());
};

bool disp_stream_global::separate_columns_sparse() const
{
    return true;
};
bool disp_stream_global::disp_header_only() const
{
    return options().get_option<bool>(matcl::opt::disp::header_only());
};

std::string disp_stream_global::get_col_name(Integer c) const
{
    return boost::lexical_cast<std::string>(c+1);
};
std::string disp_stream_global::get_row_name(Integer r) const
{
    return boost::lexical_cast<std::string>(r+1);
};
std::string disp_stream_global::get_rows_label() const
{
    return "";
};
disp_stream_ptr matcl::default_disp_stream()
{
    disp_stream_ptr ds = disp_stream_ptr(new disp_stream_global(global_output_stream()));
    return ds;
};

//------------------------------------------------------------------------
//                          forwarding_disp_stream
//------------------------------------------------------------------------
static output_stream_ptr get_output_stream(const disp_stream_ptr& other_stream)
{
    if (!other_stream)
        throw std::runtime_error("uninitialized output_stream_ptr used");

    return other_stream->impl()->get_output_stream();
};
forwarding_disp_stream::forwarding_disp_stream(const disp_stream_ptr& other_stream, const output_stream_ptr& os)
    :disp_stream(os), m_impl(other_stream)
{
    if (!other_stream)
        throw std::runtime_error("uninitialized output_stream_ptr used");
    if (!other_stream)
        throw std::runtime_error("uninitialized disp_stream_ptr used");
};
forwarding_disp_stream::forwarding_disp_stream(const disp_stream_ptr& other_stream, std::ostream& os)
    :disp_stream(os), m_impl(other_stream)
{
    if (!other_stream)
        throw std::runtime_error("uninitialized output_stream_ptr used");
    if (!other_stream)
        throw std::runtime_error("uninitialized disp_stream_ptr used");
};
forwarding_disp_stream::forwarding_disp_stream(const disp_stream_ptr& other_stream)
    :disp_stream(get_output_stream(other_stream)), m_impl(other_stream)
{
    if (!other_stream)
        throw std::runtime_error("uninitialized disp_stream_ptr used");
}

forwarding_disp_stream::~forwarding_disp_stream()
{};

std::string forwarding_disp_stream::get_value_name(matcl::value_code vt) const
{
    return m_impl->get_value_name(vt);
};

bool forwarding_disp_stream::show_matrix_header() const
{
    return m_impl->show_matrix_header();
};

void forwarding_disp_stream::general_info_string(line_printer& p, bool header_only) const
{
    return m_impl->general_info_string(p,header_only);
};

void forwarding_disp_stream::general_info_scalar(line_printer& p, matcl::value_code vt, bool as_matrix,
                                                 const std::string& type_name) const
{
    return m_impl->general_info_scalar(p, vt, as_matrix, type_name);
};

void forwarding_disp_stream::general_info_dense_matrix(line_printer& p, Integer r, Integer c, matcl::value_code vt,
                                          const std::string& struct_name) const
{
    return m_impl->general_info_dense_matrix(p, r,c,vt, struct_name);
};
void forwarding_disp_stream::general_info_banded_matrix(line_printer& p, Integer r, Integer c, Integer fd, Integer ld, 
                                           matcl::value_code vt, const std::string& struct_name) const
{
    return m_impl->general_info_banded_matrix(p, r,c,fd,ld,vt, struct_name);
};
void forwarding_disp_stream::general_info_sparse_matrix(line_printer& p, Integer r, Integer c, Integer nz, matcl::value_code vt,
                                           const std::string& struct_name) const
{
    return m_impl->general_info_sparse_matrix(p, r,c,nz,vt, struct_name);
};

void forwarding_disp_stream::display_empty_matrix(line_printer& p, Integer r, Integer c) const
{
    return m_impl->display_empty_matrix(p, r,c);
};

bool forwarding_disp_stream::show_final_continuation() const		
{ 
    return m_impl->show_final_continuation();
};
std::string forwarding_disp_stream::get_final_continuation() const
{ 
    return m_impl->get_final_continuation();
};
Integer forwarding_disp_stream::get_number_subvalues() const
{
    return m_impl->get_number_subvalues();
};
align_type forwarding_disp_stream::get_align_row_header() const
{
    return m_impl->get_align_row_header();
};
align_type forwarding_disp_stream::get_align_col(Integer c) const
{
    return m_impl->get_align_col(c);
};

Integer	forwarding_disp_stream::get_terminal_width() const
{
    return m_impl->get_terminal_width();
};

std::string forwarding_disp_stream::get_col_separator() const
{
    return m_impl->get_col_separator();
};
std::string forwarding_disp_stream::get_col_values_separator() const
{
    return m_impl->get_col_values_separator();
};

Integer	forwarding_disp_stream::get_max_value_width(Integer value_pos) const
{
    return m_impl->get_max_value_width(value_pos);
};

Integer forwarding_disp_stream::get_sparse_separator_width() const
{
    return m_impl->get_sparse_separator_width();
};
std::string forwarding_disp_stream::get_sparse_separator() const
{
    return m_impl->get_sparse_separator();
};
bool forwarding_disp_stream::display_zero() const
{
    return m_impl->display_zero();
};
bool forwarding_disp_stream::ignore_lower_triangle() const
{
    return m_impl->ignore_lower_triangle();
};
disp_mode forwarding_disp_stream::get_disp_mode() const
{
    return m_impl->get_disp_mode();
};

void forwarding_disp_stream::start_display(line_printer& p) const
{
    return m_impl->start_display(p);
};
void forwarding_disp_stream::end_display(line_printer& p) const
{
    return m_impl->end_display(p);
}

std::string forwarding_disp_stream::end_line(Integer line) const
{
    return m_impl->end_line(line);
};
std::string forwarding_disp_stream::new_line(Integer line) const
{
    return m_impl->new_line(line);
};
Integer forwarding_disp_stream::new_line_width(Integer terminal_width) const
{
    return m_impl->new_line_width(terminal_width);
};
Integer forwarding_disp_stream::end_line_width(Integer terminal_width) const
{
    return m_impl->end_line_width(terminal_width);
};
std::string forwarding_disp_stream::new_line_special(line_type type) const
{
    return m_impl->new_line_special(type);
};
std::string forwarding_disp_stream::end_line_special(line_type type) const
{
    return m_impl->end_line_special(type);
};
Integer forwarding_disp_stream::new_line_special_width(Integer terminal_width, line_type type) const
{
    return m_impl->new_line_special_width(terminal_width, type);
};
Integer forwarding_disp_stream::end_line_special_width(Integer terminal_width, line_type type) const
{
    return m_impl->end_line_special_width(terminal_width, type);
};

forwarding_disp_stream::string forwarding_disp_stream::get_labels_row_id() const
{
    return m_impl->get_labels_row_id();
};
bool forwarding_disp_stream::show_row_headers() const
{
    return m_impl->show_row_headers();
};
std::string forwarding_disp_stream::get_first_col_separator() const
{
    return m_impl->get_first_col_separator();
};
bool forwarding_disp_stream::show_column_header() const
{
    return m_impl->show_column_header();
};
bool forwarding_disp_stream::show_first_row_separator() const
{
    return m_impl->show_first_row_separator();
};
std::string forwarding_disp_stream::get_subvalue_label(Integer value_pos) const
{
    return m_impl->get_subvalue_label(value_pos);
};
bool forwarding_disp_stream::show_values_labels() const
{
    return m_impl->show_values_labels();
};

void forwarding_disp_stream::display_matrix_name(line_printer& p) const
{
    return m_impl->display_matrix_name(p);
};
void forwarding_disp_stream::start_display_matrix_block(line_printer& p, Integer block_width, Integer first_col,
                            Integer last_col) const
{
    return m_impl->start_display_matrix_block(p, block_width, first_col, last_col);
}
void forwarding_disp_stream::end_display_matrix_block(line_printer& p, Integer block_width) const
{
    return m_impl->end_display_matrix_block(p, block_width);
};
void forwarding_disp_stream::start_display_matrix_block_sparse(line_printer& p, Integer block_width) const
{
    return m_impl->start_display_matrix_block_sparse(p, block_width);
};
void forwarding_disp_stream::end_display_matrix_block_sparse(line_printer& p, Integer block_width) const
{
    return m_impl->end_display_matrix_block_sparse(p, block_width);
};



Integer forwarding_disp_stream::get_max_cols() const
{
    return m_impl->get_max_cols();
};
Integer forwarding_disp_stream::get_max_rows() const
{
    return m_impl->get_max_rows();
};
Integer forwarding_disp_stream::get_max_nnz() const
{
    return m_impl->get_max_nnz();
}
bool forwarding_disp_stream::restrict_sparse_matrix_size() const
{
    return m_impl->restrict_sparse_matrix_size();
};
Integer forwarding_disp_stream::get_precision() const
{
    return m_impl->get_precision();
};
const disp_stream* forwarding_disp_stream::get_orginal_stream() const
{
    return m_impl.get();
};
bool forwarding_disp_stream::separate_columns_sparse() const
{
    return m_impl->separate_columns_sparse();
};
bool forwarding_disp_stream::disp_header_only() const
{
    return m_impl->disp_header_only();
};
Integer forwarding_disp_stream::get_row_block() const
{
    return m_impl->get_row_block();
};
std::string forwarding_disp_stream::get_col_name(Integer c) const
{
    return m_impl->get_col_name(c);
};
std::string forwarding_disp_stream::get_row_name(Integer r) const
{
    return m_impl->get_row_name(r);
};
std::string forwarding_disp_stream::get_rows_label() const
{
    return m_impl->get_rows_label();
};

disp_stream_default::disp_stream_default(const output_stream_ptr& os)
    :forwarding_disp_stream(default_disp_stream(), os)
{};
disp_stream_default::disp_stream_default(std::ostream& os)
    :forwarding_disp_stream(default_disp_stream(), os)
{};

disp_stream_default::disp_stream_default()
    :forwarding_disp_stream(default_disp_stream())
{};

disp_stream_default::~disp_stream_default()
{};

//-----------------------------------------------------------------------------------------
//                          disp_stream_labels
//-----------------------------------------------------------------------------------------
disp_labels_provider::disp_labels_provider()
    :m_given_by_functor(false)
{};
disp_labels_provider::~disp_labels_provider()
{};

disp_labels_provider::disp_labels_provider(const std::vector<std::string>& labels)
    :m_labels(labels), m_given_by_functor(false)
{};
disp_labels_provider::disp_labels_provider(std::vector<std::string>&& labels)
    :m_labels(std::move(labels)), m_given_by_functor(false)
{};

disp_labels_provider::disp_labels_provider(const functor_type& labels_functor)
    :m_functor(labels_functor), m_given_by_functor(true)
{};

disp_labels_provider::disp_labels_provider(std::initializer_list<std::string> labels)
    :m_given_by_functor(false)
{
    m_labels.insert(m_labels.begin(), labels.begin(), labels.end());
};

std::string disp_labels_provider::get(Integer pos) const
{
    if (m_given_by_functor == true)
        return m_functor(pos);
            
    if (pos < (Integer)m_labels.size())
        return m_labels[pos];
    else
        return boost::lexical_cast<std::string>(pos+1);
};

disp_stream_label::disp_stream_label(const disp_stream_ptr& other_stream, const disp_labels_provider& rows, 
                    const disp_labels_provider& cols, const std::string& row_label)
    :forwarding_disp_stream(other_stream), m_labels_rows(rows), m_labels_cols(cols), m_row_label(row_label)
{};

disp_stream_label::disp_stream_label(const disp_labels_provider& rows, const disp_labels_provider& cols)
    :forwarding_disp_stream(default_disp_stream()), m_labels_rows(rows), m_labels_cols(cols)
{};

matcl::disp_mode disp_stream_label::get_disp_mode() const
{
    return disp_mode::scalar_dense;
};

std::string disp_stream_label::get_row_name(Integer r) const
{
    return m_labels_rows.get(r);
};
std::string disp_stream_label::get_rows_label() const
{
    return m_row_label;
};
std::string disp_stream_label::get_col_name(Integer c) const
{
    return m_labels_cols.get(c);
};

}