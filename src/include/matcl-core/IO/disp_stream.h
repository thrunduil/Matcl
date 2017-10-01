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
#include "matcl-core/matrix/enums.h"
#include "matcl-core/IO/output_stream.h"

#include <string>
#include <iosfwd>
#include <vector>
#include <memory>

namespace matcl 
{

#pragma warning(push)
#pragma warning(disable:4251)   //warning C4251: needs to have dll-interface

class MATCL_CORE_EXPORT line_printer
{
    public:
        virtual ~line_printer(){};

        // print empty line without prefix and sufix
        virtual void    disp_empty_line() = 0;

        // print string without prefix and sufix
        virtual void    disp_raw_string(const std::string& str, align_type at) = 0;

        // print string with prefix and sufix
        virtual void    disp_string(const std::string& str, align_type at) = 0;
};
//-----------------------------------------------------------------------------------------
//                          abstract disp_stream
//-----------------------------------------------------------------------------------------

using disp_stream_ptr = std::shared_ptr<class disp_stream>;

// this class controls behavior of matrix printers
class MATCL_CORE_EXPORT disp_stream
{
    public:
        enum class line_type
        {
            // sequence of events
            start_display,      // begin preambule
            matrix_name,        // display matrix name
            matrix_type,        // display matrix type
            end_display,        // final line

            // special events
            empty_matrix,       // displaying empty matrix, after matrix_type
            empty_line,         // displaying empty line
        };

    protected:
        using impl_type = std::shared_ptr<details::disp_stream_impl>;
        using string    = std::string;

    private:
        impl_type           m_impl;

    public:
        explicit disp_stream(const output_stream_ptr& os);
        explicit disp_stream(std::ostream& os);

        virtual ~disp_stream();

    public:
        //-------------------------------------------------------------------
        //              initialization, finalization
        //-------------------------------------------------------------------
        virtual void        start_display(line_printer& p) const = 0;
        virtual void        end_display(line_printer& p) const = 0;

        //-------------------------------------------------------------------
        //                          preambule
        //-------------------------------------------------------------------
        virtual void        general_info_string(line_printer& p, bool header_only) const = 0;
        virtual void        general_info_scalar(line_printer& p, matcl::value_code vt, bool as_matrix,
                                    const std::string& type_name) const = 0;
        virtual void        general_info_dense_matrix(line_printer& p, Integer r, Integer c, matcl::value_code vt,
                                    const std::string& struct_name) const = 0;
        virtual void        general_info_banded_matrix(line_printer& p, Integer r, Integer c, Integer fd,
                                    Integer ld, matcl::value_code vt, const std::string& struct_name) const = 0;
        virtual void        general_info_sparse_matrix(line_printer& p, Integer r, Integer c, Integer nz,
                                    matcl::value_code vt, const std::string& struct_name) const = 0;
        virtual void        display_empty_matrix(line_printer& p, Integer r, Integer c) const = 0;     
        virtual void        display_matrix_name(line_printer& p) const = 0;

        // a new block of dense matrix will be displayed for columns [first_col, last_col].
        virtual void        start_display_matrix_block(line_printer& p, Integer block_width, Integer first_col,
                                    Integer last_col) const = 0;

        // finalize current block for dense matrix
        virtual void        end_display_matrix_block(line_printer& p, Integer block_width) const  = 0;
        
        // a new block of sparse matrix will be displayed
        virtual void        start_display_matrix_block_sparse(line_printer& p, Integer block_width) const  = 0;
        
        // finalize current block for sparse matrix
        virtual void        end_display_matrix_block_sparse(line_printer& p, Integer block_width) const  = 0;

        //-------------------------------------------------------------------
        //                  printing options
        //-------------------------------------------------------------------
        virtual bool        show_matrix_header()  const = 0;
        virtual bool        show_column_header() const = 0;
        virtual bool        show_row_headers()  const = 0;
        virtual bool        show_values_labels() const = 0;
        virtual bool        show_first_row_separator() const = 0;        
        virtual bool        show_final_continuation() const = 0;
        virtual bool        separate_columns_sparse() const = 0;        
        virtual string      get_col_separator() const = 0;        
        virtual string      get_first_col_separator() const = 0;                
        virtual string      get_col_values_separator() const = 0;
        virtual string      get_sparse_separator() const = 0;
        virtual Integer     get_sparse_separator_width() const = 0;
        virtual string      get_labels_row_id() const = 0;        
        virtual std::string get_final_continuation() const = 0;
        virtual align_type  get_align_row_header() const = 0;
        virtual align_type  get_align_col(Integer c) const = 0;

        //-------------------------------------------------------------------
        //                  line prefix and sufix
        //-------------------------------------------------------------------
        // string at the beginning of new row of data, matrix_line is row number
        // for dense matrices and element number for sparse matrices (1-based),
        // matrix_line = 0 mean, that special matrix line is printed (i.e. header,
        // separator). New line character is automatically added. Strings are 
        // automatically adjusted to new line or end line width.
        virtual std::string new_line(Integer matrix_line) const = 0;

        //string at the end of new row of data
        virtual std::string end_line(Integer matrix_line) const = 0;

        // total size of new line and end line string cannot exceed 50% 
        // of terminal width
        virtual Integer     new_line_width(Integer terminal_width) const = 0;
        virtual Integer     end_line_width(Integer terminal_width) const = 0;

        // new lines and end lines in for preambule, header, finalized, etc.
        virtual std::string new_line_special(line_type type) const = 0;
        virtual std::string end_line_special(line_type type) const = 0;

        // total size of new line and end line string cannot exceed 50% 
        // of terminal width
        virtual Integer     new_line_special_width(Integer terminal_width, line_type type) const = 0;
        virtual Integer     end_line_special_width(Integer terminal_width, line_type type) const = 0;

        //-------------------------------------------------------------------
        //                  object value properties
        //-------------------------------------------------------------------
        // one object may contain many independent values

        // number of subvalues for one object
        virtual Integer     get_number_subvalues() const = 0;

        // maximum size for given subvalue, (1-based), -1 if unknown
        virtual Integer     get_max_value_width(Integer subvalue_index) const = 0;

        // label for given subvalue
        virtual std::string get_subvalue_label(Integer subvalue_index) const = 0;

        // row and col labels, indices are 0-based
        virtual string      get_row_name(Integer r) const = 0;
        virtual string      get_col_name(Integer c) const = 0;

        // Label of the row names column
        virtual string      get_rows_label() const = 0;

        //-------------------------------------------------------------------
        //                 general matrix properties
        //-------------------------------------------------------------------
        virtual std::string get_value_name(matcl::value_code vt) const = 0;        

        //-------------------------------------------------------------------
        //                 general printing properties
        //-------------------------------------------------------------------
        virtual Integer     get_terminal_width() const = 0;
        
        //maximum number of columns printed
        virtual Integer     get_max_cols() const = 0;
        
        //maximum number of rows printed
        virtual Integer     get_max_rows() const = 0;
        
        //maximum number of nonzero elements printed for sparse matrices
        virtual Integer     get_max_nnz() const = 0;

        //add limit on number of rows and columns printed for sparse matrices
        virtual bool        restrict_sparse_matrix_size() const = 0;

        //do display zero elements?
        virtual bool        display_zero() const = 0;

        //display lower traingle for symmetric and hermitian matrices?
        virtual bool        ignore_lower_triangle() const = 0;

        //number of significant digits printed
        virtual Integer     get_precision() const = 0;

        //how to display scalars and sparse matrices?
        virtual disp_mode   get_disp_mode() const = 0;

        //if true then only header will be displayed
        virtual bool        disp_header_only() const = 0;

        //distance between row separators, if <= 0 then row separators
        //are not displayed
        virtual Integer     get_row_block() const = 0;

    //internal use
    public:
        impl_type           impl() const;
};

//-----------------------------------------------------------------------------------------
//                          forwarding_disp_stream
//-----------------------------------------------------------------------------------------
// this class implements all disp_stream interface by forwarding calls
// to other disp_stream. In derived clas one can redefine some of functions.
// Other disp streams should inherit from this class.
class MATCL_CORE_EXPORT forwarding_disp_stream : public disp_stream
{
    private:
        disp_stream_ptr m_impl;

    public:
        // take implementation from other_stream, but use own output stream
        forwarding_disp_stream(const disp_stream_ptr& other_stream, const output_stream_ptr& os);
        forwarding_disp_stream(const disp_stream_ptr& other_stream, std::ostream& os);

        // take implementation and output stream from other_stream
        explicit forwarding_disp_stream(const disp_stream_ptr& other_stream);

        virtual ~forwarding_disp_stream();

    public:
        virtual void        start_display(line_printer& p) const override;
        virtual void        end_display(line_printer& p) const override;
        virtual void        general_info_string(line_printer& p, bool header_only) const override;
        virtual void        general_info_scalar(line_printer& p, matcl::value_code vt, bool as_matrix,
                                    const std::string& type_name) const override;
        virtual void        general_info_dense_matrix(line_printer& p, Integer r, Integer c, matcl::value_code vt, 
                                    const std::string& struct_name) const override;
        virtual void        general_info_banded_matrix(line_printer& p, Integer r, Integer c, 
                                    Integer fd, Integer ld, matcl::value_code vt, 
                                    const std::string& struct_name) const override;
        virtual void        general_info_sparse_matrix(line_printer& p, Integer r, Integer c, Integer nz,
                                    matcl::value_code vt, const std::string& struct_name) const override;
        virtual void        display_empty_matrix(line_printer& p, Integer r, Integer c) const override;
        virtual bool        show_matrix_header()  const override;
        virtual bool        show_column_header() const override;        
        virtual string      get_col_separator() const override;
        virtual bool        show_row_headers()  const override;
        virtual string      get_first_col_separator() const override;                
        virtual bool        show_values_labels() const override;
        virtual bool        show_first_row_separator() const override;
        virtual string      get_col_values_separator() const override;
        virtual bool        separate_columns_sparse() const override;
        virtual bool        disp_header_only() const override;
        virtual Integer     get_row_block() const override;
        virtual string      get_sparse_separator() const override;
        virtual Integer     get_sparse_separator_width() const override;
        virtual string      get_labels_row_id() const override;
        virtual bool        show_final_continuation() const override;
        virtual Integer     get_terminal_width() const override;
        virtual Integer     get_max_cols() const override;
        virtual Integer     get_max_rows() const override;
        virtual Integer     get_max_nnz() const override;
        virtual bool        restrict_sparse_matrix_size() const override;
        virtual bool        display_zero() const override;
        virtual bool        ignore_lower_triangle() const override;
        virtual Integer     get_precision() const override;
        virtual disp_mode   get_disp_mode() const override;
        virtual std::string new_line(Integer matrix_line) const override;
        virtual std::string end_line(Integer matrix_line) const override;
        virtual Integer     new_line_width(Integer terminal_width) const override;
        virtual Integer     end_line_width(Integer terminal_width) const override;
        virtual std::string new_line_special(line_type type) const override;
        virtual std::string end_line_special(line_type type) const override;
        virtual Integer     new_line_special_width(Integer terminal_width, line_type type) const override;
        virtual Integer     end_line_special_width(Integer terminal_width, line_type type) const override;
        virtual std::string get_final_continuation() const override;
        virtual std::string get_value_name(matcl::value_code vt) const override;
        virtual void        display_matrix_name(line_printer& p)  const override;
        virtual void        start_display_matrix_block(line_printer& p, Integer block_width, Integer first_col,
                                    Integer last_col) const  override;
        virtual void        end_display_matrix_block(line_printer& p, Integer block_width) const  override;
        virtual void        start_display_matrix_block_sparse(line_printer& p, Integer block_width) const  override;
                            // finalize current block for sparse matrix
        virtual void        end_display_matrix_block_sparse(line_printer& p, Integer block_width) const  override;
        virtual std::string get_subvalue_label(Integer value_pos) const override;
        virtual Integer     get_number_subvalues() const override;
        virtual Integer     get_max_value_width(Integer value_pos) const override;
        virtual string      get_row_name(Integer r) const override;
        virtual string      get_col_name(Integer c) const override;
        virtual string      get_rows_label() const override;
        virtual align_type  get_align_row_header() const override;
        virtual align_type  get_align_col(Integer c) const override;

    public:
        // Get orginal disp stream
        const disp_stream*  get_orginal_stream() const;
};

//-----------------------------------------------------------------------------------------
//                          disp_stream_default
//-----------------------------------------------------------------------------------------
// this stream takes all options from default_disp_stream but allows to
// change output stream; this class is not thread safe and instances of this
// and derived classes should not be shared between threads
class MATCL_CORE_EXPORT disp_stream_default : public forwarding_disp_stream
{
    public:
        // use global output_stream
        disp_stream_default();

        explicit disp_stream_default(const output_stream_ptr& os);
        explicit disp_stream_default(std::ostream& os);                 

        virtual ~disp_stream_default();
};

// default disp stream, which contains default printing options and prints on
// global output stream. Global printing options can be changed using matcl_options
// class, and changing default values of options. See options_disp and 
// matcl_options for details. Default disp stream is thread safe and should not be
// shared between threads
MATCL_CORE_EXPORT disp_stream_ptr    default_disp_stream();

//-----------------------------------------------------------------------------------------
//                          disp_stream_labels
//-----------------------------------------------------------------------------------------
// this class is storing labels
class MATCL_CORE_EXPORT disp_labels_provider
{
    private:
        using functor_type  = std::function<std::string (Integer pos)>;
        using vector_type   = std::vector<std::string>;

    private:
        vector_type         m_labels;
        functor_type        m_functor;
        bool                m_given_by_functor;

    public:
        // no labels
        disp_labels_provider();
        ~disp_labels_provider();

        //labels size can be lower or higher than matrix size
        disp_labels_provider(const std::vector<std::string>& labels);
        disp_labels_provider(std::vector<std::string>&& labels);

        // labels are provided by functor labels_functor
        disp_labels_provider(const functor_type& labels_functor);

        // labels constructed from initializer list
        disp_labels_provider(std::initializer_list<std::string> labels);

    public:
        //get label for row or column with index pos (0-based)
        std::string         get(Integer pos) const;
};

// this disp_stream allows for named row and column labels
// example: disp_stream_label({'ROW 1', 'ROW 2'}, {'COL 1', 'COL 2});
class MATCL_CORE_EXPORT disp_stream_label : public matcl::forwarding_disp_stream
{
    private:
        disp_labels_provider    m_labels_rows;
        disp_labels_provider    m_labels_cols;
        std::string             m_row_label;

    public:
        // use disp options as well as output stream from other stream
        disp_stream_label(const disp_stream_ptr& other_stream, const disp_labels_provider& rows, 
                          const disp_labels_provider& cols, const std::string& row_label = "");

        // use disp options as well as output stream from default_disp_stream()
        disp_stream_label(const disp_labels_provider& rows, const disp_labels_provider& cols);

        virtual matcl::disp_mode
                            get_disp_mode() const override;

        virtual string      get_row_name(Integer r) const override;
        virtual string      get_col_name(Integer c) const override;
        virtual string      get_rows_label() const override;
};

#pragma warning(pop)
}