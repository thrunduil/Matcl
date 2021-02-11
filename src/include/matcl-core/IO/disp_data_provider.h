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

#include "matcl-core/IO/disp_stream.h"

namespace matcl
{

// Class allows for displaying any data using formatting rules used for matrices.
// Actual display is performed by disp function.
class MATCL_CORE_EXPORT disp_data_provider : public matcl_new_delete
{
    public:
        virtual ~disp_data_provider(){};

        // Number of rows printed
        virtual Integer     rows() const = 0;

        // Number of columns printed
        virtual Integer     cols() const = 0;

        // Get element labeled with row index r and columnt index c as string.
        // This is also value of element curently pointed by iterator.
        // Row and column indices are 0-based. If width == 0, then we are in
        // measurement phase. If width > 0, then it is maximum allowed width
        // of given value. If returned string is longer, then will be truncated.
        virtual std::string get_value(const disp_stream* ds, Integer width,     
                                      align_type at, Integer r, Integer c) const = 0;

        // If given data structure does not allow for random access, then given point
        // through iterator like access. It is assumed, that data structure supports
        // forward iterator like access in rows and columns, i.e. it is assumed, that 
        // a set of possible operations contains: return to saved point, go to next row
        // and go to next column. One requires two variables that store saved points.
        virtual void        begin(){};          //point to first element of the matrix
        virtual void        hold(){};           //store current point in variable V1
        virtual void        hold_column(){};    //store current point in variable V2
                                                //V2 always points to beginning of block
                                                //currently displayed
        virtual void        next_row(){};       //go to next row in given column
        virtual void        next_column(){};    //go to next column in given row
        virtual void        restore(){};        //go to the point stored in variable V1
        virtual void        restore_column(){}; //go to the point stored in variable V2

        // If given data structure does not allow for random access iteration or forward
        // iteration, then it is rather impossible to use this class. Most of standard
        // containers satisfy this requirements. On the other hand sparse matrices do not.

        // Similar iterators, that allows for accesssing row name and column name
        // through iterator, not position number
        virtual void        begin_row_headers(){};
        virtual void        next_row_header(){};
        virtual void        begin_col_headers(){};
        virtual void        next_col_header(){};

        // Get general information about data being displayed (e.g. matrix type 
        // and size in case of matrices). If empty string is returned (default
        // implementation), then general description is not printed.
        virtual std::string get_data_type_description() const;

        // --------------------------------------------------------------------
        //                 disp stream interface
        // --------------------------------------------------------------------
        // This is part of disp stream interface, which should be modified for
        // given data set. See disp_stream class for additional information.
        // All functions accepts argument orig, which is original disp_stream
        // passed to disp function. Default implementation of the following function
        // simply forwards calls to original disp stream.

        // Data display begins
        virtual void        start_display(const disp_stream* orig, line_printer& p) const;
        // Data display ends
        virtual void        end_display(const disp_stream* orig, line_printer& p) const;
        // How to print empty data?
        virtual void        display_empty_matrix(const disp_stream* orig,line_printer& p, Integer r,
                                    Integer c) const;
        // call display_empty_matrix function for empty matrices?
        virtual bool        short_print_empty_matrix(const disp_stream* orig) const;
        // Name of data set
        virtual void        display_matrix_name(const disp_stream* orig, line_printer& p) const;
        // show matrix header?
        virtual bool        show_matrix_header(const disp_stream* orig) const;

        // A new block of dense matrix will be displayed for columns [first_col, last_col].
        virtual void        start_display_matrix_block(const disp_stream* orig, line_printer& p, 
                                    Integer block_width, Integer first_col, Integer last_col) const;

        // Finalize current block of data
        virtual void        end_display_matrix_block(const disp_stream* orig, line_printer& p, 
                                    Integer block_width) const;

        // show row and column destriptions?
        virtual bool        show_column_header_line(const disp_stream* orig) const;
        
        // show row destription? called when show_column_header_line() = true
        virtual bool        show_column_header_row(const disp_stream* orig)  const;
        
        // show columns destription? called when show_column_header_line() = true
        virtual bool        show_column_header_columns(const disp_stream* orig) const;

        // if true, then matrix can be split on blocks, when required
        virtual bool        can_split(const disp_stream* orig) const;

        // How to align row headers?
        virtual align_type  get_align_row_header(const disp_stream* orig) const;
        // How to align columns?
        virtual align_type  get_align_col(const disp_stream* orig, Integer c) const;

        //Get name of given row, index 0-based
        virtual std::string get_row_name(const disp_stream* orig, Integer r) const;
        //Get name of given column, index 0-based
        virtual std::string get_col_name(const disp_stream* orig, Integer c) const;
        // Label of the row names column
        virtual std::string get_rows_label(const disp_stream* orig) const;

        // get minimum and maximum width of column with row labels
        virtual void        get_column_width_row(const disp_stream* orig, Integer& w_min, 
                                Integer& w_max) const;

        // get minimum and maximum width of i-th column
        virtual void        get_column_width(const disp_stream* orig, Integer c, 
                                Integer& w_min, Integer& w_max) const;

    public:
        // Converting functions from Integer, Real, and Complex to string.
        // It is not strictly necessary to use these functions, but this is the 
        // only way to print Real and Complex number in the same way as for matrices
        // if width > 0, then length of resulting string is exactly width; if width < 0
        // then width is the maximum length of resulting string; if width = 0, then
        // length of resulting string is value dependent
        std::string         to_string(const disp_stream* user, Integer width, Integer val,
                                      align_type at) const;
        std::string         to_string(const disp_stream* user, Integer width, Real val,
                                      align_type at) const;
        std::string         to_string(const disp_stream* user, Integer width, const Complex& val,
                                      align_type at) const;        
        std::string         to_string(const disp_stream* user, Integer width, const dynamic::object& val,
                                      align_type at) const;        
        std::string         to_string(const disp_stream* user, Integer width, const std::string& val,
                                      align_type at) const;        
};

};