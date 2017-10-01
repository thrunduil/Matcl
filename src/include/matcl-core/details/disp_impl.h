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

#include <iostream>
#include "matcl-core/IO/disp_stream.h"
#include "matcl-core/IO/disp_data_provider.h"

namespace matcl { namespace raw { namespace details
{

namespace md = matcl :: details;

template<class V>
struct matrix_provider_base : matcl::details::label_iterators
{
    virtual ~matrix_provider_base(){};

    virtual Integer     rows() const = 0;
    virtual Integer     cols() const = 0;

    // iterator like access. It is assumed, that the only possible
    // operations are return to saved point, go to next row and go to
    // next column. We require two saved points
    virtual void        begin(){};          //point to first element of the matrix
    virtual void        hold(){};           //store current point in variable V1
    virtual void        hold_column(){};    //store current point in variable V2
    virtual void        next_row(){};       //go to next row in given column
    virtual void        next_column(){};    //go to next column in given row
    virtual void        restore(){};        //go to the point stored in variable V1
    virtual void        restore_column(){}; //go to the point stored in variable V2

    virtual void        begin_row_headers() override{};
    virtual void        next_row_header() override{};
    virtual void        begin_col_headers() override{};
    virtual void        next_col_header() override{};

    // get value at given point. If given data structure supports random access iterations
    // then required value can be accessed at row r and column c
    // if is_zero is set to true; then value is treated as zero and need not be printed
    virtual V           get_value(const disp_stream* user, Integer width, align_type at,
                            Integer r, Integer c, bool& is_zero) const = 0 ;

    virtual bool        is_symher() const = 0;
};

// this class can be specialized for T = Object
template<class T>
struct elem_handle_type
{
    using type = T;
};

// make call to printer::disp_elem for a printer stored in 
// given disp stream
template<class V, class S>
struct disp_elem_helper
{
    using handle_type = typename elem_handle_type<V>::type;

    static void eval(md::disp_stream_impl& os, const disp_stream* user, 
                    const V& elem, const std::vector<Integer>& widths, 
                     Integer col);
    
    static void eval(md::disp_stream_impl& os, const disp_stream* user, 
                    const std::string& sym, const V& elem, bool is_zero,
                    const std::vector<Integer>& widths, Integer col);
};

//--------------------------------------------------------------------------
//                         impl disp matrix
//--------------------------------------------------------------------------

template<class V, class S>
class disp_matrix
{};

template<class V>
class disp_matrix<V,struct_dense>
{   
    private:
        using handle_type = typename elem_handle_type<V>::type;

        struct measures
        {
            Integer chw;    //column header width;
            Integer csw;    //col separator width;
            Integer rhw;    //row headers width;
            Integer fcw;    //first col separator width;
        };

    public:
        static void eval_matrix_body(md::disp_stream_impl& os, const disp_stream* user, 
                                    matrix_provider_base<V>& m);

    private:
        static void measure_column_header(md::disp_stream_impl& os, const disp_stream* user, 
                        measures& ms, Integer lc, Integer cc, bool need_cont, bool add_cont);

        static void make_column_header(md::disp_stream_impl& os, const disp_stream* user, 
                                    measures& ms, Integer lc, Integer cc, bool need_cont, 
                                    bool add_cont);

        static void make_row_separator(md::disp_stream_impl& os, const disp_stream* user, 
                                    measures& ms);

        static void make_row(md::disp_stream_impl& os, const disp_stream* user, measures& ms, 
                            Integer i, Integer lc, Integer cc, matrix_provider_base<V>& mat,
                            bool value, bool need_cont, bool add_cont, bool is_sym);
};

}}}

namespace matcl { namespace raw
{

// configure matrix display (alignment, row and column labels, ect.)
class MATCL_CORE_EXPORT disp_stream_data_provider : public forwarding_disp_stream
{
    private:
        const disp_data_provider&   m_dp;
        const disp_stream*          m_owner;

    public:
        disp_stream_data_provider(const disp_stream_ptr& impl, const disp_data_provider &dp);
            
        // part of disp_stream interface
        virtual void start_display(line_printer& p) const override;

        virtual void end_display(line_printer& p) const override;
        
        virtual void display_empty_matrix(line_printer& p, Integer r,
                                Integer c) const override;

        virtual void display_matrix_name(line_printer& p) const;

        // a new block of dense matrix will be displayed 
        // for columns [first_col, last_col].
        virtual void start_display_matrix_block(line_printer& p, Integer block_width, 
                                Integer first_col, Integer last_col) const override;

        // finalize current block for dense matrix
        virtual void end_display_matrix_block(line_printer& p, Integer 
                        block_width) const override;

        virtual bool        show_column_header() const override;
        virtual bool        show_row_headers()  const override;
        virtual align_type  get_align_row_header() const override;
        virtual align_type  get_align_col(Integer c) const;
        virtual string      get_row_name(Integer r) const override;
        virtual string      get_col_name(Integer c) const override;
        virtual string      get_rows_label() const override;
};

MATCL_CORE_EXPORT void  disp(const disp_stream_ptr &os, disp_data_provider& v);

}};

#include "matcl-core/details/disp_impl.inl"