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

#include "matcl-scalar/config.h"
#include "matcl-scalar/general/fwd_decls.h"
#include "matcl-core/IO/disp_stream.h"
#include "matcl-core/options/matcl_options.h"

#include <string>

#pragma warning(push)
#pragma warning(disable : 4251)  //needs to have dll-interface to be used by clients

namespace matcl
{

// format of a column
class MATCL_SCALAR_EXPORT column_info
{
    private:
        std::string         m_name;
        align_type          m_align;
        int                 m_width;

    public:
        // create uninitialized object
        column_info();

        // define column name, alignment type and preferred width 
        // in given column
        column_info(const std::string& name, align_type align, int width);

    public:
        // return column name
        const std::string&  name() const;

        // return alignment type
        align_type          align() const;

        // return preferred column width
        int                 width() const;
};

// class preforming formatted display of a vector of data
// according to specified format
class MATCL_SCALAR_EXPORT formatted_disp
{
    private:
        using impl_type     = std::shared_ptr<details::formatted_disp_stream>;

    private:
        impl_type           m_impl;

    public:
        // intialize with given disp stream and printing options
        formatted_disp(const disp_stream_ptr& os = default_disp_stream(),
                        const options& opts = options());

    public:
        // set format of row labels; if not set, then row labels are not printed
        void                set_row_label(const std::string& name, align_type al,
                                int width);
        void                set_row_label(const column_info& format);

        // set format of next column
        void                add_column(const std::string& name, align_type al,
                                int width);
        void                add_column(const column_info& format);

        // display column header
        void                disp_header();

        // display next row; row label 'label' is displayed, when set_row_label
        // function is called, otherwise is ignored; number of remaining arguments
        // must be equal to number of columns
        template<class ... Args>
        void                disp_row(const std::string& label, Args&& ... args);        

    public:
        // return true if row labels will be displayed
        bool                display_row_label() const;

        // return number of columns to be displayed
        Integer             number_columns() const;

        // return format of column 'col'
        const column_info&  get_column_format(Integer col) const;

        // return format of row labels
        const column_info&  get_row_label_format() const;

    private:
        void                disp_row_object(const std::string& label, 
                                std::initializer_list<Object> vals);

        template<class Type>
        Object              make_object(Type&& t);        
};

template<class ... Args>
void formatted_disp::disp_row(const std::string& label, Args&& ... args)
{
    disp_row_object(label, {make_object(std::forward<Args>(args))...});
};

template<class Type>
Object formatted_disp::make_object(Type&& t)
{
    return Object(std::forward<Type>(t));
}

};

#pragma warning(pop)