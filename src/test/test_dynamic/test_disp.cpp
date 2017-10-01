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

#include "test_disp.h"

#include "matcl-core/IO/disp_data_provider.h"
#include "matcl-core/IO/base_io.h"
#include "matcl-scalar/object.h"
#include "matcl-scalar/IO/scalar_io.h"

#include <map>
#include <iostream>

namespace matcl { namespace test
{

void test::test_disp()
{
    disp_tester dt;
    dt.make();
};

class column_info
{
    private:
        std::string         m_name;
        align_type          m_align;
        int                 m_width;

    public:
        column_info();
        column_info(const std::string& name, align_type align, int width);

    public:
        const std::string&  name() const;
        align_type          align() const;
        int                 width() const;
};

column_info::column_info()
    :m_width(0)
{}

column_info::column_info(const std::string& name, align_type align, int width)
    :m_name(name), m_align(align), m_width(width)
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

/*
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
        // Name of data set
        virtual void        display_matrix_name(const disp_stream* orig, line_printer& p) const;

        // A new block of dense matrix will be displayed for columns [first_col, last_col].
        virtual void        start_display_matrix_block(const disp_stream* orig, line_printer& p, 
                                    Integer block_width, Integer first_col, Integer last_col) const;

        // Finalize current block of data
        virtual void        end_display_matrix_block(const disp_stream* orig, line_printer& p, 
                                    Integer block_width) const;

        // Display information about columns?
        virtual bool        show_column_header(const disp_stream* orig) const;
        // Display information about rows?
        virtual bool        show_row_headers(const disp_stream* orig)  const;
*/

class formatted_disp : public disp_data_provider
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
        formatted_disp(const disp_stream_ptr& os = default_disp_stream(),
                        const options& opts = options())
            :m_header_printed(false), m_disp_stream(os), m_options(opts)
        {};

        formatted_disp(const formatted_disp&) = delete;
        formatted_disp& operator=(const formatted_disp&) = delete;

    public:
        void                set_row_label(const std::string& name, align_type al,
                                int width);
        void                add_column(const std::string& name, align_type al,
                                int width);

        void                disp_header();

        template<class ... Args>
        void                disp_row(const std::string& label, Args&& ... args);        

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

        virtual void        begin_row_headers() override;
        virtual void        next_row_header() override;

        virtual std::string get_value(const disp_stream* ds, Integer width, align_type at, 
                                      Integer r, Integer c) const override;

        virtual std::string get_row_name(const disp_stream* orig, Integer r) const;

        virtual std::string get_rows_label(const disp_stream* orig) const override;
        virtual align_type  get_align_row_header(const disp_stream* orig) const override;
        virtual std::string get_col_name(const disp_stream* orig, Integer c) const override;                
        virtual align_type  get_align_col(const disp_stream* orig, Integer c) const override;

    private:
        template<class Type>
        Object              make_object(Type&& t);

        void                disp_row_object(const std::string& label, 
                                std::initializer_list<Object> vals);
};

void formatted_disp::set_row_label(const std::string& name, align_type al, int width)
{
    m_row_label = column_info(name, al, width);
};

void formatted_disp::add_column(const std::string& name, align_type al,
                        int width)
{
    m_cols.push_back(column_info(name, al, width));
};

std::string formatted_disp::get_rows_label(const disp_stream* orig) const
{
    (void)orig;
    return m_row_label.name();
};

align_type formatted_disp::get_align_row_header(const disp_stream* orig) const
{
    (void)orig;
    return m_row_label.align();
};

std::string formatted_disp::get_col_name(const disp_stream* orig, Integer c) const
{
    (void)orig;
    return m_cols[c].name();
};

align_type formatted_disp::get_align_col(const disp_stream* orig, Integer c) const
{
    (void)orig;
    return m_cols[c].align();
};

std::string formatted_disp::get_row_name(const disp_stream* orig, Integer r) const
{
    (void)orig;
    (void)r;
    return m_row_name;
}

Integer formatted_disp::cols() const 
{ 
    return (Integer)m_cols.size();
};

Integer formatted_disp::rows() const
{
    return m_header_printed == false ? 0 : 1;
}

void formatted_disp::begin()
{
    m_col   = 0;
};

void formatted_disp::hold()
{
    m_col_V1    = m_col;
};

void formatted_disp::hold_column()
{
    m_col_V2    = m_col;
};

void formatted_disp::next_row()
{};

void formatted_disp::next_column()
{
    ++m_col;
};

void formatted_disp::restore()
{
    m_col       = m_col_V1;
};

void formatted_disp::restore_column()
{
    m_col       = m_col_V2;
};

void formatted_disp::begin_row_headers()
{}

void formatted_disp::next_row_header()
{};

std::string formatted_disp::get_value(const disp_stream* ds, Integer width, align_type at, 
                                Integer r, Integer c) const
{ 
    (void)r;

    const Object& val   = m_data[c];
    align_type align    = m_cols[c].align();
    width               = m_cols[c].width();
    return this->to_string(ds, width, val, align);
};

void formatted_disp::disp_header()
{
    if (m_header_printed == true)
        return;

    matcl::disp(*this, m_disp_stream, m_options);
    m_header_printed    = true;
}

template<class Type>
Object formatted_disp::make_object(Type&& t)
{
    return Object(std::forward<Type>(t));
}

template<class ... Args>
void formatted_disp::disp_row(const std::string& label, Args&& ... args)
{
    disp_row_object(label, {make_object(std::forward<Args>(args))...});
};

void formatted_disp::disp_row_object(const std::string& label, 
                        std::initializer_list<Object> vals)
{
    if (vals.size() != m_cols.size())
        throw std::runtime_error("TODO");

    m_data.resize(m_cols.size());

    for (const Object& elem : vals)
        m_data.push_back(elem);

    disp_header();
    matcl::disp(*this, m_disp_stream, m_options);
};

void disp_tester::make()
{
    Real v = 1.0/13;
    Object obj(v);

    disp(v);
    disp(obj);    

    std::cout << v << "\n";
    std::cout << obj << "\n";

    std::cout << to_string(v) << "\n";        
    std::cout << to_string(obj) << "\n";

    formatted_disp dm;

    dm.set_row_label("it", align_type::left, 5);
    dm.add_column("value", align_type::left, 19);
    dm.add_column("diff", align_type::left, 10);
    dm.add_column("norm", align_type::left, 10);

    dm.disp_header();

    dm.disp_row("1", 1.0, "abc1", 1.0/13);
    dm.disp_row("2", 2.0, "abc2", 1.0/14);
    dm.disp_row("3", 2.0, "abc4", 1.0/15);
};

}};