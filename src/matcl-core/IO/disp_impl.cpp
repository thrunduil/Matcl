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

#include "matcl-core/details/IO/disp_impl.h"
#include "matcl-core/memory/global_objects.h"
#include "matcl-core/details/IO/disp_stream_impl.h"
#include "matcl-core/matrix/matrix_traits.h"
#include "matcl-core/IO/disp_data_provider.h"

#include <algorithm>

#pragma warning(push) 
#pragma warning(disable:4702) //unreachable code

namespace mr = matcl::raw;
namespace mrd = matcl::raw::details;
namespace md = matcl::details;

namespace matcl { namespace raw { namespace details
{

struct matrix_provider_data : public matrix_provider_base<std::string>
{
    disp_data_provider&   m_mat;

    matrix_provider_data(disp_data_provider& mat)
        :m_mat(mat)
    {};

    virtual ~matrix_provider_data(){};

    virtual Integer rows() const override
    {
        return m_mat.rows();
    };
    virtual Integer cols() const override
    {
        return m_mat.cols();
    };

    virtual std::string get_value(const disp_stream* user, Integer width, align_type at,
                                    Integer r, Integer c, bool& is_zero) const override
    {
        std::string ret = m_mat.get_value(user,width,at,r,c);
        is_zero = (ret.size() == 0);
        return std::move(ret);
    };

    virtual bool is_symher() const override
    {
        return false;
    };

    virtual void begin() override               { return m_mat.begin();};
    virtual void hold() override                { return m_mat.hold();};
    virtual void hold_column() override         { return m_mat.hold_column(); };
    virtual void next_row() override            { return m_mat.next_row(); };
    virtual void next_column() override         { return m_mat.next_column(); };
    virtual void restore() override             { return m_mat.restore(); };
    virtual void restore_column() override      { return m_mat.restore_column(); };
    virtual void begin_row_headers() override   { return m_mat.begin_row_headers(); };
    virtual void next_row_header() override     { return m_mat.next_row_header(); };
    virtual void begin_col_headers() override   { return m_mat.begin_col_headers(); };
    virtual void next_col_header() override     { return m_mat.next_col_header(); };
};

struct struct_unknown{};
    
template<>
struct disp_matrix<disp_data_provider,struct_unknown>
{
    static void eval(md::disp_stream_impl& os, const disp_stream* user, 
                        disp_data_provider& mat)
    {            
        os.do_init_display(user,value_code::v_object);
        os.do_start_display(user);	
            
        std::string data_type = mat.get_data_type_description();
        os.do_displaying_data(user,mat.rows(),mat.cols(), data_type, align_type::left);

        if (os.do_display_data(user) == false)
        {
            os.do_end_display(user);
            return;
        };

        matrix_provider_data mdp(mat);
        return disp_matrix<std::string,struct_dense>::eval_matrix_body(os,user,mdp);
    };
};

}}};

namespace matcl { namespace raw
{

//--------------------------------------------------------------------------
//                         disp_stream_data_provider 
//--------------------------------------------------------------------------
disp_stream_data_provider::disp_stream_data_provider(const disp_stream_ptr& impl, 
                                                     const disp_data_provider &dp)
    :forwarding_disp_stream(impl), m_dp(dp)
{
    m_owner = forwarding_disp_stream::get_orginal_stream();
};

void disp_stream_data_provider::start_display(line_printer& p) const
{
    return m_dp.start_display(m_owner,p);
}

void disp_stream_data_provider::end_display(line_printer& p) const
{
    return m_dp.end_display(m_owner,p);
};

void disp_stream_data_provider::display_empty_matrix(line_printer& p, Integer r,
                        Integer c) const
{
    return m_dp.display_empty_matrix(m_owner,p,r,c);
};
bool disp_stream_data_provider::short_print_empty_matrix() const
{
    return m_dp.short_print_empty_matrix(m_owner);
};
bool disp_stream_data_provider::show_matrix_header() const
{
    return m_dp.show_matrix_header(m_owner);
};

void disp_stream_data_provider::display_matrix_name(line_printer& p) const
{
    return m_dp.display_matrix_name(m_owner,p);
};

void disp_stream_data_provider::start_display_matrix_block(line_printer& p, 
                 Integer block_width, Integer first_col, Integer last_col) const
{
    return m_dp.start_display_matrix_block(m_owner,p, block_width, first_col, last_col);
};

void disp_stream_data_provider::end_display_matrix_block(line_printer& p,
                                            Integer block_width) const
{
    return m_dp.end_display_matrix_block(m_owner,p, block_width);
};

bool disp_stream_data_provider::show_column_header_line() const
{
    return m_dp.show_column_header_line(m_owner);
}

bool disp_stream_data_provider::show_column_header_row() const
{
    return m_dp.show_column_header_row(m_owner);
}

bool disp_stream_data_provider::show_column_header_columns() const
{
    return m_dp.show_column_header_columns(m_owner);
}
bool disp_stream_data_provider::can_split() const
{
    return m_dp.can_split(m_owner);
};

align_type disp_stream_data_provider::get_align_row_header() const
{
    return m_dp.get_align_row_header(m_owner);
};

align_type disp_stream_data_provider::get_align_col(Integer c) const
{
    return m_dp.get_align_col(m_owner, c);
};

std::string disp_stream_data_provider::get_row_name(Integer r) const
{
    return m_dp.get_row_name(m_owner,r);
}

std::string disp_stream_data_provider::get_col_name(Integer c) const
{
    return m_dp.get_col_name(m_owner,c);
}

std::string disp_stream_data_provider::get_rows_label() const
{
    return m_dp.get_rows_label(m_owner);
};

void disp_stream_data_provider::get_column_width_row(Integer& w_min, Integer& w_max) const
{
    return m_dp.get_column_width_row(m_owner, w_min, w_max);
};

void disp_stream_data_provider::get_column_width(Integer c, Integer& w_min, Integer& w_max) const
{
    return m_dp.get_column_width(m_owner, c, w_min, w_max);
};

//--------------------------------------------------------------------------
//                       public functions impl
//--------------------------------------------------------------------------

void raw::disp(const disp_stream_ptr& os, disp_data_provider &c)
{	
    disp_stream_data_provider local(os, c);
    mrd::disp_matrix<disp_data_provider, mrd::struct_unknown>::eval(*(local.impl()), &local, c);
};

};};

#pragma warning(pop) 