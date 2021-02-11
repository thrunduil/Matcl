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

#include "matcl-core/details/IO/disp_stream_options.h"
#include "matcl-core/memory/global_objects.h"
#include "matcl-core/options/options_disp.h"

namespace matcl
{

disp_stream_options::disp_stream_options(const disp_stream_ptr& other_stream, const options& opts)
    :forwarding_disp_stream(other_stream), m_options(opts)
{};
disp_stream_options::disp_stream_options(const options& opts)
    :forwarding_disp_stream(default_disp_stream()), m_options(opts)
{};

Integer disp_stream_options::get_terminal_width() const
{
    optional<Integer> opt = m_options.get<Integer>(matcl::opt::disp::terminal_width());
    if (opt)
        return opt.value();
    else
        return forwarding_disp_stream::get_terminal_width();
};
Integer disp_stream_options::get_max_cols() const
{
    optional<Integer> opt = m_options.get<Integer>(matcl::opt::disp::max_cols());
    if (opt)
        return opt.value();
    else
        return forwarding_disp_stream::get_max_cols();
};
Integer disp_stream_options::get_max_rows() const
{
    optional<Integer> opt = m_options.get<Integer>(matcl::opt::disp::max_rows());
    if (opt)
        return opt.value();
    else
        return forwarding_disp_stream::get_max_rows();
}
Integer disp_stream_options::get_max_nnz() const
{
    optional<Integer> opt = m_options.get<Integer>(matcl::opt::disp::max_nnz());
    if (opt)
        return opt.value();
    else
        return forwarding_disp_stream::get_max_nnz();
};
Integer disp_stream_options::get_precision() const
{
    optional<Integer> opt = m_options.get<Integer>(matcl::opt::disp::precision());
    if (opt)
        return opt.value();
    else
        return forwarding_disp_stream::get_precision();
};
bool disp_stream_options::restrict_sparse_matrix_size() const
{
    optional<bool> opt = m_options.get<bool>(matcl::opt::disp::restrict_sparse_matrix_size());
    if (opt)
        return opt.value();
    else
        return forwarding_disp_stream::restrict_sparse_matrix_size();
};
bool disp_stream_options::display_zero() const
{
    optional<bool> opt = m_options.get<bool>(matcl::opt::disp::display_zero());
    if (opt)
        return opt.value();
    else
        return forwarding_disp_stream::display_zero();
};
bool disp_stream_options::ignore_lower_triangle() const
{
    optional<bool> opt = m_options.get<bool>(matcl::opt::disp::ignore_lower_triangle());
    if (opt)
        return opt.value();
    else
        return forwarding_disp_stream::ignore_lower_triangle();
};
bool disp_stream_options::disp_header_only() const
{
    optional<bool> opt = m_options.get<bool>(matcl::opt::disp::header_only());
    if (opt)
        return opt.value();
    else
        return forwarding_disp_stream::disp_header_only();
};
Integer disp_stream_options::get_row_block() const
{
    optional<Integer> opt = m_options.get<Integer>(matcl::opt::disp::row_block());
    if (opt)
        return opt.value();
    else
        return forwarding_disp_stream::get_row_block();
};
disp_mode disp_stream_options::get_disp_mode() const
{
    optional<Integer> opt = m_options.get<Integer>(matcl::opt::disp::disp_mode());
    if (opt)
        return static_cast<disp_mode>(opt.value());
    else
        return forwarding_disp_stream::get_disp_mode();
};

};
