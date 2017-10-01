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
#include "matcl-core/IO/disp_stream.h"
#include "matcl-core/options/matcl_options.h"

namespace matcl
{

// this disp_stream allows for setting printing options locally
class MATCL_CORE_EXPORT disp_stream_options : public matcl::forwarding_disp_stream
{
    private:
        options           m_options;

    public:
        // use most of disp options as well as output stream from other stream
        // general printing options are taken from opts argument. See options_disp
        // for a list of available options. If given option is not set in opts
        // argument, then appropriate method from other_stream is called, not the
        // value taken from default options or predefined options.
        disp_stream_options(const disp_stream_ptr& other_stream, const options& opts);

        // use global disp stream as base disp_stream
        explicit disp_stream_options(const options& opts);

        // See disp_stream for description of the following methods.
        virtual Integer     get_terminal_width() const override;
        virtual Integer     get_max_cols() const override;
        virtual Integer     get_max_rows() const override;
        virtual Integer     get_max_nnz() const override;
        virtual bool        restrict_sparse_matrix_size() const override;
        virtual bool        display_zero() const override;
        virtual bool        ignore_lower_triangle() const override;
        virtual bool        disp_header_only() const override;
        virtual Integer     get_row_block() const override;
        virtual Integer     get_precision() const override;
        virtual disp_mode   get_disp_mode() const override;
};

};
