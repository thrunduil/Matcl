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

#include "matcl-core/IO/base_io.h"
#include "matcl-core/options/options_disp.h"
#include "matcl-core/details/disp_stream_options.h"
#include "matcl-core/details/disp_impl.h"

namespace matcl
{

void matcl::disp_header(disp_data_provider& m, const disp_stream_ptr& os, const options& opts)
{
    options opts2 = opts;
    opts2.set(matcl::opt::disp::header_only(true));
    return disp(m, os, opts2);
};

void matcl::disp(disp_data_provider& m, const disp_stream_ptr& os, const options& opts)
{
    if (opts.size() == 0)
        return raw::disp(os,m);
    else
        return raw::disp(std::make_shared<disp_stream_options>(os, opts), m);
};

};
