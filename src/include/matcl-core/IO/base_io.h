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
#include "matcl-core/IO/disp_data_provider.h"
#include "matcl-core/options/matcl_options.h"

namespace matcl
{

namespace md = matcl::details;

//--------------------------------------------------------------------
//                      PRETTY PRINTING
//--------------------------------------------------------------------

// Display any data using the same formatting rules as for matrices. Display
// data using global disp stream or local disp stream. Data are represented by
// abstract class disp_data_provider defined elsewhere. Options controls how printing
// is performed, see options_disp for details
MATCL_CORE_EXPORT 
void            disp(disp_data_provider& data, const disp_stream_ptr& os 
                    = default_disp_stream(), const options& opts = options());

// Display heder of any data using the same formatting rules as for matrices. Display
// data using global disp stream or local disp stream. Data are represented by
// abstract class disp_data_provider defined elsewhere. Options controls how printing
// is performed, see options_disp for details
MATCL_CORE_EXPORT
void            disp_header(disp_data_provider& data, const disp_stream_ptr& os 
                    = default_disp_stream(), const options& opts = options());

};
