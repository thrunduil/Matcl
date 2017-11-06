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
#include <fstream>
#include <memory>

namespace matcl 
{

// set log file; every messages printed to global output stream
// (for example when disp function is called on global output stream)
// will also be printed to the log file
MATCL_CORE_EXPORT 
void                set_logger(const std::shared_ptr<std::ofstream>& log_file);

// return pointer to previously set log file; return empty pointer
// if log file was not set
MATCL_CORE_EXPORT
std::shared_ptr<std::ofstream>
                    get_logger();   

// return previously set logger as output stream; this stream is synchronized; 
// return empty pointer if log file was not set
MATCL_CORE_EXPORT 
output_stream_ptr   get_logger_output_stream();

}