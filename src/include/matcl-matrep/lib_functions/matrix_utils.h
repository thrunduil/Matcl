/*
 *  This file is a part of Matrix Computation Library (MATCL)
 *
 *  Copyright (c) Pawe³ Kowal 2017 - 2018
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

#include "matcl-scalar/lib_functions/utils.h"
#include "matcl-matrep/matrix/matrix.h"

namespace matcl
{

// start global timer
// redeclaration of function defined in matcl-scalar
MATCL_SCALAR_EXPORT 
void            tic();

// stop global timer and return time elapsed from last tic in seconds
// redeclaration of function defined in matcl-scalar
MATCL_SCALAR_EXPORT
Real            toc();

// stop global timer and return time elapsed from last tic in seconds
// represented as string
// redeclaration of function defined in matcl-scalar
MATCL_SCALAR_EXPORT
std::string     tocstr();

// stop global timer; calculatr time elapsed from last tic in seconds
// and print elapsed time on global_output_stream
// redeclaration of function defined in matcl-scalar
MATCL_SCALAR_EXPORT
void            tocdisp();

};