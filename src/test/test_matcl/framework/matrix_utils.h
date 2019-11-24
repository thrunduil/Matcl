/*
 *  This file is a part of Matrix Computation Library (MATCL)
 *
 *  Copyright (c) Pawe³ Kowal 2018
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

#include "matcl-matrep/matrix/matrix.h"
#include "matcl-core/profile/timer.h"

namespace matcl
{

Real norm_1(const Matrix& m);

// check struct flags
bool    has_struct_diag(const Matrix& m);
bool    has_struct_tril(const Matrix& m);
bool    has_struct_triu(const Matrix& m);
bool    has_struct_qtril(const Matrix& m);
bool    has_struct_qtriu(const Matrix& m);
bool    has_struct_hessl(const Matrix& m);
bool    has_struct_hessu(const Matrix& m);

bool    has_struct(const Matrix& m, struct_flag sf);

};