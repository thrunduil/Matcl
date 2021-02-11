/*
 *  This file is a part of Matrix Computation Library (MATCL)
 *
 *  Copyright (c) Pawe³ Kowal 2018 - 2021
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

#include "matcl-matrep/matcl_matrep.h"

namespace matcl { namespace test
{

Matrix	randn_str(Integer m, Integer n);
Matrix	frandn_str(Integer m, Integer n);
Matrix	crandn_str(Integer m, Integer n);
Matrix	fcrandn_str(Integer m, Integer n);

Matrix	sprandn_str(Integer m, Integer n, Real d);
Matrix	fsprandn_str(Integer m, Integer n, Real d);
Matrix	csprandn_str(Integer m, Integer n, Real d);
Matrix	fcsprandn_str(Integer m, Integer n, Real d);

Matrix	randn_band_str(Integer m, Integer n, Integer fd, Integer ld);
Matrix	frandn_band_str(Integer m, Integer n, Integer fd, Integer ld);
Matrix	crandn_band_str(Integer m, Integer n, Integer fd, Integer ld);
Matrix	fcrandn_band_str(Integer m, Integer n, Integer fd, Integer ld);

};};