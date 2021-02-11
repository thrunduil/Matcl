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
Integer irandn_sp(bool with_nan);
Real	randn_sp(bool with_nan, bool with_inf = true);
Float   frandn_sp(bool with_nan, bool with_inf = true);

Matrix	randn_sp(Integer m, Integer n,bool with_nan, bool with_inf = true);
Matrix	frandn_sp(Integer m, Integer n,bool with_nan, bool with_inf = true);
Matrix	crandn_sp(Integer m, Integer n,bool with_nan);
Matrix	fcrandn_sp(Integer m, Integer n,bool with_nan);

Matrix	sprandn_sp(Integer m, Integer n, Real d,bool with_nan, bool with_inf = true);
Matrix	fsprandn_sp(Integer m, Integer n, Real d,bool with_nan, bool with_inf = true);
Matrix	csprandn_sp(Integer m, Integer n, Real d,bool with_nan);
Matrix	fcsprandn_sp(Integer m, Integer n, Real d,bool with_nan);

Matrix	randn_band_sp(Integer m, Integer n, Integer fd, Integer ld, bool with_nan, bool with_inf = true);
Matrix	frandn_band_sp(Integer m, Integer n, Integer fd, Integer ld, bool with_nan, bool with_inf = true);
Matrix	crandn_band_sp(Integer m, Integer n, Integer fd, Integer ld, bool with_nan);
Matrix	fcrandn_band_sp(Integer m, Integer n, Integer fd, Integer ld, bool with_nan);

};};