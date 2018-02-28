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

#include "matcl-matrep/matrix/matrix.h"

namespace matcl { namespace details { namespace scal_func
{

Float_complex   impl_erf(const Float_complex &z);
Complex         impl_erf(const Complex &x);

Float_complex   impl_erfc(const Float_complex &z);
Complex         impl_erfc(const Complex &x);

Float_complex   impl_gamma(const Float_complex &z);
Complex         impl_gamma(const Complex &x);

Float_complex   impl_gammaln(const Float_complex &z);
Complex         impl_gammaln(const Complex &x);

}}}