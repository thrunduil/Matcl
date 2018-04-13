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

#include "matcl-linalg/norms_error/cond.h"
#include "matcl-matrep/matcl_matrep.h"

#include "matcl-linalg/general/linalg_exception.h"
#include "matcl-linalg/decompositions/svd.h"

namespace matcl
{

Real cond(const Matrix& A)
{
    if (!A.is_square())
        throw error::square_matrix_required(A.rows(), A.cols());
    if (A.numel() <= 1) return 1.;
    Matrix s1 = svd1(A,svd_algorithm::dc);
    return (s1(1) / s1(end)).get_scalar<Real>();
}

}
