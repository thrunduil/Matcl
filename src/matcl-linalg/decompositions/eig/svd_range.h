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

#include "matcl-matrep/matrix/matrix.h"
#include "matcl-linalg/general/config_linalg.h"
#include "matcl-linalg/linear_eq/linsolve_object.h"
#include "matcl-linalg/special_matrices/unitary_matrix.h"
#include "matcl-matrep/details/enablers.h"
#include "matcl-linalg/decompositions/schur.h"

namespace matcl
{

struct svd_range
{
    static void eval_range(Matrix& ret, const Matrix& A, Real VL, Real VU);
    static void eval_index(Matrix& ret, const Matrix& A, Integer IF, Integer IL);
    static void eval2_range(mat_tup_3& ret, const Matrix& A, Real VL, Real VU);
    static void eval2_index(mat_tup_3& ret, const Matrix& A, Integer IF, Integer IL);
};

}

