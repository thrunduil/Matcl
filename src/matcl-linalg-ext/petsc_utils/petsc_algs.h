/*
 *  This file is a part of Morfa Matrix Lib.
 *
 *  Copyright (c) Pawe³ Kowal 2011-2016
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

#if 0
TODO

#include "matcl-linalg/general/config_linalg.h"
#include "matcl-linalg/details/linalg_fwd.h"
#include "matcl-matrep/matrix/matrix.h"
#include "matcl-core/options/matcl_options.h"
#include "matcl-linalg/linear_eq/petsc_enums.h"
#include "matcl-matrep/matrix/permvec.h"

namespace matcl { namespace petsc
{

enum class coarsen_type
{
    hem, mis
};

permvec reorder_sym(const Matrix& mat, petsc_ordering ord_type);

tuple<permvec,permvec>
reorder_nsym(const Matrix& mat, petsc_ordering ord_type);

Matrix coarsen(const Matrix& A, const permvec& perm, coarsen_type ct);

}};

#endif
