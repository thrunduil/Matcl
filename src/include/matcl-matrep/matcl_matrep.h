/*
 *  This file is a part of Matrix Computation Library (MATCL)
 *
 *  Copyright (c) Pawe³ Kowal 2017 - 2021
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

// main header of matcl-matrep

#include "matcl-core/matrix/matrix_traits.h"
#include "matcl-core/lib_functions/constants.h"

#include "matcl-dynamic/matcl_function_names.h"
#include "matcl-scalar/matcl_scalar.h"

#include "matcl-matrep/lib_functions/matrix_utils.h"
#include "matcl-matrep/matrix/matrix_concat.h"
#include "matcl-matrep/lib_functions/matrix_gen.h"

#include "matcl-matrep/IO/matrix_io.h"
#include "matcl-matrep/lib_functions/manip.h"
#include "matcl-matrep/lib_functions/vecfunc.h"
#include "matcl-matrep/lib_functions/func_unary.h"
#include "matcl-matrep/lib_functions/func_binary.h"
#include "matcl-matrep/lib_functions/func_matrix.h"
#include "matcl-matrep/lib_functions/eval_functors.h"
#include "matcl-matrep/matrix/colon.h"
#include "matcl-matrep/details/matrix.inl"
#include "matcl-matrep/matrix/permvec.h"
#include "matcl-matrep/matrix/matrix_rep_sparse.h"
#include "matcl-matrep/matrix/matrix_rep_dense.h"
#include "matcl-matrep/matrix/matrix_rep_band.h"
#include "matcl-matrep/matrix/matrix_rep_functions.h"
#include "matcl-matrep/matrix/struct_flag_ext.h"

#include "matcl-scalar/object.h"
#include "matcl-matrep/objects/object_matrix.h"
#include "matcl-matrep/matrix/unique_matrix.h"
#include "matcl-matrep/objects/object_matrix.h"
