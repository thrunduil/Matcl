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

#include "matcl-linalg/decompositions/chol.h"
#include "matcl-linalg/decompositions/cholmod.h"
#include "matcl-linalg/decompositions/eig_functions.h"
#include "matcl-linalg/decompositions/givens.h"
#include "matcl-linalg/decompositions/gschur.h"
#include "matcl-linalg/decompositions/gschur_sym.h"
#include "matcl-linalg/decompositions/hess.h"
#include "matcl-linalg/decompositions/householder.h"
#include "matcl-linalg/decompositions/ldl.h"
#include "matcl-linalg/decompositions/lu.h"
#include "matcl-linalg/decompositions/qr.h"
#include "matcl-linalg/decompositions/schur.h"
#include "matcl-linalg/decompositions/svd.h"
#include "matcl-linalg/decompositions/balancing.h"

#include "matcl-linalg/krylov/krylov.h"
#include "matcl-linalg/krylov/arnoldi.h"
#include "matcl-linalg/krylov/block_arnoldi.h"

#include "matcl-linalg/linear_eq/linsolve.h"

#include "matcl-linalg/norms_error/norm.h"
#include "matcl-linalg/norms_error/norm.h"
#include "matcl-linalg/norms_error/cond.h"

//#include "matcl-linalg/matrix_eq/sylvester_equation.h"
//#include "matcl-linalg/matrix_eq/riccati.h"
//#include "matcl-linalg/matrix_eq/lyapunov.h"

#include "matcl-linalg/utils/utils.h"
#include "matcl-linalg/general/control_exceptions.h"
#include "matcl-linalg/special_matrices/struct_flag_linalg.h"
#include "matcl-linalg/graph/matcl_graph.h"

