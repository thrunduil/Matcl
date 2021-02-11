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

#include "matcl-linalg/graph/matcl_graph.h"
//#include "matcl-linalg/petsc_utils/petsc_algs.h"
#include "matcl-linalg/general/linalg_exception.h"

namespace matcl
{

Matrix matcl::coarser_hem(const Matrix& A)
{
    (void)A;
    //TODO
    throw;
    //petsc::coarsen_type ct  = petsc::coarsen_type::hem;
    //permvec perm            = permvec::identity(A.rows());
    //return petsc::coarsen(A, perm, ct);
};

Matrix matcl::coarser_mis(const Matrix& A)
{
    //TODO
    (void)A;
    throw;

    //petsc::coarsen_type ct  = petsc::coarsen_type::mis;
    //permvec perm            = permvec::identity(A.rows());
    //return petsc::coarsen(A, perm, ct);
};

};