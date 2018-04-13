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
#include "matcl-linalg/general/config_linalg.h"

//TODO: cleanup

namespace matcl
{
/**
 * Solve Sylvester equation AX + XB = C, for A MxM, B NxN, C MxN matrices.
 *
 * @param   A   The MxM left factor
 * @param   B   The NxN right factor
 * @param   C   The MxN right hand side
 * @param   fast_exit_linsolve_allowed      If false, will never fallback to linsolve
 * 							  
 * @return  Solution matrix X
 */   
MATCL_LINALG_EXPORT 
Matrix solve_sylvester(const Matrix& A, const Matrix& B, const Matrix& C,
                       const bool fast_exit_linsolve_allowed = true);


/**
 * Solve discrete time Sylvester equation AXB + X = C, for A MxM, B NxN, C MxN matrices.
 * @param   A   The MxM left factor
 * @param   B   The NxN right factor
 * @param   C   The MxN right hand side
 * @return  Solution matrix X
 */   
MATCL_LINALG_EXPORT 
Matrix solve_dsylvester(const Matrix& A, const Matrix& B, const Matrix& C);


/**
* Solve generalized Sylvester equation
*
*       A1 X + Y A2 = C1,
*       B1 X + Y B2 = C2
*
* A1, B1 are NxN matrices
* A2, B2 are MxM matrices
* C1, C2 are NxM matrices

* returns [X,Y]
 */   

MATCL_LINALG_EXPORT 
mat_tup_2 solve_gsylvester(const Matrix& A1, const Matrix& A2, 
                           const Matrix& B1, const Matrix& B2, 
                           const Matrix& C1, const Matrix& C2);

};
