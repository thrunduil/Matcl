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
 * Solve discrete algebraic Lyapunov equations, DALE (Stein Equation)
 *
 * Solves the DALE (discrete time algebraic Lyapunov equation) of the form:
 * A'XA - X = C
 * where C is symmetric
 *
 * @param   A           NxN matrix
 * @param   C           NxN symmetric matrix
 *
 * @throw   error::error_size_lyapunov                  if arguments nonconformant
 * @throw   error::error_lyapunov                       if method failed to converge 
 * @throw   error::error_complex_value_type_not_allowed if given Copmlex matrix
 * 							    
 * @return  2-tuple with: 
 *          Solution NxN matrix X,
 *          Forward error estimate
 */   
MATCL_LINALG_EXPORT mat_tup_2 solve_dale    (const Matrix& A, const Matrix& C);

/**
 * Solve continuous algebraic Lyapunov equations, CALE
 *
 * Solves the CALE (continuous time algebraic Lyapunov equation) of the form:
 * A'X + XA = C
 * where C is symmetric
 *
 * @param   A           NxN matrix
 * @param   C           NxN symmetric matrix
 *
 * @throw   error::error_size_lyapunov                  if arguments nonconformant
 * @throw   error::error_lyapunov                       if method failed to converge 
 * @throw   error::object_value_type_not_allowed  if given Object matrix
 * @throw   error::error_complex_value_type_not_allowed if given Copmlex matrix
 * 							    
 * @return  2-tuple with: 
 *          Solution NxN matrix X,
 *          Forward error estimate
 */   
MATCL_LINALG_EXPORT mat_tup_2 solve_cale (const Matrix& A, const Matrix& C);

/**
 * Solve sparse discrete algebraic Lyapunov equations, DALE
 *
 * Solves the DALE (discrete time algebraic Lyapunov equation) of the form:
 * A'XA - X = - GG'
 * where A large and sparse, G thin. A must be stable, i.e. abs(eig(A)) < 1
 * X is returned factored as X = ZZ^H
 * for details, see Benner, P., Fassbender, H.: On the numerical solution of large-scale sparse discrete-time Riccati equations. Adv. Comput. Math. 35, 119–147 (2011)
 *
 * @param   A               NxN sparse large matrix
 * @param   G               NxM matrix
 * @param   l_zero_cap      Maximum number of ADI iterations
 *
 * @throw   error::error_size_lyapunov                  if arguments nonconformant
 * @throw   error::error_lyapunov                       if underlying lyapunov method failed to converge 
 * @throw   error::object_value_type_not_allowed  if given Object matrix
 * @throw   error::error_complex_value_type_not_allowed if given Copmlex matrix
 * 							    
 * @return  Solution factor Z, X = ZZ^H
 */   
MATCL_LINALG_EXPORT Matrix sparse_dale(const Matrix& A, const Matrix& G,
                                       const Integer l_zero_cap = 50);

/**
 * Solve sparse continuous algebraic Lyapunov equations, CALE
 *
 * Solves the CALE (continuous time algebraic Lyapunov equation) of the form:
 * A'X + XA = - GG'
 * where A large and sparse, G thin. A must be stable, i.e. real(eig(A)) < 0
 * X is returned factored as X = ZZ^H
 * for details, see  A MATLAB Toolbox for Large Lyapunov and Riccati Equations, Model Reduction Problems, and Linear Quadratic Optimal Control Problems. Users’ Guide (Version 1.0 (1999) by T Penzl, LYAPACK 
 *
 * @param   A               NxN sparse large matrix
 * @param   G               NxM matrix
 * @param   l_zero_cap      Maximum number of ADI iterations
 * @param   V_update_tol    Tolerance to early-stop ADI iterations in case Z updates become small
 *
 * @throw   error::error_size_lyapunov                  if arguments nonconformant
 * @throw   error::error_lyapunov                       if underlying lyapunov method failed to converge 
 * @throw   error::object_value_type_not_allowed  if given Object matrix
 * @throw   error::error_complex_value_type_not_allowed if given Copmlex matrix
 * 							    
 * @return  Solution factor Z, X = ZZ^H
 */   
MATCL_LINALG_EXPORT Matrix sparse_cale(const Matrix& A, const Matrix& G, 
                                       const Integer l_zero_cap = 50, const Real V_update_tol = 1e-9);

}