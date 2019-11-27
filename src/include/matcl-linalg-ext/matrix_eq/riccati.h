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
 * Solve discrete algebraic riccati equations, DARE
 *
 * Solves the DARE (discrete time algebraic Riccati equation) of the form:
 * X = A'XA - (L+A'XB)(R + B'XB)^(-1)(L+A'XB)' + Q
 * where Q = C'C, R = D'D and L = C'D for some C and D
 * WARNING: R + B'XB may be singular for the X returned
 *
 * @author  Piotr Dobaczewski
 * @date    16/01/2013
 *
 * @param   A           NxN matrix
 * @param   B           NxM matrix
 * @param   Q           NxN matrix, Q = C'C
 * @param   R           MxM matrix, R = D'D
 * @param   L           NxM matrix, L = C'D
 *
 * @throw   error::error_size_riccati                   if arguments nonconformant
 * @throw   error::error_riccati                        if method failed to converge 
 * @throw   error::object_value_type_not_allowed  if given Object matrix
 * @throw   error::error_complex_value_type_not_allowed if given Copmlex matrix
 * 							    
 * @return  3-tuple with: 
 *          Solution NxN symmetric matrix X,
 *          An estimate of the reciprocal of the condition number (in the 1-norm) of the N-th order system of algebraic equations from which the solution matrix X is obtained.
 *          Spectrum of the 2N-by-2N matrix pair
 */   
MATCL_LINALG_EXPORT mat_tup_3 solve_dare    (const Matrix& A, const Matrix& B, const Matrix& Q,
                                             const Matrix& R, const Matrix& L);

/**
 * Solve continous algebraic riccati equations, CARE
 *
 * Solves the CARE (continous time algebraic Riccati equation) of the form:
 * Q + A'X + XA - (L+XB)R^(-1) (L+XB)' = 0    
 * where Q = C'C, R = D'D and L = C'D for some C and D
 * WARNING: No checks for singularity of R performed
 *
 * @param   A           NxN matrix
 * @param   B           NxM matrix
 * @param   Q           NxN matrix, Q = C'C
 * @param   R           MxM matrix, R = D'D
 * @param   L           NxM matrix, L = C'D
 *
 * @throw   error::error_size_riccati                   if arguments nonconformant
 * @throw   error::error_riccati                        if method failed to converge 
 * @throw   error::object_value_type_not_allowed  if given Object matrix
 * @throw   error::error_complex_value_type_not_allowed if given Copmlex matrix
 * 							  
 * @return  3-tuple with: 
 *          Solution NxN symmetric matrix X,
 *          An estimate of the reciprocal of the condition number (in the 1-norm) of the N-th order system of algebraic equations from which the solution matrix X is obtained.
 *          Spectrum of the 2N-by-2N matrix pair
 */   
MATCL_LINALG_EXPORT mat_tup_3 solve_care (const Matrix& A, const Matrix& B, const Matrix& Q,
                                          const Matrix& R, const Matrix& L);

/**
 * Solve sparse discrete algebraic Riccati equations, DARE
 *
 * Solves the DARE (discrete time algebraic Riccati equation) of the form:
 * C'QC + A'XA - X - (A'XB + C'S')(R + B'XB)^{-1}(A'XB + C'S')' = 0
 * where Q > 0 symmetric, R > 0 symmetric, A large and sparse, S, C and B thin.
 * [ Q  S ]
 * [ S' R ] must be also positive definite.
 * A must be stable, i.e. abs(eig(A)) < 1.
 * X is returned factored as X = ZZ^H
 * for details, see Benner, P., Fassbender, H.: On the numerical solution of large-scale sparse discrete-time Riccati equations. Adv. Comput. Math. 35, 119–147 (2011)
 *
 * @param   A               NxN sparse large matrix
 * @param   B               NxM matrix
 * @param   C               QxN matrix
 * @param   Q               QxQ symmetric positive definite matrix
 * @param   R               MxM symmetric positive definite matrix
 * @param   S               MxQ matrix
 * @param   lyapunov_tol    Tolerance to early-stop Newton iterations if Lyapunov equation matrices converge
 * @param   kmax            Maximum number of Newton iterations
 * @param   l_zero_cap      Maximum number of ADI iterations for Lyapunov equation
 *
 * @throw   error::error_size_riccati                   if arguments nonconformant
 * @throw   error::error_lyapunov                       if underlying lyapunov method failed to converge 
 * @throw   error::error::error_nonposdef          if Q, R or [Q S; S' R] is not positive definite
 * @throw   error::object_value_type_not_allowed  if given Object matrix
 * @throw   error::error_complex_value_type_not_allowed if given Copmlex matrix
 * 							    
 * @return  Solution factor Z, X = ZZ^H
 */   
MATCL_LINALG_EXPORT Matrix sparse_dare(const Matrix& A, const Matrix& B, const Matrix& C, const Matrix& Q,
                                       const Matrix& R, const Matrix& S, const Real lyapunov_tol = 1e-14, 
                                       const Integer kmax = 20, 
                                       const Integer l_zero_cap = 50);

/**
 * Solve sparse continuous algebraic Riccati equations, CARE
 *
 * Solves the CARE (continuous time algebraic Riccati equation) of the form:
 * C'QC + A'X + XA - XBR^{-1}B'X' = 0
 * where Q > 0 symmetric, R > 0 symmetric, A large and sparse, C and B thin.
 * A must be stable, i.e. real(eig(A)) < 0.
 * X is returned factored as X = ZZ^H
 * for details, see  A MATLAB Toolbox for Large Lyapunov and Riccati Equations, Model Reduction Problems, and Linear Quadratic Optimal Control Problems. Users’ Guide (Version 1.0 (1999) by T Penzl, LYAPACK 
 *
 * @param   A               NxN sparse large matrix
 * @param   B               NxM matrix
 * @param   C               QxN matrix
 * @param   Q               QxQ symmetric positive definite matrix
 * @param   R               MxM symmetric positive definite matrix
 * @param   lyapunov_tol    Tolerance to early-stop Newton iterations if Lyapunov equation matrices converge
 * @param   kmax            Maximum number of Newton iterations
 * @param   l_zero_cap      Maximum number of ADI iterations for Lyapunov equation
 * @param   V_update_tol    Tolerance to early-stop ADI iterations for Lyapunov Equation
 *
 * @throw   error::error_size_riccati                   if arguments nonconformant
 * @throw   error::error_lyapunov                       if underlying lyapunov method failed to converge 
 * @throw   error::error::error_nonposdef          if Q or R is not positive definite
 * @throw   error::object_value_type_not_allowed  if given Object matrix
 * @throw   error::error_complex_value_type_not_allowed if given Copmlex matrix
 * 							    
 * @return  Solution factor Z, X = ZZ^H
 */   
MATCL_LINALG_EXPORT Matrix sparse_care(const Matrix& A, const Matrix& B, const Matrix& C, const Matrix& Q,
                                       const Matrix& R, const Real lyapunov_tol = 1e-14, 
                                       const Integer kmax = 20, 
                                       const Integer l_zero_cap = 50, const Real V_update_tol = 1e-9);

}