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
#include "matcl-core/profile/timer.h"

namespace matcl
{

Real norm_1(const Matrix& m);

// return true if m is a scalar representing value 1.0
bool    is_one(const Matrix& m);

// check struct flags
bool    has_struct_diag(const Matrix& m);
bool    has_struct_tril(const Matrix& m);
bool    has_struct_triu(const Matrix& m);
bool    has_struct_qtril(const Matrix& m);
bool    has_struct_qtriu(const Matrix& m);
bool    has_struct_hessl(const Matrix& m);
bool    has_struct_hessu(const Matrix& m);
bool    has_struct_unitary(const Matrix& m);

bool    has_struct(const Matrix& m, struct_flag sf);

// return eps * norm_1(mat)
Real    epsilon_mat(const Matrix& mat);

// return mult * min(mat.rows(), mat.cols()) * epsilon_mat(mat)
Real    error_tolerance(Real mult, const Matrix& mat);

// return mult * min(mat1.rows(), mat1.cols()) 
//             * (epsilon_mat(mat1) + epsilon_mat(mat2))
Real    error_tolerance2(Real mult, const Matrix& mat1, const Matrix& mat2);

// return normwise forward error bound for matrix multiplication, i.e.
// ||C - Cbar|| <= error_mult(1.0, A, B)
// where Cbar = A * B, and C is true value of this matrix product
Real    error_mult(Real mult, const Matrix& A, const Matrix& B);

// return normwise forward error bound for matrix multiplication, i.e.
// ||X -Xbar|| <= error_mult(1.0, A, B, C)
// where Xbar = A * B * C, and X is true value of this matrix product
Real    error_mult(Real mult, const Matrix& A, const Matrix& B, const Matrix& C);

// make scalar of given type; warnings are not printed
template<class V>
Matrix   make_scalar(const V& v, value_code vc);

};