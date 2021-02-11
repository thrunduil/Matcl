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

#include "unary_function.h"

namespace matcl { namespace test
{

Real unary_function::eval(const Matrix& mat,bool show_res, int code)
{
	m_is_error = false;
	m_error.clear();

	return eval_mat(mat,show_res, code);
};
Real unary_function::eval_scal(const Scalar& s,bool show_res, int code)
{
	m_is_error = false;
	m_error.clear();

	return eval_scalar(s,show_res,code);
};

};};