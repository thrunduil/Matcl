/*
 *  This file is a part of Matrix Computation Library (MATCL)
 *
 *  Copyright (c) Pawe³ Kowal 2019
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

#include "mkgen/matrix/scalar.h"
#include "mkgen/matrix/predefined_scalars.h"
#include "mkgen/matrix/matrix.h"
#include "mkgen/details/matrix/matrix_printer.h"
#include "mkgen/expression/expressions.h"

#include "mkgen/TODO/expression/ct_matrix_expr.inl"
#include "mkgen/TODO/utils/utils.h"
#include "mkgen/TODO/expression/mat_assign.h"
#include "mkgen/TODO/expression/mat_temporary.h"
#include "mkgen/TODO/expression/mat_other.h"
#include "mkgen/TODO/evaler/temp_storage.h"
#include "mkgen/TODO/evaler/expr_evaler.h"
#include "mkgen/TODO/expression/for_expr.h"
#include "mkgen/TODO/expression/mat_call.h"
#include "mkgen/TODO/matrix/computation.h"
#include "mkgen/TODO/matrix/printer.h"