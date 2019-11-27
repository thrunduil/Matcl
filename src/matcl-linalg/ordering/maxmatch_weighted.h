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

#include "matcl-core/matrix/scalar_types.h"

namespace matcl { namespace details
{

// Provides the core Hungarian Algorithm implementation for solving the
// minimum or maximum sum assignment problem as per Duff and Koster.

// M        number of rows
// N        number of columns
// ptr_c    column pointers, size N + 1
// ptr_r    row pointers, size nnz
// val      value of the entry that corresponds to row, all values val(k) must 
//          be non-negative, size nnz
// iperm    matching itself: row i is matched to column iperm(i); size M
// jperm    matching itself: column i is matched to row jperm(i); size N
// num      cardinality of the matching   
// dualu    is the reduced weight for rows, size M
// dualv    is the reduced weight for cols, size N
// iwork    size 3 * N + 2 * M
// work     size max(M,N)
template<class V>
void hungarian_match(bool minimum, Integer M, Integer N, const Integer* ptr_c, const Integer* ptr_r, 
        const V* val, Integer* iperm, Integer* jperm, Integer& num, V* dualu, V* dualv, Integer* iwork,
        V* work);

}};
