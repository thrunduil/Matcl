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
#include "matcl-linalg/special_matrices/unitary_matrix.h"

namespace matcl { namespace details
{

template<class Val>
struct quern_unitary_mat_mult
{
    public:
        Integer     m_rows;
        Integer     m_cols;
        const int*  row_start;
        const int*  column_index;
        const Val*  value;

    public:
        quern_unitary_mat_mult(Integer rows, Integer cols, const int* row_ind, 
                               const int* col_ind, const Val* values);

        void mult_right(matcl::Matrix& ret, const matcl::Matrix& X, trans_type t_unitary) const;
        void mult_left(matcl::Matrix& ret, const matcl::Matrix& X, trans_type t_unitary) const;
};


}};