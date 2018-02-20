/*
 *  This file is a part of Matrix Computation Library (MATCL)
 *
 *  Copyright (c) Pawe³ Kowal 2017 - 2018
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

#include "matcl-internals/container/mat_d.h"
#include "matcl-internals/error/error_check_basic.h"

namespace matcl { namespace raw 
{
    
template <class value_type>
force_inline const value_type& 
Matrix<value_type,struct_dense>::operator()(Integer i, Integer j) const
{
    Integer r = Matrix<value_type,struct_dense>::base_type::m_data.m_rows;
    Integer c = Matrix<value_type,struct_dense>::base_type::m_data.m_cols;
    error::check_index(i, j, r, c);

    Integer pos = i-1+(j-1)*Matrix<value_type,struct_dense>::ld();
    return Matrix<value_type,struct_dense>::base_type::ptr()[pos];
}

};};
