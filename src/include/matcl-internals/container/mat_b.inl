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

#include "matcl-internals/error/error_check_basic.h"
#include "matcl-internals/container/mat_b.h"
#include "matcl-core/details/integer.h"

namespace matcl { namespace raw 
{

template <class value_type>
inline const value_type& Matrix<value_type,struct_banded>::operator()(Integer i, Integer j) const
{
    error::check_index_band(i, j, m_rows, m_cols, m_ldiags, m_udiags);
    return base_type::ptr()[m_udiags + i - j + (j-1)*ld()];
}

template <class value_type>
inline Integer Matrix<value_type,struct_banded>::rows() const
{
    return m_rows;
}

template <class value_type>
inline Integer Matrix<value_type,struct_banded>::cols() const
{
    return m_cols;
}

template <class value_type>
inline Integer Matrix<value_type,struct_banded>::max_udiags() const
{
    auto off        = Matrix<value_type,struct_banded>::base_type::m_data.m_root_ptr
                        .offset(Matrix<value_type,struct_banded>::base_type::m_data.m_ptr);

    return m_udiags + Integer(off) % ld();
};

template <class value_type>
inline Integer Matrix<value_type,struct_banded>::max_ldiags() const
{
    return Matrix<value_type,struct_banded>::max_rows() - max_udiags() - 1;
};

template <class value_type>
inline Integer Matrix<value_type,struct_banded>::length() const
{ 
    return (rows() == 0 || cols() == 0)? 0 : std::max(rows(), cols()); 
};

template <class value_type>
inline const value_type* Matrix<value_type,struct_banded>::rep_ptr() const
{
    return base_type::ptr();
};

template <class value_type>
inline value_type* Matrix<value_type,struct_banded>::rep_ptr()
{
    return base_type::ptr();
};

template <class value_type>
inline Integer Matrix<value_type,struct_banded>::first_diag() const
{
    return -m_ldiags;
}

template <class value_type>
inline Integer Matrix<value_type,struct_banded>::last_diag() const
{
    return m_udiags;
}

template <class value_type>
inline Integer Matrix<value_type,struct_banded>::number_superdiagonals() const
{
    return m_udiags;
}

template <class value_type>
inline Integer Matrix<value_type,struct_banded>::number_subdiagonals() const
{
    return m_ldiags;
}

template <class value_type>
inline Integer Matrix<value_type,struct_banded>::first_row(Integer col) const
{
    return (col - m_udiags < 0) ? 0 : col - m_udiags;
};

template <class value_type>
inline Integer Matrix<value_type,struct_banded>::last_row(Integer col) const
{
    Integer r = rows() - 1;
    return (r < col + m_ldiags) ? r : col + m_ldiags;
};

template <class value_type>
inline Integer Matrix<value_type,struct_banded>::first_nonzero_column() const
{
    return std::max(first_diag(), 0);
};

template <class value_type>
inline Integer Matrix<value_type,struct_banded>::last_nonzero_column() const
{
    return std::min(rows() - 1 + last_diag(), cols() - 1);
};

template <class value_type>
inline Integer Matrix<value_type,struct_banded>::first_nonzero_row() const
{
    return std::max(-last_diag(), 0);
};

template <class value_type>
inline Integer Matrix<value_type,struct_banded>::last_nonzero_row() const
{
    return std::min(cols() - 1 - first_diag(), rows() - 1);
};

template <class value_type>
inline Integer Matrix<value_type,struct_banded>::first_col(Integer r) const
{
    return std::max(0, r - m_ldiags);
};

template <class value_type>
inline Integer Matrix<value_type,struct_banded>::last_col(Integer r) const
{
    return std::min(m_cols - 1, r + m_udiags);
};


template <class value_type>
inline Integer Matrix<value_type,struct_banded>::first_elem_pos(Integer col) const
{
    return (m_udiags < col ) ? 0 : m_udiags - col;
};

template <class value_type>
inline Integer Matrix<value_type,struct_banded>::first_elem_pos_row(Integer row) const
{
    Integer col = first_col(row);
    return m_udiags + row - col + imult(col,ld());
};

template <class value_type>
inline Integer Matrix<value_type,struct_banded>::first_elem_diag(Integer d) const
{
    return (d <= 0) ? m_udiags - d : m_udiags-d + imult(d,ld());
};

template <class value_type>
inline bool Matrix<value_type,struct_banded>::has_diag(Integer d) const
{
    return d <= last_diag() && d >= first_diag();
};

template <class value_type>
inline Integer Matrix<value_type,struct_banded>::diag_length(Integer d) const
{
    return (d >= 0) ? std::min(m_rows, m_cols - d) : std::min(m_rows + d, m_cols);
};

template <class value_type>
inline Integer Matrix<value_type,struct_banded>::first_row_on_diag(Integer d) const
{
    return (d >= 0) ? 0 : -d;
};

template <class value_type>
inline Integer Matrix<value_type,struct_banded>::first_col_on_diag(Integer d) const
{
    return (d >= 0) ? d : 0;
};

template <class value_type>
inline Integer Matrix<value_type,struct_banded>::element_pos(Integer r, Integer c) const
{
    return m_udiags + r + c * (ld() - 1);
};

template <class value_type>
inline Integer Matrix<value_type,struct_banded>::ld() const
{
    return base_type::ld();
};

};};