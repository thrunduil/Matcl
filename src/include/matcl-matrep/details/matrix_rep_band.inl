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

#include "matcl-matrep/matrix/matrix_rep_band.h"

namespace matcl
{

template<class T>
const T& band_matrix<T,true>::elem(Integer r, Integer c) const
{
    return m_base_ptr[m_udiags + r - c + (c-1) * m_base_ld];
};

template<class T>
const T& band_matrix<T,true>::elem(Integer r, Integer c)
{
    return m_base_ptr[m_udiags + r - c + (c-1) * m_base_ld];
};

template<class T>
T& band_matrix<T,false>::elem(Integer r, Integer c)
{
    return m_base_ptr[m_udiags + r - c + (c-1) * m_base_ld];
};

template<class T>
Integer	band_matrix<T,true>::first_row(Integer col) const
{
    return (col - m_udiags < 0) ? 0 : col - m_udiags;
};

template<class T>
Integer	band_matrix<T,true>::last_row(Integer col) const
{
    Integer r = rows() - 1;
    return (r < col + m_ldiags) ? r : col + m_ldiags;
};

template<class T>
Integer	band_matrix<T,true>::first_col(Integer r) const
{
    return std::max(0, r - m_ldiags);
};

template<class T>
Integer	band_matrix<T,true>::last_col(Integer r) const
{
    return std::min(m_cols - 1, r + m_udiags);
};

template<class T>
Integer	band_matrix<T,true>::first_nonzero_column() const
{
    return std::max(first_diag(), 0);
};

template<class T>
Integer	band_matrix<T,true>::last_nonzero_column() const
{
    return std::min(rows() - 1 + last_diag(), cols() - 1);
};

template<class T>
Integer	band_matrix<T,true>::first_nonzero_row() const
{
    return std::max(-last_diag(), 0);
};

template<class T>
Integer	band_matrix<T,true>::last_nonzero_row() const
{
    return std::min(cols() - 1 - first_diag(), rows() - 1);
};

template<class T>
Integer	band_matrix<T,true>::first_elem_pos(Integer col) const
{
    return (m_udiags < col ) ? 0 : m_udiags - col;
};

template<class T>
Integer	band_matrix<T,true>::first_elem_pos_row(Integer row) const
{
    Integer col = first_col(row);
    return m_udiags + row - col + imult(col,ld());
};

template<class T>
Integer	band_matrix<T,true>::first_elem_diag(Integer d) const
{
    return (d <= 0) ? m_udiags - d : m_udiags-d + imult(d,ld());
};

template<class T>
Integer band_matrix<T,true>::element_pos(Integer r, Integer c) const
{
    return m_udiags + r + c * (m_base_ld - 1);
};

template<class T>
bool band_matrix<T,true>::has_diag(Integer d) const
{
    return (d >= 0) ? d <= m_cols - 1 : -d <= m_rows - 1;
};

template<class T>
Integer band_matrix<T,true>::diag_length(Integer d) const
{
    return (d >= 0) ? std::min(m_rows, m_cols - d) : std::min(m_rows + d, m_cols);
};

template<class T>
Integer band_matrix<T,true>::first_row_on_diag(Integer d) const
{
    return (d >= 0) ? 0 : -d;
};

template<class T>
Integer band_matrix<T,true>::first_col_on_diag(Integer d) const
{
    return (d >= 0) ? d : 0;
};

template<class T>
inline sub_band_matrix<T>::sub_band_matrix(matrix_type* m, const colon& c1)
    :m_matrix(m), m_colon_1(&c1),m_colon_2(nullptr),m_d(0)
{}

template<class T>
inline sub_band_matrix<T>::sub_band_matrix(matrix_type* m, const colon& c1, const colon& c2)
    :m_matrix(m), m_colon_1(&c1),m_colon_2(&c2),m_d(0)
{};

template<class T>
inline sub_band_matrix<T>::sub_band_matrix(Integer diag, matrix_type* m)
    :m_matrix(m), m_colon_1(nullptr), m_colon_2(nullptr), m_d(diag)
{};

template<class T>
inline sub_band_matrix<T>::sub_band_matrix(const sub_band_matrix& m)
    :m_matrix(m.m_matrix),m_colon_1(m.m_colon_1),m_colon_2(m.m_colon_2),m_d(m.m_d)
{};

template<class T>
inline sub_band_matrix<T>::sub_band_matrix(sub_band_matrix&& m)
    :m_matrix(m.m_matrix),m_colon_1(m.m_colon_1),m_colon_2(m.m_colon_2),m_d(m.m_d)
{};

template<class T>
inline sub_band_matrix_1<T>::sub_band_matrix_1(matrix_type* m, Integer r)
    :m_matrix(m), m_ind_1(r)
{}

template<class T>
inline sub_band_matrix_1<T>::sub_band_matrix_1(const sub_band_matrix_1& m)
    :m_matrix(m.m_matrix),m_ind_1(m.m_ind_1)
{};

template<class T>
inline sub_band_matrix_1<T>::sub_band_matrix_1(sub_band_matrix_1&& m)
    :m_matrix(m.m_matrix),m_ind_1(m.m_ind_1)
{};

template<class T>
inline sub_band_matrix_2<T>::sub_band_matrix_2(matrix_type* m, Integer r, Integer c)
    :m_matrix(m), m_ind_1(r), m_ind_2(c)
{}

template<class T>
inline sub_band_matrix_2<T>::sub_band_matrix_2(const sub_band_matrix_2& m)
    :m_matrix(m.m_matrix),m_ind_1(m.m_ind_1),m_ind_2(m.m_ind_2)
{};

template<class T>
inline sub_band_matrix_2<T>::sub_band_matrix_2(sub_band_matrix_2&& m)
    :m_matrix(m.m_matrix),m_ind_1(m.m_ind_1),m_ind_2(m.m_ind_2)
{};

};