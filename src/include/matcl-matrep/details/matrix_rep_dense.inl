/*
 *  This file is a part of Matrix Computation Library (MATCL)
 *
 *  Copyright (c) Pawe� Kowal 2017 - 2018
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

#include "matcl-matrep/matrix/matrix_rep_dense.h"
#include "matcl-matrep/IO/matrix_io.h"

namespace matcl
{

//-----------------------------------------------------------------
//                  dense_matrix
//-----------------------------------------------------------------
template<class T>
inline const T& dense_matrix<T, true>::operator()(Integer r) const
{
    if (r < 1 || r > m_size)
        throw_error_single_index(r,m_size);

    return m_ptr[r-1];
};

template<class T>
inline const T& dense_matrix<T, true>::operator()(Integer r)
{
    if (r < 1 || r > m_size)
        throw_error_single_index(r,m_size);

    return m_ptr[r-1];
};

template<class T>
inline T& dense_matrix<T, false>::operator()(Integer r)
{
    if (r < 1 || r > m_size)
        throw_error_single_index(r,m_size);

    return m_ptr[r-1];
};

template<class T>
inline const T& dense_matrix<T, true>::operator()(Integer r, Integer c) const
{
    if ((r < 1) || (r > m_rows) || (c < 1) || (c > m_cols))
        throw_error_double_index(r, c, m_rows, m_cols);

    return m_ptr[r-1 + (c-1) * m_ld];
};

template<class T>
inline const T& dense_matrix<T, true>::operator()(Integer r, Integer c)
{
    if ((r < 1) || (r > m_rows) || (c < 1) || (c > m_cols))
        throw_error_double_index(r, c, m_rows, m_cols);

    return m_ptr[r-1 + (c-1) * m_ld];
};

template<class T>
inline T& dense_matrix<T, false>::operator()(Integer r, Integer c)
{
    if ((r < 1) || (r > m_rows) || (c < 1) || (c > m_cols))
        throw_error_double_index(r, c, m_rows, m_cols);

    return m_ptr[r-1 + (c-1) * m_ld];
};

//-----------------------------------------------------------------
//                  sub_dense_matrix
//-----------------------------------------------------------------
template<class T>
inline sub_dense_matrix<T>::sub_dense_matrix(matrix_type* m, const colon& c1)
    :m_matrix(m), m_colon_1(&c1),m_colon_2(nullptr),m_d(0)
{}

template<class T>
inline sub_dense_matrix<T>::sub_dense_matrix(matrix_type* m, const colon& c1, const colon& c2)
    :m_matrix(m), m_colon_1(&c1),m_colon_2(&c2),m_d(0)
{};

template<class T>
inline sub_dense_matrix<T>::sub_dense_matrix(Integer diag, matrix_type* m)
    :m_matrix(m), m_colon_1(nullptr), m_colon_2(nullptr), m_d(diag)
{};

template<class T>
inline sub_dense_matrix<T>::sub_dense_matrix(const sub_dense_matrix& m)
    :m_matrix(m.m_matrix),m_colon_1(m.m_colon_1),m_colon_2(m.m_colon_2),m_d(m.m_d)
{};

template<class T>
inline sub_dense_matrix<T>::sub_dense_matrix(sub_dense_matrix&& m)
    :m_matrix(m.m_matrix),m_colon_1(m.m_colon_1),m_colon_2(m.m_colon_2),m_d(m.m_d)
{};

//-----------------------------------------------------------------
//                  dense_row
//-----------------------------------------------------------------
template<class T>
inline dense_row<T>::dense_row()
{};

template<class T>
inline dense_row<T>::~dense_row()
{};

template<class T>
inline dense_row<T>::dense_row(const dense_row& other)
    :base_type(other)
{};

template<class T>
inline dense_row<T>::dense_row(dense_row&& other)
    :base_type(std::move(other))
{};

template<class T>
inline dense_row<T>& dense_row<T>::operator=(const dense_row& rhs)
{
    base_type::operator=(rhs);
    return *this;
};

template<class T>
inline dense_row<T>& dense_row<T>::operator=(dense_row&& rhs)
{
    base_type::operator=(std::move(rhs));
    return *this;
};

template<class T>
inline dense_row<T>& dense_row<T>::add(const T& rhs)
{
    base_type::add(rhs);
    return *this;
};

template<class T>
inline dense_row<T>& dense_row<T>::add(T&& rhs)
{
    base_type::add(std::move(rhs));
    return *this;
};

template<class T>
inline dense_row<T>& dense_row<T>::add(const dense_matrix<T>& rhs)
{
    base_type::add(Matrix(rhs));
    return *this;
};

template<class T>
inline dense_row<T>& dense_row<T>::add(dense_matrix<T>&& rhs)
{
    base_type::add(Matrix(std::move(rhs)));
    return *this;
};

template<class T>
inline dense_row<T>& dense_row<T>::add(const dense_row& rhs)
{
    base_type::add(rhs);
    return *this;
};

template<class T>
inline dense_row<T>& dense_row<T>::add(dense_row&& rhs)
{
    base_type::add(std::move(rhs));
    return *this;
};

template<class T>
inline dense_row<T>& dense_row<T>::add(const dense_col<T>& rhs)
{
    base_type::add(rhs);
    return *this;
};

template<class T>
inline dense_row<T>& dense_row<T>::add(dense_col<T>&& rhs)
{
    base_type::add(std::move(rhs));
    return *this;
};

template<class T>
inline const dense_matrix<T> dense_row<T>::to_matrix() const
{
    return dense_matrix<T>(base_type::to_matrix());
}

//-----------------------------------------------------------------
//                  dense_col
//-----------------------------------------------------------------
template<class T>
inline dense_col<T>::dense_col()
{};

template<class T>
inline dense_col<T>::dense_col(const dense_col& other)
    :base_type(other)
{};

template<class T>
inline dense_col<T>::dense_col(dense_col&& other)
    :base_type(std::move(other))
{};

template<class T>
inline dense_col<T>& dense_col<T>::operator=(const dense_col& rhs)
{
    base_type::operator=(rhs);
    return *this;
};

template<class T>
inline dense_col<T>& dense_col<T>::operator=(dense_col&& rhs)
{
    base_type::operator=(std::move(rhs));
    return *this;
};

template<class T>
inline dense_col<T>::~dense_col()
{};

template<class T>
inline dense_col<T>& dense_col<T>::add(const T& rhs)
{
    base_type::add(rhs);
    return *this;
};

template<class T>
inline dense_col<T>& dense_col<T>::add(T&& rhs)
{
    base_type::add(std::move(rhs));
    return *this;
};

template<class T>
inline dense_col<T>& dense_col<T>::add(const dense_matrix<T>& rhs)
{
    base_type::add(Matrix(rhs));
    return *this;
};

template<class T>
inline dense_col<T>& dense_col<T>::add(dense_matrix<T>&& rhs)
{
    base_type::add(Matrix(std::move(rhs)));
    return *this;
};

template<class T>
inline dense_col<T>& dense_col<T>::add(const dense_row<T>& rhs)
{
    base_type::add(rhs);
    return *this;
};

template<class T>
inline dense_col<T>& dense_col<T>::add(dense_row<T>&& rhs)
{
    base_type::add(std::move(rhs));
    return *this;
};

template<class T>
inline dense_col<T>& dense_col<T>::add(const dense_col& rhs)
{
    base_type::add(rhs);
    return *this;
};

template<class T>
inline dense_col<T>& dense_col<T>::add(dense_col&& rhs)
{
    base_type::add(std::move(rhs));
    return *this;
};

template<class T>
inline const dense_matrix<T> dense_col<T>::to_matrix() const
{
    return dense_matrix<T>(base_type::to_matrix());
};

};
