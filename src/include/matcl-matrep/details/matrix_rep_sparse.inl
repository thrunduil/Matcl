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

#include "matcl-matrep/matrix/matrix_rep_sparse.h"

namespace matcl
{

//-----------------------------------------------------------------
//                  sub_sparse_matrix
//-----------------------------------------------------------------
template<class T>
inline sub_sparse_matrix<T>::sub_sparse_matrix(matrix_type* m, const colon& c1)
    :m_matrix(m), m_colon_1(&c1),m_colon_2(nullptr),m_d(0)
{}

template<class T>
inline sub_sparse_matrix<T>::sub_sparse_matrix(matrix_type* m, const colon& c1, const colon& c2)
    :m_matrix(m), m_colon_1(&c1),m_colon_2(&c2),m_d(0)
{};

template<class T>
inline sub_sparse_matrix<T>::sub_sparse_matrix(Integer diag, matrix_type* m)
    :m_matrix(m), m_colon_1(nullptr), m_colon_2(nullptr), m_d(diag)
{};


template<class T>
inline sub_sparse_matrix<T>::sub_sparse_matrix(const sub_sparse_matrix& m)
    :m_matrix(m.m_matrix),m_colon_1(m.m_colon_1),m_colon_2(m.m_colon_2),m_d(m.m_d)
{};

template<class T>
inline sub_sparse_matrix<T>::sub_sparse_matrix(sub_sparse_matrix&& m)
    :m_matrix(m.m_matrix),m_colon_1(m.m_colon_1),m_colon_2(m.m_colon_2),m_d(m.m_d)
{};

template<class T>
inline sub_sparse_matrix_1<T>::sub_sparse_matrix_1(matrix_type* m, Integer r)
    :m_matrix(m), m_ind_1(r)
{}

template<class T>
inline sub_sparse_matrix_1<T>::sub_sparse_matrix_1(const sub_sparse_matrix_1& m)
    :m_matrix(m.m_matrix),m_ind_1(m.m_ind_1)
{};

template<class T>
inline sub_sparse_matrix_1<T>::sub_sparse_matrix_1(sub_sparse_matrix_1&& m)
    :m_matrix(m.m_matrix),m_ind_1(m.m_ind_1)
{};

template<class T>
inline sub_sparse_matrix_2<T>::sub_sparse_matrix_2(matrix_type* m, Integer r, Integer c)
    :m_matrix(m), m_ind_1(r), m_ind_2(c)
{}

template<class T>
inline sub_sparse_matrix_2<T>::sub_sparse_matrix_2(const sub_sparse_matrix_2& m)
    :m_matrix(m.m_matrix),m_ind_1(m.m_ind_1),m_ind_2(m.m_ind_2)
{};

template<class T>
inline sub_sparse_matrix_2<T>::sub_sparse_matrix_2(sub_sparse_matrix_2&& m)
    :m_matrix(m.m_matrix),m_ind_1(m.m_ind_1),m_ind_2(m.m_ind_2)
{};

//-----------------------------------------------------------------
//                  sparse_row
//-----------------------------------------------------------------
template<class T>
inline sparse_row<T>::sparse_row()
{};

template<class T>
inline sparse_row<T>::~sparse_row()
{};

template<class T>
inline sparse_row<T>::sparse_row(const sparse_row& other)
    :base_type(other)
{};

template<class T>
inline sparse_row<T>::sparse_row(sparse_row&& other)
    :base_type(std::move(other))
{};

template<class T>
inline sparse_row<T>& sparse_row<T>::operator=(const sparse_row& rhs)
{
    base_type::operator=(rhs);
    return *this;
};

template<class T>
inline sparse_row<T>& sparse_row<T>::operator=(sparse_row&& rhs)
{
    base_type::operator=(std::move(rhs));
    return *this;
};

template<class T>
inline sparse_row<T>& sparse_row<T>::add(const T& rhs)
{
    base_type::add(rhs);
    return *this;
};

template<class T>
inline sparse_row<T>& sparse_row<T>::add(T&& rhs)
{
    base_type::add(std::move(rhs));
    return *this;
};

template<class T>
inline sparse_row<T>& sparse_row<T>::add(const sparse_matrix<T>& rhs)
{
    base_type::add(Matrix(rhs));
    return *this;
};

template<class T>
inline sparse_row<T>& sparse_row<T>::add(sparse_matrix<T>&& rhs)
{
    base_type::add(Matrix(std::move(rhs)));
    return *this;
};

template<class T>
inline sparse_row<T>& sparse_row<T>::add(const band_matrix<T>& rhs)
{
    base_type::add(Matrix(rhs));
    return *this;
};

template<class T>
inline sparse_row<T>& sparse_row<T>::add(band_matrix<T>&& rhs)
{
    base_type::add(Matrix(std::move(rhs)));
    return *this;
};

template<class T>
inline sparse_row<T>& sparse_row<T>::add(const sparse_row& rhs)
{
    base_type::add(rhs);
    return *this;
};

template<class T>
inline sparse_row<T>& sparse_row<T>::add(sparse_row&& rhs)
{
    base_type::add(std::move(rhs));
    return *this;
};

template<class T>
inline sparse_row<T>& sparse_row<T>::add(const sparse_col<T>& rhs)
{
    base_type::add(rhs);
    return *this;
};

template<class T>
inline sparse_row<T>& sparse_row<T>::add(sparse_col<T>&& rhs)
{
    base_type::add(std::move(rhs));
    return *this;
};

template<class T>
inline const sparse_matrix<T> sparse_row<T>::to_matrix() const
{
    return sparse_matrix<T>(base_type::to_matrix());
}

//-----------------------------------------------------------------
//                  sparse_col
//-----------------------------------------------------------------
template<class T>
inline sparse_col<T>::sparse_col()
{};

template<class T>
inline sparse_col<T>::sparse_col(const sparse_col& other)
    :base_type(other)
{};

template<class T>
inline sparse_col<T>::sparse_col(sparse_col&& other)
    :base_type(std::move(other))
{};

template<class T>
inline sparse_col<T>& sparse_col<T>::operator=(const sparse_col& rhs)
{
    base_type::operator=(rhs);
    return *this;
};

template<class T>
inline sparse_col<T>& sparse_col<T>::operator=(sparse_col&& rhs)
{
    base_type::operator=(std::move(rhs));
    return *this;
};

template<class T>
inline sparse_col<T>::~sparse_col()
{};

template<class T>
inline sparse_col<T>& sparse_col<T>::add(const T& rhs)
{
    base_type::add(rhs);
    return *this;
};

template<class T>
inline sparse_col<T>& sparse_col<T>::add(T&& rhs)
{
    base_type::add(std::move(rhs));
    return *this;
};

template<class T>
inline sparse_col<T>& sparse_col<T>::add(const sparse_matrix<T>& rhs)
{
    base_type::add(Matrix(rhs));
    return *this;
};

template<class T>
inline sparse_col<T>& sparse_col<T>::add(sparse_matrix<T>&& rhs)
{
    base_type::add(Matrix(std::move(rhs)));
    return *this;
};

template<class T>
inline sparse_col<T>& sparse_col<T>::add(const band_matrix<T>& rhs)
{
    base_type::add(Matrix(rhs));
    return *this;
};

template<class T>
inline sparse_col<T>& sparse_col<T>::add(band_matrix<T>&& rhs)
{
    base_type::add(Matrix(std::move(rhs)));
    return *this;
};

template<class T>
inline sparse_col<T>& sparse_col<T>::add(const sparse_row<T>& rhs)
{
    base_type::add(Matrix(rhs));
    return *this;
};

template<class T>
inline sparse_col<T>& sparse_col<T>::add(sparse_row<T>&& rhs)
{
    base_type::add(std::move(rhs));
    return *this;
};

template<class T>
inline sparse_col<T>& sparse_col<T>::add(const sparse_col& rhs)
{
    base_type::add(rhs);
    return *this;
};

template<class T>
inline sparse_col<T>& sparse_col<T>::add(sparse_col&& rhs)
{
    base_type::add(std::move(rhs));
    return *this;
};

template<class T>
inline const sparse_matrix<T> sparse_col<T>::to_matrix() const
{
    return sparse_matrix<T>(base_type::to_matrix());
};

//-----------------------------------------------------------------
//                  sparse_matrix
//-----------------------------------------------------------------

template<class T>
inline const Integer* sparse_matrix<T, true>::ptr_c() const
{ 
    return m_c; 
}

template<class T>
inline const Integer* sparse_matrix<T, true>::ptr_r() const
{ 
    return *m_r; 
}

template<class T>
inline const typename sparse_matrix<T, true>::value_type*
sparse_matrix<T, true>::ptr_x() const
{ 
    return *m_x; 
}

template<class T>
inline Integer sparse_matrix<T, true>::offset() const
{ 
    return m_offset; 
};

template<class T>
inline Integer sparse_matrix<T, true>::max_cols() const
{ 
    return m_max_cols; 
}

template<class T>
inline Integer sparse_matrix<T, true>::nnz() const
{ 
    return m_c? m_c[m_cols]-m_offset : 0; 
}

template<class T>
inline const Integer* sparse_matrix<T, true>::ptr_c()
{ 
    return m_c; 
}

template<class T>
inline const Integer* sparse_matrix<T, true>::ptr_r()
{ 
    return *m_r; 
}

template<class T>
inline const typename sparse_matrix<T, true>::value_type*
sparse_matrix<T, true>::ptr_x()
{ 
    return *m_x; 
}

template<class T>
inline const Matrix& sparse_matrix<T, true>::to_matrix() const &
{ 
    return m_matrix; 
};

template<class T>
inline Matrix&& sparse_matrix<T, true>::to_matrix() &&
{ 
    return std::move(m_matrix); 
};

template<class T>
inline Integer sparse_matrix<T, true>::rows() const
{ 
    return m_rows; 
}; 

template<class T>
inline Integer sparse_matrix<T, true>::cols() const
{ 
    return m_cols; 
};

template<class T>
inline const struct_flag sparse_matrix<T, true>::get_struct() const
{ 
    return *m_flag; 
};   

template<class T>
inline sparse_matrix<T, true>::operator bool() const
{ 
    return (bool)m_matrix; 
};

template<class T>
inline bool sparse_matrix<T, true>::is_empty() const
{ 
    return rows() == 0 || cols() == 0; 
};

template<class T>
inline bool sparse_matrix<T, true>::is_scalar() const
{ 
    return rows() == 1 && cols() == 1; 
};

template<class T>
inline bool sparse_matrix<T, true>::is_square() const
{ 
    return rows() == cols(); 
};

template<class T>
inline bool sparse_matrix<T, true>::is_vector() const
{ 
    return rows() == 1 || cols() == 1 || is_empty() == true; 
};

template<class T>
inline bool sparse_matrix<T, true>::is_matrix_type() const
{ 
    return true; 
}

template<class T>
inline bool sparse_matrix<T, true>::is_scalar_type() const
{ 
    return false; 
};

template<class T>
inline bool sparse_matrix<T, true>::is_unique() const
{ 
    return m_matrix.is_unique(); 
}

template<class T>
inline value_code sparse_matrix<T, true>::get_value_code() const
{ 
    return m_matrix.get_value_code(); 
}

template<class T>
inline struct_code sparse_matrix<T, true>::get_struct_code() const
{ 
    return m_matrix.get_struct_code(); 
}

template<class T>
inline mat_code sparse_matrix<T, true>::get_matrix_code() const
{ 
    return m_matrix.get_matrix_code(); 
}

template<class T>
inline Integer* sparse_matrix<T,false>::ptr_c()
{ 
    return m_c; 
}

template<class T>
inline Integer* sparse_matrix<T,false>::ptr_r()
{ 
    return *m_r; 
}

template<class T>
inline typename sparse_matrix<T,false>::value_type*
sparse_matrix<T,false>::ptr_x()
{ 
    return *m_x; 
}

};
