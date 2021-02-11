/*
 *  This file is a part of Matrix Computation Library (MATCL)
 *
 *  Copyright (c) Pawe³ Kowal 2017 - 2021
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

#include "matcl-internals/container/mat_s.h"
#include "matcl-internals/error/error_check_basic.h"

namespace matcl { namespace raw 
{

template <class value_type>
inline value_type sparse_matrix_base<value_type>::operator()(Integer i, Integer j) const
{
    error::check_index(i, j, m_data.rows(), m_data.cols());
    return m_data.get_elem(i-1, j-1);
}

template <class value_type>
inline Integer sparse_matrix_base<value_type>::rows() const
{
    return m_data.rows();
}

template <class value_type>
inline Integer sparse_matrix_base<value_type>::length() const
{ 
    return (rows() == 0 || cols() == 0)? 0 : std::max(rows(), cols()); 
};

template <class value_type>
inline Integer sparse_matrix_base<value_type>::cols() const
{
    return m_data.cols();
}

template <class value_type>
inline Integer sparse_matrix_base<value_type>::max_cols() const
{
    return m_data.max_cols();
}

template <class value_type>
inline Integer sparse_matrix_base<value_type>::total_cols() const
{
    return m_data.total_cols();
}

template <class value_type>
inline Integer sparse_matrix_base<value_type>::nnz() const
{
    return m_data.nnz();
}

template <class value_type>
inline Integer sparse_matrix_base<value_type>::nzmax() const
{
    return m_data.nzmax();
}

template <class value_type>
inline sparse_matrix_base<value_type>
sparse_matrix_base<value_type>::make_unique(bool keep_maxcol) const
{
    return sparse_matrix_base(m_data.make_unique(keep_maxcol));
};

template <class value_type>
template<class Real_T>
sparse_matrix_base<value_type>::sparse_matrix_base(type_info ti, const Integer* ri, const Integer* ci, 
                    const Real_T* xr, const Real_T* xi, Integer r, Integer c, Integer nnz,
                    typename matcl::details::enable_if_val_complex<value_type,Real_T>::type)
:m_data(ti)
{
    construct2(ri, ci, xr,xi, r, c, nnz, nnz);
};

template <class value_type>
template<class Real_T>
sparse_matrix_base<value_type>::sparse_matrix_base(type_info ti, const Integer* ri, const Integer* ci, 
                    const Real_T* xr, const Real_T* xi, 
                    Integer r, Integer c, Integer nnz, Integer nzmax,
                    typename matcl::details::enable_if_val_complex<value_type,Real_T>::type)
:m_data(ti)
{
    construct2(ri, ci, xr,xi, r, c, nnz, nzmax);
};

template<class value_type>
inline Matrix<value_type,struct_sparse>
Matrix<value_type,struct_sparse>::copy(bool keep_maxcol) const
{
    return base_type::copy(keep_maxcol);
}

template<class value_type>
inline Matrix<value_type,struct_sparse>
Matrix<value_type,struct_sparse>::clone(bool keep_maxcol) const
{
    return base_type::clone(keep_maxcol);
};
        
template<class value_type>
inline Matrix<value_type,struct_sparse>
Matrix<value_type,struct_sparse>::make_unique(bool keep_bufor) const
{
    return base_type::make_unique(keep_bufor);
}
                
template<class value_type>
inline Matrix<value_type,struct_sparse>
Matrix<value_type,struct_sparse>::resize(Integer r, Integer c)
{
    return base_type::resize(r, c);
}

template<class value_type>
inline const Matrix<value_type,struct_sparse> 
Matrix<value_type,struct_sparse>::reserve(Integer r, Integer c) const
{
    return base_type::reserve(r, c);
};

template<class value_type>
inline const Matrix<value_type,struct_sparse> 
Matrix<value_type,struct_sparse>::resize(Integer r, Integer c) const
{
    return base_type::resize(r, c);
};

template<class value_type>
inline const Matrix<value_type,struct_sparse> 
Matrix<value_type,struct_sparse>::make_view(Integer c_start, Integer c_end) const
{
    return base_type::make_view(c_start, c_end);
}

};};