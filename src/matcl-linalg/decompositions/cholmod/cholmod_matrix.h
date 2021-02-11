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

#ifndef _MSC_VER
#include <stdexcept>
#endif

#include "matcl-blas-lapack/blas/details/blas_utils.h"

namespace matcl { namespace details 
{

// which part of symmetric or hermitian matrix is stored
enum class triang_part
{
    UPPER, LOWER
};

// a simple representation of a symmetric or hermitian dense matrix
template<class T> class cholmod_matrix
{
    public:
        using data_type = T;
        using self_type = cholmod_matrix;

    private:
        Integer             m_M;
        Integer             m_N;
        Integer             m_LD;
        T*                  m_data;
        triang_part         m_tr_type;

    public:
        // create Matrix from array data representing a symmetrix or hermitian
        // dense matrix of size M x N with leading dimension LD; tr_type informs
        // whether upper- or lower-triangular part of the matrix is stored in
        // array data
        cholmod_matrix(T* data, Integer M, Integer N, Integer LD, triang_part tr_type)
            :m_data(data),m_M(M),m_N(N),m_LD(LD),m_tr_type(tr_type)
        {
            if ( M < 0 || N < 0 || LD < 1)
            {
                throw std::runtime_error("invalid matrix");
            };
        };

        // get element at row i and column j; no checks are performed
        T&                  operator()(Integer i,Integer j)         { return m_data[i+j*m_LD]; };
        const T&            operator()(Integer i,Integer j) const   { return m_data[i+j*m_LD]; };

        // number of rows of the matrix
        Integer             rows() const                            { return m_M; };

        // number of columns of the matrix
        Integer             cols() const                            { return m_N; };

        // leading dimension of the matrix
        Integer             ld() const                              { return m_LD; };

        // return which part of the matrix is stored
        triang_part         triang_type() const                     { return m_tr_type; };
};

};};