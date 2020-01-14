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

namespace matcl { namespace raw 
{

template<class Raw_matrix>
class const_matrix
{
    public:
        using matrix_type   = Raw_matrix;

    private:
        matrix_type         m_matrix;

    public:
        const_matrix(const matrix_type& m)
            :m_matrix(m, Raw_matrix::copy_is_safe())
        {};

        ~const_matrix()
        {};

        const_matrix(const const_matrix& m)
            :m_matrix(m.m_matrix, Raw_matrix::copy_is_safe())
        {};

        const_matrix(const_matrix&& m)
            :m_matrix(std::move(m.m_matrix))
        {};

        const_matrix& operator=(const const_matrix& m);
        const_matrix& operator=(const_matrix&& m);

        const Raw_matrix&   get() const;

        void rebind(const matrix_type& m);
        void rebind(matrix_type&& m);
};

template<class Raw_matrix>
inline
const_matrix<Raw_matrix>& const_matrix<Raw_matrix>::operator=(const const_matrix& m)
{
    this->rebind(m.m_matrix);
    return *this;
};

template<class Raw_matrix>
inline
const_matrix<Raw_matrix>& const_matrix<Raw_matrix>::operator=(const_matrix&& m)
{
    this->rebind(std::move(m.m_matrix));
    return *this;
};

template<class Raw_matrix>
inline
const Raw_matrix& const_matrix<Raw_matrix>::get() const
{
    return m_matrix;
}

template<class Raw_matrix>
inline
void const_matrix<Raw_matrix>::rebind(const matrix_type& m)
{
    // m_matrix is not stored in a Matrix container
    m_matrix.assign_to_fresh(m);
};

template<class Raw_matrix>
inline
void const_matrix<Raw_matrix>::rebind(matrix_type&& m)
{
    // m_matrix is not stored in a Matrix container
    m_matrix.assign_to_fresh(std::move(m));
};

};};

