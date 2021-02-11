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

#include "matcl-matrep/matrix/unique_matrix.h"
#include "matcl-matrep/details/matrix_details_subs.h"

namespace matcl
{

unique_matrix::unique_matrix(Matrix& mat, bool allow_nondense)
{
    m_data.m_matrix = &mat;
    m_type          = type::mat;

    check_dense(m_data.m_matrix, allow_nondense);
    m_data.m_matrix->mark_unique(true);
};

unique_matrix::unique_matrix(matcl::sub_matrix&& sm, bool allow_nondense)
{
    m_data.m_sub    = &sm;
    m_type          = type::sub;

    check_dense(sm.m_matrix, allow_nondense);
};

unique_matrix::unique_matrix(matcl::sub_matrix_1&& sm, bool allow_nondense)
{
    m_data.m_sub_1  = &sm;
    m_type          = type::sub_1;

    check_dense(sm.m_matrix, allow_nondense);
};

unique_matrix::unique_matrix(matcl::sub_matrix_2&& sm, bool allow_nondense)
{
    m_data.m_sub_2  = &sm;
    m_type          = type::sub_2;

    check_dense(sm.m_matrix, allow_nondense);
};

unique_matrix::~unique_matrix()
{
    if (m_type == type::mat)
        m_data.m_matrix->mark_unique(false);

    m_conv_sub.mark_unique(false);
};

Matrix unique_matrix::to_matrix() &&
{
    Matrix ret;
    switch(m_type)
    {
        case type::mat:
            return *m_data.m_matrix;
        case type::sub:
            ret = *m_data.m_sub;
            break;
        case type::sub_1:
            ret = *m_data.m_sub_1;
            break;
        case type::sub_2:
            ret = *m_data.m_sub_2;
            break;
        default:
            //invalid case
            matcl_assert(0,"invalid case");
            return Matrix();
    }

    m_conv_sub.mark_unique(false);
    m_conv_sub = ret;
    m_conv_sub.mark_unique(true);
    return ret;
};

Matrix& unique_matrix::operator=(const matcl::sub_matrix& sm) &&
{
    switch(m_type)
    {
        case type::mat:
        {
            m_data.m_matrix->mark_unique(false);
            return (*m_data.m_matrix = sm);
        }
        case type::sub:
            return (std::move(*m_data.m_sub).assign_unique(sm));
        case type::sub_1:
            return (std::move(*m_data.m_sub_1).assign_unique(sm));
        case type::sub_2:
            return (std::move(*m_data.m_sub_2).assign_unique(sm));
        default:
            //invalid case
            matcl_assert(0,"invalid case");
            throw error::error_general("invalid case");
    }
};

Matrix& unique_matrix::operator=(const matcl::sub_matrix_1& sm) &&
{
    switch(m_type)
    {
        case type::mat:
        {
            m_data.m_matrix->mark_unique(false);
            return (*m_data.m_matrix = sm);
        }
        case type::sub:
            return (std::move(*m_data.m_sub).assign_unique(sm));
        case type::sub_1:
            return (std::move(*m_data.m_sub_1).assign_unique(sm));
        case type::sub_2:
            return (std::move(*m_data.m_sub_2).assign_unique(sm));
        default:
            //invalid case
            matcl_assert(0,"invalid case");
            throw error::error_general("invalid case");
    }
};

Matrix& unique_matrix::operator=(const matcl::sub_matrix_2& sm) &&
{
    switch(m_type)
    {
        case type::mat:
        {
            m_data.m_matrix->mark_unique(false);
            return (*m_data.m_matrix = sm);
        }
        case type::sub:
            return (std::move(*m_data.m_sub).assign_unique(sm));
        case type::sub_1:
            return (std::move(*m_data.m_sub_1).assign_unique(sm));
        case type::sub_2:
            return (std::move(*m_data.m_sub_2).assign_unique(sm));
        default:
            //invalid case
            matcl_assert(0,"invalid case");
            throw error::error_general("invalid case");
    }
};

Matrix& unique_matrix::operator=(const Matrix& sm) &&
{
    switch(m_type)
    {
        case type::mat:
        {
            m_data.m_matrix->mark_unique(false);
            return (*m_data.m_matrix = sm);
        }
        case type::sub:
            return (std::move(*m_data.m_sub).assign_unique(sm));
        case type::sub_1:
            return (std::move(*m_data.m_sub_1).assign_unique(sm));
        case type::sub_2:
            return (std::move(*m_data.m_sub_2).assign_unique(sm));
        default:
            //invalid case
            matcl_assert(0,"invalid case");
            throw error::error_general("invalid case");
    }
};

Matrix& unique_matrix::operator=(Matrix&& sm) &&
{
    switch(m_type)
    {
        case type::mat:
        {
            m_data.m_matrix->mark_unique(false);
            return (*m_data.m_matrix = std::move(sm));
        }
        case type::sub:
            return std::move(*m_data.m_sub).assign_unique(std::move(sm));
        case type::sub_1:
            return std::move(*m_data.m_sub_1).assign_unique(std::move(sm));
        case type::sub_2:
            return std::move(*m_data.m_sub_2).assign_unique(std::move(sm));
        default:
            //invalid case
            matcl_assert(0,"invalid case");
            throw error::error_general("invalid case");
    }
};

Matrix& unique_matrix::operator=(unique_matrix&& sm) &&
{
    switch(m_type)
    {
        case type::mat:
        {
            m_data.m_matrix->mark_unique(false);
            return (*m_data.m_matrix = std::move(sm));
        }
        case type::sub:
            return (std::move(*m_data.m_sub).assign_unique(std::move(sm)));
        case type::sub_1:
            return (std::move(*m_data.m_sub_1).assign_unique(std::move(sm)));
        case type::sub_2:
            return (std::move(*m_data.m_sub_2).assign_unique(std::move(sm)));
        default:
            //invalid case
            matcl_assert(0,"invalid case");
            throw error::error_general("invalid case");
    }
}

void unique_matrix::check_dense(const Matrix* mat, bool allow_nondense) const
{
    if (allow_nondense == false && mat->get_struct_code() != struct_code::struct_dense)
        throw error::dense_matrix_required(mat->get_struct_code());
};

bool unique_matrix::is_matrix() const
{
    return m_type == type::mat;
}

bool unique_matrix::is_submatrix() const
{
    return m_type == type::sub;
}

bool unique_matrix::is_submatrix_1() const
{
    return m_type == type::sub_1;
}

bool unique_matrix::is_submatrix_2() const
{
    return m_type == type::sub_2;
}

Matrix& unique_matrix::get_matrix() &&
{
    if (is_matrix() == false)
        throw error::invalid_extract_from_unique_matrix();

    return *m_data.m_matrix;
}

sub_matrix&& unique_matrix::get_submatrix() &&
{
    if (is_submatrix() == false)
        throw error::invalid_extract_from_unique_matrix();

    return std::move(*m_data.m_sub);
}

sub_matrix_1&& unique_matrix::get_submatrix_1() &&
{
    if (is_submatrix_1() == false)
        throw error::invalid_extract_from_unique_matrix();

    return std::move(*m_data.m_sub_1);
}

sub_matrix_2&& unique_matrix::get_submatrix_2() &&
{
    if (is_submatrix_2() == false)
        throw error::invalid_extract_from_unique_matrix();

    return std::move(*m_data.m_sub_2);
}

};
