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

#include "matcl-matrep/matrix/colon.h"
#include "matcl-matrep/details/fwd_decls.h"
#include "matcl-matrep/matrix/matrix_concat.h"
#include "matcl-matrep/details/matrix.inl"
#include "matcl-matrep/lib_functions/manip.h"
#include "matcl-internals/container/mat_d.h"

namespace matcl
{

colon::colon(const Matrix& mat)
{
    create(mat);
};

colon::colon(Matrix&& mat)
{
    create(std::move(mat));
};

colon::colon(const permvec& mat)
{
    create(mat.to_matrix());
};

colon::colon(permvec&& mat)
{
    create(mat.to_matrix());
};

void colon::create(Matrix&& mat)
{
    Matrix tmp(std::move(mat));
    create(tmp);
};

void colon::create(const Matrix& mat)
{
    m_s		= 0;
    m_i		= 0;
    m_e		= 0;
    m_flag	= t_matrix1;
    m_mat_size.m_mat_12 = nullptr;

    Integer mr  = mat.rows();
    Integer mc  = mat.cols();

    const raw::integer_dense& imat = mat.impl<raw::integer_dense>();

    if (mr == 0 || mc == 0)
    {
        m_mat_size.m_mat_12 = new Matrix(imat.make_explicit(),false);
        return;
    }    
    
    if (mr == 1 && mc == 1)
    {
        Integer p1  = *imat.ptr();
        m_s         = p1;
        m_i         = 1;
        m_e         = p1;
        m_flag      = t_range_simple;
        return;
    };

    const Integer* ptr = imat.ptr();

    Integer s   = ptr[0];
    Integer l   = ptr[1];
    Integer i   = l - s;

    if (i == 0)
    {
        //conversion to raw colon is not possible
        m_mat_size.m_mat_12 = new Matrix(imat.make_explicit(),false);    
        return;
    };

    for(Integer k = 2; k < imat.size(); ++k)
    {
        Integer c   = ptr[k];
        Integer di  = c - l;
        l           = c;

        if (di != i)
        {
            //conversion to raw colon is not possible
            m_mat_size.m_mat_12 = new Matrix(imat.make_explicit(),false);
            return;
        };
    };

    m_s         = s;
    m_i         = i;
    m_e         = l;
    m_flag      = t_range_mat;

    matrix_size_type* tmp   = new matrix_size_type{imat.rows(),imat.cols()};
    m_mat_size.m_mat_size   = tmp;
};

void colon::create(const Matrix& mat_1, const Matrix& mat_2)
{
    m_s		= 0;
    m_i		= 0;
    m_e		= 0;
    m_flag	= t_matrix2;
    m_mat_size.m_mat_12 = nullptr;

    if (mat_2.rows() != mat_1.rows() || mat_2.cols() != mat_1.cols())
        throw error::invalid_size2(mat_2.rows(), mat_2.cols(), mat_1.rows(), mat_2.cols());

    const raw::integer_dense& imat1 = mat_1.impl<raw::integer_dense>();
    const raw::integer_dense& imat2 = mat_2.impl<raw::integer_dense>();

    m_mat_size.m_mat_12     = new Matrix[2];
    m_mat_size.m_mat_12[0]  = Matrix(imat1.make_explicit(),false);
    m_mat_size.m_mat_12[1]  = Matrix(imat2.make_explicit(),false);
    return;
};

Integer colon::get_matrix_rows() const
{
    return m_mat_size.m_mat_size->rows;
};

Integer colon::get_matrix_cols() const
{
    return m_mat_size.m_mat_size->cols;
};

colon::~colon()
{
    if (m_flag == t_matrix1)
    {
        if (m_mat_size.m_mat_12)
            delete m_mat_size.m_mat_12;
    }
    else if (m_flag == t_matrix2)
    {
        if (m_mat_size.m_mat_12)
            delete[] m_mat_size.m_mat_12;
    }
    else
    {
        if(m_mat_size.m_mat_size)
            delete m_mat_size.m_mat_size;
    }
};

colon::colon(colon&& other)
    : m_s(other.m_s), m_i(other.m_i), m_e(other.m_e), m_flag(other.m_flag)
    , m_mat_size(other.m_mat_size)
{
    other.m_mat_size.m_mat_12 = nullptr;
};

colon colon::copy() const
{
    colon ret;

    ret.m_s     = this->m_s;
    ret.m_i     = this->m_i;
    ret.m_e     = this->m_e;
    ret.m_flag  = this->m_flag;

    if (this->m_mat_size.m_mat_12 != nullptr)
    {
        if (this->m_flag == t_matrix1)
        {
            ret.m_mat_size.m_mat_12 = new Matrix(*this->m_mat_size.m_mat_12);
        }
        else if (this->m_flag == t_matrix2)
        {
            ret.m_mat_size.m_mat_12     = new Matrix[2];
            ret.m_mat_size.m_mat_12[0]  = this->m_mat_size.m_mat_12[0];
            ret.m_mat_size.m_mat_12[1]  = this->m_mat_size.m_mat_12[1];
        }
        else
        {
            ret.m_mat_size.m_mat_size = new matrix_size_type(*this->m_mat_size.m_mat_size);
        }
    };

    return ret;
};

colon& colon::operator=(colon&& other)
{
    m_s     = other.m_s;
    m_i     = other.m_i;
    m_e     = other.m_e;
    m_flag  = other.m_flag;

    std::swap(m_mat_size,other.m_mat_size);

    return *this;
};

};
