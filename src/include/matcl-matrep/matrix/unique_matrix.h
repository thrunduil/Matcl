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

#include "matcl-matrep/matrix/matrix.h"

//TODO: unique matrix for typed matrices

namespace matcl
{

// mark given matrix or submatrix as effectively unique, i.e. it is safe to modify
// given matrix or submatrix even if is not unique; for example:
//     1. assignment:     Y(sel_1)    = Y(sel_2)
//     2. gemm function:  gemm(alpha, Y(sel_2), Y(sel_3), t1, t2, beta, Y(sel_1))
//     3. inplace func:   func(Y)
// where sel_1, sel_2, sel_3 are appropriate colon operations, t1, t2 are trans_type;
// only unique matrices can be modified, therefore in case (1, 2) copy will always take
// place and in case 3 inplace version will be used only if Y is unique and temporaty
//
// unique_matrix class allows for relaxing this assumption by calling:
//     1. unique_matrix(Y(sel_1))  = Y(sel_2)
//     2. gemm(alpha, Y(sel_2), Y(sel_3), t1, t2, beta, unique_matrix(Y(sel_1)))
//     3. func(unique_matrix(Y))
// then Y(sel_1) will be assumed to be effectively unique and no copy will take place
//
// note that this is unsafe assumption and make break down the code if modification of
// Y(sel_1) will change other matrices in the expression as well or change matrices
// used in other threads
//
// unique_matrix class is short lived and should be used only in expressions;
// uniqueness assumption affects only given submatrix, other submatrices of given matrix
// as well as the whole matrix are not assumed to be unique; this uniqueness assumption
// is removed when destructor of unique_matrix is called, 
class MATCL_MATREP_EXPORT unique_matrix
{
    private:
        enum class type
        {
            mat, sub, sub_1, sub_2
        };

        union data
        {
            Matrix*         m_matrix;
            sub_matrix*     m_sub;
            sub_matrix_1*   m_sub_1;
            sub_matrix_2*   m_sub_2;
        };

        data                m_data;
        type                m_type;
        Matrix              m_conv_sub;

    public:
        // mark given matrix mat as effectively unique, this is especially
        // dangerous assumption if mat is not dense; on default exception is
        // thrown if mat is no dense unless allow_nondense is set to true
        // notice that mat can be a submatrix of other matrix if view was 
        // created
        unique_matrix(Matrix& mat, bool allow_nondense = false);

        // mark given submatrix submat as effectively unique, see also 
        // constructor from a Matrix
        unique_matrix(matcl::sub_matrix&& submat, bool allow_nondense = false);
        unique_matrix(matcl::sub_matrix_1&& submat, bool allow_nondense = false);
        unique_matrix(matcl::sub_matrix_2&& submat, bool allow_nondense = false);

        // destructor, remove uniqueness assumption
        ~unique_matrix();

        // assignments from a Matrix; call assignment operator defined in Matrix 
        // or one of submatrix classes
        Matrix&                 operator=(const Matrix& mat) &&;
        Matrix&                 operator=(Matrix&& mat) &&;

        // assignments from a Matrix, umat is converted to a Matrix first; 
        // call assignment operator defined in Matrix or one of submatrix classes
        Matrix&                 operator=(unique_matrix&& umat) &&;

        // assignments from submatrices; submatrices are converted to matrices
        // and then standard assignment take place defined in Matrix or one of 
        // submatrix classes
        Matrix&                 operator=(const matcl::sub_matrix&) &&;
        Matrix&                 operator=(const matcl::sub_matrix_1&) &&;
        Matrix&                 operator=(const matcl::sub_matrix_2&) &&;

        // return true if unique_matrix was contructored from a matrix
        bool                    is_matrix() const;

        // return true if unique_matrix was contructored from a sub_matrix
        bool                    is_submatrix() const;

        // return true if unique_matrix was contructored from a sub_matrix_1
        bool                    is_submatrix_1() const;

        // return true if unique_matrix was contructored from a sub_matrix_2
        bool                    is_submatrix_2() const;

        // return matrix used to construct this unique_matrix, throw exception
        // if is_matrix() == false
        Matrix&                 get_matrix() &&;

        // return sub_matrix used to construct this unique_matrix, throw exception
        // if is_submatrix() == false
        sub_matrix&&            get_submatrix() &&;

        // return sub_matrix_1 used to construct this unique_matrix, throw exception
        // if is_submatrix_1() == false
        sub_matrix_1&&          get_submatrix_1() &&;

        // return sub_matrix_2 used to construct this unique_matrix, throw exception
        // if is_submatrix_2() == false
        sub_matrix_2&&          get_submatrix_2() &&;

    private:
        unique_matrix() = delete;
        unique_matrix(const unique_matrix& m) = delete;
        unique_matrix(unique_matrix&& m) = delete;

        Matrix                  to_matrix() &&;
        void                    check_dense(const Matrix* mat, bool allow_nondense) const;

        friend Matrix;
};

};