/*
 *  This file is a part of Morfa Matrix Lib.
 *
 *  Copyright (c) Pawe³ Kowal 2011-2016
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

#if 0
TODO

#include "matcl-linalg/petsc_utils/petsc_matrix.h"
#include "matcl-matrep/matcl_matrep.h"
#include "matcl-linalg/general/mmlib_petsc_exception.h"

namespace matcl { namespace details
{

//---------------------------------------------------------------
//                  petsc_matrix_str
//---------------------------------------------------------------
petsc_matrix_str<struct_dense>::petsc_matrix_str(const Matrix& A_in)
    : petsc_matrix(A_in), Araw(ti::ti_empty()) 
{
    if (A_in.get_matrix_code() != mat_code::real_dense)
    {
        Matrix Ac   = convert(A_in, mat_code::real_dense);
        Araw.assign_to_fresh(Ac.get_impl_unique<raw::real_dense>());    
    }
    else
    {
        Araw.assign_to_fresh(A_in.get_impl<raw::real_dense>().make_unique());
    };

    m_Aksp = mp::smart_mat::mat_create_seq_dense(PETSC_COMM_SELF, Araw.rows(), Araw.cols(), Araw.ld(), 
                                                Araw.ptr());

    m_Aksp.set_struct_flags(A_in);
}

Matrix petsc_matrix_str<struct_dense>::mmul_right(const Matrix& X, trans_type tA) const
{
    return mmul(Matrix(Araw,false), X, tA);
}

petsc_matrix_str<struct_sparse>::petsc_matrix_str(const Matrix& A_in)
    : petsc_matrix(A_in), m_A(A_in), m_At(ti::ti_empty()), m_column_indices(ti::ti_empty())
{
    // NOTE: Petsc uses CSR not CSC
    // we need to transpose the matrix

    Matrix Ac   = convert(trans(A_in), mat_code::real_sparse);
    Mat Ac_rep  = Ac.get_impl_unique<raw::real_sparse>();

    m_At.assign_to_fresh(Ac_rep);

    Mat::sparse_ccs& sparse_ccs = m_At.rep();    

    if (sparse_ccs.offset() == 0)
    {            
        m_Aksp  = mp::smart_mat::mat_create_seq_aij_with_arrays(PETSC_COMM_SELF, A_in.rows(), A_in.cols(),
                        sparse_ccs.ptr_c(), sparse_ccs.ptr_r(), sparse_ccs.ptr_x());
    }
    else
    {
        //Petsc requires that column indices starts from 0
        //offset can be nonzero if one creates column view of
        //sparse matrix and assign sym structure

        Integer M       = sparse_ccs.cols();
        Integer offset  = sparse_ccs.offset();

        m_column_indices.assign_to_fresh(Mat_I(ti::ti_empty(), M, 1));

        const Integer* ptr_mat  = sparse_ccs.ptr_c();
        Integer* ptr_new        = m_column_indices.ptr();

        for (Integer i = 0; i < M; ++i)
            ptr_new[i]  = ptr_mat[i] - offset;

        m_Aksp      = mp::smart_mat::mat_create_seq_aij_with_arrays(PETSC_COMM_SELF, A_in.rows(), A_in.cols(),
                        ptr_new, sparse_ccs.ptr_r() + offset, sparse_ccs.ptr_x() + offset);
    };

    m_Aksp.set_struct_flags(A_in);
}

linear_operator petsc_matrix_str<struct_sparse>::get_linop() const
{
    return linear_operator(m_A);
};
Matrix petsc_matrix_str<struct_sparse>::get_matrix() const
{
    return m_A;
};

Matrix petsc_matrix_str<struct_sparse>::mmul_right(const Matrix& X, trans_type tA) const
{
    return mmul(m_A, X, tA);
}

petsc_matrix_linop::petsc_matrix_linop(const linear_operator& A)
    : petsc_matrix(A)
{
    m_Aksp = mp::smart_mat::mat_create_shell(PETSC_COMM_SELF, A);
}

Matrix petsc_matrix_linop::mmul_right(const Matrix& X, trans_type tA) const
{
    return m_Araw.mmul_right(X, tA);
}
Matrix petsc_matrix_linop::get_matrix() const
{
    throw error::unable_convert_linear_operator_to_matrix();
};

petsc_matrix::petsc_matrix(const Matrix& A_in)
    : m_is_finite( A_in.all_finite() ), m_rows(A_in.rows()), m_cols(A_in.cols())
    , m_value_code(A_in.get_value_code())
{}

petsc_matrix::petsc_matrix(const linear_operator& A_in)
    : m_is_finite(true), m_rows(A_in.rows()), m_cols(A_in.cols())
    , m_value_code(A_in.get_value_code())
{}

petsc_matrix::ptr_type petsc_matrix::create(const Matrix& A)
{
    //TODO:
    if (A.get_value_code() == value_code::v_complex)
        throw error::error_complex_value_type_not_allowed();

    if (A.get_value_code() == value_code::v_object)
        throw error::object_value_type_not_allowed("petsc_matrix");

    if(A.numel() <= 1)
        return ptr_type(new petsc_matrix_str<struct_dense>(A));
    
    switch (A.get_struct_code())
    {
        case struct_code::struct_banded:
            return ptr_type(new petsc_matrix_str<struct_sparse>(sparse(A)));
        case struct_code::struct_dense:
        case struct_code::struct_scalar:
            return ptr_type(new petsc_matrix_str<struct_dense>(A));
        case struct_code::struct_sparse:
            return ptr_type(new petsc_matrix_str<struct_sparse>(A));
        default:
            assert(0 && "Impossible Matrix struct type encountered");
            return ptr_type(new petsc_matrix_str<struct_dense>(A));
    }
}
petsc_matrix::ptr_type petsc_matrix::create(const linear_operator& A)
{
    //TODO:
    if (A.get_value_code() == value_code::v_complex)
        throw error::error_complex_value_type_not_allowed();

    if (A.get_value_code() == value_code::v_object)
        throw error::object_value_type_not_allowed("petsc_matrix");
    
    return ptr_type(new petsc_matrix_linop(A.convert(value_code::v_real)));
}

}};

#endif