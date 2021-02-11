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
#include "matcl-core/error/exception_classes.h"
#include "matcl-matrep/details/isa.h"
#include "matcl-matrep/details/extract_type_switch.h" 
#include "matcl-matrep/details/assign_visitor.h" 

namespace matcl { namespace details
{
    //  Handle superfluous cloning on self-assignment by appriopriate modification
    //  of RHS matrix 
    //  
    //  In case when LHS in assignment is unique we force make_unique on RHS
    //  when dealing with self-assignment, i.e. V[...] = V[...]
    void MATCL_MATREP_EXPORT prepare_for_assign(const Matrix& lhs, const Matrix& rhs);
    void MATCL_MATREP_EXPORT prepare_for_assign(const Matrix& lhs, const Matrix& rhs1, 
                                                const Matrix& rhs2);

};};

namespace matcl
{

// class representing a submatrix A(c) or A(c1, c2), where c, c1, c2
// are colons; only Matrix can create instance of this class; this class
// is intended to live as short as possible
class MATCL_MATREP_EXPORT sub_matrix
{
    private:
        Matrix*                 m_matrix;
        const colon*            m_colon_1;
        const colon*            m_colon_2;
        Integer                 m_d;
        
    private:
        sub_matrix(Matrix* m, const colon& c1);
        sub_matrix(Matrix* m, const colon& c1, const colon& c2);
        sub_matrix(Integer diag, Matrix* m);

        template<class S>
        Matrix& assign_scalar(const S& val) const;

    public:
        // move constructor; should not be used explicitly
        sub_matrix(sub_matrix&& m);

        // create submatrix
        const Matrix            to_matrix() const;

        // assign scalar
        template<class S>
        typename details::enable_if<details::is_scalar<S>::value,Matrix&>::type
        operator=(const S&) const &&;

        // assign matrix or submatrix
        Matrix&                 operator=(const Matrix&) const &&;
        Matrix&                 operator=(const sub_matrix&) const &&;
        Matrix&                 operator=(const sub_matrix_1&) const &&;
        Matrix&                 operator=(const sub_matrix_2&) const &&;

        // for sparse matrices remove elements v satisfying abs(v) <= tol;
        // do nothing for dense and band matrices
        Matrix&                 drop_sparse(Real tol = 0.0) const &&;

        // for sparse matrices add selected elements to the structure 
        // (i.e. do nothing if element exists or add new element with zero
        // value if element does not exist); 
        // do nothing for dense and band matrices
        Matrix&                 add_sparse() const &&;

        // represent given submatrix A(c1, c2) as:
        //     A(first_row:1:first_row+rows-1, first_col:1:first_col+cols-1)
        // throw exception is this is not possible
        Matrix&                 get_view_info(Integer& first_row, Integer& rows, 
                                    Integer& first_col, Integer& cols) const;

        //--------------------------------------------------------------------
        //      selected functions defined for Matrix
        //--------------------------------------------------------------------
        // get scalar of given type; matrix must store scalars or 1x1 matrices;
        // conversions are allowed
        template<class V> V     get_scalar() const;
 
        // get number of rows; return 1 for scalars
        Integer                 rows() const;

        // get number of columns; return 1 for scalars
        Integer                 cols() const;
        
        // larger or rows() and cols(); return zero for empty matrices
        Integer                 length() const;
                
        // number of rows times number of columns
        Real                    numel() const;

        // check if all elements are finite
        bool                    all_finite() const;

        // conversion to bool, throw error if Matrix is not convertible to scalar and
        // return true for nonzero scalars
        explicit                operator bool() const;

    private:
        sub_matrix(){};
        sub_matrix(const sub_matrix& m);

        // assign matrix or submatrix, assume that this submatrix is effectively
        // unique
        Matrix&                 assign_unique(const Matrix&) const &&;
        Matrix&                 assign_unique(const sub_matrix&) const &&;
        Matrix&                 assign_unique(const sub_matrix_1&) const &&;
        Matrix&                 assign_unique(const sub_matrix_2&) const &&;

        friend Matrix;
        friend sub_matrix_1;
        friend sub_matrix_2;
        friend unique_matrix;
};

// class representing a submatrix A(ind), where ind is Integer
// only Matrix can create instance of this class; this class
// is intended to live as short as possible
class MATCL_MATREP_EXPORT sub_matrix_1
{
    private:
        Matrix*                 m_matrix;
        Integer                 m_ind_1;

    private:
        sub_matrix_1(Matrix* m, Integer i);

        template<class S>
        Matrix& assign_scalar(const S& val) const;

    public:        
        // move constructor; should not be used explicitly
        sub_matrix_1(sub_matrix_1&& m);

        // create submatrix
        const Matrix            to_matrix() const;

        // assign scalar
        template<class S>
        typename details::enable_if<details::is_scalar<S>::value,Matrix&>::type
        operator=(const S&) const &&;

        // assign matrix or submatrix
        Matrix&                 operator=(const Matrix&) const &&;
        Matrix&                 operator=(const sub_matrix&) const &&;
        Matrix&                 operator=(const sub_matrix_1&) const &&;
        Matrix&                 operator=(const sub_matrix_2&) const &&;

        // for sparse matrices remove elements v satisfying abs(v) <= tol;
        // do nothing for dense and band matrices
        Matrix&                 drop_sparse(Real tol = 0.0) const &&;

        // for sparse matrices add selected elements to the structure 
        // (i.e. do nothing if element exists or add new element with zero
        // value if element does not exist); 
        // do nothing for dense and band matrices
        Matrix&                 add_sparse() const &&;

        // represent given submatrix A(i) as:
        //     A(first_row:1:first_row+rows-1, first_col:1:first_col+cols-1)
        // throw exception is this is not possible
        Matrix&                 get_view_info(Integer& first_row, Integer& rows, 
                                    Integer& first_col, Integer& cols) const;

        //--------------------------------------------------------------------
        //      selected functions defined for Matrix
        //--------------------------------------------------------------------
        // get scalar of given type; matrix must store scalars or 1x1 matrices;
        // conversions are allowed
        template<class V> V     get_scalar() const;
 
        // get number of rows; return 1 for scalars
        Integer                 rows() const;

        // get number of columns; return 1 for scalars
        Integer                 cols() const;
        
        // larger or rows() and cols(); return zero for empty matrices
        Integer                 length() const;
                
        // number of rows times number of columns
        Real                    numel() const;

        // check if all elements are finite
        bool                    all_finite() const;

        // conversion to bool, throw error if Matrix is not convertible to scalar and
        // return true for nonzero scalars
        explicit                operator bool() const;

    private:
        sub_matrix_1(){};
        sub_matrix_1(const sub_matrix_1& m);

        // assign matrix or submatrix, assume that this submatrix is effectively
        // unique
        Matrix&                 assign_unique(const Matrix&) const &&;
        Matrix&                 assign_unique(const sub_matrix&) const &&;
        Matrix&                 assign_unique(const sub_matrix_1&) const &&;
        Matrix&                 assign_unique(const sub_matrix_2&) const &&;

        friend Matrix;
        friend sub_matrix;
        friend sub_matrix_2;
        friend unique_matrix;
};

// class representing a submatrix A(ind1, ind2), where ind1, ind2 are Integer
// only Matrix can create instance of this class; this class
// is intended to live as short as possible
class MATCL_MATREP_EXPORT sub_matrix_2
{
    private:
        Matrix*                 m_matrix;
        Integer                 m_ind_1;
        Integer                 m_ind_2;

    private:
        sub_matrix_2(Matrix* m, Integer i, Integer j);

        template<class S>
        Matrix& assign_scalar(const S& val) const;

    public:        
        // move constructor; should not be used explicitly
        sub_matrix_2(sub_matrix_2&& m);

        // create submatrix
        const Matrix            to_matrix() const;

        // assign scalar
        template<class S>
        typename details::enable_if<details::is_scalar<S>::value,Matrix&>::type
        operator=(const S&) const &&;

        // assign matrix or submatrix
        Matrix&                 operator=(const Matrix&) const &&;
        Matrix&                 operator=(const sub_matrix&) const &&;
        Matrix&                 operator=(const sub_matrix_1&) const &&;
        Matrix&                 operator=(const sub_matrix_2&) const &&;

        // for sparse matrices remove elements v satisfying abs(v) <= tol;
        // do nothing for dense and band matrices
        Matrix&                 drop_sparse(Real tol = 0.0) const &&;

        // for sparse matrices add selected elements to the structure 
        // (i.e. do nothing if element exists or add new element with zero
        // value if element does not exist); 
        // do nothing for dense and band matrices
        Matrix&                 add_sparse() const &&;

        // represent given submatrix A(i, j) as:
        //     A(first_row:1:first_row+rows-1, first_col:1:first_col+cols-1)
        // throw exception is this is not possible
        Matrix&                 get_view_info(Integer& first_row, Integer& rows, 
                                    Integer& first_col, Integer& cols) const;

        //--------------------------------------------------------------------
        //      selected functions defined for Matrix
        //--------------------------------------------------------------------
        // get scalar of given type; matrix must store scalars or 1x1 matrices;
        // conversions are allowed
        template<class V> V     get_scalar() const;
 
        // get number of rows; return 1 for scalars
        Integer                 rows() const;

        // get number of columns; return 1 for scalars
        Integer                 cols() const;
        
        // larger or rows() and cols(); return zero for empty matrices
        Integer                 length() const;
                
        // number of rows times number of columns
        Real                    numel() const;

        // check if all elements are finite
        bool                    all_finite() const;

        // conversion to bool, throw error if Matrix is not convertible to scalar and
        // return true for nonzero scalars
        explicit                operator bool() const;

    private:
        sub_matrix_2();
        sub_matrix_2(const sub_matrix_2& m);

        // assign matrix or submatrix, assume that this submatrix is effectively
        // unique
        Matrix&                 assign_unique(const Matrix&) const &&;
        Matrix&                 assign_unique(const sub_matrix&) const &&;
        Matrix&                 assign_unique(const sub_matrix_1&) const &&;
        Matrix&                 assign_unique(const sub_matrix_2&) const &&;

        friend Matrix;
        friend sub_matrix_1;
        friend sub_matrix;
        friend unique_matrix;
};

};

#include "matcl-matrep/details/matrix_details_subs.h"
