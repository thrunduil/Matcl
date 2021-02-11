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

#include "matcl-matrep/details/fwd_decls.h"
#include "matcl-matrep/matrix/struct_flag.h"
#include "matcl-matrep/matrix/matrix.h"

namespace matcl
{

// class representing a submatrix A(c) or A(c1, c2), where c, c1, c2
// are colons; only Matrix can create instance of this class; this class
// is intended to live as short as possible
template<class T>
class MATCL_MATREP_EXPORT sub_dense_matrix
{
    private:
        using matrix_type   = dense_matrix<T, true>;

    private:
        matrix_type*        m_matrix;
        const colon*        m_colon_1;
        const colon*        m_colon_2;
        Integer             m_d;
        
    private:
        inline sub_dense_matrix(matrix_type* m, const colon& c1);
        inline sub_dense_matrix(matrix_type* m, const colon& c1, const colon& c2);
        inline sub_dense_matrix(Integer diag, matrix_type* m);

    public:        
        // move constructor; should not be used explicitly
        sub_dense_matrix(sub_dense_matrix&& m);

        // create submatrix
        const matrix_type   to_matrix() const;

        // create submatrix as a raw Matrix
        const Matrix        to_raw_matrix() const;

        // assign matrix or submatrix
        matrix_type&        operator=(const matrix_type&) const &&;
        matrix_type&        operator=(const T&) const &&;
        matrix_type&        operator=(const sub_dense_matrix<T>&) const &&;

        //--------------------------------------------------------------------
        //      selected functions defined for Matrix
        //--------------------------------------------------------------------
        // get scalar of given type; matrix must store scalars or 1x1 matrices;
        // conversions are allowed
        template<class V> V get_scalar() const;
 
        // get number of rows; return 1 for scalars
        Integer             rows() const;

        // get number of columns; return 1 for scalars
        Integer             cols() const;
        
        // larger or rows() and cols(); return zero for empty matrices
        Integer             length() const;
                
        // number of rows times number of columns
        Real                numel() const;

        // check if all elements are finite
        bool                all_finite() const;

        // conversion to bool, throw error if Matrix is not convertible to scalar and
        // return true for nonzero scalars
        explicit            operator bool() const;

    private:
        sub_dense_matrix(){};
        inline sub_dense_matrix(const sub_dense_matrix& m);

        friend matrix_type;
};

template<class T>
class sub_band_matrix;

template<class T>
class sub_band_matrix_1;

template<class T>
class sub_band_matrix_2;

// class representing a submatrix A(c) or A(c1, c2), where c, c1, c2
// are colons; only Matrix can create instance of this class; this class
// is intended to live as short as possible
template<class T>
class MATCL_MATREP_EXPORT sub_band_matrix
{
    private:
        using matrix_type   = band_matrix<T, true>;

    private:
        matrix_type*        m_matrix;
        const colon*        m_colon_1;
        const colon*        m_colon_2;
        Integer             m_d;
        
    private:
        inline sub_band_matrix(matrix_type* m, const colon& c1);
        inline sub_band_matrix(matrix_type* m, const colon& c1, const colon& c2);
        inline sub_band_matrix(Integer diag, matrix_type* m);

    public:        
        // move constructor; should not be used explicitly
        sub_band_matrix(sub_band_matrix&& m);

        // create submatrix
        const matrix_type   to_matrix() const;

        // create submatrix as a raw Matrix
        const Matrix        to_raw_matrix() const;

        // assign matrix or submatrix
        matrix_type&        operator=(const matrix_type&) const &&;
        matrix_type&        operator=(const T&) const &&;
        matrix_type&        operator=(const sub_band_matrix<T>&) const &&;
        matrix_type&        operator=(const sub_band_matrix_1<T>&) const &&;
        matrix_type&        operator=(const sub_band_matrix_2<T>&) const &&;

        //--------------------------------------------------------------------
        //      selected functions defined for Matrix
        //--------------------------------------------------------------------
        // get scalar of given type; matrix must store scalars or 1x1 matrices;
        // conversions are allowed
        template<class V> V get_scalar() const;
 
        // get number of rows; return 1 for scalars
        Integer             rows() const;

        // get number of columns; return 1 for scalars
        Integer             cols() const;
        
        // larger or rows() and cols(); return zero for empty matrices
        Integer             length() const;
                
        // number of rows times number of columns
        Real                numel() const;

        // check if all elements are finite
        bool                all_finite() const;

        // conversion to bool, throw error if Matrix is not convertible to scalar and
        // return true for nonzero scalars
        explicit            operator bool() const;
    private:
        sub_band_matrix(){};
        inline sub_band_matrix(const sub_band_matrix& m);

        friend matrix_type;
};

// class representing a submatrix A(ind), where ind is Integer
// only Matrix can create instance of this class; this class
// is intended to live as short as possible
template<class T>
class MATCL_MATREP_EXPORT sub_band_matrix_1
{
    private:
        using matrix_type   = band_matrix<T, true>;

    private:
        matrix_type*        m_matrix;
        Integer             m_ind_1;
        
    private:
        inline sub_band_matrix_1(matrix_type* m, Integer r);

    public:        
        // move constructor; should not be used explicitly
        sub_band_matrix_1(sub_band_matrix_1&& m);

        // create submatrix
        T                   to_matrix() const;

        // assign matrix or submatrix
        matrix_type&        operator=(const matrix_type&) const &&;
        matrix_type&        operator=(const T&) const &&;
        matrix_type&        operator=(const sub_band_matrix<T>&) const &&;
        matrix_type&        operator=(const sub_band_matrix_1<T>&) const &&;
        matrix_type&        operator=(const sub_band_matrix_2<T>&) const &&;

        //--------------------------------------------------------------------
        //      selected functions defined for Matrix
        //--------------------------------------------------------------------
        // get scalar of given type; matrix must store scalars or 1x1 matrices;
        // conversions are allowed
        template<class V> V get_scalar() const;
 
        // get number of rows; return 1 for scalars
        Integer             rows() const;

        // get number of columns; return 1 for scalars
        Integer             cols() const;
        
        // larger or rows() and cols(); return zero for empty matrices
        Integer             length() const;
                
        // number of rows times number of columns
        Real                numel() const;

        // check if all elements are finite
        bool                all_finite() const;

        // conversion to bool, throw error if Matrix is not convertible to scalar and
        // return true for nonzero scalars
        explicit            operator bool() const;

    private:
        sub_band_matrix_1(){};
        inline sub_band_matrix_1(const sub_band_matrix_1& m);

        friend matrix_type;
};

// class representing a submatrix A(ind1, ind2), where ind1, ind2 are Integer
// only Matrix can create instance of this class; this class
// is intended to live as short as possible
template<class T>
class MATCL_MATREP_EXPORT sub_band_matrix_2
{
    private:
        using matrix_type   = band_matrix<T, true>;

    private:
        matrix_type*        m_matrix;
        Integer             m_ind_1;
        Integer             m_ind_2;
        
    private:
        inline sub_band_matrix_2(matrix_type* m, Integer r, Integer c);

    public:        
        // move constructor; should not be used explicitly
        sub_band_matrix_2(sub_band_matrix_2&& m);

        // create submatrix
        T                   to_matrix() const;

        // assign matrix or submatrix
        matrix_type&        operator=(const matrix_type&) const &&;
        matrix_type&        operator=(const T&) const &&;
        matrix_type&        operator=(const sub_band_matrix<T>&) const &&;
        matrix_type&        operator=(const sub_band_matrix_1<T>&) const &&;
        matrix_type&        operator=(const sub_band_matrix_2<T>&) const &&;

        //--------------------------------------------------------------------
        //      selected functions defined for Matrix
        //--------------------------------------------------------------------
        // get scalar of given type; matrix must store scalars or 1x1 matrices;
        // conversions are allowed
        template<class V> V get_scalar() const;
 
        // get number of rows; return 1 for scalars
        Integer             rows() const;

        // get number of columns; return 1 for scalars
        Integer             cols() const;
        
        // larger or rows() and cols(); return zero for empty matrices
        Integer             length() const;
                
        // number of rows times number of columns
        Real                numel() const;

        // check if all elements are finite
        bool                all_finite() const;

        // conversion to bool, throw error if Matrix is not convertible to scalar and
        // return true for nonzero scalars
        explicit            operator bool() const;

    private:
        sub_band_matrix_2(){};
        inline sub_band_matrix_2(const sub_band_matrix_2& m);

        friend matrix_type;
};

template<class T>
class sub_sparse_matrix;

template<class T>
class sub_sparse_matrix_1;

template<class T>
class sub_sparse_matrix_2;

// class representing a submatrix A(c) or A(c1, c2), where c, c1, c2
// are colons; only Matrix can create instance of this class; this class
// is intended to live as short as possible
template<class T>
class MATCL_MATREP_EXPORT sub_sparse_matrix
{
    private:
        using matrix_type   = sparse_matrix<T, true>;

    private:
        matrix_type*        m_matrix;
        const colon*        m_colon_1;
        const colon*        m_colon_2;
        Integer             m_d;
        
    private:
        inline sub_sparse_matrix(matrix_type* m, const colon& c1);
        inline sub_sparse_matrix(matrix_type* m, const colon& c1, const colon& c2);
        inline sub_sparse_matrix(Integer diag, matrix_type* m);

    public:        
        // move constructor; should not be used explicitly
        sub_sparse_matrix(sub_sparse_matrix&& m);

        // create submatrix
        const matrix_type   to_matrix() const;

        // create submatrix as a raw Matrix
        const Matrix        to_raw_matrix() const;

        // assign matrix or submatrix
        matrix_type&        operator=(const matrix_type&) const &&;
        matrix_type&        operator=(const T&) const &&;
        matrix_type&        operator=(const sub_sparse_matrix<T>&) const &&;
        matrix_type&        operator=(const sub_sparse_matrix_1<T>&) const &&;
        matrix_type&        operator=(const sub_sparse_matrix_2<T>&) const &&;

        // remove elements v satisfying abs(v) <= tol;
        matrix_type&        drop_sparse(Real tol = 0.0) const &&;

        // add selected elements to the structure (i.e. do nothing if element
        // exists or add new element with zero value if element does not exist); 
        matrix_type&        add_sparse() const &&;

        //--------------------------------------------------------------------
        //      selected functions defined for Matrix
        //--------------------------------------------------------------------
        // get scalar of given type; matrix must store scalars or 1x1 matrices;
        // conversions are allowed
        template<class V> V get_scalar() const;
 
        // get number of rows; return 1 for scalars
        Integer             rows() const;

        // get number of columns; return 1 for scalars
        Integer             cols() const;
        
        // larger or rows() and cols(); return zero for empty matrices
        Integer             length() const;
                
        // number of rows times number of columns
        Real                numel() const;

        // check if all elements are finite
        bool                all_finite() const;

        // conversion to bool, throw error if Matrix is not convertible to scalar and
        // return true for nonzero scalars
        explicit            operator bool() const;

    private:
        sub_sparse_matrix(){};
        inline sub_sparse_matrix(const sub_sparse_matrix& m);

        friend matrix_type;
};

// class representing a submatrix A(ind), where ind is Integer
// only Matrix can create instance of this class; this class
// is intended to live as short as possible
template<class T>
class MATCL_MATREP_EXPORT sub_sparse_matrix_1
{
    private:
        using matrix_type   = sparse_matrix<T, true>;

    private:
        matrix_type*        m_matrix;
        Integer             m_ind_1;
        
    private:
        inline sub_sparse_matrix_1(matrix_type* m, Integer r);

    public:        
        // move constructor; should not be used explicitly
        sub_sparse_matrix_1(sub_sparse_matrix_1&& m);

        // create submatrix
        T                   to_matrix() const;

        // assign matrix or submatrix
        matrix_type&        operator=(const matrix_type&) const &&;
        matrix_type&        operator=(const T&) const &&;
        matrix_type&        operator=(const sub_sparse_matrix<T>&) const &&;
        matrix_type&        operator=(const sub_sparse_matrix_1<T>&) const &&;
        matrix_type&        operator=(const sub_sparse_matrix_2<T>&) const &&;

        // remove elements v satisfying abs(v) <= tol;
        matrix_type&        drop_sparse(Real tol = 0.0) const &&;

        // add selected elements to the structure (i.e. do nothing if element
        // exists or add new element with zero value if element does not exist); 
        matrix_type&        add_sparse() const &&;

        //--------------------------------------------------------------------
        //      selected functions defined for Matrix
        //--------------------------------------------------------------------
        // get scalar of given type; matrix must store scalars or 1x1 matrices;
        // conversions are allowed
        template<class V> V get_scalar() const;
 
        // get number of rows; return 1 for scalars
        Integer             rows() const;

        // get number of columns; return 1 for scalars
        Integer             cols() const;
        
        // larger or rows() and cols(); return zero for empty matrices
        Integer             length() const;
                
        // number of rows times number of columns
        Real                numel() const;

        // check if all elements are finite
        bool                all_finite() const;

        // conversion to bool, throw error if Matrix is not convertible to scalar and
        // return true for nonzero scalars
        explicit            operator bool() const;

    private:
        sub_sparse_matrix_1(){};
        inline sub_sparse_matrix_1(const sub_sparse_matrix_1& m);

        friend matrix_type;
};

// class representing a submatrix A(ind1, ind2), where ind1, ind2 are Integer
// only Matrix can create instance of this class; this class
// is intended to live as short as possible
template<class T>
class MATCL_MATREP_EXPORT sub_sparse_matrix_2
{
    private:
        using matrix_type   = sparse_matrix<T, true>;

    private:
        matrix_type*        m_matrix;
        Integer             m_ind_1;
        Integer             m_ind_2;
        
    private:
        inline sub_sparse_matrix_2(matrix_type* m, Integer r, Integer c);

    public:        
        // move constructor; should not be used explicitly
        sub_sparse_matrix_2(sub_sparse_matrix_2&& m);

        // create submatrix
        T                   to_matrix() const;

        // assign matrix or submatrix
        matrix_type&        operator=(const matrix_type&) const &&;
        matrix_type&        operator=(const T&) const &&;
        matrix_type&        operator=(const sub_sparse_matrix<T>&) const &&;
        matrix_type&        operator=(const sub_sparse_matrix_1<T>&) const &&;
        matrix_type&        operator=(const sub_sparse_matrix_2<T>&) const &&;

        // remove elements v satisfying abs(v) <= tol;
        matrix_type&        drop_sparse(Real tol = 0.0) const &&;

        // add selected elements to the structure (i.e. do nothing if element
        // exists or add new element with zero value if element does not exist); 
        matrix_type&        add_sparse() const &&;

        //--------------------------------------------------------------------
        //      selected functions defined for Matrix
        //--------------------------------------------------------------------
        // get scalar of given type; matrix must store scalars or 1x1 matrices;
        // conversions are allowed
        template<class V> V get_scalar() const;
 
        // get number of rows; return 1 for scalars
        Integer             rows() const;

        // get number of columns; return 1 for scalars
        Integer             cols() const;
        
        // larger or rows() and cols(); return zero for empty matrices
        Integer             length() const;
                
        // number of rows times number of columns
        Real                numel() const;

        // check if all elements are finite
        bool                all_finite() const;

        // conversion to bool, throw error if Matrix is not convertible to scalar and
        // return true for nonzero scalars
        explicit            operator bool() const;

    private:
        sub_sparse_matrix_2(){};
        inline sub_sparse_matrix_2(const sub_sparse_matrix_2& m);

        friend matrix_type;
};

};
