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
#include "matcl-matrep/matrix/matrix_rep_sub.h"
#include "matcl-matrep/lib_functions/eval_functors.h"

namespace matcl
{

// band matrix with (extended) Lapack's representation; an m-by-n band matrix with
// first diagonal indef fd and last diagonal index ld may be stored as dense matrix
// with ld-fd+1 rows and n columns; in opposite to Lapack band matrices need not have
// main diagonal; if main diagonal is present, then number of subdiagonal is -fd and
// number of superdiagonals is ld; columns of the matrix are stored in corresponding 
// columns of the matrix, and diagonals of the matrix are stored in rows of the array;
// therefore nonzero element aij is stored in B(ld+1+i-j,j) for max(1,j-ld) <= i
// <= min(m,j-fd), where B is the representation of band matrix as a dense matrix.
// band_matrix bahaves according to copy-on-write scheme; band_matrix is 
// thread-safe, but given instance cannot be shared between threads unless all
// accesses are read only.
template<class T>
class MATCL_MATREP_EXPORT band_matrix<T, true>
{
    private:
        using mat_type              = raw::Matrix<T,struct_banded>;
        using sub_band_matrix       = sub_band_matrix<T>;
        using sub_band_matrix_1     = sub_band_matrix_1<T>;
        using sub_band_matrix_2     = sub_band_matrix_2<T>;
        using dense_matrix          = dense_matrix<T,true>;
        using sparse_matrix         = sparse_matrix<T,true>;

    public:
        using value_type    = T;

    private:
        Integer             m_ldiags;
        Integer             m_udiags;
        Integer             m_rows;
        Integer             m_cols;
        Integer             m_base_ld;
        Integer             m_base_size;
        value_type*         m_base_ptr;
        struct_flag*        m_flag;
        matcl::Matrix       m_matrix;

        friend band_matrix<T,false>;
        friend matcl::sub_band_matrix<T>;
        friend matcl::sub_band_matrix_1<T>;
        friend matcl::sub_band_matrix_2<T>;

    public:    
        //--------------------------------------------------------------------
        //      constructors, destructor, assignments
        //--------------------------------------------------------------------

        // default constructor create 1x1 matrix with zero scalar
        band_matrix();

        // conversion from bool to band_matrix is explicit
        explicit band_matrix(bool val);

        // conversion from scalar to band_matrix is explicit
        explicit band_matrix(const T& val);

        // convert general matrix to band_matrix; if mat is not dense, then
        // conversions take place
        explicit band_matrix(const matcl::Matrix& mat);
        explicit band_matrix(matcl::Matrix&& mat);

        // copy constructor and move constructor
        band_matrix(const band_matrix& mat);
        band_matrix(band_matrix&& mat);

        // conversion from other typed matrices
        explicit band_matrix(const dense_matrix& mat);
        explicit band_matrix(dense_matrix&& mat);
        explicit band_matrix(const sparse_matrix& mat);
        explicit band_matrix(sparse_matrix&& mat);

        // conversion from submatrix, a new matrix is usually created
        band_matrix(const sub_band_matrix&);
        band_matrix(const sub_band_matrix_1&);
        band_matrix(const sub_band_matrix_2&);        

        // assignments 
        band_matrix&      operator=(const band_matrix&) &;
        band_matrix&      operator=(band_matrix&&) &;

        // destructor, memory is freed if refcount drops to zero
        ~band_matrix();

        //--------------------------------------------------------------------
        //      functions specific to band_matrix
        //--------------------------------------------------------------------
        // index of first nonzero diagonal (negative values for subdiagonals
        // zero for main diagonal, positive values for superdiagonals)
        inline Integer      first_diag() const;

        // index of last nonzero diagonal (negative values for subdiagonals
        // zero for main diagonal, positive values for superdiagonals)
        // last_diag >= first_diag
        inline Integer      last_diag() const;

        // equivalent to max(0, last_diag())
        inline Integer      number_superdiagonals() const;

        // equivalent to max(0, -first_diag())
        inline Integer      number_subdiagonals() const;

        // check if the matrix has main diagonal
        inline bool         has_main_diagonal() const;

        // number of elements stored in the matrix (rows x columns);
        // throw error if integer overflow occured
        Integer             size() const;

        // get the leading dimmension of the representation of a banded matrix,
        // i.e. distance between two consecutive elements in given diagional
        inline Integer      ld() const;

        // size of the dense matrix representation of this matrix
        inline Integer      base_size() const;

        // pointer to data
        inline 
        const value_type*   ptr() const;

        // index of the first row in given column; all indices are 0-based
        inline Integer      first_row(Integer col) const;

        // index of the last row in given column; all indices are 0-based
        inline Integer      last_row(Integer col) const;

        // index of the first column in given row; all indices are 0-based
        inline Integer      first_col(Integer row) const;

        // index of the last column in given row; all indices are 0-based
        inline Integer      last_col(Integer row) const;

        // index of the first column, that may contain a nonzero element; 
        // all indices are 0-based
        inline Integer      first_nonzero_column() const;

        // index of the last column, that may contain a nonzero element;
        // all indices are 0-based
        inline Integer      last_nonzero_column() const;

        // index of the first row, that may contain a nonzero element;
        // all indices are 0-based
        inline Integer      first_nonzero_row() const;

        // index of the last row, that may contain a nonzero element;
        // all indices are 0-based
        inline Integer      last_nonzero_row() const;

        // position in the array of first element in given column; in order
        // to access this element one may call ptr()[first_elem_pos(col) + col * ld()]
        inline Integer      first_elem_pos(Integer col) const;

        // position in the array of first element in given row;  in order
        // to access this element one may call ptr()[first_elem_pos_row(row)]
        inline Integer      first_elem_pos_row(Integer row) const;

        // offset of first element on the diagonal d from the root of
        // pointer to representation (d = 0 : main diagonal, d > 0: 
        // upper diagonals, d < 0: lower diagonals)
        inline Integer      first_elem_diag(Integer d) const;

        // offset of given element from the root that contains element 
        // on row r and column c if exists (r, c are 0-based)
        inline Integer      element_pos(Integer r, Integer c) const;

        // check if given diagonal d is nonzero
        inline bool         has_diag(Integer d) const;

        // length of given diagonal d
        inline Integer      diag_length(Integer d) const;

        // row index of the first element on diagonal d (0-based)
        inline Integer      first_row_on_diag(Integer d) const;

        // columns index of the first element on diagonal d (0-based)
        inline Integer      first_col_on_diag(Integer d) const;

        // convert to general matrix
        inline const Matrix& to_matrix() const &;
        inline Matrix&&      to_matrix() &&;

        // get element at given row r and column c; element must lay on 
        // one of nonzero diagonal; checks are not performed;
        // safe version of dense_matrix allows for const access only
        inline const T&     elem(Integer r, Integer c) const;
        inline const T&     elem(Integer r, Integer c);

        //--------------------------------------------------------------------
        //          const functions defined for all matrix types
        //--------------------------------------------------------------------
        
        // number of rows
        inline Integer      rows() const;

        // number of columns
        inline Integer      cols() const;
        
        // larger or rows() and cols(); return zero for empty matrices
        Integer             length() const;
        
        // estimate of number of nonzeroes elements based on representation type
        // and struct flags (without testing of stored values)
        Integer             structural_nnz() const;
        
        // estimate of number of subdiagonal based on representation type
        // and struct flags if use_flags = true, but without testing of stored values
        Integer             structural_ldiags(bool use_flags = true) const;
        
        // estimate of number of superdiagonals based on representation type
        // and struct flags if use_flags = true, but without testing of stored values
        Integer             structural_udiags(bool use_flags = true) const;
        
        // number of rows times number of columns
        Real                numel() const;        

        // check if all elements are finite
        bool                all_finite() const;

        // return additional structures assigned to this matrix
        inline 
        const struct_flag   get_struct() const;
        
        // struct flag must be valid, data are not changed
        void                set_struct(const struct_flag&) const;
        
        // struct flag must be valid, data are not changed
        void                add_struct(const struct_flag&) const;

        // return type_info of stored elements
        ti::ti_object       get_type() const;

        // conversion to bool, throw error if Matrix is not convertible to scalar and
        // return true for nonzero scalars
        inline explicit     operator bool() const;

        // true if this is 0xn or mx0 matrix
        inline bool         is_empty() const;

        // true is this is a scalar or 1x1 matrix
        inline bool         is_scalar() const;

        // true if rows() == cols()
        inline bool         is_square() const;

        // true if rows() == 1 or cols() == 1 or is_empty() == true
        inline bool         is_vector() const;

        // true if given matrix is represented as a matrix
        inline bool         is_matrix_type() const;

        // true if given matrix is represented as scalar
        inline bool         is_scalar_type() const;

        // true if reference count is 1
        inline bool         is_unique() const;

        // code of stored elements type
        inline value_code   get_value_code() const;

        // code if matrix representation (scalar, band, dense, or sparse)
        inline struct_code  get_struct_code() const;

        // matrix code (value code and struct code)
        inline mat_code     get_matrix_code() const;

        // delete rows specified by colon c.
        const sparse_matrix delrows(const colon&) const &;
        const sparse_matrix delrows(const colon&) const &&;

        // delete columns specified in colon c.
        const sparse_matrix delcols(const colon&) const &;
        const sparse_matrix delcols(const colon&) const &&;

        // delete rows specified by colon c1 and columns specified by colon c2.
        const sparse_matrix delrowscols(const colon& c1, const colon& c2) const &;
        const sparse_matrix delrowscols(const colon& c1, const colon& c2) const &&;

        // get element at position pos; checks are performed
        T                   operator()(Integer pos) const;
        // get element at given row r and column c; checks are performed
        T                   operator()(Integer r, Integer c) const;

        // get submatrix, a new matrix is usually created, unless matrix views 
        // can be created; warning: submatrix of a band matrix is generally not
        // band, if resulting matrix is not band it is better to call operator()
        // defined for Matrix (i.e. mat.to_matrix()(colon, colon) )
        const band_matrix   operator()(const colon& r) const;        
        const band_matrix   operator()(const colon& r, const colon& c) const;
        
        // get diagonal
        const dense_matrix  diag(Integer d = 0) const;

        // if matrix is unique (i.e. refcount is 1) then do nothing; otherwise
        // make copy of this matrix, note that refcount of this instance will
        // also drop to one
        band_matrix&        make_unique();

        // make independent copy
        const band_matrix   clone() const;

        //--------------------------------------------------------------------
        //          non const functions defined for all matrix types
        //--------------------------------------------------------------------

        // get or change one element
        sub_band_matrix_1   operator()(Integer r);
        sub_band_matrix_2   operator()(Integer r, Integer c);

        // get or change submatrix; warning: submatrix of a band matrix is generally not
        // band, if submatrix should be created and resulting matrix is not band it is
        // better to call operator() defined for Matrix (i.e. mat.to_matrix()(colon, colon) );
        // if assignment is performed and resulting matrix is not band, then it is better
        // to convert this matrix first to dense_matrix or sparse_matrix.
        sub_band_matrix     operator()(const colon& r);
        sub_band_matrix     operator()(const colon& r, const colon& c);

        // get or change diagonal
        sub_band_matrix     diag(Integer d = 0);

        // change matrix size; if number of rows and columns is less than capacity,
        // then no copy is done
        void                resize(Integer r, Integer c);

        // increase rows and columns capacity
        void                reserve(Integer r, Integer c);

        // change matrix size including first and last diagonal index; 
        void                resize_band(Integer r, Integer c, Integer fd, Integer ld);

        // increase rows columns and first and last diagonal capacity; 
        void                reserve_band(Integer r, Integer c, Integer fd, Integer ld);

        //--------------------------------------------------------------------
        //          static functions
        //--------------------------------------------------------------------

        // create band matrix with zeroes of size rxc, with first diagonal fd and last diagonal ld
        template<class Enable = typename details::enable_if_not_object<T, void>::type>
        static band_matrix  zeros(Integer r, Integer c,Integer fd, Integer ld);

        // create band matrix with zeroes of size rxc, with first diagonal fd and last diagonal ld
        // version for T == Object requires additionally type info of stored elements
        template<class Enable = typename details::enable_if_object<T, void>::type>
        static band_matrix  zeros(ti::ti_object ti, Integer r, Integer c,Integer fd, Integer ld);

        // create band matrix with ones of size rxc
        template<class Enable = typename details::enable_if_not_object<T, void>::type>
        static band_matrix  ones(Integer r, Integer c);

        // create band matrix with ones of size rxc
        // version for T == Object requires additionally type info of stored elements
        template<class Enable = typename details::enable_if_object<T, void>::type>
        static band_matrix  ones(ti::ti_object ti, Integer r, Integer c);

        // create band matrix with ones on diagonal of size rxc with first diagonal fd
        // and last diagonal ld
        template<class Enable = typename details::enable_if_not_object<T, void>::type>
        static band_matrix  eye(Integer r, Integer c, Integer fd, Integer ld);

        // create band matrix with ones on diagonal of size rxr with first diagonal fd
        // and last diagonal ld
        template<class Enable = typename details::enable_if_not_object<T, void>::type>
        static band_matrix  eye(Integer r, Integer fd, Integer ld);

        // create band matrix with ones on diagonal of size rxc with first diagonal fd
        // and last diagonal ld; version for T == Object requires additionally type
        // info of stored elements
        template<class Enable = typename details::enable_if_object<T, void>::type>
        static band_matrix  eye(ti::ti_object ti, Integer r, Integer c, Integer fd, Integer ld);

        // create band matrix with ones on diagonal of size rxr with first diagonal fd
        // and last diagonal ld; version for T == Object requires additionally type
        // info of stored elements
        template<class Enable = typename details::enable_if_object<T, void>::type>
        static band_matrix  eye(ti::ti_object ti, Integer r, Integer fd, Integer ld);

        // create band matrix with elements on diagonal d given by v of size kxk, where 
        // k = length(v) + abs(d)
        static band_matrix  diag(const dense_matrix& v, Integer d = 0);

        // create rxc band matrix from the columns of A and places them along the diagonals 
        // specified by d 
        static band_matrix  diags(const dense_matrix& A, const Matrix &d, Integer r, Integer c);

        // create band matrix of size rxc with first diagonal fd and last diagonal ld of uniformly
        // distributed random numbers on [0,1]
        static band_matrix  rand(Integer r, Integer c, Integer fd, Integer ld, 
                                const rand_state& rand_ptr = global_rand_state());

        // create band matrix of size rxc with first diagonal fd and last diagonal ld of normally
        // distributed random numbers
        template<class Enable = typename details::enable_if_float<T, void>::type>
        static band_matrix  randn(Integer r, Integer c, Integer fd, Integer ld, 
                                const rand_state& rand_ptr = global_rand_state());

        // create band matrix of size rows x cols with first diagonal fd and last diagonal ld,
        // with zeroes
        template<class Enable = typename details::enable_if_not_object<T, void>::type>
        static band_matrix  make(Integer rows,Integer cols, Integer fd, Integer ld);

        // create band matrix of size rows x cols with first diagonal fd and last diagonal ld,
        // with zeroes; version for T == Object requires additionally type info of stored elements
        template<class Enable = typename details::enable_if_object<T, void>::type>
        static band_matrix  make(ti::ti_object ti, Integer rows,Integer cols, Integer fd, Integer ld);

        // create band matrix of size rows x cols with first diagonal fd and last diagonal ld,
        // with elements val
        static band_matrix  make(const T& val,Integer rows,Integer cols, Integer fd, Integer ld);

        // create band matrix of size rows x cols with first diagonal fd and last diagonal ld; 
        // elements are not initialized
        template<class Enable = typename details::enable_if_not_object<T, void>::type>
        static band_matrix  make_noinit(Integer rows,Integer cols, Integer fd, Integer ld);

        // create band matrix of size rows x cols with first diagonal fd and last diagonal ld; 
        // elements are not initialized; version for T == Object requires additionally type
        // info of stored elements
        template<class Enable = typename details::enable_if_object<T, void>::type>
        static band_matrix  make_noinit(ti::ti_object ti, Integer rows,Integer cols, Integer fd, 
                                Integer ld);

    private:
        struct str_make_unique{};
        band_matrix(matcl::Matrix& mat, str_make_unique);
        band_matrix(matcl::Matrix&& mat, str_make_unique);

        void                init_from_rep(const mat_type& mat);
        void                update_rep();
};

// unsafe extension of band_matrix that allows for modification of 
// individual elements; this class is unsafe because additional requirements
// must be fulfilled, that guaratees that non const methods are called
// on unique matrices, which cannot be assured at compile time
// given instance cannot be shared between threads unless all accesses are read only
template<class T>
class MATCL_MATREP_EXPORT band_matrix<T,false> : public band_matrix<T,true>
{
    private:
        using base_type     = band_matrix<T,true>;

    public:
        using value_type    = T;

    public:
        //--------------------------------------------------------------------
        //      constructors, destructor, assignments
        //--------------------------------------------------------------------

        // default constructor create 1x1 matrix with zero scalar
        band_matrix();

        // convert general matrix to sparsr_matrix; if mat is not dense, then
        // conversions take place; this make matrix 'mat' unique, therefore copying
        // can happen; one cannot call any non const functions on 'mat' until this object
        // is destroyed
        explicit band_matrix(matcl::Matrix& mat);
        explicit band_matrix(matcl::Matrix&& mat);

        // non safe band_matrix cannot be copied, only moved
        band_matrix(const band_matrix<T,false>& mat) = delete;
        band_matrix(band_matrix<T,false>&& mat);

        // conversion from safe band_matrix; make matrix 'mat' unique, therefore
        // copying can happen; one cannot call any non const functions on 'mat' until
        // this object is destroyed
        explicit band_matrix(band_matrix<T,true>&& mat);
        explicit band_matrix(band_matrix<T,true>& mat);

        // non safe band_matrix cannot be assigned, only move assignment is allowed        
        band_matrix&        operator=(const band_matrix&) & = delete;
        band_matrix&        operator=(band_matrix&&) &;

        // destructor, memory is freed if refcount drops to zero
        ~band_matrix();

        //--------------------------------------------------------------------
        //      functions specific to band_matrix
        //--------------------------------------------------------------------
        
        // pointer to data
        inline value_type*  ptr();

        // get or modify element at given row r and column c; element must lay on 
        // one of nonzero diagonal; checks are not performed;
        inline T&           elem(Integer r, Integer c);
};

};

#include "matcl-matrep/details/matrix_rep_band.inl"
