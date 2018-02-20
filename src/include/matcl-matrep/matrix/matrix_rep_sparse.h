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

#include "matcl-matrep/details/fwd_decls.h"
#include "matcl-matrep/matrix/struct_flag.h"
#include "matcl-matrep/matrix/matrix.h"
#include "matcl-core/IO/disp_stream.h"
#include "matcl-core/options/matcl_options.h"
#include "matcl-matrep/matrix/matrix_rep_sub.h"
#include "matcl-matrep/matrix/matrix_concat.h"
#include "matcl-matrep/lib_functions/eval_functors.h"

namespace matcl
{

// matrix horizontal concatenation helper, i.e. constructor of [A, B, ...]
// instances of this class cannot be shared between threads
template<class T>
class sparse_row : public mat_row
{
    private:
        using base_type = mat_row;

    public:
        // construct empty matrix of Integer type
        sparse_row();        

        // standard destructor
        ~sparse_row();        

        // standard copy constructor and move constructor
        sparse_row(const sparse_row&);
        sparse_row(sparse_row&&);

        // standard assignment and move assignment;
        sparse_row&             operator=(const sparse_row&);
        sparse_row&             operator=(sparse_row&&);        

        // add scalar to the list of matrices
        sparse_row&             add(const T&);
        sparse_row&             add(T&&);
        sparse_row&             operator,(const T& v)               { return add(v); };
        sparse_row&             operator,(T&& v)                    { return add(std::move(v)); };

        // add matrix to the list of matrices
        sparse_row&             add(const sparse_matrix<T>&);
        sparse_row&             add(sparse_matrix<T>&&);
        sparse_row&             add(const band_matrix<T>&);
        sparse_row&             add(band_matrix<T>&&);
        sparse_row&             operator,(const sparse_matrix<T>& v){ return add(v); };
        sparse_row&             operator,(sparse_matrix<T>&& v)     { return add(std::move(v)); };
        sparse_row&             operator,(const band_matrix<T>& v)  { return add(v); };
        sparse_row&             operator,(band_matrix<T>&& v)       { return add(std::move(v)); };

        // add concatenation helpers to the list of matrices
        sparse_row&             add(const sparse_row&);
        sparse_row&             add(sparse_row&&);
        sparse_row&             add(const sparse_col<T>&);
        sparse_row&             add(sparse_col<T>&&);
        sparse_row&             operator,(const sparse_row& v)      { return add(v); };
        sparse_row&             operator,(sparse_row&& v)           { return add(std::move(v)); };
        sparse_row&             operator,(const sparse_col<T>& v)   { return add(v); };
        sparse_row&             operator,(sparse_col<T>&& v)        { return add(std::move(v)); };

        // convert to sparse_matrix
        const sparse_matrix<T>  to_matrix() const;
};

// matrix vertical concatenation helper, i.e. constructor of [A; B; ...]
// instances of this class cannot be shared between threads
template<class T>
class sparse_col : public mat_col
{
    private:
        using base_type = mat_col;

    public:
        // construct empty matrix of Integer type
        sparse_col();        

        // standard copy constructor and move constructor
        sparse_col(const sparse_col&);
        sparse_col(sparse_col&&);

        // standard assignment and move assignment;
        sparse_col&             operator=(const sparse_col&);
        sparse_col&             operator=(sparse_col&&);

        // standard destructor
        ~sparse_col();

        // add scalar to the list of matrices
        sparse_col&             add(const T&);
        sparse_col&             add(T&&);
        sparse_col&             operator,(const T& v)               { return add(v); };
        sparse_col&             operator,(T&& v)                    { return add(std::move(v)); };

        // add matrix to the list of matrices
        sparse_col&             add(const sparse_matrix<T>&);
        sparse_col&             add(sparse_matrix<T>&&);
        sparse_col&             add(const band_matrix<T>&);
        sparse_col&             add(band_matrix<T>&&);
        sparse_col&             operator,(const sparse_matrix<T>& v){ return add(v); };
        sparse_col&             operator,(sparse_matrix<T>&& v)     { return add(std::move(v)); };
        sparse_col&             operator,(const band_matrix<T>& v)  { return add(v); };
        sparse_col&             operator,(band_matrix<T>&& v)       { return add(std::move(v)); };

        // add concatenation helpers to the list of matrices
        sparse_col&             add(const sparse_row<T>&);
        sparse_col&             add(sparse_row<T>&&);
        sparse_col&             add(const sparse_col&);
        sparse_col&             add(sparse_col&&);
        sparse_col&             operator,(const sparse_row<T>& v)   { return add(v); };
        sparse_col&             operator,(sparse_row<T>&& v)        { return add(std::move(v)); };
        sparse_col&             operator,(const sparse_col& v)      { return add(v); };
        sparse_col&             operator,(sparse_col&& v)           { return add(std::move(v)); };

        // convert to sparse_matrix
        const sparse_matrix<T>  to_matrix() const;
};

// sparse matrix represented in column compressed format given by three arrays:
// column indices (of length equal to number of columns + 1), row indices and
// data array storing nnz (number of nonzero elements) elements denoted by ptr_c,
// ptr_r, and ptr_x; these arrays have 0-based indexing;
// for 0 <= c <= cols ptr_c[c] if the first position of data
// in arrays ptr_r and ptr_x that represent column c and ptr_c[c+1] -1 is the last
// position; all row indices in given column are sorted without duplications; therefore
// ptr_c[c+1] - ptr_c[c] is the number of nonzero elements in column c, ptr_r[ptr_c[c]]
// is the row index of the first element in given column and ptr_x[ptr_c[c]] is value
// of this element; additionally ptr_c[c] (called offset) need not be zero, this may
// happen if one take submatrix of type A(:, start:end).
// sparse_matrix bahaves according to copy-on-write scheme; sparse_matrix is 
// thread-safe, but given instance cannot be shared between threads unless all
// accesses are read only.
template<class T>
class MATCL_MATREP_EXPORT sparse_matrix<T, true>
{
    private:
        using mat_type              = raw::Matrix<T,struct_sparse>;
        using sub_sparse_matrix     = sub_sparse_matrix<T>;
        using sub_sparse_matrix_1   = sub_sparse_matrix_1<T>;
        using sub_sparse_matrix_2   = sub_sparse_matrix_2<T>;
        using dense_matrix          = dense_matrix<T,true>;

    public:
        using value_type    = T;

    private:
        Integer             m_rows;
        Integer             m_cols;
        Integer             m_max_cols;
        Integer             m_offset;
        Integer*            m_c;
        Integer**           m_r;
        value_type**        m_x;
        struct_flag*        m_flag;
        matcl::Matrix       m_matrix;

        friend sparse_matrix<T,false>;
        friend matcl::sub_sparse_matrix<T>;
        friend matcl::sub_sparse_matrix_1<T>;
        friend matcl::sub_sparse_matrix_2<T>;

    public:    
        //--------------------------------------------------------------------
        //      constructors, destructor, assignments
        //--------------------------------------------------------------------

        // default constructor create 1x1 matrix with zero scalar
        sparse_matrix();

        // conversion from bool to sparse_matrix is explicit
        explicit sparse_matrix(bool val);

        // conversion from scalar to sparse_matrix is explicit
        explicit sparse_matrix(const T& val);

        // convert general matrix to sparse_matrix; if mat is not dense, then
        // conversions take place
        explicit sparse_matrix(const matcl::Matrix& mat);
        explicit sparse_matrix(matcl::Matrix&& mat);

        // copy constructor and move constructor
        sparse_matrix(const sparse_matrix& mat);
        sparse_matrix(sparse_matrix&& mat);

        // conversion from other typed matrices
        explicit sparse_matrix(const dense_matrix& mat);
        explicit sparse_matrix(dense_matrix&& mat);
        explicit sparse_matrix(const band_matrix<T>& mat);
        explicit sparse_matrix(band_matrix<T>&& mat);

        // conversion from submatrix, a new matrix is usually created
        sparse_matrix(const sub_sparse_matrix&);
        sparse_matrix(const sub_sparse_matrix_1&);
        sparse_matrix(const sub_sparse_matrix_2&);        

        // conversion from dense_row (i.e. concatenation helper)
        sparse_matrix(const sparse_row<T>&);

        // conversion from dense_col (i.e. concatenation helper)
        sparse_matrix(const sparse_col<T>&);

        // assignments 
        sparse_matrix&      operator=(const sparse_matrix&) &;
        sparse_matrix&      operator=(sparse_matrix&&) &;

        // destructor, memory is freed if refcount drops to zero
        ~sparse_matrix();

        //--------------------------------------------------------------------
        //      functions specific to sparse_matrix
        //--------------------------------------------------------------------

        // pointer to the column indices; pointer is valid as long as this instance exists;
        // also calling non constant functions invalidates pointer
        const Integer*      ptr_c() const           { return m_c; }
        
        // pointer to the row indices; pointer is valid as long as this instance exists;
        // also calling non constant functions invalidates pointer
        const Integer*      ptr_r() const           { return *m_r; }
        
        // pointer to data; pointer is valid as long as this instance exists;
        // also calling non constant functions invalidates pointer
        const value_type*   ptr_x() const           { return *m_x; }

        // offset from pointer roots to the first stored elements; i.e. first
        // row index is ptr_r()[offset()], and first stored value is 
        // ptr_x()[offset()];
        Integer             offset() const          { return m_offset; };

        // maximum number of columns; if resize(r,c) is called with c <= max_cols()
        // and r >= rows(), then no memory allocations are needed
        Integer             max_cols() const        { return m_max_cols; }

        // number of nonzero elements
        Integer             nnz() const             { return m_c? m_c[m_cols]-m_offset : 0; }

        // maximum nnz allowd without reallocation; length of arrays ptr_x and
        // ptr_r is nzmax() + offset()
        Integer             nzmax() const;

        // check is row indices are sorted in all columns
        bool                is_sorted() const;

        // check if element in column c and row r exists; if element exists then
        // position 'k' in rows and values pointer is returned; otherwise
        // position of previous element is returner; r and c are 1-based, k is
        // 0-based
        bool                is_entry(Integer r, Integer c, Integer &k) const;

        // pointer to the column indices; pointer is valid as long as this instance exists;
        // also calling non constant functions invalidates pointer;
        // safe version of sparse_matrix allows for const access only
        const Integer*      ptr_c()                 { return m_c; }

        // pointer to the row indices; pointer is valid as long as this instance exists;
        // also calling non constant functions invalidates pointer;
        // safe version of sparse_matrix allows for const access only
        const Integer*      ptr_r()                 { return *m_r; }

        // pointer to data; pointer is valid as long as this instance exists;
        // also calling non constant functions invalidates pointer;
        // safe version of sparse_matrix allows for const access only
        const value_type*   ptr_x()                 { return *m_x; }

        // increase capacity; if s < 0, then unnecessary memory is released; if s == 0 then
        // capacity is doubled; if s > 0; then capacity is increased by s
        void                add_memory(Integer s = 0);

        // sort row indices in all columns
        void                sort();

        // convert to general matrix
        const Matrix&       to_matrix() const &     { return m_matrix; };
        Matrix&&            to_matrix() &&          { return std::move(m_matrix); };

        //--------------------------------------------------------------------
        //          const functions defined for all matrix types
        //--------------------------------------------------------------------
        
        // number of rows
        Integer             rows() const                { return m_rows; }; 

        // number of columns
        Integer             cols() const                { return m_cols; };
        
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
        const struct_flag   get_struct() const          { return *m_flag; };   
        
        // struct flag must be valid, data are not changed
        void                set_struct(const struct_flag&) const;
        
        // struct flag must be valid, data are not changed
        void                add_struct(const struct_flag&) const;

        // return type_info of stored elements
        ti::ti_object       get_type() const;

        // conversion to bool, throw error if Matrix is not convertible to scalar and
        // return true for nonzero scalars
        explicit            operator bool() const       { return (bool)m_matrix; };

        // true if this is 0xn or mx0 matrix
        bool                is_empty() const            { return rows() == 0 || cols() == 0; };

        // true is this is a scalar or 1x1 matrix
        bool                is_scalar() const           { return rows() == 1 && cols() == 1; };

        // true if rows() == cols()
        bool                is_square() const           { return rows() == cols(); };

        // true if rows() == 1 or cols() == 1 or is_empty() == true
        bool                is_vector() const           { return rows() == 1 || cols() == 1 || is_empty() == true; };

        // true if given matrix is represented as a matrix
        bool                is_matrix_type() const      { return true; }

        // true if given matrix is represented as scalar
        bool                is_scalar_type() const      { return false; };

        // true if reference count is 1
        bool                is_unique() const           { return m_matrix.is_unique(); }

        // code of stored elements type
        value_code          get_value_code() const      { return m_matrix.get_value_code(); }

        // code if matrix representation (scalar, band, dense, or sparse)
        struct_code         get_struct_code() const     { return m_matrix.get_struct_code(); }

        // matrix code (value code and struct code)
        mat_code            get_matrix_code() const     { return m_matrix.get_matrix_code(); }

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
        // can be created
        const sparse_matrix operator()(const colon& r) const;        
        const sparse_matrix operator()(const colon& r, const colon& c) const;
        
        // get diagonal
        const dense_matrix  diag(Integer d = 0) const;

        // if matrix is unique (i.e. refcount is 1) then do nothing; otherwise
        // make copy of this matrix, note that refcount of this instance will
        // also drop to one
        const sparse_matrix& make_unique() const;

        // make independent copy
        const sparse_matrix clone() const;

        //--------------------------------------------------------------------
        //          non const functions defined for all matrix types
        //--------------------------------------------------------------------

        // get or change submatrix
        sub_sparse_matrix_1 operator()(Integer r);
        sub_sparse_matrix_2 operator()(Integer r, Integer c);
        sub_sparse_matrix   operator()(const colon& r);
        sub_sparse_matrix   operator()(const colon& r, const colon& c);

        // get or change diagonal
        sub_sparse_matrix   diag(Integer d = 0);

        // change matrix size; if number of rows and columns is less than capacity,
        // then no copy is done
        void                resize(Integer r, Integer c);

        // increase rows and columns capacity
        void                reserve(Integer r, Integer c);

        // change matrix size including number of sub- and superdiagonals; 
        // this is equivalent to resize(r,c);
        void                resize_band(Integer r, Integer c, Integer ld, Integer ud);

        // increase rows columns and sub- and superdiagonals capacity; 
        // this is equivalent to reserve(r,c) 
        void                reserve_band(Integer r, Integer c, Integer ld, Integer ud);

        //--------------------------------------------------------------------
        //          static functions
        //--------------------------------------------------------------------

        // create sparse matrix with zeroes of size rxc, with allocated space for nnz elements
        template<class Enable = typename details::enable_if_not_object<T, void>::type>
        static sparse_matrix    zeros(Integer r, Integer c,Integer nnz = 0);

        // create sparse matrix with zeroes of size rxc, with allocated space for nnz elements
        // version for T == Object requires additionally type info of stored elements
        template<class Enable = typename details::enable_if_object<T, void>::type>
        static sparse_matrix    zeros(ti::ti_object ti, Integer r, Integer c,Integer nnz = 0);

        // create sparse matrix with ones of size rxc
        template<class Enable = typename details::enable_if_not_object<T, void>::type>
        static sparse_matrix    ones(Integer r, Integer c);

        // create sparse matrix with ones of size rxc
        // version for T == Object requires additionally type info of stored elements
        template<class Enable = typename details::enable_if_object<T, void>::type>
        static sparse_matrix    ones(ti::ti_object ti, Integer r, Integer c);

        // create sparse matrix with ones on diagonal of size rxc
        template<class Enable = typename details::enable_if_not_object<T, void>::type>
        static sparse_matrix    eye(Integer r, Integer c);

        // create sparse matrix with ones on diagonal of size rxr
        template<class Enable = typename details::enable_if_not_object<T, void>::type>
        static sparse_matrix    eye(Integer r);

        // create sparse matrix with ones on diagonal of size rxc
        // version for T == Object requires additionally type info of stored elements
        template<class Enable = typename details::enable_if_object<T, void>::type>
        static sparse_matrix    eye(ti::ti_object ti, Integer r, Integer c);

        // create sparse matrix with ones on diagonal of size rxr
        // version for T == Object requires additionally type info of stored elements
        template<class Enable = typename details::enable_if_object<T, void>::type>
        static sparse_matrix    eye(ti::ti_object ti, Integer r);

        // create sparse matrix with elements on diagonal d given by v of size kxk, where 
        // k = length(v) + abs(d)
        static sparse_matrix    diag(const dense_matrix& v, Integer d = 0);

        // create rxc sparse matrix from the columns of A and places them along the diagonals 
        // specified by d 
        static sparse_matrix    diags(const dense_matrix& A, const Matrix &d, Integer r, Integer c);

        // create sparse matrix of size rxc with approximatelly d*r*c nonzero elements of uniformly
        // distributed random numbers on [0,1]
        static sparse_matrix    rand(Integer r, Integer c, Real d, const rand_state& rand_ptr 
                                    = global_rand_state());

        // create sparse matrix of size rxc with approximatelly d*r*c nonzero elements of normally
        // distributed random numbers
        template<class Enable = typename details::enable_if_float<T, void>::type>
        static sparse_matrix    randn(Integer r, Integer c, Real d, const rand_state& rand_ptr 
                                    = global_rand_state());

        // create sparse matrix of size rows x cols with nnz = 0, and capacity nzmax, 
        // with zero elements
        template<class Enable = typename details::enable_if_not_object<T, void>::type>
        static sparse_matrix    make(Integer rows,Integer cols, Integer nzmax = 0);

        // create sparse matrix of size rows x cols with nnz = 0, and capacity nzmax, 
        // with zero elements; version for T == Object requires additionally type info 
        // of stored elements
        template<class Enable = typename details::enable_if_object<T, void>::type>
        static sparse_matrix    make(ti::ti_object ti, Integer rows,Integer cols, Integer nzmax = 0);

        // create sparse matrix of size r x c from triplet representation; trip_c, trip_r, trip_x
        // are arrays of length at least nnz; allocate space for at least nzmax elements
        template<class Enable = typename details::enable_if_not_object<T, void>::type>
        static sparse_matrix    make(const Integer *trip_r, const Integer *trip_c, 
                                const T* trip_x, Integer r, Integer c, Integer nnz, Integer nzmax = 0);

        // create sparse matrix of size r x c from triplet representation; trip_c, trip_r, trip_x
        // are arrays of length at least nnz; allocate space for at least nzmax elements;
        // version for T == Object requires additionally type info of stored elements
        template<class Enable = typename details::enable_if_object<T, void>::type>
        static sparse_matrix    make(ti::ti_object ti, const Integer *trip_r, const Integer *trip_c, 
                                const T* trip_x, Integer r, Integer c, Integer nnz, Integer nzmax = 0);

        // create sparse matrix of size r x c from triplet representation; trip_c, trip_r, trip_re
        // trip_im are arrays of length at least nnz; real and imaginary elements are stored in separate 
        // arrays trip_re and trip_im
        template<class Enable = typename details::enable_if_complex<T, void>::type>
        static sparse_matrix    make_complex(const Integer *trip_r, const Integer *trip_c, 
                                    const typename details::real_type<T>::type * trip_re, 
                                    const typename details::real_type<T>::type * trip_im, 
                                    Integer r, Integer c, Integer nnz, Integer nzmax = 0);

        // create sparse matrix of size r x c from triplet representation; trip_c, trip_r, trip_x
        // are vectors of equal length nnz; allocate space for at least nzmax elements
        static sparse_matrix    make(const Matrix& trip_r, const Matrix& trip_c, 
                                    const Matrix& trip_x, Integer r, Integer c, Integer nzmax = 0);

        // create sparse matrix of size rows x cols with capacity for nzmax nonzero elements; 
        // elements are not initialized;
        template<class Enable = typename details::enable_if_not_object<T, void>::type>
        static sparse_matrix    make_noinit(Integer rows,Integer cols, Integer nzmax = 0);

        // create sparse matrix of size rows x cols with capacity for nzmax nonzero elements; 
        // elements are not initialized; version for T == Object requires additionally type
        // info of stored elements
        template<class Enable = typename details::enable_if_object<T, void>::type>
        static sparse_matrix    make_noinit(ti::ti_object ti, Integer rows,Integer cols, Integer nzmax = 0);

    private:
        struct str_make_unique{};
        sparse_matrix(matcl::Matrix& mat, str_make_unique);
        sparse_matrix(matcl::Matrix&& mat, str_make_unique);

        void                init_from_rep(const mat_type& mat);
        void                update_rep();
};

// unsafe extension of sparse_matrix that allows for modification of 
// individual elements; this class is unsafe because additional requirements
// must be fulfilled, that guaratees that non const methods are called
// on unique matrices, which cannot be assured at compile time
// given instance cannot be shared between threads unless all accesses are read only
template<class T>
class MATCL_MATREP_EXPORT sparse_matrix<T,false> : public sparse_matrix<T,true>
{
    private:
        using base_type     = sparse_matrix<T,true>;

    public:
        using value_type    = T;

    public:
        //--------------------------------------------------------------------
        //      constructors, destructor, assignments
        //--------------------------------------------------------------------

        // default constructor create 1x1 matrix with zero scalar
        sparse_matrix();

        // convert general matrix to sparsr_matrix; if mat is not dense, then
        // conversions take place; this make matrix 'mat' unique, therefore copying
        // can happen; one cannot call any non const functions on 'mat' until this object
        // is destroyed
        explicit sparse_matrix(matcl::Matrix& mat);
        explicit sparse_matrix(matcl::Matrix&& mat);

        // non safe sparse_matrix cannot be copied, only moved
        sparse_matrix(const sparse_matrix<T,false>& mat) = delete;
        sparse_matrix(sparse_matrix<T,false>&& mat);

        // conversion from safe sparse_matrix; make matrix 'mat' unique, therefore
        // copying can happen; one cannot call any non const functions on 'mat' until
        // this object is destroyed
        explicit sparse_matrix(sparse_matrix<T,true>&& mat);
        explicit sparse_matrix(sparse_matrix<T,true>& mat);

        // non safe sparse_matrix cannot be assigned, only move assignment is allowed        
        sparse_matrix&      operator=(const sparse_matrix&) & = delete;
        sparse_matrix&      operator=(sparse_matrix&&) &;

        // destructor, memory is freed if refcount drops to zero
        ~sparse_matrix();

        //--------------------------------------------------------------------
        //      functions specific to sparse_matrix
        //--------------------------------------------------------------------

        // pointer to the column indices; pointer is valid as long as this instance exists;
        // also calling non constant functions invalidates pointer;
        Integer*            ptr_c()                 { return m_c; }
        
        // pointer to the row indices; pointer is valid as long as this instance exists;
        // also calling non constant functions invalidates pointer;
        Integer*            ptr_r()                 { return *m_r; }
        
        // pointer to data; pointer is valid as long as this instance exists;
        // also calling non constant functions invalidates pointer;
        value_type*         ptr_x()                 { return *m_x; }
};

};

#include "matcl-matrep/details/matrix_rep_sparse.inl"
