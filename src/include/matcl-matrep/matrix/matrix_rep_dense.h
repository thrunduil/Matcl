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
#include "matcl-matrep/matrix/colon.h"
#include "matcl-matrep/matrix/matrix_rep_sub.h"
#include "matcl-matrep/matrix/matrix_concat.h"
#include "matcl-matrep/lib_functions/eval_functors.h"

namespace matcl
{

namespace md = matcl::details;

// matrix horizontal concatenation helper, i.e. constructor of [A, B, ...]
// instances of this class cannot be shared between threads
template<class T>
class dense_row : public mat_row
{
    private:
        using base_type = mat_row;

    public:
        // construct empty matrix of Integer type
        dense_row();        

        // standard destructor
        ~dense_row();        

        // standard copy constructor and move constructor
        dense_row(const dense_row&);
        dense_row(dense_row&&);

        // standard assignment and move assignment;
        dense_row&              operator=(const dense_row&);
        dense_row&              operator=(dense_row&&);        

        // add scalar to the list of matrices
        dense_row&              add(const T&);
        dense_row&              add(T&&);
        dense_row&              operator,(const T& v)               {return add(v);};
        dense_row&              operator,(T&& v)                    {return add(std::move(v));};

        // add matrix to the list of matrices
        dense_row&              add(const dense_matrix<T>&);
        dense_row&              add(dense_matrix<T>&&);
        dense_row&              operator,(const dense_matrix<T>& v) {return add(v);};
        dense_row&              operator,(dense_matrix<T>&& v)      {return add(std::move(v));};

        // add concatenation helpers to the list of matrices
        dense_row&              add(const dense_row&);
        dense_row&              add(dense_row&&);
        dense_row&              add(const dense_col<T>&);
        dense_row&              add(dense_col<T>&&);
        dense_row&              operator,(const dense_row& v)       {return add(v);};
        dense_row&              operator,(dense_row&& v)            {return add(std::move(v));};
        dense_row&              operator,(const dense_col<T>& v)    {return add(v);};
        dense_row&              operator,(dense_col<T>&& v)         {return add(std::move(v));};

        // convert to dense_matrix
        const dense_matrix<T>   to_matrix() const;
};

// matrix vertical concatenation helper, i.e. constructor of [A; B; ...]
// instances of this class cannot be shared between threads
template<class T>
class dense_col : public mat_col
{
    private:
        using base_type = mat_col;

    public:
        // construct empty matrix of Integer type
        dense_col();        

        // standard copy constructor and move constructor
        dense_col(const dense_col&);
        dense_col(dense_col&&);

        // standard assignment and move assignment;
        dense_col&              operator=(const dense_col&);
        dense_col&              operator=(dense_col&&);

        // standard destructor
        ~dense_col();

        // add scalar to the list of matrices
        dense_col&              add(const T&);
        dense_col&              add(T&&);
        dense_col&              operator,(const T& v)               { return add(v); };
        dense_col&              operator,(T&& v)                    { return add(std::move(v)); };

        // add matrix to the list of matrices
        dense_col&              add(const dense_matrix<T>&);
        dense_col&              add(dense_matrix<T>&&);
        dense_col&              operator,(const dense_matrix<T>& v) { return add(v); };
        dense_col&              operator,(dense_matrix<T>&& v)      { return add(std::move(v)); };

        // add concatenation helpers to the list of matrices
        dense_col&              add(const dense_row<T>&);
        dense_col&              add(dense_row<T>&&);
        dense_col&              add(const dense_col&);
        dense_col&              add(dense_col&&);
        dense_col&              operator,(const dense_row<T>& v)    { return add(v); };
        dense_col&              operator,(dense_row<T>&& v)         { return add(std::move(v)); };
        dense_col&              operator,(const dense_col& v)       { return add(v); };
        dense_col&              operator,(dense_col&& v)            { return add(std::move(v)); };

        // convert to dense_matrix
        const dense_matrix<T>   to_matrix() const;
};

// dense matrix represented in column oriented way; representation need not
// be continuous, distance from first element in columns c+1 and c is the leading
// dimension that can be larger than number of rows;
// dense_matrix bahaves according to copy-on-write scheme; dense_matrix is 
// thread-safe, but given instance cannot be shared between threads unless all
// accesses are read only.
template<class T>
class MATCL_MATREP_EXPORT dense_matrix<T, true>
{
    private:
        using mat_type          = raw::Matrix<T,struct_dense>;
        using sub_dense_matrix  = sub_dense_matrix<T>;

    public:
        using value_type    = T;

    private:
        Integer             m_rows;
        Integer             m_cols;
        Integer             m_max_rows;
        Integer             m_max_cols;
        Integer             m_ld;
        Integer             m_size;
        value_type*         m_ptr;
        struct_flag*        m_flag;
        matcl::Matrix       m_matrix;

        friend dense_matrix<T,false>;
        friend matcl::sub_dense_matrix<T>;

    public:
        //--------------------------------------------------------------------
        //      constructors, destructor, assignments
        //--------------------------------------------------------------------

        // default constructor create 1x1 matrix with zero scalar
        dense_matrix();

        // conversion from bool to dense_matrix is explicit
        explicit dense_matrix(bool val);

        // conversion from scalar to dense_matrix is explicit
        explicit dense_matrix(const T& val);

        // convert general matrix to dense_matrix; if mat is not dense, then
        // conversions take place
        explicit dense_matrix(const matcl::Matrix& mat);
        explicit dense_matrix(matcl::Matrix&& mat);

        // copy constructor and move constructor
        dense_matrix(const dense_matrix& mat);
        dense_matrix(dense_matrix&& mat);

        // conversion from other typed matrices
        explicit dense_matrix(const sparse_matrix<T>& mat);
        explicit dense_matrix(sparse_matrix<T>&& mat);
        explicit dense_matrix(const band_matrix<T>& mat);
        explicit dense_matrix(band_matrix<T>&& mat);

        // conversion from submatrix, a new matrix is usually created
        dense_matrix(const sub_dense_matrix&);
        
        // conversion from dense_row (i.e. concatenation helper)
        dense_matrix(const dense_row<T>&);

        // conversion from dense_col (i.e. concatenation helper)
        dense_matrix(const dense_col<T>&);

        // destructor, memory is freed if refcount drops to zero
        ~dense_matrix();

        // assignments 
        dense_matrix&   operator=(const dense_matrix&) &;
        dense_matrix&   operator=(dense_matrix&&) &;        

        //--------------------------------------------------------------------
        //      functions specific to dense_matrix
        //--------------------------------------------------------------------
        
        // get pointer to the first elemenent stored the matrix; return nullptr
        // if matrix is empty; pointer is valid as long as this instance exists
        const value_type*   ptr() const;

        // get pointer to the first elemenent stored the matrix; return nullptr
        // if matrix is empty; safe version of dense_matrix allows for const access only
        const value_type*   ptr();

        // number of elements stored in the matrix (rows x columns);
        Integer             size() const;

        // maximum number of rows; if one calls resize(r, c) and r <= max_rows,
        // c <= max_cols, then no memory allocations are needed
        Integer             max_rows() const;

        // maximum number of rows; if one calls resize(r, c) and r <= max_rows,
        // c <= max_cols, then no memory allocations are needed
        Integer             max_cols() const;

        // leading dimension, i.e. distance between first elements in
        // columns c+1 and c
        Integer             ld() const;

        // convert to general matrix
        const Matrix&       to_matrix() const &;
        Matrix&&            to_matrix() &&;

        // get element at position pos; checks are not performed
        // safe version of dense_matrix allows for const access only
        const T&            elem(Integer pos) const;
        const T&            elem(Integer pos);

        // get element at given row r and column c; checks are not performed
        // safe version of dense_matrix allows for const access only
        const T&            elem(Integer r, Integer c) const;
        const T&            elem(Integer r, Integer c);

        //--------------------------------------------------------------------
        //          const functions defined for all matrix types
        //--------------------------------------------------------------------
        
        // number of rows
        Integer             rows() const;

        // number of columns
        Integer             cols() const;
        
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
        const struct_flag   get_struct() const;
        
        // struct flag must be valid, data are not changed
        void                set_struct(const struct_flag&) const;
        
        // struct flag must be valid, data are not changed
        void                add_struct(const struct_flag&) const;

        // return type_info of stored elements
        ti::ti_object       get_type() const;

        // conversion to bool, throw error if Matrix is not convertible to scalar and
        // return true for nonzero scalars
        explicit            operator bool() const;

        // true if this is 0xn or mx0 matrix
        bool                is_empty() const;

        // true is this is a scalar or 1x1 matrix
        bool                is_scalar() const;

        // true if rows() == cols()
        bool                is_square() const;

        // true if rows() == 1 or cols() == 1 or is_empty() == true
        bool                is_vector() const;

        // true if given matrix is represented as a matrix
        bool                is_matrix_type() const;

        // true if given matrix is represented as scalar
        bool                is_scalar_type() const;

        // true if reference count is 1
        bool                is_unique() const;

        // code of stored elements type
        value_code          get_value_code() const;

        // code if matrix representation (scalar, band, dense, or sparse)
        struct_code         get_struct_code() const;

        // matrix code (value code and struct code)
        mat_code            get_matrix_code() const;

        // delete rows specified by colon c.
        const dense_matrix  delrows(const colon&) const &;
        const dense_matrix  delrows(const colon&) const &&;

        // delete columns specified in colon c.
        const dense_matrix  delcols(const colon&) const &;
        const dense_matrix  delcols(const colon&) const &&;

        // delete rows specified by colon c1 and columns specified by colon c2.
        const dense_matrix  delrowscols(const colon& c1, const colon& c2) const &;
        const dense_matrix  delrowscols(const colon& c1, const colon& c2) const &&;

        // get element at position pos; checks are performed
        inline const T&     operator()(Integer pos) const;
        
        // get element at given row r and column c; checks are performed
        inline const T&     operator()(Integer r, Integer c) const;
        
        // get submatrix, a new matrix is usually created, unless matrix views 
        // can be created
        const dense_matrix  operator()(const colon& r) const;        
        const dense_matrix  operator()(const colon& r, const colon& c) const;
        
        // get diagonal
        const dense_matrix  diag(Integer d = 0) const;

        // if matrix is unique (i.e. refcount is 1) then do nothing; otherwise
        // make copy of this matrix, note that refcount of this instance will
        // also drop to one
        const dense_matrix& make_unique() const;

        // make independent copy
        const dense_matrix  clone() const;

        //--------------------------------------------------------------------
        //          non const functions defined for all matrix types
        //--------------------------------------------------------------------
        // get element, safe version of dense_matrix allows for const access only
        inline const T&     operator()(Integer r);
        inline const T&     operator()(Integer r, Integer c);

        // get or change submatrix
        sub_dense_matrix    operator()(const colon& r);
        sub_dense_matrix    operator()(const colon& r, const colon& c);

        // get or change diagonal
        sub_dense_matrix    diag(Integer d = 0);

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
        // create dense matrix with zeroes of size rxc
        template<class Enable = typename details::enable_if_not_object<T, void>::type>
        static dense_matrix zeros(Integer r, Integer c);

        // create dense matrix with zeroes of size rxc, version for T == Object requires
        // additionally type info of stored elements
        template<class Enable = typename details::enable_if_object<T, void>::type>
        static dense_matrix zeros(ti::ti_object, Integer r, Integer c);

        // create dense matrix with ones of size rxc
        template<class Enable = typename details::enable_if_not_object<T, void>::type>
        static dense_matrix ones(Integer r, Integer c);

        // create dense matrix with ones of size rxc, version for T == Object requires
        // additionally type info of stored elements
        template<class Enable = typename details::enable_if_object<T, void>::type>
        static dense_matrix ones(ti::ti_object, Integer r, Integer c);

        // create dense matrix with ones on diagonal of size rxc
        template<class Enable = typename details::enable_if_not_object<T, void>::type>
        static dense_matrix eye(Integer r, Integer c);

        // create dense matrix with ones on diagonal of size rxr
        template<class Enable = typename details::enable_if_not_object<T, void>::type>
        static dense_matrix eye(Integer r);

        // create dense matrix with ones on diagonal of size rxc, version for T == Object
        // requires additionally type info of stored elements
        template<class Enable = typename details::enable_if_object<T, void>::type>
        static dense_matrix eye(ti::ti_object, Integer r, Integer c);

        // create dense matrix with ones on diagonal of size rxr, version for T == Object
        // requires additionally type info of stored elements
        template<class Enable = typename details::enable_if_object<T, void>::type>
        static dense_matrix eye(ti::ti_object, Integer r);

        // create dense matrix with elements on diagonal d given by v of size kxk, where 
        // k = length(v) + abs(d)
        static dense_matrix diag(const dense_matrix& v, Integer d = 0);

        // create rxc dense matrix from the columns of A and places them along the diagonals 
        // specified by d 
        static dense_matrix diags(const dense_matrix& A, const Matrix &d, Integer r, Integer c);

        // create dense matrix of size rxc of uniformly distributed random numbers on [0,1]
        static dense_matrix rand(Integer r, Integer c, const rand_state& rand_ptr 
                                = global_rand_state());

        // create dense matrix of size rxc of normally distributed random numbers
        template<class Enable = typename details::enable_if_float<T, void>::type>
        static dense_matrix randn(Integer r, Integer c, const rand_state& rand_ptr = global_rand_state());

        // create sequence of numbers from s to e with step 1.0
        template<class Enable = typename details::enable_if_real<T, void>::type>
        static dense_matrix range(T s, T e);

        // create sequence of numbers from s to e with step i
        template<class Enable = typename details::enable_if_real<T, void>::type>
        static dense_matrix range(T s, T i, T e);

        // create sequence of length n of equally spaced numbers from s to e
        template<class Enable = typename details::enable_if_real_float<T, void>::type>
        static dense_matrix linspace(typename details::real_type<T>::type s, 
                                typename details::real_type<T>::type e, Integer n);

        // create sequence of length n of logarithmically equally spaced numbers 
        // from 10^s to 10^e
        template<class Enable = typename details::enable_if_real_float<T, void>::type>
        static dense_matrix logspace(typename details::real_type<T>::type s, 
                                typename details::real_type<T>::type e, Integer n);

        // create dense matrix of size rows x cols with zeroes
        template<class Enable = typename details::enable_if_not_object<T, void>::type>
        static dense_matrix make(Integer rows,Integer cols);

        // create dense matrix of size rows x cols with zeroes
        // version for T == Object requires additionally type info of stored elements
        template<class Enable = typename details::enable_if_object<T, void>::type>
        static dense_matrix make(ti::ti_object ti, Integer rows,Integer cols);

        // create dense matrix of size rows x cols with elements copied from array
        // arr (column oriented storage assumed); optional argument set leading dimension 
        // (distance between consecutive columns)
        template<class Enable = typename details::enable_if_not_object<T, void>::type>
        static dense_matrix make(Integer rows,Integer cols, const T *arr);
        template<class Enable = typename details::enable_if_not_object<T, void>::type>
        static dense_matrix make(Integer rows,Integer cols, const T *arr, Integer ld);

        // create dense matrix of size rows x cols with elements copied from array
        // arr (column oriented storage assumed); version for T == Object requires
        // additionally type info of stored elements; optional argument set leading dimension 
        // (distance between consecutive columns)
        template<class Enable = typename details::enable_if_object<T, void>::type>
        static dense_matrix make(ti::ti_object ti, Integer rows,Integer cols, const T *arr);
        template<class Enable = typename details::enable_if_object<T, void>::type>
        static dense_matrix make(ti::ti_object ti, Integer rows,Integer cols, const T *arr, 
                                Integer ld);

        // create dense matrix of size rows x cols with elements stored in the array arr 
        // (column oriented storage assumed) with leading dimension ld (distance between
        // consecutive columns) array arr is not internally copied;
        //
        // this is unsafe function; one should ensure that this matrix and all matrices 
        // created from this matrix directly or indirectly are destroyed before destroying
        // the array arr
        template<class Enable = typename details::enable_if_not_object<T, void>::type>
        static dense_matrix make_foreign(Integer rows,Integer cols, T* arr, Integer ld);
        template<class Enable = typename details::enable_if_object<T, void>::type>
        static dense_matrix make_foreign(ti::ti_object ti, Integer rows,Integer cols, Object* arr, 
                                        Integer ld);

        // create dense matrix of size rows x cols with elements copied from Real arrays
        // ar_r, ar_i containing real and complex part (column oriented storage assumed)
        template<class Enable = typename details::enable_if_complex<T, void>::type>
        static dense_matrix make_complex(Integer rows,Integer cols, 
                                const typename details::real_type<T>::type *ar_r,
                                const typename details::real_type<T>::type* ar_i);

        // create dense matrix of size rows x cols with elements val
        static dense_matrix make(const T& val, Integer rows, Integer cols);

        // create dense matrix of size rows x cols; elements are not initialized; 
        // on return ptr_data is the pointer to stored data
        template<class Enable = typename details::enable_if_not_object<T, void>::type>
        static dense_matrix make_noinit(Integer rows,Integer cols, T*& ptr_data);

        // create dense matrix of size rows x cols; elements are not initialized; 
        // on return ptr_data is the pointer to stored data; version for T == Object
        // requires additionally type info of stored elements
        template<class Enable = typename details::enable_if_object<T, void>::type>
        static dense_matrix make_noinit(ti::ti_object ti, Integer rows,Integer cols, T*& ptr_data);

    private:
        struct str_make_unique{};
        dense_matrix(matcl::Matrix& mat, str_make_unique);
        dense_matrix(matcl::Matrix&& mat, str_make_unique);

        void                throw_error_single_index(Integer r, Integer size) const;
        void                throw_error_double_index(Integer r, Integer c, Integer rows, Integer cols) const;

        void                init_from_rep(const mat_type& mat);
        void                update_rep();
};

// unsafe extension of dense_matrix that allows for modification of 
// individual elements; this class is unsafe because additional requirements
// must be fulfilled, that guaratees that non const methods are called
// on unique matrices, which cannot be assured at compile time;
// given instance cannot be shared between threads unless all accesses are read only
template<class T>
class MATCL_MATREP_EXPORT dense_matrix<T,false> : public dense_matrix<T,true>
{
    private:
        using base_type     = dense_matrix<T,true>;

    public:
        using value_type    = T;

    public:
        //--------------------------------------------------------------------
        //      constructors, destructor, assignments
        //--------------------------------------------------------------------

        // default constructor create 1x1 matrix with zero scalar
        dense_matrix();

        // convert general matrix to dense_matrix; if mat is not dense, then
        // conversions take place; this make matrix 'mat' unique, therefore copying
        // can happen; one cannot call any non const functions on 'mat' until this object
        // is destroyed
        explicit dense_matrix(matcl::Matrix& mat);
        explicit dense_matrix(matcl::Matrix&& mat);

        // non safe dense_matrix cannot be copied, only moved
        dense_matrix(const dense_matrix<T,false>& mat) = delete;
        dense_matrix(dense_matrix<T,false>&& mat);

        // conversion from safe dense_matrix; make matrix 'mat' unique, therefore
        // copying can happen; one cannot call any non const functions on 'mat' until
        // this object is destroyed
        explicit dense_matrix(dense_matrix<T,true>&& mat);
        explicit dense_matrix(dense_matrix<T,true>& mat);

        // non safe dense_matrix cannot be assigned, only move assignment is allowed        
        dense_matrix&       operator=(const dense_matrix&) & = delete;
        dense_matrix&       operator=(dense_matrix&&) &;

        // destructor, memory is freed if refcount drops to zero
        ~dense_matrix();

        //--------------------------------------------------------------------
        //      functions specific to dense_matrix
        //--------------------------------------------------------------------

        // get pointer to the first elemenent stored the matrix; return nullptr
        // if matrix is empty
        value_type*         ptr();

        // get or modify element at position pos; checks are not performed
        T&                  elem(Integer pos);
        
        // get of modify element at given row r and column c; checks are not performed
        T&                  elem(Integer r, Integer c);

        //--------------------------------------------------------------------
        //          non const functions defined for all matrix types
        //--------------------------------------------------------------------
        // get or change submatrix
        inline T&           operator()(Integer r);
        inline T&           operator()(Integer r, Integer c);
};

};

#include "matcl-matrep/details/matrix_rep_dense.inl"
